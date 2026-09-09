"""Shared Sage runtime copied into each exported XHash test suite."""

import argparse
import hashlib
import json
import os
import re
from pathlib import Path
import shutil
import time

from sage.all import GF, PolynomialRing


def digest(data):
    payload = {k: data[k] for k in ("p", "variables", "equations", "relations", "free_inputs")}
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def read_case(path):
    data = json.loads(Path(path).read_text())
    assert digest(data) == data["sha256"], "Instance checksum mismatch"
    ring = PolynomialRing(GF(data["p"]), names=data["variables"], order="lex")
    polys = [ring(s) for s in data["equations"]]
    relations = [(ring(a), ring(b)) for a, b in data["relations"]]
    u = data["free_inputs"]
    assert len(polys) == len(data["variables"])
    assert len(relations) == len(polys) - u
    for i, (lhs, rhs) in enumerate(relations):
        assert polys[i] in (lhs - rhs, rhs - lhs), "Ideal relation differs from input equation"
        assert lhs.variables() == (ring.gen(u + i),)
        assert all(ring.gens().index(v) < u + i for v in rhs.variables())
    return data, ring, polys, relations


def run(case, method, suite):
    suite = Path(suite).resolve()
    parser = argparse.ArgumentParser(description=f"{case}: {method} elimination benchmark")
    parser.add_argument("--check-only", action="store_true", help="Validate input without running a solver")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--timeout", type=float, default=3600, help="Total solver-call budget in seconds")
    parser.add_argument("--drsolve", type=Path, default=None)
    args = parser.parse_args()
    if args.threads < 1 or args.timeout <= 0:
        parser.error("threads and timeout must be positive")
    data, ring, polys, relations = read_case(suite / "data" / (case + ".json"))
    print(json.dumps({"case": case, "method": method, "p": data["p"],
                      "equations": len(polys), "variables": len(ring.gens()),
                      "degrees": [int(f.total_degree()) for f in polys],
                      "sha256": data["sha256"]}, indent=2), flush=True)
    if args.check_only:
        return
    binary = args.drsolve or (Path(os.environ["DRSOLVE"]) if "DRSOLVE" in os.environ else None)
    if binary is None:
        binary = next((p / "drsolve" / "drsolve" for p in suite.parents
                       if (p / "drsolve" / "drsolve").is_file()), None)
    if binary is None or not binary.is_file():
        raise FileNotFoundError("Use --drsolve /absolute/path/to/drsolve or set DRSOLVE")
    binary = binary.resolve()
    # Some source archives lose executable mode; preserve the original binary.
    if not os.access(binary, os.X_OK):
        cached = suite / "results" / ".bin" / "drsolve"
        cached.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(binary, cached)
        cached.chmod(cached.stat().st_mode | 0o100)
        binary = cached
    namespace = {"__name__": "drsolve_interface"}
    interface = suite / "_support" / "drsolve_sage_interface.sage"
    exec(compile(interface.read_text(), str(interface), "exec"), namespace)
    out = suite / "results" / f"{case}_{method}_{time.time_ns()}"
    out.mkdir(parents=True)
    start = time.monotonic()
    record = {"case": case, "method": method, "sha256": data["sha256"],
              "threads": args.threads, "timeout_seconds": args.timeout,
              "binary": str(binary), "stages": [], "scope": "elimination workflow; no complete back-substitution"}
    original_cwd = Path.cwd()
    os.chdir(out)

    def call(label, equations, eliminate, ideal=False, subres=False):
        remaining = args.timeout - (time.monotonic() - start)
        if remaining <= 0:
            raise TimeoutError("Total elimination budget exhausted")
        kwargs = dict(field_size=data["p"], dixon_path=str(binary), threads=args.threads,
                      timeout=remaining, verbosity=2, time=True, live_output=True,
                      finput=str(out / (label + ".in")), foutput=str(out / (label + ".out")))
        if subres:
            kwargs["resultant_method"] = "subres"
        before = time.monotonic()
        print(f"\n=== {label} ===", flush=True)
        if ideal:
            answer = namespace["DixonIdeal"](equations, relations, eliminate, **kwargs)
        else:
            answer = namespace["DixonRes"](equations, eliminate, **kwargs)
        if answer is None:
            raise RuntimeError(f"Solver failed at {label}; inspect the saved input/output")
        # The bundled interface also captures the optional univariate root report.
        if isinstance(answer, str):
            answer = re.split(r"(?m)^Roots in ", answer, maxsplit=1)[0].strip()
        polynomial = ring(answer)
        if polynomial == 0:
            raise RuntimeError(f"Zero eliminant at {label}; further elimination is not valid")
        if any(v in polynomial.variables() for v in eliminate):
            raise RuntimeError(f"Eliminated variable remains after {label}")
        record["stages"].append({"stage": label, "seconds": time.monotonic() - before,
                                 "degree": int(polynomial.total_degree()),
                                 "terms": len(polynomial.monomials())})
        return polynomial

    try:
        if method == "direct":
            result = call("direct", polys, list(ring.gens()[1:]))
        elif method == "hybrid":
            reduced = []
            for j, constraint in enumerate(polys[len(relations):]):
                current = constraint
                for i in range(len(relations) - 1, -1, -1):
                    variable = ring.gen(data["free_inputs"] + i)
                    # An absent variable contributes no useful elimination step.
                    if variable not in current.variables():
                        continue
                    current = call(f"output{j}_level{i}", [current, polys[i]], [variable], ideal=True)
                reduced.append(current)
            if data["free_inputs"] == 1:
                result = reduced[0]
            else:
                result = call("final_merge", reduced, list(ring.gens()[1:data["free_inputs"]]),
                              subres=data["free_inputs"] == 2)
        else:
            raise ValueError(method)
        (out / "eliminant.txt").write_text(str(result) + "\n")
        record["status"] = "ok"
    except Exception as error:
        record.update(status="failed", error=f"{type(error).__name__}: {error}")
        raise
    finally:
        record["elapsed_seconds"] = time.monotonic() - start
        (out / "summary.json").write_text(json.dumps(record, indent=2) + "\n")
        os.chdir(original_cwd)
        print(f"\nResults: {out}\nElapsed: {record['elapsed_seconds']:.3f} seconds", flush=True)
