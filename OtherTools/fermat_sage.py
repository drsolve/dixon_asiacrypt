#!/usr/bin/env python3
"""Fermat interpolation benchmarks and summaries in one script.

Default: benchmark cases listed in time.txt (existing inputs).
--extend: generate 3x13, 3x14, ...; stop after a completed case exceeds 1800s,
          or on timeout/error. Per-case timeout defaults to 3600s.
--summary-only: merge saved JSON results into time_interp.txt; never runs Fermat.

All computation times exclude input generation. Original inputs/time.txt are kept.
"""
import argparse
import csv
from itertools import product
import json
from pathlib import Path
import re
import random
import resource
import shutil
import subprocess
import tempfile
import time

ROOT = Path(__file__).resolve().parent
FIELDS = ['case', 'status', 'wall_seconds', 'cpu_seconds', 'det_seconds',
          'matrix_size', 'result_degree', 'old_status', 'old_seconds',
          'speedup', 'returncode', 'log']


def make_input(source):
    tail = '&(D=-1);\n!!(Fill);\ns := Det[m4];\n&x;'
    if source.count(tail) != 1:
        raise ValueError('Expected exactly one final determinant driver')
    return source.replace(tail, '''&(D=-1);
&(L=1);
!!(Fill);
!!('BENCH_MATRIX ',Cols[m4]);
benchstart := &T;
s := Det[m4];
!!('BENCH_DET_MS ', &T-benchstart);
!!('BENCH_DEGREE ',Deg(s,zq0));
!!'BENCH_DONE';
&x;''')


def save_results(out, results):
    # Atomic snapshots after every case, so an interruption retains finished cases.
    for name, content in [('results.json', json.dumps(results, indent=2) + '\n'),
                          ('time.interp.txt', ''.join(
                              f"Running {r['case']}...\n  {r['status']}, time={r['wall_seconds']:.6f}s\n"
                              for r in results))]:
        tmp = out / (name + '.tmp')
        tmp.write_text(content)
        tmp.replace(out / name)
    tmp = out / 'comparison.csv.tmp'
    with tmp.open('w', newline='') as f:
        writer = csv.DictWriter(f, FIELDS)
        writer.writeheader()
        writer.writerows(results)
    tmp.replace(out / 'comparison.csv')


def run_case(executable, runtime, source, case, out, timeout,
             old_status="unavailable", old_seconds=None):
    casefile = runtime / 'case.fer'
    casefile.write_text(source)
    (out / f'{case}.interp.fer').write_text(source)
    logfile = out / f'{case}.log'
    print(f'Running {case}...', flush=True)
    cpu_start = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.perf_counter()
    status = 'ok'
    with logfile.open('wb') as log:
        proc = subprocess.Popen([str(executable)], cwd=ROOT, stdin=subprocess.PIPE,
                                stdout=log, stderr=subprocess.STDOUT)
        try:
            proc.communicate(b"&(R = '" + str(casefile).encode() + b"');\n&q\ny\n", timeout=timeout)
        except subprocess.TimeoutExpired:
            status = 'timeout'
            proc.kill()
            proc.communicate()
        except BaseException:
            proc.kill()
            proc.communicate()
            raise
    wall = time.perf_counter() - start
    cpu_end = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = cpu_end.ru_utime + cpu_end.ru_stime - cpu_start.ru_utime - cpu_start.ru_stime
    text = logfile.read_text(errors='replace')
    if status == 'ok' and (proc.returncode != 0 or 'BENCH_DONE' not in text
                           or re.search(r'\*\*\* Fermat error', text, re.I)):
        status = 'error'
    def marker(name):
        m = re.search(r'BENCH_' + name + r'\s+(-?\d+)', text)
        return int(m[1]) if m else None
    det_ms = marker('DET_MS')
    result = dict(case=case, status=status, wall_seconds=wall, cpu_seconds=cpu,
                  det_seconds=det_ms / 1000 if det_ms is not None else None,
                  matrix_size=marker('MATRIX'), result_degree=marker('DEGREE'),
                  old_status=old_status, old_seconds=float(old_seconds) if old_seconds is not None else None,
                  speedup=float(old_seconds) / wall if old_status == status == 'ok' else None,
                  returncode=proc.returncode, log=str(logfile))
    det_display = f"{result['det_seconds']:.3f}s" if det_ms is not None else "unfinished"
    print(f"  {status}, time={wall:.6f}s, det={det_display}", flush=True)
    return result


def generate_case(degree, seed=20260828, prime=536870923):
    rng = random.Random(seed)
    monomials = []
    for exponent in product(range(degree + 1), repeat=3):
        if sum(exponent) > degree:
            continue
        factors = [f'zq{i}' if power == 1 else f'zq{i}^{power}'
                   for i, power in enumerate(exponent) if power]
        monomials.append('*'.join(factors) or '1')
    equations = [' + '.join(f'{rng.randrange(1, prime)}*{m}' for m in monomials)
                 for _ in range(3)]
    template = (ROOT / 'dense_3x12_seed20260828.fer').read_text()
    tail = template[template.index(';; Dixon routines'):]
    return (f'&(p = {prime});\n&(J = zq0);\n&(J = zq1);\n&(J = zq2);\n'
            '@([d]);\nArray d[3];\n[d] := [[' + ',\n          '.join(equations)
            + ']];\n\n' + tail)


def run_batch(args, parser):
    entries = re.findall(r'Running (\d+x\d+)\.\.\.\s+(\w+), time=([\d.]+)s',
                         args.input.read_text())
    if not entries:
        parser.error('No cases found in timing input')
    if args.cases:
        unknown = set(args.cases) - {x[0] for x in entries}
        if unknown:
            parser.error(f'Unknown cases: {sorted(unknown)}')
        entries = [x for x in entries if x[0] in args.cases]
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    results_path = out / 'results.json'
    if results_path.exists() and not args.resume:
        parser.error('Results already exist; use --resume or another --output')
    results = json.loads(results_path.read_text()) if args.resume and results_path.exists() else []
    done = {r['case'] for r in results}
    # Validate all inputs before starting potentially long computations.
    inputs = {}
    for case, _, _ in entries:
        inputs[case] = make_input((ROOT / f'dense_{case}_seed20260828.fer').read_text())
    with tempfile.TemporaryDirectory(prefix='fermat-interp-') as tmp:
        runtime = Path(tmp)
        executable = runtime / 'fer64'
        shutil.copy2(args.fermat, executable)
        executable.chmod(executable.stat().st_mode | 0o100)
        (runtime / 'BACKWARD').symlink_to(args.fermat.resolve().parent / 'BACKWARD', target_is_directory=True)
        for case, old_status, old_seconds in entries:
            if case in done:
                continue
            result = run_case(executable, runtime, inputs[case], case, out, args.timeout,
                              old_status, old_seconds)
            results.append(result)
            save_results(out, results)
    print(f'Results: {out / "comparison.csv"}', flush=True)


def run_extend(args, parser):
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    config = dict(start_degree=args.start_degree, stop_seconds=args.stop_seconds,
                  timeout=args.timeout, seed=20260828, prime=536870923,
                  method='dense D=-1 L=1', timing='whole Fermat process, excludes input generation')
    config_path = out / 'config.json'
    if config_path.exists():
        if not args.resume:
            parser.error('Output already exists; use --resume or another --output')
        if json.loads(config_path.read_text()) != config:
            parser.error('Resume settings differ from saved config')
    else:
        config_path.write_text(json.dumps(config, indent=2)+'\n')
    result_path = out / 'results.json'
    results = json.loads(result_path.read_text()) if result_path.exists() else []
    def stop_reason(r):
        if r['status'] != 'ok':
            return r['status']
        if r['wall_seconds'] > args.stop_seconds:
            return 'wall_seconds_exceeded'
        return None
    if results and stop_reason(results[-1]):
        print(f"Already stopped: {results[-1]['case']}, {stop_reason(results[-1])}")
        return
    expected = [f'3x{d}' for d in range(args.start_degree, args.start_degree + len(results))]
    if [r['case'] for r in results] != expected:
        parser.error('Saved results are not a contiguous sequence from start-degree')
    with tempfile.TemporaryDirectory(prefix='fermat-extend-') as tmp:
        runtime = Path(tmp)
        exe = runtime / 'fer64'
        shutil.copy2(args.fermat, exe)
        exe.chmod(exe.stat().st_mode | 0o100)
        (runtime / 'BACKWARD').symlink_to(args.fermat.resolve().parent / 'BACKWARD', target_is_directory=True)
        degree = args.start_degree + len(results)
        while True:
            case = f'3x{degree}'
            source = generate_case(degree)
            (out / f'dense_{case}_seed20260828.fer').write_text(source)
            result = run_case(exe, runtime, make_input(source), case, out, args.timeout)
            results.append(result)
            save_results(out, results)
            reason = stop_reason(result)
            if reason:
                (out / 'stop.json').write_text(json.dumps(dict(
                    case=case, reason=reason, wall_seconds=result['wall_seconds'],
                    threshold_seconds=args.stop_seconds, timeout_seconds=args.timeout), indent=2)+'\n')
                print(f'Stopped at {case}: {reason}', flush=True)
                break
            degree += 1
    print(f'Results: {out / "time.interp.txt"}')


def summarize_results(paths, output):
    """Merge existing measurements, sorted by numeric n/degree; no computation."""
    rows = {}
    for path in paths:
        for row in json.loads(path.read_text()):
            case = row['case']
            if not re.fullmatch(r'\d+x\d+', case):
                raise ValueError(f'Invalid case name: {case}')
            if case in rows:
                raise ValueError(f'Duplicate measurement for {case}; choose one input')
            rows[case] = row
    ordered = sorted(rows.values(), key=lambda r: tuple(map(int, r['case'].split('x'))))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(''.join(
        f"Running {r['case']}...\n  {r['status']}, time={r['wall_seconds']:.6f}s\n"
        for r in ordered))
    print(f'Summary: {output} ({len(ordered)} cases; no computations run)')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--extend', action='store_true', help='Generate increasing 3xd cases')
    mode.add_argument('--summary-only', action='store_true', help='Only summarize saved results')
    parser.add_argument('--input', type=Path, default=ROOT / 'time.txt')
    parser.add_argument('--output', type=Path,
                        help='Result directory, or summary text file in --summary-only mode')
    parser.add_argument('--fermat', type=Path, default=ROOT / 'fer64')
    parser.add_argument('--timeout', type=float, default=3600)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--cases', nargs='+', help='Select cases from --input in default mode')
    parser.add_argument('--start-degree', type=int, default=13)
    parser.add_argument('--stop-seconds', type=float, default=1800)
    parser.add_argument('--results', type=Path, nargs='+',
                        default=[ROOT / 'out/interp_benchmark/results.json',
                                 ROOT / 'out/interp_3x_extend/results.json'],
                        help='Saved JSON files for --summary-only')
    args = parser.parse_args()
    if args.summary_only:
        summarize_results(args.results, args.output or ROOT / 'time_interp.txt')
        return
    if args.timeout <= 0 or args.start_degree < 1 or args.stop_seconds <= 0:
        parser.error('Degree and time limits must be positive')
    if args.extend and args.cases:
        parser.error('--cases cannot be combined with --extend')
    if args.output is None:
        args.output = ROOT / ('out/interp_3x_extend' if args.extend else 'out/interp_benchmark')
    if args.extend:
        run_extend(args, parser)
    else:
        run_batch(args, parser)


if __name__ == '__main__':
    main()
