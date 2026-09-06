#!/usr/bin/env python3
"""Export the paper's 14 numerical plots as vector PDFs.

Run with Python 3 + matplotlib + numpy from any directory. Existing notebook
models and recorded benchmark timings are reused without rerunning experiments.
Only Figure 6(d) changes its numerical model: max(best Step 1, best Step 4),
as already specified in the paper's caption. No other cost model is revised.
"""

import argparse
import ast
import csv
import io
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
TEX = HERE.parent / "tex"


def cell_source(name, index):
    notebook = json.loads((HERE / name).read_text(encoding="utf-8"))
    return "".join(notebook["cells"][index]["source"])


def load_cell(name, index, namespace=None, stop_at=None):
    """Load plotting definitions without running notebook demo/main blocks."""
    namespace = namespace if namespace is not None else {"__name__": "figure_export"}
    tree = ast.parse(cell_source(name, index))
    body = []
    for node in tree.body:
        if stop_at and isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == stop_at for t in node.targets):
            break
        if isinstance(node, ast.If) and "__name__" in ast.unparse(node.test):
            continue
        body.append(node)
    tree.body = body
    exec(compile(tree, f"{name}:cell{index}", "exec"), namespace)
    plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42})
    return namespace


def stage_models():
    return load_cell("comp_step1&4.ipynb", 0, stop_at="rows_d")


def export_fig6d(output_dir=TEX):
    """Single authoritative entry point for the corrected omega=2.37 panel."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    with plt.rc_context():
        ns = load_cell("Complexity_Comparision.ipynb", 0)
        stages = stage_models()
        rows = [stages["equal_degree_report"](m=4, d=d, s=1, omega=2.37, q=257)
                for d in range(2, 20)]
        lookup = {r["d"]: r["overall"] for r in rows}
        original_models = list(ns["MODELS"])
        label, _, color = ns["MODELS"][0]
        ns["MODELS"][0] = (label, lambda n, d, omega: lookup[d], color)
        ns["plot_vs_d"](output_dir / "complexity_vs_degree_n5w237.pdf", n=5, omega=2.37)
        columns = ["n", "d", "omega", "step1_best_method", "step1_best_log2",
                   "step4_best_method", "step4_best_log2", "dixon_max_log2",
                   "previous_step4_log2", "increase_log2"]
        columns += [label for label, _, _ in original_models[1:]]
        with (HERE / "fig6d_recomputed.csv").open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(columns)
            for r in rows:
                old = original_models[0][1](5, r["d"], 2.37)
                writer.writerow([5, r["d"], 2.37, r["step1_best_method"], r["step1_best"],
                                 r["step4_best_method"], r["step4"], r["overall"],
                                 old, r["overall"] - old] +
                                [model(5, r["d"], 2.37) for _, model, _ in original_models[1:]])
    return rows


def export_complexity():
    with plt.rc_context():
        ns = load_cell("Complexity_Comparision.ipynb", 0)
        ns["plot_vs_n"](TEX / "complexity_comparison_d2.pdf", d=2, omega=2.81)
        ns["plot_vs_n"](TEX / "complexity_comparison_d5.pdf", d=5, omega=2.81)
        ns["plot_vs_d"](TEX / "complexity_vs_degree_n3.pdf", n=3, omega=2.81)
        ns["plot_vs_d"](TEX / "complexity_vs_degree_n5.pdf", n=5, omega=2.81)
    export_fig6d()
    with plt.rc_context():
        ns = stage_models()
        rows_d = ns["degree_sweep_data"](m=2, d_min=2, d_max=50, s=1, q=257)
        rows_m = ns["variable_sweep_data"](d=3, m_min=2, m_max=20, s=1, q=257)
        for suffix, rows in [("degree", rows_d), ("variables", rows_m)]:
            fig, _ = ns["plot_step1_step4_vs_" + suffix](rows, show=False)
            fig.savefig(TEX / ("plot_step1_step4_vs_" + suffix + ".pdf"))
            plt.close(fig)


def export_ao():
    with plt.rc_context():
        ns = None
        for i in range(3):
            ns = load_cell("AO_Complexity.ipynb", i, ns)
        ns["plot_poseidon"](TEX / "poseidon_complexity.pdf", omega=2.81)
        for i, stem in [(3, "vision_complexity"), (4, "xhash_complexity")]:
            load_cell("AO_Complexity.ipynb", i, ns)
            with (TEX / (stem + ".pdf")).open("wb") as handle:
                plt.gcf().savefig(handle, format="pdf")
            plt.close()


def recorded_fermat_log():
    """Recover the exact printed timings already present in the uploaded file."""
    nb = json.loads((HERE / "fermat_sage.ipynb").read_text(encoding="utf-8"))
    text = "".join("".join(o.get("text", [])) for o in nb["cells"][1].get("outputs", []))
    if "Running 3x2..." not in text or "Running 7x2..." not in text:
        raise RuntimeError("Recorded Fermat benchmark output is missing; refusing to omit its series.")
    path = HERE / "fermat_timings_recorded.txt"
    path.write_text(text, encoding="utf-8")
    return path


def export_benchmarks():
    with plt.rc_context():
        ns = load_cell("benchmark_DixonMagmaMsolve.ipynb", 1)
        fermat_file = recorded_fermat_log()
        loader = ns["load_fermat_times"]
        ns["load_fermat_times"] = lambda: loader(fermat_file)
        data = json.loads((HERE / "benchmark_results_extended_20260514_013505.json").read_text())
        targets = {
            "standard1_n3_comparison": "GF536870923_n3_comparison",
            "standard1_n4_comparison": "GF536870923_n4_comparison",
            "standard1_n5_comparison": "GF536870923_n5_comparison",
            "standard1_d2_vary_n": "GF536870923_d2_vary_n",
        }
        def save_vector(fig, filename):
            stem = Path(filename).stem
            if stem in targets:
                fig.tight_layout()
                buffer = io.BytesIO()
                fig.savefig(buffer, format="pdf", bbox_inches="tight")
                (TEX / (targets[stem] + ".pdf")).write_bytes(buffer.getvalue())
            plt.close(fig)
        ns["save_figure"] = save_vector
        ns["plot_field_comparison"](data, "standard1", "GF(536870923)", output_dir=TEX)
        # Keep every n in this field, including n=6,7,... used by the d=2 plot.
        standard_only = dict(data)
        standard_only["results"] = {"standard1": data["results"]["standard1"]}
        ns["plot_fixed_degree_comparison"](standard_only, degree=2, output_dir=TEX)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fig6d-only", action="store_true")
    args = parser.parse_args()
    TEX.mkdir(exist_ok=True)
    if args.fig6d_only:
        export_fig6d()
    else:
        export_complexity()
        export_ao()
        export_benchmarks()
        print("Exported 14 vector PDFs; experiments were not rerun.")


if __name__ == "__main__":
    main()
