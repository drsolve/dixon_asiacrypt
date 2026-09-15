#!/usr/bin/env python3
"""Export the paper's 14 numerical plots as vector PDFs.

Run with Python 3 + matplotlib + numpy from any directory. Existing notebook
models and recorded benchmark timings are reused without rerunning experiments.
Figure 6(d) and XHash use the stage models specified in the paper.
"""

import argparse
import ast
import io
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from OtherTools.dixon_stage_models import (
        degree_sweep_data,
        plot_step1_step4_vs_degree,
        plot_step1_step4_vs_variables,
        variable_sweep_data,
    )
except ModuleNotFoundError:  # Direct execution from inside OtherTools.
    from dixon_stage_models import (
        degree_sweep_data,
        plot_step1_step4_vs_degree,
        plot_step1_step4_vs_variables,
        variable_sweep_data,
    )

HERE = Path(__file__).resolve().parent
TEX = HERE.parent / "figures"


def cell_source(name, index):
    notebook = json.loads((HERE / name).read_text(encoding="utf-8"))
    return "".join(notebook["cells"][index]["source"])


def load_cell(name, index, namespace=None):
    """Load plotting definitions without running notebook demo/main blocks."""
    namespace = namespace if namespace is not None else {"__name__": "figure_export"}
    tree = ast.parse(cell_source(name, index))
    body = []
    for node in tree.body:
        if isinstance(node, ast.If) and "__name__" in ast.unparse(node.test):
            continue
        body.append(node)
    tree.body = body
    exec(compile(tree, f"{name}:cell{index}", "exec"), namespace)
    plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42})
    return namespace


def export_fig6d(output_dir=TEX):
    """Export the Figure 6(d) implementation owned by the complexity notebook."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    with plt.rc_context():
        ns = load_cell("Complexity_Comparision.ipynb", 0)
        ns["plot_fig6d"](output_dir / "complexity_vs_degree_n5w237.pdf")


def export_complexity():
    with plt.rc_context():
        ns = load_cell("Complexity_Comparision.ipynb", 0)
        ns["plot_vs_n"](TEX / "complexity_comparison_d2.pdf", d=2, omega=2.81)
        ns["plot_vs_n"](TEX / "complexity_comparison_d5.pdf", d=5, omega=2.81)
        ns["plot_vs_d"](TEX / "complexity_vs_degree_n3.pdf", n=3, omega=2.81)
        ns["plot_vs_d"](TEX / "complexity_vs_degree_n5.pdf", n=5, omega=2.81)
    export_fig6d()
    with plt.rc_context():
        rows_d = degree_sweep_data(m=2, d_min=2, d_max=50, s=1, q=257)
        rows_m = variable_sweep_data(d=3, m_min=2, m_max=20, s=1, q=257)
        fig, _ = plot_step1_step4_vs_degree(rows_d, show=False)
        fig.savefig(TEX / "plot_step1_step4_vs_degree.pdf")
        plt.close(fig)
        fig, _ = plot_step1_step4_vs_variables(rows_m, show=False)
        fig.savefig(TEX / "plot_step1_step4_vs_variables.pdf")
        plt.close(fig)


def export_ao():
    with plt.rc_context():
        ns = None
        for i in range(3):
            ns = load_cell("AO_Complexity.ipynb", i, ns)
        ns["plot_poseidon"](TEX / "poseidon_complexity.pdf", omega=2.81)
        ns["plot_vision"](TEX / "vision_complexity.pdf")
        ns["plot_xhash"](TEX / "xhash_complexity.pdf")


def export_xhash():
    """Export the XHash figure from the AO notebook's authoritative model."""
    with plt.rc_context():
        ns = None
        for i in range(3):
            ns = load_cell("AO_Complexity.ipynb", i, ns)
        ns["plot_xhash"](TEX / "xhash_complexity.pdf")


def recorded_fermat_log():
    """Use the recorded timings directly, without overwriting benchmark data."""
    path = HERE / "fermat_timings_recorded.txt"
    text = path.read_text(encoding="utf-8")
    if "Running 3x2..." not in text or "Running 7x2..." not in text:
        raise RuntimeError("Recorded Fermat benchmark output is missing; refusing to omit its series.")
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
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--fig6d-only", action="store_true")
    group.add_argument("--benchmarks-only", action="store_true")
    group.add_argument("--xhash-only", action="store_true")
    args = parser.parse_args()
    TEX.mkdir(exist_ok=True)
    if args.xhash_only:
        export_xhash()
        print("Exported XHash vector PDF.")
    elif args.fig6d_only:
        export_fig6d()
    elif args.benchmarks_only:
        export_benchmarks()
        print("Exported 4 benchmark vector PDFs using recorded timings.")
    else:
        export_complexity()
        export_ao()
        export_benchmarks()
        print("Exported 14 vector PDFs; experiments were not rerun.")


if __name__ == "__main__":
    main()
