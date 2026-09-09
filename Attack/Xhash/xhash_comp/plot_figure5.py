"""Reproduce Figure 5 in the manuscript's Figure 2 benchmark style.

Run: python3 plot_figure5.py
Dependencies: matplotlib, PyMuPDF.
The four standalone vector PDFs are intended for the same LaTeX subfigure
layout used by Figure 2; the combined PDF/PNG is a convenient preview.
"""
from pathlib import Path
import csv
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import MultipleLocator
import fitz

DATA_DIR = Path(__file__).resolve().parent
ROOT = DATA_DIR.parent.parent / "tex"
with (DATA_DIR / "figure5_data.csv").open(newline="") as f:
    rows = [{key: float(value) for key, value in row.items()}
            for row in csv.DictReader(f)]

# Same family, sizes, colours and panel dimensions as
# benchmark_DixonMagmaMsolve.ipynb, plotting cell (Figure 2).
available_fonts = {f.name for f in font_manager.fontManager.ttflist}
serif = [name for name in ("Times New Roman", "DejaVu Serif", "Georgia")
         if name in available_fonts]
plt.rcParams.update({
    "font.family": "serif", "font.serif": serif or ["DejaVu Serif"],
    "font.size": 26, "axes.labelsize": 28, "axes.titlesize": 30,
    "xtick.labelsize": 26, "ytick.labelsize": 26, "legend.fontsize": 24,
    "axes.spines.top": True, "axes.spines.right": True,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

for k in range(1, 5):
    part = sorted((r for r in rows if r["k"] == k), key=lambda r: r["steps"])
    x = [int(r["steps"]) for r in part]
    direct, hybrid, gb = ([r[key] for r in part]
                          for key in ("direct", "hybrid", "gb"))
    for r in part:
        m = 12 * int(r["steps"]) + k
        dreg = m * (int(r["alpha"]) - 1) + 1
        assert abs(r["gb"] - r["omega"] * math.log2(math.comb(m+dreg,m))) < 1e-8

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.set_title(f"CICO-{k}", fontweight="bold")
    ax.set_xlabel("Number of steps (s)")
    ax.set_ylabel(r"Complexity ($\log_2$)")
    # Match Figure 2's orange/square, blue/circle, green/triangle palette.
    # Alternate marker positions on the close direct and GB series.
    ax.plot(x, direct, "s-", color="#ff7f0e", label="Method 1",
            linewidth=3, markersize=10, markevery=[0, 2, 4], zorder=4)
    ax.plot(x, gb, "o--", color="#1f77b4", label="Gröbner basis",
            linewidth=3, markersize=10, markerfacecolor="white",
            markeredgewidth=2, markevery=[1, 3, 5], zorder=3)
    ax.plot(x, hybrid, "^-", color="#2ca02c",
            label="Method 3" if k == 1 else "Method 1+3",
            linewidth=3, markersize=10, zorder=2)
    ymax = {1:1000, 2:1000, 3:1600, 4:2200}[k]
    ax.set_xlim(.8, 6.2)
    ax.set_ylim(0, ymax)
    ax.set_xticks(x)
    ax.yaxis.set_major_locator(MultipleLocator(200 if k <= 2 else 400))
    ax.grid(True, alpha=.3, which="both", linewidth=.8)
    ax.grid(True, alpha=.15, which="minor", linewidth=.5)
    ax.set_axisbelow(True)
    ax.legend(loc="upper left", framealpha=.95)
    fig.tight_layout()
    fig.savefig(ROOT / f"figure5_cico{k}.pdf", bbox_inches="tight",
                metadata={"Title":f"XHash12 CICO-{k} complexity"})
    plt.close(fig)

# Compose vector panels in the same 0.48-textwidth, 2x2 arrangement.
# No global title, footnotes or shared legend: these belong in the caption.
panel_width = 720
panel_height = 576
gap_x, gap_y = 28, 4
combined = fitz.open()
page = combined.new_page(width=2*panel_width+gap_x, height=2*panel_height+gap_y)
for k in range(1, 5):
    left = ((k-1) % 2) * (panel_width+gap_x)
    top = ((k-1) // 2) * (panel_height+gap_y)
    rect = fitz.Rect(left, top, left+panel_width, top+panel_height)
    with fitz.open(ROOT / f"figure5_cico{k}.pdf") as panel:
        page.show_pdf_page(rect, panel, 0)
combined.set_metadata({"title":"XHash12 complexity comparison, Figure 2 style"})
combined.save(ROOT / "xhash_complexity.pdf", garbage=4, deflate=True)
page.get_pixmap(matrix=fitz.Matrix(1.25, 1.25), alpha=False).save(
    ROOT / "xhash_complexity.png")
combined.close()
print("Created four vector panels and the combined PDF/PNG.")
