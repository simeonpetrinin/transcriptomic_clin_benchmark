#Figure 3: does the choice of clinical baseline change RNA's lift? (IPF-LASSO, rescued results)
#A lift over every clinical subset, with marked baselines incl. the published benchmarks; B strength of each marked baseline
#C strongest (partial likelihood) and weakest (inner-CV C) baselines picked in-fold, and the best-looking one picked afterwards
#Usage, from the project root: python figure_generation/fig3.py --results-dir results --out-dir figures
import argparse, os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
matplotlib.rcParams["font.family"] = ["Inter", "sans-serif", "DejaVu Sans"]  #font of the paper's figures; DejaVu Sans if Inter is not installed
from matplotlib.lines import Line2D
from common import collect

SURF, INK, INK2, MUTED, GRID, ZEBRA, BAND = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#f3f2ee", "#efeee9"
MARGIN = 0.01
MARK = {"age": ("age only", "o", "#ffffff"), "full": ("full main set", "s", "#262522"),
        "best": ("best, picked inside each training fold", "D", "#2a9d6f"), "survboard": ("SurvBoard variables (Wissel 2025)", "^", "#2a78d6"),
        "jardillier": ("Jardillier 2022 variables", "v", "#8e5bb5"), "herrmann": ("Herrmann 2021 variables", "*", "#e07b22")}

def rows_axis(ax, order, xlabel, zebra=True, band=True):  #one row per cancer
    n = len(order); ax.set_facecolor(SURF)
    for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
    ax.tick_params(colors=INK2, labelsize=8); ax.set_ylim(n - 0.5, -0.6)
    if zebra:
        for i in range(0, n, 2): ax.axhspan(i - 0.5, i + 0.5, color=ZEBRA, lw=0, zorder=0)
    if band: ax.axvspan(-MARGIN, MARGIN, color=BAND, lw=0, zorder=0); ax.axvline(0, color=MUTED, lw=0.9, zorder=0)
    ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
    ax.set_xlabel(xlabel, fontsize=9, color=INK2)

def label(ax, letter, title):
    ax.set_title(f"{letter}  {title}", fontsize=10.5, color=INK, loc="left", fontweight="bold")

def draw(D, order, out_png):
    n = len(order); rng = np.random.default_rng(1)
    fig, (a, b, c3) = plt.subplots(1, 3, figsize=(19, 0.44 * n + 2.6), dpi=150, gridspec_kw=dict(width_ratios=[2.2, 1.0, 1.6], wspace=0.06))
    fig.patch.set_facecolor(SURF)
    rows_axis(a, order, "RNA's lift in Uno's C (clinical + RNA − clinical only)")
    rows_axis(b, order, "Percentile of the baseline's clinical C\nwithin the cancer's lattice (100 = strongest)", band=False)
    rows_axis(c3, order, "RNA's lift in Uno's C")
    for i, c in enumerate(order):
        d = D[c]; L = d["s"].delta_c.values; p10, p50, p90 = np.percentile(L, [10, 50, 90]); clin = d["s"].c_ref.values
        a.scatter(L, i + rng.uniform(-0.22, 0.22, len(L)), s=3 if len(L) > 300 else 6, color="#b9b7ae", alpha=0.45, lw=0, zorder=1)
        a.plot([p10, p90], [i, i], color="#3d3c39", lw=3.6, alpha=0.8, solid_capstyle="butt", zorder=2)
        a.plot([p50, p50], [i - 0.25, i + 0.25], color=INK, lw=1.4, zorder=3)
        for k, (lab, mk, col) in MARK.items():
            if k in d and np.isfinite(d[k][0]):
                sz = 70 if mk == "*" else 34
                a.scatter([d[k][0]], [i], marker=mk, s=sz, facecolor=col, edgecolor=INK, lw=0.7, zorder=4)
                b.scatter([100 * np.mean(clin <= d[k][1])], [i], marker=mk, s=sz, facecolor=col, edgecolor=INK, lw=0.7, zorder=4)
        verdict_changes = p10 <= MARGIN < p90
        if verdict_changes: a.text(0.995, i, "depends", transform=a.get_yaxis_transform(), ha="right", va="center", fontsize=7.5, color="#2a78d6", fontweight="bold")
        bm, (wm, wlo, whi) = d["best"][0], d["weak"]; blo, bhi = d["best_ci"]
        c3.plot([bm, wm], [i, i], color="#b9b7ae", lw=2.5, zorder=1)
        c3.plot([blo, bhi], [i - 0.12] * 2, color="#2a9d6f", lw=1, zorder=2); c3.plot([wlo, whi], [i + 0.12] * 2, color="#c0304f", lw=1, zorder=2)
        c3.scatter([bm], [i], s=40, color="#2a9d6f", edgecolor=INK, lw=0.6, zorder=3, marker="D")
        c3.scatter([wm], [i], s=40, color="#c0304f", edgecolor=INK, lw=0.6, zorder=3, marker="o")
        c3.scatter([d["max"]], [i], s=55, color="white", edgecolor=INK, lw=0.9, zorder=3, marker="*")
    a.set_xlim(right=0.155)  #a few subsets beyond 0.15 are clipped
    a.set_yticks(range(n)); a.set_yticklabels(order, fontsize=9, color=INK); a.tick_params(axis="y", length=0)
    for ax in (b, c3): ax.set_yticks(range(n)); ax.set_yticklabels([]); ax.tick_params(axis="y", length=0)
    b.set_xlim(-3, 103); b.set_xticks([0, 25, 50, 75, 100])
    label(a, "A", "Lift over every clinical baseline"); label(b, "B", "How strong is each baseline?")
    label(c3, "C", "How far can baseline choice push RNA's lift?")
    hs = [Line2D([], [], marker=mk, ls="", markersize=9 if mk == "*" else 6.5, markerfacecolor=col, markeredgecolor=INK, label=lab) for lab, mk, col in MARK.values()]
    hs += [Line2D([], [], marker="o", ls="", color="#b9b7ae", label="one clinical subset"),
           Line2D([], [], color="#3d3c39", lw=3.6, alpha=0.8, label="10th–90th percentile of the lattice"),
           Line2D([], [], color=INK, marker="|", ls="", markersize=9, label="lattice median"),
           Line2D([], [], ls="", label="depends = verdict changes with the baseline (P10 ≤ 0.01 < P90)")]
    a.legend(handles=hs, loc="upper center", bbox_to_anchor=(0.62, -0.1), ncol=3, frameon=False, fontsize=7.5)
    hs = [Line2D([], [], marker="D", ls="", color="#2a9d6f", markeredgecolor=INK, label="honest: strongest baseline, picked in-fold (partial likelihood)"),
          Line2D([], [], marker="o", ls="", color="#c0304f", markeredgecolor=INK, label="honest: weakest baseline, picked in-fold (inner-CV C)"),
          Line2D([], [], marker="*", ls="", color="white", markeredgecolor=INK, markersize=9, label="cherry-picked: best-looking baseline, chosen after seeing results"),
          Line2D([], [], color="#2a9d6f", lw=1, label="90% interval, strongest (Nadeau–Bengio)"),
          Line2D([], [], color="#c0304f", lw=1, label="90% interval, weakest (Nadeau–Bengio)"),
          Line2D([], [], color="#b9b7ae", lw=2.5, label="shift in lift from strongest to weakest baseline")]
    c3.legend(handles=hs, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=1, frameon=False, fontsize=7.5)
    plt.savefig(out_png, facecolor=SURF, bbox_inches="tight", dpi=300); plt.close(fig)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default="results"); ap.add_argument("--out-dir", default="figures")
    a = ap.parse_args(); os.makedirs(a.out_dir, exist_ok=True)
    D, order = collect(a.results_dir)
    draw(D, order, os.path.join(a.out_dir, "fig3.png"))
