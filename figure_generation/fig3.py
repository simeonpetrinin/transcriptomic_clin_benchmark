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

SURF, INK, INK2, MUTED, GRID, ZEBRA, BAND = "#ffffff", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#f3f2ee", "#efeee9"
MARGIN = 0.01
MARK = {"age": ("age only", "o", "#ffffff"), "full": ("full main set", "s", "#262522"),
        "best": ("best, picked inside each training fold", "D", "#2a9d6f"), "survboard": ("SurvBoard variables (Wissel 2025)", "^", "#2a78d6"),
        "jardillier": ("Jardillier 2022 variables", "v", "#8e5bb5"), "herrmann": ("Herrmann 2021 variables", "*", "#e07b22")}

FS = dict(tick=6.5, axis=6.5, title=7.5, legend=6, flag=5.5)  #font sizes in points at print size (170 mm wide)

def rows_axis(ax, order, xlabel, zebra=True, band=True):  #one row per cancer
    n = len(order); ax.set_facecolor(SURF)
    for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
    for sp in ["left", "bottom"]: ax.spines[sp].set_linewidth(0.6)
    ax.tick_params(colors=INK2, labelsize=FS["tick"], width=0.6, length=2.5, pad=2); ax.set_ylim(n - 0.5, -0.6)
    if zebra:
        for i in range(0, n, 2): ax.axhspan(i - 0.5, i + 0.5, color=ZEBRA, lw=0, zorder=0)
    if band: ax.axvspan(-MARGIN, MARGIN, color=BAND, lw=0, zorder=0); ax.axvline(0, color=MUTED, lw=0.6, zorder=0)
    ax.grid(axis="x", color=GRID, lw=0.4); ax.set_axisbelow(True)
    ax.set_xlabel(xlabel, fontsize=FS["axis"], color=INK2, labelpad=2)

def label(ax, letter, title):
    ax.set_title(f"{letter}  {title}", fontsize=FS["title"], color=INK, loc="left", fontweight="bold", pad=4)

#drawn at print size (170 mm wide), so the font sizes are the printed sizes
def draw(D, order, out_png):
    n = len(order); rng = np.random.default_rng(1)
    fig, (a, b, c3) = plt.subplots(1, 3, figsize=(6.7, 0.29 * n + 2.4), dpi=300, gridspec_kw=dict(width_ratios=[2.2, 1.0, 1.6], wspace=0.08))
    fig.patch.set_facecolor(SURF)
    fig.subplots_adjust(left=0.065, right=0.99, top=0.94, bottom=0.235)
    rows_axis(a, order, "RNA's lift in Uno's C\n(clinical + RNA − clinical only)")
    rows_axis(b, order, "Percentile of the baseline's clinical C\nwithin the cancer's lattice\n(100 = strongest)", band=False)
    rows_axis(c3, order, "RNA's lift in Uno's C")
    for i, c in enumerate(order):
        d = D[c]; L = d["s"].delta_c.values; p10, p50, p90 = np.percentile(L, [10, 50, 90]); clin = d["s"].c_ref.values
        a.scatter(L, i + rng.uniform(-0.22, 0.22, len(L)), s=0.8 if len(L) > 300 else 1.6, color="#b9b7ae", alpha=0.5, lw=0, zorder=1)
        a.plot([p10, p90], [i, i], color="#3d3c39", lw=2.2, alpha=0.8, solid_capstyle="butt", zorder=2)
        a.plot([p50, p50], [i - 0.25, i + 0.25], color=INK, lw=0.9, zorder=3)
        for k, (lab, mk, col) in MARK.items():
            if k in d and np.isfinite(d[k][0]):
                sz = 22 if mk == "*" else 11
                a.scatter([d[k][0]], [i], marker=mk, s=sz, facecolor=col, edgecolor=INK, lw=0.4, zorder=4)
                b.scatter([100 * np.mean(clin <= d[k][1])], [i], marker=mk, s=sz, facecolor=col, edgecolor=INK, lw=0.4, zorder=4)
        verdict_changes = p10 <= MARGIN < p90
        if verdict_changes: a.text(0.995, i, "depends", transform=a.get_yaxis_transform(), ha="right", va="center", fontsize=FS["flag"], color="#2a78d6", fontweight="bold", zorder=2)
        bm, (wm, wlo, whi) = d["best"][0], d["weak"]; blo, bhi = d["best_ci"]
        c3.plot([bm, wm], [i, i], color="#b9b7ae", lw=1.6, zorder=1)
        c3.plot([blo, bhi], [i - 0.12] * 2, color="#2a9d6f", lw=0.7, zorder=2); c3.plot([wlo, whi], [i + 0.12] * 2, color="#c0304f", lw=0.7, zorder=2)
        c3.scatter([bm], [i], s=12, color="#2a9d6f", edgecolor=INK, lw=0.4, zorder=3, marker="D")
        c3.scatter([wm], [i], s=12, color="#c0304f", edgecolor=INK, lw=0.4, zorder=3, marker="o")
        c3.scatter([d["max"]], [i], s=20, color="white", edgecolor=INK, lw=0.6, zorder=3, marker="*")
    a.set_xlim(right=0.155)  #a few subsets beyond 0.15 are clipped
    a.set_yticks(range(n)); a.set_yticklabels(order, fontsize=FS["tick"], color=INK); a.tick_params(axis="y", length=0)
    for ax in (b, c3): ax.set_yticks(range(n)); ax.set_yticklabels([]); ax.tick_params(axis="y", length=0)
    b.set_xlim(-3, 103); b.set_xticks([0, 25, 50, 75, 100])
    label(a, "A", "Lift over every clinical baseline"); label(b, "B", "How strong is each\nbaseline?")
    label(c3, "C", "How far can baseline choice\npush RNA's lift?")
    hs = [Line2D([], [], marker=mk, ls="", markersize=5.5 if mk == "*" else 4, markerfacecolor=col, markeredgecolor=INK, markeredgewidth=0.4, label=lab) for lab, mk, col in MARK.values()]
    hs += [Line2D([], [], marker="o", ls="", markersize=2.5, color="#b9b7ae", label="one clinical subset"),
           Line2D([], [], color="#3d3c39", lw=2.2, alpha=0.8, label="10th–90th percentile of the lattice"),
           Line2D([], [], color=INK, marker="|", ls="", markersize=5, label="lattice median"),
           Line2D([], [], ls="", label="depends = verdict changes with\nthe baseline (P10 ≤ 0.01 < P90)")]
    fig.legend(handles=hs, loc="upper left", bbox_to_anchor=(0.01, 0.135), ncol=2, frameon=False, fontsize=FS["legend"],
               handletextpad=0.4, columnspacing=1.0, labelspacing=0.45, borderaxespad=0)
    hs = [Line2D([], [], marker="D", ls="", markersize=3.5, color="#2a9d6f", markeredgecolor=INK, markeredgewidth=0.4, label="honest: strongest baseline, picked in-fold (partial likelihood)"),
          Line2D([], [], marker="o", ls="", markersize=3.5, color="#c0304f", markeredgecolor=INK, markeredgewidth=0.4, label="honest: weakest baseline, picked in-fold (inner-CV C)"),
          Line2D([], [], marker="*", ls="", markersize=5.5, color="white", markeredgecolor=INK, markeredgewidth=0.6, label="cherry-picked: best-looking baseline, chosen after seeing results"),
          Line2D([], [], color="#2a9d6f", lw=0.7, label="90% interval, strongest (Nadeau–Bengio)"),
          Line2D([], [], color="#c0304f", lw=0.7, label="90% interval, weakest (Nadeau–Bengio)"),
          Line2D([], [], color="#b9b7ae", lw=1.6, label="shift in lift from strongest to weakest baseline")]
    fig.legend(handles=hs, loc="upper left", bbox_to_anchor=(0.555, 0.135), ncol=1, frameon=False, fontsize=FS["legend"],
               handletextpad=0.4, labelspacing=0.45, borderaxespad=0)
    plt.savefig(out_png, facecolor=SURF, bbox_inches="tight", pad_inches=0.03, dpi=920); plt.close(fig)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default="results"); ap.add_argument("--out-dir", default="figures")
    a = ap.parse_args(); os.makedirs(a.out_dir, exist_ok=True)
    D, order = collect(a.results_dir)
    draw(D, order, os.path.join(a.out_dir, "fig3.png"))
