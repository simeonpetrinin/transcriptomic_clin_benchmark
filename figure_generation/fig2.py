#Figure 2: does RNA add to the clinical model? (IPF-LASSO, rescued results)
#A mean lift over the best baseline picked inside each training fold, 90% Nadeau-Bengio interval; B the same lifts fold by fold
#C RNA vs the strongest clinical variable, each added last; D genes kept per fit by IPF-LASSO, offset and RNA alone
#Usage, from the project root: python figure_generation/fig2.py --results-dir results --out-dir figures
import argparse, os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from common import collect

SURF, INK, INK2, MUTED, BAND = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#efeee9"
MARGIN = 0.01
BOX = {"helps": "#c98500", "negligible": "#8f8172", "cant_tell": "#d6d4cc", "hurts": "#4a3aa7"}
BOX_LAB = {"helps": "RNA helps", "negligible": "negligible (inside ±0.01)", "cant_tell": "can't tell", "hurts": "RNA hurts"}
IPF, OFF, RNA, CLIN, SPREAD = "#c0304f", "#8e5bb5", "#e07b22", "#3d3c39", "#7fa7d9"
XLIM_B = (-0.15, 0.15)  #fold lifts beyond this are counted at the edge, not drawn
XMAX_D = 30
SCALE = 0.82  #font scale for the 2 x 2 print layout, 180 mm wide
FS = lambda x: round(x * SCALE, 2)
matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "DejaVu Sans"})  #editable text in the PDF

def box(lo, hi):
    if not np.isfinite(lo): return "cant_tell"
    if lo >= -MARGIN and hi <= MARGIN: return "negligible"
    if lo > 0: return "helps"
    if hi < 0: return "hurts"
    return "cant_tell"

def base(ax, n, xlabel, band=True):
    ax.set_facecolor(SURF); ax.set_ylim(n - 0.5, -0.6)
    for sp in ["top", "right", "left"]: ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color("#c3c2b7"); ax.tick_params(colors=INK2, labelsize=FS(8)); ax.tick_params(axis="y", length=0)
    for i in range(n): ax.axhline(i, color="#ecebe6", lw=0.6, zorder=0)
    if band: ax.axvspan(-MARGIN, MARGIN, color=BAND, lw=0, zorder=0); ax.axvline(0, color=MUTED, lw=0.8, zorder=0)
    ax.set_xlabel(xlabel, fontsize=FS(8.5), color=INK2); ax.set_yticks(range(n)); ax.set_yticklabels([])

def title(ax, letter, text, sub):
    ax.set_title(f"{letter}   {text}", fontsize=FS(10.5), color=INK, loc="left", fontweight="bold", pad=18)
    ax.text(0, 1.012, sub, transform=ax.transAxes, fontsize=FS(8), color=INK2, va="bottom")

def legend(ax, hs, ncol):
    ax.legend(handles=hs, loc="upper center", bbox_to_anchor=(0.5, -0.115), ncol=ncol, frameon=False, fontsize=FS(7.5), handletextpad=0.3, columnspacing=0.8)

def panel_a(ax, D, order):
    n = len(order); base(ax, n, "RNA's lift in Uno's C")
    for i, c in enumerate(order):
        bm, (lo, hi) = D[c]["best"][0], D[c]["best_ci"]; bx = box(lo, hi)
        ax.plot([lo, hi], [i, i], color=INK2 if bx == "cant_tell" else INK, lw=1.4, zorder=2, solid_capstyle="round")
        ax.scatter([bm], [i], s=50, facecolor=BOX[bx], edgecolor=INK, lw=0.7, zorder=3)
    ax.set_yticklabels(order, fontsize=FS(9), color=INK)
    title(ax, "A", "Does RNA add to the best clinical model?", "mean over folds · bar = 90% interval of the mean")
    legend(ax, [Line2D([], [], marker="o", ls="", markersize=7, markerfacecolor=BOX[k], markeredgecolor=INK, label=BOX_LAB[k])
                for k in ["helps", "negligible", "cant_tell"]], 3)

def panel_b(ax, D, order):  #spread of the fold-level lifts behind panel A
    n = len(order); base(ax, n, "RNA's lift in Uno's C, single folds"); rng = np.random.default_rng(3)
    for i, c in enumerate(order):
        f = D[c]["best_folds"]; inside = f[(f >= XLIM_B[0]) & (f <= XLIM_B[1])]
        if len(inside) > 2 and np.ptp(inside) > 0:
            v = ax.violinplot(inside, positions=[i], vert=False, widths=0.8, showextrema=False, bw_method=0.4)
            for b in v["bodies"]: b.set_facecolor(SPREAD); b.set_edgecolor("none"); b.set_alpha(0.35)
        ax.scatter(inside, i + rng.uniform(-0.18, 0.18, len(inside)), s=6, color=SPREAD, alpha=0.9, lw=0, zorder=2)
        ax.plot([np.median(f)] * 2, [i - 0.3, i + 0.3], color=INK, lw=1.2, zorder=3)
        lo_n, hi_n = int((f < XLIM_B[0]).sum()), int((f > XLIM_B[1]).sum())
        if lo_n: ax.text(XLIM_B[0], i, f"◂{lo_n}", ha="left", va="center", fontsize=FS(6.5), color=INK2)
        if hi_n: ax.text(XLIM_B[1], i, f"{hi_n}▸", ha="right", va="center", fontsize=FS(6.5), color=INK2)
    ax.set_xlim(*XLIM_B)
    title(ax, "B", "How much does one fold vary?", "one dot per fold · tick = median · ◂n = folds beyond the axis")
    legend(ax, [Line2D([], [], marker="o", ls="", markersize=4, color=SPREAD, label="one fold"),
                Line2D([], [], color=INK, marker="|", ls="", markersize=9, label="median")], 2)

def panel_c(ax, D, order):  #RNA vs the strongest clinical variable, both added last
    n = len(order); base(ax, n, "Gain in Uno's C when added to everything else")
    for i, c in enumerate(order):
        g = D[c]["last"]; rna = D[c]["full"][0]; tv, tval = g.index[0], g.iloc[0]
        ax.scatter(g.values[1:], [i] * (len(g) - 1), marker="|", s=40, color="#bdbbb3", lw=1, zorder=1)
        ax.plot([rna, tval], [i, i], color="#d9d7cf", lw=2.2, zorder=1, solid_capstyle="round")
        ax.scatter([tval], [i], s=34, color=CLIN, zorder=3); ax.scatter([rna], [i], s=46, marker="D", color=RNA, edgecolor=INK, lw=0.6, zorder=4)
        ax.annotate(tv, (max(tval, rna), i), xytext=(6, 0), textcoords="offset points", va="center", fontsize=FS(7), color=CLIN)  #label right of both marks
    ax.set_xlim(right=max(D[x]["last"].iloc[0] for x in order) + 0.035)
    title(ax, "C", "RNA vs the strongest clinical variable", "each added last to all other inputs")
    legend(ax, [Line2D([], [], marker="D", ls="", color=RNA, markeredgecolor=INK, label="RNA"),
                Line2D([], [], marker="o", ls="", color=CLIN, label="strongest clinical variable"),
                Line2D([], [], marker="|", ls="", color="#bdbbb3", markersize=9, label="other clinical variables")], 3)

def panel_d(ax, D, order):  #genes kept per fit, half violins
    n = len(order); base(ax, n, "Genes kept per fold", band=False)
    for i, c in enumerate(order):
        G = D[c]["genes"]
        for key, col, side in [("fused", IPF, -1), ("offset", OFF, 1)]:
            x = np.clip(G[key], 0, XMAX_D)
            if len(x) < 2 or np.ptp(x) == 0: ax.plot([x[0] if len(x) else 0] * 2, [i, i + side * 0.3], color=col, lw=2); continue
            v = ax.violinplot(x, positions=[i], vert=False, widths=0.85, showextrema=False, bw_method=0.35)
            for b in v["bodies"]:
                verts = b.get_paths()[0].vertices; verts[:, 1] = np.clip(verts[:, 1], i, np.inf) if side == 1 else np.clip(verts[:, 1], -np.inf, i)
                b.set_facecolor(col); b.set_edgecolor("none"); b.set_alpha(0.55)
            m = np.median(x); ax.plot([m, m], [i, i + side * 0.32], color=INK, lw=1.2)
        if len(G["anchor"]): ax.scatter([min(np.median(G["anchor"]), XMAX_D)], [i], s=30, marker="D", facecolor="white", edgecolor=INK, lw=0.7, zorder=4)
        ax.text(1.02, i, f"{100 * np.mean(G['fused'] == 0):.0f}%", transform=ax.get_yaxis_transform(), va="center", fontsize=FS(7.5), color=INK)
    ax.text(1.02, -0.9, "no\ngenes", transform=ax.get_yaxis_transform(), va="bottom", fontsize=FS(7.5), color=INK2)
    ax.set_xlim(-0.8, XMAX_D + 0.5)
    title(ax, "D", "How many genes do the models keep?", "one value per fold · tick = median")
    legend(ax, [Patch(facecolor=IPF, alpha=0.7, label="IPF-LASSO"), Patch(facecolor=OFF, alpha=0.7, label="offset model"),
                Line2D([], [], marker="D", ls="", markerfacecolor="white", markeredgecolor=INK, label="RNA alone")], 3)

#2 x 2 panels at double-column width (180 mm = 7.09 in), saved as vector PDF and 600 dpi PNG
def draw(D, order, out_base):
    fig = plt.figure(figsize=(7.09, 9.4)); fig.patch.set_facecolor(SURF)
    gs = fig.add_gridspec(2, 2, wspace=0.16, hspace=0.36, left=0.075, right=0.915, top=0.94, bottom=0.075)
    a, b, c, d = (fig.add_subplot(gs[r, k]) for r in range(2) for k in range(2))
    for ax, f in zip([a, b, c, d], [panel_a, panel_b, panel_c, panel_d]): f(ax, D, order)
    c.set_yticklabels(order, fontsize=FS(9), color=INK)  #left column carries the cancer names
    for ax in (a, b, c, d):
        for coll in ax.collections:  #smaller markers for the smaller canvas
            if hasattr(coll, "get_sizes") and len(coll.get_sizes()): coll.set_sizes(coll.get_sizes() * 0.45)
    fig.savefig(out_base + ".pdf", facecolor=SURF); fig.savefig(out_base + ".png", facecolor=SURF, dpi=600); plt.close(fig)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default="results"); ap.add_argument("--out-dir", default="figures")
    a = ap.parse_args(); os.makedirs(a.out_dir, exist_ok=True)
    D, order = collect(a.results_dir)
    draw(D, order, os.path.join(a.out_dir, "fig2"))
