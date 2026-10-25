#Figure 5: clinical C and RNA's lift by number of clinical variables (IPF-LASSO, rescued results)
#A clinical C as a share of the cancer's largest gain over 0.5, B RNA's lift; one thin line per cancer, thick line = median over cancers
#Drawn at print size (170 mm wide), so the font sizes are the printed sizes
#Usage, from the project root: python figure_generation/fig5.py --results-dir results --out-dir figures
import argparse, os
import pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
matplotlib.rcParams["font.family"] = ["Inter", "sans-serif", "DejaVu Sans"]  #font of the paper's figures; DejaVu Sans if Inter is not installed
from common import collect

SURF, INK, INK2, MUTED, GRID, BAND = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#efeee9"
MARGIN = 0.01
MIN_CANCERS = 5  #median over cancers only where at least this many cancers have that many variables

def curves(D, order):
    rows = []
    for c in order:
        s = D[c]["s"]; top = s.c_ref.max() - 0.5
        g = s.groupby("n_vars").agg(med=("c_ref", "median"), best=("c_ref", "max"), lift=("delta_c", "median"))
        g["med_share"], g["best_share"] = (g.med - 0.5) / top, (g.best - 0.5) / top
        rows += [dict(cancer=c, n_vars=k, med_share=r.med_share, best_share=r.best_share, lift=r.lift) for k, r in g.iterrows()]
    return pd.DataFrame(rows)

def figure(R, out_png):
    fig, (a, b) = plt.subplots(1, 2, figsize=(6.7, 2.9), dpi=300, gridspec_kw=dict(wspace=0.32)); fig.patch.set_facecolor(SURF)
    for ax in (a, b):
        ax.set_facecolor(SURF); ax.grid(color=GRID, lw=0.4); ax.set_axisbelow(True); ax.tick_params(colors=INK2, labelsize=7, width=0.6, length=2.5)
        for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
        for sp in ["left", "bottom"]: ax.spines[sp].set_linewidth(0.6)
        ax.set_xlabel("Number of clinical variables in the baseline", fontsize=7.5, color=INK2)
    for c, g in R.groupby("cancer"):
        a.plot(g.n_vars, g.med_share, color="#9cbfe9", lw=0.5, alpha=0.8); b.plot(g.n_vars, g.lift, color="#f2b98a", lw=0.5, alpha=0.8)
    M = R.groupby("n_vars").median(numeric_only=True); cnt = R.groupby("n_vars").cancer.nunique()
    keep = cnt[cnt >= MIN_CANCERS].index
    a.plot(keep, M.loc[keep, "med_share"], color="#1c5aa6", lw=1.6, marker="o", ms=2.5, label="median subset of that size")
    a.plot(keep, M.loc[keep, "best_share"], color="#1c5aa6", lw=1.0, ls="--", label="best subset of that size\n(picked after seeing results)")
    b.plot(keep, M.loc[keep, "lift"], color="#c4561a", lw=1.6, marker="o", ms=2.5, label="median over cancers")
    b.axhspan(-MARGIN, MARGIN, color=BAND, lw=0, zorder=0); b.axhline(0, color=MUTED, lw=0.6)
    a.set_ylabel("Share of the cancer's largest\nclinical C gain over 0.5", fontsize=7.5, color=INK2)
    b.set_ylabel("RNA's lift in Uno's C\n(median subset of that size)", fontsize=7.5, color=INK2)
    a.legend(frameon=False, fontsize=6.5, loc="lower right", handlelength=2.2); b.legend(frameon=False, fontsize=6.5, loc="upper right", handlelength=2.2)
    a.set_title("A  Clinical C grows with each added\n     variable, then levels off", fontsize=8, color=INK, loc="left", fontweight="bold", x=-0.25)
    b.set_title("B  RNA's lift shrinks as\n     variables are added", fontsize=8, color=INK, loc="left", fontweight="bold", x=-0.25)
    plt.savefig(out_png, facecolor=SURF, bbox_inches="tight", pad_inches=0.04); plt.close(fig)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default="results"); ap.add_argument("--out-dir", default="figures")
    a = ap.parse_args(); os.makedirs(a.out_dir, exist_ok=True)
    D, order = collect(a.results_dir)
    figure(curves(D, order), os.path.join(a.out_dir, "fig5.png"))
