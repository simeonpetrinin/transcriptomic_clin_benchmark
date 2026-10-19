#Figure 4: individual clinical variables (IPF-LASSO, rescued results)
#A-B Shapley values and additivity, drawn by variable_importance.py; C1 clinical C and C2 RNA's lift by number of clinical variables
#The two parts are drawn separately and stacked at 3300 px width
#Usage, from the project root: python figure_generation/fig4.py --results-dir results --out-dir figures
import argparse, os, tempfile
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
matplotlib.rcParams["font.family"] = ["Inter", "sans-serif", "DejaVu Sans"]  #font of the paper's figures; DejaVu Sans if Inter is not installed
from common import collect
import variable_importance as vi

SURF, INK, INK2, MUTED, GRID, BAND = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#efeee9"
MARGIN = 0.01

#panel C: clinical C and RNA's lift by number of clinical variables, one thin line per cancer and the median over cancers
def panel_c(D, order, out_png):
    fig, (a, b) = plt.subplots(1, 2, figsize=(13, 4.6), dpi=150, gridspec_kw=dict(wspace=0.22)); fig.patch.set_facecolor(SURF)
    for ax in (a, b):
        ax.set_facecolor(SURF); ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True); ax.tick_params(colors=INK2, labelsize=8)
        for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
    rows = []
    for c in order:
        s = D[c]["s"]; top = s.c_ref.max() - 0.5
        g = s.groupby("n_vars").agg(med=("c_ref", "median"), best=("c_ref", "max"), lift=("delta_c", "median"))
        g["med_share"], g["best_share"] = (g.med - 0.5) / top, (g.best - 0.5) / top
        a.plot(g.index, g.med_share, color="#9cbfe9", lw=0.8, alpha=0.8); b.plot(g.index, g.lift, color="#f2b98a", lw=0.8, alpha=0.8)
        rows += [dict(cancer=c, n_vars=k, med_share=r.med_share, best_share=r.best_share, lift=r.lift) for k, r in g.iterrows()]
    R = pd.DataFrame(rows); M = R.groupby("n_vars").median(numeric_only=True); cnt = R.groupby("n_vars").cancer.nunique()
    keep = cnt[cnt >= 5].index  #median over cancers only where at least 5 cancers have that many variables
    a.plot(keep, M.loc[keep, "med_share"], color="#1c5aa6", lw=2.6, marker="o", ms=4, label="median subset of that size (median over cancers)")
    a.plot(keep, M.loc[keep, "best_share"], color="#1c5aa6", lw=1.6, ls="--", label="best subset of that size (picked after seeing results)")
    b.plot(keep, M.loc[keep, "lift"], color="#c4561a", lw=2.6, marker="o", ms=4, label="median over cancers")
    b.axhspan(-MARGIN, MARGIN, color=BAND, lw=0, zorder=0); b.axhline(0, color=MUTED, lw=0.9)
    a.set_xlabel("Number of clinical variables in the baseline", fontsize=9, color=INK2); b.set_xlabel("Number of clinical variables in the baseline", fontsize=9, color=INK2)
    a.set_ylabel("Share of the cancer's largest\nclinical C gain over 0.5", fontsize=9, color=INK2); b.set_ylabel("RNA's lift in Uno's C\n(median subset of that size)", fontsize=9, color=INK2)
    a.legend(frameon=False, fontsize=8, loc="lower right"); b.legend(frameon=False, fontsize=8, loc="upper right")
    a.set_title("C1  Clinical C grows with each added variable, then levels off", fontsize=10, color=INK, loc="left", fontweight="bold")
    b.set_title("C2  RNA's lift shrinks as variables are added", fontsize=10, color=INK, loc="left", fontweight="bold")
    plt.savefig(out_png, facecolor=SURF, bbox_inches="tight"); plt.close(fig)
    return R

#stack PNGs at the same width, with a panel letter where given
def stitch(pngs, letters, out_png, width=3300):
    from PIL import Image, ImageDraw, ImageFont
    ims = []
    for p in pngs:
        im = Image.open(p).convert("RGB"); ims.append(im.resize((width, int(im.height * width / im.width))))
    out = Image.new("RGB", (width, sum(i.height for i in ims)), (252, 252, 251)); y = 0; dr = ImageDraw.Draw(out)
    font = ImageFont.truetype(os.path.join(matplotlib.get_data_path(), "fonts", "ttf", "DejaVuSans-Bold.ttf"), 60)
    for im, L in zip(ims, letters):
        out.paste(im, (0, y))
        if L: dr.text((20, y + 10), L, fill=(11, 11, 11), font=font)
        y += im.height
    out.save(out_png)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default="results"); ap.add_argument("--out-dir", default="figures")
    a = ap.parse_args(); os.makedirs(a.out_dir, exist_ok=True)
    R, curves, T, P, A, B = vi.collect(a.results_dir)
    D, order = collect(a.results_dir)
    with tempfile.TemporaryDirectory() as tmp:
        ab, c = os.path.join(tmp, "ab.png"), os.path.join(tmp, "c.png")
        vi.main_figure(R, P, A, B, ab)
        panel_c(D, order, c)
        stitch([ab, c], [None, "C"], os.path.join(a.out_dir, "fig4.png"))
