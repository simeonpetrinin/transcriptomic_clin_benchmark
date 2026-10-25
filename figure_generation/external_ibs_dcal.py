#Supplementary figures: clinical-only against fused IBS and D-calibration per clinical subset, external cohorts
#D-calibration panels switch to log axes where the range spans more than a factor of 8
#Usage, from the project root: python figure_generation/external_ibs_dcal.py --results-dir results --out-dir figures
import argparse, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
matplotlib.rcParams["font.family"] = ["Inter", "sans-serif", "DejaVu Sans"]  #font of the paper's figures; DejaVu Sans if Inter is not installed
from matplotlib.ticker import MaxNLocator, NullFormatter, LogLocator, ScalarFormatter

ap = argparse.ArgumentParser()
ap.add_argument("--results-dir", default="results"); ap.add_argument("--out-dir", default="figures")
args = ap.parse_args(); os.makedirs(args.out_dir, exist_ok=True)
RES, OUT = os.path.join(args.results_dir, "external"), os.path.join(args.out_dir, "external")
COHORTS = [("GSE39582", "COAD"), ("GSE68465", "LUAD"), ("GSE65858", "HNSC"), ("GSE13507", "BLCA"),
           ("GSE84437", "STAD"), ("BEATAML", "LAML"), ("GSE14520", "LIHC"), ("CGGA_LGG", "LGG"), ("CGGA_GBM", "GBM"),
           ("METABRIC", "BRCA")]
fl = lambda x: x.astype(str).str.upper().eq("TRUE")
INK, LINE, PT = "#222222", "#8A8A8A", "#2B6CB0"

def subsets(c):
    d = pd.read_csv(os.path.join(RES, f"{c}_results.csv"), low_memory=False)
    bad = fl(d.fused_path_short) | fl(d.fused_solver_error) | d.delta_c.isna() | fl(d.is_anchor)
    return d[~bad].groupby("subset_id")[["ibs_clin", "ibs_fused", "dcal_clin", "dcal_fused"]].mean()

data = {c: subsets(c) for c, _ in COHORTS}
for metric, label, lower in [("ibs", "IBS", True), ("dcal", "D-calibration (χ²)", True)]:
    nrow = -(-len(COHORTS) // 3)
    fig, axes = plt.subplots(nrow, 3, figsize=(11, 11 * nrow / 3), dpi=150)
    for ax in axes.flat[len(COHORTS):]: ax.axis("off")
    for ax, (c, t) in zip(axes.flat, COHORTS):
        s = data[c][[f"{metric}_clin", f"{metric}_fused"]].dropna()
        x, y = s.iloc[:, 0].values, s.iloc[:, 1].values
        lo, hi = min(x.min(), y.min()), max(x.max(), y.max()); pad = (hi - lo) * 0.06 or 0.01
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color=LINE, lw=1, ls="--", zorder=1)
        ax.scatter(x, y, s=12 if len(x) > 100 else 22, color=PT, alpha=0.6, lw=0, zorder=2)
        if len(x) > 2 and np.ptp(x) > 0:
            b = np.polyfit(x, y, 1); xx = np.array([x.min(), x.max()])
            ax.plot(xx, b[1] + b[0] * xx, color=INK, lw=1.6, zorder=3)
            ax.text(0.04, 0.95, f"slope {b[0]:.2f}\n{len(x)} subsets", transform=ax.transAxes, va="top", fontsize=8.5, color=INK)
        ax.set_xlim(lo - pad, hi + pad); ax.set_ylim(lo - pad, hi + pad)
        if metric == "dcal" and hi / lo > 8:
            ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lo * 0.9, hi * 1.1); ax.set_ylim(lo * 0.9, hi * 1.1)
            for a_ in (ax.xaxis, ax.yaxis): a_.set_minor_formatter(NullFormatter()); a_.set_major_formatter(ScalarFormatter())
        if ax.get_xscale() == "linear": ax.xaxis.set_major_locator(MaxNLocator(4)); ax.yaxis.set_major_locator(MaxNLocator(5))
        ax.set_title(f"{c} ({t})", fontsize=10, color=INK)
        ax.tick_params(labelsize=8, colors="#555555")
        for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
    fig.supxlabel(f"Clinical model {label}", fontsize=11); fig.supylabel(f"Fused model {label}", fontsize=11)
    fig.tight_layout(rect=(0.02, 0.02, 1, 1))
    fig.savefig(f"{OUT}_{metric}.png", facecolor="white"); plt.close(fig)
