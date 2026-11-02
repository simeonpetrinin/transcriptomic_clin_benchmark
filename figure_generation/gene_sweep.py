#Supplementary gene-count sweep: A dC, -dIBS and -dD-calibration, B genes selected, per robustness subset, 100 / 500 / all genes
#Only cells that ran in all three settings; each dot is one subset (fold -> repeat -> mean)
#500 genes = the main lattice with the rescue overlaid; rescued cells have no recomputed D-calibration and are left out of that row
#Usage, from the project root: python figure_generation/gene_sweep.py --results-dir results --out-dir figures
import argparse, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

CODES = ["BLCA", "BRCA", "CESC", "COAD", "ESCA", "GBM", "HNSC", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
         "LUAD", "LUSC", "OV", "PAAD", "READ", "SARC", "SKCM", "STAD", "UCEC"]
KEYS = ["variables", "repeat_id", "fold"]  #main-lattice and robustness subset ids differ, the variable strings do not
SETTINGS = [("g100", "100 genes", "_robust_g100"), ("g500", "500 genes", "MAIN"), ("gall", "All genes", "_robust_gall")]
COLOURS = {"g100": "#2a78d6", "g500": "#eb6834", "gall": "#1baf7a"}
METRICS = [("delta_c", "ΔC-index (Uno)", "higher = RNA helps"), ("delta_ibs", "−ΔIBS", "higher = RNA helps"),
           ("delta_dcal", "−ΔD-calibration (χ², 5 bins)", "higher = RNA helps"),
           ("fused_n_genes", "Genes selected\n(mean per subset)", "RNA genes kept by IPF-LASSO")]
SIGN = {"delta_ibs": -1, "delta_dcal": -1}  #flipped so that up means RNA helps, as for delta C
PANELS = ["A", "A", "A", "B"]  #panel letter per row: A = the gene sweep, B = genes selected

#500-gene cells: main lattice, rescued cells replaced by the rescue run (D-calibration set missing there)
def main_with_rescue(results_dir, code):
    cols = ["subset_id"] + KEYS + [k for k, _, _ in METRICS]
    m = pd.read_csv(os.path.join(results_dir, f"{code}_results.csv"), usecols=cols)
    f = os.path.join(results_dir, f"{code}_rescue.csv")
    if os.path.exists(f):
        r = pd.read_csv(f, usecols=["subset_id", "repeat_id", "fold", "delta_c", "delta_ibs", "n_genes"])
        r = r.rename(columns={"delta_c": "rc", "delta_ibs": "ri", "n_genes": "rg"})
        m = m.merge(r, on=["subset_id", "repeat_id", "fold"], how="left", indicator=True)
        hit = m["_merge"] == "both"
        m.loc[hit, "delta_c"] = m.loc[hit, "rc"]; m.loc[hit, "delta_ibs"] = m.loc[hit, "ri"]
        m.loc[hit, "fused_n_genes"] = m.loc[hit, "rg"]; m.loc[hit, "delta_dcal"] = np.nan  #not recomputed by the rescue
        m = m.drop(columns=["rc", "ri", "rg", "_merge"])
    return m.drop(columns="subset_id")

#per-subset values for every cancer, metric and setting, plus the cells dropped because a setting failed
def collect(results_dir):
    rows, dropped = [], {}
    for code in CODES:
        sweep = [os.path.join(results_dir, f"{code}{suf}_results.csv") for _, _, suf in SETTINGS if suf != "MAIN"]
        if not all(os.path.exists(f) for f in sweep + [os.path.join(results_dir, f"{code}_results.csv")]): continue
        m = None
        for lab, _, suf in SETTINGS:
            if suf == "MAIN": x = main_with_rescue(results_dir, code)
            else: x = pd.read_csv(os.path.join(results_dir, f"{code}{suf}_results.csv"), usecols=KEYS + [k for k, _, _ in METRICS])
            x = x.rename(columns={k: f"{k}|{lab}" for k, _, _ in METRICS})
            m = x if m is None else m.merge(x, on=KEYS)  #inner join keeps only the robustness cells
        have = [lab for lab, _, _ in SETTINGS]
        ok = m[[f"delta_c|{lab}" for lab in have]].notna().all(axis=1)  #cells that ran in every setting
        dropped[code] = (int((~ok).sum()), len(m))
        m = m[ok]
        for k, _, _ in METRICS:
            cc = [f"{k}|{lab}" for lab in have]; mk = m[m[cc].notna().all(axis=1)]  #D-calibration also drops rescued cells
            per_sub = mk.groupby(["variables", "repeat_id"])[cc].mean().groupby("variables").mean()
            for lab in have:
                for sid, v in per_sub[f"{k}|{lab}"].items():
                    rows.append(dict(cancer=code, metric=k, setting=lab, subset=sid, value=SIGN.get(k, 1) * v,
                                     sign_flipped=k in SIGN))  #stored values already carry the flip
    return pd.DataFrame(rows), dropped

#four stacked rows, cancers ordered by the 500-gene median delta C; drawn at print size (170 mm wide), so the font sizes are the printed sizes
def draw(S, dropped, out_png):
    order = S[(S.setting == "g500") & (S.metric == "delta_c")].groupby("cancer").value.median().sort_values(ascending=False).index.tolist()
    fig, axes = plt.subplots(len(METRICS), 1, figsize=(6.7, 7.6), dpi=300, sharex=True); fig.patch.set_facecolor("#ffffff")
    w = 0.26; rng = np.random.default_rng(1)
    for r, (ax, (k, ylab, note)) in enumerate(zip(axes, METRICS)):
        ax.set_facecolor("#ffffff")
        if r == 0 or PANELS[r] != PANELS[r - 1]:  #letter on the first row of each panel
            ax.text(-0.075, 1.06, PANELS[r], transform=ax.transAxes, fontsize=9, fontweight="bold", color="#0b0b0b", va="bottom")
        for i, code in enumerate(order):
            for j, (lab, _, _) in enumerate(SETTINGS):
                v = S[(S.cancer == code) & (S.setting == lab) & (S.metric == k)].value.values
                if not len(v): continue
                p = i + (j - 1) * w
                ax.boxplot(v, positions=[p], widths=w * 0.8, patch_artist=True, showfliers=False,
                           boxprops=dict(facecolor=COLOURS[lab] + "40", edgecolor=COLOURS[lab], lw=0.5),
                           medianprops=dict(color=COLOURS[lab], lw=0.9), whiskerprops=dict(color=COLOURS[lab], lw=0.5),
                           capprops=dict(color=COLOURS[lab], lw=0.5))
                ax.scatter(p + rng.uniform(-w * 0.22, w * 0.22, len(v)), v, s=1.2, color=COLOURS[lab], alpha=.85, zorder=3, linewidths=0)
        if k != "fused_n_genes": ax.axhline(0, color="#898781", lw=0.5, zorder=1)  #no zero line for gene counts
        ax.set_ylabel(ylab, color="#0b0b0b", fontsize=6.5, labelpad=2)
        ax.text(0.0, 1.015, note, transform=ax.transAxes, va="bottom", fontsize=5.5, color="#52514e")
        ax.grid(axis="y", color="#e1e0d9", lw=.4); ax.set_axisbelow(True)
        for s in ["top", "right"]: ax.spines[s].set_visible(False)
        for s in ["left", "bottom"]: ax.spines[s].set_color("#c3c2b7"); ax.spines[s].set_linewidth(0.6)
        ax.tick_params(colors="#52514e", labelsize=6, width=0.5, length=2, pad=1.5)
        for i in range(len(order) - 1): ax.axvline(i + .5, color="#e1e0d9", lw=.4)
    axes[-1].set_xticks(range(len(order))); axes[-1].set_xticklabels(order, color="#0b0b0b", fontsize=6, rotation=90)
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=COLOURS[lab] + "40", edgecolor=COLOURS[lab], lw=0.5) for lab, _, _ in SETTINGS]
    axes[0].legend(handles, [n for _, n, _ in SETTINGS], frameon=False, ncol=3, loc="lower right", bbox_to_anchor=(1.0, 0.98), fontsize=6, handlelength=1.2, columnspacing=1.0)
    fig.suptitle(f"RNA lift over the clinical model across the robustness subsets, {len(order)} cancers,\nby candidate gene count (IPF-LASSO)",
                 x=0.01, ha="left", fontsize=7.5, color="#0b0b0b", fontweight="bold")
    big = [f"{c} {d} of {n}" for c, (d, n) in dropped.items() if n and d / n >= 0.2]  #cancers losing at least 20% of cells
    foot = ("Δ = clinical + RNA minus clinical; IBS and D-calibration are shown as −Δ so that up means RNA helps in every row.\nOne dot = one subset (fold → repeat → mean); cells that ran in every setting only. "
            "Ordered by 500-gene ΔC median.\n500 genes = main lattice with rescued cells; rescued cells are left out of the −ΔD-calibration row (not recomputed).\n"
            "Cells dropped because a setting failed: " + (", ".join(big) if big else "under 20% in every cancer") + ".")
    fig.text(0.01, 0.005, foot, fontsize=5.5, color="#6f6e69", linespacing=1.3)
    plt.tight_layout(rect=(0, 0.075, 1, 0.965), h_pad=1.0); plt.savefig(out_png, facecolor=fig.get_facecolor(), dpi=920)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--out-dir", default="figures")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    S, dropped = collect(a.results_dir)
    S.to_csv(os.path.join(a.out_dir, "gene_sweep_boxplot_data.csv"), index=False)
    draw(S, dropped, os.path.join(a.out_dir, "gene_sweep_boxplot.png"))
    print(S.groupby(["metric", "cancer", "setting"]).value.median().unstack().round(4).to_string())