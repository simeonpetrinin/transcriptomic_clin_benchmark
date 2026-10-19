#Supplementary runtime figure: A seconds per model fit for every model, B hours for the full lattice
#A is measured on the robustness subsets (same cells and laptop for every model)
#B is measured for the full-lattice models; for the robustness-only models it is projected from their mean seconds per fit
#at each number of clinical variables times the lattice's number of cells at that size
#Usage, from the project root: python figure_generation/runtime.py --results-dir results --out-dir figures
import argparse, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
matplotlib.rcParams["font.family"] = ["Inter", "sans-serif", "DejaVu Sans"]  #font of the paper's figures; DejaVu Sans if Inter is not installed

CODES = ["BLCA", "BRCA", "CESC", "COAD", "ESCA", "GBM", "HNSC", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
         "LUAD", "LUSC", "OV", "PAAD", "READ", "SARC", "SKCM", "STAD", "UCEC"]
#key, label, seconds column, colour, marker, run on the full lattice?
SERIES = [("clin", "Clinical Cox", "secs_clin", "#898781", "o", True),
          ("rfclin", "Clinical RSF", "secs_rfclin", "#52514e", "v", False),
          ("fused", "IPF-LASSO", "secs_fused", "#2a78d6", "o", True),
          ("offset", "Offset (Volkmann)", "secs_offset", "#eb6834", "s", True),
          ("ridge", "IPF-ridge", "secs_ridge", "#1baf7a", "D", False),
          ("cb", "CoxBoost", "secs_cb", "#eda100", "^", False),
          ("bf", "blockForest", "secs_bf", "#e87ba4", "p", False),
          ("score", "RNA score", "secs_score", "#008300", "P", False)]
MAIN_HOST = "SMONI_LENOVO"  #the main laptop; cancers run on the other machine are marked with a dagger

#spread the right-end labels so they never overlap (positions in log10, min gap in decades)
def dodge(ys, gap):
    idx = np.argsort(ys); out = np.array(ys, float)[idx]
    for i in range(1, len(out)): out[i] = max(out[i], out[i - 1] + gap)
    res = np.empty_like(out); res[idx] = out; return res

#mean seconds per fit at each number of clinical variables, filled in by straight lines where a size was not sampled
def per_size(rob, col, sizes):
    m = rob.groupby("n_vars")[col].mean().dropna()
    return pd.Series(np.interp(sizes, m.index.values, m.values), index=sizes)

def collect(results_dir):
    rows = []
    for code in CODES:
        main = pd.read_csv(os.path.join(results_dir, f"{code}_results.csv"), low_memory=False,
                           usecols=lambda c: c in ["n_vars", "is_anchor", "host", "n_train", "n_test", "secs_clin", "secs_fused", "secs_offset"])
        rob = pd.read_csv(os.path.join(results_dir, f"{code}_robust_results.csv"), low_memory=False)
        n_pat = int((main.n_train + main.n_test).iloc[0])
        cells = main.groupby("n_vars").size()  #lattice cells at each size (anchor = size 0)
        lat = rob[~rob.is_anchor.astype(str).str.upper().eq("TRUE")]  #panel A leaves out the omics-only anchor (no clinical model)
        for key, lab, col, _, _, full in SERIES:
            if col not in rob: continue
            per_fit = lat[col].mean()
            if full: total, how = main[col].sum() / 3600, "measured"
            else: total, how = float((per_size(rob, col, cells.index.values) * cells).sum()) / 3600, "projected"
            rows.append(dict(cancer=code, patients=n_pat, host=main.host.iloc[0], series=key, label=lab,
                             secs_per_fit=per_fit, lattice_hours=total, lattice_how=how, lattice_cells=int(cells.sum())))
    return pd.DataFrame(rows)

def draw(D, out_png):
    order = D.drop_duplicates("cancer").sort_values("patients")
    x = np.arange(len(order)); pos = dict(zip(order.cancer, x))
    ticks = [f"{c}{'' if h == MAIN_HOST else '†'} ({p})" for c, p, h in zip(order.cancer, order.patients, order.host)]
    fig, axes = plt.subplots(1, 2, figsize=(19, 7.2), dpi=150); fig.patch.set_facecolor("#fcfcfb")
    for ax, (val, ylab, letter) in zip(axes, [("secs_per_fit", "Mean runtime per model fit (seconds, 1 core)", "A"),
                                              ("lattice_hours", "Runtime for the full lattice (hours, 1 core)", "B")]):
        ax.set_facecolor("#fcfcfb"); ax.set_yscale("log"); ends = []
        for key, lab, _, colour, mk, full in SERIES:
            d = D[D.series == key].copy()
            if d.empty: continue
            d["x"] = d.cancer.map(pos); d = d.sort_values("x")
            proj = val == "lattice_hours" and not full  #dashed, hollow = projected, not measured
            ax.plot(d.x, d[val], color=colour, lw=2, ls="--" if proj else "-", zorder=2)
            ax.scatter(d.x, d[val], s=46, marker=mk, zorder=3, linewidths=1.6,
                       facecolor="#fcfcfb" if proj else colour, edgecolor=colour)
            ends.append((lab, d[val].iloc[-1]))
        lo, hi = np.log10(ax.get_ylim()); ys = dodge([np.log10(v) for _, v in ends], (hi - lo) * 0.035)
        for (lab, v), y in zip(ends, ys):  #direct labels at the right end, nudged apart, with a thin leader
            ax.plot([x[-1] + 0.15, x[-1] + 0.4], [v, 10 ** y], color="#c3c2b7", lw=0.8)
            ax.text(x[-1] + 0.45, 10 ** y, lab, color="#0b0b0b", fontsize=8.5, va="center")
        ax.set_xticks(x); ax.set_xticklabels(ticks, rotation=90, fontsize=8.5, color="#0b0b0b")
        ax.set_xlim(-0.5, x[-1] + 3.2)
        ax.set_xlabel("Cancer (patients), ordered by number of patients", color="#52514e")
        ax.set_ylabel(ylab, color="#0b0b0b", fontsize=10.5)
        ax.text(-0.07, 1.02, letter, transform=ax.transAxes, fontsize=16, fontweight="bold", color="#0b0b0b")
        ax.grid(axis="y", which="major", color="#e1e0d9", lw=.6); ax.set_axisbelow(True)
        for s in ["top", "right"]: ax.spines[s].set_visible(False)
        for s in ["left", "bottom"]: ax.spines[s].set_color("#c3c2b7")
        ax.tick_params(colors="#52514e")
    hs = [plt.Line2D([], [], color=c, marker=m, lw=2, markersize=7) for _, _, _, c, m, _ in SERIES]
    hs += [plt.Line2D([], [], color="#52514e", lw=2), plt.Line2D([], [], color="#52514e", lw=2, ls="--", marker="o", mfc="#fcfcfb")]
    fig.legend(hs, [l for _, l, *_ in SERIES] + ["measured", "projected"], loc="upper center", ncol=10, frameon=False, fontsize=9)
    plt.tight_layout(rect=(0, 0, 1, 0.94)); plt.savefig(out_png, facecolor=fig.get_facecolor())

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--out-dir", default="figures")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    D = collect(a.results_dir)
    D.to_csv(os.path.join(a.out_dir, "runtime_figure_data.csv"), index=False)
    draw(D, os.path.join(a.out_dir, "runtime_figure.png"))
    print(D.groupby(["label", "lattice_how"]).agg(hours=("lattice_hours", "sum"), median_secs_per_fit=("secs_per_fit", "median")).round(2).to_string())