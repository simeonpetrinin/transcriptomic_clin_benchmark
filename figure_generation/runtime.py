#Supplementary runtime figure: A median seconds per model fit on the robustness subsets (IPF-ridge and RNA score: on the full
#lattice, <CODE>_rs run), B hours for the full TCGA lattice (measured for the full-lattice models, projected from seconds per
#fit for CoxBoost, blockForest and the clinical RSF), C hours per external cohort,
#D IPF-LASSO hours with and without the rescue refits, E median seconds per fit of the RNA-only model
#Usage, from the project root: python figure_generation/runtime.py --results-dir results --out-dir figures
import argparse, glob, os, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

CODES = ["BLCA", "BRCA", "CESC", "COAD", "ESCA", "GBM", "HNSC", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
         "LUAD", "LUSC", "OV", "PAAD", "READ", "SARC", "SKCM", "STAD", "UCEC"]
#key, label, seconds column, colour, marker, run on the full lattice?
SERIES = [("clin", "Clinical Cox", "secs_clin", "#898781", "o", True),
          ("rfclin", "Clinical RSF", "secs_rfclin", "#52514e", "v", False),
          ("fused", "IPF-LASSO", "secs_fused", "#2a78d6", "o", True),
          ("offset", "Offset (Volkmann)", "secs_offset", "#eb6834", "s", True),
          ("ridge", "IPF-ridge", "secs_ridge", "#1baf7a", "D", True),
          ("cb", "CoxBoost", "secs_cb", "#eda100", "^", False),
          ("bf", "blockForest", "secs_bf", "#e87ba4", "p", False),
          ("score", "RNA score", "secs_score", "#008300", "P", True)]
RS_KEYS = {"ridge", "score"}  #read from the full-lattice <CODE>_rs run, not from the robustness run
STYLE = {k: (lab, col, mk, not full) for k, lab, _, col, mk, full in SERIES}  #label, colour, marker, projected?
STYLE.update({"glmnet": ("glmnet only", "#2a78d6", "o", False),
              "rescue": ("with rescue refits", "#4a3aa7", "D", False)})
MAIN_HOST = "SMONI_LENOVO"  #the main laptop; cancers run on the other machine are marked with a dagger
RNA_ONLY = ["fused", "ridge", "cb", "bf", "score"]  #learners with an RNA-only fit; the offset model is the same fit as IPF-LASSO without clinical variables

#spread the right-end labels so they never overlap (min gap in axis units: log10 or linear)
def dodge(ys, gap):
    idx = np.argsort(ys); out = np.array(ys, float)[idx]
    for i in range(1, len(out)): out[i] = max(out[i], out[i - 1] + gap)
    res = np.empty_like(out); res[idx] = out; return res

#mean seconds per fit at each number of clinical variables, filled in by straight lines where a size was not sampled
def per_size(rob, col, sizes):
    m = rob.groupby("n_vars")[col].mean().dropna()
    return pd.Series(np.interp(sizes, m.index.values, m.values), index=sizes)

#panels A, B, D, E: the TCGA cancers
def collect_tcga(results_dir):
    rows, resc, rna = [], [], []
    for code in CODES:
        main = pd.read_csv(os.path.join(results_dir, f"{code}_results.csv"), low_memory=False,
                           usecols=lambda c: c in ["n_vars", "is_anchor", "host", "n_train", "n_test", "secs_clin", "secs_fused", "secs_offset"])
        rob = pd.read_csv(os.path.join(results_dir, f"{code}_robust_results.csv"), low_memory=False)
        rs = pd.read_csv(os.path.join(results_dir, f"{code}_rs_results.csv"), low_memory=False,
                         usecols=lambda c: c in ["n_vars", "is_anchor", "secs_ridge", "secs_score"])
        base = dict(cancer=code, patients=int((main.n_train + main.n_test).iloc[0]), host=main.host.iloc[0])
        cells = main.groupby("n_vars").size()  #lattice cells at each size (anchor = size 0)
        lat = rob[~rob.is_anchor.astype(str).str.upper().eq("TRUE")]  #panel A leaves out the omics-only anchor (no clinical model)
        lat_rs = rs[~rs.is_anchor.astype(str).str.upper().eq("TRUE")]
        for key, lab, col, _, _, full in SERIES:
            if key in RS_KEYS:
                rows.append(dict(base, series=key, label=lab, secs_per_fit=lat_rs[col].median(), lattice_hours=rs[col].sum() / 3600,
                                 lattice_how="measured", lattice_cells=int(cells.sum())))
                continue
            if col not in rob: continue
            if full: total, how = main[col].sum() / 3600, "measured"
            else: total, how = float((per_size(rob, col, cells.index.values) * cells).sum()) / 3600, "projected"
            rows.append(dict(base, series=key, label=lab, secs_per_fit=lat[col].median(), lattice_hours=total,
                             lattice_how=how, lattice_cells=int(cells.sum())))
        #D: glmnet run, then the same plus the penalized refits of the failed cells
        f = os.path.join(results_dir, f"{code}_rescue.csv")
        extra = pd.read_csv(f).secs.sum() / 3600 if os.path.exists(f) else 0.0
        glm = main.secs_fused.sum() / 3600
        resc += [dict(base, series="glmnet", hours=glm), dict(base, series="rescue", hours=glm + extra, rescue_hours=extra)]
        #E: the omics-only anchor, RNA block only, every learner that fits one
        anc = rob[rob.is_anchor.astype(str).str.upper().eq("TRUE")]
        anc_rs = rs[rs.is_anchor.astype(str).str.upper().eq("TRUE")]
        for key in RNA_ONLY:
            col = next(c for k, _, c, *_ in SERIES if k == key)
            if key in RS_KEYS: rna.append(dict(base, series=key, secs_per_fit=anc_rs[col].median(), n_fits=int(anc_rs[col].notna().sum())))
            elif col in anc: rna.append(dict(base, series=key, secs_per_fit=anc[col].median(), n_fits=int(anc[col].notna().sum())))
    return pd.DataFrame(rows), pd.DataFrame(resc), pd.DataFrame(rna)

#panel C: the external cohorts, from their timing summaries (last row = the finished run)
def collect_external(results_dir):
    rows = []
    ext = os.path.join(results_dir, "external")
    for f in sorted(glob.glob(os.path.join(ext, "timing", "*_timing.csv"))):
        if re.search(r"_(fold_timing|timing_repeats)\.csv$", f): continue
        t = pd.read_csv(f).iloc[-1]
        man = os.path.join(ext, f"{t.cancer}_manifest.txt")
        comp = re.search(r"tcga_comparator\s*:\s*(\S+)", open(man).read()).group(1) if os.path.exists(man) else "?"
        for key, col in [("clin", "compute_secs_clin"), ("fused", "compute_secs_fused"), ("offset", "compute_secs_offset")]:
            rows.append(dict(cancer=t.cancer, comparator=comp, patients=int(t.n_patients), series=key,
                             lattice_hours=t[col] / 3600, lattice_cells=int(t.n_cells)))
    return pd.DataFrame(rows)

#one panel of lines across cohorts, labelled at the right end
def panel(ax, D, val, order, ticks, ylab, xlab, letter, log=True):
    x = np.arange(len(order)); pos = dict(zip(order, x)); ends = []
    ax.set_facecolor("#fcfcfb")
    if log: ax.set_yscale("log")
    for key in D.series.unique():
        lab, colour, mk, proj = STYLE[key]
        proj = proj and val == "lattice_hours"  #dashed, hollow = projected, not measured
        d = D[D.series == key].assign(x=lambda z: z.cancer.map(pos)).sort_values("x")
        ax.plot(d.x, d[val], color=colour, lw=2, ls="--" if proj else "-", zorder=2)
        ax.scatter(d.x, d[val], s=46, marker=mk, zorder=3, linewidths=1.6,
                   facecolor="#fcfcfb" if proj else colour, edgecolor=colour)
        ends.append((lab, d[val].iloc[-1]))
    lo, hi = ax.get_ylim()
    tf = (np.log10, lambda v: 10 ** v) if log else (lambda v: v, lambda v: v)
    ys = dodge([tf[0](v) for _, v in ends], (tf[0](hi) - tf[0](lo)) * 0.045)
    for (lab, v), y in zip(ends, ys):  #direct labels at the right end, nudged apart, with a thin leader
        ax.plot([x[-1] + 0.15, x[-1] + 0.4], [v, tf[1](y)], color="#c3c2b7", lw=0.8)
        ax.text(x[-1] + 0.45, tf[1](y), lab, color="#0b0b0b", fontsize=8.5, va="center")
    ax.set_xticks(x); ax.set_xticklabels(ticks, rotation=90, fontsize=8.5, color="#0b0b0b")
    ax.set_xlim(-0.5, x[-1] + 0.15 * len(order) + 1.2)
    ax.set_xlabel(xlab, color="#52514e"); ax.set_ylabel(ylab, color="#0b0b0b", fontsize=10.5)
    ax.text(-0.07, 1.02, letter, transform=ax.transAxes, fontsize=16, fontweight="bold", color="#0b0b0b")
    ax.grid(axis="y", which="major", color="#e1e0d9", lw=.6); ax.set_axisbelow(True)
    for s in ["top", "right"]: ax.spines[s].set_visible(False)
    for s in ["left", "bottom"]: ax.spines[s].set_color("#c3c2b7")
    ax.tick_params(colors="#52514e")

def draw(D, X, R, G, out_png):
    tc = D.drop_duplicates("cancer").sort_values("patients")
    t_order = tc.cancer.tolist()
    t_ticks = [f"{c}{'' if h == MAIN_HOST else '†'} ({p})" for c, p, h in zip(tc.cancer, tc.patients, tc.host)]
    xc = X.drop_duplicates("cancer").sort_values("patients")
    x_ticks = [f"{c} ({k}, {p})" for c, k, p in zip(xc.cancer, xc.comparator, xc.patients)]
    xlab = "Cancer (patients), ordered by number of patients"
    fig = plt.figure(figsize=(21, 15.5), dpi=150); fig.patch.set_facecolor("#fcfcfb")
    gs = fig.add_gridspec(2, 6, height_ratios=[1, 1], left=0.045, right=0.99, top=0.95, bottom=0.11, wspace=0.75, hspace=0.42)  #fixed margins: tight_layout gives up with the right-end labels
    axA, axB = fig.add_subplot(gs[0, :3]), fig.add_subplot(gs[0, 3:])
    axC, axD, axE = fig.add_subplot(gs[1, :2]), fig.add_subplot(gs[1, 2:4]), fig.add_subplot(gs[1, 4:])
    panel(axA, D, "secs_per_fit", t_order, t_ticks, "Median runtime per model fit (seconds, 1 core)", xlab, "A")
    panel(axB, D, "lattice_hours", t_order, t_ticks, "Runtime for the full lattice (hours, 1 core)", xlab, "B")
    panel(axC, X, "lattice_hours", xc.cancer.tolist(), x_ticks, "Runtime for the full lattice (hours, 1 core)",
          "External cohort (TCGA comparator, patients)", "C")
    R = R.sort_values("series", key=lambda z: z.map({"rescue": 0, "glmnet": 1}), kind="stable")  #glmnet drawn on top: unrescued cancers stay blue
    panel(axD, R, "hours", t_order, t_ticks, "IPF-LASSO runtime for the full lattice (hours, 1 core)", xlab, "D")
    panel(axE, G, "secs_per_fit", t_order, t_ticks, "RNA-only model: median runtime per fit (seconds, 1 core)", xlab, "E")
    hs = [plt.Line2D([], [], color=c, marker=m, lw=2, markersize=7) for _, _, _, c, m, _ in SERIES]
    hs += [plt.Line2D([], [], color="#52514e", lw=2), plt.Line2D([], [], color="#52514e", lw=2, ls="--", marker="o", mfc="#fcfcfb")]
    fig.legend(hs, [l for _, l, *_ in SERIES] + ["measured", "projected"], loc="upper center", ncol=10, frameon=False, fontsize=9)
    plt.savefig(out_png, facecolor=fig.get_facecolor())

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--out-dir", default="figures")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    D, R, G = collect_tcga(a.results_dir)
    X = collect_external(a.results_dir)
    for name, df in [("tcga", D), ("external", X), ("rescue", R), ("rna", G)]:
        df.to_csv(os.path.join(a.out_dir, f"runtime_figure_data_{name}.csv"), index=False)
    draw(D, X, R, G, os.path.join(a.out_dir, "runtime_figure.png"))
    print("TCGA hours:\n", D.groupby(["label", "lattice_how"]).lattice_hours.sum().round(1).to_string())
    print("\nExternal hours:\n", X.groupby("series").lattice_hours.sum().round(1).to_string())
    print("\nRescue refits, hours:", round(R.rescue_hours.sum(), 1))
    print("\nRNA-only model, median seconds per fit (median over cancers):\n", G.groupby("series").secs_per_fit.median().round(2).to_string())
