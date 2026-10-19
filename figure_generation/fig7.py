#Figure 7: agreement between models (robustness subsets) and between metrics (main lattice, rescued results)
#A median RNA lift per model; B Spearman correlation of each model's lift with IPF-LASSO's across subsets
#C split-half reliability of each metric's lift and the correlations between metrics, medians over cancers
#Usage, from the project root: python figure_generation/fig7.py --results-dir results --out-dir figures
import argparse, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = ["Inter", "sans-serif", "DejaVu Sans"]  #font of the paper's figures; DejaVu Sans if Inter is not installed
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm, Normalize

ap = argparse.ArgumentParser()
ap.add_argument("--results-dir", default="results"); ap.add_argument("--out-dir", default="figures")
args = ap.parse_args(); os.makedirs(args.out_dir, exist_ok=True)
RES, OUT = args.results_dir, os.path.join(args.out_dir, "fig7.png")
CANCERS = ["BLCA", "BRCA", "CESC", "COAD", "ESCA", "GBM", "HNSC", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
           "LUAD", "LUSC", "OV", "PAAD", "READ", "SARC", "SKCM", "STAD", "UCEC"]
NO_ICI = {"BRCA", "KIRC", "KIRP", "UCEC"}  #their ICI files hold the old truncated fits
ARMS = [("fused", "IPF-Lasso\n(main)"), ("offset", "offset"), ("ridge", "IPF-ridge"), ("score", "RNA score"),
        ("coxboost", "CoxBoost"), ("blockforest_2000", "blockForest\n(vs clinical RSF)")]
INK, MUTED = "#222222", "#666666"
fl = lambda x: x.astype(str).str.upper().eq("TRUE")


def robust(c):
    d = pd.read_csv(os.path.join(RES, f"{c}_robust_subset_summary.csv"))
    d = d[~fl(d.is_anchor)]
    return d.pivot_table(index="subset_id", columns="arm", values="delta_c")


def cells(c):
    #post-rescue per-fold lifts, same overlay and bad-cell rule as scripts/02_summarise.R
    cols = ["subset_id", "variables", "is_anchor", "repeat_id", "fold", "delta_c", "delta_ibs", "delta_dcal",
            "fused_path_short", "fused_solver_error"]
    d = pd.read_csv(os.path.join(RES, f"{c}_results.csv"), usecols=cols, low_memory=False)
    d["rescued"] = False
    rp = os.path.join(RES, f"{c}_rescue.csv")
    if os.path.exists(rp):
        r = pd.read_csv(rp, low_memory=False)
        r = r[r.error.isna() | (r.error.astype(str).str.strip() == "")]
        key = lambda df: df.variables.astype(str) + "|" + df.repeat_id.astype(str) + "|" + df.fold.astype(str)
        idx = pd.Series(d.index, index=key(d))
        i = idx.reindex(key(r)).dropna().astype(int).values
        rr = r.set_index(key(r)).reindex(key(d).iloc[i].values)
        d.loc[i, "delta_c"] = rr.delta_c.values
        d.loc[i, "delta_ibs"] = rr.delta_ibs.values
        d.loc[i, "delta_dcal"] = np.nan
        d.loc[i, ["fused_path_short", "fused_solver_error"]] = False
        d.loc[i, "rescued"] = True
    bad = fl(d.fused_path_short) | fl(d.fused_solver_error) | d.delta_c.isna()
    d = d[~bad & ~fl(d.is_anchor)].copy()
    d["C"], d["IBS"], d["Dcal"] = d.delta_c, -d.delta_ibs, -d.delta_dcal
    if c not in NO_ICI:
        ici = pd.read_csv(os.path.join(RES, f"{c}_loglik_ici.csv"),
                          usecols=["subset_id", "repeat_id", "fold", "delta_ici_q50"])
        d = d.merge(ici, on=["subset_id", "repeat_id", "fold"], how="left")
        d["ICI"] = -d.delta_ici_q50
        d.loc[d.rescued, "ICI"] = np.nan
    else:
        d["ICI"] = np.nan
    return d


def reliability(d, m):
    #split-half (odd vs even repeats) Spearman across subsets, Spearman-Brown corrected
    odd = d[d.repeat_id % 2 == 1].groupby("subset_id")[m].mean()
    even = d[d.repeat_id % 2 == 0].groupby("subset_id")[m].mean()
    pair = pd.concat([odd, even], axis=1).dropna()
    if len(pair) < 5: return np.nan
    r = pair.iloc[:, 0].corr(pair.iloc[:, 1], method="spearman")
    return 2 * r / (1 + r)


MET = ["C", "IBS", "Dcal", "ICI"]
MLAB = ["ΔC", "−ΔIBS", "−ΔD-cal", "−ΔICI"]
A, B, rel, agr = {}, {}, {m: [] for m in MET}, {(a, b): [] for a in MET for b in MET}
for c in CANCERS:
    w = robust(c)
    A[c] = {a: w[a].median() * 100 for a, _ in ARMS}
    B[c] = {a: w[["fused", a]].dropna().corr(method="spearman").iloc[0, 1] for a, _ in ARMS[1:]}
    A[c]["n"] = len(w)
    d = cells(c)
    s = d.groupby("subset_id")[MET].mean()
    for m in MET:
        if s[m].notna().sum() >= 5: rel[m].append(reliability(d, m))
    for i, a in enumerate(MET):
        for b in MET[:i]:
            p = s[[a, b]].dropna()
            if len(p) >= 5: agr[(a, b)].append(p[a].corr(p[b], method="spearman"))

A = pd.DataFrame(A).T; B = pd.DataFrame(B).T
order = A.sort_values("fused").index
summ = lambda v: (np.nanmedian(v), np.nanpercentile(v, 25), np.nanpercentile(v, 75), int(np.sum(np.isfinite(v))))

fig = plt.figure(figsize=(15, 18.5), dpi=150)
gs = fig.add_gridspec(2, 2, height_ratios=[1.45, 0.8], hspace=0.14, wspace=0.08, left=0.07, right=0.98, top=0.965, bottom=0.04)

#panel A: median lift x 100
gold_violet = LinearSegmentedColormap.from_list("gv", ["#6A51A3", "#F2F0EC", "#C8860A"])
axA = fig.add_subplot(gs[0, 0])
M = A.loc[order, [a for a, _ in ARMS]].values
lim = np.nanmax(np.abs(M))
axA.imshow(M, cmap=gold_violet, norm=TwoSlopeNorm(0, -lim, lim), aspect="auto")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        axA.text(j, i, f"{M[i, j]:.1f}", ha="center", va="center", fontsize=9, color=INK)
axA.set_xticks(range(len(ARMS))); axA.set_xticklabels([l for _, l in ARMS], fontsize=9.5)
axA.set_yticks(range(len(order))); axA.set_yticklabels([f"{c} ({int(A.loc[c, 'n'])})" for c in order], fontsize=9.5)
axA.xaxis.tick_top(); axA.tick_params(length=0)
axA.set_title("A  Median RNA lift across subsets (ΔC × 100)", loc="left", fontsize=12, fontweight="bold", pad=38)
for s_ in axA.spines.values(): s_.set_visible(False)

#panel B: Spearman with IPF-Lasso
greens = LinearSegmentedColormap.from_list("g", ["#F4FAF4", "#2E8B57", "#0B4D2C"])
axB = fig.add_subplot(gs[0, 1])
Mb = B.loc[order, [a for a, _ in ARMS[1:]]].values
axB.imshow(Mb, cmap=greens, norm=Normalize(-0.2, 1), aspect="auto")
for i in range(Mb.shape[0]):
    for j in range(Mb.shape[1]):
        v = Mb[i, j]
        axB.text(j, i, "–" if np.isnan(v) else f"{v:.2f}", ha="center", va="center", fontsize=9,
                 color="white" if v > 0.6 else INK)
axB.set_xticks(range(len(ARMS) - 1)); axB.set_xticklabels([l for _, l in ARMS[1:]], fontsize=9.5)
axB.set_yticks(range(len(order))); axB.set_yticklabels(order, fontsize=9.5)
axB.xaxis.tick_top(); axB.tick_params(length=0)
axB.set_title("B  Same ranking of subsets as IPF-Lasso? (Spearman ρ)", loc="left", fontsize=12, fontweight="bold", pad=38)
for s_ in axB.spines.values(): s_.set_visible(False)

#panel C: metric reliability (diagonal) and agreement (below)
blues = LinearSegmentedColormap.from_list("b", ["#FFFFFF", "#6FA3D8", "#1F5A99"])
axC = fig.add_subplot(gs[1, :])
k = len(MET)
cell = np.full((k, k), np.nan); txt = {}
for i, a in enumerate(MET):
    med, lo, hi, n = summ(rel[a]); cell[i, i] = med
    txt[(i, i)] = f"{med:.2f}\n[{lo:.2f} to {hi:.2f}]" + (f"\nn = {n}" if n < len(CANCERS) else "")
    for j, b in enumerate(MET[:i]):
        med, lo, hi, n = summ(agr[(a, b)]); cell[i, j] = med
        txt[(i, j)] = f"{med:.2f}\n[{lo:.2f} to {hi:.2f}]" + (f"\nn = {n}" if n < len(CANCERS) else "")
axC.imshow(np.ma.masked_invalid(cell), cmap=blues, norm=Normalize(0, 1), aspect="auto")
for (i, j), t in txt.items():
    axC.text(j, i, t, ha="center", va="center", fontsize=11, color="white" if cell[i, j] > 0.6 else INK)
    if i == j:
        axC.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, lw=2, ec=INK))
axC.set_xticks(range(k)); axC.set_xticklabels(MLAB, fontsize=12)
axC.set_yticks(range(k)); axC.set_yticklabels(MLAB, fontsize=12)
axC.tick_params(length=0)
for s_ in axC.spines.values(): s_.set_visible(False)
axC.set_title("C  Agreement between metrics across clinical subsets (Spearman ρ, median over cancers)",
              loc="left", fontsize=12, fontweight="bold", pad=10)
fig.savefig(OUT, dpi=150, facecolor="white")

#numbers quoted in the text
print("median rho with IPF-Lasso:", B.median().round(2).to_dict())
for m in MET: print("reliability", m, np.round(summ(rel[m]), 2))
for (a, b), v in agr.items():
    if v: print("agreement", a, b, np.round(summ(v), 2))
