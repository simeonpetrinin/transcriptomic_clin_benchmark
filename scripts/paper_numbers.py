#Numbers quoted in the text of the paper that no figure script prints
#Each line starts with where the number appears ([Results 3.1] etc.); numbers that the figure scripts print themselves
#(slopes and crossing points, Figure 7 correlations, Figure 8 agreement and reliabilities) are read from their metrics files
#where the text counts them, and are otherwise not repeated here
#Reads the archived results (main lattice with the rescue overlaid, as in the figures), the external results in results/external,
#data/processed for the cohort tables, and fig6_metrics.csv from the figure scripts (run figure_generation/fig6.py first)
#Usage, from the project root: python scripts/paper_numbers.py --results-dir results --processed-dir data/processed --figures-dir figures
import argparse, itertools, math, os, re, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figure_generation"))
from common import CODES, collect, fold_lifts, nb_interval, flag

MARGIN = 0.01
DROP = {"GSE13507": "invasiveness"}  #as in fig7.py: subsets with this variable are left out of the external comparisons
PAIRS = [("GSE39582", "COAD"), ("GSE68465", "LUAD"), ("GSE65858", "HNSC"), ("GSE13507", "BLCA"), ("GSE84437", "STAD"),
         ("BEATAML", "LAML"), ("GSE14520", "LIHC"), ("CGGA_LGG", "LGG"), ("CGGA_GBM", "GBM"), ("METABRIC", "BRCA")]  #external, TCGA

def say(where, text): print(f"[{where}] {text}")
def pct(x): return f"{100 * x:.1f}%"

#verdict of a 90% interval against the equivalence margin d (Section 2.10)
def box(lo, hi, d):
    if not (np.isfinite(lo) and np.isfinite(hi)): return "can't tell"
    if -d <= lo and hi <= d: return "negligible"
    if lo > 0: return "helps"
    if hi < 0: return "hurts"
    return "can't tell"

#percentile rule (Section 2.10): where the 10th-90th percentile of subset lifts sits against the margin
def pct_rule(lifts):
    p10, p90 = np.percentile(lifts, [10, 90])
    return "depends" if p10 <= MARGIN < p90 else ("above" if p10 > MARGIN else "below")

#per-subset rule (Discussion 4.5): each subset classified by its own interval, then the cohort labelled from the shares
def subset_boxes(m, tt):
    out = {}
    for v, g in m.groupby("variables"):
        mm, lo, hi = nb_interval(g.lift_c.values, tt)
        out[v] = box(lo, hi, MARGIN) if g.lift_c.notna().sum() >= len(g) / 2 else "can't tell"
    return pd.Series(out)
def subset_label(b):
    s = lambda k: (b == k).mean()
    if s("can't tell") >= 0.5: return "underpowered"
    for k in ["helps", "negligible", "hurts"]:
        if s(k) >= 0.9: return k
    if s("helps") >= 0.1 and s("negligible") + s("hurts") >= 0.1: return "depends"
    return "leans " + b.value_counts().idxmax()

#mean clinical C gain of adding a variable to every subset without it (empty model = 0.5), from the subset summary
def add_one_gains(s):
    C = dict(zip(s.vset, s.c_ref)); full = s.loc[s.n_vars.idxmax(), "vset"]; out = {}
    for v in full:
        others = sorted(full - {v}); g = []
        for r in range(len(others) + 1):
            for S in itertools.combinations(others, r):
                S = frozenset(S); g.append(C.get(S | {v}, np.nan) - (C.get(S, np.nan) if S else 0.5))
        out[v] = np.nanmean(g)
    return out

#exact Shapley value of each variable for v(S) = clinical C of S minus 0.5, as a share of v(full set)
def shapley_shares(s):
    C = dict(zip(s.vset, s.c_ref)); full = s.loc[s.n_vars.idxmax(), "vset"]; n = len(full)
    val = lambda S: (C.get(frozenset(S), np.nan) - 0.5) if S else 0.0
    out = {}
    for v in full:
        others = sorted(full - {v}); phi = 0.0
        for r in range(n):
            w = math.factorial(r) * math.factorial(n - r - 1) / math.factorial(n)
            for S in itertools.combinations(others, r): phi += w * (val(set(S) | {v}) - val(S))
        out[v] = phi / val(full)
    return out

#the six equivalence margins of one cancer (Section 2.10)
def margins(m, tt, s):
    gains = add_one_gains(s); pos = [g for g in gains.values() if g > 0]
    per = m.dropna(subset=["lift_c"]).groupby("variables").lift_c.agg(["std", "size"])
    se = np.sqrt(per["std"] ** 2 * (1 / per["size"] + tt))
    mde = ((stats.t.ppf(0.95, per["size"] - 1) + stats.t.ppf(0.80, per["size"] - 1)) * se).median()
    age = gains.get("age", np.nan)
    return {"fixed 0.01": MARGIN, "Herrmann best 0.0067": 0.0067, "Herrmann helps 0.016": 0.016,
            "typical variable": float(np.median(pos)) if len(pos) >= 2 else MARGIN,
            "age": age if age > 0 else np.nan, "CV noise": float(mde)}

#genes kept per fit with the subset it belongs to (IPF-Lasso with the rescue overlaid; offset), anchor = RNA-only model
def genes_by_subset(results_dir, code):
    cols = ["variables", "is_anchor", "repeat_id", "fold", "fused_n_genes", "offset_n_genes", "fused_path_short",
            "fused_solver_error", "offset_path_short", "offset_solver_error"]
    m = pd.read_csv(os.path.join(results_dir, f"{code}_results.csv"), usecols=cols, low_memory=False)
    fused = m.fused_n_genes.astype(float).copy(); bad = flag(m.fused_path_short) | flag(m.fused_solver_error)
    f = os.path.join(results_dir, f"{code}_rescue.csv")
    if os.path.exists(f):
        r = pd.read_csv(f); r = r[r.error.isna() | (r.error.astype(str) == "")]
        key = m.variables + "|" + m.repeat_id.astype(str) + "|" + m.fold.astype(str)
        hit = key.map(dict(zip(r.variables + "|" + r.repeat_id.astype(str) + "|" + r.fold.astype(str), r.n_genes)))
        fused[hit.notna()] = hit[hit.notna()]; bad = bad & hit.isna()
    fused[bad] = np.nan
    off = m.offset_n_genes.astype(float).where(~(flag(m.offset_path_short) | flag(m.offset_solver_error)))
    return pd.DataFrame(dict(variables=m.variables, anchor=flag(m.is_anchor), fused=fused, offset=off))

#median follow-up by reverse Kaplan-Meier (censoring counted as the event)
def km_median(time, event):
    order = np.argsort(time, kind="stable"); t, e = time[order], event[order]; s = 1.0
    for u in np.unique(t):
        d = ((t == u) & (e == 1)).sum()
        if d:
            s *= 1 - d / (t >= u).sum()
            if s <= 0.5: return u
    return np.nan

def read_rds(f):
    import rdata
    return rdata.conversion.convert(rdata.parser.parse_file(f))

#n, events, k, RNA dimensions (genes at least 80% complete), subsets and CV design of one cohort
def cohort_row(rds, results_csv, manifest):
    x = read_rds(rds); cl = x["clinical"]; e = np.asarray(x["expr"], float)
    n, ev = len(cl), int(np.asarray(cl["OS"]).sum())
    sub = pd.read_csv(results_csv, usecols=["variables", "is_anchor"], low_memory=False)
    sub = sub[~flag(sub.is_anchor)].variables.unique()
    k = max(len(v.split("+")) for v in sub)
    reps = int(re.search(r"n_repeats\s*:\s*(\d+)", open(manifest).read()).group(1))
    fu = km_median(np.asarray(cl["OS.time"], float), 1 - np.asarray(cl["OS"], int)) / 365.25
    return dict(n=n, events=ev, event_pct=round(100 * ev / n, 1), k=k, genes=int(((1 - np.isnan(e).mean(0)) >= 0.8).sum()),
                subsets=len(sub), cv=f"{reps} x 5", followup_years=round(fu, 1))

# ---------------------------------------------------------------- TCGA: verdicts against the in-fold strongest baseline
def tcga(results_dir, D):
    rows = []; F = {}
    for c in CODES:
        m, tt = fold_lifts(results_dir, c); F[c] = (m, tt)
        d = D[c]; bm, (blo, bhi) = d["best"][0], d["best_ci"]
        M = margins(m, tt, d["s"])
        rows.append(dict(cancer=c, mean=bm, lo=blo, hi=bhi, verdict=box(blo, bhi, MARGIN), rule=pct_rule(d["s"].delta_c.values),
                         **{f"m:{k}": (box(blo, bhi, v) if np.isfinite(v) else "n/a") for k, v in M.items()}))
    return pd.DataFrame(rows).set_index("cancer"), F

def report_verdicts(V):
    helps = V.index[V.verdict == "helps"].tolist()
    say("Abstract; Results 3.1", f"RNA helps (in-fold strongest baseline, margin 0.01) in {len(helps)}: {', '.join(helps)}")
    say("Results 3.1", f"negligible: {', '.join(V.index[V.verdict == 'negligible'])}; can't tell: {(V.verdict == "can't tell").sum()}; hurts: {(V.verdict == 'hurts').sum()}")
    say("Results 3.1", f"lower 90% bounds of the 'helps' cancers: {V.loc[helps, 'lo'].min():.3f} to {V.loc[helps, 'lo'].max():.3f}")
    say("Discussion 4.1", "mean lift of the 'helps' cancers: " + ", ".join(f"{c} {V.loc[c, 'mean']:.3f}" for c in V.loc[helps].sort_values('mean', ascending=False).index))
    for k in [c for c in V.columns if c.startswith("m:")]:
        say("Results 3.1; File S1", f"helps under the {k[2:]} margin: {', '.join(V.index[V[k] == 'helps'])}")
    allm = [c for c in helps if all(V.loc[c, k] in ("helps", "n/a") for k in V.columns if k.startswith("m:"))]
    say("Results 3.1; Discussion 4.1", f"helps under all six margins: {len(allm)} ({', '.join(allm)})")
    rest = V[V.verdict != "helps"].hi
    say("Abstract; Results 3.1", f"of the remaining {len(rest)}, upper 90% bound at most 0.05: {(rest <= 0.05).sum()}")
    say("Results 3.1", f"upper bound at most 0.02: {(rest <= 0.02).sum()}; 0.02-0.03: {((rest > 0.02) & (rest <= 0.03)).sum()}; "
                       f"0.03-0.05: {((rest > 0.03) & (rest <= 0.05)).sum()}; above 0.05: " + ", ".join(f"{c} {v:.3f}" for c, v in rest[rest > 0.05].items()))
    say("Results 3.1", f"probability of 4 or more false 'helps' among 21 one-sided 5% tests: {stats.binom.sf(3, 21, 0.05):.3f}")
    dep = V.index[V.rule == "depends"].tolist()
    say("Abstract; Results 3.2; Conclusions", f"verdict depends on the baseline (P10 <= 0.01 < P90 of subset lifts) in {len(dep)}: {', '.join(dep)}")

# ---------------------------------------------------------------- Results 3: fits, Figure 2 text
def report_fits(results_dir):
    cells = failed = rescued = 0; ridge = [0, 0]; cb = [0, 0]; lam = {}
    for c in CODES:
        m = pd.read_csv(os.path.join(results_dir, f"{c}_results.csv"), low_memory=False,
                        usecols=["is_anchor", "delta_c", "fused_path_short", "fused_solver_error", "fused_lambda_fraction"])
        m = m[~flag(m.is_anchor)]; bad = flag(m.fused_path_short) | flag(m.fused_solver_error) | m.delta_c.isna()
        cells += len(m); failed += bad.sum()
        f = os.path.join(results_dir, f"{c}_rescue.csv")
        if os.path.exists(f):
            r = pd.read_csv(f); rescued += (r.error.isna() | (r.error.astype(str) == "")).sum()
        lf = m.fused_lambda_fraction.astype(float)[~bad].dropna()
        lam[c] = (lf <= lf.min() + 1e-9).mean() if len(lf) and np.isclose(lf.min(), 0.05) else 0.0  #share at the 0.05 end of the path
        rs = pd.read_csv(os.path.join(results_dir, f"{c}_rs_results.csv"), usecols=["is_anchor", "ridge_error"], low_memory=False)
        rs = rs[~flag(rs.is_anchor)]; ridge[0] += len(rs); ridge[1] += rs.ridge_error.notna().sum()
        rb = pd.read_csv(os.path.join(results_dir, f"{c}_robust_results.csv"), usecols=["is_anchor", "cb_error"], low_memory=False)
        rb = rb[~flag(rb.is_anchor)]; cb[0] += len(rb); cb[1] += rb.cb_error.notna().sum()
    say("Results 3", f"IPF-Lasso fits {cells}; failed in glmnet {failed} ({pct(failed / cells)}); refitted with penalized {rescued}; "
                     f"still failed and dropped {failed - rescued} ({pct((failed - rescued) / cells)})")
    say("Results 3", f"CoxBoost fits failed: {pct(cb[1] / cb[0])} of {cb[0]}; IPF-ridge: {pct(ridge[1] / ridge[0])} of {ridge[0]}")
    lam = pd.Series(lam)
    say("Results 3.4", f"selected IPF-Lasso lambda at the end of the path (0.05 of lambda max): at most {100 * lam.max():.2f}% of fits "
                       f"({lam.idxmax()}); 0% in {(lam == 0).sum()} cancers")

def report_fig2_text(results_dir, D):
    med = {c: np.median(D[c]["best_folds"]) for c in CODES}
    say("Results 3.1", f"median in-fold fold lift inside +-0.01 in {sum(abs(v) <= MARGIN for v in med.values())} of 21")
    rf = D["READ"]["best_folds"]; lo = rf.min()
    say("Results 3.1", f"READ lowest fold {lo:.3f}; READ mean {rf.mean():.3f} with it, {rf[rf != lo].mean():.3f} without it")
    rna = {c: D[c]["full"][0] for c in CODES}; clin = {c: D[c]["last"].iloc[0] for c in CODES}
    more = [c for c in CODES if clin[c] > rna[c]]
    say("Results 3.1", f"strongest clinical variable (added last) adds more than RNA (added last) in {len(more)} of 21; RNA more in "
                       + ", ".join(f"{c} ({D[c]['last'].index[0]} {clin[c]:+.4f} vs RNA {rna[c]:+.4f})" for c in CODES if c not in more))
    zero = {}; allg = {"fused": [], "offset": [], "anchor": []}; per = []
    for c in CODES:
        g = genes_by_subset(results_dir, c); lat = g[~g.anchor]
        zero[c] = (lat.fused.dropna() == 0).mean()
        for k, x in [("fused", lat.fused), ("offset", lat.offset), ("anchor", g[g.anchor].fused)]:
            x = x.dropna(); allg[k].append(pd.DataFrame(dict(cancer=c, n=x.values)))
        sub = lat.groupby("variables").fused.mean().rename("genes").to_frame().join(D[c]["s"].set_index("variables").delta_c)
        per.append(sub)
    z = pd.Series(zero).sort_values(ascending=False)
    say("Results 3.1", "share of IPF-Lasso fits with no genes: above 50%: " + ", ".join(f"{c} {pct(v)}" for c, v in z[z > 0.5].items())
        + "; 40-50%: " + ", ".join(z[(z > 0.4) & (z <= 0.5)].index) + f"; 20-40%: {((z > 0.2) & (z <= 0.4)).sum()} cancers")
    for k, lab in [("fused", "IPF-Lasso"), ("offset", "offset"), ("anchor", "RNA-only")]:
        a = pd.concat(allg[k]); i = a.n.values.argmax()
        say("Results 3.1", f"genes per fit, {lab}: mean {a.n.mean():.1f}, median {a.n.median():.0f}, largest {a.n.max():.0f} ({a.cancer.values[i]})")
    p = pd.concat(per).dropna()
    say("Results 3.1", f"Spearman of genes selected and lift across all {len(p)} subsets: {stats.spearmanr(p.genes, p.delta_c)[0]:.2f}")

# ---------------------------------------------------------------- Results 3.2-3.3: baselines and composition
def report_baselines(results_dir, D, F):
    P = {}
    for c in CODES:
        clin = D[c]["s"].c_ref.values; P[c] = {k: 100 * np.mean(clin <= D[c][k][1]) for k in ["herrmann", "jardillier", "survboard"] if k in D[c]}
        #Herrmann's variables recoded with our preprocessing (Table S13): a separate run ("herrmann_ours") when the recoded set is not a
        #lattice subset, otherwise the lattice subset itself (listed in <cancer>_published_ours.csv)
        f = os.path.join(results_dir, f"{c}_published_subset_summary_norescue.csv")
        if os.path.exists(f):
            o = pd.read_csv(f); o = o[o.arm == "fused"].merge(pd.read_csv(os.path.join(results_dir, f"{c}_published_specs.csv"))[["subset_id", "spec"]], on="subset_id")
            o = o[o.spec == "herrmann_ours"]
            if len(o): P[c]["herrmann recoded"] = 100 * np.mean(clin <= o.c_ref.iloc[0])
        f = os.path.join(results_dir, f"{c}_published_ours.csv")
        if os.path.exists(f):
            o = pd.read_csv(f); o = o[(o.paper == "herrmann") & (o.version == "ours")]
            if len(o): P[c]["herrmann recoded"] = 100 * np.mean(clin <= o.c_ref.iloc[0])
    P = pd.DataFrame(P).T
    h = P.herrmann.dropna()
    say("Abstract; Results 3.2; Discussion 4.2", f"Herrmann et al. baseline in the weakest 10% of the lattice in {(h <= 10).sum()} of {len(h)}: "
        + ", ".join(f"{c} {v:.0f}" for c, v in h[h <= 10].sort_values().items()))
    for k in ["jardillier", "survboard", "herrmann recoded"]:
        say("Results 3.2", f"{k}: median percentile {P[k].median():.0f} over {P[k].notna().sum()} cancers")
    labels = {c: subset_label(subset_boxes(*F[c])) for c in CODES}
    lab = pd.Series(labels)
    say("Results 3.2", f"per-subset rule (margin 0.01): underpowered {(lab == 'underpowered').sum()}, depends {(lab == 'depends').sum()}; "
                       + ", ".join(f"{c} {v}" for c, v in lab[lab != 'underpowered'].items()))
    gap = pd.Series({c: D[c]["weak"][0] - D[c]["best"][0] for c in CODES})
    say("Results 3.2", f"weakest (in-fold, by C) minus strongest (in-fold, by partial likelihood) baseline lift: median {gap.median():.4f}; "
                       "largest " + ", ".join(f"{c} {v:+.3f}" for c, v in gap.sort_values(ascending=False).head(5).items())
                       + "; negative " + ", ".join(f"{c} {v:+.4f} (strongest {D[c]['best'][0]:.4f}, weakest {D[c]['weak'][0]:.4f})" for c, v in gap[gap < -0.005].sort_values().items()))
    ns = pd.Series({c: F[c][0].groupby("variables").lift_c.var().mean() / F[c][0].lift_c.var() for c in CODES})
    say("Results 3.2", f"share of fold-level lift variance within subsets: median {pct(ns.median())} (range {pct(ns.min())} to {pct(ns.max())})")

def report_composition(D):
    sh = pd.DataFrame({c: shapley_shares(D[c]["s"]) for c in CODES}).T
    med = sh.loc[:, sh.notna().sum() >= 5].median().sort_values(ascending=False)  #variables in at least 5 cancers, as in Figure 4A
    say("Results 3.3", "Shapley share of clinical C (median over cancers with the variable; cancers): "
        + ", ".join(f"{v} {100 * x:.2f}% ({sh[v].notna().sum()})" for v, x in med.items()))
    share, lift = {}, {}
    for c in CODES:
        s = D[c]["s"].copy(); s["share"] = (s.c_ref - 0.5) / (s.c_ref.max() - 0.5)
        share[c] = s.groupby("n_vars").share.median(); lift[c] = s.groupby("n_vars").delta_c.median()
    share, lift = pd.DataFrame(share), pd.DataFrame(lift)
    say("Results 3.3", "share of the strongest subset's clinical gain by number of variables (median over cancers; cancers): "
        + ", ".join(f"{k}: {pct(share.loc[k].median())} ({share.loc[k].notna().sum()})" for k in share.index if k <= 10))
    say("Results 3.3", "RNA lift by number of variables (median over cancers): " + ", ".join(f"{k}: {lift.loc[k].median():.3f}" for k in lift.index if k <= 5))
    fin = pd.Series({c: D[c]["full"][1] - D[c]["best"][1] for c in CODES})
    aft = pd.Series({c: D[c]["s"].c_ref.max() - D[c]["full"][1] for c in CODES})
    say("Results 3.3; Discussion 4.3", f"clinical C, full set minus in-fold pick: median {fin.median():+.3f} (range {fin.min():+.3f} to {fin.max():+.3f})")
    say("Results 3.3", f"clinical C, best subset after the fact minus full set: median {aft.median():.3f} (range {aft.min():.3f} to {aft.max():.3f})")
    say("Discussion 4.2", f"BRCA median clinical C over subsets: {D['BRCA']['s'].c_ref.median():.2f}")

def report_slopes(results_dir, figures_dir):
    f = os.path.join(figures_dir, "fig6_metrics.csv")
    if not os.path.exists(f): say("Results 3.4", f"skipped: {f} not found (run figure_generation/fig6.py first)"); return
    m = pd.read_csv(f).set_index("cancer")
    both = m.index[(m.slope < 1) & (m.slope_offset < 1)]
    say("Results 3.4; Conclusions", f"slope below 1 for both IPF-Lasso and offset: {len(both)} of {len(m)}")
    say("Results 3.4", "slope below 0.50 for both: " + ", ".join(f"{c} {m.slope[c]:.2f}/{m.slope_offset[c]:.2f}" for c in m.index[(m.slope < 0.5) & (m.slope_offset < 0.5)]))
    say("Results 3.4", "IPF-Lasso slopes from 0.95 up: " + ", ".join(f"{c} {v:.3f}" for c, v in m.slope[m.slope >= 0.95].sort_values().items()))
    say("Results 3.4", f"CESC slopes: IPF-Lasso {m.slope['CESC']:.3f}, offset {m.slope_offset['CESC']:.2f}")
    s = pd.read_csv(os.path.join(results_dir, "CESC_subset_summary.csv")); s = s[~flag(s.is_anchor)]
    for arm in s.arm.unique():
        x = s[s.arm == arm].delta_c; say("Results 3.4", f"CESC subsets with positive lift, {arm} arm: {(x > 0).sum()} of {len(x)}")
    corr = m.slope_noise_checked - m.slope
    say("Figure 6 legend; File S1", f"reliability-corrected minus plain slope: largest {corr.max():.2f} ({corr.idxmax()}); all others below "
                                    f"{corr.drop(corr.idxmax()).max():.3f}; {corr.idxmax()} reliability {m.reliability[corr.idxmax()]:.2f}")

# ---------------------------------------------------------------- Results 3.5: external cohorts
def report_external(results_dir, D):
    E = os.path.join(results_dir, "external"); rows = []
    for ext, tc in PAIRS:
        m, tt = fold_lifts(E, ext)
        if ext in DROP: m = m[~m.variables.str.split("+").apply(lambda v: DROP[ext] in v)]  #as in Figure 7
        b = subset_boxes(m, tt)
        sub = m.groupby("variables").agg(lift=("lift_c", "mean"), c=("c_clin", "mean"))
        t = D[tc]["s"]; tl = pd.Series(t.delta_c.values, index=t.vset)
        common = [v for v in sub.index if frozenset(v.split("+")) in tl.index]
        rows.append(dict(ext=ext, tcga=tc, label=subset_label(b), own_rule=pct_rule(sub.lift.values),
                         tcga_rule=pct_rule(t.delta_c.values), n_common=len(common),
                         ext_common=pct_rule(sub.loc[common].lift.values) if len(common) >= 3 else "n/a",
                         tcga_common=pct_rule(tl.loc[[frozenset(v.split("+")) for v in common]].values) if len(common) >= 3 else "n/a",
                         shares=b.value_counts(normalize=True), lat=sub))
    X = pd.DataFrame(rows).set_index("ext")
    say("Discussion 4.5", f"external per-subset rule: underpowered {(X.label == 'underpowered').sum()} of 10; "
                          + ", ".join(f"{e} {l}" for e, l in X.label.items() if l != "underpowered"))
    for kind, a, b_ in [("own lattices", "tcga_rule", "own_rule"), ("shared subsets", "tcga_common", "ext_common")]:
        same = X.index[X[a] == X[b_]]
        say("Results 3.5", f"percentile rule, TCGA vs external on {kind}: same in {len(same)} of 10 ("
            + ", ".join(X.tcga[same]) + "); " + ", ".join(f"{X.tcga[e]} {X[a][e]} -> {e} {X[b_][e]}" for e in X.index if e not in same))
    same_both = [e for e in X.index if (X.label[e] == "depends") == (X.own_rule[e] == "depends") and X.label[e] != "underpowered"]
    say("Discussion 4.5", f"external cohorts not underpowered whose per-subset label matches the percentile rule: {', '.join(same_both)}")
    sh = X.shares["METABRIC"].reindex(["helps", "negligible", "hurts", "can't tell"]).fillna(0)
    say("Results 3.5", "METABRIC subsets: " + ", ".join(f"{k} {100 * v:.0f}%" for k, v in sh.items()) + f" of {len(X.lat['METABRIC'])}")
    bb = subset_boxes(*fold_lifts(results_dir, "BRCA"))
    say("Results 3.5", f"TCGA BRCA subsets can't tell: {100 * (bb == "can't tell").mean():.0f}%")
    s = X.lat["METABRIC"]; q = s.c.quantile([0.25, 0.75])
    say("Results 3.5; Conclusions", f"METABRIC median lift, weakest quarter of baselines (by clinical C) {s.lift[s.c <= q[0.25]].median():+.3f}, "
                                    f"strongest quarter {s.lift[s.c >= q[0.75]].median():+.3f}")
    lg = fold_lifts(E, "CGGA_LGG")[0].groupby("variables").lift_c.mean()
    say("Results 3.1; Results 3.5", f"CGGA_LGG lift: age+sex+histology {lg['age+sex+histology']:.3f}, with IDH {lg['age+sex+histology+IDH']:.3f}")

# ---------------------------------------------------------------- Tables 1-2 and cohort notes
def report_tables(results_dir, processed_dir):
    T = {}
    for c in CODES:
        T[c] = cohort_row(os.path.join(processed_dir, f"{c}_modelling.rds"), os.path.join(results_dir, f"{c}_results.csv"),
                          os.path.join(results_dir, f"{c}_manifest.txt"))
    T = pd.DataFrame(T).T.sort_values("n", ascending=False)
    print(T.drop(columns="followup_years").to_string()); tot = T[["n", "events", "subsets"]].astype(int).sum()
    say("Table 1", f"total n {tot.n}, events {tot.events} ({100 * tot.events / tot.n:.1f}%), subsets {tot.subsets}")
    fu = T.followup_years.astype(float)
    say("Discussion 4.5", "median follow-up below 3 years: " + ", ".join(f"{c} {v:.1f}" for c, v in fu[fu < 3].sort_values().items()))
    say("Discussion 4.5", "events: " + ", ".join(f"{c} {T.events[c]}" for c in ["READ", "KIRP", "CESC"]))
    X = {}
    for ext, tc in PAIRS:
        X[ext] = dict(tcga=tc, **cohort_row(os.path.join(processed_dir, "external", f"{ext}_modelling.rds"),
                                           os.path.join(results_dir, "external", f"{ext}_results.csv"),
                                           os.path.join(results_dir, "external", f"{ext}_manifest.txt")))
    X = pd.DataFrame(X).T
    print(X.drop(columns="followup_years").to_string()); tot = X[["n", "events", "subsets"]].astype(int).sum()
    say("Table 2", f"total n {tot.n}, events {tot.events}, subsets {tot.subsets}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default="results"); ap.add_argument("--processed-dir", default=os.path.join("data", "processed"))
    ap.add_argument("--figures-dir", default="figures")
    a = ap.parse_args()
    D, order = collect(a.results_dir)
    V, F = tcga(a.results_dir, D)
    report_verdicts(V)
    report_fits(a.results_dir)
    report_fig2_text(a.results_dir, D)
    report_baselines(a.results_dir, D, F)
    report_composition(D)
    report_slopes(a.results_dir, a.figures_dir)
    report_external(a.results_dir, D)
    report_tables(a.results_dir, a.processed_dir)
