#Shared data loading for Figures 2-4
#Per cancer: fused-arm subset summaries, fold-level lifts with the rescue overlaid, the in-fold baseline picks,
#the published baselines and the number of genes kept per fit
import os
import numpy as np, pandas as pd
from scipy.stats import t as tdist

CODES = ["BLCA", "BRCA", "CESC", "COAD", "ESCA", "GBM", "HNSC", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
         "LUAD", "LUSC", "OV", "PAAD", "READ", "SARC", "SKCM", "STAD", "UCEC"]

def flag(x): return x.astype(str).str.upper().isin(["TRUE", "1"])

#fold-level lifts for one cancer, rescue overlaid as in the summaries; IBS lift sign-flipped so up = RNA helps
def fold_lifts(results_dir, code):
    cols = ["variables", "n_vars", "is_anchor", "repeat_id", "fold", "c_clin", "ibs_clin", "delta_c", "delta_ibs",
            "fused_path_short", "fused_solver_error", "n_train", "n_test"]
    m = pd.read_csv(os.path.join(results_dir, f"{code}_results.csv"), usecols=cols, low_memory=False)
    bad = flag(m.fused_path_short) | flag(m.fused_solver_error) | m.delta_c.isna()
    f = os.path.join(results_dir, f"{code}_rescue.csv")
    if os.path.exists(f):
        r = pd.read_csv(f); r = r[r.error.isna() | (r.error.astype(str) == "")]
        key = m.variables + "|" + m.repeat_id.astype(str) + "|" + m.fold.astype(str)
        rk = r.variables + "|" + r.repeat_id.astype(str) + "|" + r.fold.astype(str)
        for col in ["delta_c", "delta_ibs"]:
            hit = key.map(dict(zip(rk, r[col]))); m.loc[hit.notna(), col] = hit[hit.notna()]
        bad = bad & key.map(dict(zip(rk, r.delta_c))).isna()
    m.loc[bad, ["delta_c", "delta_ibs"]] = np.nan
    m["lift_c"], m["lift_ibs"] = m.delta_c, -m.delta_ibs
    return m[~flag(m.is_anchor)], float((m.n_test / m.n_train).mean())

#mean and 90% interval with the Nadeau-Bengio correction
def nb_interval(x, tt):
    x = np.asarray(x, float); x = x[np.isfinite(x)]; J = len(x)
    if J < 2: return np.nan, np.nan, np.nan
    h = tdist.ppf(0.95, J - 1) * np.sqrt(np.var(x, ddof=1) * (1 / J + tt))
    return x.mean(), x.mean() - h, x.mean() + h

#best and worst subsets picked inside each outer training fold (by = "pl" or "c"), with their outer clinical C and RNA lift
def picks(results_dir, code, by):
    s = pd.read_csv(os.path.join(results_dir, f"{code}_subset_summary.csv")); s = s[(s.arm == "fused") & ~s.is_anchor]
    sel = pd.read_csv(os.path.join(results_dir, f"{code}_baseline_selection.csv"))
    m, tt = fold_lifts(results_dir, code)
    key = m.variables + "|" + m.repeat_id.astype(str) + "|" + m.fold.astype(str)
    lift = dict(zip(key, m.lift_c)); cclin = dict(zip(key, m.c_clin))
    out = {}
    for end in ["best", "worst"]:
        v = sel[f"{end}_{by}_variables"]
        k = v + "|" + sel.repeat_id.astype(str) + "|" + sel.fold.astype(str)
        out[end] = pd.DataFrame(dict(variables=v, c_clin=k.map(cclin), lift=k.map(lift)))
    return s, out, tt

#fused-arm subset means, anchor dropped
def summary(results_dir, code):
    s = pd.read_csv(os.path.join(results_dir, f"{code}_subset_summary.csv"))
    s = s[(s.arm == "fused") & ~s.is_anchor.astype(str).str.upper().eq("TRUE")].copy()
    s["vset"] = s.variables.map(lambda t: frozenset(t.split("+")))
    return s

#lift and clinical C of the published benchmarks' variable sets (runs without rescue)
def published(results_dir, code):
    f = os.path.join(results_dir, f"{code}_published_subset_summary_norescue.csv"); out = {}
    if not os.path.exists(f): return out
    p = pd.read_csv(f); p = p[p.arm == "fused"]
    p = p.merge(pd.read_csv(os.path.join(results_dir, f"{code}_published_specs.csv"))[["subset_id", "spec"]], on="subset_id")
    for k in ["survboard", "jardillier", "herrmann"]:
        r = p[p.spec == k]
        if len(r): out[k] = (r.delta_c.iloc[0], r.c_ref.iloc[0])
    return out

#genes kept per fit (one outer fold of one clinical subset): IPF-LASSO with the rescue overlaid, offset, and the RNA-only anchor
def genes(results_dir, code):
    cols = ["variables", "is_anchor", "repeat_id", "fold", "fused_n_genes", "offset_n_genes", "fused_path_short", "fused_solver_error",
            "offset_path_short", "offset_solver_error"]
    m = pd.read_csv(os.path.join(results_dir, f"{code}_results.csv"), usecols=lambda c: c in cols, low_memory=False)
    anc = flag(m.is_anchor); fused = m.fused_n_genes.astype(float).copy()
    bad = flag(m.fused_path_short) | flag(m.fused_solver_error)
    f = os.path.join(results_dir, f"{code}_rescue.csv")
    if os.path.exists(f):  #rescued cells take the rescue fit's gene count, as their lift does
        r = pd.read_csv(f); r = r[r.error.isna() | (r.error.astype(str) == "")]
        key = m.variables + "|" + m.repeat_id.astype(str) + "|" + m.fold.astype(str)
        hit = key.map(dict(zip(r.variables + "|" + r.repeat_id.astype(str) + "|" + r.fold.astype(str), r.n_genes)))
        fused[hit.notna()] = hit[hit.notna()]; bad = bad & hit.isna()
    fused[bad] = np.nan
    off = m.offset_n_genes.astype(float).where(~(flag(m.offset_path_short) | flag(m.offset_solver_error))) if "offset_path_short" in m else m.offset_n_genes.astype(float)
    return dict(fused=fused[~anc].dropna().values, offset=off[~anc].dropna().values, anchor=fused[anc].dropna().values)

#everything per cancer, computed once; cancers ordered by median lattice lift, smallest on top
def collect(results_dir):
    D = {}
    for c in CODES:
        s = summary(results_dir, c); full = s.loc[s.n_vars.idxmax()]; age = s[s.variables == "age"]
        _, pl, tt = picks(results_dir, c, "pl"); _, cc, _ = picks(results_dir, c, "c")
        bm, blo, bhi = nb_interval(pl["best"].lift, tt); wm, wlo, whi = nb_interval(cc["worst"].lift, tt)
        cref = dict(zip(s.vset, s.c_ref))
        last = pd.Series({v: full.c_ref - cref.get(full.vset - {v}, np.nan) for v in full.vset}).dropna().sort_values(ascending=False)
        d = dict(s=s, full=(full.delta_c, full.c_ref), full_sd=full.delta_c_sd, age=(age.delta_c.iloc[0], age.c_ref.iloc[0]) if len(age) else (np.nan, np.nan),
                 best=(bm, nb_interval(pl["best"].c_clin, tt)[0]), best_ci=(blo, bhi), weak=(wm, wlo, whi), last=last, max=s.delta_c.max(),
                 best_folds=pl["best"].lift.dropna().values)
        d.update(published(results_dir, c)); d["genes"] = genes(results_dir, c); D[c] = d
    order = sorted(CODES, key=lambda c: np.median(D[c]["s"].delta_c))
    return D, order
