#Clinical variable importance from the full subset lattice, no refitting (Figure 4 A-B and supplementary figures)
#Exact Shapley values per cancer for clinical-only C (empty set = 0.5) and for RNA's lift (empty set = the RNA-only lift),
#Shapley interaction index for every pair, and the share of the lattice explained by single variables and by pairs
#Intervals: values per outer fold, 90% interval with the Nadeau-Bengio correction, compared with +-0.01
#Usage, from the project root: python figure_generation/variable_importance.py --results-dir results --out-dir figures/variable_importance
import argparse, os
from math import factorial
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
matplotlib.rcParams["font.family"] = ["Inter", "sans-serif", "DejaVu Sans"]  #font of the paper's figures; DejaVu Sans if Inter is not installed
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

CODES = ["BLCA", "BRCA", "CESC", "COAD", "ESCA", "GBM", "HNSC", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
         "LUAD", "LUSC", "OV", "PAAD", "READ", "SARC", "SKCM", "STAD", "UCEC"]
MIN_CANCERS = 3  #a variable in at least this many cancers gets its own row in the pan-cancer grid
SURF, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
MARGIN = 0.01  #the fixed equivalence margin, also used for interactions
FOLD = {}  #(cancer, game) -> fold-level Shapley values, variables, test/train ratio, usable; for the summed row of panel A
GREY_UNSURE, GREY_NEGL = np.array(matplotlib.colors.to_rgba("#b9b7ae")), np.array(matplotlib.colors.to_rgba("#e4e3dd"))
ADDITIVE = np.array(matplotlib.colors.to_rgba("#3f9e8f"))  #pairs with the interval inside +-MARGIN

#value-suppressing colour: full colour when the interval is clear of 0, washed towards grey when it can't tell,
#light grey (variables) or teal (pairs) when the whole interval is inside +-MARGIN
def vsup(rgba, box, pair=False):
    rgba = np.array(rgba)
    if box in ("above", "below"): return rgba
    if box == "negligible": return ADDITIVE if pair else GREY_NEGL
    return 0.35 * rgba + 0.65 * GREY_UNSURE
SEQ = LinearSegmentedColormap.from_list("seq", ["#cfe0f6", "#2a78d6", "#123a6b"])  #one hue, light to dark
INTERACT = LinearSegmentedColormap.from_list("int", ["#2e2370", "#4a3aa7", "#f3f2ee", "#c98500", "#8a5a00"])  #violet = overlap, gold = synergy

def style(ax):
    ax.set_facecolor(SURF)
    for s in ["top", "right"]: ax.spines[s].set_visible(False)
    for s in ["left", "bottom"]: ax.spines[s].set_color("#c3c2b7")
    ax.tick_params(colors=INK2, labelsize=8)

#one cancer: the fused-arm summary, the per-repeat values, and the variable order of the full set
def load(results_dir, code):
    s = pd.read_csv(os.path.join(results_dir, f"{code}_subset_summary.csv"))
    p = pd.read_csv(os.path.join(results_dir, f"{code}_per_repeat.csv"))
    s, p = s[s.arm == "fused"].copy(), p[p.arm == "fused"].copy()
    full = s.loc[s.n_vars.idxmax(), "variables"].split("+")
    bit = {v: 1 << i for i, v in enumerate(full)}
    to_mask = lambda txt: 0 if txt == "(omics_only)" else sum(bit[v] for v in txt.split("+"))
    s["mask"] = s.variables.map(to_mask); p["mask"] = p.variables.map(to_mask)
    if len(s) != 2 ** len(full): raise ValueError(f"{code}: {len(s)} subsets, expected {2 ** len(full)}")
    return s, p, full

#exact Shapley and Banzhaf values from a value for every subset (array indexed by bit mask)
def shapley(v, n):
    m = np.arange(2 ** n); size = np.array([bin(x).count("1") for x in m])
    w = np.array([factorial(k) * factorial(n - k - 1) / factorial(n) for k in range(n)])
    phi, banz = np.empty(n), np.empty(n)
    for i in range(n):
        out = m[(m >> i) & 1 == 0]  #subsets without variable i
        gain = v[out | (1 << i)] - v[out]
        phi[i], banz[i] = np.sum(w[size[out]] * gain), gain.mean()
    return phi, banz

#exact Shapley interaction index for every pair: subsets of the other k-2 variables, weight |S|!(k-|S|-2)!/(k-1)!
def interactions(v, n):
    out = np.zeros((n, n))
    if n < 2: return out
    m = np.arange(2 ** n); size = np.array([bin(x).count("1") for x in m])
    w = np.array([factorial(k) * factorial(n - k - 2) / factorial(n - 1) for k in range(n - 1)])
    for i in range(n):
        for j in range(i + 1, n):
            out_s = m[((m >> i) & 1 == 0) & ((m >> j) & 1 == 0)]  #subsets without i and j
            d = v[out_s | (1 << i) | (1 << j)] - v[out_s | (1 << i)] - v[out_s | (1 << j)] + v[out_s]
            out[i, j] = out[j, i] = np.sum(w[size[out_s]] * d)
    return out

#share of the lattice's variation explained by main effects, and by main effects + pairs
def additivity(v, n):
    m = np.arange(2 ** n); bits = [((m >> i) & 1).astype(float) for i in range(n)]
    def r2(cols):
        X = np.column_stack([np.ones(len(m))] + cols); b, *_ = np.linalg.lstsq(X, v, rcond=None)
        return 1 - np.var(v - X @ b) / np.var(v) if np.var(v) > 0 else np.nan
    pairs = [bits[i] * bits[j] for i in range(n) for j in range(i + 1, n)]
    return r2(bits), (r2(bits + pairs) if n >= 3 else 1.0)

#fold-level quantity: 90% interval with the Nadeau-Bengio correction, then four boxes
#negligible = inside +-MARGIN; above/below = whole interval above/below 0 and not inside; cant_tell = anything else
def p1_box(x, test_train):
    x = np.asarray(x, float); x = x[np.isfinite(x)]; J = len(x)
    if J < 2: return np.nan, np.nan, np.nan, "cant_tell"
    from scipy.stats import t as tdist
    se = np.sqrt(np.var(x, ddof=1) * (1 / J + test_train))  #Nadeau & Bengio 2003
    m = x.mean(); h = tdist.ppf(0.95, J - 1) * se; lo, hi = m - h, m + h
    if lo >= -MARGIN and hi <= MARGIN: box = "negligible"
    elif lo > 0: box = "above"
    elif hi < 0: box = "below"
    else: box = "cant_tell"
    return m, lo, hi, box

#fold-level values: every outer fold holds the whole lattice on one split, so Shapley values and interactions exist per fold
#the rescue is overlaid as in the summaries; failed cells leave their fold out for RNA's lift
def fold_level(results_dir, code, full):
    n = len(full); bit = {v: 1 << i for i, v in enumerate(full)}
    cols = ["variables", "repeat_id", "fold", "c_clin", "delta_c", "fused_path_short", "fused_solver_error", "n_train", "n_test"]
    m = pd.read_csv(os.path.join(results_dir, f"{code}_results.csv"), usecols=cols, low_memory=False)
    flag = lambda x: x.astype(str).str.upper().isin(["TRUE", "1"])
    bad = flag(m.fused_path_short) | flag(m.fused_solver_error) | m.delta_c.isna()
    f = os.path.join(results_dir, f"{code}_rescue.csv")
    if os.path.exists(f):  #same overlay as 02_summarise.R
        r = pd.read_csv(f); r = r[r.error.isna() | (r.error.astype(str) == "")]
        key = m.variables + "|" + m.repeat_id.astype(str) + "|" + m.fold.astype(str)
        rk = dict(zip(r.variables + "|" + r.repeat_id.astype(str) + "|" + r.fold.astype(str), r.delta_c))
        hit = key.map(rk); m.loc[hit.notna(), "delta_c"] = hit[hit.notna()]; bad = bad & hit.isna()
    m.loc[bad, "delta_c"] = np.nan
    m["mask"] = m.variables.map(lambda t: 0 if t == "(omics_only)" else sum(bit[v] for v in t.split("+")))
    test_train = float((m.n_test / m.n_train).mean())
    out = {"clinical": ([], []), "lift": ([], [])}; skipped = {"clinical": 0, "lift": 0}
    for _, g in m.groupby(["repeat_id", "fold"]):
        for game, col in [("clinical", "c_clin"), ("lift", "delta_c")]:
            v = values(g, col, n)
            if game == "clinical": v[0] = 0.5
            if np.isnan(v).any(): skipped[game] += 1; continue
            out[game][0].append(shapley(v, n)[0]); out[game][1].append(interactions(v, n))
    total = m.groupby(["repeat_id", "fold"]).ngroups
    return out, skipped, total, test_train

def values(df, col, n):
    v = np.full(2 ** n, np.nan); v[df["mask"].values] = df[col].values; return v

def collect(results_dir):
    rows, curves, tnm, pairs, add, boxes = [], {}, [], [], [], []
    for code in CODES:
        s, p, full = load(results_dir, code); n = len(full)
        anchor = s[s["mask"] == 0].iloc[0]
        for game, col in [("clinical", "c_ref"), ("lift", "delta_c")]:
            v = values(s, col, n)
            if game == "clinical": v[0] = 0.5  #the anchor's clinical model is Kaplan-Meier
            if np.isnan(v).any(): print(f"[{code}] {game}: {np.isnan(v).sum()} subsets without a value, skipped"); continue
            phi, banz = shapley(v, n)
            reps = []  #the same values within each CV repeat
            for r, g in p.groupby("repeat_id"):
                vr = values(g, col, n)
                if game == "clinical": vr[0] = 0.5
                if not np.isnan(vr).any(): reps.append(shapley(vr, n)[0])
            reps = np.array(reps) if reps else np.full((1, n), np.nan)
            total = v[-1] - v[0]  #what the values add up to
            I = interactions(v, n)
            Ir = [interactions(values(g, col, n) if game == "lift" else np.r_[0.5, values(g, col, n)[1:]], n)
                  for _, g in p.groupby("repeat_id")]  #per repeat, for the range
            Ir = np.array([x for x in Ir if np.isfinite(x).all()]) if Ir else np.full((1, n, n), np.nan)
            for i in range(n):
                for j in range(i + 1, n):
                    pairs.append(dict(cancer=code, game=game, var_a=full[i], var_b=full[j], interaction=I[i, j],
                                      rep_min=Ir[:, i, j].min(), rep_max=Ir[:, i, j].max(), n_repeats=len(Ir)))
            r_main, r_pair = additivity(v, n)
            #noise ceiling: split-half reliability of the lattice values (odd vs even repeats, Spearman-Brown corrected)
            pr = p.pivot_table(index="mask", columns="repeat_id", values=col)
            if game == "clinical": pr.loc[0] = 0.5
            pr = pr.dropna(); odd = pr[[c for c in pr.columns if c % 2 == 1]].mean(1); even = pr[[c for c in pr.columns if c % 2 == 0]].mean(1)
            rr = np.corrcoef(odd, even)[0, 1] if len(pr) > 2 and odd.std() > 0 and even.std() > 0 else np.nan
            ceiling = max(0.0, 2 * rr / (1 + rr)) if np.isfinite(rr) else np.nan
            add.append(dict(cancer=code, game=game, n_vars=n, r2_main=r_main, r2_main_pairs=r_pair, noise_ceiling=ceiling, v_full=v[-1]))
            for i, var in enumerate(full):
                rows.append(dict(cancer=code, game=game, variable=var, shapley=phi[i], banzhaf=banz[i],
                                 share=phi[i] / total if game == "clinical" else np.nan,
                                 rep_min=reps[:, i].min(), rep_max=reps[:, i].max(), rep_sd=reps[:, i].std(ddof=1) if len(reps) > 1 else np.nan,
                                 n_repeats=len(reps), total=total, v_full=v[-1], v_empty=v[0], n_vars=n))
        curves[code] = (s, full, anchor)
        fl, skipped, total, tt = fold_level(results_dir, code, full)
        for game in ["clinical", "lift"]:
            phis, ints = (np.array(x) for x in fl[game])
            usable = len(phis) >= total / 2  #fewer than half of the folds usable -> can't tell
            FOLD[(code, game)] = (phis, full, tt, usable)
            for i, var in enumerate(full):
                mm, lo, hi, box = p1_box(phis[:, i], tt) if len(phis) else (np.nan,) * 3 + ("cant_tell",)
                boxes.append(dict(cancer=code, game=game, kind="variable", var_a=var, var_b="", fold_mean=mm, ci_lo=lo, ci_hi=hi,
                                  box=box if usable else "cant_tell", n_folds=len(phis), n_folds_total=total))
            for i in range(n):
                for j in range(i + 1, n):
                    mm, lo, hi, box = p1_box(ints[:, i, j], tt) if len(ints) else (np.nan,) * 3 + ("cant_tell",)
                    boxes.append(dict(cancer=code, game=game, kind="pair", var_a=full[i], var_b=full[j], fold_mean=mm, ci_lo=lo, ci_hi=hi,
                                      box=box if usable else "cant_tell", n_folds=len(ints), n_folds_total=total))
        if skipped["lift"]: print(f"[{code}] lift: {skipped['lift']} of {total} folds left out (failed cells)")
        #T, N and M added to a subset with vs without stage (clinical C)
        if "stage" in full:
            v = values(s, "c_ref", n); v[0] = 0.5; st = 1 << full.index("stage"); m = np.arange(2 ** n)
            for var in [x for x in ["T", "N", "M"] if x in full]:
                b = 1 << full.index(var)
                for mm in m[(m & b) == 0]:
                    tnm.append(dict(cancer=code, variable=var, with_stage=bool(mm & st), gain=v[mm | b] - v[mm]))
    return pd.DataFrame(rows), curves, pd.DataFrame(tnm), pd.DataFrame(pairs), pd.DataFrame(add), pd.DataFrame(boxes)

#shared layout for the main figure: cancers by full-set clinical C, shared variables by median clinical share
def layout(R):
    clin = R[R.game == "clinical"]
    cnt = clin.groupby("variable").cancer.nunique()
    shared = cnt[cnt >= MIN_CANCERS].index.tolist()
    order_c = clin.drop_duplicates("cancer").sort_values("v_full", ascending=False).cancer.tolist()  #strongest full baseline first
    order_v = clin[clin.variable.isin(shared)].groupby("variable").share.median().sort_values(ascending=False).index.tolist()
    return shared, order_c, order_v

#one Shapley grid: variables x cancers, cancer-specific variables summed in the last row, median on the right
def draw_grid(fig, ax, R, game, title, fmt, unit, shared, order_c, order_v, Bv=None):
    rows_lab = order_v + ["cancer-specific (sum)"]
    d = R[R.game == game]
    M = np.full((len(rows_lab), len(order_c)), np.nan)
    for j, c in enumerate(order_c):
        dc = d[d.cancer == c]
        for i, v in enumerate(order_v):
            x = dc[dc.variable == v][unit]
            if len(x): M[i, j] = x.iloc[0]
        spec = dc[~dc.variable.isin(shared)][unit]
        if len(spec): M[-1, j] = spec.sum()
    lim = np.nanmax(np.abs(M)); ax.set_facecolor(SURF)
    if game == "clinical": lim = min(lim, 1.0)  #shares can pass 100% when other variables are negative; colour is capped, the numbers are not
    norm = TwoSlopeNorm(0, -lim, lim)
    box = {}  #(row, column) -> box of the fold-level interval
    if Bv is not None:
        bb = Bv[Bv.game == game]; bmap = dict(zip(zip(bb.var_a, bb.cancer), bb.box))
        for i, v in enumerate(order_v):
            for j, c in enumerate(order_c): box[(i, j)] = bmap.get((v, c))
        for j, c in enumerate(order_c):  #summed cancer-specific row: the sum per fold gets its own interval
            if (c, game) in FOLD:
                phis, full, tt, usable = FOLD[(c, game)]; idx = [k for k, v in enumerate(full) if v not in shared]
                if idx and len(phis): box[(len(order_v), j)] = p1_box(phis[:, idx].sum(axis=1), tt)[3] if usable else "cant_tell"
    RGBA = np.zeros(M.shape + (4,))
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            RGBA[i, j] = (0, 0, 0, 0) if np.isnan(M[i, j]) else vsup(INTERACT(norm(M[i, j])), box.get((i, j), "above"))
    ax.imshow(RGBA, aspect="auto")
    im = plt.cm.ScalarMappable(norm=norm, cmap=INTERACT)  #colour bar for the clear cells
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if np.isnan(M[i, j]): continue
            strong = box.get((i, j), "above") in ("above", "below")
            dark = strong and min(abs(M[i, j]) / lim, 1) > 0.55
            ax.text(j, i, fmt(M[i, j]), ha="center", va="center", fontsize=6.5, color="white" if dark else (INK if strong else INK2),
                    fontweight="bold" if strong else "normal")
    for i in range(M.shape[0]):  #pan-cancer median and the number of cancers with the variable
        if np.isfinite(M[i]).any():
            ax.text(len(order_c) - 0.3, i, f"{fmt(np.nanmedian(M[i]))}  (n={int(np.isfinite(M[i]).sum())})", va="center", fontsize=7, color=INK, ha="left")
    ax.text(len(order_c) - 0.3, -0.9, "median", fontsize=7.5, color=INK2, ha="left")
    ax.axhline(len(order_v) - 0.5, color=INK2, lw=0.8)
    ax.set_xticks(range(len(order_c))); ax.set_xticklabels(order_c, rotation=90, fontsize=8, color=INK)
    ax.set_yticks(range(len(rows_lab))); ax.set_yticklabels(rows_lab, fontsize=8.5, color=INK)
    ax.set_xlim(-0.5, len(order_c) + 2.2)
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.tick_params(length=0); ax.set_title(title, fontsize=10, color=INK, loc="left")
    cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=7, colors=INK2)

#supplementary: per cancer, every variable's fold mean and 90% interval
def dot_figure(R, B, game, xlab, out_png):
    d = R[R.game == game]; codes = d.drop_duplicates("cancer").sort_values("v_full", ascending=False).cancer.tolist()
    ncol = 4; nrow = int(np.ceil(len(codes) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(17, 3.1 * nrow), dpi=150); fig.patch.set_facecolor(SURF)
    for ax in axes.flat[len(codes):]: ax.axis("off")
    for ax, c in zip(axes.flat, codes):
        x = d[d.cancer == c].sort_values("shapley")
        b = B[(B.kind == "variable") & (B.game == game) & (B.cancer == c)].set_index("var_a").loc[x.variable]
        y = np.arange(len(x)); style(ax)
        ax.axvspan(-MARGIN, MARGIN, color=GREY_NEGL, zorder=0, lw=0)  #the +-0.01 margin
        ax.axvline(0, color=MUTED, lw=0.8)
        col = [{"above": "#c98500", "below": "#4a3aa7"}.get(bx, "#9a988f") for bx in b.box]  #gold/violet = clear, grey = not
        ax.hlines(y, b.ci_lo, b.ci_hi, color=col, lw=2.2, alpha=0.55)
        ax.scatter(b.fold_mean, y, s=26, color=col, zorder=3)
        ax.set_yticks(y); ax.set_yticklabels(x.variable, fontsize=7.5, color=INK)
        t = x.total.iloc[0]
        ax.set_title(f"{c}  ({len(x)} vars, sum = {t:+.3f})", fontsize=9, color=INK, loc="left")
        ax.grid(axis="x", color=GRID, lw=.5); ax.set_axisbelow(True)
        ax.xaxis.set_major_locator(plt.MaxNLocator(4))  #fewer ticks, no overlap
    fig.supxlabel(xlab, fontsize=9, color=INK2)
    plt.tight_layout(rect=(0, 0.02, 1, 1)); plt.savefig(out_png, facecolor=SURF); plt.close(fig)

#supplementary, one per cancer: specification curve of RNA's lift with the variables of each subset underneath
def spec_curve(code, s, full, anchor, out_png):
    d = s[s["mask"] != 0].sort_values("delta_c").reset_index(drop=True); x = np.arange(len(d))
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(12, 3.2 + 0.24 * len(full)), dpi=150, sharex=True,
                                 gridspec_kw=dict(height_ratios=[2.2, 0.28 * len(full) + 0.4], hspace=0.05))
    fig.patch.set_facecolor(SURF); style(a1); style(a2)
    sz = 14 if len(d) < 300 else (6 if len(d) < 1500 else 2.5)
    a1.vlines(x, d.delta_c - d.delta_c_sd, d.delta_c + d.delta_c_sd, color="#e1e0d9", lw=0.6 if len(d) < 1500 else 0.25, zorder=1)
    sc = a1.scatter(x, d.delta_c, c=d.c_ref, cmap=SEQ, s=sz, zorder=3, linewidths=0)
    a1.axhline(0, color=MUTED, lw=0.9)
    a1.axhline(anchor.delta_c, color=INK2, lw=1, ls="--")
    a1.annotate("RNA only (no clinical variables)", (0, anchor.delta_c), xytext=(4, 4), textcoords="offset points", fontsize=7.5, color=INK2,
                bbox=dict(facecolor=SURF, edgecolor="none", pad=0.5))  #label sits above the line
    a1.set_ylabel("ΔC (clinical + RNA − clinical)", fontsize=9, color=INK)
    a1.set_title(f"{code}: RNA's lift over all {len(d)} clinical baselines, ranked", fontsize=10.5, color=INK, loc="left")
    cb = fig.colorbar(sc, ax=[a1, a2], fraction=0.018, pad=0.01); cb.set_label("clinical-only C of the baseline", fontsize=8, color=INK2)
    cb.outline.set_visible(False); cb.ax.tick_params(labelsize=7, colors=INK2)
    M = np.array([[(m >> i) & 1 for m in d["mask"]] for i in range(len(full))])
    a2.imshow(M, aspect="auto", cmap=LinearSegmentedColormap.from_list("b", [SURF, "#262522"]), interpolation="nearest")
    a2.set_yticks(range(len(full))); a2.set_yticklabels(full, fontsize=7.5, color=INK)
    a2.set_xlabel("Clinical baselines, ranked by ΔC", fontsize=9, color=INK2); a2.tick_params(length=0)
    plt.savefig(out_png, facecolor=SURF, bbox_inches="tight"); plt.close(fig)

#supplementary: what T, N and M add to the clinical C, in subsets without vs with stage
def tnm_figure(T, out_png):
    codes = T.cancer.unique().tolist(); ncol = 5; nrow = int(np.ceil(len(codes) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(17, 2.9 * nrow), dpi=150); fig.patch.set_facecolor(SURF)
    for ax in axes.flat[len(codes):]: ax.axis("off")
    for ax, c in zip(axes.flat, codes):
        d = T[T.cancer == c]; vs = [v for v in ["T", "N", "M"] if v in d.variable.values]; style(ax)
        for k, v in enumerate(vs):
            for j, (ws, col) in enumerate([(False, "#9cbfe9"), (True, "#1c5aa6")]):
                g = d[(d.variable == v) & (d.with_stage == ws)].gain.values
                ax.boxplot(g, positions=[k + (j - 0.5) * 0.36], widths=0.3, patch_artist=True, showfliers=False,
                           boxprops=dict(facecolor=col, edgecolor=col), medianprops=dict(color=INK, lw=1.2),
                           whiskerprops=dict(color=col), capprops=dict(color=col))
        ax.axhline(0, color=MUTED, lw=0.8); ax.set_xticks(range(len(vs))); ax.set_xticklabels(vs, color=INK)
        ax.set_title(c, fontsize=9, color=INK, loc="left"); ax.grid(axis="y", color=GRID, lw=.5); ax.set_axisbelow(True)
    h = [plt.Rectangle((0, 0), 1, 1, color="#9cbfe9"), plt.Rectangle((0, 0), 1, 1, color="#1c5aa6")]
    fig.legend(h, ["added to a subset without stage", "added to a subset with stage"], loc="upper center", ncol=2, frameon=False, fontsize=9)
    fig.supylabel("Gain in clinical-only C from adding the variable", fontsize=9, color=INK2)
    plt.tight_layout(rect=(0, 0, 1, 0.96)); plt.savefig(out_png, facecolor=SURF); plt.close(fig)

#additivity per cancer: share of the lattice's variation explained by single variables, and with every pair added
def draw_additivity(ax, A, order_c):
    style(ax); y = np.arange(len(order_c))
    for game, col, dy in [("clinical", "#2a78d6", -0.17), ("lift", "#eb6834", 0.17)]:
        d = A[A.game == game].set_index("cancer").loc[order_c]
        ax.hlines(y + dy, d.r2_main, d.r2_main_pairs, color=col, lw=1.6, alpha=0.6)
        ax.scatter(d.r2_main, y + dy, s=30, color=col, zorder=3)
        ax.scatter(d.r2_main_pairs, y + dy, s=30, facecolor=SURF, edgecolor=col, lw=1.6, zorder=3)
        ax.scatter(d.noise_ceiling, y + dy, marker="|", s=110, color=INK, lw=1.6, zorder=4)  #the most any model could explain
    ax.set_yticks(y); ax.set_yticklabels(order_c, fontsize=8, color=INK); ax.invert_yaxis(); ax.set_xlim(0, 1.02)
    ax.set_xlabel("Share of the variation across all subsets explained", fontsize=9, color=INK2)
    ax.grid(axis="x", color=GRID, lw=.5); ax.set_axisbelow(True)
    hs = [plt.Line2D([], [], color="#2a78d6", marker="o", lw=0), plt.Line2D([], [], color="#eb6834", marker="o", lw=0),
          plt.Line2D([], [], color=INK2, marker="o", lw=0), plt.Line2D([], [], color=INK2, marker="o", lw=0, mfc=SURF),
          plt.Line2D([], [], color=INK, marker="|", lw=0, markersize=10, mew=1.6)]
    ax.legend(hs, ["clinical-only C", "RNA's lift (ΔC)", "variables one by one", "plus every pair", "noise ceiling"], frameon=False, fontsize=8, loc="lower left")
    ax.set_title("How additive is each cancer's lattice?", fontsize=10, color=INK, loc="left")

#pairs across cancers: share of cancers in each box; RNA's lift is sign-flipped so that overlap means the same in both triangles
def draw_triangle(fig, ax, P, shared_order):
    #P: fold-level pair table, one row per cancer, game and pair; each cell is a bar split into overlap | additive | works together | can't tell
    k = len(shared_order); ax.set_facecolor(SURF)
    seg_col = [INTERACT(0.08), ADDITIVE, INTERACT(0.92), GREY_NEGL]
    for game, upper, sign in [("clinical", True, 1), ("lift", False, -1)]:
        d = P[P.game == game]
        for i in range(k):
            for j in range(i + 1, k):
                a, b = shared_order[i], shared_order[j]
                x = d[((d.var_a == a) & (d.var_b == b)) | ((d.var_a == b) & (d.var_b == a))]
                if len(x) < 2: continue
                over_lab, syn_lab = ("below", "above") if sign == 1 else ("above", "below")  #lift flipped: overlap = less taken together
                parts = [(x.box == over_lab).sum(), (x.box == "negligible").sum(), (x.box == syn_lab).sum(), (x.box == "cant_tell").sum()]
                r, c = (i, j) if upper else (j, i)
                left = c - 0.43
                for n_part, col in zip(parts, seg_col):
                    w = 0.86 * n_part / len(x)
                    if w > 0: ax.add_patch(plt.Rectangle((left, r - 0.32), w, 0.36, facecolor=col, edgecolor=SURF, lw=0.6))
                    left += w
                ax.text(c, r + 0.27, f"of {len(x)}", ha="center", va="center", fontsize=6, color=INK2)
    ax.plot([-0.5, k - 0.5], [-0.5, k - 0.5], color=INK2, lw=0.8)
    ax.set_xlim(-0.5, k - 0.5); ax.set_ylim(k - 0.5, -0.5); ax.set_aspect("equal")
    ax.set_xticks(range(k)); ax.set_xticklabels(shared_order, rotation=90, fontsize=8, color=INK)
    ax.set_yticks(range(k)); ax.set_yticklabels(shared_order, fontsize=8, color=INK)
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.tick_params(length=0)
    ax.set_title("Pairs across cancers (90% corrected interval vs ±0.01)\n"
                 "upper right: clinical-only C   ·   lower left: what the pair takes from RNA's lift", fontsize=10, color=INK, loc="left")
    hs = [plt.Rectangle((0, 0), 1, 1, facecolor=c, edgecolor="none") for c in seg_col]
    ax.legend(hs, ["overlap", "additive (within ±0.01)", "work together", "can't tell"], loc="upper center", bbox_to_anchor=(0.5, -0.12),
              ncol=4, frameon=False, fontsize=8.5)

#Figure 4 A-B: Shapley values (clinical share | change in RNA's lift) and additivity (per cancer | pairwise interactions)
def main_figure(R, P, A, B, out_png):
    shared, order_c, order_v = layout(R)
    fig = plt.figure(figsize=(22, 16.5), dpi=150); fig.patch.set_facecolor(SURF)
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], width_ratios=[1, 1], hspace=0.28, wspace=0.18,
                          left=0.07, right=0.97, top=0.95, bottom=0.07)
    a1, a2 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    b1, b2 = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
    Bv, Bp = B[B.kind == "variable"], B[B.kind == "pair"]
    draw_grid(fig, a1, R, "clinical", "Share of the full clinical baseline (% of C above 0.5)", lambda x: f"{100 * x:.0f}", "share", shared, order_c, order_v, Bv)
    draw_grid(fig, a2, R, "lift", "Change in RNA's lift (ΔC × 100; negative = takes over RNA's lift)", lambda x: f"{100 * x:.1f}", "shapley", shared, order_c, order_v, Bv)
    draw_additivity(b1, A, order_c)
    draw_triangle(fig, b2, Bp, order_v)
    for ax, letter in [(a1, "A"), (b1, "B")]:
        ax.text(-0.09, 1.04, letter, transform=ax.transAxes, fontsize=18, fontweight="bold", color=INK)
    fig.text(0.07, 0.975, "Shapley values", fontsize=12, color=INK, fontweight="bold")
    fig.text(0.07, 0.505, "Are the variables additive?", fontsize=12, color=INK, fontweight="bold")
    plt.savefig(out_png, facecolor=SURF); plt.close(fig)

#supplementary: every cancer's interaction matrix, upper triangle clinical C, lower triangle RNA's lift
def interaction_by_cancer(R, P, B, out_png):
    codes = R[R.game == "clinical"].drop_duplicates("cancer").sort_values("v_full", ascending=False).cancer.tolist()
    lim = P.interaction.abs().quantile(0.99) * 100
    ncol = 4; nrow = int(np.ceil(len(codes) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(17, 4.3 * nrow), dpi=150); fig.patch.set_facecolor(SURF)
    for ax in axes.flat[len(codes):]: ax.axis("off")
    for ax, c in zip(axes.flat, codes):
        vs = R[(R.game == "clinical") & (R.cancer == c)].sort_values("shapley", ascending=False).variable.tolist(); k = len(vs)
        M = np.full((k, k), np.nan); BX = {}
        bp = B[(B.kind == "pair") & (B.cancer == c)]; bmap = dict(zip(zip(bp.game, bp.var_a, bp.var_b), bp.box))
        for game, upper, sign in [("clinical", True, 1), ("lift", False, -1)]:  #lift flipped: negative = overlap in both triangles
            d = P[(P.game == game) & (P.cancer == c)]
            for row in d.itertuples():
                i, j = sorted([vs.index(row.var_a), vs.index(row.var_b)]); rc = (i, j) if upper else (j, i)
                M[rc] = sign * row.interaction * 100; BX[rc] = bmap.get((game, row.var_a, row.var_b), "cant_tell")
        ax.set_facecolor(SURF); norm = TwoSlopeNorm(0, -lim, lim)
        RGBA = np.zeros((k, k, 4))
        for (i, j), bx in BX.items(): RGBA[i, j] = vsup(INTERACT(norm(M[i, j])), bx, pair=True)  #faded = can't tell, teal = additive
        ax.imshow(RGBA); im = plt.cm.ScalarMappable(norm=norm, cmap=INTERACT)
        if k <= 9:
            for (i, j), bx in BX.items():
                strong = bx in ("above", "below")
                ax.text(j, i, f"{M[i, j]:.1f}", ha="center", va="center", fontsize=5.5, fontweight="bold" if strong else "normal",
                        color="white" if strong and abs(M[i, j]) / lim > 0.55 else (INK if strong else INK2))
        ax.plot([-0.5, k - 0.5], [-0.5, k - 0.5], color=INK2, lw=0.6)
        ax.set_xticks(range(k)); ax.set_xticklabels(vs, rotation=90, fontsize=6.5, color=INK)
        ax.set_yticks(range(k)); ax.set_yticklabels(vs, fontsize=6.5, color=INK)
        for sp in ax.spines.values(): sp.set_visible(False)
        ax.tick_params(length=0); ax.set_title(c, fontsize=9, color=INK, loc="left")
    cb = fig.colorbar(im, ax=axes, fraction=0.015, pad=0.01); cb.outline.set_visible(False)
    cb.set_label("Shapley interaction × 100 (violet = overlap, gold = works together; faded = can't tell, teal = additive)\nupper: clinical-only C · lower: what the pair takes from RNA's lift (sign flipped)", fontsize=8, color=INK2)
    plt.savefig(out_png, facecolor=SURF, bbox_inches="tight"); plt.close(fig)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--out-dir", default="figures/variable_importance")
    a = ap.parse_args()
    os.makedirs(os.path.join(a.out_dir, "spec_curves"), exist_ok=True)
    R, curves, T, P, A, B = collect(a.results_dir)
    B.to_csv(os.path.join(a.out_dir, "fold_level_boxes.csv"), index=False)
    R.to_csv(os.path.join(a.out_dir, "shapley_values.csv"), index=False)
    P.to_csv(os.path.join(a.out_dir, "shapley_interactions.csv"), index=False)
    A.to_csv(os.path.join(a.out_dir, "lattice_additivity.csv"), index=False)
    T.to_csv(os.path.join(a.out_dir, "tnm_with_without_stage.csv"), index=False)
    main_figure(R, P, A, B, os.path.join(a.out_dir, "fig_variable_importance_main.png"))
    dot_figure(R, B, "clinical", "Shapley value for the clinical-only C (fold mean, 90% corrected interval; shaded = ±0.01)", os.path.join(a.out_dir, "fig_shapley_clinical_by_cancer.png"))
    dot_figure(R, B, "lift", "Shapley value for RNA's lift, ΔC (fold mean, 90% corrected interval; shaded = ±0.01)", os.path.join(a.out_dir, "fig_shapley_lift_by_cancer.png"))
    for code, (s, full, anchor) in curves.items():
        spec_curve(code, s, full, anchor, os.path.join(a.out_dir, "spec_curves", f"{code}_spec_curve.png"))
    tnm_figure(T, os.path.join(a.out_dir, "fig_tnm_with_without_stage.png"))
    interaction_by_cancer(R, P, B, os.path.join(a.out_dir, "fig_interactions_by_cancer.png"))
    print(R.pivot_table(index="variable", columns="game", values=["share", "shapley"], aggfunc="median").round(3).to_string())
