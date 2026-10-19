#Figure 6: the 9 external cohorts
#A clinical-only C against clinical + RNA C for every clinical subset, as in Figure 5
#B correlation of RNA's lift between each external cohort and its TCGA cancer, over the subsets present in both
#External cells where IPF-LASSO failed are dropped, as in TCGA; there is no rescue run for the external cohorts
#Usage, from the project root: python figure_generation/fig6.py --results-dir results --out-dir figures
import argparse, os
import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
matplotlib.rcParams["font.family"] = ["Inter", "sans-serif", "DejaVu Sans"]  #font of the paper's figures; DejaVu Sans if Inter is not installed
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from scipy.stats import spearmanr, pearsonr

ap = argparse.ArgumentParser()
ap.add_argument("--results-dir", default="results"); ap.add_argument("--out-dir", default="figures")
args = ap.parse_args(); os.makedirs(args.out_dir, exist_ok=True)
R = os.path.join(args.results_dir, "")
OUT = os.path.join(args.out_dir, "")
PAIRS = [('GSE68465', 'LUAD'), ('GSE39582', 'COAD'), ('GSE84437', 'STAD'), ('GSE65858', 'HNSC'), ('GSE14520', 'LIHC'),
         ('GSE13507', 'BLCA'), ('BEATAML', 'LAML'), ('CGGA_LGG', 'LGG'), ('CGGA_GBM', 'GBM')]  #external cohort, TCGA cancer
DROP = {'GSE13507': 'invasiveness'}  #variable not in the cohort's final set: subsets containing it are left out
DELTA = 0.01  #fixed equivalence margin
B = 1000; rng = np.random.default_rng(20261011)
J, O, G = '#c0304f', '#8e5bb5', '#1a8a3c'

def crossing(a, b, lo, hi):  #x where the line meets y = x, NaN if outside [lo, hi] or parallel
    if abs(1 - a) < 1e-9: return np.nan
    x = b / (1 - a)
    return x if lo <= x <= hi else np.nan

def metrics(x, y):  #panel numbers for one set of subset means
    a, b = np.polyfit(x, y, 1); d = y - x
    return dict(slope=a, cross=crossing(a, b, x.min(), x.max()), loss=np.abs(d).mean(),
                within=(np.abs(d) <= DELTA).mean() * 100, better=(d > 0).mean() * 100)

def load_external(code):  #subset x repeat matrices (fold means per repeat) of clinical, IPF-LASSO and offset C
    d = pd.read_csv(os.path.join(R, 'external', code + '_results.csv'), low_memory=False)
    d['is_anchor'] = d.is_anchor.astype(str).str.upper().eq('TRUE')
    if code in DROP:  #keep only subsets without the dropped variable
        d = d[~d.variables.astype(str).str.split('+').apply(lambda v: DROP[code] in v)]
    d = d[d.c_fused.notna()]  #bad cells: IPF-LASSO failed
    g = d.groupby(['subset_id', 'repeat_id'])
    X = g.c_clin.mean().unstack(); Y = g.c_fused.mean().unstack(); Off = g.c_offset.mean().unstack()
    info = d.groupby('subset_id').agg(n_vars=('n_vars', 'first'), variables=('variables', 'first'), anchor=('is_anchor', 'first'))
    return X, Y, Off, info

rows = []; panels = {}; lifts = {}
for ext, tc in PAIRS:
    X, Y, Off, info = load_external(ext)
    a_ids = info.index[info.anchor]; s_ids = info.index[~info.anchor]
    anc = (X.loc[a_ids].mean(axis=1).mean(), Y.loc[a_ids].mean(axis=1).mean()) if len(a_ids) else (np.nan, np.nan)
    Xs, Ys, Os = X.loc[s_ids], Y.loc[s_ids], Off.loc[s_ids]
    x, y, off = Xs.mean(axis=1).values, Ys.mean(axis=1).values, Os.mean(axis=1).values
    m = metrics(x, y)
    se2 = np.nanvar(Xs.values, 1, ddof=1) / np.sum(np.isfinite(Xs.values), 1)  #leftover noise in each subset's clinical C
    rel = 1 - np.nanmean(se2) / np.var(x, ddof=1)  #reliability of the clinical-C spread
    Rn = Xs.shape[1]; boot = []
    for _ in range(B):  #resample repeats with replacement
        idx = rng.integers(0, Rn, Rn)
        xb, yb = np.nanmean(Xs.values[:, idx], 1), np.nanmean(Ys.values[:, idx], 1)
        ok = np.isfinite(xb) & np.isfinite(yb); boot.append(metrics(xb[ok], yb[ok]))
    bt = pd.DataFrame(boot); q = lambda k: np.nanpercentile(bt[k], [2.5, 97.5]); cr = np.isfinite(bt.cross).mean() * 100
    po = np.polyfit(x, off, 1)
    r = dict(cohort=ext, tcga=tc, n_subsets=len(x), repeats=Rn, pct_better=m['better'], loss=m['loss'], loss_lo=q('loss')[0],
             loss_hi=q('loss')[1], pct_within=m['within'], slope=m['slope'], slope_lo=q('slope')[0], slope_hi=q('slope')[1],
             reliability=rel, slope_noise_checked=m['slope'] / rel, cross=m['cross'],
             cross_lo=q('cross')[0] if cr > 0 else np.nan, cross_hi=q('cross')[1] if cr > 0 else np.nan, pct_resamples_crossing=cr,
             slope_offset=po[0])
    panels[ext] = (x, y, off, info.loc[s_ids, 'n_vars'].values, anc, r)
    #panel B: lift per subset in both, matched on the set of variable names
    ext_lift = pd.Series(y - x, index=[frozenset(v.split('+')) for v in info.loc[s_ids, 'variables']])
    t = pd.read_csv(R + tc + '_subset_summary.csv'); t = t[(t.arm == 'fused') & ~t.is_anchor.astype(str).str.upper().eq('TRUE')]
    tcga_lift = pd.Series(t.delta_c.values, index=[frozenset(v.split('+')) for v in t.variables])
    common = [k for k in ext_lift.index if k in tcga_lift.index]
    a_, b_ = tcga_lift.loc[common].values, ext_lift.loc[common].values; ok = np.isfinite(a_) & np.isfinite(b_)
    a_, b_ = a_[ok], b_[ok]; nv = np.array([len(k) for k, o in zip(common, ok) if o])
    sp = spearmanr(a_, b_)[0] if len(a_) > 2 else np.nan; pe = pearsonr(a_, b_)[0] if len(a_) > 2 else np.nan
    shared = sorted(set().union(*common)) if common else []
    r.update(matched_subsets=len(a_), shared_vars='+'.join(shared), spearman_tcga_vs_ext=sp, pearson_tcga_vs_ext=pe)
    lifts[ext] = (a_, b_, nv, sp, pe, shared)
    rows.append(r)

T = pd.DataFrame(rows); T.round(4).to_csv(OUT + 'fig6_metrics.csv', index=False)
vmax = max(p[3].max() for p in panels.values()); cmap = plt.cm.viridis
fig = plt.figure(figsize=(14.5, 16.5))
gs = fig.add_gridspec(4, 4, height_ratios=[1, 1, 1, 0.75], hspace=0.42, wspace=0.32)
#Panel A
for k, (ext, tc) in enumerate(PAIRS):
    ax = fig.add_subplot(gs[k // 3, k % 3])
    x, y, off, nv, anc, r = panels[ext]
    vals = np.r_[x, y, off[np.isfinite(off)], [v for v in anc if np.isfinite(v)]]
    pad = 0.03 * (np.nanmax(vals) - np.nanmin(vals)); lo, hi = np.nanmin(vals) - pad, np.nanmax(vals) + pad
    ax.fill_between([lo, hi], [lo, hi], hi, color='#e8f4ee', zorder=0)
    ax.plot([lo, hi], [lo, hi], color='#888', lw=0.9, ls=(0, (4, 3)), zorder=1)
    sz = 5 if len(x) > 500 else 12
    ax.scatter(x, off, s=sz * 0.6, marker='s', color='#bbb', alpha=0.5, lw=0, zorder=2)
    ax.scatter(x, y, s=sz, c=nv, cmap=cmap, vmin=0, vmax=vmax, alpha=0.85, lw=0, zorder=3)
    if np.all(np.isfinite(anc)): ax.scatter([anc[0]], [anc[1]], marker='D', s=28, color='#e07b22', zorder=4)
    xs = np.array([x.min(), x.max()]); pj = np.polyfit(x, y, 1); po = np.polyfit(x, off, 1)
    ax.plot(xs, np.polyval(pj, xs), color=J, lw=1.5, zorder=5); ax.plot(xs, np.polyval(po, xs), color=O, lw=1.5, ls='--', zorder=5)
    ax.plot(xs, y.mean() + r['slope_noise_checked'] * (xs - x.mean()), color=G, lw=1.4, ls=(0, (1.5, 1.5)), zorder=5)
    if np.isfinite(r['cross']):
        ax.scatter([r['cross']], [r['cross']], marker='x', s=40, color=J, lw=1.6, zorder=6)
        ax.plot([r['cross_lo'], r['cross_hi']], [r['cross_lo'], r['cross_hi']], color=J, lw=3, alpha=0.35, zorder=5)
        c1, c2 = f"crosses at {r['cross']:.3f} [{r['cross_lo']:.3f}–{r['cross_hi']:.3f}]", f"in {r['pct_resamples_crossing']:.0f}% of resamples"
    else:
        c1, c2 = 'no crossing', f"({r['pct_resamples_crossing']:.0f}% of resamples cross)"
    lines = [(f"offset slope {r['slope_offset']:.2f}", O), (f"noise-checked slope {r['slope_noise_checked']:.2f}", G),
             (f"slope {r['slope']:.2f} [{r['slope_lo']:.2f}–{r['slope_hi']:.2f}]", J),
             (f"{r['pct_better']:.0f}% above diagonal", '#2a9d6f' if r['pct_better'] >= 50 else '#c0304f'),
             (f"{r['pct_within']:.0f}% within ±0.01", '#333'),
             (f"mean |ΔC| {r['loss']:.3f} [{r['loss_lo']:.3f}–{r['loss_hi']:.3f}]", '#333'), (c1, J), (c2, J)]
    n = len(lines)
    for i, (t, col) in enumerate(lines):
        ax.text(0.97, 0.03 + (n - 1 - i) * 0.052, t, transform=ax.transAxes, fontsize=5.6, color=col, ha='right', va='bottom')
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect('equal')
    ax.set_title(f'{ext} ({tc}; {int(nv.max())} vars, {r["repeats"]} repeats)', fontsize=8.5, pad=3)
    ax.set_xlabel("clinical-only model (Uno's C)", fontsize=7); ax.set_ylabel("clinical + omics (Uno's C)", fontsize=7)
    ax.tick_params(labelsize=6); ax.xaxis.set_major_locator(MaxNLocator(5)); ax.yaxis.set_major_locator(MaxNLocator(6))
    if k == 0: ax.text(-0.28, 1.12, 'A', transform=ax.transAxes, fontsize=16, fontweight='bold')
lg = fig.add_subplot(gs[0:3, 3]); lg.axis('off')
h = [Line2D([], [], marker='o', ls='', color='#666', label='subset, joint penalised fit (IPF-LASSO)'),
     Line2D([], [], marker='s', ls='', color='#bbb', label='subset, frozen-offset fit'),
     Line2D([], [], marker='D', ls='', color='#e07b22', label='omics-only anchor'),
     Line2D([], [], color='#888', ls=(0, (4, 3)), label='identity: omics adds nothing'),
     Line2D([], [], color=J, label='OLS fit, joint penalised'), Line2D([], [], color=O, ls='--', label='OLS fit, frozen-offset'),
     Line2D([], [], color=G, ls=(0, (1.5, 1.5)), lw=1.4, label='noise-checked fit, joint penalised'),
     Line2D([], [], marker='x', ls='', color=J, label='joint fit crosses the identity line'),
     Line2D([], [], color=J, lw=3, alpha=0.35, label='95% interval of the crossing')]
lg.legend(handles=h, loc='upper left', frameon=False, fontsize=7.5, title='Panel A', title_fontsize=8.5)
cax = lg.inset_axes([0.05, 0.42, 0.9, 0.025])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0, vmax), cmap=cmap), cax=cax, orientation='horizontal')
cb.set_label('number of clinical variables in the baseline', fontsize=7.5); cb.ax.tick_params(labelsize=6.5)
#Panel B: one row per TCGA/external pair
ax = fig.add_subplot(gs[3, 0:3]); ax.axis('off')
ax.text(-0.09, 1.02, 'B', transform=ax.transAxes, fontsize=16, fontweight='bold')
cells = [[tc, ext, ', '.join(lifts[ext][5]), str(len(lifts[ext][0])), f'{lifts[ext][3]:.2f}', f'{lifts[ext][4]:.2f}'] for ext, tc in PAIRS]
tb = ax.table(cellText=cells, colLabels=['TCGA', 'external cohort', 'shared variables', 'matched\nsubsets', 'Spearman ρ', 'Pearson r'],
              colWidths=[0.07, 0.14, 0.36, 0.09, 0.1, 0.1], loc='upper center', cellLoc='center')
tb.auto_set_font_size(False); tb.set_fontsize(8); tb.scale(1, 1.45)
for (i, j), c in tb.get_celld().items():
    c.set_edgecolor('#ccc')
    if i == 0: c.set_facecolor('#eeeeee'); c.set_text_props(fontweight='bold'); c.set_height(c.get_height() * 1.4)
    if j == 2 and i > 0: c.set_text_props(ha='left'); c._loc = 'left'
fig.savefig(OUT + 'fig6.png', dpi=150, bbox_inches='tight')
print(T[['cohort', 'tcga', 'n_subsets', 'pct_better', 'loss', 'slope', 'slope_noise_checked', 'cross', 'pct_resamples_crossing',
         'matched_subsets', 'spearman_tcga_vs_ext', 'pearson_tcga_vs_ext']].round(3).to_string(index=False))
