#Figure 6: clinical-only C against clinical + RNA C for every clinical subset, 21 TCGA cancers (IPF-LASSO and offset)
#fig6_metrics.csv (supplementary table): share above the diagonal, mean |dC|, crossing point and slope, with 95% intervals
#from resampling repeats, and the noise-checked slope (slope divided by the reliability of clinical C); the figure shows the slope
#Usage, from the project root: python figure_generation/fig6.py --results-dir results --out-dir figures
import argparse, glob, os
import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
matplotlib.rcParams["font.family"] = ["Inter", "sans-serif", "DejaVu Sans"]  #font of the paper's figures; DejaVu Sans if Inter is not installed
from matplotlib.lines import Line2D

ap = argparse.ArgumentParser()
ap.add_argument("--results-dir", default="results"); ap.add_argument("--out-dir", default="figures")
args = ap.parse_args(); os.makedirs(args.out_dir, exist_ok=True)
R = os.path.join(args.results_dir, "")
OUT = os.path.join(args.out_dir, "")
DELTA = 0.01  #fixed equivalence margin
B = 1000  #resamples of repeats
rng = np.random.default_rng(20261011)  #fixed seed so the intervals are reproducible

def fit(x, y):  #OLS line y = a*x + b, None if too few points
    m = np.isfinite(x) & np.isfinite(y)
    return np.polyfit(x[m], y[m], 1) if m.sum() > 2 else None

def crossing(a, b, lo, hi):  #x where the line meets y = x, NaN if outside [lo, hi] or parallel
    if abs(1 - a) < 1e-9: return np.nan
    x = b / (1 - a)
    return x if lo <= x <= hi else np.nan

def metrics(x, y):  #the panel numbers for one set of subset means
    a, b = np.polyfit(x, y, 1)
    d = y - x
    return dict(slope=a, cross=crossing(a, b, x.min(), x.max()), loss=np.abs(d).mean(),
                within=(np.abs(d) <= DELTA).mean() * 100, above=(d > 0).mean() * 100)

def load(code):  #subset x repeat matrices of clinical C and IPF-LASSO C, plus subset means for offset
    p = pd.read_csv(R + code + '_per_repeat.csv')
    p = p[~p.is_anchor.astype(str).str.upper().eq('TRUE')]
    f = p[p.arm == 'fused']; o = p[p.arm == 'offset']
    X = f.pivot(index='subset_id', columns='repeat_id', values='c_ref')
    Y = f.pivot(index='subset_id', columns='repeat_id', values='c_model').reindex_like(X)
    off = o.groupby('subset_id').c_model.mean().reindex(X.index)
    nv = f.groupby('subset_id').n_vars.first().reindex(X.index)
    s = pd.read_csv(R + code + '_subset_summary.csv')
    anc = s[(s.arm == 'fused') & s.is_anchor.astype(str).str.upper().eq('TRUE')].iloc[0]
    return X.values, Y.values, off.values, nv.values, (anc.c_ref, anc.c_model)

codes = sorted(os.path.basename(f).replace('_per_repeat.csv', '') for f in glob.glob(R + '*_per_repeat.csv'))
codes = [c for c in codes if '_' not in c]  #main lattice only, not the robustness or published runs
rows = []; panels = {}
for c in codes:
    X, Y, off, nv, anc = load(c)
    ok = np.isfinite(np.nanmean(Y, 1))  #subsets with at least one fitted repeat
    X, Y, off, nv = X[ok], Y[ok], off[ok], nv[ok]
    x, y = np.nanmean(X, 1), np.nanmean(Y, 1)
    m = metrics(x, y)
    nrep = np.sum(np.isfinite(X), 1)
    se2 = np.nanvar(X, 1, ddof=1) / nrep  #leftover noise in each subset's clinical C
    rel = 1 - np.nanmean(se2) / np.var(x, ddof=1)  #reliability: share of the clinical-C spread that is real
    Rn = X.shape[1]; boot = []
    for _ in range(B):  #resample repeats with replacement, recompute every subset mean and the metrics
        idx = rng.integers(0, Rn, Rn)
        xb, yb = np.nanmean(X[:, idx], 1), np.nanmean(Y[:, idx], 1)
        g = np.isfinite(xb) & np.isfinite(yb)
        boot.append(metrics(xb[g], yb[g]))
    bt = pd.DataFrame(boot)
    q = lambda k: np.nanpercentile(bt[k], [2.5, 97.5])
    cr_share = np.isfinite(bt.cross).mean() * 100  #share of resamples where the fit crosses inside the range
    po = fit(x, off)
    rows.append(dict(cancer=c, n_subsets=len(x), repeats=Rn, pct_above=round(m['above'], 1),
                     loss=round(m['loss'], 4), loss_lo=round(q('loss')[0], 4), loss_hi=round(q('loss')[1], 4),
                     pct_within_0_01=round(m['within'], 1),
                     slope=round(m['slope'], 3), slope_lo=round(q('slope')[0], 3), slope_hi=round(q('slope')[1], 3),
                     reliability=round(rel, 3), slope_noise_checked=round(m['slope'] / rel, 3),
                     cross=None if np.isnan(m['cross']) else round(m['cross'], 4),
                     cross_lo=None if cr_share == 0 else round(q('cross')[0], 4),
                     cross_hi=None if cr_share == 0 else round(q('cross')[1], 4),
                     pct_resamples_crossing=round(cr_share, 1), slope_offset=round(po[0], 3)))
    panels[c] = (x, y, off, nv, anc, m, po, rows[-1])

#drawn at print size (170 mm wide): each panel shows the subsets, y = x, both fits and the IPF-LASSO slope;
#every other number is in fig6_metrics.csv (supplementary table)
NC = 5; NROW = int(np.ceil((len(codes) + 1) / NC))
vmax = max(p[3].max() for p in panels.values()); cmap = plt.cm.viridis; J, O = '#c0304f', '#8e5bb5'
FH = 1.36 * NROW + 0.35
fig = plt.figure(figsize=(6.7, FH))
gs = fig.add_gridspec(NROW, NC, hspace=0.5, wspace=0.34, left=0.075, right=0.99, top=1 - 0.12 / FH, bottom=0.42 / FH)
for k, c in enumerate(codes):
    ax = fig.add_subplot(gs[k // NC, k % NC])
    x, y, off, nv, anc, m, po, r = panels[c]
    vals = np.r_[x, y, off[np.isfinite(off)], [v for v in anc if np.isfinite(v)]]
    pad = 0.03 * (np.nanmax(vals) - np.nanmin(vals)); lo, hi = np.nanmin(vals) - pad, np.nanmax(vals) + pad
    ax.fill_between([lo, hi], [lo, hi], hi, color='#e8f4ee', zorder=0)
    ax.plot([lo, hi], [lo, hi], color='#888', lw=0.6, ls=(0, (3, 2)), zorder=1)
    sz = 0.8 if len(x) > 500 else 3
    ax.scatter(x, off, s=sz * 0.6, marker='s', color='#bbb', alpha=0.5, lw=0, zorder=2)
    ax.scatter(x, y, s=sz, c=nv, cmap=cmap, vmin=0, vmax=vmax, alpha=0.85, lw=0, zorder=3)
    ax.scatter([anc[0]], [anc[1]], marker='D', s=7, color='#e07b22', zorder=4)
    xs = np.array([x.min(), x.max()]); pj = np.polyfit(x, y, 1)
    ax.plot(xs, np.polyval(pj, xs), color=J, lw=0.8, zorder=5); ax.plot(xs, np.polyval(po, xs), color=O, lw=0.8, ls='--', zorder=5)
    if not np.isnan(m['cross']): ax.scatter([m['cross']], [m['cross']], marker='x', s=12, color=J, lw=0.8, zorder=6)
    ax.text(0.96, 0.04, f"slope {m['slope']:.2f}", transform=ax.transAxes, fontsize=5.6, color=J, ha='right', va='bottom')
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect('equal')
    ax.set_title(c, fontsize=6.5, pad=2)
    ax.tick_params(labelsize=5, width=0.5, length=2, pad=1.2)
    ax.xaxis.set_major_locator(plt.MaxNLocator(3)); ax.yaxis.set_major_locator(plt.MaxNLocator(3))
    for sp in ax.spines.values(): sp.set_linewidth(0.5)
fig.text(0.53, 0.1 / FH, "clinical-only model (Uno's C)", ha='center', fontsize=7)
fig.text(0.01, 0.5, "clinical + RNA model (Uno's C)", va='center', rotation=90, fontsize=7)
#legend and colour bar in the empty slots of the last row
lg = fig.add_subplot(gs[NROW - 1, len(codes) % NC:]); lg.axis('off')
h = [Line2D([], [], marker='o', ls='', ms=3, color='#666', label='subset, IPF-Lasso'),
     Line2D([], [], marker='s', ls='', ms=3, color='#bbb', label='subset, offset model'),
     Line2D([], [], marker='D', ls='', ms=3, color='#e07b22', label='RNA-only model'),
     Line2D([], [], color='#888', lw=0.6, ls=(0, (3, 2)), label='no RNA gain (y = x)'),
     Line2D([], [], color=J, lw=0.8, label='fit, IPF-Lasso'), Line2D([], [], color=O, lw=0.8, ls='--', label='fit, offset model'),
     Line2D([], [], marker='x', ls='', ms=3.5, color=J, label='IPF-Lasso fit crosses y = x')]
lg.legend(handles=h, loc='upper left', frameon=False, fontsize=5.8, ncol=2, handlelength=1.8, columnspacing=1.0, borderaxespad=0)
cax = lg.inset_axes([0.05, 0.1, 0.9, 0.08])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0, vmax), cmap=cmap), cax=cax, orientation='horizontal')
cb.set_label('clinical variables in the baseline', fontsize=5.8, labelpad=1); cb.ax.tick_params(labelsize=5.5, width=0.5, length=2, pad=1)
cb.outline.set_linewidth(0.5)
fig.savefig(OUT + 'fig6.png', dpi=600)
pd.DataFrame(rows).to_csv(OUT + 'fig6_metrics.csv', index=False)
print(pd.DataFrame(rows)[['cancer', 'repeats', 'loss', 'loss_lo', 'loss_hi', 'slope', 'slope_lo', 'slope_hi',
                          'slope_noise_checked', 'cross', 'cross_lo', 'cross_hi', 'pct_resamples_crossing']].to_string(index=False))
