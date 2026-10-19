#Figure 5: clinical-only C against clinical + RNA C for every clinical subset, 21 TCGA cancers (IPF-LASSO and offset)
#Per panel: share above the diagonal, mean |dC|, crossing point and slope, with 95% intervals from resampling repeats,
#and the noise-checked slope (slope divided by the reliability of clinical C)
#Usage, from the project root: python figure_generation/fig5.py --results-dir results --out-dir figures
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

nc = 4; nr = int(np.ceil((len(codes) + 1) / nc))
vmax = max(p[3].max() for p in panels.values()); cmap = plt.cm.viridis; J, O, G = '#c0304f', '#8e5bb5', '#1a8a3c'
fig, axes = plt.subplots(nr, nc, figsize=(3.6 * nc, 3.7 * nr))
for ax, c in zip(axes.flat, codes):
    x, y, off, nv, anc, m, po, r = panels[c]
    vals = np.r_[x, y, off[np.isfinite(off)], anc]
    lo, hi = np.nanmin(vals) - 0.01, np.nanmax(vals) + 0.01
    ax.fill_between([lo, hi], [lo, hi], hi, color='#e8f4ee', zorder=0)
    ax.plot([lo, hi], [lo, hi], color='#888', lw=0.9, ls=(0, (4, 3)), zorder=1)
    sz = 5 if len(x) > 500 else 12
    ax.scatter(x, off, s=sz * 0.6, marker='s', color='#bbb', alpha=0.5, lw=0, zorder=2)
    ax.scatter(x, y, s=sz, c=nv, cmap=cmap, vmin=0, vmax=vmax, alpha=0.85, lw=0, zorder=3)
    ax.scatter([anc[0]], [anc[1]], marker='D', s=28, color='#e07b22', zorder=4)
    xs = np.array([x.min(), x.max()]); pj = np.polyfit(x, y, 1)
    ax.plot(xs, np.polyval(pj, xs), color=J, lw=1.5, zorder=5)
    ax.plot(xs, np.polyval(po, xs), color=O, lw=1.5, ls='--', zorder=5)
    sg = r['slope_noise_checked']  #noise-checked slope, drawn through the centre of the points
    ax.plot(xs, y.mean() + sg * (xs - x.mean()), color=G, lw=1.4, ls=(0, (1.5, 1.5)), zorder=5)
    if r['cross'] is not None:
        ax.scatter([r['cross']], [r['cross']], marker='x', s=40, color=J, lw=1.6, zorder=6)
        if r['cross_lo'] is not None:  #interval of the crossing point along the diagonal
            ax.plot([r['cross_lo'], r['cross_hi']], [r['cross_lo'], r['cross_hi']], color=J, lw=3, alpha=0.35, zorder=5)
    if r['cross'] is not None:  #crossing text in two short lines
        c1, c2 = f"crosses at {r['cross']:.3f} [{r['cross_lo']:.3f}–{r['cross_hi']:.3f}]", f"in {r['pct_resamples_crossing']:.0f}% of resamples"
    else:
        c1, c2 = "no crossing", f"({r['pct_resamples_crossing']:.0f}% of resamples cross)"
    lines = [(f"offset slope {r['slope_offset']:.2f}", O),  #short lines first, so the block stays below the diagonal
             (f"noise-checked slope {r['slope_noise_checked']:.2f}", G),
             (f"slope {r['slope']:.2f} [{r['slope_lo']:.2f}–{r['slope_hi']:.2f}]", J),
             (f"{r['pct_above']:.0f}% above diagonal", '#2a9d6f' if r['pct_above'] >= 50 else '#c0304f'),
             (f"{r['pct_within_0_01']:.0f}% within ±0.01", '#333'),
             (f"mean |ΔC| {r['loss']:.3f} [{r['loss_lo']:.3f}–{r['loss_hi']:.3f}]", '#333'),
             (c1, J), (c2, J)]  #crossing reads top to bottom
    n = len(lines)
    for i, (t, col) in enumerate(lines):  #stacked upwards from the bottom-right corner
        ax.text(0.97, 0.03 + (n - 1 - i) * 0.052, t, transform=ax.transAxes, fontsize=5.6, color=col, ha='right', va='bottom')
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect('equal')
    ax.set_title(f'{c}  ({int(nv.max())} vars, {r["repeats"]} repeats)', fontsize=9, pad=3)
    ax.set_xlabel("clinical-only model (Uno's C)", fontsize=7); ax.set_ylabel("clinical + omics (Uno's C)", fontsize=7)
    ax.tick_params(labelsize=6)
for ax in axes.flat[len(codes):]: ax.axis('off')
lg = axes.flat[len(codes)]
h = [Line2D([], [], marker='o', ls='', color='#666', label='subset, joint penalised fit (IPF-LASSO)'),
     Line2D([], [], marker='s', ls='', color='#bbb', label='subset, frozen-offset fit'),
     Line2D([], [], marker='D', ls='', color='#e07b22', label='omics-only anchor'),
     Line2D([], [], color='#888', ls=(0, (4, 3)), label='identity: omics adds nothing'),
     Line2D([], [], color=J, label='OLS fit, joint penalised'), Line2D([], [], color=O, ls='--', label='OLS fit, frozen-offset'),
     Line2D([], [], color=G, ls=(0, (1.5, 1.5)), lw=1.4, label='noise-checked fit, joint penalised'),
     Line2D([], [], marker='x', ls='', color=J, label='joint fit crosses the identity line'),
     Line2D([], [], color=J, lw=3, alpha=0.35, label='95% interval of the crossing')]
lg.legend(handles=h, loc='upper center', frameon=False, fontsize=7, title='markers / lines', title_fontsize=8)
lg.text(0.5, 0.30, '[ ] = 95% interval from resampling repeats (1,000×)\nnoise-checked = slope ÷ reliability of clinical C',
        ha='center', va='top', fontsize=6.5, transform=lg.transAxes, color='#333')
cax = lg.inset_axes([0.08, 0.06, 0.84, 0.045])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0, vmax), cmap=cmap), cax=cax, orientation='horizontal')
cb.set_label('number of clinical variables in the baseline', fontsize=7.5); cb.ax.tick_params(labelsize=6.5)
fig.tight_layout(h_pad=1.2, w_pad=1.0)
fig.savefig(OUT + 'fig5.png', dpi=160)
pd.DataFrame(rows).to_csv(OUT + 'fig5_metrics.csv', index=False)
print(pd.DataFrame(rows)[['cancer', 'repeats', 'loss', 'loss_lo', 'loss_hi', 'slope', 'slope_lo', 'slope_hi',
                          'slope_noise_checked', 'cross', 'cross_lo', 'cross_hi', 'pct_resamples_crossing']].to_string(index=False))
