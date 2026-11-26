"""
plot_wall_maps.py
=================
Figures for Paper 2 from the tables written by wall_maps.py, in the house style of pub_style.py (PNG, PDF and SVG):

    figures/wall_correction.*   reduction of the growth rate at x_m against the wall radius, one panel per Oh,
                                one curve per viscosity ratio, with the inviscid confined correction for reference
    figures/wall_distance.*     gap at which the reduction falls to 0.1% and 1%, against the dynamic-viscosity
                                length delta_mu = (mu2/(rho1 s_m))^(1/2), four density ratios, cubic fit to the two
                                lightest (wall_collapse.csv, wall_density.csv -> wall_collapse_fit.csv)

    python plot_wall_maps.py
"""
import csv, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FixedLocator, NullFormatter
from pub_style import apply_style, savefig_all, subcaption, PALETTE

apply_style()
HERE = os.path.dirname(os.path.abspath(__file__))
DATA, FIGS = os.path.join(HERE, '..', 'data'), os.path.join(HERE, '..', 'figures')
REF = '0.4'                                   # grey for reference lines (inviscid relation, threshold levels)


def read(name):
    return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(open(os.path.join(DATA, name)))]


def sci(v):
    """1.225e-3 -> 1.225 \\times 10^{-3} for a mathtext label."""
    m, e = f'{v:.3e}'.split('e')
    return rf'{float(m):g} \times 10^{{{int(e)}}}'


def log_ticks(axis, ticks):
    """Label a log axis at chosen values in Computer Modern; minor ticks stay but carry no labels."""
    axis.set_major_locator(FixedLocator(ticks))
    axis.set_major_formatter(matplotlib.ticker.FixedFormatter([f'${t:g}$' for t in ticks]))
    axis.set_minor_formatter(NullFormatter())


def fig_correction(rows):
    ohs = sorted({r['Oh'] for r in rows})
    mus = sorted({r['mu_ratio'] for r in rows})
    style = dict(zip(mus, zip(PALETTE, ('-', '--', '-.'))))
    fig, axes = plt.subplots(1, len(ohs), figsize=(7.0, 2.7), sharey=True)
    for ax, oh, tag in zip(axes, ohs, 'abc'):
        for m in mus:
            sel = [r for r in rows if r['Oh'] == oh and r['mu_ratio'] == m]
            ax.loglog([r['Rw'] for r in sel], [-100 * r['rel_change_viscous'] for r in sel], style[m][1],
                      color=style[m][0], lw=1.2, label=rf'$\mu_2/\mu_1 = {m:g}$')
        inv = [r for r in rows if r['Oh'] == oh and r['mu_ratio'] == 0.013]
        ax.loglog([r['Rw'] for r in inv], [-100 * r['rel_change_inviscid'] for r in inv], ':', color='k', lw=1.0,
                  label='inviscid, confined')
        for lev in (1, 0.1):
            ax.axhline(lev, color=REF, lw=0.8)
        subcaption(ax, rf'({tag}) $\mathrm{{Oh}} = {oh:g}$')
        ax.set_xlabel(r'$R_w / R_0$')
        ax.set_xlim(1.2, 20); ax.set_ylim(1e-4, 100)  # 100%: the largest reduction mapped is 98.3%
        log_ticks(ax.xaxis, [1.5, 2, 3, 5, 10, 20])
    axes[0].set_ylabel(r'reduction of $s$ at $x_m$ (%)')
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc='lower center', ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.0))
    savefig_all(fig, 'wall_correction', outdir=FIGS, close=True)


def fig_distance(rows):
    """Gap at which the reduction falls to 0.1% (a) and 1% (b) against delta_mu = (mu2 / (rho1 s_m))^(1/2)
    = ell_nu (rho2/rho1)^(1/2), for four density ratios. Line: cubic fit in log-log to the two lightest ratios pooled;
    coefficients and residuals (pooled and per ratio) go to data/wall_collapse_fit.csv."""
    rhos = (1.225e-4, 1.225e-3, 1.225e-2, 1.225e-1)
    col, mk = dict(zip(rhos, PALETTE)), dict(zip(rhos, ('^', 'o', 's', 'D')))
    LIGHT = rhos[:2]
    d = lambda r: r['diffusion_length'] * np.sqrt(r['rho2'])
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0), sharey=True)
    fits = []
    for ax, key, tag in ((axes[0], 'gap_0p1pct', r'(a) $0.1\%$ reduction'), (axes[1], 'gap_1pct', r'(b) $1\%$ reduction')):
        for group in [LIGHT] + [(r,) for r in rhos]:
            sel = [r for r in rows if r['rho2'] in group]
            X = np.array([d(r) for r in sel]); G = np.array([r[key] for r in sel])
            cf = np.polyfit(np.log(X), np.log(G), 3)
            res = np.abs(np.exp(np.log(G) - np.polyval(cf, np.log(X))) - 1).max()
            fits.append(['+'.join(f'{g:g}' for g in group), key] + [f'{v:.5f}' for v in cf] + [f'{X.min():.4g}', f'{X.max():.4g}', f'{res:.4f}'])
            if group == LIGHT:
                xs = np.geomspace(X.min(), X.max(), 100)
                ax.loglog(xs, np.exp(np.polyval(cf, np.log(xs))), '-', color=REF, lw=1.0, zorder=0,
                          label='fit, two lightest')
        for rho in rhos:
            sel = [r for r in rows if r['rho2'] == rho]
            ax.loglog([d(r) for r in sel], [r[key] for r in sel], mk[rho], color=col[rho], ms=3.5, mfc='none', mew=0.8,
                      label=rf'$\rho_2/\rho_1 = {sci(rho)}$')
        subcaption(ax, tag)
        ax.set_xlabel(r'$\delta_\mu / R_0$')
        ax.set_xlim(7e-3, 1.5); ax.set_ylim(0.25, 10)  # gaps span 0.32 to 8.2
        log_ticks(ax.yaxis, [0.5, 1, 2, 5, 10])
    axes[0].set_ylabel(r'gap $(R_w - R_0)/R_0$')
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc='lower center', ncol=5, frameon=False, columnspacing=1.2, bbox_to_anchor=(0.5, 0.0))
    savefig_all(fig, 'wall_distance', outdir=FIGS, close=True)
    with open(os.path.join(DATA, 'wall_collapse_fit.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['rho2', 'threshold', 'c3', 'c2', 'c1', 'c0', 'delta_min', 'delta_max', 'max_rel_residual']); w.writerows(fits)
    return fits


if __name__ == '__main__':
    fig_correction(read('wall_maps_vs_Rw.csv'))
    rows = [dict(r, rho2=1.225e-3) for r in read('wall_collapse.csv')] + read('wall_density.csv')
    for f in fig_distance(rows):
        print(f)
