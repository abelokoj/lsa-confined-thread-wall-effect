"""
wall_maps.py
============
Paper 2 parameter maps: how a coaxial no-slip wall at r = Rw changes the growth rate of a viscous liquid thread in a light viscous ambient, computed from the confined two-fluid relation (confined_dispersion.py, copied from LSA_Rayleigh_Fortran/plotting). Capillary units R0 = sigma = rho1 = 1, so mu1 = Oh and mu2 = m Oh with m = mu2/mu1.

    python wall_maps.py            # writes ../data/wall_maps_summary.csv and ../data/wall_maps_vs_Rw.csv

For every (Oh, m) at the air-like density ratio rho2/rho1 = 1.225e-3:
  - the unconfined most unstable wavenumber x_m and rate s_m;
  - the wall radius at which the confinement correction at x_m falls to 1% and 0.1% (R_1, R_01);
  - the viscous diffusion length of the ambient over one growth time, (nu2 / s_m)^(1/2);
  - against Rw: the correction at x_m, the inviscid confined correction at the same Rw, and the confined most unstable wavenumber and rate.
The script first reproduces two values from the dissertation (App. D) and stops if either fails.
"""
import csv, os, sys
import numpy as np
from scipy.optimize import brentq, minimize_scalar
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from confined_dispersion import growth_rate, inviscid_bound

RHO2 = 1.225e-3
OHS = (0.005, 0.02, 0.1)
MS = (0.013, 0.1, 1.0)
RWS = np.round(np.geomspace(1.2, 20.0, 41), 4)
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')


def most_unstable(oh, mu2, Rw=None):
    r = minimize_scalar(lambda k: -growth_rate(k, oh, RHO2, mu2, Rw), bounds=(0.3, 0.95), method='bounded', options={'xatol': 1e-6})
    return r.x, -r.fun


def wall_radius(tol, k, oh, mu2, s_unc):
    """Smallest Rw at which |s(Rw)/s_unc - 1| = tol (the correction shrinks monotonically as the wall recedes)."""
    f = lambda Rw: abs(growth_rate(k, oh, RHO2, mu2, Rw) / s_unc - 1) - tol
    return brentq(f, 1.02, 60.0, xtol=1e-6)


def checks():
    oh, mu2 = 0.02, 2.6e-4                                   # the gas of the dissertation computations, m = 0.013
    su = growth_rate(0.7, oh, RHO2, mu2)
    d15 = growth_rate(0.7, oh, RHO2, mu2, 1.5) / su - 1
    d20 = growth_rate(0.7, oh, RHO2, mu2, 20.0) / su - 1
    print(f'check: correction at Rw = 1.5, x = 0.7: {100 * d15:.3f}% (dissertation App. D: -3.159%); at Rw = 20: {d20:.1e}')
    assert abs(100 * d15 + 3.159) < 5e-4, 'does not reproduce the dissertation value at Rw = 1.5'
    assert abs(d20) < 1e-9, 'wall at 20 R0 should not change the rate'


def main():
    checks()
    os.makedirs(DATA, exist_ok=True)
    summary, rows = [], []
    for oh in OHS:
        for m in MS:
            mu2 = m * oh
            km, sm = most_unstable(oh, mu2)
            r1 = wall_radius(1e-2, km, oh, mu2, sm)
            r01 = wall_radius(1e-3, km, oh, mu2, sm)
            ldiff = np.sqrt(mu2 / RHO2 / sm)
            s_inv = inviscid_bound(km, RHO2)
            summary.append([oh, m, f'{km:.5f}', f'{sm:.6e}', f'{r1:.4f}', f'{r01:.4f}', f'{ldiff:.4f}'])
            print(f'Oh={oh:<6} m={m:<6} x_m={km:.4f} s_m={sm:.5f}  R_1%={r1:.3f}  R_0.1%={r01:.3f}  diffusion length={ldiff:.3f}', flush=True)
            for Rw in RWS:
                s = growth_rate(km, oh, RHO2, mu2, Rw)
                inv = inviscid_bound(km, RHO2, Rw) / s_inv - 1
                kc, sc = most_unstable(oh, mu2, Rw)
                rows.append([oh, m, f'{Rw:.4f}', f'{km:.5f}', f'{s / sm - 1:.6e}', f'{inv:.6e}', f'{kc:.5f}', f'{sc:.6e}'])
    with open(os.path.join(DATA, 'wall_maps_summary.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['Oh', 'mu_ratio', 'x_m', 's_m', 'R_1pct', 'R_0p1pct', 'diffusion_length']); w.writerows(summary)
    with open(os.path.join(DATA, 'wall_maps_vs_Rw.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['Oh', 'mu_ratio', 'Rw', 'x_m_unconfined', 'rel_change_viscous', 'rel_change_inviscid', 'x_m_confined', 's_m_confined']); w.writerows(rows)
    print('written to', os.path.abspath(DATA))


if __name__ == '__main__':
    main()
