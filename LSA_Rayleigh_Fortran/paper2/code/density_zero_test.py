"""
density_zero_test.py
====================
Direct test that the near-wall effect does not come from the inertia of the ambient: the wall effect at x = 0.69 with the
ambient density of the maps (rho2/rho1 = 1.225e-3) and with rho2 -> 1e-8, at fixed ambient viscosity. Section "Maps of the
wall effect" quotes these values. Writes ../data/density_zero_test.csv.

    python density_zero_test.py
"""
import csv, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from confined_dispersion import growth_rate

X = 0.69
CASES = ((0.005, 0.013, 1.5), (0.02, 0.013, 1.5), (0.02, 0.013, 3.0), (0.02, 1.0, 2.0), (0.1, 0.1, 2.0), (0.1, 1.0, 4.0))
rows = []
for oh, m, Rw in CASES:
    e = [100 * (growth_rate(X, oh, rho, m * oh, Rw) / growth_rate(X, oh, rho, m * oh) - 1) for rho in (1.225e-3, 1e-8)]
    rows.append([oh, m, Rw, f'{e[0]:.4f}', f'{e[1]:.4f}', f'{abs(e[1] / e[0] - 1):.4f}', f'{e[1] - e[0]:+.4f}'])
    print(rows[-1])
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'density_zero_test.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['Oh', 'mu_ratio', 'Rw', 'wall_effect_pct_rho2_1.225e-3', 'wall_effect_pct_rho2_1e-8', 'relative_change', 'difference_pp'])
    w.writerows(rows)
