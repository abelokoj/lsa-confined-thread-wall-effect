"""
wall_density.py
===============
Repeats the 40-case sweep of wall_collapse.py (rho2/rho1 = 1.225e-3) at three other density ratios, so that
plot_wall_maps.py can test how the wall gap depends on the density of the ambient. Writes ../data/wall_density.csv.

    python wall_density.py
"""
import csv, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wall_maps
from wall_maps import DATA, most_unstable, wall_radius

RHOS = (1.225e-4, 1.225e-2, 1.225e-1)
OHS = (0.005, 0.01, 0.02, 0.05, 0.1)
MS = np.round(np.geomspace(0.005, 2.0, 8), 4)

rows = []
for rho2 in RHOS:
    wall_maps.RHO2 = rho2                                   # the helpers read the module-level density ratio
    for oh in OHS:
        for m in MS:
            mu2 = m * oh
            km, sm = most_unstable(oh, mu2)
            g1, g01 = wall_radius(1e-2, km, oh, mu2, sm) - 1, wall_radius(1e-3, km, oh, mu2, sm) - 1
            rows.append([rho2, oh, m, f'{km:.5f}', f'{sm:.6e}', f'{np.sqrt(mu2 / rho2 / sm):.5g}', f'{g1:.5g}', f'{g01:.5g}'])
    print(f'rho2/rho1 = {rho2:g} done', flush=True)
with open(os.path.join(DATA, 'wall_density.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['rho2', 'Oh', 'mu_ratio', 'x_m', 's_m', 'diffusion_length', 'gap_1pct', 'gap_0p1pct'])
    w.writerows(rows)
