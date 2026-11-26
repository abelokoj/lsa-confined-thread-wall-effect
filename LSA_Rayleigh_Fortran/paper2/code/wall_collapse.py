"""
wall_collapse.py
================
Denser sweep to test whether the gap needed for a given correction scales with the radial decay length of the ambient
disturbance, 1/l2 = (x_m^2 + s_m/nu2)^(-1/2), which tends to the diffusion length (nu2/s_m)^(1/2) when that is short
and to 1/x_m when it is long. Writes ../data/wall_collapse.csv and prints the range of gap * l2 for 1% and 0.1%.

    python wall_collapse.py
"""
import csv, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wall_maps import RHO2, DATA, most_unstable, wall_radius

OHS = (0.005, 0.01, 0.02, 0.05, 0.1)
MS = np.round(np.geomspace(0.005, 2.0, 8), 4)

rows = []
for oh in OHS:
    for m in MS:
        mu2 = m * oh
        km, sm = most_unstable(oh, mu2)
        l2 = np.sqrt(km ** 2 + sm * RHO2 / mu2)
        g1, g01 = wall_radius(1e-2, km, oh, mu2, sm) - 1, wall_radius(1e-3, km, oh, mu2, sm) - 1
        rows.append((oh, m, km, sm, np.sqrt(mu2 / RHO2 / sm), 1 / l2, g1, g01))
        print(f'Oh={oh:<6} m={m:<7} gap_1%*l2={g1 * l2:.2f}  gap_0.1%*l2={g01 * l2:.2f}', flush=True)
with open(os.path.join(DATA, 'wall_collapse.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['Oh', 'mu_ratio', 'x_m', 's_m', 'diffusion_length', 'decay_length', 'gap_1pct', 'gap_0p1pct'])
    w.writerows([[f'{v:.6g}' for v in r] for r in rows])
a = np.array(rows)
for col, name in ((6, '1%'), (7, '0.1%')):
    r = a[:, col] / a[:, 5]
    print(f'gap/decay length for {name}: {r.min():.2f} to {r.max():.2f} (median {np.median(r):.2f}); gap/diffusion length: {(a[:, col] / a[:, 4]).min():.2f} to {(a[:, col] / a[:, 4]).max():.2f}')
