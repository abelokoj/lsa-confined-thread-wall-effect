"""
compare_simulations.py
======================
Compares the growth rates of the October 2026 follow-up simulations (twin-run method of Paper 1, CLS-CMOF solver,
16 or 32 cells per radius; rows copied from LSA_Rayleigh_Fortran/results/followup_2026-10-22/analysis/summary.csv into
../data/simulations_followup_2026-10-22.csv) with the confined relation for each run's own wall and fluids.

    python compare_simulations.py      # prints the table and writes ../data/simulations_vs_relation.csv
    python compare_simulations.py SUMMARY OUT   # the same for new runs: SUMMARY is a runs/analysis/summary.csv, OUT the table to write

Fluids: Oh = 0.02, rho2/rho1 = 1.225e-3, mu2 = 2.6e-4 (mu2/mu1 = 0.013) and wall at Rw = 4 unless the run changes them.
The script checks that its unconfined rate reproduces the visc_two_fluid column of the table (Fortran code) to 1e-8.
"""
import csv, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from confined_dispersion import growth_rate

OH, RHO2, MU2 = 0.02, 1.225e-3, 2.6e-4
# label: (Rw, mu2, rho2, grid cells per radius, group)
RUNS = {'sw070_16': (4, MU2, RHO2, 16, 'gas'), 'mug10_16': (4, MU2 / 10, RHO2, 16, 'gas'), 'rhog10_16': (4, MU2, RHO2 / 10, 16, 'gas'),
        'visc010_16': (4, 0.1 * OH, RHO2, 16, 'viscous'), 'visc030_16': (4, 0.3 * OH, RHO2, 16, 'viscous'),
        'visc100_16': (4, 1.0 * OH, RHO2, 16, 'viscous'), 'visc100r6_16': (6, 1.0 * OH, RHO2, 16, 'viscous'),
        'visc100r8_16': (8, 1.0 * OH, RHO2, 16, 'viscous'),
        'conf30_16': (3, MU2, RHO2, 16, 'wall'), 'conf20_16': (2, MU2, RHO2, 16, 'wall'), 'conf15_16': (1.5, MU2, RHO2, 16, 'wall'),
        'conf15_32': (1.5, MU2, RHO2, 32, 'wall'), 'conf15k050_16': (1.5, MU2, RHO2, 16, 'wall'),
        'conf15k060_16': (1.5, MU2, RHO2, 16, 'wall'), 'conf15k080_16': (1.5, MU2, RHO2, 16, 'wall')}
HERE = os.path.dirname(os.path.abspath(__file__))
SUMMARY = sys.argv[1] if len(sys.argv) > 2 else os.path.join(HERE, '..', 'data', 'simulations_followup_2026-10-22.csv')
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, '..', 'data', 'simulations_vs_relation.csv')
rows = {r['label'].strip(): r for r in csv.DictReader(open(SUMMARY))}
out = []
for lab, (Rw, mu2, rho2, n, grp) in RUNS.items():
    r = rows[lab]; k = float(r['kr0'])
    su, sc = growth_rate(k, OH, rho2, mu2), growth_rate(k, OH, rho2, mu2, Rw)
    assert abs(su / float(r['visc_two_fluid']) - 1) < 1e-8, (lab, su, r['visc_two_fluid'])
    si = float(r['rate_interface']); sl = r['dmd_levelset_median']
    ls = f'{100 * (float(sl) / sc - 1):+.2f}' if sl else ''
    out.append([grp, lab, k, Rw, f'{mu2 / OH:.4g}', f'{rho2:.4g}', n, r['efolds'], f'{100 * (sc / su - 1):+.3f}',
                f'{100 * (si / su - 1):+.3f}', f'{100 * (si / sc - 1):+.3f}', ls])
hdr = ['group', 'run', 'x', 'Rw', 'mu_ratio', 'rho2', 'cells', 'efolds', 'theory_wall_effect_pct', 'measured_vs_unconfined_pct',
       'measured_vs_confined_pct', 'levelset_vs_confined_pct']
with open(OUT, 'w', newline='') as f:
    w = csv.writer(f); w.writerow(hdr); w.writerows(out)
for o in [hdr] + out:
    print(' | '.join(str(v) for v in o))
