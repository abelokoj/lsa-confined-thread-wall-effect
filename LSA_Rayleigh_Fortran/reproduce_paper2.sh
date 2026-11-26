#!/usr/bin/env bash
# Reproduces the results of Paper 2 (wall effects on the Rayleigh-Plateau instability of a thread in a light viscous ambient) at three levels.
#   ./reproduce_paper2.sh 1          build the analysis code, run its tests and print the confined growth rates
#   ./reproduce_paper2.sh 2          recompute every table, figure and check of the paper (no solver)
#   ./reproduce_paper2.sh 3          submit the fifteen simulations (this folder must be inside amrex_implicit_interfaces/NS_fluids_lib)
#   ./reproduce_paper2.sh 3-analyse  compare the new runs with the confined relation, once every job has ended
# Run it from this folder. Level 2 writes the figures to paper2/figures/, rewrites the tables in paper2/data/ and compares them with the stored ones.
set -euo pipefail
cd "$(dirname "$0")"
PY=${PYTHON:-python3}

level1 () {
  make -j8
  make test | tail -12
  # the air-like ambient of the simulations at k R0 = 0.7: the wall at 4 R0 and at 1.5 R0 (Table 1)
  ./bin/run_lsa theory --kr0 0.7 --r-wall 4 | tail -4
  ./bin/run_lsa theory --kr0 0.7 --r-wall 1.5 | tail -4
}

level2 () {
  [[ -x bin/run_lsa ]] || level1
  "$PY" -c "import numpy, scipy, matplotlib" 2>/dev/null || { echo "install the Python packages first: $PY -m pip install -r plotting/requirements.txt"; exit 1; }
  mkdir -p figures
  rm -rf figures/stored_data && cp -r paper2/data figures/stored_data
  {
    echo "== Section 5. Checks of the evaluation"
    echo "-- unbounded limit, against the independent Fortran tables of the two-fluid relation"
    rm -rf figures/theory && cp -r results/ch6_figures_2026-10-22/theory figures/
    (cd plotting && "$PY" make_confined_tables.py ../figures/theory)
    echo "-- the wall moving out, at k R0 = 0.7, Oh = 0.02 and the air-like ambient"
    (cd paper2/code && "$PY" -c "
from confined_dispersion import growth_rate
s0 = growth_rate(0.7, 0.02, 1.225e-3, 2.6e-4)
for rw in (1.5, 2, 3, 4, 6, 20):
    print(f'R_w = {rw:4}: {100 * (growth_rate(0.7, 0.02, 1.225e-3, 2.6e-4, rw) / s0 - 1):+.3e} %')")
    echo "-- inviscid confined limit"
    (cd plotting && "$PY" check_inviscid_limit.py)
    echo
    echo "== Section 6 and Table 1. The simulations: interface rates refitted from their records, then compared with the relation"
    scripts/refit_records.sh | grep -E "^(case|conf|mug10|rhog10|sw070|visc)"
    (cd paper2/code && "$PY" compare_simulations.py)
    echo
    echo "== Section 7 and Table 2. Maps of the wall effect, the collapse on delta_mu and the fitted gap"
    for s in wall_maps wall_collapse wall_density fit_check fit_bound density_zero_test; do
      echo "-- $s.py"; (cd paper2/code && "$PY" $s.py)
    done
    echo "-- the cases of Ooms et al. (Section 8)"
    (cd paper2/code && "$PY" ooms_cases.py)
    (cd paper2/code && "$PY" plot_wall_maps.py)                                   # Figs. 2 and 3
    echo
    echo "== The tables in paper2/data/ against the stored ones (numbers compared to a relative 1e-4)"
    compare_tables figures/stored_data paper2/data
  } 2>&1 | tee figures/paper2_tables.txt
  echo "Output: figures/paper2_tables.txt; figures: paper2/figures/ (see README.md)."
}

compare_tables () {   # stored new
  "$PY" - "$1" "$2" <<'PY'
import os, re, sys
old, new = sys.argv[1:3]
num = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")
for name in sorted(os.listdir(old)):
    a = open(os.path.join(old, name), encoding="utf-8").read().replace("\r", "").splitlines()
    b = open(os.path.join(new, name), encoding="utf-8").read().replace("\r", "").splitlines()
    if a == b:
        print(f"identical            {name}"); continue
    worst, other = 0.0, []
    for i, (x, y) in enumerate(zip(a, b)):
        if x == y or x.startswith("Runtime"): continue   # fit_bound.txt reports its own run time
        nx, ny = num.findall(x), num.findall(y)
        if num.sub("#", x) != num.sub("#", y) or len(nx) != len(ny):
            other.append(i + 1); continue
        for p, q in zip(nx, ny):
            p, q = float(p), float(q)
            if abs(p - q) > 1e-12:   # values below 1e-12 are round-off (a wall effect of order 1e-15 far from the thread)
                worst = max(worst, abs(p - q) / max(abs(p), abs(q)))
    if len(a) != len(b): other.append("length")
    if other:
        print(f"DIFFERENT            {name}: text differs at lines {other[:5]}")
    else:
        print(f"{'same to 1e-4' if worst <= 1e-4 else 'DIFFERENT':20} {name} (largest relative difference {worst:.1e}, ignoring differences below 1e-12)")
PY
}

level3 () {
  [[ -x ../amr2d.gnu.MPI.ex ]] || { echo "build the solver first: ../amr2d.gnu.MPI.ex not found (see README.md, Level 3)"; exit 1; }
  make -j8 && make test | tail -1
  ACCOUNT="${ACCOUNT:?set ACCOUNT to your Slurm account}" scripts/submit_campaign.sh campaign_p2_16.txt campaign_check.txt
  echo "When every job has ended (squeue --me), run: ./reproduce_paper2.sh 3-analyse"
}

level3_analyse () {
  scripts/refit_records.sh runs
  "$PY" paper2/code/compare_simulations.py runs/analysis/summary.csv runs/analysis/simulations_vs_relation.csv
}

case "${1:-}" in
  1) level1 ;;
  2) level2 ;;
  3) level3 ;;
  3-analyse) level3_analyse ;;
  *) sed -n '2,8p' "$0"; exit 1 ;;
esac
