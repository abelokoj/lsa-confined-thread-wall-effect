# Code and data for Paper 2: wall effects on the Rayleigh–Plateau instability

This repository contains the code and the data used in the paper *Wall effects on the Rayleigh–Plateau instability of a liquid thread in a light viscous ambient* (A. Okojunu, M. Sussman and U. S. Nair). Every table, figure and check of the paper can be recomputed from it, and the fifteen simulations can be repeated with it.

The confined relation, the maps and the fits are computed in Python and need neither the solver nor a cluster. The simulations were run with the two-phase solver `amrex_implicit_interfaces` of M. Sussman, at a fixed commit, and analyzed with the code in the folder `LSA_Rayleigh_Fortran/`, which is placed inside the solver tree to run them. Levels 1 and 2 below need only this repository.

## Contents

| Path | Contents |
|---|---|
| `LSA_Rayleigh_Fortran/paper2/code` | the evaluation of the confined relation (`confined_dispersion.py`), the maps, the fits of the gap, the comparison with the simulations and the figures |
| `LSA_Rayleigh_Fortran/paper2/data` | the tables these scripts write, as used in the paper |
| `LSA_Rayleigh_Fortran/src`, `tests`, `Makefile` | the analysis code of the simulations (Fortran 2008, LAPACK): the input decks, the interface-amplitude fit and the dispersion relations, with its test suite |
| `LSA_Rayleigh_Fortran/scripts`, `plotting` | the scripts that run and analyze the simulations, and the independent checks of the relation |
| `LSA_Rayleigh_Fortran/results` | the amplitude records and run records of the simulations, and the Fortran tables of the unbounded relation |
| `LSA_Rayleigh_Fortran/campaign_p2_16.txt`, `campaign_check.txt` | the case files of the simulations, one line per run |
| `LSA_Rayleigh_Fortran/reproduce_paper2.sh` | the commands of the three levels below |

## Level 1. Check the code (about one minute)

On Linux, or WSL on Windows, with gfortran, make and LAPACK:

```bash
git clone https://github.com/abelokoj/lsa-confined-thread-wall-effect.git
cd lsa-confined-thread-wall-effect/LSA_Rayleigh_Fortran
./reproduce_paper2.sh 1
```

The test suite must end with `ALL TEST PROGRAMS PASSED`. The script then prints the growth rate at k R0 = 0.7 for the air-like ambient of the simulations: 0.328617 with the wall at 4 R0 and 0.318305 with the wall at 1.5 R0, a reduction of 3.16%.

## Level 2. Recompute the tables, figures and checks (about three minutes)

```bash
python3 -m pip install -r plotting/requirements.txt
./reproduce_paper2.sh 2
```

The output is written to `figures/paper2_tables.txt` and the figures to `paper2/figures/`. Level 2 rewrites every table in `paper2/data/` and compares it with the stored one; on a different platform or version of SciPy, a few numbers may differ in the last printed digit.

| In the paper | Printed or drawn by Level 2 as |
|---|---|
| Section 5, checks of the evaluation | "Section 5. Checks of the evaluation" |
| Table 1 and Section 6, the simulations | the refitted interface rates of the fifteen runs, then the output of `compare_simulations.py` |
| Table 2 and Section 7, the maps | `wall_maps.py` (`paper2/data/wall_maps_summary.csv`) |
| Section 7, the collapse on δμ, the fitted gap and its test | `wall_collapse.py`, `wall_density.py`, `fit_check.py`, `fit_bound.py` |
| Section 7, the test at zero ambient density | `density_zero_test.py` |
| Section 8, the cases of Ooms et al. | `ooms_cases.py` |
| Figs. 2 and 3 | `paper2/figures/wall_correction`, `paper2/figures/wall_distance` |

Fig. 1 is a drawing of the evaluation and is not computed from data. The interface growth rates of the simulations are refitted from their amplitude records with the Fortran code that produced them.

## Level 3. Run the simulations again

The runs were made on the FSU Research Computing Center cluster with `gnu/13` and `openmpi/4.1.6`, under Slurm. They used commit `dddb827c` of a fork of the solver, whose solver source is identical to commit `0b8df42` of the original repository; the run records in `results/followup_2026-10-22/run_records/` give the commit, the number of ranks and the checksum of the executable for every run. The solver is used without modification.

```bash
export AII=$HOME/amrex_implicit_interfaces
git clone https://github.com/msussman42/amrex_implicit_interfaces.git $AII
cd $AII && git checkout 0b8df42
cp -r /path/to/lsa-confined-thread-wall-effect/LSA_Rayleigh_Fortran $AII/NS_fluids_lib/
module purge; module load gnu/13; module load openmpi
cd $AII/NS_fluids_lib && make -j16              # writes amr2d.gnu.MPI.ex
cd LSA_Rayleigh_Fortran
ACCOUNT=<your Slurm account> ./reproduce_paper2.sh 3
```

Use commit `0b8df42` and not the current head: at some later commits a calculation with fluids only stops at start-up with `expecting num_FSI_outer_sweeps>=2` (`LSA_Rayleigh_Fortran/docs/solver_notes.md`).

Level 3 submits two case files as Slurm job arrays. The number of ranks is that of the original runs, and the limit per task is the longest used for them:

| Case file | Runs | Ranks | Limit per task | In the paper |
|---|---|---|---|---|
| `campaign_p2_16.txt` | the gas, viscous and near-wall series on 16 cells per radius | 16 | 12 h | Table 1 |
| `campaign_check.txt` | the wall at 1.5 R0 on 32 cells; μ2/μ1 = 1 with the wall at 6 and 8 R0 | 16 | 48 h | Table 1 |

Each perturbed run has an unperturbed twin, except the two runs with the gas density or viscosity divided by ten, whose interface rate needs none. A task that reaches its limit is not lost: submitting the same command again skips the finished cases and restarts the others from their last checkpoint. When every job has ended, `./reproduce_paper2.sh 3-analyse` refits the new records and compares the new runs with the confined relation. Bitwise reproducibility has not been tested. Another compiler, MPI library or rank count changes the round-off, so the interface rates are expected to agree closely with the stored ones but not exactly.

## Origin

The analysis code of the simulations is that of the dissertation of A. Okojunu (Florida State University, 2026), published at https://github.com/abelokoj/LSA_Rayleigh_Fortran. The evaluation of the relation, the maps and the fits were written for this paper. The confined relation itself is that of Mikami and Mason (1975), Georgiou et al. (1991) and Rykner et al. (2024).
