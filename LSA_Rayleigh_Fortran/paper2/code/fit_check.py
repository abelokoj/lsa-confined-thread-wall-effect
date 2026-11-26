"""Reduction of the growth rate when the wall is placed at the gap given by the fit (eq. fit), for the eighty cases
of the two lightest density ratios. Shows that the fit estimates the gap but does not bound the reduction.
Run from the paper folder: python code/fit_check.py  (writes data/fit_check.txt)."""
import os, sys, csv, math
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from confined_dispersion import growth_rate

DATA = os.path.join(HERE, "..", "data")
fits = {}
for r in csv.DictReader(open(os.path.join(DATA, "wall_collapse_fit.csv"))):
    if r["rho2"] == "0.0001225+0.001225":
        fits[r["threshold"]] = [float(r[c]) for c in ("c3", "c2", "c1", "c0")]
rows = [dict(r, rho2="0.001225") for r in csv.DictReader(open(os.path.join(DATA, "wall_collapse.csv")))]
rows += [r for r in csv.DictReader(open(os.path.join(DATA, "wall_density.csv"))) if float(r["rho2"]) == 1.225e-4]
assert len(rows) == 80, len(rows)

out = [f"{len(rows)} cases (rho2/rho1 = 1.225e-3 and 1.225e-4); reduction (%) of the rate at x_m with the wall at the fitted gap"]
for key, nominal in (("gap_1pct", 1.0), ("gap_0p1pct", 0.1)):
    c, red = fits[key], []
    for r in rows:
        rho, oh, mr, x = float(r["rho2"]), float(r["Oh"]), float(r["mu_ratio"]), float(r["x_m"])
        L = math.log(float(r["diffusion_length"]) * math.sqrt(rho))
        G = math.exp(((c[0] * L + c[1]) * L + c[2]) * L + c[3])
        s0 = growth_rate(x, oh, rho, mr * oh, None)
        red.append((100 * (1 - growth_rate(x, oh, rho, mr * oh, 1 + G) / s0), rho, oh, mr))
    red.sort()
    lo, hi = red[0], red[-1]
    out.append(f"{key} (nominal {nominal}%): min {lo[0]:.4f} at rho2={lo[1]}, Oh={lo[2]}, mu2/mu1={lo[3]}; "
               f"max {hi[0]:.4f} at rho2={hi[1]}, Oh={hi[2]}, mu2/mu1={hi[3]}; above nominal in {sum(v[0] > nominal for v in red)} cases")
open(os.path.join(DATA, "fit_check.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
