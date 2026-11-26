"""Conservative wall-gap rule: the pooled fit of the threshold gap (eq. fit, two lightest density ratios) times a
safety factor F, one per threshold (1% and 0.1%). F is the smallest value, to 3 significant figures, at which the
wall at F x G_fit reduces the rate at x_m by no more than the nominal tolerance in all 80 calibration cases.
The rule is then tested on held-out cases drawn log-uniformly between the grid values (fixed seed); F is not tuned
to them. Run from the paper folder: python code/fit_bound.py
Writes data/fit_bound_calibration.csv, data/fit_bound_holdout.csv and data/fit_bound.txt."""
import os, sys, csv, math, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from confined_dispersion import growth_rate
import wall_maps
from wall_maps import most_unstable, wall_radius

DATA = os.path.join(HERE, "..", "data")
SEED, N_DRAW = 20261123, 50
DELTA_LO, DELTA_HI = 0.0086, 1.07                       # validity range of the fit stated in the paper
THRESH = (("gap_1pct", 1.0, 1e-2), ("gap_0p1pct", 0.1, 1e-3))   # (column, nominal %, tol for wall_radius)
t0 = time.time()

fits = {}
for r in csv.DictReader(open(os.path.join(DATA, "wall_collapse_fit.csv"))):
    if r["rho2"] == "0.0001225+0.001225":
        fits[r["threshold"]] = [float(r[c]) for c in ("c3", "c2", "c1", "c0")]


def g_fit(key, delta):
    c, L = fits[key], math.log(delta)
    return math.exp(((c[0] * L + c[1]) * L + c[2]) * L + c[3])


def reduction(case, gap):
    """Percent reduction of the rate at x_m with the wall at r = 1 + gap (as in fit_check.py)."""
    x, oh, rho, mu2 = case["x_m"], case["Oh"], case["rho2"], case["mu_ratio"] * case["Oh"]
    return 100 * (1 - growth_rate(x, oh, rho, mu2, 1 + gap) / growth_rate(x, oh, rho, mu2, None))


# calibration: the same 80 cases and the same x_m values that fit_check.py uses
calib = [dict(r, rho2="0.001225") for r in csv.DictReader(open(os.path.join(DATA, "wall_collapse.csv")))]
calib += [r for r in csv.DictReader(open(os.path.join(DATA, "wall_density.csv"))) if float(r["rho2"]) == 1.225e-4]
calib = [{"rho2": float(r["rho2"]), "Oh": float(r["Oh"]), "mu_ratio": float(r["mu_ratio"]), "x_m": float(r["x_m"]),
          "delta_mu": float(r["diffusion_length"]) * math.sqrt(float(r["rho2"])),
          "gap_1pct": float(r["gap_1pct"]), "gap_0p1pct": float(r["gap_0p1pct"])} for r in calib]
assert len(calib) == 80, len(calib)

F, out = {}, []
for key, nom, _ in THRESH:
    worst = lambda f: max(reduction(c, f * g_fit(key, c["delta_mu"])) for c in calib)
    w1 = worst(1.0)
    ref = {"gap_1pct": 1.2894, "gap_0p1pct": 0.1245}[key]
    assert round(w1, 4) == ref, f"fit_check.py value not reproduced for {key}: {w1:.4f} vs {ref}"
    lo, hi = 1.0, 2.0
    assert worst(hi) <= nom
    while hi - lo > 1e-5:                                  # the reduction falls monotonically as the wall recedes
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if worst(mid) <= nom else (mid, hi)
    f = math.ceil(hi * 100 - 1e-9) / 100                   # 3 significant figures for 1 <= F < 10, rounded up
    assert 1 <= f < 10 and worst(f) <= nom and worst(f - 1e-2) > nom, "F is not the smallest 3-s.f. value"
    F[key] = f

for c in calib:
    for key, nom, _ in THRESH:
        g = g_fit(key, c["delta_mu"])
        c[key + "_fit"], c[key + "_FG"] = g, F[key] * g
        c[key + "_red_fit"], c[key + "_red_FG"] = reduction(c, g), reduction(c, F[key] * g)
    for key, nom, _ in THRESH:
        assert c[key + "_red_FG"] <= nom

# held-out cases: log-uniform draws inside the fitted box, kept only inside the stated delta_mu range
rng = np.random.default_rng(SEED)
draws = np.exp(rng.uniform(np.log([0.005, 0.005, 1.225e-4]), np.log([0.1, 2.0, 1.225e-3]), size=(N_DRAW, 3)))
hold, discarded = [], []
for oh, m, rho2 in draws:
    wall_maps.RHO2 = rho2                                  # the helpers read the module-level density ratio
    mu2 = m * oh
    km, sm = most_unstable(oh, mu2)
    c = {"rho2": rho2, "Oh": oh, "mu_ratio": m, "x_m": km, "s_m": sm, "delta_mu": math.sqrt(mu2 / sm)}
    if not DELTA_LO <= c["delta_mu"] <= DELTA_HI:
        discarded.append(c); continue
    for key, nom, tol in THRESH:
        g = g_fit(key, c["delta_mu"])
        c[key] = wall_radius(tol, km, oh, mu2, sm) - 1
        c[key + "_fit"], c[key + "_FG"] = g, F[key] * g
        c[key + "_red_fit"], c[key + "_red_FG"] = reduction(c, g), reduction(c, F[key] * g)
    hold.append(c)
assert len(hold) >= 30, f"only {len(hold)} held-out cases inside the delta_mu range"

cols = ["rho2", "Oh", "mu_ratio", "x_m"]
for key, _, _ in THRESH:
    cols += [key + "_exact", key + "_fit", key + "_FG", key + "_red_fit", key + "_red_FG"]
for name, cases, extra in (("calibration", calib, []), ("holdout", hold, ["delta_mu"])):
    with open(os.path.join(DATA, f"fit_bound_{name}.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols[:4] + extra + cols[4:])
        for c in cases:
            vals = [c[k] for k in cols[:4] + extra] + [c[k.replace("_exact", "")] if k.endswith("_exact") else c[k] for k in cols[4:]]
            w.writerow([f"{v:.6g}" for v in vals])


def where(c):
    return f"rho2={c['rho2']:.4g}, Oh={c['Oh']:.4g}, mu2/mu1={c['mu_ratio']:.4g}"


out.append(f"Calibration: 80 cases (rho2/rho1 = 1.225e-4 and 1.225e-3), pooled fit; reductions (%) of the rate at x_m.")
out.append(f"Held-out: seed {SEED}, {N_DRAW} log-uniform draws (Oh 0.005-0.1, mu2/mu1 0.005-2, rho2 1.225e-4-1.225e-3); "
           f"{len(hold)} used, {len(discarded)} discarded for delta_mu outside [{DELTA_LO}, {DELTA_HI}].")
for key, nom, _ in THRESH:
    out.append(f"\n{key} (nominal {nom}%): F = {F[key]:.2f}")
    for name, cases in (("calibration", calib), ("held-out", hold)):
        wf = max(cases, key=lambda c: c[key + "_red_fit"])
        wF = max(cases, key=lambda c: c[key + "_red_FG"])
        widen = np.array([c[key + "_FG"] / c[key] for c in cases])
        over = [c for c in cases if c[key + "_red_FG"] > nom]
        out.append(f"  {name}: worst at G_fit {wf[key + '_red_fit']:.4f} ({where(wf)}); worst at F*G_fit "
                   f"{wF[key + '_red_FG']:.4f} ({where(wF)}); F*G_fit/G_exact median {np.median(widen):.3f}, "
                   f"min {widen.min():.3f}, max {widen.max():.3f}; above nominal at F*G_fit: {len(over)}")
        for c in over:
            out.append(f"    exceeds: {c[key + '_red_FG']:.4f}% at {where(c)}, delta_mu={c['delta_mu']:.4g}")
out.append(f"\nRuntime {time.time() - t0:.0f} s")
open(os.path.join(DATA, "fit_bound.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
