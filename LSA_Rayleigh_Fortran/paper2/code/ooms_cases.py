"""Wall correction of the confined relation for the fluids-at-rest cases of Ooms et al. 2026 (Table 1), in capillary units."""
import sys, math
sys.path.insert(0, '.')
from confined_dispersion import growth_rate
R2, rho1, rho2, sig, nu2 = 0.00475547, 905.0, 1.0, 0.01, 1e-6
for case, R1, nu1, lam in ((1, 0.001, 6.63e-4, 0.0116), (2, 0.0015, 6.63e-4, 0.0116), (7, 0.001, 6.63e-6, 0.0116)):
    sc = math.sqrt(rho1 * sig * R1)                      # viscosity scale
    oh, mu2, k, Rw = rho1 * nu1 / sc, rho2 * nu2 / sc, 2 * math.pi * R1 / lam, R2 / R1
    s_c = growth_rate(k, oh, rho2 / rho1, mu2, Rw); s_u = growth_rate(k, oh, rho2 / rho1, mu2)
    tsig = math.sqrt(rho1 * R1**3 / sig)
    print(f"case {case}: Oh={oh:.4g} k={k:.3f} Rw/R1={Rw:.3f} s_unconf={s_u/tsig:.4f} 1/s s_conf={s_c/tsig:.4f} 1/s wall effect={100*(s_c/s_u-1):+.3g}%")
