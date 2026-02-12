#!/usr/bin/env python3
"""
Analytic solution for 1D Gaussian pulse diffusion:
    E(x,t) = (E0*sigma0)/sqrt(sigma0^2 + 2 D t) * exp(-(x-x0)^2 / (2 (sigma0^2 + 2 D t)))

with diffusion coefficient (radiation diffusion limit):
    D = c / (3 * rho * kappa_mass)   (if kappa is mass opacity [cm^2/g])
or
    D = c / (3 * kappa_vol)          (if kappa is volumetric opacity [1/cm])

This script writes 100 snapshots from t=0..30 into ./ana/
"""

import os
import numpy as np

# -----------------------
# User settings (edit here)
# -----------------------
outdir = "ana"

c = 3.0e10        # light speed (keep explicit; set to your code value)

# time
tmax = 30.0/c
nsnap = 100  # includes t=0 and t=tmax

# space grid (choose to match your simulation)
xmin, xmax = -1.0, 1.0
nx = 200

# initial Gaussian parameters
E0 = 1.0e20       # peak amplitude at t=0 (for E background=0)
x0 = 0.0          # center
sigma0 = 0.10  # initial width

# physical parameters (use ONE of the two blocks below)

# (A) kappa as mass opacity [cm^2/g]
rho = 1.0
kappa_mass = 1.0e3
D = c / (3.0 * rho * kappa_mass)

# (B) kappa as volumetric opacity [1/cm]
# c = 3.0e10
# kappa_vol = 1.0e-2
# D = c / (3.0 * kappa_vol)

# background (if your sim has non-zero background, add it here)
E_bg = 0.0

# -----------------------
# Do not edit below
# -----------------------
def EF_analytic(x: np.ndarray, t: float) -> np.ndarray, np.ndarray:
    sig2 = sigma0**2 + 2.0 * D * t
    pref = E0 * sigma0 / np.sqrt(sig2)
    E = E_bg + pref * np.exp(- (x - x0)**2 / (2.0 * sig2))
    F = D*(x - x0)/sig2 * E
    return E,F

def main():
    os.makedirs(outdir, exist_ok=True)
    x = np.linspace(xmin, xmax, nx)

    times = np.linspace(0.0, tmax, nsnap)
    for n, t in enumerate(times):
        E,F = EF_analytic(x, t)
        fn = os.path.join(outdir, f"snap{n:05d}.dat")
        # columns: x, E, t
        data = np.column_stack([x, E, F, np.full_like(x, t)])
        header = (
            f"# time= {t}\n",
            f"#   nx= {len(x)}\n",
             "# x E Fx \n"
        )
        np.savetxt(fn, data, header=header)

    print(f"Wrote {nsnap} snapshots to ./{outdir}/ (t=0..{tmax})")
    print(f"Using D = {D:g} (with c = {c:g})")

if __name__ == "__main__":
    main()
