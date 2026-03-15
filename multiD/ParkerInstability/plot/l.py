#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Half-domain (y>=0) Parker instability solver for Matsumoto+2019 background
with B-parity (midplane corrugation) AND mode tracking by vertical node count.

What it outputs:
  dispersion_m19_half_Bparity_track_nodes0_kx0to0p4_n32.txt
  columns: kx   sigma_real(selected)   nodes(vy)   sigma_imag(selected)

Selection logic at each kx:
  - Solve generalized eigenproblem A q = sigma B q (dense).
  - Keep finite eigenvalues, and |Im(sigma)| < imag_tol.
  - Compute node count from vy(y) eigenfunction (sign changes).
  - Filter by desired node count (default nodes_target=0).
  - Track branch by maximizing overlap with previous eigenvector (all variables),
    with fallback to max Re(sigma) if previous is not available.

BC and model:
  - Background: Matsumoto+2019 eq.(59)-(61) + EOS p=rho*T/gamma (gamma=1.05)
  - B0(y) from constant beta0: B0=sqrt(2 p0/beta0)
  - divB=0 enforced by Az scalar potential a(y): dBx=a', dBy=-ik a
  - Half-domain y>=0 with B-parity at y=0:
      vx(0)=0, p1(0)=0, rho1(0)=0, vy'(0)=0, a'(0)=0
  - Top boundary (y=ymax): Neumann by default

Requirements:
  numpy, scipy
"""

from __future__ import annotations
import numpy as np


# ----------------------------
# Background profiles
# ----------------------------
def g_profile(y: np.ndarray, g0: float, Hg: float) -> np.ndarray:
    return -g0 * np.tanh(y / Hg)


def T_profile_m19_half(y: np.ndarray, T0: float, T1: float, y0: float, Ht: float) -> np.ndarray:
    # half-domain: |y|=y
    return T0 + 0.5 * (T1 - T0) * (np.tanh((y - y0) / Ht) + 1.0)


def D1_matrix(N: int, dy: float) -> np.ndarray:
    """First derivative: 2nd-order centered interior; 1st-order one-sided at boundaries."""
    D = np.zeros((N, N), dtype=float)
    for i in range(1, N - 1):
        D[i, i - 1] = -0.5 / dy
        D[i, i + 1] = 0.5 / dy
    D[0, 0] = -1.0 / dy
    D[0, 1] = 1.0 / dy
    D[-1, -2] = -1.0 / dy
    D[-1, -1] = 1.0 / dy
    return D


def D2_matrix(N: int, dy: float) -> np.ndarray:
    """Second derivative: 2nd-order centered interior; 2nd-order one-sided at boundaries."""
    D2 = np.zeros((N, N), dtype=float)
    inv = 1.0 / (dy * dy)
    for i in range(1, N - 1):
        D2[i, i - 1] = inv
        D2[i, i] = -2.0 * inv
        D2[i, i + 1] = inv
    # 2nd-order one-sided at boundaries
    D2[0, 0] = 2.0 * inv
    D2[0, 1] = -5.0 * inv
    D2[0, 2] = 4.0 * inv
    D2[0, 3] = -1.0 * inv
    D2[-1, -1] = 2.0 * inv
    D2[-1, -2] = -5.0 * inv
    D2[-1, -3] = 4.0 * inv
    D2[-1, -4] = -1.0 * inv
    return D2


def build_equilibrium_half_pressure_first(
    y: np.ndarray,
    g0: float = 1.47,
    Hg: float = 5.0,
    T0: float = 1.0,
    T1: float = 25.0,
    y0: float = 10.0,
    Ht: float = 0.5,
    beta0: float = 10,
    gamma: float = 1.05,
#    gamma: float = 5.0/3.0,
#    gamma: float = 1.4,
    rho_mid: float = 1.0,
):
    """
    Integrate pressure first:
      (1+1/beta0) dp/dy = rho g, rho = gamma p / T  -> dp/dy = (gamma/fac)(g/T)p
    """
    N = len(y)
    dy = y[1] - y[0]
    D1 = D1_matrix(N, dy)

    T = T_profile_m19_half(y, T0=T0, T1=T1, y0=y0, Ht=Ht)
    g = g_profile(y, g0=g0, Hg=Hg)

    p0 = np.zeros_like(y)
    p0[0] = rho_mid * T[0] / gamma

    fac = 1.0 + 1.0 / beta0

    def rhs_p(p, Tloc, gloc):
        return (gamma / fac) * (gloc / Tloc) * p

    for i in range(0, N - 1):
        k1 = rhs_p(p0[i], T[i], g[i])
        k2 = rhs_p(p0[i] + 0.5 * dy * k1, T[i], g[i])
        k3 = rhs_p(p0[i] + 0.5 * dy * k2, T[i], g[i])
        k4 = rhs_p(p0[i] + dy * k3, T[i + 1], g[i + 1])
        p0[i + 1] = p0[i] + (dy / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

    rho0 = gamma * p0 / T
    B0 = np.sqrt(2.0 * p0 / beta0)

    dp0dy = D1 @ p0
    dB0dy = D1 @ B0

    return rho0, p0, dp0dy, T, g, B0, dB0dy


# ----------------------------
# Assemble generalized eigenproblem: A q = sigma B q (half-domain, B-parity)
# ----------------------------
def assemble_AB_half_Bparity(
    y: np.ndarray,
    kx: float,
    rho0: np.ndarray,
    p0: np.ndarray,
    dp0dy: np.ndarray,
    g: np.ndarray,
    B0: np.ndarray,
    dB0dy: np.ndarray,
    gamma: float,
    top_bc: str = "neumann",  # "neumann" or "dirichlet"
):
    N = len(y)
    dy = y[1] - y[0]
    D1 = D1_matrix(N, dy)
    D2 = D2_matrix(N, dy)

    ik = 1j * kx
    k2 = kx * kx

    # unknown blocks: [vx, vy, a, p1, rho1]
    VX, VY, AAZ, P1, R1 = 0, 1, 2, 3, 4
    nvar = 5
    dim = nvar * N

    A = np.zeros((dim, dim), dtype=complex)
    B = np.zeros((dim, dim), dtype=complex)

    def blk(i, j):
        return slice(i * N, (i + 1) * N), slice(j * N, (j + 1) * N)

    inv_rho0 = 1.0 / rho0
    drho0dy = D1 @ rho0

    # Continuity: sigma rho1 = -[ i k rho0 vx + (rho0 vy)' ]
    A[blk(R1, VX)] += -(ik) * np.diag(rho0)
    A[blk(R1, VY)] += -(np.diag(rho0) @ D1) - np.diag(drho0dy)
    B[blk(R1, R1)] += np.eye(N)

    # x-momentum: sigma vx = [ - i k p1 - i k (dB0/dy) a ] / rho0
    A[blk(VX, P1)] += -(ik) * np.diag(inv_rho0)
    A[blk(VX, AAZ)] += -(ik) * np.diag(dB0dy * inv_rho0)
    B[blk(VX, VX)] += np.eye(N)

    # y-momentum:
    # sigma vy = [ -p1' + rho1 g + B0(k^2 a - a'') - (dB0/dy) a' ] / rho0
    A[blk(VY, P1)] += -(np.diag(inv_rho0) @ D1)
    A[blk(VY, R1)] += np.diag(g * inv_rho0)
    A[blk(VY, AAZ)] += (
        np.diag(B0 * k2 * inv_rho0)
        - (np.diag(B0 * inv_rho0) @ D2)
        - (np.diag(dB0dy * inv_rho0) @ D1)
    )
    B[blk(VY, VY)] += np.eye(N)

    # Induction (Az): sigma a + B0 vy = 0
    A[blk(AAZ, VY)] += -np.diag(B0)
    B[blk(AAZ, AAZ)] += np.eye(N)

    # Pressure: sigma p1 = - vy p0' - gamma p0 (i k vx + vy')
    A[blk(P1, VY)] += -np.diag(dp0dy) - gamma * (np.diag(p0) @ D1)
    A[blk(P1, VX)] += -(gamma * ik) * np.diag(p0)
    B[blk(P1, P1)] += np.eye(N)

    # ----------------
    # Boundary conditions
    # ----------------
    def set_row_dirichlet(var, idx):
        r = var * N + idx
        A[r, :] = 0.0
        B[r, :] = 0.0
        A[r, var * N + idx] = 1.0

    def set_row_neumann_2nd(var, idx):
        r = var * N + idx
        A[r, :] = 0.0
        B[r, :] = 0.0
        if idx == 0:
            # (-3 f0 + 4 f1 - f2) = 0
            A[r, var * N + 0] = -3.0
            A[r, var * N + 1] = 4.0
            A[r, var * N + 2] = -1.0
        elif idx == N - 1:
            # (fN-3 - 4 fN-2 + 3 fN-1) = 0
            A[r, var * N + (N - 3)] = 1.0
            A[r, var * N + (N - 2)] = -4.0
            A[r, var * N + (N - 1)] = 3.0
        else:
            raise ValueError("Neumann helper only for boundaries")

    # Midplane parity for B-family:
    #   vx(0)=0, p1(0)=0, rho1(0)=0, vy'(0)=0, a'(0)=0
    set_row_dirichlet(VX, 0)
    set_row_dirichlet(P1, 0)
    set_row_dirichlet(R1, 0)
    set_row_neumann_2nd(VY, 0)
    set_row_neumann_2nd(AAZ, 0)

    # Top boundary:
    if top_bc == "neumann":
        for var in [VX, VY, AAZ, P1, R1]:
            set_row_neumann_2nd(var, N - 1)
    elif top_bc == "dirichlet":
        for var in [VX, VY, AAZ, P1, R1]:
            set_row_dirichlet(var, N - 1)
    else:
        raise ValueError("top_bc must be 'neumann' or 'dirichlet'")

    return A, B


# ----------------------------
# Node counting + tracking
# ----------------------------
def count_nodes_real(f: np.ndarray, drop: int = 2, rel_thresh: float = 1e-6) -> int:
    """
    Count sign changes of Re(f) on interior, ignoring small values.

    drop: ignore first/last 'drop' points to avoid boundary artefacts.
    rel_thresh: values smaller than rel_thresh*max(|f|) are treated as 0.
    """
    fr = np.real(f).copy()
    if fr.size <= 2 * drop + 1:
        return 0
    fr = fr[drop:-drop]

    scale = np.max(np.abs(fr))
    if scale == 0:
        return 0
    thr = rel_thresh * scale

    # map to signs with 0 for small values
    s = np.zeros_like(fr, dtype=int)
    s[fr > thr] = 1
    s[fr < -thr] = -1

    # remove zeros by forward fill then backward fill
    # (so tiny near-zeros don't create artificial nodes)
    for i in range(1, len(s)):
        if s[i] == 0:
            s[i] = s[i - 1]
    for i in range(len(s) - 2, -1, -1):
        if s[i] == 0:
            s[i] = s[i + 1]

    # count sign changes
    nodes = int(np.sum(s[1:] * s[:-1] < 0))
    return nodes


def normalize_vec(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    if n == 0:
        return v
    return v / n


def overlap_score(v_prev: np.ndarray, v: np.ndarray) -> float:
    """
    Use magnitude of complex inner product as overlap (phase-invariant).
    """
    vp = normalize_vec(v_prev)
    vc = normalize_vec(v)
    return float(np.abs(np.vdot(vp, vc)))


def solve_and_select_mode(
    A: np.ndarray,
    B: np.ndarray,
    N: int,
    nodes_target: int,
    imag_tol: float,
    prev_vec: np.ndarray | None,
    vy_rel_thresh: float = 1e-6,
):
    """
    Solve generalized eigenproblem, return selected (sigma, vec, nodes) for target nodes.
    """
    import scipy.linalg

    # eig returns right eigenvectors: A V = B V diag(w)
    w, V = scipy.linalg.eig(A, B)

    # finite only
    mask = np.isfinite(w)
    w = w[mask]
    V = V[:, mask]

    # keep (almost) real modes
    mask = np.abs(np.imag(w)) < imag_tol
    w = w[mask]
    V = V[:, mask]

    if w.size == 0:
        return None

    # compute nodes using vy block
    VX, VY, AAZ, P1, R1 = 0, 1, 2, 3, 4
    vy_block = slice(VY * N, (VY + 1) * N)

    nodes_list = np.zeros(w.size, dtype=int)
    for i in range(w.size):
        vy = V[vy_block, i]
        nodes_list[i] = count_nodes_real(vy, drop=2, rel_thresh=vy_rel_thresh)

    # filter by target nodes
    idxs = np.where(nodes_list == nodes_target)[0]
    if idxs.size == 0:
        # If no exact match, fall back to the smallest node count available
        m = int(nodes_list.min())
        idxs = np.where(nodes_list == m)[0]

    # choose by tracking or max growth
    if prev_vec is None:
        # pick max Re(sigma) within candidates
        j = idxs[np.argmax(np.real(w[idxs]))]
    else:
        # maximize overlap; tie-breaker by larger Re(sigma)
        scores = np.array([overlap_score(prev_vec, V[:, j]) for j in idxs])
        best = np.argmax(scores)
        # if multiple close, prefer larger growth
        best_idxs = idxs[np.where(scores >= scores[best] - 1e-6)[0]]
        j = best_idxs[np.argmax(np.real(w[best_idxs]))]

    return w[j], V[:, j], int(nodes_list[j])


# ----------------------------
# Main
# ----------------------------
def main():
    # ---- domain (half) ----
    # y in [0, 15*pi] as half of Matsumoto+2019 benchmark domain
    N = 201  # increase to 301 if you want tighter agreement; dense eig gets heavier
    y_max = 15.0 * np.pi
    y = np.linspace(0.0, y_max, N)

    # ---- Matsumoto+2019 parameters ----
    g0 = 1.47
    Hg = 5.0
    Ht = 0.5
    y0 = 10.0
    T0 = 1.0
    T1 = 25.0
    beta0 = 10
    gamma = 1.05
#    gamma = 5.0/3.0
#    gamma = 1.4

    rho0, p0, dp0dy, T, g, B0, dB0dy = build_equilibrium_half_pressure_first(
        y, g0=g0, Hg=Hg, T0=T0, T1=T1, y0=y0, Ht=Ht, beta0=beta0, gamma=gamma, rho_mid=1.0
    )

    # Matsumoto+2019 normalization: H0=1, so kxH0 == kx
    nk = 1
#    nk = 32*2
#    kx_list = np.linspace(0.0, 0.5, nk)
    kx_list = np.array([2*np.pi/(7.5*np.pi)])

    # boundary at top
    TOP_BC = "neumann"  # or "dirichlet"

    # tracking parameters
    nodes_target = 0       # <<< 基底モード（節0）を追跡。1にすると1節モード。
    imag_tol = 1e-7        # filter for almost-real eigenvalues
    vy_rel_thresh = 1e-6   # for node counting
    prev_vec = None

    # results
    out = np.zeros((nk, 4), dtype=float)

    for i, kx in enumerate(kx_list):
        A, Bm = assemble_AB_half_Bparity(
            y, kx, rho0, p0, dp0dy, g, B0, dB0dy, gamma=gamma, top_bc=TOP_BC
        )

        sel = solve_and_select_mode(
            A, Bm, N=N, nodes_target=nodes_target, imag_tol=imag_tol,
            prev_vec=prev_vec, vy_rel_thresh=vy_rel_thresh
        )

        if sel is None:
            sigma = 0.0 + 0.0j
            nodes = -1
            vec = None
        else:
            sigma, vec, nodes = sel

        # kx=0 clamp tiny noise
        if kx == 0.0 and np.real(sigma) < 1e-10:
            sigma = 0.0 + 0.0j

        # store
        out[i, 0] = kx
        out[i, 1] = np.real(sigma)
        out[i, 2] = float(nodes)
        out[i, 3] = np.imag(sigma)

        print(f"kx={kx:8.5f}  sigma={np.real(sigma): .6e} + i{np.imag(sigma): .2e}   nodes(vy)={nodes}")

        if vec is not None:
            prev_vec = vec  # update for tracking

    # save
#    outname = f"dispersion_m19_half_Bparity_track_nodes{nodes_target}_kx0to0p4_n32.txt"
#    outname = f"dispersion_m19_half_Bparity_track_nodes{nodes_target}_kx0to0p4_n32_gam1.4.txt"
    outname = f"test"
    header = (
        "# Matsumoto+2019 Parker dispersion (half-domain y>=0, B-parity, node-count tracking)\n"
        "# columns: kx(=kxH0)   sigma_real(selected)   nodes(vy)   sigma_imag(selected)\n"
        f"# params: g0={g0}, Hg={Hg}, T0={T0}, T1={T1}, y0={y0}, Ht={Ht}, beta0={beta0}, gamma={gamma}, N={N}, y_max={y_max}\n"
        f"# midplane parity (B-family): vx(0)=0, p1(0)=0, rho1(0)=0, vy'(0)=0, a'(0)=0\n"
        f"# top BC: {TOP_BC}\n"
        f"# tracking: nodes_target={nodes_target}, imag_tol={imag_tol}, vy_rel_thresh={vy_rel_thresh}\n"
    )
    np.savetxt(outname, out, header=header, fmt="%.10e")
    print(f"\nSaved: {outname}")


if __name__ == "__main__":
    main()
