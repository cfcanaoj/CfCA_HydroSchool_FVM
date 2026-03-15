import argparse
import csv
import sys
import numpy as np


def exact_riemann_solver(
    gam,
    rhoL,
    velL,
    preL,
    rhoR,
    velR,
    preR,
    tol=1.0e-8,
    max_iter=100,
    p_floor=1.0e-8,
):
    """
    Solve the star region pressure and velocity for the 1D Riemann problem
    in Lagrangian coordinates.

    Returns
    -------
    prest, velst : float
        Pressure and velocity in the star region.
    """
    gamp1 = gam + 1.0
    gamm1 = gam - 1.0

    # Lagrangian sound speed
    cL = np.sqrt(gam*preL*rhoL)
    cR = np.sqrt(gam*preR*rhoR)

    testP = (cR*preL + cL*preR - cR*cL*(velR - velL))/(cR + cL)
    testP = max(testP, p_floor)

    hantei = 1.0e30
    zR = 0.0
    zL = 0.0
    vR = 0.0
    vL = 0.0

    for _ in range(max_iter):
        if testP >= preL:  # shock
            wsqL = 0.5*cL**2*(gamp1*testP/preL + gamm1)/gam
            wL = np.sqrt(wsqL)
            zL = 2.0*wsqL*wL/(wsqL + cL**2)
        else:  # rarefaction
            wL = (gamm1/(2.0*gam)*(1.0 - testP/preL)/(1.0 - (testP/preL)**(gamm1/(2.0*gam)))*cL)
            zL = cL*(testP/preL)**(1.0 - gamm1/(2.0*gam))

        if testP >= preR:  # shock
            wsqR = 0.5*cR**2*(gamp1*testP/preR + gamm1)/gam
            wR = np.sqrt(wsqR)
            zR = 2.0*wsqR*wR/(wsqR + cR**2)
        else:  # rarefaction
            wR = (gamm1/(2.0*gam)*(1.0 - testP/preR)/(1.0 - (testP/preR)**(gamm1/(2.0*gam)))*cR)
            zR = cR*(testP/preR)**(1.0 - gamm1/(2.0*gam))

        vR = velR + (testP - preR)/wR
        vL = velL - (testP - preL)/wL

        hantei = zR*zL*(vR - vL)/((zR + zL)*testP)
        testP = max(testP*(1.0 - hantei), p_floor)

        if abs(hantei) <= tol:
            velst = (zL*vL + zR*vR)/(zL + zR)
            return testP, velst

    raise RuntimeError(
        f"exact_riemann_solver did not converge within {max_iter} iterations"
    )


def exact_riemann_solution(
    time,
    x0,
    gam,
    rhoL,
    velL,
    preL,
    rhoR,
    velR,
    preR,
    n_rarefac=64,
    pad_frac=0.1,
    tol=1.0e-4,
    max_iter=100,
    p_floor=1.0e-8,
):
    """
    Build piecewise points for plotting/output of the exact Riemann solution.
    """
    if gam <= 1.0:
        raise ValueError("gamma should be larger than 1")
    if rhoL <= 0.0 or rhoR <= 0.0:
        raise ValueError("density should be positive")
    if preL <= 0.0 or preR <= 0.0:
        raise ValueError("pressure should be positive")
    if time < 0.0:
        raise ValueError("time should be non-negative")
    if n_rarefac < 2:
        raise ValueError("n_rarefac should be >= 2")
    if time == 0.0:
        raise ValueError("time should be >= 0")

    gamp1 = gam + 1.0
    gamm1 = gam - 1.0

    prest, velst = exact_riemann_solver(
        gam,
        rhoL,
        velL,
        preL,
        rhoR,
        velR,
        preR,
        tol=tol,
        max_iter=max_iter,
        p_floor=p_floor,
    )

    # Eulerian sound speed for wave positions
    aL = np.sqrt(gam*preL/rhoL)
    aR = np.sqrt(gam*preR/rhoR)

    is_shock_L = prest > preL
    is_shock_R = prest > preR

    if is_shock_L:
        rhostL = rhoL*(prest/preL + gamm1/gamp1)/(gamm1/gamp1*prest/preL + 1.0)
        sL = velL - aL*np.sqrt((gamp1*prest/preL + gamm1)/(2.0*gam))
    else:
        rhostL = rhoL*(prest/preL)**(1.0/gam)
        sL_head = velL - aL
        sL_tail = velst - np.sqrt(gam*prest/rhostL)

    if is_shock_R:
        rhostR = rhoR*(prest/preR + gamm1/gamp1)/(gamm1/gamp1*prest/preR + 1.0)
        sR = velR + aR*np.sqrt((gamp1*prest/preR + gamm1)/(2.0*gam))
    else:
        rhostR = rhoR*(prest/preR)**(1.0/gam)
        sR_head = velR + aR
        sR_tail = velst + np.sqrt(gam*prest/rhostR)

    x_riem = []
    rho_riem = []
    vel_riem = []
    pre_riem = []

    if is_shock_L:
        x_riem += [x0 + sL*time, x0 + sL*time]
        rho_riem += [rhoL, rhostL]
        vel_riem += [velL, velst]
        pre_riem += [preL, prest]
    else:
        x_riem.append(x0 + sL_head*time)
        rho_riem.append(rhoL)
        vel_riem.append(velL)
        pre_riem.append(preL)

        dx = -(sL_head - sL_tail)*time/(n_rarefac - 1.0)
        for i in range(n_rarefac):
            xi = sL_head*time + i*dx
            x_riem.append(x0 + xi)
            rho_riem.append( rhoL*(2.0/gamp1 + gamm1/(gamp1*aL)*(velL - xi/time))**(2.0/gamm1))
            vel_riem.append(2.0/gamp1*(aL + gamm1/2.0*velL + xi/time))
            pre_riem.append(preL*(2.0/gamp1 + gamm1/(gamp1*aL)*(velL - xi/time))**(2.0*gam/gamm1))

    # contact discontinuity
    x_riem += [x0 + velst*time, x0 + velst*time]
    rho_riem += [rhostL, rhostR]
    vel_riem += [velst, velst]
    pre_riem += [prest, prest]

    if is_shock_R:
        x_riem += [x0 + sR*time, x0 + sR*time]
        rho_riem += [rhostR, rhoR]
        vel_riem += [velst, velR]
        pre_riem += [prest, preR]
    else:
        dx = (sR_head - sR_tail)*time/(n_rarefac - 1.0)
        for i in range(n_rarefac):
            xi = sR_tail*time + i*dx
            x_riem.append(x0 + xi)
            rho_riem.append(rhoR*(2.0/gamp1 - gamm1/(gamp1*aR)*(velR - xi/time))**(2.0/gamm1))
            vel_riem.append(2.0/gamp1*(-aR + gamm1/2.0*velR + xi/time))
            pre_riem.append(preR*(2.0/gamp1 - gamm1/(gamp1*aR)*(velR - xi/time))**(2.0*gam/gamm1))

        x_riem.append(x0 + sR_head*time)
        rho_riem.append(rhoR)
        vel_riem.append(velR)
        pre_riem.append(preR)

    # Automatically determine outer range
    xmin = min(x_riem)
    xmax = max(x_riem)
    width = xmax - xmin
    pad = pad_frac * width if width > 0.0 else 0.5

    # Append outer constant states
    x_riem = [xmin - pad] + x_riem + [xmax + pad]
    rho_riem = [rhoL] + rho_riem + [rhoR]
    vel_riem = [velL] + vel_riem + [velR]
    pre_riem = [preL] + pre_riem + [preR]   

    return x_riem, rho_riem, vel_riem, pre_riem


def get_case(case_name):
    """
    Return preset initial conditions.
    """
    case_name = case_name.lower()

    if case_name == "sod":
        return dict(
            gam=1.4,
            time=0.2,
            x0=0.0,
            rhoL=1.0,
            velL=0.0,
            preL=1.0,
            rhoR=0.125,
            velR=0.0,
            preR=0.1,
        )

    if case_name == "lax":
        return dict(
            gam=1.4,
            time=0.16,
            x0=0.0,
            rhoL=0.445,
            velL=0.698,
            preL=3.528,
            rhoR=0.5,
            velR=0.0,
            preR=0.571,
        )

    if case_name == "strong_shock":
        return dict(
            gam=1.4,
            time=0.012,
            x0=0.0,
            rhoL=1.0,
            velL=0.0,
            preL=1000.0,
            rhoR=1.0,
            velR=0.0,
            preR=0.01,
        )

    raise ValueError(f"unknown case: {case_name}")


def write_output(filename, x, rho, vel, pre, params):
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f, delimiter=" ")
        print(f"#{'rhoL':<5} = {params['rhoL']:<10.6f} | {'rhoR':<5} = {params['rhoR']:<10.6f}", file=f)
        print(f"#{'velL':<5} = {params['velL']:<10.6f} | {'velR':<5} = {params['velR']:<10.6f}", file=f)
        print(f"#{'preL':<5} = {params['preL']:<10.6f} | {'preR':<5} = {params['preR']:<10.6f}", file=f)
        print(f"# gamma = {params['gam']}", file=f)
        print(f"# time  = {params['time']}", file=f)
        print(f"# x0    = {params['x0']}", file=f)
        print("# x rho vel pre", file=f)
        writer.writerow(["#x", "rho", "vel", "pre"])
        for xi, rhoi, veli, prei in zip(x, rho, vel, pre):
            writer.writerow(
                [
                    f"{xi:.16e}",
                    f"{rhoi:.16e}",
                    f"{veli:.16e}",
                    f"{prei:.16e}",
                ]
            )

def parse_args():
    parser = argparse.ArgumentParser(
        description="Exact Riemann solver in Lagrangian coordinates",
        usage=(
        "%(prog)s --case custom --gamma GAMMA --time TIME --x0 X0 "
        "--rhoL RHOL --velL VELL --preL PREL "
        "--rhoR RHOR --velR VELR --preR PRER "
        "[--n-rarefac N_RAREFAC] [--pad_frac PAD_FRAC] "
        "[-o OUTFILE] "
        "[--tol TOL] [--max-iter MAX_ITER] [--p-floor P_FLOOR]\n"
        "       %(prog)s [--case {sod,lax,strong_shock}] "
        "[--time TIME] [--x0 X0] [-o OUTFILE]"
    ),
    )

    parser.add_argument(
        "--case",
        type=str,
        default="sod",
        choices=["sod", "lax", "strong_shock", "custom"],
        help="preset problem name or custom",
    )

    parser.add_argument("--gamma", type=float, default=None, help="ratio of specific heats")
    parser.add_argument("--time", type=float, default=None, help="time at which solution is sampled")
    parser.add_argument("--x0", type=float, default=None, help="initial discontinuity position")

    parser.add_argument("--rhoL", type=float, default=None, help="left density")
    parser.add_argument("--velL", type=float, default=None, help="left velocity")
    parser.add_argument("--preL", type=float, default=None, help="left pressure")

    parser.add_argument("--rhoR", type=float, default=None, help="right density")
    parser.add_argument("--velR", type=float, default=None, help="right velocity")
    parser.add_argument("--preR", type=float, default=None, help="right pressure")

    parser.add_argument(
        "--n-rarefac",
        type=int,
        default=64,
        help="number of sampling points inside each rarefaction fan",
    )
    parser.add_argument(
        "--pad_frac",
        type=float,
        default=0.1,
        help="extra x-range added to the left/right outside the waves",
    )

    parser.add_argument(
        "-o",
        "--outfile",
        type=str,
        default="riemann_exact.dat",
        help="output filename",
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=1.0e-4,
        help="iteration tolerance for star pressure solver",
    )
    parser.add_argument(
        "--max-iter",
        type=int,
        default=100,
        help="maximum number of iterations for star pressure solver",
    )
    parser.add_argument(
        "--p-floor",
        type=float,
        default=1.0e-8,
        help="minimum allowed pressure during iteration",
    )

    return parser.parse_args()


def merge_case_and_args(args):
    if args.case == "custom":
        required = ["gamma", "time", "x0", "rhoL", "velL", "preL", "rhoR", "velR", "preR"]
        missing = [name for name in required if getattr(args, name) is None]
        if missing:
            raise ValueError(
                "for --case custom, the following arguments are required: "
                + ", ".join("--" + m for m in missing)
            )

        params = dict(
            gam=args.gamma,
            time=args.time,
            x0=args.x0,
            rhoL=args.rhoL,
            velL=args.velL,
            preL=args.preL,
            rhoR=args.rhoR,
            velR=args.velR,
            preR=args.preR,
        )
    else:
        params = get_case(args.case)

        overrides = {
            "gam": args.gamma,
            "time": args.time,
            "x0": args.x0,
            "rhoL": args.rhoL,
            "velL": args.velL,
            "preL": args.preL,
            "rhoR": args.rhoR,
            "velR": args.velR,
            "preR": args.preR,
        }
        for key, value in overrides.items():
            if value is not None:
                params[key] = value

    return params

def print_initial_states(params, stream=sys.stdout):
    print("=== adopted initial states ===", file=stream)
    print(f"{'rhoL':<5} = {params['rhoL']:<10.6f} | {'rhoR':<5} = {params['rhoR']:<10.6f}", file=stream)
    print(f"{'velL':<5} = {params['velL']:<10.6f} | {'velR':<5} = {params['velR']:<10.6f}", file=stream)
    print(f"{'preL':<5} = {params['preL']:<10.6f} | {'preR':<5} = {params['preR']:<10.6f}", file=stream)

def main():
    args = parse_args()
    params = merge_case_and_args(args)

    print_initial_states(params)
    print(f"\ngamma = {params['gam']:.6f}")
    print(f"time  = {params['time']:.6f}")
    print(f"x0    = {params['x0']:.6f}")
    print(f"outputfile = {args.outfile}")

    x, rho, vel, pre = exact_riemann_solution(
        time=params["time"],
        x0=params["x0"],
        gam=params["gam"],
        rhoL=params["rhoL"],
        velL=params["velL"],
        preL=params["preL"],
        rhoR=params["rhoR"],
        velR=params["velR"],
        preR=params["preR"],
        n_rarefac=args.n_rarefac,
        pad_frac=args.pad_frac,
        tol=args.tol,
        max_iter=args.max_iter,
        p_floor=args.p_floor,
    )

    write_output(args.outfile, x, rho, vel, pre, params)


if __name__ == "__main__":
    main()
