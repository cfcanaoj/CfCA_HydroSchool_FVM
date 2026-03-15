import os
import numpy as np
import matplotlib.pyplot as plt
import argparse


def read_data(filetype, dirname, step):
    fileroot = os.path.join(dirname, f"snap{step:05d}")

    if filetype == "ascii":
        file = fileroot + ".dat"
        print("reading", file)

        with open(file, "r") as f:
            time = float(f.readline().split()[-1])
            nx, ny = map(int, f.readline().split()[-2:])

        arr = np.loadtxt(file, comments="#").reshape(ny, nx, -1)

        x = arr[0, :, 0]
        y = arr[:, 0, 1]

        names = ["rho", "vx", "vy", "vz", "pre", "Bx", "By", "Bz", "vorticity", "Az"]
        data_dict = {name: arr[:, :, i + 2] for i, name in enumerate(names)}

    elif filetype == "binary":
        file = fileroot + ".bin"
        print("reading", file)

        with open(file, "rb") as fp:
            time = np.fromfile(fp, np.float64, 1)[0]
            nx   = np.fromfile(fp, np.int32, 1)[0]
            ny   = np.fromfile(fp, np.int32, 1)[0]
            nhyd = np.fromfile(fp, np.int32, 1)[0]
            nbc  = np.fromfile(fp, np.int32, 1)[0]

            x = np.fromfile(fp, np.float64, nx)
            y = np.fromfile(fp, np.float64, ny)
            Q = np.fromfile(fp, np.float32, nx * ny * nhyd).reshape(ny, nx, nhyd)
            Bc = np.fromfile(fp, np.float32, nx * ny * nbc).reshape(ny, nx, nbc)
            vor = np.fromfile(fp, np.float32, nx * ny).reshape(ny, nx)
            Az  = np.fromfile(fp, np.float32, nx * ny).reshape(ny, nx)


        q_names = ["rho", "vx", "vy", "vz", "pre"]
        b_names = ["Bx", "By", "Bz"]

        data_dict = {name: Q[:, :, i] for i, name in enumerate(q_names)}
        data_dict.update({name: Bc[:, :, i] for i, name in enumerate(b_names)})
        data_dict["vorticity"] = vor
        data_dict["Az"] = Az - np.mean(Az)

    else:
        raise ValueError("filetype should be 'ascii' or 'binary'")

    return x, y, time, data_dict


def get_uniform_spacing(coord):
    d = np.diff(coord)
    d0 = d[0]
    if not np.allclose(d, d0, rtol=1e-10, atol=1e-12):
        raise ValueError("Grid is not uniform. This script assumes uniform spacing.")
    return d0


def make_kgrid(x, y):
    nx = x.size
    ny = y.size

    dx = get_uniform_spacing(x)
    dy = get_uniform_spacing(y)

    # cycles per unit length
    kx = np.fft.fftfreq(nx, d=dx)
    ky = np.fft.fftfreq(ny, d=dy)

    kx2d, ky2d = np.meshgrid(kx, ky, indexing="xy")
    kmag = np.sqrt(kx2d**2 + ky2d**2)

    return kx, ky, kx2d, ky2d, kmag, dx, dy


def shell_integrate(P2d, kmag, kbins):
    index = np.digitize(kmag.ravel(), kbins)
    Epower = np.zeros(len(kbins) - 1)

    flat = P2d.ravel()
    for n in range(1, len(kbins)):
        mask = (index == n)
        Epower[n - 1] = np.sum(flat[mask])

    return Epower


def fft2_norm(A):
    ny, nx = A.shape
    return np.fft.fft2(A) / (nx * ny)


def compute_spectrum(specname, data, kx2d, ky2d, kmag, kbins):
    if specname == "enstrophy":
        omegak = fft2_norm(data["vorticity"])
        P2d = 0.5 * np.abs(omegak)**2
        Epower = shell_integrate(P2d, kmag, kbins)

    elif specname == "kinene":
        vxk = fft2_norm(data["vx"])
        vyk = fft2_norm(data["vy"])
        P2d = 0.5 * (np.abs(vxk)**2 + np.abs(vyk)**2)
        Epower = shell_integrate(P2d, kmag, kbins)

    elif specname == "magene":
        Bxk = fft2_norm(data["Bx"])
        Byk = fft2_norm(data["By"])
        P2d = 0.5 * (np.abs(Bxk)**2 + np.abs(Byk)**2)
        Epower = shell_integrate(P2d, kmag, kbins)

    elif specname == "totene":
        vxk = fft2_norm(data["vx"])
        vyk = fft2_norm(data["vy"])
        Bxk = fft2_norm(data["Bx"])
        Byk = fft2_norm(data["By"])
        P2d = 0.5 * (np.abs(vxk)**2 + np.abs(vyk)**2 + np.abs(Bxk)**2 + np.abs(Byk)**2)
        Epower = shell_integrate(P2d, kmag, kbins)

    elif specname == "crosshelicity":
        vxk = fft2_norm(data["vx"])
        vyk = fft2_norm(data["vy"])
        Bxk = fft2_norm(data["Bx"])
        Byk = fft2_norm(data["By"])
        P2d = np.real(vxk * np.conjugate(Bxk) + vyk * np.conjugate(Byk))
        Epower = shell_integrate(P2d, kmag, kbins)

    elif specname == "magpot":
        # mean-square magnetic potential spectrum
        Azk = fft2_norm(data["Az"])
        P2d = 0.5 * np.abs(Azk)**2
        Epower = shell_integrate(P2d, kmag, kbins)

    else:
        raise KeyError(
            f"Unknown spectrum type: {specname}. "
            f"Available: enstrophy, kinene, magene, totene, crosshelicity, magpot"
        )

    return Epower


parser = argparse.ArgumentParser(
    description="Compare multiple spectrum types in each panel, one panel per directory.",
    usage="python3 PowerSpectrum_new2.py [ascii|binary] [step] [spec1] [spec2] ... --dirs [dir1] [dir2] ... [-h]",
    epilog=(
        "Example:\n"
        "  python3 PowerSpectrum_new2.py ascii 20 kinene magene totene --dirs ct hdc\n"
    ),
    formatter_class=argparse.RawTextHelpFormatter
)
parser.add_argument("filetype", choices=["ascii", "binary"], help="input file format")
parser.add_argument("step", type=int, help="step number")
parser.add_argument(
    "specnames",
    nargs="+",
    help=("spectrum types to plot "
          "(enstrophy, kinene, magene, totene, crosshelicity, magpot)")
)
parser.add_argument("--dirs", nargs="+", required=True,
                    help="directories containing snapshot files")
args = parser.parse_args()

filetype = args.filetype
step = args.step
specnames = args.specnames
dirnames = args.dirs

speclabel_dict = {
    "enstrophy": r"enstrophy$/k^2$",
    "kinene": "kinetic energy",
    "magene": "magnetic energy",
    "totene": "total energy",
    "crosshelicity": "cross helicity",
    "magpot": r"magnetic potential$/\Delta x^2$",
}

ndir = len(dirnames)

fig, axes = plt.subplots(
    1, ndir,
    figsize=(5.5 * ndir, 4.8),
    squeeze=False,
    sharex=True,
    sharey=True
)
axes = axes[0]

for idir, dirname in enumerate(dirnames):
    ax = axes[idir]

    x, y, time, data = read_data(filetype, dirname, step)

    kx, ky, kx2d, ky2d, kmag, dx, dy = make_kgrid(x, y)

    kmin_x = 1.0 / (x.size * dx)
    kmin_y = 1.0 / (y.size * dy)
    kmin = min(kmin_x, kmin_y)

    kmax_x = np.max(np.abs(kx))
    kmax_y = np.max(np.abs(ky))
    kmax = min(kmax_x, kmax_y)

    kbins = np.arange(kmin, kmax + kmin, kmin)
    kbinc = 0.5 * (kbins[1:] + kbins[:-1])

    for specname in specnames:
        Epower = compute_spectrum(specname, data, kx2d, ky2d, kmag, kbins)

        Eplot = Epower.copy()

#        if specname == "enstrophy":
#            Eplot *= dx**2
        if specname == "enstrophy": 
            Eplot = Eplot / (kbinc**2)

        if specname == "magpot":
            Eplot *= 1.0/(dx**2)

        label = speclabel_dict.get(specname, specname)

        if specname == "crosshelicity":
            valid = Eplot != 0.0
            ax.plot(kbinc[valid], np.abs(Eplot[valid]), "-", mfc="none", label=label)
        else:
            valid = Eplot > 0.0
            ax.plot(kbinc[valid], Eplot[valid], "-", mfc="none", label=label)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$k/(2\pi)$")

    if idir == 0:
        ax.set_ylabel("power spectrum")

    ax.axvline(10.0 * np.sqrt(2.0), linestyle="--", color="black")
    ax.text(10.0 * np.sqrt(2.0), 0.05, "injection scale",
            transform=ax.get_xaxis_transform(), ha="center", va="bottom")

    title = os.path.basename(os.path.normpath(dirname))
    ax.set_title(rf"{title}, time = {time:.2f}")
    ax.legend(loc="upper right")

outroot = "compare_spec_" + "_".join(specnames)
for dirname in dirnames:
    outroot += "_" + os.path.basename(os.path.normpath(dirname))

for format_fig in ["pdf", "png"]:
    outdir = format_fig + "file"
    os.makedirs(outdir, exist_ok=True)

    outputfile = os.path.join(outdir, f"{outroot}_snap{step:05d}.{format_fig}")
    print("making plot file", outputfile)
    plt.savefig(outputfile, bbox_inches="tight")

plt.show()
