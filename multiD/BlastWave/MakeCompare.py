import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.colors import Normalize, LogNorm
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

        names = ["rho", "vx", "vy", "vz", "pre", "Bx", "By", "Bz"]
        data_dict = {name: arr[:, :, i + 2] for i, name in enumerate(names)}

    elif filetype == "binary":
        file = fileroot + ".bin"
        print("reading", file)

        with open(file, "rb") as fp:
            time = np.fromfile(fp, np.float64, 1).item()
            nx, ny, nhyd, nbc = [np.fromfile(fp, np.int32, 1).item() for _ in range(4)]

            x = np.fromfile(fp, np.float64, nx)
            y = np.fromfile(fp, np.float64, ny)
            Q = np.fromfile(fp, np.float32, nx * ny * nhyd).reshape(ny, nx, nhyd)
            Bc = np.fromfile(fp, np.float32, nx * ny * nbc).reshape(ny, nx, nbc)

        q_names = ["rho", "vx", "vy", "vz", "pre"]
        b_names = ["Bx", "By", "Bz"]

        data_dict = {name: Q[:, :, i] for i, name in enumerate(q_names)}
        data_dict.update({name: Bc[:, :, i] for i, name in enumerate(b_names)})

    else:
        raise ValueError("filetype should be 'ascii' or 'binary'")

    data_dict["Bpre"] = 0.5*( data_dict["Bx"]**2 + data_dict["By"]**2 + data_dict["Bz"]**2)
    data_dict["Ekin"] = 0.5*( data_dict["vx"]**2 + data_dict["vy"]**2 + data_dict["vz"]**2)
    data_dict["beta"] = data_dict["pre"]/data_dict["Bpre"]
    return x, y, time, data_dict

def get_norm(results, varname, scale_type, vmin_manual=None, vmax_manual=None):
    vals = np.concatenate([r["data"][varname].ravel() for r in results])
    vals = vals[np.isfinite(vals)]

    if vals.size == 0:
        raise ValueError("no finite values found")

    if scale_type == "linear":
        vmin = np.nanmin(vals) if vmin_manual is None else vmin_manual
        vmax = np.nanmax(vals) if vmax_manual is None else vmax_manual

        if not np.isfinite(vmin) or not np.isfinite(vmax):
            raise ValueError("vmin/vmax must be finite")

        if vmin >= vmax:
            raise ValueError("vmin must be smaller than vmax")

        return Normalize(vmin=vmin, vmax=vmax)

    elif scale_type == "log":
        positive_vals = vals[vals > 0.0]

        if positive_vals.size == 0:
            raise ValueError(
                f"log scale cannot be used for '{varname}' because there are no positive values"
            )

        vmin = np.nanmin(positive_vals) if vmin_manual is None else vmin_manual
        vmax = np.nanmax(positive_vals) if vmax_manual is None else vmax_manual

        if not np.isfinite(vmin) or not np.isfinite(vmax):
            raise ValueError("vmin/vmax must be finite")

        if vmin <= 0.0 or vmax <= 0.0:
            raise ValueError("for log scale, vmin and vmax must be positive")

        if vmin >= vmax:
            raise ValueError("vmin must be smaller than vmax")

        return LogNorm(vmin=vmin, vmax=vmax)

    else:
        raise ValueError("scale_type should be 'linear' or 'log'")

def vecpot(x, y, Bx, By):
    Az = np.zeros_like(By)

    dy = y[1] - y[0]
    dx = x[1] - x[0]

    Az[1:, 0] = Az[0, 0] + np.cumsum(0.5 * (Bx[1:, 0] + Bx[:-1, 0]) * dy)
    Az[:, 1:] = Az[:, [0]] - np.cumsum(0.5 * (By[:, 1:] + By[:, :-1]) * dx, axis=1)

    return Az


if len(sys.argv) < 5:
    print("Usage:")
    print("  python3 MakePlot.py ascii  20 rho ct hdc")
    print("  python3 MakePlot.py binary 20 Bx  ct hdc run3")
    sys.exit(1)

if len(sys.argv) < 6:
    print("Usage:")
    print("  python3 MakeCompare.py ascii  20 rho linear ct hdc")
    print("  python3 MakeCompare.py binary 20 Bx  log    ct hdc run3")
    sys.exit(1)

parser = argparse.ArgumentParser()

parser.add_argument("filetype", choices=["ascii", "binary"])
parser.add_argument("step", type=int)
parser.add_argument("varname", type=str)
parser.add_argument("scale_type", choices=["linear", "log"])
parser.add_argument("dirnames", nargs="+")
parser.add_argument("--vmin", type=float, default=None)
parser.add_argument("--vmax", type=float, default=None)

args = parser.parse_args()

filetype = args.filetype
step = args.step
varname = args.varname
scale_type = args.scale_type
dirnames = args.dirnames
vmin_manual = args.vmin
vmax_manual = args.vmax

varlabel_dict = {
    "rho":  r"$\rho$",
    "vx":   r"$v_x$",
    "vy":   r"$v_y$",
    "vz":   r"$v_z$",
    "pre":  r"$P$",
    "Bx":   r"$B_x$",
    "By":   r"$B_y$",
    "Bz":   r"$B_z$",
    "Bpre": r"$B^2/2$",
    "Ekin": r"$v^2/2$",
    "beta": r"$\beta$",
}

varlabel = varlabel_dict.get(varname, varname)


results = []
times = []

for dirname in dirnames:
    x, y, time, data = read_data(filetype, dirname, step)

    if varname not in data:
        print(f"error: variable '{varname}' is not available in {dirname}")
        print("available variables:", ", ".join(data.keys()))
        sys.exit(1)

    Az = vecpot(x, y, data["Bx"], data["By"])
    results.append({
        "dirname": dirname,
        "x": x,
        "y": y,
        "time": time,
        "data": data,
        "Az": Az,
    })
    times.append(time)

#if not np.allclose(times, times[0]):
#    print("warning: times are not identical among directories")
for r in results: 
    print(f"  {r['dirname']}: time = {r['time']}")

time = results[0]["time"]

dx   = results[0]["x"][1] - results[0]["x"][0]
xmin = results[0]["x"][0] - 0.5*dx
xmax = results[0]["x"][-1] + 0.5*dx

dy   = results[0]["y"][1] - results[0]["y"][0]
ymin = results[0]["y"][0] - 0.5*dy
ymax = results[0]["y"][-1] + 0.5*dy

npanel = len(results)

fig_height = 4.8
panel_width = fig_height * (xmax - xmin) / (ymax - ymin)
fig_width = panel_width * npanel + 1.0

fig = plt.figure(figsize=(fig_width, fig_height))

grid = ImageGrid(
    fig, 111,
    nrows_ncols=(1, npanel),
    axes_pad=0.15,
    share_all=True,
    cbar_location="right",
    cbar_mode="single",
    cbar_size="4%",
    cbar_pad=0.10,
)

norm = get_norm(results, varname, scale_type, vmin_manual, vmax_manual)

im = None
for ax, r in zip(grid, results):
    x = r["x"]
    y = r["y"]
    Az = r["Az"]
    data = r["data"]

    im = ax.imshow( data[varname], extent=(xmin,xmax,ymin,ymax), origin="lower", norm=norm)
    ax.contour( x, y, Az, linestyles="solid", levels=20, colors="white", linewidths=0.7)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_title(r["dirname"])
    ax.set_xlabel("x axis")
    ax.minorticks_on()

grid[0].set_ylabel("y axis")

cbar = grid.cbar_axes[0].colorbar(im)
cbar.set_label(varlabel_dict.get(varname,varname))
cbar.minorticks_on()

fig.suptitle(rf"$\mathrm{{time}} = {time:.2f}$")

if len(dirnames) > 1:
    outroot = "compare"
    for dirname in dirnames:
        outroot += "_" + os.path.basename(os.path.normpath(dirname))
else:
    outroot = os.path.basename(os.path.normpath(dirnames[0]))

for format_fig in ["pdf", "png"]:
    outdir = format_fig + "file"
    os.makedirs(outdir, exist_ok=True)

    outputfile = os.path.join( outdir, varname + "_snap%05d_%s.%s" % (step, outroot, format_fig))
    print("making plot file", outputfile)
    plt.savefig(outputfile, bbox_inches="tight")

plt.show()
