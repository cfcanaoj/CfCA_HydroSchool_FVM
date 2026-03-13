import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.colors import Normalize, LogNorm


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
            Q  = np.fromfile(fp, np.float32, nx * ny * nhyd).reshape(ny, nx, nhyd)
            Bc = np.fromfile(fp, np.float32, nx * ny * nbc ).reshape(ny, nx, nbc)

        q_names = ["rho", "vx", "vy", "vz", "pre"]
        b_names = ["Bx", "By", "Bz"]

        data_dict = {name: Q[:, :, i] for i, name in enumerate(q_names)}
        data_dict.update({name: Bc[:, :, i] for i, name in enumerate(b_names)})

    else:
        raise ValueError("filetype should be ascii or binary")

    data_dict["Bpre"] = 0.5*(data_dict["Bx"]**2 + data_dict["By"]**2 + data_dict["Bz"]**2)
    data_dict["Ekin"] = 0.5*data_dict['rho']*(data_dict["vx"]**2 + data_dict["vy"]**2 + data_dict["vz"]**2)
    data_dict["beta"] = data_dict["pre"]/data_dict["Bpre"]

    return x, y, time, data_dict


def get_norm_from_values(vals, scale_type, vmin_manual=None, vmax_manual=None):
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
            raise ValueError("log scale cannot be used because there are no positive values")

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
        raise ValueError("scale_type should be linear or log")


def collect_sampled_values(filetype, step_s, step_e, dirnames, varname):
    step_m = (step_s + step_e) // 2
    steps_sample = sorted(set([step_s, step_m, step_e]))

    dirname_ref = dirnames[0]
    print("sampling steps for color range:", steps_sample)
    print("reference directory for color range:", dirname_ref)

    vals_list = []

    for step in steps_sample:
        _, _, _, data = read_data(filetype, dirname_ref, step)

        if varname not in data:
            raise ValueError(
                f"variable '{varname}' is not available in {dirname_ref} at step {step}"
            )

        vals_list.append(data[varname].ravel())

    return np.concatenate(vals_list)

parser = argparse.ArgumentParser(
    description="Create a comparison movie from multiple simulation outputs.",
    usage="python3 MakeAnime.py [ascii|binary] [step_s] [step_e] [varname] [linear|log] [dir1] [dir2] ... [-h] [--vmin VMIN] [--vmax VMAX] [--interval INTERVAL]",
    epilog=(
        "Example:\n"
        "  python3 MakeAnime.py ascii 0 20 rho linear ct hdc\n"
        "  python3 MakeAnime.py binary 10 50 beta log ct hdc --vmin 1e-2 --vmax 1e2"
    ),
    formatter_class=argparse.RawTextHelpFormatter
)
parser.add_argument("filetype", choices=["ascii", "binary"],help="input file format")
parser.add_argument("step_s", type=int, help="starting step number")
parser.add_argument("step_e", type=int, help="ending step number")
parser.add_argument("varname", type=str, help="variable to plot (rho, vx, vy, vz, P, Bx, By, Bz, Bpre, Ekin, beta)")
parser.add_argument("scale_type", choices=["linear", "log"], help="color scale type")
parser.add_argument("dirnames", nargs="+", help="one or more directories containing snapshot files")
parser.add_argument("--vmin", type=float, default=None, help="(optional) manual minimum value for the color scale")
parser.add_argument("--vmax", type=float, default=None, help="(optional) manual maximum value for the color scale")
parser.add_argument("--interval", type=int, default=200, help="(optional) frame interval in milliseconds")
args = parser.parse_args()

filetype    = args.filetype
step_s      = args.step_s
step_e      = args.step_e
varname     = args.varname
scale_type  = args.scale_type
dirnames    = args.dirnames
vmin_manual = args.vmin
vmax_manual = args.vmax
interval    = args.interval

varlabel_dict = {
    "rho":  "density",
    "vx":   "velocity x",
    "vy":   "velocity y",
    "vz":   "velocity z",
    "pre":  "gas pressure",
    "Bx":   "magnetic field x",
    "By":   "magnetic field y",
    "Bz":   "magnetic field z",
    "Bpre": "magnetic pressure",
    "Ekin": "kinetic energy",
    "beta": "plasma beta",
}


if step_e < step_s:
    raise ValueError("step_e must be >= step_s")

if (vmin_manual is not None) and (vmax_manual is not None):
    norm = get_norm_from_values(
        np.array([vmin_manual, vmax_manual]),
        scale_type,
        vmin_manual,
        vmax_manual
    )
else:
    sample_vals = collect_sampled_values(filetype, step_s, step_e, dirnames, varname)
    norm = get_norm_from_values(sample_vals, scale_type, vmin_manual, vmax_manual)

# 最初のフレームで座標範囲を決める
x0, y0, time0, data0 = read_data(filetype, dirnames[0], step_s)

dx = x0[1] - x0[0]
dy = y0[1] - y0[0]
xmin = x0[0]  - 0.5 * dx
xmax = x0[-1] + 0.5 * dx
ymin = y0[0]  - 0.5 * dy
ymax = y0[-1] + 0.5 * dy

npanel = len(dirnames)

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

for i, ax in enumerate(grid):
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel("x axis")
    ax.set_title(dirnames[i])
    ax.minorticks_on()
grid[0].set_ylabel("y axis")

graph_list = []

for istep in range(step_s, step_e + 1):
    artists = []
    times = []

    for i, dirname in enumerate(dirnames):
        x, y, time, data = read_data(filetype, dirname, istep)
        times.append(time)

        im = grid[i].imshow(
            data[varname],
            extent=(xmin, xmax, ymin, ymax),
            origin="lower",
            norm=norm,
            animated=True,
        )
        artists.append(im)

    time_text = fig.text(
        0.5, 0.96,
        rf"$\mathrm{{time}}={times[0]:.2f}$",
        ha="center", va="center"
    )
    artists.append(time_text)

    if istep == step_s:
        cbar = grid.cbar_axes[0].colorbar(artists[0])
        cbar.set_label(varlabel_dict.get(varname,varname))
        cbar.minorticks_on()

    graph_list.append(artists)

ani = animation.ArtistAnimation(fig, graph_list, interval=interval)

if len(dirnames) > 1:
    joined = "_vs_".join(os.path.basename(os.path.normpath(d)) for d in dirnames)
else:
    joined = os.path.basename(os.path.normpath(dirnames[0]))

fname_anime = f"{varname}_snap{step_s:05d}_{step_e:05d}_{joined}.mp4"
outdir = "mp4file"
os.makedirs(outdir, exist_ok=True)
outfile = os.path.join(outdir, fname_anime)

print("making animation file", outfile)
ani.save(
    outfile,
    writer="ffmpeg",
    dpi=150,
    codec="libx264",
    extra_args=[
        "-pix_fmt", "yuv420p",
        "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2"
    ]
)

plt.show()
