import os
import sys
import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.colors import Normalize, LogNorm

SCRIPT_DIR = Path(__file__).resolve().parent


def resolve_dir(dirname):
    path = Path(dirname)
    if path.exists():
        return path

    alt = SCRIPT_DIR / dirname
    if alt.exists():
        return alt

    raise FileNotFoundError(path)


def _assemble_blocks(blocks, step):
    x_starts = sorted({block[1][0] for block in blocks})
    y_starts = sorted({block[2][0] for block in blocks})

    x_offsets = {}
    y_offsets = {}
    x_parts = []
    y_parts = []

    offset = 0
    for x0 in x_starts:
        x = next(block[1] for block in blocks if block[1][0] == x0)
        x_offsets[x0] = offset
        x_parts.append(x)
        offset += x.size

    offset = 0
    for y0 in y_starts:
        y = next(block[2] for block in blocks if block[2][0] == y0)
        y_offsets[y0] = offset
        y_parts.append(y)
        offset += y.size

    x = np.concatenate(x_parts)
    y = np.concatenate(y_parts)
    fields = np.empty((y.size, x.size, blocks[0][3].shape[2]), dtype=blocks[0][3].dtype)

    time = blocks[0][0]
    for block_time, block_x, block_y, block_fields in blocks:
        if not np.isclose(block_time, time):
            raise ValueError(f"inconsistent time at step {step}")
        ix = x_offsets[block_x[0]]
        iy = y_offsets[block_y[0]]
        fields[iy:iy + block_fields.shape[0], ix:ix + block_fields.shape[1], :] = block_fields

    return x, y, time, fields


def read_data(filetype, dirname, step):
    directory = resolve_dir(dirname)

    if filetype == "ascii":
        paths = sorted(directory.glob(f"snap*-{step:05d}.dat"))
        if not paths:
            single = directory / f"snap{step:05d}.dat"
            if not single.exists():
                raise FileNotFoundError(single)
            paths = [single]

        print("reading", str(directory / f"snap*-{step:05d}.dat"))
        blocks = []
        for path in paths:
            with path.open("r", encoding="ascii") as fp:
                time = float(fp.readline().split()[-1])
                nx, ny = map(int, fp.readline().split()[-2:])
                body = fp.read()
            arr = np.fromstring(body, sep=" ")
            ncols = 11
            expected = nx * ny * ncols
            if arr.size != expected:
                raise ValueError(
                    f"unexpected number of values in {path}: expected {expected}, got {arr.size}"
                )
            arr = arr.reshape(ny, nx, ncols)
            blocks.append((time, arr[0, :, 0], arr[:, 0, 1], arr[:, :, 2:]))

        x, y, time, fields = _assemble_blocks(blocks, step)
        names = ["rho", "vx", "vy", "vz", "pre", "Bx", "By", "Bz", "psi"]
        data_dict = {name: fields[:, :, i] for i, name in enumerate(names)}

    elif filetype == "binary":
        paths = sorted(directory.glob(f"snap*-{step:05d}.bin"))
        if not paths:
            single = directory / f"snap{step:05d}.bin"
            if not single.exists():
                raise FileNotFoundError(single)
            paths = [single]

        print("reading", str(directory / f"snap*-{step:05d}.bin"))
        blocks = []
        for path in paths:
            with path.open("rb") as fp:
                time = np.fromfile(fp, np.float64, 1).item()
                nx, ny, nhyd, nbc = [np.fromfile(fp, np.int32, 1).item() for _ in range(4)]
                x = np.fromfile(fp, np.float64, nx)
                y = np.fromfile(fp, np.float64, ny)
                q = np.fromfile(fp, np.float32, nx * ny * nhyd).reshape(ny, nx, nhyd)
                b = np.fromfile(fp, np.float32, nx * ny * nbc).reshape(ny, nx, nbc)
            fields = np.concatenate((q, b), axis=2)
            blocks.append((time, x, y, fields))

        x, y, time, fields = _assemble_blocks(blocks, step)
        names = ["rho", "vx", "vy", "vz", "pre", "Bx", "By", "Bz", "psi"]
        data_dict = {name: fields[:, :, i] for i, name in enumerate(names)}

    else:
        raise ValueError("filetype should be ascii or binary")

    data_dict["Bpre"] = 0.5 * (data_dict["Bx"] ** 2 + data_dict["By"] ** 2 + data_dict["Bz"] ** 2)
    data_dict["Ekin"] = 0.5 * data_dict["rho"] * (
        data_dict["vx"] ** 2 + data_dict["vy"] ** 2 + data_dict["vz"] ** 2
    )
    data_dict["beta"] = data_dict["pre"] / data_dict["Bpre"]

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

    if scale_type == "log":
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
            raise ValueError(f"variable '{varname}' is not available in {dirname_ref} at step {step}")
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
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument("filetype", choices=["ascii", "binary"], help="input file format")
parser.add_argument("step_s", type=int, help="starting step number")
parser.add_argument("step_e", type=int, help="ending step number")
parser.add_argument("varname", type=str, help="variable to plot (rho, vx, vy, vz, pre, Bx, By, Bz, psi, Bpre, Ekin, beta)")
parser.add_argument("scale_type", choices=["linear", "log"], help="color scale type")
parser.add_argument("dirnames", nargs="+", help="one or more directories containing snapshot files")
parser.add_argument("--vmin", type=float, default=None, help="manual minimum value for the color scale")
parser.add_argument("--vmax", type=float, default=None, help="manual maximum value for the color scale")
parser.add_argument("--interval", type=int, default=200, help="frame interval in milliseconds")
args = parser.parse_args()

filetype = args.filetype
step_s = args.step_s
step_e = args.step_e
varname = args.varname
scale_type = args.scale_type
dirnames = args.dirnames
vmin_manual = args.vmin
vmax_manual = args.vmax
interval = args.interval

varlabel_dict = {
    "rho": "density",
    "vx": "velocity x",
    "vy": "velocity y",
    "vz": "velocity z",
    "pre": "gas pressure",
    "Bx": "magnetic field x",
    "By": "magnetic field y",
    "Bz": "magnetic field z",
    "psi": "divergence cleaning field",
    "Bpre": "magnetic pressure",
    "Ekin": "kinetic energy",
    "beta": "plasma beta",
}

if step_e < step_s:
    raise ValueError("step_e must be >= step_s")

if (vmin_manual is not None) and (vmax_manual is not None):
    norm = get_norm_from_values(np.array([vmin_manual, vmax_manual]), scale_type, vmin_manual, vmax_manual)
else:
    sample_vals = collect_sampled_values(filetype, step_s, step_e, dirnames, varname)
    norm = get_norm_from_values(sample_vals, scale_type, vmin_manual, vmax_manual)

x0, y0, time0, data0 = read_data(filetype, dirnames[0], step_s)
dx = x0[1] - x0[0]
dy = y0[1] - y0[0]
xmin = x0[0] - 0.5 * dx
xmax = x0[-1] + 0.5 * dx
ymin = y0[0] - 0.5 * dy
ymax = y0[-1] + 0.5 * dy

npanel = len(dirnames)
fig_height = 4.8
panel_width = fig_height * (xmax - xmin) / (ymax - ymin)
fig_width = panel_width * npanel + 1.0
fig = plt.figure(figsize=(fig_width, fig_height))

grid = ImageGrid(
    fig,
    111,
    nrows_ncols=(1, npanel),
    axes_pad=0.40,
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

    time_text = fig.text(0.5, 0.96, rf"$\mathrm{{time}}={times[0]:.2f}$", ha="center", va="center")
    artists.append(time_text)

    if istep == step_s:
        cbar = grid.cbar_axes[0].colorbar(artists[0])
        cbar.set_label(varlabel_dict.get(varname, varname))
        cbar.minorticks_on()

    graph_list.append(artists)

ani = animation.ArtistAnimation(fig, graph_list, interval=interval)

if len(dirnames) > 1:
    joined = "_".join(os.path.basename(os.path.normpath(d)) for d in dirnames)
else:
    joined = os.path.basename(os.path.normpath(dirnames[0]))

outdir = "mp4file"
os.makedirs(outdir, exist_ok=True)
if animation.writers.is_available("ffmpeg"):
    fname_anime = f"{varname}_snap{step_s:05d}_{step_e:05d}_{joined}.mp4"
    outfile = os.path.join(outdir, fname_anime)
    writer = animation.FFMpegWriter(
        codec="libx264",
        extra_args=["-pix_fmt", "yuv420p", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2"],
    )
else:
    fname_anime = f"{varname}_snap{step_s:05d}_{step_e:05d}_{joined}.gif"
    outfile = os.path.join(outdir, fname_anime)
    writer = animation.PillowWriter()

print("making animation file", outfile)
ani.save(outfile, writer=writer, dpi=150)

plt.show()
