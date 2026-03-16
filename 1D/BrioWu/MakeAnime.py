# -*- coding: utf-8 -*-
import os
import re
import argparse
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import ArtistAnimation
from matplotlib import gridspec

# =========================================================
# argument
# =========================================================
parser = argparse.ArgumentParser(
    description="Create an animation comparing snapshots from multiple directories.",
    usage="python3 MakeAnime.py [step_s] [step_e] [dir1] [dir2] ...",
    epilog=(
        "Example:\n"
        "  python3 MakeAnime.py 0 40 lax\n"
        "  python3 MakeAnime.py 0 40 lax hll\n"
    ),
    formatter_class=argparse.RawTextHelpFormatter
)
parser.add_argument("step_s", type=int, help="start step")
parser.add_argument("step_e", type=int, help="end step")
parser.add_argument("dirnames", nargs="+", help="one or more directories containing snapshot files")
args = parser.parse_args()

step_s   = args.step_s
step_e   = args.step_e
dirnames = args.dirnames

# =========================================================
# figure
# =========================================================
fig = plt.figure(figsize=(10, 10))
gs = gridspec.GridSpec(2, 2, wspace=0.30, hspace=0.2)
ax = [plt.subplot(gs[i]) for i in range(4)]

ylabel = [r"$\rho$", r"$v_x$", r"$v_y$", r"$P$"]
for i, ax0 in enumerate(ax):
    ax0.minorticks_on()
    ax0.set_xlim(-0.5, 0.5)
    ax0.set_ylabel(ylabel[i])

ax[2].set_xlabel(r"$x$")

# =========================================================
# exact solution
# =========================================================
anasol = np.loadtxt("briowu_nonregsol.dat")
tout_ana = 0.1
x_ana = anasol[:, 0]

den_ana = np.concatenate(([anasol[0, 1]], anasol[:, 1], [anasol[-1, 1]]))
vx_ana  = np.concatenate(([anasol[0, 2]], anasol[:, 2], [anasol[-1, 2]]))
vy_ana  = np.concatenate(([anasol[0, 3]], anasol[:, 3], [anasol[-1, 3]]))
By_ana  = np.concatenate(([anasol[0, 7]], anasol[:, 7], [anasol[-1, 7]]))

colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

# =========================================================
# utility
# =========================================================
def read_snapshot(dirname, istep):
    foutname = os.path.join(dirname, "snap%05d.dat" % istep)
    print("reading " + foutname)

    with open(foutname, "r") as data_file:
        line = data_file.readline()
        attributes = re.search(r'#\s*(\S+)', line)

    if attributes is None:
        raise ValueError(f"Could not read time from header: {foutname}")

    time = float(attributes.group(1))

    data = np.loadtxt(foutname)
    x   = data[:, 0]
    den = data[:, 1]
    vx  = data[:, 2]
    vy  = data[:, 3]
    By  = data[:, 7]

    return time, x, den, vx, vy, By


def make_outroot(dirnames):
    if len(dirnames) > 1:
        outroot = "compare"
        for dirname in dirnames:
            outroot += "_" + os.path.basename(os.path.normpath(dirname))
    else:
        outroot = os.path.basename(os.path.normpath(dirnames[0]))
    return outroot


# =========================================================
# animation
# =========================================================
frames = []

for istep in range(step_s, step_e + 1):
    artists = []
    time_ref = None

    for idir, dirname in enumerate(dirnames):
        time, x, den, vx, vy, By = read_snapshot(dirname, istep)

        if time_ref is None:
            time_ref = time
            x_ana1 = np.concatenate(([-0.5], x_ana * time_ref / tout_ana, [0.5]))

            pg01, = ax[0].plot(x_ana1, den_ana, '-', c="k",lw=1.5)
            pg11, = ax[1].plot(x_ana1, vx_ana, '-', c="k",lw=1.5)
            pg21, = ax[2].plot(x_ana1, vy_ana, '-', c="k",lw=1.5)
            pg31, = ax[3].plot(x_ana1, By_ana, '-', c="k",lw=1.5)
            artists.extend([pg01, pg11, pg21, pg31])

        color = colors[idir % len(colors)]

        pg00, = ax[0].plot(x, den, 'o-', c=color, mfc="none")
        pg10, = ax[1].plot(x, vx, 'o-', c=color, mfc="none")
        pg20, = ax[2].plot(x, vy, 'o-', c=color, mfc="none")
        pg30, = ax[3].plot(x, By, 'o-', c=color, mfc="none")
        artists.extend([pg00, pg10, pg20, pg30])

    if istep == step_s:
        handles = []
        labels = []

        h_exact, = ax[0].plot([], [], '-', c='k', lw=1.5)
        handles.append(h_exact)
        labels.append("exact")

        for idir, dirname in enumerate(dirnames):
            color = colors[idir % len(colors)]
            label = os.path.basename(os.path.normpath(dirname))
            h, = ax[0].plot([], [], 'o-', c=color, mfc="none")
            handles.append(h)
            labels.append(label)

        ax[0].legend(handles, labels, fontsize=9)

    pg3 = ax[0].text(
        0.5, 1.02,
        r"$\mathrm{time} = %.2f$" % time_ref,
        horizontalalignment="center",
        transform=ax[0].transAxes
    )
    artists.append(pg3)

    frames.append(artists)

ani = ArtistAnimation(fig, frames, interval=100, blit=False)

outroot = make_outroot(dirnames)
fname_anime = "animation_%s.mp4" % outroot
print("save animation file", fname_anime)

ani.save(
    fname_anime,
    writer="ffmpeg",
    dpi=150,
    codec="libx264",
    extra_args=[
        "-pix_fmt", "yuv420p",
        "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2"
    ]
)

plt.show()
