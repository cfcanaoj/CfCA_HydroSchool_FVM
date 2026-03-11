import os
import re
import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import ArtistAnimation
from matplotlib import gridspec

gam = 5.0/3.0
amp = 1e-5

if len(sys.argv) < 4:
    print("Usage: python MakeAnime.py dirname nmin nmax [save]")
    sys.exit(1)

dirname = sys.argv[1]
nmin = int(sys.argv[2])
nmax = int(sys.argv[3])
save_png = (len(sys.argv) >= 5 and sys.argv[4].lower() == "save")

fig = plt.figure(figsize=(8, 10))
gs = gridspec.GridSpec(3, 1, wspace=0.30, hspace=0.2)
ax = []
for i in range(3):
    ax.append(plt.subplot(gs[i]))

ylabel = [r"$\rho$", r"$v$", r"$P$"]
for i, ax0 in enumerate(ax):
    ax0.minorticks_on()
    ax0.set_ylim(-1.1 * amp, 1.1 * amp)
    ax0.set_xlim(-0.5, 0.5)
    ax0.set_ylabel(ylabel[i])

ax[2].set_xlabel(r"$x$")


def read_snapshot(istep):
    foutname = os.path.join(dirname, "snap%05d.dat" % istep)
    print("reading " + foutname)
    with open(foutname, 'r') as data_file:
        line = data_file.readline()
        attributes = re.search(r'#\s*(\S+)', line)
    time = float(attributes.group(1))

    data = np.loadtxt(foutname)

    x = data[:, 0]
    den = data[:, 1]
    vel = data[:, 2]
    pre = data[:, 3]
    return time, x, den, vel, pre


if save_png:
    imgdir = os.path.join(dirname, "pngfile")
    os.makedirs(imgdir, exist_ok=True)

    for istep in range(nmin, nmax + 1):
        time, x, den, vel, pre = read_snapshot(istep)

        for a in ax:
            a.cla()
            a.minorticks_on()
            a.set_ylim(-1.1 * amp, 1.1 * amp)
            a.set_xlim(-0.5, 0.5)

        for i, a in enumerate(ax):
            a.set_ylabel(ylabel[i])
        ax[2].set_xlabel(r"$x$")

        ax[0].plot(x, den - 1.0, 'o-', c="r", mfc = "none", label="numerical")
        ax[0].plot(x, 1e-5 * np.sin(2.0 * np.pi * (x - time)), '-', c="b", label="exact")
        ax[1].plot(x, vel, 'o-', c="r", mfc = "none", label="numerical")
        ax[1].plot(x, 1e-5 * np.sin(2.0 * np.pi * (x - time)), '-', c="b", label="exact")
        ax[2].plot(x, pre - 1.0 / gam, 'o-', c="r", mfc = "none", label="numerical")
        ax[2].plot(x, 1e-5 * np.sin(2.0 * np.pi * (x - time)), '-', c="b", label="exact")

        ax[0].legend()
        ax[0].text(0.0, amp * 1.15, r"$\mathrm{time} = %.2f$" % time, horizontalalignment="center")

        fname_png = os.path.join(imgdir, "snap%05d.png" % istep)
        print("save snapshot",fname_png)
        fig.savefig(fname_png, dpi=100)
else:
    frames = []
    icount = 0

    for istep in range(nmin, nmax + 1):
        time, x, den, vel, pre = read_snapshot(istep)

        pg00, = ax[0].plot(x, den - 1.0, 'o-', c="r", mfc = "none",  label="numerical")
        pg01, = ax[0].plot(x, 1e-5 * np.sin(2.0 * np.pi * (x - time)), '-', c="b", label="exact")
        pg10, = ax[1].plot(x, vel, 'o-', c="r", mfc = "none", label="numerical")
        pg11, = ax[1].plot(x, 1e-5 * np.sin(2.0 * np.pi * (x - time)), '-', c="b", label="exact")
        pg20, = ax[2].plot(x, pre - 1.0 / gam, 'o-', c="r", mfc = "none",  label="numerical")
        pg21, = ax[2].plot(x, 1e-5 * np.sin(2.0 * np.pi * (x - time)), '-', c="b", label="exact")

        pg3 = ax[0].text(0, amp * 1.15, r"$\mathrm{time} = %.2f$" % time, horizontalalignment="center")

        if icount == 0:
            ax[0].legend()

        frames.append([pg00, pg01, pg10, pg11, pg20, pg21, pg3])
        icount += 1

    ani = ArtistAnimation(fig, frames, interval=50)
    ani.save(
        os.path.join(dirname, "animation.mp4"),
        writer="ffmpeg",
        dpi=150,
        codec="libx264",
        extra_args=["-pix_fmt", "yuv420p", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2"],
    )
    plt.show()

