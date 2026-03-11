# -*- coding: utf-8 -*-
import os
import re
import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import ArtistAnimation
from matplotlib import gridspec

if len(sys.argv) < 4:
    print("Usage: python MakeAnime.py dirname nmin nmax [save]")
    sys.exit(1)

dirname = sys.argv[1]
nmin = int(sys.argv[2])
nmax = int(sys.argv[3])

save_png = (len(sys.argv) >= 5 and sys.argv[4].lower() == "save")

fig = plt.figure(figsize=(8, 10))
gs = gridspec.GridSpec(3, 1, wspace=0.30, hspace=0.2)
ax = [plt.subplot(gs[i]) for i in range(3)]

ylabel = [r"$\rho$", r"$v$", r"$P$"]
for i, ax0 in enumerate(ax):
    ax0.minorticks_on()
    ax0.set_xlim(-0.5, 0.5)
    ax0.set_ylabel(ylabel[i])

ax[2].set_xlabel(r"$x$")

anasol = np.loadtxt("sod_ana.dat")
tout_ana = 0.2
x_ana = anasol[:, 0]

den_ana = np.concatenate(([anasol[0, 1]], anasol[:, 1], [anasol[-1, 1]]))
vel_ana = np.concatenate(([anasol[0, 2]], anasol[:, 2], [anasol[-1, 2]]))
pre_ana = np.concatenate(([anasol[0, 3]], anasol[:, 3], [anasol[-1, 3]]))

# -------------------------
# save mode: PNG保存だけ
# -------------------------
if save_png:
    imgdir = os.path.join(dirname, "pngfile")
    os.makedirs(imgdir, exist_ok=True)

    for istep in range(nmin, nmax + 1):
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

        x_ana1 = np.concatenate(([-0.5], x_ana * time / tout_ana, [0.5]))

        for a in ax:
            a.cla()
            a.minorticks_on()
            a.set_xlim(-0.5, 0.5)

        for i, a in enumerate(ax):
            a.set_ylabel(ylabel[i])
        ax[2].set_xlabel(r"$x$")

        ax[0].plot(x, den, 'o-', c="r", mfc="none", label="numerical")
        ax[0].plot(x_ana1, den_ana, '-', c="b", label="exact")

        ax[1].plot(x, vel, 'o-', c="r", mfc="none", label="numerical")
        ax[1].plot(x_ana1, vel_ana, '-', c="b", label="exact")

        ax[2].plot(x, pre, 'o-', c="r", mfc="none", label="numerical")
        ax[2].plot(x_ana1, pre_ana, '-', c="b", label="exact")

        ax[0].legend()
        ax[0].text(
            0.5, 1.10,
            r"$\mathrm{time} = %.2f$" % time,
            horizontalalignment="center",
            transform=ax[0].transAxes
        )

        fname_png = os.path.join(imgdir, "snap%05d.png" % istep)
        print("save snapshot",fname_png)
        fig.savefig(fname_png, dpi=100)

# -------------------------
# normal mode: 動画生成だけ
# -------------------------
else:
    ax[0].plot([], [], 'o-', c="r", mfc="none", label="numerical")
    ax[0].plot([], [], '-', c="b", label="exact")
    ax[0].legend()

    frames = []

    for istep in range(nmin, nmax + 1):
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

        x_ana1 = np.concatenate(([-0.5], x_ana * time / tout_ana, [0.5]))

        pg00, = ax[0].plot(x, den, 'o-', c="r", mfc="none")
        pg01, = ax[0].plot(x_ana1, den_ana, '-', c="b")

        pg10, = ax[1].plot(x, vel, 'o-', c="r", mfc="none")
        pg11, = ax[1].plot(x_ana1, vel_ana, '-', c="b")

        pg20, = ax[2].plot(x, pre, 'o-', c="r", mfc="none")
        pg21, = ax[2].plot(x_ana1, pre_ana, '-', c="b")

        pg3 = ax[0].text(
            0.5, 1.10,
            r"$\mathrm{time} = %.2f$" % time,
            horizontalalignment="center",
            transform=ax[0].transAxes
        )

        frames.append([pg00, pg01, pg10, pg11, pg20, pg21, pg3])

    ani = ArtistAnimation(fig, frames, interval=50, blit=False)

    fname_anime = "animation.mp4"
    print("save animation file ",fname_anime)
    ani.save(
        os.path.join(dirname, fname_anime),
        writer="ffmpeg",
        dpi=150,
        codec="libx264",
        extra_args=["-pix_fmt", "yuv420p",
                    "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2"]
    )
    plt.show()

