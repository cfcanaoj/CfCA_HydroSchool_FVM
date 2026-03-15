import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import ArtistAnimation
from matplotlib import gridspec
import re

if len(sys.argv) < 4:
    print("Usage: python MakeAnime.py nmin nmax dirname [save]")
    sys.exit(1)

nmin = int(sys.argv[1])
nmax = int(sys.argv[2])
dirname = sys.argv[3]
save_png = (len(sys.argv) >= 5 and sys.argv[4].lower() == "save")

fig = plt.figure(figsize=(10, 10))
gs = gridspec.GridSpec(2, 2, wspace=0.30, hspace=0.2)
ax = []
for i in range(4):
    ax.append(plt.subplot(gs[i]))

ylabel = [r"$\rho$", r"$v_x$", r"$v_y$", r"$B_y$"]
for i, ax0 in enumerate(ax):
    ax0.minorticks_on()
    ax0.set_xlim(-0.5, 0.5)
    ax0.set_ylabel(ylabel[i])

ax[2].set_xlabel(r"$x$")
ax[3].set_xlabel(r"$x$")

anasol = np.loadtxt("briowu_nonregsol.dat")
tout_ana = 0.1
x_ana = anasol[:, 0]

den_ana = np.concatenate(([anasol[0, 1]], anasol[:, 1], [anasol[-1, 1]]))
vx_ana  = np.concatenate(([anasol[0, 2]], anasol[:, 2], [anasol[-1, 2]]))
vy_ana  = np.concatenate(([anasol[0, 3]], anasol[:, 3], [anasol[-1, 3]]))
by_ana  = np.concatenate(([anasol[0, 7]], anasol[:, 7], [anasol[-1, 7]]))

if save_png:
    imgdir = os.path.join(dirname, "pngfile")
    os.makedirs(imgdir, exist_ok=True)

    for istep in range(nmin, nmax + 1):
        foutname = os.path.join(dirname, "snap%05d.dat" % istep)
        print("reading " + foutname)

        with open(foutname, "r") as data_file:
            line = data_file.readline()
            attributes = re.search(r"#\s*(\S+)", line)
        time = float(attributes.group(1))

        data = np.loadtxt(foutname)

        x   = data[:, 0]
        den = data[:, 1]
        vx  = data[:, 2]
        vy  = data[:, 3]
        By  = data[:, 7]

        x_ana1 = np.concatenate(([-0.5], x_ana * time / tout_ana, [0.5]))

        for a in ax:
            a.cla()
            a.minorticks_on()
            a.set_xlim(-0.5, 0.5)

        for i, a in enumerate(ax):
            a.set_ylabel(ylabel[i])
        ax[2].set_xlabel(r"$x$")
        ax[3].set_xlabel(r"$x$")

        ax[0].plot(x, den, "o-", c="tab:blue", mfc="none", label="numerical")
        ax[0].plot(x_ana1, den_ana, "-", c="b", label="exact")

        ax[1].plot(x, vx, "o-", c="tab:blue", mfc="none", label="numerical")
        ax[1].plot(x_ana1, vx_ana, "-", c="k", label="exact")

        ax[2].plot(x, vy, "o-", c="tab:blue", mfc="none", label="numerical")
        ax[2].plot(x_ana1, vy_ana, "-", c="k", label="exact")

        ax[3].plot(x, By, "o-", c="tab:blue", mfc="none", label="numerical")
        ax[3].plot(x_ana1, by_ana, "-", c="k", label="exact")

        ax[0].legend()
        ax[0].text(0.5, 1.02, r"$\mathrm{time} = %.2f$" % time, horizontalalignment="center",\
              transform=ax[0].transAxes)

        fname_png = os.path.join(imgdir, "snap%05d_"%(istep) + dirname + ".png")
        print("save snapshot",fname_png)
        fig.savefig(fname_png, dpi=100)

else:
    frames = []
    icount = 0

    for istep in range(nmin, nmax + 1):
        foutname = os.path.join(dirname, "snap%05d.dat" % istep)
        print("reading " + foutname)

        with open(foutname, "r") as data_file:
            line = data_file.readline()
            attributes = re.search(r"#\s*(\S+)", line)
        time = float(attributes.group(1))

        data = np.loadtxt(foutname)

        x   = data[:, 0]
        den = data[:, 1]
        vx  = data[:, 2]
        vy  = data[:, 3]
        By  = data[:, 7]

        x_ana1 = np.concatenate(([-0.5], x_ana * time / tout_ana, [0.5]))

        pg00, = ax[0].plot(x, den, "o-", c="tab:blue", mfc="none", label="numerical")
        pg01, = ax[0].plot(x_ana1, den_ana, "-", c="k", label="exact")

        pg10, = ax[1].plot(x, vx, "o-", c="tab:blue", mfc="none", label="numerical")
        pg11, = ax[1].plot(x_ana1, vx_ana, "-", c="k", label="exact")

        pg20, = ax[2].plot(x, vy, "o-", c="tab:blue", mfc="none", label="numerical")
        pg21, = ax[2].plot(x_ana1, vy_ana, "-", c="k", label="exact")

        pg30, = ax[3].plot(x, By, "o-", c="tab:blue", mfc="none", label="numerical")
        pg31, = ax[3].plot(x_ana1, by_ana, "-", c="k", label="exact")

        pg3 = ax[0].text(0.5, 1.02, r"$\mathrm{time} = %.2f$" % time, horizontalalignment="center",\
              transform=ax[0].transAxes)

        if icount == 0:
            ax[0].legend()

        frames.append([pg00, pg01, pg10, pg11, pg20, pg21, pg30, pg31, pg3])
        icount += 1

    ani = ArtistAnimation(fig, frames, interval=50)

    fname_anime = "animation_" + dirname + ".mp4"
    print("save animation file ",fname_anime)
    ani.save(
        fname_anime,
        writer="ffmpeg",
        dpi=150,
        codec="libx264",
        extra_args=["-pix_fmt", "yuv420p",
                    "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2"]
    )
    plt.show()

# plt.close()
