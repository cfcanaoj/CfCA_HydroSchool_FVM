import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec


def read_textdata(file):
    with open(file, "r") as f:
        time = float(f.readline().split()[-1])

    arr = np.loadtxt(file)

    x = arr[:, 0]
    y = arr[:, 1]

    names = ["rho", "vx", "vy", "vz", "pre", "Bx", "By", "Bz"]
    data_dict = {name: arr[:, i+1] for i, name in enumerate(names)}

    return x, y, time, data_dict


#------------------------------------------------------------
# Usage
#   python MakePlot.py step dir1 [dir2 dir3 ...]
# Example
#   python MakePlot.py 10 lax roe hlld
#------------------------------------------------------------
if len(sys.argv) < 3:
    print("Usage: python MakePlot.py step dir1 [dir2 dir3 ...]")
    sys.exit(1)

step = int(sys.argv[1])
dirnames = sys.argv[2:]

fig = plt.figure(figsize=(10, 10))
gs = gridspec.GridSpec(3, 2, wspace=0.30, hspace=0.2)
ax = []
for i in range(6):
    ax.append(plt.subplot(gs[i]))

ylabel = [r"$\rho$", r"$P$", r"$v_y$", r"$v_z$",r"$B_y$",r"$B_z$"]
for i, ax0 in enumerate(ax):
    ax0.minorticks_on()
    ax0.set_xlim(-0.5, 0.5)
    ax0.set_ylabel(ylabel[i])

ax[4].set_xlabel(r"$x$")
ax[5].set_xlabel(r"$x$")

#------------------------------------------------------------
# exact solution
#------------------------------------------------------------
anasol = np.loadtxt("daiwoodward1_exact.dat")
tout_ana = 0.2
x_ana = anasol[:, 0] -0.5

den_ana = np.concatenate(([anasol[0, 1]], anasol[:, 1], [anasol[-1, 1]]))
pre_ana = np.concatenate(([anasol[0, 2]], anasol[:, 2], [anasol[-1, 2]]))
vx_ana  = np.concatenate(([anasol[0, 3]], anasol[:, 3], [anasol[-1, 3]]))
vy_ana  = np.concatenate(([anasol[0, 4]], anasol[:, 4], [anasol[-1, 4]]))
vz_ana  = np.concatenate(([anasol[0, 5]], anasol[:, 5], [anasol[-1, 5]]))
by_ana  = np.concatenate(([anasol[0, 7]], anasol[:, 7], [anasol[-1, 7]]))
bz_ana  = np.concatenate(([anasol[0, 8]], anasol[:, 8], [anasol[-1, 8]]))

#------------------------------------------------------------
# plot exact solution only once
#------------------------------------------------------------
time_ref = None
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

for idir, dirname in enumerate(dirnames):
    foutname = os.path.join(dirname, "snap%05d.dat" % step)
    print("reading " + foutname)

    x, y, time, data = read_textdata(foutname)

    if time_ref is None:
        time_ref = time
        x_ana1 = np.concatenate(([-0.5], x_ana * time_ref / tout_ana, [0.5]))

        ax[0].plot(x_ana1, den_ana, '-', c="k", lw=1.5, label="exact",zorder=99)
        ax[1].plot(x_ana1, pre_ana,  '-', c="k", lw=1.5, label="exact",zorder=99)
        ax[2].plot(x_ana1, vy_ana,  '-', c="k", lw=1.5, label="exact",zorder=99)
        ax[3].plot(x_ana1, vz_ana,  '-', c="k", lw=1.5, label="exact",zorder=99)
        ax[4].plot(x_ana1, by_ana,  '-', c="k", lw=1.5, label="exact",zorder=99)
        ax[5].plot(x_ana1, bz_ana,  '-', c="k", lw=1.5, label="exact",zorder=99)

    color = colors[idir % len(colors)]
    label = os.path.basename(os.path.normpath(dirname))

    ax[0].plot(x, data['rho'], 'o-', c=color, mfc="none", label=label)
    ax[1].plot(x, data['pre'],  'o-', c=color, mfc="none", label=label)
    ax[2].plot(x, data['vy'],  'o-', c=color, mfc="none", label=label)
    ax[3].plot(x, data['vz'],  'o-', c=color, mfc="none", label=label)
    ax[4].plot(x, data['By'],  'o-', c=color, mfc="none", label=label)
    ax[5].plot(x, data['Bz'],  'o-', c=color, mfc="none", label=label)

ax[0].text(
    0.5, 1.02, r"$\mathrm{time} = %.2f$" % (time_ref),
    horizontalalignment="center",
    transform=ax[0].transAxes
)

ax[0].legend(fontsize=9)

#------------------------------------------------------------
# output
#------------------------------------------------------------
if len(dirnames) > 1:
    outroot = "compare"
    for dirname in dirnames:
        outroot += "_" + os.path.basename(os.path.normpath(dirname))
else:
    outroot = os.path.basename(os.path.normpath(dirnames[0]))

for format_fig in ["pdf", "png"]:
    outdir = format_fig + "file"
    os.makedirs(outdir, exist_ok=True)

    outputfile = os.path.join(outdir, "snap%05d_%s.%s" % (step, outroot, format_fig))
    print("making plot file", outputfile)
    plt.savefig(outputfile, bbox_inches="tight")

plt.show()
# plt.close()
