# -*- coding: utf-8 -*-
import os
import re
import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec

if len(sys.argv) < 3:
    print("Usage: python MakePlot.py step dir1 [dir2 dir3 ...]")
    sys.exit(1)

step = int(sys.argv[1])
dirnames = sys.argv[2:]

fig = plt.figure(figsize=(8, 10))
gs = gridspec.GridSpec(3, 1, wspace=0.30, hspace=0.2)
ax = [plt.subplot(gs[i]) for i in range(3)]

ylabel = [r"$\rho$", r"$v$", r"$P$"]
for i, ax0 in enumerate(ax):
    ax0.minorticks_on()
    ax0.set_xlim(-0.5, 0.5)
    ax0.set_ylabel(ylabel[i])

ax[2].set_xlabel(r"$x$")

#------------------------------------------------------------
# exact solution
#------------------------------------------------------------
anasol = np.loadtxt("sod_ana.dat")
tout_ana = 0.2
x_ana = anasol[:, 0]

den_ana = np.concatenate(([anasol[0, 1]], anasol[:, 1], [anasol[-1, 1]]))
vel_ana = np.concatenate(([anasol[0, 2]], anasol[:, 2], [anasol[-1, 2]]))
pre_ana = np.concatenate(([anasol[0, 3]], anasol[:, 3], [anasol[-1, 3]]))

colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
time_ref = None

for idir, dirname in enumerate(dirnames):
    foutname = os.path.join(dirname, "snap%05d.dat" % step)
    print("reading " + foutname)

    with open(foutname, 'r') as data_file:
        line = data_file.readline()
        attributes = re.search(r'#\s*(\S+)', line)
    time = float(attributes.group(1))

    data = np.loadtxt(foutname)
    x   = data[:, 0]
    den = data[:, 1]
    vel = data[:, 2]
    pre = data[:, 3]

    if time_ref is None:
        time_ref = time
        x_ana1 = np.concatenate(([-0.5], x_ana * time_ref / tout_ana, [0.5]))

        ax[0].plot(x_ana1, den_ana, '-', c="k", lw=1.5, label="exact",zorder=99)
        ax[1].plot(x_ana1, vel_ana, '-', c="k", lw=1.5, label="exact",zorder=99)
        ax[2].plot(x_ana1, pre_ana, '-', c="k", lw=1.5, label="exact",zorder=99)

    color = colors[idir % len(colors)]
    label = os.path.basename(os.path.normpath(dirname))

    ax[0].plot(x, den, 'o-', c=color, mfc="none", label=label)
    ax[1].plot(x, vel, 'o-', c=color, mfc="none", label=label)
    ax[2].plot(x, pre, 'o-', c=color, mfc="none", label=label)

ax[0].legend(fontsize=9)
ax[0].text( 0.5, 1.02, r"$\mathrm{time} = %.2f$" % time_ref, horizontalalignment="center",\
    transform=ax[0].transAxes)

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
