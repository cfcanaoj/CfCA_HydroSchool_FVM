import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
from dataclasses import dataclass
import re

fsizeforfig=16
fsizeforlabel=18
plt.rcParams.update({
    # font
    "font.family": "sans-serif",
    "font.size": fsizeforfig,
    "axes.labelsize": fsizeforlabel,
    # ticks
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,

    # line / axis
    "lines.linewidth": 4.5,
    "axes.linewidth": 2.0,
    "xtick.major.width": 1.8,
    "ytick.major.width": 1.8,
    "xtick.minor.width": 1.4,
    "ytick.minor.width": 1.4,
})

from cycler import cycler
# red, blue, green, sky-blue azure, pink, orange, purple, brown
cmapudc =  cycler(color=["#66ccff","#ff2800","#0041ff" ,"#35a16B", "#ff99a0","#ff9900" ,"#9a0079", "#663300"])
plt.rcParams['axes.prop_cycle'] = cmapudc
cmap = ["#66ccff", "#ff2800","#0041ff" ,"#35a16B", "#ff99a0","#ff9900" ,"#9a0079", "#663300"]

@dataclass
class FluidSnapshot:
    modelname: str
    time: float 
    nx: int
    x:  np.ndarray
    E:  np.ndarray
    Fx: np.ndarray

def ReadSnapshot(filename,modelname) -> FluidSnapshot:
    print("reading " + filename)
    with open(filename, 'r') as data_file: 
        line = data_file.readline() 
        parts = line.split()
        time = float(parts[2])
        line = data_file.readline() 
        parts = line.split()
        nx = float(parts[2])
        data = np.genfromtxt(filename,skip_header=2)
        x, E, Fx = np.split(data,3,1)
    return FluidSnapshot(
        modelname = modelname,
        time = time,
        nx = nx,
        x = x,
        E = E,
        Fx = Fx,
        )
#dirname = sys.argv[1]
step = int(sys.argv[1])

model_dirs ={
    "ana" : "analytic",
    "sn" : r"$S_n$",
    "m1" : "M1",
   "fld" : "FLD",
}

snapshots: list[FluidSnapshot] = []
for dirname, modelname in model_dirs.items():
    filename = dirname + "/snap%05d.dat"%(step)
    snap = ReadSnapshot(filename,modelname)
    snapshots.append(snap)

fig = plt.figure(figsize=(10,8),dpi=150) 
ax  = plt.subplot(111)
ax.set_xlabel(r"$x$")
ax.set_xlim(-1.0,1.0)
ax.set_ylabel(r"$E$")
ax.set_ylim(0.0,1.0)

enorm = 1.0e20
timenorm = 3.0e10
for snap in snapshots:
# 複数のプロットからなるグラフを作成する。
    ax.plot(snap.x, snap.E/enorm,'-',label=snap.modelname)
ax.text(0.9,1.04,r"$ct = %.2f$"%(snapshots[0].time*timenorm),horizontalalignment="center",fontsize=plt.rcParams["axes.labelsize"])
ax.legend(frameon=False)

filename="cmp" + '/Ecmp%05d.png'%(step)
print("file is saved as ",filename)
plt.tight_layout()
plt.savefig(filename)
