import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from dataclasses import dataclass
import re

dirname = sys.argv[1]
step = int(sys.argv[2])

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
# red, blue, green, yellow, sky-blue azure, pink, orange, purple, brown
cmapudc =  cycler(color=["#ff2800","#0041ff" ,"#35a16B","#faf500","#66ccff", "#ff99a0","#ff9900" ,"#9a0079", "#663300"])
plt.rcParams['axes.prop_cycle'] = cmapudc
cmap = ["#ff2800","#0041ff" ,"#35a16B","#faf500","#66ccff", "#ff99a0","#ff9900" ,"#9a0079", "#663300"]

@dataclass
class FluidSnapshot:
    modelname: str
    time: float 
    nx: int
    ny: int
    x:  np.ndarray
    y:  np.ndarray
    E:  np.ndarray
    Fx: np.ndarray
    Fy: np.ndarray

def ReadSnapshot(filename,modelname) -> FluidSnapshot:
    print("reading " + filename)
    with open(filename, 'r') as data_file: 
        line = data_file.readline() # 1st line 
        parts = line.split()
        time = float(parts[2])
        print("time",time)
        line = data_file.readline() # 2nd line
        parts = line.split()
        nx = int(parts[2])
        line = data_file.readline() # 3rd line
        parts = line.split()
        ny = int(parts[2])
        print("nx,ny",nx,ny)
        data = np.genfromtxt(filename,comments="#",dtype=float,names=("x","y","E","Fx","Fy"))
        x  = data["x"].reshape(ny,nx)
        y  = data["y"].reshape(ny,nx)
        E  = data["E"].reshape(ny,nx)
        Fx = data["Fx"].reshape(ny,nx)
        Fy = data["Fy"].reshape(ny,nx)
        
        
    return FluidSnapshot(
        modelname = modelname,
        time = time,
        nx = nx,
        ny = ny,
        x = x,
        y = y,
        E = E,
        Fx = Fx,
        Fy = Fy,
        )

def makedirs(path):
    if not os.path.isdir(path):
        os.makedirs(path)

filename = dirname + "/snap%05d.dat"%(step)
snap = ReadSnapshot(filename,dirname)

fig, ax = plt.subplots(figsize=(10, 8), dpi=150)
plt.subplots_adjust(right=0.85)
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")

ax.set_aspect(1)
im = ax.pcolormesh(snap.x, snap.y, snap.E,shading="auto",cmap="afmhot",vmin=0.0,vmax=2.0)

from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="4%", pad=0.1)

cbar = fig.colorbar(im, cax=cax, orientation="vertical")

timenorm=1.0e0
ax.text(0.9,1.05,r"$ct = %.2f$"%(snap.time*timenorm),transform=ax.transAxes,horizontalalignment="center",fontsize=plt.rcParams["axes.labelsize"])

#ax.legend(frameon=False)

filename=dirname + '/E%05d.png'%(step)
print("file is saved as ",filename)
plt.savefig(filename)
