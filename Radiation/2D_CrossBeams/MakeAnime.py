import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from dataclasses import dataclass

# =====================
# 引数
# =====================
if len(sys.argv) != 4:
    print("Usage: python MakePlot.py dir nstart nend")
    sys.exit(1)

dirname = sys.argv[1]
nstart  = int(sys.argv[2])
nend    = int(sys.argv[3])

# =====================
# matplotlib 設定（そのまま）
# =====================
fsizeforfig=16
fsizeforlabel=18
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": fsizeforfig,
    "axes.labelsize": fsizeforlabel,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "lines.linewidth": 4.5,
    "axes.linewidth": 2.0,
    "xtick.major.width": 1.8,
    "ytick.major.width": 1.8,
    "xtick.minor.width": 1.4,
    "ytick.minor.width": 1.4,
})

# =====================
# データ構造
# =====================
@dataclass
class FluidSnapshot:
    modelname: str
    time: float
    nx: int
    ny: int
    x: np.ndarray
    y: np.ndarray
    E: np.ndarray
    Fx: np.ndarray
    Fy: np.ndarray

def ReadSnapshot(filename, modelname) -> FluidSnapshot:
    print("reading", filename)
    with open(filename, 'r') as f:
        time = float(f.readline().split()[2])
        nx   = int(f.readline().split()[2])
        ny   = int(f.readline().split()[2])

    data = np.genfromtxt(
        filename,
        comments="#",
        skip_header=4,
        names=("x","y","E","Fx","Fy"),
        dtype=float
    )

    return FluidSnapshot(
        modelname=modelname,
        time=time,
        nx=nx,
        ny=ny,
        x = data["x"].reshape(ny,nx),
        y = data["y"].reshape(ny,nx),
        E = data["E"].reshape(ny,nx),
        Fx= data["Fx"].reshape(ny,nx),
        Fy= data["Fy"].reshape(ny,nx),
    )

# =====================
# 初期フレーム
# =====================
snap0 = ReadSnapshot(f"{dirname}/snap{nstart:05d}.dat", dirname)

fig, ax = plt.subplots(figsize=(10, 8), dpi=150)
plt.subplots_adjust(right=0.85)
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")
ax.set_aspect(1)

im = ax.pcolormesh(
    snap0.x, snap0.y, snap0.E,
    shading="auto",
    cmap="afmhot",
    vmin=0.0,
    vmax=2.0
)

from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="4%", pad=0.1)
cbar = fig.colorbar(im, cax=cax, orientation="vertical")

timenorm = 1.0e0
time_text = ax.text(
    0.9, 1.05, "",
    transform=ax.transAxes,
    ha="center",
    fontsize=plt.rcParams["axes.labelsize"]
)

# =====================
# 更新関数
# =====================
def update(step):
    snap = ReadSnapshot(f"{dirname}/snap{step:05d}.dat", dirname)
    im.set_array(snap.E.ravel())
    time_text.set_text(rf"$ct = {snap.time*timenorm:.2f}$")
    return im, time_text

# =====================
# アニメーション作成
# =====================
frames = range(nstart, nend+1)

ani = animation.FuncAnimation(
    fig,
    update,
    frames=frames,
    blit=False
)

# =====================
# mp4 保存
# =====================
outfile = f"{dirname}/animate.mp4"
writer = animation.FFMpegWriter(fps=10, bitrate=6000)

print("writing movie:", outfile)
ani.save(outfile, writer=writer)

print("done.")
