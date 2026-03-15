# -*- coding: utf-8 -*-
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import ArtistAnimation
from matplotlib import gridspec
import re

dirname = sys.argv[1]
nmin = int(sys.argv[2])
nmax = int(sys.argv[3])

x = np.linspace(0, np.pi * 4, 100)

fig = plt.figure(figsize=(7,5)) 
gs = gridspec.GridSpec(1, 1, wspace=0.30,hspace=0.2) 
ax = []
for i in range(1):
    ax.append( plt.subplot(gs[i]) )

for i,ax0 in enumerate(ax):
    ax0.minorticks_on()
    ax0.set_xlim(-0.5,0.5)
    ax0.set_ylim(-0.15,0.15)
    ax0.grid()

ax[0].set_xlabel(r"$x$")

frames = []  # 各フレームを構成する Artist 一覧

# フレームごとの Artist を作成する。
icount = 0
for istep in range(nmin,nmax+1):
    foutname = dirname + "/snap%05d.dat"%(istep)
    print("reading " + foutname)
    with open(foutname, 'r') as data_file: 
        line = data_file.readline() 
        attributes = re.search(r'# (\S+)', line)
    time = float(attributes.group(1))

    data = np.loadtxt(foutname)

    x = data[:,0]
    By = data[:,7]
    Bz = data[:,8]

    # グラフを作成する。
    pg00, = ax[0].plot(x, By, 'o-',c="r",mfc="none",label="numerical (B_y)")
    pg01, = ax[0].plot(x, Bz, 'o-',c="b",mfc="none",label="numerical (B_z)")
    pg02, = ax[0].plot(x, 0.1*np.sin(2.0*np.pi*(x-time)), '-',c="k",label="exact sol.")
    pg03, = ax[0].plot(x, 0.1*np.cos(2.0*np.pi*(x-time)), '-',c="k")
    pg3   = ax[0].text(0.02,1.02,r"$\mathrm{time} = %.2f$"%(time),horizontalalignment="left",transform=ax[0].transAxes)

    if icount == 0: 
        ax[0].legend(ncol=3,loc="upper left", bbox_to_anchor=(0.02, 0.98), borderaxespad=0.0)

    # このフレームの Artist 一覧を追加する。
#    frames.append([pg00,pg01,pg10,pg11,pg20,pg21,pg3])
    frames.append([pg3,pg00,pg01,pg02,pg03])

    icount +=1

# アニメーションを作成する。
ani = ArtistAnimation(fig, frames, interval=200)

# mp4 画像として保存する。
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
#plt.close()
