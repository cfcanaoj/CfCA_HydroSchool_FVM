# -*- coding: utf-8 -*-
import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import ArtistAnimation
from matplotlib import gridspec
import re



#plt.rcParams['font.size'] = 20

dirname = sys.argv[1]
nmin = int(sys.argv[2])
nmax = int(sys.argv[3])

x = np.linspace(0, np.pi * 4, 100)

fig = plt.figure(figsize=(8,10)) 
gs = gridspec.GridSpec(3, 1, wspace=0.30,hspace=0.2) 
ax = []
for i in range(2):
    ax.append( plt.subplot(gs[i]) )

ylabel = [r"$E$",r"$F$"]
for i,ax0 in enumerate(ax):
    ax0.minorticks_on()
    ax0.set_xlim(0.0,1.0)
    ax0.set_ylabel(ylabel[i])

ax[1].set_xlabel(r"$x$")

frames = []  # 各フレームを構成する Artist 一覧


# フレームごとの Artist を作成する。
icount = 0
for istep in range(nmin,nmax+1):
    foutname = dirname + "/snap%05d.dat"%(istep)
    print("reading " + foutname)
    with open(foutname, 'r') as data_file: 
        line = data_file.readline() 
        parts = line.split()
    time = float(parts[2])

    data = np.genfromtxt(foutname,skip_header=2)

    x, E, F = np.split(data,3,1)

    # グラフを作成する。
    pg00, = ax[0].plot(x, E, 'o-',c="r",label="numerical")

    pg10, = ax[1].plot(x, F, 'o-',c="r",label="numerical")


    pg3 = ax[0].text(0,1.10,r"$\mathrm{time} = %.2f$"%(time),horizontalalignment="center")

    if icount == 0: 
        ax[0].legend()

    # このフレームの Artist 一覧を追加する。
    frames.append([pg00,pg10])

    icount +=1

# アニメーションを作成する。
ani = ArtistAnimation(fig, frames, interval=50)

# mp4 画像として保存する。
ani.save(dirname + "/animation.mp4", writer="imagemagick")
#plt.show()
plt.close()
