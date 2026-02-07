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

fig = plt.figure(figsize=(10,8)) 
gs = gridspec.GridSpec(1, 1, wspace=0.30,hspace=0.2) 
ax = []
for i in range(1):
    ax.append( plt.subplot(gs[i]) )

ylabel = [r"$\rho$",r"$v$",r"$P$"]
for i,ax0 in enumerate(ax):
    ax0.minorticks_on()
    ax0.set_xlim(-0.5,0.5)
    ax0.set_ylabel(ylabel[i])

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
    pg00, = ax[0].plot(x, By, 'o-',c="r",label="numerical (B_y)")
    pg01, = ax[0].plot(x, Bz, 'o-',c="r",label="numerical (B_z)")
    pg02, = ax[0].plot(x, -0.1*np.sin(2.0*np.pi*(x-time)), '-',c="k",label="exact sol. (B_y)")
    pg03, = ax[0].plot(x, -0.1*np.cos(2.0*np.pi*(x-time)), '-',c="k",label="exact sol. (B_y)")
#    pg01, = ax[0].plot(x_ana*time/tout_ana, den_ana, '-',c="b",label="exact")
#    pg10, = ax[1].plot(x, vel, 'o-',c="r",label="numerical")
#    pg11, = ax[1].plot(x_ana*time/tout_ana, vel_ana, '-',c="b",label="exact")
#    pg20, = ax[2].plot(x, pre, 'o-',c="r",label="numerical")
#    pg21, = ax[2].plot(x_ana*time/tout_ana, pre_ana, '-',c="b",label="exact")

#    pg3 = ax[0].text(0,1.10,r"$\mathrm{time} = %.2f$"%(time),horizontalalignment="center")

    if icount == 0: 
        ax[0].legend()

    # このフレームの Artist 一覧を追加する。
#    frames.append([pg00,pg01,pg10,pg11,pg20,pg21,pg3])
    frames.append([pg00,pg01,pg02,pg03])

    icount +=1

# アニメーションを作成する。
ani = ArtistAnimation(fig, frames, interval=50)

# mp4 画像として保存する。
ani.save(dirname + "/animation.mp4", writer="imagemagick")
#plt.show()
plt.close()
