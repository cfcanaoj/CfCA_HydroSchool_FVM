import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import ArtistAnimation
from matplotlib import gridspec
import re

gam = 5.0/3.0
amp = 1e-5

dirname = sys.argv[1]
nmin = int(sys.argv[2])
nmax = int(sys.argv[3])

fig = plt.figure(figsize=(8,10)) 
gs = gridspec.GridSpec(3, 1, wspace=0.30,hspace=0.2) 
ax = []
for i in range(3):
    ax.append( plt.subplot(gs[i]) )

ylabel = [r"$\rho$",r"$v$",r"$P$"]
for i,ax0 in enumerate(ax):
    ax0.minorticks_on()
    ax0.set_ylim(-1.1*amp,1.1*amp)
    ax0.set_xlim(-0.5,0.5)
    ax0.set_ylabel(ylabel[i])

ax[2].set_xlabel(r"$x$")

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
    den = data[:,1]
    vel = data[:,2]
    pre = data[:,3]

    # グラフを作成する。
    pg00, = ax[0].plot(x, den-1.0, 'o-',c="r",label="numerical")
    pg01, = ax[0].plot(x, 1e-5*np.sin(2.0*np.pi*(x-time)), '-',c="b",label="exact")
    pg10, = ax[1].plot(x, vel, 'o-',c="r",label="numerical")
    pg11, = ax[1].plot(x, 1e-5*np.sin(2.0*np.pi*(x-time)), '-',c="b",label="exact")
    pg20, = ax[2].plot(x, pre-1.0/gam, 'o-',c="r",label="numerical")
    pg21, = ax[2].plot(x, 1e-5*np.sin(2.0*np.pi*(x-time)), '-',c="b",label="exact")

    pg3 = ax[0].text(0,amp*1.15,r"$\mathrm{time} = %.2f$"%(time),horizontalalignment="center")

    if icount == 0: 
        ax[0].legend()

    # このフレームの Artist 一覧を追加する。
    frames.append([pg00,pg01,pg10,pg11,pg20,pg21,pg3])

    icount +=1

# アニメーションを作成する。
ani = ArtistAnimation(fig, frames, interval=50)

# mp4 画像として保存する。
ani.save(dirname + "/animation.mp4", writer="imagemagick")
plt.show()
#plt.close()
