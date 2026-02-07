import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
import re

dirname = sys.argv[1]
step = int(sys.argv[2])

x = np.linspace(0, np.pi * 4, 100)

fig = plt.figure(figsize=(10,10)) 
gs = gridspec.GridSpec(2, 2, wspace=0.30,hspace=0.2) 
ax = []
for i in range(2):
    ax.append( plt.subplot(gs[i]) )

ylabel = [r"$E$",r"$F_x$"]
for i,ax0 in enumerate(ax):
    ax0.minorticks_on()
    ax0.set_xlim(0.1,1.0)
    ax0.set_ylabel(ylabel[i])

ax[1].set_xlabel(r"$x$")

# フレームごとの Artist を作成する。

icount = 0
foutname = dirname + "/snap%05d.dat"%(step)
print("reading " + foutname)
with open(foutname, 'r') as data_file: 
    line = data_file.readline() 
    parts = line.split()
    time = float(parts[2])

data = np.genfromtxt(foutname,skip_header=2)
x, E, Fx = np.split(data,3,1)

# グラフを作成する。
ax[0].plot(x, E, 'o-',c="r",label="numerical")
ax[1].plot(x, Fx, 'o-',c="r",label="numerical")

ax[0].text(0,1.10,r"$\mathrm{time} = %.2f$"%(time),horizontalalignment="center")

if icount == 0: 
    ax[0].legend()


icount +=1

#plt.savefig(dirname + '/snap%05d.png'%(step),bbox_inches="tight", pat_inches=0.0,dpi=1000)
plt.savefig(dirname + '/snap%05d.pdf'%(step),bbox_inches="tight")

#plt.show()
#plt.close()
