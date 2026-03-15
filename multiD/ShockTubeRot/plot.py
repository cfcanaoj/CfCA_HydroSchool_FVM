import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import gridspec
import re

plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 15

fig = plt.figure(figsize=(10, 16)) 

ncol = 2
nrow = 3

gs = gridspec.GridSpec(nrow, ncol, wspace=0.2,hspace=0.2) 

#ana = np.loadtxt("dw_ana.dat")

ax = []
for i in range(ncol*nrow): 
    ax.append(plt.subplot(gs[i]))


istep = 40
foutname = "snap/snap%05d_ct.dat"%(istep)

data = np.loadtxt(foutname)
nx = 512
ny = 2

x = data[:,0]
y = data[:,1]
den = data[:,2]
vx = data[:,3]
vy = data[:,4]
pre = data[:,6]
bx = data[:,7]
by = data[:,8]


xi = (x + y*2.0)/np.sqrt(5.0) - 0.5*( 2.0/np.sqrt(5.0) + 2.0*ny/nx*2.0/np.sqrt(5.0))
phys = [den, pre,(vx + 2*vy)/np.sqrt(5.0) ,(-2.0*vx + vy)/np.sqrt(5.0),(bx + 2*by)/np.sqrt(5.0),(-2.0*bx + by)/np.sqrt(5.0)]
#anasol = [ana[:,1], ana[:,2], ana[:,3], ana[:,4], ana[:,6], ana[:,7], ana[:,8]]
ylabel = [r"$\rho$",r"$P$", r"$v_\xi$", r"$v_\eta$", r"$B_\xi$", r"$B_\eta$"]

for i, ax0 in enumerate(ax): 
    ax0.set_xlim(-0.5, 0.5)     
    ax0.set_xlabel(r"$\xi$")
    ax0.set_ylabel(ylabel[i])
    ax0.plot( xi, phys[i], 'o', markerfacecolor="none"  )
#    ax0.plot( ana[:,0]-0.5, anasol[i] )


plt.show()
