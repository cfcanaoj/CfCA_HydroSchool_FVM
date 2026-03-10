import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import re

def read_textdata(file):
    with open(file, "r") as f:
        time = float(f.readline().split()[-1])
        nx, ny = map(int, f.readline().split()[-2:])

    arr = np.loadtxt(file).reshape(ny, nx, -1)

    x = arr[:, :, 0]
    y = arr[:, :, 1]

    names = ["rho", "vx", "vy", "vz", "pre", "Bx", "By", "Bz"]
    data_dict = {name: arr[:, :, i+2] for i, name in enumerate(names)}

    return x, y, time, data_dict

def vecpot(x,y,Bx,By):
    Az = np.zeros_like(By)
    dy = y[1,0] - y[0,0]
    dx = x[0,1] - x[0,0]
    
    Az[1:, 0] = Az[0, 0] + np.cumsum(0.5*(Bx[1:, 0] + Bx[:-1, 0])*dy)
    Az[:, 1:] = Az[:, [0]] - np.cumsum(0.5*(By[:, 1:] + By[:, :-1])*dx, axis=1)
    
    return Az

dirname = sys.argv[1]
step = int(sys.argv[2])

foutname = dirname + "/snap%05d.dat"%(step) 
print("making plot ",foutname)

x, y, time, data = read_textdata(foutname)
Az = vecpot(x,y,data['Bx'],data['By'])

xmin = 0.0
xmax = 1.0
ymin = 0.0
ymax = 0.5

fig_width = 6.4
fig_height = fig_width*(ymax-ymin)/(xmax-xmin) + 1.9
fig = plt.figure(figsize=(fig_width, fig_height)) 

plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.xlabel("x axis") 
plt.ylabel("y axis") 

plt.text(0.5*(xmin+xmax),ymax*1.1,r"$\mathrm{time}=%.2f$"%(time),horizontalalignment="center")

im=plt.imshow(data['Bz'],extent=(xmin,xmax,ymin,ymax),origin="lower",vmin=-0.1,vmax=0.1)
cbar = plt.colorbar(im,orientation="horizontal")
cbar.set_label(r"$B_z$")

#plt.contour(x,y,Az,linestyles='solid',levels=20,colors="white")

for format_fig in ["pdf","png"]:
    outdir = dirname + "/" + format_fig + "file"
    os.makedirs(outdir, exist_ok=True)

    outputfile = outdir + "/snap%05d."%(step) + format_fig

    print("making plot file", outputfile)
    plt.savefig(outputfile,bbox_inches="tight")

plt.show()
