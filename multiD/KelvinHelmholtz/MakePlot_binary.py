import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import re

def read_bindata(file):
    with open(file, "rb") as fp: 
        time = np.fromfile(fp, np.float64, 1).item() 
        nx, ny, nhyd, nbc = [np.fromfile(fp, np.int32, 1).item() for _ in range(4)]

        xv = np.fromfile(fp, np.float64, nx)
        yv = np.fromfile(fp, np.float64, ny)
        Q  = np.fromfile(fp, np.float32, nx * ny * nhyd).reshape(ny, nx, nhyd)
        Bc = np.fromfile(fp, np.float32, nx * ny * nbc ).reshape(ny, nx, nbc)

    q_names = ['rho', 'vx', 'vy', 'vz', 'pre', 'sca']
    b_names = ['Bx', 'By', 'Bz']

    data_dict = {name: Q[:, :, i]  for i, name in enumerate(q_names)}
    data_dict.update({name: Bc[:, :, i] for i, name in enumerate(b_names)})

    return xv, yv, time, data_dict

def vecpot(x,y,Bx,By):
    Az = np.zeros_like(By)
    dy = y[1] - y[0]
    dx = x[1] - x[0]
    
    Az[1:, 0] = Az[0, 0] + np.cumsum(0.5*(Bx[1:, 0] + Bx[:-1, 0])*dy)
    Az[:, 1:] = Az[:, [0]] - np.cumsum(0.5*(By[:, 1:] + By[:, :-1])*dx, axis=1)
    
    return Az

###############################
dirname = sys.argv[1]
step = int(sys.argv[2])

foutname = dirname + "/snap%05d.bin"%(step) 
x, y, time, data = read_bindata(foutname)
Az = vecpot(x,y,data['Bx'],data['By'])

xmin = -0.5
xmax =  0.5
ymin = -1.0
ymax =  1.0

fig_height = 4.8
fig_width = fig_height*(xmax-xmin)/(ymax-ymin) + 1.9
fig = plt.figure(figsize=(fig_width, fig_height)) 

plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.xlabel("x axis") 
plt.ylabel("y axis") 

plt.text(0.5*(xmin+xmax),ymax*1.1,r"$\mathrm{time}=%.2f$"%(time),horizontalalignment="center")
im=plt.imshow(data['sca'],extent=(xmin,xmax,ymin,ymax),origin="lower",vmin=0,vmax=1)
cbar = plt.colorbar(im,orientation="vertical")
cbar.set_label("scalar field")

yy, xx = np.meshgrid(y,x,indexing="ij")
plt.contour(xx,yy,Az,linestyles='solid',levels=20,colors="white")

for format_fig in ["pdf","png"]:
    outdir = dirname + "/" + format_fig + "file"
    os.makedirs(outdir, exist_ok=True)

    outputfile = outdir + "/snap%05d."%(step) + format_fig

    print("making plot file", outputfile)
    plt.savefig(outputfile,bbox_inches="tight")

plt.show()
