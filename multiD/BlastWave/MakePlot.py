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

    return x[0,:], y[:,0], time, data_dict

def read_bindata(file):
    with open(file, "rb") as fp: 
        time = np.fromfile(fp, np.float64, 1).item() 
        nx, ny, nhyd, nbc = [np.fromfile(fp, np.int32, 1).item() for _ in range(4)]

        xv = np.fromfile(fp, np.float64, nx)
        yv = np.fromfile(fp, np.float64, ny)
        Q  = np.fromfile(fp, np.float32, nx * ny * nhyd).reshape(ny, nx, nhyd)
        Bc = np.fromfile(fp, np.float32, nx * ny * nbc ).reshape(ny, nx, nbc)

    q_names = ['rho', 'vx', 'vy', 'vz', 'pre']
    b_names = ['Bx', 'By', 'Bz']

    data_dict = {name: Q[:, :, i]  for i, name in enumerate(q_names)}
    data_dict.update({name: Bc[:, :, i] for i, name in enumerate(b_names)})

    return xv, yv, time, data_dict

def vecpot(x,y,Bx,By):
    Az = np.zeros_like(By)
    dy = y[1] - y[0]
    dx = x[0] - x[0]
    
    Az[1:, 0] = Az[0, 0] + np.cumsum(0.5*(Bx[1:, 0] + Bx[:-1, 0])*dy)
    Az[:, 1:] = Az[:, [0]] - np.cumsum(0.5*(By[:, 1:] + By[:, :-1])*dx, axis=1)
    
    return Az

filetype = sys.argv[1]
dirname = sys.argv[2]
step = int(sys.argv[3])


if filetype == "ascii": 
    foutname = dirname + "/snap%05d.dat"%(step) 
    print("making plot ",foutname)
    x, y, time, data = read_textdata(foutname)
elif filetype == "binary":
    foutname = dirname + "/snap%05d.bin"%(step) 
    print("making plot ",foutname)
    x, y, time, data = read_bindata(foutname)
else:
    print("error! filetype should be ascii or binary")
    exit(-1)

Az = vecpot(x,y,data['Bx'],data['By'])

xmin = -0.5
xmax =  0.5
ymin = -0.5
ymax =  0.5

fig_height = 4.8
fig_width = fig_height*(xmax-xmin)/(ymax-ymin) + 1.9
fig = plt.figure(figsize=(fig_width, fig_height)) 

plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.xlabel("x axis") 
plt.ylabel("y axis") 

plt.text(0.5*(xmin+xmax),ymax*1.1,r"$\mathrm{time}=%.2f$"%(time),horizontalalignment="center")

im=plt.imshow(data['rho'],extent=(xmin,xmax,ymin,ymax),origin="lower",vmin=0,vmax=3)
cbar = plt.colorbar(im,orientation="vertical")
cbar.set_label("density")

plt.contour(x,y,Az,linestyles='solid',levels=20,colors="white")

for format_fig in ["pdf","png"]:
    outdir = dirname + "/" + format_fig + "file"
    os.makedirs(outdir, exist_ok=True)

    outputfile = outdir + "/snap%05d."%(step) + format_fig

    print("making plot file", outputfile)
    plt.savefig(outputfile,bbox_inches="tight")

plt.show()
