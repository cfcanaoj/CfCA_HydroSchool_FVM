import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import re

def read_bindata(file):

    with open(foutname, "rb") as fp: 
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

dirname = sys.argv[1]
step_s = int(sys.argv[2])
step_e = int(sys.argv[3])

fig = plt.figure()  
plt.xlim(0, 1)     
plt.ylim(0, 1)
plt.xlabel("x axis") 
plt.ylabel("y axis") 

xmin = -0.5
xmax =  0.5
ymin = -1.0
ymax =  1.0

fname_anime = "animation.mp4"

graph_list = [] 
for istep in range(step_s,step_e+1):
    foutname = dirname + "/snap%05d.bin"%(istep) 
    x, y, time, data = read_bindata(foutname)

    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)

    pg00 = plt.text(0.5*(xmin+xmax),ymax*1.1,r"$\mathrm{time}=%.2f$"%(time),horizontalalignment="center")

    im=plt.imshow(data['sca'],extent=(xmin,xmax,ymin,ymax),origin="lower",vmin=0,vmax=1)

    if istep == step_s: 
        cbar = plt.colorbar(im,orientation="vertical")
        cbar.set_label("scalar field")
    graph_list.append([pg00,im])               

ani = animation.ArtistAnimation(fig, graph_list, interval=200) 
print("making animation file", fname_anime)
ani.save(dirname + "/" + fname_anime, writer="imagemagick")
plt.show()
