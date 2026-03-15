import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import re

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
ymin = -0.5
ymax =  0.5

fname_anime = "animation.mp4"

graph_list = [] 
for istep in range(step_s,step_e+1):
    foutname = dirname + "/snap%05d.bin"%(istep) 
    fp = open(foutname,"rb")

    time = np.fromfile(fp,np.float64,1)[0]
    nx   = np.fromfile(fp,np.int32,1)[0]
    ny   = np.fromfile(fp,np.int32,1)[0]
    nhyd = np.fromfile(fp,np.int32,1)[0]
    nbc  = np.fromfile(fp,np.int32,1)[0]
    xv   = np.fromfile(fp,np.float64,nx)
    yv   = np.fromfile(fp,np.float64,ny)
    Q    = np.fromfile(fp,np.float32,nx*ny*nhyd).reshape(ny,nx,nhyd)
    Bc   = np.fromfile(fp,np.float32,nx*ny*nbc).reshape(ny,nx,nbc)

    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)

    den = Q[:,:,0]

    pg00 = plt.text(0.5*(xmin+xmax),ymax*1.1,r"$\mathrm{time}=%.3f$"%(time),horizontalalignment="center")

    im=plt.imshow(den[:,:],extent=(xmin,xmax,ymin,ymax),origin="lower",vmin=0,vmax=3)


    if istep == step_s: 
        plt.colorbar(im,orientation="vertical")
    graph_list.append([pg00,im])               


ani = animation.ArtistAnimation(fig, graph_list, interval=200) 
print("making animation file", fname_anime)
ani.save(dirname + "/" + fname_anime, writer="imagemagick")
plt.show()
