import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import re

def read_textdata(file):
    with open(foutname, 'r') as data_file:
        line = data_file.readline();
        attributes1 = re.findall(r"\d+\.\d+", line)

        line = data_file.readline();
        attributes2 = re.findall(r"\d+", line)

    time = float(attributes1[0]) 
    nx = int(attributes2[0])
    ny = int(attributes2[1])

    data = np.loadtxt(foutname)

    x = data[:,0].reshape(ny,nx)
    y = data[:,1].reshape(ny,nx)

    data_dict = {}
    data_dict['rho'] = data[:,2].reshape(ny,nx)
    data_dict['vx'] = data[:,3].reshape(ny,nx)
    data_dict['vy'] = data[:,4].reshape(ny,nx)
    data_dict['vz'] = data[:,5].reshape(ny,nx)
    data_dict['pre'] = data[:,6].reshape(ny,nx)
    data_dict['sca'] = data[:,7].reshape(ny,nx)
    data_dict['Bx'] = data[:,8].reshape(ny,nx)
    data_dict['By'] = data[:,9].reshape(ny,nx)
    data_dict['Bz'] = data[:,10].reshape(ny,nx)
#    data_dict['psi'] = data[:,11].reshape(ny,nx)

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

fig = plt.figure()  
plt.xlim(0, 1)     
plt.ylim(0, 1)
plt.xlabel("x axis") 
plt.ylabel("y axis") 

xmin = -0.5
xmax =  0.5
ymin = -1.0
ymax =  1.0

plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)

pg00 = plt.text(0.5*(xmin+xmax),ymax*1.1,r"$\mathrm{time}=%.2f$"%(time),horizontalalignment="center")

im=plt.imshow(data['sca'],extent=(xmin,xmax,ymin,ymax),origin="lower",vmin=0,vmax=1)
cbar = plt.colorbar(im,orientation="vertical")
cbar.set_label("scalar field")

plt.contour(x,y,Az,linestyles='solid',levels=20,colors="white")

for format_fig in ["pdf","png"]:
    outdir = dirname + "/" + format_fig + "file"
    os.makedirs(outdir, exist_ok=True)

    outputfile = outdir + "/snap%05d."%(step) + format_fig

    print("making plot file", outputfile)
    plt.savefig(outputfile,bbox_inches="tight")

plt.show()
