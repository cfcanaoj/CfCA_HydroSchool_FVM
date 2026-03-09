import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import re

dirname = sys.argv[1]
step_s = int(sys.argv[2])

def makedirs(path):
    if not os.path.isdir(path):
        os.makedirs(path)

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
    data_dict['psi'] = data[:,11].reshape(ny,nx)

    return x, y, time, data_dict

def vecpot(x,y,Bx,By):
    Az = np.zeros_like(By)
    dy = y[1,0] - y[0,0]
    dx = x[0,1] - x[0,0]
    
    Az[1:, 0] = Az[0, 0] + np.cumsum(0.5*(Bx[1:, 0] + Bx[:-1, 0])*dy)
    Az[:, 1:] = Az[:, [0]] - np.cumsum(0.5*(By[:, 1:] + By[:, :-1])*dx, axis=1)
    
    return Az


fig = plt.figure()  
plt.xlim(0, 1)     
plt.ylim(0, 1)
plt.xlabel("x axis") 
plt.ylabel("y axis") 

xmin = -0.5
xmax =  0.5
ymin = -1.0
ymax =  1.0

for istep in range(step_s,step_s+1):
    foutname = dirname + "/snap%05d.dat"%(istep) 
    print("making plot ",foutname)

    x, y, time, data = read_textdata(foutname)

    Az = vecpot(x,y,data['Bx'],data['By'])

    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)

    pg00 = plt.text(0.5*(xmin+xmax),ymax*1.1,r"$\mathrm{time}=%.2f$"%(time),horizontalalignment="center")

    im=plt.imshow(data['sca'],extent=(xmin,xmax,ymin,ymax),origin="lower",vmin=0,vmax=1)

    plt.contour(x,y,Az,linestyles='solid',levels=20,colors="white")

    if istep == step_s: 
        plt.colorbar(im,orientation="vertical")

    makedirs(dirname + "/pdffile")
    plt.savefig(dirname + "/pdffile/snap%05d.pdf"%(istep),bbox_inches="tight")

plt.show()
