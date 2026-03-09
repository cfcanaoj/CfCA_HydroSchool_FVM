import sys
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

dirname = sys.argv[1]
step_s = int(sys.argv[2])
step_e = int(sys.argv[3])

fig = plt.figure()  
plt.xlabel("x axis") 
plt.ylabel("y axis") 

xmin = -0.5
xmax =  0.5
ymin = -1.0
ymax =  1.0

fname_anime = "animation.mp4"

graph_list = [] 
for istep in range(step_s,step_e+1):
    foutname = dirname + "/snap%05d.dat"%(istep) 

    print("making plot ",foutname)
    x, y, time, data = read_textdata(foutname)

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
