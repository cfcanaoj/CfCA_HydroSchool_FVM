import sys
import numpy as np
import matplotlib.pyplot as plt
import re

dirname = sys.argv[1]
step_s = int(sys.argv[2])
step_e = int(sys.argv[2])

xmin = -7.5*np.pi
xmax =  7.5*np.pi
ymin = -15*np.pi
ymax =  15*np.pi

for step in range(step_s,step_e+1):

    fig = plt.figure()  
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    plt.xlabel("x axis") 
    plt.ylabel("y axis") 
    
    foutname = dirname + "/snap%05d.bin"%(step) 
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
    
    Az = np.zeros([ny,nx])
    
    dy = yv[1] - yv[0]
    dx = xv[1] - xv[0]
    
    for j in range(0,ny-1):
        for i in range(0,1):
            Az[j+1,i] = Az[j,i] + 0.5*( Bc[j+1,i,0] + Bc[j,i,0] )*dy
    
    for j in range(0,ny):
        for i in range(0,nx-1):
            Az[j,i+1] = Az[j,i] - 0.5*( Bc[j,i+1,1] + Bc[j,i,1] )*dy
    
    den = Q[:,:,0]
    
    yy, xx = np.meshgrid(yv,xv,indexing="ij")
    
    print(den)
    
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    
    plt.text(0.5*(xmin+xmax),ymax*1.1,r"$\mathrm{time}=%.2f$"%(time),horizontalalignment="center")
    im=plt.imshow(np.log10(den[:,:]),extent=(xmin,xmax,ymin,ymax),origin="lower",vmin=-5,vmax=0)
    
    plt.contour(xx,yy,Az,linestyles='solid',levels=20,colors="white")
    
    plt.colorbar(im,orientation="vertical")
    
    outputfile = dirname + "/snap%05d.pdf"%(step) 
    
    print("making plot file", outputfile)
    plt.savefig(outputfile,bbox_inches="tight")

    plt.show()
