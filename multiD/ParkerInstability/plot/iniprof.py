import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import re

fig = plt.figure()  

plt.xlim(-15*np.pi, 15*np.pi)     
#plt.ylim(0, 0.25)     
#plt.ylim(0, 0.25)     
plt.yscale('log')
plt.xlabel("y",fontsize=15) 
plt.minorticks_on()

data = np.loadtxt("snap00000.dat")

iendx = 300
iendy = 600
x = data[:,0].reshape([iendy,iendx])
y = data[:,1].reshape([iendy,iendx])
den = data[:,2].reshape([iendy,iendx])
pre = data[:,6].reshape([iendy,iendx])


plt.plot( y[:,0], den[:,0], linewidth=4, label=r"$\rho$")
plt.plot( y[:,0], pre[:,0], linewidth=4, label=r"$P$")
plt.plot( y[:,0], pre[:,0]/den[:,0], linewidth=4, label=r"$T=P/\rho$")

plt.legend(ncol=1,handlelength=4)

plt.savefig("iniprof.pdf",bbox_inches="tight")
plt.show()
#    data = np.loadtxt(foutname)
#    print data
#    xv = data[:,0]
#    U  = data[:,1]


#    graph = plt.plot(xv, U, 'o-', color="black") 
#    graph_list.append(graph)               

#ani = animation.ArtistAnimation(fig, graph_list, interval=200) 
#plt.show()
