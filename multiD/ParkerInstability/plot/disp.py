import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import re

fig = plt.figure()  

plt.xlim(0, 0.5)     
plt.ylim(0, 0.25)     
plt.xlabel("k",fontsize=15) 
plt.ylabel(r"$\sigma$",fontsize=15) 
plt.minorticks_on()

data = np.loadtxt("disp_gam1.05_beta1.txt")
plt.plot( data[:,0], data[:,1], linewidth=4, label=r"$\gamma=1.05,\beta_0=1$", color="tab:red" )

data = np.loadtxt("disp_gam1.05_beta0.1.txt")
plt.plot( data[:,0], data[:,1], '--',linewidth=4, label=r"$\gamma=1.05,\beta_0=0.1$", color="tab:blue" )

data = np.loadtxt("disp_gam1.05_beta10.txt")
plt.plot( data[:,0], data[:,1], ':',linewidth=4, label=r"$\gamma=1.05,\beta_0=10$", color="tab:blue" )

data = np.loadtxt("disp_gam1.0_beta1.txt")
plt.plot( data[:,0], data[:,1], '--',linewidth=4, label=r"$\gamma=1,\beta_0=1$", color="tab:orange" )

data = np.loadtxt("disp_gam1.4_beta1.txt")
plt.plot( data[:,0], data[:,1], ':',linewidth=4, label=r"$\gamma=1.4,\beta_0=1$", color="tab:orange" )

plt.legend(ncol=2,handlelength=4)

plt.savefig("disp_parker.pdf",bbox_inches="tight")
plt.show()
#    data = np.loadtxt(foutname)
#    print data
#    xv = data[:,0]
#    U  = data[:,1]


#    graph = plt.plot(xv, U, 'o-', color="black") 
#    graph_list.append(graph)               

#ani = animation.ArtistAnimation(fig, graph_list, interval=200) 
#plt.show()
