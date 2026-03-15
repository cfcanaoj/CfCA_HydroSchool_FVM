import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import re

fig = plt.figure()  

plt.xlim(0, 50)     
plt.ylim(1e-5, 0.5)     
plt.yscale('log')
plt.xlabel("time") 
plt.minorticks_on()
plt.ylabel(r"$\sqrt{\langle B_y^2\rangle}$") 

data = np.loadtxt("ana.dat")
plt.plot( data[:,0], data[:,1], linewidth=4, label=r"$\gamma=1.05$" )

plt.plot( data[:,0], np.exp(data[:,0]*0.157)*1.9e-4,'--', color="black" )

data = np.loadtxt("../hdc_gam1.4/ana.dat")
plt.plot( data[:,0], data[:,1], linewidth=4, label=r"$\gamma=1.4$"  )

plt.plot( data[:,0], np.exp(data[:,0]*8.e-2)*4.3e-4,'--', color="black" )

plt.legend()

plt.savefig("growth_parker.pdf",bbox_inches="tight")
plt.show()
#    data = np.loadtxt(foutname)
#    print data
#    xv = data[:,0]
#    U  = data[:,1]


#    graph = plt.plot(xv, U, 'o-', color="black") 
#    graph_list.append(graph)               

#ani = animation.ArtistAnimation(fig, graph_list, interval=200) 
#plt.show()
