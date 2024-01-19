import numpy as np
import matplotlib.pyplot as plt

evoldat = np.loadtxt('evol.dat')

evoldat = np.reshape(evoldat,(1000,200,4))

evolgens = evoldat[:,:,0]
evolbest = evoldat[:,:,1]
evolalltimebest = evolbest[200::201]
print("should be best fitnesses: ",evolalltimebest[0:10])
evolavg = evoldat[:,:,2]
evolvar = evoldat[:,:,3]

# number of runs reaching some fitness
bins = [0,.05,.01,.15,.20,.25,.3,1.3,2.3,3.3,4.3,5.3,6.3,7.3,8.3,9.3,10.3,11.3,12.3,13.3,14.3,15.3,16.3,17.3,18.3,19.3,20.3]
plt.hist(evolalltimebest,bins=bins)