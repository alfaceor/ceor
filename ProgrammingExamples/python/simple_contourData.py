import numpy as np
import matplotlib.pyplot as plt
D,Rg,PotentialEnergy = np.loadtxt("longtime.dat",usecols=[5,4,3]).T
PotentialEnergy,D,Rg=np.histogram2d(D, Rg, bins=(100,100), weights=PotentialEnergy)
print len(PotentialEnergy), len(PotentialEnergy[:,0]), len(D), len(Rg)
plt.contourf(D[0:len(D)-1],Rg[0:len(Rg)-1],PotentialEnergy)
#plt.contourf(D,Rg,PotentialEnergy)
plt.show()

"""
h2d, xe, ye = np.histogram2d(D,Rg,weights=PotentialEnergy)
h2d
xe
xe[:-1]
xe[1:]
x=(xe[:-1]+xe[1:])/2
x
y = (ye[:-1] + ye[1:])/2
plt.contourf(x,y,h2d.T,origin='lower')
"""

