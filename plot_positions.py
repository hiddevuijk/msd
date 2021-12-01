import numpy as np
import matplotlib.pyplot as plt
from sys import exit

sigmaC = 1.

L = 20.
data = np.loadtxt("r.dat")
 
figure, axes = plt.subplots()
 
axes.set_aspect( 1 )

for i in range(data.shape[0]):
  #plt.scatter( data[i,0], data[i,1])
  Drawing_uncolored_circle = plt.Circle( (data[i,0], data[i,1]), sigmaC/2,color="black",fill=None)
  
  axes.add_patch( Drawing_uncolored_circle)

plt.xlim([0,L])
plt.ylim([0,L])
plt.show()
