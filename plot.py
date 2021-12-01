import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp
from sys import exit
def D(phi, k):
  A = 1- 3 *k*k
  B = 1+k*k

  D = 1 - 2 * phi * A / B
  D/= 1 + k * k

  return D

def f(x,a,b):
  return a+b*x


diff = np.loadtxt("msd.dat").T
N = diff.shape[1] - 1
t = diff[:,0]
n = 10
dt = t[-1] - t[n]

Dlist = []
D2list = []
for i in range(1,N):
  popt, pcov = sp.curve_fit(f, t[n:], diff[n:,i])
  D2 = popt[1]/2
  D2list.append(D2)

  D = (diff[-1,i] - diff[n,i]) / dt
  Dlist.append(D/2)

  plt.scatter(t,diff[:,i])

Dlist = np.asarray(Dlist)
D = np.average(Dlist)
err = np.std(Dlist)/np.sqrt(Dlist.shape[0])
print(D, err)

D2list = np.asarray(D2list)
D2 = np.average(D2list)
err = np.std(D2list)/np.sqrt(D2list.shape[0])
print(D2,err)

plt.title("D={:1.5f}  ({:1.5f})".format(D2, err))
plt.plot(t[n:],2*D*t[n:], color="black")
plt.show()

