import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp
from sys import exit

N = 320
L = 80
k = 1.
phi = N*np.pi/(4*L*L)
print(phi)

def t(k):
  A = 1- 3 *k*k
  B = 1+k*k
  return A/B

D = 1 - 2 * phi*t(k)
D /= 1 + k*k

print(D)


def f(x,a,b):
  return a+b*x

data = np.loadtxt("msd.dat")
t = data[:,0]
msd = data[:,1]

plt.plot(t, 4*D*t, color="red")

N = int(msd.shape[0]/5)

popt, pcov = sp.curve_fit(f, t[N:], msd[N:])

plt.plot(t[N:], popt[0] + popt[1]*t[N:], color='gray')

D = popt[1]/4.
print(D)
plt.plot(t, msd)
plt.show()
