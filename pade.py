import re
import sys
import numpy as np
from scipy import misc
from matplotlib import pylab as plt
import cheb_pade as cheb

#trial function
def f_x(x):
  fx = 1.0/(1.5 - np.cos(5*x)) 
  return fx

#def f_x(x):
#	fx = 1.0/((x-2)*(x+2)) 
#	return fx

def freq(a,b,N):
  delta = (b-a)/float(N)
  return np.arange(a,b,delta)

#vidberg-serene determines pade coefficients
def pade_coeff(N,z,u):
  g = np.zeros((N,N), dtype=np.complex64)
  tmp1 = np.zeros(1, dtype=np.complex64)
  tmp2 = np.zeros(1, dtype=np.complex64)
  for p in range(1, N+1):
    if p==1:
      for i in range(1,N+1):
        g[p-1, i-1]=u[i-1]
    else:
      for i in range(p,N+1):
        tmp1   = g[p-2,p-2]/g[p-2,i-1]
        tmp2   = g[p-2,i-1]/g[p-2,i-1]
        g[p-1, i-1] = (tmp1 - tmp2)/(z[i-1] - z[p-2])

  a = np.zeros(N, dtype = np.complex64)
  a = np.diagonal(g)
  return a

def pade_eval(N,z,a,w):
  acap = np.zeros(N+1, dtype=np.complex64)
  bcap = np.zeros(N+1, dtype=np.complex64)
  padapp = np.zeros(1, dtype=np.complex64)

  acap[0] = complex(0.0, 0.0)
  acap[1] = complex(a[0], 0.0)
  bcap[0] = complex(1.0, 0.0)
  bcap[1] = complex(1.0, 0.0)

  for i in range(2,N+1):
    acap[i] = acap[i-1] + (w-z[i-2]) * a[i-1] * acap[i-2]
    bcap[i] = bcap[i-1] + (w-z[i-2]) * a[i-1] * bcap[i-2]

  padapp = acap[N]/bcap[N]
  return padapp

if __name__=="__main__":
  wmax = 1
  a1 = -1.0
  b1 = 1.0
#Original function:
  xrang = np.arange(-4,4,0.01, dtype=np.complex64)
	#xrang = xrang + complex(0.0, 0.01)
  fx = np.array([f_x(x) for x in xrang])
#create list of different order interpolants:
  rh_list=[]
  sample_list=[]
  N = [20]
  rez=[]
  for n in N:
    z = np.zeros(n, dtype=np.complex64)
    u = np.zeros(n, dtype=np.complex64)
    z = np.arange(a1,b1,(b1-a1)/n)
    rez.append(np.arange(a1,b1,(b1-a1)/n))
#sample function on imaginary axis.
    z = z * complex(0,1)
#find order N interpolant
    u = [f_x(x) for x in z]
    a = pade_coeff(n, z, u) 
    rh_list.append(np.array([pade_eval(n,z,a,x) for x in xrang], dtype=np.complex64))
    sample_list.append(np.array([f_x(x) for x in z], dtype=np.complex64))

#Interpolated function:
  fig = plt.figure()
  ax1 = fig.add_subplot(111)
  ax1.set_ylim([-2,2])
  ax1.plot(xrang,fx, '-x')
  fx1 = np.array([f_x(complex(0,1.0)*x) for x in xrang])
  ax1.plot(xrang,fx1)

#Plot real and imag part of f_x (blue and green), fitted function (red), 
#and the sampling points (dots)
  for i in range(len(N)):
    ax1.plot(xrang, np.real(rh_list[i]), label = 'N = %s'%(N[i]))
    ax1.plot(rez[i], np.real(sample_list[i]),'o', label = 'sample %s'%(N[i]))
  plt.legend()
  plt.show()


