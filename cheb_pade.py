import re
import sys
import numpy as np
from scipy import misc, fftpack as dft
from matplotlib import pylab as plt

def f_x(x):
	fx = 1.0/(1.5 - np.cos(5*x)) 
	return fx

def cheb_points(m, n, wmax):
	N = m + n
	points = np.array([(np.cos((2*float(j)+1)*np.pi/(2*N+2))) for j in range(N+1)], dtype=np.complex64)
	points = points[::-1]
	nodes = points

	for i in points:
		print '0.0000', '\t%10.5f'%(i)
        print 'Affine Transformation:'
	a = -wmax
	b = wmax
	points = [(a+b)/2.0 + ((b-a)/2.0)*(np.cos((2*float(j)+1)*np.pi/(2*N+2))) for j in range(N+1)]
	points = points[::-1]
#	for i in points:
#		print '0.0000', '\t%10.5f'%(i)

	return nodes

def bary(x, fx, nodes, qw):
	y = np.zeros((x.shape), dtype=np.complex64)
	fxw = np.multiply(fx[:], qw[:])
	print qw
	for i in range(len(x)):
		dxinv = np.divide(1.0, (x[i] - nodes[:]))
		#ind   = np.nonzero((np.isfinite(dxinv[:])) 
#		ind   = np.nonzero(dxinv[:]) 
		ind   = np.nonzero([not(a) for a in np.isfinite(dxinv[:])])
		if len(ind) > 1:
			print ind
			y[i] = fx[ind]
		else:
			y[i] =  np.dot(fxw, dxinv)/np.dot(qw, dxinv)
	return y

def ratinterp(fx, nodes, m, n):
	N  = m + n
	D1 = np.diag(fx)
	#if complex might have to fft parts separately
	D  = dft.dct(D1[:,:].real, norm='ortho')
	D = np.transpose(D)
	temp = D[:n,:]
	Z  = dft.dct(temp[:n,:], norm='ortho')
	Z  = np.transpose(Z)
	#Z  = dft.dct(np.transpose(D[:n,:]), norm='ortho')
	#Z  = dft.dct(D[:n,:], norm='ortho')
	print Z
#	Z  = dft.dct(D[:n,:], norm='ortho')
	u, s, v = np.linalg.svd(Z[m+1:N,:])
	#u, s, v = np.linalg.svd(Z[:,m+1:N])
	v = v[-1,:]
	#v = v[:,-1]
	v = np.append(v, np.zeros((N+1)-len(v)))
	len(v)
	q = dft.idct(v, norm='ortho')
	w = np.array([np.power(-1.0,j)*(np.sin((2*float(j)+1)*np.pi/(2*N+2))) for j in range(N+1)], dtype=np.complex64)
	print v
	print 'q', q
	print 'w', w
	qw = np.multiply(q,w)
	return lambda x: bary(x, fx, nodes, qw)

if __name__ == '__main__':
	m = 10
	n = 10
	wmax = 1
	nodes = cheb_points(m,n,wmax)

	fx = np.array([f_x(x) for x in nodes], dtype=np.complex64)
	rh = ratinterp(fx, nodes, m, n)

	xrang = np.arange(-2,2,0.01)
	fx = np.array([f_x(x) for x in xrang])

	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ax1.set_ylim([0,3])
	ax1.plot(xrang,fx,'-x')
	ax1.plot(xrang,rh(xrang),'ro')
	plt.show()
	
