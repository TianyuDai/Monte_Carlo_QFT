import numpy as np
from numpy.linalg import inv
from scipy.optimize import curve_fit
C = []
Lt = [32, 48, 64, 80, 96, 128, 256]
LX = 8
kappa = 1./4.25
for LT in Lt: 
	M = [[0 for i in range(LX * LT)] for j in range(LX * LT)]
	for i in range(LX * LT): 
		M[i][i] = 1. 
	for i in range(LX): 
		for j in range(LT): 
			left = (LT + j - 1) % LT
			right = (j + 1) % LT
			up = (LX + i - 1) % LX
			down = (i + 1) % LX
			M[i+j*LX][i+left*LX] = -1. * kappa
			M[i+j*LX][i+right*LX] = -1. * kappa
			M[i+j*LX][up+j*LX] = -1. * kappa
			M[i+j*LX][down+j*LX] = -1. * kappa
	Minv = inv(M)
	C.append(sum(Minv[0][LT/2*LX:LT/2*LX+LX]))
print(C)
def func(x, a, b, c): 
	return a * np.exp(-b * x / 2.) + c
guess = [1, 0.1, 1]
popt, pcov = curve_fit(func, Lt, C, p0 = guess)
print(popt)
"""
k = np.arange(0.24907572, 0.24907574, 1e-09)
for kappa in k: 
	C_ = []
	for LT in Lt: 
		M = [[0 for i in range(LX * LT)] for j in range(LX * LT)]
		for i in range(LX * LT): 
			M[i][i] = 1. 
		for i in range(LX): 
			for j in range(LT): 
				left = (LT + j - 1) % LT
				right = (j + 1) % LT
				up = (LX + i - 1) % LX
				down = (i + 1) % LX
				M[i+j*LX][i+left*LX] = -1. * kappa
				M[i+j*LX][i+right*LX] = -1. * kappa
				M[i+j*LX][up+j*LX] = -1. * kappa
				M[i+j*LX][down+j*LX] = -1. * kappa
		Minv = inv(M)
		C_.append(sum(Minv[0][LT/2*LX:LT/2*LX+LX]))
	popt_, pcov_ = curve_fit(func, Lt, C_, p0 = guess)
	#if abs(popt[1]/4 - popt_[1]) < 1e-03: 
	print (kappa, abs(popt[1]/4 - popt_[1]))"""
