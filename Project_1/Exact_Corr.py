import numpy as np
from numpy.linalg import inv
LX = 16
LT = 64
kappa = 0.249376559
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
print(sum(Minv[0][LT/2*LX:LT/2*LX+LX]))
