import numpy as np
from numpy.linalg import inv
LX = 4
LT = 16
kappa = 0.124843945
M = [[0 for i in range(LX * LX * LX * LT)] for j in range(LX * LX * LX * LT)]
for i in range(LX * LX * LX * LT): 
	M[i][i] = 1. 
for i in range(LX): 
	for j in range(LX): 
		for k in range(LX): 
			for l in range(LT): 
				nb1 = (LX + i - 1) % LX
				nb2 = (i + 1) % LX
				nb3 = (LX + j - 1) % LX
				nb4 = (j + 1) % LX
				nb5 = (LX + k - 1) % LX
				nb6 = (k + 1) % LX
				nb7 = (LT + l - 1) % LT
				nb8 = (l + 1) % LT
				M[i+j*LX+k*LX*LX+l*LX*LX*LX][nb1+j*LX+k*LX*LX+l*LX*LX*LX] = -1. * kappa
				M[i+j*LX+k*LX*LX+l*LX*LX*LX][nb2+j*LX+k*LX*LX+l*LX*LX*LX] = -1. * kappa
				M[i+j*LX+k*LX*LX+l*LX*LX*LX][i+nb3*LX+k*LX*LX+l*LX*LX*LX] = -1. * kappa
				M[i+j*LX+k*LX*LX+l*LX*LX*LX][i+nb4*LX+k*LX*LX+l*LX*LX*LX] = -1. * kappa
				M[i+j*LX+k*LX*LX+l*LX*LX*LX][i+j*LX+nb5*LX*LX+l*LX*LX*LX] = -1. * kappa
				M[i+j*LX+k*LX*LX+l*LX*LX*LX][i+j*LX+nb6*LX*LX+l*LX*LX*LX] = -1. * kappa
				M[i+j*LX+k*LX*LX+l*LX*LX*LX][i+j*LX+k*LX*LX+nb7*LX*LX*LX] = -1. * kappa
				M[i+j*LX+k*LX*LX+l*LX*LX*LX][i+j*LX+k*LX*LX+nb8*LX*LX*LX] = -1. * kappa
Minv = inv(M)
print(sum(Minv[0][LT/2*LX*LX*LX:LT/2*LX*LX*LX+LX*LX*LX]))
