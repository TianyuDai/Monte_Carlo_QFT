import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#y = [3.6671260339, 1.6374308332, 0.3302833014, 0.0667231504]
#y = [3.66489, 1.63529, 0.336196, 0.0658161]
#x = [48, 64, 96, 128]
y = [0.0153983, 0.000848118, 4.77E-05, 1.62E-07]
x = [16, 24, 32, 48]
#x = [16, 24, 32, 48]
#y = [0.509644, 0.118301, 0.0272178, 0.00143473]
#y = [3.74834, 1.77304, 0.832098, 0.183132]
def func(x, a, b, c): 
	return a * np.exp(-b * x / 2.) + c
guess = [1, 0.1, 1]
popt, pcov = curve_fit(func, x, y, p0=guess, bounds = (0, [99, 1, 1]))
print(popt)
"""
x = [16, 24, 32, 48]
Lx = 8
kappa = 1./4.25
"""
