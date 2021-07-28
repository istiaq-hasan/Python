import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.special import erf
from scipy.integrate import quad
import scipy.integrate as integrate

a = 0
b = 1
n = 10

dx = (b - a) / n

x = np.linspace(a, b, n + 1)
y = 3- x + 2 * x ** 2 + x ** 3

def z(x):
    z = 3- x + 2 * x ** 2 + x ** 3
    return z

y1 = y[1 : ]
y2 = y[ : -1]

A = (dx / 2) * np.sum(y1 + y2)
B = np.trapz(y, x, dx)
C, err = integrate.quad(z, a, b)

print(A, B, C)

def integrand(t):
        return np.exp(-(t ** 2))

def damp(alpha, rho):
        return quad(integrand, alpha, rho) [0]

def norm_constant():
        return 2.0 / math.sqrt(math.pi)

a = 0

x = np.arange(-2, 2, 0.01)

e = [norm_constant() * damp(a,b) for b in x]


plt.plot(x,e,'r')
plt.grid()
plt.show()