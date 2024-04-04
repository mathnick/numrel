import numpy as np
from matplotlib import pyplot as plt

h = 0,1
n = 8
G = 6.67*(10**(-11))
r = np.zeros(n)
F = np.zeros(n)

def dis(x,y):
    for i in range (n-1):
        for j in range (n-1):
            if i!= j:
                r[i,j] = ((x[i]-x[j])**2 + (y[i]-y[j])**2)**(1/2)
                print(r[i,j])
        return r


r = dis(x, y)


def W(r):
    res = (1 / (h * (np.pi) ** (3 / 2))) * (np.exp(-(r / h)**2))
    return res
#Equação do kernel

def rho(r):
    res = m*W(r,h)
    return res
#Equação da densidade

def P(rho):
    res = rho**(5/3)
    return res
#Equação da pressão, esta coreeta ?

def gradW(r):
    res = (1 / (h * (np.pi) ** (3 / 2))) * 2 * (np.exp(-(r / h)**2)) * (x/h**2)
    return res
#Equação do gradiente de W

def gradP(r):
    res = np.sum(m*(p/rho)*gradW)
    return res

def Fg(r,m):
    for i in range(n-1):
        for j in range(n-1):
            if i != j:
                F[i,j] = (G*m[i]*m[j])/(r**2)
        return F

F = Fg(r,m)
