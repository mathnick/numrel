import numpy as np
import matplotlib.pyplot as plt


n = 200

r = np.zeros(n)
m = np.zeros(n)
# rho = np.zeros(n)
# p = np.zeros(n)
r[0] = r0
m[0] = m0
# rho[0] = rho0
# p[0] = p0


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
#Equação do gradiente de p

def f(gradP,rho):
    res = -gradP/rho + gradG
    return res
#Equação da continuidade

def g(gradP,rho):
    res = v
    return res
#Equação da derivada das posições
#A partir desse momento, utilizo runge-kutta para descobrir os valores da velocidade e posição ?
#Como achar os valores de m ?
