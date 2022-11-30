import numpy as np
from matplotlib import pyplot as plt

np.set_printoptions(precision=16)

def f(t, x, y, z):
    return s*y - s*x
def g(t, x, y, z):
    return r*x -y - x*z
def m(t, x, y, z):
    return x*y - b*z

def RK4(f, x0, y0,z0, t0, h, n):

    x = np.zeros(n)
    y = np.zeros(n)
    z = np.zeros(n)
    t = np.zeros(n)
    x[0] = x0
    y[0] = y0
    z[0] = z0
    t[0] = t0

    for k in range(n-1):
        k1 = h*f(t[k], x[k], y[k], z[k])
        k2 = h*f(t[k]+h/2, x[k] + k1/2, y[k] + k1/2, z[k] + k1/2)
        k3 = h*f(t[k]+h/2, x[k] + k2/2, y[k] + k2/2, z[k] + k2/2)
        k4 = h*f(t[k]+h, x[k] + k3, y[k] + k3, z[k] + k3)
        kk1 = h * g(t[k], x[k], y[k], z[k])
        kk2 = h * g(t[k] + h / 2, x[k] + kk1 / 2, y[k] + kk1/2, z[k] + kk1 / 2)
        kk3 = h * g(t[k] + h / 2, x[k] + kk2 / 2, y[k] + kk2/2, z[k] + kk2 / 2)
        kk4 = h * g(t[k] + h, x[k] + kk3, y[k] + kk3, z[k] + kk3)
        kkk1 = h * m(t[k], x[k], y[k], z[k])
        kkk2 = h * m(t[k] + h / 2, x[k] + kkk1 / 2, y[k] + kkk1 / 2, z[k] + kkk1 / 2)
        kkk3 = h * m(t[k] + h / 2, x[k] + kkk2 / 2, y[k] + kkk2 / 2, z[k] + kkk2 / 2)
        kkk4 = h * m(t[k] + h, x[k] + kkk3, y[k] + kkk3, z[k] + kkk3)
        x[k+1] = x[k] + (1 / 6) * (k1 + 2*k2 + 2*k3 +k4)
        y[k+1] = y[k] + (1 / 6) * (kk1 + 2*kk2 + 2*kk3 + kk4)
        z[k+1] = z[k] + (1 / 6) * (kkk1 + 2 * kkk2 + 2 * kkk3 + kkk4)
        t[k+1] = t[k] + h
        print(x[k+1], y[k+1], z[k+1], t[k+1])
    return x, y, z, t


x0=1.0
y0=1.0
z0=20.01
t0=0
T=14.0
n=1000
s=10.0
b=8.0/3.0
r=28.0
h=T/n

x, y, z, t = RK4(f, x0, y0, z0, t0, h, n)

plt.figure()
plt.plot(t,x)
plt.title('Runge Kutta 4-1 - x')
#plt.plot(t,y)
#plt.plot(t,z)

plt.show()
