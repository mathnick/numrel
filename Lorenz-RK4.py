import numpy as np
from matplotlib import pyplot as plt

def f(t, x, y, z):
    return s*y - s*x
def g(t, y, x, z):
    return r*x -y - x*z
def m(t, z, y, x):
    return x*y - b*z

def RK4(f, g, m, x0, y0,z0, t0, h, n):

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
        k4 = h*f(t[k]+h, x[k] + k3, y[k] + k3, z[k] + k2/2)
        kk1 = h * g(t[k], y[k], x[k], z[k])
        kk2 = h * g(t[k] + h / 2, y[k] + kk1 / 2, x[k] + kk1/2, z[k] + kk1 / 2)
        kk3 = h * g(t[k] + h / 2, y[k] + kk2 / 2, x[k] + kk2/2, z[k] + kk2 / 2)
        kk4 = h * g(t[k] + h, y[k] + kk3, x[k] + kk3, z[k] + kk3)
        kkk1 = h * m(t[k], z[k], y[k], x[k])
        kkk2 = h * m(t[k] + h / 2, z[k] + kkk1 / 2, y[k] + kkk1 / 2, x[k] + kkk1 / 2)
        kkk3 = h * m(t[k] + h / 2, z[k] + kkk2 / 2, y[k] + kkk2 / 2, x[k] + kkk2 / 2)
        kkk4 = h * m(t[k] + h, z[k] + kkk3, y[k] + kkk3, x[k] + kkk3)
        x[k+1] = x[k] + (1/6)*(k1 + 2*k2 + 2*k3 +k4)
        y[k+1] = y[k] + (1/6)*(kk1 + 2*kk2 + 2*kk3 + kk4)
        z[k+1] = z[k] + (1 / 6) * (kkk1 + 2 * kkk2 + 2 * kkk3 + kkk4)
        t[k+1] = t[k] + h
        print(x[k+1], y[k+1], z[k+1], t[k+1])
    return x, y, z, t


x0=1
y0=1
z0=1
t0=0
h=0.01
n=2000
s=10
b=8/3
r=28

x, y, z, t = RK4(f, g, m, x0, y0, z0, t0, h, n)

plt.figure()
plt.plot(t,x)
plt.plot(t,y)
plt.plot(t,z)

plt.show()
