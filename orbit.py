import numpy as np
from matplotlib import pyplot as plt



def f(t, x, y, vx, vy):
    return vx

def g(t, y, x, vy, vx):
    return vy

def i(t, x, y, vx, vy):
    return (G*M*x)/(((x ** 2) + (y ** 2)) ** (3 / 2))

def j(t, y, x, vy, vx):
    return (G*M*y)/(((x ** 2) + (y ** 2)) ** (3 / 2))


def RK4(i, j, x0, y0, vx0, vy0, t0, h, n):
    x = np.zeros(n)
    y = np.zeros(n)
    vx = np.zeros(n)
    vy = np.zeros(n)
    t = np.zeros(n)
    x[0] = x0
    y[0] = y0
    vx[0] = vx0
    vy[0] = vy0
    t[0] = t0

    for k in range(n-1):
        kvx1 = h * f(t[k], x[k], y[k], vx[k], vy[k])
        kvx2 = h * f(t[k] + h / 2, x[k] + kvx1 / 2, y[k] + kvx1 / 2, vx[k] + kvx1 / 2, vy[k] + kvx1 / 2)
        kvx3 = h * f(t[k] + h / 2, x[k] + kvx2 / 2, y[k] + kvx2 / 2, vx[k] + kvx2 / 2, vy[k] + kvx2 / 2)
        kvx4 = h * f(t[k] + h, x[k] + kvx3, y[k] + kvx3, vx[k] + kvx3 / 2, vy[k] + kvx3)

        kvy1 = h * j(t[k], y[k], x[k], vy[k], vx[k])
        kvy2 = h * j(t[k] + h / 2, y[k] + kvy1 / 2, x[k] + kvy1 / 2, vy[k] + kvy1 / 2, vx[k] + kvy1 / 2)
        kvy3 = h * j(t[k] + h / 2, y[k] + kvy2 / 2, x[k] + kvy2 / 2, vy[k] + kvy2 / 2, vx[k] + kvy2 / 2)
        kvy4 = h * j(t[k] + h, y[k] + kvy3, x[k] + kvy3, vy[k] + kvy3, vx[k] + kvy3)

        kx1 = h * i(t[k], x[k], y[k], vx[k], vy[k])
        kx2 = h * i(t[k] + h / 2, x[k] + kx1 / 2, y[k] + kx1 / 2, vx[k] + kx1 / 2, vy[k] + kx1 / 2)
        kx3 = h * i(t[k] + h / 2, x[k] + kx2 / 2, y[k] + kx2 / 2, vx[k] + kx2 / 2, vy[k] + kx2 / 2)
        kx4 = h * i(t[k] + h, x[k] + kx3, y[k] + kx3, vx[k] + kx3/2, vy[k] + kx3)

        ky1 = h * j(t[k], y[k], x[k], vy[k], vx[k])
        ky2 = h * j(t[k] + h / 2, y[k] + ky1 / 2, x[k] + ky1 / 2, vy[k] + ky1 / 2, vx[k] + ky1 / 2)
        ky3 = h * j(t[k] + h / 2, y[k] + ky2 / 2, x[k] + ky2 / 2, vy[k] + ky2 / 2, vx[k] + ky2 / 2)
        ky4 = h * j(t[k] + h, y[k] + ky3, x[k] + ky3, vy[k] + ky3, vx[k] +ky3)

        x[k+1] = x[k] + (1/6)*(kx1 + 2*kx2 + 2*kx3 +kx4)
        y[k+1] = y[k] + (1/6)*(ky1 + 2*ky2 + 2*ky3 + ky4)
        vx[k+1] = vx[k] + (1/6) * (kvx1 + 2 * kvx2 + 2 * kvx3 + kvx4)
        vy[k+1] = vy[k] + (1/6) * (kvy1 + 2 * kvy2 + 2 * kvy3 + kvy4)
        t[k+1] = t[k] + h
        print(x[k+1], y[k+1], vx[k+1], vy[k+1], t[k+1])
    return x, y, vx, vy, t


x0 = 15
y0 = 1
vx0 = 1000
vy0 = 1000
t0 = 0
h = 0.01
n = 2000
G = 6.67408*(10^(-11))
M = 2*(10^30)


x, y, vx, vy, t = RK4(i, j, x0, y0, vx0, vy0, t0, h, n)

plt.plot(x,t)
plt.title('x x t')
plt.figure()

plt.plot(x,y)
plt.title('órbita')
plt.figure()
