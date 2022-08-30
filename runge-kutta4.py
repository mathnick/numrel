import numpy as np
from matplotlib import pyplot as plt

def f(x,y):
    return x*(a-b*y)
def g(x,y):
    return y*(-c+d*x)

def rk4(f, g, x0, y0, t0, h, n):
    #f-população de coelhos, g-população de lobos
    x = np.zeros(n)
    y = np.zeros(n)
    t = np.zeros(n)
    x[0] = x0
    y[0] = y0
    t[0] = t0

    for k in range(n-1):
        k1 = h*f(x[k],y[k])
        k2 = h*f(x[k]+h/2,y[k] + k1/2)
        k3 = h*f(x[k]+h/2,y[k] + k2/2)
        k4 = h*f(x[k]+h,y[k] + k3)
        kk1 = h * g(x[k], y[k])
        kk2 = h * g(x[k] + h / 2, y[k] + kk1 / 2)
        kk3 = h * g(x[k] + h / 2, y[k] + kk2 / 2)
        kk4 = h * g(x[k] + h, y[k] + kk3)
        x[k+1] = x[k] + (1/6)*(k1 + 2*k2 + 2*k3 +k4)
        y[k+1] = y[k] + (1/6)*(kk1 + 2*kk2 + 2*kk3 + kk4)
        t[k+1] = t[k] + h
        print(x[k+1],y[k+1],t[k+1])
    return x,y,t


x0=1
y0=1
t0=0
h=0.01
n=2000
a=1
b=1
c=1
d=2

x, y, t = rk4(f, g, x0, y0, t0, h, n)

plt.figure()
plt.plot(t,x,)
plt.plot(t,y,)
plt.title('Evolução das popoulções ao longo do tempo')
plt.xlabel('tempo')
plt.ylabel('população')


plt.figure()
plt.plot(x,y)
plt.title('Relação entre as populções de lobos e coelhos')
plt.xlabel('coelhos')
plt.ylabel('lobos')

plt.show()
