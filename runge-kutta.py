import numpy as np
from matplotlib import pyplot as plt

def f(x,y):
    return x*(a-b*y)
def g(x,y):
    return y*(-c+d*x)

def euler(f, g, x0, y0, t0, h, n):
    #f-população de coelhos, g-população de lobos
    x = np.zeros(n)
    y = np.zeros(n)
    t = np.zeros(n)
    x[0] = x0
    y[0] = y0
    t[0] = t0

    for k in range(n-1):
        x[k+1] = x[k] + h*f(x[k]+h/2,y[k]+(h/2)*f(x[k],y[k]))
        y[k+1] = y[k] + h*g(x[k]+h/2,y[k]+(h/2)*g(x[k],y[k]))
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

x, y, t = euler(f, g, x0, y0, t0, h, n)

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
