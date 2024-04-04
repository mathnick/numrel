import numpy as np
from matplotlib import pyplot as plt

n = 9
a = 1
x = np.zeros(n)
y = np.zeros(n)
#r = np.zeros(n)
m = np.zeros(n)
x[0] = 1
y[0] = 0
m[0] = 1


print(x[0], y[0], m[0])


def dados(x, y, m):
    for k in range (n-1):
         x[k] = a*np.cos(k*(np.pi/180)*45)
         y[k] = a*np.sin(k*(np.pi/180)*45)
         m[k] = 1
         #r[k] = ((x[k])**2 + (y[k])**2)**(1/2)
         print( x[k], y[k], m[k])
    return  x , y, m

x, y, m = dados(x, y, m)

plt.scatter(x, y,)

plt.xlabel('Eixo X')
plt.ylabel('Eixo Y')


plt.title('Gráfico das posições iniciais')

plt.show()

