import numpy as np
from matplotlib import pyplot as plt

h = 0.1
n = 9
m = 20
G = 6.67*(10**(-11))
a = 1
x_0 = np.zeros((n, 1))
y_0 = np.zeros((n, 1))
z_0 = np.zeros((n, 1))
vx_0 = np.zeros((n, 1))
vy_0 = np.zeros((n, 1))
vz_0 = np.zeros((n, 1))
x_0[0,0] = 1
y_0[0,0] = 0
z_0[0,0] = 0
vx_0[0,0] = 0
vy_0[0,0] = 0
vz_0[0,0] = 0

# pos = np.zeros((n,3))
# vel = np.zeros((n,3))


R = np.zeros((n,m))
Fx = np.zeros((n,m))
Fy = np.zeros((n,m))
Fz = np.zeros((n,m))
M = np.zeros(n)
M[0] = 1

# x = np.zeros(n,m)
# y = np.zeros(n,m)
# z = np.zeros(n,m)
# vx = np.zeros(n,m)
# vy = np.zeros(n,m)
# vz = np.zeros(n,m)
# m = np.zeros(n,m)
# x[0] = x_0[0, 0]
# y[0] = y_0[0, 0]
# z[0] = z_0[0, 0]
# vx[0] = vx_0[0, 0]
# vy[0] = vy_0[0, 0]
# vz[0] = vz_0[0, 0]
# m[0] = 1


# print(x_0[0], y_0[0], z_0[0], m[0])


def dados_in(x_0, y_0, z_0, vx_0, vy_0, vz_, m):
    for k in range (n):
        x_0[k] = a*np.cos(k*(np.pi/180)*45)
        y_0[k] = a*np.sin(k*(np.pi/180)*45)
        z_0[k] = 0
        vx_0[k] = 0
        vy_0[k] = 0
        vz_0[k] = 0
        m[k] = 1
        print(x_0[k], y_0[k], z_0[k], vx_0[k], vy_0[k], vz_0[k], m[k])
    return  x_0 , y_0, z_0, vx_0, vy_0, vz_0, m

x_0, y_0, z_0, vx_0, vy_0, vz_0, m = dados_in(x_0, y_0, z_0, vx_0, vy_0, vz_0, m)

x = np.zeros((n,m))
y = np.zeros((n,m))
z = np.zeros((n,m))
vx = np.zeros((n,m))
vy = np.zeros((n,m))
vz = np.zeros((n,m))
x[n,0] = x_0[n, 0]
y[n,0] = y_0[n, 0]
z[n,0] = z_0[n, 0]
vx[n,0] = vx_0[n, 0]
vy[n,0] = vy_0[n, 0]
vz[n,0] = vz_0[n, 0]

pos = np.zeros((n,3))
vel = np.zeros((n,3))


# def pos(x,y):
#     for k in range (n-1):
#         r[k] = ((x[k]))

def dis(x ,y ,z):
    for i in range (n-1):
        for j in range (n-1):
            if i!= j:
                R[i,j] = ((x[i]-x[j])**2 + (y[i]-y[j])**2 + (z[i] - z[j])**2)**(1/2)
                print(R[i,j])
    return R


R = dis(x, y, z)


def W(x, y, z, h ):
    res = (1 / (h * (np.pi) ** (3 / 2))) * (np.exp(-(R / h)**2))
    return W
#Equação do kernel

def rho(R):
    res = m*W(R,h)
    return rho
#Equação da densidade

def P(rho):
    res = rho**(5/3)
    return P
#Equação da pressão, esta coreeta ?

def gradW(R):
    W = (1 / (h * (np.pi) ** (3 / 2))) * 2 * (np.exp(-(R / h)**2)) * (1/h**2)
    Wx = 2*x*W
    Wy = 2*y*W
    Wz = 2*z*W
    return Wx, Wy, Wz
Wx, Wy, Wz = gradW(R)

def gradP(R):
    res = np.sum(m*(P/rho)*(Wx + Wy + Wz))
    return res

def Fgx(R,m):
    for i in range(n-1):
        for j in range(n-1):
            if i != j:
                Fx[i,j] = ((G*m[i]*m[j]*(x[i] - x[j]))/(R**(3/2)))
        return Fx

Fx = Fgx(R,m)

def Fgy(R,m):
    for i in range(n-1):
        for j in range(n-1):
            if i != j:
                Fy[i,j] = ((G*m[i]*m[j]*(y[i] - y[j]))/(R**(3/2)))
        return Fy

Fy = Fgy(R,m)

def Fgz(R,m):
    for i in range(n-1):
        for j in range(n-1):
            if i != j:
                Fz[i,j] = ((G*m[i]*m[j]*(z[i] - z[j]))/(R**(3/2)))
        return Fz

Fz = Fgz(R,m)

def f(t,u):
    x, y, z, vx, vy, vz = u

    return np.array([vx, vy, vz, Fx, Fy, Fz])

def RK4_(f, t, u, h):
    k1 = h * f(t, u)
    k2 = h * f(t + h/2, u + k1/2)
    k3 = h * f(t + h/2, u + k2/2)
    k4 = h * f(t + h, u + k3)
    return u + (1/6) * (k1 + 2*k2 + 2*k3 + k4)

def RK4(f, t0, tf, u0, h):
    n = int((tf - t0) / h) + 1
    t = np.linspace(t0, tf, n)
    u = np.zeros((n, len(u0)))
    u[0] = u0

    for i in range(n-1):
        u[i+1] = RK4_(f, t[i], u[i], h)
        pos[i] = (u[:, 0], u[:, 1], u[:, 2])
        vel[i] = (u[:, 3], u[:, 4], u[:, 5])

    return t, u, pos, vel

t0, tf = 0, 20
u0 = np.array([x_0, y_0, z_0, vx_0, vy_0, vz_0])
h=0.1
t, u, pos, vel = RK4(f, t0, tf, u0, h)

import matplotlib.pyplot as plt
plt.scatter(pos[:, 0], pos[:, 1])
plt.xlabel('Posição X')
plt.ylabel('Posição Y')
plt.title('Simulação SPH')
plt.show()
