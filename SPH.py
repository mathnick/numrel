import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

h = 0.1
n = 15
G = 6.67e-11

def dados_in(n, r):
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    x_0 = r * np.cos(theta)
    y_0 = r * np.sin(theta)
    z_0 = np.zeros(n)
    vx_0 = np.zeros(n)
    vy_0 = np.zeros(n)
    vz_0 = np.zeros(n)
    M = np.ones(n) * 10 
    return x_0, y_0, z_0, vx_0, vy_0, vz_0, M

x_0, y_0, z_0, vx_0, vy_0, vz_0, M = dados_in(n, r=1.0)

x = np.zeros((n,))
y = np.zeros((n,))
z = np.zeros((n,))
vx = np.zeros((n,))
vy = np.zeros((n,))
vz = np.zeros((n,))
x[...] = x_0
y[...] = y_0
z[...] = z_0

def distancia(x, y, z):
    dx = x[:, np.newaxis] - x[np.newaxis, :]
    dy = y[:, np.newaxis] - y[np.newaxis, :]
    dz = z[:, np.newaxis] - z[np.newaxis, :]
    return dx, dy, dz

def kernel(r, h):
    return (1 / (h * np.sqrt(np.pi))) * np.exp(-(r / h)**2)

def densidade(x, y, z, M, h):
    n = len(x)
    densidades = np.zeros(n)
    distancias = distancia(x, y, z)  # Calcula as distâncias entre as partículas
    for i in range(n):
        W = kernel(distancias[0][:, i], h)  # Correção aqui
        densidades[i] = np.sum(M * W)
    return densidades

def forca_pressao(x, y, z, M, h):
    n = len(x)
    forcas_px = np.zeros(n)
    forcas_py = np.zeros(n)
    forcas_pz = np.zeros(n)
    densidades = densidade(x, y, z, M, h)  
    distancias = distancia(x, y, z) 
    for i in range(n):
        W = kernel(distancias[0][:, i], h) 
        for j in range(n):
            if i != j:
                forcas_px[i] -= (M[j] / densidades[j]) * (x[i] - x[j]) * W[j]
                forcas_py[i] -= (M[j] / densidades[j]) * (y[i] - y[j]) * W[j]
                forcas_pz[i] -= (M[j] / densidades[j]) * (z[i] - z[j]) * W[j]
    return forcas_px, forcas_py, forcas_pz

def forca_gravitacional(x, y, z, M):
    n = len(x)
    forcas_x = np.zeros(n)
    forcas_y = np.zeros(n)
    forcas_z = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if i != j:
                r = np.sqrt((x[i] - x[j])**2 + (y[i] - y[j])**2 + (z[i] - z[j])**2)
                F = G * M[i] * M[j] / r**2
                theta = np.arctan2(y[j] - y[i], x[j] - x[i])
                phi = np.arccos((z[j] - z[i]) / r)
                forcas_x[i] += F * np.sin(phi) * np.cos(theta)
                forcas_y[i] += F * np.sin(phi) * np.sin(theta)
                forcas_z[i] += F * np.cos(phi)
    return forcas_x, forcas_y, forcas_z

def atualizar(x, y, z, vx, vy, vz, M, h, dt):
    forcas_px, forcas_py, forcas_pz = forca_pressao(x, y, z, M, h)
    forcas_grav_x, forcas_grav_y, forcas_grav_z = forca_gravitacional(x, y, z, M)
    vx += (dt / M) * (forcas_px + forcas_grav_x)
    vy += (dt / M) * (forcas_py + forcas_grav_y)
    vz += (dt / M) * (forcas_pz + forcas_grav_z)
    x += dt * vx
    y += dt * vy
    z += dt * vz
    return x, y, z, vx, vy, vz

dt = 0.1  # Intervalo de tempo
num_frames = 100 

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
ax.set_zlim(-2, 2)
particles, = ax.plot([], [], [], 'bo', ms=8)

def init():
    particles.set_data([], [])
    particles.set_3d_properties([])
    return particles,

def update(frame):
    global x, y, z, vx, vy, vz
    x, y, z, vx, vy, vz = atualizar(x, y, z, vx, vy, vz, M, h, dt)
    particles.set_data(x, y)
    particles.set_3d_properties(z)
    return particles,

ani = FuncAnimation(fig, update, frames=num_frames, init_func=init, blit=True)

plt.show()
