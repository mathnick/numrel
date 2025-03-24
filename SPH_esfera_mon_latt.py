import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from mpl_toolkits.mplot3d import Axes3D

h = 5
n = 500
G = 6.67e-11
k = 0.1
u = 2
interval = 150


# Função para gerar dados iniciais
def dados_in(n, r):
    phi = np.random.uniform(0, np.pi, n)
    theta = np.random.uniform(0, 2 * np.pi, n)
    u = np.random.uniform(0, 1, n)

    r_ = r * u ** (1 / 3)
    x_0 = r_ * np.sin(phi) * np.cos(theta)
    y_0 = r_ * np.sin(phi) * np.sin(theta)
    z_0 = r_ * np.cos(phi)

    vx_0 = np.zeros(n)
    vy_0 = np.zeros(n)
    vz_0 = np.zeros(n)
    M = np.ones(n) * (0.005)
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


def monaghan_lattanzio_kernel(r, h, dim=3):
    q = r / h
    sigma = {1: 2 / 3, 2: 10 / (7 * np.pi), 3: 1 / np.pi}[dim]  # Constante de normalização

    W = np.zeros_like(q)
    mask1 = (q >= 0) & (q <= 1)
    mask2 = (q > 1) & (q <= 2)

    W[mask1] = sigma / (h ** dim) * (1 - 1.5 * q[mask1] ** 2 + 0.75 * q[mask1] ** 3)
    W[mask2] = sigma / (h ** dim) * 0.25 * (2 - q[mask2]) ** 3

    return W


def grad_monaghan_lattanzio_kernel(r, h, dim=3):
    q = r / h
    sigma = {1: 2 / 3, 2: 10 / (7 * np.pi), 3: 1 / np.pi}[dim]  # Constante de normalização

    grad_W = np.zeros_like(q)
    mask1 = (q >= 0) & (q <= 1)
    mask2 = (q > 1) & (q <= 2)

    grad_W[mask1] = sigma / (h ** (dim + 1)) * (-3 * q[mask1] + 2.25 * q[mask1] ** 2)
    grad_W[mask2] = sigma / (h ** (dim + 1)) * (-0.75 * (2 - q[mask2]) ** 2)

    return grad_W


def densidade(x, y, z, M, h):
    n = len(x)
    densidades = np.zeros(n)
    distancias = distancia(x, y, z)
    for i in range(n):
        r = np.sqrt(distancias[0][:, i] ** 2 + distancias[1][:, i] ** 2 + distancias[2][:, i] ** 2)
        W = monaghan_lattanzio_kernel(r, h)
        densidades[i] = np.sum(M * W)

    return densidades

print(densidade(x,y,z,M,h))

def forca_pressao(x, y, z, M, h, densidades):
    n = len(x)
    forcas_px = np.zeros(n)
    forcas_py = np.zeros(n)
    forcas_pz = np.zeros(n)
    distancias = distancia(x, y, z)
    pressao = k * (densidades ** ((u + 1) / u))
    for i in range(n):
        for j in range(n):
            if i != j:
                dx, dy, dz = distancias[0][i, j], distancias[1][i, j], distancias[2][i, j]
                r = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
                grad_W = grad_monaghan_lattanzio_kernel(r, h)
                forcas_px[i] += M[j] * pressao[j] * grad_W * dx / r
                forcas_py[i] += M[j] * pressao[j] * grad_W * dy / r
                forcas_pz[i] += M[j] * pressao[j] * grad_W * dz / r
    forcas_px /= densidades
    forcas_py /= densidades
    forcas_pz /= densidades
    return forcas_px, forcas_py, forcas_pz


def forca_gravitacional(x, y, z, M):
    n = len(x)
    forcas_x = np.zeros(n)
    forcas_y = np.zeros(n)
    forcas_z = np.zeros(n)
    F = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if i != j:
                r = np.sqrt((x[i] - x[j]) ** 2 + (y[i] - y[j]) ** 2 + (z[i] - z[j]) ** 2)
                F = G * M[i] * M[j] / r ** 2
                theta = np.arctan2(y[j] - y[i], x[j] - x[i])
                phi = np.arccos((z[j] - z[i]) / r)
                forcas_x[i] += F * np.sin(phi) * np.cos(theta)
                forcas_y[i] += F * np.sin(phi) * np.sin(theta)
                forcas_z[i] += F * np.cos(phi)
    return forcas_x, forcas_y, forcas_z


def derivadas(x, y, z, vx, vy, vz, M, h):
    densidades = densidade(x, y, z, M, h)
    forcas_px, forcas_py, forcas_pz = forca_pressao(x, y, z, M, h, densidades)
    forcas_grav_x, forcas_grav_y, forcas_grav_z = forca_gravitacional(x, y, z, M)

    ax = (forcas_px + forcas_grav_x) / M
    ay = (forcas_py + forcas_grav_y) / M
    az = (forcas_pz + forcas_grav_z) / M

    return vx, vy, vz, ax, ay, az


def atualizar_rk4(x, y, z, vx, vy, vz, M, h, dt):
    k1_vx, k1_vy, k1_vz, k1_ax, k1_ay, k1_az = derivadas(x, y, z, vx, vy, vz, M, h)
    k1_rx, k1_ry, k1_rz = vx, vy, vz

    k2_vx, k2_vy, k2_vz, k2_ax, k2_ay, k2_az = derivadas(x + 0.5 * dt * k1_rx, y + 0.5 * dt * k1_ry,
                                                         z + 0.5 * dt * k1_rz, vx + 0.5 * dt * k1_ax,
                                                         vy + 0.5 * dt * k1_ay, vz + 0.5 * dt * k1_az, M, h)
    k2_rx, k2_ry, k2_rz = vx + 0.5 * dt * k1_ax, vy + 0.5 * dt * k1_ay, vz + 0.5 * dt * k1_az

    k3_vx, k3_vy, k3_vz, k3_ax, k3_ay, k3_az = derivadas(x + 0.5 * dt * k2_rx, y + 0.5 * dt * k2_ry,
                                                         z + 0.5 * dt * k2_rz, vx + 0.5 * dt * k2_ax,
                                                         vy + 0.5 * dt * k2_ay, vz + 0.5 * dt * k2_az, M, h)
    k3_rx, k3_ry, k3_rz = vx + 0.5 * dt * k2_ax, vy + 0.5 * dt * k2_ay, vz + 0.5 * dt * k2_az

    k4_vx, k4_vy, k4_vz, k4_ax, k4_ay, k4_az = derivadas(x + dt * k3_rx, y + dt * k3_ry, z + dt * k3_rz,
                                                         vx + dt * k3_ax, vy + dt * k3_ay, vz + dt * k3_az, M, h)
    k4_rx, k4_ry, k4_rz = vx + dt * k3_ax, vy + dt * k3_ay, vz + dt * k3_az

    x += (dt / 6) * (k1_rx + 2 * k2_rx + 2 * k3_rx + k4_rx)
    y += (dt / 6) * (k1_ry + 2 * k2_ry + 2 * k3_ry + k4_ry)
    z += (dt / 6) * (k1_rz + 2 * k2_rz + 2 * k3_rz + k4_rz)

    vx += (dt / 6) * (k1_ax + 2 * k2_ax + 2 * k3_ax + k4_ax)
    vy += (dt / 6) * (k1_ay + 2 * k2_ay + 2 * k3_ay + k4_ay)
    vz += (dt / 6) * (k1_az + 2 * k2_az + 2 * k3_az + k4_az)

    return x, y, z, vx, vy, vz



dt = 8
num_frames = 300


fig = plt.figure(figsize=(12, 6))
ax1 = fig.add_subplot(121, projection='3d')
ax1.set_xlim(-2, 2)
ax1.set_ylim(-2, 2)
ax1.set_zlim(-2, 2)
particulas, = ax1.plot([], [], [], 'bo', ms=8)

ax2 = fig.add_subplot(122)
ax2.set_xlabel('Raio')
ax2.set_ylabel('Densidade')
ax2.set_xlim(0, 1.1)
ax2.set_ylim(0.0055, 0.0065)
densidade_plot, = ax2.plot([], [], 'r-')

plt.tight_layout()

def init():
    particulas.set_data([], [])
    particulas.set_3d_properties([])
    densidade_plot.set_data([], [])
    return particulas, densidade_plot


def update(frame):
    global x, y, z, vx, vy, vz
    x, y, z, vx, vy, vz = atualizar_rk4(x, y, z, vx, vy, vz, M, h, dt)
    particulas.set_data(x, y)
    particulas.set_3d_properties(z)

    raio = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    densidades = densidade(x, y, z, M, h)
    sorted_indices = np.argsort(raio)
    raio_sorted = raio[sorted_indices]
    densidades_sorted = densidades[sorted_indices]

    densidade_plot.set_data(raio_sorted, densidades_sorted)



    return particulas, densidade_plot


ani = FuncAnimation(fig, update, frames=num_frames, init_func=init, blit=True, interval = interval)

#ani.save('animacao_esfera_monaghan_lattanzio2.gif', writer=PillowWriter(fps=30))

plt.show()
