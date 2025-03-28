import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from mpl_toolkits.mplot3d import Axes3D

h = 5
n = 30
G = 6.67e-11
k = 0.1
u = 2


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
    return np.exp(-(r / h) ** 2) / (h ** 3 * np.pi ** (3 / 2))


def grad_kernel(x, y, z, h):
    w_ = -2 / (h ** 5 * np.pi ** (3 / 2))
    Wx = w_ * x * np.exp(-(x ** 2 + y ** 2 + z ** 2) / h ** 2)
    Wy = w_ * y * np.exp(-(x ** 2 + y ** 2 + z ** 2) / h ** 2)
    Wz = w_ * z * np.exp(-(x ** 2 + y ** 2 + z ** 2) / h ** 2)
    return Wx, Wy, Wz


def densidade(x, y, z, M, h):
    n = len(x)
    densidades = np.zeros(n)
    distancias = distancia(x, y, z)
    for i in range(n):
        r = np.sqrt(distancias[0][:, i] ** 2 + distancias[1][:, i] ** 2 + distancias[2][:, i] ** 2)
        W = kernel(r, h)
        densidades[i] = np.sum(M * W)
    return densidades


def forca_pressao(x, y, z, M, h, densidades):
    n = len(x)
    forcas_px = np.zeros(n)
    forcas_py = np.zeros(n)
    forcas_pz = np.zeros(n)
    distancias = distancia(x, y, z)
    pressao = k * (densidades ** ((u+1)/u))
    for i in range(n):
        for j in range(n):
            if i != j:
                dx, dy, dz = distancias[0][i, j], distancias[1][i, j], distancias[2][i, j]
                Wx, Wy, Wz = grad_kernel(dx, dy, dz, h)
                forcas_px[i] += M[j] * pressao[j] * Wx / densidades[j]
                forcas_py[i] += M[j] * pressao[j] * Wy / densidades[j]
                forcas_pz[i] += M[j] * pressao[j] * Wz / densidades[j]
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


dt = 1
num_frames = 500

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
ax.set_zlim(-2, 2)
particulas, = ax.plot([], [], [], 'bo', ms=8)


def init():
    particulas.set_data([], [])
    particulas.set_3d_properties([])
    return particulas,


def update(frame):
    global x, y, z, vx, vy, vz
    x, y, z, vx, vy, vz = atualizar_rk4(x, y, z, vx, vy, vz, M, h, dt)
    particulas.set_data(x, y)
    particulas.set_3d_properties(z)
    return particulas,


ani = FuncAnimation(fig, update, frames=num_frames, init_func=init, blit=True)

ani.save('animacao_circulo.gif', writer=PillowWriter(fps=30))

plt.show()
