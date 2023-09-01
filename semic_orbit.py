import numpy as np
import bpy

bpy.ops.mesh.primitive_uv_sphere_add(radius=0.5, location=(0.0, 0.0, 0.0))

def f_(t, x, y, vx, vy):
    res = vx
    return res

def g_(t, x, y, vx, vy):
    res = vy
    return res

def i_(t, x, y, vx, vy):
    res = (-G*M*x)/(((x ** 2) + (y ** 2)) ** (3 / 2))
    return res

def j_(t, x, y, vx, vy):
    res = (-G*M*y)/(((x ** 2) + (y ** 2)) ** (3 / 2))
    return res

def RK4(f, g, i, j, x0, y0, vx0, vy0, t0, h, n):
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
    print(x[0], y[0], vx[0], vy[0], t[0])
    for k in range(n-1):
        kx1 = h * f(t[k], x[k], y[k], vx[k], vy[k])
        ky1 = h * g(t[k], x[k], y[k], vx[k], vy[k])
        kvx1 = h * i(t[k], x[k], y[k], vx[k], vy[k])
        kvy1 = h * j(t[k], x[k], y[k], vx[k], vy[k])
        
        kx2 = h * f(t[k] + h / 2, x[k] + kx1 / 2, y[k] + ky1 / 2, vx[k] + kvx1 / 2, vy[k] + kvy1 / 2)
        ky2 = h * g(t[k] + h / 2, x[k] + kx1 / 2, y[k] + ky1 / 2, vx[k] + kvx1 / 2, vy[k] + kvy1 / 2)
        kvx2 = h * i(t[k] + h / 2, x[k] + kx1 / 2, y[k] + ky1 / 2, vx[k] + kvx1 / 2, vy[k] + kvy1 / 2)
        kvy2 = h * j(t[k] + h / 2, x[k] + kx1 / 2, y[k] + ky1 / 2, vx[k] + kvx1 / 2, vy[k] + kvy1 / 2)
        
        kx3 = h * f(t[k] + h / 2, x[k] + kx2 / 2, y[k] + ky2 / 2, vx[k] + kvx2 / 2, vy[k] + kvy2 / 2)
        ky3 = h * g(t[k] + h / 2, x[k] + kx2 / 2, y[k] + ky2 / 2, vx[k] + kvx2 / 2, vy[k] + kvy2 / 2)
        kvx3 = h * i(t[k] + h / 2, x[k] + kx2 / 2, y[k] + ky2 / 2, vx[k] + kvx2 / 2, vy[k] + kvy2 / 2)
        kvy3 = h * j(t[k] + h / 2, x[k] + kx2 / 2, y[k] + ky2 / 2, vx[k] + kvx2 / 2, vy[k] + kvy2 / 2)
        
        kx4  = h * f(t[k] + h, x[k] + kx3, y[k] + ky3, vx[k] + kvx3, vy[k] + kvy3)
        ky4  = h * g(t[k] + h, x[k] + kx3, y[k] + ky3, vx[k] + kvx3, vy[k] + kvy3)
        kvx4 = h * i(t[k] + h, x[k] + kx3, y[k] + ky3, vx[k] + kvx3, vy[k] + kvy3)
        kvy4 = h * j(t[k] + h, x[k] + kx3, y[k] + ky3, vx[k] + kvx3, vy[k] + kvy3)

        x[k+1] = x[k] + (1/6)*(kx1 + 2*kx2 + 2*kx3 + kx4)
        y[k+1] = y[k] + (1/6)*(ky1 + 2*ky2 + 2*ky3 + ky4)
        vx[k+1] = vx[k] + (1/6) * (kvx1 + 2 * kvx2 + 2 * kvx3 + kvx4)
        vy[k+1] = vy[k] + (1/6) * (kvy1 + 2 * kvy2 + 2 * kvy3 + kvy4)
        t[k+1] = t[k] + h

        obj = bpy.context.object
        obj.location[0] = x[k+1]
        obj.keyframe_insert(data_path="location", frame=k, index=0)
        obj.location[2] = y[k+1]
        obj.keyframe_insert(data_path="location", frame=k, index=1)
        
    return x, y, vx, vy, t


x0_ = 4*10**(12)
y0_ = 0
vx0_ = 0
vy0_ = 500
t0_ = 0
h_ = 5000
n_ = 200000
G = 6.67*(10**(-11))
M = 1.99*(10**30)


x, y, vx, vy, t = RK4(f_, g_, i_, j_, x0_, y0_, vx0_, vy0_, t0_, h_, n_)
