import bpy
import numpy as np

bpy.ops.mesh.primitive_uv_sphere_add(radius=0.5, location=(0.0, 0.0, 0.0))

def z(t, vini, angulo, gravidade):
    return vini * np.sin(angulo) * t - (gravidade/2)*(t**2)

def x(t, vini, angulo):
   return vini * np.cos(angulo)*t

g = 9.8
theta = (np.pi / 180) * 80
v0 = 90
tR = 2*v0*np.sin(theta)/g
h = 0.01
dt = h * tR

for t in range (100) :
      obj = bpy.context.object
      obj.location[0] = x(t*dt, v0, theta)
      obj.keyframe_insert(data_path="location", frame=t, index=0)
      obj.location[2] = z(t*dt, v0, theta, g)
      obj.keyframe_insert(data_path="location", frame=t, index=2)

# import bpy
# bpy.ops.mesh.primitive_uv_sphere_add(radius=0.25, location=(0.0, 0.0, 10.0))
# def z(t):
#     return 10  - 0.001*(t**2)
# def x(t):
#    return u*t
# u = 0.2

# for t in range (100) :
#       obj = bpy.context.object
#       obj.location[0] = x(t)
#       obj.keyframe_insert(data_path="location", frame=t, index=0)
#       obj.location[2] = z(t)
#       obj.keyframe_insert(data_path="location", frame=t, index=2)
