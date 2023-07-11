import bpy
def x(t):
   return v*t
v = 2.2
for t in range (5) :
      bpy.ops.mesh.primitive_uv_sphere_add(radius=0.25, location=(x(t), 0.0, 0.0))
