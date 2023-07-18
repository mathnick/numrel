import bpy
bpy.ops.mesh.primitive_uv_sphere_add(radius=0.25, location=(0.0, 0.0, 10.0))
def z(t):
    return 10  - 0.001*(t**2)
def x(t):
   return u*t
u = 0.2

for t in range (100) :
      obj = bpy.context.object
      obj.location[0] = x(t)
      obj.keyframe_insert(data_path="location", frame=t, index=0)
      obj.location[2] = z(t)
      obj.keyframe_insert(data_path="location", frame=t, index=2)
