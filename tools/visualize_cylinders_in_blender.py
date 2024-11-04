import bpy
import math


def point_at(x, y, z):
    bpy.ops.mesh.primitive_ico_sphere_add(size=0.5, location=(x, y, z))


def cylinder_between(x1, y1, z1, x2, y2, z2, r):
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1
    dist = math.sqrt(dx**2 + dy**2 + dz**2)

    bpy.ops.mesh.primitive_cylinder_add(radius=r, depth=dist, location=(0, 0, dist / 2))

    for obj in bpy.context.selected_objects:
        obj.name = "Origin" + obj.name

    bpy.ops.mesh.primitive_cylinder_add(
        radius=r, depth=dist, location=(dx / 2 + x1, dy / 2 + y1, dz / 2 + z1)
    )

    phi = math.atan2(dy, dx)
    theta = math.acos(dz / dist)

    bpy.context.object.rotation_euler[1] = theta
    bpy.context.object.rotation_euler[2] = phi


# cylinder_between(-24.927, 1647.6, 80.4342, -24.5942, 1648.11, 80.3892, 0.035)
# cylinder_between(-25.5692, 1647.49, 79.4002, -26.5383, 1652.38, 79.3819, 2.955)
# point_at(0.849938, -0.891369, -0.0307044)
