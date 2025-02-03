import bpy

# File path to the text file containing sphere data
file_path = bpy.path.abspath("/home/localadmin/Documents/CATERPillar/arthurs_analysis/voxel_1.swc")


# Clear all objects in the current scene (optional, to start fresh)
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

# Dictionary to store objects by axon ID
axons = {}

# Function to create a sphere at a given location with a specified radius
def create_sphere(location, radius):
    bpy.ops.mesh.primitive_uv_sphere_add(radius=radius, location=location, segments=4, ring_count=2)
    return bpy.context.object

# Open and read the file
with open(file_path, 'r') as file:
    for e, line in enumerate(file):
        if (e == 0):
            continue
        if ((e-1)%2 != 0):
            continue
        if (e > 1000):
            break
        # Parse each line into values
        id_ax, id_sph, id_branch, Type, x, y, z, Rin, Rout, P = line.split()

        # Convert necessary values to floats
        id_ax = int(id_ax)
        x, y, z, Rout = float(x), float(y), float(z), float(Rout)
        
        # Create a sphere at the given location with the specified outer radius (Rout)
        sphere = create_sphere(location=(x, y, z), radius=Rout)

        # Add the sphere to the appropriate axon group
        if id_ax not in axons:
            axons[id_ax] = []
        axons[id_ax].append(sphere)

# Combine all spheres belonging to the same axon
for axon_id, spheres in axons.items():
    bpy.ops.object.select_all(action='DESELECT')
    # Select all spheres for this axon
    for obj in spheres:
        obj.select_set(True)
    # Set one active object (needed for join)
    bpy.context.view_layer.objects.active = spheres[0]
    # Join them into a single object
    bpy.ops.object.join()
    # Optionally rename the new object to identify the axon
    bpy.context.object.name = f"Axon_{axon_id}"
