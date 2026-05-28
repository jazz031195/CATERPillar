import bpy
import bmesh
from mathutils import Vector
import math
import random

# --------------------------- CONFIG ---------------------------
FILE_PATH = bpy.path.abspath("/home/localadmin/Documents/MCDS/Permeable_MCDS/output/SMI_pred/complex_axons/f_0.5.swc")

SKIP_HEADER_LINES = 1      # lines to skip at top of file
TAKE_EVERY_NTH = 2         # 2 means take every other line (to thin data); set 1 to take all
ROW_LIMIT = 50000           # e.g. 10000 to preview, or None for all

USE_ROUT = True            # True: use Rout column; False: use Rin
RADIUS_SCALE = 1.0         # global scale for radii (scene units)
COORD_SCALE = 1        # scale xyz from µm->mm or similar; 1.0 means raw units

ICOSPHERE_SUBDIV = 1       # 0..2 (keep low for speed)
SHADE_SMOOTH = True

MAKE_MATERIAL = True
MATERIAL_NAME = "Spheres_Mat"

REALIZE_INSTANCES = False  # keep False for speed; True if you need a real mesh

COLL_NAME_POINTS = "Spheres_Points_Collection"
COLL_NAME_INST   = "Spheres_Instance_Collection"
GN_GROUP_NAME    = "GN_InstanceSpheresGroup"
POINTS_OBJ_NAME  = "SpheresPoints"
INSTANCE_OBJ_NAME= "InstanceSphere"
POINTS_MESH_NAME = "SpheresPointsMesh"

# Clean scene (optional). Set to True to wipe current scene.
CLEAN_SCENE = True
# --------------------------------------------------------------


# ----------------------- UTILITIES ----------------------------
def ensure_collection(name):
    coll = bpy.data.collections.get(name)
    if not coll:
        coll = bpy.data.collections.new(name)
        bpy.context.scene.collection.children.link(coll)
    return coll

def clear_scene():
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete(use_global=False)
    for block in (bpy.data.meshes, bpy.data.objects, bpy.data.materials, bpy.data.node_groups):
        for datablock in list(block):
            try:
                block.remove(datablock)
            except:
                pass

def set_obj_collection(obj, coll):
    # unlink from all, link to desired
    for c in obj.users_collection:
        c.objects.unlink(obj)
    coll.objects.link(obj)

def set_object_info_target(n_obj, obj):
    """
    Blender ≤3.6: set via n_obj.object
    Blender 4.x : set via input socket "Object"
    """
    try:
        n_obj.object = obj
    except AttributeError:
        # Blender 4.x
        try:
            n_obj.inputs["Object"].default_value = obj
        except (KeyError, AttributeError):
            # fallback to first input
            n_obj.inputs[0].default_value = obj

def set_object_info_as_instance(n_obj, flag=True):
    """
    Blender ≤3.6: property; Blender 4.x: input socket "As Instance"
    """
    try:
        # 3.x property path varies; try inputs first (4.x)
        if "As Instance" in {s.name for s in n_obj.inputs}:
            n_obj.inputs["As Instance"].default_value = flag
            return
    except Exception:
        pass
    # try attributes
    for attr in ("as_instance", "AsInstance"):
        if hasattr(n_obj, attr):
            setattr(n_obj, attr, flag)
            return
    # last resort: ignore (node may default to instancing)

def add_group_io_sockets(node_group):
    """
    Add Geometry input/output to the node group in a way that works on 3.6 and 4.x
    """
    # Blender 4.x uses interface
    if hasattr(node_group, "interface"):
        iface = node_group.interface
        names_in = [s.name for s in iface.items_tree if s.in_out == 'INPUT']
        names_out= [s.name for s in iface.items_tree if s.in_out == 'OUTPUT']
        if "Geometry" not in names_in:
            iface.new_socket(name="Geometry", in_out='INPUT', socket_type='NodeSocketGeometry')
        if "Geometry" not in names_out:
            iface.new_socket(name="Geometry", in_out='OUTPUT', socket_type='NodeSocketGeometry')
    else:
        # Blender 3.6
        if "Geometry" not in [s.name for s in node_group.inputs]:
            node_group.inputs.new("NodeSocketGeometry", "Geometry")
        if "Geometry" not in [s.name for s in node_group.outputs]:
            node_group.outputs.new("NodeSocketGeometry", "Geometry")

def try_new_node(nodes, bl_idname):
    try:
        return nodes.new(bl_idname)
    except RuntimeError:
        return None

# ---------------------- LOAD THE SWC --------------------------
def load_swc_points(path, skip_header=1, nth=1, limit=None, use_rout=True,
                    coord_scale=1.0, radius_scale=1.0):
    """
    Returns:
        points: list[(x,y,z)]
        radii : list[float]
    """
    points, radii = [], []
    with open(path, 'r') as f:
        for i, line in enumerate(f):
            if not line.strip():
                continue
            if i < skip_header:
                continue
            if nth > 1 and ((i - skip_header) % nth != 0):
                continue
            parts = line.split()
            if len(parts) < 9:
                continue
            # id_ax, id_sph, id_branch, Type, x, y, z, Rin, Rout, P
            x = float(parts[4]) * coord_scale
            y = float(parts[5]) * coord_scale
            z = float(parts[6]) * coord_scale
            Rin  = float(parts[7]) * radius_scale
            Rout = float(parts[8]) * radius_scale
            r = Rout if use_rout else Rin
            points.append((x, y, z))
            radii.append(max(r, 0.0))  # clamp negatives
            if limit is not None and len(points) >= limit:
                break
    return points, radii

# ---------------------- BUILD POINTS MESH ---------------------
def build_points_object(points, radii, obj_name=POINTS_OBJ_NAME, mesh_name=POINTS_MESH_NAME):
    mesh = bpy.data.meshes.new(mesh_name)
    obj  = bpy.data.objects.new(obj_name, mesh)
    bpy.context.scene.collection.objects.link(obj)

    # Build vertices via bmesh
    bm = bmesh.new()
    for p in points:
        bm.verts.new(p)
    bm.verts.ensure_lookup_table()
    bm.to_mesh(mesh)
    bm.free()

    # Add a float attribute "radius" on points domain
    attr = mesh.attributes.get("radius")
    if not attr:
        attr = mesh.attributes.new(name="radius", type='FLOAT', domain='POINT')
    # Fill attribute
    for i, r in enumerate(radii):
        attr.data[i].value = float(r)

    return obj

# ---------------------- INSTANCE SPHERE -----------------------
def build_instance_sphere(obj_name=INSTANCE_OBJ_NAME, subdivisions=1, radius=1.0, shade_smooth=True):
    bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=subdivisions, radius=radius, location=(0,0,0))
    inst = bpy.context.object
    inst.name = obj_name
    if shade_smooth:
        bpy.ops.object.shade_smooth()
    return inst

def attach_material(obj, mat_name=MATERIAL_NAME):
    mat = bpy.data.materials.get(mat_name)
    if mat is None:
        mat = bpy.data.materials.new(mat_name)
        mat.use_nodes = True
        nt = mat.node_tree
        # Simple principled color setup (optional)
        principled = nt.nodes.get("Principled BSDF")
        if principled:
            principled.inputs["Roughness"].default_value = 0.6
            principled.inputs["Metallic"].default_value = 0.0
    if obj.data and (mat.name not in [m.name for m in obj.data.materials]):
        obj.data.materials.append(mat)
    return mat

# --------------------- GEOMETRY NODES -------------------------
def build_gn_instancer(points_obj, instance_obj, group_name=GN_GROUP_NAME,
                       realize_instances=False):
    """
    Create a GN modifier on points_obj that:
        - reads attribute 'radius'
        - instances instance_obj on points
        - scales instances by 'radius'
        - (optionally) realizes instances
    Works on Blender 3.6 and 4.x.
    """
    gn = points_obj.modifiers.new(name="GN_InstanceSpheres", type='NODES')
    node_group = bpy.data.node_groups.new(group_name, 'GeometryNodeTree')
    gn.node_group = node_group

    nodes = node_group.nodes
    links = node_group.links
    nodes.clear()

    # I/O sockets
    add_group_io_sockets(node_group)

    n_in  = nodes.new("NodeGroupInput");  n_in.location  = (-1000, 0)
    n_out = nodes.new("NodeGroupOutput"); n_out.location = ( 1000, 0)

    # Try Named Attribute (preferred). If not available, fall back to Capture Attribute.
    n_named = try_new_node(nodes, "GeometryNodeInputNamedAttribute")

    if n_named is not None:
        # Named Attribute path
        n_named.location = (-800, -200)
        n_named.data_type = 'FLOAT'
        n_named.inputs["Name"].default_value = "radius"
        attr_socket = n_named.outputs.get("Attribute")

        # Instance on Points
        n_iop = nodes.new("GeometryNodeInstanceOnPoints"); n_iop.location = (0, 0)

        # Object Info
        n_obj = nodes.new("GeometryNodeObjectInfo"); n_obj.location = (-500, 200)
        set_object_info_target(n_obj, instance_obj)
        set_object_info_as_instance(n_obj, True)

        # Links
        links.new(n_in.outputs["Geometry"], n_iop.inputs["Points"])
        links.new(n_obj.outputs["Geometry"], n_iop.inputs["Instance"])
        if attr_socket is not None:
            links.new(attr_socket, n_iop.inputs["Scale"])

        last_geo = n_iop.outputs["Instances"]

    else:
        # Capture Attribute fallback
        n_capture = nodes.new("GeometryNodeCaptureAttribute"); n_capture.location = (-600, 0)
        n_capture.data_type = 'FLOAT'
        # feed geometry from input
        links.new(n_in.outputs["Geometry"], n_capture.inputs["Geometry"])

        # A Named Attribute node just to read 'radius' as a field (may still exist here)
        n_read = try_new_node(nodes, "GeometryNodeInputNamedAttribute")
        if n_read is not None:
            n_read.location = (-900, -200)
            n_read.data_type = 'FLOAT'
            n_read.inputs["Name"].default_value = "radius"
            links.new(n_read.outputs["Attribute"], n_capture.inputs["Value"])
        else:
            # As a last resort, set default value (no scaling)
            pass

        # Instance on Points
        n_iop = nodes.new("GeometryNodeInstanceOnPoints"); n_iop.location = (200, 0)

        # Object Info
        n_obj = nodes.new("GeometryNodeObjectInfo"); n_obj.location = (-200, 200)
        set_object_info_target(n_obj, instance_obj)
        set_object_info_as_instance(n_obj, True)

        links.new(n_capture.outputs["Geometry"], n_iop.inputs["Points"])
        links.new(n_obj.outputs["Geometry"], n_iop.inputs["Instance"])
        # Capture Attribute outputs the captured value as a field
        try:
            links.new(n_capture.outputs["Attribute"], n_iop.inputs["Scale"])
        except KeyError:
            # Some builds label the socket differently; ignore scaling if missing
            pass

        last_geo = n_iop.outputs["Instances"]

    if realize_instances:
        n_real = nodes.new("GeometryNodeRealizeInstances"); n_real.location = (600, 0)
        links.new(last_geo, n_real.inputs["Geometry"])
        links.new(n_real.outputs["Geometry"], n_out.inputs["Geometry"])
    else:
        links.new(last_geo, n_out.inputs["Geometry"])

    return gn

# -------------------------- MAIN ------------------------------
def main():
    if CLEAN_SCENE:
        clear_scene()

    # Collections
    coll_points = ensure_collection(COLL_NAME_POINTS)
    coll_inst   = ensure_collection(COLL_NAME_INST)

    # Load data
    points, radii = load_swc_points(
        FILE_PATH,
        skip_header=SKIP_HEADER_LINES,
        nth=TAKE_EVERY_NTH,
        limit=ROW_LIMIT,
        use_rout=USE_ROUT,
        coord_scale=COORD_SCALE,
        radius_scale=RADIUS_SCALE
    )
    print(f"Loaded {len(points)} spheres from {FILE_PATH}")

    if len(points) == 0:
        raise RuntimeError("No points loaded. Check file path and parsing options.")

    # Build points object with 'radius' attribute
    obj_points = build_points_object(points, radii)
    set_obj_collection(obj_points, coll_points)

    # Build instance sphere (low poly)
    inst = build_instance_sphere(subdivisions=ICOSPHERE_SUBDIV, radius=1.0, shade_smooth=SHADE_SMOOTH)
    set_obj_collection(inst, coll_inst)

    if MAKE_MATERIAL:
        attach_material(inst, MATERIAL_NAME)

    # Add GN instancer modifier to points object
    build_gn_instancer(obj_points, inst, group_name=GN_GROUP_NAME, realize_instances=REALIZE_INSTANCES)

    # Focus view on points
    def focus_selected_safe():
        # Only try if we have a screen and at least one 3D View
        win = bpy.context.window_manager.windows[0] if bpy.context.window_manager.windows else None
        if not win or not bpy.context.screen:
            return  # likely headless or no UI; skip

        area = next((a for a in bpy.context.screen.areas if a.type == 'VIEW_3D'), None)
        if not area:
            return  # no 3D View open; skip

        region = next((r for r in area.regions if r.type == 'WINDOW'), None)
        if not region:
            return

        with bpy.context.temp_override(window=win, area=area, region=region):
            try:
                bpy.ops.view3d.view_selected(use_all_regions=False)
            except Exception:
                pass

    # … after creating/ selecting your points object:
    bpy.context.view_layer.objects.active = obj_points
    for o in bpy.context.selected_objects:
        o.select_set(False)
    obj_points.select_set(True)

    focus_selected_safe()  # call the safe helper


    print("Geometry Nodes instancing set up. If you need a real mesh, set REALIZE_INSTANCES=True (slower).")

# Run
if __name__ == "__main__":
    main()
