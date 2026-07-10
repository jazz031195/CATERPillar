import numpy as np
import pyvista as pv
from simulationgraphs import read_swc_file
def spheres_to_surface_mesh(sphere_centers, radii, resolution=100):
    """
    Transforms a list of overlapping spheres into a single continuous surface mesh.
    """
    # 1. Define the 3D bounding box for our grid
    min_bounds = np.min(sphere_centers - radii[:, None], axis=0)
    max_bounds = np.max(sphere_centers + radii[:, None], axis=0)
    
    # Add a small padding
    padding = np.max(radii)
    min_bounds -= padding
    max_bounds += padding

    # 2. Create a 3D uniform grid (Voxel space)
    grid = pv.ImageData(
        dimensions=(resolution, resolution, resolution),
        spacing=((max_bounds - min_bounds) / (resolution - 1)),
        origin=min_bounds
    )

    # 3. Get the coordinates of all points in the grid
    points = grid.points

    # 4. Calculate the Scalar Field (Signed Distance or inside/outside boolean)
    # We will compute the distance from every grid point to every sphere center.
    # To avoid memory errors with huge arrays, we do this iteratively or use a KDTree.
    # For simplicity, here is a vectorized distance field:
    
    # Initialize the field with a large positive number
    scalar_field = np.full(grid.n_points, np.inf)

    for center, r in zip(sphere_centers, radii):
        # Distance squared from all grid points to this sphere center
        dist_sq = np.sum((points - center)**2, axis=1)
        # We want values <= 0 to be "inside" the sphere
        surface_dist = dist_sq - (r**2)
        # The union of all spheres is the minimum of their individual distance fields
        scalar_field = np.minimum(scalar_field, surface_dist)

    # Add the scalar field to the grid
    grid["distance_field"] = scalar_field

    # 5. Run Marching Cubes (contouring) at the threshold of 0 (the surface)
    # PyVista uses the incredibly fast 'Flying Edges' algorithm under the hood
    mesh = grid.contour([0.0], scalars="distance_field")

    return mesh

def create_mesh(cell_type, cell_id=None, max_cell_id=None):
    df_ = df[df['cell_type'] == cell_type].copy()
    if cell_id is not None:
        df_ = df_[df_['cell_id'] == cell_id]
    if max_cell_id is not None:
        random_ids = np.random.choice(df_['cell_id'].unique(), size=min(max_cell_id, len(df_['cell_id'].unique())), replace=False)
        df_ = df_[df_['cell_id'].isin(random_ids)]

    x = df_['X'].values
    y = df_['Y'].values
    z = df_['Z'].values

    centers = np.column_stack((x, y, z))
    radii = df_['outer_radius'].values

    print(f"Creating mesh for {cell_type} with {len(centers)} spheres.")


    max_number_spheres = 300000  # Adjust based on your system's memory capacity
    if len(centers) > max_number_spheres:
        print(f"Warning: Too many spheres ({len(centers)}). Consider subsampling to {max_number_spheres} for performance.")
        indices = np.random.choice(len(centers), max_number_spheres, replace=False)
        centers = centers[indices]
        radii = radii[indices]

    # Generate the unified mesh
    surface_mesh = spheres_to_surface_mesh(centers, radii, resolution=200)
    return surface_mesh

if __name__ == "__main__":
    # --- Example Usage ---
    file_path = "/home/localadmin/Documents/Santi/Healthy_Voxel.csv"
    
    # Assuming read_swc_file returns a pandas DataFrame
    df = read_swc_file(file_path)

    # 1. Create the meshes and save them to SEPARATE variables
    mesh_vessels = create_mesh(cell_type="glial_cell", cell_id=None, max_cell_id=2)
    mesh_axons = create_mesh(cell_type="axon", cell_id=None, max_cell_id=25)

    # 2. ADDING DATA ELEMENTS TO THE MESH
    # If you want to tag the meshes with an ID array before saving or merging
    # Let's say: 1 for vessels, 2 for axons
    if mesh_vessels.n_points > 0:
        mesh_vessels["structure_type"] = np.ones(mesh_vessels.n_points)
    
    if mesh_axons.n_points > 0:
        mesh_axons["structure_type"] = np.full(mesh_axons.n_points, 2.0)

    # 3. MERGING INTO ONE MESH (Optional)
    # If you need a single mesh object for export, you can merge them:
    # unified_mesh = mesh_vessels.merge(mesh_axons)
    # unified_mesh.save("combined_structures.vtk")

    # --- Plotting ---
    plotter = pv.Plotter()
    
    # Add the blood vessels to the scene (colored red)
    if mesh_vessels.n_points > 0:
        plotter.add_mesh(
            mesh_vessels, 
            color="green", 
            smooth_shading=True, 
            specular=0.5, 
            label="Glial Cells"
        )
        
    # Add the axons to the scene (colored green)
    if mesh_axons.n_points > 0:
        plotter.add_mesh(
            mesh_axons, 
            color="blue", 
            smooth_shading=True, 
            specular=0.5,
            label="Axons"
        )

    # Add a legend so you know which is which
    plotter.add_legend()
    plotter.add_title("Unified Cell Surfaces (Marching Cubes)")
    plotter.show()