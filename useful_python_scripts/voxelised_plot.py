import numpy as np
import pyvista as pv
from simulationgraphs import read_swc_file
import numpy as np
import pyvista as pv
import time

def spheres_to_pure_voxels_fast(sphere_centers, radii, resolution=150):
    """
    Highly optimized voxelization using localized bounding boxes (splatting).
    'resolution' is the number of solid VOXELS along each axis.
    """
    print("Starting voxelization...")
    start_time = time.time()

    # 1. Define the global bounding box
    min_bounds = np.min(sphere_centers - radii[:, None], axis=0)
    max_bounds = np.max(sphere_centers + radii[:, None], axis=0)
    padding = np.max(radii)
    min_bounds -= padding
    max_bounds += padding

    # 2. Calculate voxel spacing
    # Spacing is total length divided by the number of voxels
    spacing = (max_bounds - min_bounds) / resolution

    # 3. Create 1D coordinate arrays for the CELL CENTERS
    # A cell center is half a spacing step inward from the edges
    x_coords = np.linspace(min_bounds[0] + spacing[0]/2, max_bounds[0] - spacing[0]/2, resolution)
    y_coords = np.linspace(min_bounds[1] + spacing[1]/2, max_bounds[1] - spacing[1]/2, resolution)
    z_coords = np.linspace(min_bounds[2] + spacing[2]/2, max_bounds[2] - spacing[2]/2, resolution)

    # 4. Initialize the master 3D boolean mask (False = empty, True = solid)
    inside_mask = np.zeros((resolution, resolution, resolution), dtype=bool)

    # 5. FAST LOOP: Localized Splatting
    for center, r in zip(sphere_centers, radii):
        
        # Physical bounding box of THIS sphere
        s_min = center - r
        s_max = center + r

        # Convert physical bounds to Voxel Grid Indices
        idx_min = np.maximum(0, np.floor((s_min - min_bounds) / spacing)).astype(int)
        idx_max = np.minimum(resolution - 1, np.floor((s_max - min_bounds) / spacing)).astype(int)

        # Extract only the local 1D coordinates for this specific bounding box
        xi = x_coords[idx_min[0] : idx_max[0] + 1].reshape(-1, 1, 1)
        yi = y_coords[idx_min[1] : idx_max[1] + 1].reshape(1, -1, 1)
        zi = z_coords[idx_min[2] : idx_max[2] + 1].reshape(1, 1, -1)

        # Calculate distances ONLY for the voxels inside this tiny sub-grid
        dist_sq = (xi - center[0])**2 + (yi - center[1])**2 + (zi - center[2])**2

        # Update ONLY the tiny slice of the master mask
        inside_mask[idx_min[0] : idx_max[0] + 1,
                    idx_min[1] : idx_max[1] + 1,
                    idx_min[2] : idx_max[2] + 1] |= (dist_sq <= r**2)

    print(f"Math finished in {time.time() - start_time:.2f} seconds.")
    
    # 6. Push the completed mask into PyVista
    print("Building PyVista mesh...")
    
    # THE FIX: PyVista dimensions require Points. Points = Voxels + 1.
    grid = pv.ImageData(
        dimensions=(resolution + 1, resolution + 1, resolution + 1),
        spacing=spacing,
        origin=min_bounds
    )
    
    # Flatten the 3D mask using Fortran order ('F') to match VTK/PyVista's memory layout
    grid.cell_data["is_inside"] = inside_mask.flatten(order="F").astype(int)

    # 7. Extract the solid cubes
    voxel_mesh = grid.threshold(0.5, scalars="is_inside")
    
    print("Done!")
    return voxel_mesh

if __name__ == "__main__":
    # --- Example Usage ---
    file_path = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/incoherent_blood_flow/Voxel.csv"
    df = read_swc_file(file_path)
    cell_type = "blood_vessel"
    cell_id = None

    df = df[df['cell_type'] == cell_type]
    if cell_id is not None:
        df = df[df['cell_id'] == cell_id]
    x = df['X'].values
    y = df['Y'].values
    z = df['Z'].values

    centers = np.column_stack((x, y, z))
    radii = df['outer_radius'].values

    # Generate the blocky voxel mesh (Lower resolution makes blocks more obvious)
    voxels = spheres_to_pure_voxels_fast(centers, radii, resolution=1000)

    # Plot it
    plotter = pv.Plotter()

    # 'show_edges=True' draws the black lines around each individual cubic voxel
    plotter.add_mesh(voxels, color="coral", show_edges=True, ambient=0.2, diffuse=0.8)
    plotter.add_title("Blocky Voxel Representation")
    plotter.show()