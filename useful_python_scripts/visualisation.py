import pyvista as pv
import pandas as pd
import numpy as np
from simulationgraphs import get_random_element
import matplotlib.colors as mcolors
from tqdm import tqdm

def load_data(file_path):
    """Load tube data from a text file."""
    columns = ["id_ax", "id_sph", "id_branch", "Type", "X", "Y", "Z", "Rin", "Rout", "P"]
    df = pd.read_csv(file_path, delim_whitespace=True, names=columns, comment='#')
    # delete first row
    df = df.iloc[1:] 
    return df


def load_spheres_to_dataframe(file_path):
    """
    Reads the spheres data from a text file and filters rows where x < 50 and y < 50.

    Args:
        file_path (str): Path to the spheres file.

    Returns:
        pd.DataFrame: Filtered DataFrame containing only relevant spheres.
    """
    # Define column names based on the provided data structure
    col_names = ["id_ax", "id_sph", "id_branch", "Type", "X", "Y", "Z", "Rin", "Rout", "P"]
    
    # Load data into a DataFrame
    df = pd.read_csv(file_path, delim_whitespace=True, names=col_names, skiprows=1)

    # Apply filtering (only keep spheres within x < 50 and y < 50)
    df_filtered = df[(df["X"] < 50) & (df["Y"] < 50)].copy()

    return df_filtered


def plot_points(file_path):
    """
    Reads a file, extracts X, Y, Z, and Rout values, and plots them as a point cloud.
    
    Parameters:
    - file_path (str): Path to the input file.

    Returns:
    - None (Displays a 3D scatter plot)
    """

    # Load data into a Pandas DataFrame
    df = load_spheres_to_dataframe(file_path)
    colors = mcolors.CSS4_COLORS
    df["color"] = df["id_ax"].apply(lambda x: get_random_element(colors, seed = x)[1])
    df["color_rgb"] = df["color"].apply(lambda c: mcolors.to_rgb(c)) 

    # Extract X, Y, Z coordinates and radius (Rout)
    x_values = df["X"].values
    y_values = df["Y"].values
    z_values = df["Z"].values
    R_values = df["Rout"].values  # Sphere radii
    cols = df["color_rgb"].values

    # Create a MultiBlock container (to store multiple sphere meshes)
    spheres = pv.MultiBlock()

    # Generate a sphere for each point
    for (x, y, z, r,color) in zip(x_values, y_values, z_values, R_values, cols):
        sphere = pv.Sphere(radius=r, center=(x, y, z),theta_resolution = 1 ,phi_resolution = 1)  # Create sphere at (x, y, z) with radius r
        spheres.append(sphere)

    # Create a PyVista plotter object
    plotter = pv.Plotter()
    
    # Add all spheres to the plot
    plotter.add_mesh(
        spheres,
        opacity=1,
        show_edges=False,
        color="blue"
    )

    # Set camera position for better visualization
    plotter.view_isometric()

    # Show the plot
    plotter.show()


if __name__ == "__main__":
    file_path = "/home/localadmin/Documents/CATERPillar/c2/voxel6.swc"
    plot_points(file_path)