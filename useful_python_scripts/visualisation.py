import numpy as np
import pandas as pd
import pyvista as pv
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor

def load_sphere_data(file_path):
    """
    Load the data from a text file and return a pandas DataFrame.
    Ensure the numerical columns are correctly typed.
    """
    # Load the data from the text file into a DataFrame
    column_names = ["id_ax", "id_sph", "id_branch", "Type", "X", "Y", "Z", "Rin", "Rout", "P"]

    # Skip the header if it exists (if there is no header, remove the skiprows argument)
    data = pd.read_csv(file_path, delim_whitespace=True, names=column_names, skiprows=1, dtype={
        "id_ax": int,
        "id_sph": int,
        "id_branch": int,
        "Type": str,
        "X": float,
        "Y": float,
        "Z": float,
        "Rin": float,
        "Rout": float,
        "P": float
    })

    return data

def assign_colors_by_ax_id(data):
    unique_ax_ids = data['id_ax'].unique()
    cmap = plt.get_cmap('tab20')
    colors = cmap(np.linspace(0, 1, len(unique_ax_ids)))
    ax_id_to_color = {ax_id: colors[i] for i, ax_id in enumerate(unique_ax_ids)}
    data['color'] = data['id_ax'].apply(lambda ax_id: ax_id_to_color[ax_id])
    return data


def plot_spheres(data):
    """
    Plot the spheres using PyVista, where each sphere is placed at (X, Y, Z) with radius Rout
    and colored by 'id_ax'.
    """
    # Initialize the PyVista plotter
    plotter = pv.Plotter()

    # Add progress bar using tqdm to track progress
    for _, row in tqdm(data.iterrows(), total=len(data), desc="Plotting Spheres", unit="sphere"):
        
        center = np.array([row['X'], row['Y'], row['Z']])
        radius = float(row['Rout'])  # Ensure radius is a float
        color = row['color']

        # Add a sphere to the plotter
        plotter.add_mesh(pv.Sphere(radius=radius, center=center), color=mcolors.to_hex(color), opacity=1)

    # Show the plot
    plotter.show()

def create_sphere(row):
    """
    Create a PyVista sphere mesh for the given row.
    """
    center = np.array([row['X'], row['Y'], row['Z']])
    radius = float(row['Rout'])  # Ensure radius is a float
    color = 'white'  # Spheres will be white
    opacity = 0.5    # Transparent white

    # Return the sphere mesh and color/opacity
    return pv.Sphere(radius=radius, center=center), mcolors.to_hex(color), opacity

def plot_spheres_black_white(data):
    """
    Plot the spheres using PyVista, where each sphere is placed at (X, Y, Z) with radius Rout
    and spheres appear as transparent white with a black background.
    """
    # Initialize the PyVista plotter with a black background
    plotter = pv.Plotter()
    plotter.set_background('black')  # Set the background to black

    # Prepare spheres in parallel
    spheres = []
    with ThreadPoolExecutor() as executor:
        results = list(tqdm(executor.map(create_sphere, [row for _, row in data.iterrows()]), 
                            total=len(data), desc="Creating Spheres", unit="sphere"))
        spheres.extend(results)

    # Add the spheres to the plotter sequentially
    for sphere, color, opacity in tqdm(spheres, desc="Plotting Spheres", unit="sphere"):
        plotter.add_mesh(sphere, color=color, opacity=opacity)

    # Show the plot
    plotter.show()

if __name__ == "__main__":
    # Load the data from the file
    file_path = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/overlapping_factor/factor_2.swc"  # Replace with your file path
    data = load_sphere_data(file_path)

    # Assign colors based on 'id_ax'
    #data = assign_colors_by_ax_id(data)
    data["color"] = ["orange"]*len(data)

    data = data.loc[data["Type"] == "axon"]

    data = data.loc[data["id_ax"] == 1]

    # Plot the spheres
    plot_spheres(data)
