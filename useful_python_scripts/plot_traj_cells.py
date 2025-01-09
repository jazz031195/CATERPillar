import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
import plotly.graph_objects as go
from random import seed, randint
from simulationgraphs import read_swc_file, get_random_element
import matplotlib.colors as mcolors
import plotly.colors as colors

def read_bin_file(traj_file_path, swc_file_path, glial_only = False):

    N = 20  # Size scaling factor
    colors_ = mcolors.CSS4_COLORS  # Use CSS colors for better compatibility

    # Read and preprocess data
    df = read_swc_file(swc_file_path)

    # Color assignment with a fallback in case of errors
    df["color"] = df["ax_id"].apply(lambda x: get_random_element(colors_, seed=x)[1])

    df["color"] = list(map(lambda x, y, z: "white" if x == "axon" else y, list(df["type"]), list(df["color"]), list(df["ax_id"])))

    # Filter for axons in a specific region
    df_ = df.loc[df["type"] == "axon"]
    df_ = df_.groupby("ax_id").first().reset_index()
    ax_ids = df_["ax_id"].unique()

    # Filter axons by selected ax_ids
    df_axons = df.loc[df["type"] == "axon"]
    df_axons = df_axons.loc[df_axons["ax_id"].isin(ax_ids)]

    # Glial cells selection
    df_glial = df.loc[df["type"] != "axon"]
    # df_glial = df_glial.loc[df_glial["ax_id"]==1]
    df = pd.concat([df_axons, df_glial])

    # Calculate distance from (0,0,0) for sorting purposes
    df["distance_to_point"] = np.linalg.norm(df[["x", "y", "z"]], axis=1)
    df = df.sort_values(by="distance_to_point")

    if glial_only:
        df = df.loc[df["type"] != "axon"]

    # Ensure marker size is valid
    df["Rout"] = df["Rout"].fillna(0)  # Replace NaN or missing sizes with a default value

    # Create a 3D scatter plot
    scatter_cell = go.Scatter3d(
        x=df["x"],
        y=df["y"],
        z=df["z"],
        mode="markers",
        name="Axons",
        marker=dict(
            sizemode="diameter",
            size=df["Rout"] * N,  # Scale the marker size
            color=df["color"],  # Set the color to white
            opacity=0.5,  # Set transparency to 50%
            line=dict(
                color="rgba(0, 0, 0, 0)",
                width=0
            )
        )
    )



    traj_part = np.fromfile(traj_file_path, dtype="float32")
    chunks = 1
    nbr_steps=115500
    valid_traj_part = traj_part[np.isfinite(traj_part)]  # Ensure all entries are valid floats
    
    scatters = []
    for i in range(chunks):
        start = (nbr_steps+1) * i
        end = nbr_steps * (i + 1)-1
        
        chunk = valid_traj_part[3*start:3*end]
        xs, ys, zs = chunk[0::3], chunk[1::3], chunk[2::3]

        print("done")
        # Example scatter plot
        colours = colors.qualitative.Plotly[:10]
        c = colours[0]  # Assuming e is defined somewhere in your code

        scatter = go.Scatter3d(
            x=[i for i in xs],
            y=[i for i in ys],
            z=[i for i in zs],
            mode="markers+lines",  # Include both markers and lines
            name=f"Axon",
            marker=dict(
                sizemode="diameter",
                size=1,  # Set the size of scatter points for the axon
                color=c,  # Set the color of scatter points for the axon
                line=dict(
                    color="rgba(0, 0, 0, 0.6)",  # Set color to semi-transparent black
                    width=2  # Set the width of the lines
                )
            ),
            line=dict(
                color="rgba(0, 0, 0, 0.6)",  # Set color to semi-transparent black
                width=2  # Set the width of the lines
            )
        )

        layout = go.Layout(
            scene=dict(
                xaxis=dict(title='X [mm]'),
                yaxis=dict(title='Y [mm]'),
                zaxis=dict(title='Z [mm]')
            )
        )
        scatters.append(scatter)

    scatters.append(scatter_cell)

    # Create the figure
    fig = go.Figure(data=scatters, layout=layout)

    # Show the figure
    fig.show()



if __name__ == "__main__":
    swc_file = "/home/localadmin/Documents/CATERPillar/astrocytes_0.swc"
    traj_file = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/test_0.traj"
    read_bin_file(traj_file, swc_file, glial_only=True)