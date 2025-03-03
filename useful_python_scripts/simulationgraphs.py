import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.colors as colors
import os
import glob
from matplotlib.patches import Circle
import matplotlib.colors as mcolors
import math
import random
import pyvista as pv
from tqdm import tqdm
import pandas as pd
from matplotlib.patches import Circle
from multiprocessing import Pool
import copy
from scipy.signal import find_peaks




def read_swc_file(file_path):
    columns = ["id_ax","sph_id", "branch_id", "type", "x", "y", "z", "Rin","Rout", "P"]
    df = pd.read_csv(file_path, sep=' ', names=columns)

    df = df.iloc[1:]
    df["x"] = [float(i) for i in list(df["x"])]
    df["y"] = [float(i) for i in list(df["y"])]
    df["z"] = [float(i) for i in list(df["z"])]
    df["Rout"] = [float(i) for i in list(df["Rout"])]
    df["Rin"] = [float(i) for i in list(df["Rin"])]
    df["id_ax"] = [float(i) for i in list(df["id_ax"])]
    df["sph_id"] = [float(i) for i in list(df["sph_id"])]
    df["P"] = [float(i) for i in list(df["P"])]
    df["branch_id"] = [float(i) for i in list(df["branch_id"])]
    df["type"] = [str(i) for i in list(df["type"])]
    return df



def radius_histogram(df):
    df['R'] = pd.to_numeric(df['R'], errors='coerce')
    df["Diameter"] = df["R"]*2
    df = df.groupby(by = "id_ax").mean()
    sns.histplot(data=df, x="Diameter", color='blue', bins=30, kde=True)
    plt.xlabel('Diameter')
    plt.ylabel('Frequency')
    plt.title('Diameter Histogram')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

def mean_dist_between_maxima(file_path):
    df = read_swc_file(file_path)

    axon_ids = df['id_ax'].unique()
    mean_distances = []

    # Step 1: Process each axon to find mean distances between peaks
    for axon_id in axon_ids:
        # Extract radius and z values for the current axon
        radius = df.loc[df['id_ax'] == axon_id, 'Rin'].values
        z_values = df.loc[df['id_ax'] == axon_id, 'z'].values

        # Find local maxima in the radius values
        peaks, _ = find_peaks(radius, distance = 1, prominence=0.5)
        
        # Calculate the distances between consecutive peaks based on z-values
        if len(peaks) > 1:
            z_peaks = z_values[peaks]  # Get the corresponding z-values for the peaks
            z_distances = np.diff(z_peaks)  # Calculate the distance between consecutive peaks
            mean_distance = np.mean(z_distances)  # Calculate the mean distance
            mean_distances.append(mean_distance)

    # Step 2: Plot a random radius variation with z and mark the peaks
    random_axon_id = random.choice(axon_ids)  # Pick a random axon ID
    radius_random = df.loc[df['id_ax'] == random_axon_id, 'Rout'].values
    z_random = df.loc[df['id_ax'] == random_axon_id, 'z'].values

    # Find peaks for the random axon
    peaks_random, _ = find_peaks(radius_random)

    # Plot the radius variation with respect to z
    plt.figure(figsize=(10, 6))
    plt.plot(z_random, radius_random, label=f'Axon ID: {random_axon_id}', color='b')
    
    # Mark the peaks on the plot
    plt.plot(z_random[peaks_random], radius_random[peaks_random], "ro", label='Peaks')

    # Add labels and legend
    plt.xlabel('z (Position)')
    plt.ylabel('Radius')
    plt.title('Radius Variation with z for a Random Axon with Local Maxima')
    plt.legend()
    plt.show()

    # Step 3: Create a histogram of mean distances
    plt.hist(mean_distances, bins=100, edgecolor='black')
    plt.xlabel('Mean Distance Between Local Maxima (z)')
    plt.ylabel('Frequency')
    plt.title('Histogram of Mean Distances Between Local Maxima (z) Across Radii')
    plt.show()


def diameter_variation(file_path, num_axons=10, max_z=None):
    df = read_swc_file(file_path)
    print(df)
    plt.figure(figsize=(12, 6))  # Adjust the figure size as needed

    # Get the first 'num_axons' axon IDs
    first_n_axon_ids = df['id_ax'].unique()[:num_axons]
    df_subset = df[df['id_ax'].isin(first_n_axon_ids)]

    # Filter data until the specific 'z' value if provided
    if max_z is not None:
        df_subset = df_subset[df_subset['z'] <= max_z]

    # Convert 'sph_id' to categorical to ensure proper x-axis alignment
    df_subset['z'] = pd.Categorical(df_subset['z'])

    # Calculate diameter from radius
    df_subset['Diameter'] = df_subset['Rin'] * 2

    # Plot smooth curves without individual data points
    sns.lineplot(data=df_subset, x="z", y="Diameter", hue="id_ax", ci=None, legend=False, markers=True)
    
    plt.xlabel("z")
    plt.ylabel("2r")
    if max_z is not None:
        plt.title(f"Sphere Diameter for First {num_axons-1} Axons (until z={max_z})")
    else:
        plt.title(f"Sphere Diameter for First {num_axons-1} Axons")  # Updated title with the number of axons

    plt.show()

def create_subplots(file_path):
    df = read_swc_file(file_path)
    df["Outer Diameter"] = df["Rout"]*2
    df["Inner Diameter"] = df["Rin"]*2
    axons = df.loc[df["type"] == "axon"]
    # add column "myelinatin" with true if Rin = Rout, otherwise false
    axons["myelinated"] = axons["Rin"] != axons["Rout"]
    #axons= axons.loc[axons["myelinated"] == False]
    # calculate covariance of Rout of each axon
    cov_out_mean = axons[["Outer Diameter", "id_ax"]].groupby(by="id_ax").mean()
    cov_out_std = axons[["Outer Diameter", "id_ax"]].groupby(by="id_ax").std()
    cov_out = pd.DataFrame()
    cov_out["Coefficient of Variation"] = cov_out_std/cov_out_mean
    cov_in_mean = axons[["Inner Diameter", "id_ax"]].groupby(by="id_ax").mean()
    cov_in_std = axons[["Inner Diameter", "id_ax"]].groupby(by="id_ax").std()
    cov_in = pd.DataFrame()
    cov_in["Coefficient of Variation"] = cov_in_std/cov_in_mean
    # delete nans
    cov_in = cov_in.dropna()
    cov_out = cov_out.dropna()

    print("length cov_out_mean", len(cov_out_mean))
    print("length cov_out_std", len(cov_out_std))
    print("length cov_in_mean", len(cov_in_mean))
    print("length cov_in_std", len(cov_in_std))
    print("length cov_out", len(cov_out))
    print("length cov_in", len(cov_in))

    # add g-ratio column = Rin/Rout
    axons["g_ratio"] = axons["Rin"]/axons["Rout"]


    # define figure with 4 subplots
    fig, axs = plt.subplots(2, 2, figsize=(10, 10))
    # plot "Rout" histogram with seaborn
    sns.histplot(data=cov_out_mean, x="Outer Diameter", color='blue',  kde=True, ax=axs[0, 0], label="Outer Diameter", bins=30, binrange = (0, 5))
    sns.histplot(data=cov_in_mean, x="Inner Diameter", color='red',  kde=True, ax=axs[0, 0], label="Inner Diameter", bins=30, binrange = (0, 5))
    axs[0, 0].legend()
    axs[0, 0].set_title('Diameter Histogram')
    axs[0, 0].set_xlabel('Diameter [µm]')
    axs[0, 0].set_ylabel('Count')

    
    # plot Rin wrt g-ratio with seaborn
    sns.scatterplot(data=axons, x="Inner Diameter", y="g_ratio", color = "red", ax=axs[0, 1], label="Inner Diameter")
    sns.scatterplot(data=axons, x="Outer Diameter", y="g_ratio", color = "blue", ax=axs[0, 1], label="Outer Diameter")
    axs[0, 1].set_title('Diameter vs. g-ratio')
    axs[0, 1].set_xlabel('Diameter [µm]')
    axs[0, 1].set_ylabel('g-ratio')
    axs[0, 1].legend()


    # plot covariance of Rout
    sns.histplot(data=cov_out, x= "Coefficient of Variation", color='blue',  kde=True, ax=axs[1, 0], label="Outer Diameter", bins=30, binrange = (0, 1))
    sns.histplot(data=cov_in, x= "Coefficient of Variation", color='red',  kde=True, ax=axs[1, 0], label="Inner Diameter", bins=30, binrange = (0, 1))
    axs[1, 0].legend()
    axs[1, 0].set_title('Coefficient of Variation of Diameter')
    axs[1, 0].set_xlabel('Coefficient of Variation')
    axs[1, 0].set_ylabel('Count')
    

    # plot tortuosity per axon with seaborn = total length / distance between first and last sphere
    tortuosities, radii = tortuosity(axons)

    tort = pd.DataFrame()
    tort["Tortuosity"] = tortuosities
    sns.histplot(data = tort, x= "Tortuosity", ax=axs[1, 1], color= "purple", kde = True)
    # set maximum and minimum values for x-axis
    axs[1, 1].set_xlim([1, 1.4])
    axs[1, 1].set_title('Tortuosity Histogram')
    axs[1, 1].set_xlabel('Tortuosity')
    axs[1, 1].set_ylabel('Count')

    plt.tight_layout()
    plt.show()



    
    
    

def all_tortuosity(ondulation_factors, std_deviations):
    all_dfs = []
    for o in ondulation_factors:
        for std in std_deviations:
            new_df = pd.DataFrame(columns=["Tortuosity", "Radius", "Std", "Ondulation_factor"])
            file_path = f"/home/localadmin/Documents/Melina_branch/Sim_Growth/data_std/std_dev_{std}/ond_factor_{o}/growth_icvf_0.10_cap_24_vox_50_factor_2_0.swc"
            df = read_swc_file(file_path)
            tort, rad = tortuosity(df)
            new_df["Tortuosity"] = tort
            new_df["Radius"] = rad
            new_df["Std"] = [std]*len(new_df)
            new_df["Ondulation_factor"] = [o]*len(new_df)
            all_dfs.append(new_df)
    
    df = pd.concat(all_dfs).reset_index()

    print(df)

    sns.boxplot(data =df, y = "Tortuosity", x = "Std", hue = "Ondulation_factor")
    plt.show()


def tortuosity(df):
    if len(df) == 0:
        print("No axons found")
        return [], []

    nbr_axons = int(df.iloc[len(df)-1 ]["id_ax"])
    tortuosities = []
    radii = []
    for axon in range(nbr_axons):
 
        df_ = df.loc[df["id_ax"]== axon].reset_index()
        if(len(df_)== 0):
            continue
        
        # Calculate Euclidean distances between consecutive spheres
        distances = []

        for i in range(1, len(df_)):

            distance = math.sqrt((df_.at[i, 'x'] - df_.at[i-1, 'x'])**2 +
                                (df_.at[i, 'y'] - df_.at[i-1, 'y'])**2 +
                                (df_.at[i, 'z'] - df_.at[i-1, 'z'])**2)
    
            distances.append(distance)
            

        # Calculate total length of the axon
        total_length = sum(distances)

        # Calculate distance between the first and last sphere
        first_last_distance = math.sqrt((df_.at[0, 'x'] - df_.at[len(df_)-1, 'x'])**2 +
                                        (df_.at[0, 'y'] - df_.at[len(df_)-1, 'y'])**2 +
                                        (df_.at[0, 'z'] - df_.at[len(df_)-1, 'z'])**2)

        if first_last_distance == 0:
            continue
        # Calculate tortuosity
        tortuosity = total_length / first_last_distance
        tortuosities.append(float(tortuosity))
        radii.append(float(df_.at[0,"Rout"]))
   
    return tortuosities, radii


def draw_axons(file_path, glial_only = False):

    N = 5
    colors = mcolors._colors_full_map #dictionary of all colors
    df = read_swc_file(file_path)
    df["color"] = df["id_ax"].apply(lambda x: get_random_element(colors, seed = x)[1])

    # distance from x,y,z to point(0,0,0)
    df["distance_to_point"] = np.linalg.norm(df[["x", "y", "z"]], axis=1)

    # sort by distance to point
    df = df.sort_values(by="distance_to_point")

    if glial_only:
        df = df.loc[df["type"] != "axon"]

    # Create a scatter plot for the axons
    scatter = go.Scatter3d(
        x=df["x"],
        y=df["y"],
        z=df["z"],
        type="scatter3d",
        mode="markers",
        name="Axons",
        marker=dict(
            sizemode="diameter",
            size=df["Rout"]*N,
            color=df["color"],
            line=dict(
                color="rgba(0, 0, 0, 0)",
                width=0
            )
        )
    )


    layout = go.Layout(
        scene=dict(
            xaxis=dict(title='X [µm]'),
            yaxis=dict(title='Y [µm]'),
            zaxis=dict(title='Z [µm]')
        )
    )

    # Create the figure
    fig = go.Figure(data=scatter, layout=layout)
    # Show the figure
    fig.show()

def get_random_element_list(list, seed = 0):
    random.seed(seed)
    return random.choice(list)

def draw_cells(file_path, plot_type="all", axon_indices=None, astrocyte_indices=None):
    """
    Draw structures (axons or astrocytes) in 3D with customizable options.

    Parameters:
    - file_path: str, path to the data file.
    - plot_type: str, "all" to plot everything, "axons" to plot only axons, 
                 or "astrocytes" to plot only astrocytes.
    - axon_indices: list of int, indices of specific axons to plot. If None, all axons are plotted.
    """
    N = 10  # Size scaling factor
    colors = mcolors.CSS4_COLORS

    # Step 1: Read data
    df = read_swc_file(file_path)

    # Step 2: Assign colors to axons
    df["color"] = df["id_ax"].apply(lambda x: get_random_element(colors, seed=x)[1])

    # Step 3: Prepare data for plotting
    if plot_type == "axons" or (plot_type == "all" and axon_indices):
        df = df[df["type"] == "axon"]
        if axon_indices is not None:
            id_axs = df["id_ax"].unique()
            axon_indices = [id_axs[i] for i in axon_indices]
            df = df[df["id_ax"].isin(axon_indices)]
    elif plot_type == "astrocytes":
        df = df[df["type"] != "axon"]
        if astrocyte_indices is not None:
            id_axs = df["id_ax"].unique()
            astrocyte_indices = [id_axs[i] for i in astrocyte_indices]
            df = df[df["id_ax"].isin(astrocyte_indices)]
    elif plot_type == "all":
        pass
    else:
        raise ValueError("Invalid plot_type. Choose 'all', 'axons', or 'astrocytes'.")

    if (len(df) == 0):
        raise ValueError("No data to plot.")
    # Step 4: Clean and sort data
    df["distance_to_point"] = np.linalg.norm(df[["x", "y", "z"]], axis=1)
    df = df.sort_values(by="distance_to_point")
    df["Rout"] = df["Rout"].fillna(0)  # Handle missing sizes
    print(df)
    # Step 5: Create 3D scatter plot
    scatter = go.Scatter3d(
        x=df["x"],
        y=df["y"],
        z=df["z"],
        mode="markers",
        name="Structures",
        marker=dict(
            sizemode="diameter",
            size=df["Rout"] * N,
            color=df["color"],
            opacity=0.85,
            line=dict(color=df["color"], width=0)
        )
    )

    # Step 6: Define layout
    layout = go.Layout(
        scene=dict(
            xaxis=dict(title='X [µm]', showbackground=False, showgrid=False, showline=False, showticklabels=False),
            yaxis=dict(title='Y [µm]', showbackground=False, showgrid=False, showline=False, showticklabels=False),
            zaxis=dict(title='Z [µm]', showbackground=False, showgrid=False, showline=False, showticklabels=False),
            aspectmode="auto"
        ),
        paper_bgcolor="black",
        plot_bgcolor="black",
        showlegend=False
    )

    # Step 7: Show the plot
    fig = go.Figure(data=[scatter], layout=layout)
    fig.show()


def draw_axons_black_white(file_path, glial_only=False):
    N = 20  # Size scaling factor
    colors = mcolors.CSS4_COLORS  # Use CSS colors for better compatibility

    # Read and preprocess data
    df = read_swc_file(file_path)
    
    # Color assignment (for spheres we will use white)
    df["color"] = "white"
    
    grey_colors = ["#d3d3d3", "#a9a9a9", "#696969", "#808080", "#778899", "#708090", "#2f4f4f", "#708090", "#778899", "#808080", "#696969", "#a9a9a9", "#d3d3d3"]
    df["color"] = list(map(lambda x, y, z: "white" if x == "axon" else y, list(df["type"]), list(df["color"]), list(df["id_ax"])))

    # Filter for axons in a specific region
    df_ = df.loc[df["type"] == "axon"]
    df_ = df_.groupby("id_ax").first().reset_index()
    id_axs = df_["id_ax"].unique()

    # Filter axons by selected id_axs
    df_axons = df.loc[df["type"] == "axon"]
    df_axons = df_axons.loc[df_axons["id_ax"].isin(id_axs)]

    # Glial cells selection
    df_glial = df.loc[df["type"] != "axon"]
    # df_glial = df_glial.loc[df_glial["id_ax"]==1]
    df = pd.concat([df_axons, df_glial])

    # Calculate distance from (0,0,0) for sorting purposes
    df["distance_to_point"] = np.linalg.norm(df[["x", "y", "z"]], axis=1)
    df = df.sort_values(by="distance_to_point")

    if glial_only:
        df = df.loc[df["type"] != "axon"]

    # Ensure marker size is valid
    df["Rout"] = df["Rout"].fillna(0)  # Replace NaN or missing sizes with a default value

    # Create a 3D scatter plot
    scatter = go.Scatter3d(
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

    # Layout configuration with black background
    layout = go.Layout(
        scene=dict(
            xaxis=dict(title='X [µm]', showbackground=False, showgrid=False, showline=False, showticklabels=False),
            yaxis=dict(title='Y [µm]', showbackground=False, showgrid=False, showline=False, showticklabels=False),
            zaxis=dict(title='Z [µm]', showbackground=False, showgrid=False, showline=False, showticklabels=False),
            aspectmode="auto"  # Ensure the aspect ratio is set correctly
        ),
        paper_bgcolor='black',  # Black background
        plot_bgcolor='black',   # Black background for the plot
        showlegend=False
    )

    # Create the figure
    fig = go.Figure(data=[scatter], layout=layout)

    # Show the figure
    fig.show()
def get_spheres_array(df):
    axons = []
    current_axon_id = None
    current_axon = []
    
    df['x'] = pd.to_numeric(df['x'], errors='coerce')
    df['y'] = pd.to_numeric(df['y'], errors='coerce')
    df['z'] = pd.to_numeric(df['z'], errors='coerce')
    df['Rout'] = pd.to_numeric(df['Rout'], errors='coerce')

    for _, row in df.iterrows(): # loops over each row of the df
        id_ax = row["id_ax"]
        x = row["x"]
        y = row["y"]
        z = row["z"]
        r = row["Rout"]  
        type_ = row["type"]

        if id_ax != current_axon_id: # passing to next axon 
            if current_axon_id is not None:
                axons.append(current_axon) # full list
                current_axon = [] # emptying list
            current_axon_id = id_ax # update axon number

        current_axon.append([x, y, z, r, type_])

    if current_axon: # if the list is not empty
        axons.append(current_axon) # last axon

    return axons



def draw_circles(center_radii):
    """
    Draw circles in a 2D plot.

    Parameters:
        center_radii (list of tuples): List of tuples, where each tuple contains (x, y, radius).
        ax (matplotlib.axes._axes.Axes, optional): Axes object to draw the circles on. If not provided, a new plot will be created.
        **kwargs: Additional keyword arguments to customize the appearance of circles.
    """

    fig, ax = plt.subplots()

    for x, y, radius, color in center_radii:
        circle = Circle((x, y), radius, color=color, alpha=0.5)
        ax.add_patch(circle)

    #ax.set_aspect('equal', adjustable='datalim')  # Equal aspect ratio

    return ax

def find_closest_to_value(list_of_lists, value):
    closest_diff = float('inf')
    closest_list = None
    
    for e, sublist in enumerate(list_of_lists):
        if len(sublist) >= 3:

            diff = abs(sublist[2] - value)
            if diff < closest_diff:
                closest_diff = diff
                index = e
                closest_value = sublist[2]
    
    return index, closest_value

def draw_spheres(file_path, limit, z_slice):
    df = read_swc_file(file_path)
    colors = mcolors._colors_full_map #dictionary of all colors
    circles = []

    axons = df.loc[df["type"] == "axon"]
    for axon in axons["id_ax"].unique():
        axon_i = axons.loc[axons["id_ax"] == axon]


        axon_in_slice = axon_i.loc[(axon_i["z"]-axon_i["Rin"] < z_slice) & (z_slice < axon_i["z"]+axon_i["Rin"])]
        random_key, random_value = get_random_element(colors, seed = axon)
        c = random_value
        for i, row in axon_in_slice.iterrows():
            x = row["x"]
            y = row["y"]
            z = row["z"]
            r = row["Rin"]
            
            Rnew = math.sqrt(r*r - (z - z_slice)**2)
            circles.append((x, y, Rnew, c))

        

    glial_df = df.loc[df["type"] != "axon"]
    for glial in glial_df["id_ax"].unique():
        glial_i = glial_df.loc[glial_df["id_ax"] == glial]
        glial_in_slice = glial_i.loc[(glial_i["z"]-glial_i["Rout"] < z_slice) & (z_slice < glial_i["z"]+glial_i["Rout"])]
        random_key, random_value = get_random_element(colors, seed = axon)
        c = random_value
        for i, row in glial_in_slice.iterrows():
            x = row["x"]
            y = row["y"]
            z = row["z"]
            r = row["Rout"]
            
            Rnew = math.sqrt(r*r - (z - z_slice)**2)
            circles.append((x, y, Rnew, c))

    draw_circles(circles)

    # Set plot limits
    plt.xlim(0, limit)
    plt.ylim(0, limit)

    plt.title('2D Circles')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.grid()
    plt.show()

def read_data(filename):
    with open(filename, 'r') as file:
        data = {}
        for line in file:
            key, value = line.strip().split()
            data[key] = float(value)
        return data

def combine_files(file_list):
    combined_data = {}
    for file in file_list:
        data = read_data(file)
        for key, value in data.items():
            if key in combined_data:
                combined_data[key].append(value) # adds element to key list
            else:
                combined_data[key] = [value] # creates new key and adds element
    return combined_data

def vox_time_plot(file_list):
    data = combine_files(file_list)
    # sns.lineplot(x=data['Voxel'], y=data['Duration'])
    # plt.xlabel('Voxel size (µm)')
    # plt.ylabel('Time (s)')
    # plt.title('Time vs. Voxel Size')
    # plt.legend(title='icvf = ' + str(data['icvf'][0]) + '\n' + 'capacity = ' + str(data['Capacity'][0]), loc='upper left')
    # plt.show()
    df = pd.DataFrame(data)


    unique_cap_values = df['Capacity'].unique()

    
    for capacity in unique_cap_values:
        subset_data = df[df['Capacity'] == capacity]
        
        sns.lineplot(x='Voxel', y='Duration', data=subset_data, label=f'capacity = {capacity}')

    plt.xlabel('Voxel size (µm)')
    plt.ylabel('Time (s)')
    plt.title('Time vs. Voxel Size for Different Capacity Values')
    plt.legend(loc='upper left')
    plt.show()

def cap_time_plot(file_list):
    data = combine_files(file_list)
    keys = ['Duration', 'Capacity', 'Voxel']
    data = {key: data[key] for key in keys}

    df = pd.DataFrame(data)
    df['Duration'] =  df['Duration']/60

    df_pivot = df.pivot_table(index='Voxel', columns='Capacity', values='Duration', aggfunc='mean')
    

    # Set the heatmap parameters
    sns.heatmap(df_pivot,
                annot=True,
                fmt=".2f",  # Format for the annotations (optional, adjust as needed)
                cmap='viridis')

    plt.title('Simulation Heat Map')
    plt.xlabel('Capacity')
    plt.ylabel('Voxel size (um)')
    plt.show()

def get_text_from_folder(folder_path,  straight = False):

    txt_files = glob.glob(os.path.join(folder_path, f"*.txt"))
    txt_files.sort()
    list_txt = []
    for txt_file in txt_files:
        if straight and "straight" in txt_file:
            list_txt.append(txt_file)
        elif not straight and "straight" not in txt_file:
            list_txt.append(txt_file)
    return list_txt
    
def get_swc_from_folder(folder_path, icvf):

    txt_files = glob.glob(os.path.join(folder_path, f"growth*"))
    txt_files.sort()
    list_txt = []
    if icvf != None:
        for txt_file in txt_files:
            if (str(icvf) in txt_file):
                list_txt.append(txt_file)
        return list_txt
    else:
        return txt_files

def read_swc_file_np(file_path, axon_nbr):
# Define column names
    columns = ["id_ax", "sph_id", "branch_id", "type", "x", "y", "z", "Rin", "Rout", "P"]
    
    # Read the file using numpy.genfromtxt to handle non-numeric values
    data_ = np.genfromtxt(file_path, skip_header=1, dtype=None, encoding=None)

    # Convert selected columns to float
    data = np.ones((len(data_), len(columns)))*np.nan
    for i in range(len(data_)) :
        if(data_[i][3] != "axon" and int(data_[i][0]) == axon_nbr):
    
            for j in range(len(columns)):
                if (columns[j] != "type"):
                    data[i][j] = float(data_[i][j])

    # Find rows containing NaN values
    nan_mask = np.isnan(data).all(axis=1)

    # Filter out rows with NaN values
    data = data[~nan_mask]


    # Create a dictionary mapping column names to their corresponding data
    data_dict = {column: data[:, i] for i, column in enumerate(columns)}
    
    return data_dict

def get_random_element(dictionary, seed):
    random.seed(seed)
    key, value = random.choice(list(dictionary.items()))
    return key, value

def draw_spheres_pyvista(file_path, axons = True, processes = True, chosen_id = None):
    df = read_swc_file(file_path)
    colors = mcolors._colors_full_map #dictionary of all colors
    plotter = pv.Plotter()

    scale = 500
    if axons:
        df_axons = df[df["type"] == "axon"]
        if (chosen_id is None):
            for id in tqdm(df_axons["id_ax"].unique(), desc="Processing"):
                axon_i = df_axons.loc[df_axons["id_ax"] == id]
                N = np.array(axon_i[["x", "y", "z"]].astype(float))
                Rout = np.array(axon_i["Rout"])
                Rin = np.array(axon_i["Rin"])
                random_key, random_value = get_random_element(colors, seed = id)
                c = random_value

                for i, p in enumerate(N):
                    #plotter.add_mesh(pv.PolyData(p), point_size=R[i]*scale, color=c, render_points_as_spheres=True)
                    if Rin[i] != Rout[i]:
                        plotter.add_points(p, render_points_as_spheres=True, point_size=Rout[i]*scale, color=c, opacity = 0.4)
                        plotter.add_points(p, render_points_as_spheres=True, point_size=Rin[i]*scale, color=c, opacity = 1)
                    else:
                        plotter.add_points(p, render_points_as_spheres=True, point_size=Rout[i]*scale, color=c, opacity = 1)
        else:
            axon_i = df_axons.loc[df_axons["id_ax"] == chosen_id]
            N = np.array(axon_i[["x", "y", "z"]].astype(float))
            Rout = np.array(axon_i["Rout"])
            Rin = np.array(axon_i["Rin"])
            if chosen_id is not None:
                random_key, random_value = get_random_element(colors, seed = chosen_id)
            else:
                random_key, random_value = get_random_element(colors, seed = 0)
            c = random_value

            for i, p in enumerate(N):
                #plotter.add_mesh(pv.PolyData(p), point_size=R[i]*scale, color=c, render_points_as_spheres=True)
                if Rin[i] != Rout[i]:
                    plotter.add_points(p, render_points_as_spheres=True, point_size=Rout[i]*scale, color=c, opacity = 0.4)
                    plotter.add_points(p, render_points_as_spheres=True, point_size=Rin[i]*scale, color=c, opacity = 1)
                else:
                    plotter.add_points(p, render_points_as_spheres=True, point_size=Rout[i]*scale, color=c, opacity = 1)
    
    df_glial = df[df["type"] == "glial"]
    for id in tqdm(df_glial["id_ax"].unique(), desc="Processing"):
        if chosen_id is not None:
            if id != chosen_id:
                continue
        glial_i = df_glial.loc[df_glial["id_ax"] == id]
        N = np.array(glial_i[["x", "y", "z"]].astype(float))
        R = np.array(glial_i["Rout"])
        if chosen_id is not None:
            random_key, random_value = get_random_element(colors, seed = chosen_id)
        else:
            random_key, random_value = get_random_element(colors, seed = id)
        c = random_value
        for i, p in enumerate(N):
            plotter.add_points(p, render_points_as_spheres=True, point_size=R[i]*scale, color=c)
                
            #plotter.add_mesh(pv.PolyData(p), point_size=R[i]*scale, color=c, render_points_as_spheres=True)
    if processes:
        df_glial_ramification = df[df["type"] == "glialRamification"]
        for id in tqdm(df_glial_ramification["id_ax"].unique(), desc="Processing"):
            if chosen_id is not None:
                if id != chosen_id:
                    continue
     
            glial_ramification_i = df_glial_ramification.loc[df_glial_ramification["id_ax"] == id]
            N = np.array(glial_ramification_i[["x", "y", "z"]].astype(float))
            R = np.array(glial_ramification_i["Rout"])
            if chosen_id is not None:
                random_key, random_value = get_random_element(colors, seed = chosen_id)
            else:
                random_key, random_value = get_random_element(colors, seed = id)
            c = random_value
            print(len(N))
            for i, p in enumerate(N):
                plotter.add_points(p, render_points_as_spheres=True, point_size=R[i]*scale, color=c)
                #plotter.add_mesh(pv.PolyData(p), point_size=R[i]*scale, color=c, render_points_as_spheres=True)

    plotter.show()

def draw_one_axon_pyvista(file_path, chosen_id = None):
    df = read_swc_file(file_path)
    colors = mcolors._colors_full_map #dictionary of all colors
    plotter = pv.Plotter()
    df_axons = df[df["type"] == "axon"]

    scale = 150

    axon_i = df_axons.loc[df_axons["id_ax"] == chosen_id]

    N = np.array(axon_i[["x", "y", "z"]].astype(float))
    Rout = np.array(axon_i["Rout"])
    Rin = np.array(axon_i["Rin"])
    random_key, random_value = get_random_element(colors, seed = chosen_id)
    c = random_value

    for i, p in enumerate(N):
        #plotter.add_mesh(pv.PolyData(p), point_size=R[i]*scale, color=c, render_points_as_spheres=True)
        if Rin[i] != Rout[i]:
            plotter.add_points(p, render_points_as_spheres=True, point_size=Rout[i]*scale, color="blue", opacity = 0.3)
            plotter.add_points(p, render_points_as_spheres=True, point_size=Rin[i]*scale, color="red", opacity = 1)
        else:
            plotter.add_points(p, render_points_as_spheres=True, point_size=Rout[i]*scale, color="purple", opacity = 0.8)

    plotter.show()

def draw_one_glial_pyvista(file_path, chosen_id = None):
    df = read_swc_file(file_path)
    colors = mcolors._colors_full_map #dictionary of all colors
    plotter = pv.Plotter()
    df_axons = df[df["type"] != "axon"]

    scale = 20

    axon_i = df_axons.loc[df_axons["id_ax"] == chosen_id]

    N = np.array(axon_i[["x", "y", "z"]].astype(float))
    Rout = np.array(axon_i["Rout"])
    Rin = np.array(axon_i["Rin"])
    random_key, random_value = get_random_element(colors, seed = chosen_id)
    c = random_value

    for i, p in enumerate(N):

        plotter.add_points(p, render_points_as_spheres=True, point_size=Rout[i]*scale, color="green", opacity = 0.8)

    plotter.show()


def sholl_intersection(file_path):
    df = read_swc_file(file_path)
    sphere_around_soma_radii = [5, 7, 10, 15, 20 , 25, 30, 40, 50, 60, 80]
    intersections_list_all = []
    glial_df = df.loc[df["type"] != "axon"]
    for glial in glial_df["id_ax"].unique():
        intersections_list = []
        glial_i = glial_df.loc[glial_df["id_ax"] == glial]
        soma_glial = glial_i.loc[glial_i["type"] == "glialSoma"]
        print(soma_glial)
        soma_position = np.array(soma_glial[["x", "y", "z"]].astype(float))[0]
        print(soma_position)
        processes_glial = glial_i.loc[glial_i["type"] == "glialRamification"]
        processes_glial["distance_to_soma"] = np.linalg.norm(np.array(processes_glial[["x", "y", "z"]].astype(float)) - soma_position, axis=1)
        # sort by distance to soma
        processes_glial = processes_glial.sort_values(by="distance_to_soma")
        # keep only first and last element of each branch_id
        processes_glial_first = processes_glial.drop_duplicates(subset=["branch_id"], keep="first")
        processes_glial_last = processes_glial.drop_duplicates(subset=["branch_id"], keep="last")

        distance_first = np.array(processes_glial_first[["distance_to_soma"]].astype(float))
        distance_last = np.array(processes_glial_last[["distance_to_soma"]].astype(float))
        # find how many times glial cross each sphere
        for r in sphere_around_soma_radii:
            intersections = 0
            for i, p in enumerate(distance_first):
                if distance_first[i] < r  and distance_last[i] > r :
                    intersections += 1
            intersections_list.append(intersections)
        intersections_list_all.append(intersections_list)
    return intersections_list_all, sphere_around_soma_radii

def sholl_intersections(file_path_WM, file_path_GM = None):

    intersections_list_WM, sphere_around_soma_radii = sholl_intersection(file_path_WM)
    if file_path_GM is not None:
        intersections_list_GM, sphere_around_soma_radii = sholl_intersection(file_path_GM)

        # plot sholl intersections wrt distance from soma
        plt.figure()
        for intersections, tissue in zip([intersections_list_WM, intersections_list_GM], ["With axons", "Without axons"]):
            if (tissue == "With axons"):
                plt.plot(sphere_around_soma_radii, np.mean(intersections, axis=0), color='blue', label=tissue)
            else:
                plt.plot(sphere_around_soma_radii, np.mean(intersections, axis=0), color='red', label=tissue)
            for intersections in intersections:
                if (tissue == "With axons"):
                    plt.plot(sphere_around_soma_radii, intersections, alpha=0.1, color='blue')
                else:
                    plt.plot(sphere_around_soma_radii, intersections, alpha=0.1, color='red')
        plt.xlabel("Distance from soma (µm)")
        plt.ylabel("Intersections")
        plt.title("Sholl intersections")
        plt.legend()
        plt.show()
    else:
        plt.figure()
        plt.plot(sphere_around_soma_radii, np.mean(intersections_list_WM, axis=0), color='blue')
        for intersections in intersections_list_WM:
            plt.plot(sphere_around_soma_radii, intersections, alpha=0.1, color='blue')
        plt.xlabel("Distance from soma (µm)")
        plt.ylabel("Intersections")
        plt.title("Sholl intersections")
        plt.show()

def varying_sholl(files, stds, lengths):
    all_intersections = []
    #colors = ["blue", "red", "purple", "orange", "black", "yellow", "pink", "brown", "grey", "cyan"]
    colors = ["blue", "black", "purple", "grey",  "cyan"]
    for i, file in enumerate(files):
        intersections_list, sphere_around_soma_radii = sholl_intersection(file)
        all_intersections.append(intersections_list)

    plt.figure()
    for i, intersections in enumerate(all_intersections):
        c = colors[i]
        plt.plot(sphere_around_soma_radii, np.mean(intersections, axis=0),  label = f"Mean process length = {lengths[i]}, std = {stds[i]}", color= c)
        for intersections in intersections:
            plt.plot(sphere_around_soma_radii, intersections, alpha=0.1, color= c)
    plt.xlabel("Distance from soma (µm)")
    plt.ylabel("Intersections")
    plt.title("Sholl intersections")
    plt.legend()
    plt.show()

def vox_size_analysis(files):
    df = pd.DataFrame()
    for file in files:
        data = read_data(file)
        df = pd.concat([df, pd.DataFrame(data, index=[0])])
    df["Log(Duration)"] = np.log(df["Duration"])
    print(df)
    #plot
    sns.boxplot(data = df, y = "Log(Duration)", x = "Voxel")
    sns.swarmplot(data = df, y = "Log(Duration)", x = "Voxel", color="black")
    plt.ylabel("Log (Duration (s))")
    plt.xlabel("Voxel cube length (µm)")
    plt.show()


def tortuosity_plot(folder_path):
    
    folders =[ "std_0.1", "std_0.2", "std_0.3", "std_0.4", "std_0.5"] 
    tortuosities = []
    stds =[] 
    
    for folder in folders:
        print(folder)
        all_dfs = []
        path = os.path.join(folder_path, folder)
        files =  glob.glob(os.path.join(path, f"growth*"))
        for file in files:
            print(file)
            df = read_swc_file(file)
            tort, radii = tortuosity(df)
            tortuosities.extend(tort)
            stds.extend([float(folder.split("_")[1])]*len(tort))
            df_copy = copy.deepcopy(df)
            all_dfs.append(df_copy)


    df_final = pd.DataFrame()
    df_final["Tortuosity"] = tortuosities
    df_final["Std"] = stds
    # Your lmplot with customized colors for dots and regression line
    sns.lmplot(data=df_final, y="Tortuosity", x="Std", order=2,
            scatter_kws={'color': 'black'},   # Set dots to black
            line_kws={'color': 'red'})        # Set line to red

    # Display the plot
    plt.show()

    

if __name__ == "__main__":

    file_path_GM= "/home/localadmin/Documents/CATERPillar/growth_vox_50_factor_4_0.swc"
    file_path_WM = "/home/localadmin/Documents/CATERPillar/growth_vox_50_factor_4_0.swc"
    file_path = "/home/localadmin/Documents/CATERPillar/test.swc"
    # draw_spheres(file_path, 50, 49)
    #tortuosity_plot("/home/localadmin/Documents/CATERPillar/tortuosities")
    # draw_one_glial_pyvista(file_path_GM, chosen_id = 1)
    #draw_spheres_pyvista(file_path, axons = True ,chosen_id = 1)
    #create_subplots(file_path)
    # create_subplots(file_path)
    #draw_axons_with_shading(file_path, glial_only = True)
    draw_cells(file_path, plot_type="axons")
    # draw_spheres(file_path, 150, 30)
    #diameter_variation(file_path, num_axons = 10, max_z = 15)
    #mean_dist_between_maxima(file_path)
    
    #sholl_intersections(file_path_WM, file_path_GM)
    #files = []
    #for v in [30, 50, 100, 200]:
    #    for cap in [20]:
    #        for i in range(10):
    #            files.append(f"/home/localadmin/Documents/CATERPillar/analysis_time/simulation_vox_{v}_cap_{cap}_{i}.swc")
    #vox_size_analysis(files)

    #files = []
    #stds = []
    #lengths = []
    #for std in [10]:
    #    for length in [20,10, 30,40]:
    #        lengths.append(length)
    #        stds.append(std)
    #        files.append(f"/home/localadmin/Documents/CATERPillar/growth_vox_100_factor_2_0_length_{length}_std_{std}.swc")
    #varying_sholl(files, stds, lengths)
    
    #files = []
    #stds = []
    #lengths = []
    #for std in [10, 1, 5, 15, 20]:
    #    for length in [20]:
    #        lengths.append(length)
    #        stds.append(std)
    #        files.append(f"/home/localadmin/Documents/CATERPillar/growth_vox_100_factor_2_0_length_{length}_std_{std}.swc")
    
    
    #varying_sholl(files, stds, lengths)
    #folder_path = "/home/localadmin/Documents/CATERPillar/tortuosities/"
    #tortuosity_plot(folder_path)

    


    