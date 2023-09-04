import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.colors as colors
import os
import chardet 
import glob
from matplotlib.patches import Circle

def read_swc_file(file_path):
    columns = ["ax_id","sph_id", "type", "x", "y", "z", "R", "P"]
    df = pd.read_csv(file_path, sep=' ', names=columns)
    df = df.iloc[1:]

    return df

def radius_file(file_path): 
    columns = ["Type1", "ax_id", "Type2", "x", "y", "z", "Distance", "R", "R0", "Tortuosity"]
    df = pd.read_csv(file_path, sep='\s+', names=columns)
    return df

def radius_histogram(df):
    df['R'] = pd.to_numeric(df['R'], errors='coerce')
    sns.histplot(data=df, x="R", color='blue', bins=30, kde=True)
    plt.xlabel('Radius')
    plt.ylabel('Frequency')
    plt.title('Radius Histogram')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

def diameter_variation(file_path, num_axons=10, max_z=None):
    df = radius_file(file_path)
    plt.figure(figsize=(12, 6))  # Adjust the figure size as needed
    df['R'] = pd.to_numeric(df['R'], errors='coerce')

    # Get the first 'num_axons' axon IDs
    first_n_axon_ids = df['ax_id'].unique()[:num_axons]
    df_subset = df[df['ax_id'].isin(first_n_axon_ids)]

    # Filter data until the specific 'z' value if provided
    if max_z is not None:
        max_z_str = str(max_z)
        df_subset = df_subset[df_subset['Distance'] <= max_z_str]

    # Convert 'sph_id' to categorical to ensure proper x-axis alignment
    df_subset['Distance'] = pd.Categorical(df_subset['Distance'])

    # Calculate diameter from radius
    df_subset['Diameter'] = df_subset['R'] * 2

    # Plot smooth curves without individual data points
    sns.lineplot(data=df_subset, x="Distance", y="Diameter", hue="ax_id", ci=None, legend=False)
    
    plt.xlabel("z")
    plt.ylabel("2r")
    if max_z is not None:
        plt.title(f"Sphere Diameter for First {num_axons-1} Axons (until z={max_z})")
    else:
        plt.title(f"Sphere Diameter for First {num_axons-1} Axons")  # Updated title with the number of axons

    plt.xticks([])  # Remove the x-axis labels
    plt.show()

def create_subplots(df, num_axons=10, max_z=None):
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(12, 6))

    # Diameter variation with the first subplot
    plt.sca(axes[0])
    df['R'] = pd.to_numeric(df['R'], errors='coerce')

    first_n_axon_ids = df['ax_id'].unique()[-num_axons:] # get the first 'num_axons' axon IDs
    df_subset = df[df['ax_id'].isin(first_n_axon_ids)]

    # Filter data until the specific 'z' value if provided
    if max_z is not None:
        max_z_str = str(max_z)
        df_subset = df_subset[df_subset['Distance'] <= max_z_str]

    df_subset['Distance'] = pd.Categorical(df_subset['Distance']) # ensures proper x-axis alignment

    df_subset['Diameter'] = df_subset['R'] * 2

    sns.lineplot(data=df_subset, x="Distance", y="Diameter", hue="ax_id", ci=None, legend=False)
    plt.xlabel("z (µm)")
    plt.ylabel("2r (µm)")
    if max_z is not None:
        plt.title(f"Sphere Diameter for First {num_axons-1} Axons (until z={max_z})")
    else:
        plt.title(f"Sphere Diameter for First {num_axons-1} Axons") 

    plt.xticks([])  # Remove the x-axis labels

    # Radius histogram with the second subplot
    plt.sca(axes[1])

    df['R'] = pd.to_numeric(df['R'], errors='coerce')
    df['Diameter'] = df['R'] * 2
    sns.histplot(data=df, x="Diameter", color='blue', bins=30, kde=True)
    plt.xlabel('Diameter (µm)')
    plt.ylabel('Frequency')
    plt.title('Diameter Histogram')
    plt.xticks(rotation=45)
    
    # Coeff variation for the third subplot
    plt.sca(axes[2])
    df['R'] = pd.to_numeric(df['R'], errors='coerce')
    df['Diameter'] = df['R'] * 2
    # Calculate the coefficient of variation for each 'ax_id'
    cv_data = df.groupby('ax_id')['Diameter'].agg(lambda x: (x.std() / x.mean()))

    # Plot the histogram of coefficient of variation
    sns.histplot(data=cv_data, color='blue', bins=30, kde=True)
    plt.xlabel('CV of diameter')
    plt.ylabel('Frequency')
    plt.title('CV Histogram')
    plt.xticks(rotation=45)

    # Add an overall title for the entire figure
    plt.suptitle('Substrate Analysis', fontsize=16)

    # Adjust the layout to prevent overlapping of titles and labels
    plt.tight_layout()

    # Display the figure
    plt.show()

def tortuosity_plot(df):
    df['radius'] = pd.to_numeric(df['R'], errors='coerce')
    df['Tortuosity'] = pd.to_numeric(df['Tortuosity'], errors='coerce')

    g = sns.JointGrid(data=df, x="radius", y="Tortuosity", space=0)
    g.plot_joint(sns.kdeplot, fill=True, cmap="rocket")
    g.plot_marginals(sns.histplot, color='blue', alpha=1, bins=50)

    plt.xlabel('radius')
    plt.ylabel('Tortuosity')
    plt.title('Tortuosity plot')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

def draw_axons(df):
    axons = get_spheres_array(df)[1:]
    scatters = []
    colours = colors.qualitative.Plotly[:10]
    for e, axon in enumerate(axons):
        c = colours[e % 10]
        # Create a scatter plot for the axon
        scatter = go.Scatter3d(
            x=[s[0] for s in axon],
            y=[s[1] for s in axon],
            z=[s[2] for s in axon],
            mode="markers",
            name=f"Axon {e}",
            marker=dict(
                sizemode="diameter",
                size=[s[3] * 10 for s in axon],  # Set the size of scatter points for the axon
                color=c,  # Set the color of scatter points for the axon
                line=dict(
                    color="rgba(0, 0, 0, 0)",  # Set color to transparent (alpha=0)
                    width=0  # Set width to 0 to remove the contour
                )
            )
        )
        scatters.append(scatter)

    layout = go.Layout(
        scene=dict(
            xaxis=dict(title='X [µm]'),
            yaxis=dict(title='Y [µm]'),
            zaxis=dict(title='Z [µm]')
        )
    )

    # Create the figure
    fig = go.Figure(data=scatters, layout=layout)
    # Show the figure
    fig.show()

def get_spheres_array(df):
    axons = []
    current_axon_id = None
    current_axon = []
    
    df['x'] = pd.to_numeric(df['x'], errors='coerce')
    df['y'] = pd.to_numeric(df['y'], errors='coerce')
    df['z'] = pd.to_numeric(df['z'], errors='coerce')
    df['R'] = pd.to_numeric(df['R'], errors='coerce')

    for _, row in df.iterrows(): # loops over each row of the df
        ax_id = row["ax_id"]
        x = row["x"]
        y = row["y"]
        z = row["z"]
        r = row["R"]  

        if ax_id != current_axon_id: # passing to next axon 
            if current_axon_id is not None:
                axons.append(current_axon) # full list
                current_axon = [] # emptying list
            current_axon_id = ax_id # update axon number

        current_axon.append([x, y, z, r])

    if current_axon: # if the list is not empty
        axons.append(current_axon) # last axon

    return axons

def draw_circles(center_radii, ax=None, **kwargs):
    """
    Draw circles in a 2D plot.

    Parameters:
        center_radii (list of tuples): List of tuples, where each tuple contains (x, y, radius).
        ax (matplotlib.axes._axes.Axes, optional): Axes object to draw the circles on. If not provided, a new plot will be created.
        **kwargs: Additional keyword arguments to customize the appearance of circles.
    """
    if ax is None:
        fig, ax = plt.subplots()

    for x,y,z, radius in center_radii:
        circle = Circle((x,y,z), radius, **kwargs)
        ax.add_patch(circle)

    ax.set_aspect('equal', adjustable='datalim')  # Equal aspect ratio

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

def draw_spheres(file_path, limit, z):
    df = read_swc_file(file_path)

    axons = get_spheres_array(df)
    circles = []

    for axon in axons:

        if (len(axon) > 1):
            index, value =  find_closest_to_value(axon, z)
            distance = np.abs(z-value)
            R = axon[index][3]
            if distance == 0:
                new_r = R
            else:
                if distance < R :
                    new_r = np.sqrt(np.abs(R*R-distance*distance))
                else:
                    continue
            circles.append([axon[index][0],axon[index][1], axon[index][2], new_r])
        
    draw_circles(circles, color='blue', alpha=0.5, linewidth=2)

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


if __name__ == "__main__":


    # create_subplots(radius_file(file2), 51)
    folder = "/home/localadmin/Documents/Melina_branch/Sim_Growth/"
    file_list = get_text_from_folder(folder, straight=False)

    vox_time_plot(file_list)

    cap_time_plot(file_list)

    file_list = get_text_from_folder(folder, straight=True)

    vox_time_plot(file_list)

    #print(file_list[2])
    #draw_spheres(file_list[1], 30, z = 15)
    cap_time_plot(file_list)
    
    #file = "/home/localadmin/Documents/Melina_branch/Sim_Growth/axons_icvf_0.30_cap_1_vox_50.swc"
    #df = radius_file(file)
    #create_subplots(df, num_axons=10, max_z=None)

