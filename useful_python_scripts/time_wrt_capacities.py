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
from simulationgraphs import read_swc_file, draw_circles

def extract_values_cap(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

        values = {}
        for line in lines:
     
            key= line.strip().split(' ')[0]
            value = line.strip().split(' ')[1] 
            if key in ['Duration', 'Capacity']:
                values[key] = int(value)
        return values
    
def extract_values_vox(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        values = {}
        for line in lines:
            key, value = line.strip().split(' ')
            if key in ['Duration', 'Voxel']:
                values[key] = int(value)
        return values
    
def extract_values_icvf(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        values = {}
        for line in lines:
            key, value = line.strip().split(' ')
            if key in ['Duration', 'icvf']:
                values[key] = float(value)
        return values
    
def normalise_column(col, data):
    # Calculate the average "adc [um²/ms]" for type = "cylinders" for each combination
    average = data[data['Capacity'] == 1][col].mean()
    # Divide "adc [um²/ms]" by the average value
    title = 'normalized_'+col
    data[title] = data[col] / average
    return data

def varying_cap(icvfs, vox_size, repetitions):
    capacities = [1, 3, 6, 12, 18, 24]

    datas = []
    for icvf in icvfs:
        df = pd.DataFrame(columns=['Duration', 'Capacity', 'ICVF'])
        for rep in range(repetitions):
            for cap in capacities:
                path = f"/home/localadmin/Documents/Melina_branch/Sim_Growth/data/simulation_icvf_{icvf}_cap_{cap}_vox_{vox_size}_{rep}.swc"
                values = extract_values_cap(path)
                values["ICVF"] = icvf
                df.loc[len(df)] = values

        df = normalise_column("Duration", df)
        datas.append(df)

    df = pd.concat(datas)

    df["Duration"] = df["Duration"]/60

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))

    sns.set(style="whitegrid")

    sns.lineplot(x="Capacity", y="normalized_Duration", hue = "ICVF",data=df.reset_index(), ax=axes[0], ci="sd")


    axes[0].set_xlabel("Capacity")
    axes[0].set_ylabel("Normalized Run-time")
    axes[0].set_title("Normalized Run-time vs Capacity")

    df_ = df.loc[df["Capacity"] == 24].reset_index()

    sns.set(style="whitegrid")
    sns.barplot(y="Duration", x="ICVF", hue = "ICVF", data=df_, ax=axes[1])

    axes[1].set_xlabel("ICVF")
    axes[1].set_ylabel("Run-time [min]")
    axes[1].set_title("Run-time vs ICVF")

    plt.suptitle('Run-time vs ICVF and Capacity', fontsize=16)

    plt.tight_layout()
    plt.show()


def varying_icvf():
    icvfs = ["0.30","0.50", "0.70"]

    df = pd.DataFrame(columns = ['Duration', 'icvf'])

    for rep in range(3):
        for icvf  in icvfs:

            path = f"/home/localadmin/Documents/Melina_branch/Sim_Growth/data/simulation_icvf_{icvf}_cap_24_vox_50_{rep}.swc"
            values = extract_values_icvf(path)
            # Add new line using loc method
            df.loc[len(df)] = values

    print(df)
    df["icvf"] = [round(i,1) for i in list(df["icvf"])]
    # Create a bar plot using Seaborn
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 6))
    sns.swarmplot(x="icvf", y="Duration", data=df, palette="Blues_d")

    # Add labels and title
    plt.xlabel("icvf")
    plt.ylabel("Duration [s]")
    plt.title("Duration vs icvf")

    # Show the plot
    plt.show()

def varying_vox():
    voxs = [30,50,100]

    df = pd.DataFrame(columns = ['Duration', 'Voxel'])

    for rep in range(2):
        for vox  in voxs:

            path = f"/home/localadmin/Documents/Melina_branch/Sim_Growth/data/simulation_icvf_0.30_cap_24_vox_{vox}_{rep}.swc"
            values = extract_values_vox(path)
            # Add new line using loc method
            df.loc[len(df)] = values


    # Create a bar plot using Seaborn
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 6))
    sns.swarmplot(x="Voxel", y="Duration", data=df, palette="Blues_d")

    # Add labels and title
    plt.xlabel("Voxel Size")
    plt.ylabel("Duration [s]")
    plt.title("Duration vs Voxel size")

    # Show the plot
    plt.show()



def draw_spheres(file_paths, z_slice, icvfs):
    fig, axes = plt.subplots(nrows=1, ncols=len(file_paths), figsize=(12, 12))
    colors = sns.color_palette()  # Get default Seaborn color palette
    for e, file_path in enumerate(file_paths):
        df = read_swc_file(file_path)
        df["In_slice"] = list(map(lambda r, z: True if z_slice > (z - r) and z_slice < (z + r) else False, list(df["R"]), list(df["z"])))
        sliced_df = df.loc[df["In_slice"] == True]
        sliced_df["R_eff"] = list(map(lambda r, z: np.sqrt(np.abs(r * r - (z - z_slice) * (z - z_slice))), list(sliced_df["R"]), list(sliced_df["z"])))

        circles = []
        for i in range(len(sliced_df)):
            circles.append([list(sliced_df["x"])[i], list(sliced_df["y"])[i], list(sliced_df["R_eff"])[i]])
        
        color = colors[e % len(colors)]  # Cycle through default Seaborn palette colors
        draw_circles(circles, color=color, linewidth=2, ax=axes[e])
        axes[e].set_xlabel("X [um]")
        axes[e].set_ylabel("Y [um]")
        axes[e].set_title(f"ICVF = {icvfs[e]}%")
        
        # Set explicit axis limits
        axes[e].set_xlim(0, 50)  # Update with appropriate limits for X
        axes[e].set_ylim(0, 50)  # Update with appropriate limits for Y

    plt.grid()
    plt.tight_layout()  # Ensures tight layout to prevent overlap
    plt.show()

def cap_wrt_run_time(folder_path):
    # get all the files in the folder
    folder_path_50 = folder_path + "/packing_0.5/"
    folder_path_30 = folder_path + "/packing_0.3/"
    folder_path_70 = folder_path + "/packing_0.7/"
    all_folders = [folder_path_50, folder_path_30, folder_path_70]
    all_packings = [50, 30, 70]
    all_dfs = []
    for folder, packing in zip(all_folders, all_packings):
        all_values = []
        files= glob.glob(folder + "/*.swc")
        for file in files:
            values = extract_values_cap(file)
            all_values.append(values)
        df = pd.DataFrame(all_values, columns = ['Duration', 'Capacity'])
        df["packing"] = [packing]*len(df)
        df["Normalized Duration"] = df["Duration"]/np.mean(df.loc[df["Capacity"] == 1, "Duration"].values)
        all_dfs.append(df)
    df = pd.concat(all_dfs)
    df = df.reset_index()
    sns.lineplot(x="Capacity", y="Normalized Duration", hue = "packing", data=df)
    plt.show()
    sns.lineplot(x="packing", y="Duration",  data=df.loc[df["Capacity"]==200 ] )
    plt.show()
   


if __name__ == "__main__":
    folder_path = "/home/localadmin/Documents/CATERPillar/capacity_runtime"
    cap_wrt_run_time(folder_path)

