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

def extract_values_cap(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        values = {}
        for line in lines:
            key, value = line.strip().split(' ')
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
    
def varying_cap(icvfs, vox_size, repetitions):
    capacities = [3,6,12,18, 24]

    df = pd.DataFrame(columns = ['Duration', 'Capacity', 'ICVF'])
    for icvf in icvfs:
        for rep in range(repetitions):
            for cap  in capacities:

                path = f"/home/localadmin/Documents/Melina_branch/Sim_Growth/data/simulation_icvf_{icvf}_cap_{cap}_vox_{vox_size}_{rep}.swc"
                values = extract_values_cap(path)
                # Add new line using loc method
                values["ICVF"]= icvf
                df.loc[len(df)] = values

    print(df)
    # Create a bar plot using Seaborn
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 6))
    sns.violinplot(x="Capacity", y="Duration", hue = "ICVF", data=df, palette="Blues_d")

    # Add labels and title
    plt.xlabel("Capacity")
    plt.ylabel("Duration [s]")
    plt.title("Duration vs Capacity")

    # Show the plot
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

varying_cap(["0.30", "0.50"], 50, 3)
#varying_icvf()