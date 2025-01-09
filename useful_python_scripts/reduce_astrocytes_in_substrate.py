import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
from pathlib import Path
from scipy import stats
import random
from simulationgraphs import get_spheres_array, read_swc_file


def compute_icvf_astrocyte_processes(df_astrocyte, limit):
    branches = df_astrocyte['branch_id'].unique()
    all_volume = 0
    soma = df_astrocyte[df_astrocyte['branch_id'] == -1]
    count_discarded = 0

    if soma.empty:
        raise ValueError("Soma is not defined in the input data.")

    soma_center = soma.iloc[0][['x', 'y', 'z']].values
    soma_radius = float(soma.iloc[0]['Rout'])

    for branch in branches:

        if branch == -1:  # Skip the soma itself
            continue

        df_branch = df_astrocyte[df_astrocyte['branch_id'] == branch]

        i = -1

        for sphere1, sphere2 in zip(df_branch.iloc[:-1].itertuples(), df_branch.iloc[1:].itertuples()):

            i = i + 1
            # Calculate distance from sphere1 to the soma
            distance_to_soma = np.linalg.norm(np.array([sphere1.x, sphere1.y, sphere1.z]) - soma_center)

            if distance_to_soma < soma_radius:
                count_discarded += 1
                continue

            # Calculate distance between two consecutive spheres
            distance_between_spheres = np.linalg.norm(
                np.array([sphere1.x, sphere1.y, sphere1.z]) - np.array([sphere2.x, sphere2.y, sphere2.z])
            )
            volume = np.pi * (sphere1.Rout**2 + sphere2.Rout**2 + sphere2.Rout * sphere1.Rout) * distance_between_spheres / 3
            # Add volume using the formula for frustum volume
            all_volume += (volume)

            
    print(f"Total volume /tot_vol: {all_volume/(limit**3)}")
    print(f"Total volume: {all_volume}")
    print(f"Discarded {count_discarded} spheres")

    return all_volume

def compute_icvf_astrocyte_somas(df_astrocyte) :

    df_soma = df_astrocyte[df_astrocyte['branch_id'] == -1]
    volume = float(df_soma['Rout'])**3 * np.pi * 4/3
    return volume

def write_new_swc (path, df):
    columns =  list(df.columns)
    first_line = (" ").join(columns) + "\n"
    with open(path, 'w') as file:
        file.write(first_line)
        for i in range(len(df)):
            line = df.iloc[i]
            file.write(f"{line['ax_id']} {line['sph_id']} {line['branch_id']} {line['type']} {line['x']} {line['y']} {line['z']} {line['Rin']} {line['Rout']} {line['P']}\n")
    

if __name__ == "__main__":

    # Path to the sphere file
    desired_astrocyte_processes_icvf = 1
    desired_soma_icvf = desired_astrocyte_processes_icvf/10
    limit = 150

    swc_file = "/home/localadmin/Documents/CATERPillar/example_morphology/astrocytes.swc"
    new_swc_file = f"/home/localadmin/Documents/CATERPillar/test_{desired_astrocyte_processes_icvf}.swc"

    total_volume = (limit)**3

    print("total volume: ", total_volume)

    df = read_swc_file(swc_file)

    df_astrocytes_to_keep = pd.DataFrame(columns=df.columns)

    df_astrocytes = df[df['type'] != "axon"]
    df_axons = df[df['type'] == "axon"]
    astrocyte_nbrs = df_astrocytes['ax_id'].unique()
    stop = False
    while (stop == False):
        print("Attempting to reach desired ICVF")
        reached_icvf = 0
        reached_soma_icvf = 0
        #shuffle the astrocytes
        random.shuffle(astrocyte_nbrs)
        for astrocyte_nbr in astrocyte_nbrs:
            print(f"Computing ICVF for astrocyte {astrocyte_nbr}")
            df_astrocyte = df_astrocytes[df_astrocytes['ax_id'] == astrocyte_nbr]
            soma = df_astrocyte[df_astrocyte['branch_id'] == -1]
            icvf_soma = compute_icvf_astrocyte_somas(df_astrocyte) / total_volume
            icvf = compute_icvf_astrocyte_processes(df_astrocyte, limit)/total_volume
            #if reached_icvf < desired_astrocyte_processes_icvf:
            reached_icvf += icvf
            reached_soma_icvf += icvf_soma
            print(f"ICVF for processes of astrocyte {astrocyte_nbr}: {icvf}, soma radius: {float(soma['Rout'])}")
            df_astrocytes_to_keep = pd.concat([df_astrocytes_to_keep, df_astrocyte])
            #else:
                #break
        
        # error margin of 0.0001
        #if reached_soma_icvf < desired_soma_icvf - 0.0001:
        #    stop = False
        #else:
        #    stop = True

        stop = True
            
    print(f"Intravolume Fraction for astrocyte somas: {reached_soma_icvf}")
    print(f"Intravolume Fraction for astrocyte processes: {reached_icvf}")
    new_df = pd.concat([df_astrocytes_to_keep, df_axons])

    #write_new_swc (new_swc_file, new_df)
