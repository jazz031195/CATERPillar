import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.colors as colors
import os
import chardet 
import glob
from simulationgraphs import read_swc_file



def interpolate_axons_swc(df, factor):

    num = factor
    nbr_spheres_to_add_inbetween = num -1

    print("Axon")
    df = df.loc[df["type"]== "axon"]

    # Create a DataFrame containing NaN values with the same columns as the original DataFrame
    nan_df = pd.DataFrame(np.nan, index=np.arange(len(df)*nbr_spheres_to_add_inbetween), columns=df.columns)
    nan_df["index"] = nan_df.index/nbr_spheres_to_add_inbetween
    nan_df.index = nan_df["index"]
    # nan_df["old"] =["false"]*len(nan_df)      
    del nan_df["index"]

    # df["old"] =["true"]*len(df)

    # Concatenate the original DataFrame with the NaN DataFrame, with NaN DataFrame inserted between each row
    interpolated_df = pd.concat([df, nan_df], axis=0).sort_index(kind='merge')

    interpolated_df  = interpolated_df.interpolate()

    # Filter rows where "id_ax" value is a round number
    interpolated_df = interpolated_df[interpolated_df['ax_id'] == interpolated_df['ax_id'].round()]


    interpolated_df["sph_id"] = list(map(lambda x: int(x*num), list(interpolated_df.sph_id)))
    
    interpolated_df["ax_id"] = list(map(lambda x: int(x), list(interpolated_df.ax_id)))
    interpolated_df["type"] = ["axon"]* len(interpolated_df)
    interpolated_df["P"] = list(map(lambda x: int(x*num)+nbr_spheres_to_add_inbetween, list(interpolated_df.P)))
    
    # Remove duplicate rows
    interpolated_df = interpolated_df.drop_duplicates()

    # Round all numerical values in the DataFrame to 4 digits
    interpolated_df = interpolated_df.round(6)

    return interpolated_df.reset_index(drop=True)

def interpolate_branch_glial_cell(df_processes_id_branch, branch, factor):

    num = factor
    nbr_spheres_to_add_inbetween = num -1
    # Create a DataFrame containing NaN values with the same columns as the original DataFrame
    nan_df = pd.DataFrame(np.nan, index=np.arange(len(df_processes_id_branch)*nbr_spheres_to_add_inbetween), columns=df_processes_id_branch.columns)
    nan_df["index"] = nan_df.index/3
    nan_df.index = nan_df["index"]
    # nan_df["old"] =["false"]*len(nan_df)      
    del nan_df["index"]
    # Concatenate the original DataFrame with the NaN DataFrame, with NaN DataFrame inserted between each row
    interpolated_df = pd.concat([df_processes_id_branch, nan_df], axis=0).sort_index(kind='merge')
    columns_to_convert = interpolated_df.columns.difference(['type'])
    interpolated_df[columns_to_convert] = interpolated_df[columns_to_convert].astype(float)
    interpolated_df  = interpolated_df.interpolate(method='linear')
    interpolated_df["sph_id"] = list(map(lambda x: int(x*num), list(interpolated_df.sph_id)))
    interpolated_df["P"] = [list(interpolated_df.P)[-1]]*len(interpolated_df)
    interpolated_df["branch_id"] = list(map(lambda x: int(branch), list(interpolated_df.branch_id)))
    interpolated_df = interpolated_df.iloc[1:]
    # replace nan in type with "glialRamification"
    interpolated_df["type"] = interpolated_df["type"].fillna("glialRamification")
    # round all numerical values in the DataFrame to 4 digits
    interpolated_df = interpolated_df.round(6)
    

    return interpolated_df


def interpolate_glial_cells(df):

    print("Glial Cells")
    columns_to_convert = df.columns.difference(['type'])
    df[columns_to_convert] = df[columns_to_convert].astype(float)
    df_somas = df.loc[df["type"] == "glial"]
    df_processes = df.loc[df["type"] == "glialRamification"]

    all_branches = []
    print("Branches")
    for _, soma in df_somas.iterrows():  # Iterate over soma cells
        id = int(soma["ax_id"])
        df_processes_id = df_processes.loc[(df_processes["ax_id"] == id) & (df_processes["P"] <= 0)]
        branches = df_processes_id["branch_id"].unique()
        for branch in branches:
            df_processes_id_branch = df_processes_id.loc[df_processes_id["branch_id"] == branch] 
            if (len(df_processes_id_branch) == 1):
                all_branches.append(df_processes_id_branch)
                continue
            # Add soma to the beginning of the dataframe
            df_processes_id_branch = pd.concat([soma.to_frame().T, df_processes_id_branch], ignore_index=True)
            interpolated_df = interpolate_branch_glial_cell(df_processes_id_branch, branch)
            all_branches.append(interpolated_df)

    
    df_processes_not_connected_to_soma = df_processes.loc[df_processes["P"] > 0]
    object_ids = df_processes_not_connected_to_soma["ax_id"].unique()
    print("Extra Branches")
    for ax_id in object_ids:
        df_processes_not_connected_to_soma_object = df_processes_not_connected_to_soma.loc[df_processes_not_connected_to_soma["ax_id"] == ax_id]
        #print(" ax id ", ax_id)
        #print(df_processes_not_connected_to_soma_object)
        parents = df_processes_not_connected_to_soma["P"].unique()
        
        for parent in parents:
            #print(" parent ", parent)
            # find parent sphere 
            parent_df = df_processes.loc[(df_processes["sph_id"]== parent) & (df_processes["ax_id"]== ax_id)]
            #print(" parent df ", parent_df)
            # find connected processes
            process_connected = df_processes_not_connected_to_soma_object.loc[df_processes_not_connected_to_soma_object["P"]== parent]
            #print(" connected ", process_connected)
            branches = process_connected["branch_id"].unique()
            for branch in branches:
                #print(" branch ", branch)
                process_connected_branch_ = process_connected.loc[process_connected["branch_id"] == branch]
        
                process_connected_branch = pd.concat([parent_df, process_connected_branch_], ignore_index=True)
                #print(" process connected branch ", process_connected_branch)
                interpolated_df = interpolate_branch_glial_cell(process_connected_branch, branch)
                #print(" interpolated ", interpolated_df)
                all_branches.append(interpolated_df)


    interpolated_df = pd.concat(all_branches, axis=0).sort_index(kind='merge')
    # round all numerical values in the DataFrame to 4 digits
    interpolated_df = interpolated_df.round(6)
    # merge with somas
    interpolated_df = pd.concat([df_somas, interpolated_df], axis=0).sort_index(kind='merge')
    interpolated_df = interpolated_df.drop_duplicates()
    return interpolated_df.reset_index(drop=True)


# Load original data from text file into DataFrame
folder= '/home/localadmin/Documents/MCDS/Permeable_MCDS/output/overlapping_factor'
original_df = read_swc_file(f'/{folder}/factor_1.swc')


#glial_df = interpolate_glial_cells(original_df)
factor = 8
axons_df = interpolate_axons_swc(original_df, factor)

# Create a new DataFrame containing the interpolated axons and glial cells
#interpolated_df = pd.concat([axons_df, glial_df], axis=0).sort_index(kind='merge')

interpolated_df = axons_df
# organise by "axi_id" and "sph_id"
interpolated_df = interpolated_df.sort_values(by=["type", "ax_id", "branch_id", "sph_id"])
print(interpolated_df)
# Save the new DataFrame to a text file
interpolated_df.to_csv(f'{folder}/growth_vox_100_factor_{factor}.swc', sep=' ', index=False, header=True)


    