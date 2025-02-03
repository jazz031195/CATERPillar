import numpy as np
import os
from compute_icvf import compute_intravolume_fraction
import pandas as pd
import sys
import seaborn as sns
import matplotlib.pyplot as plt
sys.path.append('/home/localadmin/Documents/MCDS/Permeable_MCDS/useful_functions')
from useful_functions import get_scheme_info_iso, get_scheme_info, read_binary_file

def calculate_DWI(file, compartment, scheme_path):

    if "waveform" in scheme_path:      
        b_values, directions = get_scheme_info_iso(scheme_path)
    else:
        b_values, directions = get_scheme_info(scheme_path)

    DWIs = read_binary_file(file)

    DWIs =[float(i) for i in DWIs] 
    # calculate angle between direction and vector (0,0,1)
    angles = [ np.arccos(np.dot(x, [0,0,1])/np.linalg.norm(x)) for x in directions]
    # angle must be from 0 to pi
    angles = [ x if x <= np.pi/2 else np.pi-x for x in angles]

    df = pd.DataFrame(columns=['angle (rad)', 'b_value', 'DWI'])
    df["DWI"] = DWIs
    df["angle (rad)"] = angles
    df["b_value"] = b_values
    df["compartment"] =[compartment]*len(df)
    print(df.loc[df["b_value"] == 1])
    print(df.loc[df["b_value"] == 0.2])
    print(df)


    return df

def ADC(b0, b1, df):

    df_b0 = df.loc[df["b_value"] == b0]
    df_b1 = df.loc[df["b_value"] == b1]

    print(df_b0)
    print(df_b1)
    new_df = pd.DataFrame(columns=['angle (rad)', 'ADC'])
    angles = df_b0["angle (rad)"].unique()
    all_angles = []
    all_adcs = []
    for angle in angles:
        adc = np.log(df_b0.loc[df_b0["angle (rad)"] == angle, "DWI"].values/df_b1.loc[df_b1["angle (rad)"] == angle, "DWI"].values)/(b1-b0)
        all_angles.append(angle)
        all_adcs.append(adc[0])
    new_df["angle (rad)"] = all_angles
    new_df["ADC"] = all_adcs
    return new_df

def calculate_ADC(b0, b1, file, compartment, scheme_path):
    
    df = calculate_DWI(file, compartment, scheme_path)
    df_b0 = ADC(b0, b1, df)

    return df_b0

def compute_all_icvfs(directory):

    swc_files = os.listdir(directory)  
    swc_filenames =[f for f in swc_files if "swc" in f]
    swc_files = [f"{directory}/{f}" for f in swc_filenames]

    bounds = (0, 100, 0, 100, 0, 100)
    factor = 4

    df = pd.DataFrame(columns=['filename', 'icvf'])
    for i, swc_file in enumerate(swc_files):
        print(f"Computing ICVF for {swc_file}")
        icvf = compute_intravolume_fraction(bounds, swc_file, factor)
        print(f"ICVF: {icvf}")
        filename = swc_filenames[i].split(".swc")[0]
        df = df._append({'filename': filename, 'icvf': icvf}, ignore_index=True)

    df.to_csv(f"{directory}/icvf.csv", index=False)


def compute_df_no_compartments(directory):

    icvf_df = pd.read_csv(f"{directory}/icvf.csv")
    all_dfs= [] 

    folder = f"{directory}/"
    # find all files in the folder
    DWI_files = os.listdir(folder)  
    DWI_filenames =[f for f in DWI_files if "DWI" in f and "img" not in f] 
    DWI_files = [f"{folder}/{f}" for f in DWI_filenames]

    for i, file in enumerate(DWI_files):
        print(f"Calculating ADC for {file}")
        if "iso" in file:
            scheme_path = "/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/iso_waveform_vec_b.txt"
            sequence = "iso"
        else:
            scheme_path = "/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/PGSE_21_dir_12_b.scheme"
            sequence = "pgse"
        filename = DWI_filenames[i]
        df = calculate_DWI(file, "", scheme_path).copy()
        print (df)

        df["sequence"] = [sequence]*len(df)
        df["swelling"] = [filename.split("voxel_")[1] ]*len(df)
        all_dfs.append(df)
    
    df_total = pd.concat(all_dfs)
    # add DWI of df with same direction, angle, nbr,  sequence but different compartment 
    df_total = df_total.groupby(["swelling", 'angle (rad)', 'b_value', "sequence"], as_index=False)['DWI'].sum()
    df_total["compartment"] = ["all"]*len(df_total)
    df_total = ADC(0.2, 1, df_total)
    return df_total



            
def compute_df_compartment(directory, compartment):

    folder = f"{directory}/{compartment}"
    # find all files in the folder
    DWI_files = os.listdir(folder)  
    DWI_filenames =[f for f in DWI_files if "DWI" in f and "img" not in f] 
    DWI_files = [f"{folder}/{f}" for f in DWI_filenames]

    all_dfs= [] 
    for i, file in enumerate(DWI_files):
        print(f"Calculating ADC for {file}")
        if "iso" in file:
            scheme_path = "/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/iso_waveform_vec_b.txt"
            sequence = "iso"
        
        else:
            scheme_path = "/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/PGSE_21_dir_12_b.scheme"
            sequence = "pgse"
            
        
        b0 = 0.2
        b1 = 1
        df = calculate_ADC(b0, b1, file, compartment, scheme_path)
        filename = DWI_filenames[i].split(compartment)[0]
        df["filename"] = [filename]*len(df)
        df["sequence"] = [sequence]*len(df)
        all_dfs.append(df)

    return pd.concat(all_dfs)

def calculate_relative(group):
    base_adc = group.loc[group['swelling'] == "0", 'ADC'].values[0]
    group['ADC_relative_change'] = (group['ADC']-base_adc)*100 / base_adc
    return group


def analysis(directory):
    compartment = "extra"
    df = compute_df_compartment(directory, compartment)
    
    df["swelling"] = df["filename"].apply(lambda x: x.split("_")[1])

    df_pgse = df.loc[df["sequence"] == "pgse"]
    df_pgse = df_pgse.groupby('angle (rad)').apply(calculate_relative)
    df_pgse.to_csv(f"{directory}/{compartment}_pgse.csv", index=False) 
    sns.lmplot(x='angle (rad)', y='ADC_relative_change', data=df_pgse, hue='swelling')
    plt.title(f"{compartment} ADC relative change, PGSE")
    plt.show()

    sns.lmplot(x='angle (rad)', y='ADC', data=df_pgse, hue='swelling')
    plt.title(f"{compartment} ADC , linear")
    plt.show()


    df_iso = df.loc[df["sequence"] == "iso"]
    df_iso = df_iso.groupby('angle (rad)').apply(calculate_relative)
    df_iso.to_csv(f"{directory}/{compartment}_iso.csv", index=False)
    sns.lmplot(x='angle (rad)', y='ADC_relative_change', data=df_iso, hue='swelling')
    plt.title(f"{compartment} ADC relative change, ISO")
    plt.show()

    sns.lmplot(x='angle (rad)', y='ADC', data=df_iso, hue='swelling')
    plt.title(f"{compartment} ADC , ISO")
    plt.show()
    print(df_iso)


    

def all_analysis(directory):
    df = compute_df_no_compartments(directory)
    print(df)
    df_pgse = df.loc[df["sequence"] == "pgse"]
    df_pgse = df_pgse.groupby('angle (rad)').apply(calculate_relative)
    # to csv
    df_pgse.to_csv(f"{directory}/all_pgse.csv", index=False)
    df_iso = df.loc[df["sequence"] == "iso"]
    df_iso = df_iso.groupby('angle (rad)').apply(calculate_relative)
    # to csv
    df_iso.to_csv(f"{directory}/all_iso.csv", index=False)

    # plot
    sns.lmplot(x='angle (rad)', y='ADC_relative_change', data=df_iso, hue='swelling')
    plt.title("ADC relative change, ISO")
    plt.show()

    sns.lmplot(x='angle (rad)', y='ADC_relative_change', data=df_pgse, hue='swelling')
    plt.title("ADC relative change, PGSE")
    plt.show()

    sns.lmplot(x='angle (rad)', y='ADC', data=df_iso, hue='swelling')
    plt.title("ADC, ISO")
    plt.show()

    sns.lmplot(x='angle (rad)', y='ADC', data=df_pgse, hue='swelling')
    plt.title("ADC, PGSE")
    plt.show()

if __name__ == '__main__':

    directory = "/home/localadmin/Documents/CATERPillar/arthurs_analysis"
    #directory = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/verif_arthurs_analysis"
    analysis(directory)