import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
from pathlib import Path
from scipy import stats
import random
from simulationgraphs import get_spheres_array, read_swc_file
from scipy.optimize import least_squares
from reduce_astrocytes_in_substrate import write_new_swc

def compute_icvf_axons(df_axon, limit,factor):

    total_volume = 0
    for i in np.arange(factor,len(df_axon),factor):
        sphere1 = df_axon.iloc[i-factor]
        sphere2 = df_axon.iloc[i]
        distance_between_spheres = np.linalg.norm(
            np.array([sphere1.x, sphere1.y, sphere1.z]) - np.array([sphere2.x, sphere2.y, sphere2.z])
        )
        volume = np.pi * (sphere1.Rout**2 + sphere2.Rout**2 + sphere2.Rout * sphere1.Rout) * distance_between_spheres / 3
        total_volume += volume

    return total_volume/(limit**3)


def myelin_thickness(inner_radius):
    """
    Compute the myelin thickness based on the inner radius.
    """
    return 0.35 + 0.006 * 2.0 * inner_radius + 0.024 * np.log(2.0 * inner_radius)

def compute_outer_radius(inner_radius):
    """
    Compute the outer radius based on the inner radius.
    """
    return inner_radius + myelin_thickness(inner_radius)

def compute_inner_radius(outer_radius):
    """
    Compute the inner radius given an outer radius using least squares optimization.

    Parameters:
        outer_radius (float): The outer radius value.

    Returns:
        float: The computed inner radius.
    """
    def residual(inner_radius):
        # Residual function: difference between given and computed outer radius
        return compute_outer_radius(inner_radius) - outer_radius

    # Initial guess for the inner radius (can be tuned based on the expected range)
    initial_guess = outer_radius / 2.0

    # Solve using least squares
    result = least_squares(residual, initial_guess, bounds=(0, outer_radius))  # Non-negative inner radius

    if result.success:
        return result.x[0]
    else:
        raise ValueError("Optimization failed to converge.")
    
def add_myelin_to_one_axon(df_axon):

    df_axon["Rin"] = list(map(lambda x: compute_inner_radius(x), list(df_axon["Rout"])))
    
    return df_axon

def add_myelin(f_myelinated_axons_to_reach, path_swc, factor):
    f_myelinated_axons_reached = 0

    df = read_swc_file(path_swc)
    df_axons = df[df['type'] == "axon"]
    limit = 150

    new_df = pd.DataFrame(columns=df.columns)
    ax_ids = df_axons['ax_id'].unique()
    #shuffle ax_ids
    random.shuffle(ax_ids)
    for axon_id in ax_ids:
        print(f"Processing axon {axon_id}")
        df_axon = df_axons[df_axons['ax_id'] == axon_id]
        if (f_myelinated_axons_reached >= f_myelinated_axons_to_reach):
            new_df = pd.concat([new_df, df_axon])
        else:
            icvf_axon = compute_icvf_axons(df_axon, limit,factor)
            if (icvf_axon <= (f_myelinated_axons_to_reach-f_myelinated_axons_reached)+0.05):
                print(f"myelinated")
                df_axon = add_myelin_to_one_axon(df_axon)
                f_myelinated_axons_reached += icvf_axon
            new_df = pd.concat([new_df, df_axon])

    print(f"Reached ICVF: {f_myelinated_axons_reached}")
    return new_df

def compute_total_icvf(path_swc, factor):
    df = read_swc_file(path_swc)
    df_axons = df[df['type'] == "axon"]
    limit = 150

    ax_ids = df_axons['ax_id'].unique()
    total_icvf = 0
    for axon_id in ax_ids:
        print(f"Processing axon {axon_id}")
        df_axon = df_axons[df_axons['ax_id'] == axon_id]
        icvf_axon = compute_icvf_axons(df_axon, limit,factor)
        total_icvf += icvf_axon

    return total_icvf

if __name__ == "__main__":
    path_swc = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/SMI_pred/complex_axons/f_0.2.swc"
    f_myelinated_axons_to_reach = 0.1
    factor = 4
    total_icvf = compute_total_icvf(path_swc, factor)
    print   (f"Total ICVF: {total_icvf}")
    #new_df = add_myelin(f_myelinated_axons_to_reach, path_swc,factor)
    #new_path_swc = f"/home/localadmin/Documents/MCDS/Permeable_MCDS/output/SMI_pred/permeable_axons/axons_{f_myelinated_axons_to_reach}.swc"
    #write_new_swc(new_path_swc, new_df)