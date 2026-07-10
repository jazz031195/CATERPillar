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
from compute_icvf import intra_volume

def reduce_myelin_volume(input_path_swc, output_path, reduction_fraction=0.5):
    """
    Reduces the myelin volume of axon segments by a specified fraction 
    by adjusting only the outer radius.
    
    Parameters:
    df (pd.DataFrame): DataFrame containing 'inner_radius', 'outer_radius', and 'cell_id'
    reduction_fraction (float): The fraction by which to reduce the volume (default 0.5 for 50%)
    
    Returns:
    pd.DataFrame: A copy of the DataFrame with updated 'outer_radius' values.
    """
    # Create a copy to avoid SettingWithCopy warnings on the original dataframe
    df = read_swc_file(input_path_swc)
    df_axons = df[df['cell_type'] == "axon"]
    
    # Extract the radii as numpy arrays for clean calculations
    r_in = df_axons['inner_radius'].values
    r_out_old = df_axons['outer_radius'].values
    
    # Calculate the new outer radius using the derived geometric formula
    # Factor is the amount of volume to KEEP (e.g., 1.0 - 0.50 = 0.50)
    retention_factor = 1.0 - reduction_fraction
    
    r_out_new_squared = (retention_factor * (r_out_old**2 - r_in**2)) + r_in**2
    r_out_new = np.sqrt(r_out_new_squared)
    
    # Update the dataframe
    df_axons['outer_radius'] = r_out_new

    #save the modified dataframe to a new SWC file
    pd.concat([ df_axons, df[df['cell_type'] != "axon"]]).to_csv(output_path, index=False, sep=" ")
    
    print(f"Myelin volume reduced by {reduction_fraction*100:.1f}%. Updated SWC saved to {output_path}")



def simulate_axon_loss(input_path_swc, output_path_swc, volume_func, factor, square_bounds, 
                       cell_type="axon", reduction_fraction=0.3):
    """
    Randomly removes entire cells of a specific type until the total 
    intracellular volume is reduced by the specified fraction.
    """
    df = read_swc_file(input_path_swc)
    
    # 1. Isolate the specific cell type to calculate the baseline
    df_type = df[df["cell_type"] == cell_type].copy()
    
    # 2. Calculate the initial total volume and our target reduction
    total_vol = volume_func(df_type, factor, square_bounds, "intra")
    target_reduction_vol = total_vol * reduction_fraction
    
    # 3. Get unique cell IDs and shuffle them for random deletion
    cell_ids = df_type["cell_id"].unique()
    np.random.shuffle(cell_ids)
    
    removed_vol = 0.0
    ids_to_remove = []
    
    # 4. Iteratively select cells to remove until the target volume is met
    for cid in cell_ids:
        if removed_vol >= target_reduction_vol:
            break
            
        df_single_cell = df_type[df_type["cell_id"] == cid]
        cell_vol = volume_func(df_single_cell, factor, square_bounds, "intra")
        
        removed_vol += cell_vol
        ids_to_remove.append(cid)
        
    # 5. Clean filtering of the original dataframe
    # Find rows that are BOTH the target cell_type AND in the removal list
    mask_to_drop = (df["cell_type"] == cell_type) & (df["cell_id"].isin(ids_to_remove))
    
    # Keep everything that is NOT in the drop mask
    df_final = df[~mask_to_drop].copy()
    
    # Print a quick summary
    actual_reduction_pct = (removed_vol / total_vol) * 100 if total_vol > 0 else 0
    print(f"--- {cell_type.capitalize()} Loss Summary ---")
    print(f"Initial Volume:   {total_vol:.2f}")
    print(f"Target Reduction: {target_reduction_vol:.2f} ({(reduction_fraction*100):.1f}%)")
    print(f"Actual Reduction: {removed_vol:.2f} ({actual_reduction_pct:.1f}%)")
    print(f"Cells Removed:    {len(ids_to_remove)} out of {len(cell_ids)} {cell_type}s")

    # 6. Save the modified dataframe to a new SWC file
    df_final.to_csv(output_path_swc, index=False, sep=" ")
    print(f"Updated SWC with reduced {cell_type} volume saved to {output_path_swc}")
    
    return df_final # Returning it just in case you need to use it immediately

if __name__ == "__main__":
    path_swc = "/home/localadmin/Documents/Santi/Inflammation.csv"
    output_path1 = "/home/localadmin/Documents/Santi/Inflammation_reduced_myelin.csv"
    factor = 4
    square_bounds = (0, 100, 0, 100, 0, 100)
    reduce_myelin_volume(path_swc, output_path1, reduction_fraction=0.5)
    output_path2 = "/home/localadmin/Documents/Santi/Neuroinfl.csv"
    simulate_axon_loss(output_path1, output_path2, intra_volume, factor, square_bounds, cell_type="axon", reduction_fraction=0.3)

