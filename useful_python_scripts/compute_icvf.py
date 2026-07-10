import numpy as np
from tqdm import tqdm
from simulationgraphs import read_swc_file
import pandas as pd

def intra_volume(df, factor, bounds, compartment="intra"):
    R_col = 'inner_radius' if compartment == "intra" else 'outer_radius'
    
    # Helper function for boundary checking moved to the top
    def in_bounds(x, y, z, b):
        return (x >= b[0]) & (x <= b[1]) & (y >= b[2]) & (y <= b[3]) & (z >= b[4]) & (z <= b[5])
    
    # 1. Handle the downsampling factor
    if factor > 1:
        # Take every 'factor'-th row within each branch, plus the last row to close the segment
        df_sub = df.groupby(['cell_id', 'component', 'component_id']).apply(
            lambda x: pd.concat([x.iloc[::factor], x.iloc[[-1]]]) if len(x) % factor != 1 else x.iloc[::factor]
        ).drop_duplicates().reset_index(drop=True)
    else:
        df_sub = df.copy()

    group_cols = ['cell_id', 'component', 'component_id']
    
    # --- NEW: Catch Single-Sphere Components ---
    # Find components that only have exactly 1 point 
    group_sizes = df_sub.groupby(group_cols)['X'].transform('size')
    single_spheres_df = df_sub[group_sizes == 1]
    
    # 2. Align "Point 1" and "Point 2" on the same row using shift
    df_sub['X2'] = df_sub.groupby(group_cols)['X'].shift(-1)
    df_sub['Y2'] = df_sub.groupby(group_cols)['Y'].shift(-1)
    df_sub['Z2'] = df_sub.groupby(group_cols)['Z'].shift(-1)
    df_sub['R2'] = df_sub.groupby(group_cols)[R_col].shift(-1)
    
    # Drop the last point of every branch (it has no "next" point to form a segment with)
    segments = df_sub.dropna(subset=['X2']).copy()
    
    # --- SEGMENT VOLUME CALCULATION ---
    total_volume = 0.0
    if len(segments) > 0:
        X1, Y1, Z1, R1 = segments['X'].values, segments['Y'].values, segments['Z'].values, segments[R_col].values
        X2, Y2, Z2, R2 = segments['X2'].values, segments['Y2'].values, segments['Z2'].values, segments['R2'].values

        # Euclidean distance & Conical Frustum volume
        dist = np.sqrt((X1 - X2)**2 + (Y1 - Y2)**2 + (Z1 - Z2)**2)
        vol = np.pi * (R1**2 + R2**2 + R1*R2) * dist / 3.0

        # Boundary checks
        b1 = in_bounds(X1, Y1, Z1, bounds)
        b2 = in_bounds(X2, Y2, Z2, bounds)
        both_in = b1 & b2
        one_in = b1 ^ b2  

        # Sum segment volumes
        total_volume += np.sum(vol[both_in]) + np.sum(vol[one_in] / 2.0)
    
    # --- SINGLE SPHERE VOLUME CALCULATION ---
    if len(single_spheres_df) > 0:
        s_X = single_spheres_df['X'].values
        s_Y = single_spheres_df['Y'].values
        s_Z = single_spheres_df['Z'].values
        s_R = single_spheres_df[R_col].values
        
        # Check if these specific spheres fall within the bounding box
        s_in = in_bounds(s_X, s_Y, s_Z, bounds)
        
        # Calculate volume of a perfect sphere: V = 4/3 * pi * r^3
        s_vol = (4.0 / 3.0) * np.pi * (s_R**3)
        
        # Add the volumes of the single spheres that are inside the bounds
        total_volume += np.sum(s_vol[s_in])
    
    return total_volume

def compute_volume_fractions(square_bounds, file, factor):

    df =  read_swc_file(file)

    cell_types = df['cell_type'].unique()

    results_fraction = {"IAS": 0.0, "IAS+MYELIN":0.0, "EAS": 0.0, "CBV": 0.0, "IGS": 0.0}

    results_nbr_walkers = {"IAS": 0, "IAS+MYELIN":0, "EAS": 0, "CBV": 0, "IGS": 0}

    total_vol = (square_bounds[1] - square_bounds[0]) * (square_bounds[3] - square_bounds[2]) * (square_bounds[5] - square_bounds[4])

    print(f"Total volume of the bounding box: {total_vol:.2f} cubic microns")

    density = 1.0

    total_walkers = int(total_vol * density)

    for cell_type in cell_types:
        vol = 0.0
        print(f"Found cell type: {cell_type}")
        df_type = df.loc[df["cell_type"] == cell_type].copy()
        vol = intra_volume(df_type, factor, square_bounds, "intra")
        if cell_type == "axon":
            results_fraction["IAS"] += vol/total_vol
            results_nbr_walkers["IAS"] += int(vol * density)
            vol_total_axon = intra_volume(df_type, factor, square_bounds, "extra")
            results_fraction["IAS+MYELIN"] += vol_total_axon/total_vol
            results_nbr_walkers["IAS+MYELIN"] += int(vol_total_axon * density)
        elif cell_type == "glial_cell":
            results_fraction["IGS"] += vol/total_vol
            results_nbr_walkers["IGS"] += int(vol * density)
        elif cell_type == "blood_vessel":
            results_fraction["CBV"] += vol/total_vol
            results_nbr_walkers["CBV"] += int(vol * density)
    results_fraction["EAS"] = 1 - results_fraction["IAS+MYELIN"] - results_fraction["IGS"] - results_fraction["CBV"]
    results_nbr_walkers["EAS"] = total_walkers - results_nbr_walkers["IAS+MYELIN"] - results_nbr_walkers["IGS"] - results_nbr_walkers["CBV"]

    return results_fraction, results_nbr_walkers

def swap_inner_radius_outer_radius_columns(input_file, output_file):
    df =  read_swc_file(input_file)
    if 'inner_radius' in df.columns and 'outer_radius' in df.columns:
        df = df.rename(columns={'inner_radius': 'temp_radius'})
        df = df.rename(columns={'outer_radius': 'inner_radius'})
        df = df.rename(columns={'temp_radius': 'outer_radius'})
        # swap order of columns to match original format
        df = df[['cell_type', 'cell_id',  'component', 'component_id', 'X', 'Y', 'Z', 'inner_radius', 'outer_radius']]
        df.to_csv(output_file, sep=" ", index=False)
    else:
        print("Error: The input file must contain 'inner_radius' and 'outer_radius' columns.")

# Example usage
if __name__ == "__main__":

    # Path to the sphere file
    sphere_file = "/home/localadmin/Documents/Santi/Neuroinfl.csv"

    output_file = "/home/localadmin/Documents/Santi/Inflammation_swapped.csv"

    #swap_inner_radius_outer_radius_columns(sphere_file, output_file)

    # Define square bounds (xmin, xmax, ymin, ymax)
    square_bounds = (0, 100, 0, 100, 0, 100)

    factor = 4

    # Compute the intravolume fraction
    results_fraction, results_nbr_walkers = compute_volume_fractions(square_bounds, sphere_file, factor)
    print(results_nbr_walkers)
    print(results_fraction)