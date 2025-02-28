import numpy as np
from tqdm import tqdm
from simulationgraphs import read_swc_file


def sphere_within_bounds(sphere, bounds):
    return (bounds[0] <= sphere['x'] <= bounds[1] and
            bounds[2] <= sphere['y'] <= bounds[3] and
            bounds[4] <= sphere['z'] <= bounds[5])

def volume_cylinders(axon_df, factor, bounds, compartment = "intra"):

    if compartment == "intra":
        R = 'Rin'
    else:
        R = 'Rout'
    volume = 0
    for i in np.arange(factor, len(axon_df), factor):
        sphere1 = axon_df.iloc[i-factor]
        sphere2 = axon_df.iloc[i]
        
        if (sphere_within_bounds(sphere1, bounds) and sphere_within_bounds(sphere2, bounds)):
            volume += np.pi * (sphere1[R]**2 + sphere2[R]**2 + sphere2[R]*sphere1[R]) * np.linalg.norm(sphere1[['x', 'y', 'z']] - sphere2[['x', 'y', 'z']]) / 3
        elif (sphere_within_bounds(sphere1, bounds) or sphere_within_bounds(sphere2, bounds)):
            volume += np.pi * (sphere1[R]**2 + sphere2[R]**2 + sphere2[R]*sphere1[R]) * np.linalg.norm(sphere1[['x', 'y', 'z']] - sphere2[['x', 'y', 'z']]) / 6
    return volume

def intra_volume(df, factor, bounds, compartment):
    axons_ids = df['id_ax'].unique()
    volume_tot = 0
    for axon_id in axons_ids:
        print(f"Computing volume for axon {axon_id}")
        axon_df = df[df['id_ax'] == axon_id]
        volume_tot += volume_cylinders(axon_df, factor, bounds, compartment)

    return volume_tot

def compute_intravolume_fraction(square_bounds, file, factor, compartment, cells = None):

    df =  read_swc_file(file)
    if cells == "axons":
        df = df.loc[df["type"] == "axon"]
    elif cells == "astrocytes":
        df = df.loc[df["type"] != "axon"]

    intra_vol = intra_volume(df, factor, square_bounds, compartment)

    total_vol = (square_bounds[1] - square_bounds[0]) * (square_bounds[3] - square_bounds[2]) * (square_bounds[5] - square_bounds[4])
    
    return intra_vol / total_vol

# Example usage
if __name__ == "__main__":

    # Path to the sphere file
    sphere_file = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/SMI_pred/axons_astrocytes/astrocytes_0.06.swc"

    # Define square bounds (xmin, xmax, ymin, ymax)
    square_bounds = (0, 150, 0, 150, 0, 150)

    factor = 4

    compartment ="intra"

    # Compute the intravolume fraction
    icvf = compute_intravolume_fraction(square_bounds, sphere_file, factor, compartment, cells = "axons")
    print(f"Intravolume Fraction: {icvf}")