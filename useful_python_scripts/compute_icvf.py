import numpy as np
from tqdm import tqdm
from simulationgraphs import read_swc_file


def sphere_volume(radius):
    return (4.0 / 3.0) * np.pi * radius**3

def overlap_volume(r1, r2, d):
    if d >= r1 + r2:
        return 0.0  # No overlap

    part1 = (r1 + r2 - d)**2
    part2 = d**2 + 2 * d * (r1 + r2) - 3 * (r1 - r2)**2
    return (np.pi * part1 * part2) / (12.0 * d)

def sphere_within_bounds(sphere, bounds):
    return (bounds[0] <= sphere['X'] <= bounds[1] and
            bounds[2] <= sphere['Y'] <= bounds[3] and
            bounds[4] <= sphere['Z'] <= bounds[5])

def volume_cylinders(axon_df, factor, bounds, compartment="intra"):
    R_col = 'inner_radius' if compartment == "intra" else 'outer_radius'
    volume = 0.0

    for i in np.arange(0, len(axon_df), factor):
        sphere1 = axon_df.iloc[i]
        R1 = sphere1[R_col] 
        X1, Y1, Z1 = sphere1['X'], sphere1['Y'], sphere1['Z']
        sphere2 = axon_df.iloc[i+factor] if (i+factor) < len(axon_df) else axon_df.iloc[-1]
        R2 = sphere2[R_col]
        X2, Y2, Z2 = sphere2['X'], sphere2['Y'], sphere2['Z']
        if (sphere_within_bounds(sphere1, bounds) and sphere_within_bounds(sphere2, bounds)):
            volume += np.pi * (R1**2 + R2**2 + R2*R1) * np.linalg.norm(np.array([X1,Y1,Z1] ) - np.array([X2,Y2,Z2] )) / 3
        elif (sphere_within_bounds(sphere1, bounds) or sphere_within_bounds(sphere2, bounds)):
            volume += np.pi * (R1**2 + R2**2 + R2*R1) * np.linalg.norm(np.array([X1,Y1,Z1] ) - np.array([X2,Y2,Z2] )) / 6
    return volume

def intra_volume(df, factor, bounds, compartment):
    axons_ids = df['cell_id'].unique()
    volume_tot = 0
    for axon_id in axons_ids:
        print(f"Computing volume for cell {axon_id}")
        axon_df = df[df['cell_id'] == axon_id]
        volume_tot += volume_cylinders(axon_df, factor, bounds, compartment)

    return volume_tot

def compute_intravolume_fraction(square_bounds, file, factor, compartment, cells = None):

    df =  read_swc_file(file)
    if cells == "axon":
        df = df.loc[df["cell_type"] == "axon"]
    elif cells == "glial cell":
        df = df.loc[df["cell_type"] == "glial_cell"]
    elif cells == "blood vessel":
        print(df)
        df = df.loc[df["cell_type"] == "blood_vessel"]

    intra_vol = intra_volume(df, factor, square_bounds, compartment)

    total_vol = (square_bounds[1] - square_bounds[0]) * (square_bounds[3] - square_bounds[2]) * (square_bounds[5] - square_bounds[4])
    
    return intra_vol / total_vol

# Example usage
if __name__ == "__main__":

    # Path to the sphere file
    sphere_file = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/incoherent_blood_flow/State4.csv"

    # Define square bounds (xmin, xmax, ymin, ymax)
    square_bounds = (0, 150, 0, 150, 0, 150)

    factor = 4

    compartment ="intra"

    # Compute the intravolume fraction
    icvf = compute_intravolume_fraction(square_bounds, sphere_file, factor, compartment, cells = "axon")
    print(f"Intravolume Fraction: {icvf}")