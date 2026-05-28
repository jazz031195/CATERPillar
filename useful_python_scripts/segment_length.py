import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import detrend
from simulationgraphs import read_swc_file
import matplotlib.pyplot as plt

import numpy as np
from scipy.optimize import curve_fit

def segment_length(df):

    angle_threshold =70  # degrees;

    segment_lengths = []

    for (cell_id, component_id), vessel in df.groupby(["cell_id", "component_id"]):

        coords = vessel[["X", "Y", "Z"]].values

        # vectors between consecutive points
        v = np.diff(coords, axis=0)   

        # angle between consecutive vectors
        v1 = v[:-1]
        v2 = v[1:]

        cosang = np.sum(v1 * v2, axis=1) / (
            np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1)
        )
        cosang = np.clip(cosang, -1, 1)

        angles = np.degrees(np.arccos(cosang))

        # keep endpoints + strong turning points
        turning_idx = np.where(angles > angle_threshold)[0] +1
        #keep_idx = np.r_[0, turning_idx, len(coords)-1]
        key_coords = coords[turning_idx]

        # distance between main turning points
        main_lengths = np.linalg.norm(np.diff(key_coords, axis=0), axis=1)

        # exclude below 12
        main_lengths = main_lengths[main_lengths > 12]

        segment_lengths.extend(main_lengths)
    plt.figure(figsize=(8, 5))
    plt.hist(segment_lengths, bins=50, edgecolor='black')
    plt.title('Distribution of Segment Lengths')
    plt.xlabel('Segment Length (microns)')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.show()
    average_branch_length = np.mean(segment_lengths)
    return average_branch_length


if __name__ == "__main__":
    # Example usage with dummy data
    # Replace this with actual 3D points from an axon skeleton
    path = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/incoherent_blood_flow/State1.csv"
    cells = "blood vessel"
    df =  read_swc_file(path)
    if cells == "axon":
        df = df.loc[df["cell_type"] == "axon"]
    elif cells == "glial cell":
        df = df.loc[df["cell_type"] == "glial_cell"]
    elif cells == "blood vessel":
        df = df.loc[df["cell_type"] == "blood_vessel"]

    segment = segment_length(df)
    print(f"Average segment length: {segment:.2f} microns")