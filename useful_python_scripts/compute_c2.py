import numpy as np
from simulationgraphs import get_spheres_array, read_swc_file
import matplotlib.pyplot as plt

def cos_angle(v1, v2):
    """
    Compute the cosine of the angle between two 3D vectors.

    Parameters:
        v1 (np.ndarray): First vector (3D position).
        v2 (np.ndarray): Second vector (3D position).

    Returns:
        float: Cosine of the angle between the two vectors.
    """
    # Ensure inputs are numpy arrays
    v1 = np.array(v1)
    v2 = np.array(v2)

    # Compute the dot product
    dot_product = np.dot(v1, v2)
    
    # Compute the magnitudes (Euclidean norms)
    magnitude_v1 = np.linalg.norm(v1)
    magnitude_v2 = np.linalg.norm(v2)

    # Avoid division by zero
    if magnitude_v1 == 0 or magnitude_v2 == 0:
        raise ValueError("One of the vectors has zero magnitude; cannot compute angle.")

    # Compute cosine of the angle
    cos_theta = dot_product / (magnitude_v1 * magnitude_v2)

    return cos_theta

def compute_c2_length(axon_df):

    first_sphere = axon_df.iloc[0]
    last_sphere = axon_df.iloc[-1]
    z_axis = np.array([0, 0, 1])
    v1 = np.array([last_sphere['x'] - first_sphere['x'], last_sphere['y'] - first_sphere['y'], last_sphere['z'] - first_sphere['z']])
    v1 = v1 / np.linalg.norm(v1)
    cos2 = np.dot(v1, z_axis)**2 
    return cos2


def compute_c2(axon_df, factor, limit=16):
    """
    Compute the C2 metric for a given axon.

    Parameters:
        axon_df (pd.DataFrame): DataFrame containing the axon data with 'x', 'y', 'z' columns.
        factor (int): Step size for sampling spheres along the axon.
        limit (float): Distance threshold for valid sphere pairs (default: 15).

    Returns:
        float: The computed C2 metric value.
    """
    c2_values = []
    z_axis = np.array([0, 0, 1])

    # Extract position data for easier access
    positions = axon_df[['x', 'y', 'z']].values

    # Iterate over sampled spheres based on the 'factor'
    for i in range(factor, len(positions), factor):

        sphere1 = positions[i - factor]
        distances = np.linalg.norm(positions[i:] - sphere1, axis=1)

        # Find the first valid sphere (distance > limit)
        valid_indices = np.where(distances > limit)[0]
        if not len(valid_indices):
            continue
        
        # Take the first valid sphere (short-circuit)
        sphere2 = positions[i + valid_indices[0]]
        
        # Compute cosine squared of the angle between vectors
        v1 = sphere2 - sphere1
        cos_theta_sq = np.dot(v1, z_axis)**2 / (np.linalg.norm(v1)**2 * np.linalg.norm(z_axis)**2)
        c2_values.append(cos_theta_sq)

    # Compute the average C2 value
    if len(c2_values) > 0:
        return np.mean(c2_values)
    return np.nan



def compute_c2_axons(df, factor, limit):

    axons_ids = df['ax_id'].unique()
    c2_tot = 0
    tot_nbr = 0
    for axon_id in axons_ids:
        # print(f"Computing C2 for axon {axon_id}")
        axon_df = df[df['ax_id'] == axon_id]
        cos2= compute_c2(axon_df, factor, limit)
        print("C2 : ", cos2)
        if cos2 is not np.nan:
            c2_tot += cos2
            tot_nbr += 1

    return c2_tot / tot_nbr

def plot_angle_histogram(df, factor, limit):
    axons_ids = df['ax_id'].unique()
    all_c2 = []
    for axon_id in axons_ids:
        print(f"Computing C2 for axon {axon_id}")
        axon_df = df[df['ax_id'] == axon_id]
        cos2= compute_c2(axon_df, factor, limit)
        print("C2 : ", cos2)
        all_c2.append(cos2)
    
    #histogram
    plt.hist(all_c2, bins=20)
    plt.xlabel('C2')
    plt.ylabel('Frequency')
    plt.title('C2 distribution')
    plt.show()

# Example usage
if __name__ == "__main__":

    # Path to the sphere file
    sphere_file = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/SMI_pred/axons_astrocytes/astrocytes_0.0.swc"

    factor = 4
    limit = 30

    # Load sphere data
    df = read_swc_file(sphere_file)

    # Compute the C2 metric
    c2 = compute_c2_axons(df, factor, limit)
    print(f"C2 Metric: {c2}")
