import numpy as np
from simulationgraphs import get_spheres_array, read_swc_file
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from sklearn.decomposition import PCA
from scipy.signal import periodogram

import numpy as np
from scipy.fft import fft, fftfreq
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

def estimate_undulation_fft(points_3d, plot = False):
    """
    Estimate the dominant wavelength and amplitude of axonal undulation
    using FFT after projecting the skeleton onto its main axis.

    Parameters:
        points_3d (ndarray): Nx3 array of 3D points (x, y, z) along the axon.

    Returns:
        dominant_wavelength (float): Dominant undulation wavelength.
        undulation_amplitude (float): Root-mean-square undulation amplitude.
    """

    # Step 1: PCA to align axon with main axis
    pca = PCA(n_components=3)
    pca.fit(points_3d)
    transformed = pca.transform(points_3d)

    # Use the main axis as the new 'z'
    z_aligned = transformed[:, 0]
    lateral_displacement = np.sqrt(transformed[:, 1]**2 + transformed[:, 2]**2)

    # Step 2: Sort by z and resample to uniform spacing
    sort_idx = np.argsort(z_aligned)
    z_sorted = z_aligned[sort_idx]
    w_sorted = lateral_displacement[sort_idx]

    # Interpolate to uniform grid
    z_uniform = np.linspace(z_sorted.min(), z_sorted.max(), len(z_sorted))
    w_uniform = np.interp(z_uniform, z_sorted, w_sorted)

    # Step 3: FFT
    n = len(z_uniform)
    dz = np.mean(np.diff(z_uniform))
    freqs = fftfreq(n, d=dz)
    fft_vals = fft(w_uniform)
    power = np.abs(fft_vals)**2

    # Use only positive frequencies
    mask = freqs > 0
    freqs = freqs[mask]
    power = power[mask]

    # Find dominant frequency
    dominant_idx = np.argmax(power)
    dominant_freq = freqs[dominant_idx]
    dominant_wavelength = 1 / dominant_freq

    # Amplitude estimate (RMS)
    undulation_amplitude = np.std(w_uniform)

    if plot:

        # Optional plot for diagnostics
        plt.figure(figsize=(10, 4))
        plt.plot(z_uniform, w_uniform)
        plt.xlabel("Distance along axon (μm)")
        plt.ylabel("Lateral displacement (μm)")
        plt.title("Undulation Profile (after PCA)")
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    return dominant_wavelength, undulation_amplitude


import numpy as np

def compute_axon_tortuosity(df):
    """
    Calculates the Arc-Chord tortuosity ratio for each axon.
    Tortuosity = (Total Arc Length) / (Distance between endpoints)
    
    Parameters:
    df (pd.DataFrame): DataFrame containing 'cell_id', 'component_id', 'X', 'Y', 'Z'
    
    Returns:
    dict: Mapping of (cell_id, component_id) to its individual tortuosity value.
    float: The average tortuosity across all valid axons.
    """
    tortuosity_results = {}
    valid_tortuosity_values = []
    
    for (cell_id, component_id), axon in df.groupby(["cell_id", "component_id"]):
        coords = axon[["X", "Y", "Z"]].values
        
        # An axon must have at least 2 points to have a length
        if len(coords) < 2:
            continue
            
        # 1. Calculate Arc Length (the actual winding path)
        # np.diff gets the vectors between consecutive points, norm gets their lengths
        arc_length = np.sum(np.linalg.norm(np.diff(coords, axis=0), axis=1))
        
        # 2. Calculate Chord Length (straight-line distance from first to last point)
        chord_length = np.linalg.norm(coords[-1] - coords[0])
        
        # 3. Calculate Tortuosity (must be >= 1.0)
        if chord_length > 0:
            tortuosity = arc_length / chord_length
            tortuosity_results[(cell_id, component_id)] = tortuosity
            valid_tortuosity_values.append(tortuosity)
        else:
            # If start and end points are exactly the same (a perfect loop), 
            # the ratio is mathematically undefined.
            tortuosity_results[(cell_id, component_id)] = np.inf
            
    # Calculate the mean across all valid axons
    mean_tortuosity = np.mean(valid_tortuosity_values) if valid_tortuosity_values else np.nan

    print(f"Computed tortuosity for {len(tortuosity_results)} axons. Average tortuosity: {mean_tortuosity:.3f}")
    
    return tortuosity_results, mean_tortuosity

if __name__ == "__main__":

    # Path to the sphere file
    sphere_file = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/incoherent_blood_flow/State1.csv"

    factor = 4
    limit = 150

    # Load sphere data
    df = read_swc_file(sphere_file)

    df = df[df['cell_type'] == 'blood_vessel']  # Filter for axons only

    compute_axon_tortuosity(df)

    assert(0)

    unique_cells = df['cell_id'].unique()

    amplitudes = []
    wavelengths = []

    for cell_id in unique_cells:
        print(f"Processing cell ID: {cell_id}/ {len(unique_cells)}")
        df_cell = df[df['cell_id'] == cell_id]
        skeleton_points = df_cell[['X', 'Y', 'Z']].values 
        dominant_wavelength, undulation_amplitude = estimate_undulation_fft(skeleton_points, plot = False)
        amplitudes.append(undulation_amplitude)
        wavelengths.append(dominant_wavelength)

    #plot 
    plt.figure(figsize=(12, 6))
    plt.plot(wavelengths, amplitudes, 'o', markersize=5)
    plt.xlabel("Wavelength (μm)")
    plt.ylabel("Amplitude (μm)")
    plt.title("Undulation Amplitude vs Wavelength")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

