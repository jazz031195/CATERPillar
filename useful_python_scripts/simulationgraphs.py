import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.colors as colors
import os
import glob
from matplotlib.patches import Circle
import matplotlib.colors as mcolors
import math
import random
import pyvista as pv
from tqdm import tqdm
import pandas as pd
from matplotlib.patches import Circle
from multiprocessing import Pool
import copy
from scipy.signal import find_peaks
from scipy.interpolate import interp1d
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import gaussian_kde
     

def read_swc_file(file_path):

    if "swc" in file_path:
        print(file_path, " is old format")
        columns = ["id_ax","sph_id", "branch_id", "type", "x", "y", "z", "Rin","Rout", "P"]
        df = pd.read_csv(file_path, sep=' ', names=columns)

        df = df.iloc[1:]
        df["x"] = [float(i) for i in list(df["x"])]
        df["y"] = [float(i) for i in list(df["y"])]
        df["z"] = [float(i) for i in list(df["z"])]
        df["Rout"] = [float(i) for i in list(df["Rout"])]
        df["Rin"] = [float(i) for i in list(df["Rin"])]
        df["id_ax"] = [float(i) for i in list(df["id_ax"])]
        df["sph_id"] = [float(i) for i in list(df["sph_id"])]
        df["P"] = [float(i) for i in list(df["P"])]
        df["branch_id"] = [float(i) for i in list(df["branch_id"])]
        df["type"] = [str(i) for i in list(df["type"])]

        # rename columns
        df.rename(columns={"Rout": "outer_radius", "Rin": "inner_radius", "x": "X", "y": "Y", "z": "Z", 
                           "id_ax": "cell_id", "sph_id": "component_id", "type": "cell_type", "branch_id": "component_id"}, inplace=True)
        
        df["component"] = list(df["cell_type"])
        df["component"] = df["component"].apply(lambda x: "branch" if x == "Process" else x)
        df["component"] = df["component"].apply(lambda x: "soma" if x == "CellSoma" else x)

        df["cell_type"] = df["cell_type"].apply(lambda x: "glial_cell" if x != "axon" else x)

        del df["P"]

        return df
    
    elif "csv" in file_path:
        print(file_path, " is new format")
        return read_swc_file_new(file_path)

def read_swc_file_new(file_path):
    # Matches CaterpillarGrowth::create_SWC_file's current header exactly:
    # "cell_type cell_id component component_id parent_component_id X Y Z inner_radius outer_radius".
    # parent_component_id was previously missing here, silently shifting every
    # column after component_id one slot to the left (X read as component_id,
    # Y as X, etc.) for any file written after that column was added.
    columns = ["cell_type","cell_id", "component", "component_id", "parent_component_id", "X", "Y", "Z", "inner_radius","outer_radius"]
    df = pd.read_csv(file_path, sep=' ', names=columns)

    df = df.iloc[1:]
    df["X"] = [float(i) for i in list(df["X"])]
    df["Y"] = [float(i) for i in list(df["Y"])]
    df["Z"] = [float(i) for i in list(df["Z"])]
    df["outer_radius"] = [float(i) for i in list(df["outer_radius"])]
    df["inner_radius"] = [float(i) for i in list(df["inner_radius"])]
    df["cell_id"] = [float(i) for i in list(df["cell_id"])]
    df["component_id"] = [float(i) for i in list(df["component_id"])]
    df["parent_component_id"] = [float(i) for i in list(df["parent_component_id"])]
    df["cell_type"] = [str(i) for i in list(df["cell_type"])]
    df["component"] = [str(i) for i in list(df["component"])]
    return df


def plot_radius_power_spectrum(
    file_path,
    which_radius="outer_radius",   # "outer_radius" or "inner_radius"
    detrend="linear",              # "linear" | "mean" | "none"
    step_um=None,                  # None => auto (median center spacing per axon)
    min_points=16,                 # skip very short axons
    mode="mean",                   # "mean" | "each" | "both"
):

    # --- read SWC (your format) ---
    cols = ["cell_type","cell_id","component","component_id","X","Y","Z","inner_radius","outer_radius"]
    df = pd.read_csv(file_path, sep=" ", names=cols).iloc[1:].copy()
    for c in ["X","Y","Z","inner_radius","outer_radius","cell_id","component_id"]:
        df[c] = df[c].astype(float)

    # --- compute PSD per axon ---
    psds = []        # list of (cell_id, freqs, Sxx)
    steps = []       # keep sampling step for info
    for cid, g in df.groupby("cell_id"):
        g = g.sort_values("component_id")
        P = g[["X","Y","Z"]].to_numpy()
        r = g[which_radius].to_numpy()
        if len(P) < 2: 
            continue
        d = np.linalg.norm(P[1:] - P[:-1], axis=1)
        s = np.concatenate([[0.0], np.cumsum(d)])
        L = s[-1]
        if not np.isfinite(L) or L <= 0: 
            continue

        # uniform resampling
        step = (np.median(d[d > 0]) if step_um is None and np.any(d > 0)
                else (L / max(len(P)-1, 1) if step_um is None else float(step_um)))
        s_u = np.arange(0.0, L, step)
        if len(s_u) < min_points: 
            continue
        r_u = np.interp(s_u, s, r)

        # detrend
        if detrend == "mean":
            r_d = r_u - r_u.mean()
        elif detrend == "linear":
            p = np.polyfit(s_u, r_u, 1)
            r_d = r_u - (p[0]*s_u + p[1])
        else:
            r_d = r_u

        # window + PSD (one-sided, Parseval-consistent)
        w = np.hanning(len(r_d))
        xw = r_d * w
        Fs = 1.0 / step                 # samples per µm
        X = np.fft.rfft(xw)
        freqs = np.fft.rfftfreq(len(xw), d=step)  # cycles/µm
        Sxx = (np.abs(X)**2) / (Fs * np.sum(w**2))
        if len(xw) % 2 == 0: Sxx[1:-1] *= 2.0
        else:                Sxx[1:]    *= 2.0

        psds.append((cid, freqs, Sxx))
        steps.append(step)

    if not psds:
        raise ValueError("No valid axons found to compute PSD.")

    # --- plotting ---
    if mode in ("each", "both"):
        plt.figure()
        for cid, f, sxx in psds:
            plt.plot(f, sxx, label=f"cell_id={cid}")
        plt.xlabel("Frequency (cycles / µm)")
        plt.ylabel("PSD")
        plt.title("Power Spectrum of Radius Variation (per axon)")
        plt.grid(True)
        plt.legend()

    if mode in ("mean", "both"):
        # build a common frequency grid and average
        # conservative grid: min max-freq, max min-df
        max_fmax = min(f.max() for _, f, _ in psds if len(f) > 1)
        min_df   = max(np.diff(f).min() for _, f, _ in psds if len(f) > 2)
        fgrid = np.arange(0.0, max_fmax, min_df)
        stack = []
        for _, f, sxx in psds:
            stack.append(np.interp(fgrid, f, sxx))
        mean_psd = np.mean(np.vstack(stack), axis=0)

        plt.figure()
        plt.plot(fgrid, mean_psd)
        plt.xlabel("Frequency (cycles / µm)")
        plt.ylabel("PSD")
        plt.title("Mean Power Spectrum of Radius Variation (across axons)")
        plt.grid(True)

    plt.show()
    
def radius_histogram(df):
    df['R'] = pd.to_numeric(df['R'], errors='coerce')
    df["Diameter"] = df["R"]*2
    df = df.groupby(by = "cell_id").mean()
    sns.histplot(data=df, x="Diameter", color='blue', bins=30, kde=True)
    plt.xlabel('Diameter')
    plt.ylabel('Frequency')
    plt.title('Diameter Histogram')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

def mean_dist_between_maxima(file_path):

    print(file_path)

    df = read_swc_file(file_path)

    axon_ids = df['cell_id'].unique()
    mean_distances_Rin = []
    mean_distances_Rout = []

    # Step 1: Compute mean distances between peaks for Rin and Rout
    for axon_id in axon_ids:
        for radius_key, results_list in zip(['inner_radius', 'outer_radius'], [mean_distances_Rin, mean_distances_Rout]):
            radius = df.loc[df['cell_id'] == axon_id, radius_key].values
            z_values = df.loc[df['cell_id'] == axon_id, 'Z'].values

            peaks, _ = find_peaks(radius,width=0.25, prominence=0.1)

            if len(peaks) > 1:
                z_peaks = z_values[peaks]
                z_distances = np.diff(z_peaks)
                mean_distance = np.mean(z_distances)
                results_list.append(mean_distance)

        # Step 2: Plot example of Rin and Rout profiles with peaks
        random_axon_id = random.choice(axon_ids)
        z_random = df.loc[df['cell_id'] == random_axon_id, 'Z'].values
        radius_rin = df.loc[df['cell_id'] == random_axon_id, 'inner_radius'].values
        radius_rout = df.loc[df['cell_id'] == random_axon_id, 'outer_radius'].values

        peaks_rin, _ = find_peaks(radius_rin, prominence=0.5)
        peaks_rout, _ = find_peaks(radius_rout, prominence=0.5)

        plt.figure(figsize=(10, 6))
        plt.plot(z_random, radius_rin, label='inner_radius', color='blue')
        plt.plot(z_random[peaks_rin], radius_rin[peaks_rin], "bo", label='Peaks Rin')

        plt.plot(z_random, radius_rout, label='outer_radius', color='red')
        plt.plot(z_random[peaks_rout], radius_rout[peaks_rout], "ro", label='Peaks Rout')

        plt.xlabel('z (Position) [μm]')
        plt.ylabel('Radius [μm]')
        plt.title(f'Radius Profiles with Local Maxima (Axon ID: {random_axon_id})')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

        # 95% confidence interval of mean_distances_Rin
        p = np.percentile(mean_distances_Rin, 95)
        print(f"95% percentile for Rin: {p:.2f} μm")

        # Step 3: Histogram with sns.histplot
        plt.figure(figsize=(10, 6))
        sns.histplot(mean_distances_Rin, bins=50, kde=True, color='red', label='Rin')
        sns.histplot(mean_distances_Rout, bins=50, kde=True, color='blue', label='Rout')
        plt.grid(False)
        plt.xlabel('Mean Distance Between Local Maxima [μm]')
        plt.ylabel('Count')
        plt.title('Mean Distance Between Local Maxima')
        plt.tight_layout()
        plt.show()

def diameter_variation(file_path, num_axons=10, max_z=None, num_interp_points=200):
    df = read_swc_file(file_path)

    # Get unique axons and select a subset
    unique_axons = df['cell_id'].unique()
    selected_ids = random.sample(list(unique_axons), min(num_axons, len(unique_axons)))
    selected_ids = [int(i) for i in selected_ids]
    print(f"Selected axon IDs to plot: {selected_ids}")
    
    df_subset = df[df['cell_id'].isin(selected_ids)]

    # Optional z filtering
    if max_z is not None:
        df_subset = df_subset[df_subset['z'] <= max_z]

    # Convert and compute diameter
    df_subset['Z'] = df_subset['Z'].astype(float)
    df_subset['Diameter'] = df_subset['inner_radius'].astype(float) * 2

    plt.figure(figsize=(12, 6))

    # Interpolate and plot each axon
    for axon_id in selected_ids:
        axon_data = df_subset[df_subset['cell_id'] == axon_id].sort_values(by='Z')
        z_vals = axon_data['Z'].values
        diameters = axon_data['Diameter'].values

        if len(z_vals) < 2:
            continue  # skip if not enough points to interpolate

        # Create interpolation function
        interp_func = interp1d(z_vals, diameters, kind='quadratic', bounds_error=False, fill_value="interpolate")

        # Generate interpolated z-values
        z_interp = np.linspace(z_vals.min(), z_vals.max(), num_interp_points)
        d_interp = interp_func(z_interp)

        plt.plot(z_interp, d_interp, label=f'Axon {axon_id}', alpha=0.8)

    plt.xlabel("Z")
    plt.ylabel("Diameter (2 × Rin)")
    plt.title(f"Interpolated Sphere Diameter for {len(selected_ids)} Axons" + (f" (until z={max_z})" if max_z else ""))
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def coefficient_of_variation(file_path):
    df = read_swc_file(file_path)

    axons = df.loc[df["cell_type"] == "axon"]
    # add column "myelinatin" with true if Rin = Rout, otherwise false
    cov_out_mean = axons[["outer_radius", "cell_id"]].groupby(by="cell_id").mean()
    if len(cov_out_mean) == 0:
        return 0.0
    cov_out_std = axons[["outer_radius", "cell_id"]].groupby(by="cell_id").std()
    cov_out = pd.DataFrame()
    cov_out["Coefficient of Variation"] = cov_out_std/cov_out_mean

    cov_out = cov_out.dropna()

    mean_value = cov_out["Coefficient of Variation"].mean()

    return mean_value
    


def create_subplots(file_path):
    df = read_swc_file(file_path)

    df["Outer Diameter"] = df["outer_radius"]*2
    df["Inner Diameter"] = df["inner_radius"]*2
    axons = df.loc[df["cell_type"] == "axon"]
    # add column "myelinatin" with true if Rin = Rout, otherwise false
    axons["myelinated"] = axons["inner_radius"] != axons["outer_radius"]
    #axons= axons.loc[axons["myelinated"] == False]
    # calculate covariance of Rout of each axon
    cov_out_mean = axons[["Outer Diameter", "cell_id"]].groupby(by="cell_id").mean()
    cov_out_std = axons[["Outer Diameter", "cell_id"]].groupby(by="cell_id").std()
    cov_out = pd.DataFrame()
    cov_out["Coefficient of Variation"] = cov_out_std/cov_out_mean
    cov_in_mean = axons[["Inner Diameter", "cell_id"]].groupby(by="cell_id").mean()
    cov_in_std = axons[["Inner Diameter", "cell_id"]].groupby(by="cell_id").std()
    cov_in = pd.DataFrame()
    cov_in["Coefficient of Variation"] = cov_in_std/cov_in_mean
    # delete nans
    cov_in = cov_in.dropna()
    cov_out = cov_out.dropna()
    # add g-ratio column = Rin/Rout
    axons["g_ratio"] = axons["inner_radius"]/axons["outer_radius"]

    print("length cov_out_mean", len(cov_out_mean))
    print("length cov_out_std", len(cov_out_std))
    print("length cov_in_mean", len(cov_in_mean))
    print("length cov_in_std", len(cov_in_std))
    print("length cov_out", len(cov_out))
    print("length cov_in", len(cov_in))




    # define figure with 4 subplots
    fig, axs = plt.subplots(2, 2, figsize=(10, 10))
    # plot "outer_radius" histogram with seaborn
    sns.histplot(data=cov_out_mean, x="Outer Diameter", color='blue',  kde=True, ax=axs[0, 0], label="Outer Diameter", bins=30, binrange = (0, 5))
    sns.histplot(data=cov_in_mean, x="Inner Diameter", color='red',  kde=True, ax=axs[0, 0], label="Inner Diameter", bins=30, binrange = (0, 5))
    axs[0, 0].legend()
    axs[0, 0].set_title('Diameter Histogram')
    axs[0, 0].set_xlabel('Diameter [µm]')
    axs[0, 0].set_ylabel('Count')

    
    # plot Rin wrt g-ratio with seaborn
    sns.scatterplot(data=axons, x="Inner Diameter", y="g_ratio", color = "red", ax=axs[0, 1], label="Inner Diameter")
    sns.scatterplot(data=axons, x="Outer Diameter", y="g_ratio", color = "blue", ax=axs[0, 1], label="Outer Diameter")
    axs[0, 1].set_title('Diameter vs. g-ratio')
    axs[0, 1].set_xlabel('Diameter [µm]')
    axs[0, 1].set_ylabel('g-ratio')
    axs[0, 1].legend()


    # plot covariance of Rout
    sns.histplot(data=cov_out, x= "Coefficient of Variation", color='blue',  kde=True, ax=axs[1, 0], label="Outer Diameter", bins=30, binrange = (0, 1))
    sns.histplot(data=cov_in, x= "Coefficient of Variation", color='red',  kde=True, ax=axs[1, 0], label="Inner Diameter", bins=30, binrange = (0, 1))
    axs[1, 0].legend()
    axs[1, 0].set_title('Coefficient of Variation of Diameter')
    axs[1, 0].set_xlabel('Coefficient of Variation')
    axs[1, 0].set_ylabel('Count')
    

    tortuosities, radii = tortuosity(axons)


    tort = pd.DataFrame()
    tort["Tortuosity"] = tortuosities
    sns.histplot(data = tort, x= "Tortuosity", ax=axs[1, 1], color= "purple", kde = True)
    # set maximum and minimum values for x-axis
    axs[1, 1].set_xlim([1, 1.4])
    axs[1, 1].set_title('Tortuosity Histogram')
    axs[1, 1].set_xlabel('Tortuosity')
    axs[1, 1].set_ylabel('Count')

    plt.tight_layout()
    plt.show()



    
    
    

def all_tortuosity(ondulation_factors, std_deviations):
    all_dfs = []
    for o in ondulation_factors:
        for std in std_deviations:
            new_df = pd.DataFrame(columns=["Tortuosity", "Radius", "Std", "Ondulation_factor"])
            file_path = f"/home/localadmin/Documents/Melina_branch/Sim_Growth/data_std/std_dev_{std}/ond_factor_{o}/growth_icvf_0.10_cap_24_vox_50_factor_2_0.swc"
            df = read_swc_file(file_path)
            tort, rad = tortuosity(df)
            new_df["Tortuosity"] = tort
            new_df["Radius"] = rad
            new_df["Std"] = [std]*len(new_df)
            new_df["Ondulation_factor"] = [o]*len(new_df)
            all_dfs.append(new_df)
    
    df = pd.concat(all_dfs).reset_index()

    print(df)

    sns.boxplot(data =df, y = "Tortuosity", x = "Std", hue = "Ondulation_factor")
    plt.show()


def tortuosity(df):


    str_id = "cell_id"
    str_x = "X"
    str_y = "Y"
    str_z = "Z"
    str_rout = "outer_radius"

    if len(df) == 0:
        print("No axons found")
        return [], []

    nbr_axons = int(df.iloc[len(df)-1 ][str_id])
    tortuosities = []
    radii = []
    for axon in range(nbr_axons):
 
        df_ = df.loc[df[str_id]== axon].reset_index()
        if(len(df_)== 0):
            continue
        
        # Calculate Euclidean distances between consecutive spheres
        distances = []

        for i in range(1, len(df_)):

            distance = math.sqrt((df_.at[i, str_x] - df_.at[i-1, str_x])**2 +
                                (df_.at[i, str_y] - df_.at[i-1, str_y])**2 +
                                (df_.at[i, str_z] - df_.at[i-1, str_z])**2)
    
            distances.append(distance)
            

        # Calculate total length of the axon
        total_length = sum(distances)

        # Calculate distance between the first and last sphere
        first_last_distance = math.sqrt((df_.at[0, str_x] - df_.at[len(df_)-1, str_x])**2 +
                                        (df_.at[0, str_y] - df_.at[len(df_)-1, str_y])**2 +
                                        (df_.at[0, str_z] - df_.at[len(df_)-1, str_z])**2)

        if first_last_distance == 0:
            continue
        # Calculate tortuosity
        tortuosity = total_length / first_last_distance
        tortuosities.append(float(tortuosity))
        radii.append(float(df_.at[0,str_rout]))
   
    return tortuosities, radii


def draw_axons(file_path, glial_only = False):

    N = 5
    colors = mcolors._colors_full_map #dictionary of all colors
    df = read_swc_file(file_path)
    df["color"] = df["cell_id"].apply(lambda x: get_random_element(colors, seed = x)[1])

    # distance from x,y,z to point(0,0,0)
    df["distance_to_point"] = np.linalg.norm(df[["X", "Y", "Z"]], axis=1)

    # sort by distance to point
    df = df.sort_values(by="distance_to_point")

    if glial_only:
        df = df.loc[df["cell_type"] != "axon"]

    # Create a scatter plot for the axons
    scatter = go.Scatter3d(
        x=df["X"],
        y=df["Y"],
        z=df["Z"],
        type="scatter3d",
        mode="markers",
        name="Axons",
        marker=dict(
            sizemode="diameter",
            size=df["outer_radius"]*N,
            color=df["color"],
            line=dict(
                color="rgba(0, 0, 0, 0)",
                width=0
            )
        )
    )


    layout = go.Layout(
        scene=dict(
            xaxis=dict(title='X [µm]'),
            yaxis=dict(title='Y [µm]'),
            zaxis=dict(title='Z [µm]')
        )
    )

    # Create the figure
    fig = go.Figure(data=scatter, layout=layout)
    # Show the figure
    fig.show()

def get_random_element_list(list, seed = 0):
    random.seed(seed)
    return random.choice(list)

def draw_cells(file_path, plot_type="all", axon_indices=None, astrocyte_indices=None):
    """
    Draw structures (axons or astrocytes) in 3D with customizable options.

    Parameters:
    - file_path: str, path to the data file.
    - plot_type: str, "all" to plot everything, "axons" to plot only axons, 
                 or "astrocytes" to plot only astrocytes.
    - axon_indices: list of int, indices of specific axons to plot. If None, all axons are plotted.
    """
    N = 10  # Size scaling factor
    colors = mcolors.CSS4_COLORS

    # Step 1: Read data
    df = read_swc_file(file_path)

    # Step 2: Assign colors to axons
    df["color"] = df["cell_id"].apply(lambda x: get_random_element(colors, seed=x)[1])

    # Step 3: Prepare data for plotting
    if plot_type == "axons" or (plot_type == "all" and axon_indices):
        df = df[df["cell_type"] == "axon"]
        if axon_indices is not None:
            id_axs = df["cell_id"].unique()
            axon_indices = [id_axs[i] for i in axon_indices]
            df = df[df["cell_id"].isin(axon_indices)]
    elif plot_type == "astrocytes":
        df = df[df["cell_type"] != "axon"]
        if astrocyte_indices is not None:
            id_axs = df["cell_id"].unique()
            astrocyte_indices = [id_axs[i] for i in astrocyte_indices]
            df = df[df["cell_id"].isin(astrocyte_indices)]
    elif plot_type == "all":
        pass
    else:
        raise ValueError("Invalid plot_type. Choose 'all', 'axons', or 'astrocytes'.")

    if (len(df) == 0):
        raise ValueError("No data to plot.")
    # Step 4: Clean and sort data
    df["distance_to_point"] = np.linalg.norm(df[["X", "Y", "Z"]], axis=1)
    df = df.sort_values(by="distance_to_point")
    df["outer_radius"] = df["outer_radius"].fillna(0)  # Handle missing sizes
    print(df)
    # Step 5: Create 3D scatter plot
    scatter = go.Scatter3d(
        x=df["X"],
        y=df["Y"],
        z=df["Z"],
        mode="markers",
        name="Structures",
        marker=dict(
            sizemode="diameter",
            size=df["outer_radius"] * N,
            color=["black"]*len(df),
            opacity=0.85,
            line=dict(color=df["color"], width=0)
        )
    )

    # Step 6: Define layout
    layout = go.Layout(
        scene=dict(
            xaxis=dict(title='X [µm]', showbackground=False, showgrid=False, showline=False, showticklabels=False),
            yaxis=dict(title='Y [µm]', showbackground=False, showgrid=False, showline=False, showticklabels=False),
            zaxis=dict(title='Z [µm]', showbackground=False, showgrid=False, showline=False, showticklabels=False),
            aspectmode="auto"
        ),
        paper_bgcolor="white",
        plot_bgcolor="white",
        showlegend=False
    )

    # Step 7: Show the plot
    fig = go.Figure(data=[scatter], layout=layout)
    fig.show()


def draw_axons_black_white(file_path, glial_only=False):
    N = 20  # Size scaling factor
    colors = mcolors.CSS4_COLORS  # Use CSS colors for better compatibility

    # Read and preprocess data
    df = read_swc_file(file_path)
    
    # Color assignment (for spheres we will use white)
    df["color"] = "white"
    
    grey_colors = ["#d3d3d3", "#a9a9a9", "#696969", "#808080", "#778899", "#708090", "#2f4f4f", "#708090", "#778899", "#808080", "#696969", "#a9a9a9", "#d3d3d3"]
    df["color"] = list(map(lambda x, y, z: "white" if x == "axon" else y, list(df["cell_type"]), list(df["color"]), list(df["cell_id"])))

    # Filter for axons in a specific region
    df_ = df.loc[df["cell_type"] == "axon"]
    df_ = df_.groupby("cell_id").first().reset_index()
    id_axs = df_["cell_id"].unique()

    # Filter axons by selected id_axs
    df_axons = df.loc[df["cell_type"] == "axon"]
    df_axons = df_axons.loc[df_axons["cell_id"].isin(id_axs)]

    # Glial cells selection
    df_glial = df.loc[df["cell_type"] != "axon"]
    # df_glial = df_glial.loc[df_glial["cell_id"]==1]
    df = pd.concat([df_axons, df_glial])

    # Calculate distance from (0,0,0) for sorting purposes
    df["distance_to_point"] = np.linalg.norm(df[["X", "Y", "Z"]], axis=1)
    df = df.sort_values(by="distance_to_point")

    if glial_only:
        df = df.loc[df["cell_type"] != "axon"]

    # Ensure marker size is valid
    df["outer_radius"] = df["outer_radius"].fillna(0)  # Replace NaN or missing sizes with a default value

    # Create a 3D scatter plot
    scatter = go.Scatter3d(
        x=df["X"],
        y=df["Y"],
        z=df["Z"],
        mode="markers",
        name="Axons",
        marker=dict(
            sizemode="diameter",
            size=df["outer_radius"] * N,  # Scale the marker size
            color=df["color"],  # Set the color to white
            opacity=0.5,  # Set transparency to 50%
            line=dict(
                color="rgba(0, 0, 0, 0)",
                width=0
            )
        )
    )

    # Layout configuration with black background
    layout = go.Layout(
        scene=dict(
            xaxis=dict(title='X [µm]', showbackground=False, showgrid=False, showline=False, showticklabels=False),
            yaxis=dict(title='Y [µm]', showbackground=False, showgrid=False, showline=False, showticklabels=False),
            zaxis=dict(title='Z [µm]', showbackground=False, showgrid=False, showline=False, showticklabels=False),
            aspectmode="auto"  # Ensure the aspect ratio is set correctly
        ),
        paper_bgcolor='black',  # Black background
        plot_bgcolor='black',   # Black background for the plot
        showlegend=False
    )

    # Create the figure
    fig = go.Figure(data=[scatter], layout=layout)

    # Show the figure
    fig.show()
def get_spheres_array(df):
    axons = []
    current_axon_id = None
    current_axon = []
    
    df['x'] = pd.to_numeric(df['x'], errors='coerce')
    df['y'] = pd.to_numeric(df['y'], errors='coerce')
    df['z'] = pd.to_numeric(df['z'], errors='coerce')
    df['Rout'] = pd.to_numeric(df['Rout'], errors='coerce')

    for _, row in df.iterrows(): # loops over each row of the df
        id_ax = row["cell_id"]
        x = row["X"]
        y = row["Y"]
        z = row["Z"]
        r = row["outer_radius"]  
        type_ = row["cell_type"]

        if id_ax != current_axon_id: # passing to next axon 
            if current_axon_id is not None:
                axons.append(current_axon) # full list
                current_axon = [] # emptying list
            current_axon_id = id_ax # update axon number

        current_axon.append([x, y, z, r, type_])

    if current_axon: # if the list is not empty
        axons.append(current_axon) # last axon

    return axons



def draw_circles(center_radii):
    """
    Draw circles in a 2D plot.

    Parameters:
        center_radii (list of tuples): List of tuples, where each tuple contains (x, y, radius).
        ax (matplotlib.axes._axes.Axes, optional): Axes object to draw the circles on. If not provided, a new plot will be created.
        **kwargs: Additional keyword arguments to customize the appearance of circles.
    """

    fig, ax = plt.subplots()

    for x, y, radius, color in center_radii:
        circle = Circle((x, y), radius, color=color, alpha=0.5)
        ax.add_patch(circle)

    #ax.set_aspect('equal', adjustable='datalim')  # Equal aspect ratio

    return ax

def find_closest_to_value(list_of_lists, value):
    closest_diff = float('inf')
    closest_list = None
    
    for e, sublist in enumerate(list_of_lists):
        if len(sublist) >= 3:

            diff = abs(sublist[2] - value)
            if diff < closest_diff:
                closest_diff = diff
                index = e
                closest_value = sublist[2]
    
    return index, closest_value

def draw_spheres(file_path, limit, z_slice):
    df = read_swc_file(file_path)
    colors = mcolors._colors_full_map #dictionary of all colors
    circles = []

    axons = df.loc[df["cell_type"] == "axon"]
    for axon in axons["cell_id"].unique():
        axon_i = axons.loc[axons["cell_id"] == axon]


        axon_in_slice = axon_i.loc[(axon_i["Z"]-axon_i["inner_radius"] < z_slice) & (z_slice < axon_i["Z"]+axon_i["inner_radius"])]
        random_key, random_value = get_random_element(colors, seed = axon)
        c = random_value
        for i, row in axon_in_slice.iterrows():
            x = row["X"]
            y = row["Y"]
            z = row["Z"]
            r = row["inner_radius"]
            
            Rnew = math.sqrt(r*r - (z - z_slice)**2)
            circles.append((x, y, Rnew, c))

        

    glial_df = df.loc[df["cell_type"] != "axon"]
    for glial in glial_df["cell_id"].unique():
        glial_i = glial_df.loc[glial_df["cell_id"] == glial]
        glial_in_slice = glial_i.loc[(glial_i["Z"]-glial_i["outer_radius"] < z_slice) & (z_slice < glial_i["Z"]+glial_i["outer_radius"])]
        random_key, random_value = get_random_element(colors, seed = glial)
        c = random_value
        for i, row in glial_in_slice.iterrows():
            x = row["X"]
            y = row["Y"]
            z = row["Z"]
            r = row["outer_radius"]
            
            Rnew = math.sqrt(r*r - (z - z_slice)**2)
            circles.append((x, y, Rnew, c))

    draw_circles(circles)

    # Set plot limits
    plt.xlim(0, limit)
    plt.ylim(0, limit)

    plt.title('2D Circles')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.grid()
    plt.show()

def read_data(filename):
    with open(filename, 'r') as file:
        data = {}
        for line in file:
            key, value = line.strip().split()
            data[key] = float(value)
        return data

def combine_files(file_list):
    combined_data = {}
    for file in file_list:
        data = read_data(file)
        for key, value in data.items():
            if key in combined_data:
                combined_data[key].append(value) # adds element to key list
            else:
                combined_data[key] = [value] # creates new key and adds element
    return combined_data

def vox_time_plot(file_list):
    data = combine_files(file_list)
    # sns.lineplot(x=data['Voxel'], y=data['Duration'])
    # plt.xlabel('Voxel size (µm)')
    # plt.ylabel('Time (s)')
    # plt.title('Time vs. Voxel Size')
    # plt.legend(title='icvf = ' + str(data['icvf'][0]) + '\n' + 'capacity = ' + str(data['Capacity'][0]), loc='upper left')
    # plt.show()
    df = pd.DataFrame(data)


    unique_cap_values = df['Capacity'].unique()

    
    for capacity in unique_cap_values:
        subset_data = df[df['Capacity'] == capacity]
        
        sns.lineplot(x='Voxel', y='Duration', data=subset_data, label=f'capacity = {capacity}')

    plt.xlabel('Voxel size (µm)')
    plt.ylabel('Time (s)')
    plt.title('Time vs. Voxel Size for Different Capacity Values')
    plt.legend(loc='upper left')
    plt.show()

def cap_time_plot(file_list):
    data = combine_files(file_list)
    keys = ['Duration', 'Capacity', 'Voxel']
    data = {key: data[key] for key in keys}

    df = pd.DataFrame(data)
    df['Duration'] =  df['Duration']/60

    df_pivot = df.pivot_table(index='Voxel', columns='Capacity', values='Duration', aggfunc='mean')
    

    # Set the heatmap parameters
    sns.heatmap(df_pivot,
                annot=True,
                fmt=".2f",  # Format for the annotations (optional, adjust as needed)
                cmap='viridis')

    plt.title('Simulation Heat Map')
    plt.xlabel('Capacity')
    plt.ylabel('Voxel size (um)')
    plt.show()

def get_text_from_folder(folder_path,  straight = False):

    txt_files = glob.glob(os.path.join(folder_path, f"*.txt"))
    txt_files.sort()
    list_txt = []
    for txt_file in txt_files:
        if straight and "straight" in txt_file:
            list_txt.append(txt_file)
        elif not straight and "straight" not in txt_file:
            list_txt.append(txt_file)
    return list_txt
    
def get_swc_from_folder(folder_path, icvf):

    txt_files = glob.glob(os.path.join(folder_path, f"growth*"))
    txt_files.sort()
    list_txt = []
    if icvf != None:
        for txt_file in txt_files:
            if (str(icvf) in txt_file):
                list_txt.append(txt_file)
        return list_txt
    else:
        return txt_files

def read_swc_file_np(file_path, axon_nbr):
# Define column names
    columns = ["cell_id", "sph_id", "branch_id", "cell_type", "X", "Y", "Z", "inner_radius", "outer_radius", "P"]
    
    # Read the file using numpy.genfromtxt to handle non-numeric values
    data_ = np.genfromtxt(file_path, skip_header=1, dtype=None, encoding=None)

    # Convert selected columns to float
    data = np.ones((len(data_), len(columns)))*np.nan
    for i in range(len(data_)) :
        if(data_[i][3] != "axon" and int(data_[i][0]) == axon_nbr):
    
            for j in range(len(columns)):
                if (columns[j] != "cell_type"):
                    data[i][j] = float(data_[i][j])

    # Find rows containing NaN values
    nan_mask = np.isnan(data).all(axis=1)

    # Filter out rows with NaN values
    data = data[~nan_mask]


    # Create a dictionary mapping column names to their corresponding data
    data_dict = {column: data[:, i] for i, column in enumerate(columns)}
    
    return data_dict

def get_random_element(dictionary, seed):
    random.seed(seed)
    key, value = random.choice(list(dictionary.items()))
    return key, value

def draw_spheres_pyvista(file_path, cell_types=None, chosen_id=None):
    # Default to these three if no list is provided
    if cell_types is None:
        cell_types = ["axon", "glial", "glialRamification", "blood_vessel" ]
        
    df = read_swc_file(file_path)
    colors = mcolors._colors_full_map # dictionary of all colors
    plotter = pv.Plotter()
    scale = 1
    
    # Loop through only the cell types requested
    for c_type in cell_types:
        df_filtered = df[df["cell_type"] == c_type]
        
        if df_filtered.empty:
            print(f"No data found for cell type: {c_type}")
            continue
            
        # Process each unique cell ID for the current cell type
        for cell_id in tqdm(df_filtered["cell_id"].unique(), desc=f"Processing {c_type}"):
            # Skip if a specific ID is chosen and this isn't it
            if chosen_id is not None and cell_id != chosen_id:
                continue
                
            cell_data = df_filtered[df_filtered["cell_id"] == cell_id]
            N = np.array(cell_data[["X", "Y", "Z"]].astype(float))
            Rout = np.array(cell_data["outer_radius"])
            
            # Determine color seed
            seed_val = chosen_id if chosen_id is not None else cell_id
            _, c = get_random_element(colors, seed=seed_val)
            
            # Axons have specific inner/outer radius opacity logic
            if c_type == "axon":
                Rin = np.array(cell_data["inner_radius"])
                for i, p in enumerate(N):
                    if Rin[i] != Rout[i]:
                        plotter.add_points(p, render_points_as_spheres=True, point_size=Rout[i]*scale, color=c, opacity=0.4)
                        plotter.add_points(p, render_points_as_spheres=True, point_size=Rin[i]*scale, color=c, opacity=1.0)
                    else:
                        plotter.add_points(p, render_points_as_spheres=True, point_size=Rout[i]*scale, color=c, opacity=1.0)
                        
            # All other cell types just use the outer radius
            else:
                for i, p in enumerate(N):
                    plotter.add_points(p, render_points_as_spheres=True, point_size=Rout[i]*scale, color=c)

    plotter.show()

def draw_one_axon_pyvista(file_path, chosen_id = None):
    df = read_swc_file(file_path)
    colors = mcolors._colors_full_map #dictionary of all colors
    plotter = pv.Plotter()
    df_axons = df[df["cell_type"] == "axon"]

    scale = 150

    axon_i = df_axons.loc[df_axons["cell_id"] == chosen_id]

    N = np.array(axon_i[["X", "Y", "Z"]].astype(float))
    Rout = np.array(axon_i["outer_radius"])
    Rin = np.array(axon_i["inner_radius"])
    random_key, random_value = get_random_element(colors, seed = chosen_id)
    c = random_value

    for i, p in enumerate(N):
        #plotter.add_mesh(pv.PolyData(p), point_size=R[i]*scale, color=c, render_points_as_spheres=True)
        if Rin[i] != Rout[i]:
            plotter.add_points(p, render_points_as_spheres=True, point_size=Rout[i]*scale, color="blue", opacity = 0.1)
            plotter.add_points(p, render_points_as_spheres=True, point_size=Rin[i]*scale, color="red", opacity = 1)
        else:
            plotter.add_points(p, render_points_as_spheres=True, point_size=Rout[i]*scale, color="purple", opacity = 0.8)

    plotter.show()

def draw_one_glial_pyvista(file_path, chosen_id = None):
    df = read_swc_file(file_path)
    colors = mcolors._colors_full_map #dictionary of all colors
    plotter = pv.Plotter()
    df_axons = df[df["cell_type"] != "axon"]

    scale = 20

    axon_i = df_axons.loc[df_axons["cell_id"] == chosen_id]

    N = np.array(axon_i[["X", "Y", "Z"]].astype(float))
    Rout = np.array(axon_i["outer_radius"])
    Rin = np.array(axon_i["inner_radius"])
    random_key, random_value = get_random_element(colors, seed = chosen_id)
    c = random_value

    for i, p in enumerate(N):

        plotter.add_points(p, render_points_as_spheres=True, point_size=Rout[i]*scale, color="black", opacity = 0.8)

    plotter.show()


def sholl_intersection(file_path):
    df = read_swc_file(file_path)
    sphere_around_soma_radii = [5, 7, 10, 15, 20 , 25, 30, 40, 50, 60, 80]
    intersections_list_all = []
    glial_df = df.loc[df["cell_type"] != "axon"]
    for glial in glial_df["cell_id"].unique():
        intersections_list = []
        glial_i = glial_df.loc[glial_df["cell_id"] == glial]
        soma_glial = glial_i.loc[glial_i["cell_type"] == "CellSoma"]
        print(soma_glial)
        soma_position = np.array(soma_glial[["X", "Y", "Z"]].astype(float))[0]
        print(soma_position)
        processes_glial = glial_i.loc[glial_i["cell_type"] == "Process"]
        processes_glial["distance_to_soma"] = np.linalg.norm(np.array(processes_glial[["X", "Y", "Z"]].astype(float)) - soma_position, axis=1)
        # sort by distance to soma
        processes_glial = processes_glial.sort_values(by="distance_to_soma")
        # keep only first and last element of each branch_id
        processes_glial_first = processes_glial.drop_duplicates(subset=["branch_id"], keep="first")
        processes_glial_last = processes_glial.drop_duplicates(subset=["branch_id"], keep="last")

        distance_first = np.array(processes_glial_first[["distance_to_soma"]].astype(float))
        distance_last = np.array(processes_glial_last[["distance_to_soma"]].astype(float))
        # find how many times glial cross each sphere
        for r in sphere_around_soma_radii:
            intersections = 0
            for i, p in enumerate(distance_first):
                if distance_first[i] < r  and distance_last[i] > r :
                    intersections += 1
            intersections_list.append(intersections)
        intersections_list_all.append(intersections_list)
    return intersections_list_all, sphere_around_soma_radii

def sholl_intersections(file_path_WM, file_path_GM = None):

    intersections_list_WM, sphere_around_soma_radii = sholl_intersection(file_path_WM)
    if file_path_GM is not None:
        intersections_list_GM, sphere_around_soma_radii = sholl_intersection(file_path_GM)

        # plot sholl intersections wrt distance from soma
        plt.figure()
        for intersections, tissue in zip([intersections_list_WM, intersections_list_GM], ["With axons", "Without axons"]):
            if (tissue == "With axons"):
                plt.plot(sphere_around_soma_radii, np.mean(intersections, axis=0), color='blue', label=tissue)
            else:
                plt.plot(sphere_around_soma_radii, np.mean(intersections, axis=0), color='red', label=tissue)
            for intersections in intersections:
                if (tissue == "With axons"):
                    plt.plot(sphere_around_soma_radii, intersections, alpha=0.1, color='blue')
                else:
                    plt.plot(sphere_around_soma_radii, intersections, alpha=0.1, color='red')
        plt.xlabel("Distance from soma (µm)")
        plt.ylabel("Intersections")
        plt.title("Sholl intersections")
        plt.legend()
        plt.show()
    else:
        plt.figure()
        plt.plot(sphere_around_soma_radii, np.mean(intersections_list_WM, axis=0), color='blue')
        for intersections in intersections_list_WM:
            plt.plot(sphere_around_soma_radii, intersections, alpha=0.1, color='blue')
        plt.xlabel("Distance from soma (µm)")
        plt.ylabel("Intersections")
        plt.title("Sholl intersections")
        plt.show()

def varying_sholl(files, stds, lengths):
    all_intersections = []
    #colors = ["blue", "red", "purple", "orange", "black", "yellow", "pink", "brown", "grey", "cyan"]
    colors = ["blue", "black", "purple", "grey",  "cyan"]
    for i, file in enumerate(files):
        intersections_list, sphere_around_soma_radii = sholl_intersection(file)
        all_intersections.append(intersections_list)

    plt.figure()
    for i, intersections in enumerate(all_intersections):
        c = colors[i]
        plt.plot(sphere_around_soma_radii, np.mean(intersections, axis=0),  label = f"Mean process length = {lengths[i]}, std = {stds[i]}", color= c)
        for intersections in intersections:
            plt.plot(sphere_around_soma_radii, intersections, alpha=0.1, color= c)
    plt.xlabel("Distance from soma (µm)")
    plt.ylabel("Intersections")
    plt.title("Sholl intersections")
    plt.legend()
    plt.show()

def vox_size_analysis(files):
    df = pd.DataFrame()
    for file in files:
        data = read_data(file)
        df = pd.concat([df, pd.DataFrame(data, index=[0])])
    df["Log(Duration)"] = np.log(df["Duration"])
    print(df)
    #plot
    sns.boxplot(data = df, y = "Log(Duration)", x = "Voxel")
    sns.swarmplot(data = df, y = "Log(Duration)", x = "Voxel", color="black")
    plt.ylabel("Log (Duration (s))")
    plt.xlabel("Voxel cube length (µm)")
    plt.show()


def tortuosity_ring_plot(path, subsample=10000, tau=10.0, axial=True):
    df = read_swc_file(path)
    diff_vectors = []
    for axon_id in df["cell_id"].unique():
        coords = df.loc[df["cell_id"]==axon_id, ["X","Y","Z"]].to_numpy()
        for i in range(1, len(coords)):
            v = coords[i] - coords[i-1]
            n = np.linalg.norm(v)
            if n > 0:
                diff_vectors.append(v / n)

    V = np.array(diff_vectors)
    if V.size == 0:
        print("No vectors found."); return

    # Axial Watson: identify antipodal directions if desired
    if axial:
        # Map all vectors to the upper hemisphere to respect axial symmetry
        V = np.where(V[:,2:3] >= 0, V, -V)

    if subsample and V.shape[0] > subsample:
        rng = np.random.default_rng(0)
        V = V[rng.choice(V.shape[0], subsample, replace=False)]

    # Rotation-invariant spherical KDE (Watson/vMF kernel)
    density = watson_kde_density(V, tau=tau, axial=axial)

    # Plot colored points on the unit sphere
    xs, ys, zs = V[:,0], V[:,1], V[:,2]
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111, projection='3d')
    sc = ax.scatter(xs, ys, zs, c=density, cmap="viridis", s=12, alpha=0.85)

    # wireframe sphere
    u = np.linspace(0, 2*np.pi, 60)
    v = np.linspace(0, np.pi, 30)
    ax.plot_wireframe(np.outer(np.cos(u), np.sin(v)),
                      np.outer(np.sin(u), np.sin(v)),
                      np.outer(np.ones_like(u), np.cos(v)),
                      color="lightgray", linewidth=0.4, alpha=0.35)

    fig.colorbar(sc, ax=ax, shrink=0.6, label="Watson-KDE (relative)")
    ax.set_box_aspect([1,1,1])
    ax.set_title(f"Spherical density (Watson KDE, tau={tau}, axial={axial})")
    plt.show()



def watson_kde_density(X, tau=50.0, axial=True):
    """
    X: (N,3) unit vectors.
    tau: concentration parameter (larger -> tighter kernel). Try 10..200.
    axial: if True, uses axial Watson kernel exp(tau * (x·y)^2),
           else uses vMF-like exp(tau * (x·y)).
    Returns: (N,) unnormalized densities at each sample (leave as relative colors).
    """
    # Optionally fold antipodes for axial symmetry
    if axial:
        # work with |dot| by squaring later; no fold needed
        pass
    else:
        # directional: ensure X is unit
        X = X / np.linalg.norm(X, axis=1, keepdims=True)

    # Compute all pairwise dot products efficiently in blocks if N large
    # Here simple version; for large N use chunking.
    D = X @ X.T                       # (N,N) with dot products
    if axial:
        K = np.exp(tau * (D**2))      # Watson kernel
    else:
        K = np.exp(tau * D)           # vMF-like kernel

    # Leave out self-term bias if you like; it barely matters for coloring
    np.fill_diagonal(K, 0.0)
    dens = K.mean(axis=1)
    return dens

def p2_plot(path, subsample=10000, tau=80.0, axial=True):
    df = read_swc_file(path)
    diff_vectors = []
    for axon_id in df["cell_id"].unique():
        pts = df.loc[df["cell_id"] == axon_id, ["X", "Y", "Z"]].to_numpy()
        if pts.shape[0] < 2: 
            continue
        v = pts[-1] - pts[0]
        n = np.linalg.norm(v)
        if n > 0:
            diff_vectors.append(v / n)

    V = np.array(diff_vectors)
    if V.size == 0:
        print("No vectors found."); return

    # Axial Watson: identify antipodal directions if desired
    if axial:
        # Map all vectors to the upper hemisphere to respect axial symmetry
        V = np.where(V[:,2:3] >= 0, V, -V)

    if subsample and V.shape[0] > subsample:
        rng = np.random.default_rng(0)
        V = V[rng.choice(V.shape[0], subsample, replace=False)]

    # Rotation-invariant spherical KDE (Watson/vMF kernel)
    density = watson_kde_density(V, tau=tau, axial=axial)

    # Plot colored points on the unit sphere
    xs, ys, zs = V[:,0], V[:,1], V[:,2]
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111, projection='3d')
    sc = ax.scatter(xs, ys, zs, c=density, cmap="viridis", s=12, alpha=0.85)

    # wireframe sphere
    u = np.linspace(0, 2*np.pi, 60)
    v = np.linspace(0, np.pi, 30)
    ax.plot_wireframe(np.outer(np.cos(u), np.sin(v)),
                      np.outer(np.sin(u), np.sin(v)),
                      np.outer(np.ones_like(u), np.cos(v)),
                      color="lightgray", linewidth=0.4, alpha=0.35)

    fig.colorbar(sc, ax=ax, shrink=0.6, label="Watson-KDE (relative)")
    ax.set_box_aspect([1,1,1])
    ax.set_title(f"Spherical density (Watson KDE, tau={tau}, axial={axial})")
    plt.show()
        
def tortuosity_plot(folder_path):
    
    folders =[ "std_0.1", "std_0.2", "std_0.3", "std_0.4", "std_0.5"] 
    tortuosities = []
    stds =[] 
    
    for folder in folders:
        print(folder)
        all_dfs = []
        path = os.path.join(folder_path, folder)
        files =  glob.glob(os.path.join(path, f"growth*"))
        for file in files:
            print(file)
            df = read_swc_file(file)
            tort, radii = tortuosity(df)
            tortuosities.extend(tort)
            stds.extend([float(folder.split("_")[1])]*len(tort))
            df_copy = copy.deepcopy(df)
            all_dfs.append(df_copy)


    df_final = pd.DataFrame()
    df_final["Tortuosity"] = tortuosities
    df_final["Std"] = stds
    # Your lmplot with customized colors for dots and regression line
    sns.lmplot(data=df_final, y="Tortuosity", x="Std", order=2,
            scatter_kws={'color': 'black'},   # Set dots to black
            line_kws={'color': 'red'})        # Set line to red

    # Display the plot
    plt.show()

def plot_radii(file_path, cell_type):
    swc_df = read_swc_file(file_path)

    swc_df = swc_df[swc_df["cell_type"] == cell_type].copy()

    if cell_type =="glial_cell":
        swc_df = swc_df[swc_df["component"] == "soma"].copy()

    # mean inner and outer radius per cell
    mean_radii = swc_df.groupby("cell_id")[["outer_radius", "inner_radius"]].mean().reset_index()
    # plot histogram of mean inner and outer radius
    plt.figure(figsize=(12, 6))
    plt.subplot(1, 2, 1)
    plt.hist(mean_radii["outer_radius"], bins=20, color="blue", alpha=0.7, label="Outer Radius")
    plt.hist(mean_radii["inner_radius"], bins=20, color="#d62728", alpha=0.7, label="Inner Radius")
    plt.xlabel("Mean Radius (µm)")
    plt.ylabel("Frequency")
    plt.title(f"Histogram of Mean Radii for {cell_type.capitalize()}")
    plt.legend()
    plt.subplot(1, 2, 2)
    plt.scatter(mean_radii["outer_radius"], mean_radii["inner_radius"], color="purple", alpha=0.7)
    plt.xlabel("Mean Outer Radius (µm)")
    plt.ylabel("Mean Inner Radius (µm)")
    plt.title(f"Mean Inner vs Outer Radius for {cell_type.capitalize()}")
    plt.grid(True)
    plt.show()


def plot_radii_with_length(file_path, axon_id, z_step=0.005):
    swc_df = read_swc_file(file_path)

    axon_df = swc_df[swc_df["cell_id"] == axon_id].copy()
    
    # --- Interpolation of effective radius along z-axis ---
    # Extract spheres
    spheres = axon_df[["X", "Y", "Z", "outer_radius", "inner_radius"]].values

    # Define interpolation range
    z_min = spheres[:, 2].min()
    z_max = spheres[:, 2].max()
    z_values = np.arange(z_min, z_max + z_step, z_step)
    interpolated_outer_radii = []
    
    for z in z_values:
        r_effective = 0
        for x, y, z0, r, _ in spheres:
            dz = z - z0
            if abs(dz) <= r:
                r_z = np.sqrt(r**2 - dz**2)
                r_effective = max(r_effective, r_z)
        interpolated_outer_radii.append(r_effective)

    interpolated_inner_radii = []
    
    for z in z_values:
        r_effective = 0
        for x, y, z0, _, r in spheres:
            dz = z - z0
            if abs(dz) <= r:
                r_z = np.sqrt(r**2 - dz**2)
                r_effective = max(r_effective, r_z)
        interpolated_inner_radii.append(r_effective)

    # Plot interpolated radius profile
    plt.plot(z_values, interpolated_outer_radii, color="blue", linewidth=5, alpha = 0.7)
    plt.plot(z_values, interpolated_inner_radii,  color = "#d62728", linewidth=5, alpha = 0.7)
    plt.xlabel("Z")
    plt.ylabel("Effective radius at z")
    plt.title(f"Interpolated Cross-Sectional Radius for Axon {axon_id}")
    plt.grid(True)
    plt.show()


if __name__ == "__main__":

    file_path = "/home/localadmin/Documents/CATERPillar/tests/run_75pct_eps03_uf5_4thread_100um_swf100_rerun.csv"

    create_subplots(file_path)

    


    