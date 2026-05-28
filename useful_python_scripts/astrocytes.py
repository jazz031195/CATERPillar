import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from simulationgraphs import read_swc_file
from scipy.interpolate import interp1d
from scipy.stats import ttest_ind
import matplotlib.patches as patches
from statsmodels.stats.multitest import multipletests

def sholl_intersection(file_path):
    df = read_swc_file(file_path)
    sphere_around_soma_radii = [5, 7, 10, 15, 20 , 25, 30, 40, 50, 60, 80]
    intersections_list_all = []
    glial_df = df.loc[df["cell_type"] != "axon"]
    for glial in glial_df["cell_id"].unique():
        intersections_list = []
        glial_i = glial_df.loc[glial_df["cell_id"] == glial]
        soma_glial = glial_i.loc[glial_i["component"] == "soma"]
        print(soma_glial)
        soma_position = np.array(soma_glial[["X", "Y", "Z"]].astype(float))[0]
        print(soma_position)
        processes_glial = glial_i.loc[glial_i["component"] == "branch"]
        processes_glial["distance_to_soma"] = np.linalg.norm(np.array(processes_glial[["X", "Y", "Z"]].astype(float)) - soma_position, axis=1)
        # sort by distance to soma
        processes_glial = processes_glial.sort_values(by="distance_to_soma")
        # keep only first and last element of each branch_id
        processes_glial_first = processes_glial.drop_duplicates(subset=["component_id"], keep="first")
        processes_glial_last = processes_glial.drop_duplicates(subset=["component_id"], keep="last")

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
                plt.plot(sphere_around_soma_radii, np.mean(intersections, axis=0), color='blue')
            else:
                plt.plot(sphere_around_soma_radii, np.mean(intersections, axis=0), color='red')
            for intersections in intersections:
                if (tissue == "With axons"):
                    plt.plot(sphere_around_soma_radii, intersections, alpha=0.05, color='blue')
                else:
                    plt.plot(sphere_around_soma_radii, intersections, alpha=0.05, color='red')
        plt.xlabel("Distance from soma (µm)", fontsize=14)
        plt.ylabel("Intersections", fontsize=14)
        plt.title("Sholl intersections", fontsize=14)
        plt.tick_params(axis='both', labelsize=14)
        plt.legend()
        plt.show()
        plt.legend().set_visible(False)
    else:
        plt.figure()
        plt.plot(sphere_around_soma_radii, np.mean(intersections_list_WM, axis=0), color='blue')
        for intersections in intersections_list_WM:
            plt.plot(sphere_around_soma_radii, intersections, alpha=0.05, color='blue')
        plt.xlabel("Distance from soma (µm)")
        plt.ylabel("Intersections")
        plt.title("Sholl intersections")
        plt.show()


def total_process_length(file_path_WM, file_path_GM):
    df = read_swc_file(file_path_WM)
    glial_df = df.loc[df["cell_type"] != "axon"]
    all_process_lengths_wm = []
    all_mean_process_length_wm = []
    nbr_processes_wm = []
    all_std_process_length_wm = []
    all_radii_wm = []
    for glial in glial_df["cell_id"].unique():
        glial_i = glial_df.loc[glial_df["cell_id"] == glial]
        processes_glial = glial_i.loc[glial_i["component"] == "branch"]
        # keep only first and last element of each branch_id
        processes_glial_first = processes_glial.drop_duplicates(subset=["component_id"], keep="first")
        processes_glial_last = processes_glial.drop_duplicates(subset=["component_id"], keep="last")
        # calculate length of each process
        process_lengths = []
        mean_radii = []
        nbr_process = len(processes_glial_first)
        for i in range(len(processes_glial_first)):
            coord_first = np.array(processes_glial_first.iloc[i][["X", "Y", "Z"]].astype(float))
            coord_last = np.array(processes_glial_last.iloc[i][["X", "Y", "Z"]].astype(float))
            length = np.linalg.norm(coord_first - coord_last)
            process_lengths.append(length)
            mean_radii.append(float(processes_glial_first.iloc[i]["outer_radius"]))

        if (nbr_process == 0):
            continue
        all_std_process_length_wm.append(np.std(process_lengths))
        all_mean_process_length_wm.append(np.mean(process_lengths))
        all_process_lengths_wm.append(np.sum(process_lengths))
        nbr_processes_wm.append(nbr_process)
        all_radii_wm.append(np.mean(mean_radii))

    df = read_swc_file(file_path_GM)
    glial_df = df.loc[df["cell_type"] != "axon"]
    all_process_lengths_gm = []
    nbr_processes_gm = []
    all_mean_process_length_gm = []
    all_std_process_length_gm = []
    all_radii_gm = []

    for glial in glial_df["cell_id"].unique():
        glial_i = glial_df.loc[glial_df["cell_id"] == glial]
        processes_glial = glial_i.loc[glial_i["component"] == "branch"]
        # keep only first and last element of each branch_id
        processes_glial_first = processes_glial.drop_duplicates(subset=["component_id"], keep="first")
        processes_glial_last = processes_glial.drop_duplicates(subset=["component_id"], keep="last")
        # calculate length of each process
        process_lengths = []
        mean_radii = []
        nbr_process = len(processes_glial_first)
        for i in range(len(processes_glial_first)):
            coord_first = np.array(processes_glial_first.iloc[i][["X", "Y", "Z"]].astype(float))
            coord_last = np.array(processes_glial_last.iloc[i][["X", "Y", "Z"]].astype(float))
            length = np.linalg.norm(coord_first - coord_last)
            process_lengths.append(length)
            mean_radii.append( float(processes_glial_first.iloc[i]["outer_radius"]))
        if (nbr_process == 0):
            continue
        all_mean_process_length_gm.append(np.mean(process_lengths))
        all_process_lengths_gm.append(np.sum(process_lengths))
        all_std_process_length_gm.append(np.std(process_lengths))
        nbr_processes_gm.append(nbr_process)
        all_radii_gm.append(np.mean(mean_radii))

    data = pd.DataFrame({
        "Glial cells": ["GM"] * len(all_process_lengths_gm) + ["WM"] * len(all_process_lengths_wm),
        "Length (µm)": all_process_lengths_gm + all_process_lengths_wm,
        "Nbr processes": nbr_processes_gm + nbr_processes_wm,
        "Mean process length (µm)": all_mean_process_length_gm + all_mean_process_length_wm,
        "Std process length (µm)": all_std_process_length_gm + all_std_process_length_wm,
        "Mean radius (µm)": all_radii_gm + all_radii_wm
    })
    plot_(data, primary = False)


def total_central_process_length(file_path_WM, file_path_GM):
    df = read_swc_file(file_path_WM)
    glial_df = df.loc[df["cell_type"] != "axon"]
    all_total_process_lengths_wm = []
    nbr_processes_wm = []
    all_mean_process_length_wm = []
    all_std_process_length_wm = []
    all_radii_wm = []
    
    for glial in glial_df["cell_id"].unique():
        glial_i = glial_df.loc[glial_df["cell_id"] == glial]
        processes_glial = glial_i.loc[glial_i["component"] == "branch"]
        # keep only first and last element of each branch_id
        processes_glial_first = processes_glial.drop_duplicates(subset=["component_id"], keep="first")
        processes_glial_last = processes_glial.drop_duplicates(subset=["component_id"], keep="last")
        coord_soma = np.array(glial_i.loc[glial_i["component"] == "soma"][["X", "Y", "Z"]].astype(float))
        radius_soma = float(glial_i.loc[glial_i["component"] == "soma"]["outer_radius"])
        # calculate length of each process
        process_lengths = []
        nbr_process = 0
        mean_radii = []
        for i in range(len(processes_glial_first)):
            coord_first = np.array(processes_glial_first.iloc[i][["X", "Y", "Z"]].astype(float))
            distance = np.linalg.norm(coord_first - coord_soma)
            if (distance > radius_soma/2):
                continue
            coord_last = np.array(processes_glial_last.iloc[i][["X", "Y", "Z"]].astype(float))
            length = np.linalg.norm(coord_first - coord_last)
            process_lengths.append(length)
            nbr_process += 1
            mean_radii.append(float(processes_glial_first.iloc[i]["outer_radius"]))

        if (nbr_process == 0):
            continue
        all_mean_process_length_wm.append(np.mean(process_lengths))
        all_std_process_length_wm.append(np.std(process_lengths))
        all_total_process_lengths_wm.append(np.sum(process_lengths))
        nbr_processes_wm.append(nbr_process)
        all_radii_wm.append(np.mean(mean_radii))

    df = read_swc_file(file_path_GM)
    glial_df = df.loc[df["cell_type"] != "axon"]
    all_total_process_lengths_gm = []
    nbr_processes_gm = []
    all_mean_process_length_gm = []
    all_std_process_length_gm = []
    all_radii_gm = []
    for glial in glial_df["cell_id"].unique():
        glial_i = glial_df.loc[glial_df["cell_id"] == glial]
        processes_glial = glial_i.loc[glial_i["component"] == "branch"]
        # keep only first and last element of each branch_id
        processes_glial_first = processes_glial.drop_duplicates(subset=["component_id"], keep="first")
        processes_glial_last = processes_glial.drop_duplicates(subset=["component_id"], keep="last")
        coord_soma = np.array(glial_i.loc[glial_i["component"] == "soma"][["X", "Y", "Z"]].astype(float))
        radius_soma = float(glial_i.loc[glial_i["component"] == "soma"]["outer_radius"])
        # calculate length of each process
        total_process_lengths = []
        mean_radii = []
        nbr_process = 0
        for i in range(len(processes_glial_first)):
            coord_first = np.array(processes_glial_first.iloc[i][["X", "Y", "Z"]].astype(float))
            distance = np.linalg.norm(coord_first - coord_soma)
            if (distance > radius_soma/2):
                continue
            coord_last = np.array(processes_glial_last.iloc[i][["X", "Y", "Z"]].astype(float))
            length = np.linalg.norm(coord_first - coord_last)
            total_process_lengths.append(length)
            nbr_process += 1
            mean_radii.append(float(processes_glial_first.iloc[i]["outer_radius"]))

        if (nbr_process == 0):
            continue

        all_mean_process_length_gm.append(np.mean(total_process_lengths))
        all_total_process_lengths_gm.append(np.sum(total_process_lengths))
        all_std_process_length_gm.append(np.std(total_process_lengths))
        nbr_processes_gm.append(nbr_process)
        all_radii_gm.append(np.mean(mean_radii))


    data = pd.DataFrame({
        "Glial cells": ["GM"] * len(all_total_process_lengths_gm) + ["WM"] * len(all_total_process_lengths_wm),
        "Length (µm)": all_total_process_lengths_gm + all_total_process_lengths_wm,
        "Nbr processes": nbr_processes_gm + nbr_processes_wm,
        "Mean process length (µm)": all_mean_process_length_gm + all_mean_process_length_wm,
        "Std process length (µm)": all_std_process_length_gm + all_std_process_length_wm,
        "Mean radius (µm)": all_radii_gm + all_radii_wm
    })

    plot_(data, primary = True)
        
def hist_range(ax, unique_groups, expected_range):

    x_position = 2  # place it after WM (0) and GM (1)

    # Add a dummy tick label to x-axis
    xticks = list(unique_groups) + ["Hist. range (GM)"]
    ax.set_xticks(range(len(xticks)))
    ax.set_xticklabels(xticks)

    # Draw the expected range box
    box = patches.Rectangle((x_position - 0.2, expected_range[0]),
                            width=0.4,
                            height=expected_range[1] - expected_range[0],
                            linewidth=1.5,
                            edgecolor='gray',
                            facecolor='pink',
                            alpha=0.5)
    ax.add_patch(box)
    return ax

def plot_(data, primary = False):
    # Create subplots
    fig, axes = plt.subplots(1, 4, figsize=(18, 5), sharex=False)

    # Define color palette as a dictionary
    unique_groups = data["Glial cells"].unique()
    palette_dict = dict(zip(unique_groups, ["red", "blue"]))

    # --- Plot 1: Number of primary processes ---
    sns.stripplot(ax=axes[0], x="Glial cells", y="Nbr processes", data=data, hue = "Glial cells", palette=palette_dict, alpha = 0.5, s = 10 )
    axes[0].set_title("Number of processes",fontsize=14)
    axes[0].set_xlabel("Glial cells",fontsize=14)
    max_value = data["Nbr processes"].max()
    axes[0].set_ylim(0, max_value + max_value/5) 
    if primary:
        axes[0].set_ylabel("Nbr primary processes",fontsize=14)
    else:
        axes[0].set_ylabel("Nbr processes",fontsize=14)
    axes[0].tick_params(axis='both', labelsize=14)


    # Add median lines
    for i, group in enumerate(unique_groups):
        median = data[data["Glial cells"] == group]["Nbr processes"].median()
        axes[0].hlines(median, i - 0.2, i + 0.2, colors='black', linestyles='-', linewidth=5, alpha = 0.5)

    # --- Add "Expected" range as a third box ---
    if primary :
        expected_range = (5, 12) 
    else:
        expected_range = (50, 140)
    axes[0] = hist_range(axes[0], unique_groups, expected_range)

    # --- Plot 2: Total primary process length ---
    sns.stripplot(ax=axes[1], x="Glial cells", y="Length (µm)", hue = "Glial cells", data=data, palette=palette_dict, alpha = 0.5, s = 10 )
    axes[1].set_title("Total process length",fontsize=14)
    axes[1].set_xlabel("Glial cells",fontsize=14)
    max_value = data["Length (µm)"].max()
    axes[1].set_ylim(0, max_value + max_value/5) 
    
    if primary:
        axes[1].set_ylabel("Total primary process length (µm)",fontsize=14)
    else:
        axes[1].set_ylabel("Total process length (µm)",fontsize=14)
    axes[1].tick_params(axis='both', labelsize=14)


    for i, group in enumerate(unique_groups):
        median = data[data["Glial cells"] == group]["Length (µm)"].median()
        axes[1].hlines(median, i - 0.2, i + 0.2, colors='black', linestyles='-', linewidth=5, alpha = 0.5)

    if primary :
        expected_range = (200, 400) 
    else:
        expected_range = (600, 1800)
    axes[1] = hist_range(axes[1], unique_groups, expected_range)

    # --- Plot 3: Mean process length ---
    sns.stripplot(ax=axes[2], x="Glial cells", y="Mean process length (µm)", hue = "Glial cells", data=data, palette=palette_dict, alpha = 0.5, s = 10 )
    axes[2].set_title("Mean process length",fontsize=14) 
    axes[2].set_xlabel("Glial cells", fontsize=14)
    if primary:
        axes[2].set_ylabel("Mean primary process length (µm)", fontsize=14)
    else:
        axes[2].set_ylabel("Mean process length (µm)", fontsize=14)
    axes[2].tick_params(axis='both', labelsize=14)
    max_value = data["Mean process length (µm)"].max()
    axes[2].set_ylim(0, max_value + max_value/5) 


    for i, group in enumerate(unique_groups):
        median = data[data["Glial cells"] == group]["Mean process length (µm)"].median()
        axes[2].hlines(median, i - 0.2, i + 0.2, colors='black', linestyles='-', linewidth=5, alpha = 0.5)

    if primary :
        expected_range = (20, 35) 
    else:
        expected_range = (10, 16)
    axes[2] = hist_range(axes[2], unique_groups, expected_range)

    
    # --- Plot 3: Mean process length ---
    sns.stripplot(ax=axes[3], x="Glial cells", y="Mean radius (µm)", hue = "Glial cells", data=data, palette=palette_dict, alpha = 0.5, s = 10 )
    axes[3].set_title("Mean radius (µm)",fontsize=14) 
    axes[3].set_xlabel("Glial cells", fontsize=14)
    max_value = data["Mean radius (µm)"].max()
    axes[3].set_ylim(0, max_value + max_value/5) 
    if primary:
        axes[3].set_ylabel("Mean primary process radius (µm)", fontsize=14)
    else:
        axes[3].set_ylabel("Mean process radius (µm)", fontsize=14)
    axes[3].tick_params(axis='both', labelsize=14)

    for i, group in enumerate(unique_groups):
        median = data[data["Glial cells"] == group]["Mean radius (µm)"].median()
        axes[3].hlines(median, i - 0.2, i + 0.2, colors='black', linestyles='-', linewidth=5, alpha = 0.5)

    if not primary :
        expected_range = (0.1, 1.7) 
        axes[3] = hist_range(axes[3], unique_groups, expected_range)

    plt.tight_layout()
    plt.show()

    gm_data = data.loc[data["Glial cells"] == "GM"]
    wm_data = data.loc[data["Glial cells"] == "WM"]

    metrics = [
        "Nbr processes",
        "Length (µm)",
        "Mean process length (µm)",
        "Mean radius (µm)"
    ]


    # 1) Run unpaired t-tests for each metric
    for metric in metrics:
        x = gm_data[metric].dropna().values
        y = wm_data[metric].dropna().values

        t_stat, p = ttest_ind(x, y, equal_var=False)
        print("Metric:", metric, " p-value:", p)


    
def radius_with_length(file_path):
    df = read_swc_file(file_path)
    glial_df = df.loc[df["cell_type"] != "axon"]

    interpolated_profiles = []
    distance_samples = np.linspace(5, 50, 200)  # You can adjust max distance and resolution

    for glial in glial_df["cell_id"].unique():
        glial_i = glial_df.loc[glial_df["cell_id"] == glial]
        soma_rows = glial_i.loc[glial_i["component"] == "soma"]
        processes_glial = glial_i.loc[glial_i["component"] == "branch"]

        if soma_rows.empty or processes_glial.empty:
            continue

        coord_soma = np.array(soma_rows[["X", "Y", "Z"]].astype(float)).mean(axis=0)

        distances = []
        radii = []

        for _, sphere in processes_glial.iterrows():
            radius = float(sphere["outer_radius"])
            distance = np.linalg.norm(
                np.array([sphere["X"], sphere["Y"], sphere["Z"]]).astype(float) - coord_soma
            )
            distances.append(distance)
            radii.append(radius)

        if len(distances) < 5:
            continue  # Skip under-sampled cells

        # Sort by distance
        distances, radii = zip(*sorted(zip(distances, radii)))

        # Interpolate for this glial cell
        try:
            interp_func = interp1d(distances, radii, bounds_error=False, fill_value=np.nan)
            interpolated_radii = interp_func(distance_samples)
            interpolated_profiles.append(interpolated_radii)
        except Exception as e:
            print(f"Skipping glial cell {glial} due to interpolation error: {e}")
            continue

    if not interpolated_profiles:
        print("No valid glial cells for interpolation.")
        return

    # Compute mean and std across cells
    interpolated_profiles = np.array(interpolated_profiles)
    mean_radii = np.nanmean(interpolated_profiles, axis=0)
    std_radii = np.nanstd(interpolated_profiles, axis=0)

    return distance_samples, mean_radii, std_radii

def radius_with_length_comparison(file_path_wm, file_path_gm):

    distance_sampled_wm, mean_radii_wm, std_radii_wm = radius_with_length(file_path_wm)
    distance_sampled_gm, mean_radii_gm, std_radii_gm = radius_with_length(file_path_gm)

    # Plot
    plt.figure(figsize=(8, 5))
    plt.plot(distance_sampled_wm, mean_radii_wm, label='Mean radius', color = "blue")
    plt.fill_between(distance_sampled_wm, mean_radii_wm - std_radii_wm, mean_radii_wm + std_radii_wm, alpha=0.3, label='±1 std')
    plt.plot(distance_sampled_gm, mean_radii_gm, label='GM Mean radius', color = "red")
    plt.fill_between(distance_sampled_gm, mean_radii_gm - std_radii_gm, mean_radii_gm + std_radii_gm, alpha=0.3, label='GM ±1 std')
    plt.xlabel("Distance from soma (µm)", fontsize = 14)
    plt.ylabel("Radius (µm)", fontsize = 14)
    plt.tick_params(axis='both', labelsize=14)
    plt.title("Mean radius variation as a function of distance to soma", fontsize = 14)
    plt.grid(True)
    plt.legend().set_visible(False)
    plt.tight_layout()

    plt.show()

def compute_fa_for_glial_cell(file_path):
    from numpy.linalg import eigvalsh

    df = read_swc_file(file_path)
    glial_df = df[df["cell_type"] != "axon"]
    results = []

    for glial in glial_df["cell_id"].unique():
        glial_i = glial_df[glial_df["cell_id"] == glial]
        processes = glial_i[glial_i["component"] == "branch"]
        soma = glial_i[glial_i["component"] == "soma"]

        if soma.empty or processes.empty:
            continue

        # Soma coordinate (mean in case multiple points)
        soma_coord = soma[["X", "Y", "Z"]].astype(float).mean().to_numpy()

        # Use the last point of each branch to represent the end of a process
        process_endpoints = processes.drop_duplicates(subset=["component_id"], keep="last")
        process_firstpoints = processes.drop_duplicates(subset=["component_id"], keep="first")
        vectors = (process_endpoints[["X", "Y", "Z"]].astype(float).to_numpy() - process_firstpoints[["X", "Y", "Z"]].astype(float).to_numpy())
        #vectors = process_endpoints[["X", "Y", "Z"]].astype(float).to_numpy() - soma_coord
        radii = (process_firstpoints["outer_radius"].astype(float).to_numpy() + process_endpoints["outer_radius"].astype(float).to_numpy())/2

        if vectors.shape[0] < 3:
            continue  # Not enough data to compute tensor
        
        # Normalize each vector (unit direction), then scale by its radius
        weighted_vectors = vectors * radii[:, None] * radii[:, None] * np.pi

        # Compute 3x3 covariance matrix of the vectors
        cov = np.cov(weighted_vectors.T)

        # Compute eigenvalues
        eigenvalues = eigvalsh(cov)  # Sorted ascending
        if np.sum(eigenvalues) == 0:
            continue

        # Compute FA
        mean_lambda = np.mean(eigenvalues)
        numerator = np.sqrt(np.sum((eigenvalues - mean_lambda) ** 2)) * np.sqrt(3/2)
        denominator = np.sqrt(np.sum(eigenvalues ** 2))
        FA = numerator / denominator

        results.append({
            "cell_id": glial,
            "FA": FA,
            "N_processes": len(process_endpoints)
        })

    return pd.DataFrame(results)

def FA_comparison(file_path_WM, file_path_GM):
    """
    Compare the Fractional Anisotropy (FA) of glial cells in white matter (WM) and gray matter (GM).
    """

    fa_df_wm = compute_fa_for_glial_cell(file_path_WM)
    fa_df_gm = compute_fa_for_glial_cell(file_path_GM)

    fa_df_wm["Region"] = "WM"
    fa_df_gm["Region"] = "GM"

    fa_df = pd.concat([fa_df_wm, fa_df_gm], ignore_index=True)


    # Extract FA values as numpy arrays
    wm_fa = fa_df_wm["FA"].values
    gm_fa = fa_df_gm["FA"].values

    print(f"Number of glial cells in WM: {len(wm_fa)}, GM: {len(gm_fa)}")

    # Trim to same length
    min_len = min(len(wm_fa), len(gm_fa))
    wm_fa = wm_fa[:min_len].flatten()
    gm_fa = gm_fa[:min_len].flatten()

    # Remove pairs with NaN or Inf
    valid_mask = np.isfinite(wm_fa) & np.isfinite(gm_fa)
    wm_fa = wm_fa[valid_mask]
    gm_fa = gm_fa[valid_mask]

    stat, p_value = ttest_ind(wm_fa, gm_fa)
    print(f"Wilcoxon paired test statistic: {stat}, p-value: {p_value}")

    # Define color palette as a dictionary
    palette_dict = dict(zip(["GM", "WM"], ["red", "blue"]))

    # --- Plot 1: Number of primary processes ---
    sns.stripplot( x="Region", y="FA", data=fa_df, hue = "Region", palette=palette_dict, alpha = 0.5, s = 10 )
    plt.title("FA",fontsize=14)
    plt.xlabel("Glial cells",fontsize=14)
    plt.ylabel("FA",fontsize=14)
    plt.tick_params(axis='both', labelsize=14)


    # Add median lines
    for i, group in enumerate(["WM", "GM"]):
        median = fa_df[fa_df["Region"] == group]["FA"].median()
        plt.hlines(median, i - 0.2, i + 0.2, colors='black', linestyles='-', linewidth=5, alpha = 0.5)

    plt.ylim([0.1,1.2])
    plt.show()


if __name__ == "__main__":

    file_path_GM= "/home/localadmin/Documents/CATERPillar/astrocytes/GM_2.csv"
    file_path_WM = "/home/localadmin/Documents/CATERPillar/astrocytes/WM_8.csv"

    FA_comparison(file_path_WM, file_path_GM)
    total_central_process_length(file_path_WM, file_path_GM)
    total_process_length(file_path_WM, file_path_GM)
    sholl_intersections(file_path_WM, file_path_GM)
    # radius_with_length_comparison(file_path_WM, file_path_GM)