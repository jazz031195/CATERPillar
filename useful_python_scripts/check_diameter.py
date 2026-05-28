import numpy as np
import os
from simulationgraphs import read_swc_file
import matplotlib.pyplot as plt
import matplotlib.cm as cm
# Example usage
def plot_diameters_from_folder():

    directory = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/SMI_pred/"

    folders = ["bead_freq", "complex_axons", "p2", "undulations"]
    variables = ["COV\n", "icvf", "p2", "Undulation"]

    # Create a 2x2 grid of subplots since there are 4 folders
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(14, 10))
    axes = axes.flatten() # Flatten the 2D array of axes for easy zipping

    for ax, folder, variable in zip(axes, folders, variables):
        path = directory + folder + "/combined/"

        # find all folders
        if not os.path.exists(path):
            print(f"Directory not found: {path}")
            continue
            
        subfolders = [f.path for f in os.scandir(path) if f.is_dir()]

        values = []
        all_radii = []
        
        for subfolder in subfolders:
            smi_output_path = subfolder + "/SMI_output_SNR_None_None.txt"
            
            if not os.path.exists(smi_output_path):
                continue
                
            with open(smi_output_path, "r") as f:
                lines = f.readlines()
                for i, line in enumerate(lines):
                    
                    if i == 0:
                        # find index of variable
                        try:
                            index = line.split(" ").index(variable)
                            print(f"Found {variable} at index {index}")
                        except ValueError:
                            index = -1
                    elif i == 1 and index != -1:
                        value = line.split(" ")[index]
                        values.append(float(value))
                    else:
                        break

            # Use basename so we don't append a full path to another full path
            subfolder_name = os.path.basename(subfolder)
            voxel_path = directory + folder + f"/{subfolder_name}.csv"

            # check if the file exists
            if not os.path.exists(voxel_path):
                voxel_path = directory + folder + f"/{subfolder_name}.swc"
                if not os.path.exists(voxel_path):
                    print(f"File not found for {subfolder_name}")
                    continue

            print(f"Processing {voxel_path} with {variable.strip()} = {values[-1]}")

            # Assuming read_swc_file_new is defined elsewhere in your script
            df = read_swc_file(voxel_path)
            
            df = df.groupby(by="cell_id").mean(numeric_only=True)
            radii = df["outer_radius"].values

            print(radii)
            all_radii.append(radii)

        # Sort values and all_radii together based on values (lowest to highest)
        values, all_radii = zip(*sorted(zip(values, all_radii)))

        # Generate a smooth sequence of colors based on the number of distributions
        colors = cm.viridis(np.linspace(0, 0.9, len(all_radii)))

        # Plot all histograms for this specific folder on its assigned subplot
        for r, value, color in zip(all_radii, values, colors):
            
            if variable == "icvf":
                variable = "f"
            # Remove the weights calculation and just add density=True
            ax.hist(r, bins=50, density=True, histtype='step', linewidth=2.5, 
                    alpha=0.9, color=color, label=f"{variable.strip()} = {value:.2f}")
            
        if folder == "bead_freq":
            title = "Beading Strength (COV)"
        elif folder == "complex_axons":
            title = "Axonal Water Fraction (f)"
        elif folder == "p2":
            title = "Fibre Dispersion (p2)"
        elif folder == "undulations":
            title = "Undulation Strength"
        # Format the subplot
        ax.set_title(f"Variation of: {title}", fontsize=14, fontweight='bold')
        ax.set_xlabel("Radius (μm)", fontsize=12)
        ax.set_ylabel("Probability Density", fontsize=12)
        ax.legend(loc='upper right')

        del values, all_radii


    # Adjust layout so labels and titles don't overlap, then display
    plt.tight_layout()
    plt.show()


def plot_radii_distribution(file_path, cell_type = None):
    """
    Reads a single CSV or SWC file, calculates the mean outer radius per cell, 
    and plots the probability density distribution.
    """
    if not os.path.exists(file_path):
        print(f"Error: File not found at {file_path}")
        return

    print(f"Processing: {os.path.basename(file_path)}")

    # 1. Load the data
    df = read_swc_file(file_path)

    if (cell_type is not None):
        df = df.loc[df["cell_type"]== cell_type]

    # 2. Extract the mean outer radius per cell
    df_grouped = df.groupby(by="cell_id").mean(numeric_only=True)
    radii = df_grouped["outer_radius"].values

    if len(radii) == 0:
        print("Warning: No radii data found after grouping.")
        return

    mean = np.mean(radii)
    print("Mean radius:" ,mean)
    # 3. Plot the distribution
    plt.figure(figsize=(8, 6))
    plt.hist(radii, bins=50, density=True, histtype='step', 
             linewidth=2.5, color='#2980b9', alpha=0.9)

    # 4. Format the plot
    plt.title("Distribution of Mean Cell Radii", fontsize=15, fontweight='bold')
    plt.xlabel("Radius (μm)", fontsize=13)
    plt.ylabel("Probability Density", fontsize=13)
    
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.show()

# Example usage:
if __name__ == "__main__":
    target_csv = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/incoherent_blood_flow/State4.csv"
    plot_radii_distribution(target_csv, cell_type = "blood_vessel")