import os
import matplotlib.pyplot as plt
import numpy as np
def extract_data_from_txt(file_path):
    """Extract Duration and Axon icvf from a single .txt file."""
    duration = None
    axon_icvf = None
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("Duration"):
                duration = float(line.split()[1])
            elif line.startswith("Axon icvf"):
                axon_icvf = float(line.split()[2])
    return duration, axon_icvf

def plot_icvf_vs_duration(folder_path):
    durations = []
    icvfs = []

    for filename in os.listdir(folder_path):
        if filename.endswith("info.txt"):
            file_path = os.path.join(folder_path, filename)
            duration, icvf = extract_data_from_txt(file_path)
            if duration is not None and icvf is not None:
                durations.append(duration / 60)  # Convert seconds to minutes
                icvfs.append(icvf)

    fig, ax = plt.subplots(figsize=(9, 6))
    fig.patch.set_facecolor('black')          # Full figure background
    ax.set_facecolor('black')                 # Plot background

    # Scatter plot
    ax.scatter(icvfs, durations,
               color='white',                 # White points
               edgecolor='white',
               s=80, alpha=0.9)

    # Labels and title
    ax.set_xlabel('Axon f', fontsize=14, weight='bold', color='white')
    ax.set_ylabel('Run-time [min]', fontsize=14, weight='bold', color='white')
    ax.set_title('Axon f vs Run-time', fontsize=16, weight='bold', color='white')

    # Ticks color
    ax.tick_params(axis='x', colors='white', labelsize=12)
    ax.tick_params(axis='y', colors='white', labelsize=12)

    # Grid
    ax.grid(True, linestyle='--', linewidth=0.5, color='gray', alpha=0.3)

    #ax.set_ylim(-0.5, 3)
    plt.tight_layout()
    plt.show()



if __name__ == "__main__":
    # Replace with the path to your folder containing .txt files
    #folder_path = "/home/localadmin/Documents/CATERPillar/runtime-packing/try2_30/"
    folder_path = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/SMI_pred/complex_axons/"
    
    # Call the function to plot ICVF vs Duration
    plot_icvf_vs_duration(folder_path)
