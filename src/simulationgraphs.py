import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def read_swc_file(file_path):
    columns = ["id", "type", "x", "y", "z", "radius", "parent"]
    df = pd.read_csv(file_path, sep=' ', names=columns)

    return df

def radius_file(file_path):
    columns = ["Type1", "ax_id", "Type2", "z",  "R", "R0", "Tortuosity"]
    df = pd.read_csv(file_path, sep='\s+', names=columns)
    return df

def radius_histogram(df):
    df['R'] = pd.to_numeric(df['R'], errors='coerce')
    sns.histplot(data=df, x="R", color='blue', bins=30, kde=True)
    plt.xlabel('Radius')
    plt.ylabel('Frequency')
    plt.title('Radius Histogram')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()


def diameter_variation(df, num_axons=10, max_z=None):
    plt.figure(figsize=(12, 6))  # Adjust the figure size as needed
    df['R'] = pd.to_numeric(df['R'], errors='coerce')

    # Get the first 'num_axons' axon IDs
    first_n_axon_ids = df['ax_id'].unique()[:num_axons]
    df_subset = df[df['ax_id'].isin(first_n_axon_ids)]

    # Filter data until the specific 'z' value if provided
    if max_z is not None:
        max_z_str = str(max_z)
        df_subset = df_subset[df_subset['z'] <= max_z_str]

    # Convert 'sph_id' to categorical to ensure proper x-axis alignment
    df_subset['z'] = pd.Categorical(df_subset['z'])

    # Calculate diameter from radius
    df_subset['Diameter'] = df_subset['R'] * 2

    # Plot smooth curves without individual data points
    sns.lineplot(data=df_subset, x="z", y="Diameter", hue="ax_id", ci=None, legend=False)
    
    plt.xlabel("z")
    plt.ylabel("2r")
    if max_z is not None:
        plt.title(f"Sphere Diameter for First {num_axons-1} Axons (until z={max_z})")
    else:
        plt.title(f"Sphere Diameter for First {num_axons-1} Axons")  # Updated title with the number of axons

    plt.xticks([])  # Remove the x-axis labels
    plt.show()


def create_subplots(df, num_axons=10, max_z=None):
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(12, 6))

    # Diameter variation with the first subplot
    plt.sca(axes[0])
    df['R'] = pd.to_numeric(df['R'], errors='coerce')

    first_n_axon_ids = df['ax_id'].unique()[:num_axons] # get the first 'num_axons' axon IDs
    df_subset = df[df['ax_id'].isin(first_n_axon_ids)]

    # Filter data until the specific 'z' value if provided
    if max_z is not None:
        max_z_str = str(max_z)
        df_subset = df_subset[df_subset['z'] <= max_z_str]

    df_subset['z'] = pd.Categorical(df_subset['z']) # ensures proper x-axis alignment

    df_subset['Diameter'] = df_subset['R'] * 2

    sns.lineplot(data=df_subset, x="z", y="Diameter", hue="ax_id", ci=None, legend=False)
    plt.xlabel("z (µm)")
    plt.ylabel("2r (µm)")
    if max_z is not None:
        plt.title(f"Sphere Diameter for First {num_axons-1} Axons (until z={max_z})")
    else:
        plt.title(f"Sphere Diameter for First {num_axons-1} Axons") 

    plt.xticks([])  # Remove the x-axis labels

    # Radius histogram with the second subplot
    plt.sca(axes[1])

    df['R'] = pd.to_numeric(df['R'], errors='coerce')
    df['Diameter'] = df['R'] * 2
    sns.histplot(data=df, x="Diameter", color='blue', bins=30, kde=True)
    plt.xlabel('Diameter (µm)')
    plt.ylabel('Frequency')
    plt.title('Diameter Histogram')
    plt.xticks(rotation=45)
    
    # Coeff variation for the third subplot
    plt.sca(axes[2])
    df['R'] = pd.to_numeric(df['R'], errors='coerce')
    # Calculate the coefficient of variation for each 'ax_id'
    cv_data = df.groupby('ax_id')['R'].agg(lambda x: (x.std() / x.mean()))

    # Plot the histogram of coefficient of variation
    sns.histplot(data=cv_data, color='blue', bins=30, kde=True)
    plt.xlabel('CV of radius')
    plt.ylabel('Frequency')
    plt.title('CV Histogram')
    plt.xticks(rotation=45)

    # Add an overall title for the entire figure
    plt.suptitle('Substrate Analysis', fontsize=16)

    # Adjust the layout to prevent overlapping of titles and labels
    plt.tight_layout()

    # Display the figure
    plt.show()

def plot_tortuosity_(df):
    df['Tortuosity'] = pd.to_numeric(df['Tortuosity'], errors='coerce')
    g = sns.JointGrid(data=df, x="R0", y="Tortuosity", space=0)
    g.plot_joint(sns.kdeplot,
             fill=True,
              cmap="rocket")
    g.plot_marginals(sns.histplot, color="#03051A", alpha=1, bins=25)
    plt.show()

def tortuosity_plot(df):
    df['radius'] = pd.to_numeric(df['R'], errors='coerce')
    df['Tortuosity'] = pd.to_numeric(df['Tortuosity'], errors='coerce')

    g = sns.JointGrid(data=df, x="radius", y="Tortuosity", space=0)
    g.plot_joint(sns.kdeplot, fill=True, cmap="rocket")
    g.plot_marginals(sns.histplot, color='blue', alpha=1, bins=25)

    plt.xlabel('radius')
    plt.ylabel('Tortuosity')
    plt.title('Tortuosity plot')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()



if __name__ == "__main__":
    swc_file_path = "/Users/melina/Desktop/EPFL/BachelorProject/Sim_Growth/axon_simulation.swc"
    radius_file_path = "/Users/melina/Desktop/EPFL/BachelorProject/Sim_Growth/radius.swc"
    # graph = read_swc_file(swc_file_path)
    graph2 = radius_file(radius_file_path)
    tortuosity_plot(graph2)
    create_subplots(graph2, 149)