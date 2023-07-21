import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def read_swc_file(file_path):
    columns = ["id", "type", "x", "y", "z", "radius", "parent"]
    df = pd.read_csv(file_path, sep=' ', names=columns)
    return df


def radius_histogram(df):
    sns.histplot(data=df, x="radius", color='blue', bins=10, kde=False)  # Use sns.histplot() instead of histoplot
    plt.xlabel('Radius')
    plt.ylabel('Frequency')
    plt.title('Radius Histogram')
    plt.xticks(rotation=45)  # Rotates the x-axis labels for better readability
    plt.tight_layout()  # Ensures the plot components fit within the figure area
    plt.show()


def radius_file(file_path):
    columns = ["ax_id", "Type", "sph_id", "Type2",  "R"]
    df = pd.read_csv(file_path, sep='\s+', comment='#', names=columns)
    return df



def radius_variation(df):
    plt.figure(figsize=(12, 6))  # Adjust the figure size as needed

    sns.lineplot(data=df, x="sph_id", y="R", hue="ax_id", marker="o")

    plt.xlabel("Sphere ID")
    plt.ylabel("Sphere Radius")
    plt.title("Sphere Radius for Each Axon")
    plt.legend(title="Axon ID")
    plt.show()

if __name__ == "__main__":
    swc_file_path = "/Users/melina/Desktop/EPFL/BachelorProject/Sim_Growth/axon_simulation.swc"
    # radius_file_path = "/Users/melina/Desktop/EPFL/BachelorProject/Sim_Growth/radius.swc"
    graph = read_swc_file(swc_file_path)
    # graph2 = radius_file(radius_file_path)
    radius_histogram(graph)
    # radius_variation(graph2)