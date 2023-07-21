import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def read_swc_file(file_path):
    columns = ["id", "type", "x", "y", "z", "radius", "parent"]
    df = pd.read_csv(file_path, sep=' ', comment='#', names=columns)
    return df


def radius_file(file_path):
    columns = ["ax_id", "Type", "sph_id", "Type2",  "R"]
    df = pd.read_csv(file_path, sep='\s+', comment='#', names=columns)
    return df


def radius_histogram(df):
    radius_data = df["radius"]
    print(df)
    sns.histplot(radius_data, kde=True, color='blue')
    plt.xlabel('Radius')
    plt.ylabel('Frequency')
    plt.title('Radius Histogram')
    plt.show()


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
    radius_file_path = "/Users/melina/Desktop/EPFL/BachelorProject/Sim_Growth/radius.swc"
    graph = read_swc_file(swc_file_path)
    # graph2 = radius_file(radius_file_path)
    radius_histogram(graph)
    # radius_variation(graph2)
