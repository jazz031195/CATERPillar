import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import binned_statistic_2d

def rotation_matrix_from_vectors(vec1, vec2):
    """
    Find the rotation matrix that aligns vec1 to vec2.
    :param vec1: A 3D "source" vector.
    :param vec2: A 3D "destination" vector.
    :return: A 3x3 rotation matrix that aligns vec1 with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]],
                     [v[2], 0, -v[0]],
                     [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2 + 1e-10))
    return rotation_matrix

def generate_random_points_on_sphere(std, num_points=1000):
    """
    Generate uniformly distributed points on the surface of a sphere.
    """
    points = []
    for _ in range(num_points):
        # Generate random theta and phi values for spherical coordinates
        theta = np.random.normal(0, std)
        cos_phi = np.random.normal(0, std)
        phi = np.arccos(cos_phi)

        # Ensure theta is within the range [0, 2*pi] and phi is within the range [0, pi]
        theta = theta % (2 * np.pi)
        phi = np.clip(phi, 0, np.pi)


        # Spherical to Cartesian conversion
        x = np.sin(phi) * np.cos(theta)
        y = np.sin(phi) * np.sin(theta)
        z = np.cos(phi)

        points.append([x, y, z])

    return np.array(points)

def apply_bias_toward_target(points, target):
    """
    Apply a bias to the points to cluster them toward the target direction.
    This is done by rotating random points on a sphere toward the desired direction.
    """
    # Step 1: Compute the rotation matrix from [0, 0, 1] to the target
    R = rotation_matrix_from_vectors(np.array([1, 0, 0]), target)

    # Step 2: Rotate all generated points using the computed rotation matrix
    rotated_points = points.dot(R.T)

    return rotated_points

def plot_points_on_sphere(points, target):
    """
    Plot the generated points on a 3D sphere with specific visual adjustments.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Set the background color to black
    ax.set_facecolor('black')
    fig.patch.set_facecolor('black')

    # Remove the axes and grid
    ax.set_axis_off()

    # Extract x, y, z coordinates from points
    x_vals = points[:, 0]
    y_vals = points[:, 1]
    z_vals = points[:, 2]

    # Scatter plot of the points on the sphere in purple
    ax.scatter(x_vals, y_vals, z_vals, c='purple', marker='o', alpha = 0.3)

    # Plot the target point in gold
    ax.scatter(target[0], target[1], target[2], c='gold', marker='o', s=100)

    # Set the aspect ratio to be equal to make the sphere look correct
    ax.set_box_aspect([1, 1, 1])

    # Plot the sphere
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, color='white', alpha=0.3)

    plt.show()

def plot_density_on_sphere_with_target(points, target):
    """
    Plot the density of points on a sphere as a heatmap, with a target point highlighted.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Set the background color to black
    ax.set_facecolor('black')
    fig.patch.set_facecolor('black')

    # Remove the axes and grid
    ax.set_axis_off()

    # Map the points to spherical coordinates
    x_vals = points[:, 0]
    y_vals = points[:, 1]
    z_vals = points[:, 2]

    phi = np.arctan2(y_vals, x_vals)  # Longitude
    theta = np.arccos(z_vals / np.linalg.norm(points, axis=1))  # Colatitude

    # Create a grid for the sphere
    u = np.linspace(0, 2 * np.pi, 200)  # Longitude
    v = np.linspace(0, np.pi, 100)  # Colatitude
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones(np.size(u)), np.cos(v))

    # rid phi on nans
    phi = np.nan_to_num(phi)
    # rid theta on nans
    theta = np.nan_to_num(theta)

    print(len(phi), len(theta)) 
    # Calculate the density of points in spherical bins
    density, _, _, _ = binned_statistic_2d(phi, theta, None, statistic='count', bins=[200, 100])

    # Normalize density to [0, 1] for coloring
    density_normalized = density / np.max(density)
    colors = plt.cm.plasma(density_normalized)

    # Reshape colors to match the grid dimensions
    facecolors = colors[: x.shape[0], : x.shape[1], :]

    # Plot the sphere with heatmap colors
    ax.plot_surface(x, y, z, rstride=1, cstride=1, facecolors=facecolors, linewidth=0, alpha=1)

    # Highlight the target point
    ax.scatter(target[0], target[1], target[2], c='lime', marker='o', s=100)

    plt.show()

# Main function to generate and plot points
if __name__ == "__main__":
    desired_target = np.array([1, 0, 1])  # Change this to any desired point in space

    # Step 1: Generate random points on the surface of a sphere
    points = generate_random_points_on_sphere(std= 0.2, num_points=1000)

    # Step 2: Apply bias toward the desired target
    biased_points = apply_bias_toward_target(points,  target=desired_target)

    # Step 3: Plot the biased points with the desired target
    plot_points_on_sphere(biased_points, desired_target)
