"""
Plots mean/std of glial cell process lengths, primary and secondary
processes shown separately.

"Length" here is the real (arc) length along a process -- the sum of
distances between its consecutive spheres, in growth order -- not the
straight-line distance between its first and last sphere. Glial processes
wander, so the two can differ meaningfully; arc length is also what
mean_process_length/std_process_length actually parameterize internally
(see growPrimaryBranch/growSecondaryBranch in src/grow_glial_cells.cpp), so
it's the metric that's directly comparable to those configured values.

A process is classified as "primary" if its first sphere sits on (or very
near) its own cell's soma surface, "secondary" otherwise (it attaches
partway along another process instead). This mirrors the geometric
convention total_central_process_length() in astrocytes.py already uses,
and is more robust here than trusting the CSV's parent_component_id alone:
branch index 0 doubles as both "this branch's own id" and the sentinel
CaterpillarGrowth::create_SWC_file's resolveParentBranches uses for
"attaches directly to the soma", so a secondary process that happens to
attach to branch 0 would be indistinguishable from a true primary process
by parent_component_id alone.

Note: the CSV doesn't record which glial population (pop1/pop2/pop3) a cell
belongs to, so this pools every glial cell in the file together -- it can't
report primary/secondary stats per population.
"""

import numpy as np
import matplotlib.pyplot as plt
from simulationgraphs import read_swc_file


def compute_process_lengths(file_path, primary_soma_multiplier=1.5):
    """
    Returns (primary_lengths, secondary_lengths): lists of each glial
    process's own real (arc) length in microns, split by whether it
    attaches directly to its cell's soma (primary) or partway along
    another process (secondary).

    primary_soma_multiplier: a process's first sphere is considered "on the
    soma" if it's within this many soma radii of the soma center. Primary
    processes start essentially on the soma surface (~1 soma radius away);
    secondary processes attach much farther out along another process, so
    this cleanly separates the two even with some margin for slack.
    """
    df = read_swc_file(file_path)
    glial_df = df.loc[df["cell_type"] == "glial_cell"]

    primary_lengths = []
    secondary_lengths = []

    for cell_id, cell_df in glial_df.groupby("cell_id"):
        soma_rows = cell_df.loc[cell_df["component"] == "soma"]
        if soma_rows.empty:
            continue
        soma_center = soma_rows.iloc[0][["X", "Y", "Z"]].to_numpy(dtype=float)
        soma_radius = float(soma_rows.iloc[0]["outer_radius"])

        branches = cell_df.loc[cell_df["component"] == "branch"]
        for component_id, branch_df in branches.groupby("component_id"):
            coords = branch_df[["X", "Y", "Z"]].to_numpy(dtype=float)
            if len(coords) < 2:
                continue

            segment_lengths = np.linalg.norm(np.diff(coords, axis=0), axis=1)
            length = float(np.sum(segment_lengths))

            distance_first_to_soma = np.linalg.norm(coords[0] - soma_center)
            is_primary = distance_first_to_soma <= primary_soma_multiplier * soma_radius

            if is_primary:
                primary_lengths.append(length)
            else:
                secondary_lengths.append(length)

    return primary_lengths, secondary_lengths


def plot_process_lengths(primary_lengths, secondary_lengths, title=None):
    """
    Plots primary and secondary process length distributions side by side,
    each annotated with mean +/- std, and prints the same to stdout.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for ax, lengths, label, color in zip(
        axes,
        [primary_lengths, secondary_lengths],
        ["Primary", "Secondary"],
        ["tab:blue", "tab:orange"],
    ):
        if len(lengths) == 0:
            ax.set_title(f"{label} processes (n=0)")
            print(f"{label} processes: n=0")
            continue

        mean = np.mean(lengths)
        std = np.std(lengths)

        ax.hist(lengths, bins=30, color=color, edgecolor="black", alpha=0.8)
        ax.axvline(mean, color="black", linewidth=2, label=f"mean = {mean:.2f} µm")
        ax.axvline(mean - std, color="black", linestyle="--", linewidth=1, label=f"std = {std:.2f} µm")
        ax.axvline(mean + std, color="black", linestyle="--", linewidth=1)
        ax.set_title(f"{label} processes (n={len(lengths)})")
        ax.set_xlabel("Process length (µm)")
        ax.set_ylabel("Count")
        ax.legend()

        print(f"{label} processes: n={len(lengths)}  mean={mean:.3f} µm  std={std:.3f} µm")

    if title:
        fig.suptitle(title)
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    file_path = "/home/localadmin/Documents/CATERPillar/tests/Voxel.csv"

    primary_lengths, secondary_lengths = compute_process_lengths(file_path)
    plot_process_lengths(primary_lengths, secondary_lengths, title=file_path)
