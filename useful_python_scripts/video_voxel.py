import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from matplotlib import animation
from simulationgraphs import read_swc_file
import math
def create_video_with_soma_axon_process(
    input_swc, output_mp4,
    fps=5, resolution=(12, 8),
    group_size=10
):
    """
    Animation:
      Frame 0: show all CellSoma
      Next frames: add axon spheres in groups
      Then frames: add Process spheres in groups
    """
    # 1) Load & split
    df = read_swc_file(input_swc)
    df_soma    = df[(df['type'] == 'CellSoma') & (df["id_ax"] == 0)].reset_index(drop=True)
    df_axon    = df[(df['type'] == 'axon') & (df["id_ax"] <10)].reset_index(drop=True)
    df_process = df[(df['type'] == 'Process') & (df["id_ax"] == 0)].reset_index(drop=True)

    # 2) Unit-sphere mesh
    u = np.linspace(0, 2*np.pi, resolution[0])
    v = np.linspace(0,   np.pi,  resolution[1])
    X0 = np.outer(np.cos(u), np.sin(v))
    Y0 = np.outer(np.sin(u), np.sin(v))
    Z0 = np.outer(np.ones_like(u), np.cos(v))

    # 3) Figure setup
    fig = plt.figure(figsize=(8,6))
    ax  = fig.add_subplot(111, projection='3d')
    ax.set_box_aspect([1,1,1])

    # 4) Global bounds
    all_coords = np.vstack([
        df_soma[['x','y','z']].values,
        df_axon[['x','y','z']].values,
        df_process[['x','y','z']].values
    ])
    all_radii = np.hstack([
        df_soma['Rout'].values,
        df_axon['Rout'].values,
        df_process['Rout'].values
    ])
    mins = all_coords - all_radii[:,None]
    maxs = all_coords + all_radii[:,None]
    mn, mx = mins.min(), maxs.max()
    ax.set_xlim(mn, mx); ax.set_ylim(mn, mx); ax.set_zlim(mn, mx)

    # 5) Create & hide artists
    soma_artists = []
    for _, r in df_soma.iterrows():
        X = X0*r.Rout + r.x
        Y = Y0*r.Rout + r.y
        Z = Z0*r.Rout + r.z
        surf = ax.plot_surface(X,Y,Z, color='orange', alpha=1, linewidth=0, shade=True)
        surf.set_visible(False)
        soma_artists.append(surf)

    axon_artists = []
    for _, r in df_axon.iterrows():
        X = X0*r.Rout + r.x
        Y = Y0*r.Rout + r.y
        Z = Z0*r.Rout + r.z
        surf = ax.plot_surface(X,Y,Z, color='blue', alpha=0.5, linewidth=0, shade=True)
        surf.set_visible(False)
        axon_artists.append(surf)

    process_artists = []
    for _, r in df_process.iterrows():
        X = X0*r.Rout + r.x
        Y = Y0*r.Rout + r.y
        Z = Z0*r.Rout + r.z
        surf = ax.plot_surface(X,Y,Z, color='orange', alpha=0.5, linewidth=0, shade=True)
        surf.set_visible(False)
        process_artists.append(surf)

    # 6) Frame counts
    n_axon_frames    = math.ceil(len(axon_artists)    / group_size)
    n_process_frames = math.ceil(len(process_artists) / group_size)
    total_frames = 1 + n_axon_frames + n_process_frames

    # 7) Update function
    def update(f):
        changed = []
        if f == 0:
            # show all somas
            for art in soma_artists:
                art.set_visible(True)
                changed.append(art)
        elif 1 <= f <= n_axon_frames:
            # add next group of axons
            i = f - 1
            start = i*group_size
            end   = min(start+group_size, len(axon_artists))
            for art in axon_artists[start:end]:
                art.set_visible(True)
                changed.append(art)
        else:
            # add next group of processes
            j = f - 1 - n_axon_frames
            start = j*group_size
            end   = min(start+group_size, len(process_artists))
            for art in process_artists[start:end]:
                art.set_visible(True)
                changed.append(art)

        #ax.set_title(f"Frame {f+1}/{total_frames}", fontsize=12)
        return changed

    # 8) Animate & save
    ani = animation.FuncAnimation(
        fig, update,
        frames=total_frames,
        blit=True,
        interval=1000/fps,
        repeat=False
    )
    ani.save(output_mp4, writer='ffmpeg', fps=fps)
    plt.close(fig)
    print("Video saved to", output_mp4)

# Example usage:
create_video_with_soma_axon_process('/home/localadmin/Documents/MCDS/Permeable_MCDS/output/incoherent_blood_flow/Voxel.csv', 'spheres_animation.mp4', fps=10, group_size=100)
