import numpy as np
import pandas as pd
from simulationgraphs import read_swc_file


def _orientation_tensor_from_pairs(positions, target_sep, tol=2.0, stride=1, weight="chord"):
    """
    Build an orientation tensor from pairs of points along a 3D curve.

    positions : (N,3) array
    target_sep : float
        Desired *arc-length* separation (in same units as positions). If you prefer
        Euclidean, set use_arc_length=False below.
    tol : float
        Acceptable tolerance around target_sep (arc-length).
    stride : int
        Step between starting indices to reduce computation.
    weight : {'chord','uniform','arc'}
        'chord'   -> w = |p_j - p_i|
        'uniform' -> w = 1
        'arc'     -> w = arc-length separation

    Returns
    -------
    M : (3,3) ndarray
    W : float
    """

    # Precompute arc-length along the polyline
    diffs = np.diff(positions, axis=0)
    seglen = np.linalg.norm(diffs, axis=1)
    s = np.zeros(len(positions))
    s[1:] = np.cumsum(seglen)

    M = np.zeros((3,3), dtype=float)
    W = 0.0

    # For each start point, find an end point with arc-length close to target
    for i in range(0, len(positions) - 1, stride):
        si = s[i]
        # binary search (or argmin) over j>i to get closest arc-length
        # restrict candidates by tolerance window to avoid excessive curvature bias
        lo = np.searchsorted(s, si + max(target_sep - tol, si), side='left')
        hi = np.searchsorted(s, si + (target_sep + tol), side='right')
        if lo >= len(s):
            continue
        j_candidates = np.arange(lo, min(hi, len(s)))
        if j_candidates.size == 0:
            continue

        # choose the j with arc-length closest to target
        j = j_candidates[np.argmin(np.abs(s[j_candidates] - (si + target_sep)))]
        v = positions[j] - positions[i]
        Lchord = np.linalg.norm(v)
        if Lchord == 0:
            continue
        u = v / Lchord

        if weight == "chord":
            w = Lchord
        elif weight == "uniform":
            w = 1.0
        elif weight == "arc":
            w = float(s[j] - s[i])
        else:
            raise ValueError("weight must be 'chord','uniform', or 'arc'")

        M += w * np.outer(u, u)
        W += w

    return M, W


def compute_p2(positions, target_sep, tol=2.0, stride=1, weight="chord"):
    """
    Compute orientation tensor for one axon, then return (M, W).

    axon_df : DataFrame with columns ['x','y','z']
    target_sep : desired *arc-length* separation between paired points
    tol : tolerance around target_sep (arc-length)
    stride : step between starting indices
    weight : 'chord' | 'uniform' | 'arc'
    """
    
    return _orientation_tensor_from_pairs(positions, target_sep, tol, stride, weight)


def compute_p2_axons(path_df, target_sep, tol=2.0, stride=1, weight="chord"):
    """
    Aggregate over all axons in an SWC-derived DataFrame and return p2.

    Returns np.nan if no valid pairs found.
    """

    str_id = "cell_id"
    str_type = "cell_type"
    str_x = "X"
    str_y = "Y"
    str_z = "Z"

    df = read_swc_file(path_df)
    df = df.loc[df[str_type] == "axon"]  # axon only
    axon_ids = df[str_id].unique()
    if axon_ids.size == 0:
        return np.nan

    M_tot = np.zeros((3,3), dtype=float)
    W_tot = 0.0

    for ax_id in axon_ids:
        ax_df = df.loc[df[str_id] == ax_id]
        positions = ax_df[[str_x, str_y, str_z]].to_numpy(dtype=float)
        M, W = compute_p2(positions, target_sep, tol, stride, weight)
        if W > 0:
            M_tot += M
            W_tot += W

    if W_tot == 0:
        return np.nan

    # Normalize to the orientation tensor
    M_tot /= W_tot

    # Symmetrize defensively (numerical noise)
    M_tot = 0.5 * (M_tot + M_tot.T)

    # Eigen-decomp (ascending order)
    vals, vecs = np.linalg.eigh(M_tot)
    c2 = float(vals[-1])  # principal eigenvalue

    # Map to p2, clamp small numerical drift
    p2 = 1.5 * c2 - 0.5
    p2 = float(np.clip(p2, -0.5, 1.0))
    return p2

# Example usage
if __name__ == "__main__":

    # Path to the sphere file
    sphere_file = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/incoherent_blood_flow/State1.csv"

    factor = 4
    limit = 15

    # Compute the C2 metric
    p2 = compute_p2_axons(sphere_file, target_sep = limit)
    print(f"p2 Metric: {p2}")
