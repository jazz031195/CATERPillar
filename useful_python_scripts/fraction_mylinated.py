import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
from pathlib import Path
from scipy import stats
import random
from simulationgraphs import get_spheres_array, read_swc_file
from reduce_myelinated_axons_in_substrate import compute_icvf_axons

if __name__ == "__main__":

    # Load the astrocyte data
    file = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/SMI_pred/permeable_axons/axons_0.5.swc"
    df= read_swc_file(file)
    cell_ids = df['ax_id'].unique()
    factor = 4
    limit = 150

    myelinated_icvf = 0
    unmyelinated_icvf = 0
    for cell_id in cell_ids:
        print(f"Computing ICVF for cell {cell_id}")
        df_cell = df[df['ax_id'] == cell_id]
        if float(df_cell.iloc[0]['Rout']) !=  float(df_cell.iloc[0]['Rin']):
            myelinated_icvf += compute_icvf_axons(df_cell, limit,factor)
        else:
            unmyelinated_icvf += compute_icvf_axons(df_cell, limit,factor)
    
    print(f"Myelinated ICVF: {myelinated_icvf}")
    print(f"Unmyelinated ICVF: {unmyelinated_icvf}")

