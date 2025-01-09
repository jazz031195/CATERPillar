import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
from pathlib import Path
from scipy import stats
import random
from simulationgraphs import get_spheres_array, read_swc_file

def shrink_swc(path, shrink_perc):
    # create an array with dwi values
    new_lines = []
    i = 0
    with open(path) as file:
        for line in file:
     
            if i != 0:
                # read each line in the file
                # convert the line to a float and append it to the signal array
                line_values =  line.strip().split()
                Rout = float(line_values[8])
                Rin = float(line_values[7])
                new_Rout = shrink_radius(shrink_perc, Rout)
                new_Rin = shrink_radius(shrink_perc, Rin)
                new_line = line_values
                new_line[7] = new_Rin
                new_line[8] = new_Rout 
                new_lines.append(new_line)
            else:
                new_lines.append(line.strip().split())
            i = i+ 1
    return new_lines

def swell_swc(path, shrink_perc):
    # create an array with dwi values
    new_lines = []
    i = 0
    with open(path) as file:
        for line in file:
     
            if i != 0:
                # read each line in the file
                # convert the line to a float and append it to the signal array
                line_values =  line.strip().split()
                Rout = float(line_values[8])
                Rin = float(line_values[7])
                new_Rout = swell_radius(shrink_perc, Rout)
                new_Rin = swell_radius(shrink_perc, Rin)
                new_line = line_values
                new_line[7] = new_Rin
                new_line[8] = new_Rout 
                new_lines.append(new_line)
            else:
                new_lines.append(line.strip().split())
            i = i+ 1
    return new_lines

def write_new_swc_file(path, new_lines):
    with open(path, 'w') as file:
        for line in new_lines:
            file.write(f"{line[0]} {line[1]} {line[2]} {line[3]} {line[4]} {line[5]} {line[6]} {line[7]} {line[8]} {line[9]}\n")
    
    
def shrink_radius(perc, swollen_radius):
    return swollen_radius/np.sqrt(1+perc)

def swell_radius(perc, initial_radius):
    return initial_radius*np.sqrt(1+perc)



if __name__ == "__main__":

    # shrink by 1%
    path = "/home/localadmin/Documents/CATERPillar/arthurs_analysis/voxel.swc"
    shrink_perc = 0.01
    new_lines = shrink_swc(path, shrink_perc)
    new_path = "/home/localadmin/Documents/CATERPillar/arthurs_analysis/voxel_0.swc"
    write_new_swc_file(new_path, new_lines)

    path = "/home/localadmin/Documents/CATERPillar/arthurs_analysis/voxel_0.swc"

    # swell by 0.25%
    new_lines = swell_swc(path, 0.0025)
    new_path = "/home/localadmin/Documents/CATERPillar/arthurs_analysis/voxel_0.25.swc"
    write_new_swc_file(new_path, new_lines)

    # swell by 0.5%
    new_lines = swell_swc(path, 0.005)
    new_path = "/home/localadmin/Documents/CATERPillar/arthurs_analysis/voxel_0.5.swc"
    write_new_swc_file(new_path, new_lines)

    # swell by 0.75%
    new_lines = swell_swc(path, 0.0075)
    new_path = "/home/localadmin/Documents/CATERPillar/arthurs_analysis/voxel_0.75.swc"
    write_new_swc_file(new_path, new_lines)
