import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
from pathlib import Path
from scipy import stats
import random
from simulationgraphs import get_spheres_array, read_swc_file

def get_all_spheres(path):
    # create an array with dwi values
    spheres = []
    i = 0
    with open(path) as file:
        for line in file:
            if i != 0:
            # read each line in the file
            # convert the line to a float and append it to the signal array
                line_values =  line.strip().split()
                spheres.append(line_values)
          
            i = i+ 1
    return np.array(spheres)

def get_all_axons(path):
    spheres = get_all_spheres(path)
    axons = []
    spheres_of_axon=[]
    last_ax_id = -1 # id of last sphere to be appended
    for s, sphere in enumerate(spheres):

        ax_id = sphere[0]
        if s!=0 and last_ax_id != ax_id:
            #print ("new axon")
            axons.append(spheres_of_axon)
            spheres_of_axon=[]

        spheres_of_axon.append(sphere)
        last_ax_id = ax_id

    axons.append(spheres_of_axon)
    return axons

def shrink_radius(perc, swollen_radius):
    return swollen_radius/np.sqrt(1+perc)

def swell_radius(perc, initial_radius):
    return initial_radius*np.sqrt(1+perc)

def create_boolean_list(size, true_percentage):
    if not (0 <= true_percentage <= 100):
        raise ValueError("Percentage must be between 0 and 100.")

    true_count = int(size * true_percentage / 100)
    false_count = size - true_count

    boolean_list = [True] * true_count + [False] * false_count
    random.shuffle(boolean_list)

    return boolean_list

def get_first_line_from_file(file_path):
    with open(file_path, 'r') as file:
        first_line = file.readline()
    return first_line  

def write_list_of_lists_to_file(file_path, lines):
    with open(file_path, 'w') as file:
        for e,line in enumerate(lines):
            if e == 0 :
                str_line = ''.join(line) 
                file.write(str_line)
            else:
                str_line = ' '.join(line) + '\n'
                file.write(str_line)

def adjust_percentage_of_true(input_list,  target_percentage):

    if not all(isinstance(item, bool) for item in input_list):
        raise ValueError("The input list should contain only boolean values.")

    if not (0 <= target_percentage <= 100):
        raise ValueError("Percentages must be between 0 and 100.")

    true_indices = [i for i, val in enumerate(input_list) if val]
    true_count = len(true_indices)
    target_true_count = round(len(input_list) * target_percentage / 100)
     #print("target_true_count :",target_true_count)
    #print("true_count :",true_count)
    if target_true_count > true_count:
        
        raise ValueError("The target percentage must be smaller than the initial percentage.")

    random.shuffle(true_indices)
    new_true_indices = true_indices[:target_true_count]

    result_list = [True if i in new_true_indices else False for i, val in enumerate(input_list)]
    return result_list

def swelling_axons_vector(path, perc_swollen_axons):
    
    axons = get_all_axons(path)
    
    # set all swellings to true 
    to_swell = [True]*len(axons)
    # 70, 50, 30
    perc_swollen_axons = sorted(perc_swollen_axons, reverse=True)
    
    swelling_axons = []

    for perc in perc_swollen_axons:

        print("Percentage ", perc)
        # list of bool, if true -> the axon does not shrink
        to_swell = adjust_percentage_of_true(to_swell, perc)
        swelling_axons.append(to_swell)

    return swelling_axons

def shrink_axons(path, perc_swelling, not_shrinking_axons, perc_swollen_axons, output_paths):

    header = get_first_line_from_file(path)

    axons = get_all_axons(path)

    perc_swollen_axons = sorted(perc_swollen_axons, reverse=True)
        
    # for every swelling percentage : 100, 70, etc. 
    for e, to_not_shrink in enumerate(not_shrinking_axons):

        new_lines = []

        new_lines.append(header)

        for i,axon in enumerate(axons):

            old_radius = float(axon[0][6])

         
            if (to_not_shrink[i]):
                for sphere in axon:
                    new_lines.append(sphere)
                    #print(sphere)
            else:
                # to shrink
                for sphere in axon:
                    sphere_ = sphere
                    old_radius = float(sphere_[6])
                    new_radius = shrink_radius(perc_swelling, old_radius)
                    sphere_[6] = str(new_radius)
                    #print(sphere)
                    new_lines.append(sphere_)
        axons = get_all_axons(path)
        
        write_list_of_lists_to_file(output_paths[e], new_lines)


def swell_axons(path, perc_swelling, not_shrinking_axons, perc_swollen_axons, output_paths):

    header = get_first_line_from_file(path)

    axons = get_all_axons(path)

    perc_swollen_axons = sorted(perc_swollen_axons, reverse=True)
        
    # for every swelling percentage : 100, 70, etc. 
    for e, to_swell in enumerate(not_shrinking_axons):

        new_lines = []

        new_lines.append(header)

        for i,axon in enumerate(axons):

            old_radius = float(axon[0][6])

         
            if (to_swell[i]):
                for sphere in axon:
                    sphere_ = sphere
                    old_radius = float(sphere_[6])
                    new_radius = swell_radius(perc_swelling, old_radius)
                    sphere_[6] = str(new_radius)
                    new_lines.append(sphere)
                    #print(sphere)
            else:
                # to shrink
                for sphere in axon:
                    new_lines.append(sphere)

        axons = get_all_axons(path)
        
        write_list_of_lists_to_file(output_paths[e], new_lines)

def get_total_area(df):
    axons = get_spheres_array(df)
    area = 0
    for axon in axons : 
        area = area + axon[0][3]*axon[0][3] # add area of first sphere
    return area 

def vary_perc_shrinking_of_each_axon(path):

    # shrink all by 1% 
    perc_swollen_axons = [0]
    perc_swelling = 0.01
    swelling_axons = swelling_axons_vector(path, perc_swollen_axons) # 100% False
    output_path1 = path[:-4]+"_swell_0.swc"
    shrink_axons(path, perc_swelling, swelling_axons, perc_swollen_axons, [output_path1])

    # swell by 0.25% , 0.5%, 0.75%, 1%
    path2 = output_path1
    perc_swollen_axons = [100]
    swelling_axons = swelling_axons_vector(path, perc_swollen_axons) # 100% True
    perc_swellings = [0.0025, 0.0050, 0.0075]
    output_paths = [path[:-4]+"_swell_"+str(perc_swelling)+path[-4:] for perc_swelling in perc_swellings]
    
    for i, perc_swelling in enumerate(perc_swellings):
        swell_axons(path2, perc_swelling, swelling_axons, perc_swollen_axons, [output_paths[i]])
        df = read_swc_file(output_paths[i])
        area = get_total_area(df)
        print(f"Area for {perc_swelling}% : {area}")




def vary_perc_swelling_axons():
    path = "/home/localadmin/Documents/Melina_branch/Sim_Growth/growth_icvf_0.50_cap_1_vox_50_factor_8_shrink_0.swc"

    #perc_swollen_axons = [0, 30, 50, 70, 100]
    perc_swollen_axons = [0]
    perc_swelling = 0.0025
    not_shrinking_axons = swelling_axons_vector(path, perc_swollen_axons)

    output_paths = [path[:-4]+"_swell_"+str(swell)+path[-4:] for swell in perc_swollen_axons]
        
    shrink_axons(path, perc_swelling, not_shrinking_axons, perc_swollen_axons, output_paths)
    #path1 = "/home/localadmin/Documents/Melina_branch/Sim_Growth/growth_icvf_0.50_cap_1_vox_50_factor_4.swc"
    #shrink_axons(path1, perc_swelling, not_shrinking_axons, perc_swollen_axons)
    #path2 = "/home/localadmin/Documents/Melina_branch/Sim_Growth/growth_icvf_0.50_cap_1_vox_50_factor_8.swc"
    #shrink_axons(path2, perc_swelling, not_shrinking_axons, perc_swollen_axons)

    for i in perc_swollen_axons:
        file_path = path[:-4]+"_swell_"+str(i)+path[-4:]
        df = read_swc_file(file_path)
        area = get_total_area(df)
        print(f"Area for {i}% : {area}")


vary_perc_shrinking_of_each_axon("/home/localadmin/Documents/Melina_branch/Sim_Growth/growth_icvf_0.50_cap_1_vox_50_factor_8.swc")