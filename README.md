# CATERPillar  :Computational Axonal Threading Engine for Realistic Proliferation

![GitHub Logo](https://github.com/jazz031195/CATERPillar/blob/updated/caterpillar.jpg)

## Table of Contents
- [Installation](#section-1)
- [Run the program 2](#section-2)
  - [Linux](#subsection-21)
  - [MacOs](#subsection-22)
- [Parameters](#section-3)
- [Outputs of program](#section-4)
  
## Installations
### MacOs
Type the following in a terminal :
- /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
- brew install freeglut


## Run the program
### Linux 
- In the compiler_linux.sh, modify the following line : cd "/home/localadmin/Documents/Sim_Growth/src/", to the path that leads to the src folder on your local machine.
- Adapt the args.conf file to your preference. 
- In the terminal, use "cd" to go to the folder.
- Type : "./run_linux.sh"

### MacOs
- In the compiler_mac.sh, modify the following line : cd "/home/localadmin/Documents/Sim_Growth/src/", to the path that leads to the src folder on your local machine.
- Adapt the args.conf file to your preference. 
- In the terminal, use "cd" to go to the folder.
- Type : "./run_mac.sh"

## Parameters
The parameters are organised in the args.conf file. Some parameters may have multiple values (those written as <parameter>) over which the program will iterate. The parameter value will then appear in the name of the file.
The parameters that can be modified are the following :

- data_directory (string): path to the folder on your local computer
- alpha (double) : alpha parameter for a Gamma Distribution. The radii of the axons arew taken from this distribution.
- beta (double) : beta parameter for a Gamma Distribution. The radii of the axons arew taken from this distribution.
- vox_sizes (list of int) : The maximum distance (in micrometers) in space to define the corner of the 3D box to work in.
- min_radius (double) : minimum radius of the axons.
- tortuous (0 or 1) : if 1, the axons will be tortuouss, overwise they will grow straight.
- beading_variation (double between 0 and 1) : The radius of each axon will vary periodically between the original radius and radius/beading_variation. If set to 1, there will be no beading.
- draw (0 or 1) : If set to 1, the voxel will be drawn with all the growing axons. If so, no output with the axons' information will be created.
- regrow_thr (int) : If the amount of regrowing attempts during which no new axons are made is above the chosen threshold, the growth process stops.
- spheres_overlap_factors (list of int, must be above 2) : The distance between spheres in each axon is the radius/spheres_overlap_factor.
- icvf (double between 0 and 1) : Intra Compartement Volume Fraction
- capacities (list of int, must be bellow the number of processors available on computer) : number of parallel threads to grow axons simultaniously
- repetitions (int, must be 1 or above) : Number of repetitions for the growing process

## Outputs of program
- name of file : directory + "/simulation_icvf_" + icvf + "_cap_" + capacity + "_vox_" + vox_size + "_"+ repetition + ".txt"
  -  Duration (in seconds)
  -  Num_axons (number of axons)
  -  Capacity (number of processors)
  -  Voxel (voxel size in micrometers)
  -  icvf (Intra Compartement Volume Fraction reached)
- name of file : directory + "/growth_icvf_" + icvf + "_cap_" + capacity + "_vox_" + vox_size + "_factor_"+ factor + "_" + repetition + ".swc"
  - Position of each sphere of each axon placed
  - id_ax: ID of axon
  - id_sph: ID of sphere
  - Type : type of cell (axon)
  - X Y Z : coordinates of center of sphere in 3D space in micrometers
  - R : radius in micrometers
  - P : ID of parent sphere 
## Developers 
- Jasmine Nguyen-Duc
- Melina Cherchali
- Ines deRiedmatten
