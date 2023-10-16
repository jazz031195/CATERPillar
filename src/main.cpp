#include <iostream>
#include "Eigen/Core"
#include "Eigen/Dense"
#include <fstream>
#include <SFML/Graphics.hpp>
#include "axongammadistribution.h"
#include <chrono>

typedef unsigned int uint;

using namespace std;
using namespace Eigen;

int main()
{
    for (int rep = 0; rep < 1; rep++){
        std::vector<int> vox_sizes = {20};
        for (int i = 0; i < vox_sizes.size(); i++){
            std::vector<int> capacities = {25};
            for (int c = 0; c < capacities.size(); c++){
                // number of axons
                unsigned int number_axons = 100;
                int axon_capacity = capacities[c];

                // constants for gamma distribution, mean = 0.5 um
                double alpha = 5.0;
                double beta = 0.1;
                double vox_size= vox_sizes[i];

                // min and max limits of voxel
                Eigen::Vector3d min_l = {0, 0, 0};
                Eigen::Vector3d max_l = {vox_size, vox_size, vox_size}; // um

                // minimum radius
                double min_radius = 0.15; // um

                bool can_shrink = false;

                // simulation parameters
                bool tortuous = true;
                bool draw = false;
                //double beading_variation = 0.55;
                double beading_variation = 1;

                // density parameters
                double icvf = 0.6;

                // number of regrowth batches allowed
                int regrow_thr = 5;

                // distance between spheres = radius /spheres_overlap_factor
                std::vector<int> spheres_overlap_factors = {2};

                auto startTime = std::chrono::high_resolution_clock::now();

                // create distribution of axons
                AxonGammaDistribution *AxonDistribution = new AxonGammaDistribution(number_axons, axon_capacity, alpha, beta, min_l, max_l, min_radius, tortuous, draw, regrow_thr, can_shrink, beading_variation);
                AxonDistribution->set_icvf(icvf, max_l[0], max_l[1]);

                cout << "AxonDistribution created" << endl;


                AxonDistribution->createSubstrate();


                // Open the output file stream to the desired file path
                //std::ofstream axons_file("/Users/melina/Desktop/EPFL/BachelorProject/Sim_Growth/files/axons_icvf_" + std::to_string(icvf).substr(0, 4) + " _cap_" + std::to_string(axon_capacity) + "_vox_" + std::to_string(max_l[0]).substr(0, 3) + ".swc");
                //std::ofstream simulation_file("/Users/melina/Desktop/EPFL/BachelorProject/Sim_Growth/files/simulation_icvf_" + std::to_string(icvf).substr(0, 4) + " _cap_" + std::to_string(axon_capacity) + "_vox_" + std::to_string(max_l[0]).substr(0, 3) + ".txt");
                //std::ofstream swc_file("/Users/melina/Desktop/EPFL/BachelorProject/Sim_Growth/files/growth_icvf_" + std::to_string(icvf).substr(0, 4) + " _cap_" + std::to_string(axon_capacity) + "_vox_" + std::to_string(max_l[0]).substr(0, 3) + ".swc");


                std::string icvf_str = std::to_string(icvf);
                std::string vox_size_str = std::to_string(vox_size);
                

                std::string axons_file_name;
                std::string simulation_file_name;
                std::string swc_file_name;
                

                if (!draw){
                    std::string directory = "/home/localadmin/Documents/Melina_branch/Sim_Growth";
                    std::string filename;
                    if (tortuous) {
                        filename = directory + "/axons_icvf_" + std::to_string(icvf).substr(0, 4) + "_cap_" + std::to_string(axon_capacity) + "_vox_" + std::to_string(vox_size).substr(0, 2) + ".swc"; // Replace with your file name
                    }
                    else{
                        filename = directory + "/axons_icvf_" + std::to_string(icvf).substr(0, 4) + "_cap_" + std::to_string(axon_capacity) + "_vox_" + std::to_string(vox_size).substr(0, 2) +  "_straight.swc"; // Replace with your file name
                    }
                    std::ifstream file(filename);
                    if (file.good()) {
                        for (int n = 0; n < 10; n++) {
                            std::string new_filename = filename.substr(0, filename.size() - 4) + "_" + std::to_string(n) + ".swc";
                            std::ifstream new_file(new_filename);
                            if (new_file.good()) {
                                continue;
                            } else {
                                if (tortuous){
                                    axons_file_name = (directory + "axons_icvf_" + icvf_str.substr(0, 4) + "_cap_" + std::to_string(axon_capacity) + "_vox_" + vox_size_str.substr(0, 2) +   "_" + std::to_string(n) +".swc");
                                    simulation_file_name = (directory + "/simulation_icvf_" + icvf_str.substr(0, 4) + "_cap_" + std::to_string(axon_capacity) + "_vox_" + vox_size_str.substr(0, 2) +   "_" + std::to_string(n) +".txt");
                                    swc_file_name = (directory + "/growth_icvf_" + icvf_str.substr(0, 4) + "_cap_" + std::to_string(axon_capacity) + "_vox_" + vox_size_str.substr(0, 2) +   "_" + std::to_string(n) +".swc");
                                }
                                else{
                                    axons_file_name = (directory + "axons_icvf_" + icvf_str.substr(0, 4) + "_cap_" + std::to_string(axon_capacity) + "_vox_" + vox_size_str.substr(0, 2) +  "_" + std::to_string(n) +"_straight.swc");
                                    simulation_file_name = (directory + "/simulation_icvf_" + icvf_str.substr(0, 4) + "_cap_" + std::to_string(axon_capacity) + "_vox_" + vox_size_str.substr(0, 2) +  "_" + std::to_string(n) + "_straight.txt");
                                    swc_file_name = (directory + "/growth_icvf_" + icvf_str.substr(0, 4) + "_cap_" + std::to_string(axon_capacity) + "_vox_" + vox_size_str.substr(0, 2) +   "_" + std::to_string(n) + "_straight.swc");
                                }
                                break;
                            }
                        }
                    } else {
                        if (tortuous){
                            axons_file_name = (directory + "/axons_icvf_" + icvf_str.substr(0, 4) + "_cap_" + std::to_string(axon_capacity) + "_vox_" + vox_size_str.substr(0, 2) +  ".swc");
                            simulation_file_name = (directory + "/simulation_icvf_" + icvf_str.substr(0, 4) + "_cap_" + std::to_string(axon_capacity) + "_vox_" + vox_size_str.substr(0, 2) +  ".txt");
                            swc_file_name = (directory + "/growth_icvf_" + icvf_str.substr(0, 4) + "_cap_" + std::to_string(axon_capacity) + "_vox_" + vox_size_str.substr(0, 2) +  ".swc");
                        }
                        else{
                            axons_file_name = (directory + "/axons_icvf_" + icvf_str.substr(0, 4) + "_cap_" + std::to_string(axon_capacity) + "_vox_" + vox_size_str.substr(0, 2) +  "_straight.swc");
                            simulation_file_name = (directory + "/simulation_icvf_" + icvf_str.substr(0, 4) + "_cap_" + std::to_string(axon_capacity) + "_vox_" + vox_size_str.substr(0, 2) +  "_straight.txt");
                            swc_file_name = (directory + "/growth_icvf_" + icvf_str.substr(0, 4) + "_cap_" + std::to_string(axon_capacity) + "_vox_" + vox_size_str.substr(0, 2) +  "_straight.swc");
                        
                        }
                    }

                    auto endTime = std::chrono::high_resolution_clock::now();
                    auto duration = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime);


                    for (int f = 0; f < spheres_overlap_factors.size(); f++){
                        // factor for overlapping spheres
                        std::string factor_str =std::to_string(spheres_overlap_factors[f]);
                        // add to name of file
                        std::string swc_file_name_ = swc_file_name.substr(0, swc_file_name.size()-4);
                        swc_file_name_ = swc_file_name_ + "_factor_" +factor_str+"_"+std::to_string(rep)+".swc";
                        // create file
                        std::ofstream swc_file(swc_file_name_);
                        // Check if files opened successfully
                        if (!swc_file)
                        {
                            std::cerr << "Error opening output file!" << std::endl;
                            return 1;
                        }
                        // write to file
                        AxonDistribution->create_SWC_file(swc_file, spheres_overlap_factors[f]);

                        swc_file.close();
                        
                    }

                    //std::ofstream axons_file(axons_file_name);
                    std::string simulation_file_name_ = simulation_file_name.substr(0, simulation_file_name.size()-4);
                    simulation_file_name_ = simulation_file_name_ + "_"+std::to_string(rep)+".swc";
                        
                    std::ofstream simulation_file(simulation_file_name_);

                    // Check if files opened successfully
                    if (!simulation_file)
                    {
                        std::cerr << "Error opening output file!" << std::endl;
                        return 1;
                    }

                    

                    //AxonDistribution->axons_file(axons_file);
                    AxonDistribution->simulation_file(simulation_file, duration);

                    //axons_file.close();
                    simulation_file.close();
                    
                }

                cout << "End of simulation!" << endl;
                delete AxonDistribution;
            }
        }
    } 
}
