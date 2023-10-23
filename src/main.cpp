#include <iostream>
#include "Eigen/Core"
#include "Eigen/Dense"
#include <fstream>
#include <SFML/Graphics.hpp>
#include "axongammadistribution.h"
#include "parameters.h"
#include <chrono>
#include <variant>

typedef unsigned int uint;

using namespace std;
using namespace Eigen;

int main(int argn, char* argv[])
{

    string conf = "";

    if(argn == 2){
        conf = argv[1];
    }
    else{
        return -1;
    }

    // read parameters from args.conf file
    Parameters parameters;
    parameters.readConfFile(conf);

    int repetitions = parameters.repetitions;
    std::vector<double> capacities = parameters.capacities;
    std::vector<double> vox_sizes = parameters.vox_sizes;
    std::vector<double> icvfs = parameters.icvf;
    std::vector<double> spheres_overlap_factors = parameters.spheres_overlap_factors;
    int beading_variation = parameters.beading_variation;
    std::string directory = parameters.data_directory;
    bool tortuous = parameters.tortuous;
    bool draw = parameters.draw;
    double alpha = parameters.alpha;
    double beta = parameters.beta;
    int regrow_thr = parameters.regrow_thr;
    double min_radius = parameters.min_rad;

    cout << icvfs.size() << endl;
    for (int icvf_index = 0; icvf_index < icvfs.size() ; icvf_index++){
        
        for (int rep = 0; rep < repetitions; rep++){
            //std::vector<int> vox_sizes = {50};
            for (int i = 0; i < vox_sizes.size(); i++){
                //std::vector<int> capacities = {24, 18, 12, 6, 3, 1};
                for (int c = 0; c < capacities.size(); c++){

                    int axon_capacity = capacities[c];
                    double vox_size= vox_sizes[i];

                    // min and max limits of voxel
                    Eigen::Vector3d min_l = {0, 0, 0};
                    Eigen::Vector3d max_l = {vox_size, vox_size, vox_size}; // um

                    double icvf = icvfs[icvf_index];

                    // distance between spheres = radius /spheres_overlap_factor
                    //std::vector<int> spheres_overlap_factors = {2};

                    auto startTime = std::chrono::high_resolution_clock::now();

                    // create distribution of axons
                    AxonGammaDistribution *AxonDistribution = new AxonGammaDistribution(icvf, axon_capacity, alpha, beta, min_l, max_l, min_radius, tortuous, draw, regrow_thr, beading_variation);

                    cout << "AxonDistribution created" << endl;


                    AxonDistribution->createSubstrate();

                    std::string icvf_str = std::to_string(icvf);
                    std::string vox_size_str = std::to_string(vox_size);
                    

                    std::string axons_file_name;
                    std::string simulation_file_name;
                    std::string swc_file_name;
                    

                    if (!draw){
                        std::string filename;
                        std::ifstream file(filename);

                        if (tortuous){
                            simulation_file_name = (directory + "/simulation_icvf_" + icvf_str.substr(0, 4) + "_cap_" + std::to_string(axon_capacity) + "_vox_" + vox_size_str.substr(0, 2) +  ".txt");
                            swc_file_name = (directory + "/growth_icvf_" + icvf_str.substr(0, 4) + "_cap_" + std::to_string(axon_capacity) + "_vox_" + vox_size_str.substr(0, 2) +  ".swc");
                        }
                        else{
                            simulation_file_name = (directory + "/simulation_icvf_" + icvf_str.substr(0, 4) + "_cap_" + std::to_string(axon_capacity) + "_vox_" + vox_size_str.substr(0, 2) +  "_straight.txt");
                            swc_file_name = (directory + "/growth_icvf_" + icvf_str.substr(0, 4) + "_cap_" + std::to_string(axon_capacity) + "_vox_" + vox_size_str.substr(0, 2) +  "_straight.swc");
                            
                        }

                        auto endTime = std::chrono::high_resolution_clock::now();
                        auto duration = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime);


                        for (int f = 0; f < spheres_overlap_factors.size(); f++){
                            // factor for overlapping spheres
                            std::string factor_str =std::to_string(int(spheres_overlap_factors[f]));
                            // add to name of file
                            std::string swc_file_name_ = swc_file_name.substr(0, swc_file_name.size()-4);
                            swc_file_name_ = swc_file_name_ + "_factor_" +factor_str+"_"+std::to_string(rep)+".swc";
                            // create file
                            std::ofstream swc_file(swc_file_name_);
                            // Check if files opened successfully
                            if (!swc_file)
                            {
                                std::cerr << "Error opening output file : "<< swc_file_name_ << std::endl;
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
                            std::cerr << "Error opening output file : " << simulation_file_name_ <<std::endl;
                            return 1;
                        }

                        AxonDistribution->simulation_file(simulation_file, duration);
                        simulation_file.close();
                        
                    }

                    cout << "End of simulation!" << endl;
                    delete AxonDistribution;
                }
            }
        }
    } 
    return 0;
}
