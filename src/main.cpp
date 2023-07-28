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
    // number of axons
    unsigned int number_axons = 100;
    int axon_capacity = 1;

    // constants for gamma distribution, mean = 0.5 um
    double alpha = 5.0;
    double beta = 0.1;

    double vox_size= 20;

    // min and max limits of voxel
    Eigen::Vector3d min_l = {0, 0, 0};
    Eigen::Vector3d max_l = {vox_size, vox_size, vox_size}; // um

    // minimum radius
    double min_radius = 0.15; // um

    // simulation parameters
    bool tortuous = true;
    bool draw = false;

    // density parameters
    double icvf = 0.3;

    // number of regrowth batches allowed
    int regrow_thr = 1000;

    auto startTime = std::chrono::high_resolution_clock::now();

    // create distribution of axons
    AxonGammaDistribution *AxonDistribution = new AxonGammaDistribution(number_axons, axon_capacity, alpha, beta, min_l, max_l, min_radius, tortuous, draw, regrow_thr);
    AxonDistribution->set_icvf(icvf, max_l[0], max_l[1]);

    cout << "AxonDistribution created" << endl;


    AxonDistribution->createSubstrate();


    // Open the output file stream to the desired file path
    //std::ofstream axons_file("/Users/melina/Desktop/EPFL/BachelorProject/Sim_Growth/files/axons_icvf_" + std::to_string(icvf).substr(0, 4) + " _cap_" + std::to_string(axon_capacity) + "_vox_" + std::to_string(max_l[0]).substr(0, 3) + ".swc");
    //std::ofstream simulation_file("/Users/melina/Desktop/EPFL/BachelorProject/Sim_Growth/files/simulation_icvf_" + std::to_string(icvf).substr(0, 4) + " _cap_" + std::to_string(axon_capacity) + "_vox_" + std::to_string(max_l[0]).substr(0, 3) + ".txt");
    //std::ofstream swc_file("/Users/melina/Desktop/EPFL/BachelorProject/Sim_Growth/files/growth_icvf_" + std::to_string(icvf).substr(0, 4) + " _cap_" + std::to_string(axon_capacity) + "_vox_" + std::to_string(max_l[0]).substr(0, 3) + ".swc");

    std::ofstream axons_file("/home/localadmin/Documents/Melina_branch/Sim_Growth/axons_icvf_" + std::to_string(icvf).substr(0, 4) + "_cap_" + std::to_string(axon_capacity) + "_vox_" + std::to_string(vox_size).substr(0, 2) + ".swc");
    std::ofstream simulation_file("/home/localadmin/Documents/Melina_branch/Sim_Growth/simulation_icvf_" + std::to_string(icvf).substr(0, 4) + "_cap_" + std::to_string(axon_capacity) + "_vox_" + std::to_string(vox_size).substr(0, 2) + ".txt");
    std::ofstream swc_file("/home/localadmin/Documents/Melina_branch/Sim_Growth/growth_icvf_" + std::to_string(icvf).substr(0, 4) + "_cap_" + std::to_string(axon_capacity) + "_vox_" + std::to_string(vox_size).substr(0, 2) + ".swc");

    // Check if files opened successfully
    if (!axons_file || !simulation_file || !swc_file)
    {
        std::cerr << "Error opening output file!" << std::endl;
        return 1;
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime);


    AxonDistribution->axons_file(axons_file);
    AxonDistribution->simulation_file(simulation_file, duration);
    AxonDistribution->create_SWC_file(swc_file);

    axons_file.close();
    simulation_file.close();
    swc_file.close();

    cout << "End of simulation!" << endl;
    delete AxonDistribution;
}
