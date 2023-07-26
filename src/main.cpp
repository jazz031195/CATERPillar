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
    int axon_capacity = 5;

    // constants for gamma distribution, mean = 0.5 um
    double alpha = 5.0;
    double beta = 0.1;

    // min and max limits of voxel
    Eigen::Vector3d min_l = {0, 0, 0};
    Eigen::Vector3d max_l = {5, 5, 5}; // um

    // minimum radius
    double min_radius = 0.15; // um

    // simulation parameters
    bool tortuous = true;
    bool draw = false;

   // density parameters
    double icvf = 0.7;

    // create distribution of axons
    AxonGammaDistribution *AxonDistribution = new AxonGammaDistribution(number_axons, axon_capacity, alpha, beta, min_l, max_l, min_radius, tortuous, draw);
    AxonDistribution->set_icvf(icvf, max_l[0], max_l[1]);

    cout << "AxonDistribution created" << endl;

    auto startTime = std::chrono::high_resolution_clock::now();

    AxonDistribution->parallelGrowth();
    // AxonDistribution->axonDensityMap();

    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::minutes>(endTime - startTime);

    // Open the output file stream to the desired file path
    std::ofstream axons_file("/Users/melina/Desktop/EPFL/BachelorProject/Sim_Growth/axons.swc");
    std::ofstream simulation_file("/Users/melina/Desktop/EPFL/BachelorProject/Sim_Growth/simulation");

    // Check if files opened successfully
    if (!axons_file || !simulation_file)
    {
        std::cerr << "Error opening output file!" << std::endl;
        return 1;
    }

    AxonDistribution->axons_file(axons_file);
    AxonDistribution->simulation_file(simulation_file, duration);

    axons_file.close();
    simulation_file.close();

    cout << "End of simulation!" << endl;
    delete AxonDistribution;
}
