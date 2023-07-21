#include <iostream>
#include "Eigen/Core"
#include "Eigen/Dense"
#include <fstream>
#include <SFML/Graphics.hpp>
#include "axongammadistribution.h"

typedef unsigned int uint;

using namespace std;
using namespace Eigen;

int main()
{

    // number of axons
    unsigned int number_axons = 100;
    int axon_capacity = 10;

    // constants for gamma distribution, mean = 0.5 um
    double alpha = 5.0;
    double beta = 0.1;

    // min and max limits of voxel
    Eigen::Vector3d min_l = {0, 0, 0};
    Eigen::Vector3d max_l = {20, 20, 20}; // um

    // minimum radius
    double min_radius = 0.15; // um

    // simulation parameters
    bool tortuous = true;
    bool draw = false;

    // density parameters
<<<<<<< HEAD
    double icvf = 0.01;
=======
    double icvf = 0.5;
>>>>>>> 07190cf774c4d2ca0888b24bccea5cd001616ac5

    // create distribution of axons
    AxonGammaDistribution *AxonDistribution = new AxonGammaDistribution(number_axons, axon_capacity, alpha, beta, min_l, max_l, min_radius, tortuous, draw);
    AxonDistribution->set_icvf(icvf, max_l[0], max_l[1]);
    cout << "AxonDistribution created" << endl;
    AxonDistribution->parallelGrowth_();

    // Open the output file stream to the desired file path
    std::ofstream outFile("/Users/melina/Desktop/EPFL/BachelorProject/Sim_Growth/axon_simulation.swc");
    std::ofstream outFile2("/Users/melina/Desktop/EPFL/BachelorProject/Sim_Growth/radius.swc");


    // Check if the file was opened successfully
    if (!outFile)
    {
        std::cerr << "Error opening output file!" << std::endl;
        return 1;
    }

    AxonDistribution->create_SWC_file(outFile);
    outFile.close();

    cout << "End of simulation!" << endl;
    delete AxonDistribution;
}
