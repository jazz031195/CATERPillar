#include <iostream>
#include "Eigen/Core"
#include "Eigen/Dense"
#include <fstream>
#include <SFML/Graphics.hpp>
#include "axongammadistribution.h"

typedef unsigned int uint;

using namespace std;
using namespace Eigen;

int main() {
    // number of axons
    unsigned int number_axons = 100;
    // constants for gamma distribution, mean = 0.5 um 
    double alpha = 5.0;
    double beta = 0.1;
    // min and max limits of voxel
    Eigen::Vector3d min_l = {0,0,0};
    Eigen::Vector3d max_l = {50, 50, 50}; //um
    // minimum radius
    double min_radius = 0.15; // um
    // whether the substrate has tortuous axons 
    bool tortuous = true;
    // create distribution of axons
    AxonGammaDistribution* AxonDistribution = new AxonGammaDistribution(number_axons, alpha, beta, min_l, max_l, min_radius, tortuous);
    cout << "AxonDistribution created" << endl;
    AxonDistribution->createGammaSubstrate();
    //delete AxonDistribution;
}