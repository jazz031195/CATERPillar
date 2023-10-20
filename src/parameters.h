//!  Basic class to store simulation parameters =============================================================/
/*!
*   \details   Basic class to store and handle all the possible simulation parameters.
*   \author    Jasmine Nguyen-DUc
*   \date      October 2023
*   \version   0
*===========================================================================================================*/
#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <vector>
#include "Eigen/Core"
#include <utility>

using namespace std;

/*! \class Parameter
 *  \brief Class used to hold and operate all the parameters.
 */

class Parameters
{
public:
    std::string data_directory;
    std::vector<double> capacities;
    std::vector<double> vox_sizes;
    int repetitions;
    std::vector<double> icvf;
    std::vector<double> spheres_overlap_factors;
    int beading_variation;

    static int str_dist(std::string s, std::string t);
    void readConfFile(std::string conf_file_path);
    void readVector(ifstream& in , std::vector<double> &parameter);
};

#endif // PARAMETERS_H