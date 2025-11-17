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
    std::string data_directory = "";
    std::string filename = "";
    std::vector<double> vox_sizes = {50};
    int repetitions = 1;
    int spheres_overlap_factor = 2;
    double glial_pop1_icvf_soma = 0.0;
    double glial_pop1_icvf_processes = 0.0;
    double glial_pop2_icvf_soma = 0.0;
    double glial_pop2_icvf_processes = 0.0;
    int nbr_axons_populations = 1;
    int crossing_fibers_type = 0;

    double axons_wo_myelin_icvf = 0.2;
    double axons_w_myelin_icvf = 0.1;
    double blood_vessels_icvf = 0.0;
    bool tortuous = true;
    double alpha = 4.0;
    double beta = 0.25;
    int regrow_thr = 10;
    double min_rad = 0.15;
    double std_dev = 0.01;
    int ondulation_factor = 5;
    double beading_period = 10;           /*!< Frequency for the beading  */
    double beading_variation = 0.0;
    double beading_variation_std = 0.0;
    bool can_shrink = true;
    double cosPhiSquared= 1.0;
    int nbr_threads = 1;
    double std_glial_process_length;
    double mean_glial_process_length;
    static int str_dist(std::string s, std::string t);
    void readConfFile(std::string conf_file_path);

    double c1 = 0.35;
    double c2 = 0.006;
    double c3 = 0.024;

};

#endif // PARAMETERS_H