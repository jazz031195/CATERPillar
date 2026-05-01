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
    double voxel_size = 50;
    int repetitions = 1;
    int spheres_overlap_factor = 2;
    double glial_pop1_soma_icvf = 0.0;
    double glial_pop1_processes_icvf = 0.0;
    double glial_pop2_soma_icvf = 0.0;
    double glial_pop2_processes_icvf = 0.0;
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
    double epsilon = 0.01;
    int undulation_factor = 5;
    double beading_amplitude = 0.3;    
    double beading_std = 0.05;

    double epsilon_blood_vessels = 0.0;
    double mean_vessel_rad = 6;
    double std_vessel_rad = 1;

    double swelling_factor = 1.0;     

    bool axon_can_shrink = true;
    double cosPhiSquared= 1.0;
    int nbr_threads = 1;

    double c1 = 0.35;
    double c2 = 0.006;
    double c3 = 0.024;

    int glial_pop1_nbr_primary_processes = 5;
    int glial_pop2_nbr_primary_processes = 5;

    double std_glial_pop1_process_length;
    double mean_glial_pop1_process_length;
    double std_glial_pop2_process_length;
    double mean_glial_pop2_process_length;

    bool glial_pop1_branching = true;
    bool glial_pop2_branching = true;
    double glial_pop1_radius_mean = 0.5;
    double glial_pop1_radius_std = 0.1;
    double glial_pop2_radius_mean = 0.5;
    double glial_pop2_radius_std = 0.1;

};

#endif // PARAMETERS_H