#include <iostream>
#include "Eigen/Core"
#include "Eigen/Dense"
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include "axongammadistribution.h"
#include "parameters.h"
#include <chrono>
#include <variant>

typedef unsigned int uint;

using namespace std;
using namespace Eigen;

// Function to check if a file exists
bool fileExists(const std::string& filename) {
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}

int main(int argn, char* argv[])
{
    string conf = "";

    if (argn == 2) {
        conf = argv[1];
    } else {
        return -1;
    }

    // read parameters from args.conf file
    Parameters parameters;
    parameters.readConfFile(conf);

    int repetitions = parameters.repetitions;
    std::vector<double> vox_sizes = parameters.vox_sizes;

    double glial_pop1_icvf_soma = parameters.glial_pop1_icvf_soma;
    double glial_pop1_icvf_processes = parameters.glial_pop1_icvf_processes;
    double glial_pop2_icvf_soma = parameters.glial_pop2_icvf_soma;
    double glial_pop2_icvf_processes = parameters.glial_pop2_icvf_processes;
    double axons_wo_myelin_icvf = parameters.axons_wo_myelin_icvf;
    double axons_w_myelin_icvf = parameters.axons_w_myelin_icvf;
    double blood_vessels_icvf = parameters.blood_vessels_icvf;
    int spheres_overlap_factor = parameters.spheres_overlap_factor;
    double beading_variation = parameters.beading_variation;
    double beading_variation_std = parameters.beading_variation_std;
    std::string directory = parameters.data_directory;
    std::string filename = parameters.filename;
    bool tortuous = parameters.tortuous;
    double alpha = parameters.alpha;
    double beta = parameters.beta;
    int regrow_thr = parameters.regrow_thr;
    double min_radius = parameters.min_rad;
    double std_dev = parameters.std_dev;
    int ondulation_factor = parameters.ondulation_factor;
    double beading_period = parameters.beading_period;
    bool can_shrink = parameters.can_shrink;
    double cosPhiSquared = parameters.cosPhiSquared;
    int nbr_threads = parameters.nbr_threads;
    int nbr_axons_populations = parameters.nbr_axons_populations;
    int crossing_fibers_type = parameters.crossing_fibers_type;
    double mean_glial_process_length = parameters.mean_glial_process_length;
    double std_glial_process_length = parameters.std_glial_process_length;

    double glial_pop1_radius_mean = 5.0;
    double glial_pop1_radius_std = 0.5;
    double glial_pop2_radius_mean = 5.0;
    double glial_pop2_radius_std = 0.5;
    bool glial_pop1_branching = true;
    bool glial_pop2_branching = true;
    int nbr_primary_processes_pop1 = 5;
    int nbr_primary_processes_pop2 = 5;

    double c1 = 0.35;
    double c2 = 0.006;
    double c3 = 0.024;

    for (int rep = 0; rep < repetitions; rep++) {
        for (int i = 0; i < vox_sizes.size(); i++) {
            double vox_size = vox_sizes[i];

            // min and max limits of voxel
            Eigen::Vector3d min_l = {0, 0, 0};
            Eigen::Vector3d max_l = {vox_size, vox_size, vox_size}; // um

            auto startTime = std::chrono::high_resolution_clock::now();

            AxonGammaDistribution* AxonDistribution = new AxonGammaDistribution(
                axons_wo_myelin_icvf, axons_w_myelin_icvf, glial_pop1_icvf_soma, glial_pop1_icvf_processes, glial_pop2_icvf_soma, glial_pop2_icvf_processes, blood_vessels_icvf, alpha, beta,
                                             min_l, max_l, min_radius, regrow_thr, beading_variation, beading_variation_std, std_dev, ondulation_factor, spheres_overlap_factor, can_shrink, cosPhiSquared, nbr_threads, nbr_axons_populations, crossing_fibers_type, 
                                              mean_glial_process_length, std_glial_process_length, mean_glial_process_length, std_glial_process_length,
                                              glial_pop1_radius_mean, glial_pop1_radius_std, glial_pop2_radius_mean, glial_pop2_radius_std, glial_pop1_branching, glial_pop2_branching, nbr_primary_processes_pop1, nbr_primary_processes_pop2, c1, c2, c3);

            // Generate unique filenames
            int file_suffix = 0;
            std::string info_filename;
            std::string swc_filename;

            do {
                std::ostringstream suffix;
                if (file_suffix > 0) {
                    suffix << "_rep" << (file_suffix < 10 ? "0" : "") << file_suffix;
                }
                std::string suffix_str = suffix.str();
                info_filename = directory + "/" + filename + suffix_str + "_info.txt";
                swc_filename = directory + "/" + filename + suffix_str + ".swc";
                file_suffix++;
            } while (fileExists(swc_filename) || fileExists(info_filename));

            // Create and open files
            std::ofstream swc_file(swc_filename);
            std::ofstream info_file(info_filename);

            // Check if files opened successfully
            if (!swc_file) {
                std::cerr << "Error opening output file: " << swc_filename << std::endl;
                delete AxonDistribution;
                return 1;
            } else {
                cout << "Creating file: " << swc_filename << endl;
            }

            if (!info_file) {
                std::cerr << "Error opening output file: " << info_filename << std::endl;
                delete AxonDistribution;
                return 1;
            } else {
                cout << "Creating file: " << info_filename << endl;
            }

            AxonDistribution->createSubstrate();

            // Write to SWC file
            AxonDistribution->create_SWC_file(swc_file);
            swc_file.close();

            auto endTime = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime);

            // Write to info file
            AxonDistribution->simulation_file(info_file, duration);
            info_file.close();

            cout << "End of simulation!" << endl;
            delete AxonDistribution;
        }
    }

    return 0;
}
