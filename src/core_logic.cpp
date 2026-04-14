#include "core_logic.h"
#include "axongammadistribution.h"
#include <iostream>
#include <chrono>
#include <Eigen/Dense>
#include <nlohmann/json.hpp>
#include <fstream>
#include <QDir>
#include <QString>


void CoreLogic::runSimulation(const Parameters& params) {
    std::cout << "\n========================================" << std::endl;
    std::cout << " Starting CATERPillar Simulation Engine " << std::endl;
    std::cout << " Output Directory: " << params.data_directory << std::endl;
    std::cout << "========================================" << std::endl;

    // Loop through the number of requested repetitions
    for (int rep = 0; rep < params.repetitions; rep++) {
        std::cout << "\n[Repetition " << (rep + 1) << " of " << params.repetitions << "]" << std::endl;

        std::cout << "  -> Generating voxel size: " << params.voxel_size << " um" << std::endl;

        // Define min and max limits of the voxel
        Eigen::Vector3d min_l = {0.0, 0.0, 0.0};
        Eigen::Vector3d max_l = {params.voxel_size, params.voxel_size, params.voxel_size};

        // Start the execution timer
        auto startTime = std::chrono::high_resolution_clock::now();

        // ==========================================
        // THE CORE MATH
        // ==========================================
        
        // 1. Initialize the distribution engine using our clean parameters object
        AxonGammaDistribution AxonDistribution(params, min_l, max_l);

        // 2. Grow the Substrate!
        AxonDistribution.createSubstrate();

        // ==========================================


        std::string simulation_file_name;
        std::string swc_file_name;

        std::ifstream file(params.filename);

        if (rep ==0){
            simulation_file_name = (params.data_directory + "/" + params.filename + "_growth_info.txt" );
            swc_file_name = (params.data_directory + "/" + params.filename + ".csv");
        }
        else{
            simulation_file_name = (params.data_directory + "/" + params.filename + "/growth_info_" + std::to_string(rep) + ".txt");
            swc_file_name = (params.data_directory + "/" + params.filename + "_" + std::to_string(rep) + ".csv");
        }
        std::ofstream swc_file(swc_file_name);
        std::ofstream simulation_file(simulation_file_name);

        // Check if files opened successfully
        if (!swc_file)
        {
            std::cerr << "Error opening output file : "<< swc_file_name << std::endl;

        }
        cout << "Creating file: " << swc_file_name << endl;
        // write to file
        AxonDistribution.create_SWC_file(swc_file);
        swc_file.close();

        // Check if files opened successfully

        if (!simulation_file)
        {
            std::cerr << "Error opening output file : " << simulation_file_name <<std::endl;
        }

        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime);
        AxonDistribution.simulation_file(simulation_file, duration);
        simulation_file.close();

    }


    std::cout << "\n========================================" << std::endl;
    std::cout << " All simulations completed successfully! " << std::endl;
    std::cout << "========================================" << std::endl;
}

void CoreLogic::runSimulationFromJson(const std::string& jsonFilePath) {
    std::ifstream file(jsonFilePath);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open config file " << jsonFilePath << std::endl;
        return;
    }

    nlohmann::json data = nlohmann::json::parse(file);
    Parameters params;

    // ==========================================
    // 1. General Parameters
    // ==========================================
    params.data_directory = data["GeneralParameters"]["OutputDirectory"];
    // --- DIRECTORY CHECK & CREATION (C++14 / Qt Way) ---
    QString dirPath = QString::fromStdString(params.data_directory);
    QDir dir(dirPath);
    
    if (!dir.exists()) {
        std::cout << "Warning: Output directory does not exist. Creating it now: " 
                  << params.data_directory << std::endl;
        
        // mkpath creates the folder and any missing parent folders
        if (!dir.mkpath(".")) {
            std::cerr << "CRITICAL ERROR: Failed to create directory. Check permissions." << std::endl;
            return; // Abort the simulation
        }
    }


    params.filename = data["GeneralParameters"].value("Filename", "simulation_output");
    
    // Updated to match your new Parameters class variable
    params.voxel_size = double(data["GeneralParameters"]["VoxelEdgeLength"]);
    
    params.repetitions = data["GeneralParameters"]["Repetitions"];
    params.spheres_overlap_factor = data["GeneralParameters"]["OverlappingFactor"];

    // ==========================================
    // 2. Axon Parameters 
    // ==========================================
    params.axons_wo_myelin_icvf = double(data["AxonParameters"]["AxonsICVF"]) / 100.0;
    params.axons_w_myelin_icvf = double(data["AxonParameters"]["AxonsWithMyelinICVF"]) / 100.0;
    params.blood_vessels_icvf = double(data["AxonParameters"]["BloodVesselsICVF"]) / 100.0;
    
    params.nbr_axons_populations = data["AxonParameters"]["NumberOfPopulations"];
    params.crossing_fibers_type = data["AxonParameters"]["CrossingFibersType"];
    params.nbr_threads = data["AxonParameters"]["NumberOfThreads"];
    
    // Axon morphology
    params.alpha = data["AxonParameters"]["Alpha"];
    params.beta = data["AxonParameters"]["Beta"];
    params.min_rad = data["AxonParameters"]["MinRadius"];
    params.std_dev = data["AxonParameters"]["Tortuosity_Epsilon"];
    params.cosPhiSquared = data["AxonParameters"]["FODF_c2"];
    
    // Updated beading variables
    params.beading_amplitude = data["AxonParameters"]["BeadingAmplitude"];
    params.beading_std = data["AxonParameters"]["BeadingStd"];
    
    // Myelin constants
    params.c1 = data["AxonParameters"]["K1"];
    params.c2 = data["AxonParameters"]["K2"];
    params.c3 = data["AxonParameters"]["K3"];

    // Simulation engine flags
    params.tortuous = data["AxonParameters"].value("Tortuous", true);
    params.can_shrink = data["AxonParameters"].value("CanShrink", true);
    params.regrow_thr = data["AxonParameters"].value("RegrowThreshold", 10);
    params.ondulation_factor = data["AxonParameters"].value("OndulationFactor", 5);

    // ==========================================
    // 3. Glial Parameters
    // ==========================================
    
    // Population 1
    params.glial_pop1_icvf_soma = double(data["GlialParameters"]["Pop1SomaICVF"]) / 100.0;
    params.glial_pop1_icvf_processes = double(data["GlialParameters"]["Pop1ProcessesICVF"]) / 100.0;
    params.glial_pop1_radius_mean = data["GlialParameters"]["Pop1SomaRadiusMean"];
    params.glial_pop1_radius_std = data["GlialParameters"]["Pop1SomaRadiusStd"];
    params.mean_glial_pop1_process_length = data["GlialParameters"]["Pop1MeanProcessLength"];
    params.std_glial_pop1_process_length = data["GlialParameters"]["Pop1StdProcessLength"];
    params.glial_pop1_nbr_primary_processes = data["GlialParameters"]["Pop1NbrPrimaryProcesses"];
    params.glial_pop1_branching = data["GlialParameters"]["Pop1Branching"];

    // Population 2
    params.glial_pop2_icvf_soma = double(data["GlialParameters"]["Pop2SomaICVF"]) / 100.0;
    params.glial_pop2_icvf_processes = double(data["GlialParameters"]["Pop2ProcessesICVF"]) / 100.0;
    params.glial_pop2_radius_mean = data["GlialParameters"]["Pop2SomaRadiusMean"];
    params.glial_pop2_radius_std = data["GlialParameters"]["Pop2SomaRadiusStd"];
    params.mean_glial_pop2_process_length = data["GlialParameters"]["Pop2MeanProcessLength"];
    params.std_glial_pop2_process_length = data["GlialParameters"]["Pop2StdProcessLength"];
    params.glial_pop2_nbr_primary_processes = data["GlialParameters"]["Pop2NbrPrimaryProcesses"];
    params.glial_pop2_branching = data["GlialParameters"]["Pop2Branching"];

    // ==========================================
    // Run the math!
    // ==========================================
    runSimulation(params);


}