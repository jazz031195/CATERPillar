#include "core_logic.h"
#include "CaterpillarGrowth.h"
#include <iostream>
#include <chrono>
#include <Eigen/Dense>
#include <nlohmann/json.hpp>
#include <fstream>
#include <filesystem>


CoreLogic::SimResult CoreLogic::runSimulation(const Parameters& params,
                                               ProgressCallback on_growth_progress,
                                               ProgressCallback on_swelling_progress) {
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
        CaterpillarGrowth Sim(params, min_l, max_l);
        Sim.on_growth_progress = on_growth_progress;
        Sim.on_swelling_progress = on_swelling_progress;

        // 2. Grow the Substrate!
        Sim.createSubstrate();

        // ==========================================


        std::string simulation_file_name;
        std::string swc_file_name;

        std::ifstream file(params.filename);

        if (rep ==0){
            simulation_file_name = (params.data_directory + "/" + params.filename + "_growth_info.txt" );
            swc_file_name = (params.data_directory + "/" + params.filename + ".csv");
        }
        else{
            simulation_file_name = (params.data_directory + "/" + params.filename + "_growth_info_" + std::to_string(rep) + ".txt");
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
        Sim.create_SWC_file(swc_file);
        swc_file.close();

        // Check if files opened successfully

        if (!simulation_file)
        {
            std::cerr << "Error opening output file : " << simulation_file_name <<std::endl;
        }

        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime);
        Sim.simulation_file(simulation_file, duration);
        simulation_file.close();

        if (rep == params.repetitions -1){
            std::cout << "\n========================================" << std::endl;
            std::cout << " All simulations completed successfully! " << std::endl;
            std::cout << "========================================" << std::endl;

            return std::make_tuple(Sim.axons,
                                Sim.blood_vessels,
                                Sim.glial_pop1,
                                Sim.glial_pop2,
                                Sim.glial_pop3,
                                Sim.min_limits,
                                Sim.max_limits);
                } 

    }


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
    // --- DIRECTORY CHECK & CREATION ---
    std::filesystem::path dirPath(params.data_directory);

    if (!std::filesystem::exists(dirPath)) {
        std::cout << "Warning: Output directory does not exist. Creating it now: "
                  << params.data_directory << std::endl;

        // create_directories creates the folder and any missing parent folders
        std::error_code ec;
        std::filesystem::create_directories(dirPath, ec);
        if (ec) {
            std::cerr << "CRITICAL ERROR: Failed to create directory. Check permissions." << std::endl;
            return; // Abort the simulation
        }
    }


    params.filename = data["GeneralParameters"].value("Filename", "simulation_output");
    
    // Updated to match your new Parameters class variable
    params.voxel_size = double(data["GeneralParameters"]["VoxelEdgeLength"]);
    
    params.repetitions = data["GeneralParameters"]["Repetitions"];
    params.spheres_overlap_factor = data["GeneralParameters"]["OverlappingFactor"];
    params.nbr_threads = data["GeneralParameters"]["NumberOfThreads"];

    // ==========================================
    // 2. Axon Parameters
    // ==========================================
    params.axons_wo_myelin_icvf = double(data["AxonParameters"]["AxonsICVF"]) / 100.0;
    params.axons_w_myelin_icvf = double(data["AxonParameters"]["AxonsWithMyelinICVF"]) / 100.0;
    // Arteriole = the main vessel (formerly "blood vessels"/"trunk"); Capillaries =
    // the branches that emerge from it. Renamed to match the new two-tier model.
    params.blood_vessels_icvf = double(data["AxonParameters"]["ArterioleICVF"]) / 100.0;
    params.blood_vessels_processes_icvf = data["AxonParameters"].value("CapillariesICVF", 0.0) / 100.0;
    params.capillary_radius = data["AxonParameters"].value("CapillaryRadius", 1.0);
    params.max_generations = data["AxonParameters"].value("MaxGenerations", 7);
    params.blood_vessels_voxel_size = data["AxonParameters"].value("BloodVesselsVoxelEdgeLength", 0.0);

    params.nbr_axons_populations = data["AxonParameters"]["NumberOfPopulations"];
    params.crossing_fibers_type = data["AxonParameters"]["CrossingFibersType"];
    
    // Axon morphology
    params.alpha = data["AxonParameters"]["Alpha"];
    params.beta = data["AxonParameters"]["Beta"];
    params.min_rad = data["AxonParameters"]["MinRadius"];
    params.epsilon = data["AxonParameters"]["Tortuosity_Epsilon"];
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
    params.axon_can_shrink = data["AxonParameters"].value("CanShrink", true);
    params.regrow_thr = data["AxonParameters"].value("RegrowThreshold", 10);
    params.undulation_factor = data["AxonParameters"].value("UndulationFactor", 5);
    params.swelling_factor = data["AxonParameters"].value("SwellingFactor", 1.0);

    // ==========================================
    // 3. Glial Parameters
    // ==========================================
    
    // Population 1
    params.glial_pop1_soma_icvf = double(data["GlialParameters"]["Pop1SomaICVF"]) / 100.0;
    params.glial_pop1_processes_icvf = double(data["GlialParameters"]["Pop1ProcessesICVF"]) / 100.0;
    params.glial_pop1_radius_mean = data["GlialParameters"]["Pop1SomaRadiusMean"];
    params.glial_pop1_radius_std = data["GlialParameters"]["Pop1SomaRadiusStd"];
    params.mean_glial_pop1_process_length = data["GlialParameters"]["Pop1MeanProcessLength"];
    params.std_glial_pop1_process_length = data["GlialParameters"]["Pop1StdProcessLength"];
    params.glial_pop1_nbr_primary_processes = data["GlialParameters"]["Pop1NbrPrimaryProcesses"];
    params.glial_pop1_branching = data["GlialParameters"]["Pop1Branching"];
    params.glial_pop1_minimum_process_radius = data["GlialParameters"].value("Pop1MinimumProcessRadius", 0.15);

    // Population 2
    params.glial_pop2_soma_icvf = double(data["GlialParameters"]["Pop2SomaICVF"]) / 100.0;
    params.glial_pop2_processes_icvf = double(data["GlialParameters"]["Pop2ProcessesICVF"]) / 100.0;
    params.glial_pop2_radius_mean = data["GlialParameters"]["Pop2SomaRadiusMean"];
    params.glial_pop2_radius_std = data["GlialParameters"]["Pop2SomaRadiusStd"];
    params.mean_glial_pop2_process_length = data["GlialParameters"]["Pop2MeanProcessLength"];
    params.std_glial_pop2_process_length = data["GlialParameters"]["Pop2StdProcessLength"];
    params.glial_pop2_nbr_primary_processes = data["GlialParameters"]["Pop2NbrPrimaryProcesses"];
    params.glial_pop2_branching = data["GlialParameters"]["Pop2Branching"];
    params.glial_pop2_minimum_process_radius = data["GlialParameters"].value("Pop2MinimumProcessRadius", 0.15);

    // Population 3 (optional -- defaults keep older config files without it working unchanged)
    params.glial_pop3_soma_icvf = data["GlialParameters"].value("Pop3SomaICVF", 0.0) / 100.0;
    params.glial_pop3_processes_icvf = data["GlialParameters"].value("Pop3ProcessesICVF", 0.0) / 100.0;
    params.glial_pop3_radius_mean = data["GlialParameters"].value("Pop3SomaRadiusMean", 3.0);
    params.glial_pop3_radius_std = data["GlialParameters"].value("Pop3SomaRadiusStd", 0.1);
    params.mean_glial_pop3_process_length = data["GlialParameters"].value("Pop3MeanProcessLength", 30.0);
    params.std_glial_pop3_process_length = data["GlialParameters"].value("Pop3StdProcessLength", 15.0);
    params.glial_pop3_nbr_primary_processes = data["GlialParameters"].value("Pop3NbrPrimaryProcesses", 5);
    params.glial_pop3_branching = data["GlialParameters"].value("Pop3Branching", true);
    params.glial_pop3_minimum_process_radius = data["GlialParameters"].value("Pop3MinimumProcessRadius", 0.15);

    // ==========================================
    // Run the math!
    // ==========================================
    runSimulation(params);


}

void CoreLogic::writeParametersToJson(const Parameters& params, const std::string& jsonFilePath) {
    nlohmann::json data;

    data["GeneralParameters"]["OutputDirectory"] = params.data_directory;
    data["GeneralParameters"]["Filename"] = params.filename;
    data["GeneralParameters"]["VoxelEdgeLength"] = params.voxel_size;
    data["GeneralParameters"]["Repetitions"] = params.repetitions;
    data["GeneralParameters"]["OverlappingFactor"] = params.spheres_overlap_factor;
    data["GeneralParameters"]["NumberOfThreads"] = params.nbr_threads;

    data["AxonParameters"]["AxonsICVF"] = params.axons_wo_myelin_icvf * 100.0;
    data["AxonParameters"]["AxonsWithMyelinICVF"] = params.axons_w_myelin_icvf * 100.0;
    data["AxonParameters"]["ArterioleICVF"] = params.blood_vessels_icvf * 100.0;
    data["AxonParameters"]["CapillariesICVF"] = params.blood_vessels_processes_icvf * 100.0;
    data["AxonParameters"]["CapillaryRadius"] = params.capillary_radius;
    data["AxonParameters"]["MaxGenerations"] = params.max_generations;
    data["AxonParameters"]["BloodVesselsVoxelEdgeLength"] = params.blood_vessels_voxel_size;
    data["AxonParameters"]["NumberOfPopulations"] = params.nbr_axons_populations;
    data["AxonParameters"]["CrossingFibersType"] = params.crossing_fibers_type;
    data["AxonParameters"]["Alpha"] = params.alpha;
    data["AxonParameters"]["Beta"] = params.beta;
    data["AxonParameters"]["MinRadius"] = params.min_rad;
    data["AxonParameters"]["Tortuosity_Epsilon"] = params.epsilon;
    data["AxonParameters"]["FODF_c2"] = params.cosPhiSquared;
    data["AxonParameters"]["BeadingAmplitude"] = params.beading_amplitude;
    data["AxonParameters"]["BeadingStd"] = params.beading_std;
    data["AxonParameters"]["K1"] = params.c1;
    data["AxonParameters"]["K2"] = params.c2;
    data["AxonParameters"]["K3"] = params.c3;
    data["AxonParameters"]["Tortuous"] = params.tortuous;
    data["AxonParameters"]["CanShrink"] = params.axon_can_shrink;
    data["AxonParameters"]["RegrowThreshold"] = params.regrow_thr;
    data["AxonParameters"]["UndulationFactor"] = params.undulation_factor;
    data["AxonParameters"]["SwellingFactor"] = params.swelling_factor;

    data["GlialParameters"]["Pop1SomaICVF"] = params.glial_pop1_soma_icvf * 100.0;
    data["GlialParameters"]["Pop1ProcessesICVF"] = params.glial_pop1_processes_icvf * 100.0;
    data["GlialParameters"]["Pop1SomaRadiusMean"] = params.glial_pop1_radius_mean;
    data["GlialParameters"]["Pop1SomaRadiusStd"] = params.glial_pop1_radius_std;
    data["GlialParameters"]["Pop1MeanProcessLength"] = params.mean_glial_pop1_process_length;
    data["GlialParameters"]["Pop1StdProcessLength"] = params.std_glial_pop1_process_length;
    data["GlialParameters"]["Pop1NbrPrimaryProcesses"] = params.glial_pop1_nbr_primary_processes;
    data["GlialParameters"]["Pop1Branching"] = params.glial_pop1_branching;
    data["GlialParameters"]["Pop1MinimumProcessRadius"] = params.glial_pop1_minimum_process_radius;

    data["GlialParameters"]["Pop2SomaICVF"] = params.glial_pop2_soma_icvf * 100.0;
    data["GlialParameters"]["Pop2ProcessesICVF"] = params.glial_pop2_processes_icvf * 100.0;
    data["GlialParameters"]["Pop2SomaRadiusMean"] = params.glial_pop2_radius_mean;
    data["GlialParameters"]["Pop2SomaRadiusStd"] = params.glial_pop2_radius_std;
    data["GlialParameters"]["Pop2MeanProcessLength"] = params.mean_glial_pop2_process_length;
    data["GlialParameters"]["Pop2StdProcessLength"] = params.std_glial_pop2_process_length;
    data["GlialParameters"]["Pop2NbrPrimaryProcesses"] = params.glial_pop2_nbr_primary_processes;
    data["GlialParameters"]["Pop2Branching"] = params.glial_pop2_branching;
    data["GlialParameters"]["Pop2MinimumProcessRadius"] = params.glial_pop2_minimum_process_radius;

    data["GlialParameters"]["Pop3SomaICVF"] = params.glial_pop3_soma_icvf * 100.0;
    data["GlialParameters"]["Pop3ProcessesICVF"] = params.glial_pop3_processes_icvf * 100.0;
    data["GlialParameters"]["Pop3SomaRadiusMean"] = params.glial_pop3_radius_mean;
    data["GlialParameters"]["Pop3SomaRadiusStd"] = params.glial_pop3_radius_std;
    data["GlialParameters"]["Pop3MeanProcessLength"] = params.mean_glial_pop3_process_length;
    data["GlialParameters"]["Pop3StdProcessLength"] = params.std_glial_pop3_process_length;
    data["GlialParameters"]["Pop3NbrPrimaryProcesses"] = params.glial_pop3_nbr_primary_processes;
    data["GlialParameters"]["Pop3Branching"] = params.glial_pop3_branching;
    data["GlialParameters"]["Pop3MinimumProcessRadius"] = params.glial_pop3_minimum_process_radius;

    std::ofstream out(jsonFilePath);
    if (!out) {
        std::cerr << "Error opening output file : " << jsonFilePath << std::endl;
        return;
    }
    out << data.dump(4);
}