#pragma once
#include "src/parameters.h"
#include "../src/Axon.h"
#include "../src/Blood_Vessel.h"
#include "../src/Glial.h"
#include <functional>
#include <string>
#include <tuple>
#include <vector>

class CoreLogic {
public:
    using SimResult = std::tuple<std::vector<Axon>,
                                    std::vector<Blood_Vessel>,
                                    std::vector<Glial>,
                                    std::vector<Glial>>;
    using ProgressCallback = std::function<void(double, double)>;

    // This function reads the JSON, fills the struct, and calls runSimulation()
    static void runSimulationFromJson(const std::string& jsonFilePath);

    // on_growth_progress: (depth travelled, total box depth), fired once per
    // depth layer. on_swelling_progress: (current axon ICVF, target axon
    // ICVF), fired once per swelling round. Both optional -- left unset
    // (default-constructed empty std::function), CLI callers are unaffected.
    static SimResult runSimulation(const Parameters& params,
                                    ProgressCallback on_growth_progress = nullptr,
                                    ProgressCallback on_swelling_progress = nullptr);
};