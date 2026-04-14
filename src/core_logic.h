#pragma once
#include "src/parameters.h"
#include <string>

class CoreLogic {
public:

    // This function reads the JSON, fills the struct, and calls runSimulation()
    static void runSimulationFromJson(const std::string& jsonFilePath);
    static void runSimulation(const Parameters& params);
};