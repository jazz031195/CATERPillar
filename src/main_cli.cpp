#include <iostream>
#include <string>
#include "core_logic.h"

void printUsage() {
    std::cout << "Usage: CATERPillar-cli --config <config.json>\n";
    std::cout << "  Runs the white matter growth simulation headlessly, without the GUI.\n";
}

int main(int argc, char *argv[]) {
    if (argc < 3 || std::string(argv[1]) != "--config") {
        printUsage();
        return -1;
    }

    std::string jsonFilePath = argv[2];
    std::cout << "Starting in headless mode using: " << jsonFilePath << std::endl;

    CoreLogic::runSimulationFromJson(jsonFilePath);

    return 0;
}
