#pragma once
#include "src/parameters.h"
#include "../src/Axon.h"
#include "../src/Blood_Vessel.h"
#include "../src/Glial.h"
#include "Eigen/Core"
#include <functional>
#include <string>
#include <tuple>
#include <vector>

class CoreLogic {
public:
    // Trailing Eigen::Vector3d pair is the real (small) voxel's actual final
    // placement -- min_limits/max_limits, see CaterpillarGrowth::PlaceSmallVoxel
    // -- which can land anywhere inside the (padded) blood-vessel voxel, not
    // necessarily at the origin. Callers that only care about cell/vessel
    // geometry (e.g. drawing a bounding wireframe) must use these rather than
    // assuming an origin-anchored box of edge Parameters::voxel_size.
    using SimResult = std::tuple<std::vector<Axon>,
                                    std::vector<Blood_Vessel>,
                                    std::vector<Glial>,
                                    std::vector<Glial>,
                                    std::vector<Glial>,
                                    Eigen::Vector3d,
                                    Eigen::Vector3d>;
    using ProgressCallback = std::function<void(double, double)>;

    // This function reads the JSON, fills the struct, and calls runSimulation()
    static void runSimulationFromJson(const std::string& jsonFilePath);

    // Inverse of runSimulationFromJson's parsing: writes params back out as a
    // config JSON with the exact same keys/sections, so the file this
    // produces can itself be fed back in via --config. Used by the GUI to
    // save a record of the exact parameters behind each run, alongside its
    // .log (see Window::StartSimulation).
    static void writeParametersToJson(const Parameters& params, const std::string& jsonFilePath);

    // on_growth_progress: (depth travelled, total box depth), fired once per
    // depth layer. on_swelling_progress: (current axon ICVF, target axon
    // ICVF), fired once per swelling round. Both optional -- left unset
    // (default-constructed empty std::function), CLI callers are unaffected.
    static SimResult runSimulation(const Parameters& params,
                                    ProgressCallback on_growth_progress = nullptr,
                                    ProgressCallback on_swelling_progress = nullptr);
};