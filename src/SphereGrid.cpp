#include "SphereGrid.h"
#include "constants.h"
#include <algorithm>
#include <cmath>
#include <limits>

SphereGrid::SphereGrid()
    : min_limits(Eigen::Vector3d::Zero()), max_limits(Eigen::Vector3d::Zero()), voxel_size(1.0), nx(1), ny(1), nz(1)
{
    voxels.resize(1);
}

SphereGrid::SphereGrid(const Eigen::Vector3d &min_limits_, const Eigen::Vector3d &max_limits_, const double &voxel_size_)
    : min_limits(min_limits_), max_limits(max_limits_), voxel_size(voxel_size_)
{
    Eigen::Vector3d extent = max_limits - min_limits;
    nx = std::max(1, static_cast<int>(std::ceil(extent[0] / voxel_size)));
    ny = std::max(1, static_cast<int>(std::ceil(extent[1] / voxel_size)));
    nz = std::max(1, static_cast<int>(std::ceil(extent[2] / voxel_size)));
    voxels.resize(static_cast<size_t>(nx) * static_cast<size_t>(ny) * static_cast<size_t>(nz));
}

void SphereGrid::clear()
{
    for (auto &voxel : voxels) {
        voxel.clear();
    }
}

Eigen::Vector3i SphereGrid::clampedVoxelIndex(const Eigen::Vector3d &position) const
{
    int ix = static_cast<int>(std::floor((position[0] - min_limits[0]) / voxel_size));
    int iy = static_cast<int>(std::floor((position[1] - min_limits[1]) / voxel_size));
    int iz = static_cast<int>(std::floor((position[2] - min_limits[2]) / voxel_size));

    ix = std::min(std::max(ix, 0), nx - 1);
    iy = std::min(std::max(iy, 0), ny - 1);
    iz = std::min(std::max(iz, 0), nz - 1);

    return Eigen::Vector3i(ix, iy, iz);
}

int SphereGrid::linearIndex(int ix, int iy, int iz) const
{
    return ix + iy * nx + iz * nx * ny;
}

void SphereGrid::insert(const Sphere &sph)
{
    Eigen::Vector3d radius_vec(sph.radius, sph.radius, sph.radius);
    Eigen::Vector3i min_idx = clampedVoxelIndex(sph.center - radius_vec);
    Eigen::Vector3i max_idx = clampedVoxelIndex(sph.center + radius_vec);

    SphereGridEntry entry{sph.id, sph.object_id, sph.object_type, sph.center, sph.radius, sph.branch_id};

    for (int ix = min_idx[0]; ix <= max_idx[0]; ++ix) {
        for (int iy = min_idx[1]; iy <= max_idx[1]; ++iy) {
            for (int iz = min_idx[2]; iz <= max_idx[2]; ++iz) {
                voxels[linearIndex(ix, iy, iz)].push_back(entry);
            }
        }
    }
}

void SphereGrid::remove(const Sphere &sph)
{
    Eigen::Vector3d radius_vec(sph.radius, sph.radius, sph.radius);
    Eigen::Vector3i min_idx = clampedVoxelIndex(sph.center - radius_vec);
    Eigen::Vector3i max_idx = clampedVoxelIndex(sph.center + radius_vec);

    for (int ix = min_idx[0]; ix <= max_idx[0]; ++ix) {
        for (int iy = min_idx[1]; iy <= max_idx[1]; ++iy) {
            for (int iz = min_idx[2]; iz <= max_idx[2]; ++iz) {
                auto &voxel = voxels[linearIndex(ix, iy, iz)];
                voxel.erase(
                    std::remove_if(voxel.begin(), voxel.end(),
                        [&](const SphereGridEntry &entry) {
                            return entry.sphere_id == sph.id &&
                                   entry.cell_id == sph.object_id &&
                                   entry.object_id == sph.object_type;
                        }),
                    voxel.end());
            }
        }
    }
}

std::vector<SphereGridEntry> SphereGrid::query(const Eigen::Vector3d &center, const double &radius) const
{
    std::vector<SphereGridEntry> results;

    Eigen::Vector3d radius_vec(radius, radius, radius);
    Eigen::Vector3i min_idx = clampedVoxelIndex(center - radius_vec);
    Eigen::Vector3i max_idx = clampedVoxelIndex(center + radius_vec);

    for (int ix = min_idx[0]; ix <= max_idx[0]; ++ix) {
        for (int iy = min_idx[1]; iy <= max_idx[1]; ++iy) {
            for (int iz = min_idx[2]; iz <= max_idx[2]; ++iz) {
                const auto &voxel = voxels[linearIndex(ix, iy, iz)];
                results.insert(results.end(), voxel.begin(), voxel.end());
            }
        }
    }

    return results;
}

// Cheap per-axis interval overlap test (no sqrt).
// Lets canSpherebePlaced reject most candidates before paying for the exact distance check.
static bool axisIntervalsOverlap(double c1, double r1, double c2, double r2, double tolerance)
{
    double min1 = c1 - r1;
    double max1 = c1 + r1 + tolerance;
    double min2 = c2 - r2;
    double max2 = c2 + r2 + tolerance;
    return !(max1 < min2 || max2 < min1);
}

bool SphereGrid::canSpherebePlaced(const Sphere &sph, bool check_collision_with_branches, int extra_excluded_branch_id) const
{
    // Deliberately does not go through query(): this is the single most
    // frequently called function in the whole codebase (every growth
    // attempt, every push, every shrink-fallback and swelling bisection
    // iteration, seedAllAxons, checkNoCollisions), and query() would
    // heap-allocate a fresh std::vector and copy in every candidate before
    // any of them are even checked. Iterating voxels/entries directly here
    // instead means the common case -- a real collision found early, or no
    // entries at all in range -- never touches the heap, and returns the
    // instant a genuine overlap is confirmed rather than after collecting
    // every candidate first.
    double tolerance = barrier_tickness;

    Eigen::Vector3d radius_vec(sph.radius, sph.radius, sph.radius);
    Eigen::Vector3i min_idx = clampedVoxelIndex(sph.center - radius_vec);
    Eigen::Vector3i max_idx = clampedVoxelIndex(sph.center + radius_vec);

    for (int ix = min_idx[0]; ix <= max_idx[0]; ++ix) {
        for (int iy = min_idx[1]; iy <= max_idx[1]; ++iy) {
            for (int iz = min_idx[2]; iz <= max_idx[2]; ++iz) {
                const auto &voxel = voxels[linearIndex(ix, iy, iz)];
                for (const auto &entry : voxel) {
                    // skip spheres belonging to the same object (self)
                    if (entry.object_id == sph.object_type &&
                        entry.cell_id == sph.object_id &&
                        (!check_collision_with_branches || entry.branch_id == sph.branch_id || entry.branch_id == extra_excluded_branch_id)) {
                        continue;
                    }

                    bool boxes_overlap = true;
                    for (int axis = 0; axis < 3 && boxes_overlap; ++axis) {
                        if (!axisIntervalsOverlap(entry.center[axis], entry.radius, sph.center[axis], sph.radius, tolerance)) {
                            boxes_overlap = false;
                        }
                    }
                    if (!boxes_overlap) {
                        continue;
                    }

                    double dist = (entry.center - sph.center).norm();
                    if (dist <= entry.radius + sph.radius + tolerance) {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

bool SphereGrid::findWorstOverlap(const Eigen::Vector3d &center, double radius, double search_radius,
                                   int self_object_type, int self_object_id,
                                   Eigen::Vector3d &blocker_center, double &blocker_radius) const
{
    Eigen::Vector3d search_vec(search_radius, search_radius, search_radius);
    Eigen::Vector3i min_idx = clampedVoxelIndex(center - search_vec);
    Eigen::Vector3i max_idx = clampedVoxelIndex(center + search_vec);

    double worst_overlap = 0.0;
    bool found = false;

    for (int ix = min_idx[0]; ix <= max_idx[0]; ++ix) {
        for (int iy = min_idx[1]; iy <= max_idx[1]; ++iy) {
            for (int iz = min_idx[2]; iz <= max_idx[2]; ++iz) {
                const auto &voxel = voxels[linearIndex(ix, iy, iz)];
                for (const auto &entry : voxel) {
                    if (entry.object_id == self_object_type && entry.cell_id == self_object_id) {
                        continue; // this object's own spheres
                    }
                    double d = (entry.center - center).norm();
                    double overlap = (radius + entry.radius) - d;
                    if (overlap > worst_overlap) {
                        worst_overlap = overlap;
                        blocker_center = entry.center;
                        blocker_radius = entry.radius;
                        found = true;
                    }
                }
            }
        }
    }
    return found;
}

bool SphereGrid::findNearest(const Eigen::Vector3d &center, double search_radius,
                              int self_object_type, int self_object_id,
                              Eigen::Vector3d &nearest_center, double &nearest_radius) const
{
    Eigen::Vector3d search_vec(search_radius, search_radius, search_radius);
    Eigen::Vector3i min_idx = clampedVoxelIndex(center - search_vec);
    Eigen::Vector3i max_idx = clampedVoxelIndex(center + search_vec);

    double best_dist = std::numeric_limits<double>::max();
    bool found = false;

    for (int ix = min_idx[0]; ix <= max_idx[0]; ++ix) {
        for (int iy = min_idx[1]; iy <= max_idx[1]; ++iy) {
            for (int iz = min_idx[2]; iz <= max_idx[2]; ++iz) {
                const auto &voxel = voxels[linearIndex(ix, iy, iz)];
                for (const auto &entry : voxel) {
                    if (entry.object_id == self_object_type && entry.cell_id == self_object_id) {
                        continue; // this object's own spheres
                    }
                    double d = (entry.center - center).norm();
                    if (d < best_dist) {
                        best_dist = d;
                        nearest_center = entry.center;
                        nearest_radius = entry.radius;
                        found = true;
                    }
                }
            }
        }
    }
    return found;
}
