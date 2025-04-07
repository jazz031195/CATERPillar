#ifndef CELLGROWTH_H
#define CELLGROWTH_H

#include <vector>
#include <Eigen/Dense>
#include "Axon.h"
#include "Glial.h"
#include "sphere.h"
#include "threads.h"

/**
 * @brief Abstract base class for cell growth (axon or glial).
 * Handles logic for growth within spatial constraints and collision checks.
 */
class CellGrowth
{
public:
    virtual ~CellGrowth(); 
    CellGrowth();
    CellGrowth(const CellGrowth& other);

    int finished = 0;

    bool check_borders(const Eigen::Vector3d& min_l, const Eigen::Vector3d& max_l, const Eigen::Vector3d& pos, const double& distance_to_border);
    double clamp(double value, double lower, double upper);
    Eigen::Vector3d generate_random_point_on_sphere(double std);
    Eigen::Matrix3d rotation_matrix_from_vectors(const Eigen::Vector3d& vec1, const Eigen::Vector3d& vec2);
    Eigen::Vector3d apply_bias_toward_target(const Eigen::Vector3d& point, const Eigen::Vector3d& target);
    bool canSpherebePlaced(const Sphere& sph);

protected:
    std::vector<Glial> glial_cells;
    const std::vector<Axon>& axons;  

    Eigen::Vector3d min_limits;
    Eigen::Vector3d max_limits;
    Eigen::Vector3d extended_min_limits;
    Eigen::Vector3d extended_max_limits;

    double std_dev;
    double min_radius;
    int grow_straight = 0;

    CellGrowth(const std::vector<Axon>& axons_,
               const std::vector<Glial>& astrocytes,
               const std::vector<Glial>& oligodendrocytes,
               const Eigen::Vector3d& extended_min_limits_,
               const Eigen::Vector3d& extended_max_limits_,
               const Eigen::Vector3d& min_limits_,
               const Eigen::Vector3d& max_limits_,
               const double& std_dev_,
               const double& min_radius_)
        : axons(axons_),
          min_limits(min_limits_),
          max_limits(max_limits_),
          extended_min_limits(extended_min_limits_),
          extended_max_limits(extended_max_limits_),
          std_dev(std_dev_),
          min_radius(min_radius_),
          grow_straight(std_dev_ == 0.0 ? 0 : 1),
          finished(0)
    {
        glial_cells = astrocytes;
        glial_cells.insert(glial_cells.end(), oligodendrocytes.begin(), oligodendrocytes.end());
    }
};

#endif // CELLGROWTH_H
