#ifndef CELLGROWTH_H
#define CELLGROWTH_H

#include <vector>
#include <atomic>
#include <Eigen/Dense>
#include "Axon.h"
#include "Glial.h"
#include "sphere.h"
#include "threads.h"

/**
 * @brief The Growth class handles the logic for growing axons or glial cells
 *        within specified bounding limits, checking collisions against
 *        other axons and glial cells.
 */
class CellGrowth
{
public:

    std::vector<Axon> axons;
    std::vector<Glial> glial_cells;
    Eigen::Vector3d min_limits;
    Eigen::Vector3d max_limits;
    Eigen::Vector3d extended_min_limits;
    Eigen::Vector3d extended_max_limits;
    double std_dev;
    double min_radius;
    int grow_straight;
    int finished;

    CellGrowth();
  
    CellGrowth(const std::vector<Axon>& axons_, const std::vector<Glial>& astrocytes, const std::vector<Glial>& oligodendrocytes, const Eigen::Vector3d &extended_min_limits_, const Eigen::Vector3d& extended_max_limits_, const Eigen::Vector3d& min_limits_, const Eigen::Vector3d& max_limits_, const double& std_dev_, const double& min_radius_)
    : axons(axons_), glial_cells(astrocytes), min_limits(min_limits_), max_limits(max_limits_), extended_min_limits(extended_min_limits_), extended_max_limits(extended_max_limits_), std_dev(std_dev_), min_radius(min_radius_), finished(false)
    {
        glial_cells.insert(glial_cells.end(), oligodendrocytes.begin(), oligodendrocytes.end());

        if (std_dev == 0.0)
        {
            grow_straight = 0;
        }

    }

    virtual ~CellGrowth();

    CellGrowth(const CellGrowth &other);

    bool check_borders(const Eigen::Vector3d&  min_l, const Eigen::Vector3d&  max_l, const Eigen::Vector3d& pos, const double& distance_to_border);

    double clamp(double value, double lower, double upper);

    Eigen::Vector3d generate_random_point_on_sphere(double std);

    Eigen::Matrix3d rotation_matrix_from_vectors(const Eigen::Vector3d &vec1, const Eigen::Vector3d &vec2);

    Eigen::Vector3d apply_bias_toward_target(const Eigen::Vector3d &point, const Eigen::Vector3d &target);

    bool canSpherebePlaced(Sphere sph);
        

};

#endif // CELLGROWTH_H
