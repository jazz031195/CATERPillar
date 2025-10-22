#ifndef CELLGROWTH_H
#define CELLGROWTH_H

#include <vector>
#include <Eigen/Dense>
#include "Axon.h"
#include "Glial.h"
#include "sphere.h"
#include "threads.h"
#include <random>

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
    void update_environment(const std::vector<Axon>* axons_,
                                const std::vector<Glial>* glial_pop1_,
                                const std::vector<Glial>* glial_pop2_) noexcept;
    bool check_borders(const Eigen::Vector3d& min_l, const Eigen::Vector3d& max_l, const Eigen::Vector3d& pos, const double& distance_to_border);
    double clamp(double value, double lower, double upper);
    Eigen::Vector3d generate_random_point_on_sphere(double std);
    Eigen::Matrix3d rotation_matrix_from_vectors(const Eigen::Vector3d& vec1, const Eigen::Vector3d& vec2);
    Eigen::Vector3d apply_bias_toward_target(const Eigen::Vector3d& point, const Eigen::Vector3d& target);
    bool canSpherebePlaced(Sphere& sph);
    bool checkAxonsOverlap(Sphere &sph);
    const std::vector<Axon>& AX() const;
    const std::vector<Axon>* AXptr() const;

    bool finished = false;

protected:

    const std::vector<Glial>* glial_pop1;
    const std::vector<Glial>* glial_pop2;
    const std::vector<Axon>* axons;  

    Eigen::Vector3d min_limits;
    Eigen::Vector3d max_limits;
    Eigen::Vector3d extended_min_limits;
    Eigen::Vector3d extended_max_limits;

    double std_dev;
    double min_radius;
    int grow_straight;

    std::mt19937 gen;

    CellGrowth(const std::vector<Axon>* axons_,
               const std::vector<Glial>* glial_pop1_,
               const std::vector<Glial>* glial_pop2_,
               const Eigen::Vector3d& extended_min_limits_,
               const Eigen::Vector3d& extended_max_limits_,
               const Eigen::Vector3d& min_limits_,
               const Eigen::Vector3d& max_limits_,
               const double& std_dev_,
               const double& min_radius_)
        : min_limits(min_limits_),
          max_limits(max_limits_),
          extended_min_limits(extended_min_limits_),
          extended_max_limits(extended_max_limits_),
          std_dev(std_dev_),
          min_radius(min_radius_),
          grow_straight(std_dev_ == 0.0 ? 0 : 1),
          finished(false),
          axons(axons_),
          glial_pop1(glial_pop1_),
          glial_pop2(glial_pop2_) {
        std::random_device rd;
        gen.seed(rd());
    }
};

#endif // CELLGROWTH_H
