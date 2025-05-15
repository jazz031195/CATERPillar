#ifndef AXONGROWTH_H
#define AXONGROWTH_H

#include <vector>
#include <Eigen/Dense>
#include "Axon.h"
#include "Glial.h"
#include "grow_cells.h"
#include "sphere.h"
#include "threads.h"

/**
 * @brief Handles the logic for axon growth including sphere placement and 
 *        optional myelination logic. Inherits from CellGrowth.
 */
class AxonGrowth : public CellGrowth
{
public:
    Axon& axon_to_grow;
    std::vector<Axon> axons;  

    AxonGrowth() = delete;

    AxonGrowth(Axon& axon_to_grow_,
               const std::vector<Glial>& astrocytes,
               const std::vector<Glial>& oligodendrocytes,
               const std::vector<Axon>& axons_,
               const Eigen::Vector3d& extended_min_limits_,
               const Eigen::Vector3d& extended_max_limits_,
               const Eigen::Vector3d& min_limits_,
               const Eigen::Vector3d& max_limits_,
               const double& std_dev_,
               const double& min_radius_);

    ~AxonGrowth();
    AxonGrowth(const AxonGrowth& other);

    // Growth and placement
    bool AddOneSphere(double radius_,
                      bool create_sphere,
                      int grow_straight,
                      const int& factor);

    void add_spheres(Sphere& sph,
                     const Sphere& last_sphere,
                     const int& factor);

    // Positioning
    Eigen::Vector3d find_next_center_straight(const double distance,
                                   const Sphere& s,
                                   const std::vector<Sphere>& spheres);

    Eigen::Vector3d find_next_center(const Sphere& s,
                          const double dist_,
                          const std::vector<Sphere>& spheres,
                          const Eigen::Vector3d& target,
                          const bool& is_axon);

    // Collision detection
    bool isSphereCollidingwithAxon(Axon ax, Sphere sph);

    // Override if needed
    bool pushAxonSpheres(Axon &axon, const Sphere &sph);

    bool canSpherebePlaced(Sphere &sph);

    bool checkAxonsOverlap(Sphere &sph);
    Eigen::Vector3d readapt_sphere_position(const Sphere &s, const Axon &neighbour_axon, bool can_readapt);
    Eigen::Vector3d find_closest_neighbour(const Sphere &sphere);
};

#endif // AXONGROWTH_H
