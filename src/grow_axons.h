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
    std::vector<Axon> modified_axons;  
    std::vector<int> ind_axons_pushed;

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
    void find_next_center_straight(double distance,
                                   Sphere& s,
                                   const std::vector<Sphere>& spheres);

    void find_next_center(Sphere& s,
                          double dist_,
                          const std::vector<Sphere>& spheres,
                          const Eigen::Vector3d& target,
                          const bool& is_axon);

    // Collision detection
    bool isSphereCollidingwithAxon(Axon ax, Sphere sph);

    // Override if needed
    bool pushAxonSpheres(Axon &axon, const Sphere &sph);

    bool canSpherebePlaced(const Sphere &sph, const bool &with_push);
};

#endif // AXONGROWTH_H
