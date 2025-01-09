#ifndef AXONGROWTH_H
#define AXONGROWTH_H

#include <vector>
#include <atomic>
#include <Eigen/Dense>
#include "Axon.h"
#include "Glial.h"
#include "grow_cells.h"
#include "sphere.h"
#include "threads.h"

/**
 * @brief The Growth class handles the logic for growing axons or glial cells
 *        within specified bounding limits, checking collisions against
 *        other axons and glial cells.
 */
class AxonGrowth : public CellGrowth
{
public:

    Axon &axon_to_grow;

    /**
     * @brief Constructor for growing an axon.
     * @param axon_to_grow_         The axon to grow.
     * @param astrocytes            Existing astrocytes.
     * @param oligodendrocytes      Existing oligodendrocytes.
     * @param axons_                Reference to the global axons container 
     *                              (modified in place).
     * @param extended_min_limits_  Extended bounding box minimum.
     * @param extended_max_limits_  Extended bounding box maximum.
     * @param min_limits_           Tight bounding box minimum.
     * @param max_limits_           Tight bounding box maximum.
     * @param std_dev_              Standard deviation for random growth angle.
     * @param min_radius_           Minimum allowed radius for growth.
     */

    AxonGrowth();

    AxonGrowth(Axon &axon_to_grow_,
               const std::vector<Glial>& astrocytes,
               const std::vector<Glial>& oligodendrocytes,
               const std::vector<Axon>& axons_,
               const Eigen::Vector3d &extended_min_limits_,
               const Eigen::Vector3d &extended_max_limits_,
               const Eigen::Vector3d &min_limits_,
               const Eigen::Vector3d &max_limits_,
               const double &std_dev_,
               const double &min_radius_);

    ~AxonGrowth();

    AxonGrowth(const AxonGrowth &other);


    /**
     * @brief Attempts to push other axons away if the newly placed sphere
     *        collides with them.
     * @param sph The sphere that caused collision.
     * @return true if successful, false otherwise.
     */
    bool PushOtherAxons(const Sphere &sph);

    /**
     * @brief Pushes the spheres of a specific axon to resolve collisions.
     * @return true if push succeeded, false if collision remains.
     */
    bool pushAxonSpheres(Axon &axon, const Sphere &sph);

    /**
     * @brief Checks if the sphere collides with an axon.
     */
    bool isSphereCollidingwithAxon(Axon ax, Sphere sph);


    /**
     * @brief Attempts to place a new sphere for the axon.
     * @param radius_       The desired radius for the new sphere.
     * @param create_sphere If true, actually add the sphere once placed.
     * @param grow_straight Flag controlling whether to grow in a straight line.
     * @param factor        Possibly used for intermediate spheres.
     * @return true if successful, false otherwise.
     */
    bool AddOneSphere(double radius_,
                      bool create_sphere,
                      int grow_straight,
                      const int &factor);


    /**
     * @brief Adds intermediate spheres between last_sphere and sph for an axon.
     */
    void add_spheres(Sphere &sph,
                         const Sphere &last_sphere,
                         const int &factor);

    void find_next_center_straight(double distance, Sphere &s, const std::vector<Sphere> &spheres);
    double myelin_thickness_(const double &inner_radius);
    double originalFunction_(double x, double outerRadius);
    double derivative_(double x);
    double findInnerRadius_(const double &outerRadius);
    void find_next_center(Sphere &s,  double dist_, const std::vector<Sphere> &spheres, const Eigen::Vector3d &target, const bool &is_axon);

};

#endif // AXONGROWTH_H
