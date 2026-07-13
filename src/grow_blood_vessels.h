#ifndef BLOODVESSELGROWTH_H
#define BLOODVESSELGROWTH_H

#include <vector>
#include <functional>
#include <Eigen/Dense>
#include "Axon.h"
#include "Blood_Vessel.h"
#include "grow_cells.h"
#include "sphere.h"
#include "threads.h"

/**
 * @brief Handles the logic for axon growth including sphere placement and
 *        optional myelination logic. Inherits from CellGrowth.
 */
class BloodVesselGrowth : public CellGrowth
{
public:
    Blood_Vessel& bv_to_grow;  /*!< Reference to the Blood_Vessel being grown */

    BloodVesselGrowth() = delete;

    BloodVesselGrowth(Blood_Vessel &bv_to_grow_,
                       const SphereGrid* sphere_grid_,
                       const Eigen::Vector3d &extended_min_limits_,
                       const Eigen::Vector3d &extended_max_limits_,
                       const Eigen::Vector3d &min_limits_,
                       const Eigen::Vector3d &max_limits_,
                       const double &epsilon_,
                       const double &min_radius_);


    ~BloodVesselGrowth();

    BloodVesselGrowth(const BloodVesselGrowth& other);

        // Growth and placement (main vessel, branch 0)
    bool AddOneSphere(double radius_,
                      bool create_sphere,
                      int grow_straight,
                      const int& factor);

    void add_spheres(Sphere& sph,
                     const Sphere& last_sphere,
                     const int& factor);

    // Positioning (main vessel, branch 0)
    Eigen::Vector3d find_next_center_straight(const double distance,
                                   const std::vector<Sphere>& spheres);

    Eigen::Vector3d find_next_center(const double dist_,
                          const std::vector<Sphere>& spheres,
                          const Eigen::Vector3d& target);

    void growthThread(double& stuck_radius, int& stuck_index, int factor, bool axon_can_shrink);

    double RandomradiusVariation();
    bool shrinkRadius(const double &radius_to_shrink, const bool& axon_can_shrink, const int &factor);
    void update_straight(bool can_grow_, int &grow_straight, int &straight_growths);

    // Branch growth (branch_id >= 1), mirrors GlialCellGrowth's secondary-branch machinery
    bool AddOneSphere(const double &radius_, const bool &create_sphere, int &grow_straight, const int &i, const bool &check_collision_with_branches, const int &parent, const int &factor);
    void add_spheres(Sphere &sph, const Sphere &last_sphere, const bool &check_collision_with_branches, const int &factor, const int &index_ram_spheres);
    void find_next_center(Sphere &s, double dist_, const std::vector<Sphere> &spheres, const Eigen::Vector3d &target);
    bool GenerateFirstSphereinProcess(Sphere &first_sphere, Eigen::Vector3d &attractor, const double &radius, const Sphere &sphere_to_emerge_from, const Eigen::Vector3d &vector_to_prev_center, const int &nbr_spheres, const int &nbr_spheres_between, const int &vessel_id, const int &branch_id);
    std::vector<Sphere> addIntermediateSpheres(const Sphere &random_sphere, const Sphere &first_sphere, const int &branch_nbr, const int &nbr_spheres, const int &nbr_spheres_between, const std::function<double(double)> &compute_radius, const double &t_start, const double &t_end);
    bool growBranch(int &nbr_spheres, const int &factor);

};

#endif // BLOODVESSELGROWTH_H
