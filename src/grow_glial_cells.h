#ifndef GLIALCELLGROWTH_H
#define GLIALCELLGROWTH_H

#include <vector>
#include <atomic>
#include <functional>
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
class GlialCellGrowth : public CellGrowth

{
    public:

        Glial &glial_cell_to_grow;

        GlialCellGrowth();

        ~GlialCellGrowth();

        GlialCellGrowth(const GlialCellGrowth &other);

        GlialCellGrowth(Glial &glial_cell_to_grow_,
            const SphereGrid* sphere_grid_,
            const Eigen::Vector3d &extended_min_limits_,
            const Eigen::Vector3d &extended_max_limits_,
            const Eigen::Vector3d &min_limits_,
            const Eigen::Vector3d &max_limits_,
            const double &min_radius_, const double &epsilon_= 0.4);

        void add_spheres(Sphere &sph, const Sphere &last_sphere, const bool &check_collision_with_branches, const int &factor, const int &index_ram_spheres);
        bool AddOneSphere(const double &radius_, const bool &create_sphere, int &grow_straight, const int &i, const bool &check_collision_with_branches, const int &parent, const int &factor);
        void find_next_center_straight(double distance, Sphere &s, const std::vector<Sphere> &spheres);
        /*!
         *  \param std_override When >= 0, used in place of the class's own
         *         epsilon as the angular spread for this call's direction
         *         sample -- lets AddOneSphere widen the search on later
         *         retries (see its own doc comment) instead of resampling
         *         the same narrow, attractor-biased cone every time.
         */
        void find_next_center(Sphere &s,  double dist_, const std::vector<Sphere> &spheres, const Eigen::Vector3d &target, double std_override = -1.0);
        void growFirstPrimaryBranches(const int &number_ramification_points, int &nbr_spheres, const double &mean_process_length, const double &std_process_length, const int &factor);
        bool growPrimaryBranch(int &nbr_spheres, const double &mean_primary_process_length, const double &std_primary_process_length, const int &factor);
        bool growSecondaryBranch(int &nbr_spheres, const double &mean_process_length, const double &std_process_length, const int &factor);
        bool GenerateFirstSphereinProcess(Sphere &first_sphere, Eigen::Vector3d &attractor, const double &radius, const Sphere &sphere_to_emerge_from, const Eigen::Vector3d &vector_to_prev_center, const int &nbr_spheres, const int &nbr_spheres_between, const int &cell_id, const int &branch_id, const bool &primary_process);
        std::vector<Sphere> addIntermediateSpheres(const Sphere &random_sphere, const Sphere &first_sphere,  const int &branch_nbr, const int &nbr_spheres, const int &nbr_spheres_between, const std::function<double(double)> &compute_radius, const double &t_start, const double &t_end);

};

#endif // GLIALCELLGROWTH_H


    

