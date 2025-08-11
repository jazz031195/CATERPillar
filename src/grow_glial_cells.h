#ifndef GLIALCELLGROWTH_H
#define GLIALCELLGROWTH_H

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
class GlialCellGrowth : public CellGrowth

{
    public:

        Glial &glial_cell_to_grow;

        GlialCellGrowth();

        ~GlialCellGrowth();

        GlialCellGrowth(const GlialCellGrowth &other);

        GlialCellGrowth(Glial &glial_cell_to_grow_,
            const std::vector<Glial> &astrocytes,
            const std::vector<Glial> &oligodendrocytes,
            const std::vector<Axon> &axons_,
            const Eigen::Vector3d &extended_min_limits_,
            const Eigen::Vector3d &extended_max_limits_,
            const Eigen::Vector3d &min_limits_,
            const Eigen::Vector3d &max_limits_,
            const double &std_dev_,
            const double &min_radius_);

        void add_spheres(Sphere &sph, const Sphere &last_sphere, const bool &check_collision_with_branches, const int &factor, const int &index_ram_spheres);
        bool AddOneSphere(const double &radius_, const bool &create_sphere, int &grow_straight, const int &i, const bool &check_collision_with_branches, const int &parent, const int &factor);
        bool collideswithItself(Sphere &sph);
        void find_next_center_straight(double distance, Sphere &s, const std::vector<Sphere> &spheres);
        void find_next_center(Sphere &s,  double dist_, const std::vector<Sphere> &spheres, const Eigen::Vector3d &target);
        bool canSpherebePlaced(Sphere &sph);
        Eigen::Vector3d readapt_sphere_position(const Sphere &s, const Axon &neighbour_axon, bool can_readapt);
        bool checkAxonsOverlap(Sphere &sph);

};

#endif // GLIALCELLGROWTH_H


    

