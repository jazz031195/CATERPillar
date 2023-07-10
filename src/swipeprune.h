//!  Projections Class =============================================================/
/*!
*   \details   Class to construct projections of all objects on each of the axes. This speeds up the collision check 
*   \author    Jasmine Nguyen-Duc
*   \date      January 2023
*   \version   0.2
=================================================================================================*/

#ifndef PROJECTIONS_H
#define PROJECTIONS_H

#include "Eigen/Core"
#include <vector>
#include "constants.h"
#include <iostream>


class Projections
{
public:

    struct projection_pt{
        double position;
        int axon_id;
        int sph_id;
    };

    // projection of spheres in x, y, and z axes
    std::vector<projection_pt> sph_projections_x;
    std::vector<projection_pt> sph_projections_y;
    std::vector<projection_pt> sph_projections_z;
    
    Projections();
    void clear_projections();
    std::vector<projection_pt> find_collisions(projection_pt proj_on_axis_min, projection_pt proj_on_axis_max,std::vector<projection_pt> projections_on_axis, double distance_to_be_inside);
    std::vector<std::vector<projection_pt>> find_collisions_all_axes(Eigen::Vector3d &position, double rad, int ax_id, double distance_to_be_inside);
    bool isProjInside(std::vector<projection_pt> projs, projection_pt p);
    void append_right_place(projection_pt p1, projection_pt p2, int axis);
};

#endif // PROJECTIONS_H
