
#ifndef GROWAXONS_H
#define GROWAXONS_H

#include "Eigen/Core"
#include <vector>
#include "constants.h"
#include <iostream>
#include "Axon.h"
#include "dynamic_sphere.h"





class Growth
{
public:

    std::vector<Axon> env_axons; 
    Axon* axon_to_grow;                       /*!< Axon vector                                                            */
    Eigen::Vector3d voxel_size;
    bool tortuous;
    Dynamic_Sphere sphere_to_add;
    std::vector<Eigen::Vector3d> centers;
    bool finished;
    double max_radius;
    bool grow_straight;                     /*Sometimes, you want the axon to grow staright and sometimes you want it to grow in a different random direction*/


    Growth(){};
    ~Growth();

    Growth (Axon*,  std::vector<Axon>, Eigen::Vector3d, bool, double, bool);
    
    bool GrowAxon();
    bool GrowFirstSphere();
    void add_next_sphere(Dynamic_Sphere added_sphere, std::vector<Eigen::Vector3d>& centers, std::vector<double>& sph_radii);
    void find_next_center( Dynamic_Sphere& s, vector<Eigen::Vector3d> centers, double dist_, double rad);
    bool isSphereColliding(Dynamic_Sphere sph);
    bool check_borders(Eigen::Vector3d pos, double distance_to_border, Eigen::Vector2d& twin_delta_pos);
    void find_next_center_straight(vector<Eigen::Vector3d> centers, double distance,Dynamic_Sphere& s);
};

#endif 
