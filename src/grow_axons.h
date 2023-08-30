
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
    std::vector<Axon> axons;
    std::vector<Axon> axons_to_regrow;

    Axon axon_to_grow; /* Axon vector */
    Eigen::Vector3d voxel_size;
    bool tortuous;
    Dynamic_Sphere sphere_to_add;
    bool finished;
    double max_radius; /* initial radius */ 
    bool grow_straight; /* Sometimes, you want the axon to grow staright and sometimes you want it to grow in a different random direction */
    double radius;
    Growth(){};
    ~Growth(){};

    Growth(Axon, std::vector<Axon>, std::vector<Axon>,  Eigen::Vector3d, bool, double, double, int);

    bool GrowAxon();
    bool GrowFirstSphere();
    void find_next_center(Dynamic_Sphere &s,double dist_);
    //bool isSphereColliding(Dynamic_Sphere sph);
    bool isSphereColliding_(Dynamic_Sphere sph);
    std::vector<int> checkAxisForCollision(Dynamic_Sphere sph, int axis);
    bool check_borders(Eigen::Vector3d pos, double distance_to_border, Eigen::Vector2d &twin_delta_pos);
    void find_next_center_straight(double distance, Dynamic_Sphere &s);
    bool TestGrowAxonAtPos(Eigen::Vector3d position_to_test);
    bool TestGrowAxon(Eigen::Vector3d &position_that_worked);
    void initialise_env_axons();


};

#endif
