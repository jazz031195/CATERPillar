
#ifndef GROWAXONS_H
#define GROWAXONS_H

#include "Eigen/Core"
#include <vector>
#include "constants.h"
#include <iostream>
#include "Axon.h"
#include "sphere.h"

class Growth
{
public:

    std::vector<Axon> env_axons;            /* Only axosn that are close to the one growing */
    std::vector<Axon> axons;                /* All axons in environment */
    Axon axon_to_grow;                      /* Axon to grow */
    Eigen::Vector3d voxel_size;             /* Voxel size in um */ 
    bool tortuous;                          /* True if the axon is in general tortuous */ 
    Sphere sphere_to_add;                   /* Sphere to add to axon */ 
    bool finished;                          /* True if the axon has finished growing */ 
    double max_radius;                      /* Maximum radius */ 
    bool grow_straight;                     /* True if the next sphere should grow straight */


    Growth(){};
    ~Growth(){};

    Growth(Axon&, std::vector<Axon>,  Eigen::Vector3d, bool, double, int);

    /*!
     *  \brief Adds sphere to axon
     */
    void initialise_env_axons();
    bool GrowAxon(double radius, bool create_sphere);
    bool GrowFirstSphere();
    void find_next_center(Sphere &s,double dist_);
    void find_next_center_straight(double distance, Sphere &s);
    bool isSphereColliding(Sphere sph);
    std::vector<int> checkAxisForCollision(Sphere sph, int axis);
    bool check_borders(Eigen::Vector3d pos, double distance_to_border);
    

};

#endif
