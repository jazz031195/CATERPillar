
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
     *  \brief Initialises environment axons that are close by
     */
    void initialise_env_axons();

    /*!
     *  \param radius Radius of sphere to add
     *  \param create_sphere Whether to add the sphere to the list of spheres or not
     *  \brief Grows axon by adding another sphere to it
     */
    bool GrowAxon(double radius, bool create_sphere);
    /*!
     *  \param dist_ distance to add the sphere to (from last sphere)
     *  \param s sphere to add
     *  \brief Finds position of center of sphere to add
     */
    void find_next_center(Sphere &s, double dist_);
    /*!
     *  \param dist_ distance to add the sphere to (from last sphere)
     *  \param s sphere to add
     *  \brief Finds position of center of sphere to add if it must be added in a straight line
     */
    void find_next_center_straight(double distance, Sphere &s);

    /*!
     *  \param sh sphere to check
     *  \brief Checks if sphere is colliding with environment
     */
    bool isSphereColliding(Sphere sph);
    /*!
     *  \param sh sphere to check
     *  \param axis axis number
     *  \brief Checks if sphere is colliding with environment with respect to an axis
     */
    std::vector<int> checkAxisForCollision(Sphere sph, int axis);
    /*!
     *  \param distance_to_border dictance to be inside voxel
     *  \param pos position to check
     *  \brief Checks if position is within the voxel
     */
    bool check_borders(Eigen::Vector3d pos, double distance_to_border);

    /*!
     *  \param sph sphere to check
     *  \param ax axon to check
     *  \brief Checks if sph collides with specific axon
     */
    bool isSphereCollidingwithAxon(Axon ax, Sphere sph);
    

};

#endif
