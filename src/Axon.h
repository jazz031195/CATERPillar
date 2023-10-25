//!  Axon Obstacle Derived Class =============================================================/
/*!
*   \details   Axon class derived from an Obstacle
*              in the direction set by begin, end.
*   \author    Jasmine Nguyen-Duc
*   \date      September 2023
*   \version   1.42
=================================================================================================*/

#ifndef AXON_H
#define AXON_H

#include "sphere.h"
#include "obstacle.h"
#include <vector>

using namespace std;

/// @brief
class Axon : public Obstacle
{
public:
    int id;                                         /*!< ID of axon */
    std::vector<Sphere> spheres;            /*!< spheres in axon */
    double radius;                                  /*!< radius of axon */
    Eigen::Vector3d begin;                          /*!< position of first sphere */
    Eigen::Vector3d end;                            /*!< target position to grow towards */
    std::vector<Eigen::Vector2d> Box;               /*!< Box with <min, max> for each axis (x,y,z) */
    int growth_attempts;

    /*!
     *  \brief Default constructor. Does nothing
     */
    Axon(){};

    ~Axon(){
        spheres.clear();
        Box.clear();
    };


    Axon(int id_, Eigen::Vector3d begin_, Eigen::Vector3d end_, double radius_)
    {
        id = id_;
        begin = begin_;
        end = end_;
        radius = radius_;
        spheres.clear();
        Box.clear();
        growth_attempts = 0;
    }

    Axon(Axon const &ax);

    void keep_one_sphere(Eigen::Vector3d begin_, Eigen::Vector3d end_);
    /*!
     *  \param sphere_to_add sphere to add
     *  \brief Adds sphere to axon
     */
    void add_sphere(Sphere sphere_to_add);

    /*!
     *  \param sph sphere to check
     *  \brief Checks if a sphere collides with this axon
     */
    bool isSphereInsideAxon(Sphere sph);

    /*!
     *  \param sph sphere to check
        \param axis axis (0, 1 or 2)
     *  \brief Checks if a sphere collides with this axon along a specific axis
     */
    std::vector<int> checkAxisForCollision(Sphere sph, int axis);

    /*!
     *  \param position position in voxel
        \param distance_to_be_inside distance to be inside Box
     *  \brief Checks if a position is near this axon (near the Box)
     */
    bool isNearAxon(Eigen::Vector3d position, double distance_to_be_inside);

    /*!
     *  \brief Deletes all spheres in axon except the first one. The Box is updated.
     */
    void keep_only_first_sphere();

    /*!
     *  \brief Deletes all spheres in axon. The Box is updated.
     */
    void destroy();

};

#endif // AXON_H
