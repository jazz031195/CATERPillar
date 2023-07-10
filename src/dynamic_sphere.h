//!  Sphere Obstacle Derived Class =============================================================/
/*!
*   \details   Sphere obstacle class derived from an Obstacle. Defines a analyitical sphere of radius
*   r and centered in center
*   \author    Jonathan Rafael
*   \date      2020
*   \version   1.5
=================================================================================================*/

#ifndef DYN_SPHERE_H
#define DYN_SPHERE_H

#include "obstacle.h"

class Dynamic_Sphere : public Obstacle
{
public:

    int id;
    Eigen::Vector3d center;    /*!< Cilinder Axis reference Points, P should be the "center"      */
    double radius;             /*!< Radius of the Sphere                                          */
    int ax_id;


    /*!
     *  \brief Default constructor. Does nothing
     */
    Dynamic_Sphere(){}
    /*!
     *  \brief Default destructor. Does nothing
     */
    ~Dynamic_Sphere(){}

    /*!
     *  \param center Sphere origin
     *  \param radius Sphere's radius
     *  \param scale  overall scale for when reading files.
     *  \brief Initialize everything.
     */
    Dynamic_Sphere(int id_, int ax_id_, Eigen::Vector3d center_, double radius_){

        id =id_;
        ax_id = ax_id_;
        center = center_;
        radius = radius_;
    }

    /*!
     *  \brief constrcutor by copy
     */
    Dynamic_Sphere(Dynamic_Sphere const &sph);
    bool isInside(Eigen::Vector3d pos, double distance_to_be_inside);
    bool distSmallerThan(Eigen::Vector3d pos, double distance);
    double minDistance(Eigen::Vector3d O);
    bool CollideswithSphere(Dynamic_Sphere sphere_, double distance_to_be_inside);
 

};  

#endif // Sphere_H
