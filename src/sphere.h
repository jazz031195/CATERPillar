//!  Sphere Obstacle Derived Class =============================================================/
/*!
*   \details   Sphere obstacle class derived from an Obstacle. 
*   \author    Jasmine Nguyen-Duc
*   \date      September 2023
*   \version   1.5
=================================================================================================*/

#ifndef DYN_SPHERE_H
#define DYN_SPHERE_H

#include "obstacle.h"

class Sphere : public Obstacle
{
public:

    int id;                    /*!< ID of sphere     */
    Eigen::Vector3d center;    /*!< Position of center of the sphere      */
    double radius;             /*!< Radius of the Sphere                  */
    int ax_id;                 /*!< ID of axon to which the sphere belongs to     */


    /*!
     *  \brief Default constructor. Does nothing
     */
    Sphere(){}
    /*!
     *  \brief Default destructor. Does nothing
     */
    ~Sphere(){}

    /*!
     *  \param center_ Sphere origin
     *  \param radius_ Sphere's radius
     *  \param id_ ID of sphere
     *  \param ax_id_ ID of axon to which the sphere belongs to 
     *  \brief Initialize everything.
     */
    Sphere(int id_, int ax_id_, Eigen::Vector3d center_, double radius_){

        id =id_;
        ax_id = ax_id_;
        center = center_;
        radius = radius_;
    }

    /*!
     *  \brief constrcutor by copy
     */
    Sphere(Sphere const &sph);

    /*!
     *  \param pos Position to check
     *  \param distance_to_be_inside Distance to be inside sphere
     *  \brief Returns true if position is inside sphere
     */
    bool isInside(Eigen::Vector3d pos, double distance_to_be_inside);

    /*!
     *  \param O Position to check
     *  \brief Distance from sphere to position
     */
    double minDistance(Eigen::Vector3d O);

    /*!
     *  \param sphere_ sphere to check
     *  \param distance_to_be_inside Distance to be inside sphere
     *  \brief checks if this sphere collides with another sohere
     */
    bool CollideswithSphere(Sphere sphere_, double distance_to_be_inside);
 

};  

#endif // Sphere_H
