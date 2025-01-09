//!  Sphere Obstacle Derived Class =============================================================/
/*!
*   \details   Sphere obstacle class derived from an Obstacle. 
*   \author    Jasmine Nguyen-Duc
*   \date      September 2023
*   \version   1.5
=================================================================================================*/

#ifndef DYN_SPHERE_H
#define DYN_SPHERE_H

#include "Eigen/Core"
#include "obstacle.h"
#include <vector>

class Sphere : public Obstacle
{
public:

    int id;                    /*!< ID of sphere     */
    Eigen::Vector3d center;    /*!< Position of center of the sphere      */
    double radius;             /*!< Radius of the Sphere                  */
    int object_id;              /*!< ID of object to which the sphere belongs to     */
    int object_type;         /*!< Type of object to which the sphere belongs to (0 : outer axon, 1:glial, 2: inner axon)    */
    int branch_id;           /*!< ID of branch of glial cell to which the sphere belongs to     */
    int parent_id;           /*!< ID of parent sphere of the sphere     */
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
    Sphere(int id_, int object_id_, int object_type_, Eigen::Vector3d center_, double radius_, int branch_id_ = -1, int parent_id_ = -1){

        id =id_;
        object_id = object_id_;
        object_type = object_type_;
        center = center_;
        radius = radius_;
        branch_id = branch_id_;
        parent_id = parent_id_;
    }


    /*!
     *  \brief constrcutor by copy
     */

    Sphere& operator=(const Sphere& sph) {
        if (this != &sph) {
            center = sph.center;
            radius = sph.radius;
            object_id = sph.object_id;
            object_type = sph.object_type;
            id = sph.id;
            branch_id = sph.branch_id;
            parent_id = sph.parent_id;
        }
        return *this;
    }

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
    double minDistance(const Eigen::Vector3d &O) const;

    /*!
     *  \param sphere_ sphere to check
     *  \param distance_to_be_inside Distance to be inside sphere
     *  \brief checks if this sphere collides with another sohere
     */
    bool CollideswithSphere(const Sphere &sphere_,const double &distance_to_be_inside) const;
    
    void getPointOnSphereSurface(Eigen::Vector3d &point, Eigen::Vector3d &vector, const Eigen::Vector3d &vector2) const;

    Eigen::Vector3d getNormalOnSpherePoint(Eigen::Vector3d point);


};  

#endif // Sphere_H
