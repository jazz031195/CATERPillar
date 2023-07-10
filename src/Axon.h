//!  Cylinder Obstacle Derived Class =============================================================/
/*!
*   \details   Cylinder class derived from an Obstacle. Defines infinite long cylinders
*              in the direction set by P,Q.
*   \author    Jasmine Nguyen-Duc
*   \date      September 2022 
*   \version   1.42
=================================================================================================*/


#ifndef AXON_H
#define AXON_H

#include "dynamic_sphere.h"
#include "swipeprune.h"
#include "obstacle.h"
#include <vector>

using namespace std;

/// @brief 
class Axon : public Obstacle
{
public:
    int id;
    std::vector<Dynamic_Sphere> spheres; 
    double radius;
    Eigen::Vector3d begin;
    Eigen::Vector3d end;
    // Box with <min, max> for each axis (x,y,z)
    std::vector<Eigen::Vector2d> Box;
    Projections projections;


    /*!
     *  \brief Default constructor. Does nothing
     */
    Axon(){};

    ~Axon(){};


    Axon(int id_,  Eigen::Vector3d begin_,Eigen::Vector3d end_ ,  double radius_){

        id = id_;
        begin = begin_;
        end = end_;
        radius =radius_;
        spheres.clear();

    }
    Axon(Axon const &ax);

    void add_sphere(Dynamic_Sphere sphere_to_add);
    bool isNearAxon(Eigen::Vector3d position,  double distance_to_be_inside);
    bool isPosInsideAxon(Eigen::Vector3d &position,  double distance_to_be_inside, std::vector<int>& sphere_ids,double max_radius );
    void add_projection(Dynamic_Sphere sphere_to_add);
    void threadGrowth();
};


#endif // AXON_H
