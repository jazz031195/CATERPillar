#include "dynamic_sphere.h"
#include "constants.h"
#include <Eigen/Dense>
#include <iostream>


using namespace Eigen;
using namespace std;



Dynamic_Sphere::Dynamic_Sphere(const Dynamic_Sphere &sph)
{
    center = sph.center;
    radius = sph.radius;
    ax_id = sph.ax_id;
    id = sph.id;

}



double Dynamic_Sphere::minDistance(Eigen::Vector3d O){

    Vector3d m = O - this->center;
    // minimum distance to the cylinder axis.
    double distance_to_sphere = m.norm();

    //Minimum distance to the shpere wall.
    double d_ = (distance_to_sphere - this->radius);
   // return d_>0.0?d_:0.0;
    return d_;
}


bool Dynamic_Sphere::isInside(Eigen::Vector3d pos, double distance_to_be_inside){
    double d_ = (pos - this->center).norm();
 
    d_ = d_-this->radius;
    
   // return d_>0.0?d_:0.0;
    return d_ <= distance_to_be_inside;
}

bool Dynamic_Sphere::distSmallerThan(Eigen::Vector3d pos, double distance){
    double d_ = (pos - this->center).norm();
    
    return d_ <= distance;
}

bool Dynamic_Sphere::CollideswithSphere(Dynamic_Sphere sphere_, double distance_to_be_inside){
    double d_ = (sphere_.center - this->center).norm();
    if (d_ <= sphere_.radius + this->radius + distance_to_be_inside){
        return true;
    }
    else{
        return false;
    }
}
