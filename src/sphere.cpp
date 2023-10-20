#include "sphere.h"
#include "constants.h"
#include <Eigen/Dense>
#include <iostream>


using namespace Eigen;
using namespace std;



Sphere::Sphere(const Sphere &sph)
{
    center = sph.center;
    radius = sph.radius;
    ax_id = sph.ax_id;
    id = sph.id;

}


double Sphere::minDistance(Eigen::Vector3d O){

    Vector3d m = O - this->center;
    // minimum distance to the cylinder axis.
    double distance_to_sphere = m.norm();

    //Minimum distance to the shpere wall.
    double d_ = (distance_to_sphere - this->radius);
   // return d_>0.0?d_:0.0;
    return d_;
}


bool Sphere::isInside(Eigen::Vector3d pos, double distance_to_be_inside){
    double d_ = (pos - this->center).norm();
 
    d_ = d_-this->radius;
    
   // return d_>0.0?d_:0.0;
    return d_ <= distance_to_be_inside;
}


bool Sphere::CollideswithSphere(Sphere sphere_, double distance_to_be_inside){
    double d_ = (sphere_.center - this->center).norm();
    if (d_ <= sphere_.radius + this->radius + distance_to_be_inside){
        return true;
    }
    else{
        return false;
    }
}
