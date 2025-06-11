#include "sphere.h"
#include "constants.h"
#include <iostream>
#include <random>

using namespace Eigen;
using namespace std;


double Sphere::minDistance(const Eigen::Vector3d &O) const{

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


bool Sphere::CollideswithSphere(const Sphere &sphere_,const double &distance_to_be_inside) const{
    double d_ = (sphere_.center - this->center).norm();
    if (d_ <= sphere_.radius + this->radius + distance_to_be_inside){
        return true;
    }
    else{
        return false;
    }
}


void Sphere::getPointOnSphereSurface(Eigen::Vector3d &point, Eigen::Vector3d &vector, const Eigen::Vector3d &vector2, const bool& primary_process) const {
    // random device
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<> disTheta(0, 2 * M_PI);  // Uniform distribution for azimuthal angle
    std::uniform_real_distribution<> disCosPhi(-1, 1);       // Uniform distribution for cosine of polar angle

    bool valid = false;
    int tries = 0;
    while (!valid && tries < 1000) {
        // Generate random spherical coordinates
        double theta = disTheta(gen);     // Azimuthal angle (0 to 2 * pi)
        double cosPhi = disCosPhi(gen);   // Cosine of polar angle (-1 to 1)
        double sinPhi = sqrt(1 - cosPhi * cosPhi); // Sin of the polar angle

        // Convert spherical coordinates to Cartesian coordinates (uniformly on sphere)
        vector = Eigen::Vector3d(sinPhi * cos(theta), sinPhi * sin(theta), cosPhi);
        vector.normalize();  // Ensure vector is unit length
        point = this->center + (this->radius * vector);  // Scale by the radius to get point on surface

        if (!primary_process) {
            double angle = acos(vector.dot(vector2));
            // Check if the angle is <= pi / 2
            if (angle <= M_PI / 4.0) {
                valid = true;
            }
            else{
                tries++;
                valid = false;
            }
        }
        else{
            valid = true;
        }
    }
}




Eigen::Vector3d Sphere::getNormalOnSpherePoint(Eigen::Vector3d point){
    return (point - this->center).normalized();
}
