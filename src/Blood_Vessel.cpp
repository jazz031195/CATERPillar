#include "Blood_Vessel.h"
#include "constants.h"
#include "Eigen/Dense"
#include <iostream>


using namespace Eigen;
using namespace std;


Blood_Vessel::Blood_Vessel(){}

Blood_Vessel::~Blood_Vessel()
{
    spheres.clear();
}


void Blood_Vessel::keep_one_sphere(){
    Sphere s (this->spheres[0]);
    spheres.clear();
    add_sphere(s);
    growth_attempts += 1;


}

void Blood_Vessel::destroy(){
    spheres.clear();
    growth_attempts = 0;

}

void Blood_Vessel::addToGrid(SphereGrid &grid) const{

    for (const auto &sph : spheres){
        grid.insert(sph);
    }
}

void Blood_Vessel::add_first_sphere(const Sphere &s){

    this->spheres.clear();
    this->spheres.push_back(s);

}


void Blood_Vessel::add_sphere(const Sphere &sphere_to_add){
    // add sphere to list of spheres
    this->spheres.push_back(sphere_to_add);
}


void Blood_Vessel::update_Volume(const int &factor, const Eigen::Vector3d &min_limits, const Eigen::Vector3d &max_limits){
    double new_volume = 0.0;

    for (size_t i = 1; i < spheres.size(); ++i) {
        const Sphere &last_sphere = spheres[i - 1];
        const Sphere &current_sphere = spheres[i];

        double distance = (current_sphere.center - last_sphere.center).norm();

        // Volume of truncated cone
        double segment_volume = M_PI * distance * (
            current_sphere.radius * current_sphere.radius +
            last_sphere.radius * last_sphere.radius +
            current_sphere.radius * last_sphere.radius) / 3.0;

        bool current_in_bounds = Obstacle::check_borders(min_limits, max_limits, current_sphere.center, barrier_tickness);
        bool last_in_bounds = Obstacle::check_borders(min_limits, max_limits, last_sphere.center, barrier_tickness);

        if (current_in_bounds && last_in_bounds) {
            new_volume += segment_volume;
        } else if (current_in_bounds || last_in_bounds) {
            new_volume += segment_volume / 2.0;
        }
        // If both are out of bounds, add nothing
    }

    volume = new_volume;
}
