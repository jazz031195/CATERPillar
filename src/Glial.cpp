#include "Glial.h"
#include "Eigen/Dense"
#include "constants.h"
#include <iostream>

using namespace Eigen;
using namespace std;

Glial::Glial()
{}

Glial::~Glial()
{}

void Glial::initialise_boxes(){


    big_box = Box {0,0,0,0,0,0};
    soma_box = Box {0,0,0,0,0,0};

    Sphere sph = this->soma;
    double sph_highest_x_val = sph.center[0] + sph.radius;
    double sph_lowest_x_val = sph.center[0] - sph.radius;
    double sph_highest_y_val = sph.center[1] + sph.radius;
    double sph_lowest_y_val = sph.center[1] - sph.radius;
    double sph_highest_z_val = sph.center[2] + sph.radius;
    double sph_lowest_z_val = sph.center[2] - sph.radius;

    soma_box.x_max = sph_highest_x_val;
    soma_box.x_min = sph_lowest_x_val;
    soma_box.y_max = sph_highest_y_val;
    soma_box.y_min = sph_lowest_y_val;
    soma_box.z_max = sph_highest_z_val;
    soma_box.z_min= sph_lowest_z_val;

    big_box = soma_box;
}

void Glial::update_Box(Box& box, const Sphere &sph){

    double sph_highest_x_val = sph.center[0] + sph.radius;
    double sph_lowest_x_val = sph.center[0] - sph.radius;
    double sph_highest_y_val = sph.center[1] + sph.radius;
    double sph_lowest_y_val = sph.center[1] - sph.radius;
    double sph_highest_z_val = sph.center[2] + sph.radius;
    double sph_lowest_z_val = sph.center[2] - sph.radius;

    if (sph_highest_x_val > box.x_max){
        box.x_max = sph_highest_x_val;
    }
    if (sph_lowest_x_val < box.x_min){
        box.x_min = sph_lowest_x_val;
    }
    if (sph_highest_y_val > box.y_max){
        box.y_max = sph_highest_y_val;
    }
    if (sph_lowest_y_val < box.y_min){
        box.y_min = sph_lowest_y_val;
    }
    if (sph_highest_z_val > box.z_max){
        box.z_max = sph_highest_z_val;
    }
    if (sph_lowest_z_val < box.z_min){
        box.z_min = sph_lowest_z_val;
    }

}
void Glial::update_all_boxes(const Sphere &sph){

    int branch_id = sph.branch_id;
    if (ramification_boxes.size() <= branch_id) {
        ramification_boxes.resize(branch_id + 1);
        Box new_box = {0,0,0,0,0,0};
        ramification_boxes[branch_id] = new_box;
    }
    update_Box(ramification_boxes[branch_id], sph);
    update_Box(big_box, sph);
    
}

bool Glial::isNearBox(const Box& box, const Eigen::Vector3d &position, const double &distance_to_be_inside) const {

    if ((position[0] >= box.x_min - distance_to_be_inside) && (position[0] <= box.x_max + distance_to_be_inside)) {
        if ((position[1] >= box.y_min - distance_to_be_inside) && (position[1] <= box.y_max + distance_to_be_inside)) {
            if ((position[2] >= box.z_min - distance_to_be_inside) && (position[2] <= box.z_max + distance_to_be_inside)) {
                return true; // point is near the glial cell
            }
        }
    }
    return false;
}
bool Glial::isNearGlialCell(const Eigen::Vector3d &position, const double &distance_to_be_inside) const{
    
    if (isNearBox(big_box, position, distance_to_be_inside)){
        return true;
    }
    else{
        return false;
    }

}

bool Glial::collidesWithRamification(const Eigen::Vector3d &position, const double &distance_to_be_inside) const{

    if (ramification_boxes.empty()){
        return false;
    }

    std::vector<int> indices_close_boxes;

    int i = 0;
    for (const auto& box : ramification_boxes) {
        if (isNearBox(box, position, distance_to_be_inside)) {
            indices_close_boxes.push_back(i);
        }
        i += 1;
    }

    if (indices_close_boxes.empty()){
        return false; // no ramification close to the point
    }

    for (const auto& index : indices_close_boxes) {
        const std::vector<Sphere>& branch = ramification_spheres[index];
        for (const auto& sphere : branch) {
            if ((sphere.center - position).norm() <= sphere.radius + distance_to_be_inside) {
                return true; // collision with ramification
            }
        }
    }
    return false; // no collision with any ramification

}

bool Glial::collidesWithItsOwnRamification(const Sphere &sph, const double &distance_to_be_inside) const{


    std::vector<int> indices_close_boxes;
    Eigen::Vector3d position = sph.center;

    int i = 0;
    for (const auto& box : ramification_boxes) {
        if (isNearBox(box, position, distance_to_be_inside)) {
            if(i == sph.branch_id){
                continue; // do not check the branch the sphere belongs to
            }
            indices_close_boxes.push_back(i);
        }
        i += 1;
    }

    if (indices_close_boxes.empty()){
        return false; // no ramification close to the point
    }

    for (const auto& index : indices_close_boxes) {
        const std::vector<Sphere>& branch = ramification_spheres[index];
        for (const auto& sphere : branch) {
            if ((sphere.center - position).norm() <= sphere.radius + sph.radius +  distance_to_be_inside) {
                return true; // collision with ramification
            }
        }
    }
    return false; // no collision with any ramification

}
bool Glial::collides_with_GlialCell(const Sphere &sph, const double &distance_to_be_inside) const{


    if (sph.object_type == 1 && sph.object_id == this->id){
        return false; // do not collide with itself
    }  
    else if (this->soma.CollideswithSphere(sph, distance_to_be_inside)) // check if sphere collides with soma
    {
        return true; // overlap
    }
    
    else if (!isNearGlialCell(sph.center, sph.radius + distance_to_be_inside)) // check if sphere is near glial cell
    {
        return false; 
    } 
    else if (collidesWithRamification(sph.center, sph.radius + distance_to_be_inside)) // overlap
    {
        return true;
    } 
    else{
        return false; // no overlap
    } 

}

void Glial::compute_processes_icvf(const int &factor, const Eigen::Vector3d &min_limits, const Eigen::Vector3d &max_limits) {
    volume_processes = 0.0;
    for (int b = 0; b < ramification_spheres.size(); ++b) {
        for (int i = factor; i < ramification_spheres[b].size(); i += factor) {
            if ((soma.center-ramification_spheres[b][i-factor].center).norm() < soma.radius) {
                continue;
            }
            Sphere &s = ramification_spheres[b][i-factor];
            Sphere &s_next = ramification_spheres[b][i];
            if (s_next.center [0] + s_next.radius < min_limits[0] || s_next.center [0] - s_next.radius > max_limits[0] ||
                s_next.center [1] + s_next.radius < min_limits[1] || s_next.center [1] - s_next.radius > max_limits[1] ||
                s_next.center [2] + s_next.radius < min_limits[2] || s_next.center [2] - s_next.radius > max_limits[2]) {
                continue;
            }
            
            double distance = (s_next.center - s.center).norm();
            double v = M_PI * (s.radius * s.radius + s_next.radius * s_next.radius + s.radius * s_next.radius) * distance / 3.0;
            volume_processes += v;

        }
    }
}