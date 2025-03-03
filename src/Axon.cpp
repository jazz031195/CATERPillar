#include "Axon.h"
#include "constants.h"
#include "Eigen/Dense"
#include <iostream>


using namespace Eigen;
using namespace std;


Axon::Axon()
{}

Axon::~Axon()
{
    inner_spheres.clear();
    outer_spheres.clear();
    Box.clear();
}

void Axon::ModifyRadiusFirstSphere(const double &new_radius){
    outer_spheres[0].radius = new_radius;

}

void Axon::keep_one_sphere(){
    Sphere s (this->outer_spheres[0]);
    outer_spheres.clear();
    Box.clear();
    add_sphere(s);
    growth_attempts += 1;

}

void Axon::destroy(){
    outer_spheres.clear();
    Box.clear();
    growth_attempts = 0;
}

void Axon::updateBox(){
    Box.clear();
    for (int i = 0; i < outer_spheres.size(); i++){
        Sphere s = outer_spheres[i];
        double sph_highest_x_val = s.center[0]+ s.radius;
        double sph_lowest_x_val = s.center[0] -s.radius;
        double sph_highest_y_val = s.center[1] +s.radius;
        double sph_lowest_y_val = s.center[1] -s.radius;
        double sph_highest_z_val = s.center[2] +s.radius;
        double sph_lowest_z_val = s.center[2] -s.radius;

        if (i == 0){
            Box.push_back({sph_lowest_x_val, sph_highest_x_val});
            Box.push_back({sph_lowest_y_val, sph_highest_y_val});
            Box.push_back({sph_lowest_z_val, sph_highest_z_val});
        }
        else{
            if (Box[0][0] > sph_lowest_x_val){
                Box[0][0] = sph_lowest_x_val;
            }
            if (Box[0][1] < sph_highest_x_val){
                Box[0][1] = sph_highest_x_val;
            }
            if (Box[1][0] > sph_lowest_y_val){
                Box[1][0] = sph_lowest_y_val;
            }
            if (Box[1][1] < sph_highest_y_val){
                Box[1][1] = sph_highest_y_val;
            }
            if (Box[2][0] > sph_lowest_z_val){
                Box[2][0] = sph_lowest_z_val;
            }
            if (Box[2][1] < sph_highest_z_val){
                Box[2][1] = sph_highest_z_val;
            }
        }
    }
} 
void Axon::add_sphere(const Sphere &sphere_to_add){

    // value of center of sphere at x that has the highest x center value
    double sph_highest_x_val;
    // value of center of sphere at y that has the highest y center value
    double sph_highest_y_val;
    // value of center of sphere at x that has the lowest x center value
    double sph_lowest_x_val;
    // value of center of sphere at y has the lowest y center value
    double sph_lowest_y_val;
    // value of center of sphere at z that has the highest z center value
    double sph_highest_z_val;
    // value of center of sphere at z that has the lowest z center value
    double sph_lowest_z_val;

    // add sphere to list of spheres
    this->outer_spheres.push_back(sphere_to_add);


    // if there is only one sphere in list
    if (this->outer_spheres.size() == 1){
        // create box around that one sphere
        sph_highest_x_val = sphere_to_add.center[0]+ sphere_to_add.radius;
        sph_lowest_x_val = sphere_to_add.center[0] -sphere_to_add.radius;
        sph_highest_y_val = sphere_to_add.center[1] +sphere_to_add.radius;
        sph_lowest_y_val = sphere_to_add.center[1] -sphere_to_add.radius;
        sph_highest_z_val = sphere_to_add.center[2] +sphere_to_add.radius;
        sph_lowest_z_val = sphere_to_add.center[2] -sphere_to_add.radius;

        Box.clear();
        //x
        Box.push_back({sph_lowest_x_val, sph_highest_x_val});
        //y
        Box.push_back({sph_lowest_y_val, sph_highest_y_val});
        //z
        Box.push_back({sph_lowest_z_val, sph_highest_z_val});

        this->begin = sphere_to_add.center;

    }
    else{

        // take values in the box
        sph_highest_x_val = Box[0][1];
        sph_lowest_x_val = Box[0][0];
        sph_highest_y_val = Box[1][1];
        sph_lowest_y_val = Box[1][0];
        sph_highest_z_val = Box[2][1];
        sph_lowest_z_val = Box[2][0];

        //cout << "---------" << endl;
        //cout << "Box: " << Box[0][0] << " " << Box[0][1] << " " << Box[1][0] << " " << Box[1][1] << " " << Box[2][0] << " " << Box[2][1] << endl;

        // if the extemity of sphere is higher than extermity of Box (x)
        if (sph_highest_x_val < sphere_to_add.center[0]+ sphere_to_add.radius){
            Box[0][1] = sphere_to_add.center[0]+ sphere_to_add.radius;
        }
        // if the extemity of sphere is lower than extermity of Box (x)
        if(sph_lowest_x_val > sphere_to_add.center[0]- sphere_to_add.radius){
            Box[0][0] = sphere_to_add.center[0]- sphere_to_add.radius;
        }
        // if the extemity of sphere is higher than extermity of Box (y)
        if (sph_highest_y_val < sphere_to_add.center[1]+ sphere_to_add.radius){
            Box[1][1] = sphere_to_add.center[1]+ sphere_to_add.radius;
        }
        // if the extemity of sphere is lower than extermity of Box (y)
        if(sph_lowest_y_val > sphere_to_add.center[1]- sphere_to_add.radius){
            Box[1][0]= sphere_to_add.center[1]- sphere_to_add.radius;
        }
        // if the extemity of sphere is higher than extermity of Box (z)
        if (sph_highest_z_val < sphere_to_add.center[2]+ sphere_to_add.radius){
            Box[2][1] = sphere_to_add.center[2]+ sphere_to_add.radius;
        }
        // if the extemity of sphere is lower than extermity of Box (z)
        if(sph_lowest_z_val > sphere_to_add.center[2]- sphere_to_add.radius){
            Box[2][0] = sphere_to_add.center[2]- sphere_to_add.radius;
        }
    //cout << "Box: " << Box[0][0] << " " << Box[0][1] << " " << Box[1][0] << " " << Box[1][1] << " " << Box[2][0] << " " << Box[2][1] << endl;
    }
}

bool check_position_is_in_Box(Eigen::Vector3d position, std::vector<Eigen::Vector2d> Box, double buffer_distance){
    if ((position[0] >=  Box[0][0]- buffer_distance)  && (position[0] <= Box[0][1] + buffer_distance)){
        if ((position[1] >= Box[1][0] - buffer_distance) && (position[1] <=  Box[1][1] + buffer_distance)){
            return true;
        }
    }
    return false;
}
 

bool Axon::isNearAxon(const Eigen::Vector3d &position, const double &distance_to_be_inside) const{

    if (outer_spheres.size()==0){
        return false;
    }
    else{
        if (check_position_is_in_Box(position, Box, distance_to_be_inside)){
            return true;
        }
    }
    return false;
}

std::vector<int> findCommonIntegers(const std::vector<int>& vec1, const std::vector<int>& vec2, const std::vector<int>& vec3) {
    std::vector<int> result;
    std::set_intersection(vec1.begin(), vec1.end(),
                          vec2.begin(), vec2.end(),
                          std::back_inserter(result));
    std::vector<int> commonIntegers;
    std::set_intersection(result.begin(), result.end(),
                          vec3.begin(), vec3.end(),
                          std::back_inserter(commonIntegers));
    return commonIntegers;
}
std::vector<int> Axon::checkAxisForCollision(const Sphere &sph, const int &axis) const{

	std::vector<int> spheres_id_to_check;
	for (auto i = 0; i < outer_spheres.size(); ++i) {

            double min_i = outer_spheres[i].center[axis] - outer_spheres[i].radius - barrier_tickness;
            double max_j = sph.center[axis] + sph.radius + barrier_tickness;

			if (min_i> max_j) {
				continue;
			}
			else {
				double max_i = outer_spheres[i].center[axis] + outer_spheres[i].radius;
                double min_j = sph.center[axis] - sph.radius;
                if (min_j> max_i) {
                    continue;
                }
                else{
                    spheres_id_to_check.push_back(i);
                }
			}
	}

    return spheres_id_to_check;
}

std::vector<int> Axon::checkAxisForInnerCollision(const Sphere &sph, const int &axis) const{

    if (inner_spheres.size() == 0){
        return {};
    }

	std::vector<int> spheres_id_to_check;
	for (auto i = 0; i < inner_spheres.size(); ++i) {

            double min_i = inner_spheres[i].center[axis] - inner_spheres[i].radius - barrier_tickness;
            double max_j = sph.center[axis] + sph.radius + barrier_tickness;

			if (min_i> max_j) {
				continue;
			}
			else {
				double max_i = inner_spheres[i].center[axis] + inner_spheres[i].radius;
                double min_j = sph.center[axis] - sph.radius;
                if (min_j> max_i) {
                    continue;
                }
                else{
                    spheres_id_to_check.push_back(i);
                }
			}
	}

    return spheres_id_to_check;
}

bool Axon::isSphereInsideAxon(const Sphere &sph) const{

    if (isNearAxon(sph.center, 2*sph.radius + barrier_tickness)){ // if near axon
        if(!(sph.object_type == 0 && sph.object_id == id)){ 
            std::vector<std::vector<int>> spheres_id_to_check;
            for (auto axis = 0; axis < 3; ++axis) {
                spheres_id_to_check.push_back(checkAxisForCollision(sph, axis)); // check for collision along 1 axis
                if (spheres_id_to_check[axis].size() == 0){
                    return false;
                }
            }
            // find common ids in all 3 axes
            std::vector<int> spheres_to_check_all_axes = findCommonIntegers(spheres_id_to_check[0], spheres_id_to_check[1], spheres_id_to_check[2]);
            for (auto i = 0; i < spheres_to_check_all_axes.size(); ++i) {
                Sphere sphere_to_check = outer_spheres[spheres_to_check_all_axes[i]];
                if (sph.minDistance(sphere_to_check.center) < sphere_to_check.radius){
                    return true;
                }
            }
            spheres_id_to_check.clear();
            spheres_to_check_all_axes.clear();
        }
    }

    return false;
}

std::vector<int> Axon::isSphereInsideAxon_(const Sphere &sph) {

    std::vector<int> spheres_coliding;
    if (isNearAxon(sph.center, 2*sph.radius + barrier_tickness)){ // if near axon
        if(!(sph.object_type == 0 && sph.object_id == id)){ 
            std::vector<std::vector<int>> spheres_id_to_check;
            for (auto axis = 0; axis < 3; ++axis) {
                spheres_id_to_check.push_back(checkAxisForCollision(sph, axis)); // check for collision along 1 axis
                //if (spheres_id_to_check[axis].size() == 0){
                //    return spheres_coliding;
                //}
            }
            // find common ids in all 3 axes
            std::vector<int> spheres_to_check_all_axes = findCommonIntegers(spheres_id_to_check[0], spheres_id_to_check[1], spheres_id_to_check[2]);
            for (auto i = 0; i < spheres_to_check_all_axes.size(); ++i) {
                Sphere sphere_to_check = outer_spheres[spheres_to_check_all_axes[i]];
                if (sph.minDistance(sphere_to_check.center) < sphere_to_check.radius){
                    spheres_coliding.push_back(spheres_to_check_all_axes[i]);
                }
            }
            spheres_id_to_check.clear();
            spheres_to_check_all_axes.clear();
        }
    }

    return spheres_coliding;
}

bool Axon::isSphereInsideInnerAxon(const Sphere &sph) const{

    if (isNearAxon(sph.center, 2*sph.radius + barrier_tickness)){ // if near axon
        if(!(sph.object_type == 0 && sph.object_id == id)){ 
            std::vector<std::vector<int>> spheres_id_to_check;
            for (auto axis = 0; axis < 3; ++axis) {
                spheres_id_to_check.push_back(checkAxisForInnerCollision(sph, axis)); // check for collision along 1 axis
                if (spheres_id_to_check[axis].size() == 0){
                    return false;
                }
            }
            // find common ids in all 3 axes
            std::vector<int> spheres_to_check_all_axes = findCommonIntegers(spheres_id_to_check[0], spheres_id_to_check[1], spheres_id_to_check[2]);
            for (auto i = 0; i < spheres_to_check_all_axes.size(); ++i) {
                Sphere sphere_to_check = inner_spheres[spheres_to_check_all_axes[i]];
                if (sph.minDistance(sphere_to_check.center) < sphere_to_check.radius){
                    return true;
                }
            }
            spheres_id_to_check.clear();
            spheres_to_check_all_axes.clear();
        }
    }

    return false;
}

// Function to check if a point is inside a dilated box
bool check_borders(const Eigen::Vector3d&  min_l, const Eigen::Vector3d&  max_l, const Eigen::Vector3d& pos, const double& distance_to_border) {

    
    // Check if the point is inside the dilated box
    for (int i = 0; i < 3; ++i) {
        double min_bound = min_l[0] - distance_to_border;
        double max_bound = max_l[1] + distance_to_border;
        if (pos[i] < min_bound || pos[i] > max_bound) {
            return false; // Point is outside the dilated box
        }
    }
    
    return true; // Point is inside the dilated box
}


void Axon::update_Volume(const int &factor, const Eigen::Vector3d &min_limits, const Eigen::Vector3d &max_limits){

    double new_volume = 0.0;
    for (auto i = factor; i < outer_spheres.size(); i+= factor) {
        Sphere current_sphere = outer_spheres[i];
        Sphere last_sphere = outer_spheres[i-factor];

        bool current_sphere_within_boundaries = check_borders(min_limits, max_limits, current_sphere.center, barrier_tickness);
        bool last_sphere_within_boundaries = check_borders(min_limits, max_limits, last_sphere.center, barrier_tickness);
        double distance = (current_sphere.center - last_sphere.center).norm();
        double volume_cone = M_PI * distance *
                            (current_sphere.radius * current_sphere.radius +
                             last_sphere.radius * last_sphere.radius +
                             current_sphere.radius * last_sphere.radius) / 3;
        if (current_sphere_within_boundaries && last_sphere_within_boundaries){
            new_volume += volume_cone;
        }
        else if (current_sphere_within_boundaries || last_sphere_within_boundaries){ 
            new_volume += volume_cone/2;
        } 
    }
    volume = new_volume;

} 

