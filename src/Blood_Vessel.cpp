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
    Box.clear();
}


void Blood_Vessel::keep_one_sphere(){
    Sphere s (this->spheres[0]);
    spheres.clear();
    Box.clear();
    add_sphere(s);
    growth_attempts += 1;


}

void Blood_Vessel::destroy(){
    spheres.clear();
    Box.clear();
    growth_attempts = 0;

}

void Blood_Vessel::add_first_sphere(const Sphere &s){

    this->spheres.clear();
    this->spheres.push_back(s);
    updateBox();
}

void Blood_Vessel::updateBox(){
    Box.clear();
    for (int i = 0; i < spheres.size(); i++){
        Sphere s = spheres[i];
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

void Blood_Vessel::add_sphere(const Sphere &sphere_to_add){

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
    this->spheres.push_back(sphere_to_add);


    // if there is only one sphere in list
    if (this->spheres.size() == 1){
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


bool Blood_Vessel::isNearBlood_Vessel(const Eigen::Vector3d &position, const double &distance_to_be_inside) const{

    if (spheres.size()==0){
        return false;
    }
    else{
        if (Obstacle::check_position_is_in_Box(position, Box, distance_to_be_inside)){
            return true;
        }
    }
    return false;
}




std::vector<int> Blood_Vessel::checkAxisForCollision(const Sphere &sph, const int &axis) const{

	std::vector<int> spheres_id_to_check;
	for (auto i = 0; i < spheres.size(); ++i) {

            double min_i = spheres[i].center[axis] - spheres[i].radius - barrier_tickness;
            double max_j = sph.center[axis] + sph.radius + barrier_tickness;

			if (min_i> max_j) {
				continue;
			}
			else {
				double max_i = spheres[i].center[axis] + spheres[i].radius;
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

bool Blood_Vessel::isSphereInsideBlood_Vessel(const Sphere &sph) const{

    if (sph.object_type == 3 && sph.object_id == id){
        return false;
    }
    
    if (isNearBlood_Vessel(sph.center, 2*sph.radius + barrier_tickness)){ // if near axon
        if(!(sph.object_type == 0 && sph.object_id == id)){ 
            std::vector<std::vector<int>> spheres_id_to_check;
            for (auto axis = 0; axis < 3; ++axis) {
                spheres_id_to_check.push_back(checkAxisForCollision(sph, axis)); // check for collision along 1 axis
                if (spheres_id_to_check[axis].size() == 0){
                    return false;
                }
            }
            // find common ids in all 3 axes
            auto spheres_to_check_all_axes = Obstacle::findCommonIntegers(spheres_id_to_check[0], spheres_id_to_check[1], spheres_id_to_check[2]);
            for (auto i = 0; i < spheres_to_check_all_axes.size(); ++i) {
                Sphere sphere_to_check = spheres[spheres_to_check_all_axes[i]];
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
