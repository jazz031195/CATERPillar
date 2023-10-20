#include "Axon.h"
#include "constants.h"
#include "Eigen/Dense"
#include <iostream>


using namespace Eigen;
using namespace std;


Axon::Axon(const Axon &ax)
{
    id = ax.id;
    spheres = ax.spheres;
    radius = ax.radius;
    begin = ax.begin;
    end = ax.end;
    Box = ax.Box;

};



void Axon::keep_only_first_sphere(){
    Sphere s = spheres[0];
    spheres.clear();
    Box.clear();
    add_sphere(s);

}

void Axon::destroy(){
    spheres.clear();
    Box.clear();
}


void Axon::add_sphere(Sphere sphere_to_add){

    // value of center of sphere at x that has the highest x center value
    double sph_highest_x_val;
    // value of center of sphere at y that has the highest y center value
    double sph_highest_y_val;
    // value of center of sphere at x that has the lowest x center value
    double sph_lowest_x_val;
    // value of center of sphere at y has the lowest y center value
    double sph_lowest_y_val;

    // add sphere to list of spheres
    this->spheres.push_back(sphere_to_add);


    // if there is only one sphere in list
    if (spheres.size() == 1){
        // create box around that one sphere
        sph_highest_x_val = sphere_to_add.center[0]+ sphere_to_add.radius;
        sph_lowest_x_val = sphere_to_add.center[0] -sphere_to_add.radius;
        sph_highest_y_val = sphere_to_add.center[1] +sphere_to_add.radius;
        sph_lowest_y_val = sphere_to_add.center[1] -sphere_to_add.radius;

        Box.clear();
        //x
        Box.push_back({sph_lowest_x_val, sph_highest_x_val});
        //y
        Box.push_back({sph_lowest_y_val, sph_highest_y_val});
        //z
        Box.push_back({sphere_to_add.center[2] - sphere_to_add.radius, sphere_to_add.center[2] + sphere_to_add.radius});

        this->begin = sphere_to_add.center;

    }
    else{

        // take values in the box
        sph_highest_x_val = Box[0][1];
        sph_lowest_x_val = Box[0][0];
        sph_highest_y_val = Box[1][1];
        sph_lowest_y_val = Box[1][0];

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
        // modify box in z direction
        Box[2][1] =  sphere_to_add.center[2] + sphere_to_add.radius;

    }

}


bool check_with_edge(Eigen::Vector3d position, Vector2d x_limits, Vector2d y_limits){ 
    if ((position[0] >=  x_limits[0])  && (position[0] <= x_limits[1])){
        if ((position[1] >= y_limits[0]) && (position[1] <=  y_limits[1])){
            return true;
        }
    }
    return false;
} 

bool Axon::isNearAxon(Eigen::Vector3d position, double distance_to_be_inside){

    if (spheres.size()==0){
        return false;
    }
    Eigen::Vector2d x_limits = Box[0];
    Eigen::Vector2d y_limits = Box[1];

    if(check_with_edge(position, x_limits+ Eigen::Vector2d{-distance_to_be_inside,distance_to_be_inside}, y_limits+ Eigen::Vector2d{-distance_to_be_inside,distance_to_be_inside})){

        return true;
    }  
    else{
        return false;
    }

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
std::vector<int> Axon::checkAxisForCollision(Sphere sph, int axis){

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

bool Axon::isSphereInsideAxon(Sphere sph){

    if(isNearAxon(sph.center, 2*sph.radius + barrier_tickness)){ // if near axon
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
            Sphere sphere_to_check = spheres[spheres_to_check_all_axes[i]];
            if (sph.minDistance(sphere_to_check.center) <= sphere_to_check.radius+ barrier_tickness){
                return true;
            }
        }
        spheres_id_to_check.clear();
        spheres_to_check_all_axes.clear();
    }

    return false;
}


