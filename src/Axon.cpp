#include "Axon.h"
#include "constants.h"
#include "swipeprune.h"
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
    projections = ax.projections;

};


void Axon::add_sphere(Dynamic_Sphere sphere_to_add){

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

        this->end = sphere_to_add.center;
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

        this->end = sphere_to_add.center;
    }


    add_projection(sphere_to_add);


}

void Axon::add_projection(Dynamic_Sphere sphere_to_add){

    // 2d limits of axon
    Vector2d x_limits;
    Vector2d y_limits;
    // projections are in descending order. When added, it is added at the right position.
    for (unsigned axis = 0; axis < 3; ++axis) {

        double position1, position2;
        int sphere_id = spheres.size()-1;

        // projections
        // center + radius
        position1 = sphere_to_add.center[axis] + sphere_to_add.radius;
        Projections::projection_pt p1 {position1, id, sphere_id};
        // center - radius
        position2 = sphere_to_add.center[axis] - sphere_to_add.radius;
        Projections::projection_pt p2 {position2, id, sphere_id};

        projections.append_right_place(p1, p2, axis);
        
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

    Eigen::Vector2d x_limits = Box[0];
    Eigen::Vector2d y_limits = Box[1];

    if(check_with_edge(position, x_limits+ Eigen::Vector2d{-distance_to_be_inside,distance_to_be_inside}, y_limits+ Eigen::Vector2d{-distance_to_be_inside,distance_to_be_inside})){

        return true;
    }  
    else{
        return false;
    }

}



bool Axon::isPosInsideAxon(Eigen::Vector3d &position,  double distance_to_be_inside, double max_radius){
    // when checking collision with walker -> check with normal radius
    // when checking with collisions of other axons -> check with max_radius so there is room for swelling
    std::vector<std::vector<Projections::projection_pt>> coliding_projs;
    bool colliding_all_axes;
    Dynamic_Sphere sphere_ ;
    double rad;

    // if position is in box with axon inside
    if(isNearAxon(position, distance_to_be_inside)){
        // find all projections in between the two projections of the edges
        
        coliding_projs = projections.find_collisions_all_axes(position, max_radius + barrier_tickness, id, distance_to_be_inside);

        if (coliding_projs.size() == 3){ 
          
            // for all coliding objects in x 
            for(unsigned j = 0; j < coliding_projs[0].size() ; j++){ 

                const Projections::projection_pt coliding_proj = coliding_projs[0][j];
                // if the same coliding objects are also in y and z but are not from same objects

                colliding_all_axes = (projections.isProjInside(coliding_projs[1], coliding_proj) && projections.isProjInside(coliding_projs[2], coliding_proj));

                if (colliding_all_axes){

                    sphere_ = spheres[coliding_proj.sph_id];
                    
                    if (sphere_.minDistance(position) <= distance_to_be_inside){ 

                        return true;
                        
                        /*
                        cout << " Axon : "<< id <<"Position : [" << position[0] << ", "<< position[1] << ", " << position[2] << "]" << endl; 
                        cout << "           distance to sphere :" << coliding_proj.sph_id << " : " <<  sphere_.minDistance(position) ;
                        cout << ", Sphere position : [" << sphere_.center[0] << ", "<< sphere_.center[1] << ", " << sphere_.center[2] << "]" << endl;  
                        */
                    }  
                }
            }
        }

        return false;
        
    }
    else{
        return false;
    }
    return false;
} 


void Axon::threadGrowth() {
    std::cout << "Growing Axon number: \n " << id << std::endl;
}