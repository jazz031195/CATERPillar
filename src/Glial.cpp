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

void Glial::update_ramifications_boxes(const Sphere &sph){

    int branch_id = sph.branch_id;
    double sph_highest_x_val = sph.center[0] + sph.radius;
    double sph_lowest_x_val = sph.center[0] - sph.radius;
    double sph_highest_y_val = sph.center[1] + sph.radius;
    double sph_lowest_y_val = sph.center[1] - sph.radius;
    double sph_highest_z_val = sph.center[2] + sph.radius;
    double sph_lowest_z_val = sph.center[2] - sph.radius;

    if (sph_highest_x_val > ramification_boxes[branch_id][0][1]){
        ramification_boxes[branch_id][0][1] = sph_highest_x_val;
    }
    if (sph_lowest_x_val < ramification_boxes[branch_id][0][0]){
        ramification_boxes[branch_id][0][0] = sph_lowest_x_val;
    }
    if (sph_highest_y_val > ramification_boxes[branch_id][1][1]){
        ramification_boxes[branch_id][1][1] = sph_highest_y_val;
    }
    if (sph_lowest_y_val < ramification_boxes[branch_id][1][0]){
        ramification_boxes[branch_id][1][0] = sph_lowest_y_val;
    }
    if (sph_highest_z_val > ramification_boxes[branch_id][2][1]){
        ramification_boxes[branch_id][2][1] = sph_highest_z_val;
    }
    if (sph_lowest_z_val < ramification_boxes[branch_id][2][0]){
        ramification_boxes[branch_id][2][0] = sph_lowest_z_val;
    }
    
}


void Glial::update_Box(){


    std::vector<Sphere> spheres_to_add;

    for (int i = 0; i < ramification_spheres.size(); i++){
        for (int j = 0; j < ramification_spheres[i].size(); j++){
            spheres_to_add.push_back(ramification_spheres[i][j]);
        }
    }

    spheres_to_add.push_back(soma);

    for (int i = 0; i < spheres_to_add.size(); i++){

        Sphere sphere_to_add = spheres_to_add[i];
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

        if (this->Box.empty()){
            // intialise box to soma
            // create box around that one sphere
            sph_highest_x_val = soma.center[0]+ soma.radius;
            sph_lowest_x_val = soma.center[0] -soma.radius;
            sph_highest_y_val = soma.center[1] +soma.radius;
            sph_lowest_y_val = soma.center[1] -soma.radius;
            sph_highest_z_val = soma.center[2] +soma.radius;
            sph_lowest_z_val = soma.center[2] -soma.radius;
            //x
            Box.push_back({sph_lowest_x_val, sph_highest_x_val});
            //y
            Box.push_back({sph_lowest_y_val, sph_highest_y_val});
            //z
            Box.push_back({sph_lowest_z_val, sph_highest_z_val});

        }
        else{
            // take values in the box
            sph_highest_x_val = Box[0][1];
            sph_lowest_x_val = Box[0][0];
            sph_highest_y_val = Box[1][1];
            sph_lowest_y_val = Box[1][0];
            sph_highest_z_val = Box[2][1];
            sph_lowest_z_val = Box[2][0];
        }
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

void Glial::update_Box(const Sphere &sph){

    if (Box.empty()){
        cerr << "Error: Box is empty. Cannot update Box with sphere." << endl;
        assert(0);
    }

    update_ramifications_boxes(sph);

    double sph_highest_x_val = sph.center[0] + sph.radius;
    double sph_lowest_x_val = sph.center[0] - sph.radius;
    double sph_highest_y_val = sph.center[1] + sph.radius;
    double sph_lowest_y_val = sph.center[1] - sph.radius;
    double sph_highest_z_val = sph.center[2] + sph.radius;
    double sph_lowest_z_val = sph.center[2] - sph.radius;

    // Update Box values
    if (sph_highest_x_val > Box[0][1]) Box[0][1] = sph_highest_x_val;
    if (sph_lowest_x_val < Box[0][0]) Box[0][0] = sph_lowest_x_val;
    if (sph_highest_y_val > Box[1][1]) Box[1][1] = sph_highest_y_val;
    if (sph_lowest_y_val < Box[1][0]) Box[1][0] = sph_lowest_y_val;
    if (sph_highest_z_val > Box[2][1]) Box[2][1] = sph_highest_z_val;
    if (sph_lowest_z_val < Box[2][0]) Box[2][0] = sph_lowest_z_val;
}


bool Glial::isNearGlialCell(const Eigen::Vector3d &position, const double &distance_to_be_inside) const{
    

    if (Box.empty()){
        cerr << "Error: Box is empty. Cannot check proximity to glial cell." << endl;
        assert(0);
        return false;
    }
    if ((position[0] >=  Box[0][0]- distance_to_be_inside)  && (position[0] <= Box[0][1] + distance_to_be_inside)){
        if ((position[1] >= Box[1][0] - distance_to_be_inside) && (position[1] <=  Box[1][1] + distance_to_be_inside)){
            if ((position[2] >= Box[2][0] - distance_to_be_inside) && (position[2] <= Box[2][1] + distance_to_be_inside)){
                return true; // point is inside the glial cell
            }
        }
    }
    return false;
}

bool Glial::collidesWithRamification(const Eigen::Vector3d &position, const double &distance_to_be_inside) const{

    if (ramification_boxes.empty()){
        return false;
    }

    std::vector<int> indices_close_boxes;

    int i = 0;
    for (const auto& box : ramification_boxes) {
        if ((position[0] >= box[0][0] - distance_to_be_inside) && (position[0] <= box[0][1] + distance_to_be_inside)) {
            if ((position[1] >= box[1][0] - distance_to_be_inside) && (position[1] <= box[1][1] + distance_to_be_inside)) {
                if ((position[2] >= box[2][0] - distance_to_be_inside) && (position[2] <= box[2][1] + distance_to_be_inside)) {
                    indices_close_boxes.push_back(i); // store index of box that is close to the point
                }
            }
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

bool Glial::collides_with_GlialCell(const Sphere &sph, const double &distance_to_be_inside) const{


    if (sph.object_type == 1 && sph.object_id == this->id || sph.object_type == 0){
        return false; // do not collide with itself
    }  
    else if (!isNearGlialCell(sph.center, sph.radius + distance_to_be_inside)) // check if sphere is near glial cell
    {
        return false; 
    } 
    else if (soma.CollideswithSphere(sph, distance_to_be_inside)) // check if sphere collides with soma
    {
        return true; // overlap
    }
    else if (collidesWithRamification(sph.center, sph.radius + distance_to_be_inside)) // overlap
    {
        return true;
    } 
    else{
        return false; // no overlap
    } 

}

void Glial::compute_processes_icvf(const int &factor) {
    volume_processes = 0.0;
    for (int b = 0; b < ramification_spheres.size(); ++b) {
        for (int i = factor; i < ramification_spheres[b].size(); i += factor) {
            if ((soma.center-ramification_spheres[b][i-factor].center).norm() < soma.radius) {
                continue;
            }
            Sphere &s = ramification_spheres[b][i-factor];
            Sphere &s_next = ramification_spheres[b][i];
            double distance = (s_next.center - s.center).norm();
            double v = M_PI * (s.radius * s.radius + s_next.radius * s_next.radius + s.radius * s_next.radius) * distance / 3.0;
            volume_processes += v;

        }
    }
}