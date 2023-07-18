#include "swipeprune.h"
#include <algorithm> // std::sort
#include <random>
#include <iostream>
#include <list>

using namespace std;
using namespace Eigen;

Projections::Projections(){
    sph_projections_x.clear();
    sph_projections_y.clear();
    sph_projections_z.clear();
}

void Projections::clear_projections(){
    sph_projections_x.clear();
    sph_projections_y.clear();
    sph_projections_z.clear();
}

bool isSmallerProj(const Projections::projection_pt& obj1, const Projections::projection_pt& obj2) {
    return obj1.position < obj2.position;
}

void Projections::append_right_place(Projections::projection_pt p1, Projections::projection_pt p2, int axis){
    

    std::vector<projection_pt> axon_projection_on_axis;

    if (axis == 0){
        axon_projection_on_axis = sph_projections_x;
    }
    else if (axis == 1){
        axon_projection_on_axis = sph_projections_y;
    }
    else if (axis == 2){
        axon_projection_on_axis = sph_projections_z;
    }

    // Find the position to insert the new object
    std::vector<Projections::projection_pt>::iterator index_min_it = std::lower_bound(axon_projection_on_axis.begin(), axon_projection_on_axis.end(), p1, isSmallerProj);
    
    // Insert the new object at the found position
    axon_projection_on_axis.insert(index_min_it, p1);

    // Find the position to insert the new object
    std::vector<Projections::projection_pt>::iterator index_max_it = std::lower_bound(axon_projection_on_axis.begin(), axon_projection_on_axis.end(), p2, isSmallerProj);
    
    // Insert the new object at the found position
    axon_projection_on_axis.insert(index_max_it, p2);
                

    
    if (axis == 0){
        sph_projections_x = axon_projection_on_axis;
    }
    else if (axis == 1){
        sph_projections_y = axon_projection_on_axis;
    }
    else if (axis == 2){
        sph_projections_z = axon_projection_on_axis;
    }
}


std::vector<Projections::projection_pt> Projections::find_collisions(projection_pt proj_on_axis_min, projection_pt proj_on_axis_max,std::vector<projection_pt> projections_on_axis){
    
    std::vector<projection_pt> closest_spheres;

    if (projections_on_axis.size()== 0){
        return closest_spheres; 
    } 

    // Find the position to insert the new object
    std::vector<Projections::projection_pt>::iterator index_min_it = std::lower_bound(projections_on_axis.begin(), projections_on_axis.end(), proj_on_axis_min, isSmallerProj);
    int index_min = std::distance(projections_on_axis.begin(), index_min_it);
    std::vector<Projections::projection_pt>::iterator index_max_it = std::lower_bound(projections_on_axis.begin(), projections_on_axis.end(), proj_on_axis_max, isSmallerProj);
    int index_max = std::distance(projections_on_axis.begin(), index_max_it);

    if (index_min == index_max){
        return closest_spheres; 
    } 

    for (unsigned i = index_min ; i < index_max; i++){
        
        projection_pt s{projections_on_axis[i].position, projections_on_axis[i].axon_id, projections_on_axis[i].sph_id};  
        closest_spheres.push_back(s);
        
    } 
    return closest_spheres;
}

bool Projections::isProjInside(std::vector<Projections::projection_pt> projs, Projections::projection_pt p){
    
    // search for s in spheres_
    if (projs.size() == 0){
        return false;
    } 
    for (unsigned i = 0; i < projs.size(); i++){
        if (projs[i].axon_id == p.axon_id && projs[i].sph_id == p.sph_id){
            return true;

        }  
    } 
    return false;
}  


std::vector<std::vector<Projections::projection_pt>> Projections::find_collisions_all_axes(Vector3d &position, double rad, int ax_id, double distance_to_be_inside){
    std::vector<std::vector<projection_pt>> coliding_projs;
    std::vector<projection_pt> colisions_axis_projs;

    // on all axes
    for (unsigned x = 0; x < 3; x++){

        colisions_axis_projs.clear();
        
        projection_pt proj_on_axis_min {position[x]- rad -distance_to_be_inside, ax_id, 1000};
        // get max projection
        projection_pt proj_on_axis_max {position[x] + rad + distance_to_be_inside, ax_id, 1000};

        if (x== 0){
            
            colisions_axis_projs = find_collisions(proj_on_axis_min, proj_on_axis_max, sph_projections_x);

        }  
        else if (x == 1) {

            colisions_axis_projs = find_collisions(proj_on_axis_min, proj_on_axis_max, sph_projections_y);
        }  
        else{

            colisions_axis_projs = find_collisions(proj_on_axis_min, proj_on_axis_max, sph_projections_z);
            
        } 
        if (colisions_axis_projs.size()== 0){
            return coliding_projs;

        } 
        else{ 
            coliding_projs.push_back(colisions_axis_projs);

        } 
    }
    return coliding_projs;
}


