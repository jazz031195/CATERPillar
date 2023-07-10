#include "swipeprune.h"
#include <algorithm> // std::sort
#include <random>


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

void Projections::append_right_place(Projections::projection_pt p1, Projections::projection_pt p2, int axis){
    
    double position1 = p1.position;
    double position2 = p2.position;

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

    //  find the first element bigger than the new position
    auto pos1 = std::find_if(axon_projection_on_axis.begin(), axon_projection_on_axis.end(), [position1](projection_pt s) {
        return s.position >= position1 ;
    });
    // And then insert the new element at this position
    axon_projection_on_axis.insert(pos1, p1);

    //  find the first element bigger than the new position
    auto pos2 = std::find_if(axon_projection_on_axis.begin(), axon_projection_on_axis.end(), [position2](projection_pt s) {
        return s.position >= position2 ;
    });
                
    // And then insert the new element at this position
    axon_projection_on_axis.insert(pos2, p2);
    
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


std::vector<Projections::projection_pt> Projections::find_collisions(projection_pt proj_on_axis_min, projection_pt proj_on_axis_max,std::vector<projection_pt> projections_on_axis, double distance_to_be_inside){
    
    std::vector<projection_pt> closest_spheres;

    if (projections_on_axis.size()== 0){
        return closest_spheres; 
    } 

    // projection after which projections are bigger than min
    auto pos_min = std::find_if(projections_on_axis.begin(), projections_on_axis.end(), [proj_on_axis_min, distance_to_be_inside](projection_pt s) {
        return s.position >= proj_on_axis_min.position + distance_to_be_inside ;
    });

    // projection index after which projections are bigger than min
    unsigned index_min = std::distance(std::begin(projections_on_axis), pos_min);  

    // projection index after which projections are bigger than max
    auto pos_max = std::find_if(projections_on_axis.begin(), projections_on_axis.end(), [proj_on_axis_max, distance_to_be_inside](projection_pt s) {
        return s.position >= proj_on_axis_max.position - distance_to_be_inside;
    });

    // projection index after which projections are bigger than max
    unsigned index_max = std::distance(std::begin(projections_on_axis), pos_max);

    if (index_min == index_max){
        return closest_spheres; 
    } 

    for (unsigned i = index_min ; i < index_max; i++){
        
        projection_pt s{projections_on_axis[i].position, projections_on_axis[i].axon_id, projections_on_axis[i].sph_id};  
        // if s not already in closest_spheres
        if (!isProjInside(closest_spheres,s)){ 
            closest_spheres.push_back(s);
        } 
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
            break;
        }  
    } 
    return false;
}  


std::vector<std::vector<Projections::projection_pt>> Projections::find_collisions_all_axes(Vector3d &position, double rad, int ax_id, double distance_to_be_inside){
    std::vector<std::vector<projection_pt>> coliding_projs;
    std::vector<projection_pt> colisions_axis_projs;
    coliding_projs.clear();
    // on all axes
    for (unsigned x = 0; x < 3; x++){

        colisions_axis_projs.clear();
        
        projection_pt proj_on_axis_min {position[x]- rad, ax_id, 1000};
        // get max projection
        projection_pt proj_on_axis_max {position[x] + rad, ax_id, 1000};

        if (x== 0){

            colisions_axis_projs = find_collisions(proj_on_axis_min, proj_on_axis_max, sph_projections_x, distance_to_be_inside);

        }  
        else if (x == 1) {

            colisions_axis_projs = find_collisions(proj_on_axis_min, proj_on_axis_max, sph_projections_y, distance_to_be_inside);
        }  
        else{

            colisions_axis_projs = find_collisions(proj_on_axis_min, proj_on_axis_max, sph_projections_z, distance_to_be_inside);
            
        } 
        if (colisions_axis_projs.size()== 0){
            return coliding_projs;
            break;
        } 
        else{ 
            coliding_projs.push_back(colisions_axis_projs);

        } 
    }
    return coliding_projs;
}


