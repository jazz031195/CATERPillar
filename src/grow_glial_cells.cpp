#include "axongammadistribution.h"
#include "grow_glial_cells.h"
#include <algorithm> // std::sort
#include <random>
#include <chrono>
#include <future>
#include "Eigen/Dense"
#include <thread>
#include "threads.h"

using namespace std;
using namespace Eigen;
using namespace std::chrono;

GlialCellGrowth::GlialCellGrowth() : CellGrowth(), glial_cell_to_grow(*(Glial*)nullptr)  {}

GlialCellGrowth::~GlialCellGrowth() {}

GlialCellGrowth::GlialCellGrowth(Glial &glial_cell_to_grow_,
            const std::vector<Glial> &astrocytes,
            const std::vector<Glial> &oligodendrocytes,
            const std::vector<Axon> &axons_,
            const Eigen::Vector3d &extended_min_limits_,
            const Eigen::Vector3d &extended_max_limits_,
            const Eigen::Vector3d &min_limits_,
            const Eigen::Vector3d &max_limits_,
            const double &std_dev_,
            const double &min_radius_) : CellGrowth(axons_, astrocytes, oligodendrocytes, extended_min_limits_, extended_max_limits_, min_limits_, max_limits_, std_dev_, min_radius_), glial_cell_to_grow(glial_cell_to_grow_){}

GlialCellGrowth::GlialCellGrowth(const GlialCellGrowth &other)
  : CellGrowth(other),             // Call the base-class copy constructor
    glial_cell_to_grow(other.glial_cell_to_grow)  // references must be bound here
{}


void GlialCellGrowth::add_spheres(Sphere &sph, const Sphere &last_sphere, const bool &check_collision_with_branches, const int &factor, const int &index_ram_spheres){
    
    // nbr of spheres to add in between
    int nbr_spheres = factor - 1;
    if (factor > 1){
        
        // distance between two consecutive spheres
        double distance = (sph.center - last_sphere.center).norm();
        Eigen::Vector3d vector = (sph.center - last_sphere.center).normalized();
        double distance_between_spheres = distance/(nbr_spheres+1);
        int id_;
        for (int i = 0 ; i < nbr_spheres; i++){
            Eigen::Vector3d position = last_sphere.center + vector*distance_between_spheres*(i+1);
            double rad = last_sphere.radius + (sph.radius - last_sphere.radius)*(i+1)/(nbr_spheres+1);
            id_ = last_sphere.id + i + 1;
            Sphere s(id_, last_sphere.object_id, last_sphere.object_type, position, rad, last_sphere.branch_id, sph.parent_id);
            bool can_grow_;
            if (check_collision_with_branches){
                can_grow_ = canSpherebePlaced(s) && !collideswithOtherBranches(s);
            }
            else{
                can_grow_ = canSpherebePlaced(s);
            }
            if(can_grow_){
                glial_cell_to_grow.ramification_spheres[index_ram_spheres].push_back(s);
            }

        }
    }
    sph.id = last_sphere.id + nbr_spheres + 1;
    glial_cell_to_grow.ramification_spheres[index_ram_spheres].push_back(sph); 
    glial_cell_to_grow.update_Box();
}


bool GlialCellGrowth::AddOneSphere(const double &radius_, const bool &create_sphere, int &grow_straight, const int &i, const bool &check_collision_with_branches, const int &parent, const int &factor)
{
    grow_straight = 0;
    // if sphere collides with environment
    bool can_grow_;
 
    if (glial_cell_to_grow.ramification_spheres[i].size() > 100000)
    {
        return false;
    }

    //cout << "ramification spheres : " << glial_cell_to_grow.ramification_spheres[i].size() << endl;
    double distance;
    Eigen::Vector3d destination = glial_cell_to_grow.attractors[i];
    //cout << "destination : " << destination << endl;
    if (glial_cell_to_grow.ramification_spheres.size() == 0 ){
        assert(0);
    }
    Sphere last_sphere = glial_cell_to_grow.ramification_spheres[i][glial_cell_to_grow.ramification_spheres[i].size() - 1];
    // the distance between two consecutive spheres is the maximum radius / 2
    double max_radius_ = max(radius_, last_sphere.radius);

    distance = (max_radius_);

    can_grow_ = false;
    int tries = 0;
    int id_ = last_sphere.id + factor;
    Sphere s(id_, glial_cell_to_grow.id, 1, {0,0,0}, radius_, i);
    while (!can_grow_ && tries < 10){

        if (std_dev != 0.0)
        {
            if (grow_straight==1)
            {
                find_next_center_straight(distance, s, glial_cell_to_grow.ramification_spheres[i]);
            }
            else
            {
                find_next_center(s, distance, glial_cell_to_grow.ramification_spheres[i], destination, false);
            }
            // check if there is a collision
            if (check_collision_with_branches){
                can_grow_ = canSpherebePlaced(s) && !collideswithOtherBranches(s);
            }
            else{
                can_grow_ = canSpherebePlaced(s);
            }
        }
        else
        {
            cout << "no std dev" << endl;
            assert(0);  
            find_next_center(s, distance, glial_cell_to_grow.ramification_spheres[i], destination, false);
            if (check_collision_with_branches){
                can_grow_ = canSpherebePlaced(s) && !collideswithOtherBranches(s);
            }
            else{
                can_grow_ = canSpherebePlaced(s);
            }
        }
        if (!can_grow_){
            tries += 1;
        }
    }
    
    if (can_grow_ )
    {
        if (create_sphere){
            //cout << "adding sphere to glial, radius "<< s.radius << ", psotion : " << s.center << endl;
            s.parent_id = parent;
            add_spheres(s, last_sphere, check_collision_with_branches, factor, i);
        }
        // if is not in inside voxel
        if (!check_borders(min_limits, max_limits, s.center, s.radius)){ 
            finished = true;
            //cout << "hits end of voxel" << endl;
            return false;
        }
        else{
            //cout << "can grow" << endl;

            return true;
        }
    }
    else // collides and >1000 tries
    {
        //cout << "could not grow glial" << endl;
        return false;
    }
        
}

bool GlialCellGrowth::collideswithOtherBranches(Sphere &sph){

    // if collides with own soma
    if (glial_cell_to_grow.soma.CollideswithSphere(sph, barrier_tickness)){
        return true;
    }

    // check other branches of same glial cell
    for (long unsigned int i = 0; i < glial_cell_to_grow.ramification_spheres.size(); i++){
        if(glial_cell_to_grow.ramification_spheres[i].size() > 0){ 
            std::vector<Sphere> branch = glial_cell_to_grow.ramification_spheres[i];
            if (branch[0].branch_id != sph.branch_id){
                for (long unsigned int k = 0; k < branch.size(); k++){
                    Sphere sph_ = branch[k];
                    if (sph_.CollideswithSphere(sph, barrier_tickness)){
                        return true;
                    }
                }
            }
        } 
    } 
    


    return false;
}

void GlialCellGrowth::find_next_center_straight(double distance, Sphere &s, const std::vector<Sphere> &spheres)
{

    Eigen::Vector3d normal = (spheres[spheres.size() - 1].center-glial_cell_to_grow.soma.center).normalized();
    Eigen::Vector3d new_center = spheres[spheres.size() - 1].center + normal*distance;
    s.center = new_center;
}

void GlialCellGrowth::find_next_center(Sphere &s,  double dist_, const std::vector<Sphere> &spheres, const Eigen::Vector3d &target, const bool &is_axon)
{

    Eigen::Vector3d target_ = target;
    Eigen::Vector3d vector_to_target = target_ - spheres[spheres.size() - 1].center;
    vector_to_target = vector_to_target.normalized();
    Eigen::Vector3d vector = generate_random_point_on_sphere(std_dev);
    vector = apply_bias_toward_target(vector, vector_to_target);
    Eigen::Vector3d position = spheres[spheres.size() - 1].center + dist_ * vector.normalized();
    s.center = position;
}