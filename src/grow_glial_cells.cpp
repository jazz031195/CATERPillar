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


bool GlialCellGrowth::canSpherebePlaced(Sphere &sph){


    bool axons_check = checkAxonsOverlap(sph);

    if (!axons_check){
        return false;
    }
    // check collision other glial cells 
    for (auto &glial : glial_cells)
    {
        if (glial.collides_with_GlialCell(sph, barrier_tickness)){
            return false;
        } 
        
    }


    return true;
} 
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
                can_grow_ = canSpherebePlaced(s) && !collideswithItself(s);
            }
            else{
                can_grow_ = canSpherebePlaced(s);
            }
            if(can_grow_){
                glial_cell_to_grow.ramification_spheres[index_ram_spheres].push_back(s);
                glial_cell_to_grow.update_Box(s);
            }


        }
    }
    sph.id = last_sphere.id + nbr_spheres + 1;
    glial_cell_to_grow.ramification_spheres[index_ram_spheres].push_back(sph); 
    glial_cell_to_grow.update_Box(sph);
}


bool GlialCellGrowth::AddOneSphere(const double &radius_, const bool &create_sphere, int &grow_straight, const int &i, const bool &check_collision_with_branches, const int &parent, const int &factor)
{
    if (glial_cell_to_grow.ramification_boxes.size() != glial_cell_to_grow.ramification_spheres.size())
    {
        std::cout << "ramification boxes and spheres are not the same size" << std::endl;
        cout << "ramification boxes size: " << glial_cell_to_grow.ramification_boxes.size() << endl;
        cout << "ramification spheres size: " << glial_cell_to_grow.ramification_spheres.size() << endl;
        assert(0);
    }
 
    if (glial_cell_to_grow.ramification_spheres[i].size() > 1e10)
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

    bool can_grow_ = false;

    int id_ = last_sphere.id + factor;
    Sphere s(id_, glial_cell_to_grow.id, 1, {0,0,0}, radius_, i);

    int tries = 0;
    int threshold_tries = 100;

    while (!can_grow_ && tries < threshold_tries){

        find_next_center(s, distance, glial_cell_to_grow.ramification_spheres[i], destination);
        // check if there is a collision
        if (check_collision_with_branches){
            can_grow_ = canSpherebePlaced(s) && !collideswithItself(s);
        }
        else{
            can_grow_ = canSpherebePlaced(s);
        }
        tries += 1;
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

Eigen::Vector3d GlialCellGrowth::readapt_sphere_position(const Sphere &s, const Axon &neighbour_axon, bool can_readapt){

    Eigen::Vector3d sphere_center = s.center;
    for (int i = 0; i < neighbour_axon.outer_spheres.size(); i++) {
        const Sphere& other_sphere = neighbour_axon.outer_spheres[i];
        double distance = (sphere_center - other_sphere.center).norm();
        double overlap = s.radius + other_sphere.radius + 1e-6 - distance;
        if (overlap > 0){
            Eigen::Vector3d vector = (sphere_center - other_sphere.center).normalized()*(overlap);
            sphere_center += vector;
        }
    }
    int branch_id = glial_cell_to_grow.ramification_spheres.size()-1;
    const Sphere& prev_sphere = glial_cell_to_grow.ramification_spheres[branch_id].back();
    double new_distance = (sphere_center - prev_sphere.center).norm();
    if (new_distance > max(prev_sphere.radius, s.radius)){
        can_readapt = true;
        return sphere_center;
    }
    else{
        can_readapt = false;
        return {0, 0, 0}; // Return a zero vector if the sphere cannot be placed
    }
}


bool GlialCellGrowth::checkAxonsOverlap(Sphere &sph){

    for (const auto& axon : axons) {
        if (axon.isSphereInsideAxon(sph)) {
            return false;  
            
        }
    }
    return true; 
    
}


bool GlialCellGrowth::collideswithItself(Sphere &sph){

    // if collides with own soma
    if (glial_cell_to_grow.soma.CollideswithSphere(sph, barrier_tickness)){
        return true;
    }

    double soma_radius = glial_cell_to_grow.soma.radius;
    double x_min = glial_cell_to_grow.soma.center[0] - soma_radius - barrier_tickness;
    double x_max = glial_cell_to_grow.soma.center[0] + soma_radius + barrier_tickness;
    double y_min = glial_cell_to_grow.soma.center[1] - soma_radius - barrier_tickness;
    double y_max = glial_cell_to_grow.soma.center[1] + soma_radius + barrier_tickness;
    double z_min = glial_cell_to_grow.soma.center[2] - soma_radius - barrier_tickness;
    double z_max = glial_cell_to_grow.soma.center[2] + soma_radius + barrier_tickness;
    int branch_id = sph.branch_id;

    if (glial_cell_to_grow.ramification_boxes.size() != glial_cell_to_grow.ramification_spheres.size())
    {
        std::cout << "ramification boxes and spheres are not the same size" << std::endl;
        cout << "ramification boxes size: " << glial_cell_to_grow.ramification_boxes.size() << endl;
        cout << "ramification spheres size: " << glial_cell_to_grow.ramification_spheres.size() << endl;
        assert(0);
    }

    // check other branches of same glial cell
    for (long unsigned int i = 0; i < glial_cell_to_grow.ramification_boxes.size(); i++)
    {
        if (glial_cell_to_grow.ramification_boxes[i].empty())
        {
            continue; // skip empty boxes
        }

        if (i == branch_id)
        {
            continue; // skip own branch
        }

        // if inside box 
        if (x_min < glial_cell_to_grow.ramification_boxes[i][0][1] &&
            x_max > glial_cell_to_grow.ramification_boxes[i][0][0] &&
            y_min < glial_cell_to_grow.ramification_boxes[i][1][1] &&
            y_max > glial_cell_to_grow.ramification_boxes[i][1][0] &&
            z_min < glial_cell_to_grow.ramification_boxes[i][2][1] &&
            z_max > glial_cell_to_grow.ramification_boxes[i][2][0])
        {
            // check if collides with spheres in branch
            if (glial_cell_to_grow.ramification_spheres[i].size() > 0)
            {
                std::vector<Sphere> branch = glial_cell_to_grow.ramification_spheres[i];
                for (long unsigned int k = 0; k < branch.size(); k++)
                {
                    Sphere sph_ = branch[k];
                    if (sph_.CollideswithSphere(sph, barrier_tickness))
                    {
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

void GlialCellGrowth::find_next_center(Sphere &s,  double dist_, const std::vector<Sphere> &spheres, const Eigen::Vector3d &target)
{

    Eigen::Vector3d target_ = target;
    Eigen::Vector3d vector_to_target = target_ - spheres[spheres.size() - 1].center;
    vector_to_target = vector_to_target.normalized();
    double std_ = 0.4;
    Eigen::Vector3d vector = generate_random_point_on_sphere(std_);
    vector = apply_bias_toward_target(vector, vector_to_target);
    Eigen::Vector3d position = spheres[spheres.size() - 1].center + dist_ * vector.normalized();
    s.center = position;
}