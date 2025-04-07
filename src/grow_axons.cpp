#include "axongammadistribution.h"
#include "grow_axons.h"
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


AxonGrowth::~AxonGrowth() {}

AxonGrowth::AxonGrowth(Axon &axon_to_grow_,
                       const std::vector<Glial> &astrocytes,
                       const std::vector<Glial> &oligodendrocytes,
                       const std::vector<Axon> &axons_,
                       const Eigen::Vector3d &extended_min_limits_,
                       const Eigen::Vector3d &extended_max_limits_,
                       const Eigen::Vector3d &min_limits_,
                       const Eigen::Vector3d &max_limits_,
                       const double &std_dev_,
                       const double &min_radius_)
    : CellGrowth(axons_, astrocytes, oligodendrocytes,
                 extended_min_limits_, extended_max_limits_,
                 min_limits_, max_limits_,
                 std_dev_, min_radius_),
      axon_to_grow(axon_to_grow_), modified_axons(axons_),
      axons(axons_)
{ ind_axons_pushed.clear(); }

AxonGrowth::AxonGrowth(const AxonGrowth &other)
    : CellGrowth(other), // call base copy constructor
      axon_to_grow(other.axon_to_grow),
      axons(other.axons) {}

/*
bool AxonGrowth::PushOtherAxons(const Sphere &sph) {
    
    if (canSpherebePlaced(sph)){
        cout <<"sphere : " << sph.id << " axon : "<< sph.object_id << " can be placed without push" << endl;
        return true;
    } 
    
    //cout <<"------" << endl;
    std::vector<Axon> old_axons = axons;

    for (int j = 0; j < axons.size(); j++) {
        int index = j;
        //cout << "colliding axon id: " << axons[index].id << endl;
        if (!pushAxonSpheres(axons[index], sph)) {
            axons = std::move(old_axons); // Restore the original axons
            //cout <<"cannot place sphere : " << sph.id << " axon : "<< sph.object_id << endl;
            return false;
        }
    }

    if (!canSpherebePlaced(sph)){
        cout <<"cannot place sphere : " << sph.id << " axon : "<< sph.object_id << endl;
        assert(0);
    }
    
    else{
        cout <<"sphere : " << sph.id << " axon : "<< sph.object_id << " can be placed with push" << endl;
    }
    

    return true;
}
*/

bool AxonGrowth::pushAxonSpheres(Axon &axon, const Sphere &sph) {

    int count = 0;
    for (int i = 0; i < axon.outer_spheres.size(); i++) {
        Sphere &other_sphere = axon.outer_spheres[i];
        Sphere sph_before;
        Sphere sph_after;
        if (i > 0){
            sph_before = axon.outer_spheres[i - 1];
        }
        if (i < axon.outer_spheres.size() - 1){
            sph_after = axon.outer_spheres[i + 1];
        }
        double distance = (sph.center - other_sphere.center).norm();
        double overlap = sph.radius + other_sphere.radius + 1e-6 - distance;

        if (distance < sph.radius + other_sphere.radius) {
            count += 1;
            Eigen::Vector3d push_vector = (other_sphere.center - sph.center).normalized() * overlap;
            Sphere new_sphere = Sphere(other_sphere.id, other_sphere.object_id, other_sphere.object_type,
                                       other_sphere.center + push_vector, other_sphere.radius);
            
            if (i > 0) {
                double distance_to_sphere_before = (new_sphere.center - sph_before.center).norm();
                if (distance_to_sphere_before > max(sph_before.radius, new_sphere.radius)) {
                    return false;
                }
            }

            if (i < axon.outer_spheres.size() - 1) {
                double distance_to_sphere_after = (new_sphere.center - sph_after.center).norm();
                if (distance_to_sphere_after > max(sph_before.radius, new_sphere.radius)) {
                    return false;
                }
            }

            if (canSpherebePlaced(new_sphere, false) && check_borders(extended_min_limits, extended_max_limits, new_sphere.center, new_sphere.radius)) {
                axon.outer_spheres[i] = std::move(new_sphere); // Update the sphere
                //cout <<"sphere : " << new_sphere.id << " axon : "<< new_sphere.object_id << " can be pushed" << endl;
            } 
            else {
                return false; // Cannot place the sphere
            }
        }
    }
    axon.updateBox();

    return true;
}

bool AxonGrowth::canSpherebePlaced(const Sphere &sph, const bool &with_push){

    std::vector<int> axons_pushed_by_sphere;


    for (int i = 0 ; i < axons.size(); i++) {

        if (!(axons[i].id == sph.object_id && sph.object_type == 0))
        {
            // Check overlap
            if (axons[i].isSphereInsideAxon(sph)) 
            {
                
                if (with_push) {
                    Axon axon_to_push = axons[i];
                    
                    bool can_axon_be_pushed = pushAxonSpheres(axon_to_push, sph);
                    if (!can_axon_be_pushed){
                        return false;
                    }
                    else{
                        modified_axons[i] = std::move(axon_to_push);
                        axons_pushed_by_sphere.push_back(i);
                    }
                }
                else{
                    return false;
                }
                
                return false;
            }
        }
    }

    // check collision other glial cells 
    for (auto &glial : glial_cells)
    {
        if (sph.object_type == 0 || (sph.object_type != 0 && sph.object_id != glial.id) ){
            if (glial.isNearGlialCell(sph.center, 2*sph.radius+1e-6)){
                if (glial.collides_with_GlialCell(sph)){
                    return false;
                }
                // check with branches of other glial cells
                for (long unsigned int i = 0; i < glial.ramification_spheres.size(); i++){
                    std::vector<Sphere> branch = glial.ramification_spheres[i];
                    for (long unsigned int k = 0; k < branch.size(); k++){
                        Sphere sph_ = branch[k];
                        if (sph_.CollideswithSphere(sph, barrier_tickness)){
                            return false;
                        }
                    }
                }
            }
        }
    }

    if (axons_pushed_by_sphere.size() > 0){
        for (auto &index : axons_pushed_by_sphere){
            // if index is not already in the list
            if (std::find(ind_axons_pushed.begin(), ind_axons_pushed.end(), index) == ind_axons_pushed.end()){
                ind_axons_pushed.push_back(index);
            }
        }    
    }

    return true;
} 

bool AxonGrowth::isSphereCollidingwithAxon(Axon ax, Sphere sph){

    return ax.isSphereInsideAxon(sph);
}


void AxonGrowth::find_next_center_straight(double distance, Sphere &s, const std::vector<Sphere> &spheres)
{
    if (spheres.size() < 2){
        assert(0);
        cout << "spheres size : " << spheres.size() << endl;
        return;
    }
    Eigen::Vector3d last_center = spheres[spheres.size() - 1].center;
    Eigen::Vector3d before_last_center = spheres[spheres.size() - 2].center;
    Eigen::Vector3d straight_vector = (last_center - before_last_center).normalized();
    straight_vector *= distance;
    Eigen::Vector3d new_center = last_center + straight_vector;
    s.center = new_center;

}


bool AxonGrowth::AddOneSphere(double radius_, bool create_sphere, int grow_straight, const int &factor)
{
    // Basic validation
    if (axon_to_grow.outer_spheres.empty()) {
        std::cerr << "EMPTY AXON!" << std::endl;
        assert(false); // or return false;
    }

    assert(axon_to_grow.growth_axis >= 0 && axon_to_grow.growth_axis < 3);
    
    bool is_allowed_to_stop_early = axon_to_grow.outside_voxel;
    

    // If the last sphere's center is beyond extended_max_limits, the axon is considered fully grown
    Sphere last_sphere = axon_to_grow.outer_spheres.back();
    if (last_sphere.center[axon_to_grow.growth_axis] >= extended_max_limits[axon_to_grow.growth_axis] && !is_allowed_to_stop_early) {
        finished = true;
        return true; // Axon is done
    }
    else if (!check_borders(min_limits, max_limits, last_sphere.center, last_sphere.radius) && is_allowed_to_stop_early) {
        finished = true;
        return true; // Axon is done
    }
    

    // Prepare
    double max_radius_ = std::max(radius_, last_sphere.radius);
    if (axon_to_grow.myelin_sheath && axon_to_grow.inner_radius < max_radius_) {
        max_radius_ = axon_to_grow.inner_radius;
    }
    double distance = max_radius_;
    int threshold_tries = (std_dev == 0.0) ? 1 : 1000;

    // New sphere to attempt placing
    Sphere s(axon_to_grow.outer_spheres.size() + factor,
             axon_to_grow.id,
             /*object_type=*/0,
             axon_to_grow.begin,
             radius_);

    bool can_grow_ = false;
    int tries = 0;

    // Helper lambda: tries to place sphere (find next center, check if inside, collision-check).
    // Returns true if placed successfully, false otherwise.
    auto attemptPlacement = [&](Sphere &candidate, int triesCount, bool with_push) -> bool {
        // 1) Find next center
        if (std_dev != 0.0) {
            // If grow_straight is set, use the straight function, else the normal function
            if (grow_straight == 1) {
                find_next_center_straight(distance, candidate, axon_to_grow.outer_spheres);
            } else {
                find_next_center(candidate, distance, axon_to_grow.outer_spheres, axon_to_grow.end, true);
            }
        } else {
            // If std_dev == 0, we skip the "straight vs. random" logic and always use find_next_center
            find_next_center(candidate, distance, axon_to_grow.outer_spheres, axon_to_grow.end, true);
        }

        // 2) Check if inside
        if (!check_borders(extended_min_limits, extended_max_limits, candidate.center, candidate.radius) && !is_allowed_to_stop_early) {
            return false;
        }


        bool canPlace = canSpherebePlaced(candidate, with_push);

        return canPlace;
        

    };
 

    while (!can_grow_ && tries < threshold_tries)
    {
        bool with_push = false;
        if (tries > threshold_tries/2) {
            with_push = true;
        }
        // Attempt the main placement
        bool success = attemptPlacement(s, tries, with_push);
        if (success) {
            can_grow_ = true;
            break;
        }
        
        // If still not placed, and we used "grow_straight == 1", then attempt a fallback (non-straight) approach
        if (!can_grow_ && (std_dev != 0.0) && (grow_straight == 1)) {
            grow_straight == 0; // Fallback to normal growth
            if (attemptPlacement(s, tries, with_push)) {
                can_grow_ = true;
                break;
            }
        }
        
        tries++;
    }

    // Evaluate the result
    if (!can_grow_) {
        // Collides or max tries reached
        return false;
    }

    // If can_grow_ == true
    if (create_sphere) {
        s.parent_id = last_sphere.id;
        add_spheres(s, last_sphere, factor);
        // Update volume if within the stricter [min_limits, max_limits]
        Sphere newly_added = axon_to_grow.outer_spheres.back(); // s with final coords
        if (check_borders(min_limits, max_limits, newly_added.center, newly_added.radius)) {
            axon_to_grow.update_Volume(factor, min_limits, max_limits);
        }
    }

    // If we reach the edge of voxel after adding the new sphere
    Sphere current_last = axon_to_grow.outer_spheres.back();
    
    if (current_last.center[axon_to_grow.growth_axis] + current_last.radius > extended_max_limits[axon_to_grow.growth_axis]) {
        finished = true;
    }
    else if (!check_borders(min_limits, max_limits, current_last.center, current_last.radius) && is_allowed_to_stop_early) {
        finished = true;
    }
    return true;
}


void AxonGrowth::add_spheres(Sphere &sph, const Sphere &last_sphere, const int &factor){
    
    // nbr of spheres to add in between
    int nbr_spheres = factor - 1;
    int last_id = last_sphere.id;

    if (factor > 1){
        // distance between two consecutive spheres
        double distance = (sph.center - last_sphere.center).norm();
        Eigen::Vector3d vector = (sph.center - last_sphere.center).normalized();
        double distance_between_spheres = distance/(nbr_spheres+1);
        int id_;
        for (int i = 0 ; i < nbr_spheres; i++){
            Eigen::Vector3d position = last_sphere.center + vector*distance_between_spheres*(i+1);
            //double length_axon = axon_length(axon_to_grow);
            //double rad = radius_variation(axon_to_grow, length_axon, factor, beading_period, min_radius);
            double rad = last_sphere.radius + (sph.radius - last_sphere.radius)*(i+1)/(nbr_spheres+1);
            id_ = last_id+ 1;
            Sphere s(id_, sph.object_id, sph.object_type, position, rad, sph.branch_id, sph.parent_id);
            
            bool can_grow_ = canSpherebePlaced(s, true);
            if(can_grow_){
                last_id = s.id;
                axon_to_grow.add_sphere(s);
                //cout <<"sphere : " << s.id << " axon : "<< s.object_id << " can be placed as interpolated" << endl;
            }
        }
    }
    sph.id = last_id + 1;
    axon_to_grow.add_sphere(sph);
    //cout <<"sphere : " << sph.id << " axon : "<< sph.object_id << " can be placed as last" << endl;
    
}
void AxonGrowth::find_next_center(Sphere &s,  double dist_, const std::vector<Sphere> &spheres, const Eigen::Vector3d &target, const bool &is_axon)
{

    Eigen::Vector3d target_ = target;
    Eigen::Vector3d additional_vector = axon_to_grow.end - axon_to_grow.begin;
    additional_vector.normalize();
    additional_vector = additional_vector * 10;

    if ((target_-spheres[spheres.size() - 1].center).norm()< 10){
        target_ = target_ + additional_vector;
    }
    Eigen::Vector3d vector_to_target = target_ - spheres[spheres.size() - 1].center;
    vector_to_target = vector_to_target.normalized();
    Eigen::Vector3d vector = generate_random_point_on_sphere(std_dev);
    vector = apply_bias_toward_target(vector, vector_to_target);
    Eigen::Vector3d position = spheres[spheres.size() - 1].center + dist_ * vector.normalized();
    s.center = position;
}



