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
      axon_to_grow(axon_to_grow_), 
      axons(axons_)
{}

AxonGrowth::AxonGrowth(const AxonGrowth &other)
    : CellGrowth(other), // call base copy constructor
      axon_to_grow(other.axon_to_grow),
      axons(other.axons) {}


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

            if (canSpherebePlaced(new_sphere) && check_borders(extended_min_limits, extended_max_limits, new_sphere.center, new_sphere.radius)) {
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

Eigen::Vector3d AxonGrowth::readapt_sphere_position(const Sphere &s, const Axon &neighbour_axon, bool can_readapt){

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
    const Sphere& prev_sphere = axon_to_grow.outer_spheres[axon_to_grow.outer_spheres.size() - 1];
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

bool AxonGrowth::checkAxonsOverlap(Sphere &sph){
    const int max_iterations = 10;
    int iteration = 0;
    bool moved = true;

    while (moved && iteration < max_iterations) {
        moved = false;
        for (const auto& axon : axons) {
            if (!(axon.id == sph.object_id && sph.object_type == 0)) {
                if (axon.isSphereInsideAxon(sph)) {
                    bool can_readapt = false;
                    Eigen::Vector3d new_center = readapt_sphere_position(sph, axon, can_readapt);
                    if (can_readapt) {
                        sph.center = new_center;
                        moved = true;  // we moved, need to recheck all
                        break;         // break to restart loop from beginning
                    } else {
                        return false;  // can't resolve
                    }
                }
            }
        }
        iteration++;
    }
    if (iteration >= max_iterations){
        return false; // if we moved, it means we are maybe colliding
    }
    else{
        return true;  // if we get here, it's placed safely
    }
}

bool AxonGrowth::canSpherebePlaced(Sphere &sph){

    bool axons_check = checkAxonsOverlap(sph);

    if (!axons_check){
        return false;
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


    return true;
} 

bool AxonGrowth::isSphereCollidingwithAxon(Axon ax, Sphere sph){

    return ax.isSphereInsideAxon(sph);
}


Eigen::Vector3d AxonGrowth::find_next_center_straight(const double distance, const Sphere &s, const std::vector<Sphere> &spheres)
{
    if (spheres.size() < 2){
        assert(0);
    }
    Eigen::Vector3d last_center = spheres[spheres.size() - 1].center;
    Eigen::Vector3d before_last_center = spheres[spheres.size() - 2].center;
    Eigen::Vector3d straight_vector = (last_center - before_last_center).normalized();
    straight_vector *= distance;
    Eigen::Vector3d new_center = last_center + straight_vector;
    return new_center;

}


bool AxonGrowth::AddOneSphere(double radius_, bool create_sphere, int grow_straight, const int &factor)
{
    // Basic validation
    if (axon_to_grow.outer_spheres.empty()) {
        std::cerr << "EMPTY AXON!" << std::endl;
        assert(false); // or return false;
    }

    if (axon_to_grow.outer_spheres.size() > (max_limits-min_limits).norm()*factor/(axon_to_grow.radius)*100) {
        finished = true;
        return false; // Axon has grown too long
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
    int threshold_tries = (std_dev == 0.0) ? 1 : 100; // tries 100 times to place sphere

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
    auto attemptPlacement = [&](Sphere &candidate, int triesCount) -> bool {
        // 1) Find next center
        if (std_dev != 0.0) {
            // If grow_straight is set, use the straight function, else the normal function
            if (grow_straight == 1) {
                candidate.center= find_next_center_straight(distance, candidate, axon_to_grow.outer_spheres);
            } else {
                candidate.center= find_next_center(candidate, distance, axon_to_grow.outer_spheres, axon_to_grow.end, true);
            }
        } else {
            // If std_dev == 0, we skip the "straight vs. random" logic and always use find_next_center
            candidate.center= find_next_center(candidate, distance, axon_to_grow.outer_spheres, axon_to_grow.end, true);
        }

        // 2) Check if inside
        if (!check_borders(extended_min_limits, extended_max_limits, candidate.center, candidate.radius) && !is_allowed_to_stop_early) {
            return false;
        }


        bool canPlace = canSpherebePlaced(candidate);

        return canPlace;
        

    };
 

    while (!can_grow_ && tries < threshold_tries)
    {

        /*
        if (tries > 0.75*threshold_tries) {
            with_push = true;
        }
        */
        //cout <<"axon : "<< axon_to_grow.id <<" with push : " << with_push << endl;
        // Attempt the main placement
        bool success = attemptPlacement(s, tries);
        if (success) {
            can_grow_ = true;
            break;
        }
        
        // If still not placed, and we used "grow_straight == 1", then attempt a fallback (non-straight) approach
        if (!can_grow_ && (std_dev != 0.0) && (grow_straight == 1)) {
            grow_straight == 0; // Fallback to normal growth
            if (attemptPlacement(s, tries)) {
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
            
            bool can_grow_ = canSpherebePlaced(s);
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

Eigen::Vector3d AxonGrowth::find_closest_neighbour(const Sphere &sphere){
    double min_distance = std::numeric_limits<double>::max();
    Eigen::Vector3d centre_closest_sphere = Eigen::Vector3d::Zero();
    bool found = false;

    for (const auto &axon : axons) {
        if (axon.isNearAxon(sphere.center, 2 * sphere.radius + barrier_tickness)) {
            if (axon.id != sphere.object_id) {
                for (const auto &sph : axon.outer_spheres) {
                    double distance = (sphere.center - sph.center).norm() - (sphere.radius + sph.radius);
                    if (distance < min_distance && distance < 2 * sphere.radius + barrier_tickness) {
                        min_distance = distance;
                        centre_closest_sphere = sph.center;
                        found = true;
                    }
                }
            }
        }
    }

    if (found) {
        return (centre_closest_sphere - sphere.center);
    } else {
        return Eigen::Vector3d::Zero();  // explicitly handle no-neighbor-found case
    }
}

Eigen::Vector3d AxonGrowth::find_next_center(const Sphere &s, const double dist_, 
                                             const std::vector<Sphere> &spheres, 
                                             const Eigen::Vector3d &target, const bool &is_axon)
{
    Eigen::Vector3d last_center = spheres.back().center;
    Eigen::Vector3d target_direction = (target - last_center).normalized();

    Eigen::Vector3d additional_vector = (axon_to_grow.end - axon_to_grow.begin).normalized() * 10.0;

    if ((target - last_center).norm() < 10) {
        target_direction = (target + additional_vector - last_center).normalized();
    }

    Eigen::Vector3d biased_random_vector = apply_bias_toward_target(
        generate_random_point_on_sphere(std_dev), target_direction
    );

    Eigen::Vector3d neighbor_vector = find_closest_neighbour(spheres.back());
    bool has_valid_neighbor = neighbor_vector.norm() > 1e-6;

    Eigen::Vector3d combined_vector;
    if (has_valid_neighbor) {
        // Adjust weighting explicitly: stronger attraction to target
        combined_vector = (0.9 * biased_random_vector + 0.1 * neighbor_vector.normalized()).normalized();
    } else {
        combined_vector = biased_random_vector.normalized();
    }

    if (spheres.size() > 2) {
        Eigen::Vector3d previous_vector = (last_center - spheres[spheres.size() - 2].center).normalized();
        double cos_angle = previous_vector.dot(combined_vector.normalized());
        cos_angle = std::max(-1.0, std::min(1.0, cos_angle)); // clamp to [-1, 1]

        double angle = acos(cos_angle);
        int nbr_tries = 0;

        while(angle > M_PI / 4 && nbr_tries < 10) {
            biased_random_vector = apply_bias_toward_target(
                generate_random_point_on_sphere(std_dev), target_direction
            );
            if (has_valid_neighbor) {
                combined_vector = (0.9 * biased_random_vector + 0.1 * neighbor_vector.normalized()).normalized();
            } else {
                combined_vector = biased_random_vector.normalized();
            }
            cos_angle = previous_vector.dot(combined_vector.normalized());
            cos_angle = std::max(-1.0, std::min(1.0, cos_angle)); // clamp to [-1, 1]
            angle = acos(cos_angle);
            nbr_tries += 1;
        }

        if (angle > M_PI / 4) {
            combined_vector = previous_vector;
        }
    }

    return last_center + dist_ * combined_vector;
}