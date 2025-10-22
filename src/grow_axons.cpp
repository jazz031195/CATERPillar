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
                       const std::vector<Glial>* glial_pop1_,
                       const std::vector<Glial>* glial_pop2_,
                       const std::vector<Axon>* axons_,
                       const Eigen::Vector3d &extended_min_limits_,
                       const Eigen::Vector3d &extended_max_limits_,
                       const Eigen::Vector3d &min_limits_,
                       const Eigen::Vector3d &max_limits_,
                       const double &std_dev_,
                       const double &min_radius_)
    : CellGrowth(axons_, glial_pop1_, glial_pop2_,
                 extended_min_limits_, extended_max_limits_,
                 min_limits_, max_limits_,
                 std_dev_, min_radius_),
      axon_to_grow(axon_to_grow_)
{}

AxonGrowth::AxonGrowth(const AxonGrowth &other)
    : CellGrowth(other), // call base copy constructor
      axon_to_grow(other.axon_to_grow) {}



Eigen::Vector3d AxonGrowth::find_next_center_straight(const double distance, const std::vector<Sphere> &spheres)
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
                candidate.center= find_next_center_straight(distance, axon_to_grow.outer_spheres);
            } else {
                candidate.center= find_next_center(distance, axon_to_grow.outer_spheres, axon_to_grow.end);
            }
        } else {
            // If std_dev == 0, we skip the "straight vs. random" logic and always use find_next_center
            candidate.center= find_next_center(distance, axon_to_grow.outer_spheres, axon_to_grow.end);
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


Eigen::Vector3d AxonGrowth::find_next_center(const double dist_, 
                                             const std::vector<Sphere> &spheres, 
                                             const Eigen::Vector3d &target)
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

    
    if (spheres.size() > 2) {
        Eigen::Vector3d previous_vector = (last_center - spheres[spheres.size() - 2].center).normalized();
        double cos_angle = previous_vector.dot(biased_random_vector.normalized());
        cos_angle = std::max(-1.0, std::min(1.0, cos_angle)); // clamp to [-1, 1]

        double angle = acos(cos_angle);
        int nbr_tries = 0;
        double angle_limit = M_PI / 2; 
        
        while(angle > angle_limit && nbr_tries < 10) {
            biased_random_vector = apply_bias_toward_target(
                generate_random_point_on_sphere(std_dev), target_direction
            );

            cos_angle = previous_vector.dot(biased_random_vector.normalized());
            cos_angle = std::max(-1.0, std::min(1.0, cos_angle)); // clamp to [-1, 1]
            angle = acos(cos_angle);
            nbr_tries += 1;

        }
        
        if (angle > angle_limit) {
            biased_random_vector = previous_vector;
        }
        
        
    }
    

    return last_center + dist_ * biased_random_vector;
}

double AxonGrowth::RandomradiusVariation() 
{ 
    double prev_radius = axon_to_grow.outer_spheres[axon_to_grow.outer_spheres.size()-1].radius; 

    double prev_radius_clamped = prev_radius;

    if (prev_radius_clamped < axon_to_grow.radius-axon_to_grow.beading_amplitude*axon_to_grow.radius) 
    { 
        prev_radius_clamped = axon_to_grow.radius-axon_to_grow.beading_amplitude*axon_to_grow.radius; 
    } 
    else if (prev_radius_clamped > axon_to_grow.radius+axon_to_grow.beading_amplitude*axon_to_grow.radius) 
    { 
        prev_radius_clamped = axon_to_grow.radius+axon_to_grow.beading_amplitude*axon_to_grow.radius; 
    } 

    double standard_deviation = axon_to_grow.radius * axon_to_grow.beading_std; 

    std::normal_distribution<double> dis(prev_radius_clamped, standard_deviation); 

    double random_radius = dis(gen); 

    while (random_radius > prev_radius_clamped+3* standard_deviation || random_radius < prev_radius_clamped-3*standard_deviation) 
    { 
        random_radius = dis(gen); 
    } 
    if (random_radius < min_radius) 
    { 
        random_radius = min_radius; 
    } 
    return random_radius; 
}

bool AxonGrowth::shrinkRadius(const double &radius_to_shrink, const bool& axon_can_shrink, const int &factor)
{

    bool can_grow_;
    // find position that works for smallest radius
    can_grow_ = false;
    double rad = radius_to_shrink;
    double initial_rad = radius_to_shrink;
    bool can_grow_min_rad = AddOneSphere(min_radius, false, 0, factor);
    double intervals = (initial_rad) / 10;

    if(!axon_can_shrink){
        return false;
    }
    if (!can_grow_min_rad)
    {
        // cout << "minimum radius doesn't fit, min rad :"<< min_radius << endl;
        return false;
    }
    else
    {
        while (!can_grow_ && rad > min_radius+intervals)
        {
            rad -= intervals;
            can_grow_ = AddOneSphere(rad, true, 0, factor);
        }
        if (can_grow_)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}

void AxonGrowth::update_straight(bool can_grow_, int &grow_straight, int &straight_growths)
{

    if (can_grow_)
    {

        if (grow_straight == 1)
        {
            if (straight_growths >= axon_to_grow.undulation_factor) // if axon has been growing straight for a number of spheres in a row
            {
                grow_straight = 0; // set to false so that next step doesn't go straight
                straight_growths = 0;
            }
            else
            {
                straight_growths += 1;
            }
        }
        else
        {
            // if the sphere hadn't grown straight previously . set to straight for next "undulation_factor" spheres
            grow_straight = 1; // set to true
        }
    }
    else
    {
        if (grow_straight == 1) // if when growing straight it collides with environment
        {
            grow_straight = 0; // set to false so that next step doesn't go straight
            straight_growths = 0;
        }
    }
}

void AxonGrowth::growthThread(
    double& stuck_radius,
    int& stuck_index,
    int factor,
    bool axon_can_shrink
) {
    auto set_stuck = [&](bool done) {
        finished = done ? true : false;                  // explicit member
        stuck_radius   = done ? axon_to_grow.radius : -1.0;
        stuck_index    = done ? axon_to_grow.id     : -1;
    };
    int  grow_straight     = 0;
    int  straight_growths  = 0;

    std::size_t tries = 0;
    const std::size_t max_tries = 10000;  // large enough; growth logic should end earlier

    while (!finished && tries < max_tries) {
        ++tries;

        // Radius beading
        double varied_radius = axon_to_grow.radius;
        if (axon_to_grow.beading_amplitude > 0) {
            varied_radius = RandomradiusVariation();
        }

        // Try to place a sphere
        const bool grew = AddOneSphere(varied_radius, /*create_sphere=*/true, grow_straight, factor);

        // Growth process may have marked itself finished
        if (finished) {
            set_stuck(false);
            update_straight(grew, grow_straight, straight_growths);
            break;
        }

        if (grew) {
            set_stuck(false);
            update_straight(true, grow_straight, straight_growths);
            tries = 0;
            continue;
        }

        // Could not grow: try shrinking if allowed
        if (axon_can_shrink) {
            const bool shrank_and_grew = shrinkRadius(varied_radius, axon_can_shrink, factor);
            if (shrank_and_grew) {
                set_stuck(false);
                update_straight(true, grow_straight, straight_growths);
                continue;
            }
            // fall through to retry/destroy
        }

        // No shrink or shrink failed → retry a few times, then give up
        if (axon_to_grow.growth_attempts < 3) {
            axon_to_grow.keep_one_sphere();  // retry from same position
            set_stuck(false);
        } else {
            axon_to_grow.destroy();
            set_stuck(true);
            break;
        }

        update_straight(false, grow_straight, straight_growths);
    }

    if (tries >= max_tries && !finished) {
        // safety stop
        axon_to_grow.destroy();
        set_stuck(true);
    }

    if (finished && !axon_to_grow.outer_spheres.empty()) {
        if (axon_to_grow.outer_spheres.size() < 10) {

            axon_to_grow.destroy();
            set_stuck(true);
        }
    }

}

