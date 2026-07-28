#include "CaterpillarGrowth.h"
#include "grow_blood_vessels.h"
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


BloodVesselGrowth::~BloodVesselGrowth() {}

BloodVesselGrowth::BloodVesselGrowth(Blood_Vessel &bv_to_grow_,
                       const SphereGrid* sphere_grid_,
                       const Eigen::Vector3d &extended_min_limits_,
                       const Eigen::Vector3d &extended_max_limits_,
                       const Eigen::Vector3d &min_limits_,
                       const Eigen::Vector3d &max_limits_,
                       const double &epsilon_,
                       const double &min_radius_,
                       const double &capillary_radius_,
                       const int &max_generations_)
    : CellGrowth(sphere_grid_,
                 extended_min_limits_, extended_max_limits_,
                 min_limits_, max_limits_,
                 epsilon_, min_radius_),
      bv_to_grow(bv_to_grow_),
      capillary_radius(capillary_radius_),
      max_generations(max_generations_)
{}

BloodVesselGrowth::BloodVesselGrowth(const BloodVesselGrowth &other)
    : CellGrowth(other), // call base copy constructor
      bv_to_grow(other.bv_to_grow),
      capillary_radius(other.capillary_radius),
      max_generations(other.max_generations) {}


// =====================================================================================
// Main vessel growth (branch 0)
// =====================================================================================

Eigen::Vector3d BloodVesselGrowth::find_next_center_straight(const double distance, const std::vector<Sphere> &spheres)
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



bool BloodVesselGrowth::AddOneSphere(double radius_, bool create_sphere, int grow_straight, const int &factor)
{
    std::vector<Sphere> &trunk = bv_to_grow.ramification_spheres[0];

    // Basic validation
    if (trunk.empty()) {
        std::cerr << "EMPTY BLOOD VESSEL!" << std::endl;
        assert(false); // or return false;
    }

    // Sized off extended_min/max_limits (the blood-vessel growth voxel, which may be
    // padded larger than the real one) rather than min/max_limits, so this safety cap
    // can't truncate a trunk before it reaches its actual (possibly farther) target face.
    if (trunk.size() > (extended_max_limits-extended_min_limits).norm()*factor/(bv_to_grow.radius)*100) {
        finished = true;
        return false; // Vessel has grown too long
    }

    assert(bv_to_grow.growth_axis >= 0 && bv_to_grow.growth_axis < 3);

    bool is_allowed_to_stop_early = false;

    // If the last sphere's center is beyond extended_max_limits, the vessel is considered fully grown
    Sphere last_sphere = trunk.back();
    if (last_sphere.center[bv_to_grow.growth_axis] >= extended_max_limits[bv_to_grow.growth_axis] && !is_allowed_to_stop_early) {
        finished = true;

        return true; // Vessel is done
    }
    else if (!check_borders(min_limits, max_limits, last_sphere.center, last_sphere.radius) && is_allowed_to_stop_early) {
        finished = true;
        return true; // Vessel is done
    }


    // Prepare
    double max_radius_ = std::max(radius_, last_sphere.radius);

    double distance = max_radius_;
    int threshold_tries = (epsilon == 0.0) ? 1 : 100; // tries 100 times to place sphere

    // New sphere to attempt placing
    Sphere s(trunk.size() + factor,
             bv_to_grow.id,
             /*object_type=*/blood_constant,
             bv_to_grow.begin,
             radius_,
             /*branch_id=*/0);

    bool can_grow_ = false;
    int tries = 0;

    // Helper lambda: tries to place sphere (find next center, check if inside, collision-check).
    // Returns true if placed successfully, false otherwise.
    auto attemptPlacement = [&](Sphere &candidate, int triesCount) -> bool {
        // 1) Find next center
        if (epsilon != 0.0) {

            // If grow_straight is set, use the straight function, else the normal function
            if (grow_straight == 1) {
                candidate.center= find_next_center_straight(distance, trunk);
            } else {
                candidate.center= find_next_center(distance, trunk, bv_to_grow.end);

            }
        } else {
            // If epsilon == 0, we skip the "straight vs. random" logic and always use find_next_center
            candidate.center= find_next_center(distance, trunk, bv_to_grow.end);
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
        // Attempt the main placement
        bool success = attemptPlacement(s, tries);
        if (success) {
            can_grow_ = true;
            break;
        }

        // If still not placed, and we used "grow_straight == 1", then attempt a fallback (non-straight) approach
        if (!can_grow_ && (epsilon != 0.0) && (grow_straight == 1)) {
            grow_straight = 0; // Fallback to normal growth
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
        // Refresh volume (both real-voxel and blood-vessel-voxel figures) as long as
        // the newest sphere is still within the growth region -- gated on the wider
        // extended_min/max_limits (not min/max_limits) so this keeps running for
        // spheres grown past the real voxel into the padded region too.
        Sphere newly_added = trunk.back(); // s with final coords
        if (check_borders(extended_min_limits, extended_max_limits, newly_added.center, newly_added.radius)) {
            bv_to_grow.update_Volume(factor, min_limits, max_limits, extended_min_limits, extended_max_limits);
        }
    }

    // If we reach the edge of voxel after adding the new sphere
    Sphere current_last = trunk.back();

    if (current_last.center[bv_to_grow.growth_axis] + current_last.radius > extended_max_limits[bv_to_grow.growth_axis]) {

        finished = true;
    }
    else if (!check_borders(min_limits, max_limits, current_last.center, current_last.radius) && is_allowed_to_stop_early) {

        finished = true;
    }
    return true;
}


void BloodVesselGrowth::add_spheres(Sphere &sph, const Sphere &last_sphere, const int &factor){

    // nbr of spheres to add in between
    int nbr_spheres = factor - 1;

    if (factor > 1){
        // distance between two consecutive spheres
        double distance = (sph.center - last_sphere.center).norm();
        Eigen::Vector3d vector = (sph.center - last_sphere.center).normalized();
        double distance_between_spheres = distance/(nbr_spheres+1);
        for (int i = 0 ; i < nbr_spheres; i++){
            Eigen::Vector3d position = last_sphere.center + vector*distance_between_spheres*(i+1);
            double rad = last_sphere.radius + (sph.radius - last_sphere.radius)*(i+1)/(nbr_spheres+1);
            Sphere s(bv_to_grow.next_sphere_id++, sph.object_id, sph.object_type, position, rad, sph.branch_id, sph.parent_id);

            bool can_grow_ = canSpherebePlaced(s);
            if(can_grow_){
                bv_to_grow.add_sphere(s);
            }
        }
    }
    sph.id = bv_to_grow.next_sphere_id++;
    bv_to_grow.add_sphere(sph);

}


Eigen::Vector3d BloodVesselGrowth::find_next_center(const double dist_,
                                             const std::vector<Sphere> &spheres,
                                             const Eigen::Vector3d &target)
{
    Eigen::Vector3d last_center = spheres.back().center;
    Eigen::Vector3d target_direction = (target - last_center).normalized();

    Eigen::Vector3d additional_vector = (bv_to_grow.end - bv_to_grow.begin).normalized() * 10.0;

    if ((target - last_center).norm() < 10) {
        target_direction = (target + additional_vector - last_center).normalized();
    }

    Eigen::Vector3d biased_random_vector = apply_bias_toward_target(
        generate_random_point_on_sphere(epsilon), target_direction
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
                generate_random_point_on_sphere(epsilon), target_direction
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

double BloodVesselGrowth::RandomradiusVariation()
{
    std::vector<Sphere> &trunk = bv_to_grow.ramification_spheres[0];
    double prev_radius = trunk.back().radius;

    double prev_radius_clamped = prev_radius;

    if (prev_radius_clamped < bv_to_grow.radius-bv_to_grow.beading_amplitude*bv_to_grow.radius)
    {
        prev_radius_clamped = bv_to_grow.radius-bv_to_grow.beading_amplitude*bv_to_grow.radius;
    }
    else if (prev_radius_clamped > bv_to_grow.radius+bv_to_grow.beading_amplitude*bv_to_grow.radius)
    {
        prev_radius_clamped = bv_to_grow.radius+bv_to_grow.beading_amplitude*bv_to_grow.radius;
    }

    double standard_deviation = bv_to_grow.radius * bv_to_grow.beading_std;

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

bool BloodVesselGrowth::shrinkRadius(const double &radius_to_shrink, const bool& bv_can_shrink, const int &factor)
{

    bool can_grow_;
    // find position that works for smallest radius
    can_grow_ = false;
    double rad = radius_to_shrink;
    double initial_rad = radius_to_shrink;
    bool can_grow_min_rad = AddOneSphere(min_radius, false, 0, factor);
    double intervals = (initial_rad) / 10;

    if(!bv_can_shrink){
        return false;
    }
    if (!can_grow_min_rad)
    {
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

void BloodVesselGrowth::update_straight(bool can_grow_, int &grow_straight, int &straight_growths)
{

    if (can_grow_)
    {

        if (grow_straight == 1)
        {
            if (straight_growths >= bv_to_grow.undulation_factor) // if vessel has been growing straight for a number of spheres in a row
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

void BloodVesselGrowth::growthThread(
    double& stuck_radius,
    int& stuck_index,
    int factor,
    bool bv_can_shrink
) {
    auto set_stuck = [&](bool done) {
        finished = done ? true : false;                  // explicit member
        stuck_radius   = done ? bv_to_grow.radius : -1.0;
        stuck_index    = done ? bv_to_grow.id     : -1;
    };
    int  grow_straight     = 0;
    int  straight_growths  = 0;

    std::size_t tries = 0;
    const std::size_t max_tries = 1000000;  // large enough; growth logic should end earlier
    int total_nbr_growth_attempts = 10;

    while (!finished && tries < max_tries) {
        ++tries;

        // Radius beading
        double varied_radius = bv_to_grow.radius;
        if (bv_to_grow.beading_amplitude > 0) {
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
        if (bv_can_shrink) {
            const bool shrank_and_grew = shrinkRadius(varied_radius, bv_can_shrink, factor);
            if (shrank_and_grew) {
                set_stuck(false);
                update_straight(true, grow_straight, straight_growths);
                continue;
            }
            // fall through to retry/destroy
        }

        // No shrink or shrink failed → retry a few times, then give up
        if (bv_to_grow.growth_attempts < total_nbr_growth_attempts) {
            bv_to_grow.keep_one_sphere();  // retry from same position
            set_stuck(false);
        } else {
            bv_to_grow.destroy();
            set_stuck(true);
            break;
        }

        update_straight(false, grow_straight, straight_growths);
    }

    if (tries >= max_tries && !finished) {
        // safety stop
        bv_to_grow.destroy();
        set_stuck(true);
    }

    if (finished && !bv_to_grow.ramification_spheres.empty() && !bv_to_grow.ramification_spheres[0].empty()) {
        if (bv_to_grow.ramification_spheres[0].size() < 10) {

            bv_to_grow.destroy();
            set_stuck(true);
        }
    }

}


// =====================================================================================
// Branch growth (branch_id >= 1) - mirrors GlialCellGrowth's secondary-branch machinery
// =====================================================================================

void BloodVesselGrowth::add_spheres(Sphere &sph, const Sphere &last_sphere, const bool &check_collision_with_branches, const int &factor, const int &index_ram_spheres, int extra_excluded_branch_id){

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
            id_ = bv_to_grow.next_sphere_id++;
            Sphere s(id_, last_sphere.object_id, last_sphere.object_type, position, rad, last_sphere.branch_id, sph.parent_id);
            bool can_grow_ = canSpherebePlaced(s, check_collision_with_branches, extra_excluded_branch_id);
            if(can_grow_){
                bv_to_grow.ramification_spheres[index_ram_spheres].push_back(s);
            }
        }
    }
    sph.id = bv_to_grow.next_sphere_id++;
    bv_to_grow.ramification_spheres[index_ram_spheres].push_back(sph);
}


bool BloodVesselGrowth::AddOneSphere(const double &radius_, const bool &create_sphere, int &grow_straight, const int &i, const bool &check_collision_with_branches, const int &parent, const int &factor, int extra_excluded_branch_id)
{

    if (bv_to_grow.ramification_spheres[i].size() > 1e10)
    {
        return false;
    }

    double distance;
    Eigen::Vector3d destination = bv_to_grow.attractors[i];
    if (bv_to_grow.ramification_spheres.size() == 0 ){
        assert(0);
    }
    Sphere last_sphere = bv_to_grow.ramification_spheres[i][bv_to_grow.ramification_spheres[i].size() - 1];
    // the distance between two consecutive spheres is the maximum radius / 2
    double max_radius_ = max(radius_, last_sphere.radius);

    distance = (max_radius_);

    bool can_grow_ = false;

    int id_ = last_sphere.id + factor;
    Sphere s(id_, bv_to_grow.id, blood_constant, {0,0,0}, radius_, i);

    int tries = 0;
    int threshold_tries = 100;

    while (!can_grow_ && tries < threshold_tries){

        find_next_center(s, distance, bv_to_grow.ramification_spheres[i], destination);
        can_grow_ = canSpherebePlaced(s, check_collision_with_branches, extra_excluded_branch_id);

        tries += 1;
    }

    if (can_grow_ )
    {
        if (create_sphere){
            s.parent_id = parent;
            add_spheres(s, last_sphere, check_collision_with_branches, factor, i, extra_excluded_branch_id);
        }
        // if is not in inside voxel
        if (!check_borders(extended_min_limits, extended_max_limits, s.center, s.radius)){
            finished = true;
            return false;
        }
        else{
            return true;
        }
    }
    else // collides and >1000 tries
    {
        return false;
    }

}

void BloodVesselGrowth::find_next_center(Sphere &s, double dist_, const std::vector<Sphere> &spheres, const Eigen::Vector3d &target)
{

    Eigen::Vector3d target_ = target;
    Eigen::Vector3d vector_to_target = target_ - spheres[spheres.size() - 1].center;
    vector_to_target = vector_to_target.normalized();
    double std_ = epsilon;
    Eigen::Vector3d vector = generate_random_point_on_sphere(std_);
    vector = apply_bias_toward_target(vector, vector_to_target);
    Eigen::Vector3d position = spheres[spheres.size() - 1].center + dist_ * vector.normalized();
    s.center = position;
}


bool BloodVesselGrowth::GenerateFirstSphereinProcess(Sphere &first_sphere, Eigen::Vector3d &attractor, const double &radius, const Sphere &sphere_to_emerge_from, const Eigen::Vector3d &vector_to_prev_center, const int &nbr_spheres, const int &nbr_spheres_between, const int &vessel_id, const int &branch_id) {

    bool stop = false;
    int tries_ = 0;
    int max_nbr_tries = 10;
    attractor = Eigen::Vector3d(0, 0, 0);

    while (!stop && tries_ < max_nbr_tries) {
        Eigen::Vector3d vector = {0, 0, 0};
        Eigen::Vector3d point = {0, 0, 0};
        // vessels have no primary/secondary distinction: always avoid folding back
        // toward the branch it emerges from, mirroring Glial's secondary-branch case.
        sphere_to_emerge_from.getPointOnSphereSurface(point, vector, vector_to_prev_center, false);
        // "+10" past extended_max_limits (not max_limits): keeps this attractor point
        // meaningfully outside the actual (possibly padded) growth voxel.
        attractor = point + (extended_max_limits[0] + 10) * vector;
        first_sphere = Sphere(bv_to_grow.next_sphere_id++, vessel_id, blood_constant, point, radius, branch_id, sphere_to_emerge_from.id);
        // Exempt only this branch's own id and its direct parent (sphere_to_emerge_from's
        // branch) -- not the whole vessel -- so a legitimate touch at the attachment point
        // doesn't get flagged, while every *other* branch is still checked from sphere 1.
        if (canSpherebePlaced(first_sphere, /*check_collision_with_branches=*/ true, sphere_to_emerge_from.branch_id)) {
            stop = true;
            // check boundaries
            bool is_inside_voxel = check_borders(extended_min_limits, extended_max_limits, first_sphere.center, first_sphere.radius);
            if (!is_inside_voxel) {
                stop = false;
                tries_++;
                if (tries_ == max_nbr_tries) return false;
            }
        }
        else {
            tries_++;
            if (tries_ == max_nbr_tries) return false;
        }
    }
    return true;
}


std::vector<Sphere> BloodVesselGrowth::addIntermediateSpheres(const Sphere &random_sphere, const Sphere &first_sphere, const int &branch_nbr, const int &nbr_spheres, const int &nbr_spheres_between, const std::function<double(double)> &compute_radius, const double &t_start, const double &t_end) {

    std::vector<Sphere> intermediate_spheres;

    // Direction vector for sphere placement
    Eigen::Vector3d direction = (first_sphere.center - random_sphere.center).normalized();
    double total_distance = (first_sphere.center - random_sphere.center).norm();
    double distance_between_spheres = total_distance / (nbr_spheres_between + 1);

    // Add intermediate spheres
    for (int i = 0; i < nbr_spheres_between; ++i) {
        Eigen::Vector3d position = random_sphere.center + direction * distance_between_spheres * (i + 1);
        double t = t_start + (t_end - t_start) * (i + 1) / (nbr_spheres_between + 1);
        double rad = compute_radius(t);

        Sphere next(
            bv_to_grow.next_sphere_id++, first_sphere.object_id, first_sphere.object_type, position, rad, branch_nbr,
            random_sphere.id);

        // Check if the sphere can be placed. Exempt only this branch's own id and its
        // direct parent (random_sphere's branch) -- not the whole vessel -- same reasoning
        // as GenerateFirstSphereinProcess.
        if (canSpherebePlaced(next, /*check_collision_with_branches=*/ true, random_sphere.branch_id)) {
            intermediate_spheres.emplace_back(next);
        }
    }

    // Add the initial branching sphere to the list
    intermediate_spheres.push_back(first_sphere);

    return intermediate_spheres;
}

bool BloodVesselGrowth::growBranch(int &nbr_spheres, const int &factor) {

    if (bv_to_grow.ramification_spheres.empty()) {
        cout << "No branches in blood vessel" << endl;
        return false;
    }
    if (nbr_spheres <= 0) {
        cout << "No spheres in blood vessel : " <<  bv_to_grow.id << endl;
        return false;
    }

    int nbr_branches = bv_to_grow.ramification_spheres.size();

    std::mt19937 rng(std::random_device{}());
    finished = false;

    // finding source of branching
    int nbr_spheres_between = factor - 1;

    // Weight each candidate branch by its own arc length, so longer branches
    // (more likely to be mature/thick vessels in practice) are more likely to
    // sprout a new sub-branch than short ones. Branches with <=1 sphere get
    // zero weight (need at least 2 spheres so vector_to_prev_sphere below is
    // well-defined) and are never selected.
    std::vector<double> branch_weights(nbr_branches, 0.0);
    bool any_candidate = false;
    for (int b = 0; b < nbr_branches; ++b) {
        const auto &branch = bv_to_grow.ramification_spheres[b];
        if (branch.size() <= 1) {
            continue;
        }
        double length = 0.0;
        for (size_t k = 1; k < branch.size(); ++k) {
            length += (branch[k].center - branch[k - 1].center).norm();
        }
        branch_weights[b] = length;
        any_candidate = true;
    }
    if (!any_candidate) {
        cout << "Could not find a branch with enough spheres for blood vessel : " << bv_to_grow.id << endl;
        return false;
    }
    std::discrete_distribution<int> branch_dist(branch_weights.begin(), branch_weights.end());

    // Cap how many generations of capillary branching are allowed from the arteriole:
    // a new branch's generation is its chosen parent's generation + 1, and that must
    // not exceed max_generations. Resample the candidate parent up to 20 times looking
    // for one that keeps the new branch within the limit; if none is found, fall back
    // to sprouting directly off the arteriole (branch 0, generation 0) so growth still
    // makes progress instead of stalling.
    const int max_branch_selection_tries = 20;
    int random_branch = -1;
    for (int try_i = 0; try_i < max_branch_selection_tries; ++try_i) {
        int candidate = branch_dist(rng);
        int candidate_generation = bv_to_grow.branch_generation[candidate];
        if (candidate_generation + 1 <= max_generations) {
            random_branch = candidate;
            break;
        }
    }
    if (random_branch < 0) {
        random_branch = 0;
    }

    int size = bv_to_grow.ramification_spheres[random_branch].size();
    int random_sphere_ind = 1 + rand() % (size - 1);
    Sphere random_sphere = bv_to_grow.ramification_spheres[random_branch][random_sphere_ind];

    if (bv_to_grow.lengths_branches.size() <= random_branch) {
        cout << "Error: lengths_branches vector is not large enough." << endl;
        assert(0);
    }
    double old_length = bv_to_grow.lengths_branches[random_branch][random_sphere_ind];

    // Capillaries no longer target a random (Gaussian) length: they always aim to
    // cross a voxel face, so a length of 2x the voxel edge length is generous
    // enough to guarantee that (the longest possible straight-line distance
    // inside the box is its diagonal, sqrt(3)*edge < 2*edge) regardless of
    // where along the box the branch starts or which direction it grows in.
    // If a branch instead gets stuck (collision) before ever reaching a face, it
    // is discarded below.
    // Uses extended_min/max_limits (the blood-vessel growth voxel, which may be
    // padded larger than the real one) so the budget actually covers the box
    // branches are targeting -- min/max_limits would under-budget a branch aimed
    // at the farther, padded face.
    double voxel_edge_length = extended_max_limits[0] - extended_min_limits[0];
    double length_to_grow = 2.0 * voxel_edge_length;
    int nbr_non_checked_spheres = factor*3 ;

    Eigen::Vector3d vector_to_prev_sphere = (random_sphere.center - bv_to_grow.ramification_spheres[random_branch][random_sphere_ind - 1].center).normalized();
    Sphere first_sphere;
    Eigen::Vector3d attractor = Eigen::Vector3d(0, 0, 0);
    // Capillaries hold one fixed radius for their entire length, unrelated to the
    // parent (arteriole/capillary) sphere's own local radius at the attachment
    // point -- no decay, and never rescaled by Murray's law (see
    // Blood_Vessel::enforceMurraysLaw, which only ever thins the arteriole side).
    double initial_radius = capillary_radius;
    bool first_sphere_created = GenerateFirstSphereinProcess(first_sphere, attractor, initial_radius, random_sphere, vector_to_prev_sphere, nbr_spheres, nbr_spheres_between, bv_to_grow.id, nbr_branches);

    if (!first_sphere_created) {
        return false;
    }

    auto compute_radius = [&](double t) {
        return capillary_radius;
    };

    std::vector<Sphere> vector_first_spheres;
    if (factor > 1) {
        // add spheres between the first and the last
        vector_first_spheres = addIntermediateSpheres(random_sphere, first_sphere, nbr_branches, nbr_spheres, nbr_spheres_between, compute_radius, 0.0, 0.0);
    } else {
        vector_first_spheres = {first_sphere};
    }

    Blood_Vessel old_bv = bv_to_grow;
    int current_branch = bv_to_grow.ramification_spheres.size();
    bv_to_grow.ramification_spheres.resize(current_branch + 1);
    bv_to_grow.lengths_branches.resize(current_branch + 1);
    bv_to_grow.attractors.resize(current_branch + 1);
    bv_to_grow.children_branches.resize(current_branch + 1);
    bv_to_grow.branch_generation.resize(current_branch + 1);
    bv_to_grow.ramification_spheres[current_branch] = vector_first_spheres;
    bv_to_grow.lengths_branches[current_branch] = std::vector<double>(vector_first_spheres.size(), old_length);
    bv_to_grow.attractors[current_branch] = attractor;
    bv_to_grow.children_branches[random_branch].push_back(current_branch);
    bv_to_grow.branch_generation[current_branch] = bv_to_grow.branch_generation[random_branch] + 1;

    Eigen::Vector3d prev_pos = first_sphere.center;
    double distance = initial_radius;
    double total_distance = initial_radius + old_length;
    bool can_grow = true;
    int stop_criteria = 0;

    int parent = random_sphere.id;

    while (can_grow && stop_criteria < 1e5) {

        double R_ = compute_radius(distance);

        // No radius-decay stopping condition (R_ is always capillary_radius); only
        // the length budget and AddOneSphere's own collision/wall handling end growth.
        if (distance >= length_to_grow) {
            can_grow = false;

        }
        else {
            int grow_straight = 0;
            // Always check collisions against every OTHER branch of this vessel except
            // this one's own (current_branch, exempted unconditionally). The direct
            // parent (random_branch) is exempted too, but only while still within the
            // near-attachment window (first nbr_non_checked_spheres spheres) -- past
            // that, the branch must avoid its own parent same as any other branch, so it
            // can't drift back into it later in its growth. This replaces the old
            // distance-based grace period that exempted the *whole vessel* (any sibling,
            // not just the parent) for that same window.
            int extra_excluded_branch_id = (bv_to_grow.ramification_spheres.back().size() <= nbr_non_checked_spheres) ? random_branch : -1;
            can_grow = AddOneSphere(R_, true, grow_straight, current_branch, /*check_collision_with_branches=*/true, parent, factor, extra_excluded_branch_id);
            if (can_grow) {
                double segment = (prev_pos - bv_to_grow.ramification_spheres.back().back().center).norm();
                distance += segment;
                total_distance += segment;
                int nbr_spheres_in_branch_before = bv_to_grow.lengths_branches[current_branch].size();
                int nbr_spheres_in_branch_after = bv_to_grow.ramification_spheres.back().size();
                for (int i = nbr_spheres_in_branch_before; i < nbr_spheres_in_branch_after; i++) {
                    bv_to_grow.lengths_branches[current_branch].push_back(total_distance);
                }
                prev_pos = bv_to_grow.ramification_spheres.back().back().center;
                parent = bv_to_grow.ramification_spheres.back().back().id;
            }
        }

        stop_criteria++;
    }

    // finished is only set true by AddOneSphere when the branch's last placed
    // sphere fell outside the voxel box, i.e. the branch actually reached a
    // voxel plane. Any other reason the loop stopped (hit the length budget, or
    // growth got stuck against a collision) leaves finished false, and the
    // branch must be discarded.
    if (finished && bv_to_grow.ramification_spheres.back().size() > nbr_non_checked_spheres) {
        nbr_spheres += bv_to_grow.ramification_spheres.back().size();

        return true;
    }
    else{
        bv_to_grow = std::move(old_bv);
        return false;
    }
}
