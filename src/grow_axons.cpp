#include "CaterpillarGrowth.h"
#include "grow_axons.h"
#include <algorithm> // std::sort
#include <random>
#include <chrono>
#include <future>
#include <limits>
#include <cstdlib>
#include "Eigen/Dense"
#include <thread>
#include "threads.h"

using namespace std;
using namespace Eigen;
using namespace std::chrono;


// Spherical linear interpolation between two unit direction vectors, at
// parameter t in [0, 1] -- gives a constant-angular-rate rotation from v0
// to v1, used to spread a turn cycle's direction change smoothly across
// several spheres instead of jumping to the new heading in one step (see
// find_next_center/find_next_center_straight below).
static Eigen::Vector3d slerp_direction(const Eigen::Vector3d &v0, const Eigen::Vector3d &v1, double t)
{
    double dot = std::max(-1.0, std::min(1.0, v0.dot(v1)));
    double theta = std::acos(dot);
    if (theta < 1e-6) {
        // Nearly identical directions -- a linear blend is fine here and
        // avoids dividing by sin(theta) ~= 0.
        Eigen::Vector3d blended = v0 + t * (v1 - v0);
        double norm = blended.norm();
        return (norm > 1e-9) ? (blended / norm) : v0;
    }
    double sin_theta = std::sin(theta);
    double w0 = std::sin((1.0 - t) * theta) / sin_theta;
    double w1 = std::sin(t * theta) / sin_theta;
    return (w0 * v0 + w1 * v1).normalized();
}


AxonGrowth::~AxonGrowth() {}

AxonGrowth::AxonGrowth(Axon &axon_to_grow_,
                       const SphereGrid* sphere_grid_,
                       const Eigen::Vector3d &extended_min_limits_,
                       const Eigen::Vector3d &extended_max_limits_,
                       const Eigen::Vector3d &min_limits_,
                       const Eigen::Vector3d &max_limits_,
                       const double &epsilon_,
                       const double &min_radius_)
    : CellGrowth(sphere_grid_,
                 extended_min_limits_, extended_max_limits_,
                 min_limits_, max_limits_,
                 epsilon_, min_radius_),
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

    // Gradually curve from this turn cycle's starting heading toward its
    // freshly-sampled target (both captured in find_next_center, below,
    // when the target was sampled) across the whole cycle, rather than
    // having already jumped fully to the target direction and now
    // extrapolating it frozen for the rest of the cycle -- t ramps from
    // just above 0 (first step) to 1.0 (last step, fully at the target).
    double t = (axon_to_grow.straight_growths + 1.0) / (axon_to_grow.undulation_factor + 1.0);
    t = std::max(0.0, std::min(1.0, t));
    Eigen::Vector3d direction = slerp_direction(axon_to_grow.turn_start_direction, axon_to_grow.turn_target_direction, t);

    Eigen::Vector3d new_center = last_center + distance * direction;
    return new_center;
}


bool AxonGrowth::AddOneSphere(double radius_, bool create_sphere, int grow_straight, const int &factor, bool axon_can_shrink)
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
    int threshold_tries = (epsilon == 0.0) ? 1 : 100; // tries 100 times to place sphere

    // New sphere to attempt placing
    Sphere s(axon_to_grow.outer_spheres.size() + factor,
             axon_to_grow.id,
             /*object_type=*/axon_constant,
             axon_to_grow.begin,
             radius_);

    bool can_grow_ = false;
    int tries = 0;

    auto findCenter = [&]() -> Eigen::Vector3d {
        if (epsilon != 0.0 && grow_straight == 1) {
            return find_next_center_straight(distance, axon_to_grow.outer_spheres);
        }
        return find_next_center(distance, axon_to_grow.outer_spheres, axon_to_grow.end);
    };

    // Finds the neighbor causing the deepest overlap with a candidate sphere
    // (center, radius_) and returns a push vector, perpendicular to
    // growth_axis, that clears exactly that overlap (plus a small margin).
    // Self (this axon's own spheres) is never a blocker, matching
    // canSpherebePlaced's own same-object exclusion. Returns false if
    // nothing overlaps, or if the only overlap is directly along
    // growth_axis (no lateral escape direction).
    auto findPush = [&](const Eigen::Vector3d &center, Eigen::Vector3d &push) -> bool {
        auto candidates = sphere_grid->query(center, radius_ + max_radius_);
        double worst_overlap = 0.0;
        bool found = false;
        Eigen::Vector3d blocker_center = Eigen::Vector3d::Zero();
        for (const auto &entry : candidates) {
            if (entry.object_id == axon_constant && entry.cell_id == axon_to_grow.id) {
                continue; // this axon's own spheres
            }
            double d = (entry.center - center).norm();
            double overlap = (radius_ + entry.radius) - d;
            if (overlap > worst_overlap) {
                worst_overlap = overlap;
                blocker_center = entry.center;
                found = true;
            }
        }
        if (!found) {
            return false;
        }
        Eigen::Vector3d away = center - blocker_center;
        away[axon_to_grow.growth_axis] = 0.0;
        double away_norm = away.norm();
        if (away_norm < 1e-9) {
            return false;
        }
        away /= away_norm;
        const double push_margin = 1e-3 * radius_;
        push = away * (worst_overlap + push_margin);
        return true;
    };

    // Best shrink-fallback (largest achieved radius, i.e. least shrinkage)
    // found across all attempts this call -- used only if no attempt, pushed
    // or not, ever fits at the full radius_.
    bool have_fallback = false;
    Eigen::Vector3d fallback_center = Eigen::Vector3d::Zero();
    double fallback_radius = -1.0;

    // The shrink-fallback bisection is the expensive step (each iteration is
    // its own collision check): only pay for it on the last 20% of the try
    // budget, once cheap full-radius/pushed attempts have had most of their
    // chances -- always at least the final try, so a small threshold_tries
    // (e.g. 1, when epsilon == 0) still gets a shrink attempt.
    const int shrink_phase_start = std::max(0, threshold_tries - std::max(1, static_cast<int>(threshold_tries * 0.2)));

    // Push is a deterministic correction (always toward whichever neighbor
    // is nearest), not a random one -- if it fired on every failed attempt,
    // it would dominate the realized path at high packing density (where a
    // clean epsilon-random candidate rarely fits on its own), decoupling the
    // actual tortuosity from Tortuosity_Epsilon. Reserve the first half of
    // the try budget for pure, unpushed random resamples, so epsilon's
    // randomness gets a fair chance to find an unforced fit before any
    // correction kicks in; same "always at least the final try" floor as
    // shrink_phase_start for a small threshold_tries.
    const int push_phase_start = std::max(0, threshold_tries - std::max(1, static_cast<int>(threshold_tries * 0.5)));

    while (!can_grow_ && tries < threshold_tries)
    {
        // 1) Random search based on the existing attractor-biased direction.
        Eigen::Vector3d center = findCenter();
        bool inside = check_borders(extended_min_limits, extended_max_limits, center, radius_) || is_allowed_to_stop_early;

        // 2) Check for overlap at the full target radius.
        if (inside) {
            Sphere candidate(s.id, s.object_id, s.object_type, center, radius_);
            if (canSpherebePlaced(candidate)) {
                s.center = center;
                can_grow_ = true;
                break;
            }
        }

        // 3) There is overlap -- on the last 50% of tries, replace with the
        // pushed position and retry at the same full radius. Before that,
        // just move on to a fresh unpushed random resample (see
        // push_phase_start above).
        Eigen::Vector3d probe_center = center;
        bool probe_inside = inside;
        if (tries >= push_phase_start) {
            Eigen::Vector3d push;
            bool have_push = findPush(center, push);
            if (have_push) {
                probe_center = center + push;
                probe_inside = check_borders(extended_min_limits, extended_max_limits, probe_center, radius_) || is_allowed_to_stop_early;
                if (probe_inside) {
                    Sphere pushed_candidate(s.id, s.object_id, s.object_type, probe_center, radius_);
                    if (canSpherebePlaced(pushed_candidate)) {
                        s.center = probe_center;
                        can_grow_ = true;
                        break;
                    }
                }
            }
        }

        // 4) Still overlapping (pushed or not) -- on the last 20% of tries,
        // find the minimum shrinkage needed to fit at probe_center, and
        // remember it if it's the best (least-shrunk) fallback seen so far
        // this call.
        if (axon_can_shrink && probe_inside && tries >= shrink_phase_start) {
            double lo = 0.0, hi = radius_;
            Sphere floor_candidate(s.id, s.object_id, s.object_type, probe_center, min_radius);
            if (canSpherebePlaced(floor_candidate)) {
                lo = min_radius;
                for (int iter = 0; iter < 15; ++iter) {
                    double mid = (lo + hi) / 2.0;
                    Sphere trial(s.id, s.object_id, s.object_type, probe_center, mid);
                    if (canSpherebePlaced(trial)) {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                if (lo > fallback_radius) {
                    fallback_radius = lo;
                    fallback_center = probe_center;
                    have_fallback = true;
                }
            }
        }

        // 5) Try again (new random direction) until max_tries is met.
        ++tries;
    }

    // 6) If a full-radius fit was never found this way, fall back to the
    // position saved with minimum shrinkage.
    if (!can_grow_ && have_fallback) {
        s.center = fallback_center;
        s.radius = fallback_radius;
        can_grow_ = true;
    }

    // Evaluate the result
    if (!can_grow_) {
        // Collides or max tries reached, and no shrink-fallback available.
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

            // Dropping this point outright (as before) would leave its
            // neighbors twice as far apart as distance_between_spheres --
            // more than their radius sum wherever radii are locally small
            // (e.g. a beading trough), opening a real hole in the chain.
            // Shrink it via bisection instead, so *something* always sits
            // at this fixed spacing; only give up (and warn) if even a
            // near-zero radius still collides.
            if (!can_grow_) {
                double lo = 0.0, hi = rad;
                double best = 0.0;
                for (int iter = 0; iter < 20; ++iter) {
                    double mid = (lo + hi) / 2.0;
                    Sphere trial(id_, sph.object_id, sph.object_type, position, mid, sph.branch_id, sph.parent_id);
                    if (canSpherebePlaced(trial)) {
                        best = mid;
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                if (best > 0.0) {
                    s.radius = best;
                    can_grow_ = true;
                } else {
                    // Even a near-zero radius collides here: the fixed
                    // interpolated position itself sits inside another
                    // object's sphere, so no shrink can fix it. Push the
                    // position laterally away from whichever neighbor is
                    // responsible (same mechanism AddOneSphere's own findPush
                    // and the swelling phase's FindSwellPush already use),
                    // then retry the same shrink-to-fit search there --
                    // zeroing the growth_axis component keeps the push from
                    // fighting this layer's depth-cap accounting.
                    Eigen::Vector3d blocker_center;
                    double blocker_radius = 0.0;
                    double search_radius = std::max(rad, distance_between_spheres) * 2.0;
                    if (sphere_grid->findWorstOverlap(position, rad, search_radius,
                                                       sph.object_type, sph.object_id,
                                                       blocker_center, blocker_radius)) {
                        Eigen::Vector3d away = position - blocker_center;
                        away[axon_to_grow.growth_axis] = 0.0;
                        double away_norm = away.norm();
                        if (away_norm > 1e-9) {
                            away /= away_norm;
                            double overlap = (rad + blocker_radius) - (blocker_center - position).norm();
                            Eigen::Vector3d pushed_position = position + away * (overlap + 1e-3 * rad);

                            double lo2 = 0.0, hi2 = rad;
                            double best2 = 0.0;
                            for (int iter = 0; iter < 20; ++iter) {
                                double mid = (lo2 + hi2) / 2.0;
                                Sphere trial(id_, sph.object_id, sph.object_type, pushed_position, mid, sph.branch_id, sph.parent_id);
                                if (canSpherebePlaced(trial)) {
                                    best2 = mid;
                                    lo2 = mid;
                                } else {
                                    hi2 = mid;
                                }
                            }
                            if (best2 > 0.0) {
                                s.center = pushed_position;
                                s.radius = best2;
                                can_grow_ = true;
                            }
                        }
                    }
                }
            }

            if(can_grow_){
                last_id = s.id;
                axon_to_grow.add_sphere(s);
                //cout <<"sphere : " << s.id << " axon : "<< s.object_id << " can be placed as interpolated" << endl;
            } else {
                std::cerr << "Warning: axon " << axon_to_grow.id
                          << " interpolated sphere " << id_
                          << " could not be placed even at minimal radius or after a push -- possible discontinuity" << std::endl;
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
        generate_random_point_on_sphere(epsilon), target_direction
    );

    Eigen::Vector3d previous_vector = biased_random_vector; // no prior heading yet (axon's very first step)
    if (spheres.size() > 2) {
        previous_vector = (last_center - spheres[spheres.size() - 2].center).normalized();
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

    // A fresh target direction was just sampled -- rather than jumping to
    // it immediately and then holding it frozen for the rest of this turn
    // cycle (see find_next_center_straight), remember where we're turning
    // from and to, and start the same gradual slerp this step too, so the
    // whole cycle (this sampling step plus the following straight-phase
    // steps) curves smoothly toward the target instead of kinking sharply
    // on this one step alone.
    axon_to_grow.turn_start_direction = previous_vector;
    axon_to_grow.turn_target_direction = biased_random_vector;

    double t = (axon_to_grow.straight_growths + 1.0) / (axon_to_grow.undulation_factor + 1.0);
    t = std::max(0.0, std::min(1.0, t));
    Eigen::Vector3d direction = slerp_direction(previous_vector, biased_random_vector, t);

    return last_center + dist_ * direction;
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
    bool axon_can_shrink,
    std::size_t layer_start_spheres
) {
    auto set_stuck = [&](bool done) {
        finished = done ? true : false;                  // explicit member
        stuck_radius   = done ? axon_to_grow.radius : -1.0;
        stuck_index    = done ? axon_to_grow.id     : -1;
    };
    // Seeded from the axon's own persisted state (not reset to 0) so the
    // straight/random cadence carries across depth-layer boundaries instead
    // of restarting every layer -- growthThread is called on a fresh
    // AxonGrowth per axon per layer, so without this every layer would force
    // an extra epsilon-random redirection regardless of where the axon was
    // in its straight streak when the previous layer ended.
    int  grow_straight     = axon_to_grow.grow_straight;
    int  straight_growths  = axon_to_grow.straight_growths;

    std::size_t tries = 0;
    const std::size_t max_tries = 1000000;  // large enough; growth logic should end earlier
    // Swept 5/10/15/20 on the 75% packing benchmark: every value gave zero
    // abandoned axons and final ICVF within under a point of each other
    // (74.3-75.0%), but runtime scaled clearly with more attempts (374s at
    // 5 vs 707-711s at 15/20) -- the push+shrink+round-retry machinery
    // already does the heavy lifting, so extra retries mostly just cost
    // time without buying more density. TEMP: still overridable via
    // TOTAL_NBR_GROWTH_ATTEMPTS for further A/B testing.
    int total_nbr_growth_attempts = 5;
    if (const char *env_val = std::getenv("TOTAL_NBR_GROWTH_ATTEMPTS")) {
        total_nbr_growth_attempts = std::atoi(env_val);
    }
    int jostle_rounds = 0;
    int max_jostle_rounds = 5; // TEMP: overridable via MAX_JOSTLE_ROUNDS for A/B testing
    if (const char *env_val = std::getenv("MAX_JOSTLE_ROUNDS")) {
        max_jostle_rounds = std::atoi(env_val);
    }

    while (!finished && tries < max_tries) {
        ++tries;

        // Radius beading
        double varied_radius = axon_to_grow.radius;
        if (axon_to_grow.beading_amplitude > 0) {
            varied_radius = RandomradiusVariation();
        }

        // Try to place a sphere: random search based on the attractor bias,
        // pushing away from a blocker if the first try overlaps, and
        // (if axon_can_shrink) falling back to the least-shrunk position
        // found across all attempts if a full-radius fit is never found --
        // see AddOneSphere's doc comment.
        const bool grew = AddOneSphere(varied_radius, /*create_sphere=*/true, grow_straight, factor, axon_can_shrink);

        // Growth process may have marked itself finished (reached this call's
        // extended_max_limits[growth_axis] -- either this layer's depth cap or,
        // on the final layer, the true wall; the caller tells these apart)
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

        // Could not grow (no fit found even after push + shrink-fallback
        // inside AddOneSphere) → retry a few times, then give up. Both
        // the retry and the give-up path roll back only to this layer's
        // starting sphere count, not all the way to the axon's very first
        // sphere: previously-committed layers are already merged into the
        // shared sphere_grid, and fully wiping them here (as keep_one_sphere()/
        // destroy() would) would desync this axon's local state from what's
        // still sitting in that grid.
        if (axon_to_grow.growth_attempts < total_nbr_growth_attempts) {
            axon_to_grow.truncate_to(layer_start_spheres);  // retry this layer from its own start
            axon_to_grow.growth_attempts += 1;
            set_stuck(false);
        } else {
            // Retrying from the exact same spot has been exhausted. Before
            // giving up on this layer entirely, try a bounded number of
            // lateral jostles: nudge the axon's aim point (axon_to_grow.end,
            // the target find_next_center steers toward) away from whichever
            // neighbor is closest to the tip, then give it a fresh
            // growth_attempts budget. The tip itself can't be relocated like
            // a Phase-A seed circle without breaking the chain behind it, so
            // this steers future attempts around the obstruction instead of
            // moving anything already placed -- a real second chance to
            // route past a local blocker, not just resampling the same
            // biased-random directions again.
            bool jostled = false;
            if (jostle_rounds < max_jostle_rounds && !axon_to_grow.outer_spheres.empty()) {
                const Sphere &tip = axon_to_grow.outer_spheres.back();
                auto candidates = sphere_grid->query(tip.center, axon_to_grow.radius * 3.0);
                double best_dist = std::numeric_limits<double>::max();
                Eigen::Vector3d blocker_center = Eigen::Vector3d::Zero();
                bool found_blocker = false;
                for (const auto &entry : candidates) {
                    if (entry.object_id == axon_constant && entry.cell_id == axon_to_grow.id) {
                        continue; // this axon's own spheres
                    }
                    double d = (entry.center - tip.center).norm();
                    if (d < best_dist) {
                        best_dist = d;
                        blocker_center = entry.center;
                        found_blocker = true;
                    }
                }
                if (found_blocker) {
                    Eigen::Vector3d away = tip.center - blocker_center;
                    away[axon_to_grow.growth_axis] = 0.0; // steer within the cross-section only
                    double away_norm = away.norm();
                    if (away_norm > 1e-9) {
                        away /= away_norm;
                        double nudge = 0.5 * std::max(tip.radius, min_radius);
                        axon_to_grow.end += away * nudge;
                        axon_to_grow.growth_attempts = 0;
                        ++jostle_rounds;
                        jostled = true;
                    }
                }
            }

            if (jostled) {
                set_stuck(false);
            } else {
                axon_to_grow.truncate_to(layer_start_spheres);  // give up on this layer, keep prior layers
                set_stuck(true);
                break;
            }
        }

        update_straight(false, grow_straight, straight_growths);
    }

    if (tries >= max_tries && !finished) {
        // safety stop
        axon_to_grow.truncate_to(layer_start_spheres);
        set_stuck(true);
    }

    // Persist the straight/random cadence back onto the axon so the next
    // layer's growthThread call (a fresh AxonGrowth, but the same Axon)
    // picks up where this one left off.
    axon_to_grow.grow_straight = grow_straight;
    axon_to_grow.straight_growths = straight_growths;

    // Note: the previous "reached the wall but ended up with fewer than 10
    // spheres total -> discard" prune is now handled once, globally, by the
    // caller after all layers are done (growAxonsLayered) -- not here, since
    // outer_spheres.size() here only reflects progress up to and including
    // the CURRENT layer, and destroying wholesale would hit the same
    // sphere_grid desync problem for any axon with committed prior layers.
}

