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
    // Only the last sphere's center is ever read below -- direction comes from
    // the persisted turn_start_direction/turn_target_direction, not from a
    // second point -- so the real precondition is just "not empty". A stricter
    // "< 2" guard here used to assert(0) whenever a collision rollback (see
    // Axon::truncate_to) left an axon with only its seed sphere while
    // grow_straight was still 1 from before the rollback: a legitimate,
    // reachable state (now additionally guarded against by truncate_to
    // resetting grow_straight itself), not a real invariant violation.
    if (spheres.empty()){
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

    // `factor` used to belong here because outer_spheres.size() counted
    // `factor` gap-filling spheres per real hop; growth is backbone-only
    // now (see the call site below), so size() already IS the real hop
    // count and the old *factor would just make this cap `factor` times
    // more permissive than intended.
    if (axon_to_grow.outer_spheres.size() > (max_limits-min_limits).norm()/(axon_to_grow.radius)*100) {
        finished = true;
        return false; // Axon has grown too long
    }

    assert(axon_to_grow.growth_axis >= 0 && axon_to_grow.growth_axis < 3);

    bool is_allowed_to_stop_early = axon_to_grow.outside_voxel;


    // If the last sphere is beyond extended_max_limits, the axon is considered fully grown.
    // Must match hasReachedTrueWall's (and the check at line ~319 below's) center+radius
    // test, not center alone -- once extended_max_limits == the true wall (once the
    // depth-layer cap has grown past the box), a center-only test can report "finished"
    // here while hasReachedTrueWall still says "not actually at the wall" for the very
    // same sphere, since a sphere can have center < limit but center+radius > limit.
    // The caller then keeps this axon in next_active forever -- neither done nor
    // stuck -- re-entering this same immediate-return path every subsequent layer with
    // zero growth, a real livelock only bounded by growAxonsLayered's max_layers backstop.
    Sphere last_sphere = axon_to_grow.outer_spheres.back();
    if (last_sphere.center[axon_to_grow.growth_axis] + last_sphere.radius > extended_max_limits[axon_to_grow.growth_axis] && !is_allowed_to_stop_early) {
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
        // Non-allocating, same as canSpherebePlaced and FindSwellPush (the
        // swelling-phase equivalent of this exact search): query() would
        // heap-allocate a fresh std::vector and copy in every candidate
        // before any of them are even checked, and this runs inside
        // AddOneSphere's retry loop -- the hottest path in growth,
        // especially at high packing density where pushes fire often.
        Eigen::Vector3d blocker_center = Eigen::Vector3d::Zero();
        double blocker_radius = 0.0;
        bool found = sphere_grid->findWorstOverlap(center, radius_, radius_ + max_radius_,
                                                     axon_constant, axon_to_grow.id,
                                                     blocker_center, blocker_radius);
        if (!found) {
            return false;
        }
        double worst_overlap = (radius_ + blocker_radius) - (blocker_center - center).norm();
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
                    // A push is lateral (findPush zeroes the growth_axis
                    // component), so it always moves the candidate FURTHER
                    // from the previous sphere: |probe - prev| =
                    // sqrt(d^2 + |push|^2) > d. canSpherebePlaced only ever
                    // tests against *other* objects -- this axon's own
                    // spheres are deliberately excluded -- so nothing here
                    // notices the candidate drifting out of contact with its
                    // own predecessor, exactly the blind spot that lets a
                    // stranded shrink-fallback sever the chain. Require
                    // continuity explicitly, mirroring the
                    // still_touches_neighbor guard ComputeSwollenSphere
                    // already applies to the swelling-phase push.
                    Sphere pushed_candidate(s.id, s.object_id, s.object_type, probe_center, radius_);
                    const bool touches_prev =
                        (probe_center - last_sphere.center).norm() <= radius_ + last_sphere.radius;
                    if (touches_prev && canSpherebePlaced(pushed_candidate)) {
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
        // The step out to fallback_center was committed from `distance`
        // above -- i.e. sized for the radius this sphere WANTED to be, before
        // any shrinking was considered. Placing a shrunken sphere there
        // strands it: staying joined to the previous sphere needs
        // d < r_new + r_prev, and d was reserved for an r_new this sphere no
        // longer has. For a myelinated axon `distance` is the axon-level
        // inner_radius, so a sphere bisected down to a fraction of that
        // cannot reach back at all -- and a run of such spheres severs the
        // axon into disconnected pieces (measured on the 70% set: necks
        // collapsing to r=0.0007 while the stride stayed at a constant
        // 0.072, leaving 1.55% of axons in several fragments).
        // So shorten the stride to match the radius actually achieved: pull
        // the centre back along the line from the previous sphere until the
        // step is max(r_new, r_prev), which guarantees overlap because
        // max(a,b) < a+b for positive a,b. Moving toward an already-placed,
        // already-collision-free neighbour is the safe direction, but it is
        // still re-checked, and if the fully-pulled position does not fit we
        // bisect back toward the original centre and take the closest spot
        // that does.
        const Eigen::Vector3d prev_center = last_sphere.center;
        Eigen::Vector3d dir = fallback_center - prev_center;
        const double cur_d = dir.norm();
        const double want_d = std::max(fallback_radius, last_sphere.radius);
        if (cur_d > want_d && cur_d > 1e-9) {
            dir /= cur_d;
            double lo_d = want_d;   // closest (guarantees overlap), may collide
            double hi_d = cur_d;    // original (fits, but may leave a gap)
            Eigen::Vector3d best_center = fallback_center;
            Sphere pulled(s.id, s.object_id, s.object_type, prev_center + dir * lo_d, fallback_radius);
            if (canSpherebePlaced(pulled)) {
                best_center = pulled.center;
            } else {
                for (int iter = 0; iter < 12; ++iter) {
                    double mid_d = (lo_d + hi_d) / 2.0;
                    Sphere trial(s.id, s.object_id, s.object_type, prev_center + dir * mid_d, fallback_radius);
                    if (canSpherebePlaced(trial)) {
                        best_center = trial.center;
                        hi_d = mid_d;   // fits -- try to come closer still
                    } else {
                        lo_d = mid_d;
                    }
                }
            }
            fallback_center = best_center;
        }
        // Final continuity requirement. The pull above is best-effort: if the
        // fully-pulled position and every bisected position between it and the
        // original both collide, best_center stays at the ORIGINAL centre --
        // which, when that centre came from a lateral push, can sit far beyond
        // reach of the previous sphere (measured: a stride of 2.23um where
        // max(r, r_prev) was 0.98um, leaving a 0.79um hole). Placing the
        // sphere anyway is what actually severs the axon during growth.
        // There is no valid position here, so refuse the step rather than
        // commit a disconnected chain: returning false lets growthThread spend
        // its retry/jostle budget and, failing that, hand the axon to the
        // relocate-and-retry machinery, which is the designed response to
        // "this axon cannot continue from here".
        //
        // The bar here is contact (d <= r_new + r_prev), not the full growth
        // stride `distance`. Requiring the stride was tried and does keep the
        // chain tighter, but it rejects so many fallback steps that
        // abandonment roughly doubles (55 -> 114 axons on the 40um/60%/c2=0.6
        // benchmark), and abandonment is size-biased -- it preferentially
        // drops the large axons, narrowing the achieved radius distribution
        // against the target Gamma. Contact is what actually distinguishes a
        // connected axon from a severed one, so that is what is enforced.
        if ((fallback_center - last_sphere.center).norm() > fallback_radius + last_sphere.radius) {
            return false;
        }
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
        // Backbone-only during growth: no gap-filling spheres are inserted
        // here anymore (see CaterpillarGrowth::InterpolateAllAxons, run once
        // after every axon has finished growing AND swelling). Consecutive
        // backbone spheres are already spaced about one radius apart (see
        // `distance` above), close enough that other axons still growing
        // don't need infill to detect them via sphere_grid.
        axon_to_grow.add_sphere(s);
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


// add_spheres (growth-time gap-filling between consecutive backbone
// spheres) used to live here. Growth is backbone-only now -- see the
// comment at its call site in AddOneSphere -- and the same lerp-position/
// lerp-radius/bisection-shrink interpolation is done once, after every
// axon has finished growing AND swelling, by
// CaterpillarGrowth::InterpolateAllAxons (CaterpillarGrowth.cpp), which
// needs direct mutable sphere_grid access (to insert the new spheres) that
// this class's const SphereGrid* deliberately doesn't have.


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
                // Non-allocating (see findPush above for the same
                // reasoning): this fallback is rarer than findPush, but
                // still no reason to pay query()'s heap allocation for it.
                Eigen::Vector3d blocker_center = Eigen::Vector3d::Zero();
                double blocker_radius = 0.0;
                bool found_blocker = sphere_grid->findNearest(tip.center, axon_to_grow.radius * 3.0,
                                                                axon_constant, axon_to_grow.id,
                                                                blocker_center, blocker_radius);
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

