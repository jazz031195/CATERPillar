#include "CaterpillarGrowth.h"
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

// A grown process shorter than this (soma-to-tip, in um) is discarded outright
// rather than kept as a visually degenerate stub -- independent of how far
// below the population's configured mean/std a particular length draw fell.
static constexpr double kMinGlialProcessLength = 5.0;


GlialCellGrowth::~GlialCellGrowth() {}

GlialCellGrowth::GlialCellGrowth(Glial &glial_cell_to_grow_,
            const SphereGrid* sphere_grid_,
            const Eigen::Vector3d &extended_min_limits_,
            const Eigen::Vector3d &extended_max_limits_,
            const Eigen::Vector3d &min_limits_,
            const Eigen::Vector3d &max_limits_,
            const double &min_radius_, const double &epsilon_) : CellGrowth(sphere_grid_, extended_min_limits_, extended_max_limits_, min_limits_, max_limits_, epsilon_, min_radius_), glial_cell_to_grow(glial_cell_to_grow_){}

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
            id_ = glial_cell_to_grow.next_sphere_id++;
            Sphere s(id_, last_sphere.object_id, last_sphere.object_type, position, rad, last_sphere.branch_id, sph.parent_id);
            bool can_grow_ = canSpherebePlaced(s, check_collision_with_branches);
            if(can_grow_){
                glial_cell_to_grow.ramification_spheres[index_ram_spheres].push_back(s);
            }
        }
    }
    sph.id = glial_cell_to_grow.next_sphere_id++;
    glial_cell_to_grow.ramification_spheres[index_ram_spheres].push_back(sph);
}


bool GlialCellGrowth::AddOneSphere(const double &radius_, const bool &create_sphere, int &grow_straight, const int &i, const bool &check_collision_with_branches, const int &parent, const int &factor)
{

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
    Sphere s(id_, glial_cell_to_grow.id, glial_cell_constant, {0,0,0}, radius_, i);

    int tries = 0;
    // Only the first 100 tries keep AddOneSphere's original behavior --
    // epsilon-biased toward this branch's fixed attractor, unchanged, so the
    // common case (that already succeeds well within 100 tries) isn't
    // slowed down. Every one of those 100 tries resamples a fresh direction,
    // but always within the same narrow cone centered on the same fixed,
    // far-away attractor -- if a local obstacle blocks that whole cone (e.g.
    // another cell's soma sitting roughly along this branch's original
    // heading), all 100 can fail even with plenty of free space just off to
    // the side, well before the voxel is anywhere near full. Observed to
    // noticeably truncate primary glial processes short of their sampled
    // target even at ~5% overall ICVF, where genuine global crowding is an
    // unlikely explanation. Two further stages progressively widen the
    // angular spread so a genuinely blocked branch gets real chances to
    // route around the obstacle instead of just re-rolling the same cone.
    int threshold_tries = 300;

    while (!can_grow_ && tries < threshold_tries){

        double std_override = -1.0; // default: find_next_center falls back to epsilon, unchanged from before
        if (tries >= 200) {
            // cos(phi) ~ N(1, 3) clamped to [-1, 1] piles up near both poles
            // with a broad spread between -- close enough to an unbiased
            // direction over the whole sphere for a last-resort escape try.
            std_override = 3.0;
        } else if (tries >= 100) {
            std_override = epsilon + 1.0; // moderately widened first
        }

        find_next_center(s, distance, glial_cell_to_grow.ramification_spheres[i], destination, std_override);
        // check if there is a collision
        can_grow_ = canSpherebePlaced(s, check_collision_with_branches);

        tries += 1;
    }

    if (can_grow_ )
    {
        if (create_sphere){
            //cout << "adding sphere to glial, radius "<< s.radius << ", position : " << s.center << endl;
            s.parent_id = parent;
            add_spheres(s, last_sphere, check_collision_with_branches, factor, i);
        }
        // Processes are allowed to keep growing past the voxel boundary --
        // only their own sampled length budget (growPrimaryBranch/
        // growSecondaryBranch) or a genuine collision stops them. Where a
        // branch may *originate* is still gated separately (see
        // growSecondaryBranch's source-inside-voxel check).
        return true;
    }
    else // collides and >1000 tries
    {
        //cout << "could not grow glial" << endl;
        return false;
    }
        
}

void GlialCellGrowth::find_next_center_straight(double distance, Sphere &s, const std::vector<Sphere> &spheres)
{

    Eigen::Vector3d normal = (spheres[spheres.size() - 1].center-glial_cell_to_grow.soma.center).normalized();
    Eigen::Vector3d new_center = spheres[spheres.size() - 1].center + normal*distance;
    s.center = new_center;
}

void GlialCellGrowth::find_next_center(Sphere &s,  double dist_, const std::vector<Sphere> &spheres, const Eigen::Vector3d &target, double std_override)
{

    Eigen::Vector3d target_ = target;
    Eigen::Vector3d vector_to_target = target_ - spheres[spheres.size() - 1].center;
    vector_to_target = vector_to_target.normalized();
    double std_ = (std_override >= 0.0) ? std_override : epsilon;
    Eigen::Vector3d vector = generate_random_point_on_sphere(std_);
    vector = apply_bias_toward_target(vector, vector_to_target);
    Eigen::Vector3d position = spheres[spheres.size() - 1].center + dist_ * vector.normalized();
    s.center = position;
}


// Function to find intersection points between a vector and each face of the cube
Eigen::Vector3d findDistantPoint(const Eigen::Vector3d &vector, const Eigen::Vector3d &point, const double &L)
{
    double distance = L + 10; // Large initial distance
    Eigen::Vector3d distant_point = point + distance * vector;
    return distant_point;
}


bool GlialCellGrowth::GenerateFirstSphereinProcess(Sphere &first_sphere, Eigen::Vector3d &attractor, const double &radius, const Sphere &sphere_to_emerge_from, const Eigen::Vector3d &vector_to_prev_center, const int &nbr_spheres, const int &nbr_spheres_between, const int &cell_id, const int &branch_id, const bool &primary_process) {
    
    bool stop = false;
    int tries_ = 0;
    int max_nbr_tries = 10;
    attractor = Eigen::Vector3d(0, 0, 0);


    while (!stop && tries_ < max_nbr_tries) {
        Eigen::Vector3d vector = {0, 0, 0};
        Eigen::Vector3d point = {0, 0, 0};
        //cout << "Generating first sphere in process, try: " << tries_ + 1 << endl;
        sphere_to_emerge_from.getPointOnSphereSurface(point, vector, vector_to_prev_center, primary_process);
        //cout << "Point on sphere surface: " << point.transpose() << ", Direction vector: " << vector.transpose() << endl;
        attractor = findDistantPoint(vector, point, max_limits[0]);
        first_sphere = Sphere(glial_cell_to_grow.next_sphere_id++, cell_id, glial_cell_constant, point, radius, branch_id, sphere_to_emerge_from.id);
        //cout <<"check placement for first sphere at position: " << first_sphere.center.transpose() << " with radius: " << first_sphere.radius << endl;
        if (canSpherebePlaced(first_sphere, /*check_collision_with_branches=*/ false)) {
            stop = true;
        }
        else {
            //cout << "First sphere collides with existing structures." << endl;
            tries_++;
            if (tries_ == max_nbr_tries) return false;
        }
    }
    return true;
}


std::vector<Sphere> GlialCellGrowth::addIntermediateSpheres(const Sphere &random_sphere, const Sphere &first_sphere,  const int &branch_nbr, const int &nbr_spheres, const int &nbr_spheres_between, const std::function<double(double)> &compute_radius, const double &t_start, const double &t_end) {

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
            glial_cell_to_grow.next_sphere_id++, first_sphere.object_id, first_sphere.object_type, position, rad, branch_nbr,
            random_sphere.id);

        // Check if the sphere can be placed
        if (canSpherebePlaced(next, /*check_collision_with_branches=*/ false)) {
            intermediate_spheres.emplace_back(next);
        }
    }

    // Add the initial branching sphere to the list
    intermediate_spheres.push_back(first_sphere);

    return intermediate_spheres;
}

void GlialCellGrowth::growFirstPrimaryBranches(const int &number_ramification_points, int &nbr_spheres, const double &mean_process_length, const double &std_process_length, const int &factor)
{
    int nbr_tries = 0;
    int max_nbr_tries = 1000;
    for (int j = 0; j < number_ramification_points; j++)
    {
        bool has_grown = growPrimaryBranch(nbr_spheres, mean_process_length, std_process_length, factor);
        if (!has_grown && nbr_tries < max_nbr_tries){
            j = j - 1;
            nbr_tries += 1;
        }
        else if (nbr_tries >= max_nbr_tries){
            cout << "Failed to grow glial cell" << endl;
        }
        else{
            nbr_tries = 0;
        }
        
    }
}

bool GlialCellGrowth::growPrimaryBranch(int &nbr_spheres, const double &mean_primary_process_length, const double &std_primary_process_length, const int &factor) {

    //cout << "Growing primary branch for glial cell : "<<  glial_cell_to_grow.id << "/"<< glial_cells.size()<< endl;
    int j = glial_cell_to_grow.ramification_spheres.size();

    // Generate length from Gaussian distribution
    std::random_device rd;
    std::mt19937 generator(rd());
    std::normal_distribution<double> length_dist(mean_primary_process_length, std_primary_process_length);
    // Floor relative to the requested mean, not the sample itself: a Gaussian
    // draw can land small or negative (e.g. whenever std_primary_process_length
    // is comparable to or larger than the mean, as the GUI's own defaults are),
    // and a floor derived from that same draw is equally small/negative -- it
    // stops acting as a floor right when it's needed most, letting a branch
    // that grows only its first (on-the-soma) sphere slip through as "valid".
    double min_length = std::max(0.5*mean_primary_process_length, kMinGlialProcessLength);
    double length = length_dist(generator);
    int tries = 0;
    while (length < min_length && tries < 100) {
        length = length_dist(generator);
        tries++;
    }

    if (length < min_length) {
        return false;
    }

    double initial_radius = glial_cell_to_grow.soma.radius/3;
    // Cap the decay length at the mean (not the sampled) target: length here can
    // occasionally be drawn much larger than typical, and actual growth often
    // stops well short of it (collisions, running out of room, etc.) -- decay
    // calibrated to an inflated target barely moves over the distance actually
    // travelled, leaving the process looking uniformly thick instead of
    // tapering toward minimum_radius.
    double alpha = -std::log(glial_cell_to_grow.minimum_radius / initial_radius)/std::min(length, mean_primary_process_length);

    auto compute_radius = [&](double t) {
        return std::exp(-alpha * t) * initial_radius;
    };

    //cout << "Target length for primary branch: " << length << endl;


    // Find points on the surface of the glial soma sphere
    tries = 0;
    // Add the first sphere with intermediate spheres
    int parent = 0;  // parent is the soma
    int nbr_spheres_between = factor - 1;

    Eigen::Vector3d vector_to_prev_center = {0, 0, 0};
    Sphere first_sphere;
    Eigen::Vector3d attractor = Eigen::Vector3d(0, 0, 0);
    // t=0 here means "at the branch's own start", i.e. the true initial_radius
    // with no decay yet -- passing initial_radius itself as t (a radius value
    // used as if it were a distance-along-branch) was a bug.
    double first_radius = compute_radius(0.0);
    bool first_sphere_created = GenerateFirstSphereinProcess(first_sphere, attractor, first_radius, glial_cell_to_grow.soma, vector_to_prev_center, nbr_spheres, nbr_spheres_between, glial_cell_to_grow.id, j, true);

    if (!first_sphere_created) {
        return false;
    }

    // No intermediate spheres between the soma and first_sphere: unlike a
    // branch-to-branch attachment (where the parent sphere's own radius is a
    // short, branch-scale hop worth subdividing), first_sphere already sits
    // directly on the soma's surface with zero gap -- there's nothing to
    // bridge. Interpolating from the soma's *center* (as addIntermediateSpheres
    // would) placed several branch-radius-sized spheres deep inside the soma's
    // own volume, visible in the GUI as a big lump right at the base of every
    // primary process.
    std::vector<Sphere> vector_first_spheres = {first_sphere};

    //cout << "First sphere in primary branch created at position: " << first_sphere.center.transpose() << " with radius: " << first_sphere.radius << endl;

    int nbr_non_checked_spheres = factor * 3;

    // Grow spheres in the branch
    bool can_grow = true;

    Eigen::Vector3d prev_pos = first_sphere.center;

    Glial old_glial_cell = glial_cell_to_grow;

    glial_cell_to_grow.ramification_spheres.resize(j + 1);
    glial_cell_to_grow.lengths_branches.resize(j + 1);
    glial_cell_to_grow.attractors.resize(j + 1);
    glial_cell_to_grow.ramification_spheres[j] = vector_first_spheres;
    glial_cell_to_grow.attractors[j] = attractor;
    glial_cell_to_grow.lengths_branches[j] = std::vector<double>(vector_first_spheres.size(), first_sphere.radius);
    //Eigen::Vector3d expanded_space_min = min_limits - Eigen::Vector3d({expanded_for_glial_space, expanded_for_glial_space, expanded_for_glial_space});
    //Eigen::Vector3d expanded_space_max = max_limits + Eigen::Vector3d({expanded_for_glial_space, expanded_for_glial_space, expanded_for_glial_space});

    //cout << "Growing primary branch..." << endl;
    double distance = initial_radius + first_radius;
    while (can_grow) {

        // Radius at the current distance, clamped to minimum_radius so the
        // branch keeps growing in *length* toward its true sampled target
        // even once decay alone would have taken it below that floor. alpha
        // above is deliberately calibrated against min(length, mean) so the
        // radius still visibly tapers within about a mean-length's worth of
        // distance even when the true target is much longer -- but that must
        // only govern how fast the radius shrinks, not how far the branch is
        // allowed to grow. Previously R_ <= minimum_radius doubled as the
        // loop's own stop condition, so any branch with length > mean (about
        // half of all draws) silently had its actual grown length capped at
        // ~mean regardless of its sampled target -- e.g. never observed to
        // exceed mean_process_length at all, no matter how large std was.
        double R_ = std::max(compute_radius(distance), glial_cell_to_grow.minimum_radius);
        if (distance >= length) {
            break;
        } else {
            bool check_collision = glial_cell_to_grow.ramification_spheres[j].size() >= nbr_non_checked_spheres;
            int grow_straight = 0;
            //cout << "Adding sphere with radius: " << R_ << " at distance: " << distance << endl;
            can_grow = AddOneSphere(R_, true, grow_straight, j, check_collision, parent, factor);
        }

        // AddOneSphere can add spheres (via add_spheres, up to `factor` of them:
        // intermediate spheres plus the sphere itself) and *still* return false
        // right afterward (e.g. the sphere lands outside the extended voxel
        // bounds) -- so lengths_branches must be resynced against however many
        // spheres actually landed in ramification_spheres[j], regardless of
        // AddOneSphere's return value, or it silently falls behind. Once that
        // happens it stays wrong for the rest of this branch's life, and any
        // later secondary branch that attaches partway along it (via
        // lengths_branches[branch][index] in growSecondaryBranch) reads
        // past the end of a too-short lengths_branches[j] -- undefined
        // behavior, observed in practice as garbage old_length values.
        {
            int nbr_spheres_in_branch_before = glial_cell_to_grow.lengths_branches[j].size();
            int nbr_spheres_in_branch_after = glial_cell_to_grow.ramification_spheres[j].size();
            if (nbr_spheres_in_branch_after > nbr_spheres_in_branch_before) {
                const auto &new_sphere = glial_cell_to_grow.ramification_spheres[j].back();
                distance += (prev_pos - new_sphere.center).norm();
                prev_pos = new_sphere.center;
                for (int i = nbr_spheres_in_branch_before; i < nbr_spheres_in_branch_after; i++) {
                    glial_cell_to_grow.lengths_branches[j].push_back(distance);
                }
            }
        }
    }

    //cout << "Finished growing primary branch. Total length grown: " << distance << endl;


    // Check if the branch has grown sufficiently
    if (distance < min_length) {
        glial_cell_to_grow = std::move(old_glial_cell);
        return false;
    }
    else if (glial_cell_to_grow.ramification_spheres[j].size() != 0) {
        //cout << "Branch grown" << endl;
        // Update the number of spheres
        nbr_spheres += glial_cell_to_grow.ramification_spheres.back().size();

        return true;
    }
    return true;


}

bool GlialCellGrowth::growSecondaryBranch(int &nbr_spheres, const double &mean_process_length, const double &std_process_length, const int &factor) {

    //cout << "Growing secondary branch for glial cell : "<<  glial_cell.id<< endl;
    if (glial_cell_to_grow.ramification_spheres.empty()) {
        cout << "No branches in glial cell" << endl;
        return false;
    }
    if (nbr_spheres <= 0) {
        cout << "No spheres in glial cell : " <<  glial_cell_to_grow.id << endl;
        return false;
    }

    int nbr_branches = glial_cell_to_grow.ramification_spheres.size();

    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> branch_dist(0, nbr_branches - 1);
    finished = false;

    // finding source of branching
    int nbr_spheres_between = factor - 1;


    int random_branch = branch_dist(rng);
    int nbr_tries = 0;
    // size <= 1 is rejected alongside empty: random_sphere_ind below must land
    // at index >= 1 (it needs a *previous* sphere to derive a direction from),
    // which a single-sphere branch can't offer.
    while (glial_cell_to_grow.ramification_spheres[random_branch].size() <= 1 && nbr_tries < 100) {
        random_branch = branch_dist(rng);
        nbr_tries++;
    }
    if (nbr_tries >= 100) {
        cout << "Could not find a non-empty branch for glial cell : " << glial_cell_to_grow.id << endl;
        return false;
    }

    int size = glial_cell_to_grow.ramification_spheres[random_branch].size();
    //int random_sphere_ind = size / 2 + rand() % (size - size / 2);
    int random_sphere_ind = 1 + rand() % (size - 1); // >= 1: needs [random_sphere_ind - 1] below
    Sphere random_sphere = glial_cell_to_grow.ramification_spheres[random_branch][random_sphere_ind];

    // A secondary branch may only originate from a point that is actually
    // inside the real voxel (min_limits/max_limits, not the extended
    // growth region) -- unlike ordinary growth, which is now allowed past
    // the boundary (see AddOneSphere), a *new* branch is not allowed to
    // spring from a source that already lies outside it. Zero margin (not
    // random_sphere.radius): the source's own radius would otherwise dilate
    // the box and let a source sitting just past the true edge still count
    // as "inside".
    if (!check_borders(min_limits, max_limits, random_sphere.center, 0.0)) {
        return false;
    }

    if (glial_cell_to_grow.lengths_branches.size() <= random_branch) {
        cout << "glial_cell.lengths_branches.size() : " << glial_cell_to_grow.lengths_branches.size() << endl;
        cout << "random_branch : " << random_branch << endl;
        cout << "Error: lengths_branches vector is not large enough." << endl;
        assert(0);
    }
    if (glial_cell_to_grow.lengths_branches[random_branch].size() != glial_cell_to_grow.ramification_spheres[random_branch].size()) {
        // Would previously read random_sphere_ind out of bounds on
        // lengths_branches[random_branch] (undefined behavior, observed as
        // garbage old_length values) whenever a branch's growth loop added
        // spheres via AddOneSphere without a matching lengths_branches
        // update -- see the resync fix in growPrimaryBranch/growSecondaryBranch's
        // growth loops. Kept as a loud assertion rather than silently
        // trusting the sizes match, since a future change reintroducing that
        // desync should fail immediately here instead of corrupting length
        // sampling downstream.
        cout << "glial_cell.ramification_spheres[" << random_branch << "].size() : "
             << glial_cell_to_grow.ramification_spheres[random_branch].size()
             << " != lengths_branches[" << random_branch << "].size() : "
             << glial_cell_to_grow.lengths_branches[random_branch].size() << endl;
        assert(0);
    }
    double old_length = glial_cell_to_grow.lengths_branches[random_branch][random_sphere_ind];
    std::random_device rd;
    std::mt19937 generator(rd());

    // Cap this branch's own total real (travelled arc-length, not straight-line
    // Euclidean distance -- glial processes wander, so "length" here means the
    // actual path length, matching how process length is specified/measured
    // elsewhere in this model) reach from the soma at its parent branch's own
    // total real length: lengths_branches[random_branch].back() is the parent's
    // own cumulative arc length at its tip, and old_length (already one of the
    // parent's own lengths_branches entries) is how far along the parent this
    // branch is attaching, so the parent's own remaining reach from here is
    // their difference.
    double parent_total_length = glial_cell_to_grow.lengths_branches[random_branch].back();
    double max_length_to_grow = std::max(0.0, parent_total_length - old_length);
    if (max_length_to_grow < kMinGlialProcessLength) {
        // Attaching this close to (or past) the parent's own tip leaves no
        // room for a worthwhile child branch without it overshooting the
        // parent -- discard the attempt outright, same as growPrimaryBranch's
        // own "not enough room" rejection.
        return false;
    }

    // Every process's own length -- primary or secondary, regardless of
    // nesting depth or how far along its parent it attaches -- is drawn from
    // the same N(mean_process_length, std_process_length), matching
    // growPrimaryBranch. This used to be sampled from
    // N(mean_process_length - old_length, std_process_length) instead, aiming
    // to make the *cumulative* soma-to-tip length hit N(mean, std): but
    // combined with max_length_to_grow (the parent-reach cap above), that
    // shrank the *effective* sampling window down to a sliver for any branch
    // attaching more than a little way along its parent, badly compressing
    // the realized standard deviation of secondary lengths well below
    // std_process_length -- measured at less than half the configured value
    // in practice. Sampling the branch's own length independently of
    // old_length, and letting max_length_to_grow act purely as a rejection
    // filter (redraw/discard, never shift the distribution), keeps every
    // accepted sample's underlying distribution faithful to std_process_length.
    std::normal_distribution<double> length_dist(mean_process_length, std_process_length);
    double length_to_grow = length_dist(generator);
    // Just the absolute floor here (unlike growPrimaryBranch's 0.5*mean, which
    // has no upper bound to conflict with): max_length_to_grow is often well
    // under mean_process_length for a branch attaching partway along its
    // parent, so a floor derived from the mean (e.g. 0.75*mean) would frequently
    // sit *above* max_length_to_grow, making the two bounds contradictory and
    // rejecting nearly every secondary branch attempt outright. The upfront
    // "max_length_to_grow < kMinGlialProcessLength" rejection above already
    // ensures a non-empty window remains.
    double min_length_to_grow = kMinGlialProcessLength;
    int nbr_non_checked_spheres = factor*3 ;

    // Redraw (not clamp) until the sample actually respects both bounds:
    // clamping every over-cap draw to exactly max_length_to_grow would pile up
    // an artificial spike of branches at exactly the parent's own length,
    // rather than a natural distribution of shorter-than-parent lengths.
    int count = 0;
    while (count < 100 && (length_to_grow < min_length_to_grow || length_to_grow > max_length_to_grow)) {
        length_to_grow = length_dist(generator);
        count++;
    }

    if (length_to_grow < min_length_to_grow || length_to_grow > max_length_to_grow) {
        return false;
    }

    Eigen::Vector3d vector_to_prev_sphere = (random_sphere.center - glial_cell_to_grow.ramification_spheres[random_branch][random_sphere_ind - 1].center).normalized();
    Sphere first_sphere;
    Eigen::Vector3d attractor = Eigen::Vector3d(0, 0, 0);
    double initial_radius = random_sphere.radius;
    //cout << "Initial radius for new branch: " << initial_radius << endl;
    bool first_sphere_created = GenerateFirstSphereinProcess(first_sphere, attractor, initial_radius, random_sphere, vector_to_prev_sphere, nbr_spheres, nbr_spheres_between, glial_cell_to_grow.id, nbr_branches, false);
    //cout << "First sphere created at position: " << first_sphere.center.transpose() << " with radius: " << first_sphere.radius << endl;
    // Cap the decay length at the mean (not the sampled) target, same reasoning
    // as growPrimaryBranch: length_to_grow (a remaining-budget sample) can end
    // up larger than typical, and actual growth often stops well short of it,
    // leaving radius barely decayed if alpha is calibrated to the inflated value.
    double alpha = -std::log(glial_cell_to_grow.minimum_radius / initial_radius)/std::max(length_to_grow, 20.0);

    auto compute_radius = [&](double t) {
        return std::exp(-alpha * t) * initial_radius;
    };
    

    if (!first_sphere_created) {
        return false;
    }

    std::vector<Sphere> vector_first_spheres;
    if (factor > 1) {
        // add spheres between the first and the last
        vector_first_spheres = addIntermediateSpheres(random_sphere, first_sphere, nbr_branches, nbr_spheres, nbr_spheres_between, compute_radius, 0.0, 0.0);
    } else {
        vector_first_spheres = {first_sphere};
    }

    Glial old_glial_cell = glial_cell_to_grow;
    int current_branch = glial_cell_to_grow.ramification_spheres.size();
    if (glial_cell_to_grow.ramification_spheres.size() <= current_branch) {
        glial_cell_to_grow.ramification_spheres.resize(current_branch + 1);
        glial_cell_to_grow.lengths_branches.resize(current_branch + 1);
        glial_cell_to_grow.lengths_branches[current_branch].resize(vector_first_spheres.size());
    }
    glial_cell_to_grow.ramification_spheres[current_branch]= vector_first_spheres;
    glial_cell_to_grow.lengths_branches[current_branch] = std::vector<double>(vector_first_spheres.size(), old_length);
    glial_cell_to_grow.attractors.push_back(attractor);

    // grow glial cell process
    //Eigen::Vector3d expanded_min_limits = min_limits - Eigen::Vector3d({expanded_for_glial_space, expanded_for_glial_space, expanded_for_glial_space});
    //Eigen::Vector3d expanded_max_limits = max_limits + Eigen::Vector3d({expanded_for_glial_space, expanded_for_glial_space, expanded_for_glial_space});

    Eigen::Vector3d prev_pos = first_sphere.center;
    double distance = initial_radius;
    double total_distance = initial_radius + old_length;
    bool can_grow = true;
    int stop_criteria = 0;

    int parent = random_sphere.id;

    while (can_grow && stop_criteria < 1e5) {

        double R_ = compute_radius(distance);

        if (R_ <= glial_cell_to_grow.minimum_radius || distance >= length_to_grow) {
            can_grow = false;

        } 
        else {
            int grow_straight = 0;
            can_grow = AddOneSphere(R_, true, grow_straight, current_branch, glial_cell_to_grow.ramification_spheres.back().size() > nbr_non_checked_spheres, parent, factor);

            // Same resync as growPrimaryBranch: AddOneSphere can add spheres
            // and still return false right after (e.g. lands outside the
            // extended voxel bounds), so lengths_branches must stay synced
            // with ramification_spheres[current_branch] regardless of the
            // return value -- otherwise it silently falls behind, and any
            // later branch attaching partway along *this* one reads past the
            // end of a too-short lengths_branches[current_branch].
            int nbr_spheres_in_branch_before = glial_cell_to_grow.lengths_branches[current_branch].size();
            int nbr_spheres_in_branch_after = glial_cell_to_grow.ramification_spheres[current_branch].size();
            if (nbr_spheres_in_branch_after > nbr_spheres_in_branch_before) {
                double segment = (prev_pos - glial_cell_to_grow.ramification_spheres.back().back().center).norm();
                distance += segment;
                total_distance += segment;
                for (int i = nbr_spheres_in_branch_before; i < nbr_spheres_in_branch_after; i++) {
                    glial_cell_to_grow.lengths_branches[current_branch].push_back(total_distance);
                }
                prev_pos = glial_cell_to_grow.ramification_spheres.back().back().center;
                parent = glial_cell_to_grow.ramification_spheres.back().back().id;
            }
        }

        stop_criteria++;
    }

    if (distance > min_length_to_grow && glial_cell_to_grow.ramification_spheres.back().size() > nbr_non_checked_spheres) {
        finished = true;
        nbr_spheres += glial_cell_to_grow.ramification_spheres.back().size();

        return true;
    }
    else{
        glial_cell_to_grow = std::move(old_glial_cell);
        //cout << "Could not grow secondary branch for glial cell : "<<  glial_cell.id<< endl;
        return false;
    }
}