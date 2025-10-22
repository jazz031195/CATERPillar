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
            const std::vector<Glial>* glial_pop1_,
            const std::vector<Glial>* glial_pop2_,
            const std::vector<Axon>* axons_,
            const Eigen::Vector3d &extended_min_limits_,
            const Eigen::Vector3d &extended_max_limits_,
            const Eigen::Vector3d &min_limits_,
            const Eigen::Vector3d &max_limits_,
            const double &min_radius_, const double &std_dev_) : CellGrowth(axons_, glial_pop1_, glial_pop2_, extended_min_limits_, extended_max_limits_, min_limits_, max_limits_, std_dev_, min_radius_), glial_cell_to_grow(glial_cell_to_grow_){}

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
                can_grow_ = canSpherebePlaced(s) && !collideswithItself(s);
            }
            else{
                can_grow_ = canSpherebePlaced(s);
            }
            if(can_grow_){
                glial_cell_to_grow.ramification_spheres[index_ram_spheres].push_back(s);
                glial_cell_to_grow.update_all_boxes(s);
            }


        }
    }
    sph.id = last_sphere.id + nbr_spheres + 1;
    glial_cell_to_grow.ramification_spheres[index_ram_spheres].push_back(sph); 
    glial_cell_to_grow.update_all_boxes(sph);
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
            //cout << "adding sphere to glial, radius "<< s.radius << ", position : " << s.center << endl;
            s.parent_id = parent;
            add_spheres(s, last_sphere, check_collision_with_branches, factor, i);
        }
        // if is not in inside voxel
        if (!check_borders(extended_min_limits, extended_max_limits, s.center, s.radius)){ 
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



bool GlialCellGrowth::collideswithItself(Sphere &sph){


    // if collides with own soma
    if (glial_cell_to_grow.soma.CollideswithSphere(sph, barrier_tickness)){
        return true;
    }

    // check other branches of same glial cell
    bool collides_with_branches = glial_cell_to_grow.collidesWithItsOwnRamification(sph, barrier_tickness);
    if (collides_with_branches){
        return true;
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
    double std_ = std_dev;
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
        first_sphere = Sphere(nbr_spheres + nbr_spheres_between + 1, cell_id, 1, point, radius, branch_id, sphere_to_emerge_from.id);
        //cout <<"check placement for first sphere at position: " << first_sphere.center.transpose() << " with radius: " << first_sphere.radius << endl;
        if (canSpherebePlaced(first_sphere)) {
            stop = true;
            // check boundaries
            bool is_inside_voxel = check_borders(extended_min_limits, extended_max_limits, first_sphere.center, first_sphere.radius);
            if (!is_inside_voxel) {
                stop = false;
                tries_++;
                //cout << "First sphere is outside voxel boundaries." << endl;
                if (tries_ == max_nbr_tries) return false;
            } 
        } 
        else {
            //cout << "First sphere collides with existing structures." << endl;
            tries_++;
            if (tries_ == max_nbr_tries) return false;
        }
    }
    return true;
}


std::vector<Sphere> GlialCellGrowth::addIntermediateSpheres(const Sphere &random_sphere, const Sphere &first_sphere,  const int &branch_nbr, const int &nbr_spheres, const int &nbr_spheres_between) {
    
    std::vector<Sphere> intermediate_spheres;

    // Direction vector for sphere placement
    Eigen::Vector3d direction = (first_sphere.center - random_sphere.center).normalized();
    double total_distance = (first_sphere.center - random_sphere.center).norm();
    double distance_between_spheres = total_distance / (nbr_spheres_between + 1);

    // Add intermediate spheres
    for (int i = 0; i < nbr_spheres_between; ++i) {
        Eigen::Vector3d position = random_sphere.center + direction * distance_between_spheres * (i + 1);
        double rad = random_sphere.radius +
                        (first_sphere.radius - random_sphere.radius) * (i + 1) / (nbr_spheres_between + 1);

        Sphere next(
            nbr_spheres + i + 1, first_sphere.object_id, first_sphere.object_type, position, rad, branch_nbr,
            random_sphere.id);

        // Check if the sphere can be placed
        if (canSpherebePlaced(next)) {
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
    double length = length_dist(generator);
    double min_length = 0.5*length; // Minimum length to grow a branch
    int tries = 0;
    while (length < min_length && tries < 100) {
        length = length_dist(generator);
    }

    if (length < min_length) {
        return false;
    }

    double initial_radius = glial_cell_to_grow.soma.radius/3;
    double alpha = -std::log(glial_cell_to_grow.minimum_radius / initial_radius)/length;

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
    double first_radius = compute_radius(initial_radius);
    bool first_sphere_created = GenerateFirstSphereinProcess(first_sphere, attractor, first_radius, glial_cell_to_grow.soma, vector_to_prev_center, nbr_spheres, nbr_spheres_between, glial_cell_to_grow.id, j, true);

    if (!first_sphere_created) {
        return false;
    }

    std::vector<Sphere> vector_first_spheres;
    if (factor >1)
    {
        // add spheres between the first and the last
        vector_first_spheres = addIntermediateSpheres(glial_cell_to_grow.soma, first_sphere, j, nbr_spheres, nbr_spheres_between);
    }
    else{
        vector_first_spheres = {first_sphere};
    }

    //cout << "First sphere in primary branch created at position: " << first_sphere.center.transpose() << " with radius: " << first_sphere.radius << endl;

    int nbr_non_checked_spheres = factor * 3;

    // Grow spheres in the branch
    bool can_grow = true;

    Eigen::Vector3d prev_pos = first_sphere.center;

    Glial old_glial_cell = glial_cell_to_grow;

    glial_cell_to_grow.ramification_spheres.resize(j + 1);
    glial_cell_to_grow.lengths_branches.resize(j + 1);
    glial_cell_to_grow.ramification_boxes.resize(j + 1);
    glial_cell_to_grow.attractors.resize(j + 1);
    glial_cell_to_grow.ramification_spheres[j] = vector_first_spheres;
    glial_cell_to_grow.attractors[j] = attractor;
    glial_cell_to_grow.lengths_branches[j] = std::vector<double>(vector_first_spheres.size(), first_sphere.radius);
    //Eigen::Vector3d expanded_space_min = min_limits - Eigen::Vector3d({expanded_for_glial_space, expanded_for_glial_space, expanded_for_glial_space});
    //Eigen::Vector3d expanded_space_max = max_limits + Eigen::Vector3d({expanded_for_glial_space, expanded_for_glial_space, expanded_for_glial_space});

    //cout << "Growing primary branch..." << endl;
    double distance = initial_radius + first_radius;
    while (can_grow) {
 
        // Calculate the radius at the current distance
        double R_ = compute_radius(distance);
        if (R_ <= glial_cell_to_grow.minimum_radius or distance >= length) {
            break;
        } else {
            bool check_collision = glial_cell_to_grow.ramification_spheres[j].size() >= nbr_non_checked_spheres;
            int grow_straight = 0;
            //cout << "Adding sphere with radius: " << R_ << " at distance: " << distance << endl;
            can_grow = AddOneSphere(R_, true, grow_straight, j, check_collision, parent, factor);
        }

        if (can_grow) {
            const auto &new_sphere = glial_cell_to_grow.ramification_spheres[j].back();
            distance += (prev_pos - new_sphere.center).norm();
            prev_pos = new_sphere.center;
            int nbr_spheres_in_branch_before = glial_cell_to_grow.lengths_branches[j].size();
            int nbr_spheres_in_branch_after = glial_cell_to_grow.ramification_spheres.back().size();
            for (int i = nbr_spheres_in_branch_before; i < nbr_spheres_in_branch_after; i++) {
                glial_cell_to_grow.lengths_branches[j].push_back(distance);
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
    while (glial_cell_to_grow.ramification_spheres[random_branch].empty() && nbr_tries < 100) {
        random_branch = branch_dist(rng);
        nbr_tries++;
    }
    if (nbr_tries >= 100) {
        cout << "Could not find a non-empty branch for glial cell : " << glial_cell_to_grow.id << endl;
        return false;
    }

    int size = glial_cell_to_grow.ramification_spheres[random_branch].size();
    //int random_sphere_ind = size / 2 + rand() % (size - size / 2);
    int random_sphere_ind = rand() % size ;
    Sphere random_sphere = glial_cell_to_grow.ramification_spheres[random_branch][random_sphere_ind];
    
    if (glial_cell_to_grow.lengths_branches.size() <= random_branch) {
        cout << "glial_cell.lengths_branches.size() : " << glial_cell_to_grow.lengths_branches.size() << endl;
        cout << "random_branch : " << random_branch << endl;
        cout << "Error: lengths_branches vector is not large enough." << endl;
        assert(0);
    }
    double old_length = glial_cell_to_grow.lengths_branches[random_branch][random_sphere_ind];
    std::random_device rd;
    std::mt19937 generator(rd());

    std::normal_distribution<double> length_dist(mean_process_length - old_length, std_process_length);
    double length_to_grow = length_dist(generator);
    double min_length_to_grow = 0.75*length_to_grow;
    int nbr_non_checked_spheres = factor*3 ;

    int count = 0;
    while (count < 100 && length_to_grow < min_length_to_grow) {
        length_to_grow = length_dist(generator);
        count++;
    }

    if (length_to_grow < min_length_to_grow) {
        return false;
    }

    Eigen::Vector3d vector_to_prev_sphere = (random_sphere.center - glial_cell_to_grow.ramification_spheres[random_branch][random_sphere_ind - 1].center).normalized();
    Sphere first_sphere;
    Eigen::Vector3d attractor = Eigen::Vector3d(0, 0, 0);
    double initial_radius = random_sphere.radius;
    //cout << "Initial radius for new branch: " << initial_radius << endl;
    bool first_sphere_created = GenerateFirstSphereinProcess(first_sphere, attractor, initial_radius, random_sphere, vector_to_prev_sphere, nbr_spheres, nbr_spheres_between, glial_cell_to_grow.id, nbr_branches, false);
    //cout << "First sphere created at position: " << first_sphere.center.transpose() << " with radius: " << first_sphere.radius << endl;
    double alpha = -std::log(glial_cell_to_grow.minimum_radius / initial_radius)/length_to_grow;

    auto compute_radius = [&](double t) {
        return std::exp(-alpha * t) * initial_radius;
    };
    

    if (!first_sphere_created) {
        return false;
    }

    std::vector<Sphere> vector_first_spheres;
    if (factor > 1) {
        // add spheres between the first and the last
        vector_first_spheres = addIntermediateSpheres(random_sphere, first_sphere, nbr_branches, nbr_spheres, nbr_spheres_between);
    } else {
        vector_first_spheres = {first_sphere};
    }

    Glial old_glial_cell = glial_cell_to_grow;
    int current_branch = glial_cell_to_grow.ramification_spheres.size();
    if (glial_cell_to_grow.ramification_spheres.size() <= current_branch) {
        glial_cell_to_grow.ramification_spheres.resize(current_branch + 1);
        glial_cell_to_grow.ramification_boxes.resize(current_branch + 1);
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
            if (can_grow) {
                double segment = (prev_pos - glial_cell_to_grow.ramification_spheres.back().back().center).norm();
                distance += segment;
                total_distance += segment;
                int nbr_spheres_in_branch_before = glial_cell_to_grow.lengths_branches[current_branch].size();
                int nbr_spheres_in_branch_after = glial_cell_to_grow.ramification_spheres.back().size();
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