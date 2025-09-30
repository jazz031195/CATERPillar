#include "axongammadistribution.h"
#include "Glial.h"
#include "grow_axons.h"
#include "grow_glial_cells.h"
#include <algorithm> // std::sort
#include <random>
#include <chrono>
#include <thread>
#include <mutex>
#include <cmath>
#include "threads.h"
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <cassert>
#include <atomic>
#include <iomanip>
#include <functional>


using namespace std;
using namespace Eigen;
using namespace std::chrono;

//* Auxiliare method to split words in a line using the spaces*//
template <typename Out>
void _split_(const std::string &s, char delim, Out result)
{
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim))
    {
        *(result++) = item;
    }
}
std::vector<std::string> _split_line(const std::string &s, char delim)
{
    std::vector<std::string> elems;
    _split_(s, delim, std::back_inserter(elems));
    return elems;
}

AxonGammaDistribution::AxonGammaDistribution(const double &axons_wo_myelin_icvf_, const double &axons_w_myelin_icvf_, const double &glial_pop1_icvf_soma_, const double &glial_pop1_icvf_branches_, const double &glial_pop2_icvf_soma_, const double &glial_pop2_icvf_branches_, const double &a, const double &b,
                                             Eigen::Vector3d &min_l, Eigen::Vector3d &max_l, const double &min_radius_,
                                              const int &regrow_thr_, const double &beading_variation_, const double &std_dev_, const int &ondulation_factor_, const int &factor_, const bool &can_shrink_, const double &cosPhiSquared_, const double &nbr_threads_, const int &nbr_axons_populations_, const int &crossing_fibers_type_, 
                                              const double &mean_glial_pop1_process_length_, const double &std_glial_pop1_process_length_, const double &mean_glial_pop2_process_length_, const double &std_glial_pop2_process_length_,
                                              const double &glial_pop1_radius_mean_, const double &glial_pop1_radius_std_, const double &glial_pop2_radius_mean_, const double &glial_pop2_radius_std_, const bool &glial_pop1_branching_, const bool &glial_pop2_branching_, const int &nbr_primary_processes_pop1_, const int &nbr_primary_processes_pop2_,
                                              const double &c1_, const double &c2_, const double &c3_)
{
    
    std::random_device rd;
    gen.seed(rd());
    target_axons_wo_myelin_icvf = axons_wo_myelin_icvf_;
    target_axons_w_myelin_icvf = axons_w_myelin_icvf_;
    target_axons_icvf = axons_wo_myelin_icvf_ + axons_w_myelin_icvf_;
    target_glial_pop1_soma_icvf = glial_pop1_icvf_soma_;
    target_glial_pop1_branches_icvf = glial_pop1_icvf_branches_;
    target_glial_pop2_soma_icvf = glial_pop2_icvf_soma_;
    target_glial_pop2_branches_icvf = glial_pop2_icvf_branches_;
    nbr_axons_populations = nbr_axons_populations_;

    glial_pop1_soma_icvf = 0.0;
    glial_pop1_branches_icvf = 0.0;
    glial_pop2_soma_icvf = 0.0;
    glial_pop2_branches_icvf = 0.0;
    myelin_icvf = 0.0;
    axons_icvf = 0.0;
    extracellular_icvf = 0.0;

    crossing_fibers_type = crossing_fibers_type_;
    nbr_primary_processes_pop1 = nbr_primary_processes_pop1_;
    nbr_primary_processes_pop2 = nbr_primary_processes_pop2_;

    alpha = a;
    beta = b;
    min_limits = min_l;
    max_limits = max_l;

    cosPhiSquared = cosPhiSquared_;

    min_radius = min_radius_;
    regrow_thr = regrow_thr_;
    beading_variation = beading_variation_;
    std_dev = std_dev_;
    ondulation_factor = ondulation_factor_;

    glial_pop1_radius_mean = glial_pop1_radius_mean_;
    glial_pop1_radius_std = glial_pop1_radius_std_;
    glial_pop2_radius_mean = glial_pop2_radius_mean_;
    glial_pop2_radius_std = glial_pop2_radius_std_;
    glial_pop1_branching = glial_pop1_branching_;
    glial_pop2_branching = glial_pop2_branching_;

    factor = factor_;
    axon_can_shrink = can_shrink_;

    axons.clear();
    glial_pop1.clear();
    glial_pop2.clear();

    total_volume = (max_l[0] - min_l[0]) * (max_l[1] - min_l[1]) * (max_l[2] - min_l[2]);
    nbr_threads = nbr_threads_;

    mean_glial_pop1_process_length = mean_glial_pop1_process_length_;
    std_glial_pop1_process_length = std_glial_pop1_process_length_;
    mean_glial_pop2_process_length = mean_glial_pop2_process_length_;
    std_glial_pop2_process_length = std_glial_pop2_process_length_;

    c1 = c1_;
    c2 = c2_;
    c3 = c3_;
    /*
    cdf = {
        {4, 8, 16, 32, 64, 128}, // Kappas
        {5, 10, 15, 30, 45, 60, 75}, // Angles
        {   // CDF values
            {0.025, 0.095, 0.20, 0.56, 0.80, 0.91, 0.97},
            {0.055, 0.200, 0.39, 0.84, 0.97, 0.99, 1.00},
            {0.110, 0.370, 0.65, 0.98, 1.00, 1.00, 1.00},
            {0.210, 0.610, 0.88, 1.00, 1.00, 1.00, 1.00},
            {0.380, 0.850, 0.99, 1.00, 1.00, 1.00, 1.00},
            {0.620, 0.980, 1.00, 1.00, 1.00, 1.00, 1.00}
        }
    }; // https://www.sciencedirect.com/science/article/pii/S1053811911001376?via%3Dihub
    */

}
void display_progress(double nbr_axons, double number_obstacles)
{
    int cTotalLength = 50;
    double lProgress = nbr_axons / number_obstacles;
    if (lProgress > 1)
    {
        lProgress = 1;
    }
    std::cout << "\r[" <<                                     //'\r' aka carriage return should move printer's cursor back at the beginning of the current line
        string(int(cTotalLength * lProgress), '*') <<         // printing filled part
        string(int(cTotalLength * (1 - lProgress)), '-') <<   // printing empty part
        "] " << nbr_axons << "/" << number_obstacles << endl; // printing percentage
}


// Function to find intersection points between a vector and each face of the cube
Eigen::Vector3d findDistantPoint(const Eigen::Vector3d &vector, const Eigen::Vector3d &point, const double &L)
{
    double distance = L + 10; // Large initial distance
    Eigen::Vector3d distant_point = point + distance * vector;
    return distant_point;
}

std::vector<Sphere> AxonGammaDistribution::addIntermediateSpheres(const Sphere &random_sphere, const Sphere &first_sphere,  const int &branch_nbr, const int &nbr_spheres, const int &nbr_spheres_between) {
    
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
        if (canSpherebePlaced(next, axons, glial_pop1, glial_pop2)) {
            intermediate_spheres.emplace_back(next);
        }
    }

    // Add the initial branching sphere to the list
    intermediate_spheres.push_back(first_sphere);

    return intermediate_spheres;
}

bool AxonGammaDistribution::growSecondaryBranch(Glial &glial_cell, int &nbr_spheres, const double &mean_process_length, const double &std_process_length) {

    //cout << "Growing seondary branch for glial cell : "<<  glial_cell.id<< endl;
    if (glial_cell.ramification_spheres.empty()) {
        cout << "No branches in glial cell" << endl;
        return false;
    }
    if (nbr_spheres <= 0) {
        cout << "No spheres in glial cell : " <<  glial_cell.id << endl;
        return false;
    }

    int nbr_branches = glial_cell.ramification_spheres.size();

    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> branch_dist(0, nbr_branches - 1);
    bool finished = false;

    // finding source of branching
    int nbr_spheres_between = factor - 1;


    int random_branch = branch_dist(rng);
    while (glial_cell.ramification_spheres[random_branch].empty()) {
        random_branch = branch_dist(rng);
    }

    int size = glial_cell.ramification_spheres[random_branch].size();
    //int random_sphere_ind = size / 2 + rand() % (size - size / 2);
    int random_sphere_ind = rand() % size ;
    Sphere random_sphere = glial_cell.ramification_spheres[random_branch][random_sphere_ind];
    if (glial_cell.lengths_branches.size() <= random_branch) {
        cout << "glial_cell.lengths_branches.size() : " << glial_cell.lengths_branches.size() << endl;
        cout << "random_branch : " << random_branch << endl;
        cout << "Error: lengths_branches vector is not large enough." << endl;
        assert(0);
    }
    double old_length = glial_cell.lengths_branches[random_branch][random_sphere_ind];
    std::random_device rd;
    std::mt19937 generator(rd());

    std::normal_distribution<double> length_dist(mean_process_length - old_length, std_process_length);
    double length_to_grow = length_dist(generator);
    double min_length_to_grow = 0.5*length_to_grow;
    int nbr_non_checked_spheres = factor*3 ;


    int count = 0;
    while (count < 100 && length_to_grow < min_length_to_grow) {
        length_to_grow = length_dist(generator);
        count++;
    }

    if (length_to_grow < min_length_to_grow) {
        return false;
    }

    Eigen::Vector3d vector_to_prev_sphere = (random_sphere.center - glial_cell.ramification_spheres[random_branch][random_sphere_ind - 1].center).normalized();
    Sphere first_sphere;
    Eigen::Vector3d attractor = Eigen::Vector3d(0, 0, 0);
    double initial_radius = random_sphere.radius;
    bool first_sphere_created = GenerateFirstSphereinProcess(first_sphere, attractor, initial_radius, random_sphere, vector_to_prev_sphere, nbr_spheres, nbr_spheres_between, glial_cell.id, nbr_branches, false);

    double alpha = -std::log(glial_cell.minimum_radius / initial_radius)/length_to_grow;

    auto compute_radius = [&](double t) {
        return std::exp(-alpha * t) * initial_radius;
    };
    
    //length_to_grow = dist_non_checked_spheres + length_to_grow;
    //cout <<" length to grow : " << length_to_grow << endl;


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

    Glial glial_cell_ = glial_cell;
    int current_branch = glial_cell_.ramification_spheres.size();
    if (glial_cell_.ramification_spheres.size() <= current_branch) {
        glial_cell_.ramification_spheres.resize(current_branch + 1);
        glial_cell_.ramification_boxes.resize(current_branch + 1);
        glial_cell_.lengths_branches.resize(current_branch + 1);
        glial_cell_.lengths_branches[current_branch].resize(vector_first_spheres.size());
    }
    glial_cell_.ramification_spheres[current_branch]= vector_first_spheres;
    glial_cell_.lengths_branches[current_branch] = std::vector<double>(vector_first_spheres.size(), old_length);
    glial_cell_.attractors.push_back(attractor);

    // grow glial cell process

    GlialCellGrowth growth(glial_cell_,  glial_pop1, glial_pop2, axons, min_limits, max_limits, min_limits, max_limits, std_dev, glial_cell_.minimum_radius);

    Eigen::Vector3d prev_pos = first_sphere.center;
    double distance = initial_radius;
    double total_distance = initial_radius + old_length;
    bool can_grow = true;
    int stop_criteria = 0;

    int parent = random_sphere.id;

    while (can_grow && stop_criteria < 1e5) {

        double R_ = compute_radius(distance);

        if (R_ <= glial_cell_.minimum_radius || distance >= length_to_grow) {
            can_grow = false;

        } 
        else {
            int grow_straight = 0;
            can_grow = growth.AddOneSphere(R_, true, grow_straight, current_branch, glial_cell_.ramification_spheres.back().size() > nbr_non_checked_spheres, parent, factor);
            if (can_grow) {
                glial_cell_ = std::move(growth.glial_cell_to_grow);
                double segment = (prev_pos - glial_cell_.ramification_spheres.back().back().center).norm();
                distance += segment;
                total_distance += segment;
                int nbr_spheres_in_branch_before = glial_cell_.lengths_branches[current_branch].size();
                int nbr_spheres_in_branch_after = glial_cell_.ramification_spheres.back().size();
                for (int i = nbr_spheres_in_branch_before; i < nbr_spheres_in_branch_after; i++) {
                    glial_cell_.lengths_branches[current_branch].push_back(total_distance);
                }
                prev_pos = glial_cell_.ramification_spheres.back().back().center;
                parent = glial_cell_.ramification_spheres.back().back().id;
            }
        }

        stop_criteria++;
    }

    if (distance > min_length_to_grow && glial_cell_.ramification_spheres.back().size() > nbr_non_checked_spheres) {
        finished = true;
        nbr_spheres += glial_cell_.ramification_spheres.back().size();
        glial_cell = std::move(glial_cell_);
        return true;
    }
    else{
        return false;
    }
}

bool AxonGammaDistribution::GenerateFirstSphereinProcess(Sphere &first_sphere, Eigen::Vector3d &attractor, const double &radius, const Sphere &sphere_to_emerge_from, const Eigen::Vector3d &vector_to_prev_center, const int &nbr_spheres, const int &nbr_spheres_between, const int &cell_id, const int &branch_id, const bool &primary_process) {
    
    bool stop = false;
    int tries_ = 0;
    int max_nbr_tries = 100;
    attractor = Eigen::Vector3d(0, 0, 0);


    while (!stop && tries_ < max_nbr_tries) {
        Eigen::Vector3d vector = {0, 0, 0};
        Eigen::Vector3d point = {0, 0, 0};
        sphere_to_emerge_from.getPointOnSphereSurface(point, vector, vector_to_prev_center, primary_process);
        attractor = findDistantPoint(vector, point, max_limits[0]);
        first_sphere = Sphere(nbr_spheres + nbr_spheres_between + 1, cell_id, 1, point, radius, branch_id, sphere_to_emerge_from.id);
        if (canSpherebePlaced(first_sphere, axons, glial_pop1, glial_pop2)) {
            stop = true;
            // check boundaries
            if (!withinBounds(first_sphere.center, glial_pop1_radius_mean)) {
                stop = false;
            }  
        } else {
            tries_++;
            if (tries_ == max_nbr_tries) return false;
        }
    }
    return true;
}

bool AxonGammaDistribution::growPrimaryBranch(Glial &glial_cell, int &nbr_spheres, const double &mean_primary_process_length, const double &std_primary_process_length) {

    //cout << "Growing primary branch for glial cell : "<<  glial_cell.id << "/"<< glial_pop1.size()<< endl;
    int j = glial_cell.ramification_spheres.size();
    glial_cell.ramification_spheres.resize(j + 1);
    glial_cell.lengths_branches.resize(j + 1);
    glial_cell.ramification_boxes.resize(j + 1);
    glial_cell.attractors.resize(j + 1);

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

    double initial_radius = glial_cell.soma.radius/3;
    double alpha = -std::log(glial_cell.minimum_radius / initial_radius)/length;

    auto compute_radius = [&](double t) {
        return std::exp(-alpha * t) * initial_radius;
    };


    // Find points on the surface of the glial soma sphere
    bool stop = false;
    tries = 0;
    // Add the first sphere with intermediate spheres
    int parent = 0;  // parent is the soma
    int nbr_spheres_between = factor - 1;

    Eigen::Vector3d vector_to_prev_center = {0, 0, 0};
    Sphere first_sphere;
    Eigen::Vector3d attractor = Eigen::Vector3d(0, 0, 0);
    double first_radius = compute_radius(initial_radius);
    bool first_sphere_created = GenerateFirstSphereinProcess(first_sphere, attractor, first_radius, glial_cell.soma, vector_to_prev_center, nbr_spheres, nbr_spheres_between, glial_cell.id, j, true);

    if (!first_sphere_created) {
        return false;
    }

    std::vector<Sphere> vector_first_spheres;
    if (factor >1)
    {
        // add spheres between the first and the last
        vector_first_spheres = addIntermediateSpheres(glial_cell.soma, first_sphere, j, nbr_spheres, nbr_spheres_between);
    }
    else{
        vector_first_spheres = {first_sphere};
    }

    int nbr_non_checked_spheres = factor * 3;

    // Grow spheres in the branch
    bool can_grow = true;

    Eigen::Vector3d prev_pos = first_sphere.center;

    Glial glial_cell_ = glial_cell;
    glial_cell_.ramification_spheres[j] = vector_first_spheres;
    glial_cell_.attractors[j] = attractor;
    glial_cell_.lengths_branches[j] = std::vector<double>(vector_first_spheres.size(), first_sphere.radius);

    GlialCellGrowth growth(glial_cell_, glial_pop1, glial_pop2, axons, min_limits, max_limits, min_limits, max_limits, std_dev, glial_cell_.minimum_radius);
    
    double distance = initial_radius + first_radius;
    while (can_grow) {
 
        // Calculate the radius at the current distance
        double R_ = compute_radius(distance);
        if (R_ <= glial_cell_.minimum_radius or distance >= length) {
            break;
        } else {
            bool check_collision = glial_cell.ramification_spheres[j].size() >= nbr_non_checked_spheres;
            int grow_straight = 0;
            can_grow = growth.AddOneSphere(R_, true, grow_straight, j, check_collision, parent, factor);
        }

        if (can_grow) {
            glial_cell_ = growth.glial_cell_to_grow;
            const auto &new_sphere = glial_cell_.ramification_spheres[j].back();
            distance += (prev_pos - new_sphere.center).norm();
            prev_pos = new_sphere.center;
            int nbr_spheres_in_branch_before = glial_cell_.lengths_branches[j].size();
            int nbr_spheres_in_branch_after = glial_cell_.ramification_spheres.back().size();
            for (int i = nbr_spheres_in_branch_before; i < nbr_spheres_in_branch_after; i++) {
                glial_cell_.lengths_branches[j].push_back(distance);
            }
        }
    }


    // Check if the branch has grown sufficiently
    if (distance < min_length) {
        return false;
    }
    else {
        //cout << "Branch grown" << endl;
        glial_cell = std::move(glial_cell_);
        // Update the number of spheres
        nbr_spheres += glial_cell.ramification_spheres.back().size();
        
        return true;
    }

}


void AxonGammaDistribution::growFirstPrimaryBranches(Glial &glial_cell, const int &number_ramification_points, int &nbr_spheres, const double &mean_process_length, const double &std_process_length)
{


    int nbr_tries = 0;
    for (int j = 0; j < number_ramification_points; j++)
    {
        bool has_grown = growPrimaryBranch(glial_cell, nbr_spheres, mean_process_length, std_process_length);
        if (!has_grown && nbr_tries < 100){
            j = j - 1;
            nbr_tries += 1;
        }
        else if (nbr_tries >= 1000){
            cout << "Failed to grow glial cell" << endl;
        }
        else{
            nbr_tries = 0;
        }
        
    }


}

bool AxonGammaDistribution::withinBounds(const Eigen::Vector3d &pos, const double &distance)
{
    // Check if the point is inside the dilated box
    for (int i = 0; i < 3; ++i) {
        double min_bound = min_limits[i] - distance;
        double max_bound = max_limits[i] + distance;
        if (pos[i] < min_bound || pos[i] > max_bound) {
            return false; // Point is outside the dilated box
        }
    }
    
    return true; // Point is inside the dilated box
}

bool AxonGammaDistribution::get_begin_end_point(Eigen::Vector3d &Q, Eigen::Vector3d &D, double &angle)
{

    std::uniform_real_distribution<double> udist(0, 1);


    int axis1, axis2, axis3, choice;
    Eigen::Vector3d min_limits_, max_limits_;

    if (nbr_axons_populations == 1){
        axis1 = 0;
        axis2 = 1;
        axis3 = 2;

    }
    else if (nbr_axons_populations == 2){
        choice = rand() % 2;
        if (choice == 0){
            axis1 = 0;
            axis2 = 1;
            axis3 = 2;
        }
        else{
            axis1 = 0;
            axis2 = 2;
            axis3 = 1;
        }
    } 
    else{
        choice = rand() % 3;
        if (choice == 0){
            axis1 = 0;
            axis2 = 1;
            axis3 = 2;
        }
        else if (choice == 1){
            axis1 = 0;
            axis2 = 2;
            axis3 = 1;
        }
        else{
            axis1 = 1;
            axis2 = 2;
            axis3 = 0;
        }
    } 

    if (crossing_fibers_type == 0 && nbr_axons_populations == 2){

        if (nbr_axons_populations == 2){ 
            Eigen::Vector3d min_limits_pop1 = min_limits;
            Eigen::Vector3d max_limits_pop1 = max_limits;
            max_limits_pop1[axis1] = max_limits[axis1]/2;
            Eigen::Vector3d min_limits_pop2 = min_limits;
            min_limits_pop2[axis1] = max_limits[axis1]/2;
            Eigen::Vector3d max_limits_pop2 = max_limits;


            if (choice == 0){
                min_limits_ = min_limits_pop1;
                max_limits_ = max_limits_pop1;

            }
            else{
                min_limits_ = min_limits_pop2;
                max_limits_ = max_limits_pop2;
            }
        }

    }
    else{
        min_limits_ = min_limits;
        max_limits_ = max_limits;
    }  

    double t = udist(gen);
    double axis1_pos = (t * (max_limits_[axis1])) + (1 - t) * (min_limits_[axis1]);
    t = udist(gen);
    double axis2_pos = (t * (max_limits_[axis2])) + (1 - t) * (min_limits_[axis2]);

    Q = {min_limits[axis3], min_limits[axis3], min_limits[axis3]};
    Q[axis1] = axis1_pos;
    Q[axis2] = axis2_pos;

    D = {max_limits[axis3], max_limits[axis3], max_limits[axis3]};
    D[axis1] = axis1_pos;
    D[axis2] = axis2_pos;
    bool outside_voxel = false;
    if (cosPhiSquared != 1.0){
        D = randomPointOnPlane(Q, D, axis1, axis2, axis3, angle, outside_voxel);
        // angle between the axon and the z-axis
        Eigen::Vector3d z_axis = {0, 0, 1};
        Eigen::Vector3d v = D - Q;
        double cosPhi = v.dot(z_axis) / (v.norm() * z_axis.norm());
        double a = acos(cosPhi);
        if (abs(a-angle)>0.02) {
            cout << "Error in angle calculation" << endl;
            cout << "Angle: " << angle << " a: " << a << endl;
            assert(0);
        }
        return !outside_voxel;
        
    }

    return true;
    
}

double AxonGammaDistribution::myelin_thickness(const double &inner_radius){
    return (c1 + c2 * 2.0 * inner_radius + c3 * log(2.0 * inner_radius));
}  
// Set substrate attributes
void AxonGammaDistribution::generate_radii(std::vector<double> &radii_, std::vector<bool> &has_myelin)
{
    axons_w_myelin_icvf = 0.0;
    axons_wo_myelin_icvf = 0.0;
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::gamma_distribution<double> distribution(alpha, beta);

    double icvf_to_reach = target_axons_wo_myelin_icvf + target_axons_w_myelin_icvf;

    if (icvf_to_reach > 0)
    {

        int tried = 0;

        double icvf_ = 0;

        double VolIntra = 0;

        double height = max_limits[2] - min_limits[2];

        while (icvf_ < icvf_to_reach + 0.05)
        {
            if (tried > 1000)
            {
                std::string message = "Radii distribution cannot be sampled [Min. radius Error]\n";
                std::cout << message << std::endl;
                assert(0);
            }
            double jkr = distribution(generator);

            // generates the radii in a list
            if (jkr > min_radius)
            {
                if (target_axons_w_myelin_icvf > 0){   

                    if (axons_w_myelin_icvf < target_axons_w_myelin_icvf){
                        double thickness = myelin_thickness(jkr);
                        jkr = jkr + thickness;
                        axons_w_myelin_icvf += (jkr * jkr * M_PI * height)/total_volume;
                        has_myelin.push_back(true);
                    } 
                    else if (axons_wo_myelin_icvf < target_axons_wo_myelin_icvf){
                        axons_wo_myelin_icvf += (jkr * jkr * M_PI * height)/total_volume;
                        has_myelin.push_back(false);
                    }
                    else{
                        break;
                    }
                }
                else{
                    if (axons_wo_myelin_icvf < target_axons_wo_myelin_icvf){
                        axons_wo_myelin_icvf += (jkr * jkr * M_PI * height)/total_volume;
                        has_myelin.push_back(false);
                    }
                    else{
                        break;
                    }
                }  
                radii_.push_back(jkr);
                tried = 0;
                icvf_ = axons_wo_myelin_icvf+ axons_w_myelin_icvf;

            }
            else
            {
                tried += 1;
            }
        }

            // Create a vector of indices
        std::vector<size_t> indices(radii_.size());
        for (size_t i = 0; i < indices.size(); ++i) {
            indices[i] = i;
        }

        // Sort indices based on the values in radii_
        std::sort(indices.begin(), indices.end(), [&radii_](size_t i1, size_t i2) {
            return radii_[i1] > radii_[i2];
        });

        // Create temporary vectors to hold the sorted values
        std::vector<double> sorted_radii_(radii_.size());
        std::vector<bool> sorted_bools(has_myelin.size());

        for (size_t i = 0; i < indices.size(); ++i) {
            sorted_radii_[i] = radii_[indices[i]];
            sorted_bools[i] = has_myelin[indices[i]];
        }

        // Assign the sorted values back to the original vectors
        radii_ = sorted_radii_;
        has_myelin = sorted_bools;

        if (!radii_.empty()) {
            max_radius = radii_[0];
        } else {
            max_radius = 0; // Or some default value
        }
        // cout << "Maximum radius :" << max_radius << endl;

        std::cout << "Number of axons :" << radii_.size() << endl;
        std::cout <<"initial icvf axons without myelin :" << axons_wo_myelin_icvf << endl;
        std::cout <<"initial icvf axons with myelin :" << axons_w_myelin_icvf << endl;
        axons_wo_myelin_icvf = 0.0;
        axons_w_myelin_icvf = 0.0;

    }

}

bool AxonGammaDistribution::PlaceAxon(const int &axon_id, const double &radius_for_axon, const Eigen::Vector3d &Q, const Eigen::Vector3d &D, std::vector<Axon> &new_axons, const bool &has_myelin, const double &angle_, const bool &outside_voxel)
{

    Axon ax = Axon(axon_id, Q, D, radius_for_axon, beading_variation, has_myelin, angle_, outside_voxel); // axons for regrow batch

    if (has_myelin) {
        double inner_radius = findInnerRadius(radius_for_axon);
        ax.inner_radius = inner_radius;
    }

    Sphere sphere = Sphere(0, ax.id, 0, Q, radius_for_axon);

    bool no_overlap_axons = canSpherebePlaced(sphere, axons, glial_pop1, glial_pop2);
    bool no_overlap_new_axons = canSpherebePlaced(sphere, new_axons, glial_pop1, glial_pop2);

    if (no_overlap_new_axons && no_overlap_axons)
    {
        ax.add_sphere(sphere);
        new_axons.push_back(ax);
        return true;
    }
    else // after comparing with all axons
    {
        return false;
    }
}

bool AxonGammaDistribution::collideswithOtherBranches(const Sphere &sph, const Glial &glial_cell_to_grow)
{
    // if collides with own soma
    if (glial_cell_to_grow.soma.CollideswithSphere(sph, barrier_tickness))
    {
        // cout << "       collides with own soma" << endl;
        return true;
    }

    // check other branches of same glial cell
    for (long unsigned int i = 0; i < glial_cell_to_grow.ramification_spheres.size(); i++)
    {
        if (glial_cell_to_grow.ramification_spheres[i].size() > 0)
        {
            if (glial_cell_to_grow.ramification_spheres[i][0].branch_id != sph.branch_id)
            {
                std::vector<Sphere> branch = glial_cell_to_grow.ramification_spheres[i];
                for (long unsigned int k = 0; k < branch.size(); k++)
                {
                    Sphere sph_ = branch[k];
                    if (sph_.CollideswithSphere(sph, barrier_tickness))
                    {
                        if (k > 20)
                        {
                            // cout << "       collides with sphere k :" << k<< " branch : " << sph_.branch_id<< endl;
                            return true;
                        }
                    }
                }
            }
        }
    }

    return false;
}

// Function to check if a point is inside a dilated box
bool AxonGammaDistribution::check_borders(const Eigen::Vector3d&  min_l, const Eigen::Vector3d&  max_l, const Eigen::Vector3d& pos, const double& distance_to_border) {

    
    // Check if the point is inside the dilated box
    for (int i = 0; i < 3; ++i) {
        double min_bound = min_l[i] - distance_to_border;
        double max_bound = max_l[i] + distance_to_border;
        if (pos[i] < min_bound || pos[i] > max_bound) {
            return false; // Point is outside the dilated box
        }
    }
    
    return true; // Point is inside the dilated box
}


bool AxonGammaDistribution::canSpherebePlaced(const Sphere &sph, const std::vector<Axon> &axs, const std::vector<Glial> &astros, const std::vector<Glial> &oligos
) 
{

    for (int i = 0 ; i < axs.size(); i++) {

        if (!(axs[i].id == sph.object_id && sph.object_type == 0))
        {
            // Check overlap
            if (axs[i].isSphereInsideAxon(sph)) 
            {
                
                return false;
            }
        }
    }

    std::vector<Glial> glial_cells = astros;    
    glial_cells.insert(glial_cells.end(), oligos.begin(), oligos.end());
    // check collision other glial cells 
    for (auto &glial : glial_cells)
    {
        if (glial.collides_with_GlialCell(sph, barrier_tickness)){
            return false;
        }
    }
    return true;
}


// Main function to perform the analysis with parallel threads
void AxonGammaDistribution::ICVF(const std::vector<Axon> &axs, const std::vector<Glial> &glial_pop1, const std::vector<Glial> &oligos) {


    axons_w_myelin_icvf = 0.0;
    axons_wo_myelin_icvf = 0.0;
    for (const auto &axon : axs) {
        if (axon.myelin_sheath) {
            axons_w_myelin_icvf += axon.volume;
        } else {
            axons_wo_myelin_icvf += axon.volume;
        }
    }
    glial_pop1_branches_icvf = 0.0;
    for (const auto &glial : glial_pop1) {
        glial_pop1_branches_icvf += glial.volume_processes;
    }

    glial_pop1_soma_icvf = 0.0;
    for (const auto &glial : glial_pop1) {
        glial_pop1_soma_icvf += glial.volume_soma;
    }
    glial_pop2_branches_icvf = 0.0;
    for (const auto &glial : oligos) {
        glial_pop2_branches_icvf += glial.volume_processes;
    }

    glial_pop2_soma_icvf = 0.0;
    for (const auto &glial : oligos) {
        glial_pop2_soma_icvf += glial.volume_soma;
    }


    axons_w_myelin_icvf = axons_w_myelin_icvf / total_volume;
    axons_wo_myelin_icvf = axons_wo_myelin_icvf / total_volume;
    axons_icvf = axons_w_myelin_icvf + axons_wo_myelin_icvf;
    glial_pop1_branches_icvf = glial_pop1_branches_icvf / total_volume;
    glial_pop1_soma_icvf = glial_pop1_soma_icvf / total_volume;
    glial_pop2_branches_icvf = glial_pop2_branches_icvf / total_volume;
    glial_pop2_soma_icvf = glial_pop2_soma_icvf / total_volume;
    extracellular_icvf = 1 - (axons_w_myelin_icvf + axons_wo_myelin_icvf + glial_pop1_branches_icvf + glial_pop1_soma_icvf + glial_pop2_branches_icvf + glial_pop2_soma_icvf);

}

// Interpolation helper function
double interpolate(double x, const std::vector<double>& xs, const std::vector<double>& ys) {
    if (x <= xs.front()) {
        return ys.front();
    }
    if (x >= xs.back()) {
        return ys.back();
    }

    auto it = std::lower_bound(xs.begin(), xs.end(), x);
    size_t idx = std::distance(xs.begin(), it) - 1;

    double x0 = xs[idx];
    double x1 = xs[idx + 1];
    double y0 = ys[idx];
    double y1 = ys[idx + 1];

    return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}

// Numerical integration using the trapezoidal rule
double integrate(std::function<double(double)> f, double a, double b, int n = 1000) {
    double h = (b - a) / n; // Step size
    double sum = 0.5 * (f(a) + f(b)); // Endpoints contribution
    for (int i = 1; i < n; ++i) {
        sum += f(a + i * h);
    }
    return sum * h;
}

// Helper function to compute erfi (imaginary error function)
double erfi(double x) {
    auto erf_integrand = [](double t) { return std::exp(t * t); };
    double integral_value = integrate(erf_integrand, 0, x, 1000); // Numerical integration
    return (2 / std::sqrt(M_PI)) * integral_value;
}


// Computes c2(kappa) for the axial Watson distribution (kappa >= 0).
// Safe across kappa=0 and large kappa.
static inline double c2_of_kappa(double kappa) {
    if (kappa <= 0.0) return 1.0/3.0;
    const double rt = std::sqrt(kappa);

    // F = (sqrt(pi)/2) * e^{-kappa} * erfi(sqrt(kappa))
    // so that 1/(2*sqrt(kappa)*F) = e^{kappa}/(sqrt(pi)*erfi(sqrt(kappa))*sqrt(kappa))
    const double erfi_rt = erfi(rt);               // assume you have a stable erfi
    const double emk     = std::exp(-kappa);
    const double F       = 0.5 * std::sqrt(M_PI) * emk * erfi_rt;

    // c2 = 1/(2*sqrt(kappa)*F) - 1/(2*kappa)
    const double term1 = 1.0 / (2.0 * rt * F);
    const double term2 = 1.0 / (2.0 * kappa);
    return term1 - term2;
}

// Invert c2 to kappa by monotone bisection on [0, kappa_max].
double AxonGammaDistribution::c2toKappa(double c2_target,
                                        double c2_tol =1e-6,
                                        double kappa_max=64) {
    // Clamp noisy inputs
    if (c2_target <= 1.0/3.0) return 0.0;

    // For near-perfect alignment, use large finite kappa (or caller-specific cap)
    if (c2_target >= 1.0 - 1e-12) {
        // Large-kappa asymptotic: c2 ≈ 1 - 1/(2*kappa)  =>  kappa ≈ 1/(2*(1-c2))
        double guess = 1.0 / std::max(2.0 * (1.0 - c2_target), 1e-12);
        return std::min(guess, kappa_max);
    }

    // Ensure the bracket [lo, hi] contains the solution
    double lo = 0.0;
    double hi = std::min(kappa_max, 1.0 / std::max(2.0 * (1.0 - c2_target), 1e-8)); // asymptotic-based hi
    double c2_lo = 1.0/3.0;           // c2(0) = 1/3
    double c2_hi = c2_of_kappa(hi);

    // If hi is not high enough, expand exponentially until c2_hi >= target or we hit kappa_max
    while (c2_hi < c2_target && hi < kappa_max) {
        lo = hi; c2_lo = c2_hi;
        hi = std::min(hi * 2.0, kappa_max);
        c2_hi = c2_of_kappa(hi);
        if (hi >= kappa_max && c2_hi < c2_target) return kappa_max; // saturated
    }

    // Bisection
    for (int it = 0; it < 100; ++it) {
        double mid   = 0.5 * (lo + hi);
        double c2mid = c2_of_kappa(mid);

        // Check tolerance in c2-space (more meaningful than kappa-space)
        if (std::abs(c2mid - c2_target) <= c2_tol) return mid;

        if (c2mid < c2_target) { lo = mid; c2_lo = c2mid; }
        else                   { hi = mid; c2_hi = c2mid; }
    }
    return 0.5 * (lo + hi); // fallback
}

// Generic clamp: keeps x within [low, high].
template <typename T>
constexpr const T& clamp(const T& x, const T& low, const T& high) {
    return (x < low) ? low : (x > high) ? high : x;
}

// Sample theta in [0, pi/2] from axial Watson(kappa >= 0)
double AxonGammaDistribution::draw_angle(double kappa) {
    std::uniform_real_distribution<double> U(0.0, 1.0);
    double u = U(gen);

    if (kappa <= 0.0) {
        // isotropic axial => mu ~ U(0,1)
        double mu = u;
        return std::acos(mu);
    }

    const double rt = std::sqrt(kappa);
    const double erfi_rt = erfi(rt);                   // same erfi you used before
    const double norm_c = 2.0*rt / (std::sqrt(M_PI)*erfi_rt); // F'(mu) coefficient

    // Initial guess: use small/large-kappa heuristics
    double mu = (kappa < 1.0)
              ? u                       // near-uniform
              : std::min(1.0, std::sqrt(std::max(0.0, std::log1p((erfi_rt*u)/(erfi_rt*(1.0-u))) / kappa))); // crude

    // Newton solve F(mu)=u, clamp mu to [0,1]
    for (int it = 0; it < 20; ++it) {
        double Fmu = erfi(rt*mu) / erfi_rt;
        double fmu = norm_c * std::exp(kappa*mu*mu);
        double step = (Fmu - u) / std::max(1e-16, fmu);
        mu = clamp(mu - step, 0.0, 1.0);
        if (std::abs(step) < 1e-12) break;
    }
    return std::acos(mu);
}


Eigen::Vector3d AxonGammaDistribution::randomPointOnPlane(const Eigen::Vector3d &begin, const Eigen::Vector3d &end, const int &axis1, const int &axis2, const int &axis3, double &angle, bool &outside_voxel) {
 
    double L = (end - begin).norm();

    double phi = angle;
    outside_voxel = false;

    // Seed the random number generator with a random device
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_real_distribution<double> dist(0, 2*M_PI);

    double d = L * std::tan(phi);
    Eigen::Vector3d new_end = end;

    if (angle > 0.2*M_PI) {
        Eigen::Vector3d center_plane = {0,0,0};
        center_plane[axis1] = (max_limits[axis1]/2+3*min_limits[axis1]/2);
        center_plane[axis2] = (max_limits[axis2]/2+3*min_limits[axis2]/2);
        Eigen::Vector3d vector_to_center = center_plane - begin;
        vector_to_center.normalize();
        double theta = dist(rng);
        new_end[axis1] = end[axis1] + d * std::cos(theta);
        new_end[axis2] = end[axis2] + d * std::sin(theta);
        double cos = (new_end-begin).dot(vector_to_center)/(new_end-begin).norm();
        int tries = 0;
        // try to make direction towards center of the voxel
        while (cos < 0) {
            theta = dist(rng);
            new_end[axis1] = end[axis1] + d * std::cos(theta);
            new_end[axis2] = end[axis2] + d * std::sin(theta);
            cos = (new_end-begin).dot(vector_to_center)/(new_end-begin).norm();
            if (tries > 10){
                new_end = end + d * vector_to_center;
                break;
            }
            tries += 1;
        }
        if(new_end[axis1] < min_limits[axis1] || new_end[axis1] > max_limits[axis1] || new_end[axis2] < min_limits[axis2] || new_end[axis2] > max_limits[axis2]) {
            outside_voxel = true;
        }
    }
    else{
        // Generate a random angle theta in the plane
        double theta = dist(rng);
        new_end[axis1] = end[axis1] + d * std::cos(theta);
        new_end[axis2] = end[axis2] + d * std::sin(theta);
        int tries = 0;
        while(new_end[axis1] < min_limits[axis1] || new_end[axis1] > max_limits[axis1] || new_end[axis2] < min_limits[axis2] || new_end[axis2] > max_limits[axis2]) {
            theta = dist(rng);
            new_end[axis1] = end[axis1] + d * std::cos(theta);
            new_end[axis2] = end[axis2] + d * std::sin(theta);
            tries += 1;
            if (tries > 10){
                outside_voxel = true;
                return new_end;
            }
        }

    }

    return new_end;
}

void calculate_c2(std::vector<double> &angles) {
    double cos_2 = 0.0;
    for (auto angle : angles) {
        cos_2 += std::cos(angle) * std::cos(angle);
    }
    cos_2 = cos_2 / angles.size();
    cout << "cos_2 :" << cos_2 << endl;
}



void AxonGammaDistribution::createBatches(std::vector<double> &radii_, std::vector<Axon> &new_axons, std::vector<bool> &has_myelin, std::vector<double> &angles)
{
    const long int tries_threshold = static_cast<long int>((max_limits[0] - min_limits[0]) * (max_limits[1] - min_limits[1]) * 10);
    new_axons.clear();
    new_axons.reserve(radii_.size());

    std::vector<int> indices_to_erase;
    int axon_index = 0;

    for (size_t i = 0; i < radii_.size(); ++i) {
        bool placed = false;
        int tries = 0;

        const double radius = radii_[i];
        const bool has_myelin_ = has_myelin[i];
        double angle = angles[i];

        while (!placed && tries < tries_threshold) {
            Vector3d Q, D;
            bool outside_voxel = !get_begin_end_point(Q, D, angle);

            placed = PlaceAxon(axon_index, radius, Q, D, new_axons, has_myelin_, angle, outside_voxel);

            if (!placed) {
                ++tries;
            }
        }

        if (tries >= tries_threshold) {
            std::cout << "Could not place an initial sphere on plane, axon discarded" << std::endl;
            indices_to_erase.push_back(static_cast<int>(i));
        } else {
            ++axon_index;
        }
    }

    // Erase in reverse order to avoid index shifting
    for (auto it = indices_to_erase.rbegin(); it != indices_to_erase.rend(); ++it) {
        radii_.erase(radii_.begin() + *it);
        has_myelin.erase(has_myelin.begin() + *it);
        angles.erase(angles.begin() + *it);
    }
}



void AxonGammaDistribution::setBatches(const int &num_axons, std::vector<int> &subsets)
{

    if (num_axons <= nbr_threads)
    {
        num_batches = 1;
        subsets = std::vector<int>(1, num_axons);
    }
    else
    {
        if (num_axons % nbr_threads == 0)
        {
            num_batches = num_axons / nbr_threads;
            subsets = std::vector<int>(num_batches, nbr_threads);
        }
        else
        {
            int left = num_axons % nbr_threads;
            subsets = std::vector<int>(int(num_axons / nbr_threads), nbr_threads); // max capacity axons in all batches except last
            subsets.push_back(left);                                                   // last batch has the extra axons left
            num_batches = subsets.size();                                              // number of batches
        }
    }
}
void AxonGammaDistribution::createSubstrate()
{
    parallelGrowth();
}

std::vector<double> AxonGammaDistribution::generate_angles(const int &num_samples){

    std::vector<double> angles;
    for (int i = 0; i < num_samples; i++)
    {
        double phi = draw_angle(kappa);
        angles.push_back(phi);
    }
    return angles;
}

void AxonGammaDistribution::GrowAllAxons(){
    
    if (cosPhiSquared != 1.0){

        kappa = c2toKappa(cosPhiSquared);
        cout <<"kappa :" << kappa << endl;
    }
    
    // generate radii with gamma distribution
    std::vector<bool> has_myelin = {};
    generate_radii(radii, has_myelin);
    std::vector<double> angles(radii.size(), 0.0);
    
    if (cosPhiSquared != 1.0){
        angles = generate_angles(radii.size());
    } 

    cout <<"number of threads :" << nbr_threads << endl;

    calculate_c2(angles);

    if (radii.size() != 0){
        // std::cout << "Parallel growth simulation" << endl;

        bool stop = false;

        while (!stop)
        {

            if (radii.size() == 0)
            {
                stop = true;
                break;
            }
            std::vector<int> subsets; // depending on which batch

            createBatches(radii, axons, has_myelin, angles);
            int num_obstacles = radii.size();

            //  divides obstacles into subsets
            setBatches(num_obstacles, subsets);


            // cout << " length of subsets :" << subsets.size() << endl;
            growBatches(radii, subsets, has_myelin, angles); // grows all axons, fills list of stuck radii, deletes empty axons
            std::vector<double> radii_to_regrow = stuck_radii;
            std::vector<int> indices_to_regrow = stuck_indices;

            ICVF(axons, glial_pop1, glial_pop2);
            cout << "ICVF axons :" << axons_icvf << endl;

            
            int nbr_attempts = 0;
            double percentage_swelling = 0.1;
            double minimum_percentage_swelling = 1e-6;

            while (axons_icvf < target_axons_icvf && nbr_attempts < 1000)
            {
                double old_icvf = axons_icvf;
                SwellAxons(percentage_swelling);
                cout << "new ICVF " << axons_icvf  << " old icvf : "<< old_icvf <<" percentage swelling :"<< percentage_swelling<< endl;
                if ((axons_icvf - old_icvf) < 1e-5 && percentage_swelling >= minimum_percentage_swelling)
                {
                    percentage_swelling = percentage_swelling/3;
                }
                else if (percentage_swelling < minimum_percentage_swelling)
                {
                    break;
                }
                nbr_attempts += 1;
            }
            
            
            stop = true;
            
        }
    }

}


bool AxonGammaDistribution::SwellSphere(Sphere &sph, const double &percentage){
    
    double R = sph.radius;
    double R_ = R + percentage*R;
    Sphere sph_ = Sphere(sph.id, sph.object_id, sph.object_type, sph.center, R_);
    // if can be placed in the substrate
    if (canSpherebePlaced(sph_, axons, glial_pop1, glial_pop2))
    {
        sph = sph_;
        return true;
    }
    else{
        return false;
    }
}

void AxonGammaDistribution::SwellAxon(Axon &ax, const double &percentage) {

    // Enqueue tasks for each sphere
    std::vector<Sphere> new_spheres;
    for (size_t i = 0; i < ax.outer_spheres.size(); i++) {
        
        Sphere sph = ax.outer_spheres[i];
        bool can_swell = SwellSphere(sph, percentage);
        if (can_swell){
            new_spheres.push_back(sph);
        }
        else{
            // if cannot swell, push the sphere back to the axon
            new_spheres.push_back(ax.outer_spheres[i]);
        }
    }
    ax.outer_spheres = new_spheres;
}


void AxonGammaDistribution::SwellAxons(const double &percentage){
    
    int nbr_axons_without_spheres = 0;
    for (auto &axon : axons)
    {
        if (axon.outer_spheres.size()== 0)
        {
            nbr_axons_without_spheres += 1;
            continue;
        }
        // compute mean radius of the axon
        double mean_radius = 0.0;
        for (size_t i = 0; i < axon.outer_spheres.size(); i++) {
            mean_radius += axon.outer_spheres[i].radius;
        }
        mean_radius /= axon.outer_spheres.size();
        if (mean_radius + percentage*mean_radius > axon.radius){
            continue;
        }
        SwellAxon(axon, percentage);
        axon.update_Volume(factor, min_limits, max_limits);
        axon.updateBox();
    }
    // check ICVF 
    ICVF(axons, glial_pop1, glial_pop2);

}

// Growing substrate
void AxonGammaDistribution::parallelGrowth()
{
    // place glial cells
    cout << "nbr_threads :" << nbr_threads << endl;
    cout <<"overlapping_factor :" << factor << endl;
    cout << "Place Glial Cells" << endl;
    PlaceGlialCells();
    cout << "Grow all Axons" << endl;
    GrowAllAxons();
    cout << "Grow Myelin" << endl;
    add_Myelin();

    ICVF(axons, glial_pop1, glial_pop2);

    cout << "GrowAllGlialCells" << endl;

    GrowAllGlialCells();
    
    std::vector<double> stuck_radii_; 
    std::vector<int> stuck_indices_;
    //bool cells_ok = FinalCheck(axons, stuck_radii_, stuck_indices_);
    bool cells_ok = true;

    if (!cells_ok)
    {
        cout << "Axons Final check failed" << endl;
    }
    else
    {
        cout << "Axons Final check passed" << endl;
    }
    
    ICVF(axons, glial_pop1, glial_pop2);
}

void AxonGammaDistribution::GrowAllGlialCells() {

    if (glial_pop1.size() <= 0.0 && glial_pop2.size() <= 0.0)
    {
        return;
    } 
    if (target_glial_pop1_branches_icvf <= 0.0 && target_glial_pop2_branches_icvf <= 0.0)
    {
        return;
    }

    if (target_glial_pop1_branches_icvf >0.0 && glial_pop1.size() > 0)
    {   
        // Growing extra branches for glial_pop1
        growBranches(glial_pop1, 1);
    } 
    if (target_glial_pop2_branches_icvf >0.0 &&  glial_pop2.size() > 0)
    {
        // Growing extra branches for glial_pop2
        growBranches(glial_pop2, 2);
    }
}

void AxonGammaDistribution::ThreadGrowthGlialCells(std::vector<Glial>& glial_cell_list, const double &mean_process_length, const double &std_process_length, std::vector<int> &nbr_spheres, const int &max_nbr_processes) {
    
    //ThreadPool pool(nbr_threads);
    //std::vector<std::future<void>> futures; // Store futures for synchronization

    for (size_t i = 0; i < glial_cell_list.size(); ++i) {

        //if (glial_cell_list[i].allow_branching && glial_cell_list[i].ramification_spheres.size() < max_nbr_processes){
        if (glial_cell_list[i].allow_branching){
            //futures.emplace_back(pool.enqueueTask(
            //    [this, &glial_cell_list, &nbr_spheres, &mean_process_length, &std_process_length, i]() {
            //        growSecondaryBranch(glial_cell_list[i], nbr_spheres[i], mean_process_length, std_process_length);
            //    }
            //));
            //cout <<"growing secondary branch for glial cell " << glial_cell_list[i].id << endl;
            growSecondaryBranch(glial_cell_list[i], nbr_spheres[i], mean_process_length, std_process_length);
        }
        //else if (glial_cell_list[i].ramification_spheres.size() < max_nbr_processes){
        else{
            //futures.emplace_back(pool.enqueueTask(
            //    [this, &glial_cell_list, &nbr_spheres, &mean_process_length, &std_process_length, i]() {
            //        growPrimaryBranch(glial_cell_list[i], nbr_spheres[i], mean_process_length, std_process_length);
            //    }
            //));
            //cout <<"growing primary branch for glial cell " << glial_cell_list[i].id << endl;
            growPrimaryBranch(glial_cell_list[i], nbr_spheres[i], mean_process_length, std_process_length);
        }
    }
    /*
    for (auto& future : futures) {
        future.get();
    }
    
    for (size_t i = 0; i < glial_cell_list.size(); ++i) {
        auto& glial_cell1 = glial_cell_list[i];
        auto& spheres1 = glial_cell1.ramification_spheres.back();

        for (size_t j = i + 1; j < glial_cell_list.size(); ++j) {
            auto& glial_cell2 = glial_cell_list[j];
            auto& spheres2 = glial_cell2.ramification_spheres.back();

            for (auto& sphere1 : spheres1) {
                for (auto& sphere2 : spheres2) {
                    if (sphere1.CollideswithSphere(sphere2, barrier_tickness)) {

                        glial_cell1.ramification_spheres.back().clear();
                        glial_cell2.ramification_spheres.back().clear();
                        glial_cell1.lengths_branches.pop_back();
                        glial_cell2.lengths_branches.pop_back();
                        glial_cell1.attractors.pop_back();
                        glial_cell2.attractors.pop_back();
                        glial_cell1.update_Box();
                        glial_cell2.update_Box();
                        //cout <<"Collision detected between glial cells " << glial_cell1.id << " and " << glial_cell2.id << endl;

                        // Once collision is handled, skip to the next pair
                        goto next_pair;
                    }
                }
            }
        next_pair:;
        }
    }
    */
    

}
void AxonGammaDistribution::growBranches(std::vector<Glial>& glial_cell_list, const int &population_nbr) {

    double target_branches_icvf = 0.0;
    double current_branches_icvf = 0.0;

    int nbr_primary_processes = 1;

    double mean_process_length = 0.0;
    double std_process_length = 0.0;

    double mean_primary_process_length = 0.0;
    double std_primary_process_length = 0.0;

    int max_nbr_processes = 0.0;
    
    if (population_nbr == 1){
        current_branches_icvf = glial_pop1_branches_icvf;
        target_branches_icvf = target_glial_pop1_branches_icvf;
        nbr_primary_processes = nbr_primary_processes_pop1;
        mean_process_length = mean_glial_pop1_process_length;
        std_process_length = std_glial_pop1_process_length;
    }
    else if (population_nbr == 2){
        current_branches_icvf = glial_pop2_branches_icvf;
        target_branches_icvf = target_glial_pop2_branches_icvf;
        nbr_primary_processes = nbr_primary_processes_pop2;
        mean_process_length = mean_glial_pop2_process_length;
        std_process_length = std_glial_pop2_process_length;
    }

    if (target_branches_icvf == 0.0){
        return;
    }

    // Grow first primary branches
    std::vector<int> nbr_spheres(glial_cell_list.size(), 0);


    if (target_branches_icvf >0.0){
        cout << "Growing first primary branches for glial cells" << endl;
        for (int i = 0; i < glial_cell_list.size(); i++) 
        {
            growFirstPrimaryBranches(glial_cell_list[i], nbr_primary_processes, nbr_spheres[i], mean_process_length, std_process_length);
            glial_cell_list[i].compute_processes_icvf(factor, min_limits, max_limits);
        }
    }

    if (population_nbr == 1) {
        ICVF(axons, glial_cell_list, glial_pop2);
        current_branches_icvf = glial_pop1_branches_icvf;
    }
    else if (population_nbr == 2) {
        ICVF(axons, glial_pop1, glial_cell_list);
        current_branches_icvf = glial_pop2_branches_icvf;
    }

    display_progress(current_branches_icvf, target_branches_icvf);

    if (current_branches_icvf >= target_branches_icvf) {
        return;
    }

    int nbr_tries = 0;
    const int max_tries = 1e5;
    double previous_branches_icvf = current_branches_icvf;
    int tolerance = 0;
    
    while (current_branches_icvf < target_branches_icvf) {

        ThreadGrowthGlialCells(glial_cell_list, mean_process_length, std_process_length, nbr_spheres, max_nbr_processes);

        if (nbr_tries%10 == 0 && nbr_tries >0) {
            // Recompute processes and update ICVF for all glial_pop1 in parallel
            for (size_t i = 0; i < glial_cell_list.size(); ++i) {
                glial_cell_list[i].compute_processes_icvf(factor, min_limits, max_limits);
            }
            // Update global ICVF and display progress
            ICVF(axons, glial_pop1, glial_cell_list);
            if (population_nbr == 1) {
                current_branches_icvf = glial_pop1_branches_icvf;
            }
            else if (population_nbr == 2) {
                current_branches_icvf = glial_pop2_branches_icvf;
            }

            display_progress(current_branches_icvf, target_branches_icvf);

            if (current_branches_icvf >= target_branches_icvf) {
                break;
            }
 
            if (std::abs(current_branches_icvf - previous_branches_icvf) < 1e-10) {
                tolerance++;
            }
            if (tolerance > 100) {
                std::cerr << "Stuck in a loop, stopping growth!" << std::endl;
                target_branches_icvf= current_branches_icvf;
                break;
            }
            
            previous_branches_icvf = current_branches_icvf;
            
        }

        // Break if too many attempts
        if (nbr_tries > max_tries) {
            std::cerr << "Max attempts reached while growing branches!" << std::endl;
            break;
        }
        nbr_tries++;

    }

}
void AxonGammaDistribution::PlaceGlialCells() {
    // --- helpers (local lambdas) ---------------------------------------------
    auto draw_positive_radius = [&](std::normal_distribution<>& dis) {
        double r;
        do { r = dis(gen); } while (r <= 0.0);
        return r;
    };

    // Signed clearance to the box (>= 0 => sphere fully inside by that margin).
    auto signed_box_margin = [&](const Sphere& s) {
        double m = std::numeric_limits<double>::infinity();
        for (int i = 0; i < 3; ++i) {
            double left  = (s.center[i] - min_limits[i]) - s.radius;
            double right = (max_limits[i] - s.center[i]) - s.radius;
            m = std::min(m, std::min(left, right));
        }
        return m;
    };

    auto soma_volume = [](double r) {
        return (4.0 * M_PI * r * r * r) / 3.0;
    };

    // Expand the sampling window to allow placements outside the bounds.
    // If you want a different expansion per population, split this.
    const double expand = std::max(0.0, glial_pop1_radius_mean);

    // --- Population 1 (astrocytes) ------------------------------------------
    std::normal_distribution<> dis_radius_pop1(glial_pop1_radius_mean, glial_pop1_radius_std);

    double glial_icvf_1 = 0.0;
    const int MAX_GLOBAL_FAILS_1 = 20000;
    int global_fail_count_1 = 0;

    std::uniform_real_distribution<> base_x1(min_limits[0] - expand, max_limits[0] + expand);
    std::uniform_real_distribution<> base_y1(min_limits[1] - expand, max_limits[1] + expand);
    std::uniform_real_distribution<> base_z1(min_limits[2] - expand, max_limits[2] + expand);

    while (glial_icvf_1 < target_glial_pop1_soma_icvf && global_fail_count_1 < MAX_GLOBAL_FAILS_1) {
        Sphere s;
        bool placed = false;
        bool fully_inside = false;

        for (int attempt = 0; attempt < 10000; ++attempt) {
            const double rad = draw_positive_radius(dis_radius_pop1);

            s = Sphere(
                /*object_id*/ 0,                        // keep your ID scheme if needed
                /*index*/ static_cast<int>(glial_pop1.size()),
                /*type*/ 1,
                { base_x1(gen), base_y1(gen), base_z1(gen) },
                rad
            );

            if (canSpherebePlaced(s, axons, glial_pop1, glial_pop2)) {
                fully_inside = (signed_box_margin(s) >= 0.0);
                placed = true;
                break;
            }
        }

        if (!placed) {
            ++global_fail_count_1;
            continue; // try placing another cell; loop terminates via MAX_GLOBAL_FAILS_1
        }

        // Always add to the population
        Glial glial_cell = Glial(s.object_id, s, glial_pop1_branching);
        glial_pop1.push_back(glial_cell);

        // Count ICVF only if fully inside
        if (fully_inside) {
            glial_icvf_1 += soma_volume(glial_cell.soma.radius) / total_volume;
        }
    }

    glial_pop1_soma_icvf = glial_icvf_1;

    // --- Population 2 (oligodendrocytes) ------------------------------------
    std::normal_distribution<> dis_radius_pop2(glial_pop2_radius_mean, glial_pop2_radius_std);

    double glial_icvf_2 = 0.0;
    const int MAX_GLOBAL_FAILS_2 = 20000;
    int global_fail_count_2 = 0;

    std::uniform_real_distribution<> base_x2(min_limits[0] - expand, max_limits[0] + expand);
    std::uniform_real_distribution<> base_y2(min_limits[1] - expand, max_limits[1] + expand);
    std::uniform_real_distribution<> base_z2(min_limits[2] - expand, max_limits[2] + expand);

    while (glial_icvf_2 < target_glial_pop2_soma_icvf && global_fail_count_2 < MAX_GLOBAL_FAILS_2) {
        Sphere s;
        bool placed = false;
        bool fully_inside = false;

        for (int attempt = 0; attempt < 10000; ++attempt) {
            const double rad = draw_positive_radius(dis_radius_pop2);

            s = Sphere(
                /*object_id*/ 0,
                /*index*/ static_cast<int>(glial_pop2.size()),
                /*type*/ 1,
                { base_x2(gen), base_y2(gen), base_z2(gen) },
                rad
            );

            if (canSpherebePlaced(s, axons, glial_pop1, glial_pop2)) {
                fully_inside = (signed_box_margin(s) >= 0.0);
                placed = true;
                break;
            }
        }

        if (!placed) {
            ++global_fail_count_2;
            continue;
        }

        // Always add to the population
        Glial glial_cell = Glial(s.object_id, s); // pop2 ctor (no branching arg)
        glial_pop2.push_back(glial_cell);

        // Count ICVF only if fully inside
        if (fully_inside) {
            glial_icvf_2 += soma_volume(glial_cell.soma.radius) / total_volume;
        }
    }

    glial_pop2_soma_icvf = glial_icvf_2;
}

bool AxonGammaDistribution::FinalCheck(std::vector<Axon> &axs, std::vector<double> &stuck_radii_, std::vector<int> &stuck_indices_)
{
    // std::cout << "--- Final Check ---" << endl;
    std::vector<Axon> final_axons;
    bool not_collide;

    for (long unsigned int j = 0; j < axs.size(); j++)
    {
        bool all_spheres_can_be_placed = true;
        for (long unsigned int i = 0; i < axs[j].outer_spheres.size(); i++)
        { // for all spheres
            if (!canSpherebePlaced(axs[j].outer_spheres[i], axs, glial_pop1, glial_pop2))
            {
                std::cout << " Axon :" << axs[j].id << ", sphere : " << axs[j].outer_spheres[i].id << " collides with environment !" << endl;
                all_spheres_can_be_placed = false;
                break;
            }
        }
        if (all_spheres_can_be_placed)
        {
            final_axons.push_back(axs[j]);
        }
        else
        {
            stuck_radii_.push_back(axs[j].radius);
            stuck_indices_.push_back(axs[j].id);
        }
    }
    if (final_axons.size() == axs.size())
    {
        //std::cout << " No Axon collides with environment !" << endl;
        not_collide = true;
    }
    else
    {
        //std::cout << " Axon COLLIDE with environment !" << endl;
        not_collide = false;
        axs.clear();
        axs = final_axons;
    }

    return not_collide;
}
bool AxonGammaDistribution::SanityCheck(std::vector<Axon>& growing_axons,
                                        std::vector<double>& stuck_radii_,
                                        std::vector<int>& stuck_indices_) {
    std::unordered_set<int> collided_ids;
    std::vector<Axon> axons_to_check_collision_with = growing_axons;

    if (growing_axons.size() <=1) {
        return true;  // No axons to check
    }

    int nbr_empty_spheres = 0;
    for (auto& axon : growing_axons) {
        if (axon.outer_spheres.empty()) {
            nbr_empty_spheres++;
        }
    }
    if (nbr_empty_spheres == growing_axons.size()) {
        return true;  // No axons to check
    }

    // generate random number between 0 and growing_axons.size()
    std::uniform_int_distribution<> dis(0, growing_axons.size() );
    int random_index = dis(gen);

    // check if this axon has spheres
    while (growing_axons[random_index].outer_spheres.size() <= 1) {
        random_index = dis(gen);
    }

    for (size_t i = 0; i < growing_axons.size(); ++i) {  // Start from second axon

        if (i == random_index) continue;  // Skip the randomly selected axon

        const auto& axon = growing_axons[i];

        if (axon.outer_spheres.empty()) continue;

        for (const auto& sphere : axon.outer_spheres) {
            if (!canSpherebePlaced(sphere, axons_to_check_collision_with, glial_pop1, glial_pop2)) {
                collided_ids.insert(axon.id);
                axons_to_check_collision_with.erase(
                    std::remove_if(axons_to_check_collision_with.begin(), axons_to_check_collision_with.end(),
                                   [&axon](const Axon& a) { return a.id == axon.id; }),
                    axons_to_check_collision_with.end()
                );
                break;
            }
        }
    }

    bool all_clear = collided_ids.empty();

    for (auto& axon : growing_axons) {
        if (collided_ids.count(axon.id)) {
            stuck_radii_.push_back(axon.radius);
            stuck_indices_.push_back(axon.id);
            axon.destroy();
            axon.update_Volume(factor, min_limits, max_limits);
            axon.updateBox();
        }
    }

    return all_clear;
}


void update_straight(bool can_grow_, int &grow_straight, int &straight_growths, int ondulation_factor)
{

    if (can_grow_)
    {

        if (grow_straight == 1)
        {
            if (straight_growths >= ondulation_factor) // if axon has been growing straight for a number of spheres in a row
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
            // if the sphere hadn't grown straight previously . set to straight for next "ondulation_factor" spheres
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

void AxonGammaDistribution::growthThread(
    std::vector<Axon>& axs,
    Axon& axon,
    AxonGrowth& growth,
    int& finished,
    int& grow_straight,
    int& straight_growths,
    double& stuck_radius,
    int& stuck_index
) {


    double varied_radius = axon.radius;
    // 1) Possibly vary the radius (beading)
    if (beading_variation > 0){
        varied_radius =RandomradiusVariation(axon);
        
    }


    // 2) Attempt to add a sphere (returns true if successful)
    bool can_grow = growth.AddOneSphere(varied_radius, /*create_sphere=*/true, grow_straight, factor);
    // 3) Sync axon data back
    //    growth.axon_to_grow updated by AddOneSphere
    growth.axon_to_grow.growth_attempts = axon.growth_attempts; 

    // 4) Check if axon is finished
    if (growth.finished) {
        //cout << "axon " << axon.id << " is finished" << endl;
        finished = 1; 
        stuck_radius = -1;
        stuck_index = -1;
    }
    else { 
        // The axon is still growing
        if (!can_grow && axon_can_shrink) {
            // We tried adding a sphere and failed to grow. Attempt to shrink radius.
            bool shrinkOk = shrinkRadius(growth, varied_radius, axon); 
            if (!shrinkOk) {
                // If shrinking didn't help either
                if (axon.growth_attempts < 3) {
                    //cout <<"axon will try again on same spot" << endl;
                    // Retry with a new attempt
                    axon.keep_one_sphere();  // Revert the last sphere or do whatever keep_one_sphere does
                    growth.axon_to_grow = axon;
                    finished = 0;
                    stuck_radius = -1;
                    stuck_index = -1;
                    //cout << "axon " << axon.id << " is stuck, growth attempts : "<<  axon.growth_attempts<< " center : (" << axon.begin[0] << ", " << axon.begin[1]<< ", "<< axon.begin[2] << ")" << endl;
                }
                else {
                    //cout <<"axon is stuck and will try again somewhere else" << endl;
                    // The axon is stuck
                    axon.destroy(); 
                    finished = 1;
                    stuck_radius = axon.radius;
                    stuck_index = axon.id;
                    //cout <<"try axon in different position" << endl;
                }
            }
            else {
                //cout <<"axon " << axon.id << " shrank successfully" << endl;
                // We shrank successfully and placed a sphere 
                if (growth.finished) {
                    finished = 1;
                }
            }
        }
        else if (!can_grow && !axon_can_shrink) {
            if (axon.growth_attempts < 3) {
                // Retry at same positions 10 times
                axon.keep_one_sphere();  // Revert the last sphere or do whatever keep_one_sphere does
                growth.axon_to_grow = axon;
                finished = 0;
              
            }
            else {
                //cout <<"axon " << axon.id << " is stuck" << endl;
                // The axon is stuck
                axon.destroy(); 
                finished = 1;
                stuck_radius = axon.radius;
                stuck_index = axon.id;

            }
        }
        else if (can_grow) {
            finished = 0;
        }
    }

    // 5) Update "straightness" logic
    update_straight(can_grow, grow_straight, straight_growths, ondulation_factor);
}


void AxonGammaDistribution::growAxon(Axon& axon_to_grow, int &index, double& stuck_radius, int& stuck_index) {
    int finished = 0;
    int grow_straight = 0;
    int straight_growths = 0;

    // adjust growth axis
    /* 
    int i;
    double maxVal = axon_to_grow.begin.cwiseAbs().minCoeff(&i);
    if (i != axon_to_grow.growth_axis) {
        axon_to_grow.growth_axis = i; // set growth axis to the one with the smallest absolute value
        cout << "axon " << axon_to_grow.id << " growth axis changed to " << axon_to_grow.growth_axis << endl;

    }
    */


    // Possibly alter the beading amplitude if beading_variation > 0
    if (beading_variation > 0) {

        axon_to_grow.beading_amplitude = beading_variation;
        // Clamp beading amplitude between 0 and 1
        if (axon_to_grow.beading_amplitude < 0.0) {
            axon_to_grow.beading_amplitude = 0.0;
        }
        else if (axon_to_grow.beading_amplitude > 1.0) {
            axon_to_grow.beading_amplitude = 1.0;
        }
    }

    // Initialize a Growth object for this axon
    AxonGrowth growth(axon_to_grow, glial_pop1, glial_pop2, axons, min_limits, max_limits, min_limits, max_limits, std_dev, min_radius);

    int tries_  = 0;
    int max_nbr_tries = 1e10; // growth class stops itself anyway
    // Keep trying to grow until finished
    while (finished == 0 && tries_ < max_nbr_tries) {

        if (!axon_to_grow.outer_spheres.empty()) {
            // Attempt a single growth step
            growthThread(axons, axon_to_grow, growth, finished, grow_straight, straight_growths, stuck_radius, stuck_index);

            if (stuck_radius>0) {
                tries_ += 1;
            }
        }
        else {
            //cout <<"finished = 1" << endl;
            cout <<"discarded axon " << axon_to_grow.id << " with radius " << axon_to_grow.radius << endl;
            // Axon is empty => it was discarded or destroyed
            finished = 1;
        }

    }

    if (axon_to_grow.outer_spheres.size()==1) {
        assert(0);
    }
    else if (tries_ == max_nbr_tries) {
        //cout << "axon " << axon_to_grow.id << " is stuck" << endl;
        axon_to_grow.destroy();
        stuck_radius = axon_to_grow.radius;
        stuck_index = axon_to_grow.id;
        //ind_axons_pushed = {};
    }
    if (finished == 1 && !axon_to_grow.outer_spheres.empty()){
        Sphere last_sphere = axon_to_grow.outer_spheres.back();
        if (!(last_sphere.center[axon_to_grow.growth_axis] + last_sphere.radius > max_limits[axon_to_grow.growth_axis])) {
            cout << " not close to target limits, axon " << axon_to_grow.id << " is finished" << endl;
            cout <<"axon_to_grow.growth_axis : " << axon_to_grow.growth_axis << endl;
            //assert(0);
        }
        if (axon_to_grow.outer_spheres.size() < 10){
            cout << "Axon " << axon_to_grow.id << " has too few spheres: " << axon_to_grow.outer_spheres.size() << endl;
            cout <<"growth axis : " << axon_to_grow.growth_axis << endl;
            cout << "starting point : (" << axon_to_grow.begin[0] << ", " << axon_to_grow.begin[1] << ", " << axon_to_grow.begin[2] << ")" << endl;
            int i;
            double maxVal = axon_to_grow.begin.cwiseAbs().minCoeff(&i);
            cout << " growth axis should be : " << i << endl;
            axon_to_grow.destroy();
            stuck_radius = axon_to_grow.radius;
            stuck_index = axon_to_grow.id;
        }
    }
    /*
    else{
        ind_axons_pushed =  growth.ind_axons_pushed;
        for (int i : ind_axons_pushed) {
            axons_pushed.push_back(growth.modified_axons[i]);
        }
    }
    */
}
void AxonGammaDistribution::processBatchWithThreadPool(
    std::vector<Axon>& axons_to_grow,
    std::vector<int>& indices,
    std::vector<double>& stuck_radii,
    std::vector<int>& stuck_indices) 
{
    ThreadPool pool(nbr_threads);

    std::vector<std::future<void>> futures; // Store futures for synchronization

    for (size_t i = 0; i < axons_to_grow.size(); ++i) {
        futures.emplace_back(pool.enqueueTask(
            [this, &axons_to_grow, &stuck_radii, &stuck_indices, &indices, i]() {
                growAxon(axons_to_grow[i], indices[i], stuck_radii[i], stuck_indices[i]);
            }
        ));
    }

    // Wait for all tasks to complete
    for (auto& future : futures) {
        future.get();
    }

}

// Function to check if two vectors share any common element
bool hasCommonElement(const std::vector<int>& a, const std::vector<int>& b) {
    std::unordered_set<int> set_a(a.begin(), a.end());
    for (int val : b) {
        if (set_a.count(val)) return true;
    }
    return false;
}

std::vector<int> removeOverlappingVectors(
    std::vector<std::vector<int>>& intGroups,
    std::vector<std::vector<Axon>>& axonGroups)
{
    std::unordered_set<int> indicesToRemove;

    for (size_t i = 0; i < intGroups.size(); ++i) {
        for (size_t j = i + 1; j < intGroups.size(); ++j) {
            if (indicesToRemove.count(i) || indicesToRemove.count(j)) continue;
            if (hasCommonElement(intGroups[i], intGroups[j])) {
                indicesToRemove.insert(i);
                indicesToRemove.insert(j);
            }
        }
    }

    // Convert to vector and sort in reverse so we can erase safely
    std::vector<int> toErase(indicesToRemove.begin(), indicesToRemove.end());
    std::sort(toErase.rbegin(), toErase.rend());

    for (int idx : toErase) {
        intGroups.erase(intGroups.begin() + idx);
        axonGroups.erase(axonGroups.begin() + idx);
    }

    return toErase;
}

void AxonGammaDistribution::ModifyAxonsStartingPoint(std::vector<int> &stuck_indices){

    std::vector<int> new_stuck_indices = {};
    for (int i = 0; i < axons.size(); i++) {

        if (std::find(stuck_indices.begin(), stuck_indices.end(), axons[i].id) != stuck_indices.end()) {
            //cout <<"Place axon in new position" << endl;
            Vector3d Q, D;
            double angle_ = axons[i].angle;
            Eigen::Vector3d initial_begin = axons[i].begin;
            bool outside_voxel = !get_begin_end_point(Q, D, angle_);
            Sphere sphere = Sphere(0, axons[i].id, 0, Q, axons[i].radius);
            bool no_overlap_axons = canSpherebePlaced(sphere, axons, glial_pop1, glial_pop2);
            int nbr_tries = 0;
            while(!no_overlap_axons && nbr_tries < 1000){
                outside_voxel = !get_begin_end_point(Q, D, angle_);
                sphere = Sphere(0, axons[i].id, 0, Q, axons[i].radius);
                no_overlap_axons = canSpherebePlaced(sphere, axons, glial_pop1, glial_pop2);
                nbr_tries += 1;
            }
            if (no_overlap_axons) {
                axons[i].outer_spheres.clear();
                axons[i].outer_spheres.push_back(sphere);
                axons[i].update_Volume(factor, min_limits, max_limits);
                axons[i].updateBox();
                axons[i].end = D;
                axons[i].outside_voxel = outside_voxel;
                axons[i].angle = angle_;
                axons[i].begin = Q;
                new_stuck_indices.push_back(axons[i].id);
                //cout <<"Placed axon in new position" << endl;
            }
            else{
                new_stuck_indices.push_back(axons[i].id);
                //cout <<"Couldnt place axon in new position" << endl;
            }
        }
    }
    stuck_indices = new_stuck_indices;
}

void AxonGammaDistribution::createBatch(const std::vector<double> &radii_, const std::vector<int> &indices, const int &num_subset, const int &first_index_batch, std::vector<Axon> &new_axons, const std::vector<bool> &has_myelin, std::vector<double> &angles)
{
    long int tries_threshold = (max_limits[0]- min_limits[0]) * (max_limits[1]-min_limits[1]) * 5;
    new_axons.clear();

    std::vector<double>::const_iterator startIterator = radii_.begin() + first_index_batch;             // Start from index
    std::vector<double>::const_iterator stopIterator = radii_.begin() + first_index_batch + num_subset; // Stop when batch is finished
    // radii in batch
    std::vector<double> batch_radii(startIterator, stopIterator);
    double angle_;
    bool outside_voxel = false;

    cout <<"Placing axon" << endl;

    // cout << "batch_radii.size :" << batch_radii.size() << endl;
    for (long unsigned int i = 0; i < batch_radii.size(); ++i) // create axon for each radius
    {
        // std::cout << "radii created :" << i << "/" << batch_radii.size() << endl;
        bool next = false;
        int tries = 0;

        while (!next && tries < tries_threshold)
        {
            
            Vector3d Q, D;
            
            int index_of_axon = indices[i + first_index_batch];
            angle_ = angles[indices[i + first_index_batch]];
            outside_voxel = !get_begin_end_point(Q, D, angle_);
            bool has_myelin_ = has_myelin[indices[i + first_index_batch]];
            if (tries> tries_threshold/2) {
                next = PlaceAxon(index_of_axon, batch_radii[i], Q, D, new_axons, has_myelin_, angle_, outside_voxel);

            }
            else{
                next = PlaceAxon(index_of_axon, batch_radii[i], Q, D, new_axons, has_myelin_, angle_, outside_voxel);
            }
            //cout << "Q : " << Q << endl;
            if (!next)
            {
                tries += 1;
                
            }
        }

    }
    // cout << "new_axons.size :" << new_axons.size() << endl;
}

void AxonGammaDistribution::growBatch(int &number_axons_to_grow, std::vector<double> &radii_,std::vector<int> &indices, std::vector<double> &stuck_radii_, std::vector<int> &stuck_indices_, const int &first_index_batch, std::vector<Axon> &growing_axons, std::vector <bool> &has_myelin, std::vector<double> &angles)
{

    if (first_index_batch+number_axons_to_grow > radii_.size()){
        number_axons_to_grow = radii_.size() - first_index_batch;
    }

    std::vector<Axon> batch_growing_axons(growing_axons.begin() + first_index_batch,
                                      growing_axons.begin() + first_index_batch + number_axons_to_grow);
                                    

    //vector full of -1 with size of number of axons
    std::vector<double> all_stuck_radii(number_axons_to_grow, -1);  

    std::vector<int> all_stuck_indices(number_axons_to_grow, -1);   

    //std::vector<std::vector<int>> ind_axons_pushed;  

    //std::vector<std::vector<Axon>> axons_pushed;   

    //ind_axons_pushed.resize(number_axons_to_grow);
    //axons_pushed.resize(number_axons_to_grow); 

    processBatchWithThreadPool(batch_growing_axons, indices, all_stuck_radii, all_stuck_indices);


    for (int i = 0; i < number_axons_to_grow; i++) // for each axon
    {
        
        if (all_stuck_radii[i] > 0)
        {
            stuck_radii_.push_back(all_stuck_radii[i]);
            stuck_indices_.push_back(all_stuck_indices[i]);
        }
    }

    /*
    
    std::vector<int> ind_growing_axons_to_remove = removeOverlappingVectors(ind_axons_pushed, axons_pushed);
    
    cout <<"number fo duplicates : " << ind_growing_axons_to_remove.size() << endl;

    for (int i = 0; i < ind_growing_axons_to_remove.size(); i++) // for each axon
    {
        // growing axons that pushed the same static axons are labelled as stuck
        stuck_radii_.push_back(all_stuck_radii[ind_growing_axons_to_remove[i]]);
        stuck_indices_.push_back(all_stuck_indices[ind_growing_axons_to_remove[i]]);
        batch_growing_axons[ind_growing_axons_to_remove[i]].destroy();
    }

    std::vector<int> all_indices;
    for (const auto& vec : ind_axons_pushed) {
        all_indices.insert(all_indices.end(), vec.begin(), vec.end());
    }

    // Sort and remove duplicates
    std::sort(all_indices.begin(), all_indices.end());
    all_indices.erase(std::unique(all_indices.begin(), all_indices.end()), all_indices.end());

    cout <<" number of axons pushed : " << all_indices.size() << endl;

    for (const auto& group : axons_pushed) {
        batch_growing_axons.insert(batch_growing_axons.end(), group.begin(), group.end());
    }
    */
    
    //cout <<"number of axons to check in batch_growing_axons : " << batch_growing_axons.size() << endl;


    SanityCheck(batch_growing_axons, stuck_radii_, stuck_indices_);
        
    //cout << stuck_radii_.size() << " axons are stuck" << endl;

    // add axons in batch that worked in axons
    for (int i = 0; i < batch_growing_axons.size(); i++)
    {
        
        int index = batch_growing_axons[i].id;
        // if index not in stuck indices
        if (std::find(stuck_indices_.begin(), stuck_indices_.end(), index) == stuck_indices_.end())
        {
            if (batch_growing_axons[i].outer_spheres.size() == 0)
            {
                //cout << "batch_growing_axons[i].id :"<< batch_growing_axons[i].id<< " has 0 spheres but is not in stuck_indices_" << endl;
                continue;
            }
            for (int j = 0; j < axons.size(); j++)
            {
                if (axons[j].id == index)
                {
                    axons[j] = std::move(batch_growing_axons[i]);
                }
            }     
        }  
    }
    //cout << "Grown axons added to list" << endl;

    /*

    bool no_collision = FinalCheck(axons, stuck_radii_, stuck_indices_);
    if (!no_collision)
    {
        cout << "Final check failed" << endl;
        assert(0);
    }
    else{
        cout << "Final check passed" << endl;
    }
    */
    
    
    
    //cout << "axons.size() after check: " << axons.size() << endl;
}

void AxonGammaDistribution::growBatches(std::vector<double> &radii_, std::vector<int> &subsets_, std::vector<bool> &has_myelin, std::vector<double> &angles)
{

    // indices for axons
    std::vector<int> indices(radii_.size());
    std::iota(indices.begin(), indices.end(), 0);
    
    stuck_radii.clear();

    int nbr_tries = 0;

    for (int j = 0; j < num_batches; j++) // batches of axon growth
    {
        auto startTime = std::chrono::high_resolution_clock::now();
        std::cout << "---   Batch " << j << "   --- " << endl;

        display_progress(j, num_batches);
        // index of first axon in batch
        int first_index_batch = j * nbr_threads;

        vector<int> finished(subsets_[j], 0);         // 0 for false
        vector<int> grow_straight(subsets_[j], 1);    // 1 for true
        vector<int> straight_growths(subsets_[j], 0); // for each axon
        vector<int> shrink_tries(subsets_[j], 0);     // for each axon
        vector<int> restart_tries(subsets_[j], 0);    // for each axon

        std::vector<double> stuck_radii_;
        std::vector<int> stuck_indices_;
        growBatch(subsets_[j], radii_, indices, stuck_radii_, stuck_indices_, first_index_batch, axons, has_myelin, angles);
        int tries_threshold = regrow_thr;
        // tries in 10 times to regrow the stuck axons
        std::vector<double> new_stuck_radii_;
        std::vector<int> new_stuck_indices_;
        double old_stuck_radii_size =  0;
        std::vector<Axon> stuck_axons;
   
        while (stuck_radii_.size() > 0 && nbr_tries < tries_threshold)
        {
            if (nbr_tries > tries_threshold/2) {
                axon_can_shrink = true;
            }
            new_stuck_radii_.clear();
            new_stuck_indices_.clear();
            //cout << " Regrow " << stuck_radii_.size() << " axons" << endl;
            int number_axons_to_grow = stuck_radii_.size();
            ModifyAxonsStartingPoint(stuck_indices_); 

            if (stuck_indices_.size() == 0){
                stuck_radii_.clear();
                break;
            }

            stuck_axons.reserve(stuck_indices_.size());
            for (unsigned i = 0; i < stuck_indices_.size(); i++)
            {
                stuck_axons.push_back(axons[stuck_indices_[i]]);
            }
            growBatch(number_axons_to_grow, stuck_radii_, stuck_indices_, new_stuck_radii_, new_stuck_indices_, 0, stuck_axons, has_myelin, angles);
            stuck_radii_ = new_stuck_radii_;
            stuck_indices_ = new_stuck_indices_;
            if(stuck_radii_.size() == old_stuck_radii_size){
                nbr_tries += 1;
            }
            else{
                nbr_tries = 0;
            }
            old_stuck_radii_size = stuck_radii_.size();
            stuck_axons.clear();
        }

        stuck_radii.clear();
        stuck_indices.clear();
        if (nbr_tries >= tries_threshold)
        {
            for (unsigned i = 0; i < stuck_radii_.size(); i++)
            {
                stuck_radii.push_back(stuck_radii_[i]);
                stuck_indices.push_back(stuck_indices_[i]);
            }
            nbr_tries = 0;
        }

        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime);


    } // end for batches

    int counts = 0;
    for (int i = axons.size() - 1; i >= 0; --i) {
        if (axons[i].outer_spheres.size() <= 1) {
            axons.erase(axons.begin() + i);
            counts++;
        }
    }

    cout <<"Number of axons abandoned : " << counts << endl;

}

double get_axonal_length(Axon axon)
{
    double l = 0;
    if (axon.outer_spheres.size() > 1)
    {
        for (long unsigned int i = 1; i < axon.outer_spheres.size(); i++)
        {
            double dist = (axon.outer_spheres[i - 1].center - axon.outer_spheres[i].center).norm();
            l += dist;
        }
        return l;
    }
    else
    {
        return 0;
    }
}



double AxonGammaDistribution::RandomradiusVariation(Axon &axon) {
    double prev_radius = axon.outer_spheres[axon.outer_spheres.size()-1].radius;


    if (prev_radius < axon.radius-axon.beading_amplitude*axon.radius)
    {
        prev_radius = axon.radius-axon.beading_amplitude*axon.radius;
    }
    else if (prev_radius > axon.radius+axon.beading_amplitude*axon.radius)
    {
        prev_radius = axon.radius+axon.beading_amplitude*axon.radius;
    }

    double standard_deviation = axon.beading_amplitude * axon.radius / 3;

    std::normal_distribution<double> dis(prev_radius, standard_deviation);

    double random_radius = dis(gen);

    while (random_radius > prev_radius+3* standard_deviation || random_radius < prev_radius-3*standard_deviation)
    {
        random_radius = dis(gen);
    }
    
    if (random_radius < min_radius)
    {
        random_radius = min_radius;
    }

    return random_radius;
}
/*
double AxonGammaDistribution::RandomradiusVariation(Axon& axon)   // you pass an id or index for the axon
{

    std::mt19937_64 rng{std::random_device{}()};
    std::normal_distribution<double> N01{0.0, 1.0};
    double corr_len = 5.0;

    double prev_radius = axon.outer_spheres[axon.outer_spheres.size()-1].radius;
    const double R_mean = axon.radius;
    const double A      = axon.beading_amplitude; 
    const double Rmin   = R_mean * (1.0 - A);
    const double Rmax   = R_mean * (1.0 + A);

    // step length
    double ds = prev_radius/factor;
    const double ds_min = 1e-6 * std::max(1.0, R_mean);
    ds = std::max(ds, ds_min);

    // OU/AR(1) update
    const double alpha     = std::exp(-ds / std::max(1e-12, corr_len));
    const double sigma_eps = std::sqrt(std::max(0.0, 1.0 - alpha*alpha));

    double z = axon.z_stat_rad_variation;
    z = alpha * z + sigma_eps * N01(rng);

    // map latent Gaussian into [1-A, 1+A]
    double m = 1.0 + A * std::tanh(z);
    double R = R_mean * m;

    // enforce bounds
    if (R < min_radius) R = min_radius;
    if (R < Rmin) R = Rmin;
    if (R > Rmax) R = Rmax;

    axon.z_stat_rad_variation = z; // store updated latent state

    return R;

}
*/
// Axon growth
double AxonGammaDistribution::radiusVariation(const Axon &axon)
{
    double mean_radius = axon.radius;
    double length = get_axonal_length(axon);

    double amplitude = mean_radius * axon.beading_amplitude;
    double beading_period = 1;

    double omega = 2 * M_PI / (beading_period*mean_radius);

    double lambda = - omega *axon.phase_shift*beading_period*mean_radius;

    double r = amplitude * sin(omega * length + lambda) + mean_radius;

    
    if (r < min_radius)
    {
        r = min_radius;
    }


    return r;
}

bool AxonGammaDistribution::shrinkRadius(AxonGrowth &growth, const double &radius_to_shrink, Axon &axon)
{

    bool can_grow_;
    // find position that works for smallest radius
    can_grow_ = false;
    double rad = radius_to_shrink;
    double initial_rad = radius_to_shrink;
    bool can_grow_min_rad = growth.AddOneSphere(min_radius, false, 0, factor);
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
            can_grow_ = growth.AddOneSphere(rad, true, 0, factor);
            axons = growth.axons;
        }
        if (can_grow_)
        {
            axon = growth.axon_to_grow;
            return true;
        }
        else
        {
            // cout << "cannot shrink after dichotomy" << endl;
            return false;
        }
    }
}


double volumeFrustumCone(double r1, double r2, double h)
{
    return M_PI * h * (r1 * r1 + r2 * r2 + r1 * r2) / 3;
}


void AxonGammaDistribution ::create_SWC_file(std::ostream &out)
{
    std::vector<Axon> final_axons;

    final_axons = axons;


    out << "id_cell id_sph id_branch Type X Y Z Rin Rout P" << endl;
    std::sort(final_axons.begin(), final_axons.end(), [](const Axon a, Axon b) -> bool
              { return a.radius < b.radius; }); // sort by size

    double cos_2_all = 0.0;
    for (long unsigned int i = 0; i < final_axons.size(); i++)
    {
        if (final_axons[i].outer_spheres.size() == 1){
            cout << "Axon " << final_axons[i].id << " has only one sphere" << endl;
        } 

        Eigen::Vector3d v = final_axons[i].outer_spheres[final_axons[i].outer_spheres.size()-1].center - final_axons[i].outer_spheres[0].center;
        Eigen::Vector3d z = Eigen::Vector3d(0, 0, 1);
        double cos_angle_2 =  v.dot(z)/(v.norm());
        cos_angle_2 = cos_angle_2*cos_angle_2;
        cos_2_all += cos_angle_2;

        for (long unsigned int j = 0; j < final_axons[i].outer_spheres.size(); j++)
        {
                if (final_axons[i].inner_spheres.size() > 0)
            {
                out << i << " " << j << " -1 axon " << final_axons[i].outer_spheres[j].center[0] << " " << final_axons[i].outer_spheres[j].center[1] << " " << final_axons[i].outer_spheres[j].center[2] << " " << final_axons[i].inner_spheres[j].radius << " " << final_axons[i].outer_spheres[j].radius << " " << final_axons[i].outer_spheres[j].parent_id << endl;
            }
            else
            {
                out << i << " " << j << " -1 axon " << final_axons[i].outer_spheres[j].center[0] << " " << final_axons[i].outer_spheres[j].center[1] << " " << final_axons[i].outer_spheres[j].center[2] << " " << final_axons[i].outer_spheres[j].radius << " " << final_axons[i].outer_spheres[j].radius << " " << final_axons[i].outer_spheres[j].parent_id << endl;
            }
            
        }


    }
    cosPhiSquared = cos_2_all/final_axons.size();


    for (auto &glial_cell : glial_pop1)
    {
        out << glial_cell.id << " " << glial_cell.soma.id << " -1 CellSoma " << glial_cell.soma.center[0] << " " << glial_cell.soma.center[1] << " " << glial_cell.soma.center[2] << " " << glial_cell.soma.radius << " " << glial_cell.soma.radius << " " << -1 << endl;
        for (long unsigned int j = 0; j < glial_cell.ramification_spheres.size(); j++)
        {
            for (long unsigned int k = 0; k < glial_cell.ramification_spheres[j].size(); k++)
            {
                out << glial_cell.id << " " << glial_cell.ramification_spheres[j][k].id << " " << glial_cell.ramification_spheres[j][k].branch_id << " Process " << glial_cell.ramification_spheres[j][k].center[0] << " " << glial_cell.ramification_spheres[j][k].center[1] << " " << glial_cell.ramification_spheres[j][k].center[2] << " " << glial_cell.ramification_spheres[j][k].radius << " " << glial_cell.ramification_spheres[j][k].radius << " " << glial_cell.ramification_spheres[j][k].parent_id << endl;
            }
        }
    }

    for (auto &glial_cell : glial_pop2)
    {
        out << glial_cell.id << " " << glial_cell.soma.id << " -1 CellSoma " << glial_cell.soma.center[0] << " " << glial_cell.soma.center[1] << " " << glial_cell.soma.center[2] << " " << glial_cell.soma.radius << " " << glial_cell.soma.radius << " " << -1 << endl;
        for (long unsigned int j = 0; j < glial_cell.ramification_spheres.size(); j++)
        {
            for (long unsigned int k = 0; k < glial_cell.ramification_spheres[j].size(); k++)
            {
                out << glial_cell.id << " " << glial_cell.ramification_spheres[j][k].id << " " << glial_cell.ramification_spheres[j][k].branch_id << " Process " << glial_cell.ramification_spheres[j][k].center[0] << " " << glial_cell.ramification_spheres[j][k].center[1] << " " << glial_cell.ramification_spheres[j][k].center[2] << " " << glial_cell.ramification_spheres[j][k].radius << " " << glial_cell.ramification_spheres[j][k].radius << " " << glial_cell.ramification_spheres[j][k].parent_id << endl;
            }
        }
    }
}

void AxonGammaDistribution::simulation_file(std::ostream &out, const std::chrono::seconds &duration)
{
    out << "Duration " << duration.count() << endl;
    out << "Num_axons " << axons.size() << endl;
    out << "Voxel " << max_limits[0] << endl;
    out << "Axon icvf " << axons_w_myelin_icvf + axons_wo_myelin_icvf << endl;
    out << "Axon without myelin icvf " <<  axons_wo_myelin_icvf << endl;
    out << "Axon with myelin icvf " << axons_w_myelin_icvf << endl;
    out << "Myelin icvf " << myelin_icvf << endl;
    out << "Glial cell population 1 icvf soma " << glial_pop1_soma_icvf << endl;
    out << "Glial cell population 1 icvf branches " << glial_pop1_branches_icvf << endl;
    out << "Glial cell population 2 icvf soma " << glial_pop2_soma_icvf << endl;
    out << "Glial cell population 2 icvf branches " << glial_pop2_branches_icvf << endl;
    out << "Tortuosity (std of gaussians) " << std_dev << endl;
    out << "Beading amplitude " << beading_variation << endl;
    out << "Overlapping factor " << factor << endl;
    out << "Total icvf " << axons_w_myelin_icvf + axons_wo_myelin_icvf + glial_pop1_soma_icvf + glial_pop1_branches_icvf + glial_pop2_soma_icvf + glial_pop2_branches_icvf << endl;
    out << "C2 " << cosPhiSquared << endl;
    out << "Number of threads " << nbr_threads << endl;
}

std::vector<Eigen::Vector3d> equallySpacedPoints(const Eigen::Vector3d &point1, const Eigen::Vector3d &point2, int n)
{
    std::vector<Eigen::Vector3d> result;

    // Calculate the step size for each dimension
    double stepX = (point2[0] - point1[0]) / (n + 1);
    double stepY = (point2[1] - point1[1]) / (n + 1);
    double stepZ = (point2[2] - point1[2]) / (n + 1);

    // Generate the equally spaced points
    for (int i = 1; i <= n; ++i)
    {
        Eigen::Vector3d newPoint;
        newPoint[0] = point1[0] + i * stepX;
        newPoint[1] = point1[1] + i * stepY;
        newPoint[2] = point1[2] + i * stepZ;
        result.push_back(newPoint);
    }

    return result;
}

std::vector<double> equallySpacedValues(double start, double end, int n)
{
    std::vector<double> result;

    // Calculate the step size
    double step = (end - start) / (n + 1);

    // Generate the equally spaced values
    for (int i = 1; i <= n; ++i)
    {
        double newValue = start + i * step;
        result.push_back(newValue);
    }

    return result;
}


double AxonGammaDistribution::originalFunction(const double &x, const double &outerRadius) {
    return outerRadius - (myelin_thickness(x) + x);
}

double AxonGammaDistribution::derivative(const double &x) {
    if (x <= 0.0001) return -1.0; // Avoid division by zero
    return -(c2 * 2.0 + c3 / (2.0 * x) + 1);
}

double AxonGammaDistribution::findInnerRadius(const double &outerRadius) {
    if (outerRadius <= 0.0) return 0.0;  // Handle invalid inputs

    double guess = outerRadius * 0.7;  // Better initial guess
    double tolerance = 1e-3;
    double step_limit = 0.1 * outerRadius;  // Prevent huge jumps

    int max_iterations = 100;  // Avoid infinite loops
    int iterations = 0;

    while (fabs(originalFunction(guess, outerRadius)) > tolerance) {
        double step = originalFunction(guess, outerRadius) / derivative(guess);

        // Clamp step size to prevent excessive jumps
        if (fabs(step) > step_limit) {
            step = (step > 0 ? step_limit : -step_limit);
        }

        double new_guess = guess - step;

        // Ensure it does not go negative
        if (new_guess <= 0.0) {
            new_guess = 0.01; // Prevent invalid radius
        }

        // Stop if change is very small
        if (fabs(new_guess - guess) < tolerance) {
            break;
        }

        guess = new_guess;

        if (++iterations >= max_iterations) {
            break; // Prevent infinite loops
        }
    }

    if (guess/outerRadius < 0.2) {
        return 0.2*outerRadius;
    }
    else{
        return guess;
    }

}

void AxonGammaDistribution::add_Myelin()
{

    if (target_axons_w_myelin_icvf == 0.0){
        for (long unsigned int k = 0; k < axons.size(); k++){
            axons[k].inner_spheres = axons[k].outer_spheres;
        }
        return;
    }

    double innerRadius;
    Sphere inner_sphere;
    int index;
    // create vector from 0 to axons.size()
    std::vector<int> indices(axons.size());
    std::iota(indices.begin(), indices.end(), 0);
    // shuffle
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(indices.begin(), indices.end(), g);


    for (long unsigned int k = 0; k < indices.size(); k++)
    {
        index = indices[k];
        // ranvier node probability
        int nbr_spheres_to_delete = 1 / axons[index].radius; // a ranvier node is approx 1 um
        nbr_spheres_to_delete = nbr_spheres_to_delete*factor;
        double prob_ranvier_node = 600*factor/axons[index].radius;

        int ranvier_count = 0;
        bool ranvier = false;

        std::uniform_int_distribution<> dis(0, prob_ranvier_node);

        for (long unsigned int i = 0; i < axons[index].outer_spheres.size(); ++i)
        {


            if (axons[index].myelin_sheath ){
                innerRadius = findInnerRadius(axons[index].outer_spheres[i].radius);
            } 
            else {
                innerRadius = axons[index].outer_spheres[i].radius;
            }

            inner_sphere = Sphere(axons[index].outer_spheres[i].id, 2, axons[index].outer_spheres[i].object_id, axons[index].outer_spheres[i].center, innerRadius);
            axons[index].inner_spheres.push_back(inner_sphere);
            if (axons[index].myelin_sheath ){
                // get a random number between 0 and nbr_spheres_for_ranvier
                int ranvier_node = dis(gen);
                if (ranvier_node == 1){
                    ranvier = true;
                }
                // delete some spheres to create ranvier nodes
                if (ranvier && ranvier_count < nbr_spheres_to_delete){
                    axons[index].outer_spheres[i].radius = innerRadius;
                    ranvier_count += 1;
                    //cout << "Ranvier node at sphere " << i << " of axon " << axons[index].id << " position : "<< axons[index].outer_spheres[i].center  << endl;
                    
                }
                else if (ranvier && ranvier_count >= nbr_spheres_to_delete){
                    ranvier = false;
                    ranvier_count = 0;
                }
            } 
        }

    }
    for (long unsigned int k = 0; k < axons.size(); k++){
        if (axons[k].inner_spheres.size() == 0){
            axons[k].inner_spheres = axons[k].outer_spheres;
        }
    }
    
}
