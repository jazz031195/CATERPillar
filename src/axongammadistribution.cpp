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
    
    num_obstacles = 0;
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
    };

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
Eigen::Vector3d findDistantPoint(const Eigen::Vector3d &vector, const Eigen::Vector3d &point)
{
    double distance = 1000;
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
        if (canSpherebePlaced(next, axons, glial_pop1, glial_pop2, false)) {
            intermediate_spheres.emplace_back(next);
        }
    }

    // Add the initial branching sphere to the list
    intermediate_spheres.push_back(first_sphere);

    return intermediate_spheres;
}

bool AxonGammaDistribution::growSecondaryBranch(Glial &glial_cell, int &nbr_spheres, const double &mean_process_length, const double &std_process_length) {

    if (glial_cell.ramification_spheres.empty()) {
        cout << "No branches in glial cell" << endl;
        return false;
    }
    if (nbr_spheres <= 0) {
        cout << "No spheres in glial cell" << endl;
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

    int random_sphere_ind = rand() %  glial_cell.ramification_spheres[random_branch].size() ;
    Sphere random_sphere = glial_cell.ramification_spheres[random_branch][random_sphere_ind];
    double old_length = (random_sphere.center - glial_cell.soma.center).norm();
    std::random_device rd;
    std::mt19937 generator(rd());
    std::normal_distribution<double> length_dist(mean_process_length, std_process_length);
    double length = 0;
    while (length <= old_length) {
        length = length_dist(generator);
    }
    double length_to_grow = length - old_length;

    Eigen::Vector3d vector_to_prev_sphere = (random_sphere.center - glial_cell.ramification_spheres[random_branch][random_sphere_ind - 1].center).normalized();
    Sphere first_sphere;
    Eigen::Vector3d attractor = Eigen::Vector3d(0, 0, 0);
    double radius = random_sphere.radius;
    bool first_sphere_created = GenerateFirstSphereinProcess(first_sphere, attractor, radius, random_sphere, vector_to_prev_sphere, nbr_spheres, nbr_spheres_between, glial_cell.id, nbr_branches);

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
    glial_cell_.ramification_spheres.push_back(vector_first_spheres);
    glial_cell_.principle_processes_lengths.push_back(length);
    glial_cell_.attractors.push_back(attractor);

    // grow glial cell process

    GlialCellGrowth growth(glial_cell_,  glial_pop1, glial_pop2, axons, min_limits, max_limits, min_limits, max_limits, std_dev, glial_cell_.minimum_radius);

    double alpha = -std::log(glial_cell_.minimum_radius / first_sphere.radius);
    Eigen::Vector3d prev_pos = first_sphere.center;
    double distance = radius;
    bool can_grow = true;
    int stop_criteria = 0;
    int nbr_non_checked_spheres = factor * 3;

    int parent = random_sphere.id;

    while (can_grow && stop_criteria < 1000) {
        double R_ = std::exp(-alpha * distance / length_to_grow) * first_sphere.radius;

        if (R_ < glial_cell_.minimum_radius) {
            can_grow = false;

        } else {
            int grow_straight = 0;
            can_grow = growth.AddOneSphere(R_, true, grow_straight, glial_cell_.ramification_spheres.size() - 1, glial_cell_.ramification_spheres.back().size() > nbr_non_checked_spheres, parent, factor);
            if (can_grow) {
                glial_cell_ = std::move(growth.glial_cell_to_grow);
                distance += (prev_pos - glial_cell_.ramification_spheres.back().back().center).norm();
                prev_pos = glial_cell_.ramification_spheres.back().back().center;
                parent = glial_cell_.ramification_spheres.back().back().id;
                
            }
        }

        stop_criteria++;
    }

    if (glial_cell_.ramification_spheres.back().size() > nbr_non_checked_spheres + factor) {
        finished = true;
        nbr_spheres += glial_cell_.ramification_spheres.back().size();
        glial_cell = std::move(glial_cell_);


        return true;
        
    }
         
    return false;
}

bool AxonGammaDistribution::GenerateFirstSphereinProcess(Sphere &first_sphere, Eigen::Vector3d &attractor, const double &radius, const Sphere &sphere_to_emerge_from, const Eigen::Vector3d &vector_to_prev_center, const int &nbr_spheres, const int &nbr_spheres_between, const int &cell_id, const int &branch_id) {
    
    bool stop = false;
    int tries_ = 0;
    int max_nbr_tries = 100;
    attractor = Eigen::Vector3d(0, 0, 0);


    while (!stop && tries_ < max_nbr_tries) {
        Eigen::Vector3d vector = {0, 0, 0};
        Eigen::Vector3d point = {0, 0, 0};
        sphere_to_emerge_from.getPointOnSphereSurface(point, vector, vector_to_prev_center);
        attractor = findDistantPoint(vector, point);
        first_sphere = Sphere(nbr_spheres + nbr_spheres_between + 1, cell_id, 1, point, radius, branch_id, sphere_to_emerge_from.id);
        if (canSpherebePlaced(first_sphere, axons, glial_pop1, glial_pop2, false)) {
            stop = true;
        } else {
            tries_++;
            if (tries_ == max_nbr_tries) return false;
        }
    }
    return true;
}

bool AxonGammaDistribution::growPrimaryBranch(Glial &glial_cell, const int &nbr_spheres, const double &mean_process_length, const double &std_process_length) {


    int j = glial_cell.ramification_spheres.size();
    
    // Ensure vectors have space for branch `j`
    if (glial_cell.principle_processes_lengths.size() <= j) {
        glial_cell.principle_processes_lengths.resize(j + 1, 0);
    }
    if (glial_cell.ramification_spheres.size() <= j) {
        glial_cell.ramification_spheres.resize(j + 1);
    }
    if (glial_cell.attractors.size() <= j) {
        glial_cell.attractors.resize(j + 1);
    }

    // Generate length from Gaussian distribution
    std::random_device rd;
    std::mt19937 generator(rd());
    std::normal_distribution<double> length_dist(mean_process_length, std_process_length);
    double length = 0;
    while (length < 3) {
        length = length_dist(generator);
    }
    glial_cell.principle_processes_lengths[j] = length;

    // Find points on the surface of the glial soma sphere
    bool stop = false;
    int tries = 0;
    // Add the first sphere with intermediate spheres
    int parent = 0;  // parent is the soma
    int nbr_spheres_between = factor - 1;

    Eigen::Vector3d vector_to_prev_center = {0, 0, 0};
    Sphere first_sphere;
    Eigen::Vector3d attractor = Eigen::Vector3d(0, 0, 0);
    double initial_radius = glial_cell.soma.radius/3;
    bool first_sphere_created = GenerateFirstSphereinProcess(first_sphere, attractor, initial_radius, glial_cell.soma, vector_to_prev_center, nbr_spheres, nbr_spheres_between, glial_cell.id, j);

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

    double alpha = -std::log(glial_cell.minimum_radius / initial_radius);

    Eigen::Vector3d prev_pos = first_sphere.center;

    Glial glial_cell_ = glial_cell;
    glial_cell_.ramification_spheres[j] = vector_first_spheres;
    glial_cell_.attractors[j] = attractor;


    GlialCellGrowth growth(glial_cell_, glial_pop1, glial_pop2, axons, min_limits, max_limits, min_limits, max_limits, std_dev, glial_cell_.minimum_radius);
    
    double distance = initial_radius;
    while (can_grow) {
 
        // Calculate the radius at the current distance
        double R_ = std::exp(-alpha * distance / length) * initial_radius;
        if (R_ < glial_cell_.minimum_radius) {
            can_grow = false;
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

        }
    }

    // Check if the branch has grown sufficiently
    if (glial_cell_.ramification_spheres[j].size() < nbr_non_checked_spheres + factor) {
        glial_cell_.ramification_spheres.pop_back();
        glial_cell_.principle_processes_lengths.pop_back();
        glial_cell_.attractors.pop_back();
        //cout << "Branch did not grow enough" << endl;
        return false;
    }
    else {
        //cout << "Branch grown" << endl;
        glial_cell = glial_cell_;
        return true;
    }

}


void AxonGammaDistribution::growFirstPrimaryBranches(Glial &glial_cell, const int &number_ramification_points, int &nbr_spheres, const double &mean_process_length, const double &std_process_length)
{

    glial_cell.ramification_spheres = std::vector<std::vector<Sphere>>(number_ramification_points, std::vector<Sphere>());
    glial_cell.principle_processes_lengths = std::vector<double>(number_ramification_points);


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
            nbr_spheres += glial_cell.ramification_spheres.back().size();
        }
        
    }


}

bool AxonGammaDistribution::withinBounds(const Eigen::Vector3d &pos, const double &distance)
{
    // Check if the point is inside the dilated box
    for (int i = 0; i < 3; ++i) {
        double min_bound = min_limits[i];
        double max_bound = max_limits[i];
        if (pos[i] < min_bound || pos[i] > max_bound) {
            return false; // Point is outside the dilated box
        }
    }
    
    return true; // Point is inside the dilated box
}

bool AxonGammaDistribution::get_begin_end_point(Eigen::Vector3d &Q, Eigen::Vector3d &D, double &angle)
{

    std::random_device rd;
    std::mt19937 gen(rd());
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

        while (icvf_ < icvf_to_reach)
        {
            if (tried > 1000)
            {
                std::string message = "Radii distribution cannot be sampled [Min. radius Error]\n";
                std::cout << message << std::endl;
                assert(0);
            }
            double jkr = distribution(generator);

            // generates the radii in a list
            if (jkr > min_radius && jkr < 2)
            {
                if (target_axons_w_myelin_icvf > 0){   

                    if (axons_w_myelin_icvf < target_axons_w_myelin_icvf){
                        has_myelin.push_back(true);
                        jkr += myelin_thickness(jkr);
                        axons_w_myelin_icvf += (jkr * jkr * M_PI * (max_limits[2] - min_limits[2]))/total_volume;
                    } 
                    else if (axons_wo_myelin_icvf < target_axons_wo_myelin_icvf){
                        axons_wo_myelin_icvf += (jkr * jkr * M_PI * (max_limits[2] - min_limits[2]))/total_volume;
                        has_myelin.push_back(false);
                    }
                    else{
                        break;
                    } 
                }
                else{
                    if (axons_wo_myelin_icvf < target_axons_wo_myelin_icvf){
                        axons_wo_myelin_icvf += (jkr * jkr * M_PI * (max_limits[2] - min_limits[2]))/total_volume;
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

        num_obstacles = radii_.size();

        std::cout << "Number of axons :" << num_obstacles << endl;
        std::cout <<"initial icvf axons without myelin :" << axons_wo_myelin_icvf << endl;
        std::cout <<"initial icvf axons with myelin :" << axons_w_myelin_icvf << endl;
        axons_wo_myelin_icvf = 0.0;
        axons_w_myelin_icvf = 0.0;

    }
    else
    {
        num_obstacles = 0;
    }
}

bool AxonGammaDistribution::PlaceAxon(const int &axon_id, const double &radius_for_axon, const Eigen::Vector3d &Q, const Eigen::Vector3d &D, std::vector<Axon> &new_axons, const bool &has_myelin, const double &angle_, const bool &outside_voxel, const bool &push_axons)
{

    Axon ax = Axon(axon_id, Q, D, radius_for_axon, beading_variation, has_myelin, angle_, outside_voxel); // axons for regrow batch
    
    if (has_myelin) {
        double inner_radius = findInnerRadius(radius_for_axon);
        ax.inner_radius = inner_radius;
    }

    Sphere sphere = Sphere(0, ax.id, 0, Q, radius_for_axon);

    bool no_overlap_axons = canSpherebePlaced(sphere, axons, glial_pop1, glial_pop2, push_axons);
    bool no_overlap_new_axons = canSpherebePlaced(sphere, new_axons, glial_pop1, glial_pop2, false);

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

bool AxonGammaDistribution::pushAxonSpheres(std::vector<Axon> &axs, const std::vector<Glial> &astros, const std::vector<Glial> &oligos, Axon &axon, const Sphere &sph) {

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

            if (canSpherebePlaced(new_sphere, axs, astros, oligos, false) && check_borders(min_limits, max_limits, new_sphere.center, new_sphere.radius)) {
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

bool AxonGammaDistribution::canSpherebePlaced(const Sphere &sph, std::vector<Axon> &axs, const std::vector<Glial> &astros, const std::vector<Glial> &oligos, const bool &with_push) 
{

    for (int i = 0 ; i < axs.size(); i++) {

        if (!(axs[i].id == sph.object_id && sph.object_type == 0))
        {
            // Check overlap
            if (axs[i].isSphereInsideAxon(sph)) 
            {
                
                if (with_push) {
                    Axon axon_to_push = axs[i];
                    
                    bool can_axon_be_pushed = pushAxonSpheres(axs, astros, oligos, axon_to_push, sph);
                    if (!can_axon_be_pushed){
                        return false;
                    }
                    else{
                        //cout <<" pushed an axon before growth! " << endl;
                        axs[i] = std::move(axon_to_push);
                    }
                }
                else{
                    return false;
                }
                
                return false;
            }
        }
    }

    std::vector<Glial> glial_cells = astros;    
    glial_cells.insert(glial_cells.end(), oligos.begin(), oligos.end());
    // check collision other glial cells 
    for (auto &glial : glial_cells)
    {
        if (sph.object_type == 0 || (sph.object_type == 1 && sph.object_id != glial.id)){
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

void AxonGammaDistribution::process_point(const Eigen::Vector3d &point, const std::vector<Axon> &axs, const std::vector<Glial> &glial_pop1, const std::vector<Glial> &oligos, int &axons_count, int &myelin_count, int &glial_pop1_somas_count, int &glial_pop1_process_count, int &oligos_somas_count, int &oligos_process_count, int &extracellular_count) {
    
    Sphere sph(-1, -1, -1, point, 0.001);
    bool collision_detected = false;

    // Check collision with axons
    for (const auto &axon : axs) {
        if (axon.myelin_sheath){

            if (axon.isSphereInsideInnerAxon(sph)){
                axons_count += 1;
                collision_detected = true;
            } 
            else if (axon.isSphereInsideAxon(sph)) {
                myelin_count += 1;
                collision_detected = true;
            }
        }
        else{
            if (axon.isSphereInsideAxon(sph)) {
                axons_count += 1;
                collision_detected = true;
            }
        }  
        if (collision_detected) return;
    }

    // Check collision with glial_pop1
    for (const auto &glial : glial_pop1) {
        if (glial.isNearGlialCell(sph.center, sph.radius )) {
            if (glial.collides_with_GlialCell(sph)) {
                glial_pop1_somas_count += 1;
                collision_detected = true;
                break;
            }
            // Check with branches of glial_pop1
            for (const auto &branch : glial.ramification_spheres) {
                for (const auto &sph_ : branch) {
                    if (sph_.CollideswithSphere(sph,0)) {
                        glial_pop1_process_count += 1;
                        collision_detected = true;
                        break;
                    }
                }
                if (collision_detected) break;
            }
        }
        if (collision_detected) break;
    }

    if (collision_detected) return;

    // Check collision with glial_pop2
    for (const auto &glial : oligos) {
        if (glial.isNearGlialCell(sph.center,  sph.radius )) {
            if (glial.collides_with_GlialCell(sph)) {
                oligos_somas_count += 1;
                collision_detected = true;
                break;
            }
            // Check with branches of glial_pop2
            for (const auto &branch : glial.ramification_spheres) {
                for (const auto &sph_ : branch) {
                    if (sph_.CollideswithSphere(sph, 0)) {
                        oligos_process_count += 1;
                        collision_detected = true;
                        break;
                    }
                }
                if (collision_detected) break;
            }
        }
        if (collision_detected) break;
    }

    if (!collision_detected) {
        extracellular_count += 1;  // If no collisions, increment extracellular count
    }
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


// Function to find the kappa value corresponding to a given c2
double AxonGammaDistribution::c2toKappa(double c2, double tol = 1e-3, double kappa_min = 0, double kappa_max = 64) {
    if (c2 == 1.0) {
        return std::numeric_limits<double>::infinity(); // Return infinity for perfect alignment
    } else if (c2 < 1.0 / 3.0) {
        return 0; // Return 0 for isotropic distribution
    }

    // Generate a range of kappa values
    std::vector<double> kappas;
    for (double k = kappa_min; k <= kappa_max; k += tol) {
        kappas.push_back(k);
    }

    // Compute Fs and c2s for each kappa
    std::vector<double> Fs(kappas.size());
    std::vector<double> c2s(kappas.size());

    for (size_t i = 0; i < kappas.size(); ++i) {
        double sqrt_kappa = std::sqrt(kappas[i]);
        double exp_neg_kappa = std::exp(-kappas[i]);
        Fs[i] = (std::sqrt(M_PI) / 2.0) * exp_neg_kappa * erfi(sqrt_kappa);

        if (kappas[i] > 0) {
            c2s[i] = 1.0 / (2.0 * sqrt_kappa * Fs[i]) - 1.0 / (2.0 * kappas[i]);
        } else {
            c2s[i] = 1.0 / 3.0; // Edge case for kappa = 0
        }
    }

    // Set the last c2 value explicitly to 1
    c2s.back() = 1.0;

    // Find the kappa value closest to the target c2
    double min_diff = std::numeric_limits<double>::max();
    size_t idx_c2 = 0;
    for (size_t i = 0; i < c2s.size(); ++i) {
        double diff = std::abs(c2 - c2s[i]);
        if (diff < min_diff) {
            min_diff = diff;
            idx_c2 = i;
        }
    }

    double kappa = kappas[idx_c2];
    return std::max(0.0, kappa); // Ensure kappa is non-negative
}

// Function to draw angles from the Watson distribution given c2
double AxonGammaDistribution::draw_angle(double kappa) {

    // Find kappa bounds
    double max_kappa = *std::max_element(cdf.kappas.begin(), cdf.kappas.end());
    double min_kappa = *std::min_element(cdf.kappas.begin(), cdf.kappas.end());

    if (kappa > max_kappa) {
        return 0.0;
    }
    if (kappa < min_kappa) {
        
        cout << "Kappa value is out of bounds "<< kappa << ", new kappa value : " << min_kappa << endl;
        kappa = min_kappa;

    }

    // Interpolate for the kappa row
    std::vector<double> interpolated_row(cdf.angles.size());
    for (size_t i = 0; i < cdf.angles.size(); ++i) {
        std::vector<double> cdf_column;
        for (const auto& row : cdf.data) {
            cdf_column.push_back(row[i]);
        }
        interpolated_row[i] = interpolate(kappa, cdf.kappas, cdf_column);
    }
    
    // Generate random values and map them to angles
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    double random_value = dis(gen);
    double sampled_angle = interpolate(random_value, interpolated_row, cdf.angles);
    sampled_angle = sampled_angle*M_PI/180;

    return sampled_angle;
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



void AxonGammaDistribution::createBatches(std::vector<double> &radii_, std::vector<int> &indices, std::vector<Axon> &new_axons, std::vector<bool> &has_myelin, std::vector<double> &angles)
{
    long int tries_threshold = (max_limits[0]- min_limits[0]) * (max_limits[1]-min_limits[1]) * 500;
    new_axons.clear();

    double angle_;
    bool outside_voxel = false;

    std::vector<int> indices_to_erase = {};

    for (long unsigned int i = 0; i < indices.size(); ++i) // create axon for each radius
    {
        // std::cout << "radii created :" << i << "/" << batch_radii.size() << endl;
        bool next = false;
        int tries = 0;

        while (!next && tries < tries_threshold)
        {
            
            Vector3d Q, D;
            
            int index_of_axon = indices[i];
            angle_ = angles[indices[i]];
            outside_voxel = !get_begin_end_point(Q, D, angle_);
            bool has_myelin_ = has_myelin[indices[i]];
            if (tries> tries_threshold/2) {
                next = PlaceAxon(index_of_axon, radii_[indices[i]], Q, D, new_axons, has_myelin_, angle_, outside_voxel, true);

            }
            else{
                next = PlaceAxon(index_of_axon, radii_[indices[i]], Q, D, new_axons, has_myelin_, angle_, outside_voxel, false);
            }
            //cout << "Q : " << Q << endl;
            if (!next)
            {
                tries += 1;
                
            }
        }
        if (tries >= tries_threshold)
        {
            cout << "Tries exceeded threshold" << endl;
            indices_to_erase.push_back(i);
        }
    }

    for (long unsigned int i = 0; i < indices_to_erase.size(); i++)
    {
        indices.erase(indices.begin() + indices_to_erase[i]);
        radii.erase(radii.begin() + indices_to_erase[i]);
        has_myelin.erase(has_myelin.begin() + indices_to_erase[i]);
        angles.erase(angles.begin() + indices_to_erase[i]);
    }

    // check if number of new axons is same as number of indices
    if (indices.size() != new_axons.size())
    {
        cout << "Number of axons created is not equal to number of indices" << endl;
        assert(indices.size() == new_axons.size());
    }
    /*
    for (long unsigned int i = 0; i < new_axons.size(); i++)
    {
        //cout <<" (" << new_axons[i].radius << "," << new_axons[i].begin[0] << "," << new_axons[i].begin[1] << ")," << endl;
        cout <<"axon.id :  "<< new_axons[i].id << endl;
    }
    assert(0);
    */
    

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
    std::vector<double> angles(num_obstacles, 0.0);
    
    if (cosPhiSquared != 1.0){
        angles = generate_angles(num_obstacles);
    } 

    cout <<"number of threads :" << nbr_threads << endl;

    calculate_c2(angles);

    if (num_obstacles != 0){
        // indices for axons
        std::vector<int> indices;
        for (unsigned i = 0; i < num_obstacles; i++)
        {
            indices.push_back(i);
        }
        // std::cout << "Parallel growth simulation" << endl;

        bool stop = false;

        while (!stop)
        {

            num_obstacles = radii.size();

            if (num_obstacles == 0)
            {
                stop = true;
                break;
            }
            std::vector<int> subsets; // depending on which batch
            // std::cout << "radii size =" << num_obstacles << endl;
            //  divides obstacles into subsets
            setBatches(num_obstacles, subsets);

            createBatches(radii, indices, axons, has_myelin, angles);

            // cout << " length of subsets :" << subsets.size() << endl;
            growBatches(radii, indices, subsets, has_myelin, angles); // grows all axons, fills list of stuck radii, deletes empty axons
            std::vector<double> radii_to_regrow = stuck_radii;
            std::vector<int> indices_to_regrow = stuck_indices;

            ICVF(axons, glial_pop1, glial_pop2);
            cout << "ICVF axons :" << axons_icvf << endl;

            
            int nbr_attempts = 0;
            double percentage_swelling = 0.5;
            double minimum_percentage_swelling = 1e-6;

            while (axons_icvf < target_axons_icvf && nbr_attempts < 1000)
            {
                double old_icvf = axons_icvf;
                SwellAxons(percentage_swelling);
                cout << "new ICVF " << axons_icvf  << " old icvf : "<< old_icvf <<" percentage swelling :"<< percentage_swelling<< endl;
                if ((axons_icvf - old_icvf) < 1e-6 && percentage_swelling >= minimum_percentage_swelling)
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
    if (canSpherebePlaced(sph_, axons, glial_pop1, glial_pop2, true))
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
    for (size_t i = 0; i < ax.outer_spheres.size(); i++) {
        
        Sphere sph = ax.outer_spheres[i];
        bool can_swell = SwellSphere(sph, percentage);
        if (can_swell){
            ax.outer_spheres[i] = sph;
        }
    }
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
        SwellAxon(axon, percentage);
        axon.update_Volume(factor, min_limits, max_limits);
        axon.updateBox();
        // check ICVF 
        ICVF(axons, glial_pop1, glial_pop2);
        
        if (axons_icvf > target_axons_icvf)
        {
            break;
        }
    }

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
    
    if (glial_pop1.size() > 0.0 || glial_pop2.size() > 0.0)
    {
        if (target_glial_pop1_branches_icvf > 0.0 || target_glial_pop2_branches_icvf > 0.0)
        {
            GrowAllGlialCells();
        }
    }
    
    std::vector<double> stuck_radii_; 
    std::vector<int> stuck_indices_;
    bool cells_ok = FinalCheck(axons, stuck_radii_, stuck_indices_);

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

void AxonGammaDistribution::growBranches(std::vector<Glial>& glial_cell_list, const int &population_nbr) {

    double target_branches_icvf = 0.0;
    double current_branches_icvf = 0.0;

    int nbr_primary_processes = 1;

    double mean_process_length = 0.0;
    double std_process_length = 0.0;

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
        for (int i = 0; i < glial_cell_list.size(); i++) 
        {
            growFirstPrimaryBranches(glial_cell_list[i], nbr_primary_processes, nbr_spheres[i], mean_process_length, std_process_length);
            glial_cell_list[i].compute_processes_icvf(factor);
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
    const int max_tries = 1000000;

    while (current_branches_icvf < target_branches_icvf) {

       for (size_t i = 0; i < glial_cell_list.size(); ++i) {
            if (glial_cell_list[i].allow_branching){
                growSecondaryBranch(glial_cell_list[i], nbr_spheres[i], mean_process_length, std_process_length);
            }
            else{
                growPrimaryBranch(glial_cell_list[i], nbr_spheres[i], mean_process_length, std_process_length);
            }
        }
        if (nbr_tries%20 == 0 && nbr_tries >0) {
            // Recompute processes and update ICVF for all glial_pop1 in parallel
            for (size_t i = 0; i < glial_cell_list.size(); ++i) {
                glial_cell_list[i].compute_processes_icvf(factor);
            }
            // Update global ICVF and display progress
            if (population_nbr == 1) {
                ICVF(axons, glial_cell_list, glial_pop2);
                current_branches_icvf = glial_pop1_branches_icvf;
            }
            else if (population_nbr == 2) {
                ICVF(axons, glial_pop1, glial_cell_list);
                current_branches_icvf = glial_pop2_branches_icvf;
            }

            display_progress(current_branches_icvf, target_branches_icvf);
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
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_radius(glial_pop1_radius_mean, glial_pop1_radius_std);
    
    //Astrocytes
    double glial_icvf_ = 0.0;

    while (glial_icvf_ < target_glial_pop1_soma_icvf) {
        Sphere s;
        bool can_be_placed = false;

        for (int attempt = 0; attempt < 10000; ++attempt) {
            // Generate random position

            double rad = dis_radius(gen);

            if (rad < 0.5){
                rad = 0.5;
            }
            else if (rad > 6){
                rad = 6;
            }

            double min_distance_to_border = rad;
            std::uniform_real_distribution<> dis_x(min_limits[0] + min_distance_to_border, max_limits[0] - min_distance_to_border);
            std::uniform_real_distribution<> dis_y(min_limits[1] + min_distance_to_border, max_limits[1] - min_distance_to_border);
            std::uniform_real_distribution<> dis_z(min_limits[2] + min_distance_to_border, max_limits[2] - min_distance_to_border);

            s = Sphere(0, glial_pop1.size(), 1, {
                dis_x(gen), dis_y(gen), dis_z(gen)
            }, rad);

            // Check if sphere can be placed
            if (canSpherebePlaced(s, axons, glial_pop1, glial_pop2, false)) {
                // Check distance to voxel border
                bool within_bounds = true;
                for (int i = 0; i < 3; ++i) {
                    if (s.center[i] <= min_limits[i] + s.radius || s.center[i] >= max_limits[i] - s.radius) {
                        within_bounds = false;
                        break;
                    }
                }
                if (within_bounds) {
                    can_be_placed = true;
                    break;  // Exit loop if placed
                }
            }
        }

        if (can_be_placed) {
            Glial glial_cell = Glial(s.object_id, s, glial_pop1_branching);
            glial_pop1.push_back(glial_cell);
            double glial_volume = 4 * M_PI * glial_cell.soma.radius * glial_cell.soma.radius * glial_cell.soma.radius / 3;
            glial_icvf_ += glial_volume / total_volume;

        } else {
            break;  // Exit while loop if unable to place more glial cells
        }
    }
    glial_pop1_soma_icvf = glial_icvf_;

    //glial_pop2

    glial_icvf_ = 0.0;
    std::uniform_real_distribution<> dis_radius_oligo(glial_pop2_radius_mean, glial_pop2_radius_std);

    while (glial_icvf_ < target_glial_pop2_soma_icvf) {
        Sphere s;
        bool can_be_placed = false;

        for (int attempt = 0; attempt < 10000; ++attempt) {
            // Generate random position

            double rad = dis_radius_oligo(gen);

            double min_distance_to_border = rad;
            std::uniform_real_distribution<> dis_x(min_limits[0] + min_distance_to_border, max_limits[0] - min_distance_to_border);
            std::uniform_real_distribution<> dis_y(min_limits[1] + min_distance_to_border, max_limits[1] - min_distance_to_border);
            std::uniform_real_distribution<> dis_z(min_limits[2] + min_distance_to_border, max_limits[2] - min_distance_to_border);

            s = Sphere(0, glial_pop2.size(), 1, {
                dis_x(gen), dis_y(gen), dis_z(gen)
            }, rad);

            // Check if sphere can be placed
            if (canSpherebePlaced(s, axons, glial_pop1, glial_pop2, false)) {
                // Check distance to voxel border
                bool within_bounds = true;
                for (int i = 0; i < 3; ++i) {
                    if (s.center[i] <= min_limits[i] + s.radius || s.center[i] >= max_limits[i] - s.radius) {
                        within_bounds = false;
                        break;
                    }
                }
                if (within_bounds) {
                    can_be_placed = true;
                    break;  // Exit loop if placed
                }
            }
        }

        if (can_be_placed) {
            Glial glial_cell = Glial(s.object_id, s);
            glial_pop2.push_back(glial_cell);
            double glial_volume = 4 * M_PI * glial_cell.soma.radius * glial_cell.soma.radius * glial_cell.soma.radius / 3;
            glial_icvf_ += glial_volume / total_volume;

        } else {
            break;  // Exit while loop if unable to place more glial cells
        }
    }
    glial_pop2_soma_icvf = glial_icvf_;

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
            if (!canSpherebePlaced(axs[j].outer_spheres[i], axs, glial_pop1, glial_pop2, false))
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

bool AxonGammaDistribution::SanityCheck(std::vector<Axon>& growing_axons, std::vector<double> &stuck_radii_, std::vector<int> &stuck_indices_) {
    
    std::vector<double> results(growing_axons.size(), -1);
    std::vector<Axon> axons_to_check_collision_with = growing_axons;

    for (long unsigned int j = 0; j < growing_axons.size(); ++j) {

        results[j] = -1;
        if (growing_axons[j].outer_spheres.size() > 0) {
            for (unsigned k = 0; k < growing_axons[j].outer_spheres.size(); ++k) {
                Sphere last_sphere = growing_axons[j].outer_spheres[k];

                if (!canSpherebePlaced(last_sphere, axons_to_check_collision_with, glial_pop1, glial_pop2, false)) {
                    results[j] = growing_axons[j].id;

                    auto foundObject = std::find_if(axons_to_check_collision_with.begin(), axons_to_check_collision_with.end(),
                        [j, &growing_axons](const Axon& ax) {
                            return ax.id == growing_axons[j].id;
                        });

                    if (foundObject != axons_to_check_collision_with.end()) {
                        axons_to_check_collision_with.erase(
                            std::remove_if(axons_to_check_collision_with.begin(), axons_to_check_collision_with.end(),
                                [j, &growing_axons](const Axon& ax) {
                                    return ax.id == growing_axons[j].id;
                                }),
                            axons_to_check_collision_with.end()
                        );
                    }

                    //std::cout << "Sanity check: axon " << growing_axons[j].id << " collides with environment" << std::endl;
                }
            }
        }
    }

    bool not_collide = true;
    std::vector<int> positions_to_destroy;
    for (double result : results) {
        if (result >= 0) {
            not_collide = false;

            auto foundObject = std::find_if(growing_axons.begin(), growing_axons.end(), [result](const Axon& ax) {
                return ax.id == result;
            });

            if (foundObject != growing_axons.end()) {
                positions_to_destroy.push_back(std::distance(growing_axons.begin(), foundObject));
            } else {
                assert(0);
            }
        }
    }

    for (int pos : positions_to_destroy) {
        growing_axons[pos].destroy();
        stuck_radii_.push_back(growing_axons[pos].radius);
        stuck_indices_.push_back(growing_axons[pos].id);
    }

    return not_collide;
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
    double& stuck_radii_,
    int& stuck_indices_
) {
   
    // 1) Possibly vary the radius (beading)
    double varied_radius = (beading_variation > 0)
        ? RandomradiusVariation(axon)
        //? radiusVariation(axon)
        : axon.radius;

    // 2) Attempt to add a sphere (returns true if successful)
    bool can_grow = growth.AddOneSphere(varied_radius, /*create_sphere=*/true, grow_straight, factor);
    // 3) Sync axon data back
    //    growth.axon_to_grow updated by AddOneSphere
    growth.axon_to_grow.growth_attempts = axon.growth_attempts; 

    // 4) Check if axon is finished
    if (growth.finished) {
        //cout << "axon " << axon.id << " is finished" << endl;
        finished = 1; 
        stuck_radii_ = -1;
    }
    else { 
        // The axon is still growing
        if (!can_grow && axon_can_shrink) {
            // We tried adding a sphere and failed to grow. Attempt to shrink radius.
            bool shrinkOk = shrinkRadius(growth, varied_radius, axon); 
            if (!shrinkOk) {
                // If shrinking didn't help either
                if (axon.growth_attempts < 10) {
                    //cout <<"axon " << axon.id << " only has one sphere" << endl;
                    // Retry with a new attempt
                    axon.keep_one_sphere();  // Revert the last sphere or do whatever keep_one_sphere does
                    growth.ind_axons_pushed.clear();
                    growth.axon_to_grow = axon;
                    finished = 0;
                    stuck_radii_ = -1;
                }
                else {
                    //cout <<"axon " << axon.id << " is stuck" << endl;
                    // The axon is stuck
                    axon.destroy(); 
                    growth.ind_axons_pushed.clear();
                    finished = 1;
                    stuck_radii_ = axon.radius;
                    stuck_indices_ = axon.id;
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
            if (axon.growth_attempts < 10) {
                //cout <<"axon " << axon.id << " only has one sphere" << endl;
                // Retry with a new attempt
                axon.keep_one_sphere();  // Revert the last sphere or do whatever keep_one_sphere does
                growth.axon_to_grow = axon;
                finished = 0;
              
            }
            else {
                //cout <<"axon " << axon.id << " is stuck" << endl;
                // The axon is stuck
                axon.destroy(); 
                finished = 1;
                stuck_radii_ = axon.radius;
                stuck_indices_ = axon.id;
            }
        }
        else if (can_grow) {
            finished = 0;
        }
    }

    // 5) Update "straightness" logic
    update_straight(can_grow, grow_straight, straight_growths, ondulation_factor);
}


void AxonGammaDistribution::growAxon(Axon& axon_to_grow, double& stuck_radii_, int& stuck_indices_, std::vector<int>& ind_axons_pushed, std::vector<Axon> &axons_pushed) {
    int finished = 0;
    int grow_straight = 0;
    int straight_growths = 0;


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
    int max_nbr_tries = 1e10;
    // Keep trying to grow until finished
    while (finished == 0 && tries_ < max_nbr_tries) {

        if (!axon_to_grow.outer_spheres.empty()) {
            // Attempt a single growth step
            growthThread(axons, axon_to_grow, growth, finished, grow_straight, straight_growths, stuck_radii_, stuck_indices_);

            if (stuck_radii_>0) {
                tries_ += 1;
            }
        }
        else {
            //cout <<"finished = 1" << endl;
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
        stuck_radii_ = axon_to_grow.radius;
        stuck_indices_ = axon_to_grow.id;
        ind_axons_pushed = {};
    }
    else{
        ind_axons_pushed =  growth.ind_axons_pushed;
        for (int i : ind_axons_pushed) {
            axons_pushed.push_back(growth.modified_axons[i]);
        }
    }
}
void AxonGammaDistribution::processBatchWithThreadPool(
    std::vector<Axon>& axons_to_grow,
    std::vector<double>& stuck_radii,
    std::vector<int>& stuck_indices, std::vector<std::vector<int>> &ind_axons_pushed, std::vector<std::vector<Axon>> &axons_pushed) 
{
    ThreadPool pool(nbr_threads);

    std::vector<std::future<void>> futures; // Store futures for synchronization

    for (size_t i = 0; i < axons_to_grow.size(); ++i) {
        futures.emplace_back(pool.enqueueTask(
            [this, &axons_to_grow, &stuck_radii, &stuck_indices, &ind_axons_pushed, &axons_pushed, i]() {
                growAxon(axons_to_grow[i], stuck_radii[i], stuck_indices[i], ind_axons_pushed[i], axons_pushed[i]);
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

void AxonGammaDistribution::ModifyAxonsStartingPoint(const std::vector<int> &stuck_indices){
    for (int i = 0; i < axons.size(); i++) {

        if (std::find(stuck_indices.begin(), stuck_indices.end(), axons[i].id) != stuck_indices.end()) {
            Vector3d Q, D;
            int index_of_axon = axons[i].id;
            double angle_ = axons[i].angle;
            bool outside_voxel = !get_begin_end_point(Q, D, angle_);
            Sphere sphere = Sphere(0, axons[i].id, 0, Q, axons[i].radius);
            bool no_overlap_axons = canSpherebePlaced(sphere, axons, glial_pop1, glial_pop2, false);
            int nbr_tries = 0;
            while(!no_overlap_axons && nbr_tries < 1000){
                outside_voxel = !get_begin_end_point(Q, D, angle_);
                sphere = Sphere(0, axons[i].id, 0, Q, axons[i].radius);
                no_overlap_axons = canSpherebePlaced(sphere, axons, glial_pop1, glial_pop2, false);
            }
            if (no_overlap_axons) {
                axons[i].outer_spheres.clear();
                axons[i].outer_spheres.push_back(sphere);
                axons[i].update_Volume(factor, min_limits, max_limits);
                axons[i].updateBox();
                axons[i].end = D;
                axons[i].outside_voxel = outside_voxel;
                axons[i].angle = angle_;
                break;
            }
        }
    }
}

void AxonGammaDistribution::createBatch(const std::vector<double> &radii_, const std::vector<int> &indices, const int &num_subset, const int &first_index_batch, std::vector<Axon> &new_axons, const std::vector<bool> &has_myelin, std::vector<double> &angles)
{
    long int tries_threshold = (max_limits[0]- min_limits[0]) * (max_limits[1]-min_limits[1]) * 50;
    new_axons.clear();

    std::vector<double>::const_iterator startIterator = radii_.begin() + first_index_batch;             // Start from index
    std::vector<double>::const_iterator stopIterator = radii_.begin() + first_index_batch + num_subset; // Stop when batch is finished
    // radii in batch
    std::vector<double> batch_radii(startIterator, stopIterator);
    double angle_;
    bool outside_voxel = false;

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
                next = PlaceAxon(index_of_axon, batch_radii[i], Q, D, new_axons, has_myelin_, angle_, outside_voxel, false);

            }
            else{
                next = PlaceAxon(index_of_axon, batch_radii[i], Q, D, new_axons, has_myelin_, angle_, outside_voxel, false);
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

    cout <<"growing axons from : " << batch_growing_axons[0].id << "to " << batch_growing_axons[number_axons_to_grow-1].id << endl;
    //vector full of -1 with size of number of axons
    std::vector<double> all_stuck_radii(number_axons_to_grow, -1);  

    std::vector<int> all_stuck_indices(number_axons_to_grow, -1);   

    std::vector<std::vector<int>> ind_axons_pushed;  

    std::vector<std::vector<Axon>> axons_pushed;   

    ind_axons_pushed.resize(number_axons_to_grow);
    axons_pushed.resize(number_axons_to_grow); 

    processBatchWithThreadPool(batch_growing_axons, all_stuck_radii, all_stuck_indices, ind_axons_pushed, axons_pushed);


    for (int i = 0; i < number_axons_to_grow; i++) // for each axon
    {
        if (all_stuck_radii[i] > 0)
        {
            stuck_radii_.push_back(all_stuck_radii[i]);
            stuck_indices_.push_back(all_stuck_indices[i]);
        }
    }
    
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
    
    cout <<"number of axons to check in sanity check : " << batch_growing_axons.size() << endl;

    SanityCheck(batch_growing_axons, stuck_radii_, stuck_indices_);
    
    cout << stuck_radii_.size() << " axons are stuck" << endl;

    // add axons in batch that worked in axons
    for (long unsigned int i = 0; i < batch_growing_axons.size(); i++)
    {
        //cout << "growing_axons["<<growing_axons[i].id <<"].outer_spheres.size(): " << growing_axons[i].outer_spheres.size() << endl;
        if (batch_growing_axons[i].outer_spheres.size() > 0)
        {
            for (long unsigned int j = 0; j < axons.size(); j++) {
                if (axons[j].id == batch_growing_axons[i].id) {
                    // erase axon from axons
                    axons[j] = std::move(batch_growing_axons[i]);
                }
            }
        }
    }

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

void AxonGammaDistribution::growBatches(std::vector<double> &radii_, std::vector<int> &indices, std::vector<int> &subsets_, std::vector<bool> &has_myelin, std::vector<double> &angles)
{

    
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
        if (axons[i].outer_spheres.size() == 1) {
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

double AxonGammaDistribution::RandomradiusVariation(Axon &axon)
{

    double prev_radius = axon.outer_spheres[axon.outer_spheres.size()-1].radius;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dis(prev_radius, 0.5);

    double random_radius = dis(gen);
    
    if (random_radius < axon.radius-axon.beading_amplitude*axon.radius)
    {
        random_radius = axon.radius-axon.beading_amplitude*axon.radius;
    }
    else if (random_radius > axon.radius+axon.beading_amplitude*axon.radius)
    {
        random_radius = axon.radius+axon.beading_amplitude*axon.radius;
    }
    if (random_radius < min_radius)
    {
        random_radius = min_radius;
    }
    return random_radius;
}


// Axon growth
double AxonGammaDistribution::radiusVariation(Axon &axon)
{
    double mean_radius = axon.radius;
    double length = get_axonal_length(axon);

    double amplitude = mean_radius * axon.beading_amplitude;
    double beading_period = 1;

    double omega = 2 * M_PI / (beading_period*mean_radius);

    double lambda = - omega *axon.phase_shift*beading_period*mean_radius;

    double r = amplitude * sin(omega * length + lambda) + mean_radius;

    if (axon.outer_spheres.size() <= factor){
        /*
        cout << "lambda :" << lambda << endl;
        cout << "omega :" << omega << endl;
        cout << "phase shift :" << axon.phase_shift << endl;
        cout << "axon.outer_spheres.size() :" << axon.outer_spheres.size() << endl;
        cout << " angle : " << (omega * length + lambda)/M_PI << endl;
        */
        
        double new_radius = amplitude * sin(lambda) + mean_radius; // adjust first sphere
        if (new_radius < min_radius)
        {
            axon.ModifyRadiusFirstSphere(min_radius);
        }
        else if (new_radius < axon.outer_spheres[0].radius){
            axon.ModifyRadiusFirstSphere(new_radius);
        }
        
    } 
    

    if (r < min_radius)
    {
        r = min_radius;
    }


    return r;
}

bool AxonGammaDistribution::shrinkRadius(AxonGrowth &growth, const double &radius_to_shrink, Axon &axon)
{

    bool can_grow_;
    Eigen::Vector3d position_that_worked;
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

        while (!can_grow_ && rad > min_radius)
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
        //double angle = acos(v.dot(z)/(v.norm()));
        //cout <<"axon :" << final_axons[i].id << " angle : " << angle << " was supposed to be :"<< final_axons[i].angle << endl;

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
        std::random_device rd;
        std::mt19937 gen(rd());
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
