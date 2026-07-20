#include "CaterpillarGrowth.h"
#include "Glial.h"
#include "grow_axons.h"
#include "grow_glial_cells.h"
#include "grow_blood_vessels.h"
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
#include <limits>
#include <cstdlib>



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

CaterpillarGrowth::CaterpillarGrowth(const Parameters &params, const Eigen::Vector3d &min_l, const Eigen::Vector3d &max_l)
{
    std::random_device rd;
    gen.seed(rd());

    // Axon ICVFs
    target_axons_wo_myelin_icvf = params.axons_wo_myelin_icvf;
    target_axons_w_myelin_icvf = params.axons_w_myelin_icvf;
    target_axons_icvf = params.axons_wo_myelin_icvf + params.axons_w_myelin_icvf;
    
    // Glial ICVFs
    target_glial_pop1_soma_icvf = params.glial_pop1_soma_icvf;
    target_glial_pop1_processes_icvf = params.glial_pop1_processes_icvf; 
    target_glial_pop2_soma_icvf = params.glial_pop2_soma_icvf;
    target_glial_pop2_processes_icvf = params.glial_pop2_processes_icvf;
    
    target_blood_vessels_icvf = params.blood_vessels_icvf;
    
    nbr_axons_populations = params.nbr_axons_populations;
    expanded_for_glial_space = 0;

    // Initialize trackers to 0
    glial_pop1_soma_icvf = 0.0;
    glial_pop1_processes_icvf = 0.0;
    glial_pop2_soma_icvf = 0.0;
    glial_pop2_processes_icvf = 0.0;
    myelin_icvf = 0.0;
    axons_icvf = 0.0;
    extracellular_icvf = 0.0;
    blood_vessels_icvf = 0.0;

    crossing_fibers_type = params.crossing_fibers_type;
    swelling_factor = params.swelling_factor;
    
    // Glial morphology
    glial_pop1_nbr_primary_processes = params.glial_pop1_nbr_primary_processes; // Assuming you add these to Parameters.h
    glial_pop2_nbr_primary_processes = params.glial_pop2_nbr_primary_processes;
    mean_glial_pop1_process_length = params.mean_glial_pop1_process_length;
    std_glial_pop1_process_length = params.std_glial_pop1_process_length;
    mean_glial_pop2_process_length = params.mean_glial_pop2_process_length;
    std_glial_pop2_process_length = params.std_glial_pop2_process_length;
    // Note: You may need to add pop2 mean/std lengths to your Parameters.h!
    
    glial_pop1_radius_mean = params.glial_pop1_radius_mean;
    glial_pop1_radius_std = params.glial_pop1_radius_std;
    glial_pop2_radius_mean = params.glial_pop2_radius_mean;
    glial_pop2_radius_std = params.glial_pop2_radius_std;
    glial_pop1_branching = params.glial_pop1_branching;
    glial_pop2_branching = params.glial_pop2_branching;

    // blood vessels
    epsilon_blood_vessels = params.epsilon_blood_vessels;
    mean_vessel_rad = params.mean_vessel_rad;
    std_vessel_rad = params.std_vessel_rad;
    target_blood_vessels_processes_icvf = params.blood_vessels_processes_icvf;
    blood_vessels_processes_icvf = 0.0;

    // Axon morphology
    alpha = params.alpha;
    beta = params.beta;
    cosPhiSquared = params.cosPhiSquared;
    min_radius = params.min_rad;
    regrow_thr = params.regrow_thr;
    beading_amplitude = params.beading_amplitude;
    beading_std = params.beading_std;
    epsilon = params.epsilon;
    undulation_factor = params.undulation_factor;
    spheres_overlap_factor = params.spheres_overlap_factor;
    axon_can_shrink = params.axon_can_shrink;

    // Myelin parameters
    c1 = params.c1;
    c2 = params.c2;
    c3 = params.c3;

    // Spatial boundaries
    min_limits = min_l;
    max_limits = max_l;
    total_volume = (max_l[0] - min_l[0]) * (max_l[1] - min_l[1]) * (max_l[2] - min_l[2]);
    
    nbr_threads = params.nbr_threads;

    // Clear vectors
    axons.clear();
    glial_pop1.clear();
    glial_pop2.clear();
    blood_vessels.clear();

    // Size grid voxels off the largest sphere radius among populations actually
    // being grown in this simulation. A disabled population's mean radius has
    // no bearing on how densely this grid gets packed, so it's excluded rather
    // than folded into the max (which would make cells needlessly coarse for
    // whatever IS being grown) -- e.g. a large but disabled glial population
    // shouldn't oversize the grid for an axon-only run.
    std::vector<double> active_radii;
    if (target_axons_wo_myelin_icvf > 0.0 || target_axons_w_myelin_icvf > 0.0) {
        // min_radius is only ever used as a floor/clamp on drawn axon radii
        // (see e.g. the "if (r < min_radius) r = min_radius;" clamp below), not
        // a representative size -- axon radii are actually drawn from
        // Gamma(alpha, beta), whose mean is alpha*beta.
        active_radii.push_back(alpha * beta);
    }
    if (target_glial_pop1_soma_icvf > 0.0) {
        active_radii.push_back(glial_pop1_radius_mean);
    }
    if (target_glial_pop2_soma_icvf > 0.0) {
        active_radii.push_back(glial_pop2_radius_mean);
    }
    if (target_blood_vessels_icvf > 0.0) {
        active_radii.push_back(mean_vessel_rad);
    }

    if (active_radii.empty()) {
        grid_voxel_size = 1.0;
    } else {
        grid_voxel_size = 2.0 * *std::max_element(active_radii.begin(), active_radii.end());
    }
    if (grid_voxel_size <= 0.0) {
        grid_voxel_size = 1.0;
    }
    sphere_grid = SphereGrid(min_limits, max_limits, grid_voxel_size);
    // Persistent scratch grid for SanityCheck's per-sub-batch reconciliation
    // (see below): built once here and .clear()'d before each use, rather
    // than reconstructing (and re-allocating its full nx*ny*nz voxel array)
    // every single sub-batch/round during growth.
    batch_scratch_grid = SphereGrid(min_limits, max_limits, grid_voxel_size);

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



bool CaterpillarGrowth::get_begin_end_point(Eigen::Vector3d &Q, Eigen::Vector3d &D, double &angle)
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

double CaterpillarGrowth::myelin_thickness(const double &inner_radius){
    return (c1 + c2 * 2.0 * inner_radius + c3 * log(2.0 * inner_radius));
}  
// Set substrate attributes
void CaterpillarGrowth::generate_radii(std::vector<double> &radii_, std::vector<bool> &has_myelin)
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
            // True target radii, matching the sampled Gamma(alpha,beta)
            // distribution -- seedAllAxons applies swelling_factor itself
            // when actually placing each axon, so this sampling/sorting
            // pass stays based on real target sizes (right axon count and
            // shape for the requested ICVF), independent of how much
            // smaller they're grown before the post-growth swelling pass
            // brings them back up.
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

bool CaterpillarGrowth::PlaceAxon(const int &axon_id, const double &seed_radius, const double &growth_radius, const Eigen::Vector3d &Q, const Eigen::Vector3d &D, std::vector<Axon> &new_axons, const bool &has_myelin, const double &angle_, const bool &outside_voxel)
{

    // The axon's own radius field drives every subsequent AddOneSphere call
    // (step size, beading center) -- growth_radius (swelling_factor * true
    // target, < seed_radius when swelling_factor < 1) gives real, sustained
    // lateral room to wander through the whole depth. The seed sphere
    // itself uses seed_radius (the true, un-shrunk target) so Phase A's
    // circle packing reserves genuine target-density spacing between
    // neighbors from the very first cross-section, rather than seeding
    // (and so starting every axon's depth growth from) an already
    // under-packed layout.
    Axon ax = Axon(axon_id, Q, D, growth_radius, beading_amplitude, beading_std, undulation_factor, has_myelin, angle_, outside_voxel); // axons for regrow batch

    if (has_myelin) {
        double inner_radius = findInnerRadius(growth_radius);
        ax.inner_radius = inner_radius;
    }

    Sphere sphere = Sphere(0, ax.id, axon_constant, Q, seed_radius);

    bool no_overlap = sphere_grid.canSpherebePlaced(sphere);

    if (no_overlap)
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

bool CaterpillarGrowth::collideswithOtherBranches(const Sphere &sph, const Glial &glial_cell_to_grow)
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
bool CaterpillarGrowth::check_borders(const Eigen::Vector3d&  min_l, const Eigen::Vector3d&  max_l, const Eigen::Vector3d& pos, const double& distance_to_border) {

    
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



// Main function to perform the analysis with parallel threads
void CaterpillarGrowth::ICVF(const std::vector<Axon> &axs, const std::vector<Glial> &glial_pop1, const std::vector<Glial> &oligos, const std::vector<Blood_Vessel> &blood_vessels) {


    axons_w_myelin_icvf = 0.0;
    axons_wo_myelin_icvf = 0.0;
    for (const auto &axon : axs) {
        if (axon.myelin_sheath) {
            axons_w_myelin_icvf += axon.volume;
        } else {
            axons_wo_myelin_icvf += axon.volume;
        }
    }
    glial_pop1_processes_icvf = 0.0;
    for (const auto &glial : glial_pop1) {
        glial_pop1_processes_icvf += glial.volume_processes;
    }

    glial_pop1_soma_icvf = 0.0;
    for (const auto &glial : glial_pop1) {
        glial_pop1_soma_icvf += glial.soma.sphereBoxIntersectionVolume(min_limits, max_limits, /*eps_rel=*/1e-6);
    }
    glial_pop2_processes_icvf = 0.0;
    for (const auto &glial : oligos) {
        glial_pop2_processes_icvf += glial.volume_processes;
    }

    glial_pop2_soma_icvf = 0.0;
    for (const auto &glial : oligos) {
        glial_pop2_soma_icvf += glial.soma.sphereBoxIntersectionVolume(min_limits, max_limits, /*eps_rel=*/1e-6);
    }

    blood_vessels_icvf = 0.0;
    for (const auto &bv : blood_vessels) {
        blood_vessels_icvf += bv.volume;
    }

    blood_vessels_processes_icvf = 0.0;
    for (const auto &bv : blood_vessels) {
        blood_vessels_processes_icvf += bv.volume_processes;
    }

    axons_w_myelin_icvf = axons_w_myelin_icvf / total_volume;
    axons_wo_myelin_icvf = axons_wo_myelin_icvf / total_volume;
    axons_icvf = axons_w_myelin_icvf + axons_wo_myelin_icvf;
    glial_pop1_processes_icvf = glial_pop1_processes_icvf / total_volume;
    glial_pop1_soma_icvf = glial_pop1_soma_icvf / total_volume;
    glial_pop2_processes_icvf = glial_pop2_processes_icvf / total_volume;
    glial_pop2_soma_icvf = glial_pop2_soma_icvf / total_volume;
    blood_vessels_processes_icvf = blood_vessels_processes_icvf / total_volume;
    extracellular_icvf = 1 - (axons_w_myelin_icvf + axons_wo_myelin_icvf + glial_pop1_processes_icvf + glial_pop1_soma_icvf + glial_pop2_processes_icvf + glial_pop2_soma_icvf + blood_vessels_processes_icvf);
    blood_vessels_icvf = blood_vessels_icvf / total_volume;
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
double CaterpillarGrowth::c2toKappa(double c2_target,
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

// Builds a monotonic (mu, CDF) table for the axial Watson(kappa) distribution
// via ONE cumulative trapezoidal pass over exp(t^2) from 0 to sqrt(kappa),
// instead of a fresh 1000-point erfi() integration per query. kappa is fixed
// for an entire axon population, so this only needs to run once per
// GrowAllAxons() call -- see sample_from_cdf_table, which then inverts it by
// binary search + linear interpolation for each individual axon, replacing
// what used to be a 20-iteration Newton solve (each iteration itself paying
// for a fresh erfi() integral) per axon.
static void build_watson_cdf_table(double kappa, int n_points,
                                    std::vector<double> &mu_grid,
                                    std::vector<double> &cdf_grid) {
    mu_grid.assign(n_points + 1, 0.0);
    cdf_grid.assign(n_points + 1, 0.0);

    if (kappa <= 0.0) {
        // isotropic axial: F(mu) = mu
        for (int j = 0; j <= n_points; ++j) {
            double mu = double(j) / n_points;
            mu_grid[j] = mu;
            cdf_grid[j] = mu;
        }
        return;
    }

    const double rt = std::sqrt(kappa);
    // Cumulative trapezoidal integral of exp(t^2) from 0 to rt, sampled at
    // the same grid used for mu (t = rt*mu): one pass yields I(rt*mu_j) for
    // every grid point at once (the (2/sqrt(pi)) erfi prefactor cancels in
    // the CDF ratio below, so it's never needed).
    std::vector<double> I(n_points + 1, 0.0);
    double h = rt / n_points;
    double prev_f = 1.0; // exp(0^2)
    for (int j = 1; j <= n_points; ++j) {
        double t = j * h;
        double f = std::exp(t * t);
        I[j] = I[j - 1] + 0.5 * (prev_f + f) * h;
        prev_f = f;
    }

    double total = I[n_points];
    for (int j = 0; j <= n_points; ++j) {
        mu_grid[j] = double(j) / n_points;
        cdf_grid[j] = (total > 0.0) ? (I[j] / total) : mu_grid[j];
    }
}

// Inverts a monotonic CDF table (see build_watson_cdf_table) for a single
// uniform draw u via binary search + linear interpolation -- O(log
// n_points), no exp()/erfi() calls at sampling time.
static double sample_from_cdf_table(const std::vector<double> &mu_grid,
                                     const std::vector<double> &cdf_grid,
                                     double u) {
    auto it = std::lower_bound(cdf_grid.begin(), cdf_grid.end(), u);
    if (it == cdf_grid.begin()) return mu_grid.front();
    if (it == cdf_grid.end()) return mu_grid.back();
    size_t idx = static_cast<size_t>(it - cdf_grid.begin());
    double c0 = cdf_grid[idx - 1], c1 = cdf_grid[idx];
    double m0 = mu_grid[idx - 1], m1 = mu_grid[idx];
    double t = (c1 > c0) ? (u - c0) / (c1 - c0) : 0.0;
    return m0 + t * (m1 - m0);
}


Eigen::Vector3d CaterpillarGrowth::randomPointOnPlane(const Eigen::Vector3d &begin, const Eigen::Vector3d &end, const int &axis1, const int &axis2, const int &axis3, double &angle, bool &outside_voxel) {
 
    double L = (end - begin).norm();

    double phi = angle;
    outside_voxel = false;

    // Reuse the shared, seedable generator (this function is only ever
    // called from seedAllAxons's sequential per-axon loop, never in
    // parallel) instead of a fresh std::random_device-seeded engine per
    // call: that broke run-to-run reproducibility under a fixed seed (every
    // other random decision in a run goes through gen) and paid for a full
    // Mersenne Twister reseed on every single axon for no benefit.
    std::uniform_real_distribution<double> dist(0, 2*M_PI);

    double d = L * std::tan(phi);
    Eigen::Vector3d new_end = end;

    if (angle > 0.2*M_PI) {
        Eigen::Vector3d center_plane = {0,0,0};
        center_plane[axis1] = (max_limits[axis1]/2+3*min_limits[axis1]/2);
        center_plane[axis2] = (max_limits[axis2]/2+3*min_limits[axis2]/2);
        Eigen::Vector3d vector_to_center = center_plane - begin;
        vector_to_center.normalize();
        double theta = dist(gen);
        new_end[axis1] = end[axis1] + d * std::cos(theta);
        new_end[axis2] = end[axis2] + d * std::sin(theta);
        double cos = (new_end-begin).dot(vector_to_center)/(new_end-begin).norm();
        int tries = 0;
        // try to make direction towards center of the voxel
        while (cos < 0) {
            theta = dist(gen);
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
        double theta = dist(gen);
        new_end[axis1] = end[axis1] + d * std::cos(theta);
        new_end[axis2] = end[axis2] + d * std::sin(theta);
        int tries = 0;
        while(new_end[axis1] < min_limits[axis1] || new_end[axis1] > max_limits[axis1] || new_end[axis2] < min_limits[axis2] || new_end[axis2] > max_limits[axis2]) {
            theta = dist(gen);
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



void CaterpillarGrowth::seedAllAxons(std::vector<double> &radii_, std::vector<Axon> &new_axons, std::vector<bool> &has_myelin, std::vector<double> &angles)
{
    // Place every axon's seed directly at its full target radius, one at a
    // time, largest-first (generate_radii's existing sort order). Each
    // axon's search budget scales with (mean_radius/target_radius)^2 -- a
    // fixed budget (the old plane_area*10 for every axon regardless of size)
    // under-searches for a mid-size axon squeezing into an already-crowded
    // plane late in the sequence, so smaller-than-average axons get
    // proportionally more tries, larger ones fewer.
    //
    // If an axon still can't find a valid spot after that budget, we do NOT
    // discard it outright (discarding preferentially loses whichever sizes
    // are hardest to place late in a crowded pass -- almost always the
    // smaller ones, since we go largest-first -- which would systematically
    // skew the realized population away from the sampled Gamma(alpha,beta)
    // shape) and we do NOT scrap the whole pass and restart everything
    // either (restarting the full batch on any single late failure scales
    // exponentially badly with axon count, since the odds every single axon
    // succeeds in one pass shrinks like p^N). Instead we redraw a fresh
    // radius from the same Gamma(alpha,beta) distribution for that slot and
    // retry -- the realized population stays an honest i.i.d. sample from
    // the intended distribution (just resampled when a given draw didn't
    // fit), and only a slot that fails even after max_resamples redraws is
    // finally discarded.
    const double plane_area = (max_limits[0] - min_limits[0]) * (max_limits[1] - min_limits[1]);
    const double mean_radius = alpha * beta;
    const int max_resamples = 50;

    std::gamma_distribution<double> gamma_dist(alpha, beta);

    new_axons.clear();
    new_axons.reserve(radii_.size());

    std::vector<int> indices_to_erase;
    int axon_index = 0;
    int axons_resampled = 0; // how many axons needed >=1 resample (diagnostic only)
    long long total_resamples = 0;

    for (size_t i = 0; i < radii_.size(); ++i) {
        // target_radius (true target, un-shrunk) sizes the seed sphere and
        // Phase A's placement/overlap check, so the cross-section is seeded
        // at genuine target density -- growth_radius (swelling_factor *
        // target_radius; 1.0 = no change, backward compatible) is what the
        // Axon actually grows at from then on, giving a tortuous,
        // undulating path real, sustained lateral room for its *entire*
        // depth. The final swelling pass (SwellAxons) then expands every
        // grown sphere back toward its own true target size afterward,
        // non-uniformly -- each sphere grows as far as its own local room
        // allows, rather than every sphere being offered the same
        // percentage and failing outright if that doesn't fit. See
        // PlaceAxon.
        double target_radius = radii_[i];
        double growth_radius = target_radius * swelling_factor;
        const bool has_myelin_ = has_myelin[i];
        double angle = angles[i];

        bool placed = false;
        bool this_axon_resampled = false;
        for (int resample = 0; resample <= max_resamples && !placed; ++resample) {
            const long int tries_threshold = std::max<long int>(
                100,
                static_cast<long int>(plane_area * 10.0 * (mean_radius * mean_radius) / (target_radius * target_radius)));

            int tries = 0;
            while (!placed && tries < tries_threshold) {
                Vector3d Q, D;
                bool outside_voxel = !get_begin_end_point(Q, D, angle);
                placed = PlaceAxon(axon_index, target_radius, growth_radius, Q, D, new_axons, has_myelin_, angle, outside_voxel);
                if (!placed) {
                    ++tries;
                }
            }

            if (!placed && resample < max_resamples) {
                double resampled;
                int reject_tries = 0;
                do {
                    resampled = gamma_dist(gen);
                    ++reject_tries;
                } while (resampled <= min_radius && reject_tries < 1000);
                if (resampled <= min_radius) {
                    resampled = min_radius * 1.01;
                }
                if (has_myelin_) {
                    resampled += myelin_thickness(resampled);
                }
                target_radius = resampled;
                growth_radius = target_radius * swelling_factor;
                this_axon_resampled = true;
                ++total_resamples;
            }
        }

        if (this_axon_resampled) {
            ++axons_resampled;
        }

        if (!placed) {
            std::cout << "Could not place axon after " << max_resamples
                       << " resamples, discarding" << std::endl;
            indices_to_erase.push_back(static_cast<int>(i));
            continue;
        }

        // Keep whatever radius was actually placed (may be a resampled
        // replacement, not the original draw) so downstream bookkeeping
        // (ICVF target tracking, myelin thickness) matches reality.
        radii_[i] = target_radius;

        // Register the seed in sphere_grid immediately so later seeds in
        // this same pass actually see it.
        sphere_grid.insert(new_axons.back().outer_spheres[0]);
        ++axon_index;
    }

    for (auto it = indices_to_erase.rbegin(); it != indices_to_erase.rend(); ++it) {
        radii_.erase(radii_.begin() + *it);
        has_myelin.erase(has_myelin.begin() + *it);
        angles.erase(angles.begin() + *it);
    }

    if (std::getenv("VERIFY_RADII_DIST") != nullptr) {
        cout << "DEBUG seedAllAxons: axons_resampled=" << axons_resampled
             << " / " << (radii_.size() + indices_to_erase.size())
             << " total_resample_events=" << total_resamples
             << " discarded=" << indices_to_erase.size() << endl;
    }
}



std::vector<double> CaterpillarGrowth::generate_angles(const int &num_samples){

    std::vector<double> angles;
    angles.reserve(num_samples);

    // kappa is fixed for the whole population, so the inverse-CDF table
    // (see build_watson_cdf_table) is built once here rather than per axon.
    std::vector<double> mu_grid, cdf_grid;
    build_watson_cdf_table(kappa, 2000, mu_grid, cdf_grid);

    std::uniform_real_distribution<double> U(0.0, 1.0);
    for (int i = 0; i < num_samples; i++)
    {
        double mu = sample_from_cdf_table(mu_grid, cdf_grid, U(gen));
        angles.push_back(std::acos(mu));
    }
    return angles;
}

void CaterpillarGrowth::GrowAllAxons(){
    
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
        double pre_seed_mean = 0.0, pre_seed_std = 0.0;
        if (std::getenv("VERIFY_RADII_DIST") != nullptr) {
            for (double r : radii) pre_seed_mean += r;
            pre_seed_mean /= radii.size();
            for (double r : radii) pre_seed_std += (r - pre_seed_mean) * (r - pre_seed_mean);
            pre_seed_std = std::sqrt(pre_seed_std / radii.size());
            cout << "DEBUG pre-seed radii: n=" << radii.size()
                 << " mean=" << pre_seed_mean << " std=" << pre_seed_std
                 << " (target gamma mean=" << (alpha * beta)
                 << " std=" << (std::sqrt(alpha) * beta) << ")" << endl;
        }

        // Phase A: place every axon's seed directly at its (possibly
        // swelling_factor-shrunk) radius via 2D circle packing, fully
        // resolved before any 3D threaded growth starts.
        seedAllAxons(radii, axons, has_myelin, angles);

        if (std::getenv("VERIFY_RADII_DIST") != nullptr) {
            double post_mean = 0.0, post_std = 0.0;
            for (double r : radii) post_mean += r;
            post_mean /= radii.size();
            for (double r : radii) post_std += (r - post_mean) * (r - post_mean);
            post_std = std::sqrt(post_std / radii.size());
            cout << "DEBUG post-seed radii: n=" << radii.size()
                 << " mean=" << post_mean << " std=" << post_std << endl;
        }

        {
            double seed_area = 0.0;
            for (const auto &ax : axons) {
                if (ax.outer_spheres.empty()) continue;
                seed_area += M_PI * ax.outer_spheres[0].radius * ax.outer_spheres[0].radius;
            }
            int ga = axons.empty() ? 2 : axons[0].growth_axis;
            int d1 = (ga + 1) % 3, d2 = (ga + 2) % 3;
            double plane_area = (max_limits[d1] - min_limits[d1]) * (max_limits[d2] - min_limits[d2]);
            cout << "DEBUG Phase A seed area fraction: " << (seed_area / plane_area)
                 << " (target icvf " << target_axons_icvf << ")" << endl;
        }

        if (std::getenv("PHASE_A_ONLY") != nullptr) {
            std::exit(0);
        }

        // Phase B: thread-pool growth of the pre-seeded axons, in shallow
        // depth-layers so no axon's tortuous wandering can claim uncontested
        // territory far ahead of axons that haven't started growing yet. No
        // relocate-and-retry of stuck axons -- growAxonsLayered discards any
        // axon that never manages to grow beyond its seed instead.
        growAxonsLayered(radii, has_myelin, angles);

        ICVF(axons, glial_pop1, glial_pop2, blood_vessels);

        if (beading_amplitude == 0 || beading_std == 0){
            return;
        }

        cout << "ICVF axons :" << axons_icvf << endl;

        // Non-uniform final swelling: each sphere grows to its own local
        // maximum (via bisection, bounded by its true target radius) in
        // one shot per round, instead of every sphere being offered the
        // same percentage and failing outright if that doesn't fit -- a
        // tightly-pinched sphere no longer has to wait for the whole
        // population's growth rate to shrink down to its own limit before
        // it gets any growth at all.
        SwellAxons();
        cout << "new ICVF " << axons_icvf << endl;

        // Post-swelling cleanup: independent per-sphere growth (especially
        // the uncapped/push-assisted paths) can leave a sphere entirely
        // inside a neighbor along the same axon's own chain. checkNoCollisions
        // never catches this -- it deliberately skips same-object comparisons,
        // since a branch is supposed to touch its own trunk -- but such a
        // sphere presents no obstacle surface of its own (anything that could
        // touch it would already have touched the bigger neighbor containing
        // it first), so it's just dead weight in the output. Drop it and
        // recompute volume/ICVF for any axon that changed.
        int total_engulfed = 0;
        bool any_axon_changed = false;
        for (auto &ax : axons) {
            std::vector<Sphere> removed;
            ax.removeEngulfedSpheres(removed);
            if (!removed.empty()) {
                for (const auto &sph : removed) {
                    sphere_grid.remove(sph);
                }
                ax.update_Volume(spheres_overlap_factor, min_limits, max_limits);
                total_engulfed += static_cast<int>(removed.size());
                any_axon_changed = true;
            }
        }
        if (any_axon_changed) {
            ICVF(axons, glial_pop1, glial_pop2, blood_vessels);
        }
        cout << "Removed " << total_engulfed << " fully-engulfed sphere(s) after swelling; ICVF now " << axons_icvf << endl;
    }

}


// Read-only: finds the largest radius, up to cap_radius, that sph could
// grow to this round (offered up to `percentage` of its current radius,
// clipped to whatever actually fits via a short bisection) -- doesn't
// mutate sph or sphere_grid, only queries it, so it's safe to call from
// multiple threads at once as long as none of those spheres collide-check
// against each other (SphereGrid::canSpherebePlaced already excludes a
// sphere's own axon from collision, so spheres belonging to the same axon
// satisfy this -- see SwellAxon, the only caller). Returns sph.radius
// unchanged if it can't grow at all this round.
double CaterpillarGrowth::ComputeSwollenRadius(const Sphere &sph, const double &percentage, const double &cap_radius) const {

    if (cap_radius <= sph.radius) {
        return sph.radius;
    }

    double target = std::min(sph.radius + percentage * sph.radius, cap_radius);
    Sphere target_candidate(sph.id, sph.object_id, sph.object_type, sph.center, target);
    if (sphere_grid.canSpherebePlaced(target_candidate)) {
        return target;
    }

    double lo = sph.radius, hi = target;
    for (int iter = 0; iter < 10; ++iter) {
        double mid = (lo + hi) / 2.0;
        Sphere trial(sph.id, sph.object_id, sph.object_type, sph.center, mid);
        if (sphere_grid.canSpherebePlaced(trial)) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return lo;
}

bool CaterpillarGrowth::SwellAxon(Axon &ax, const double &percentage, const std::vector<double> &caps, ThreadPool &pool) {

    size_t n = ax.outer_spheres.size();

    // Compute every sphere's achievable new radius in parallel -- spheres
    // within this one axon never collide-check against each other (see
    // ComputeSwollenRadius), and no other axon is being modified while
    // this axon is being processed, so the read-only queries below don't
    // race against each other. sphere_grid itself is only ever mutated
    // afterward, sequentially.
    std::vector<std::future<double>> futures;
    futures.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        futures.emplace_back(pool.enqueueTask(
            [this, &ax, &percentage, &caps, i]() {
                return ComputeSwollenRadius(ax.outer_spheres[i], percentage, caps[i]);
            }
        ));
    }

    // Join every future BEFORE touching sphere_grid below: .get() only
    // blocks for its own future, so mutating the grid for sphere i inside
    // the same loop that's still waiting on sphere i+1's future would let
    // that mutation race against a still-in-flight query on another
    // thread. Collecting all results first, then mutating only once
    // everything has finished, avoids that entirely.
    std::vector<double> new_radii(n);
    for (size_t i = 0; i < n; ++i) {
        new_radii[i] = futures[i].get();
    }

    std::vector<Sphere> new_spheres;
    new_spheres.reserve(n);
    bool any_changed = false;
    for (size_t i = 0; i < n; ++i) {
        double new_radius = new_radii[i];
        if (new_radius > ax.outer_spheres[i].radius) {
            Sphere sph = ax.outer_spheres[i];
            sph.radius = new_radius;
            // swap the grid entry for this sphere
            sphere_grid.remove(ax.outer_spheres[i]);
            sphere_grid.insert(sph);
            new_spheres.push_back(sph);
            any_changed = true;
        } else {
            // if it couldn't grow, push the sphere back unchanged (grid entry unchanged)
            new_spheres.push_back(ax.outer_spheres[i]);
        }
    }
    ax.outer_spheres = new_spheres;
    return any_changed;
}


// Read-only: finds the neighbor causing the deepest overlap with sph at
// trial_radius (own axon's spheres excluded, matching canSpherebePlaced's
// same-object exclusion), and returns a push vector -- perpendicular to
// growth_axis, so it never changes how far along the axon this sphere sits
// -- that moves sph away from that neighbor. Mirrors AddOneSphere's
// findPush (grow_axons.cpp), applied here during the final swelling pass
// instead of initial growth. Returns false (push/blocker_radius untouched)
// if nothing overlaps at trial_radius, or the only overlap has no lateral
// escape direction.
bool CaterpillarGrowth::FindSwellPush(const Sphere &sph, const double &trial_radius, int growth_axis, Eigen::Vector3d &push, double &blocker_radius) const {

    // Deliberately avoids sphere_grid.query() -- same reasoning as
    // canSpherebePlaced (SphereGrid.cpp): iterating voxels/entries directly
    // via findWorstOverlap skips the heap-allocated intermediate vector.
    // Unlike canSpherebePlaced this can't early-exit (every candidate must
    // be compared to find the *worst* overlap), but this is still called
    // often enough during the push-assisted swelling fallback that
    // avoiding the allocation matters.
    Eigen::Vector3d blocker_center;
    bool found = sphere_grid.findWorstOverlap(sph.center, trial_radius, trial_radius * 2.0,
                                               sph.object_type, sph.object_id,
                                               blocker_center, blocker_radius);
    if (!found) {
        return false;
    }
    double worst_overlap = (trial_radius + blocker_radius) - (blocker_center - sph.center).norm();

    Eigen::Vector3d away = sph.center - blocker_center;
    away[growth_axis] = 0.0;
    double away_norm = away.norm();
    if (away_norm < 1e-9) {
        return false;
    }
    away /= away_norm;
    const double push_margin = 1e-3 * trial_radius;
    push = away * (worst_overlap + push_margin);
    return true;
}

// Read-only: like ComputeSwollenRadius, but if growing in place at the
// current center can't reach this round's target (something specific is
// blocking it), also tries pushing sph away from whichever neighbor is the
// worst blocker at that target radius -- capped so the resulting
// center-to-center distance to that neighbor never exceeds
// max(sph.radius, blocker_radius)/swelling_factor, keeping the push
// bounded to roughly one true-target-diameter's worth of clearance rather
// than an unbounded displacement. Only used as a fallback phase (see
// SwellAxons) once the plain in-place pass has already converged without
// reaching target_axons_icvf. Never mutates sph or sphere_grid -- safe to
// call concurrently across spheres of the same axon (see SwellAxonWithPush).
Sphere CaterpillarGrowth::ComputeSwollenSphere(const Sphere &sph, const double &percentage, const double &cap_radius, int growth_axis,
                                                const Sphere *prev_sphere, const Sphere *next_sphere) const {

    double target = std::min(sph.radius + percentage * sph.radius, cap_radius);
    if (target <= sph.radius) {
        return sph;
    }

    Sphere target_candidate(sph.id, sph.object_id, sph.object_type, sph.center, target);
    if (sphere_grid.canSpherebePlaced(target_candidate)) {
        return Sphere(sph.id, sph.object_id, sph.object_type, sph.center, target);
    }

    double lo = sph.radius, hi = target;
    for (int iter = 0; iter < 10; ++iter) {
        double mid = (lo + hi) / 2.0;
        Sphere trial(sph.id, sph.object_id, sph.object_type, sph.center, mid);
        if (sphere_grid.canSpherebePlaced(trial)) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    double best_radius = lo;
    if (best_radius >= target - 1e-9) {
        return Sphere(sph.id, sph.object_id, sph.object_type, sph.center, best_radius);
    }

    // In-place growth was capped short of this round's target -- try
    // pushing away from whichever neighbor is blocking it at that target.
    Eigen::Vector3d push;
    double blocker_radius = 0.0;
    if (!FindSwellPush(sph, target, growth_axis, push, blocker_radius)) {
        return Sphere(sph.id, sph.object_id, sph.object_type, sph.center, best_radius);
    }

    double max_gap = std::max(sph.radius, blocker_radius) / std::max(swelling_factor, 1e-6);
    if (push.norm() > max_gap) {
        push = push.normalized() * max_gap;
    }
    Eigen::Vector3d pushed_center = sph.center + push;

    double lo2 = sph.radius, hi2 = target;
    Sphere pushed_target(sph.id, sph.object_id, sph.object_type, pushed_center, target);
    bool full_target_fits = sphere_grid.canSpherebePlaced(pushed_target);
    if (!full_target_fits) {
        for (int iter = 0; iter < 10; ++iter) {
            double mid = (lo2 + hi2) / 2.0;
            Sphere trial(sph.id, sph.object_id, sph.object_type, pushed_center, mid);
            if (sphere_grid.canSpherebePlaced(trial)) {
                lo2 = mid;
            } else {
                hi2 = mid;
            }
        }
    } else {
        lo2 = target;
    }

    // A push chases external density, not at the cost of separating this
    // sphere from its own immediate chain neighbors -- canSpherebePlaced
    // only ever checks against *other* objects (self is deliberately
    // excluded), so nothing above would catch that on its own. Reject the
    // push outright (fall back to staying in place) rather than accepting a
    // radius that collision-avoids everyone else but opens a hole in this
    // axon's own continuity; checked against each neighbor's current,
    // pre-round state, since every sphere in the axon is being computed
    // concurrently this round (see SwellAxonWithPush).
    auto still_touches_neighbor = [](const Eigen::Vector3d &c, double r, const Sphere *nb) {
        if (!nb) return true;
        return (c - nb->center).norm() <= r + nb->radius;
    };
    bool keeps_continuity = still_touches_neighbor(pushed_center, lo2, prev_sphere) &&
                            still_touches_neighbor(pushed_center, lo2, next_sphere);

    // Only actually move the sphere if the pushed position beat staying put
    // and didn't break continuity with either neighbor.
    if (lo2 > best_radius && keeps_continuity) {
        return Sphere(sph.id, sph.object_id, sph.object_type, pushed_center, lo2);
    }
    return Sphere(sph.id, sph.object_id, sph.object_type, sph.center, best_radius);
}

bool CaterpillarGrowth::SwellAxonWithPush(Axon &ax, const double &percentage, const std::vector<double> &caps, ThreadPool &pool) {

    size_t n = ax.outer_spheres.size();
    int growth_axis = ax.growth_axis;

    std::vector<std::future<Sphere>> futures;
    futures.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        futures.emplace_back(pool.enqueueTask(
            [this, &ax, &percentage, &caps, growth_axis, i, n]() {
                // Read-only neighbor lookups against ax.outer_spheres are
                // safe here: nothing mutates that vector until every future
                // in this round has been joined (see below).
                const Sphere *prev = (i > 0) ? &ax.outer_spheres[i - 1] : nullptr;
                const Sphere *next = (i + 1 < n) ? &ax.outer_spheres[i + 1] : nullptr;
                return ComputeSwollenSphere(ax.outer_spheres[i], percentage, caps[i], growth_axis, prev, next);
            }
        ));
    }

    // Join every future BEFORE touching sphere_grid, same reasoning as
    // SwellAxon: mutating the grid for an earlier sphere while a later
    // one's query is still in flight on another thread would race.
    std::vector<Sphere> computed(n);
    for (size_t i = 0; i < n; ++i) {
        computed[i] = futures[i].get();
    }

    // Post-hoc consistency repair: ComputeSwollenSphere's own continuity
    // check (see its prev_sphere/next_sphere parameters) only ever compares
    // against a neighbor's *pre-round* state, because every sphere in this
    // axon was computed independently and concurrently above -- it has no
    // way to know whether that neighbor ALSO moved this round. Two adjacent
    // pushes can therefore each individually look fine against the other's
    // stale position and still end up not touching once both are applied.
    // Repair that here: walk the chain and, for any consecutive pair that
    // ends up separated, revert whichever of the two actually changed this
    // round (both, if both did) back to its safe, definitely-touching
    // pre-round state -- iterated to a fixed point since a revert can
    // occasionally reopen an already-checked neighboring pair.
    auto changed_this_round = [&](size_t idx) {
        return (computed[idx].center - ax.outer_spheres[idx].center).norm() > 1e-9 ||
               computed[idx].radius > ax.outer_spheres[idx].radius;
    };
    for (int pass = 0; pass < 5 && n > 1; ++pass) {
        bool any_reverted = false;
        for (size_t i = 1; i < n; ++i) {
            double d = (computed[i].center - computed[i - 1].center).norm();
            if (d > computed[i].radius + computed[i - 1].radius) {
                if (changed_this_round(i)) {
                    computed[i] = ax.outer_spheres[i];
                    any_reverted = true;
                }
                if (changed_this_round(i - 1)) {
                    computed[i - 1] = ax.outer_spheres[i - 1];
                    any_reverted = true;
                }
            }
        }
        if (!any_reverted) break;
    }

    std::vector<Sphere> new_spheres;
    new_spheres.reserve(n);
    bool any_changed = false;
    for (size_t i = 0; i < n; ++i) {
        const Sphere &result = computed[i];
        bool changed = (result.radius > ax.outer_spheres[i].radius) ||
                        ((result.center - ax.outer_spheres[i].center).norm() > 1e-9);
        if (changed) {
            sphere_grid.remove(ax.outer_spheres[i]);
            sphere_grid.insert(result);
            new_spheres.push_back(result);
            any_changed = true;
        } else {
            new_spheres.push_back(ax.outer_spheres[i]);
        }
    }
    ax.outer_spheres = new_spheres;
    return any_changed;
}


void CaterpillarGrowth::SwellAxons(){

    // Uncapped for every swelling_factor value: whatever room a sphere has
    // locally, it swells into, chasing the global target_axons_icvf with no
    // per-sphere ceiling -- including for swelling_factor < 1.0, where the
    // sphere had been deliberately grown thin (see seedAllAxons/PlaceAxon)
    // and could in principle be capped back at its own true un-shrunk
    // target. Uncapping that case too trades the shrink-then-regrow
    // approach's tighter radius-distribution guarantee for extra reachable
    // ICVF, at the cost of potentially more push-driven path distortion.
    const bool uncapped = true;

    // Fix every sphere's true target cap ONCE, from its as-grown (i.e.
    // pre-swell) radius: recomputing radius/swelling_factor from an
    // already-partially-swollen radius on a later round would overshoot
    // the true target. The seed sphere (index 0 of every axon) is already
    // at its true, un-shrunk size (see PlaceAxon/seedAllAxons), so its cap
    // is itself -- every later sphere's cap is radius/swelling_factor.
    std::vector<std::vector<double>> caps(axons.size());
    for (size_t a = 0; a < axons.size(); ++a) {
        const auto &spheres = axons[a].outer_spheres;
        caps[a].resize(spheres.size());
        for (size_t i = 0; i < spheres.size(); ++i) {
            caps[a][i] = uncapped
                ? std::numeric_limits<double>::infinity()
                : ((i == 0) ? spheres[i].radius : (spheres[i].radius / swelling_factor));
        }
    }

    // One pool, reused for every axon and every round (spinning up
    // nbr_threads OS threads per axon per round would dwarf the actual
    // work for axons with few spheres). Axons themselves are still
    // processed sequentially -- only the spheres within a single axon are
    // computed in parallel, see SwellAxon.
    ThreadPool pool(nbr_threads);

    // Same gradual, many-round structure as the original uniform scheme
    // (start at a generous percentage, shrink it whenever a round makes
    // negligible progress, stop once it's negligibly small or after 1000
    // rounds): letting every axon grow a little each round, rather than
    // greedily jumping straight to each sphere's own local maximum in
    // whatever order axons happen to be processed, is what keeps this fair
    // across the whole population instead of letting the first-processed
    // axons claim shared space before later ones get a turn. What's new
    // here is that each sphere is also capped at its own true target (via
    // ComputeSwollenRadius/SwellAxon above) and a round's growth step is
    // clipped to whatever actually fits rather than being all-or-nothing,
    // so a tightly-pinched sphere gets partial credit instead of zero
    // progress every round until the population-wide percentage happens to
    // shrink below its own limit.
    int nbr_attempts = 0;
    double percentage_swelling = 0.1;
    const double minimum_percentage_swelling = 1e-6;

    while (axons_icvf < target_axons_icvf && nbr_attempts < 1000) {
        double old_icvf = axons_icvf;
        // Only recompute an axon's volume (O(its sphere count), often
        // hundreds to 1000+) if this round actually changed one of its
        // spheres, and only re-sum the global ICVF if at least one axon
        // anywhere changed -- most axons stop changing well before the
        // loop as a whole converges, so in the later rounds this skips the
        // large majority of otherwise-guaranteed-no-op recomputation.
        bool any_axon_changed = false;
        for (size_t a = 0; a < axons.size(); ++a) {
            if (axons[a].outer_spheres.empty()) {
                continue;
            }
            if (SwellAxon(axons[a], percentage_swelling, caps[a], pool)) {
                axons[a].update_Volume(spheres_overlap_factor, min_limits, max_limits);
                any_axon_changed = true;
            }
        }
        if (any_axon_changed) {
            ICVF(axons, glial_pop1, glial_pop2, blood_vessels);
        }
        cout << "new ICVF " << axons_icvf << " old icvf : " << old_icvf
             << " percentage swelling :" << percentage_swelling << endl;
        if (on_swelling_progress) {
            on_swelling_progress(axons_icvf, target_axons_icvf);
        }
        if ((axons_icvf - old_icvf) < 1e-5 && percentage_swelling >= minimum_percentage_swelling) {
            percentage_swelling = percentage_swelling / 3;
        } else if (percentage_swelling < minimum_percentage_swelling) {
            break;
        }
        ++nbr_attempts;
    }

    // Fallback: only reached if the plain in-place pass above has already
    // converged (percentage shrunk below minimum, or 1000 rounds) without
    // reaching target_axons_icvf. Re-run the same gradual, many-round
    // structure, but now each sphere that's locally blocked may also push
    // itself away from whichever neighbor is causing the tightest pinch
    // (see ComputeSwollenSphere/FindSwellPush) before settling for a
    // smaller radius -- since that's a real repositioning, not just a
    // radius search, it's worth the extra cost only once simple growth in
    // place has already been exhausted.
    if (axons_icvf < target_axons_icvf) {
        cout << "In-place swelling converged short of target (" << axons_icvf
             << " < " << target_axons_icvf << ") -- trying push-assisted swelling" << endl;
        nbr_attempts = 0;
        percentage_swelling = 0.1;
        while (axons_icvf < target_axons_icvf && nbr_attempts < 1000) {
            double old_icvf = axons_icvf;
            bool any_axon_changed = false;
            for (size_t a = 0; a < axons.size(); ++a) {
                if (axons[a].outer_spheres.empty()) {
                    continue;
                }
                if (SwellAxonWithPush(axons[a], percentage_swelling, caps[a], pool)) {
                    axons[a].update_Volume(spheres_overlap_factor, min_limits, max_limits);
                    any_axon_changed = true;
                }
            }
            if (any_axon_changed) {
                ICVF(axons, glial_pop1, glial_pop2, blood_vessels);
            }
            cout << "new ICVF (push) " << axons_icvf << " old icvf : " << old_icvf
                 << " percentage swelling :" << percentage_swelling << endl;
            if (on_swelling_progress) {
                on_swelling_progress(axons_icvf, target_axons_icvf);
            }
            if ((axons_icvf - old_icvf) < 1e-5 && percentage_swelling >= minimum_percentage_swelling) {
                percentage_swelling = percentage_swelling / 3;
            } else if (percentage_swelling < minimum_percentage_swelling) {
                break;
            }
            ++nbr_attempts;
        }
    }
}

// Growing substrate
void CaterpillarGrowth::createSubstrate()
{
    // place glial cells
    cout << "Place Blood Vessels" << endl;
    PlaceBloodVessels();
    GrowBloodVessels();
    cout << "Grow Blood Vessel Branches" << endl;
    GrowBloodVesselBranches();
    ApplyMurraysLawToBloodVessels();
    cout << "Place Glial Cells" << endl;
    PlaceGlialCells();
    cout << "Grow all Axons" << endl;
    GrowAllAxons();
    cout << "Grow Myelin" << endl;
    add_Myelin();
    ICVF(axons, glial_pop1, glial_pop2, blood_vessels);
    cout << "GrowAllGlialCells" << endl;
    GrowAllGlialCells();

    bool cells_ok = checkNoCollisions();

    if (!cells_ok)
    {
        cout << "Final collision check failed" << endl;
    }
    else
    {
        cout << "Final collision check passed" << endl;
    }

    ICVF(axons, glial_pop1, glial_pop2, blood_vessels);
}

void CaterpillarGrowth::PlaceBloodVessels(){

    
    double achieved_icvf = 0.0;
    const int MAX_GLOBAL_FAILS = 10000;
    bool placed = false;

    blood_vessels.clear();

    std::normal_distribution<> dis_blood_vessel(mean_vessel_rad, std_vessel_rad);

    std::uniform_real_distribution<> base_x1(min_limits[0], max_limits[0]);
    std::uniform_real_distribution<> base_y1(min_limits[1], max_limits[1]);

    int global_fail_count = 0;

    double occupancy_many_branches = 0.5;

    while(achieved_icvf < target_blood_vessels_icvf*occupancy_many_branches && global_fail_count < MAX_GLOBAL_FAILS){
        Sphere s;
        for (int attempt = 0; attempt < 100; ++attempt) {
            double rad = dis_blood_vessel(gen);
            s = Sphere(
                /*object_id*/ 0,                        // keep your ID scheme if needed
                /*object_index*/ static_cast<int>(blood_vessels.size()),
                /*type*/ blood_constant,
                { base_x1(gen), base_y1(gen), 0},
                rad
            );
            if (sphere_grid.canSpherebePlaced(s)) {
                placed = true;
                break;
            }
        }
        if (!placed) {
            ++global_fail_count;
            continue; // try placing another cell; loop terminates via MAX_GLOBAL_FAILS_1
        }
        else{
            achieved_icvf += M_PI*std::pow(s.radius,2)*(max_limits[2]-min_limits[2])/total_volume;
            int id = blood_vessels.size();
            Eigen::Vector3d begin = s.center;
            Eigen::Vector3d end = s.center + Eigen::Vector3d(0,0,(max_limits[2]-min_limits[2]));
            Blood_Vessel bv (id, begin, end, s.radius, /*beading_amplitude=*/0.0, /*beading_std=*/0.0, /*undulation_factor=*/1);
            bv.add_first_sphere(s);
            blood_vessels.push_back(bv);
        }
    }

    blood_vessels_icvf = achieved_icvf;
}

void CaterpillarGrowth::GrowBloodVessels() {
    bool   bv_can_shrink = false;
    double stuck_radius  = 0.0;
    int    stuck_index   = 0;

    const std::size_t initial_n = blood_vessels.size();

    for (std::size_t j = 0; j < initial_n; ++j) {
        Blood_Vessel bv = blood_vessels[j];

        // Prefer references instead of raw pointers where possible.
        BloodVesselGrowth grow(
            bv, &sphere_grid,
            min_limits, max_limits, min_limits, max_limits, epsilon_blood_vessels, barrier_tickness
        );

        // NOTE: name says "Thread"—ensure it is synchronous here (blocking).
        // If it spawns async work, you must join before using results.
        grow.growthThread(stuck_radius, stuck_index, spheres_overlap_factor, bv_can_shrink);

        blood_vessels[j] = std::move(bv);
        blood_vessels[j].addToGrid(sphere_grid);

        display_progress(static_cast<int>(j), static_cast<int>(initial_n));
    }

    blood_vessels.erase(
        std::remove_if(blood_vessels.begin() + initial_n, blood_vessels.end(),
                    [](const Blood_Vessel& bv) { return bv.ramification_spheres.empty() || bv.ramification_spheres[0].size() <= 1; }),
        blood_vessels.end());

    ICVF(axons, glial_pop1, glial_pop2, blood_vessels);
}

void CaterpillarGrowth::GrowBloodVesselBranches() {

    if (blood_vessels.empty()) {
        return;
    }
    if (target_blood_vessels_processes_icvf <= 0.0) {
        return;
    }

    std::vector<BloodVesselGrowth> growths;
    growths.reserve(blood_vessels.size());
    for (size_t i = 0; i < blood_vessels.size(); ++i) {
        growths.emplace_back(blood_vessels[i], &sphere_grid, min_limits, max_limits,
                             min_limits, max_limits, epsilon_blood_vessels, barrier_tickness);
    }

    // Each vessel starts with only its main vessel (branch 0). Seed a length-tracking
    // vector for it so a branch can emerge off any of its spheres with old_length = 0.
    std::vector<int> nbr_spheres(blood_vessels.size(), 0);
    for (size_t i = 0; i < blood_vessels.size(); ++i) {
        int trunk_size = blood_vessels[i].ramification_spheres[0].size();
        blood_vessels[i].lengths_branches.resize(1);
        blood_vessels[i].lengths_branches[0] = std::vector<double>(trunk_size, 0.0);
        blood_vessels[i].attractors.resize(1); // placeholder for branch 0, never read
        blood_vessels[i].children_branches.resize(1); // branch 0 (trunk) has no children yet
        nbr_spheres[i] = trunk_size;
    }

    ICVF(axons, glial_pop1, glial_pop2, blood_vessels);
    display_progress(blood_vessels_processes_icvf, target_blood_vessels_processes_icvf);
    if (blood_vessels_processes_icvf >= target_blood_vessels_processes_icvf) return;

    int nbr_tries = 0;
    const int max_tries = 100000;
    double prev_icvf = blood_vessels_processes_icvf;
    int stall_count = 0;
    const int stall_limit = 100;
    const double rel_eps = 1e-6;

    while (blood_vessels_processes_icvf < target_blood_vessels_processes_icvf && nbr_tries <= max_tries) {
        for (size_t i = 0; i < growths.size(); ++i) {
            bool grew = growths[i].growBranch(nbr_spheres[i], spheres_overlap_factor);
            blood_vessels[i] = growths[i].bv_to_grow;
            if (grew) {
                // Register the newly grown branch immediately, so other vessels'
                // branches (grown later in this same loop, or in later iterations)
                // can see it instead of growing blind to it.
                for (const auto &sph : blood_vessels[i].ramification_spheres.back()) {
                    sphere_grid.insert(sph);
                }
            }
        }

        if (nbr_tries % 10 == 0) {
            for (auto& bv : blood_vessels) {
                bv.compute_processes_icvf(spheres_overlap_factor, min_limits, max_limits);
            }
            ICVF(axons, glial_pop1, glial_pop2, blood_vessels);
            display_progress(blood_vessels_processes_icvf, target_blood_vessels_processes_icvf);

            if (blood_vessels_processes_icvf >= target_blood_vessels_processes_icvf) break;

            double denom = std::max(1.0, std::abs(prev_icvf));
            if (std::abs(blood_vessels_processes_icvf - prev_icvf) / denom < rel_eps) {
                if (++stall_count > stall_limit) {
                    std::cerr << "Stuck in a loop, stopping blood vessel branch growth!\n";
                    break;
                }
            } else {
                stall_count = 0;
            }
            prev_icvf = blood_vessels_processes_icvf;
        }
        ++nbr_tries;
    }

    if (nbr_tries > max_tries) {
        std::cerr << "Max attempts reached while growing blood vessel branches!\n";
    }

    // No final grid sync needed here: the trunk was added to the grid when it
    // finished growing (in GrowBloodVessels), and each branch was added as soon
    // as it was grown, above.
}

void CaterpillarGrowth::ApplyMurraysLawToBloodVessels() {

    for (auto &bv : blood_vessels) {
        bv.enforceMurraysLaw();
        // some branches may have shrunk below the minimum radius they're meant to
        // decay toward: delete them (and their own sub-branches) rather than keep
        // a vanishingly thin sliver.
        bv.pruneUndersizedBranches(min_limits, max_limits);
        // enforceMurraysLaw only rescales radii (centers never move), so a branch's
        // attachment point to its parent can end up with too little combined radius
        // for the (fixed) distance between them once upstream junctions have
        // cascaded shrinks onto that parent sphere; close any such junction gap
        // before checking/fixing intra-branch spacing below.
        bv.bridgeJunctionGaps(spheres_overlap_factor);
        // shrinking radii above doesn't touch sphere spacing, so a chain that
        // overlapped fine at its original radius can now have visible gaps;
        // reinsert spheres so consecutive spheres never drift too far apart.
        bv.reinterpolateAfterShrink(spheres_overlap_factor);

        // radii/spheres changed: refresh this vessel's own volume figures
        bv.update_Volume(spheres_overlap_factor, min_limits, max_limits);
        bv.compute_processes_icvf(spheres_overlap_factor, min_limits, max_limits);
    }

    // sphere_grid still holds the pre-correction radii; checkNoCollisions() (called
    // later in createSubstrate) rebuilds the grid from scratch, so no sync needed here.
    ICVF(axons, glial_pop1, glial_pop2, blood_vessels);
}

void CaterpillarGrowth::GrowAllGlialCells() {

    if (glial_pop1.size() <= 0.0 && glial_pop2.size() <= 0.0)
    {
        return;
    } 
    if (target_glial_pop1_processes_icvf <= 0.0 && target_glial_pop2_processes_icvf <= 0.0)
    {
        return;
    }

    if (target_glial_pop1_processes_icvf >0.0 && glial_pop1.size() > 0)
    {   
        // Growing extra branches for glial_pop1
        growBranches(1);
    } 
    if (target_glial_pop2_processes_icvf >0.0 &&  glial_pop2.size() > 0)
    {
        // Growing extra branches for glial_pop2
        growBranches(2);
    }
}


void CaterpillarGrowth::growBranches(const int &population_nbr) {

    // Pick the population once
    auto& pop = (population_nbr == 1 ? glial_pop1 : glial_pop2);

    // Per-pop parameters
    double target_icvf = (population_nbr == 1) ? target_glial_pop1_processes_icvf
                                               : target_glial_pop2_processes_icvf;
    double current_icvf = (population_nbr == 1) ? glial_pop1_processes_icvf
                                                : glial_pop2_processes_icvf;

    int    nbr_primary_processes = (population_nbr == 1) ? glial_pop1_nbr_primary_processes
                                                         : glial_pop2_nbr_primary_processes;
    double mean_len  = (population_nbr == 1) ? mean_glial_pop1_process_length
                                             : mean_glial_pop2_process_length;
    double std_len   = (population_nbr == 1) ? std_glial_pop1_process_length
                                             : std_glial_pop2_process_length;

    if (target_icvf <= 0.0) return;

    Eigen::Vector3d extended_min_limits = min_limits - Eigen::Vector3d::Constant(expanded_for_glial_space);
    Eigen::Vector3d extended_max_limits = max_limits + Eigen::Vector3d::Constant(expanded_for_glial_space);
    std::vector<GlialCellGrowth> growths;
    growths.reserve(pop.size());
    for (size_t i = 0; i < pop.size(); ++i) {
        growths.emplace_back(pop[i], &sphere_grid, extended_min_limits, extended_max_limits,
                             min_limits, max_limits, min_radius);
    }

    // Grow first primary branches
    std::vector<int> nbr_spheres(pop.size(), 0);


    // Grow primary branches round-by-round: every cell attempts its j-th
    // primary branch in the same round, in nbr_threads-sized sub-batches of
    // mutually-blind private copies (mirroring growAxonsLayered's
    // layer/sub-batch structure for axons), then reconciled against both
    // the shared sphere_grid and this sub-batch's own siblings before
    // committing -- rather than growing one cell's entire set of primary
    // branches sequentially before moving to the next. Different cells'
    // primary branches emerge from different somas and radiate outward,
    // so collisions between them are the exception rather than the rule;
    // when one does happen, the offending attempt is simply retried (fresh
    // random direction) rather than the whole batch being penalized.
    for (int j = 0; j < nbr_primary_processes; ++j) {
        for (size_t batch_start = 0; batch_start < pop.size(); batch_start += static_cast<size_t>(nbr_threads)) {
            size_t batch_end = std::min(pop.size(), batch_start + static_cast<size_t>(nbr_threads));
            size_t batch_size = batch_end - batch_start;

            const int max_primary_tries = 1000;
            std::vector<bool> done(batch_size, false);
            std::vector<int> tries(batch_size, 0);

            bool any_pending = true;
            while (any_pending) {
                std::vector<Glial> attempt_copies;
                std::vector<int> attempt_nbr_spheres;
                std::vector<size_t> attempt_k; // index within [0, batch_size)
                attempt_copies.reserve(batch_size);
                attempt_nbr_spheres.reserve(batch_size);
                attempt_k.reserve(batch_size);
                for (size_t k = 0; k < batch_size; ++k) {
                    if (done[k]) continue;
                    size_t i = batch_start + k;
                    attempt_copies.push_back(pop[i]);
                    attempt_nbr_spheres.push_back(nbr_spheres[i]);
                    attempt_k.push_back(k);
                }

                std::vector<GlialCellGrowth> temp_growths;
                temp_growths.reserve(attempt_copies.size());
                for (size_t a = 0; a < attempt_copies.size(); ++a) {
                    temp_growths.emplace_back(attempt_copies[a], &sphere_grid, extended_min_limits,
                                              extended_max_limits, min_limits, max_limits, min_radius);
                }

                ThreadPool pool(nbr_threads);
                std::vector<std::future<bool>> futures;
                futures.reserve(temp_growths.size());
                for (size_t a = 0; a < temp_growths.size(); ++a) {
                    futures.emplace_back(pool.enqueueTask(
                        [this, &temp_growths, &attempt_nbr_spheres, a, mean_len, std_len]() {
                            return temp_growths[a].growPrimaryBranch(attempt_nbr_spheres[a], mean_len, std_len, spheres_overlap_factor);
                        }
                    ));
                }
                std::vector<bool> grew(temp_growths.size());
                for (size_t a = 0; a < futures.size(); ++a) {
                    grew[a] = futures[a].get();
                }

                // Reconcile: this round's new branches were grown as private
                // copies invisible to each other, so check every one of them
                // against both the already-committed sphere_grid and a
                // scratch grid of this round's own siblings before
                // committing any of them.
                batch_scratch_grid.clear();
                for (size_t a = 0; a < attempt_copies.size(); ++a) {
                    if (!grew[a]) continue;
                    for (const auto &sph : attempt_copies[a].ramification_spheres.back()) {
                        batch_scratch_grid.insert(sph);
                    }
                }
                for (size_t a = 0; a < attempt_copies.size(); ++a) {
                    size_t k = attempt_k[a];
                    if (!grew[a]) {
                        if (++tries[k] >= max_primary_tries) {
                            cout << "Failed to grow glial cell" << endl;
                            done[k] = true;
                        }
                        continue;
                    }
                    // check_collision_with_branches=false: exclude ANY of
                    // this cell's own spheres (any branch), not just this
                    // exact branch. growPrimaryBranch's own internal
                    // AddOneSphere calls already correctly handled
                    // branch-vs-branch collision for this cell (including
                    // deliberately allowing a new branch's first several
                    // spheres to sit close to sibling branches near the
                    // soma) -- reconciliation here only needs to catch
                    // genuinely new risk from parallel growth: collisions
                    // against a *different* cell's branch grown blind to
                    // this one in the same round.
                    bool collided = false;
                    for (const auto &sph : attempt_copies[a].ramification_spheres.back()) {
                        if (!sphere_grid.canSpherebePlaced(sph, false) || !batch_scratch_grid.canSpherebePlaced(sph, false)) {
                            collided = true;
                            break;
                        }
                    }
                    if (collided) {
                        if (++tries[k] >= max_primary_tries) {
                            cout << "Failed to grow glial cell" << endl;
                            done[k] = true;
                        }
                        continue;
                    }
                    size_t i = batch_start + k;
                    pop[i] = attempt_copies[a];
                    nbr_spheres[i] = attempt_nbr_spheres[a];
                    for (const auto &sph : pop[i].ramification_spheres.back()) {
                        sphere_grid.insert(sph);
                    }
                    done[k] = true;
                }

                any_pending = false;
                for (size_t k = 0; k < batch_size; ++k) {
                    if (!done[k]) { any_pending = true; break; }
                }
            }
        }
    }

    for (size_t i = 0; i < pop.size(); ++i) {
        pop[i].compute_processes_icvf(spheres_overlap_factor, min_limits, max_limits);
    }

    ICVF(axons, glial_pop1, glial_pop2, blood_vessels); // updates glial_pop1_processes_icvf and glial_pop2_processes_icvf
    current_icvf = (population_nbr == 1) ? glial_pop1_processes_icvf
                                         : glial_pop2_processes_icvf;
    display_progress(current_icvf, target_icvf);
    if (current_icvf >= target_icvf) return;

    // Iterative growth
    int nbr_tries = 0;
    const int   max_tries   = 100000;
    double      prev_icvf   = current_icvf;
    int         stall_count = 0;
    const int   stall_limit = 100;
    const double rel_eps    = 1e-6;

    while (current_icvf < target_icvf && nbr_tries <= max_tries) {
        for (size_t i = 0; i < growths.size(); ++i) {
            bool grew;
            if (growths[i].glial_cell_to_grow.allow_branching) {
                grew = growths[i].growSecondaryBranch(nbr_spheres[i], mean_len, std_len, spheres_overlap_factor);
            } else {
                grew = growths[i].growPrimaryBranch(nbr_spheres[i], mean_len, std_len, spheres_overlap_factor);
            }
            pop[i] = growths[i].glial_cell_to_grow;
            if (grew) {
                // Register the newly grown branch immediately, so other cells'
                // branches (grown later in this same loop, or in later
                // iterations) can see it instead of growing blind to it.
                for (const auto &sph : pop[i].ramification_spheres.back()) {
                    sphere_grid.insert(sph);
                }
            }
        }

        if (nbr_tries % 10 == 0) {
            for (auto& g : pop) {
                g.compute_processes_icvf(spheres_overlap_factor, min_limits, max_limits);
            }
            ICVF(axons, glial_pop1, glial_pop2, blood_vessels);
            current_icvf = (population_nbr == 1) ? glial_pop1_processes_icvf
                                                 : glial_pop2_processes_icvf;
            display_progress(current_icvf, target_icvf);

            if (current_icvf >= target_icvf) break;

            double denom = std::max(1.0, std::abs(prev_icvf));
            if (std::abs(current_icvf - prev_icvf) / denom < rel_eps) {
                if (++stall_count > stall_limit) {
                    std::cerr << "Stuck in a loop, stopping growth!\n";
                    break;
                }
            } else {
                stall_count = 0;
            }
            prev_icvf = current_icvf;
        }
        ++nbr_tries;
    }

    if (nbr_tries > max_tries) {
        std::cerr << "Max attempts reached while growing branches!\n";
    }

    // No final grid sync needed here: the soma was added to the grid when the
    // cell was placed, and each branch was added as soon as it was grown, above.
}

void CaterpillarGrowth::PlaceGlialCells() {
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

    // fully outside the box
    auto fully_outside_box = [&](const Sphere& s) {
        for (int i = 0; i < 3; ++i) {
            if (s.center[i] + s.radius < min_limits[i] || s.center[i] - s.radius > max_limits[i]) {
                return true;
            }
        }
        return false;
    };

    auto soma_volume = [](double r) {
        return (4.0 * M_PI * r * r * r) / 3.0;
    };

    // Expand the sampling window to allow placements outside the bounds.
    // If you want a different expansion per population, split this.
    expanded_for_glial_space = std::max(0.0, glial_pop1_radius_mean*3);

    // --- Population 1 (astrocytes) ------------------------------------------
    std::normal_distribution<> dis_radius_pop1(glial_pop1_radius_mean, glial_pop1_radius_std);

    double glial_icvf_1 = 0.0;
    const int MAX_GLOBAL_FAILS_1 = 20000;
    int global_fail_count_1 = 0;

    std::uniform_real_distribution<> base_x1(min_limits[0] - expanded_for_glial_space, max_limits[0] + expanded_for_glial_space);
    std::uniform_real_distribution<> base_y1(min_limits[1] - expanded_for_glial_space, max_limits[1] + expanded_for_glial_space);
    std::uniform_real_distribution<> base_z1(min_limits[2] - expanded_for_glial_space, max_limits[2] + expanded_for_glial_space);

    while (glial_icvf_1 < target_glial_pop1_soma_icvf && global_fail_count_1 < MAX_GLOBAL_FAILS_1) {
        Sphere s;
        bool placed = false;
        bool fully_inside = false;
        bool fully_outside = false;

        for (int attempt = 0; attempt < 10000; ++attempt) {
            const double rad = draw_positive_radius(dis_radius_pop1);

            s = Sphere(
                /*object_id*/ 0,                        // keep your ID scheme if needed
                /*index*/ static_cast<int>(glial_pop1.size()),
                /*type*/ glial_cell_constant,
                { base_x1(gen), base_y1(gen), base_z1(gen) },
                rad
            );

            if (sphere_grid.canSpherebePlaced(s)) {
                fully_inside = (signed_box_margin(s) >= 0.0);
                fully_outside = fully_outside_box(s);
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
        glial_cell.addToGrid(sphere_grid);
        glial_pop1.push_back(glial_cell);

        // Count ICVF only if fully inside
        if (fully_inside) {
            glial_icvf_1 += soma_volume(glial_cell.soma.radius) / total_volume;
        }
        else if (!fully_outside) {
            glial_icvf_1 += glial_cell.soma.sphereBoxIntersectionVolume(min_limits, max_limits, /*eps_rel=*/1e-6) / total_volume;
        }
    }

    auto transform = [&](double x, double r) {
            if (x > expanded_for_glial_space){
                x = max_limits[0] + (x - (expanded_for_glial_space-r));
            }
            else{
                x = min_limits[0] - ((expanded_for_glial_space-r) - x);
            }
            return x;
        };

    if (target_glial_pop1_soma_icvf > 0.0) {

        Sphere s;

        // add some extra in extra space
        for (int i = 0; i < 10; i++) {

            const double rad = draw_positive_radius(dis_radius_pop1);

            std::uniform_real_distribution<> base_z1(0, 2*(expanded_for_glial_space-rad));

            double x = base_x1(gen);
            double y = base_y1(gen);
            double z = base_z1(gen);
            z = transform(z, rad);

            s = Sphere(
                        /*object_id*/ 0,                        // keep your ID scheme if needed
                        /*index*/ static_cast<int>(glial_pop1.size()),
                        /*type*/ glial_cell_constant,
                        {x,y,z},
                        rad
                    );
            if (sphere_grid.canSpherebePlaced(s)) {
                Glial glial_cell = Glial(s.object_id, s, glial_pop1_branching);
                glial_cell.addToGrid(sphere_grid);
                glial_pop1.push_back(glial_cell);
            }
        }
    }
    

    glial_pop1_soma_icvf = glial_icvf_1;

    // --- Population 2 (oligodendrocytes) ------------------------------------
    std::normal_distribution<> dis_radius_pop2(glial_pop2_radius_mean, glial_pop2_radius_std);

    double glial_icvf_2 = 0.0;
    const int MAX_GLOBAL_FAILS_2 = 20000;
    int global_fail_count_2 = 0;

    std::uniform_real_distribution<> base_x2(min_limits[0] - expanded_for_glial_space, max_limits[0] + expanded_for_glial_space);
    std::uniform_real_distribution<> base_y2(min_limits[1] - expanded_for_glial_space, max_limits[1] + expanded_for_glial_space);
    std::uniform_real_distribution<> base_z2(min_limits[2] - expanded_for_glial_space, max_limits[2] + expanded_for_glial_space);

    while (glial_icvf_2 < target_glial_pop2_soma_icvf && global_fail_count_2 < MAX_GLOBAL_FAILS_2) {
        Sphere s;
        bool placed = false;
        bool fully_inside = false;
        bool fully_outside = false;

        for (int attempt = 0; attempt < 10000; ++attempt) {
            const double rad = draw_positive_radius(dis_radius_pop2);

            s = Sphere(
                /*object_id*/ 0,
                /*index*/ static_cast<int>(glial_pop2.size()),
                /*type*/ glial_cell_constant,
                { base_x2(gen), base_y2(gen), base_z2(gen) },
                rad
            );

            if (sphere_grid.canSpherebePlaced(s)) {
                fully_inside = (signed_box_margin(s) >= 0.0);
                fully_outside = fully_outside_box(s);
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
        glial_cell.addToGrid(sphere_grid);
        glial_pop2.push_back(glial_cell);

        // Count ICVF only if fully inside
        if (fully_inside) {
            glial_icvf_2 += soma_volume(glial_cell.soma.radius) / total_volume;
        }
        else if (!fully_outside) {
            glial_icvf_2 += glial_cell.soma.sphereBoxIntersectionVolume(min_limits, max_limits, /*eps_rel=*/1e-6) / total_volume;
        }
    }

    glial_pop2_soma_icvf = glial_icvf_2;

    if (target_glial_pop2_soma_icvf > 0.0) {

        Sphere s;
        // add some extra in extra space
        for (int i = 0; i < 10; i++) {

            const double rad = draw_positive_radius(dis_radius_pop2);

            std::uniform_real_distribution<> base_z2(0, 2*(expanded_for_glial_space-rad));

            double x = base_x2(gen);
            double y = base_y2(gen);
            double z = base_z2(gen);
            z = transform(z, rad);

            s = Sphere(
                        /*object_id*/ 0,                        // keep your ID scheme if needed
                        /*index*/ static_cast<int>(glial_pop2.size()),
                        /*type*/ glial_cell_constant,
                        {x,y,z},
                        rad
                    );
            if (sphere_grid.canSpherebePlaced(s)) {
                Glial glial_cell = Glial(s.object_id, s, glial_pop2_branching);
                glial_cell.addToGrid(sphere_grid);
                glial_pop2.push_back(glial_cell);
            }
        }
    }
}


bool CaterpillarGrowth::checkNoCollisions()
{
    // Rebuild from the current, authoritative state rather than trusting whatever
    // incremental grid updates happened during growth/regrow/swelling.
    sphere_grid.clear();
    for (auto &axon : axons) {
        axon.addToGrid(sphere_grid);
    }
    for (auto &bv : blood_vessels) {
        bv.addToGrid(sphere_grid);
    }
    for (auto &g : glial_pop1) {
        g.addToGrid(sphere_grid);
    }
    for (auto &g : glial_pop2) {
        g.addToGrid(sphere_grid);
    }

    int nbr_collisions = 0;

    // check_collision_with_branches=false everywhere below: this pass checks each sphere
    // against *other objects*, not against its own object's other branches. A branch's
    // spheres near its point of origin are supposed to touch the trunk/soma/parent branch
    // they emerged from (that's how a branch attaches) -- true self-intersection between
    // a cell's own branches is already prevented during growth (see collideswithItself).

    for (const auto &axon : axons) {
        for (const auto &sph : axon.outer_spheres) {
            if (!sphere_grid.canSpherebePlaced(sph, false)) {
                cout << "Axon : " << axon.id << ", sphere : " << sph.id << " collides with environment !" << endl;
                nbr_collisions++;
            }
        }
    }
    for (const auto &bv : blood_vessels) {
        for (const auto &branch : bv.ramification_spheres) {
            for (const auto &sph : branch) {
                if (!sphere_grid.canSpherebePlaced(sph, false)) {
                    cout << "Blood vessel : " << bv.id << ", sphere : " << sph.id << " (branch_id=" << sph.branch_id << ") collides with environment !" << endl;
                    nbr_collisions++;
                }
            }
        }
    }
    for (const auto &g : glial_pop1) {
        if (!sphere_grid.canSpherebePlaced(g.soma, false)) {
            cout << "Glial pop1 : " << g.id << ", soma collides with environment !" << endl;
            nbr_collisions++;
        }
        for (const auto &branch : g.ramification_spheres) {
            for (const auto &sph : branch) {
                if (!sphere_grid.canSpherebePlaced(sph, false)) {
                    cout << "Glial pop1 : " << g.id << ", sphere : " << sph.id << " collides with environment !" << endl;
                    nbr_collisions++;
                }
            }
        }
    }
    for (const auto &g : glial_pop2) {
        if (!sphere_grid.canSpherebePlaced(g.soma, false)) {
            cout << "Glial pop2 : " << g.id << ", soma collides with environment !" << endl;
            nbr_collisions++;
        }
        for (const auto &branch : g.ramification_spheres) {
            for (const auto &sph : branch) {
                if (!sphere_grid.canSpherebePlaced(sph, false)) {
                    cout << "Glial pop2 : " << g.id << ", sphere : " << sph.id << " collides with environment !" << endl;
                    nbr_collisions++;
                }
            }
        }
    }

    if (nbr_collisions > 0) {
        cout << "checkNoCollisions: " << nbr_collisions << " colliding sphere(s) found." << endl;
        return false;
    }
    return true;
}

bool CaterpillarGrowth::SanityCheck(std::vector<Axon>& growing_axons,
                                        std::vector<double>& stuck_radii_,
                                        std::vector<int>& stuck_indices_,
                                        const std::vector<std::size_t> &layer_start_spheres) {

    if (growing_axons.size() <= 1) {
        return true;  // No axons to check against each other
    }

    bool any_new_spheres = false;
    for (size_t i = 0; i < growing_axons.size(); ++i) {
        if (growing_axons[i].outer_spheres.size() > layer_start_spheres[i]) {
            any_new_spheres = true;
            break;
        }
    }
    if (!any_new_spheres) {
        return true;  // No new growth this layer to check
    }

    // growing_axons are private copies grown in parallel by processBatchWithThreadPool:
    // none of this layer's new spheres are in sphere_grid yet, so batch-mates growing
    // in the same layer are invisible to each other during growth. Reuse the persistent
    // batch_scratch_grid (cleared, not reconstructed -- this runs every sub-batch/round)
    // over just this layer's new segments -- everything before layer_start_spheres[i]
    // was already checked and committed in an earlier layer, so re-inserting/re-checking
    // it here would be pure waste.
    batch_scratch_grid.clear();
    for (size_t i = 0; i < growing_axons.size(); ++i) {
        auto &spheres = growing_axons[i].outer_spheres;
        for (size_t k = layer_start_spheres[i]; k < spheres.size(); ++k) {
            batch_scratch_grid.insert(spheres[k]);
        }
    }

    std::unordered_set<int> collided_ids;
    for (size_t i = 0; i < growing_axons.size(); ++i) {
        const auto &axon = growing_axons[i];
        const auto &spheres = axon.outer_spheres;
        for (size_t k = layer_start_spheres[i]; k < spheres.size(); ++k) {
            if (!sphere_grid.canSpherebePlaced(spheres[k]) || !batch_scratch_grid.canSpherebePlaced(spheres[k])) {
                collided_ids.insert(axon.id);
                break;
            }
        }
    }

    bool all_clear = collided_ids.empty();

    for (size_t i = 0; i < growing_axons.size(); ++i) {
        Axon &axon = growing_axons[i];
        if (collided_ids.count(axon.id)) {
            stuck_radii_.push_back(axon.radius);
            stuck_indices_.push_back(axon.id);
            // Roll back only this layer's new (colliding) segment, not
            // previously-committed layers already merged into sphere_grid.
            axon.truncate_to(layer_start_spheres[i]);
            axon.update_Volume(spheres_overlap_factor, min_limits, max_limits);
        }
    }

    return all_clear;
}



bool CaterpillarGrowth::hasReachedTrueWall(const Axon& ax) const {
    if (ax.outer_spheres.empty()) return false;
    const Sphere &last = ax.outer_spheres.back();
    // Mirrors AddOneSphere's own "reached the wall" check (grow_axons.cpp),
    // but always against the true max_limits, never a layer-capped
    // extended_max_limits, so the caller can tell "reached this layer's cap
    // only" apart from "actually reached the box's real wall."
    return last.center[ax.growth_axis] + last.radius > max_limits[ax.growth_axis];
}

void CaterpillarGrowth::growAxon(Axon& axon_to_grow, int &index, double& stuck_radius, int& stuck_index, double layer_depth_cap, std::size_t layer_start_spheres, bool can_shrink_this_round) {

    // Cap only this axon's own growth_axis component of extended_max_limits
    // at this layer's depth (or its true wall, whichever is nearer); true
    // min_limits/max_limits (used for all other border checks) are untouched.
    Eigen::Vector3d capped_max_limits = max_limits;
    capped_max_limits[axon_to_grow.growth_axis] = std::min(max_limits[axon_to_grow.growth_axis], layer_depth_cap);

    AxonGrowth growth(axon_to_grow, &sphere_grid, min_limits, capped_max_limits, min_limits, max_limits, epsilon, min_radius);

    growth.growthThread(stuck_radius, stuck_index, spheres_overlap_factor, can_shrink_this_round, layer_start_spheres);

}

void CaterpillarGrowth::processBatchWithThreadPool(
    std::vector<Axon>& axons_to_grow,
    std::vector<int>& indices,
    std::vector<double>& stuck_radii,
    std::vector<int>& stuck_indices,
    double layer_depth_cap,
    const std::vector<std::size_t> &layer_start_spheres,
    bool can_shrink_this_round)
{
    ThreadPool pool(nbr_threads);

    std::vector<std::future<void>> futures; // Store futures for synchronization

    for (size_t i = 0; i < axons_to_grow.size(); ++i) {
        futures.emplace_back(pool.enqueueTask(
            [this, &axons_to_grow, &stuck_radii, &stuck_indices, &indices, &layer_start_spheres, layer_depth_cap, can_shrink_this_round, i]() {
                growAxon(axons_to_grow[i], indices[i], stuck_radii[i], stuck_indices[i], layer_depth_cap, layer_start_spheres[i], can_shrink_this_round);
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

void CaterpillarGrowth::growAxonsLayered(std::vector<double> &radii_, std::vector<bool> &has_myelin, std::vector<double> &angles)
{
    stuck_radii.clear();
    stuck_indices.clear();

    if (axons.empty()) {
        return;
    }

    // Layer thickness scales with the already-computed grid_voxel_size
    // (itself sized off the active populations' characteristic radius)
    // rather than a fixed absolute constant, which would be meaningless
    // across configs at very different scales. Swept empirically across
    // 1/2/4/6/10/16x grid_voxel_size on the 70% packing benchmark: 1x won on
    // every axis simultaneously (0 abandoned, highest final ICVF, fastest
    // runtime) -- thinner layers mean more, smaller depth-fairness
    // checkpoints, so no axon races far ahead of its neighbors before they
    // get a turn, which also means fewer sub-batch collisions and fewer
    // costly retry rounds overall. TEMP: still overridable via
    // LAYER_THICKNESS_FACTOR for further A/B testing.
    double layer_thickness_factor = 1.0;
    if (const char *env_val = std::getenv("LAYER_THICKNESS_FACTOR")) {
        layer_thickness_factor = std::atof(env_val);
    }
    const double layer_thickness = std::max(1.0, layer_thickness_factor * grid_voxel_size);

    // How many push-only (no shrink) rounds a sub-batch gets per layer
    // before falling back to a shrink-enabled round -- see the sub-batch
    // retry loop below. TEMP: overridable via MAX_NO_SHRINK_ROUNDS for A/B
    // testing; defaults to 3.
    int max_no_shrink_rounds = 3;
    if (const char *env_val = std::getenv("MAX_NO_SHRINK_ROUNDS")) {
        max_no_shrink_rounds = std::atoi(env_val);
    }

    // active[k] = index into axons (and, in parallel, radii_/has_myelin/angles,
    // which stay index-aligned with axons throughout seeding).
    std::vector<int> active(axons.size());
    std::iota(active.begin(), active.end(), 0);

    int layer = 0;
    const int max_layers = 100000; // defensive backstop against runaway iteration
    // All boxes in this codebase are cubes with min_limits = {0,0,0}, so one
    // absolute depth cap works regardless of which axis a given axon's own
    // growth_axis happens to be.
    const double box_depth = max_limits[2] - min_limits[2];

    while (!active.empty() && layer < max_layers) {
        std::cout << "---   Depth layer " << layer << " (" << active.size() << " active axons)   --- " << endl;
        double layer_depth_cap = (layer + 1) * layer_thickness;
        display_progress(std::min(layer_depth_cap, box_depth), box_depth);
        if (on_growth_progress) {
            on_growth_progress(std::min(layer_depth_cap, box_depth), box_depth);
        }

        // Within this layer, still process active axons in nbr_threads-sized
        // sub-batches, sequentially -- not all of them in one ThreadPool pass.
        // Growing hundreds of axons at once as mutually-blind private copies
        // makes SanityCheck's cross-axon collisions far more frequent than
        // before (many more axons invisible to each other simultaneously),
        // and a SanityCheck collision gets zero retries (unlike a
        // growth_attempts exhaustion, which gets 10) -- so at that scale it
        // was permanently killing a large fraction of axons every layer.
        // Sub-batching restores the original nbr_threads-scale blindness
        // (each sub-batch's new spheres are committed to sphere_grid before
        // the next sub-batch in the SAME layer starts), while still bounding
        // how far ahead any axon can get to just this one layer's thickness.
        std::vector<int> next_active;
        next_active.reserve(active.size());

        for (size_t batch_start = 0; batch_start < active.size(); batch_start += static_cast<size_t>(nbr_threads)) {
            size_t batch_end = std::min(active.size(), batch_start + static_cast<size_t>(nbr_threads));

            std::vector<int> layer_indices;
            std::vector<std::size_t> layer_start_spheres;
            std::vector<Eigen::Vector3d> layer_start_end;
            std::vector<int> layer_start_grow_straight;
            std::vector<int> layer_start_straight_growths;
            layer_indices.reserve(batch_end - batch_start);
            layer_start_spheres.reserve(batch_end - batch_start);
            layer_start_end.reserve(batch_end - batch_start);
            layer_start_grow_straight.reserve(batch_end - batch_start);
            layer_start_straight_growths.reserve(batch_end - batch_start);
            for (size_t b = batch_start; b < batch_end; ++b) {
                int idx = active[b];
                layer_start_spheres.push_back(axons[idx].outer_spheres.size());
                layer_indices.push_back(axons[idx].id);
                // Captured here (axons[idx] itself is never mutated during
                // the retry rounds below) so a failed round's jostle
                // (which nudges .end) or straight/random state can be
                // reset back to this layer's true starting point on the
                // next retry, without needing a full re-copy from
                // axons[idx] -- see the round > 0 branch below.
                layer_start_end.push_back(axons[idx].end);
                layer_start_grow_straight.push_back(axons[idx].grow_straight);
                layer_start_straight_growths.push_back(axons[idx].straight_growths);
            }

            // Retry the whole sub-batch (push-based placement only, no
            // shrinking -- see AddOneSphere) for up to max_no_shrink_rounds
            // rounds before allowing any axon to shrink at all. Each round
            // re-grows every axon in the sub-batch fresh from the
            // already-committed state (undoing any partial progress from a
            // previous failed round), and SanityCheck reconciles collisions
            // among this sub-batch's mutually-blind private copies exactly
            // as before. Only once that budget is exhausted does one final
            // round enable shrinking (and, built into AddOneSphere, its
            // "remember the least-shrunk position seen so far" fallback) --
            // giving the cheap, non-destructive push-only strategy several
            // honest chances first, since a sub-batch-mate that only
            // recently collided may free up room again once others finish
            // moving.
            std::vector<Axon> layer_axons;
            std::vector<double> stuck_radii_;
            std::vector<int> stuck_indices_;

            for (int round = 0; ; ++round) {
                bool can_shrink_this_round = axon_can_shrink && (round >= max_no_shrink_rounds);

                if (round == 0) {
                    layer_axons.clear();
                    layer_axons.reserve(batch_end - batch_start);
                    for (size_t b = batch_start; b < batch_end; ++b) {
                        int idx = active[b];
                        axons[idx].growth_attempts = 0; // reset-per-round retry budget
                        layer_axons.push_back(axons[idx]); // private copy, grown invisibly to its sub-batch-mates until reconciled below
                    }
                } else {
                    // axons[idx] itself is never touched during these retry
                    // rounds (only merged back after the loop exits), so
                    // re-copying its entire growth history -- every sphere
                    // from every previously-committed layer, hundreds to
                    // 1000+ for a deep axon -- would just reproduce
                    // byte-for-byte what's already sitting in layer_axons[i]
                    // up to layer_start_spheres[i]. Truncating the existing
                    // private copy back to that point (a plain vector
                    // resize-down, no reallocation) undoes the failed
                    // round's partial growth for free instead. .end and the
                    // straight/random state must also be reset explicitly,
                    // not just outer_spheres -- a failed round's jostling
                    // (growthThread) mutates .end, and grow_straight/
                    // straight_growths persist across calls by design (see
                    // AxonGrowth::growthThread), so without this a failed
                    // round's leftover jostle/straight-state would carry
                    // into the next retry instead of starting fresh, unlike
                    // a real re-copy from the untouched axons[idx].
                    for (size_t i = 0; i < layer_axons.size(); ++i) {
                        layer_axons[i].truncate_to(layer_start_spheres[i]);
                        layer_axons[i].growth_attempts = 0;
                        layer_axons[i].end = layer_start_end[i];
                        layer_axons[i].grow_straight = layer_start_grow_straight[i];
                        layer_axons[i].straight_growths = layer_start_straight_growths[i];
                    }
                }

                std::vector<double> raw_stuck_radii(layer_axons.size(), -1);
                std::vector<int> raw_stuck_indices(layer_axons.size(), -1);
                processBatchWithThreadPool(layer_axons, layer_indices, raw_stuck_radii, raw_stuck_indices, layer_depth_cap, layer_start_spheres, can_shrink_this_round);

                stuck_radii_.clear();
                stuck_indices_.clear();
                for (size_t i = 0; i < layer_axons.size(); ++i) {
                    if (raw_stuck_radii[i] > 0) {
                        stuck_radii_.push_back(raw_stuck_radii[i]);
                        stuck_indices_.push_back(raw_stuck_indices[i]);
                    }
                }

                // layer_axons grew as private copies, invisible to each other
                // during the ThreadPool pass; reconcile collisions among this
                // sub-batch's new segments (and against the already-committed
                // environment) here.
                SanityCheck(layer_axons, stuck_radii_, stuck_indices_, layer_start_spheres);

                if (stuck_indices_.empty() || can_shrink_this_round) {
                    // Either the whole sub-batch succeeded cleanly this
                    // round, or this was already the final (shrink-enabled)
                    // round -- nothing further to try either way.
                    break;
                }
            }

            // Merge this sub-batch's results back into the main axons list
            // and sphere_grid immediately (so the NEXT sub-batch in this same
            // layer sees them), inserting only the spheres actually added
            // this layer, and decide which axons are still active for the
            // next layer.
            for (size_t i = 0; i < layer_axons.size(); ++i) {
                int idx = active[batch_start + i];
                Axon &grown = layer_axons[i];

                for (size_t k = layer_start_spheres[i]; k < grown.outer_spheres.size(); ++k) {
                    sphere_grid.insert(grown.outer_spheres[k]);
                }
                axons[idx] = std::move(grown);

                bool is_stuck = std::find(stuck_indices_.begin(), stuck_indices_.end(), axons[idx].id) != stuck_indices_.end();
                if (is_stuck) {
                    // Stuck within this layer (retries exhausted, or a
                    // collision SanityCheck couldn't resolve): keep whatever
                    // was committed through the previous layer (already
                    // truncated back to that point by growthThread/
                    // SanityCheck), and exclude from future layers. Not
                    // relocated/retried at a different position --
                    // seedAllAxons already chose this axon's position
                    // carefully, so there's no reason to expect an
                    // arbitrary new spot would do better.
                    stuck_radii.push_back(axons[idx].radius);
                    stuck_indices.push_back(axons[idx].id);
                    continue;
                }

                if (hasReachedTrueWall(axons[idx])) {
                    continue; // done, exclude from future layers
                }

                next_active.push_back(idx); // still short of its true wall, keep going next layer
            }
        }

        active = std::move(next_active);
        ++layer;
    }

    if (layer >= max_layers) {
        std::cerr << "Warning: depth-layered axon growth hit the max layer count safety backstop\n";
    }

    // Discard any axon that got permanently stuck (retries and jostling
    // exhausted without ever reaching the true wall -- tracked in
    // stuck_indices as it happens, above), or that ended up too trivial
    // overall (fewer than 10 spheres total) even if it did reach the wall.
    // An axon is never kept at a partial depth: it either makes it all the
    // way, or it is discarded entirely, matching this function's documented
    // intent (see the "No relocate-and-retry" comment at its call site) --
    // both cases need every already-committed sphere removed from
    // sphere_grid, not just the first, since layered growth can leave
    // several committed spheres behind before an axon is deemed not worth
    // keeping.
    std::unordered_set<int> stuck_id_set(stuck_indices.begin(), stuck_indices.end());
    int counts = 0;
    for (int i = static_cast<int>(axons.size()) - 1; i >= 0; --i) {
        bool too_short = axons[i].outer_spheres.size() < 10;
        bool got_stuck = stuck_id_set.count(axons[i].id) > 0;
        if (too_short || got_stuck) {
            for (const auto &sph : axons[i].outer_spheres) {
                sphere_grid.remove(sph);
            }
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





// Axon growth
double CaterpillarGrowth::radiusVariation(const Axon &axon)
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



double volumeFrustumCone(double r1, double r2, double h)
{
    return M_PI * h * (r1 * r1 + r2 * r2 + r1 * r2) / 3;
}


void CaterpillarGrowth ::create_SWC_file(std::ostream &out)
{
    std::vector<Axon> final_axons;

    final_axons = axons;


    out << "cell_type cell_id component component_id X Y Z inner_radius outer_radius" << endl;
    std::sort(final_axons.begin(), final_axons.end(), [](const Axon a, Axon b) -> bool
              { return a.radius > b.radius; }); // sort by size

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
            int cell_id = final_axons[i].outer_spheres[j].object_id;
            int component_id = 0;
            std::string cell_type = "axon";
            std::string component = "axon";
            double x = final_axons[i].outer_spheres[j].center[0];
            double y = final_axons[i].outer_spheres[j].center[1];
            double z = final_axons[i].outer_spheres[j].center[2];
            double outer_radius = final_axons[i].outer_spheres[j].radius;
            
            if (final_axons[i].inner_spheres.size() > 0)
            {
                double inner_radius = final_axons[i].inner_spheres[j].radius;
                out << cell_type << " " <<  cell_id << " " << component << " " << component_id << " " << x << " " << y << " " << z << " " << inner_radius << " " << outer_radius << endl;
            }
            else
            {
                out << cell_type << " " <<  cell_id << " " << component << " " << component_id << " " << x << " " << y << " " << z << " " << outer_radius << " " << outer_radius << endl;
            }
            
        }


    }
    cosPhiSquared = cos_2_all/final_axons.size();


    for (auto &glial_cell : glial_pop1)
    {

        int cell_id = glial_cell.id;
        int component_id = 0;
        std::string cell_type = "glial_cell";
        std::string component = "soma";
        double x = glial_cell.soma.center[0];
        double y = glial_cell.soma.center[1];
        double z = glial_cell.soma.center[2];
        double outer_radius = glial_cell.soma.radius;

        out << cell_type << " " <<  cell_id << " " << component << " " << component_id << " " << x << " " << y << " " << z << " " << outer_radius << " " << outer_radius << endl;
        
        for (long unsigned int j = 0; j < glial_cell.ramification_spheres.size(); j++)
        {
            cell_id = glial_cell.id;
            for (long unsigned int k = 0; k < glial_cell.ramification_spheres[j].size(); k++)
            {
                component_id = glial_cell.ramification_spheres[j][k].branch_id;
                x = glial_cell.ramification_spheres[j][k].center[0];
                y = glial_cell.ramification_spheres[j][k].center[1];
                z = glial_cell.ramification_spheres[j][k].center[2];
                outer_radius = glial_cell.ramification_spheres[j][k].radius;
                component = "branch";
                out << cell_type << " " <<  cell_id << " " << component << " " << component_id << " " << x << " " << y << " " << z << " " << outer_radius << " " << outer_radius << endl;
            }
        }
    }

    for (auto &glial_cell : glial_pop2)
    {
        int cell_id = glial_cell.id;
        int component_id = 0;
        std::string cell_type = "glial_cell";
        std::string component = "soma";
        double x = glial_cell.soma.center[0];
        double y = glial_cell.soma.center[1];
        double z = glial_cell.soma.center[2];
        double outer_radius = glial_cell.soma.radius;

        out << cell_type << " " <<  cell_id << " " << component << " " << component_id << " " << x << " " << y << " " << z << " " << outer_radius << " " << outer_radius << endl;
        
        for (long unsigned int j = 0; j < glial_cell.ramification_spheres.size(); j++)
        {
            cell_id = glial_cell.id;
            for (long unsigned int k = 0; k < glial_cell.ramification_spheres[j].size(); k++)
            {
                component_id = glial_cell.ramification_spheres[j][k].branch_id;
                x = glial_cell.ramification_spheres[j][k].center[0];
                y = glial_cell.ramification_spheres[j][k].center[1];
                z = glial_cell.ramification_spheres[j][k].center[2];
                outer_radius = glial_cell.ramification_spheres[j][k].radius;
                component = "branch";
                out << cell_type << " " <<  cell_id << " " << component << " " << component_id << " " << x << " " << y << " " << z << " " << outer_radius << " " << outer_radius << endl;
            }
        }
    }

    for (auto &bv : blood_vessels)
    {
        for (auto &branch : bv.ramification_spheres) {
        for (auto &s : branch) {
            int cell_id = bv.id;
            int component_id = s.branch_id;
            std::string cell_type = "blood_vessel";
            std::string component = (component_id == 0) ? "blood_vessel" : "branch";
            double x = s.center[0];
            double y = s.center[1];
            double z = s.center[2];
            double outer_radius = s.radius;

            out << cell_type << " " <<  cell_id << " " << component << " " << component_id << " " << x << " " << y << " " << z << " " << outer_radius << " " << outer_radius << endl;

        }
        }

    }
}

void CaterpillarGrowth::simulation_file(std::ostream &out, const std::chrono::seconds &duration)
{
    // Ensure booleans print as "true"/"false" instead of "1"/"0"
    out << std::boolalpha;

    // --- Original / Calculated Metrics ---
    out << "Duration " << duration.count() << std::endl;
    out << "Num_axons " << axons.size() << std::endl;
    out << "Voxel Size " << max_limits[0] << std::endl;
    out << "Axon icvf " << axons_w_myelin_icvf + axons_wo_myelin_icvf << std::endl;
    out << "Axon without myelin icvf " <<  axons_wo_myelin_icvf << std::endl;
    out << "Axon with myelin icvf " << axons_w_myelin_icvf << std::endl;
    out << "Myelin icvf " << myelin_icvf << std::endl;
    out << "Glial cell population 1 icvf soma " << glial_pop1_soma_icvf << std::endl;
    out << "Glial cell population 1 icvf branches " << glial_pop1_processes_icvf << std::endl;
    out << "Glial cell population 2 icvf soma " << glial_pop2_soma_icvf << std::endl;
    out << "Glial cell population 2 icvf branches " << glial_pop2_processes_icvf << std::endl;
    out << "Blood vessel icvf " << blood_vessels_icvf << std::endl;
    out << "Total icvf " << axons_w_myelin_icvf + axons_wo_myelin_icvf + glial_pop1_soma_icvf + glial_pop1_processes_icvf + glial_pop2_soma_icvf + glial_pop2_processes_icvf + blood_vessels_icvf << std::endl;
    
    // --- General Simulation Parameters ---
    out << "nbr_axons_populations " << nbr_axons_populations << std::endl;
    out << "crossing_fibers_type " << crossing_fibers_type << std::endl;
    out << "Number of threads " << nbr_threads << std::endl;

    // --- Growth & Morphology Parameters ---
    out << "alpha " << alpha << std::endl;
    out << "beta " << beta << std::endl;
    out << "regrowth threshold " << regrow_thr << std::endl;
    out << "minimum sphere radius " << min_radius << std::endl;
    out << "epsilon (Tortuosity) " << epsilon << std::endl;
    out << "undulation_factor " << undulation_factor << std::endl;
    out << "beading_amplitude " << beading_amplitude << std::endl;
    out << "beading_std " << beading_std << std::endl;
    out << "swelling_factor " << swelling_factor << std::endl;
    out << "axons can shrink " << axon_can_shrink << std::endl;
    out << "spheres overlap factor " << spheres_overlap_factor << std::endl;

    // --- Blood Vessels ---
    out << "epsilon_blood_vessels " << epsilon_blood_vessels << std::endl;
    out << "mean_vessel_rad " << mean_vessel_rad << std::endl;
    out << "std_vessel_rad " << std_vessel_rad << std::endl;

    // --- Constants ---
    out << "cosPhiSquared " << cosPhiSquared << std::endl;
    out << "c1 " << c1 << std::endl;
    out << "c2 " << c2 << std::endl;
    out << "c3 " << c3 << std::endl;

    // --- Glial Population 1 Parameters ---
    out << "glial_pop1_nbr_primary_processes " << glial_pop1_nbr_primary_processes << std::endl;
    out << "glial_pop1_branching " << glial_pop1_branching << std::endl;
    out << "mean_glial_pop1_process_length " << mean_glial_pop1_process_length << std::endl;
    out << "std_glial_pop1_process_length " << std_glial_pop1_process_length << std::endl;
    out << "glial_pop1_radius_mean " << glial_pop1_radius_mean << std::endl;
    out << "glial_pop1_radius_std " << glial_pop1_radius_std << std::endl;

    // --- Glial Population 2 Parameters ---
    out << "glial_pop2_nbr_primary_processes " << glial_pop2_nbr_primary_processes << std::endl;
    out << "glial_pop2_branching " << glial_pop2_branching << std::endl;
    out << "mean_glial_pop2_process_length " << mean_glial_pop2_process_length << std::endl;
    out << "std_glial_pop2_process_length " << std_glial_pop2_process_length << std::endl;
    out << "glial_pop2_radius_mean " << glial_pop2_radius_mean << std::endl;
    out << "glial_pop2_radius_std " << glial_pop2_radius_std << std::endl;
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


double CaterpillarGrowth::originalFunction(const double &x, const double &outerRadius) {
    return outerRadius - (myelin_thickness(x) + x);
}

double CaterpillarGrowth::derivative(const double &x) {
    if (x <= 0.0001) return -1.0; // Avoid division by zero
    return -(c2 * 2.0 + c3 / (2.0 * x) + 1);
}

double CaterpillarGrowth::findInnerRadius(const double &outerRadius) {
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

void CaterpillarGrowth::add_Myelin()
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
        nbr_spheres_to_delete = nbr_spheres_to_delete*spheres_overlap_factor;
        double prob_ranvier_node = 600*spheres_overlap_factor/axons[index].radius;

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

            inner_sphere = Sphere(axons[index].outer_spheres[i].id, inner_axon_constant, axons[index].outer_spheres[i].object_id, axons[index].outer_spheres[i].center, innerRadius);
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
