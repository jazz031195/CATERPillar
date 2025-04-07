#include "axongammadistribution.h"
#include "grow_cells.h"
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


CellGrowth::~CellGrowth() {}

CellGrowth::CellGrowth(const CellGrowth &other)
  : axons(other.axons),
    glial_cells(other.glial_cells),
    min_limits(other.min_limits),
    max_limits(other.max_limits),
    extended_min_limits(other.extended_min_limits),
    extended_max_limits(other.extended_max_limits),
    std_dev(other.std_dev),
    min_radius(other.min_radius),
    grow_straight(other.grow_straight),
    finished(other.finished)
{}


// Function to check if a point is inside a dilated box
bool CellGrowth::check_borders(const Eigen::Vector3d&  min_l, const Eigen::Vector3d&  max_l, const Eigen::Vector3d& pos, const double& distance_to_border) {

    
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



double CellGrowth::clamp(double value, double lower, double upper) {
    return std::max(lower, std::min(value, upper));
}

Eigen::Vector3d CellGrowth::generate_random_point_on_sphere(double std) {

    // Random number generation
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dist(0.0, std); // Normal distribution with mean=0 and std deviation=std

    // Generate random theta and phi values for spherical coordinates
    double theta = dist(gen);  // Random theta value
    double cos_phi = dist(gen);  // Random cos_phi value
    double phi = acos(cos_phi);  // Convert cos_phi to phi

    // Ensure theta is within the range [0, 2*pi] and phi is within the range [0, pi]
    theta = fmod(theta, 2 * M_PI);
    phi = clamp(phi, 0.0, M_PI/2);  // Ensure phi is within [0, pi]

    // Spherical to Cartesian conversion
    double x = sin(phi) * cos(theta);
    double y = sin(phi) * sin(theta);
    double z = cos(phi);

    return Eigen::Vector3d(x, y, z);

}

// Function to compute the rotation matrix from vector1 to vector2
Eigen::Matrix3d CellGrowth::rotation_matrix_from_vectors(const Eigen::Vector3d &vec1, const Eigen::Vector3d &vec2) {
    Eigen::Vector3d a = vec1.normalized();
    Eigen::Vector3d b = vec2.normalized();
    
    Eigen::Vector3d v = a.cross(b);
    double c = a.dot(b);
    double s = v.norm();
    
    Eigen::Matrix3d kmat;
    kmat <<  0, -v[2], v[1],
             v[2], 0, -v[0],
            -v[1], v[0], 0;
    
    Eigen::Matrix3d rotation_matrix = Eigen::Matrix3d::Identity() + kmat + kmat * kmat * ((1 - c) / (s * s + 1e-10));
    
    return rotation_matrix;
}

// Function to apply bias toward the target for a single point
Eigen::Vector3d CellGrowth::apply_bias_toward_target(const Eigen::Vector3d &point, const Eigen::Vector3d &target) {
    // Step 1: Compute the rotation matrix from [1, 0, 0] to the target
    Eigen::Vector3d reference(1, 0, 0);  // Reference vector [1, 0, 0]
    Eigen::Matrix3d R = rotation_matrix_from_vectors(reference, target);

    // Step 2: Rotate the point using the computed rotation matrix
    Eigen::Vector3d rotated_point = R * point;

    return rotated_point;
}





bool CellGrowth::canSpherebePlaced(const Sphere &sph){

    for (auto &axon : axons)
    {
        if (!(axon.id == sph.object_id && sph.object_type == 0))
        {
            // Check overlap
            if (axon.isSphereInsideAxon(sph)) 
            {

                return false;
                
            }
        }
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



