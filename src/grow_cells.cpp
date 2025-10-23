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
    glial_pop1(other.glial_pop1),
    glial_pop2(other.glial_pop2),
    min_limits(other.min_limits),
    max_limits(other.max_limits),
    extended_min_limits(other.extended_min_limits),
    extended_max_limits(other.extended_max_limits),
    std_dev(other.std_dev),
    min_radius(other.min_radius),
    grow_straight(other.grow_straight),
    finished(other.finished)
{}


const std::vector<Axon>& CellGrowth::AX() const { return *axons; }

const std::vector<Axon>* CellGrowth::AXptr() const { return axons; }

void CellGrowth::update_environment(const std::vector<Axon>* axons_,
                                const std::vector<Glial>* glial_pop1_,
                                const std::vector<Glial>* glial_pop2_) noexcept {
    axons      = axons_;
    glial_pop1 = glial_pop1_;
    glial_pop2 = glial_pop2_;
}

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
    // azimuth: uniform (no azimuthal bias)
    std::uniform_real_distribution<double> U(0.0, 2.0*M_PI);
    double theta = U(gen);

    // polar: bias toward the pole via cos(phi) ~ N(1, std)
    std::normal_distribution<double> N(1.0, std);
    double cphi = clamp(N(gen), -1.0, 1.0);  // cos(phi)
    double sphi = std::sqrt(std::max(0.0, 1.0 - cphi*cphi));

    double x = sphi * std::cos(theta);
    double y = sphi * std::sin(theta);
    double z = cphi;                // near 1 when std is small
    return Eigen::Vector3d(x, y, z);  // unit by construction
}


// Function to compute the rotation matrix from vector1 to vector2
Eigen::Matrix3d CellGrowth::rotation_matrix_from_vectors(const Eigen::Vector3d &vec1, const Eigen::Vector3d &vec2) {
    Eigen::Vector3d a = vec1.normalized();
    Eigen::Vector3d b = vec2.normalized();
    
    double c = a.dot(b);

    if (c > 1.0 - 1e-12) return Eigen::Matrix3d::Identity();       // already aligned
    if (c < -1.0 + 1e-12) {                                         // 180°
        Eigen::Vector3d axis = a.unitOrthogonal();
        return Eigen::AngleAxisd(M_PI, axis).toRotationMatrix();
    }

    Eigen::Vector3d v = a.cross(b);
    double s = v.norm();
    Eigen::Matrix3d K;
    K <<   0,   -v.z(),  v.y(),
         v.z(),     0,  -v.x(),
        -v.y(),  v.x(),    0;
    return Eigen::Matrix3d::Identity() + K + K*K * ((1 - c)/(s*s));
}

// Function to apply bias toward the target for a single point
Eigen::Vector3d CellGrowth::apply_bias_toward_target(const Eigen::Vector3d &point, const Eigen::Vector3d &target) {
    // Step 1: Compute the rotation matrix from [0, 0, 0] to the target
    Eigen::Vector3d reference(0, 0, 1);  // Reference vector [0, 0, 1]
    Eigen::Matrix3d R = rotation_matrix_from_vectors(reference, target);

    // Step 2: Rotate the point using the computed rotation matrix
    Eigen::Vector3d rotated_point = (R * point).normalized();

    return rotated_point;
}


bool CellGrowth::checkAxonsOverlap(Sphere &sph){

    for (const auto& axon : *axons) {
        if (axon.isSphereInsideAxon(sph)) {
            return false;  
        }
    }
    return true; 
    
}


bool CellGrowth::canSpherebePlaced(Sphere &sph){


    bool axons_check = checkAxonsOverlap(sph);

    if (!axons_check){
        return false;
    }


    // check collision other glial cells 
    for (const Glial& glial : *glial_pop1) {

        if (glial.collides_with_GlialCell(sph, barrier_tickness)) {

            return false;
        }

    }

    for (const Glial& glial : *glial_pop2) {
        if (glial.collides_with_GlialCell(sph, barrier_tickness)) {
            return false;
        }
    }


    return true;
} 
