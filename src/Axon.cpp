#include "Axon.h"
#include "constants.h"
#include "Eigen/Dense"
#include <iostream>


using namespace Eigen;
using namespace std;


Axon::Axon()
{}

Axon::~Axon()
{
    inner_spheres.clear();
    outer_spheres.clear();
}


void Axon::keep_one_sphere(){
    Sphere s (this->outer_spheres[0]);
    outer_spheres.clear();
    add_sphere(s);
    growth_attempts += 1;


}

void Axon::truncate_to(std::size_t n){
    if (outer_spheres.size() > n) {
        outer_spheres.resize(n);
    }
    // growth_attempts is deliberately left untouched here -- callers decide
    // whether/how to adjust it (e.g. growthThread's retry path increments it
    // itself, matching keep_one_sphere()'s existing convention).
}

void Axon::destroy(){
    outer_spheres.clear();
    growth_attempts = 0;

}

void Axon::addToGrid(SphereGrid &grid) const{

    for (const auto &sph : outer_spheres){
        grid.insert(sph);
    }
    for (const auto &sph : inner_spheres){
        grid.insert(sph);
    }
}


void Axon::add_sphere(const Sphere &sphere_to_add){

    // add sphere to list of spheres
    this->outer_spheres.push_back(sphere_to_add);
}

 

void Axon::removeEngulfedSpheres(std::vector<Sphere> &removed){

    if (outer_spheres.size() < 2) {
        return;
    }

    std::vector<Sphere> kept;
    kept.reserve(outer_spheres.size());

    for (const Sphere &candidate : outer_spheres) {
        bool discard_candidate = false;

        // Collapse backward: pop any kept sphere the candidate now fully
        // engulfs, or -- if the candidate itself is fully engulfed by the
        // current last kept sphere -- drop the candidate and stop.
        while (!kept.empty()) {
            const Sphere &top = kept.back();
            double dist = (candidate.center - top.center).norm();

            if (dist + candidate.radius <= top.radius) {
                // candidate entirely inside top: useless, don't keep it
                discard_candidate = true;
                break;
            }
            if (dist + top.radius <= candidate.radius) {
                // top entirely inside candidate: top was useless all along
                removed.push_back(top);
                kept.pop_back();
                continue;
            }
            break; // neither engulfs the other
        }

        if (discard_candidate) {
            removed.push_back(candidate);
        } else {
            kept.push_back(candidate);
        }
    }

    outer_spheres = std::move(kept);
}

void Axon::update_Volume(const int &factor, const Eigen::Vector3d &min_limits, const Eigen::Vector3d &max_limits){
    double new_volume = 0.0;

    for (size_t i = 1; i < outer_spheres.size(); ++i) {
        const Sphere &last_sphere = outer_spheres[i - 1];
        const Sphere &current_sphere = outer_spheres[i];

        double distance = (current_sphere.center - last_sphere.center).norm();

        // Volume of truncated cone
        double segment_volume = M_PI * distance * (
            current_sphere.radius * current_sphere.radius +
            last_sphere.radius * last_sphere.radius +
            current_sphere.radius * last_sphere.radius) / 3.0;

        bool current_in_bounds = Obstacle::check_borders(min_limits, max_limits, current_sphere.center, barrier_tickness);
        bool last_in_bounds = Obstacle::check_borders(min_limits, max_limits, last_sphere.center, barrier_tickness);

        if (current_in_bounds && last_in_bounds) {
            new_volume += segment_volume;
        } else if (current_in_bounds || last_in_bounds) {
            new_volume += segment_volume / 2.0;
        }
        // If both are out of bounds, add nothing
    }

    volume = new_volume;
}


/*
double sphereVolume(double radius) {
    return (4.0 / 3.0) * M_PI * std::pow(radius, 3);
}

double overlapVolume(double r1, double r2, double d) {
    if (d >= r1 + r2) return 0.0; // No overlap

    double part1 = (r1 + r2 - d) * (r1 + r2 - d);
    double part2 = d * d + 2 * d * (r1 + r2) - 3 * (r1 - r2) * (r1 - r2);
    return (M_PI * part1 * part2) / (12.0 * d);
}

void Axon::update_Volume(const int &factor, const Eigen::Vector3d &min_limits, const Eigen::Vector3d &max_limits) {
    double new_volume = 0.0;

    if (outer_spheres.empty())
        return;

    for (size_t i = 0; i < outer_spheres.size(); ++i) {
        const Sphere &current_sphere = outer_spheres[i];

        bool current_in_bounds = check_borders(min_limits, max_limits, current_sphere.center, barrier_tickness);

        if (!current_in_bounds) {
            continue; // Skip this sphere if it's out of bounds
        }
        
        // Always add volume of current sphere
        new_volume += sphereVolume(current_sphere.radius);

        // If not the first, subtract overlapping volume with previous sphere (if both in bounds)
        if (i > 0) {
            const Sphere &prev_sphere = outer_spheres[i - 1];

            bool prev_in_bounds = check_borders(min_limits, max_limits, prev_sphere.center, barrier_tickness);

            if (current_in_bounds && prev_in_bounds) {
                double d = (current_sphere.center - prev_sphere.center).norm();
                double ov = overlapVolume(current_sphere.radius, prev_sphere.radius, d);
                new_volume -= ov;
            }
        }
    }

    volume = new_volume;
}
*/