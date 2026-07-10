#include "Glial.h"
#include "Eigen/Dense"
#include "constants.h"
#include <iostream>

using namespace Eigen;
using namespace std;

Glial::Glial()
{}

Glial::~Glial()
{}



void Glial::addToGrid(SphereGrid &grid) const{

    grid.insert(soma);
    for (const auto &branch : ramification_spheres){
        for (const auto &sph : branch){
            grid.insert(sph);
        }
    }
}


void Glial::compute_processes_icvf(const int &factor, const Eigen::Vector3d &min_limits, const Eigen::Vector3d &max_limits) {
    volume_processes = 0.0;
    for (int b = 0; b < ramification_spheres.size(); ++b) {
        for (int i = factor; i < ramification_spheres[b].size(); i += factor) {
            if ((soma.center-ramification_spheres[b][i-factor].center).norm() < soma.radius) {
                continue;
            }
            Sphere &s = ramification_spheres[b][i-factor];
            Sphere &s_next = ramification_spheres[b][i];
            if (s_next.center [0] + s_next.radius < min_limits[0] || s_next.center [0] - s_next.radius > max_limits[0] ||
                s_next.center [1] + s_next.radius < min_limits[1] || s_next.center [1] - s_next.radius > max_limits[1] ||
                s_next.center [2] + s_next.radius < min_limits[2] || s_next.center [2] - s_next.radius > max_limits[2]) {
                continue;
            }
            
            double distance = (s_next.center - s.center).norm();
            double v = M_PI * (s.radius * s.radius + s_next.radius * s_next.radius + s.radius * s_next.radius) * distance / 3.0;
            volume_processes += v;

        }
    }
}