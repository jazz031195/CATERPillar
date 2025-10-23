#include "sphere.h"
#include "constants.h"
#include <iostream>
#include <random>

using namespace Eigen;
using namespace std;


double Sphere::minDistance(const Eigen::Vector3d &O) const{

    Vector3d m = O - this->center;
    // minimum distance to the cylinder axis.
    double distance_to_sphere = m.norm();

    //Minimum distance to the shpere wall.
    double d_ = (distance_to_sphere - this->radius);
   // return d_>0.0?d_:0.0;
    return d_;
}



bool Sphere::isInside(Eigen::Vector3d pos, double distance_to_be_inside){
    double d_ = (pos - this->center).norm();
 
    d_ = d_-this->radius;
    
   // return d_>0.0?d_:0.0;
    return d_ <= distance_to_be_inside;
}

bool Sphere::CollideswithSphere(const Sphere &other, const double &tolerance) const {
    double dist = (other.center - center).norm();
    return dist <= (other.radius + radius + tolerance);
}

void Sphere::getPointOnSphereSurface(Eigen::Vector3d &point, Eigen::Vector3d &vector, const Eigen::Vector3d &vector2, const bool& primary_process) const {
    // random device
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<> disTheta(0, 2 * M_PI);  // Uniform distribution for azimuthal angle
    std::uniform_real_distribution<> disCosPhi(-1, 1);       // Uniform distribution for cosine of polar angle

    bool valid = false;
    int tries = 0;
    while (!valid && tries < 1000) {
        // Generate random spherical coordinates
        double theta = disTheta(gen);     // Azimuthal angle (0 to 2 * pi)
        double cosPhi = disCosPhi(gen);   // Cosine of polar angle (-1 to 1)
        double sinPhi = sqrt(1 - cosPhi * cosPhi); // Sin of the polar angle

        // Convert spherical coordinates to Cartesian coordinates (uniformly on sphere)
        vector = Eigen::Vector3d(sinPhi * cos(theta), sinPhi * sin(theta), cosPhi);
        vector.normalize();  // Ensure vector is unit length
        point = this->center + (this->radius * vector);  // Scale by the radius to get point on surface

        if (!primary_process) {
            double angle = acos(vector.dot(vector2));
            // Check if the angle is <= pi / 2
            if (angle <= M_PI / 2.0) {
                valid = true;
            }
            else{
                tries++;
                valid = false;
            }
        }
        else{
            valid = true;
        }
    }
}




Eigen::Vector3d Sphere::getNormalOnSpherePoint(Eigen::Vector3d point){
    return (point - this->center).normalized();
}



// Fast helpers
inline double sqr(double x) { return x * x; }

// Squared distance from point p to AABB (0 if inside)
inline double minDist2PointAABB(const Eigen::Vector3d& p,
                                const Eigen::Vector3d& bmin,
                                const Eigen::Vector3d& bmax)
{
    double d2 = 0.0;
    for (int i = 0; i < 3; ++i) {
        if (p[i] < bmin[i]) { double d = bmin[i] - p[i]; d2 += d*d; }
        else if (p[i] > bmax[i]) { double d = p[i] - bmax[i]; d2 += d*d; }
    }
    return d2;
}

// Squared farthest distance from point p to any point in AABB
inline double maxDist2PointAABB(const Eigen::Vector3d& p,
                                const Eigen::Vector3d& bmin,
                                const Eigen::Vector3d& bmax)
{
    double d2 = 0.0;
    for (int i = 0; i < 3; ++i) {
        const double dmin = p[i] - bmin[i];
        const double dmax = p[i] - bmax[i];
        const double d = (std::abs(dmin) > std::abs(dmax)) ? dmin : dmax;
        d2 += d * d;
    }
    return d2;
}

// Estimate fraction of cell inside sphere by sampling 8 corners + center (cheap and decent).
inline double Sphere::estimateFraction9(const Eigen::Vector3d& bmin,
                                const Eigen::Vector3d& bmax,
                                const double R2)
{
    int inside = 0;
    // 8 corners
    for (int mask = 0; mask < 8; ++mask) {
        Eigen::Vector3d q{
            (mask & 1) ? bmax.x() : bmin.x(),
            (mask & 2) ? bmax.y() : bmin.y(),
            (mask & 4) ? bmax.z() : bmin.z()
        };
        const Eigen::Vector3d d = q - center;
        if (d.squaredNorm() <= R2) ++inside;
    }
    // center
    Eigen::Vector3d c = 0.5 * (bmin + bmax);
    if ( (c - center).squaredNorm() <= R2 ) ++inside;

    return inside / 9.0; // in [0,1]
}

// Main function.
// eps_rel is relative to the box volume (default 1e-6).
// max_depth caps refinement (safety); 32 is plenty for doubles.
double Sphere::sphereBoxIntersectionVolume(const Eigen::Vector3d& min_limits,
                                   const Eigen::Vector3d& max_limits,
                                   double eps_rel,
                                   int    max_depth)
{
    const Eigen::Vector3d& C = center;
    const double R  = radius;
    const double R2 = R * R;

    // Early rejects / accepts on the whole box.
    const double boxVol =
        (max_limits.x() - min_limits.x()) *
        (max_limits.y() - min_limits.y()) *
        (max_limits.z() - min_limits.z());

    if (boxVol <= 0.0) return 0.0;

    // No intersection if closest distance from center to box > R
    if (minDist2PointAABB(C, min_limits, max_limits) >= R2) return 0.0;

    // Sphere fully inside the box?
    auto inside_box = [&](int axis) -> bool {
        return (C[axis] - R >= min_limits[axis]) && (C[axis] + R <= max_limits[axis]);
    };
    if (inside_box(0) && inside_box(1) && inside_box(2)) {
        return (4.0 / 3.0) * M_PI * R2 * R;
    }

    // Box fully inside the sphere?
    if (maxDist2PointAABB(C, min_limits, max_limits) <= R2) {
        return boxVol;
    }

    // Adaptive octree (iterative)
    struct Cell {
        Eigen::Vector3d bmin, bmax;
        int depth;
    };
    std::vector<Cell> stack;
    stack.reserve(1 << 12);
    stack.push_back({min_limits, max_limits, 0});

    const double epsAbs = std::max(1e-18, eps_rel) * boxVol;
    double volume = 0.0;

    while (!stack.empty()) {
        Cell cell = stack.back();
        stack.pop_back();

        const Eigen::Vector3d& a = cell.bmin;
        const Eigen::Vector3d& b = cell.bmax;

        // Quick tests
        const double dmin2 = minDist2PointAABB(C, a, b);
        if (dmin2 >= R2) continue; // entirely outside

        const double dmax2 = maxDist2PointAABB(C, a, b);
        const double vol = (b.x() - a.x()) * (b.y() - a.y()) * (b.z() - a.z());

        if (dmax2 <= R2) {         // entirely inside
            volume += vol;
            continue;
        }

        // Stop if tiny cell or depth limit reached -> estimate fraction and accumulate
        const bool tiny = vol <= epsAbs || cell.depth >= max_depth;
        if (tiny) {
            const double f = estimateFraction9(a, b, R2);
            volume += f * vol;
            continue;
        }

        // Subdivide into 8 octants
        Eigen::Vector3d mid = 0.5 * (a + b);
        // Precompute children mins/maxs
        Eigen::Vector3d xs[2] = { {a.x(), a.y(), a.z()}, {mid.x(), mid.y(), mid.z()} };
        Eigen::Vector3d xe[2] = { {mid.x(), mid.y(), mid.z()}, {b.x(),  b.y(),  b.z()} };

        // Push 8 children (order doesn't matter). Unrolled for speed.
        // child index bits: xbit(1), ybit(2), zbit(4)
        stack.push_back({ {a.x(), a.y(), a.z()}, {mid.x(), mid.y(), mid.z()}, cell.depth+1 });
        stack.push_back({ {mid.x(), a.y(), a.z()}, {b.x(),   mid.y(), mid.z()}, cell.depth+1 });
        stack.push_back({ {a.x(), mid.y(), a.z()}, {mid.x(), b.y(),   mid.z()}, cell.depth+1 });
        stack.push_back({ {mid.x(), mid.y(), a.z()}, {b.x(),   b.y(),   mid.z()}, cell.depth+1 });
        stack.push_back({ {a.x(), a.y(), mid.z()}, {mid.x(), mid.y(), b.z()}, cell.depth+1 });
        stack.push_back({ {mid.x(), a.y(), mid.z()}, {b.x(),   mid.y(), b.z()}, cell.depth+1 });
        stack.push_back({ {a.x(), mid.y(), mid.z()}, {mid.x(), b.y(),   b.z()}, cell.depth+1 });
        stack.push_back({ {mid.x(), mid.y(), mid.z()}, {b.x(),   b.y(),   b.z()}, cell.depth+1 });
    }

    return volume;
}