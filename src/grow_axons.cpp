#include "axongammadistribution.h"
#include "grow_axons.h"
#include <algorithm> // std::sort
#include <random>
#include <chrono>
#include <future>

using namespace std;
using namespace Eigen;
using namespace std::chrono;

Growth::Growth(Axon& axon_to_grow_, std::vector<Axon> axons_, Eigen::Vector3d voxel_size_, bool tortuous_, double max_radius_, int grow_straight_)
{

    axon_to_grow = axon_to_grow_;
    voxel_size = voxel_size_;
    tortuous = tortuous_;
    finished = false;
    max_radius = max_radius_;
    grow_straight = grow_straight_;
    axons = axons_;

    if (!tortuous)
    {
        grow_straight = 0;
    }
    initialise_env_axons();

}

bool areBoxesCloserThanDistance(const std::vector<Eigen::Vector2d> & box1, const std::vector<Eigen::Vector2d> & box2, double distanceThreshold) {
    for (int i = 0; i < 3; ++i) { // Loop over dimensions (x, y, z)
        double minPosBox1 = box1[i][0];
        double maxPosBox1 = box1[i][1];
        double minPosBox2 = box2[i][0];
        double maxPosBox2 = box2[i][1];
        
        // Calculate the distance between the closest points of the two intervals
        double distance = std::max(0.0, std::max(minPosBox1 - maxPosBox2, minPosBox2 - maxPosBox1));
        
        // Check if the distance is smaller than the specified threshold
        if (distance >= distanceThreshold) {
            return false; // Boxes are not closer in this dimension
        }
    }
    return true; // Boxes are closer in all dimensions
}

void Growth::initialise_env_axons(){

    for (unsigned i = 0; i < axons.size(); i++){
        // if id is in axon_ids
        if (axons[i].id != axon_to_grow.id && (axons[i].Box, axon_to_grow.Box, 5)){
            env_axons.push_back(axons[i]);
        }
    }

    axons.clear();
}

bool Growth::check_borders(Eigen::Vector3d pos, double distance_to_border)
{

    Eigen::Vector2d new_min_limits = {distance_to_border, distance_to_border};
    Eigen::Vector2d new_max_limits = {voxel_size[0], voxel_size[1]};


    if (pos[0] < new_min_limits[0] -distance_to_border )
    {
        return false;
    }
    if (pos[1] < new_min_limits[1]-distance_to_border)
    {
        return false;
    }
    if (pos[0] > new_max_limits[0]+distance_to_border)
    {
        return false;
    }
    if (pos[1] > new_max_limits[1] + distance_to_border)
    {
        return false;
    }

    return true;

}


std::tuple<double, double> phi_theta_to_target(Eigen::Vector3d new_pos, Eigen::Vector3d end)
{
    /*
    Find phi and theta angles in spherical coordinates
    see : http://douglashalse.com/index.php/2019/09/25/spherical-coordinates/

    */

    Eigen::Vector3d vector_to_target;
    if (end[2] < new_pos[2])
    {
        vector_to_target = {end[0] - new_pos[0], end[1] - new_pos[1], end[2] - 0.01 - new_pos[2]};
    }
    else
    {
        vector_to_target = {end[0] - new_pos[0], end[1] - new_pos[1], end[2] + 0.01 - new_pos[2]};
    }
    vector_to_target = vector_to_target.normalized();

    double phi_to_target;
    double theta_to_target;

    // theta is angle between (1,0) and (x,y)
    // varies between 0 and 2pi
    theta_to_target = atan2(vector_to_target[1], vector_to_target[0]);
    if (theta_to_target < 0)
    {
        theta_to_target += 2 * M_PI;
    }

    // phi is angle between (0,0,1) and (x,y,z)
    // varies between 0 and pi
    if (vector_to_target[2] == 0)
    {
        phi_to_target = M_PI / 2;
    }
    else if (vector_to_target == Eigen::Vector3d({0, 0, -1}))
    {
        phi_to_target = M_PI;
    }
    else if (vector_to_target == Eigen::Vector3d({0, 0, 1}))
    {
        phi_to_target = 0;
    }
    else
    {
        // varies between -pi/2 and pi/2
        phi_to_target = atan((sqrt(vector_to_target[0] * vector_to_target[0] + vector_to_target[1] * vector_to_target[1])) / vector_to_target[2]);
        if (phi_to_target < 0)
        {
            phi_to_target += M_PI;
        }
    }

    return std::make_tuple(phi_to_target, theta_to_target);
}

/*
bool Growth::isSphereColliding(Sphere sph)
{

    Vector3d position = sph.center;

    for (unsigned i = 0; i < env_axons.size(); i++)
    {
        if (sph.ax_id != env_axons[i].id)
        {

            double distance_to_be_inside = sph.radius;
            bool isinside = env_axons[i].isPosInsideAxon(position, distance_to_be_inside, max_radius);

            if (isinside)
            {
                return true;
                break;
            }
        }
    }
    return false;
}
*/

/*
bool Growth::isSphereColliding(Sphere sph) {

    int num_threads;
    if (env_axons.size() < thread_capacity){
        num_threads = env_axons.size(); // Set the desired number of threads
    }
    else{
        num_threads = thread_capacity;
    }
    std::vector<std::future<bool>> futures;

    for (unsigned i = 0; i < env_axons.size(); ++i) {
        futures.emplace_back(std::async(std::launch::async, [&](unsigned index) {
            return isSphereCollidingwithAxon(env_axons[index], sph);
        }, i));
    }

    for (auto& future : futures) {
        if (future.get()) {
            return true;
        }
    }

    return false;
}
*/


bool Growth::isSphereColliding(Sphere sph){

    for (unsigned j = 0; j < env_axons.size(); j ++) {

        if (isSphereCollidingwithAxon(env_axons[j], sph)){
            return true;
        }
    }
    return false;
           
}


bool Growth::isSphereCollidingwithAxon(Axon ax, Sphere sph){

    return ax.isSphereInsideAxon(sph);
}



void Growth::find_next_center(Sphere &s,  double dist_)
{

    // phi and theta angles in spherical coordinates
    // to update position with
    double phi, theta;
    // phi and theta angles in spherical coordinates
    // that links the current position to the target position
    double phi_to_target, theta_to_target;
    // phi and theta angles in spherical coordinates
    // that links previous position with current one
    double prev_phi, prev_theta;
    // number of tries
    int tries = 0;
    // random
    std::random_device rd;
    std::mt19937 gen(rd());

    // differences in position from intial to final position
    double delta_x;
    double delta_y;
    double delta_z;

    // position of new sphere to add
    Eigen::Vector3d new_pos;
    // current position (position of last sphere)
    Eigen::Vector3d curr_pos = {axon_to_grow.spheres[axon_to_grow.spheres.size() - 1].center[0], axon_to_grow.spheres[axon_to_grow.spheres.size() - 1].center[1], axon_to_grow.spheres[axon_to_grow.spheres.size() - 1].center[2]};
    // current position (position of 4 spheres before last )
    Eigen::Vector3d prev_pos;
    if (axon_to_grow.spheres.size() > 2)
    {
        prev_pos = {axon_to_grow.spheres[axon_to_grow.spheres.size() - 2].center[0], axon_to_grow.spheres[axon_to_grow.spheres.size() - 2].center[1], axon_to_grow.spheres[axon_to_grow.spheres.size() - 2].center[2]};
    }
    else
    {
        prev_pos = {axon_to_grow.spheres[axon_to_grow.spheres.size() - 1].center[0], axon_to_grow.spheres[axon_to_grow.spheres.size() - 1].center[1], axon_to_grow.spheres[axon_to_grow.spheres.size() - 1].center[2] - dist_};
    }

    tie(phi_to_target, theta_to_target) = phi_theta_to_target(curr_pos, axon_to_grow.end + Eigen::Vector3d{0, 0, 10});

    // if not tortuous, phi and theta are those that lead to target
    if (!tortuous)
    {
        phi = phi_to_target;
        theta = theta_to_target;
    }
    // if tortuous, phi and theta comme from distribution
    else
    {
        bool inside_voxel = false;
        int restart = 0;

        std::normal_distribution<float> phi_dist(phi_to_target / M_PI, 0.12);
        std::normal_distribution<float> theta_dist(theta_to_target / M_PI, 0.12);
            

        phi = phi_dist(gen) * M_PI;
        theta = theta_dist(gen) * M_PI;

    }
    // spherical coordinates to cartesian
    delta_x = dist_ * cos(theta) * sin(phi);
    delta_y = dist_ * sin(theta) * sin(phi);
    delta_z = dist_ * cos(phi);

    // update new_pos
    new_pos = curr_pos;
    new_pos[0] += delta_x;
    new_pos[1] += delta_y;
    new_pos[2] += delta_z;
    // sphere to return
    s.center = new_pos;

}



void Growth::find_next_center_straight(double distance, Sphere &s)
{
    Eigen::Vector3d last_center = axon_to_grow.spheres[axon_to_grow.spheres.size() - 1].center;
    Eigen::Vector3d before_last_center = axon_to_grow.spheres[axon_to_grow.spheres.size() - 2].center;
    Eigen::Vector3d straight_vector = (last_center - before_last_center).normalized();
    straight_vector *= distance;
    Eigen::Vector3d new_center = last_center + straight_vector;
    s.center = new_center;
}


bool Growth::GrowAxon(double radius_, bool create_sphere)
{

    // if sphere collides with environment
    bool collides;
    bool can_grow_;
    // the distance between two consecutive spheres is the maximum radius / 2
    double max_radius_ = max(radius_, axon_to_grow.spheres[axon_to_grow.spheres.size()-1].radius);

    double distance = (max_radius_)/2;

    if (axon_to_grow.spheres[axon_to_grow.spheres.size() -1].center[2] < voxel_size[2]) // still growing
    {

            collides = true;
            can_grow_ = false;
            int tries = 0;
            Sphere s(axon_to_grow.spheres.size(), axon_to_grow.id, axon_to_grow.begin, radius_);
            while (!can_grow_ and tries < 1000)
            {
                // find the center of next sphere by taking a random position

                if (tortuous)
                {
                    if (grow_straight==1)
                    {
                        find_next_center_straight(distance, s);
                    }
                    else
                    {
                        find_next_center(s, distance);
                    }
                    // check if there is a collision
                    collides = isSphereColliding(s);
                    // check if the sphere is inside the voxel
                    bool isinside = check_borders(s.center, s.radius);

                    // can grow if there is no collision and if sphere is inside voxel
                    if (!collides && isinside){
                        can_grow_ = true;
                    }
                    else{
                        // if it was straight, try not straight
                        if (grow_straight==1){
                            find_next_center(s, distance);
                            // check if there is a collision
                            collides = isSphereColliding(s);
                            if (radius_ != s.radius){
                                assert(0);
                            }
                            // check if the sphere is inside the voxel
                            bool isinside = check_borders(s.center, s.radius);
                            if (!collides && isinside){
                                can_grow_ = true;
                            }
                            else{
                                can_grow_ = false;
                            }
                        }
                        else{
                            can_grow_ = false;
                        }
                    }

                    tries += 1;
                }
                else
                {
                    find_next_center(s, distance);

                    collides = isSphereColliding(s);
                    // check if the sphere is inside the voxel
                    bool isinside = check_borders(s.center, s.radius);
                    if (!collides && isinside){
                        can_grow_ = true;
                    }
                    else{
                        can_grow_ = false;
                    }
                }
            }

            if (can_grow_ )
            {
                if (create_sphere){
                    axon_to_grow.add_sphere(s); // we want to change the radius of the sphere
                }
                //cout << "sphere added to axon, new size :" << axon_to_grow.spheres.size() << endl;
                // if we reach edge of voxel
                if (axon_to_grow.spheres[axon_to_grow.spheres.size() -1].center[2] > voxel_size[2])
                {
                    //cout << "Axon " << axon_to_grow.id << " done!!" << endl;
                    //cout << axon_to_grow.spheres[axon_to_grow.spheres.size() -1].center[2]  << endl;
                    //cout << voxel_size[2]  << endl;
                    finished = true;
                    
                }
                return true;
                
            }
            else // collides and >1000 tries
            {
                return false;
            }


        
    }
    else // axon is fully grown
    {
        //cout << "Axon " << axon_to_grow.id << " done!" << endl;
        finished = true;
        return true;
    }
}
