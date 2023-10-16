#include "axongammadistribution.h"
#include "grow_axons.h"
#include <algorithm> // std::sort
#include <random>
#include <chrono>
#include "swipeprune.h"

using namespace std;
using namespace Eigen;
using namespace std::chrono;

Growth::Growth(Axon& axon_to_grow_, std::vector<Axon> axons_, std::vector<Axon> axons_to_regrow_, Eigen::Vector3d voxel_size_, bool tortuous_, double max_radius_, int grow_straight_)
{

    axon_to_grow = axon_to_grow_;
    voxel_size = voxel_size_;
    tortuous = tortuous_;
    finished = false;
    max_radius = max_radius_;
    grow_straight = grow_straight_;
    axons_to_regrow = axons_to_regrow_;
    axons = axons_;

    if (!tortuous)
    {
        grow_straight = 0;
    }
    initialise_env_axons();

}

void Growth::initialise_env_axons(){

    for (unsigned i = 0; i < axons.size(); i++){
        // if id is in axon_ids
        auto it = std::find(axon_to_grow.nearby_axons.begin(), axon_to_grow.nearby_axons.end(), axons[i].id);
        if (it != axon_to_grow.nearby_axons.end()){
            env_axons.push_back(axons[i]);
        }
    }

    for (unsigned i = 0; i < axons_to_regrow.size(); i++){
        auto it = std::find(axon_to_grow.nearby_axons.begin(), axon_to_grow.nearby_axons.end(), axons_to_regrow[i].id);
        if (it != axon_to_grow.nearby_axons.end()){
            env_axons.push_back(axons_to_regrow[i]);
        }
    }

    //cout << "env_axons :" << env_axons.size() << endl; 
    axons.clear();
    axons_to_regrow.clear();

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

/*
bool Growth::check_borders(Eigen::Vector3d pos, double distance_to_border, Eigen::Vector2d &twin_delta_pos)
{

    Eigen::Vector2d new_min_limits = {distance_to_border, distance_to_border};
    Eigen::Vector2d new_max_limits = {voxel_size[0] - distance_to_border, voxel_size[1] - distance_to_border};

    twin_delta_pos = {0.0, 0.0};

    if ((pos[0] - new_min_limits[0]) < 0)
    {
        // min plane of x
        twin_delta_pos[0] = voxel_size[0];
    }
    if ((pos[1] - new_min_limits[1]) < 0)
    {
        // min plane of y
        twin_delta_pos[1] = voxel_size[1];
    }
    if ((pos[0] - new_max_limits[0]) > 0)
    {
        // max plane of x
        twin_delta_pos[0] = -voxel_size[0];
    }
    if ((pos[1] - new_max_limits[1]) > 0)
    {
        // max plane of y
        twin_delta_pos[1] = -voxel_size[1];
    }
    if (twin_delta_pos != Eigen::Vector2d({0.0, 0.0}))
    {
        return true;
    }
    else
    {
        return false;
    }
}
*/
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
bool Growth::isSphereColliding(Dynamic_Sphere sph)
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

bool Growth::isSphereColliding_(Dynamic_Sphere sph){
    // for all axons
    for (unsigned i = 0; i < env_axons.size(); i++)
    {
        int target_id = env_axons[i].id;

        if (sph.ax_id != target_id){

            if (env_axons[i].isSphereInsideAxon_(sph)){
                return true;
            }
            
        }
    }
    return false;
}


bool Growth::isSphereColliding_long(Dynamic_Sphere sph){
    // for all axons
    for (unsigned i = 0; i < env_axons.size(); i++)
    {
        int target_id = env_axons[i].id;

        if (sph.ax_id != target_id){

            if (env_axons[i].isSphereInsideAxon_long(sph)){
                return true;
            }
            
        }
    }
    return false;
}

void Growth::find_next_center(Dynamic_Sphere &s,  double dist_)
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

bool Growth::GrowFirstSphere()
{

    // if sphere collides with environment
    bool collides;
    // first sphere to add
    Dynamic_Sphere s1(0, axon_to_grow.id, axon_to_grow.begin, axon_to_grow.radius);
    // check if first sphere collides

    collides = isSphereColliding_(s1);

    bool isinside = check_borders(s1.center, s1.radius);

    if (!collides && isinside)
    {
        sphere_to_add = s1;
        return true;
    }
    else
    {
        return false;
    }
}

void Growth::find_next_center_straight(double distance, Dynamic_Sphere &s)
{
    Eigen::Vector3d last_center = axon_to_grow.spheres[axon_to_grow.spheres.size() - 1].center;
    Eigen::Vector3d before_last_center = axon_to_grow.spheres[axon_to_grow.spheres.size() - 2].center;
    Eigen::Vector3d straight_vector = (last_center - before_last_center).normalized();
    straight_vector *= distance;
    Eigen::Vector3d new_center = last_center + straight_vector;
    s.center = new_center;
}

bool Growth::TestGrowAxonAtPos(Eigen::Vector3d position_to_test, double radius_to_test)
{

    Dynamic_Sphere s(axon_to_grow.spheres.size(), axon_to_grow.id, axon_to_grow.begin, radius_to_test);
    bool iscolliding = isSphereColliding_(s);
    bool isinside = check_borders(s.center, s.radius);
    if (!iscolliding && isinside){
        return true;
    }
    else{
        return false;
    }
}

bool Growth::TestGrowAxon(Eigen::Vector3d &position_that_worked, double radius_to_test)
{

    // if sphere collides with environment
    bool collides;
    bool can_grow_;

    double distance = (radius_to_test)/2;

    if (axon_to_grow.spheres[axon_to_grow.spheres.size() -1].center[2] < voxel_size[2]) // still growing
    {

        if (axon_to_grow.spheres.size() == 0)
        {
            bool first_sphere_grown = GrowFirstSphere();

            if (first_sphere_grown)
            {
                position_that_worked = sphere_to_add.center;
                return true;
            }
            else
            {
                cout << " Cannot place axon at this location" << endl;
                return false;
            }
        }
        else
        {
            collides = true;
            can_grow_ = false;
            int tries = 0;
            Dynamic_Sphere s(axon_to_grow.spheres.size(), axon_to_grow.id, axon_to_grow.begin, radius_to_test);
            while (!can_grow_ and tries < 10000)
            {
                // find the center of next sphere by taking a random position

                if (tortuous)
                {
                    if (grow_straight==1)
                    {
                        find_next_center_straight( distance, s);
                    }
                    else
                    {
                        find_next_center(s, distance);
                    }
                    // check if there is a collision
                    collides = isSphereColliding_(s);
                    // check if the sphere is inside the voxel
                    bool isinside = check_borders(s.center, s.radius);
                    // can grow if there is no collision and if sphere is inside voxel
                    if (!collides && isinside){
                        can_grow_ = true;
                    }
                    else{
                        can_grow_ = false;
                    }
                    if (grow_straight == 1 && !can_grow_)
                    {
                        grow_straight = 0;
                    }
                    tries += 1;
                }
                else
                {
                    find_next_center(s, distance);

                    collides = isSphereColliding_(s);
                }
                
            }

            if (can_grow_)
            {

                // if we reach edge of voxel
                if (axon_to_grow.spheres[axon_to_grow.spheres.size() -1].center[2] > voxel_size[2])
                {
                    finished = true;
                }
                position_that_worked = s.center;
                return true;
                
            }
            else // collides and >1000 tries
            {
                return false;
            }
        }
    }
    else // axon is fully grown
    {
        finished = true;
        return true;
    }
}

bool Growth::GrowAxon(double radius)
{

    // if sphere collides with environment
    bool collides;
    bool can_grow_;

    double distance = (radius)/2;

    if (axon_to_grow.spheres[axon_to_grow.spheres.size() -1].center[2] < voxel_size[2]) // still growing
    {

        if (axon_to_grow.spheres.size() == 0)
        {
            bool first_sphere_grown = GrowFirstSphere();

            if (first_sphere_grown)
            {
                axon_to_grow.add_sphere(sphere_to_add);
                return true;
            }
            else
            {
                cout << " Cannot place axon at this location" << endl;
                return false;
            }
        }
        else
        {
            collides = true;
            can_grow_ = false;
            int tries = 0;
            Dynamic_Sphere s(axon_to_grow.spheres.size(), axon_to_grow.id, axon_to_grow.begin, radius);
            while (!can_grow_ and tries < 10000)
            {
                // find the center of next sphere by taking a random position

                if (tortuous)
                {
                    if (grow_straight==1)
                    {
                        find_next_center_straight( distance, s);
                    }
                    else
                    {
                        find_next_center(s, distance);
                    }
                    // check if there is a collision
                    collides = isSphereColliding_(s);
                    // check if the sphere is inside the voxel
                    bool isinside = check_borders(s.center, s.radius);
                    // can grow if there is no collision and if sphere is inside voxel
                    if (!collides && isinside){
                        can_grow_ = true;
                    }
                    else{
                        can_grow_ = false;
                    }
                    /*
                    if (!collides){
                        if (isSphereColliding_long(s)){
                            cout << "something wrong in GrowAxon (it should collide)" << endl;
                            assert(0);
                        }
                    }
                    else{
                        if (!isSphereColliding_long(s)){
                            cout << "something wrong in GrowAxon (it shouldn't collide)" << endl;
                            assert(0);
                        }
                    }
                    */

                    if (grow_straight == 1 && !can_grow_)
                    {
                        grow_straight = 0;
                    }
                    tries += 1;
                }
                else
                {
                    find_next_center(s, distance);

                    collides = isSphereColliding_(s);
                }
                
            }

            if (can_grow_)
            {

                axon_to_grow.add_sphere(s); // we want to change the radius of the sphere
                //cout << "sphere added to axon, new size :" << axon_to_grow.spheres.size() << endl;
                // if we reach edge of voxel
                if (axon_to_grow.spheres[axon_to_grow.spheres.size() -1].center[2] > voxel_size[2])
                {
                    cout << "Axon " << axon_to_grow.id << " done!!" << endl;
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
    }
    else // axon is fully grown
    {
        cout << "Axon " << axon_to_grow.id << " done!" << endl;
        finished = true;
        return true;
    }
}
