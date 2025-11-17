#include "obstacle.h"
#include <math.h>
#include <unordered_set>
#include <algorithm>

using namespace Eigen;
using namespace std;

Obstacle::Obstacle():percolation(0),T2(0),id(-1)
{
}

// Function to check if a point is inside a dilated box
bool Obstacle::check_borders(const Eigen::Vector3d&  min_l, const Eigen::Vector3d&  max_l, const Eigen::Vector3d& pos, const double& distance_to_border) {

    
    // Check if the point is inside the dilated box
    for (int i = 0; i < 3; ++i) {
        double min_bound = min_l[0] - distance_to_border;
        double max_bound = max_l[1] + distance_to_border;
        if (pos[i] < min_bound || pos[i] > max_bound) {
            return false; // Point is outside the dilated box
        }
    }
    
    return true; // Point is inside the dilated box
}


std::vector<int> Obstacle::findCommonIntegers(const std::vector<int>& vec1, const std::vector<int>& vec2, const std::vector<int>& vec3) {
    std::vector<int> result;
    std::set_intersection(vec1.begin(), vec1.end(),
                          vec2.begin(), vec2.end(),
                          std::back_inserter(result));
    std::vector<int> commonIntegers;
    std::set_intersection(result.begin(), result.end(),
                          vec3.begin(), vec3.end(),
                          std::back_inserter(commonIntegers));
    return commonIntegers;
}


bool Obstacle::check_position_is_in_Box(Eigen::Vector3d position, std::vector<Eigen::Vector2d> Box, double buffer_distance){
    if ((position[0] >=  Box[0][0]- buffer_distance)  && (position[0] <= Box[0][1] + buffer_distance)){
        if ((position[1] >= Box[1][0] - buffer_distance) && (position[1] <=  Box[1][1] + buffer_distance)){
            return true;
        }
    }
    return false;
}