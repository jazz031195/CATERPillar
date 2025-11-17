//!  Obstacle Base Class ==============================================================================/
/*!
*   \details   Father class to define the base of any other obstacle (wall or substrate)
*   \author    Jonathan Rafael
*   \date      November 2016
*   \version   1.42
 =====================================================================================================*/

#ifndef OBSTACLE_H
#define OBSTACLE_H

#include "Eigen/Dense"
#include <vector>
#include "Eigen/Core"

class Obstacle
{
public:

    int id;                         /*!< Unique id of the simulation                                                */
    int count_perc_crossings;       /*!< Auxiliar value to count the number of percolatin crossings in a simulation */
    double percolation;             /*!< Percolation value between 0 and 1.                                         */
    double T2;                      /*!< T2 decay, not used by default                                              */

    /*! \fn  Obstacle
     *  \brief Default constructor. Does nothing.
     */
    Obstacle();
    static std::vector<int> findCommonIntegers(const std::vector<int>& vec1, const std::vector<int>& vec2, const std::vector<int>& vec3);
    static bool check_borders(const Eigen::Vector3d&  min_l, const Eigen::Vector3d&  max_l, const Eigen::Vector3d& pos, const double& distance_to_border);
    static bool check_position_is_in_Box(Eigen::Vector3d position, std::vector<Eigen::Vector2d> Box, double buffer_distance);


};

#endif // OBSTACLE_H
