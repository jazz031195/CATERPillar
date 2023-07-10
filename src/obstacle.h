//!  Obstacle Base Class ==============================================================================/
/*!
*   \details   Father class to define the base of any other obstacle (wall or substrate)
*   \author    Jonathan Rafael
*   \date      November 2016
*   \version   1.42
 =====================================================================================================*/

#ifndef OBSTACLE_H
#define OBSTACLE_H

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



};

#endif // OBSTACLE_H
