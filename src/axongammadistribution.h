//!  CylinderGammaDistribution Class =============================================================/
/*!
*   \details   Class to construct a substrate taken from a Gamma distribution of radiis placed in
*              a single voxel structure.
*   \author    Jasmine Nguyen-Duc
*   \date      Septembre 2022
*   \version   0.2
=================================================================================================*/

#ifndef AXONGAMMADISTRIBUTION_H
#define AXONGAMMADISTRIBUTION_H

#include "Eigen/Core"
#include <vector>
#include "constants.h"
#include <iostream>
#include "Axon.h"
#include "grow_axons.h"
#include "dynamic_sphere.h"
#include <thread>
#include <GL/gl.h>
#include <GL/glu.h>
#include <SFML/Graphics.hpp>

class AxonGammaDistribution
{
public:
    std::vector<Axon> axons; /*!< Axon vector                                                            */
    unsigned num_obstacles;  /*!< number of cylnders fit inside the substrate                                */
    double alpha;            /*!< alpha coefficient of the Gamma distribution                                */
    double beta;             /*!< beta coefficient of the gamma distribution                                 */
    double icvf;

    Eigen::Vector3d min_limits; /*!< voxel min limits (if any) (bottom left corner)                             */
    Eigen::Vector3d max_limits; /*!< voxel max limits (if any)                                                  */

    bool tortuous;

    std::vector<double> tortuosities; /*!< ODF                                               */
    double min_radius;
    std::vector<GLfloat> colours;

    double max_radius;

    /*!
     *  \param P_ Cylinder origin
     *  \param Q_ cylinder direction.
     *  \param radius_ cylinder's radius
     *  \param scale scale factor for the values passed. Useful when reading a file.
     *  \brief Initialize everything.
     */
    AxonGammaDistribution() {}

    /*!
     *  \param P_ Cylinder origin
     *  \param Q_ cylinder direction.
     *  \param radius_ cylinder's radius
     *  \param scale scale factor for the values passed. Useful when reading a file.
     *  \brief Initialize everything.
     */
    AxonGammaDistribution(unsigned &, double, double, Eigen::Vector3d &, Eigen::Vector3d &, double, bool);

    /*!
     *  \brief Shows a small histogram of the gamma distribution
     */
    void displayGammaDistribution();

    /*!
     *  \brief Samples and constructs a Gamma distribution
     */
    void createGammaSubstrate();

    /*!
     *  \param radii list of radii associated to axons
     *  \brief Creates a list with all the axons to grow
     */
    void createAxons(std::vector<double> radii);

    /*!
     *  \brief Creates and displays a parallel growth of all axons 
     */
    void parallelGrowth(); 

    /*!
     *  \param ax axon to grow
     *  \param can_grow assesses if growth possible
     *  \param finished assesses if growth finished
     *  \brief Grows a single sphere for each axon
     */
    void growthThread(Axon &ax, bool &can_grow, bool &finished);

    void drawWorld(sf::Window &window, GLfloat colour);

    void generate_radii(std::vector<double> &radiis);
    /*!
     *  \brief Prints the cylinders positions in a file or output stream.
     *  \param out ostream where to write the info.
     */
    void printSubstrate(std::ostream &out);
    bool check_borders(Eigen::Vector3d pos, double distance_to_border, Eigen::Vector2d &twin_delta_pos);
    void get_begin_end_point(Eigen::Vector3d &Q, Eigen::Vector3d &D);

private:
    /*!
     *  \brief Computes Intra Celular Volum Fraction given the voxel limits and the list of added cylinders.
     *  \param cylinders List of included cylinders.
     *  \param min_limits voxel min limits.
     *  \param max_limits voxel max limits.
     */
    double computeICVF();

    void computeMinimalSize(std::vector<double> radiis, double icvf_, Eigen::Vector3d &l);
};

#endif // AXONGAMMADISTRIBUTION_H
