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
    std::vector<Axon> axons; /*!< Axon vector                                                           */
    unsigned num_obstacles;  /*!< number of cylnders fit inside the substrate */
    int num_batches = 5;
    std::mutex axonsMutex;
    double alpha; /*!< alpha coefficient of the Gamma distribution                           */
    double beta;  /*!< beta coefficient of the gamma distribution                            */
    double icvf;

    Eigen::Vector3d min_limits; /*!< voxel min limits (if any) (bottom left corner)                     */
    Eigen::Vector3d max_limits; /*!< voxel max limits (if any)                                          */

    bool tortuous;
    bool draw;

    std::vector<double> tortuosities; /*!< ODF                                                          */
    double min_radius;
    std::vector<GLfloat> colours;

    double max_radius;

    /*!
     *  \brief Initialize everything.
     */
    AxonGammaDistribution() {}

    /*!
     *  \brief Initialize everything.
     */
    AxonGammaDistribution(unsigned &, int &, double, double, Eigen::Vector3d &, Eigen::Vector3d &, double, bool, bool);

    /*!
     *  \brief Sets icvf and overwrites num_obstacles and num_batches
     */
    void set_icvf(double icvf_, double x, double y);

    /*!
     *  \brief Shows a small histogram of the gamma distribution
     */
    void displayGammaDistribution();

    /*!
     *  \brief Creates and displays sequential axonal growth
     */
    void createGammaSubstrate();

    /*!
     *  \param radii list of radii associated to axons
     *  \brief Creates a list with all the axons to grow
     */
    void createAxons(std::vector<double> radii);

    /*!
     *  \param index position in the 2d form of the axons vector
     *  \param can_grow assesses if growth possible
     *  \param finished assesses if growth finished, true = 1 and false = 0
     *  \param grow_straight determines growth direction, true = 1 and false = 0
     *  \param stuck number of times growth was not possible
     *  \param stuck number of straight growths
     *  \brief Grows a single sphere for each axon
     */
    void growthThread(int index, bool &can_grow, int &finished, int &grow_straight, int &stuck, int &straight_growths, int time, double radius);

    void radiusVariation(Axon &axon, int time, double radius);

    /*!
     *  \brief Creates and displays a parallel growth of all axons
     */
    void parallelGrowth();

    /*!
     *  \param ax axon to draw
     *  \param window drawing window
     *  \param colour axon colour
     *  \brief draws sequential axon growth
     */
    void drawWorld_seq(Axon *ax, sf::Window &window, GLfloat colour);

    /*!
     *  \param row row of the current batch's axon vector
     *  \param window drawing window
     *  \brief draws sequential axon growth
     */
    void drawWorld(unsigned int row, sf::Window &window);

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
