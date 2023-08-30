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
#include <mutex>

class AxonGammaDistribution
{
public:
    std::vector<Axon> axons;   /*!< Axon vector  */
    std::vector<double> radii; /*!< axons radii  */

    std::vector<Axon> axons_to_regrow;
    std::vector<double> stuck_radii; /*!< radii of stuck axons */

    unsigned num_obstacles; /*!< number of cylnders fit inside the substrate */
    int num_batches;
    int axon_capacity; /* Safe number of axon per batch to avoid crash */
    std::mutex axonsMutex;
    std::mutex stuckMutex;

    double alpha; /*!< alpha coefficient of the Gamma distribution                           */
    double beta;  /*!< beta coefficient of the gamma distribution                            */
    double icvf;
    int regrow_count = 0; /*!< number of axons to regrow */
    int regrow_thr;       /*!< Number of regrowth batches allowed*/

    Eigen::Vector3d min_limits; /*!< voxel min limits (if any) (bottom left corner)                     */
    Eigen::Vector3d max_limits; /*!< voxel max limits (if any)                                          */

    bool tortuous;
    bool draw;

    std::vector<double> tortuosities; /*!< ODF                                                          */
    double min_radius;
    std::vector<GLfloat> colours;

    double max_radius; // between all axons in env

    double variation_perc = 0.55;
    //double variation_perc = 1;

    bool can_shrink;

    /*!
     *  \brief Initialize everything.
     */
    AxonGammaDistribution() {}

    /*!
     *  \brief Initialize everything.
     */
    AxonGammaDistribution(unsigned &, int &, double, double, Eigen::Vector3d &, Eigen::Vector3d &, double, bool, bool, int, bool);

    /*!
     *  \brief Sets icvf and overwrites num_obstacles and num_batches
     */
    void set_icvf(double icvf_, double x, double y);

    /*!
     *  \brief Shows a small histogram of the gamma distribution
     */
    void displayGammaDistribution();

     /*!

     *  \brief Generates a radius from gamma distribution
     */

    double generate_radius();

    double computeAreaICVF(std::vector<Axon> new_axons);

    /*!
     *  \param radii list of radii associated to axons
     *  \brief Creates a list with all the axons to grow
     */
    void createAxons(std::vector<double> &radii, std::vector<Axon> &new_axons, bool regrowth);

    /*!
     *  \param index position in the 2d form of the axons vector
     *  \param can_grow assesses if growth possible
     *  \param finished assesses if growth finished, true = 1 and false = 0
     *  \param grow_straight determines growth direction, true = 1 and false = 0
     *  \param stuck number of times growth was not possible
     *  \param stuck number of straight growths
     *  \brief Grows a single sphere for each axon
     */
    void growthThread(Axon &axon, int &finished, int &grow_straight, int &straight_growths, int time, int &shrink_tries, int &restart_tries);

    /*!
     *  \brief Causes sinusoidal fluctuation of the radii
     */
    double radiusVariation(Axon &axon, int time);

    /*!
     *  \brief Creates a parallel growth of all axons
     */
    void parallelGrowth();

    /*!
     *  \brief Creates and displays a parallel growth of all axons
     */
    void growthVisualisation();
    void growBatches(std::vector<Axon> &ax_list, std::vector<double> &radii_, std::vector<int> &num_subsets_);
    void setBatches(int num_axons, std::vector<int> &num_subsets);
    void drawBatches(sf::Window &window, std::vector<Axon> &ax_list, std::vector<double> &radii_, std::vector<int> &num_subsets_,float zoomLevel,bool isDragging, sf::Vector2i lastMousePos, bool isRightDragging, sf::Vector2i lastRightMousePos, sf::Vector2i currentMousePos, sf::Vector2i mouseDelta, sf::Vector2i prevousDisplacement, sf::Vector2i rightMouseDelta, sf::Vector2i prevousRotation, float rotationFactor, float displacementFactor);
    double computeICVF_();
    void createSubstrate();

    /*!
     *  \brief Shrinks the radius to allow passage between axons
     */
    bool shrinkRadius(double radius_to_shrink, Axon &axon, int grow_straight);

    /*!
     *  \brief Finds a radius for which shrinkage allows passage
     */
    void dichotomy(Eigen::Vector3d position_that_worked, Axon axon, double initial_rad, double &last_rad, int grow_straight);

    /*!
     *  \param row row of the current batch's axon vector
     *  \param window drawing window
     *  \brief draws sequential axon growth
     */
    void drawWorld(std::vector<Axon> ax_list, unsigned int row, sf::Window &window, int num_ax);

    /*!
     *  \brief Generates list of radii following gamma distribution
     */
    void generate_radii(std::vector<double> &radiis);

    /*!
     *  \brief Draws 2d map of axon's initial placement
     */
    void axonDensityMap();

    void create_SWC_file(std::ostream &out);

    void axons_file(std::ostream &out);
    void simulation_file(std::ostream &out, std::chrono::seconds duration);

    /*!
     *  \brief Prints the cylinders positions in a file or output stream.
     *  \param out ostream where to write the info.
     */
    void printSubstrate(std::ostream &out);

    bool check_borders(Eigen::Vector3d pos, double distance_to_border, Eigen::Vector2d &twin_delta_pos);

    void get_begin_end_point(Eigen::Vector3d &Q, Eigen::Vector3d &D);
    void PlaceAxon(double radius, bool regrowth, int i, int &tries, bool &next, Eigen::Vector3d Q, Eigen::Vector3d D, std::vector<Axon> &all_axons, std::vector<Axon> &new_axons);
    
    void PlaceTwinAxons(double radius, bool regrowth, int i, int &tries, bool &next, std::vector<Eigen::Vector3d> Qs, std::vector<Axon> &all_axons, std::vector<Axon> &new_axons);
    std::vector<Eigen::Vector3d> FindTwins(Eigen::Vector3d Q, double rad);
    bool withinBounds(Eigen::Vector3d pos, double distance);

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
