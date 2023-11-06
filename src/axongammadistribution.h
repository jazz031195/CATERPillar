//! AxonGammaDistribution Class ================================================================ /
/*!
*   \details   This class constructs a substrate taken from a Gamma distribution of radii placed in
*              a single voxel structure.
*   \author    Jasmine Nguyen-Duc
*   \date      September 2023
=============================================================================================== */

#ifndef AXONGAMMADISTRIBUTION_H
#define AXONGAMMADISTRIBUTION_H

#include "Eigen/Core"
#include <vector>
#include "constants.h"
#include <iostream>
#include "Axon.h"
#include "grow_axons.h"
#include "sphere.h"
#include <thread>
#include <GL/gl.h>
#include <GL/glu.h>
#include <SFML/Graphics.hpp>
#include <mutex>

class AxonGammaDistribution
{
public:
    std::vector<Axon> axons;            /*!< Vector of axons */
    std::vector<double> radii;          /*!< Axon radii */
    std::vector<Axon> growing_axons;    /*!< Axons that are growing in threads */
    std::vector<double> stuck_radii;    /*!< Radii of stuck axons, to regrow */
    unsigned num_obstacles;             /*!< Number of axons in the substrate */
    int num_batches;                    /*!< Number of batches of axons */
    int axon_capacity;                  /*!< Number of axons per batch (do not choose more than the maximum number of cores to avoid a crash) */
    double icvf;                        /*!< Intracellular Compartment Volume Fraction */
    double target_icvf;                 /*!< Target Intracellular Compartment Volume Fraction */
    bool tortuous;                      /*!< If true, the axons are tortuous; otherwise, they are straight */

    std::vector<double> tortuosities;   /*!< Orientation dispersion fraction of each axon */
    
    double alpha;                       /*!< Alpha coefficient of the Gamma distribution of the radii */
    double beta;                        /*!< Beta coefficient of the gamma distribution of the radii */

    int regrow_count = 0;               /*!< Number of axons to regrow */
    int regrow_thr;                     /*!< Number of regrowth batches allowed */

    Eigen::Vector3d min_limits;         /*!< Voxel min limits (if any) (bottom left corner) */
    Eigen::Vector3d max_limits;         /*!< Voxel max limits (if any) */

    bool draw;                          /*!< If true, the voxel is drawn in a window */

    double min_radius;                  /*!< Minimum radius value of all axons */
    double max_radius;                  /*!< Maximum radius value of all axons */

    std::vector<GLfloat> colours;       /*!< List of colors for the visualization of axons */

    double beading_variation;              /*!< For beading: percentage of variation between the maximum radius and minimum radius in each axon. If set to 1, there is no "beading" */
    double std_dev;                     /*< Standard deviation for gaussian distribution in generation of directions to grow in */ 
    int ondulation_factor;              /*!< axon ondulation factor : the number of spheres during whoch the axon grows straight before picking a direction from gaussian distribution  */
    double beading_frequency;           /*!< Frequency for the beading  */
    /*!
     *  \brief Initialize everything.
     */
    AxonGammaDistribution() {}

    /*!
     *  \brief Initialize everything.
     */
    AxonGammaDistribution(double, int &, double, double, Eigen::Vector3d &, Eigen::Vector3d &, double, bool, bool, int, double, double, int, double);
    
    /*!
     *  \param window Window to visualize growing axons
     *  \param stuck_radii_ radii of axons that got stuck
     *  \param radii_ radii of all axons to grow
     *  \param number_axons_to_grow number of axons to grow
     *  \param first_index_batch first axon in radii_ to grow 
     *  \param num_subsets_ List of batches with number of axons to grow in each 
     *  \param regrowth If true, the axons are regrowing after failing in previous batches
     *  \param growing_axons List of axons that have grown
     *  \param other_parameters For visualization purposes
     *  \brief Draw the batches of growing axons
     */
    void drawBatch(sf::Window &window, int number_axons_to_grow, std::vector<double> &radii_, std::vector<double> &stuck_radii_, bool regrowth, int first_index_batch, std::vector<Axon> &growing_axons, float &zoomLevel,bool &isDragging, sf::Vector2i &lastMousePos, bool &isRightDragging, sf::Vector2i &lastRightMousePos, sf::Vector2i &currentMousePos, sf::Vector2i &mouseDelta, sf::Vector2i &prevousDisplacement, sf::Vector2i &rightMouseDelta, sf::Vector2i &prevousRotation, float &rotationFactor, float &displacementFactor);
    /*!
     *  \param window Window to visualize growing axons
     *  \param ax_list List of axons
     *  \param num_subsets_ List of batches with number of axons to grow in each 
     *  \param regrowth If true, the axons are regrowing after failing in previous batches
     *  \param other_parameters For visualization purposes
     *  \brief Draw the batches of growing axons
     */
    void drawBatches(sf::Window &window, std::vector<Axon> &ax_list, std::vector<double> radii_, std::vector<int> num_subsets_, float zoomLevel, bool isDragging, sf::Vector2i lastMousePos, bool isRightDragging, sf::Vector2i lastRightMousePos, sf::Vector2i currentMousePos, sf::Vector2i mouseDelta, sf::Vector2i prevousDisplacement, sf::Vector2i rightMouseDelta, sf::Vector2i prevousRotation, float rotationFactor, float displacementFactor, bool regrowth);
    
     /*!
     *  \param ax_list List of axons to draw.
     *  \brief Draws sequential axon growth.
     */
    void drawWorld(std::vector<Axon> ax_list);

    /*!
     *  \brief Updates the window for visualization with respect to zooming, dragging, etc.
     */
    void update_window(sf::Window &window, float zoomLevel, sf::Vector2i prevousDisplacement, sf::Vector2i prevousRotation, float rotationFactor, float displacementFactor);


    /*!
     *  \param new_axons Vector of axons to use
     *  \brief Calculate the intracellular area / total area while using only the first sphere of each axon
     */
    double computeAreaICVF(std::vector<Axon> new_axons);
    
    /*!
     *  \param axon Axon to grow one extra sphere on
     *  \param growth Grwoth object with environment in memory
     *  \param finished Assesses if growth finished, true = 1, and false = 0
     *  \param grow_straight Determines whether the growth is straight or in a random direction, true = 1, and false = 0
     *  \param straight_growths Last number of spheres in a row that grow in a straight line (4 max)
     *  \param regrowth If true, the axons are regrowing after failing in previous batches
     *  \brief Grows a single sphere for each axon
     */
    void growthThread(Axon &axon, Growth &growth, int &finished, int &grow_straight, int &straight_growths, bool regrowth, std::vector<double> &stuck_radii_);

    /*!
     *  \param axon Axon to modify
     *  \brief Causes sinusoidal fluctuation of the radii
     */
    double radiusVariation(Axon axon);

    /*!
     *  \brief Creates a parallel growth of all axons
     */
    void parallelGrowth();

    /*!
     *  \brief Creates and displays a parallel growth of all axons
     */
    void growthVisualisation();

    /*!
     *  \brief Generates list of radii following gamma distribution
     */
    void generate_radii(std::vector<double> &radiis);

    /*!
     *  \param axon Axon to grow one extra sphere on
     *  \param finished Assesses if growth finished, true = 1, and false = 0
     *  \param grow_straight Determines whether the growth is straight or in a random direction, true = 1, and false = 0
     *  \param straight_growths Last number of spheres in a row that grow in a straight line (4 max)
     *  \param regrowth If true, the axons are regrowing after failing in previous batches
     *  \brief Creates a parallel growth of a batch of axons
     */
    void growBatches(std::vector<double> &radii_, std::vector<int> &num_subsets_, bool regrowth);
    
    /*!
     *  \param Q Starting point of axon
     *  \param D Ending point of axon
     *  \brief Gets a position in the plane for the starting point of an axon
     */
    void get_begin_end_point(Eigen::Vector3d &Q, Eigen::Vector3d &D);

    /*!
     *  \param num_axon Number of axons
     *  \param num_subsets Number of batches
     *  \brief Sets the number of batches based on capacity, icvf, and voxel size
     */
    void setBatches(int num_axons, std::vector<int> &num_subsets);

    
    /*!
     *  \param window Window to visualize growing axons
     *  \param stuck_radii_ radii of axons that got stuck
     *  \param radii_ radii of all axons to grow
     *  \param number_axons_to_grow number of axons to grow
     *  \param first_index_batch first axon in radii_ to grow 
     *  \param num_subsets_ List of batches with number of axons to grow in each 
     *  \param regrowth If true, the axons are regrowing after failing in previous batches
     *  \param growing_axons List of axons that have grown
     *  \brief Draw the batches of growing axons
     */

    void growBatch(int number_axons_to_grow, std::vector<double> &radii_, std::vector<double> &stuck_radii_, bool regrowth, int first_index_batch, std::vector<Axon> &growing_axons);
    
    /*!
     *  \brief Creates the entire substrate
     */
    void createSubstrate();

    
    /*!
     *  \param growth Growth object with knowledge of the environment
     *  \param radius_to_shrink Radius to shrink
     *  \param axon Axon to shrink
     *  \brief Shrinks the radius to allow passage between axons
     */
    bool shrinkRadius(Growth& growth, double radius_to_shrink, Axon &axon);

    /*!
     *  \brief Swells axons
     */
    bool swellAxons();

    /*!
     *  \param new_sphere Sphere to swell
     *  \param swelling_perc Percentage of volume to gain during swelling
     *  \brief Finds the best swelling percentage so that the sphere doesn't overlap with the environment
     */
    double dichotomy_swelling(Sphere new_sphere, double swelling_perc);

    /*!
     *  \param overlapping_factor Distance between spheres is radius/overlapping_factor
     *  \param final_axons Resulting axons
     *  \brief Adds spheres in between existing ones
     */
    // Function to calculate intermediate points and insert N lines
    void fill_with_overlapping_spheres(int overlapping_factor, std::vector<Axon> &final_axons);

    /*!
     *  \param sph Sphere
     *  \param axs Axons to check overlapping with
     *  \brief Checks if a sphere overlaps with any of the axons in axs
     */
    bool canSpherebePlaced(Sphere sph, std::vector<Axon> axs);

    /*!
     *  \param axs Axons to check 
     *  \brief Checks if any axon overlaps with another
     */
    bool FinalCheck(std::vector<Axon>& axs);
    /*!
     *  \param radii_ List of axon radii 
        \param num_subset Number of batches
        \param regrowth If true, the axons are regrowing after failing in previous batches
        \param first_index_batch Index of batch
        \param new_axons Axons created in batch
     *  \brief Creates a batch of axons by creating their first sphere on a plane
     */
    void createBatch(std::vector<double> radii_, int num_subset, bool regrowth, int first_index_batch, std::vector<Axon> &new_axons);
    
    /*!
     *  \param growing_axons Axons that have grown in batch
     *  \brief Checks if the axons that have grown in batch collide with each other
     */
    bool SanityCheck(std::vector<Axon>& growing_axons, std::vector<double>& stuck_radii_);

    /*!
     *  \param axon_to_grow Axon that is grown
     *  \param regrowth If true, the axons are regrowing after failing in previous batches
     *  \brief Grows the entire axon
     */
    void growAxon(Axon& axon_to_grow, bool regrowth, std::vector<double>& stuck_radii_);

    /*!
     *  \param pos Position
     *  \param distance Distance to be inside voxel
     *  \brief Check if the position is inside the voxel
     */
    bool withinBounds(Eigen::Vector3d pos, double distance);
    
    /*!
     *  \param radius Radius of axon.
     *  \param regrowth If true, the axons are regrowing after failing in previous batches.
     *  \param i Index value.
     *  \param tries Number of tries before axon is placed.
     *  \param Q Starting point of axon.
     *  \param D Ending point of axon.
     *  \param new_axons Newly placed axons.
     *  \brief Places an axon in voxel.
     */
    bool PlaceAxon(double radius_for_axon, double radius_for_first_sphere, bool regrowth, int i, int &tries, Eigen::Vector3d Q, Eigen::Vector3d D, std::vector<Axon> &new_axons);

    /*!
     *  \param out Output stream to write SWC data.
     *  \param overlapping_factor Distance between spheres is radius / overlapping_factor.
     *  \brief Writes a SWC file to save the position of axons.
     */
    void create_SWC_file(std::ostream &out, int overlapping_factor);

    /*!
     *  \param out Output stream to write simulation details.
     *  \param duration Duration of the growth in seconds.
     *  \brief Writes to a file some details on the simulation (duration, etc.).
     */
    void simulation_file(std::ostream &out, std::chrono::seconds duration);

    private:

    /*!
     *  \brief Computes Intra Cellular Volum Fraction given the voxel limits and the list of added cylinders.
     */
    double computeICVF();

};

#endif // AXONGAMMADISTRIBUTION_H

