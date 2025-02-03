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
#include "Glial.h"
#include "grow_axons.h"
#include "sphere.h"
#include <thread>
#include <mutex>

class AxonGammaDistribution
{
public:
    std::vector<Axon> axons;            /*!< Vector of axons */
    std::vector<double> radii;          /*!< Axon radii */
    std::vector<Axon> growing_axons;    /*!< Axons that are growing in threads */
    std::vector<double> stuck_radii;    /*!< Radii of stuck axons, to regrow */
    std::vector<int> stuck_indices;    /*!< Indices of stuck axons, to regrow */
    std::vector<Glial> astrocytes;     /*!< Vector of astrocytes */
    std::vector<Glial> oligodendrocytes;    /*!< Vector of oligodendrocytes */
    int nbr_axons_populations;                /*!< Number of populations of axons (1-3) */
    int crossing_fibers_type;                /*!< Type of crossing fibers (0 : sheet crossing, 1 : interwoven crossing) */

    unsigned num_obstacles;             /*!< Number of axons in the substrate */
    int num_batches;                    /*!< Number of batches of axons */
    int nbr_threads;                   /*!< Number of threads to grow axons */

    double target_axons_w_myelin_icvf;        /*!< Intracellular Compartment Volume Fraction of axons with myelin */
    double target_axons_wo_myelin_icvf;         /*!< Intracellular Compartment Volume Fraction of axons without myelin */
    double target_astrocytes_soma_icvf;        /*!< Astrocyte Intracellular Compartment Volume Fraction */
    double target_astrocytes_branches_icvf;         /*!< Astrocyte Intracellular Compartment Volume Fraction */
    double target_oligodendrocytes_soma_icvf;        /*!< oligodendrocytes Intracellular Compartment Volume Fraction */
    double target_oligodendrocytes_branches_icvf;    /*!< oligodendrocytes Intracellular Compartment Volume Fraction */
    double total_volume; /*!< Total volume of the voxel */
    double target_axons_icvf;        

    double axons_w_myelin_icvf;        /*!< Intracellular Compartment Volume Fraction of axons with myelin */
    double axons_wo_myelin_icvf;         /*!< Intracellular Compartment Volume Fraction of axons without myelin */
    double astrocytes_soma_icvf;        /*!< Astrocyte Intracellular Compartment Volume Fraction */
    double astrocytes_branches_icvf;         /*!< Astrocyte Intracellular Compartment Volume Fraction */
    double oligodendrocytes_soma_icvf;        /*!< oligodendrocytes Intracellular Compartment Volume Fraction */
    double oligodendrocytes_branches_icvf;    /*!< oligodendrocytes Intracellular Compartment Volume Fraction */
    double axons_icvf;        /*!< Intracellular Compartment Volume Fraction of axons without myelin */
    double myelin_icvf;         /*!< Intracellular Compartment Volume Fraction of axons with myelin */
    double extracellular_icvf;        /*!< Extracellular Compartment Volume Fraction */

    int factor;                         /*!< Factor to divide the radii by */
    bool axon_can_shrink;               /*!< If true, the axons can shrink to allow passage between them */
    double cosPhiSquared;              /*!< Cosine of the angle between the axon and the plane squared */

    double alpha;                       /*!< Alpha coefficient of the Gamma distribution of the radii */
    double beta;                        /*!< Beta coefficient of the gamma distribution of the radii */

    int regrow_count = 0;               /*!< Number of axons to regrow */
    int regrow_thr;                     /*!< Number of regrowth batches allowed */

    Eigen::Vector3d min_limits;         /*!< Voxel min limits (if any) (bottom left corner) */
    Eigen::Vector3d max_limits;         /*!< Voxel max limits (if any) */


    double min_radius;                  /*!< Minimum radius value of all axons */
    double max_radius;                  /*!< Maximum radius value of all axons */

    double beading_variation;              /*!< For beading: percentage of variation between the maximum radius and minimum radius in each axon. If set to 1, there is no "beading" */
    double std_dev;                     /*< Standard deviation for gaussian distribution in generation of directions to grow in */ 
    int ondulation_factor;              /*!< axon ondulation factor : the number of spheres during whoch the axon grows straight before picking a direction from gaussian distribution  */
    double beading_period;           /*!< Period for the beading  */

    struct CDF {
      std::vector<double> kappas;              // Row indices (kappas)
      std::vector<double> angles;              // Column indices (angles)
      std::vector<std::vector<double>> data;   // 2D data (CDF values)
    };

    CDF cdf;
    double kappa;
    /*!
     *  \brief Initialize everything.
     */
    AxonGammaDistribution() {}

    /*!
     *  \brief Initialize everything.
     */
    AxonGammaDistribution(const double &axons_wo_myelin_icvf_, const double &axons_w_myelin_icvf_, const double &astrocytes_icvf_soma_, const double &astrocytes_icvf_branches_, const double &oligodendrocytes_icvf_soma_, const double &oligodendrocytes_icvf_branches_, const double &a, const double &b,
                                             Eigen::Vector3d &min_l, Eigen::Vector3d &max_l, const double &min_radius_,
                                              const int &regrow_thr_, const double &beading_variation_, const double &std_dev_, const int &ondulation_factor_, const double &beading_period_, const int &factor_, const bool &can_shrink_, const double &cosPhiSquared_, const double &nbr_threads_, const int &nbr_axons_populations_, const int &crossing_fibers_type_);
    
    
    /*!
       *  \param axon axon to grow
         *  \param growth Growth object with knowledge of the environment
         * \param finished 1 if the axon finiched growing, otherwise 0
         * \param grow_straight 1 if the axon is growing straight, otherwise 0
         * \param straight_growths Number of straight growths
         * \param stuck_radii_ radii of axons that got stuck
         * \param stuck_indices_ indices of axons that got stuck
     *  \brief Grows a single sphere for each axon
     */
    void growthThread(std::vector <Axon> &axs, Axon &axon, AxonGrowth &growth, int &finished, int &grow_straight, int &straight_growths, double &stuck_radii_, int &stuck_indices_);

    /*!
     *  \param axon Axon to modify
     *  \brief Causes sinusoidal fluctuation of the radii
     */
    double radiusVariation(Axon &axon);

    /*!
     *  \brief Creates a parallel growth of all axons
     */
    void parallelGrowth();


    /*!
       *  \param radiis List of axon radii
     *  \brief Generates list of radii following gamma distribution
     */
    void generate_radii(std::vector<double> &radiis, std::vector<bool> &has_myelin);

    /*!
       *  \param radii_ List of axon radii
      *  \param indices List of axon indices
      *  \param num_subsets_ List of number of axons in each batch
     *  \brief Creates a parallel growth of a batch of axons
     */
    void growBatches(std::vector<double> &radii_, std::vector<int> &indices, std::vector<int> &num_subsets_, std::vector<bool> &has_myelin, std::vector<double> &angles);
    
    /*!
     *  \param Q Starting point of axon
     *  \param D Ending point of axon
     *  \brief Gets a position in the plane for the starting point of an axon
     */
    bool get_begin_end_point(Eigen::Vector3d &Q, Eigen::Vector3d &D, double &angle);

    /*!
     *  \param num_axon Number of axons
     *  \param num_subsets Number of batches
     *  \brief Sets the number of batches based on capacity, icvf, and voxel size
     */
    void setBatches(const int &num_axons, std::vector<int> &num_subsets);

    
    /*!
     *  \param number_axons_to_grow Number of axons to grow
       *  \param radii_ List of axon radii
       * \param indices List of axon indices
       * \param stuck_radii_ radii of axons that got stuck
       * \param stuck_indices_ indices of axons that got stuck
       * \param first_index_batch First index of batch
       * \param growing_axons Axons that are growing
     *  \brief Draw the batches of growing axons
     */

    void growBatch(int &number_axons_to_grow, std::vector<double> &radii_,std::vector<int> &indices, std::vector<double> &stuck_radii_, std::vector<int> &stuck_indices_, const int &first_index_batch, std::vector<Axon> &growing_axons, std::vector<bool> &has_myelin, std::vector<double> &angles);
    
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
    bool shrinkRadius(AxonGrowth &growth, const double &radius_to_shrink, Axon &axon);

    /*!
     *  \param sph Sphere
     *  \param axs Axons to check overlapping with
     * \param gls Glial cells to check overlapping with
     *  \brief Checks if a sphere overlaps with any of the axons in axs
     */
    bool canSpherebePlaced(const Sphere &sph, const std::vector<Axon> &axs, const std::vector<Glial> &astros, const std::vector<Glial> &oligos) const;

    /*!
     *  \param axs Axons to check overlapping with
      * \param stuck_radii_ radii of axons that got stuck
     *  \brief Checks if any axon overlaps with another
     */
    bool FinalCheck(std::vector<Axon> &axs, std::vector<double> &stuck_radii_, std::vector<int> &stuck_indices_);
    /*!
     *  \param radii_ List of axon radii 
        \param indices List of axon indices
        \param first_index_batch First index of batch
        \param new_axons Axons created in batch
     *  \brief Creates a batch of axons by creating their first sphere on a plane
     */
    void createBatch(const std::vector<double> &radii_, const std::vector<int> &indices, const int &num_subset, const int &first_index_batch, std::vector<Axon> &new_axons, const std::vector<bool> &has_myelin, std::vector<double> &angles);
    
    /*!
     *  \param growing_axons Axons that have grown in batch
       *\param stuck_radii_ radii of axons that got stuck
       *\param stuck_indices_ indices of axons that got stuck
     *  \brief Checks if the axons that have grown in batch collide with each other
     */
    bool SanityCheck(std::vector<Axon>& growing_axons, std::vector<double> &stuck_radii_, std::vector<int> &stuck_indices_);

    /*!
     *  \param axon_to_grow Axon that is grown
     *  \param stuck_radii_ radii of axons that got stuck
     *  \param stuck_indices_ indices of axons that got stuck
     *  \brief Grows the entire axon
     */
    void growAxon(std::vector <Axon> &axs, Axon& axon_to_grow, double &stuck_radii_, int &stuck_indices_);

    /*!
     *  \param pos Position
     *  \param distance Distance to be inside voxel
     *  \brief Check if the position is inside the voxel
     */
    bool withinBounds(const Eigen::Vector3d &pos, const double &distance);
    
    /*!
     *  \param radius_for_axon Radius of axon.
     *  \param Q Starting point of axon.
     *  \param D Ending point of axon.
     *  \param new_axons Newly placed axons.
     *  \brief Places an axon in voxel.
     */
    bool PlaceAxon(const int &axon_id, const double &radius_for_axon, const Eigen::Vector3d &Q, const Eigen::Vector3d &D, std::vector<Axon> &new_axons, const bool &has_myelin, const double &angle, const bool &outside_voxel);

    /*!
     *  \param out Output stream to write SWC data.
     *  \param overlapping_factor Distance between spheres is radius / overlapping_factor.
     *  \brief Writes a SWC file to save the position of axons.
     */
    void create_SWC_file(std::ostream &out);

    /*!
     *  \param out Output stream to write simulation details.
     *  \param duration Duration of the growth in seconds.
     *  \brief Writes to a file some details on the simulation (duration, etc.).
     */
    void simulation_file(std::ostream &out, const std::chrono::seconds &duration);

    private:


    /*!
     *  \brief Places the glial cells in the voxel.
     */
    void PlaceGlialCells();


    void growGlialCell(Glial &glial_cell, const int &number_ramification_points, int &nbr_spheres);


    bool growExtraBranchesinGlialCells(Glial &glial_cell, int &nbr_spheres);

    /*!
     *  \brief Add myelin sheath by creating an inner_axonal membrane
     */

    

    void add_Myelin();


    /*!
     *  \brief Check that processes of a sphere do not overlap with each other
         \param sph Sphere in a process to check
         \param glial_cell_to_grow Glial cell growing
     
     */
    bool collideswithOtherBranches(const Sphere &sph, const Glial &glial_cell_to_grow);

    /*!
     *  \brief Grow all glial cells
     */


    void GrowAllGlialCells();

    /*!
     *  \brief Grow all axons
     */

    void GrowAllAxons();

    /*!
     *  \brief Swells all axons and updates the ICVF
    */

    void SwellAxons(const double &percentage);

    /*!
     *  \brief Swells all spheres in axons
         \param ax Axon to swell
    */

    void SwellAxon(Axon &ax, const double &percentage);

    /*!
     *  \brief Swell a sphere if possible by increasing its radius by 1%
    */

    bool SwellSphere(Sphere &sph, const double &percentage);

    /*!
     *  \brief Finds a random point on the opposite plane based on the current point and cosPhiSquared (c2)
        \param x1 x coordinate of the point
        \param y1 y coordinate of the point
        \param L distance between the two planes
     */
    Eigen::Vector3d randomPointOnPlane(const Eigen::Vector3d &begin, const Eigen::Vector3d &end, const int &axis1, const int &axis2, const int &axis3, double &angle, bool &outside_voxel);

    double findInnerRadius(const double &outerRadius);

    bool growProcessFromSoma(Glial &glial_cell, const int &j, const int &nbr_spheres);

    double RandomradiusVariation(Axon &axon);

    void growOligodendrocytesBranches(std::vector<Glial>& oligodendrocytes, size_t astrocytes_size, std::vector<int>& nbr_spheres);

    void growExtraBranchesAstrocytes(std::vector<int>& nbr_spheres);

    void ICVF(const std::vector<Axon> &axs, const std::vector<Glial> &astrocytes, const std::vector<Glial> &oligos);
    
    void process_point(const Eigen::Vector3d &point, const std::vector<Axon> &axs, const std::vector<Glial> &astrocytes, const std::vector<Glial> &oligos, int &axons_count, int &inner_axons_count, int &astrocytes_somas_count, int &astrocytes_process_count, int &oligos_somas_count, int &oligos_process_count, int &extracellular_count);
    std::vector<Sphere> addIntermediateSpheres(const Sphere &random_sphere, const Sphere &first_sphere, const int &branch_nbr, const int &nbr_spheres, const int &nbr_spheres_between);
    bool GenerateFirstSphereinProcess(Sphere &first_sphere, Eigen::Vector3d &attractor, const double &radius, const Sphere &sphere_to_emerge_from, const Eigen::Vector3d &vector_to_prev_center, const int &nbr_spheres, const int &nbr_spheres_between, const int &cell_id, const int &branch_id);
    double draw_angle(double kappa);
    double c2toKappa(double c2, double tol, double kappa_min, double kappa_max);
    std::vector<double> generate_angles(const int &num_samples);
    void processBatchWithThreadPool(int nbr_threads, std::vector<Axon>& axons_to_grow, std::vector<double>& stuck_radii, std::vector<int>& stuck_indices);
    bool check_borders(const Eigen::Vector3d&  min_l, const Eigen::Vector3d&  max_l, const Eigen::Vector3d& pos, const double& distance_to_border);
    bool pushAxonSpheres(Axon &axon, const Sphere &sph);
    bool PushOtherAxons(const Sphere &sph);
};


#endif // AXONGAMMADISTRIBUTION_H

