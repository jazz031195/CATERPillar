//! CaterpillarGrowth Class ================================================================ /
/*!
*   \details   This class constructs a substrate taken from a Gamma distribution of radii placed in
*              a single voxel structure.
*   \author    Jasmine Nguyen-Duc
*   \date      September 2023
=============================================================================================== */

#ifndef CaterpillarGrowth_H
#define CaterpillarGrowth_H

#include "Eigen/Core"
#include <vector>
#include "constants.h"
#include "parameters.h"
#include <iostream>
#include "Axon.h"
#include "Glial.h"
#include "Blood_Vessel.h"
#include "grow_axons.h"
#include "sphere.h"
#include "SphereGrid.h"
#include <thread>
#include <mutex>
#include <random>
#include <functional>

class CaterpillarGrowth
{
public:
    std::mt19937 gen;

    /*!
     *  \brief Optional GUI hook: invoked from growAxonsLayered with (depth
     *         travelled so far, total box depth) every layer, so a caller
     *         (e.g. the Qt GUI) can drive a progress bar without CATERPillar
     *         itself depending on Qt. Left unset (empty std::function) by
     *         default -- CLI usage is entirely unaffected.
     */
    std::function<void(double completed_depth, double total_depth)> on_growth_progress;

    /*!
     *  \brief Optional GUI hook: invoked from SwellAxons with (current axon
     *         ICVF, target axon ICVF) every round of both the in-place and
     *         push-assisted swelling phases. Same opt-in/no-op-by-default
     *         design as on_growth_progress.
     */
    std::function<void(double current_icvf, double target_icvf)> on_swelling_progress;

    std::vector<Axon> axons;            /*!< Vector of axons */
    std::vector<double> radii;          /*!< Axon radii */
    std::vector<double> stuck_radii;    /*!< Radii of stuck axons, to regrow */
    std::vector<int> stuck_indices;    /*!< Indices of stuck axons, to regrow */
    std::vector<Glial> glial_pop1;     /*!< Vector of glial_pop1 */
    std::vector<Glial> glial_pop2;    /*!< Vector of glial_pop2s */
    std::vector<Glial> glial_pop3;    /*!< Vector of glial_pop3s */
    std::vector<Blood_Vessel> blood_vessels; /*!< Vector of blood vessels */

    SphereGrid sphere_grid;              /*!< Spatial grid indexing every sphere added to the environment */
    SphereGrid batch_scratch_grid;       /*!< Persistent scratch grid reused (via clear()) by SanityCheck every sub-batch, instead of rebuilding one each call */
    double grid_voxel_size;              /*!< Voxel edge length used to size sphere_grid (and any scratch grids) */

    int nbr_axons_populations;                /*!< Number of populations of axons (1-3) */
    int crossing_fibers_type;                /*!< Type of crossing fibers (0 : sheet crossing, 1 : interwoven crossing) */
    
    int nbr_threads;                   /*!< Number of threads to grow axons */
    
    double glial_pop2_radius_mean; /*!< Mean radius of glial_pop2 */
    double glial_pop1_radius_mean;       /*!< Mean radius of glial_pop1 */
    double glial_pop3_radius_mean; /*!< Mean radius of glial_pop3 */
    double glial_pop2_radius_std; /*!< Standard deviation of glial_pop2 */
    double glial_pop1_radius_std;       /*!< Standard deviation of glial_pop1 */
    double glial_pop3_radius_std; /*!< Standard deviation of glial_pop3 */
    bool glial_pop1_branching;               /*!< If true, the glial_pop1 can branch */
    bool glial_pop2_branching;               /*!< If true, the glial_pop2 can branch */
    bool glial_pop3_branching;               /*!< If true, the glial_pop3 can branch */

    double target_axons_w_myelin_icvf;        /*!< Intracellular Compartment Volume Fraction of axons with myelin */
    double target_axons_wo_myelin_icvf;         /*!< Intracellular Compartment Volume Fraction of axons without myelin */
    double target_glial_pop1_soma_icvf;        /*!< Astrocyte Intracellular Compartment Volume Fraction */
    double target_glial_pop1_processes_icvf;         /*!< Astrocyte Intracellular Compartment Volume Fraction */
    double target_glial_pop2_soma_icvf;        /*!< glial_pop2 Intracellular Compartment Volume Fraction */
    double target_glial_pop2_processes_icvf;    /*!< glial_pop2 Intracellular Compartment Volume Fraction */
    double target_glial_pop3_soma_icvf;        /*!< glial_pop3 Intracellular Compartment Volume Fraction */
    double target_glial_pop3_processes_icvf;    /*!< glial_pop3 Intracellular Compartment Volume Fraction */
    double target_blood_vessels_icvf;    /*!< Blood Vessel Compartment Volume Fraction -- target is interpreted relative to the (padded) blood-vessel voxel, see big_total_volume */
    double total_volume; /*!< Total volume of the real (small) voxel */
    double big_total_volume; /*!< Total volume of the (padded) blood-vessel voxel (bv_min_limits/bv_max_limits) */
    double target_axons_icvf;        

    double axons_w_myelin_icvf;        /*!< Intracellular Compartment Volume Fraction of axons with myelin */
    double axons_wo_myelin_icvf;         /*!< Intracellular Compartment Volume Fraction of axons without myelin */
    double glial_pop1_soma_icvf;        /*!< Astrocyte Intracellular Compartment Volume Fraction */
    double glial_pop1_processes_icvf;         /*!< Astrocyte Intracellular Compartment Volume Fraction */
    double glial_pop2_soma_icvf;        /*!< glial_pop2 Intracellular Compartment Volume Fraction */
    double glial_pop2_processes_icvf;    /*!< glial_pop2 Intracellular Compartment Volume Fraction */
    double glial_pop3_soma_icvf;        /*!< glial_pop3 Intracellular Compartment Volume Fraction */
    double glial_pop3_processes_icvf;    /*!< glial_pop3 Intracellular Compartment Volume Fraction */
    double axons_icvf;        /*!< Intracellular Compartment Volume Fraction of axons without myelin */
    double myelin_icvf;         /*!< Intracellular Compartment Volume Fraction of axons with myelin */
    double extracellular_icvf;        /*!< Extracellular Compartment Volume Fraction */
    double blood_vessels_icvf;    /*!< Arteriole Compartment Volume Fraction, relative to the real (small) voxel -- display only */
    double blood_vessels_icvf_big; /*!< Arteriole Compartment Volume Fraction, relative to the (padded) blood-vessel voxel -- what target_blood_vessels_icvf is actually compared against */
    double epsilon_blood_vessels;
    double mean_vessel_rad;
    double std_vessel_rad;
    double capillary_radius;     /*!< Fixed radius for every capillary sphere (no decay) */
    int max_generations;         /*!< Deepest allowed capillary branching depth from the arteriole */

    double target_blood_vessels_processes_icvf;  /*!< Target Intracellular Compartment Volume Fraction of capillaries -- interpreted relative to the (padded) blood-vessel voxel */
    double blood_vessels_processes_icvf;         /*!< Achieved Intracellular Compartment Volume Fraction of capillaries, relative to the real (small) voxel -- display only */
    double blood_vessels_processes_icvf_big;     /*!< Achieved Intracellular Compartment Volume Fraction of capillaries, relative to the (padded) blood-vessel voxel -- what target_blood_vessels_processes_icvf is actually compared against */

    double swelling_factor;

    int spheres_overlap_factor;                         /*!< Factor to divide the radii by */
    bool axon_can_shrink;               /*!< If true, the axons can shrink to allow passage between them */
    double cosPhiSquared;              /*!< Cosine of the angle between the axon and the plane squared */

    double alpha;                       /*!< Alpha coefficient of the Gamma distribution of the radii */
    double beta;                        /*!< Beta coefficient of the gamma distribution of the radii */

    int regrow_count = 0;               /*!< Number of axons to regrow */
    int regrow_thr;                     /*!< Number of regrowth batches allowed */

    Eigen::Vector3d min_limits;         /*!< Voxel min limits (if any) (bottom left corner) */
    Eigen::Vector3d max_limits;         /*!< Voxel max limits (if any) */

    Eigen::Vector3d bv_min_limits;       /*!< Blood vessel growth voxel min limits: min_limits symmetrically padded out to blood_vessels_voxel_size (or == min_limits if unset/not larger). Trunk/branch growth targets and stopping conditions use this box so vessels have room to actually reach an edge; seeding and ICVF still use min_limits/max_limits. */
    Eigen::Vector3d bv_max_limits;       /*!< Blood vessel growth voxel max limits, see bv_min_limits. */

    std::vector<int> intersecting_vessel_ids;              /*!< id of every blood vessel touching the (randomly placed) real voxel, set by PlaceSmallVoxel() */
    std::vector<Eigen::Vector3d> intersecting_vessel_seeds; /*!< seed (begin) coordinate of each vessel in intersecting_vessel_ids, same order -- for growth_info.txt */

    double min_radius;                  /*!< Minimum radius value of all axons */
    double max_radius;                  /*!< Maximum radius value of all axons */

    double beading_amplitude;              /*!< For beading: percentage of variation between the maximum radius and minimum radius in each axon. If set to 1, there is no "beading" */
    double beading_std;                    /*!< Standard deviation of the variation in beading */
    double epsilon;                     /*< Standard deviation for gaussian distribution in generation of directions to grow in */ 
    int undulation_factor;              /*!< axon ondulation factor : the number of spheres during whoch the axon grows straight before picking a direction from gaussian distribution  */
    double mean_glial_pop1_process_length;   /*!< Mean length of glial processes */
    double std_glial_pop1_process_length;    /*!< Standard deviation of glial processes */
    int glial_pop1_nbr_primary_processes;          /*!< Number of primary processes for glial cells */
    double mean_glial_pop2_process_length;   /*!< Mean length of glial processes */
    double std_glial_pop2_process_length;    /*!< Standard deviation of glial processes */
    int glial_pop2_nbr_primary_processes;          /*!< Number of primary processes for glial cells */
    double mean_glial_pop3_process_length;   /*!< Mean length of glial processes */
    double std_glial_pop3_process_length;    /*!< Standard deviation of glial processes */
    int glial_pop3_nbr_primary_processes;          /*!< Number of primary processes for glial cells */

    double c1;                              /*!< First coefficient of the Myelin thickness */
    double c2;                              /*!< Second coefficient of the Myelin thickness */
    double c3;                              /*!< Third coefficient of the Myelin thickness */
    double expanded_for_glial_space;         /*!< Volume expansion factor to account for glial space */

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
    CaterpillarGrowth() {}

    /*!
     *  \brief Initialize everything.
     */
    CaterpillarGrowth(const Parameters &params, const Eigen::Vector3d &min_l, const Eigen::Vector3d &max_l);
    
    
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
    void growthThread(std::vector <Axon> &axs, Axon &axon, AxonGrowth &growth, int &finished, int &grow_straight, int &straight_growths, double &stuck_radius, int &stuck_index);

    /*!
     *  \param axon Axon to modify
     *  \brief Causes sinusoidal fluctuation of the radii
     */
    double radiusVariation(const Axon &axon);

    /*!
       *  \param radiis List of axon radii
     *  \brief Generates list of radii following gamma distribution
     */
    void generate_radii(std::vector<double> &radiis, std::vector<bool> &has_myelin);

    /*!
     *  \param radii_ List of target axon radii
     *  \param has_myelin Parallel to radii_
     *  \param angles Parallel to radii_
     *  \brief Grows every (already seeded, see seedAllAxons) axon through the
     *         box in shallow depth-layers along each axon's own growth_axis: every
     *         still-active axon is advanced (largest radius first, i.e. existing vector
     *         order) only up to the current layer's depth cap before any axon moves to
     *         the next layer, so no axon's tortuous wandering can claim uncontested
     *         territory far ahead of axons that haven't started growing yet. Replaces
     *         the old setBatches/growBatches/growBatch nbr_threads-per-batch,
     *         full-depth-per-batch scheme.
     */
    void growAxonsLayered(std::vector<double> &radii_, std::vector<bool> &has_myelin, std::vector<double> &angles);

    /*!
     *  \param Q Starting point of axon
     *  \param D Ending point of axon
     *  \brief Gets a position in the plane for the starting point of an axon
     */
    bool get_begin_end_point(Eigen::Vector3d &Q, Eigen::Vector3d &D, double &angle);

    /*!
     *  \brief Creates the entire substrate
     */
    void createSubstrate();

    void PlaceBloodVessels();
    void GrowBloodVessels();
    


    /*!
     *  \brief Rebuilds sphere_grid from the current axons/glial_pop1/glial_pop2/blood_vessels
     *         and checks that no sphere collides with a sphere of another object. Meant to be
     *         run once, at the very end of the simulation; uses the grid so it stays fast even
     *         with many cells.
     */
    bool checkNoCollisions();

    /*!
     *  \param radii_ List of target axon radii (in/out: entries may be overwritten with a
     *         resampled replacement radius if the original draw never found a valid spot,
     *         and trimmed only for the rare slot that fails even after max_resamples
     *         redraws)
     *  \param new_axons Out: one Axon per successfully seeded axon, each at its actually
     *         placed radius
     *  \param has_myelin In/out, kept parallel to radii_
     *  \param angles In/out, kept parallel to radii_
     *  \brief Places every axon's seed sphere directly at its full target radius, one at a
     *         time in largest-first order (generate_radii's existing sort). Each axon's
     *         search budget (number of random position tries) scales with
     *         (mean_radius/target_radius)^2, so smaller-than-average axons get
     *         proportionally more tries than a fixed budget would give them. If an axon
     *         still can't find a valid spot, its radius is redrawn from the same
     *         Gamma(alpha,beta) distribution (not discarded, and the whole pass is not
     *         restarted) and placement retried, up to max_resamples times, so the realized
     *         population stays an honest i.i.d. sample from the intended distribution
     *         instead of systematically losing whichever sizes are hardest to place late
     *         in a crowded pass. No separate grow-in-place or jostle step is needed
     *         afterward since every placed seed is already at its true (possibly
     *         resampled) target size.
     */
    void seedAllAxons(std::vector<double> &radii_, std::vector<Axon> &new_axons, std::vector<bool> &has_myelin, std::vector<double> &angles);

    /*!
     *  \param growing_axons Axons that have grown in this layer-pass
     *  \param stuck_radii_ radii of axons that got stuck
     *  \param stuck_indices_ indices of axons that got stuck
     *  \param layer_start_spheres Per-axon sphere count as of the start of this layer
     *         (parallel to growing_axons); on collision, rolls the offending axon back
     *         to this count (its state as of the end of the previous layer) rather
     *         than destroying it outright, and the scratch collision grid is built only
     *         from spheres added since then, since everything before was already
     *         checked and committed in an earlier layer.
     *  \brief Checks if the axons that have grown in this layer-pass collide with the
     *         committed environment or with each other.
     */
    bool SanityCheck(std::vector<Axon>& growing_axons, std::vector<double> &stuck_radii_, std::vector<int> &stuck_indices_, const std::vector<std::size_t> &layer_start_spheres);

    /*!
     *  \param axon_to_grow Axon that is grown
     *  \param stuck_radii_ radii of axons that got stuck
     *  \param stuck_indices_ indices of axons that got stuck
     *  \param layer_depth_cap How far (along this axon's own growth_axis) it may grow
     *         during this layer-pass; capped at this axon's true wall if that's nearer
     *  \param layer_start_spheres This axon's outer_spheres.size() at the start of this
     *         layer-pass, so a rollback on collision/exhausted-retries undoes only
     *         spheres added during this layer
     *  \brief Grows one axon up to this layer's depth cap (or its true wall, whichever
     *         is nearer).
     */
    void growAxon(Axon& axon_to_grow, int& index, double& stuck_radius, int& stuck_index, double layer_depth_cap, std::size_t layer_start_spheres, bool can_shrink_this_round);

    /*!
     *  \param ax Axon to check
     *  \brief True if ax's last sphere has reached the box's true wall along its own
     *         growth_axis (not just an intermediate layer's depth cap) -- mirrors
     *         AddOneSphere's own wall checks, but always against the true max_limits.
     */
    bool hasReachedTrueWall(const Axon& ax) const;

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
    bool PlaceAxon(const int &axon_id, const double &seed_radius, const double &growth_radius, const Eigen::Vector3d &Q, const Eigen::Vector3d &D, std::vector<Axon> &new_axons, const bool &has_myelin, const double &angle, const bool &outside_voxel);

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
     *  \brief Gradually, non-uniformly swells every axon's spheres back
     *         toward their own true target radius (cap = radius/swelling_factor,
     *         or itself for the seed sphere), each as far as its own local
     *         room allows, and updates the ICVF.
    */

    void SwellAxons();

    /*!
     *  \brief Swells all spheres in one axon by up to `percentage` this round, each capped at its own target.
     *         Computes every sphere's new radius in parallel via pool (safe: spheres in the
     *         same axon never collide-check against each other), then applies sphere_grid
     *         updates sequentially.
         \param ax Axon to swell
         \param percentage this round's requested growth step, as a fraction of each sphere's current radius
         \param caps per-sphere target radius cap, same size as ax.outer_spheres
         \param pool thread pool shared across the whole swelling pass
     \return true if any sphere in this axon actually grew this round (so the caller
             knows whether update_Volume/ICVF need to be recomputed at all)
    */

    bool SwellAxon(Axon &ax, const double &percentage, const std::vector<double> &caps, ThreadPool &pool);

    /*!
     *  \brief Read-only: computes the largest radius (up to cap_radius) sph could grow to
     *         this round (offered up to `percentage`, clipped to whatever fits via bisection).
     *         Does not mutate sph or sphere_grid -- safe to call concurrently across spheres
     *         that don't collide-check against each other (e.g. spheres of the same axon).
    */

    double ComputeSwollenRadius(const Sphere &sph, const double &percentage, const double &cap_radius) const;

    /*!
     *  \brief Push-assisted fallback swelling for one axon, only used once plain in-place
     *         swelling (SwellAxon) has already converged short of target_axons_icvf. Same
     *         parallel-compute / sequential-apply structure as SwellAxon.
     *  \return true if any sphere in this axon actually changed (grew and/or moved) this round
    */

    bool SwellAxonWithPush(Axon &ax, const double &percentage, const std::vector<double> &caps, ThreadPool &pool);

    /*!
     *  \brief Read-only: like ComputeSwollenRadius, but if growing sph in place can't reach
     *         this round's target, also tries pushing it away from whichever neighbor is
     *         blocking it (see FindSwellPush), capped so the resulting distance to that
     *         neighbor never exceeds max(sph.radius, blocker_radius)/swelling_factor.
     *         Returns the best (center, radius) found; never mutates sph or sphere_grid.
    */

    Sphere ComputeSwollenSphere(const Sphere &sph, const double &percentage, const double &cap_radius, int growth_axis,
                                 const Sphere *prev_sphere = nullptr, const Sphere *next_sphere = nullptr) const;

    /*!
     *  \brief Read-only: finds the neighbor causing the deepest overlap with sph at
     *         trial_radius and returns a push vector (perpendicular to growth_axis) away
     *         from it, plus that neighbor's radius. Mirrors AddOneSphere's findPush.
     *         Returns false if nothing overlaps or there's no lateral escape direction.
    */

    bool FindSwellPush(const Sphere &sph, const double &trial_radius, int growth_axis, Eigen::Vector3d &push, double &blocker_radius) const;

    /*!
     *  \brief Finds a random point on the opposite plane based on the current point and cosPhiSquared (c2)
        \param x1 x coordinate of the point
        \param y1 y coordinate of the point
        \param L distance between the two planes
     */
    Eigen::Vector3d randomPointOnPlane(const Eigen::Vector3d &begin, const Eigen::Vector3d &end, const int &axis1, const int &axis2, const int &axis3, double &angle, bool &outside_voxel);

    /*!
     *  \brief Inverts myelin_thickness (outerRadius = innerRadius +
     *         myelin_thickness(innerRadius)) for a single outer radius. Uses
     *         inner_radius_lut when it's been built and outerRadius falls
     *         inside its range (the common case: add_Myelin calls this once
     *         per sphere of every myelinated axon, so a precomputed table
     *         turns what would be a full Newton's-method solve -- several
     *         iterations, two function evaluations each -- into one O(1)
     *         interpolated lookup), falling back to the exact Newton's-method
     *         solve otherwise (buildInnerRadiusLUT itself uses this fallback
     *         path to populate the table in the first place, and any query
     *         outside the table's precomputed range also lands here rather
     *         than extrapolating).
     */
    double findInnerRadius(const double &outerRadius);

    /*!
     *  \brief Precomputes inner_radius_lut: findInnerRadius's exact
     *         Newton's-method answer at n_points uniformly-spaced outer radii
     *         across [lo, hi], linearly interpolated between by
     *         findInnerRadius afterward. Meant to be called once, after the
     *         true outer-radius range in play is known (e.g. right after
     *         generate_radii sorts/sets max_radius), before any of the many
     *         per-sphere findInnerRadius calls in add_Myelin(). A no-op if
     *         hi <= lo (e.g. no myelinated axons were actually placed).
     */
    void buildInnerRadiusLUT(double lo, double hi, int n_points = 2000);

    std::vector<double> inner_radius_lut_outer; /*!< uniformly-spaced outer radii sample points, ascending */
    std::vector<double> inner_radius_lut_inner; /*!< findInnerRadius's exact answer at each inner_radius_lut_outer entry, same order */

    void growBranches(const int &population_nbr);

    /*!
     *  \brief Grows blood vessel branches (branch_id >= 1) off the already-grown main
     *         vessels, until target_blood_vessels_processes_icvf is reached. Mirrors
     *         growBranches, but there is a single population and no primary-branch phase:
     *         every vessel starts with only its main vessel (branch 0) to branch off of.
     */
    void GrowBloodVesselBranches();

    /*!
     *  \brief Enforces Murray's law at every blood vessel branching point, once all
     *         vessels (main vessel and every branch) are fully grown. See
     *         Blood_Vessel::enforceMurraysLaw for the per-vessel algorithm.
     */
    void ApplyMurraysLawToBloodVessels();

    /*!
     *  \brief Randomly repositions the real (small) voxel -- min_limits/max_limits,
     *         same edge length as before, translated -- somewhere inside the
     *         (padded) blood-vessel voxel, subject to actually intersecting at
     *         least one already-grown blood vessel (arteriole and/or capillaries).
     *         Meant to run once, after blood vessels are fully grown (including
     *         Murray's law thinning), and before any other population is placed --
     *         everything downstream (glial cells, axons, ICVF, output) uses
     *         min_limits/max_limits, so this must finalize their position first.
     *         Picks a uniformly random sphere among every grown vessel sphere,
     *         then a random small-voxel placement guaranteed to contain that
     *         sphere's center (clamped to stay inside the big voxel), rather than
     *         blindly retrying random placements until one happens to intersect --
     *         with a sparse vessel network in a large big voxel, blind placement
     *         could need many tries (or fail outright) to land near a vessel at
     *         all. If there are no blood vessels to intersect, min_limits/max_limits
     *         are left at their default (origin-anchored) position instead, since
     *         the intersection requirement would otherwise be unsatisfiable.
     *         Populates intersecting_vessel_ids/intersecting_vessel_seeds with
     *         every vessel that ends up touching the committed placement (not
     *         just the one used to seed the search), for growth_info.txt.
     */
    void PlaceSmallVoxel();

    void ICVF(const std::vector<Axon> &axs, const std::vector<Glial> &glial_pop1, const std::vector<Glial> &oligos, const std::vector<Glial> &glial_pop3, const std::vector<Blood_Vessel> &blood_vessels);
    
    double c2toKappa(double c2_target, double c2_tol, double kappa_max);
    std::vector<double> generate_angles(const int &num_samples);
    void processBatchWithThreadPool(std::vector<Axon>& axons_to_grow, std::vector<int> &indices, std::vector<double>& stuck_radii, std::vector<int>& stuck_indices, double layer_depth_cap, const std::vector<std::size_t> &layer_start_spheres, bool can_shrink_this_round);
    bool check_borders(const Eigen::Vector3d&  min_l, const Eigen::Vector3d&  max_l, const Eigen::Vector3d& pos, const double& distance_to_border);
    double originalFunction(const double &x, const double &outerRadius);
    double derivative(const double &x);
    double myelin_thickness(const double &inner_radius);

};


#endif // CaterpillarGrowth_H

