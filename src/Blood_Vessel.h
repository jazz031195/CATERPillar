//!  Axon Obstacle Derived Class =============================================================/
/*!
*   \details   Axon class derived from an Obstacle
*              in the direction set by begin, end.
*   \author    Jasmine Nguyen-Duc
*   \date      September 2023
*   \version   1.42
=================================================================================================*/

#ifndef BLOOD_VESSEL_H
#define BLOOD_VESSEL_H

#include "sphere.h"
#include "SphereGrid.h"
#include <vector>

using namespace std;

/// @brief
class Blood_Vessel : public Obstacle
{
public:
    int id;                     /*!< ID of the blood vessel */
    double radius;              /*!< Radius of the blood vessel */
    std::vector<std::vector<Sphere>> ramification_spheres; /*!< spheres of the vessel, grouped by branch_id; branch 0 is the main vessel */
    std::vector<std::vector<double>> lengths_branches;      /*!< cumulative length at each sphere of each branch */
    std::vector<Eigen::Vector3d> attractors;                /*!< per-branch target direction (unused for branch 0) */
    std::vector<std::vector<int>> children_branches;        /*!< children_branches[b] = ids of every branch that attached to branch b */
    std::vector<int> branch_generation;                     /*!< branch_generation[b] = generations of capillary branching between branch b and the arteriole (0 for the arteriole itself, 1 for a capillary attached directly to it, etc.) */
    double minimum_radius;                                  /*!< radius branches decay toward */
    double volume;                                          /*!< Volume of the main vessel (branch 0), clipped to the real (small) voxel */
    double volume_processes;                                /*!< Volume of the branches (branch_id >= 1), clipped to the real (small) voxel */
    double volume_big;                                       /*!< Volume of the main vessel (branch 0), clipped to the (padded) blood-vessel voxel */
    double volume_processes_big;                             /*!< Volume of the branches (branch_id >= 1), clipped to the (padded) blood-vessel voxel */
    Eigen::Vector3d begin;                          /*!< position of first sphere */
    Eigen::Vector3d end;                            /*!< target position to grow towards */
    int undulation_factor;                          /*!< Factor for ondulation */
    int growth_attempts;                            /*!< Number of attempts to grow axon in a row*/
    double beading_amplitude;                     /*!< Amplitude of radius beading */
    double phase_shift;                           /*!< Phase shift of radius beading */
    double beading_std;                           /*!< Standard deviation of radius beading */
    int growth_axis ;                            /*!< Axis along which the blood vessel grows (0:x, 1:y, 2:z) */
    int next_sphere_id;                           /*!< Monotonic counter handing out the id of the next sphere added anywhere in this vessel, so ids stay unique within the vessel regardless of overlap factor -- parent_id links and branch-parent lookups rely on that uniqueness. */

    Blood_Vessel();

    ~Blood_Vessel();

    Blood_Vessel(const int &id_, const Eigen::Vector3d &begin_, const Eigen::Vector3d &end_, const double &radius_, const double &beading_amplitude_, const double &beading_std_, const double &undulation_factor_){
        id = id_;
        begin = begin_;
        end = end_;
        radius = radius_;
        growth_attempts = 0;
        beading_amplitude = beading_amplitude_;
        beading_std = beading_std_;
        undulation_factor = undulation_factor_;
        growth_axis= 2;
        minimum_radius = 2.0;
        volume = 0.0;
        volume_processes = 0.0;
        volume_big = 0.0;
        volume_processes_big = 0.0;
        next_sphere_id = 0;
        ramification_spheres.clear();
        lengths_branches.clear();
        attractors.clear();
        children_branches.clear();
        branch_generation.clear();
    };

    Blood_Vessel& operator=(const Blood_Vessel &bv){
        if (this != &bv){
            id = bv.id;
            begin = bv.begin;
            end = bv.end;
            radius = bv.radius;
            ramification_spheres = bv.ramification_spheres;
            lengths_branches = bv.lengths_branches;
            attractors = bv.attractors;
            children_branches = bv.children_branches;
            branch_generation = bv.branch_generation;
            minimum_radius = bv.minimum_radius;
            undulation_factor = bv.undulation_factor;
            growth_attempts = bv.growth_attempts;
            beading_amplitude = bv.beading_amplitude;
            beading_std = bv.beading_std;
            growth_axis= bv.growth_axis;
            volume = bv.volume;
            volume_processes = bv.volume_processes;
            volume_big = bv.volume_big;
            volume_processes_big = bv.volume_processes_big;
            next_sphere_id = bv.next_sphere_id;

        }
        return *this;
    };


    void destroy();
    void keep_one_sphere();
    void add_sphere(const Sphere &sphere_to_add);
    void update_Volume(const int &factor, const Eigen::Vector3d &min_limits, const Eigen::Vector3d &max_limits, const Eigen::Vector3d &big_min_limits, const Eigen::Vector3d &big_max_limits);
    void add_first_sphere(const Sphere &s);

    /*!
     *  \brief Adds every sphere of this blood vessel (main vessel and branches) to the grid.
     */
    void addToGrid(SphereGrid &grid) const;

    /*!
     *  \brief Computes volume_processes from every branch (branch_id >= 1).
     */
    void compute_processes_icvf(const int &factor, const Eigen::Vector3d &min_limits, const Eigen::Vector3d &max_limits, const Eigen::Vector3d &big_min_limits, const Eigen::Vector3d &big_max_limits);

    /*!
     *  \brief Finds the sphere with the given id anywhere in this vessel (main vessel or
     *         any branch). Returns false if no such sphere exists.
     */
    bool findSphereById(const int &sphere_id, int &out_branch, int &out_index) const;

    /*!
     *  \brief Enforces Murray's law (R_parent^3 = r_continuation^3 + r_child^3) at every
     *         branching point, once this vessel (main vessel and all its branches) is
     *         fully grown -- for every parent branch, the arteriole (branch 0) and every
     *         capillary alike, not just junctions directly off the arteriole. Each
     *         branch's own radius is otherwise fixed (no decay along its length),
     *         assigned once at growth time from its generation (see
     *         BloodVesselGrowth::growBranch, r(g) = r0 * 2^(-g/gamma)) -- what this
     *         corrects is solely the *parent's* downstream continuation past each
     *         junction, solved for so its cube plus the (unchanged) child's cube equals
     *         the parent's own radius cube at that point. Junctions are processed from
     *         the trunk outward (branch_id increasing) and, within a parent, most
     *         upstream junction first, so a correction cascades correctly to junctions
     *         further downstream on the same branch (including a parent's later
     *         siblings, which see an already-reduced local budget from earlier ones).
     *         Before committing a junction's split, both sides are predicted: the
     *         child's own radius, and the parent's own downstream continuation (scanned
     *         for its true minimum, since it can "bead" and isn't necessarily
     *         monotonic). If either would drop below minimum_radius, the split would only
     *         end up truncating something anyway (likely losing the voxel-plane crossing a
     *         branch was grown to reach) -- so the child is discarded outright instead and
     *         the split is skipped entirely, leaving the parent's radius untouched, as
     *         though the branch had never been created, rather than needlessly thinning a
     *         parent for a branch that doesn't survive. This is also what naturally caps
     *         how many children a single parent can support before its local
     *         cross-section budget runs out.
     */
    void enforceMurraysLaw();

    /*!
     *  \brief Deletes (clears) every branch whose starting radius has dropped below
     *         minimum_radius -- and, since a branch can't meaningfully exist without its
     *         parent, its entire subtree too. A branch may also be truncated (rather than
     *         fully deleted) if the violation occurs partway along its length; but a
     *         branch is only kept if, after any such truncation, its last sphere still
     *         reaches outside [min_limits, max_limits] -- branches only exist to cross a
     *         voxel plane, so a stub left dead-ending inside the domain is deleted
     *         outright too. enforceMurraysLaw() already discards a branch upfront rather
     *         than shrink it below minimum_radius, so this is mainly a backstop for cases
     *         it doesn't cover (e.g. a branch's own continuation past a nested junction).
     *         Meant to be run after enforceMurraysLaw() and before
     *         reinterpolateAfterShrink().
     */
    void pruneUndersizedBranches(const Eigen::Vector3d &min_limits, const Eigen::Vector3d &max_limits);

    /*!
     *  \brief Clears branch (and every branch in its subtree, recursively). Does not
     *         shift any other branch's index, so parent_id / children_branches
     *         elsewhere stay valid.
     */
    void deleteSubtree(int branch);

    /*!
     *  \brief Re-densifies every branch (and the main vessel) so no two consecutive
     *         spheres are farther apart than r_max/factor, r_max being the larger of
     *         the two radii. Meant to be run after enforceMurraysLaw() (and after
     *         pruneUndersizedBranches()): shrinking a segment's radius does not touch
     *         sphere spacing, so a chain that was adequately overlapping at its
     *         original (larger) radius can end up with visible gaps once thinned out;
     *         this reinserts spheres to close them.
     *  \param sphere_grid the shared environment grid (already holds every vessel's
     *         growth-time geometry by this point): candidate spheres are checked
     *         against it, not just this vessel's own other branches, so a re-densified
     *         point can't silently land inside another vessel it happens to pass near.
     */
    void reinterpolateAfterShrink(const int &factor, const SphereGrid &sphere_grid);

    /*!
     *  \brief Closes the gap, if any, between each branch's first sphere and the
     *         parent sphere it emerged from. enforceMurraysLaw only rescales radii
     *         (sphere centers never move), so a branch's attachment point can end
     *         up with a smaller combined radius than the (fixed, original) distance
     *         between it and its parent sphere once upstream junctions have
     *         cascaded several shrinks onto that parent sphere -- reinterpolateAfterShrink
     *         cannot catch this since it only walks spheres within a single branch's
     *         own vector, never across the parent/child boundary. Meant to run after
     *         enforceMurraysLaw() and pruneUndersizedBranches(), before
     *         reinterpolateAfterShrink().
     *  \param sphere_grid see reinterpolateAfterShrink -- same reasoning.
     */
    void bridgeJunctionGaps(const int &factor, const SphereGrid &sphere_grid);

};
#endif // BLOOD_VESSEL_H
