//!  Sphere spatial grid ========================================================================/
/*!
*   \details   Uniform spatial hash grid over the simulation volume. Each voxel stores
*              the sphere id, cell id and object id of every sphere whose bounding box
*              overlaps it, so that collision checks (see canSpherebePlaced) can query a
*              small neighborhood instead of scanning every axon/glial cell/blood vessel.
*   \author    Jasmine Nguyen-Duc
*   \date      July 2026
*   \version   1.0
=================================================================================================*/

#ifndef SPHERE_GRID_H
#define SPHERE_GRID_H

#include "Eigen/Dense"
#include "sphere.h"
#include <vector>

/// @brief One sphere's identifying info and geometry, as stored in a grid voxel.
struct SphereGridEntry
{
    int sphere_id;           /*!< Sphere::id, the id of the sphere itself                              */
    int cell_id;             /*!< Sphere::object_id, the id of the axon/glial cell/blood vessel it belongs to */
    int object_id;           /*!< Sphere::object_type, the kind of object it belongs to (see sphere.h)  */
    Eigen::Vector3d center;  /*!< Sphere::center, kept here so collision checks don't need to look the sphere up elsewhere */
    double radius;           /*!< Sphere::radius     */
    int branch_id;                                                  
};

/// @brief Uniform grid used to accelerate spatial lookups of spheres.
class SphereGrid
{
public:

    SphereGrid();

    /*!
     *  \param min_limits_ lower corner of the volume covered by the grid
     *  \param max_limits_ upper corner of the volume covered by the grid
     *  \param voxel_size_ edge length of a (cubic) voxel
     */
    SphereGrid(const Eigen::Vector3d &min_limits_, const Eigen::Vector3d &max_limits_, const double &voxel_size_);

    /*!
     *  \brief Empties every voxel, keeping the grid's bounds and voxel size.
     */
    void clear();

    /*!
     *  \brief Adds a sphere to every voxel its bounding box overlaps.
     */
    void insert(const Sphere &sph);

    /*!
     *  \brief Removes a sphere's entries (matched by sphere_id/cell_id/object_id) from
     *         every voxel its bounding box overlaps. Pass the sphere's geometry as it
     *         was when inserted (e.g. before swelling changes its radius), so the same
     *         voxel range is recomputed and the stale entries are actually found.
     */
    void remove(const Sphere &sph);

    /*!
     *  \brief Returns the entries of every voxel overlapping the given query sphere.
     *         Intended as the candidate list for a subsequent narrow-phase check;
     *         the same sphere may appear more than once if it spans several voxels.
     */
    std::vector<SphereGridEntry> query(const Eigen::Vector3d &center, const double &radius) const;

    /*!
     *  \brief Checks whether a sphere can be placed without overlapping any sphere
     *         already in the grid, other than spheres belonging to the same object
     *         (matched by object_type + object_id, i.e. Sphere::object_type/object_id).
     *  \param sph the candidate sphere
     *  \param tolerance extra clearance required between surfaces (e.g. barrier_tickness)
     */
    bool canSpherebePlaced(const Sphere &sph, bool check_collision_with_branches = true) const;

    /*!
     *  \brief Finds the neighbor causing the deepest overlap with a candidate sphere of
     *         the given center/radius, excluding spheres belonging to the same object
     *         (object_type/object_id, matching canSpherebePlaced's exclusion). Iterates
     *         voxels/entries directly instead of going through query(): unlike
     *         canSpherebePlaced this can't early-exit (every candidate within
     *         search_radius must be compared to find the worst one), but it still
     *         avoids query()'s heap-allocated intermediate vector.
     *  \param center candidate sphere's center
     *  \param radius candidate sphere's radius, used for the overlap calculation
     *  \param search_radius how far out to scan for candidates (independent of radius,
     *         so callers can search wider than the candidate's own size when needed)
     *  \return true if an external neighbor overlaps; blocker_center/blocker_radius are
     *          only set when true.
     */
    bool findWorstOverlap(const Eigen::Vector3d &center, double radius, double search_radius,
                          int self_object_type, int self_object_id,
                          Eigen::Vector3d &blocker_center, double &blocker_radius) const;

private:

    Eigen::Vector3d min_limits;
    Eigen::Vector3d max_limits;
    double voxel_size;
    int nx, ny, nz;

    std::vector<std::vector<SphereGridEntry>> voxels;

    Eigen::Vector3i clampedVoxelIndex(const Eigen::Vector3d &position) const;
    int linearIndex(int ix, int iy, int iz) const;
};

#endif // SPHERE_GRID_H
