//!  Axon Obstacle Derived Class =============================================================/
/*!
*   \details   Axon class derived from an Obstacle
*              in the direction set by begin, end.
*   \author    Jasmine Nguyen-Duc
*   \date      September 2023
*   \version   1.42
=================================================================================================*/

#ifndef AXON_H
#define AXON_H

#include "sphere.h"
#include "SphereGrid.h"
#include <vector>

using namespace std;

/// @brief
class Axon : public Obstacle
{
public:


    int id;                                         /*!< ID of axon */
    std::vector<Sphere> outer_spheres;              /*!< outer spheres in axon */
    std::vector<Sphere> inner_spheres;              /*!< inner spheres in axon */
    double radius;                                  /*!< radius of axon */
    double inner_radius;                            /*!< inner radius of axon */
    double target_radius = 0.0;                     /*!< True, un-swelling_factor-shrunk target radius this axon aims to reach when swelling (set once in PlaceAxon from seed_radius) -- distinct from radius, which can be smaller (growth_radius) while the axon is still thin and growing. Used to cap the axon's own mean radius after swelling (see SwellAxons). */
    Eigen::Vector3d begin;                          /*!< position of first sphere */
    Eigen::Vector3d end;                            /*!< target position to grow towards */
    int growth_attempts;                            /*!< Number of attempts to grow axon in a row*/
    double beading_amplitude;                       /*!< Amplitude of beading */
    double beading_std;                          /*!< Standard deviation of beading */
    double phase_shift;                             /*!< Phase shift of beading */
    bool myelin_sheath;                                /*!< Axon has myelin */
    double volume;                                  /*!< Volume of axon */
    double volume_myelin;                           /*!< Volume of myelin */
    int growth_axis;                                /*!< Axis to grow along */
    double angle;
    bool outside_voxel;
    bool has_shrunk;                            /*!< Axon has shrunk */
    int undulation_factor;                          /*!< Factor for ondulation */
    int grow_straight;                              /*!< Whether the next sphere grows straight (1) or takes a fresh epsilon-random direction (0); persists across depth layers */
    int straight_growths;                           /*!< Consecutive straight growths so far in the current straight streak; persists across depth layers */
    Eigen::Vector3d turn_start_direction;            /*!< Heading at the start of the current turn cycle, interpolated away from as the cycle progresses */
    Eigen::Vector3d turn_target_direction;           /*!< Freshly-sampled heading this turn cycle is gradually curving toward, reached by the cycle's last straight step */

    /*!
     *  \brief Default constructor. Does nothing
     */
    Axon();

    ~Axon();

    Axon(const int &id_, const Eigen::Vector3d &begin_, const Eigen::Vector3d &end_, const double &radius_, const double &beading_amplitude_, const double &beading_std_, const double &undulation_factor_, const bool &myelin_sheath_, const double& angle_, const bool &outside_voxel_ = false)
    {
        id = id_;
        begin = begin_;
        end = end_;
        radius = radius_;
        inner_radius = radius_;
        target_radius = radius_; // default; PlaceAxon overrides with the true un-shrunk seed_radius
        inner_spheres.clear();
        outer_spheres.clear();
        growth_attempts = 0;
        beading_amplitude = beading_amplitude_;
        beading_std = beading_std_;
        undulation_factor = undulation_factor_;
        grow_straight = 0;
        straight_growths = 0;
        turn_start_direction = Eigen::Vector3d::Zero();
        turn_target_direction = Eigen::Vector3d::Zero();
        myelin_sheath = myelin_sheath_;
        // random phase shift between 0 and 0.5
        phase_shift = abs((double)rand() / (RAND_MAX + 1.0))/2;
        volume = 0.0;
        volume_myelin = 0.0;
        double maxVal = begin.cwiseAbs().minCoeff(&this->growth_axis); // gives growth_axis the position of the maximum value
        angle = angle_;
        outside_voxel = outside_voxel_;

    };


    Axon& operator=(const Axon& ax) {
        if (this != &ax) {
            id = ax.id;
            inner_spheres = ax.inner_spheres;
            outer_spheres = ax.outer_spheres;
            radius = ax.radius;
            begin = ax.begin;
            end = ax.end;
            growth_attempts = ax.growth_attempts;
            beading_amplitude = ax.beading_amplitude;
            undulation_factor = ax.undulation_factor;
            grow_straight = ax.grow_straight;
            straight_growths = ax.straight_growths;
            turn_start_direction = ax.turn_start_direction;
            turn_target_direction = ax.turn_target_direction;
            beading_std = ax.beading_std;
            phase_shift = ax.phase_shift;
            myelin_sheath = ax.myelin_sheath;
            volume = ax.volume;
            volume_myelin = ax.volume_myelin;
            growth_axis = ax.growth_axis;
            angle = ax.angle;
            outside_voxel = ax.outside_voxel;
            inner_radius = ax.inner_radius;
            target_radius = ax.target_radius;

        }
        return *this;
    }

    void keep_one_sphere();

    /*!
     *  \param n Number of leading spheres to keep
     *  \brief Rolls outer_spheres back to its first n elements (no-op if already <= n).
     *         Generalizes keep_one_sphere() (equivalent to truncate_to(1)) so a rollback
     *         can target any earlier point in the axon's growth, not just the very
     *         start -- needed so a mid-layer collision (depth-layered growth) only
     *         undoes spheres added during the current layer, not previously-committed
     *         layers already merged into the shared sphere grid.
     */
    void truncate_to(std::size_t n);
    /*!
     *  \param sphere_to_add sphere to add
     *  \brief Adds sphere to axon
     */
    void add_sphere(const Sphere &sphere_to_add);

    /*!
     *  \brief Adds every sphere of this axon (outer and inner) to the grid.
     */
    void addToGrid(SphereGrid &grid) const;

    /*!
     *  \brief Deletes all spheres in axon.
     */
    void destroy();

    /*!
     *  \param removed appended with every sphere dropped from outer_spheres
     *  \brief Post-swelling cleanup: a sphere fully contained inside a
     *         neighboring sphere in the chain contributes no exposed surface
     *         of its own (any obstacle it could present, a bigger neighbor
     *         already presents), so it is redundant once independent
     *         per-sphere swelling has let radii diverge enough for that to
     *         happen. Walks the chain once (monotonic-stack style) so a run
     *         of several mutually-engulfing spheres collapses correctly, not
     *         just adjacent pairs. Always keeps at least the first sphere.
     *         Caller is responsible for removing the returned spheres from
     *         the shared SphereGrid and recomputing volume/ICVF afterward.
     */
    void removeEngulfedSpheres(std::vector<Sphere> &removed);

    void update_Volume(const int &factor, const Eigen::Vector3d &min_limits, const Eigen::Vector3d &max_limits);

    /*!
     *  \brief Computes volume_myelin: the myelin sheath's own volume alone
     *         (outer boundary minus the bare axon core), by summing the same
     *         truncated-cone volume as update_Volume but over inner_spheres,
     *         then subtracting from volume (the outer boundary's volume --
     *         must already be up to date, e.g. via a prior update_Volume
     *         call). 0 for non-myelinated axons or before add_Myelin has
     *         populated inner_spheres.
     */
    void update_MyelinVolume(const Eigen::Vector3d &min_limits, const Eigen::Vector3d &max_limits);

};

#endif // AXON_H
