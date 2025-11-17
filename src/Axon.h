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
    Eigen::Vector3d begin;                          /*!< position of first sphere */
    Eigen::Vector3d end;                            /*!< target position to grow towards */
    std::vector<Eigen::Vector2d> Box;               /*!< Box with <min, max> for each axis (x,y,z) */
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
        inner_spheres.clear();
        outer_spheres.clear();
        Box.clear();
        growth_attempts = 0;
        beading_amplitude = beading_amplitude_;
        beading_std = beading_std_;
        undulation_factor = undulation_factor_;
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
            Box = ax.Box;
            growth_attempts = ax.growth_attempts;
            beading_amplitude = ax.beading_amplitude;
            undulation_factor = ax.undulation_factor;
            beading_std = ax.beading_std;
            phase_shift = ax.phase_shift;
            myelin_sheath = ax.myelin_sheath;
            volume = ax.volume;
            volume_myelin = ax.volume_myelin;
            growth_axis = ax.growth_axis;
            angle = ax.angle;
            outside_voxel = ax.outside_voxel;
            inner_radius = ax.inner_radius;

        }
        return *this;
    }

    void keep_one_sphere();
    /*!
     *  \param sphere_to_add sphere to add
     *  \brief Adds sphere to axon
     */
    void add_sphere(const Sphere &sphere_to_add);

    /*!
     *  \param sph sphere to check
     *  \brief Checks if a sphere collides with this axon
     */
    bool isSphereInsideAxon(const Sphere &sph) const;

    /*!
     *  \param sph sphere to check
        \param axis axis (0, 1 or 2)
     *  \brief Checks if a sphere collides with this axon along a specific axis
     */
    std::vector<int> checkAxisForCollision(const Sphere &sph, const int &axis) const;

    /*!
     *  \param position position in voxel
        \param distance_to_be_inside distance to be inside Box
     *  \brief Checks if a position is near this axon (inside the Box)
     */
    bool isNearAxon(const Eigen::Vector3d &position, const double &distance_to_be_inside) const;


    /*!
     *  \brief Deletes all spheres in axon.
     */
    void destroy();

    bool isSphereInsideInnerAxon(const Sphere &sph) const;

    std::vector<int> checkAxisForInnerCollision(const Sphere &sph, const int &axis) const;

    void updateBox();

    void update_Volume(const int &factor, const Eigen::Vector3d &min_limits, const Eigen::Vector3d &max_limits);

};

#endif // AXON_H
