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
#include <vector>

using namespace std;

/// @brief
class Blood_Vessel : public Obstacle
{
public:
    int id;                     /*!< ID of the blood vessel */
    double radius;              /*!< Radius of the blood vessel */
    std::vector<Sphere> spheres;              /*!< outer spheres in axon */
    Eigen::Vector3d begin;                          /*!< position of first sphere */
    Eigen::Vector3d end;                            /*!< target position to grow towards */
    int undulation_factor;                          /*!< Factor for ondulation */
    int growth_attempts;                            /*!< Number of attempts to grow axon in a row*/
    std::vector<Eigen::Vector2d> Box;               /*!< Box with <min, max> for each axis (x,y,z) */
    double beading_amplitude;                     /*!< Amplitude of radius beading */
    double phase_shift;                           /*!< Phase shift of radius beading */
    double beading_std;                           /*!< Standard deviation of radius beading */
    double volume;                             /*!< Volume of the blood vessel */
    int growth_axis ;                            /*!< Axis along which the blood vessel grows (0:x, 1:y, 2:z) */


    Blood_Vessel();

    ~Blood_Vessel();

    Blood_Vessel(const int &id_, const Eigen::Vector3d &begin_, const Eigen::Vector3d &end_, const double &radius_, const double &beading_amplitude_, const double &beading_std_, const double &undulation_factor_){
        id = id_;
        begin = begin_;
        end = end_;
        radius = radius_;
        undulation_factor = 0;
        growth_attempts = 0;
        beading_amplitude = beading_amplitude_;
        beading_std = beading_std_;
        undulation_factor = undulation_factor_;
        growth_axis= 2;
        Box.clear();
    };

    Blood_Vessel& operator=(const Blood_Vessel &bv){
        if (this != &bv){
            id = bv.id;
            begin = bv.begin;
            end = bv.end;
            radius = bv.radius;
            spheres = bv.spheres;
            undulation_factor = bv.undulation_factor;
            growth_attempts = bv.growth_attempts;
            Box = bv.Box;
            beading_amplitude = bv.beading_amplitude;
            beading_std = bv.beading_std;
            undulation_factor = bv.undulation_factor;
            growth_axis= bv.growth_axis;
        }
        return *this;
    };
    
    
    void updateBox();
    void destroy();
    void keep_one_sphere();
    void add_sphere(const Sphere &sphere_to_add);
    bool isNearBlood_Vessel(const Eigen::Vector3d &position, const double &distance_to_be_inside) const;
    bool isSphereInsideBlood_Vessel(const Sphere &sph) const;
    std::vector<int> checkAxisForCollision(const Sphere &sph, const int &axis) const;
    void update_Volume(const int &factor, const Eigen::Vector3d &min_limits, const Eigen::Vector3d &max_limits);
    void add_first_sphere(const Sphere &s);

};
#endif // BLOOD_VESSEL_H







