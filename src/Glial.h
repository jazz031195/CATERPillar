//!  Glial Obstacle Derived Class =============================================================/
/*!
*   \details   Glial class derived from an Obstacle
*              in the direction set by begin, end.
*   \author    Jasmine Nguyen-Duc
*   \date      February 2024
*   \version   1.42
=================================================================================================*/

#ifndef GLIAL_H
#define GLIAL_H

#include "sphere.h"
#include "obstacle.h"
#include <vector>


using namespace std;

class Glial : public Obstacle
{
    public : 
    int id;                                         /*!< ID of glial */
    Sphere soma;                                    /*!< soma of glial */

    std::vector<std::vector<Sphere>> ramification_spheres; /*!< ramification spheres of glial */
    std::vector<Eigen::Vector2d> Box;               /*!< Box with <min, max> for each axis (x,y,z) */
    std::vector<double> principle_processes_lengths; /*!< lengths of principle processes */
    std::vector<Eigen::Vector3d> attractors;        /*!< Attractors for the glial cell */
    double volume_soma;                                  /*!< Volume of glial */
    double volume_processes;                                  /*!< Volume of processes of glial */
    double minimum_radius;
    Glial();

    ~Glial();

    Glial(const int &id_, const Sphere &soma_)
    {
        id = id_;
        soma = soma_;
        double sph_highest_x_val = soma.center[0]+ soma.radius;
        double sph_lowest_x_val = soma.center[0] -soma.radius;
        double sph_highest_y_val = soma.center[1] +soma.radius;
        double sph_lowest_y_val = soma.center[1] -soma.radius;
        double sph_highest_z_val = soma.center[2] +soma.radius;
        double sph_lowest_z_val = soma.center[2] -soma.radius;
        volume_soma = M_PI*pow(soma.radius, 3)*4/3;
        volume_processes = 0.0;
        minimum_radius = soma.radius/20.0;;
        //x
        Box.push_back({sph_lowest_x_val, sph_highest_x_val});
        //y
        Box.push_back({sph_lowest_y_val, sph_highest_y_val});
        //z
        Box.push_back({sph_lowest_z_val, sph_highest_z_val});

        ramification_spheres.clear();
        attractors.clear();
        
    }


    Glial& operator=(const Glial& other) {
        if (this != &other) {
            id = other.id;
            soma = other.soma;
            ramification_spheres = other.ramification_spheres;
            Box = other.Box;
            principle_processes_lengths = other.principle_processes_lengths;
            volume_soma = other.volume_soma;
            volume_processes = other.volume_processes;
            minimum_radius = other.minimum_radius;
            attractors = other.attractors;
            
        }
        return *this;
    }

    /*!
    *   \brief  Add a sphere to the Box for it to be considered when checking if anothe relement is close by
    *   \param  sph sphere to add
     */

    void update_Box();

    /*!
    *   \brief  Check if a point is inside the glial cell (inside the Box)
    *   \param  point point to check
    *   \param distance_to_be_inside distance to be inside the glial cell
    *  \return true if point is inside the glial cell
     */

    bool isNearGlialCell(const Eigen::Vector3d &position, const double &distance_to_be_inside) const;
    /*!
    *   \brief  Check if a sphere collides with the glial cell
    *   \param  sph sphere to check
     */
    bool collides_with_GlialCell(const Sphere &sph) const;

    void compute_processes_icvf(const int &factor);
};

#endif // GLIAL_H