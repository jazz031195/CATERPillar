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

    struct Box {
        double x_min;
        double x_max;
        double y_min;
        double y_max;
        double z_min;
        double z_max;
    };

    int id;                                         /*!< ID of glial */
    Sphere soma;                                    /*!< soma of glial */

    std::vector<std::vector<Sphere>> ramification_spheres; /*!< ramification spheres of glial */
    Box big_box;               /*!< Box with <min, max> for each axis (x,y,z) */
    std::vector<Box> ramification_boxes;
    Box soma_box;
    std::vector<std::vector<double>> lengths_branches;               /*!< lengths of branches */
    std::vector<Eigen::Vector3d> attractors;        /*!< Attractors for the glial cell */
    double volume_soma;                                  /*!< Volume of glial */
    double volume_processes;                                  /*!< Volume of processes of glial */
    double minimum_radius;
    bool allow_branching;                          /*!< Allow branching of glial */
    Glial();

    ~Glial();

    Glial(const int &id_, const Sphere &soma_, const bool &allow_branching_ = true)
    {
        id = id_;
        soma = soma_;
        volume_soma = M_PI*pow(soma.radius, 3)*4/3;
        volume_processes = 0.0;
        minimum_radius = soma.radius/20.0;
        allow_branching = allow_branching_;

        ramification_spheres.clear();
        ramification_boxes.clear();
        attractors.clear();
        lengths_branches.clear();

        initialise_boxes();

    }


    Glial& operator=(const Glial& other) {
        if (this != &other) {
            id = other.id;
            soma = other.soma;
            ramification_spheres = other.ramification_spheres;
            ramification_boxes = other.ramification_boxes;
            soma_box = other.soma_box;
            big_box = other.big_box;
            volume_soma = other.volume_soma;
            volume_processes = other.volume_processes;
            minimum_radius = other.minimum_radius;
            attractors = other.attractors;
            lengths_branches = other.lengths_branches;
            allow_branching = other.allow_branching;
            
        }
        return *this;
    }

    /*!
    *   \brief  Add a sphere to the Box for it to be considered when checking if anothe relement is close by
    *   \param  sph sphere to add
     */

    void update_all_boxes(const Sphere &sph);
    void update_Box(Box& box, const Sphere &sph);
    bool isNearBox(const Box& box, const Eigen::Vector3d &position, const double &distance_to_be_inside) const;
    bool collidesWithItsOwnRamification(const Sphere &sph, const double &distance_to_be_inside) const;

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
    bool collides_with_GlialCell(const Sphere &sph, const double &distance_to_be_inside) const;
    bool collidesWithRamification(const Eigen::Vector3d &position, const double &distance_to_be_inside) const;

    void compute_processes_icvf(const int &factor, const Eigen::Vector3d &min_limits, const Eigen::Vector3d &max_limits);

    void initialise_boxes();
};

#endif // GLIAL_H