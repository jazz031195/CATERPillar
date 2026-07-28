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
#include "SphereGrid.h"
#include <vector>


using namespace std;



class Glial : public Obstacle
{
    public : 


    int id;                                         /*!< ID of glial */
    Sphere soma;                                    /*!< soma of glial */

    std::vector<std::vector<Sphere>> ramification_spheres; /*!< ramification spheres of glial */
    std::vector<std::vector<double>> lengths_branches;               /*!< lengths of branches */
    std::vector<Eigen::Vector3d> attractors;        /*!< Attractors for the glial cell */
    double volume_soma;                                  /*!< Volume of glial */
    double volume_processes;                                  /*!< Volume of processes of glial */
    double minimum_radius;
    bool allow_branching;                          /*!< Allow branching of glial */
    int next_sphere_id;                             /*!< Monotonic counter handing out the id of the next sphere added anywhere in this cell (soma excluded), so ids stay unique within the cell regardless of overlap factor -- parent_id links and branch-parent lookups rely on that uniqueness. */
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
        next_sphere_id = soma.id + 1; // soma occupies id soma.id; branch spheres must not collide with it

        ramification_spheres.clear();
        attractors.clear();
        lengths_branches.clear();

    }


    Glial& operator=(const Glial& other) {
        if (this != &other) {
            id = other.id;
            soma = other.soma;
            ramification_spheres = other.ramification_spheres;
            volume_soma = other.volume_soma;
            volume_processes = other.volume_processes;
            minimum_radius = other.minimum_radius;
            attractors = other.attractors;
            lengths_branches = other.lengths_branches;
            allow_branching = other.allow_branching;
            next_sphere_id = other.next_sphere_id;

        }
        return *this;
    }


    /*!
     *  \brief Adds the soma and every ramification sphere of this glial cell to the grid.
     */
    void addToGrid(SphereGrid &grid) const;

    void compute_processes_icvf(const int &factor, const Eigen::Vector3d &min_limits, const Eigen::Vector3d &max_limits);

};

#endif // GLIAL_H