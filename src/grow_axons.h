    #ifndef AXONGROWTH_H
    #define AXONGROWTH_H

    #include <vector>
    #include <Eigen/Dense>
    #include "Axon.h"
    #include "Glial.h"
    #include "grow_cells.h"
    #include "sphere.h"
    #include "threads.h"

    /**
     * @brief Handles the logic for axon growth including sphere placement and 
     *        optional myelination logic. Inherits from CellGrowth.
     */
    class AxonGrowth : public CellGrowth
    {
    public:
    Axon& axon_to_grow; 

    AxonGrowth() = delete;

    AxonGrowth(Axon& axon_to_grow_,
               const std::vector<Glial>* glial_pop1,
               const std::vector<Glial>* glial_pop2,
               const std::vector<Axon>* axons_,
               const std::vector<Blood_Vessel>* blood_vessels_,
               const Eigen::Vector3d& extended_min_limits_,
               const Eigen::Vector3d& extended_max_limits_,
               const Eigen::Vector3d& min_limits_,
               const Eigen::Vector3d& max_limits_,
               const double& std_dev_,
               const double& min_radius_);

    ~AxonGrowth();
    AxonGrowth(const AxonGrowth& other);

    // Growth and placement
    bool AddOneSphere(double radius_,
                      bool create_sphere,
                      int grow_straight,
                      const int& factor);

    void add_spheres(Sphere& sph,
                     const Sphere& last_sphere,
                     const int& factor);

    // Positioning
    Eigen::Vector3d find_next_center_straight(const double distance,
                                   const std::vector<Sphere>& spheres);

    Eigen::Vector3d find_next_center(const double dist_,
                          const std::vector<Sphere>& spheres,
                          const Eigen::Vector3d& target);

    void growthThread(double& stuck_radius, int& stuck_index, int factor, bool axon_can_shrink);
    double RandomradiusVariation();
    bool shrinkRadius(const double &radius_to_shrink, const bool& axon_can_shrink, const int &factor);
    void update_straight(bool can_grow_, int &grow_straight, int &straight_growths);

};

#endif // AXONGROWTH_H
