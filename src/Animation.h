#ifndef Animation_H
#define Animation_H

#include <string>
#include "Eigen/Core"
#include <iostream>
#include <random>
#include <SFML/Graphics.hpp>


class Animation
{
    public:

        int plane1 = 2;
        int plane2 = 1;
        sf::RenderWindow window;

        Animation(Eigen::Vector3d min_l, Eigen::Vector3d max_l);
        ~Animation();

};

#endif //Animation_H