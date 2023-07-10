#include <iostream>
#include "Eigen/Core"
#include "Eigen/Dense"
#include <fstream>
#include "Animation.h"
#include <SFML/Graphics.hpp>
#include <SFML/OpenGL.hpp>
#include <vector>


using namespace Eigen;
using namespace std;

Animation::Animation(Eigen::Vector3d min_l, Eigen::Vector3d max_l) : window(sf::VideoMode((max_l-min_l)[plane1], (max_l-min_l)[plane2]), "Animation") {}
//Animation::Animation(Eigen::Vector3d min_l, Eigen::Vector3d max_l) : window(sf::VideoMode(500, 500), "Animation") {}

Animation::~Animation() {}

