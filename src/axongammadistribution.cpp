#include "axongammadistribution.h"
#include "grow_axons.h"
#include <algorithm> // std::sort
#include <random>
#include <chrono>
#include <thread>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
// #include <SFML/Window.hpp>
// #include <SFML/Graphics.hpp>
#include <mutex>
#include <chrono>
#include <cmath>

using namespace std;
using namespace Eigen;
using namespace std::chrono;

//* Auxiliare method to split words in a line using the spaces*//
template <typename Out>
void _split_(const std::string &s, char delim, Out result)
{
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim))
    {
        *(result++) = item;
    }
}

std::vector<std::string> _split_line(const std::string &s, char delim)
{
    std::vector<std::string> elems;
    _split_(s, delim, std::back_inserter(elems));
    return elems;
}

AxonGammaDistribution::AxonGammaDistribution(unsigned &num_ax, int &axon_capacity_, double a, double b, Eigen::Vector3d &min_l, Eigen::Vector3d &max_l, double min_radius_, bool tortuous_, bool draw_)
{
    alpha = a;
    beta = b;
    min_limits = min_l;
    max_limits = max_l;
    axons.clear();
    icvf = 0;
    tortuosities.clear();
    tortuous = tortuous_;
    num_obstacles = num_ax;
    axon_capacity = axon_capacity_;
    min_radius = min_radius_;
    draw = draw_;
}

void AxonGammaDistribution::set_icvf(double icvf_, double x, double y) // overwrite num_obstacles
{
    icvf = icvf_;
    double av_radius = 0.5;
    num_obstacles = (x * y * icvf) / (M_PI * av_radius * av_radius);
    cout << "num obstacles " << num_obstacles << endl;
    if (num_obstacles > axon_capacity) // if less then capacity, only one batch
    {
        cout << "num_batches " << static_cast<int>(num_obstacles / axon_capacity) + 1 << endl;
    }
}

void AxonGammaDistribution::computeMinimalSize(std::vector<double> radiis, double icvf_, Eigen::Vector3d &l)
{
    /*
    Computes the minimal Voxel size for the chosen icvf, by assuming straight axons

    radiis : list of radii for each axon to grow
    icvf_ : Intracompartment volume fraction
    l : 3D vector with each of the 3 values being the length of the vector
    */

    double area = 0;

    for (uint i = 0; i < radiis.size(); i++)
    {
        area += radiis[i] * radiis[i] * M_PI;
    }

    double l_ = sqrt(area / icvf_);

    l = {l_, l_, l_};
}

void AxonGammaDistribution::displayGammaDistribution()
{
    /*
    Displays the Gamma Distribution of the axons
    */
    const int nrolls = 10000; // number of experiments
    const int nstars = 100;   // maximum number of stars to distribute
    string message;
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::gamma_distribution<double> distribution(alpha, beta);

    int p[11] = {};

    for (int i = 0; i < nrolls; ++i)
    {
        double number = distribution(generator);

        if (number < 10)
            ++p[int(number)];
        else
            ++p[10];
    }

    for (int i = 0; i < 9; ++i)
    {
        message = std::to_string(i) + "-" + std::to_string(i + 1) + ": " + std::string(p[i] * nstars / nrolls, '*');
    }
    message = "9-10:" + std::string(p[9] * nstars / nrolls, '*');
    std::cout << message << endl;

    message = ">10: " + std::string(p[10] * nstars / nrolls, '*') + "\n";
    std::cout << message << endl;
}

bool AxonGammaDistribution::check_borders(Eigen::Vector3d pos, double distance_to_border, Eigen::Vector2d &twin_delta_pos)
{

    Eigen::Vector2d new_min_limits = {distance_to_border, distance_to_border};
    Eigen::Vector2d new_max_limits = {max_limits[0] - distance_to_border, max_limits[1] - distance_to_border};

    twin_delta_pos = {0.0, 0.0};

    if ((pos[0] - new_min_limits[0]) < 0)
    {
        // min plane of x
        twin_delta_pos[0] = max_limits[0];
    }
    if ((pos[1] - new_min_limits[1]) < 0)
    {
        // min plane of y
        twin_delta_pos[1] = max_limits[1];
    }
    if ((pos[0] - new_max_limits[0]) > 0)
    {
        // max plane of x
        twin_delta_pos[0] = -max_limits[0];
    }
    if ((pos[1] - new_max_limits[1]) > 0)
    {
        // max plane of y
        twin_delta_pos[1] = -max_limits[1];
    }
    if (twin_delta_pos != Eigen::Vector2d({0.0, 0.0}))
    {
        return true;
    }
    else
    {
        return false;
    }
}

void display_progress(double nbr_axons, double number_obstacles)
{
    int cTotalLength = 50;
    double lProgress = nbr_axons / number_obstacles;
    std::cout << "\r[" <<                                     //'\r' aka carriage return should move printer's cursor back at the beginning of the current line
        string(int(cTotalLength * lProgress), '*') <<         // printing filled part
        string(int(cTotalLength * (1 - lProgress)), '-') <<   // printing empty part
        "] " << nbr_axons << "/" << number_obstacles << endl; // printing percentage
}

void AxonGammaDistribution::get_begin_end_point(Eigen::Vector3d &Q, Eigen::Vector3d &D)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> udist(0, 1);

    double t = udist(gen);
    // double distance_to_border = radiis[i]*sqrt(1+volume_inc_perc) + step_length;
    double x = (t * (max_limits[0])) + (1 - t) * (min_limits[0]);
    t = udist(gen);
    double y = (t * (max_limits[1])) + (1 - t) * (min_limits[1]);

    Q = {x, y, min_limits[2]};
    D = {x, y, max_limits[2]};
}

// Function to initialize the GLUT library
void initializeGLUT(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
    glutInitWindowSize(1000, 1000);
}
void drawSphere(float x, float y, float z, float radius, GLfloat colour)
{
    glPushMatrix();
    glTranslatef(x, y, z);

    GLfloat matAmbient[] = {colour};                  // Ambient color: (R, G, B) = (0.2, 0.2, 0.6)
    GLfloat matDiffuse[] = {1.0f, 1.0f, 1.0f, 1.0f};  // Diffuse color: (R, G, B) = (0.4, 0.4, 0.8)
    GLfloat matSpecular[] = {1.0f, 1.0f, 1.0f, 1.0f}; // Specular color: (R, G, B) = (1.0, 1.0, 1.0)
    GLfloat shininess = 100.0f;

    glMaterialfv(GL_FRONT, GL_AMBIENT, matAmbient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, matDiffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, matSpecular);
    glMaterialf(GL_FRONT, GL_SHININESS, shininess);
    glEnable(GL_DEPTH_TEST); // Enable depth testing

    // Draw the wireframe sphere with a white outline
    glColor3f(colour, colour, colour);
    glutWireSphere(radius, 20, 20);

    glPopMatrix();
}

// Function to initialize the OpenGL settings
void initializeOpenGL()
{
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f); // Set clear color to white

    glEnable(GL_DEPTH_TEST); // Enable depth testing
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0f, 1.0f, 0.1f, 1000.0f); // Set perspective projection parameters

    glMatrixMode(GL_MODELVIEW);

    // Enable lighting
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    // Set up light source properties
    GLfloat lightPosition[] = {100.0f, 100.0f, 100.0f, 1.0f};
    GLfloat lightAmbient[] = {0.2f, 0.2f, 0.2f, 1.0f};
    GLfloat lightDiffuse[] = {0.8f, 0.8f, 0.8f, 1.0f};
    GLfloat lightSpecular[] = {1.0f, 1.0f, 1.0f, 1.0f};

    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
    glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmbient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpecular);

    // Enable depth testing
    glEnable(GL_DEPTH_TEST);
}

GLfloat generateRandomColor()
{
    srand(static_cast<unsigned int>(time(0)));

    GLfloat color = static_cast<GLfloat>(rand()) / RAND_MAX;
    return color;
}

void AxonGammaDistribution::drawWorld(unsigned int row, sf::Window &window) // for parallel growth
{
    int num_col = num_obstacles / num_batches;

    for (unsigned int j = row * num_col; j < (row + 1) * num_col; ++j) // new axon batch
    {
        for (unsigned i = 0; i < axons[j].spheres.size(); i++)
        {
            float x = axons[j].spheres[i].center[0];
            float y = axons[j].spheres[i].center[1];
            float z = axons[j].spheres[i].center[2];
            float radius = axons[j].spheres[i].radius;
            drawSphere(x, y, z, radius, colours[j]);
        }
    }

    for (unsigned int j = 0; j < row * num_col; ++j) // already created batches
    {
        for (unsigned i = 0; i < axons[j].spheres.size(); i++)
        {
            float x = axons[j].spheres[i].center[0];
            float y = axons[j].spheres[i].center[1];
            float z = axons[j].spheres[i].center[2];
            float radius = axons[j].spheres[i].radius;
            drawSphere(x, y, z, radius, colours[j]);
        }
    }
}

void AxonGammaDistribution::generate_radii(std::vector<double> &radiis)
{

    std::random_device rd;
    std::default_random_engine generator(rd());
    std::gamma_distribution<double> distribution(alpha, beta);

    int tried = 0;

    for (unsigned i = 0; i < num_obstacles; ++i)
    {
        if (tried > 1000)
        {
            std::string message = "Radii distribution cannot be sampled [Min. radius Error]\n";
            std::cout << message << std::endl;
            assert(0);
        }
        double jkr = distribution(generator);

        // generates the radii in a list
        if (jkr < min_radius && jkr > 5) // max_radius = 5um
        {
            i--;
            tried++;
            continue;
        }
        tried = 0;

        radiis[i] = jkr;
    }
    // sorts the radii
    std::sort(radiis.begin(), radiis.end(), [](const double a, double b) -> bool
              { return a > b; });
}

void AxonGammaDistribution::createAxons(std::vector<double> radii)
{

    max_radius = radii[0];
    int tries = 0;
    for (unsigned i = 0; i < num_obstacles; ++i) // create axon for each radius
    {
        bool next = false;
        while (!next && tries < 1000000)
        {
            GLfloat new_colour = generateRandomColor();
            Vector3d Q, D;
            get_begin_end_point(Q, D);
            Axon ax = Axon(i, Q, D, radii[i]);
            Dynamic_Sphere sphere = Dynamic_Sphere(i, i, Q, radii[i]);
            ax.add_sphere(sphere);
            if (axons.empty())
            {
                axons.push_back(ax);
                colours.push_back(new_colour);
                tries = 0; // reset number of tries;
                next = true;
            }
            else // not empty
            {
                bool overlap = false;
                for (auto &axon : axons)
                {
                    if (axon.isPosInsideAxon(Q, radii[i], max_radius)) // overlap
                    {
                        overlap = true;
                        ++tries;
                        next = false;
                        break;
                    }
                }
                if (overlap == false) // after comparing with all axons
                {
                    axons.push_back(ax);
                    colours.push_back(new_colour);
                    tries = 0; // reset number of tries;
                    next = true;
                }
            }
        }
        if (tries >= 1000000)
        {
            cout << "No more space for more axons" << endl;
            num_obstacles = axons.size();
            if (num_obstacles < radii.size())
            {
                // Erase elements from 'num_obstacles' index to the end of the vector
                radii.erase(radii.begin() + num_obstacles, radii.end());
            }
            cout << "New num_obstacles = " << num_obstacles << endl;
            cout << "By default, new num_batches = " << num_batches << endl;
            break; // no more axons added due to lack of space
        }
    }
}

void AxonGammaDistribution::radiusVariation(Axon &axon, int time, double radius_)
{
    double p = 0.75;
    double frequency = 1. / radius_;
    double phase_shift = 0.;
    double amplitude = (1 - p) * radius_;
    double fluctuation = amplitude * sin(2.0 * M_PI * frequency * time + phase_shift);
    axon.radius = radius_ * (1 + p) / 2 + fluctuation;
    if (axon.radius < min_radius)
    {
        axon.radius = min_radius;
    }
}

void AxonGammaDistribution::dichotomy(Growth growth, Axon &axon, double &min_rad, double &max_rad, int radius, int &tries, double &rad)
{
    ++tries;

    double current_rad = (max_rad + min_rad) / 2;
    // Dynamic_Sphere s(axon.spheres.size() - 1, axon.id, axon.begin, current_rad);
    growth = Growth(axon, axons, max_limits, tortuous, current_rad, max_radius,false);
    bool can_grow = growth.GrowAxon();
    if (can_grow) // solution is in greater half
    {
        rad = current_rad; // update the last successful radius
        min_rad = current_rad;
    }

    else // solution is in lower half
    {
        max_rad = current_rad;
    }

    if (tries < 5)
    {
        dichotomy(growth, axon, min_rad, max_rad, radius, tries, rad);
    }

    axon.radius = rad;
}

bool AxonGammaDistribution::shrinkRadius(Growth growth, Axon &axon, int radius)
{
    double current_rad = min_radius;
    // Dynamic_Sphere s(axon.spheres.size() - 1, axon.id, axon.spheres[-1].center, current_rad);
    growth = Growth(axon, axons, max_limits, tortuous, current_rad, max_radius,false);

    double can_grow = growth.GrowAxon();
    int tries = 0;
    if (can_grow) // shrinking is useful
    {
        double rad = current_rad; // last successful radius
        double max_rad = radius;
        double min_rad = current_rad;

        dichotomy(growth, axon, min_rad, max_rad, radius, tries, rad);
        return true;
    }
    else
    {
        cout << "shrinking not useful" << endl;
        return false; // shrinking won't be useful
    }
}

void AxonGammaDistribution::growthThread(int index, bool &can_grow, int &finished, int &grow_straight, int &stuck, int &straight_growths, int time, double radius, int &shrink_tries)
{
    bool grow_straight_;
    if (grow_straight == 1)
    {
        grow_straight_ = true;
    }
    if (grow_straight == 0)
    {
        grow_straight_ = false;
    }

    radiusVariation(axons[index], time, radius);  // radius fluctuation
    Growth growth = Growth(axons[index], axons, max_limits, tortuous, axons[index].radius, max_radius, grow_straight_);

    try
    {
        can_grow = growth.GrowAxon(); // adds sphere

        if (growth.finished)
        {
            finished = 1; // 1 for true, if the growth finished
        }
        else // still growing
        {
            if (!can_grow)
            {
                if (shrink_tries < 50)
                {
                    // cout << "Shrinking axon " << axons[index].id << ", sphere " << axons[index].spheres.size();
                    {
                        std::lock_guard<std::mutex> lock(axonsMutex); // lock the radius mutex to protect radius access
                        bool shrink = shrinkRadius(growth, axons[index], radius);

                        if (!shrink) // shrinking will not help growth
                        {
                            cout << "Axon " << axons[index].id << " is stuck : stop growing !" << endl;
                            finished = 1;
                            // stuck += 1;
                        }
                        else
                        {
                            cout << "  --> done" << endl;
                        }
                    }
                    ++shrink_tries;
                }
                else
                {
                    cout << "!! Axon " << axons[index].id << " : failed growing !!" << endl;
                    finished = 1;
                }

                if (grow_straight_) // if when growing straight it collides with environment
                {
                    grow_straight = 0; // set to false so that next step doesn't go straight
                    straight_growths = 0;
                }
            }
            else
            {
                if (grow_straight_)
                {
                    if (straight_growths >= 4) // if axon has been growing straight for 4 spheres in a row
                    {
                        grow_straight = 0; // set to false so that next step doesn't go straight
                        straight_growths = 0;
                    }
                    straight_growths += 1;
                }
                else
                {
                    // if the sphere hadn't grown straight previously . set to straight for next 4 spheres
                    grow_straight = 1; // set to true
                }
                {
                    std::lock_guard<std::mutex> lock(axonsMutex); // avoid concurrent modifications and data races
                    axons[index] = growth.axon_to_grow;           // updates axon list
                }

            }
        }
    }
    catch (const std::exception &ex)
    {
        std::cerr << "Exception occurred " << std::endl;
        throw;
    }
}

void AxonGammaDistribution::growthVisualisation()
{

#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>

    // generate radii from gamma distribution
    std::vector<double> radii(num_obstacles, 0);
    generate_radii(radii);
    max_radius = radii[0];
    cout << "Parallel growth simulation" << endl;
    // threshold of tries to find a position of a sphere in axon
    int stuck_thr = 1;
    int start_overs = 0;
    bool stop = false;

    // Initialize SFML window
    sf::Window window(sf::VideoMode(800, 600), "3D Visualization");

    window.setActive();

    // Initialize OpenGL
    sf::ContextSettings settings;
    settings.depthBits = 24; // Request a 24-bit depth buffer

    // Initialize OpenGL settings
    // initializeGLUT(0, nullptr);
    initializeOpenGL();

    // Set a custom depth range
    glDepthRange(0.0f, 1000.0f);

    // Initialize zoom level and displacement
    float zoomLevel = 1.0f;
    sf::Vector2i lastMousePos;
    bool isDragging = false;
    constexpr float displacementFactor = 0.1f;
    sf::Vector2i mouseDelta;
    sf::Vector2i currentMousePos;
    sf::Vector2i prevousDisplacement;
    prevousDisplacement.x = 0;
    prevousDisplacement.y = 0;
    bool isRightDragging = false;
    sf::Vector2i lastRightMousePos;
    sf::Vector2i rightMouseDelta;
    float rotationFactor = 0.5f;
    sf::Vector2i prevousRotation;
    prevousDisplacement.x = 0;
    prevousDisplacement.y = 0;

    std::vector<std::vector<std::thread>> threads; // one for each axon

    while (!stop)
    {
        axons.clear();
        createAxons(radii);           // set the axons
        std::vector<int> num_subsets; // depending on which batch

        // Set number of batches:
        if (num_obstacles < axon_capacity)
        {
            num_batches = 1;
            num_subsets = std::vector<int>(1, num_obstacles);
        }
        if (num_obstacles % axon_capacity == 0) 
        {
            num_batches = num_obstacles / axon_capacity;
            num_subsets = std::vector<int>(num_batches, axon_capacity);
        }
        else
        {
            int left = num_obstacles % axon_capacity;
            num_batches = static_cast<int>(num_obstacles / axon_capacity) + 1;
            num_subsets = std::vector<int>(num_batches - 1, axon_capacity); // max capacity axons in all batches except last
            num_subsets.push_back(left);                                    // last batch has the extra axons left
        }

        // Start growth:
        for (unsigned j = 0; j < num_batches; j++) // batches of axon growth
        {
            cout << "---   Batch " << j << "   --- " << endl;
            std::vector<thread> row;
            bool can_grow = false;
            int stuck = 0;
            vector<int> finished(num_subsets[j], 0);         // 0 for false
            vector<int> grow_straight(num_subsets[j], 1);    // 1 for true
            vector<int> straight_growths(num_subsets[j], 0); // for each axon
            vector<int> shrink_tries(num_subsets[j], 0);     // for each axon
            bool all_finished = false;

            int time = 0; // used for radius fluctuation

            while (!all_finished && stuck < stuck_thr && window.isOpen()) // for each sphere
            {
                for (unsigned i = 0; i < num_subsets[j]; i++) // for each axon
                {
                    if (finished[i] == 0) // if the axon is not done growing
                    {
                        // Clear the color and depth buffers
                        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
                        sf::Event event;
                        while (window.pollEvent(event))
                        {
                            if (event.type == sf::Event::Closed)
                            {
                                window.close();
                            }
                            else if (event.type == sf::Event::MouseWheelScrolled)
                            {
                                // Zoom in/out based on mouse scroll
                                if (event.mouseWheelScroll.wheel == sf::Mouse::VerticalWheel)
                                {
                                    if (event.mouseWheelScroll.delta > 0)
                                    {
                                        zoomLevel *= 1.1f; // Increase zoom level
                                    }
                                    else
                                    {
                                        zoomLevel *= 0.9f; // Decrease zoom level
                                    }
                                }
                            }
                            else if (event.type == sf::Event::MouseButtonPressed)
                            {
                                if (event.mouseButton.button == sf::Mouse::Left)
                                {
                                    // Start dragging
                                    isDragging = true;
                                    lastMousePos = sf::Mouse::getPosition(window);
                                }
                                else if (event.mouseButton.button == sf::Mouse::Right)
                                {
                                    // Start right dragging
                                    isRightDragging = true;
                                    lastRightMousePos = sf::Mouse::getPosition(window);
                                }
                            }
                            else if (event.type == sf::Event::MouseButtonReleased)
                            {
                                if (event.mouseButton.button == sf::Mouse::Left)
                                {
                                    // Stop dragging
                                    isDragging = false;
                                }
                                else if (event.mouseButton.button == sf::Mouse::Right)
                                {
                                    // Stop right dragging
                                    isRightDragging = false;
                                }
                            }
                            else if (event.type == sf::Event::MouseMoved)
                            {
                                if (isDragging)
                                {
                                    // Calculate mouse displacement
                                    currentMousePos = sf::Mouse::getPosition(window);
                                    mouseDelta = currentMousePos - lastMousePos;
                                    prevousDisplacement.x += mouseDelta.x;
                                    prevousDisplacement.y += mouseDelta.y;

                                    // Update last mouse position
                                    lastMousePos = currentMousePos;
                                }
                                else if (isRightDragging)
                                {
                                    // Handle right mouse drag
                                    sf::Vector2i currentRightMousePos = sf::Mouse::getPosition(window);
                                    rightMouseDelta = currentRightMousePos - lastRightMousePos;
                                    prevousRotation.x += rightMouseDelta.x;
                                    prevousRotation.y += rightMouseDelta.y;
                                    lastRightMousePos = currentRightMousePos;
                                }
                            }
                        }

                        int index = j * num_obstacles / num_batches + i; // position within axons from all batches

                        if (axons[index].spheres.size() <= 1) // we cannot go straight if there are no first spheres as reference
                        {
                            grow_straight[i] = 0; // false
                        }

                        // add sphere at same time for all axons :
                        row.emplace_back([this, index, i, &can_grow, &finished, &grow_straight, &stuck, &straight_growths, time, radii, &shrink_tries]()
                                         { this->growthThread(index, can_grow, finished[i], grow_straight[i], stuck, straight_growths[i], time, radii[index], shrink_tries[i]); });
                    }

                } // end for axons

                ++time;

                // Set up the camera position and orientation
                glLoadIdentity();
                gluLookAt(120.0f, 120.0f, 120.0f, // Camera position
                          50.0f, 50.0f, 50.0f,    // Target position
                          0.0f, 0.0f, 100.0f);    // Up vector
                // Apply rotation
                glRotatef(prevousRotation.x * rotationFactor, 1.0f, 0.0f, 0.0f);
                glRotatef(prevousRotation.y * rotationFactor, 0.0f, 1.0f, 0.0f);
                //  Update camera position based on mouse displacement
                glTranslatef(-prevousDisplacement.x * displacementFactor, prevousDisplacement.y * displacementFactor, 0.0f);

                glScalef(zoomLevel, zoomLevel, zoomLevel); // Apply zoom transformation

                drawWorld(j, window); // draw one sphere at a time, j is the row number

                window.display(); // Display the updated window
                                  // std::this_thread::sleep_for(std::chrono::seconds(2));

                all_finished = std::all_of(finished.begin(), finished.end(), [](bool value)
                                           { return value == 1; }); // if all true, then all axons are finished growing

            } // end for spheres

            threads.emplace_back(std::move(row));

        } // end for batch

        for (auto &row : threads)
        {
            for (auto &thread : row)
            {
                thread.join();
            }
        }
        stop = true;
    }

    // messages
    std::cout << "icvf: " << icvf << " voxel size: " << max_limits[0] << std::endl;
    std::string message = "ICVF achieved: " + std::to_string(icvf * 100);
    std::cout << message << std::endl;
    std::cout << "Number axons: " << num_obstacles << std::endl;
    std::cout << "Number batches: " << num_batches << std::endl;
    std::cout << "compute icvf " << computeICVF() << "\n"
              << std::endl;
}

void AxonGammaDistribution::parallelGrowth()
{
    // generate radii from gamma distribution
    std::vector<double> radii(num_obstacles, 0);
    generate_radii(radii);
    max_radius = radii[0];
    cout << "Parallel growth simulation" << endl;
    // threshold of tries to find a position of a sphere in axon
    int stuck_thr = 1;
    int start_overs = 0;
    bool stop = false;

    std::vector<std::vector<std::thread>> threads; // one for each axon

    while (!stop)
    {
        axons.clear();
        createAxons(radii); // set the axons
        int num_subsets = num_obstacles / num_batches;

        for (unsigned j = 0; j < num_batches; j++) // batches of axon growth
        {
            cout << "---   Batch " << j << "   --- " << endl;
            std::vector<thread> row;
            bool can_grow = false;
            int stuck = 0;
            vector<int> finished(num_subsets, 0);         // 0 for false
            vector<int> grow_straight(num_subsets, 1);    // 1 for true
            vector<int> straight_growths(num_subsets, 0); // for each axon
            vector<int> shrink_tries(num_subsets, 0);     // for each axon
            bool all_finished = false;

            int time = 0; // used for radius fluctuation

            while (!all_finished && stuck < stuck_thr) // for each sphere
            {
                for (unsigned i = 0; i < num_subsets; i++) // for each axon
                {
                    if (finished[i] == 0) // if the axon is not done growing
                    {
                        int index = j * num_obstacles / num_batches + i; // position within axons from all batches
                        if (axons[index].spheres.size() <= 1)            // we cannot go straight if there are no first spheres as reference
                        {
                            grow_straight[i] = 0; // false
                        }
                        // add sphere at same time for all axons :
                        row.emplace_back([this, index, i, &can_grow, &finished, &grow_straight, &stuck, &straight_growths, time, radii, &shrink_tries]()
                                         { this->growthThread(index, can_grow, finished[i], grow_straight[i], stuck, straight_growths[i], time, radii[index], shrink_tries[i]); });
                    }

                } // end for axons

                ++time;
                all_finished = std::all_of(finished.begin(), finished.end(), [](bool value)
                                           { return value == 1; }); // if all true, then all axons are finished growing

            } // end for spheres

            threads.emplace_back(std::move(row));
            // row.clear();

        } // end for batch

        for (auto &row : threads)
        {
            for (auto &thread : row)
            {
                thread.join();
            }
        }
        stop = true;
    }

    // messages
    std::cout << "icvf: " << icvf << " voxel size: " << max_limits[0] << std::endl;
    std::string message = "ICVF achieved: " + std::to_string(icvf * 100);
    std::cout << message << std::endl;
    std::cout << "Number axons: " << num_obstacles << std::endl;
    std::cout << "Number batches: " << num_batches << std::endl;
    std::cout << "compute icvf " << computeICVF() << "\n"
              << std::endl;
}

void AxonGammaDistribution::printSubstrate(ostream &out)
{
    double scale = 1;
    out << scale << endl;
    out << icvf << endl;
    out << min_limits[2] << endl;
    out << max_limits[2] << endl;

    for (unsigned i = 0; i < axons.size(); i++)
    {
        for (unsigned s = 0; s < axons[i].spheres.size(); s++)
        {

            out << axons[i].spheres[s].center[0] / scale << " " << axons[i].spheres[s].center[1] / scale << " " << axons[i].spheres[s].center[2] / scale << " " << axons[i].spheres[s].radius / scale << endl;
        }
        out << "Axon: " << axons[i].id << " tortuosity:" << tortuosities[i] << endl;
    }
}

double AxonGammaDistribution::computeICVF()
{

    if (axons.size() == 0)
        return 0;
    double AreaV = (max_limits[0] - min_limits[0]) * (max_limits[1] - min_limits[1]) * (max_limits[2] - min_limits[2]); // total volume

    double AreaC = 0;

    double tortuosity;

    for (uint i = 0; i < axons.size(); i++) // for all axons
    {

        double ax_length = 0;
        double mean_rad = 0;
        int num_spheres = 0;

        if (axons[i].spheres.size() > 1)
        {

            for (uint j = 1; j < axons[i].spheres.size(); j++)
            {
                double l = (axons[i].spheres[j - 1].center - axons[i].spheres[j].center).norm(); // distance between centers
                ax_length += l;
                mean_rad += axons[i].spheres[j].radius;
                num_spheres += 1;
            }
            if (num_spheres != 0)
            {
                mean_rad = mean_rad / num_spheres;
            }
            else
            {
                mean_rad = 0;
            }

            tortuosity = ax_length / ((axons[i].begin - axons[i].end).norm()); // ( total distance between all centers / distance between first and last )

            tortuosities.push_back(tortuosity);

            AreaC += M_PI * mean_rad * mean_rad * ax_length; // volume of axon cylinder
        }
    }

    return AreaC / AreaV; // ( total axons volume / total volume )
}

void AxonGammaDistribution::axonDensityMap()
{
    // Create an SFML window
    sf::RenderWindow window(sf::VideoMode(1000, 1000), "2D Visualization");

    // Define a list of axons, with predetermined radii and positions
    std::vector<double> radii(num_obstacles, 0);
    generate_radii(radii);
    createAxons(radii);

    bool finished = false;
    while (window.isOpen() && !finished)
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
            {
                window.close();
            }
        }

        // Clear the window
        window.clear(sf::Color(200, 200, 200));
        // Create and draw the circles at the fixed position
        for (size_t i = 0; i < num_obstacles; ++i)
        {
            double radius = radii[i] * 50;
            Eigen::Vector3d position = {axons[i].begin[0] * 50, axons[i].begin[1] * 50, axons[i].begin[2] * 50};
            cout << position << endl;
            sf::CircleShape circle(radius);
            circle.setPosition(position[0] - radius, position[1] - radius);
            circle.setFillColor(sf::Color(0, 0, 255, 128)); // Slightly more transparent blue (alpha value 64)
            circle.setOutlineColor(sf::Color(50, 50, 50));  // Darker gray outline
            circle.setOutlineThickness(0.03f);              // Set the outline thickness

            window.draw(circle);
        }

        // Display the drawn circles
        window.display();
        std::this_thread::sleep_for(std::chrono::seconds(30));

        finished = true;
    }
}

void AxonGammaDistribution ::create_SWC_file(std::ostream &out)
{
    out << "id Type X Y Z R P" << endl;
    int previous_id;
    for (uint i = 0; i < axons.size(); i++)
    {
        previous_id = -1; // parent for first sphere of each axon
        for (uint j = 0; j < axons[i].spheres.size(); j++)
        {
            out << (i + 1) * j << " axon " << axons[i].spheres[j].center[0] << " " << axons[i].spheres[j].center[1] << " " << axons[i].spheres[j].center[2] << " " << axons[i].spheres[j].radius << " " << previous_id << endl;
            previous_id = (i + 1) * j;
        }
    }
}

void AxonGammaDistribution ::radius_file(std::ostream &out)
{
    out << "ax_id   Type   sph_id   Type2   R  " << endl;
    out << endl;

    for (uint i = 0; i < axons.size(); i++) // for each axon
    {
        for (uint j = 0; j < axons[i].spheres.size(); j++) // for each sphere
        {
            out << i << "     axon    " <<  axons[i].spheres[j].center - axons[i].begin  << "  " << "     sphere    " << axons[i].spheres[j].radius << endl;
        }
    }
}