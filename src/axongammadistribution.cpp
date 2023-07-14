#include "axongammadistribution.h"
#include "grow_axons.h"
#include <algorithm> // std::sort
#include <random>
#include <chrono>
#include <thread>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>
#include <mutex>
#include <queue>
#include <chrono>

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

AxonGammaDistribution::AxonGammaDistribution(unsigned &num_ax, double a, double b, Eigen::Vector3d &min_l, Eigen::Vector3d &max_l, double min_radius_, bool tortuous_)
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
    min_radius = min_radius_;
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

void AxonGammaDistribution::drawWorld(sf::Window &window)
{
    // Draw already created axons
    for (unsigned j = 0; j < axons.size(); j++)
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
        if (jkr < min_radius)
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
    for (unsigned i = 0; i < num_obstacles; ++i)
    {
        // cout << "loop number " << i << endl;
        bool next = false;
        while (!next)
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
                next = true;
            }
            else // not empty
            {
                for (auto &axon : axons)
                {

                    if (axon.isPosInsideAxon(Q, radii[i], max_radius) == false) // no overlapping
                    {
                        // cout << "adding axon to list" << endl;
                        axons.push_back(ax);
                        colours.push_back(new_colour);
                        next = true;
                        break;
                    }
                    else
                    {
                        // cout << "colliding" << endl;
                        next = false;
                    }
                }
            }
        }
    }
}

void AxonGammaDistribution::growthThread(Axon &ax, bool &can_grow, int &finished)
{
    Growth *growth = new Growth(new Axon(ax), axons, max_limits, tortuous, max_radius);
    can_grow = growth->GrowAxon(); // grows a single sphere
    if (growth->finished) 
    {
        finished = 1; // 1 for true, if the growth finished
    }
    ax = *(growth->axon_to_grow);  // adds a sphere to axon
}

void AxonGammaDistribution::parallelGrowth()
{

    // generate radii from gamma distribution
    std::vector<double> radii(num_obstacles, 0);
    generate_radii(radii);
    max_radius = radii[0];
    cout << "Parallel growth simulation" << endl;
    int stuck;
    // threshold of tries to find a position of a sphere in axon
    int stuck_thr = 1;
    int start_overs = 0;
    bool stop = false;

    // Initialize SFML window
    sf::Window window(sf::VideoMode(800, 600), "3D Visualization");
    window.setActive();

    // initializeGLUT(0, nullptr);
    // Initialize OpenGL
    sf::ContextSettings settings;
    settings.depthBits = 24; // Request a 24-bit depth buffer

    // Initialize OpenGL settings
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

    std::vector<thread> threads; // one for each axon

    while (!stop)
    {
        axons.clear();
        // create all axons and their growth
        // std::vector<Growth *> growths; // one for each sphere
        createAxons(radii); // set the axons

        bool can_grow = false;
        // bool finished = false;
        stuck = 0;
        int num_sphere = 0;
        vector<int> finished(num_obstacles, 0); // 0 for false
        bool all_finished = false;

        while (!all_finished && stuck < stuck_thr && window.isOpen()) // for each sphere
        {
            for (unsigned i = 0; i < num_obstacles; i++) // for each axon
            {
                sf::Event event;
                while (window.pollEvent(event))
                {
                    // Clear the color and depth buffers
                    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

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

                // add sphere at same time for all axons
                threads.push_back(std::thread([this, i, &can_grow, &finished]()
                                              { this->growthThread(axons[i], can_grow, finished[i]); }));

            } // end for axons

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

            // Clear the window
            ++num_sphere;
            drawWorld(window); // draw one sphere at a time
            window.display();  // Display the updated window
                               // std::this_thread::sleep_for(std::chrono::seconds(2));

            all_finished = std::all_of(finished.begin(), finished.end(), [](bool value)
                                       { return value == 1; }); // if all true, then all axons are finished growing

        } // end for spheres

        for (auto &thread : threads)
        {
            thread.join();
        }

        stop = true;
    }

    std::sort(axons.begin(), axons.end(), [](const Axon a, Axon b) -> bool
              { return a.id < b.id; });

    // messages
    std::cout << "icvf: " << icvf << " voxel size: " << max_limits[0] << std::endl;
    std::string message = "ICVF achieved: " + std::to_string(icvf * 100) + "\n";
    std::cout << message << std::endl;
}

void AxonGammaDistribution::createGammaSubstrate()
{
    // /*
    //     Generates the gamma distribution of axons.
    // */
    // // generate radii from gamma distribution
    // std::vector<double> radii(num_obstacles, 0);
    // generate_radii(radii);
    // max_radius = radii[0];
    // cout << "creating gamma substrate" << endl;
    // int stuck;
    // // threshold of tries to find a position of a sphere in axon
    // int stuck_thr = 1;
    // int start_overs = 0;
    // bool stop = false;

    // // Initialize SFML window
    // sf::Window window(sf::VideoMode(800, 600), "3D Visualization");
    // window.setActive();

    // initializeGLUT(0, nullptr);
    // // Initialize OpenGL
    // glClearColor(1.0f, 1.0f, 1.0f, 1.0f); // Set clear color to white
    // sf::ContextSettings settings;
    // settings.depthBits = 24; // Request a 24-bit depth buffer

    // // Initialize OpenGL settings
    // initializeOpenGL();

    // // Set a custom depth range
    // glDepthRange(0.0f, 1000.0f);

    // // Initialize zoom level and displacement
    // float zoomLevel = 1.0f;
    // sf::Vector2i lastMousePos;
    // bool isDragging = false;
    // constexpr float displacementFactor = 0.1f;
    // sf::Vector2i mouseDelta;
    // sf::Vector2i currentMousePos;
    // sf::Vector2i prevousDisplacement;
    // prevousDisplacement.x = 0;
    // prevousDisplacement.y = 0;
    // bool isRightDragging = false;
    // sf::Vector2i lastRightMousePos;
    // sf::Vector2i rightMouseDelta;
    // float rotationFactor = 0.5f;
    // sf::Vector2i prevousRotation;
    // prevousDisplacement.x = 0;
    // prevousDisplacement.y = 0;

    // while (!stop)
    // {
    //     axons.clear();
    //     // grow all axons
    //     for (unsigned i = 0; i < num_obstacles; i++)
    //     {

    //         std::cout << "Growing Axon number: " << i << std::endl;

    //         stuck = 0;
    //         bool can_grow = false;
    //         // generate random coordinates for begin and end
    //         Vector3d Q;
    //         Vector3d D;
    //         get_begin_end_point(Q, D);
    //         // initialise axon
    //         Axon *ax = new Axon(i, Q, D, radii[i]);
    //         bool finished = false;
    //         bool grow_straight = true;
    //         int straight_growths = 0;
    //         // std::cout << "Axon initialised at: " << Q << std::endl;

    //         GLfloat new_colour = generateRandomColor();
    //         // grow all spheres
    //         while (!finished && stuck < stuck_thr && window.isOpen())
    //         {

    //             sf::Event event;
    //             while (window.pollEvent(event))
    //             {
    //                 // Clear the color and depth buffers
    //                 glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    //                 if (event.type == sf::Event::Closed)
    //                 {
    //                     window.close();
    //                 }
    //                 else if (event.type == sf::Event::MouseWheelScrolled)
    //                 {
    //                     // Zoom in/out based on mouse scroll
    //                     if (event.mouseWheelScroll.wheel == sf::Mouse::VerticalWheel)
    //                     {
    //                         if (event.mouseWheelScroll.delta > 0)
    //                         {
    //                             zoomLevel *= 1.1f; // Increase zoom level
    //                         }
    //                         else
    //                         {
    //                             zoomLevel *= 0.9f; // Decrease zoom level
    //                         }
    //                     }
    //                 }
    //                 else if (event.type == sf::Event::MouseButtonPressed)
    //                 {
    //                     if (event.mouseButton.button == sf::Mouse::Left)
    //                     {
    //                         // Start dragging
    //                         isDragging = true;
    //                         lastMousePos = sf::Mouse::getPosition(window);
    //                     }
    //                     else if (event.mouseButton.button == sf::Mouse::Right)
    //                     {
    //                         // Start right dragging
    //                         isRightDragging = true;
    //                         lastRightMousePos = sf::Mouse::getPosition(window);
    //                     }
    //                 }
    //                 else if (event.type == sf::Event::MouseButtonReleased)
    //                 {
    //                     if (event.mouseButton.button == sf::Mouse::Left)
    //                     {
    //                         // Stop dragging
    //                         isDragging = false;
    //                     }
    //                     else if (event.mouseButton.button == sf::Mouse::Right)
    //                     {
    //                         // Stop right dragging
    //                         isRightDragging = false;
    //                     }
    //                 }
    //                 else if (event.type == sf::Event::MouseMoved)
    //                 {
    //                     if (isDragging)
    //                     {
    //                         // Calculate mouse displacement
    //                         currentMousePos = sf::Mouse::getPosition(window);
    //                         mouseDelta = currentMousePos - lastMousePos;
    //                         prevousDisplacement.x += mouseDelta.x;
    //                         prevousDisplacement.y += mouseDelta.y;

    //                         // Update last mouse position
    //                         lastMousePos = currentMousePos;
    //                     }
    //                     else if (isRightDragging)
    //                     {
    //                         // Handle right mouse drag
    //                         sf::Vector2i currentRightMousePos = sf::Mouse::getPosition(window);
    //                         rightMouseDelta = currentRightMousePos - lastRightMousePos;
    //                         prevousRotation.x += rightMouseDelta.x;
    //                         prevousRotation.y += rightMouseDelta.y;
    //                         lastRightMousePos = currentRightMousePos;
    //                     }
    //                 }
    //             }
    //             // we cannot go straight if there are no first spheres as reference
    //             if (ax->spheres.size() <= 1)
    //             {
    //                 grow_straight = false;
    //             }
    //             // std::cout << "Straight : " << grow_straight << std::endl;
    //             //  initialise growth
    //             // std::cout << "ax->spheres.size() :" << ax->spheres.size() << std::endl;

    //             Growth *growth = new Growth(ax, axons, max_limits, tortuous, max_radius, grow_straight);
    //             // grow next sphere
    //             can_grow = growth->GrowAxon();
    //             // cout << "can grow :" << can_grow << endl;
    //             //  is growth of axon finished
    //             finished = growth->finished;
    //             if (!can_grow)
    //             {
    //                 cout << "cannot grow" << endl;
    //                 stuck += 1;
    //                 // if when growing straight it collides with environment
    //                 if (grow_straight)
    //                 {
    //                     // set to false so that next step doesn't go straight
    //                     grow_straight = false;
    //                     straight_growths = 0;
    //                 }
    //             }
    //             else
    //             {
    //                 if (grow_straight)
    //                 {
    //                     // if axon has been growing straight for 4 spheres in a row
    //                     if (straight_growths >= 4)
    //                     {
    //                         // set to false so that next step doesn't go straight
    //                         grow_straight = false;
    //                         straight_growths = 0;
    //                     }
    //                     straight_growths += 1;
    //                 }
    //                 else
    //                 {
    //                     // if the sphere hadn't grown straight previously -> set to straight for next 4 spheres
    //                     grow_straight = true;
    //                 }

    //                 // draw sphere
    //                 ax = growth->axon_to_grow;

    //                 // Set up the camera position and orientation
    //                 glLoadIdentity();
    //                 gluLookAt(120.0f, 120.0f, 120.0f, // Camera position
    //                           50.0f, 50.0f, 50.0f,    // Target position
    //                           0.0f, 0.0f, 100.0f);    // Up vector
    //                 // Apply rotation
    //                 glRotatef(prevousRotation.x * rotationFactor, 1.0f, 0.0f, 0.0f);
    //                 glRotatef(prevousRotation.y * rotationFactor, 0.0f, 1.0f, 0.0f);
    //                 // cout << "mouseDelta.x : " << mouseDelta.x << endl;
    //                 //  Update camera position based on mouse displacement
    //                 glTranslatef(-prevousDisplacement.x * displacementFactor, prevousDisplacement.y * displacementFactor, 0.0f);

    //                 glScalef(zoomLevel, zoomLevel, zoomLevel); // Apply zoom transformation

    //                 // Clear the window
    //                 drawWorld(ax, window, new_colour);
    //                 window.display(); // Display the updated window
    //             }
    //         }
    //         if (finished)
    //         {
    //             axons.push_back(*ax);
    //             // cout << "colour added : "<< new_colour << endl;
    //             colours.push_back(new_colour);
    //         }
    //         else
    //         {
    //             start_overs += 1;
    //             // start again
    //             i--;
    //             // cout << "start over, previous  location : "<< ax->begin << endl;
    //         }
    //     } // end for axons
    //     stop = true;
    // }

    // std::sort(axons.begin(), axons.end(), [](const Axon a, Axon b) -> bool
    //           { return a.id < b.id; });

    // // messages
    // std::cout << "icvf: " << icvf << " voxel size: " << max_limits[0] << std::endl;
    // std::string message = "ICVF achieved: " + std::to_string(icvf * 100) + "\n";
    // std::cout << message << std::endl;
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
    double AreaV = (max_limits[0] - min_limits[0]) * (max_limits[1] - min_limits[1]) * (max_limits[2] - min_limits[2]);

    double AreaC = 0;

    double tortuosity;

    for (uint i = 0; i < axons.size(); i++)
    {

        double ax_length = 0;
        double mean_rad = 0;
        int num_spheres = 0;

        // if twin
        if (i > 0 && axons[i].radius == axons[i - 1].radius)
        {
            tortuosities.push_back(tortuosities[tortuosities.size() - 1]);
            continue;
        }
        else if (axons[i].spheres.size() > 1)
        {

            for (uint j = 1; j < axons[i].spheres.size(); j++)
            {
                double l = (axons[i].spheres[j - 1].center - axons[i].spheres[j].center).norm();
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

            tortuosity = ax_length / ((axons[i].begin - axons[i].end).norm());

            tortuosities.push_back(tortuosity);

            AreaC += M_PI * mean_rad * mean_rad * ax_length;
        }
    }

    return AreaC / AreaV;
}
