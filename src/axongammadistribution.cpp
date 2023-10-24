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

AxonGammaDistribution::AxonGammaDistribution(double target_icvf_, int &axon_capacity_, double a, double b,
                                             Eigen::Vector3d &min_l, Eigen::Vector3d &max_l, double min_radius_,
                                             bool tortuous_, bool draw_, int regrow_thr_,  double beading_variation)
{
    alpha = a;
    beta = b;
    min_limits = min_l;
    max_limits = max_l;
    axons.clear();
    icvf = 0;
    tortuosities.clear();
    tortuous = tortuous_;
    num_obstacles = 0;
    axon_capacity = axon_capacity_;
    min_radius = min_radius_;
    draw = draw_;
    regrow_thr = regrow_thr_;
    variation_perc = beading_variation;
    target_icvf = target_icvf_;

}


bool AxonGammaDistribution::withinBounds(Eigen::Vector3d pos, double distance)
{
    bool within;
    for (int i = 0; i < 2; i++) // check for all dimensions
    {
        if ((pos[i] < max_limits[i] + distance) && (pos[i] > min_limits[i] - distance))
        {
            within = true;
        }
        else
        {
            within = false;
            break;
        }
    }
    return within;
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

// Set substrate attributes
void AxonGammaDistribution::generate_radii(std::vector<double> &radiis)
{

    std::random_device rd;
    std::default_random_engine generator(rd());
    std::gamma_distribution<double> distribution(alpha, beta);
    radiis.clear();

    // take into account variation of radius -> increase initial icvf to compensate
    double icvf_to_reach = target_icvf/((variation_perc+1)/2);

    int tried = 0;

    double icvf_ = 0;

    double AreaIntra = 0;
    double AreaTot = max_limits[0]*max_limits[1];

    cout << "icvf_to_reach :" << icvf_to_reach << endl;

    while (icvf_< icvf_to_reach)
    {
        if (tried > 1000)
        {
            std::string message = "Radii distribution cannot be sampled [Min. radius Error]\n";
            std::cout << message << std::endl;
            assert(0);
        }
        double jkr = distribution(generator);

        // generates the radii in a list
        if (jkr > min_radius) 
        {
            radiis.push_back(jkr);
            tried = 0;
            AreaIntra += jkr*jkr*M_PI;
            icvf_ = AreaIntra/AreaTot;
        }
        else{
            tried += 1;
        }
        

    }
    // sorts the radii
    std::sort(radiis.begin(), radiis.end(), [](const double a, double b) -> bool
              { return a > b; });

    max_radius = radii[0];
    //cout << "Maximum radius :" << max_radius << endl;

    num_obstacles = radii.size();

    //std::cout << "Number of axons :" << num_obstacles << endl;
}

double AxonGammaDistribution::computeAreaICVF(std::vector<Axon> new_axons)
{
    double AreaIntra = 0;
    double AreaTot = max_limits[0]*max_limits[1];
    //double area_circle;
    for (auto &axon : new_axons)
    {
        
        if (withinBounds(axon.spheres[0].center,axon.spheres[0].radius)){
            AreaIntra += axon.radius*axon.radius*M_PI;
        }
        else if (withinBounds(axon.spheres[0].center,0)){
            AreaIntra += axon.radius*axon.radius*M_PI/2;
        }
        
    }
    return AreaIntra/AreaTot;
}

void AxonGammaDistribution::PlaceAxon(double radius, bool regrowth, int axon_index, int &tries, bool &next, Eigen::Vector3d Q, Eigen::Vector3d D, std::vector<Axon> &all_axons, std::vector<Axon> &new_axons){
    
    Axon ax;
    GLfloat new_colour = generateRandomColor();

    //cout << "max ID :" << maxId << endl;

    if (!regrowth)
    {
        ax = Axon(axon_index, Q, D, radius); // original axons
        //cout << " create axon with id :" << axon_index << endl;
    }
    else
    {
        // Initialize a variable to store the maximum id
        int maxId = -1; // Initialize to a value that is guaranteed to be less than any possible id

        // Iterate through the vector and find the maximum id
        for (const Axon& axon : axons) {
            if (axon.id > maxId) {
                maxId = axon.id;
            }
        }
        ax = Axon(axon_index + maxId + 1, Q, D, radius); // axons for regrow batch
        //cout << " create axon with id :" << axon_index + maxId + 1 << endl;
    }
    Sphere sphere = Sphere(0, ax.id, Q, radius);
    
    
    if (all_axons.empty())
    {
        ax.add_sphere(sphere);
        new_axons.push_back(ax);
        all_axons.push_back(ax);
        colours.push_back(new_colour);
        tries = 0; // reset number of tries;
        next = true;
    }
    else // not empty
    
    {
        bool no_overlap = canSpherebePlaced(sphere, all_axons);

        if (!no_overlap){
            ++tries;
            next = false;
            //cout <<"Overlap !" << endl;
        } 
        else // after comparing with all axons
        {
            ax.add_sphere(sphere);
            new_axons.push_back(ax);
            all_axons.push_back(ax);
            colours.push_back(new_colour);
            tries = 0; // reset number of tries;
            next = true;

        }
    }

}

bool AxonGammaDistribution::canSpherebePlaced(Sphere sph, std::vector<Axon> axs){

    int count = 0;
    for (auto &axon : axs)
    {

        if (count > 1){
            assert(0);
        }

        if (axon.id != sph.ax_id){ 
            if (axon.isSphereInsideAxon(sph)) // overlap
            { 
                return false;
            }

        } 
        else{
            count += 1;
        }

    }
    return true;
} 


double AxonGammaDistribution::dichotomy_swelling(Sphere new_sphere, double swelling_perc){

    double min_rad = new_sphere.radius ;
    double max_rad = new_sphere.radius + swelling_perc*new_sphere.radius;
    int restart = 0;
    bool overlap = false;
    double current_rad;

    new_sphere.radius = max_rad;
    overlap = !canSpherebePlaced(new_sphere, axons);

    if (overlap){
        while (restart < 5){
            current_rad = (min_rad + max_rad)/2;
            new_sphere.radius = current_rad;
            overlap = !canSpherebePlaced(new_sphere, axons);
            restart ++;
            if (!overlap){
                min_rad = current_rad;
            }
            else{
                max_rad = current_rad;
            }
        }
        current_rad = min_rad;
        return current_rad;
    }
    else{
        return max_rad;
    }
}
bool AxonGammaDistribution::swellAxons(){

    double obtained_icvf;
    double current_icvf = computeICVF();;
    
    if (current_icvf + 0.01 < icvf){
        for (uint i = 0; i < axons.size(); i++) // for all axons
        {
            //std::cout << " axon : " << axons[i].id << endl;
            for (uint j = 0; j < axons[i].spheres.size(); j++)
            {
                Sphere new_sphere = axons[i].spheres[j];

                double rad = axons[i].spheres[j].radius;

                double final_rad = dichotomy_swelling(new_sphere, 0.1*rad);
                axons[i].spheres[j].radius = final_rad;
            }
            obtained_icvf = computeICVF();

            if (obtained_icvf >= icvf){
                break;
            }
        }

        //std::cout << "Obtained icvf :" << obtained_icvf << endl;
        //std::cout << "Target icvf :" << icvf << endl;
        return true;
    }
    return false;

}

void AxonGammaDistribution::createBatch(std::vector<double> radii_, int num_subset, bool regrowth, int first_index_batch, std::vector<Axon> &new_axons)
{
    int tries = 0;
    std::vector<int> stuck_axons_id;
    long int tries_threshold = max_limits[0]*max_limits[1]*2;
    std::vector<Axon> all_axons = axons;
    new_axons.clear();

    std::vector<double>::iterator startIterator = radii_.begin() + first_index_batch;  // Start from index 
    std::vector<double>::iterator stopIterator = radii_.begin() + first_index_batch + num_subset;  // Stop when batch is finished
    
    // Create a new vector using the iterators
    std::vector<double> batch_radii(startIterator, stopIterator);

    for (unsigned i = 0; i < batch_radii.size(); ++i) // create axon for each radius
    {
        //std::cout << "radii created :" << i << "/" << batch_radii.size() << endl;
        bool next = false;
        
        while (!next && tries < tries_threshold)
        {
            Vector3d Q, D;
            get_begin_end_point(Q, D);
            int axon_index = first_index_batch + i;
            PlaceAxon(batch_radii[i], regrowth, axon_index, tries, next, Q, D, all_axons, new_axons);
        }
        if (tries >= tries_threshold)
        {
            std::cout << "Not enough space to add an axon -> discard it" << endl;
            tries = 0;
        }
    }
}

void AxonGammaDistribution::setBatches(int num_axons, std::vector<int> &num_subsets)
{
    if (num_axons < axon_capacity)
    {
        num_batches = 1;
        num_subsets = std::vector<int>(1, num_axons);
    }
    else
    {
        if (num_axons % axon_capacity == 0)
        {
            num_batches = num_axons / axon_capacity;
            num_subsets = std::vector<int>(num_batches, axon_capacity);
        }
        else
        {
            int left = num_axons % axon_capacity;
            num_batches = static_cast<int>(num_axons / axon_capacity) + 1;
            num_subsets = std::vector<int>(num_batches - 1, axon_capacity); // max capacity axons in all batches except last
            num_subsets.push_back(left);                                    // last batch has the extra axons left
        }
    }

}
void AxonGammaDistribution::createSubstrate()
{
    if (draw)
    {
        growthVisualisation();
    }
    else
    {
        parallelGrowth();
    }
}

// Drawing substrate
void AxonGammaDistribution::drawWorld(std::vector<Axon> ax_list) 
{
    for (unsigned int j = 0; j < ax_list.size(); ++j) // new axon batch
    {
        for (unsigned i = 0; i < ax_list[j].spheres.size(); i++)
        {
            float x = ax_list[j].spheres[i].center[0];
            float y = ax_list[j].spheres[i].center[1];
            float z = ax_list[j].spheres[i].center[2];
            float radius = ax_list[j].spheres[i].radius;
            drawSphere(x, y, z, radius, colours[ax_list[j].id]);
        }
    }
}

void AxonGammaDistribution::growthVisualisation()
{
    radii = std::vector<double>(num_obstacles, 0);
    // generate radii with gamma distribution
    generate_radii(radii);
    std::cout << "Parallel growth simulation" << endl;
    bool stop = false;

    // Initialize SFML window
    sf::Window window(sf::VideoMode(800, 600), "3D Visualization");
    window.setActive();

    // Initialize OpenGL
    sf::ContextSettings settings;
    settings.depthBits = 24; // Request a 24-bit depth buffer

    // Initialize OpenGL settings
    initializeGLUT(0, nullptr);
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

    while (!stop)
    {
        num_obstacles = radii.size();
        std::vector<int> num_subsets; // depending on which batch
        setBatches(num_obstacles, num_subsets);
        drawBatches(window, axons, radii, num_subsets,zoomLevel,isDragging, lastMousePos, isRightDragging, lastRightMousePos, currentMousePos, mouseDelta, prevousDisplacement, rightMouseDelta, prevousRotation, rotationFactor,displacementFactor, false);
        std::vector<double> radii_to_regrow = stuck_radii;
        int num_regrowth = 0;

        while (num_regrowth < regrow_thr && radii_to_regrow.size() > 0)
        {
            //std::cout << "regrowth num :" << num_regrowth << endl;
            regrow_count = radii_to_regrow.size(); // nbr of axons to regrow
            setBatches(regrow_count, num_subsets);
            double prev_stuck_radii = stuck_radii.size();
            drawBatches(window, growing_axons, radii_to_regrow, num_subsets,zoomLevel,isDragging, lastMousePos, isRightDragging, lastRightMousePos, currentMousePos, mouseDelta, prevousDisplacement, rightMouseDelta, prevousRotation, rotationFactor,displacementFactor, true);
            radii_to_regrow = stuck_radii; 
            if (prev_stuck_radii != radii_to_regrow.size()){
                num_regrowth = -1;
            }    
            ++num_regrowth;
        }

        stop = true;
    }


    //FinalCheck(axons);
    bool swellaxons = swellAxons();

    icvf = computeICVF();
    // messages
    std::cout << "icvf: " << icvf << " voxel size: " << max_limits[0] << std::endl;
    std::string message = "ICVF achieved: " + std::to_string(icvf * 100);
    std::cout << message << std::endl;
    std::cout << "Number grown axons: " << axons.size() << std::endl;
    std::cout << "Total number axons: " << num_obstacles << std::endl;
    std::cout << "Number batches: " << num_batches << std::endl;
}


void set_Visualisation_params(sf::Window &window, float &zoomLevel, bool &isDragging, sf::Vector2i &lastMousePos, bool &isRightDragging, sf::Vector2i &lastRightMousePos, sf::Vector2i &currentMousePos, sf::Vector2i &mouseDelta, sf::Vector2i &prevousDisplacement, sf::Vector2i &rightMouseDelta, sf::Vector2i &prevousRotation, float &rotationFactor, float &displacementFactor){
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
    } 
} 

void AxonGammaDistribution::update_window(sf::Window &window, float zoomLevel, sf::Vector2i prevousDisplacement,  sf::Vector2i prevousRotation, float rotationFactor, float displacementFactor){
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

    drawWorld(axons); // draw one sphere at a time, j is the row number

    window.display(); // Display the updated window
} 

void AxonGammaDistribution::drawBatches(sf::Window &window, std::vector<Axon> &ax_list, std::vector<double> radii_, std::vector<int> num_subsets_,float zoomLevel,bool isDragging, sf::Vector2i lastMousePos, bool isRightDragging, sf::Vector2i lastRightMousePos, sf::Vector2i currentMousePos, sf::Vector2i mouseDelta, sf::Vector2i prevousDisplacement, sf::Vector2i rightMouseDelta, sf::Vector2i prevousRotation, float rotationFactor, float displacementFactor, bool regrowth)
{
    stuck_radii.clear(); // start batches with no stuck radii

    for (unsigned j = 0; j < num_batches; j++) // batches of axon growth
    {
        std::cout << "---   Batch " << j << "   --- " << endl;
        display_progress(j, num_batches);
        int first_index_batch = j * axon_capacity;

        int stuck = 0;
        vector<int> finished(num_subsets_[j], 0);         // 0 for false
        vector<int> grow_straight(num_subsets_[j], 1);    // 1 for true
        vector<int> straight_growths(num_subsets_[j], 0); // for each axon
        vector<int> shrink_tries(num_subsets_[j], 0);     // for each axon
        vector<int> restart_tries(num_subsets_[j], 0);    // for each axon

        bool all_finished = false;
        bool batch_created;

        // create batch of N axons, update growing axons
        createBatch(radii_, num_subsets_[j], regrowth, first_index_batch, growing_axons);

        int nbr_growing_axons = growing_axons.size();
        std::vector<Axon> growing_axons_copy(growing_axons);

        std::vector<std::thread> all_threads;

        for (unsigned i = 0; i < nbr_growing_axons; i++) // for each axon
        {
            set_Visualisation_params(window, zoomLevel, isDragging, lastMousePos, isRightDragging, lastRightMousePos, currentMousePos, mouseDelta, prevousDisplacement, rightMouseDelta, prevousRotation, rotationFactor, displacementFactor);
            update_window(window, zoomLevel, prevousDisplacement,  prevousRotation, rotationFactor, displacementFactor);
        
            all_threads.emplace_back(&AxonGammaDistribution::growAxon, 
                                    this,
                                    std::ref(growing_axons_copy[i]),
                                    regrowth);

        } // end for axons

        for (std::thread &t : all_threads)
        {
            if (t.joinable()){
                t.join();
            } 
        }
        set_Visualisation_params(window, zoomLevel, isDragging, lastMousePos, isRightDragging, lastRightMousePos, currentMousePos, mouseDelta, prevousDisplacement, rightMouseDelta, prevousRotation, rotationFactor, displacementFactor);
        update_window(window, zoomLevel, prevousDisplacement,  prevousRotation, rotationFactor, displacementFactor);
        
        all_threads.clear(); // Vider le vecteur pour l'itération suivante

        growing_axons = growing_axons_copy;
        if (axon_capacity > 1){
            SanityCheck(growing_axons);
        }

        // add axons in batch that worked in axons
        for (unsigned i = 0; i < growing_axons.size(); i++) {
            //cout << "growing_axons["<<i <<"].spheres.size(): " << growing_axons[i].spheres.size() << endl;
            if (growing_axons[i].spheres.size()>0){
                axons.push_back(growing_axons[i]);
            }
        }


        growing_axons.clear();
        
        //std::cout << "Stuck axons: " << stuck_radii.size() << endl;

    } // end for batches


}

// Growing substrate
void AxonGammaDistribution::parallelGrowth()
{

    radii = std::vector<double>(num_obstacles, 0);
    // generate radii with gamma distribution
    generate_radii(radii);
    std::cout << "Parallel growth simulation" << endl;

    bool stop = false;

    while (!stop)
    {
   
        num_obstacles = radii.size();
        std::vector<int> num_subsets; // depending on which batch
        cout << "radii size =" << num_obstacles << endl;
        setBatches(num_obstacles, num_subsets);
        growBatches(radii, num_subsets, false); // grows all axons, fills list of stuck radii, deletes empty axons
        std::vector<double> radii_to_regrow = stuck_radii;


        int num_regrowth = 0;

        while (num_regrowth < regrow_thr && radii_to_regrow.size() > 0)
        {
 
            //std::cout << "regrowth num :" << num_regrowth << endl;
            regrow_count = radii_to_regrow.size(); // nbr of axons to regrow
            setBatches(regrow_count, num_subsets);
            double prev_stuck_radii = stuck_radii.size();
            growBatches(radii_to_regrow, num_subsets, true);                // grows stuck axons, empties and refills stuck_radii, deletes empty axons
            radii_to_regrow = stuck_radii; 
            if (prev_stuck_radii != radii_to_regrow.size()){
                num_regrowth = -1;
            }                                             
            ++num_regrowth;

        }

        stop = true;
    }

    //FinalCheck(axons);
    bool swellaxons = swellAxons();
    //if (swellaxons){
    //    FinalCheck(axons);
    //}

    // messages
    icvf = computeICVF();
    std::cout << "icvf: " << icvf << " voxel size: " << max_limits[0] << std::endl;
    std::string message = "ICVF achieved: " + std::to_string(icvf * 100);
    std::cout << message << std::endl;
    std::cout << "Number grown axons: " << axons.size() << std::endl;
    std::cout << "Total number axons: " << num_obstacles << std::endl;
    std::cout << "Number batches: " << num_batches << std::endl;

}

bool AxonGammaDistribution::FinalCheck(std::vector<Axon>& axs){

    std::vector<Axon> final_axons;
    bool not_collide;
    for (unsigned j = 0; j < axs.size(); j++) {
        bool all_spheres_can_be_placed = true;
        for (unsigned i = 0; i < axs[j].spheres.size(); i++) { // for all spheres
            if (final_axons.size() > 0 && !canSpherebePlaced(axs[j].spheres[i], final_axons)){
                //std::cout << " Axon :" << axs[j].id << ", sphere : " << axs[j].spheres[i].id << " collides with environment !" << endl;
                all_spheres_can_be_placed = false;
                break;
            }
        }
        if (all_spheres_can_be_placed){
            final_axons.push_back(axs[j]);
        }
    }
    if (final_axons.size() == axs.size()){
        std::cout << " No Axon collides with environment !" << endl;
        not_collide = true;
    }
    else{
        not_collide = false;
        axs.clear();
        axs = final_axons;
    }
                
    return not_collide;
}

bool AxonGammaDistribution::SanityCheck(std::vector<Axon>& growing_axons) {

    // Vector to store results from each thread
    std::vector<double> thread_results(axon_capacity, -1);

    // Function to be executed by each thread
    auto thread_func = [this, &growing_axons, &thread_results](unsigned thread_id, unsigned num_threads) {
        for (unsigned j = thread_id; j < growing_axons.size(); j += num_threads) {
            thread_results[thread_id] = -1;
            if (growing_axons[j].spheres.size()>0){
                for (unsigned k = 0; k < growing_axons[j].spheres.size(); k ++) {
                    Sphere last_sphere = growing_axons[j].spheres[k];
                    if (!canSpherebePlaced(last_sphere, growing_axons)) {
                        thread_results[thread_id] = growing_axons[j].id;
                        //std::cout << "Sanity check: axon " << growing_axons[j].id << " collides with environment" << std::endl;
                        //std::cout << " axon size : " << growing_axons[j].spheres.size() << endl;
      
                        break;
                    }
                }
            }
        }
    };
    
    // Create threads and run the function
    std::vector<std::thread> threads;
    for (unsigned i = 0; i < axon_capacity; ++i) {
        threads.emplace_back(thread_func, i, axon_capacity);
    }

    // Wait for all threads to finish
    for (auto& thread : threads) {
        thread.join();
    }
    
    // Check results from all threads
    bool not_collide = true; 
    bool clear_one_axon = false;
    for (double result : thread_results) {
        if (result >= 0) {
            not_collide = false;
            int position;
            auto foundObject = std::find_if(growing_axons.begin(), growing_axons.end(), [result](const Axon& ax) {
                return ax.id == result;
            });

            if (foundObject != growing_axons.end()) {
                // Calculate the position of the found object in the vector
                position = std::distance(growing_axons.begin(), foundObject);
            }
            else{
                //cout << " did not find id :" <<result<<endl; 
                assert(0);
            }
            growing_axons[position].destroy();
            stuck_radii.push_back(growing_axons[position].radius);
        }
    }
    

    return not_collide;
}

/*
void AxonGammaDistribution::growBatches(std::vector<double> &radii_, std::vector<int> &num_subsets_, bool regrowth)
{
    stuck_radii.clear();
    std::vector<std::thread> all_threads;
    
    for (unsigned j = 0; j < num_batches; j++) // batches of axon growth
    {
        auto startTime = std::chrono::high_resolution_clock::now();
        std::cout << "---   Batch " << j << "   --- " << endl;

        display_progress(j, num_batches);
        int first_index_batch = j * axon_capacity;

        int stuck = 0;
        vector<int> finished(num_subsets_[j], 0);         // 0 for false
        vector<int> grow_straight(num_subsets_[j], 1);    // 1 for true
        vector<int> straight_growths(num_subsets_[j], 0); // for each axon
        vector<int> shrink_tries(num_subsets_[j], 0);     // for each axon
        vector<int> restart_tries(num_subsets_[j], 0);    // for each axon

        bool all_finished = false;
        bool batch_created;
        // create batch of N axons, update growing axons
        createBatch(radii_, num_subsets_[j], regrowth, first_index_batch, growing_axons);

        std::cout << "growing_axons.size() : " << growing_axons.size() << endl;
        std::cout << "axons.size() : " << axons.size() << endl;

        int nbr_growing_axons = growing_axons.size();


        while (!all_finished) // for each sphere
        {

            std::vector<Axon> growing_axons_copy(growing_axons);

            for (unsigned i = 0; i < nbr_growing_axons; i++) // for each axon
            {
    
                
                if (finished[i] == 0) // if the axon is not done growing
                {
                    if (growing_axons_copy[i].spheres.size() <= 1) // we cannot go straight if there are no first spheres as reference
                    {
                        grow_straight[i] = 0; // false
                    }
                    if (growing_axons_copy[i].spheres.size() >0){
                        // add sphere at same time for all axons :
      
                        all_threads.emplace_back(&AxonGammaDistribution::growthThread, // updates stuck_radii
                                                this,
                                                std::ref(growing_axons_copy[i]),
                                                std::ref(finished[i]),
                                                std::ref(grow_straight[i]),
                                                std::ref(straight_growths[i]),
                                                std::ref(restart_tries[i]),
                                                std::ref(shrink_tries[i]),
                                                regrowth);
                    }
                    else{
                        finished[i] = 1;
                    }
                }

            } // end for axons

            for (auto &thread : all_threads)
            {
                thread.join();
            }

            all_threads.clear(); // Vider le vecteur pour l'itération suivante

            growing_axons = growing_axons_copy;
            if (axon_capacity > 1){
                SanityCheck(growing_axons);
            }

            all_finished = std::all_of(finished.begin(), finished.end(), [](bool value)
                                       { return value == 1; }); // if all true, then all axons are finished growing
        
        } // end for spheres

        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime);
        cout << " time for batch = " << duration.count() << " seconds" << endl;
        

        // add axons in batch that worked in axons
        for (unsigned i = 0; i < growing_axons.size(); i++) {
            //cout << "growing_axons["<<i <<"].spheres.size(): " << growing_axons[i].spheres.size() << endl;
            if (growing_axons[i].spheres.size()>0){
                axons.push_back(growing_axons[i]);
            }
        }
        growing_axons.clear();
        
        std::cout << "Stuck axons: " << stuck_radii.size() << endl;

    } // end for batches

}
*/

void AxonGammaDistribution::growAxon(Axon& axon_to_grow, bool regrowth){

    int finished = 0;
    int grow_straight = 0;
    int straight_growths = 0;

    Growth growth = Growth(axon_to_grow, axons, max_limits, tortuous, max_radius, grow_straight);
    

    while (finished == 0) // if the axon is not done growing
    {
        if (axon_to_grow.spheres.size() >0){

            growthThread(axon_to_grow, growth, finished, grow_straight, straight_growths, regrowth);
        }
        else{
            finished = 1;
        }
    }
}

void AxonGammaDistribution::growBatches(std::vector<double> &radii_, std::vector<int> &num_subsets_, bool regrowth)
{
    stuck_radii.clear();
    
    for (unsigned j = 0; j < num_batches; j++) // batches of axon growth
    {
        auto startTime = std::chrono::high_resolution_clock::now();
        std::cout << "---   Batch " << j << "   --- " << endl;

        display_progress(j, num_batches);
        int first_index_batch = j * axon_capacity;

        int stuck = 0;
        vector<int> finished(num_subsets_[j], 0);         // 0 for false
        vector<int> grow_straight(num_subsets_[j], 1);    // 1 for true
        vector<int> straight_growths(num_subsets_[j], 0); // for each axon
        vector<int> shrink_tries(num_subsets_[j], 0);     // for each axon
        vector<int> restart_tries(num_subsets_[j], 0);    // for each axon

        bool all_finished = false;
        bool batch_created;
        // create batch of N axons, update growing axons
        createBatch(radii_, num_subsets_[j], regrowth, first_index_batch, growing_axons);

        int nbr_growing_axons = growing_axons.size();
        std::vector<Axon> growing_axons_copy(growing_axons);

        std::vector<std::thread> all_threads;
        for (unsigned i = 0; i < nbr_growing_axons; i++) // for each axon
        {

            all_threads.emplace_back(&AxonGammaDistribution::growAxon, 
                                    this,
                                    std::ref(growing_axons_copy[i]),
                                    regrowth);

        } // end for axons

        for (std::thread &t : all_threads)
        {
            if (t.joinable()){
                t.join();
            } 
        }

        all_threads.clear(); // Vider le vecteur pour l'itération suivante

        growing_axons = growing_axons_copy;
        if (axon_capacity > 1){
            SanityCheck(growing_axons);
        }

        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime);
        cout << " time for batch = " << duration.count() << " seconds" << endl;
        

        // add axons in batch that worked in axons
        for (unsigned i = 0; i < growing_axons.size(); i++) {
            //cout << "growing_axons["<<i <<"].spheres.size(): " << growing_axons[i].spheres.size() << endl;
            if (growing_axons[i].spheres.size()>0){
                axons.push_back(growing_axons[i]);
            }
        }
        growing_axons.clear();
        
        //std::cout << "Stuck axons: " << stuck_radii.size() << endl;

    } // end for batches

}


double get_axonal_length(Axon axon){
    double l = 0;
    if (axon.spheres.size()>1){
        for (unsigned i = 1; i < axon.spheres.size(); i++){
            double dist = (axon.spheres[i-1].center-axon.spheres[i].center).norm();
            l += dist;
        } 
        return l;
    }
    else{
        return 0;
    }
} 
// Axon growth
double AxonGammaDistribution::radiusVariation(Axon axon)
{
    double initial_radius = axon.radius;

    double length = get_axonal_length(axon);
    
    double p = variation_perc;
    double amplitude = initial_radius*(1 - p)/2;
    double mean_radius = (p+1) * initial_radius/2;

    double angular_frequency = M_PI / (initial_radius*3); // Adjust the frequency to control the sinusoidal curve
    double r = amplitude * cos(angular_frequency * length) + mean_radius;
    if (r < min_radius)
    {
        r = min_radius;
    }

    return r;
}



bool AxonGammaDistribution::shrinkRadius(Growth& growth, double radius_to_shrink, Axon &axon)
{

    bool can_grow;
    Eigen::Vector3d position_that_worked;
    // find position that works for smallest radius 
    can_grow = false;
    double rad = radius_to_shrink;
    double initial_rad = radius_to_shrink;
    bool create_sphere = true;
    bool can_grow_min_rad = growth.GrowAxon(min_radius, false);
    double intervals = (initial_rad-min_radius)/5;

    if (!can_grow_min_rad){
        //cout << "minimum radius doesn't fit, min rad :"<< min_radius << endl;
        return false;
    }
    else{

        while (!can_grow && rad > min_radius){
            rad -= intervals;
            can_grow = growth.GrowAxon(rad, create_sphere);
        }
        if (can_grow){
            axon = growth.axon_to_grow;
            return true;
        }
        else{
            cout << "cannot shrink after dichotomy" << endl;
            return false;
        }
    }

}

void update_straight(bool can_grow, int &grow_straight,int &straight_growths){

    if (can_grow){

        if (grow_straight == 1)
        {
            if (straight_growths >= 4) // if axon has been growing straight for 4 spheres in a row
            {
                grow_straight = 0; // set to false so that next step doesn't go straight
                straight_growths = 0;
            }
            else{
                straight_growths += 1;
            }
        }
        else
        {
            // if the sphere hadn't grown straight previously . set to straight for next 4 spheres
            grow_straight = 1; // set to true
        }
    }
    else{
        if (grow_straight==1) // if when growing straight it collides with environment
        {
            grow_straight = 0; // set to false so that next step doesn't go straight
            straight_growths = 0;
        }
    }
}

void AxonGammaDistribution::growthThread(Axon &axon, Growth &growth,  int &finished, int &grow_straight, int &straight_growths,  bool regrowth)
{

    double varied_radius = radiusVariation(axon);
    bool can_grow;
     
    //{
        //std::lock_guard<std::mutex> lock(stuckMutex);
    can_grow = growth.GrowAxon(varied_radius, true);
    axon = growth.axon_to_grow; // updates axon list
    //}
    

    if (growth.finished)
    {
        finished = 1; // 1 for true, if the growth finished
    }
    else // still growing
    {
        if (!can_grow)
        {
            
            bool shrink;
            shrink = shrinkRadius(growth, varied_radius, axon); // adds a sphere if it works
            if (!shrink) // shrinking will not help growth
            {
                //std::cout << "Axon " << axon.id << " failed, cannot shrink !" << endl;
                axon.destroy();
                finished = 1;
                stuck_radii.push_back(axon.radius);
            }
            else{
                // successfully placed a sphere after shrinking it
                can_grow = true;
                if (growth.finished)
                {
                    finished = 1; // 1 for true, if the growth finished
                }
            }
        }
    }
    update_straight(can_grow, grow_straight, straight_growths);
}


double volumeFrustumCone(double r1,double r2, double h){
    return M_PI*h*(r1*r1+r2*r2+r1*r2)/3;
}

double AxonGammaDistribution::computeICVF()
{
    if (axons.size() == 0)
        return 0;
    double VolumeV = (max_limits[0] - min_limits[0]) * (max_limits[1] - min_limits[1]) * (max_limits[2] - min_limits[2]); // total volume
    double VolumeC = 0;
    //double tortuosity;
    for (uint i = 0; i < axons.size(); i++) // for all axons
    {
        double ax_length = 0;
        if (axons[i].spheres.size() > 1)
        {
            for (uint j = 1; j < axons[i].spheres.size(); j++)
            {
                if (axons[i].spheres[j].center[2] <= max_limits[2]){
                    double l = (axons[i].spheres[j - 1].center - axons[i].spheres[j].center).norm(); // distance between centers
                    double r1= axons[i].spheres[j - 1].radius;
                    double r2= axons[i].spheres[j].radius;
                    double area_cone =  volumeFrustumCone(r1, r2, l);
                
                    if (withinBounds(axons[i].spheres[j].center, axons[i].spheres[j].radius) && withinBounds(axons[i].spheres[j-1].center, axons[i].spheres[j-1].radius))
                    {
                        VolumeC += area_cone ;
                    }
                    else if (withinBounds(axons[i].spheres[j].center, 0) || withinBounds(axons[i].spheres[j-1].center, 0))
                    {
                        VolumeC += area_cone/2;
                    }
    
                    ax_length += l;
                }
            }
            //tortuosity = ax_length / ((axons[i].begin - axons[i].end).norm()); // ( total distance between all centers / distance between first and last )
            //tortuosities.push_back(tortuosity);
        }
    }
    return VolumeC / VolumeV; // ( total axons volume / total volume )
}


void AxonGammaDistribution ::create_SWC_file(std::ostream &out, int overlapping_factor)
{
    std::vector<Axon> final_axons;

    fill_with_overlapping_spheres(overlapping_factor,final_axons);

    out << "id_ax id_sph Type X Y Z R P" << endl;
    int previous_sph_id;


    std::sort(final_axons.begin(), final_axons.end(), [](const Axon a, Axon b) -> bool
              { return a.radius < b.radius; }); // sort by size

    for (uint i = 0; i < final_axons.size(); i++)
    {
        previous_sph_id = -1; // parent for first sphere of each axon
        for (uint j = 0; j < final_axons[i].spheres.size(); j++)
        {
            out << i << " " << j << " axon " << final_axons[i].spheres[j].center[0] << " " << final_axons[i].spheres[j].center[1] << " " << final_axons[i].spheres[j].center[2] << " " << final_axons[i].spheres[j].radius << " " << previous_sph_id << endl;
            previous_sph_id = j;
        }
    }
}

void AxonGammaDistribution::simulation_file(std::ostream &out, std::chrono::seconds duration)
{
    out << "Duration " << duration.count() << endl;
    out << "Num_axons " << axons.size() << endl;
    out << "Capacity " << axon_capacity << endl;
    out << "Voxel " << max_limits[0] << endl;
    out << "icvf " << icvf << endl;
    out << "Density " << computeICVF() << endl;
}

std::vector<Eigen::Vector3d> equallySpacedPoints(const Eigen::Vector3d& point1, const Eigen::Vector3d& point2, int n) {
    std::vector<Eigen::Vector3d> result;

    // Calculate the step size for each dimension
    double stepX = (point2[0] - point1[0]) / (n + 1);
    double stepY = (point2[1]  - point1[1] ) / (n + 1);
    double stepZ = (point2[2]  - point1[2] ) / (n + 1);

    // Generate the equally spaced points
    for (int i = 1; i <= n; ++i) {
        Eigen::Vector3d newPoint;
        newPoint[0]  = point1[0]  + i * stepX;
        newPoint[1]  = point1[1]  + i * stepY;
        newPoint[2]  = point1[2]  + i * stepZ;
        result.push_back(newPoint);
    }

    return result;
}

std::vector<double> equallySpacedValues(double start, double end, int n) {
    std::vector<double> result;

    // Calculate the step size
    double step = (end - start) / (n + 1);

    // Generate the equally spaced values
    for (int i = 1; i <= n; ++i) {
        double newValue = start + i * step;
        result.push_back(newValue);
    }

    return result;
}

// Function to calculate intermediate points and insert N lines
void AxonGammaDistribution::fill_with_overlapping_spheres(int overlapping_factor, std::vector<Axon> &final_axons) {
    if (overlapping_factor > 2){
        int N = (overlapping_factor/2)-1; 
        for (uint k = 0; k < axons.size(); k++)
        {
            std::vector<Sphere> result;
            for (size_t i = 0; i < axons[k].spheres.size() - 1; ++i) {
                result.push_back(axons[k].spheres[i]); // Add the current sphere
                
                for (int j = 1; j <= N; ++j) {
                    double factor = static_cast<double>(j) / (N + 1);
                    Sphere intermediateSphere;
                    intermediateSphere.center[0] = axons[k].spheres[i].center[0] + factor * (axons[k].spheres[i + 1].center[0] - axons[k].spheres[i].center[0]);
                    intermediateSphere.center[1] = axons[k].spheres[i].center[1] + factor * (axons[k].spheres[i + 1].center[1] - axons[k].spheres[i].center[1]);
                    intermediateSphere.center[2] = axons[k].spheres[i].center[2] + factor * (axons[k].spheres[i + 1].center[2] - axons[k].spheres[i].center[2]);
                    intermediateSphere.radius = axons[k].spheres[i].radius + factor * (axons[k].spheres[i + 1].radius - axons[k].spheres[i].radius);
                    intermediateSphere.ax_id = axons[k].spheres[i].ax_id;
                    intermediateSphere.id = i * N + j; // Adjust ID for intermediate spheres
                    if (canSpherebePlaced(intermediateSphere, axons)){
                        result.push_back(intermediateSphere);
                    }
                }
            }
            result.push_back(axons[k].spheres.back()); // Add the last sphere
            Axon final_axon = Axon (axons[k].id, axons[k].begin, axons[k].end, axons[k].radius);
            final_axon.spheres = result;
            final_axons.push_back(final_axon);
        }
    }
    else{
        final_axons = axons;
    }
}
/*
void AxonGammaDistribution::fill_with_overlapping_spheres(int overlapping_factor, std::vector<Axon> &final_axons){

    std::vector<Axon> axons_env = axons;

    int nbr_to_insert =(overlapping_factor/2)-1;

    Axon final_axon;

    if (overlapping_factor > 2){ 
        std::cout << "Increasing overlap with factor : "<< overlapping_factor << endl;

        for (uint i = 0; i < axons_env.size(); i++)
        {
            final_axon = Axon (axons_env[i].id, axons_env[i].begin, axons_env[i].end, axons_env[i].radius);
            int number_spheres = axons_env[i].spheres.size();
            // add first sphere
            final_axon.add_sphere(axons_env[i].spheres[0]);

            for (uint j = 1; j < number_spheres; j++)
            {
                std::vector<Eigen::Vector3d> centers = equallySpacedPoints(axons_env[i].spheres[j-1].center, axons_env[i].spheres[j].center, nbr_to_insert);

                std::vector<double> radii_interp = equallySpacedValues(axons_env[i].spheres[j-1].radius, axons_env[i].spheres[j].radius, nbr_to_insert);
                
                for (uint k = 0; k < radii_interp.size(); k++)
                {
                    Sphere sph_interp = Sphere(final_axon.spheres.size(), axons_env[i].id, centers[k], radii_interp[k]);

                    bool can_place = canSpherebePlaced(sph_interp, axons_env);
                    if (can_place){
                        final_axon.add_sphere(sph_interp);
                    } 

                } 
                // add last sphere of batch
                final_axon.add_sphere(axons_env[i].spheres[j]);
            }
            final_axons.push_back(final_axon);
            axons_env[i] = final_axon;
        }
    } 
    else{
        final_axons = axons;
    }
} 

*/