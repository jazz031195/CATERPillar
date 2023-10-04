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

AxonGammaDistribution::AxonGammaDistribution(unsigned &num_ax, int &axon_capacity_, double a, double b,
                                             Eigen::Vector3d &min_l, Eigen::Vector3d &max_l, double min_radius_,
                                             bool tortuous_, bool draw_, int regrow_thr_, bool can_shrink_, double beading_variation)
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
    regrow_thr = regrow_thr_;
    can_shrink=can_shrink_;
    variation_perc = beading_variation;

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
void AxonGammaDistribution::set_icvf(double icvf_, double x, double y) // overwrite num_obstacles
{
    icvf = icvf_;
    double av_radius = 0.5;
    num_obstacles = (x * y * icvf) / (M_PI * av_radius * av_radius);
    std::cout << num_obstacles << endl;
}
void AxonGammaDistribution::generate_radii(std::vector<double> &radiis)
{

    std::random_device rd;
    std::default_random_engine generator(rd());
    std::gamma_distribution<double> distribution(alpha, beta);
    radiis.clear();

    // take into account variation of radius -> increase initial icvf to compensate
    double icvf_to_reach = icvf/((variation_perc+1)/2);

    int tried = 0;

    double icvf_ = 0;

    double AreaIntra = 0;
    double AreaTot = max_limits[0]*max_limits[1];

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

    num_obstacles = radii.size();

    std::cout << "Number of axons :" << num_obstacles << endl;
}

double AxonGammaDistribution::generate_radius()

{
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::gamma_distribution<double> distribution(alpha, beta);
    double radius = 0;
    while (radius < min_radius){
        radius = distribution(generator);
    }
    return radius;
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
        
        //area_circle = squareCircleOverlap(max_limits[0] - min_limits[0], axon.spheres[0].radius, axon.spheres[0].center[0], axon.spheres[0].center[1]);
        //AreaIntra += area_circle;
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
    Dynamic_Sphere sphere = Dynamic_Sphere(0, ax.id, Q, radius);
    
    
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



bool AxonGammaDistribution::canSpherebePlaced(Dynamic_Sphere sph, std::vector<Axon> axs, bool print){

    int count = 0;
    for (auto &axon : axs)
    {

        if (count > 1){
            assert(0);
        }

        if (axon.id != sph.ax_id){ 
            if (axon.isSphereInsideAxon_(sph)) // overlap
            {
                /*
                if (!axon.isSphereInsideAxon_long(sph)){
                    cout << "isSphereInsideAxon_ not working in canSpherebePlaced !" << endl;
                    assert(0);
                }
                */
                if (print){
                    std::cout << "collides with axon : " << axon.id << endl;
                } 
                
                return false;
            }
            /*
            else{
                if (axon.isSphereInsideAxon_long(sph)){
                    cout << "isSphereInsideAxon_ not working in canSpherebePlaced !" << endl;
                    assert(0);
                }
            }
            */
        } 
        else{
            count += 1;
        }

    }
    return true;
} 



double AxonGammaDistribution::dichotomy_swelling(Dynamic_Sphere new_sphere, int index, double swelling_perc){

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
void AxonGammaDistribution::swellAxons(){
    
    for (uint i = 0; i < axons.size(); i++) // for all axons
    {
        std::cout << " axon : " << axons[i].id << endl;
        for (uint j = 0; j < axons[i].spheres.size(); j++)
        {
            Dynamic_Sphere new_sphere = axons[i].spheres[j];

            double final_rad = dichotomy_swelling(new_sphere, i, radii_swelling);
            axons[i].spheres[j].radius = final_rad;
  
        }
    }
    double obtained_icvf = computeICVF();

    std::cout << "Obtained icvf :" << obtained_icvf << endl;
    std::cout << "Target icvf :" << icvf << endl;

    // if obtained is too low
    if (obtained_icvf + 0.05 < 0.6){
        double new_icvf = obtained_icvf;
        double difference_icvf = icvf - obtained_icvf;
        int tries = 0;
        while (difference_icvf> 0.5 && tries <10){
            for (uint i = 0; i < axons.size(); i++) // for all axons
            {
                std::cout << " axon : " << axons[i].id << endl;
                for (uint j = 0; j < axons[i].spheres.size(); j++)
                {
                    Dynamic_Sphere new_sphere = axons[i].spheres[j];

                    double final_rad = dichotomy_swelling(new_sphere, i, 0.1);
                    axons[i].spheres[j].radius = final_rad;
        
                }
            }
            tries += 1;
            new_icvf = computeICVF();
            std::cout << "New icvf :" << obtained_icvf <<", try :"<< tries << endl;
            difference_icvf = icvf - new_icvf;
        }
    }
}

std::vector<Eigen::Vector3d> AxonGammaDistribution::FindTwins(Eigen::Vector3d Q, double rad){
    std::vector<Eigen::Vector3d> new_Qs;
    Eigen::Vector3d new_Q;
    new_Qs.push_back(Q);
    // overlaps with x = 0
    if (abs(Q[0]-min_limits[0])< rad){

        new_Q = Q + Eigen::Vector3d {(max_limits[0]-min_limits[0]), 0, 0};
        new_Qs.push_back(new_Q);
    }
    // overlaps with x = max
    if (abs(Q[0]-max_limits[0])< rad){
        new_Q = Q + Eigen::Vector3d {-(max_limits[0]-min_limits[0]), 0, 0};
        new_Qs.push_back(new_Q);
    }
    // overlaps with y = 0
    if (abs(Q[1]-min_limits[1])< rad){

        new_Q = Q + Eigen::Vector3d {0,(max_limits[1]-min_limits[1]), 0};
        new_Qs.push_back(new_Q);
    }
    // overlaps with y = max
    if (abs(Q[1]-max_limits[1])< rad){
        new_Q = Q + Eigen::Vector3d {0, -(max_limits[1]-min_limits[1]), 0};
        new_Qs.push_back(new_Q);
    }
    // add twin at diagonal 
    if (new_Qs.size()== 3){
        // overlaps with x = 0 and y = 0
        if (abs(Q[0]-min_limits[0])< rad && abs(Q[1]-min_limits[1])< rad){

            new_Q = Q + Eigen::Vector3d {(max_limits[1]-min_limits[1]),(max_limits[1]-min_limits[1]), 0};
            new_Qs.push_back(new_Q);
        }
        // overlaps with x = 0 and y = max
        else if (abs(Q[0]-min_limits[0])< rad && abs(Q[1]-max_limits[1])< rad){
            new_Q = Q + Eigen::Vector3d {+(max_limits[1]-min_limits[1]),-(max_limits[1]-min_limits[1]), 0};
            new_Qs.push_back(new_Q);
        }
        // overlaps with y = max and x = max
        else if (abs(Q[0]-max_limits[0])< rad && abs(Q[1]-max_limits[1])< rad){

            new_Q = Q + Eigen::Vector3d {-(max_limits[0]-min_limits[0]),-(max_limits[0]-min_limits[0]), 0};
            new_Qs.push_back(new_Q);
        }
        // overlaps with y = 0 and x = max
        else if (abs(Q[0]-max_limits[0])< rad && abs(Q[1]-min_limits[1])< rad){
            new_Q = Q + Eigen::Vector3d {-(max_limits[0]-min_limits[0]),(max_limits[0]-min_limits[0]), 0};
            new_Qs.push_back(new_Q);
        }
    }

    return new_Qs;
}

void AxonGammaDistribution::PlaceTwinAxons(double radius, bool regrowth, int i, int &tries, bool &next, std::vector<Eigen::Vector3d> Qs, std::vector<Axon> &all_axons, std::vector<Axon> &new_axons){
    
    std::vector<bool> managed;
    for (unsigned j = 0; j < Qs.size(); ++j) {
        Axon ax;
        Eigen::Vector3d D = Qs[j];
        D[2] = max_limits[2];

        if (!regrowth)
        {
            ax = Axon(i, Qs[j], D, radius); // original axons
        }
        else
        {
            ax = Axon(i + axons.size(),Qs[j], D, radius); // axons for regrow batch
        }
        Dynamic_Sphere sphere = Dynamic_Sphere(0, ax.id, Qs[j], radius);
        ax.add_sphere(sphere);

        if (all_axons.empty())
        {
            managed.push_back(true);
        }
        else // not empty
        {
            bool overlap = false;

            for (auto &axon : all_axons)
            {
                if (axon.isSphereInsideAxon_(sphere)) // overlap
                {
                    managed.push_back(false);
                    overlap = true;
                    break;

                }
            }
            if (overlap == false) // after comparing with all axons
            {
                managed.push_back(true);

            }
        }
    }
    // Check if all elements in the vector are true
    bool allTrue = std::all_of(managed.begin(), managed.end(), [](bool element) {
        return element;
    });
    if (allTrue){
        for (unsigned j = 0; j < Qs.size(); ++j) {
            Axon ax;
            GLfloat new_colour = generateRandomColor();
            Eigen::Vector3d D = Qs[j];
            D[2] = max_limits[2];

            if (!regrowth)
            {
                ax = Axon(i, Qs[j], D, radius); // original axons
            }
            else
            {
                ax = Axon(i + axons.size(),Qs[j], D, radius); // axons for regrow batch
            }
            Dynamic_Sphere sphere = Dynamic_Sphere(0, ax.id, Qs[j], radius);
            ax.add_sphere(sphere);
            new_axons.push_back(ax);
            all_axons.push_back(ax);
            colours.push_back(new_colour);
        }
        tries = 0; // reset number of tries;
        next = true;
    }
    else{
        ++tries;
        next = false;
    }
}
void AxonGammaDistribution::createAxons(std::vector<double> &radii_, std::vector<Axon> &new_axons, bool regrowth)
{
    int tries = 0;
    std::vector<int> stuck_axons_id;
    std::vector<Axon> all_axons = axons;
    new_axons.clear();
    long int tries_threshold = 100000;


    for (unsigned i = 0; i < radii_.size(); ++i) // create axon for each radius
    {
        std::cout << "radii created :" << i << "/" << radii_.size() << endl;
        //cout  << "all axons size :" <<  all_axons.size() << endl;
        bool next = false;
        while (!next && tries < tries_threshold)
        {

            Vector3d Q, D;
            get_begin_end_point(Q, D);

            PlaceAxon(radii_[i], regrowth, i, tries, next, Q, D, all_axons, new_axons);

        }
        if (tries >= tries_threshold)
        {
            std::cout << "Not enough space, restart" << endl;
            std::cout << "achieved :" << new_axons.size() << "/" << radii_.size() << endl;
            new_axons.clear();
            all_axons = axons;
            colours.clear();
            i = -1;
            tries = 0;
        }
    }
    radii_.clear();
    //cout << "All axons size :" << all_axons.size() << endl;
    //cout << "Axons size :" << axons.size() << endl;
    //cout << "New Axons to grow id :"<< endl;
    for (unsigned j = 0; j < new_axons.size(); ++j)
    {
        //cout << "   " << new_axons[j].id << endl;
        new_axons[j].nearby_axons.clear();
        if (new_axons[j].spheres.size()>0){
            for (unsigned k = 0; k < all_axons.size(); ++k) {
                if (all_axons[k].isNearAxon(new_axons[j].spheres[0].center, 10)){
                    if (all_axons[k].id != new_axons[j].id){
                        new_axons[j].add_nearby_axon(all_axons[k].id);
                    }
                    //else{
                    //    cout << "same id" << endl;
                    //    cout << "all_axons[k].radius :" << all_axons[k].radius << endl;
                    //    cout << "new_axons[j].radius :" << new_axons[j].radius << endl;
                    //}
                    //if (all_axons[k].spheres.size() == 0){
                    //    cout << "axon " << all_axons[k].id << " has no spheres" << endl;
                    //}
                }
                //else{
                //    cout << "NOT CLOSE " << endl;
                //}
            }
        }
        if (regrowth){
            std::cout << "axon " << new_axons[j].id << "has " << new_axons[j].nearby_axons.size()<< " nearby axons" << endl;
        }
        //cout << "nearby axons :" << new_axons[j].nearby_axons.size()<< endl;
        radii_.push_back(new_axons[j].radius);
    }
    
    if(!regrowth){
        double icvf_ = computeAreaICVF(new_axons);
        std::cout << "icvf achieved :" << icvf_ << endl;
    }
    
}

void AxonGammaDistribution::createBatch(std::vector<double> &radii_, int num_subset, bool regrowth, int first_index_batch, std::vector<Axon> &new_axons)
{
    int tries = 0;
    std::vector<int> stuck_axons_id;
    long int tries_threshold = 100000;
    std::vector<Axon> all_axons = axons;
    new_axons.clear();

    std::vector<double>::iterator startIterator = radii_.begin() + first_index_batch;  // Start from index 
    std::vector<double>::iterator stopIterator = radii_.begin() + first_index_batch + num_subset;  // Stop when batch is finished
    
    // Create a new vector using the iterators
    std::vector<double> batch_radii(startIterator, stopIterator);

    int placement_tries_thresh = 10;
    int placement_tries = 0;


    for (unsigned i = 0; i < batch_radii.size(); ++i) // create axon for each radius
    {
        std::cout << "radii created :" << i << "/" << batch_radii.size() << endl;
        //cout  << "all axons size :" <<  all_axons.size() << endl;
        bool next = false;
        
        while (!next && tries < tries_threshold)
        {

            Vector3d Q, D;
            get_begin_end_point(Q, D);
            int axon_index = first_index_batch + i;
            //if (withinBounds(Q, -radii_[i])){
            PlaceAxon(batch_radii[i], regrowth, axon_index, tries, next, Q, D, all_axons, new_axons);

            
            //}
            //else{
                //cout << "create twins, position : "<< Q<< endl;
                //cout << "radius : "<< radii_[i]<< endl;
                //cout << "new_axons.size before :" << new_axons.size() << endl;
                //std::vector<Vector3d> twin_Qs = FindTwins(Q,radii_[i]);
                //PlaceTwinAxons(radii_[i], regrowth, i, tries, next, twin_Qs, all_axons, new_axons);
                //cout << "new_axons.size after :" << new_axons.size() << endl;
                //cout << "next :" << next << endl;
            //}
        }
        if (tries >= tries_threshold && placement_tries < placement_tries_thresh)
        {
            //if (!regrowth) // creating all axons
            //{
            std::cout << "Not enough space, restart" << endl;
            std::cout << "achieved :" << new_axons.size() << "/" << batch_radii.size() << endl;
            new_axons.clear();
            all_axons = axons;
            colours.clear();
            i = -1;
            tries = 0;
            placement_tries += 1;
        }

            //else // regrowing some axons
            //{
            //    std::cout << "Not enough space, forget about axon " << i << endl;
            //    stuck_axons_id.push_back(i);
            //    tries = 0;
            //}
        //}
    }
    //radii_.clear();

    for (unsigned j = 0; j < new_axons.size(); ++j)
    {
        //std::cout << "   " << new_axons[j].id << endl;
        new_axons[j].nearby_axons.clear();
        if (new_axons[j].spheres.size()>0){
            for (unsigned k = 0; k < all_axons.size(); ++k) {
                if (all_axons[k].spheres.size()>0 && all_axons[k].isNearAxon(new_axons[j].spheres[0].center, 5)){
                    if (all_axons[k].id != new_axons[j].id){
                        new_axons[j].add_nearby_axon(all_axons[k].id);
                    }
                    //else{
                    //    std::cout << "same id" << endl;
                    //    std::cout << "all_axons[k].radius :" << all_axons[k].radius << endl;
                    //    std::cout << "new_axons[j].radius :" << new_axons[j].radius << endl;
                    //}
                    //if (all_axons[k].spheres.size() == 0){
                    //    std::cout << "axon " << all_axons[k].id << " has no spheres" << endl;
                    //}
                }
                //else{
                //    std::cout << "NOT CLOSE " << endl;
                //}
            }
        }
        if (regrowth){
            std::cout << "axon " << new_axons[j].id << "has " << new_axons[j].nearby_axons.size()<< " nearby axons" << endl;
        }
        //std::cout << "nearby axons :" << new_axons[j].nearby_axons.size()<< endl;
        //radii_.push_back(new_axons[j].radius);
    }

    /*
    if(!regrowth){
        double icvf_ = computeAreaICVF(new_axons);
        std::cout << "icvf achieved :" << icvf_ << endl;
    }
    */
    
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
    std::cout << "num_subsets.size : " << num_subsets.size() << endl;
    std::cout << "num_axons : " << num_axons << endl;
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
void AxonGammaDistribution::drawWorld(std::vector<Axon> ax_list, unsigned int row, sf::Window &window, int num_ax) // for parallel growth
{
    int num_col = axon_capacity;

    for (unsigned int j = row * num_col; j < row * num_col + num_ax; ++j) // new axon batch
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

    for (unsigned int j = 0; j < row * num_col; ++j) // already created batches
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

    if (ax_list.size() != axons.size()) // regrowth
    {
        for (unsigned int j = 0; j < axons.size(); ++j)
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
}
void AxonGammaDistribution::growthVisualisation()
{
    radii = std::vector<double>(num_obstacles, 0);
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
        std::vector<Axon> new_axons;
        createAxons(radii, new_axons, false);
        axons = new_axons;
        num_obstacles = axons.size();
        std::vector<int> num_subsets;
        setBatches(num_obstacles, num_subsets);
        drawBatches(window, axons, radii, num_subsets,zoomLevel,isDragging, lastMousePos, isRightDragging, lastRightMousePos, currentMousePos, mouseDelta, prevousDisplacement, rightMouseDelta, prevousRotation, rotationFactor,displacementFactor, false);
        std::vector<double> radii_to_regrow = stuck_radii;

        int num_regrowth = 0;

        while (num_regrowth < regrow_thr && radii_to_regrow.size() > 0)
        {
            createAxons(radii_to_regrow, growing_axons, true);
            regrow_count = radii_to_regrow.size(); // nbr of axons to regrow
            setBatches(regrow_count, num_subsets);
            drawBatches(window, growing_axons, radii_to_regrow, num_subsets,zoomLevel,isDragging, lastMousePos, isRightDragging, lastRightMousePos, currentMousePos, mouseDelta, prevousDisplacement, rightMouseDelta, prevousRotation, rotationFactor,displacementFactor, true);
            radii_to_regrow = stuck_radii;
            axons.insert(axons.end(), growing_axons.begin(), growing_axons.end()); // add regrown axons to the axons list
            growing_axons.clear();
            ++num_regrowth;
        }

        stop = true;
    }


    // messages
    std::cout << "icvf: " << icvf << " voxel size: " << max_limits[0] << std::endl;
    std::string message = "ICVF achieved: " + std::to_string(icvf * 100);
    std::cout << message << std::endl;
    std::cout << "Number grown axons: " << axons.size() << std::endl;
    std::cout << "Total number axons: " << num_obstacles << std::endl;
    std::cout << "Number batches: " << num_batches << std::endl;
    std::cout << "compute icvf " << computeICVF() << "\n"
              << std::endl;
}
void AxonGammaDistribution::drawBatches(sf::Window &window, std::vector<Axon> &ax_list, std::vector<double> &radii_, std::vector<int> &num_subsets,float zoomLevel,bool isDragging, sf::Vector2i lastMousePos, bool isRightDragging, sf::Vector2i lastRightMousePos, sf::Vector2i currentMousePos, sf::Vector2i mouseDelta, sf::Vector2i prevousDisplacement, sf::Vector2i rightMouseDelta, sf::Vector2i prevousRotation, float rotationFactor, float displacementFactor, bool regrowth)
{


    stuck_radii.clear(); // start batches with no stuck radii

    for (unsigned j = 0; j < num_batches; j++) // batches of axon growth
    {
        std::cout << "---   Batch " << j << "   --- " << endl;
        display_progress(j, num_batches);

        int stuck = 0;
        vector<int> finished(num_subsets[j], 0);         // 0 for false
        vector<int> grow_straight(num_subsets[j], 1);    // 1 for true
        vector<int> straight_growths(num_subsets[j], 0); // for each axon
        vector<int> shrink_tries(num_subsets[j], 0);     // for each axon
        vector<int> restart_tries(num_subsets[j], 0);    // for each axon

        bool all_finished = false;

        while (!all_finished && window.isOpen()) // for each sphere
        {
            std::vector<std::thread> all_threads;

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

                    int index = j * axon_capacity + i; // position within axons from all batches

                    if (ax_list[index].spheres.size() <= 1) // we cannot go straight if there are no first spheres as reference
                    {
                        grow_straight[i] = 0; // false
                    }

                    all_threads.emplace_back(&AxonGammaDistribution::growthThread,
                                             this,
                                             std::ref(ax_list[index]),
                                             std::ref(finished[i]),
                                             std::ref(grow_straight[i]),
                                             std::ref(straight_growths[i]),
                                             std::ref(restart_tries[i]),
                                             std::ref(shrink_tries[i]), 
                                             regrowth);

                    stuck += 1;
                }

            } // end for axons

            for (auto &thread : all_threads)
            {
                thread.join();
            }

            all_threads.clear(); // Vider le vecteur pour l'itération suivante


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

            drawWorld(ax_list, j, window, num_subsets[j]); // draw one sphere at a time, j is the row number

            window.display(); // Display the updated window
            //std::this_thread::sleep_for(std::chrono::seconds(1));

            all_finished = std::all_of(finished.begin(), finished.end(), [](bool value)
                                       { return value == 1; }); // if all true, then all axons are finished growing

        } // end for spheres

        std::cout << "Stuck axons: " << stuck_radii.size() << endl;

    } // end for batches

    num_batches = 0;
    num_subsets.clear();
    ax_list.erase(std::remove_if(ax_list.begin(), ax_list.end(), [](const Axon &ax)
                                 { return ax.spheres.empty(); }),
                  ax_list.end());
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
        //createAxons(radii, axons, false);    // set the axons
        num_obstacles = radii.size();
        std::vector<int> num_subsets; // depending on which batch
        setBatches(num_obstacles, num_subsets);
        // updates axons and stuck_radii
        growBatches(radii, num_subsets, false); // grows all axons, fills list of stuck radii, deletes empty axons
        std::vector<double> radii_to_regrow = stuck_radii;


        int num_regrowth = 0;

        while (num_regrowth < regrow_thr && radii_to_regrow.size() > 0)
        {
        
            std::cout << "regrowth num :" << num_regrowth << endl;
            //createAxons(radii_to_regrow, axons_to_regrow, true);
            regrow_count = radii_to_regrow.size(); // nbr of axons to regrow
            setBatches(regrow_count, num_subsets);
            double prev_stuck_radii = stuck_radii.size();
            growBatches(radii_to_regrow, num_subsets, true);                // grows stuck axons, empties and refills stuck_radii, deletes empty axons
            radii_to_regrow = stuck_radii; 
            if (prev_stuck_radii != radii_to_regrow.size()){
                num_regrowth = -1;
            }                                             
            //axons.insert(axons.end(), axons_to_regrow.begin(), axons_to_regrow.end()); // add regrown axons to the axons list
            
            ++num_regrowth;

        }

        stop = true;
    }
    //std::std::cout << "Swelling axons " << std::endl;
    FinalCheck();

    // messages
    std::cout << "icvf: " << icvf << " voxel size: " << max_limits[0] << std::endl;
    std::string message = "ICVF achieved: " + std::to_string(icvf * 100);
    std::cout << message << std::endl;
    std::cout << "Number grown axons: " << axons.size() << std::endl;
    std::cout << "Total number axons: " << num_obstacles << std::endl;
    std::cout << "Number batches: " << num_batches << std::endl;
    std::cout << "compute icvf " << computeICVF() << "\n"
              << std::endl;
}

bool AxonGammaDistribution::FinalCheck(){
    for (unsigned j = 0; j < axons.size(); j++) {
        for (unsigned i = 0; i < axons[j].spheres.size(); i++) {
            if (!canSpherebePlaced(axons[j].spheres[i], axons)){
                std::cout << " Axon :" << axons[j].id << ", sphere : " << axons[j].spheres[i].id << " collides with environment !" << endl;
                
                double initial_radius = axons[j].spheres[i].radius;
                bool istoobig = true;
                int tries = 0;
                while (istoobig && axons[j].spheres[i].radius > min_radius){
                    axons[j].spheres[i].radius -= initial_radius*0.1;
                    istoobig = !canSpherebePlaced(axons[j].spheres[i], axons);
                }
                if (axons[j].spheres[i].radius <= min_radius){
                    std::cout << " Axon :" << axons[j].id << ", sphere : " << axons[j].spheres[i].id << " collides with environment and cannot shrink enough !" << endl;
                    return false;
                }
                else{
                    std::cout << " Axon :" << axons[j].id << ", sphere : " << axons[j].spheres[i].id << " was shrunk from "<<initial_radius <<" to "<<axons[j].spheres[i].radius <<  " !" << endl;
                    
                }
            
            }
        }
    }
    std::cout << " No Axon collides with environment !" << endl;
                
    return true;
}
void AxonGammaDistribution::growBatches(std::vector<double> &radii_, std::vector<int> &num_subsets_, bool regrowth)
{
    stuck_radii.clear();
    
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

        std::cout << "growing_axons.size() : " << growing_axons.size() << endl;
        std::cout << "axons.size() : " << axons.size() << endl;

        int nbr_growing_axons = growing_axons.size();


        while (!all_finished) // for each sphere
        {
            std::vector<std::thread> all_threads;

            for (unsigned i = 0; i < nbr_growing_axons; i++) // for each axon
            {
    
                
                if (finished[i] == 0) // if the axon is not done growing
                {
                    if (growing_axons[i].spheres.size() <= 1) // we cannot go straight if there are no first spheres as reference
                    {
                        grow_straight[i] = 0; // false
                    }
                    if (growing_axons[i].spheres.size() >0){
                        // add sphere at same time for all axons :
      
                        all_threads.emplace_back(&AxonGammaDistribution::growthThread, // updates stuck_radii
                                                this,
                                                std::ref(growing_axons[i]),
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


            all_finished = std::all_of(finished.begin(), finished.end(), [](bool value)
                                       { return value == 1; }); // if all true, then all axons are finished growing

        } // end for spheres

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

void AxonGammaDistribution::dichotomy(Growth growth, Eigen::Vector3d position_that_worked, Axon axon, double initial_rad,  double &last_rad, int grow_straight)
{
    int tries = 0;
    double max_rad= initial_rad;
    double min_rad = min_radius;
    while (tries < 10){

        double current_rad = (max_rad + min_rad) / 2;

        bool can_grow ;

        can_grow = growth.TestGrowAxonAtPos(position_that_worked, current_rad);
        
        if (can_grow) // solution is in greater half
        {
            min_rad = current_rad;
        }

        else // solution is in lower half
        {
            max_rad = current_rad;
        }
        ++tries;
   }
   last_rad = min_rad;
   /*
   if (last_rad > initial_rad){
        std::cout << "last rad after dichtomy:" << last_rad << endl;
        std::cout << "initial_rad after dichtomy:" << initial_rad << endl;
        std::cout << "last rad after dichtomy:" << last_rad << endl;
   }
   */


}
/*
bool AxonGammaDistribution::shrinkRadius(Growth growth, double radius_to_shrink, Axon &axon, int grow_straight)
{
    
    bool can_grow;
    Eigen::Vector3d position_that_worked;
    // find position that works for smallest radius 
    can_grow = growth.TestGrowAxon(position_that_worked, min_radius);

    int tries = 0;
    if (can_grow) // shrinking is useful
    {
        double last_rad;
        dichotomy(growth, position_that_worked, axon, radius_to_shrink, last_rad, grow_straight);
        Dynamic_Sphere s(axon.spheres.size(), axon.id, position_that_worked, last_rad);
        //std::cout << "axon : "<< axon.id << " sphere : "<< axon.spheres.size()<< "shrink to :" << last_rad << ", initial at : "<< initial_radius << endl;
        axon.add_sphere(s);
        std::cout << "Axon " << axon.id << "  done after shrink from " << radius_to_shrink << " to "<<last_rad << endl;
        return true;
    }
    else
    {
        // std::cout << "shrinking not useful" << endl;
        return false; // shrinking won't be useful
    }
}
*/


bool AxonGammaDistribution::shrinkRadius(Growth growth, double radius_to_shrink, Axon &axon, int grow_straight)
{

    bool can_grow;
    Eigen::Vector3d position_that_worked;
    // find position that works for smallest radius 
    can_grow = false;
    double rad = radius_to_shrink;
    while(!can_grow && rad > min_radius){
        rad = rad - rad*0.1;
        {
            std::lock_guard<std::mutex> lock(stuckMutex);
            growth = Growth(axon, axons,growing_axons, max_limits, tortuous, rad, max_radius, grow_straight, radii_swelling);
            can_grow = growth.GrowAxon();
        }
    }
    if (rad <= min_radius){
        return false;
    }
    else{
        return true;
    }

}



void AxonGammaDistribution::growthThread(Axon &axon, int &finished, int &grow_straight, int &straight_growths, int &shrink_tries, int &restart_tries, bool regrowth)
{

    double varied_radius = radiusVariation(axon);
    Growth growth;
    bool can_grow;
     


    {
        std::lock_guard<std::mutex> lock(stuckMutex);
        growth = Growth(axon, axons,growing_axons, max_limits, tortuous, varied_radius, max_radius, grow_straight, radii_swelling);
        can_grow = growth.GrowAxon();
    }

    
    if (growth.finished)
    {
        finished = 1; // 1 for true, if the growth finished
    }
    else // still growing
    {
        if (!can_grow)
        {
            if ((regrowth || can_shrink) && shrink_tries < axon.spheres.size()/5)
            {
                {
                    bool shrink;

                    shrink = shrinkRadius(growth, varied_radius, axon, grow_straight); // adds a sphere if it works
                    
                    if (!shrink) // shrinking will not help growth
                    {
                        /*
                        if (restart_tries < 5)
                        {
                            Dynamic_Sphere sphere = axon.spheres[0];
                            axon.spheres.clear();
                            //axon.projections.clear_projections();
                            axon.add_sphere(sphere);
                            ++restart_tries;
                        }
                        else
                        {
                        */
                        std::cout << "Axon " << axon.id << " failed, cannot shrink !" << endl;
                        axon.spheres.clear();
                        finished = 1;
                        {
                            std::lock_guard<std::mutex> lock(stuckMutex);
                            stuck_radii.push_back(axon.radius);
                        }
                        //}

                    }

                    
                }
                ++shrink_tries;
            }
            else
            {
                std::cout << "!! Axon " << axon.id << " : failed growing !!" << endl;
                axon.spheres.clear();
                finished = 1;
                {
                    std::lock_guard<std::mutex> lock(stuckMutex);
                    stuck_radii.push_back(axon.radius);
                }
                
            }

            if (grow_straight==1) // if when growing straight it collides with environment
            {
                grow_straight = 0; // set to false so that next step doesn't go straight
                straight_growths = 0;
            }
        }
        else
        {
            shrink_tries = 0;

            if (grow_straight == 1)
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

            axon = growth.axon_to_grow; // updates axon list
        }
    }
}


// Relevant information
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


double AxonGammaDistribution::squareCircleOverlap(double L, double circleRadius, double circleCenterX, double circleCenterY) {
    // Calculate the half side length of the square
    double halfL = L / 2.0;

    // Calculate the distance between the center of the circle and the square's center
    double dx = std::abs(circleCenterX - halfL);
    double dy = std::abs(circleCenterY - halfL);

    // Check if the circle is entirely outside the square
    if (dx > halfL + circleRadius || dy > halfL + circleRadius) {
        return 0.0;
    }

    // Check if the circle is entirely inside the square
    if (dx <= halfL && dy <= halfL && std::hypot(dx, dy) + circleRadius <= halfL) {
        return M_PI * circleRadius * circleRadius;
    }

    // Calculate the overlap area
    double cornerDist = std::hypot(dx - halfL, dy - halfL);
    if (cornerDist <= circleRadius) {
        // Circle center is inside the square
        double sectorAngle = 2.0 * std::acos(cornerDist / circleRadius);
        double sectorArea = 0.5 * circleRadius * circleRadius * (sectorAngle - std::sin(sectorAngle));
        double squareArea = 4.0 * halfL * halfL;
        return squareArea - sectorArea;
    } else {
        // Circle center is outside the square
        double squareArea = 4.0 * halfL * halfL;
        return squareArea;
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
        if (axons[i].spheres.size() > 1)
        {
            for (uint j = 1; j < axons[i].spheres.size(); j++)
            {
                double l = (axons[i].spheres[j - 1].center - axons[i].spheres[j].center).norm(); // distance between centers
                //double area_circle1 = squareCircleOverlap(max_limits[0] - min_limits[0], axons[i].spheres[j - 1].radius, axons[i].spheres[j - 1].center[0], axons[i].spheres[j - 1].center[1]);
                //double area_circle2 = squareCircleOverlap(max_limits[0] - min_limits[0], axons[i].spheres[j].radius, axons[i].spheres[j].center[0], axons[i].spheres[j].center[1]);
                //double mean_area = (area_circle1+area_circle2)/2;
            

                double mean_r = (axons[i].spheres[j - 1].radius + axons[i].spheres[j].radius) / 2;
                
                if (withinBounds(axons[i].spheres[j].center, axons[i].spheres[j].radius) && withinBounds(axons[i].spheres[j-1].center, axons[i].spheres[j-1].radius))
                {
                    AreaC += l * M_PI * mean_r * mean_r;
                }
                else if (withinBounds(axons[i].spheres[j].center, 0) && withinBounds(axons[i].spheres[j-1].center, 0))
                {
                    AreaC += l * M_PI * mean_r * mean_r/2;
                }
                
                //AreaC += l*mean_area;
                ax_length += l;
            }
            tortuosity = ax_length / ((axons[i].begin - axons[i].end).norm()); // ( total distance between all centers / distance between first and last )
            tortuosities.push_back(tortuosity);
        }
    }
    return AreaC / AreaV; // ( total axons volume / total volume )
}


double AxonGammaDistribution::computeICVF_(){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> udist(0, 1);

    int intra = 0;
    int num_tries = max_limits[2]*max_limits[1]*max_limits[2]*100000;
    std::cout << num_tries << endl;

    for (uint i = 0; i < num_tries; i++){

        double t = udist(gen);
        // double distance_to_border = radiis[i]*sqrt(1+volume_inc_perc) + step_length;
        double x = (t * (max_limits[0])) + (1 - t) * (min_limits[0]);
        t = udist(gen);
        double y = (t * (max_limits[1])) + (1 - t) * (min_limits[1]);
        t = udist(gen);
        double z = (t * (max_limits[2])) + (1 - t) * (min_limits[2]);

        Dynamic_Sphere s = Dynamic_Sphere(-1, -1, {x,y,z}, 0.01);

        for (auto &ax : axons)
        {
            if (ax.isSphereInsideAxon_(s)) // overlap
            {
                ++intra;
                break;
            }
        }
    }

    return intra/num_tries;


    
}

// For analysis
void AxonGammaDistribution::axonDensityMap()
{
    // Create an SFML window
    sf::RenderWindow window(sf::VideoMode(1000, 1000), "2D Visualization");

    // Define a list of axons, with predetermined radii and positions
    radii = std::vector<double>(num_obstacles, 0);
    generate_radii(radii);
    createAxons(radii, axons, false);

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
void AxonGammaDistribution ::create_SWC_file(std::ostream &out, int overlapping_factor)
{
    std::vector<Axon> final_axons;

    fill_wih_overlapping_spheres(overlapping_factor,final_axons);

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
void AxonGammaDistribution::axons_file(std::ostream &out)
{
    std::sort(axons.begin(), axons.end(), [](const Axon a, Axon b) -> bool
              { return a.radius < b.radius; }); // sort by size

    out << "Type ax_id Type2 X Y Z Distance R R0 Tortuosity" << endl;
    for (uint i = 0; i < axons.size(); i++)
    {
        for (uint j = 0; j < axons[i].spheres.size(); j++)
        {
            out << " axon " << i << " sphere " << axons[i].spheres[j].center[0] << " " << axons[i].spheres[j].center[1] << " " << axons[i].spheres[j].center[2]
                << " " << axons[i].spheres[j].center[2] - axons[i].begin[2] << " " << axons[i].spheres[j].radius << " " << axons[i].radius << " " << tortuosities[i] << endl;
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

void AxonGammaDistribution::fill_wih_overlapping_spheres(int overlapping_factor, std::vector<Axon> &final_axons){

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
                    Dynamic_Sphere sph_interp = Dynamic_Sphere(final_axon.spheres.size(), axons_env[i].id, centers[k], radii_interp[k]);
                    //std::cout << "Maybe add after sphere " <<j <<" in axon :"<< axons[i].id << endl; 
                    //std::cout << "center 1 :" << axons[i].spheres[j-1].center<< " radius :" <<axons[i].spheres[j-1].radius << endl;
                    //std::cout << "center 2 :" << axons[i].spheres[j].center <<" radius :" <<axons[i].spheres[j].radius<< endl;

                    bool can_place = canSpherebePlaced(sph_interp, axons_env);
                    if (can_place){

                        //std::cout << "       Adding sphere " <<j <<" in axon :"<< axons[i].id << ", center :" <<centers[k] << ", radius :" <<radii_interp[k] << endl; 
                        // add intermediate spheres
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