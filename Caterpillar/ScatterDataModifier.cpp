#include "ScatterDataModifier.h"
#include <QtDataVisualization/Q3DTheme>
#include <QtDataVisualization/QAbstract3DGraph>
#include <QtDataVisualization/Q3DCamera>
#include <QtDataVisualization/QScatterDataProxy>
#include <QtDataVisualization/QScatter3DSeries>
#include <QtDataVisualization/QScatterDataArray>
#include <QtDataVisualization/QScatterDataItem>
#include <QComboBox>
#include <QtMath>
#include <vector>
#include <map>
#include <algorithm> // For std::max_element
#include <cstdlib> // For std::rand and std::srand
#include <random>
#include <omp.h>
#include <QApplication>



ScatterDataModifier::ScatterDataModifier(OpenGLWindow* window_)
    : window_(window_)
{
    //m_graph->activeTheme()->setType(QtDataVisualization::Q3DTheme::ThemeEbony);
    QtDataVisualization::QScatterDataProxy *proxy = new QtDataVisualization::QScatterDataProxy;

    

}

void ScatterDataModifier::setVoxelParameters(const int &voxel_size_, const std::vector<std::vector<double>> &X_axons_, const std::vector<std::vector<double>> &Y_axons_, const std::vector<std::vector<double>> &Z_axons_, const std::vector<std::vector<double>> &R_axons_, const std::vector<std::vector<double>> &X_astrocytes_, const std::vector<std::vector<double>> &Y_astrocytes_, const std::vector<std::vector<double>> &Z_astrocytes_, const std::vector<std::vector<double>> &R_astrocytes_)
{
    voxel_size = voxel_size_;
    X_axons = X_axons_;
    Y_axons = Y_axons_;
    Z_axons = Z_axons_;
    R_axons = R_axons_;
    X_astrocytes = X_astrocytes_;
    Y_astrocytes = Y_astrocytes_;
    Z_astrocytes = Z_astrocytes_;
    R_astrocytes = R_astrocytes_;
}


void ScatterDataModifier::generateRandomColors(int count, std::vector<QColor> &colors) {
    // Seed the random number generator with current time
    std::mt19937 rng(static_cast<unsigned int>(std::time(nullptr)));

    // Distribution for RGB values (0 to 255)
    std::uniform_int_distribution<int> dist(0, 255);

    // Generate random colors
    for (int i = 0; i < count; ++i) {
        int red = dist(rng);
        int green = dist(rng);
        int blue = dist(rng);
        colors.push_back(QColor(red, green, blue));
    }
}
void ScatterDataModifier::updateData() {
    if (!window_) {  // Check if window_ is a valid pointer
        qDebug() << "OpenGLWindow is not initialized!";
        return;
    }

    // Ensure that the window title and size are set before displaying
    window_->setTitle("3D Spheres Visualization");
    window_->resize(800, 600);

    // Pass the spheres data to the OpenGL window and trigger an update
    window_->setSpheres(X_axons, Y_axons, Z_axons, R_axons);

    // This will ensure the window is updated with new data
    window_->update();
}

void ScatterDataModifier::get_window(OpenGLWindow* &window)
{
    window = window_;
}

 