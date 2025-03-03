#ifndef SCATTERDATAMODIFIER_H
#define SCATTERDATAMODIFIER_H

#include <QtDataVisualization/Q3DScatter>
#include <QtDataVisualization/QScatterDataProxy>
#include <QtDataVisualization/QScatter3DSeries>
#include <QtDataVisualization/Q3DTheme>
#include <QtDataVisualization/QAbstract3DGraph>
#include <QtDataVisualization/Q3DCamera>
#include <QObject>
#include <QFont>
#include "openglwindow.h"

class ScatterDataModifier : public QObject
{
    Q_OBJECT

public:
    ScatterDataModifier(OpenGLWindow* window_);
    void updateData();
    void generateRandomColors(int count, std::vector<QColor> &colors);
    void get_window(OpenGLWindow* &window);

public slots:
    void setVoxelParameters(const int &voxel_size_, const std::vector<std::vector<double>> &X_axons_, const std::vector<std::vector<double>> &Y_axons_, const std::vector<std::vector<double>> &Z_axons_, const std::vector<std::vector<double>> &R_axons_, const std::vector<std::vector<double>> &X_astrocytes_, const std::vector<std::vector<double>> &Y_astrocytes_, const std::vector<std::vector<double>> &Z_astrocytes_, const std::vector<std::vector<double>> &R_astrocytes_);


private:
    //QtDataVisualization::Q3DScatter *m_graph;
    OpenGLWindow* window_;
    int voxel_size;
    std::vector<std::vector<double>> X_axons;
    std::vector<std::vector<double>> Y_axons;
    std::vector<std::vector<double>> Z_axons;
    std::vector<std::vector<double>> R_axons;
    std::vector<std::vector<double>> X_astrocytes;
    std::vector<std::vector<double>> Y_astrocytes;
    std::vector<std::vector<double>> Z_astrocytes;
    std::vector<std::vector<double>> R_astrocytes;

    


};

#endif // SCATTERDATAMODIFIER_H
