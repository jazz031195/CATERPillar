#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QSlider>
#include <QLabel>
#include <QString>
#include <QCheckBox>
#include <QSpinBox>
#include <QPushButton>
#include <QtDataVisualization/Q3DScatter>
#include "slidergroup.h"
#include "ScatterDataModifier.h" // Include ScatterDataModifier

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class Window : public QWidget
{
    Q_OBJECT

public:
    Window(QWidget *parent = nullptr);
    void createStatisticsMenu();
    void plotRadiusDistribution();
    void plotTortuosityDistribution();
    void ShollAnalysis();
    void resetCamera();
    
private slots:
    void onSaveButtonClicked();
    void onSelectDirectoryButtonClicked(); // Slot for selecting a directory
    void PlotCells();

private:
    void createControls(const QString &title);
    void resizeEvent(QResizeEvent *e);
    void StartSimulation();
    SlidersGroup *slidersGroup;

    QGroupBox *controlsGroup;
    QLabel *axons_icvf_qlabel;
    QLabel *axons_w_myelin_icvf_qlabel;
    QLabel *astrocyte_soma_icvf_qlabel;
    QLabel *astrocyte_processes_icvf_qlabel;
    QLabel *oligodendrocyte_soma_icvf_qlabel;
    QLabel *oligodendrocyte_processes_icvf_qlabel;
    QLabel *voxel_size_qlabel;
    QLabel *minimum_radius_qlabel;
    QLabel *nbr_threads_qlabel;
    QLabel *overlapping_factor_qlabel;
    QLabel *c2_qlabel;

    QCheckBox *Tortuous;
    QCheckBox *Beading;
    QSpinBox *axons_icvf_SpinBox;
    QSpinBox *axons_w_myelin_icvf_SpinBox;
    QSpinBox *astrocyte_soma_icvf_SpinBox;
    QSpinBox *astrocyte_processes_icvf_SpinBox;
    QSpinBox *oligodendrocyte_soma_icvf_SpinBox;
    QSpinBox *oligodendrocyte_processes_icvf_SpinBox;
    QSpinBox *voxel_size_SpinBox;
    QSpinBox *minimum_radius_SpinBox;
    QSpinBox *nbr_threads_SpinBox;
    QSpinBox *overlapping_factor_SpinBox;
    QSpinBox *c2_SpinBox;

    QBoxLayout *layout;
    QPushButton *okButton;
    QPushButton *selectDirectoryButton; // Button to select directory
    QString selectedDirectory; // String to store the selected directory

    QPushButton *statisticsButton;
    QPushButton *plotRadiusDistributionButton;
    QPushButton *plotTortuosityDistributionButton;
    QPushButton *plotShollAnalysisButton;
    QPushButton *resetCameraButton;

    // Additional member variables to store values
    int axons_icvf;
    int axons_w_myelin_icvf;
    int astrocyte_soma_icvf;
    int astrocyte_processes_icvf;
    int oligodendrocyte_soma_icvf;
    int oligodendrocyte_processes_icvf;
    int voxel_size;
    double minimum_radius;
    int nbr_threads;
    int overlapping_factor;
    double c2;

    // spheres to plot
    std::vector<std::vector<double>> X_axons;
    std::vector<std::vector<double>> Y_axons;
    std::vector<std::vector<double>> Z_axons;
    std::vector<std::vector<double>> R_axons;

    std::vector<std::vector<double>> X_astrocytes;
    std::vector<std::vector<double>> Y_astrocytes;
    std::vector<std::vector<double>> Z_astrocytes;
    std::vector<std::vector<double>> R_astrocytes;
    std::vector<std::vector<int>> Branch_astrocytes;

    // Member for the 3D scatter plot and modifier
    QtDataVisualization::Q3DScatter *graph;
    ScatterDataModifier *modifier;
    OpenGLWindow *openglWindow;
};

#endif // MAINWINDOW_H
