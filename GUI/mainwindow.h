#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QSlider>
#include <QLabel>
#include <QString>
#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QPushButton>
#include <QtDataVisualization/Q3DScatter>
#include "slidergroup.h"
#include "ScatterDataModifier.h" // Include ScatterDataModifier
#include <QComboBox>  

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
    void updateConfigurationSelectionVisibility(double value);

    
private slots:
    void onSaveButtonClicked();
    void onSelectDirectoryButtonClicked(); // Slot for selecting a directory
    void PlotCells();

private:
    void createControls(const QString &title);
    void resizeEvent(QResizeEvent *e);
    void StartSimulation();
    SlidersGroup *slidersGroup;

    QComboBox *configurationComboBox;
    QGroupBox *controlsGroup;
    QGridLayout *generalLayout;
    QGridLayout *axonsLayout;
    QGridLayout *glialLayout;
    QGroupBox *generalGroup;
    QGroupBox *axonsGroup;
    QGroupBox *glialGroup;

    QLabel *visualise_voxel_qlabel;
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
    QLabel *nbr_axons_populations_qlabel;
    QLabel *c2_qlabel;
    QLabel *epsilon_qlabel;
    QLabel *mean_process_length_qlabel;
    QLabel *std_process_length_qlabel;
    QLabel *beading_amplitude_qlabel;
    QLabel *alpha_qlabel;
    QLabel *beta_qlabel;

    QCheckBox *visualise_voxel_checkbox;
    QDoubleSpinBox *axons_icvf_SpinBox;
    QDoubleSpinBox *axons_w_myelin_icvf_SpinBox;
    QDoubleSpinBox *astrocyte_soma_icvf_SpinBox;
    QDoubleSpinBox *astrocyte_processes_icvf_SpinBox;
    QDoubleSpinBox *oligodendrocyte_soma_icvf_SpinBox;
    QDoubleSpinBox *oligodendrocyte_processes_icvf_SpinBox;
    QDoubleSpinBox *voxel_size_SpinBox;
    QDoubleSpinBox *minimum_radius_SpinBox;
    QDoubleSpinBox *nbr_threads_SpinBox;
    QDoubleSpinBox *overlapping_factor_SpinBox;
    QDoubleSpinBox *nbr_axons_populations_SpinBox;
    QDoubleSpinBox *c2_SpinBox;
    QDoubleSpinBox *epsilon_SpinBox;
    QDoubleSpinBox *mean_process_length_SpinBox;
    QDoubleSpinBox *std_process_length_SpinBox;
    QDoubleSpinBox *beading_amplitude_SpinBox;
    QDoubleSpinBox *alpha_SpinBox;
    QDoubleSpinBox *beta_SpinBox;


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
    int nbr_axons_populations;
    int crossing_fibers_type;
    double mean_process_length;
    double std_process_length;
    double beading_amplitude;
    double epsilon;
    double alpha;
    double beta;
    bool visualise_voxel;


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
