#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "Eigen/Core"
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
    void HideGlialCells();
    void HideAxons();
    void ShowAllCells();
    bool check_borders(const Eigen::Vector3d&  min_l, const Eigen::Vector3d&  max_l, const Eigen::Vector3d& pos, const double& distance_to_border);
    void StartSimulation();
    void loadConfigFromPath(const QString &fileName);
    void updateMembersFromGui();

private slots:
    void onSaveButtonClicked();
    void onSelectDirectoryButtonClicked(); // Slot for selecting a directory
    void SelectSWCFileButton();
    void PlotCells(const bool& axons_plot, const bool& glial_pop1_plot, const bool& glial_pop2_plot, const bool& blood_vessels_plot);
    void ReadGlialCellsFromSWC(const QString& filePath);
    void ReadAxonsFromSWC(const QString& filePath);
    void ReadAxonsFromCSV(const QString& fileName);
    void ReadAxonsFromFile(const QString& fileName);
    void ReadGlialCellsFromFile(const QString& fileName);
    void ReadGlialCellsFromCSV(const QString& fileName);
    void ReadBloodVesselsFromFile(const QString& fileName);
    void loadConfigFromFile();

private:
    QGroupBox* createControls(const QString &title);
    void resizeEvent(QResizeEvent *e);
    void saveConfigToPath(const QString &fileName);

    SlidersGroup *slidersGroup;

    QComboBox *configurationComboBox;
    QGroupBox *controlsGroup;
    QGridLayout *generalLayout;
    QGridLayout *axonsLayout;
    QGridLayout *glialLayout;
    QGroupBox *generalGroup;
    QGroupBox *axonsGroup;
    QGroupBox *glialGroup1;
    QGroupBox *glialGroup2;

    QLabel *nbr_repetitions_qlabel;
    QLabel *visualise_voxel_qlabel;
    QLabel *axons_icvf_qlabel;
    QLabel *axons_w_myelin_icvf_qlabel;
    QLabel *glial_pop1_soma_icvf_qlabel;
    QLabel *glial_pop1_processes_icvf_qlabel;
    QLabel *glial_pop2_soma_icvf_qlabel;
    QLabel *glial_pop2_processes_icvf_qlabel;
    QLabel *blood_vessels_icvf_qlabel;
    QLabel *voxel_size_qlabel;
    QLabel *minimum_radius_qlabel;
    QLabel *nbr_threads_qlabel;
    QLabel *overlapping_factor_qlabel;
    QLabel *nbr_axons_populations_qlabel;
    QLabel *c2_qlabel;
    QLabel *epsilon_qlabel;
    QLabel *glial_pop1_mean_process_length_qlabel;
    QLabel *glial_pop1_std_process_length_qlabel;
    QLabel *glial_pop2_mean_process_length_qlabel;
    QLabel *glial_pop2_std_process_length_qlabel;
    QLabel *beading_amplitude_qlabel;
    QLabel *beading_std_qlabel;
    QLabel *alpha_qlabel;
    QLabel *beta_qlabel;
    QLabel *glial_pop1_radius_mean_qlabel;
    QLabel *glial_pop1_radius_std_qlabel;
    QLabel *glial_pop2_radius_mean_qlabel;
    QLabel *glial_pop2_radius_std_qlabel;
    QLabel *glial_pop1_nbr_primary_processes_qlabel;
    QLabel *glial_pop2_nbr_primary_processes_qlabel;
    QLabel *glial_pop1_branching_qlabel;
    QLabel *glial_pop2_branching_qlabel;
    QLabel *k1_qlabel;
    QLabel *k2_qlabel;
    QLabel *k3_qlabel;

    QDoubleSpinBox *nbr_repetitions_SpinBox;
    QCheckBox *visualise_voxel_checkbox;
    QCheckBox *glial_pop1_branching_checkbox;
    QCheckBox *glial_pop2_branching_checkbox;
    QDoubleSpinBox *axons_icvf_SpinBox;
    QDoubleSpinBox *axons_w_myelin_icvf_SpinBox;
    QDoubleSpinBox *glial_pop1_soma_icvf_SpinBox;
    QDoubleSpinBox *glial_pop1_processes_icvf_SpinBox;
    QDoubleSpinBox *glial_pop2_soma_icvf_SpinBox;
    QDoubleSpinBox *glial_pop2_processes_icvf_SpinBox;
    QDoubleSpinBox *blood_vessels_icvf_SpinBox;
    QDoubleSpinBox *voxel_size_SpinBox;
    QDoubleSpinBox *minimum_radius_SpinBox;
    QDoubleSpinBox *nbr_threads_SpinBox;
    QDoubleSpinBox *overlapping_factor_SpinBox;
    QDoubleSpinBox *nbr_axons_populations_SpinBox;
    QDoubleSpinBox *c2_SpinBox;
    QDoubleSpinBox *epsilon_SpinBox;
    QDoubleSpinBox *glial_pop1_mean_process_length_SpinBox;
    QDoubleSpinBox *glial_pop1_std_process_length_SpinBox;
    QDoubleSpinBox *glial_pop2_mean_process_length_SpinBox;
    QDoubleSpinBox *glial_pop2_std_process_length_SpinBox;
    QDoubleSpinBox *beading_amplitude_SpinBox;
    QDoubleSpinBox *beading_std_SpinBox;
    QDoubleSpinBox *alpha_SpinBox;
    QDoubleSpinBox *beta_SpinBox;
    QDoubleSpinBox *glial_pop1_radius_mean_SpinBox;
    QDoubleSpinBox *glial_pop1_radius_std_SpinBox;
    QDoubleSpinBox *glial_pop2_radius_mean_SpinBox;
    QDoubleSpinBox *glial_pop2_radius_std_SpinBox;
    QDoubleSpinBox *glial_pop1_nbr_primary_processes_SpinBox;
    QDoubleSpinBox *glial_pop2_nbr_primary_processes_SpinBox;
    QDoubleSpinBox *k1_SpinBox;
    QDoubleSpinBox *k2_SpinBox;
    QDoubleSpinBox *k3_SpinBox;


    QBoxLayout *layout;
    QPushButton *okButton;
    QPushButton *selectDirectoryButton; // Button to select directory
    QString selectedDirectory; // String to store the selected directory
    QString SWCFile;

    QPushButton *statisticsButton;
    QPushButton *plotRadiusDistributionButton;
    QPushButton *plotTortuosityDistributionButton;
    QPushButton *plotShollAnalysisButton;
    QPushButton *resetCameraButton;
    QPushButton *visualiseButton;
    QPushButton *growButton;
    QPushButton *hideGlialCellsButton;
    QPushButton *hideAxonsButton;
    QPushButton *showAllCellsButton;

    QPushButton *loadConfigButton;    


    // Additional member variables to store values
    int nbr_repetitions;
    int axons_icvf;
    int axons_w_myelin_icvf;
    int glial_pop1_soma_icvf;
    int glial_pop1_processes_icvf;
    int glial_pop2_soma_icvf;
    int glial_pop2_processes_icvf;
    int blood_vessels_icvf;
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
    double beading_std;
    double epsilon;
    double alpha;
    double beta;
    bool visualise_voxel;
    double glial_pop1_radius_mean;
    double glial_pop1_radius_std;
    double glial_pop2_radius_mean;
    double glial_pop2_radius_std;
    int glial_pop1_nbr_primary_processes;
    int glial_pop2_nbr_primary_processes;
    bool glial_pop1_branching;
    bool glial_pop2_branching;
    double glial_pop1_mean_process_length;
    double glial_pop1_std_process_length;
    double glial_pop1_mean_primary_process_length;
    double glial_pop1_std_primary_process_length;
    double glial_pop2_mean_process_length;
    double glial_pop2_std_process_length;
    double glial_pop2_mean_primary_process_length;
    double glial_pop2_std_primary_process_length;
    int max_nbr_process_pop1;
    int max_nbr_process_pop2;
    double k1;
    double k2;
    double k3;


    // spheres to plot
    std::vector<std::vector<double>> X_axons;
    std::vector<std::vector<double>> Y_axons;
    std::vector<std::vector<double>> Z_axons;
    std::vector<std::vector<double>> R_axons;

    std::vector<std::vector<double>> X_glial_pop1;
    std::vector<std::vector<double>> Y_glial_pop1;
    std::vector<std::vector<double>> Z_glial_pop1;
    std::vector<std::vector<double>> R_glial_pop1;
    std::vector<std::vector<int>> Branch_glial_pop1;

    std::vector<std::vector<double>> X_glial_pop2;
    std::vector<std::vector<double>> Y_glial_pop2;
    std::vector<std::vector<double>> Z_glial_pop2;
    std::vector<std::vector<double>> R_glial_pop2;
    std::vector<std::vector<int>> Branch_glial_pop2;

    std::vector<std::vector<double>> X_blood_vessels;
    std::vector<std::vector<double>> Y_blood_vessels;
    std::vector<std::vector<double>> Z_blood_vessels;
    std::vector<std::vector<double>> R_blood_vessels;


    // Member for the 3D scatter plot and modifier
    QtDataVisualization::Q3DScatter *graph;
    ScatterDataModifier *modifier;
    OpenGLWindow *openglWindow;
    
};

#endif // MAINWINDOW_H
