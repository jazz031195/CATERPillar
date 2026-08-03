#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "Eigen/Core"
#include <QProcess>
#include <QWidget>
#include <QTabWidget>
#include <QStackedWidget>
#include <QToolButton>
#include <QFormLayout>
#include <QLineEdit>
#include <QMainWindow>
#include <QSlider>
#include <QLabel>
#include <QString>
#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QSpinBox>
#include <QPushButton>
#include <QProgressBar>
#include <QtDataVisualization/Q3DScatter>
#include "slidergroup.h"
#include "ScatterDataModifier.h" // Include ScatterDataModifier
#include "../src/parameters.h"
#include "../src/core_logic.h"
#include <QComboBox>
#include <thread>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class Window : public QWidget
{
    Q_OBJECT

public:
    Window(QWidget *parent = nullptr);
    ~Window();
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
    /*!
     *  \brief Parses one data line of a CATERPillar substrate CSV, tolerating both the
     *         current 10-token format (with parent_component_id) and the older 9-token
     *         format that predates it. For the 9-token case, parent is set equal to
     *         component_id, matching how the writer already fills that column for
     *         objects with no distinct parent (e.g. axons).
     *  \return false if the line doesn't have exactly 9 or 10 whitespace-separated tokens.
     */
    bool parseSubstrateCsvLine(const std::string &line, std::string &type, double &id_cell,
                                std::string &component, double &component_id, double &parent,
                                double &x, double &y, double &z, double &radius_in, double &radius_out);

private slots:
    void onSaveButtonClicked();
    void onSelectDirectoryButtonClicked(); // Slot for selecting a directory
    void SelectSWCFileButton();
    void PlotCells(const bool& axons_plot, const bool& glial_pop1_plot, const bool& glial_pop2_plot, const bool& glial_pop3_plot, const bool& blood_vessels_plot);
    void ReadGlialCellsFromSWC(const QString& filePath);
    void ReadAxonsFromSWC(const QString& filePath);
    void ReadAxonsFromCSV(const QString& fileName);
    void ReadAxonsFromFile(const QString& fileName);
    void ReadGlialCellsFromFile(const QString& fileName);
    void ReadGlialCellsFromCSV(const QString& fileName);
    void ReadBloodVesselsFromFile(const QString& fileName);
    void runMCSimulation();
    void onGrowthProgress(double completed_depth, double total_depth);
    void onSwellingProgress(double current_icvf, double target_icvf);
    void onGrowthFinished();

private:
    void initParameters();
    void buildParameterStack(QVBoxLayout *wmLayout);
    QGroupBox* createControls(const QString &title);
    void resizeEvent(QResizeEvent *e);
    void StartSimulation();

    QProgressBar *layerProgressBar = nullptr;
    QProgressBar *swellingProgressBar = nullptr;
    std::thread growthThread;
    std::vector<Axon> pendingAxons;
    std::vector<Blood_Vessel> pendingBloodVessels;
    std::vector<Glial> pendingGlialPop1;
    std::vector<Glial> pendingGlialPop2;
    std::vector<Glial> pendingGlialPop3;
    Eigen::Vector3d pendingVoxelMin;
    Eigen::Vector3d pendingVoxelMax;
    // Real (small) voxel's actual bounds for the currently-displayed data --
    // defaults to an origin-anchored box of edge parameters.voxel_size (the
    // only option when the source is a loaded CSV/SWC file with no
    // accompanying growth_info.txt to read a real placement from), but is set
    // to the true grown placement after "Grow Substrate" (see
    // CaterpillarGrowth::PlaceSmallVoxel) or a successful growth_info.txt read.
    QVector3D voxelBoundsMin = QVector3D(0.0f, 0.0f, 0.0f);
    QVector3D voxelBoundsMax = QVector3D(0.0f, 0.0f, 0.0f);
    SlidersGroup *slidersGroup;
    OpenGLWindow *openglWindow = nullptr;
    QWidget *visualizationWidget = nullptr;

    QStackedWidget *cellParamsStack;

    QLineEdit *inputN;
    QLineEdit *inputT;
    QLineEdit *inputDuration;
    QLineEdit *inputDiffIntra;
    QLineEdit *inputDiffExtra;
    QLineEdit *inputSchemeFile;
    QLineEdit *inputCsvPath;
    QLineEdit *inputExecutablePath;
    QLineEdit *inputLoadConfigPath;
    QDoubleSpinBox *inputVoxelSizeMC;
    QSpinBox *inputNumThreadsMC;
    QCheckBox *checkIncludeAxons;
    QCheckBox *checkIncludeGlial;
    QCheckBox *checkIncludeBloodVessels;
    QCheckBox *checkIniWalkersIntra;
    QCheckBox *checkIniWalkersExtra;
    QProcess *simulatorProcess; 

    QComboBox *configurationComboBox;
    QGroupBox *controlsGroup;
    QGridLayout *generalLayout;
    QGridLayout *axonsLayout;
    QGridLayout *glialLayout;
    QGroupBox *generalGroup;
    QGroupBox *axonsGroup;
    QGroupBox *glialGroup1;
    QGroupBox *glialGroup2;
    QGroupBox *glialGroup3;
    QGroupBox *myelinatedGroup;
    QGroupBox *bloodVesselGroup;

    QLabel *nbr_repetitions_qlabel;
    QLabel *visualise_voxel_qlabel;
    QLabel *axons_icvf_qlabel;
    QLabel *axons_w_myelin_icvf_qlabel;
    QLabel *glial_pop1_soma_icvf_qlabel;
    QLabel *glial_pop1_processes_icvf_qlabel;
    QLabel *glial_pop2_soma_icvf_qlabel;
    QLabel *glial_pop2_processes_icvf_qlabel;
    QLabel *glial_pop3_soma_icvf_qlabel;
    QLabel *glial_pop3_processes_icvf_qlabel;
    QLabel *blood_vessels_icvf_qlabel;
    QLabel *blood_vessels_processes_icvf_qlabel;
    QLabel *blood_vessel_voxel_size_qlabel;
    QLabel *blood_vessel_mean_radius_qlabel;
    QLabel *blood_vessel_std_radius_qlabel;
    QLabel *blood_vessel_capillary_radius_qlabel;
    QLabel *blood_vessel_max_generations_qlabel;
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
    QLabel *glial_pop3_mean_process_length_qlabel;
    QLabel *glial_pop3_std_process_length_qlabel;
    QLabel *beading_amplitude_qlabel;
    QLabel *beading_std_qlabel;
    QLabel *alpha_qlabel;
    QLabel *beta_qlabel;
    QLabel *glial_pop1_radius_mean_qlabel;
    QLabel *glial_pop1_radius_std_qlabel;
    QLabel *glial_pop2_radius_mean_qlabel;
    QLabel *glial_pop2_radius_std_qlabel;
    QLabel *glial_pop3_radius_mean_qlabel;
    QLabel *glial_pop3_radius_std_qlabel;
    QLabel *glial_pop1_minimum_process_radius_qlabel;
    QLabel *glial_pop2_minimum_process_radius_qlabel;
    QLabel *glial_pop3_minimum_process_radius_qlabel;
    QLabel *glial_pop1_nbr_primary_processes_qlabel;
    QLabel *glial_pop2_nbr_primary_processes_qlabel;
    QLabel *glial_pop3_nbr_primary_processes_qlabel;
    QLabel *glial_pop1_branching_qlabel;
    QLabel *glial_pop2_branching_qlabel;
    QLabel *glial_pop3_branching_qlabel;
    QLabel *k1_qlabel;
    QLabel *k2_qlabel;
    QLabel *k3_qlabel;

    QDoubleSpinBox *nbr_repetitions_SpinBox;
    QCheckBox *visualise_voxel_checkbox;
    QCheckBox *glial_pop1_branching_checkbox;
    QCheckBox *glial_pop2_branching_checkbox;
    QCheckBox *glial_pop3_branching_checkbox;
    QDoubleSpinBox *axons_icvf_SpinBox;
    QDoubleSpinBox *axons_w_myelin_icvf_SpinBox;
    QDoubleSpinBox *glial_pop1_soma_icvf_SpinBox;
    QDoubleSpinBox *glial_pop1_processes_icvf_SpinBox;
    QDoubleSpinBox *glial_pop2_soma_icvf_SpinBox;
    QDoubleSpinBox *glial_pop2_processes_icvf_SpinBox;
    QDoubleSpinBox *glial_pop3_soma_icvf_SpinBox;
    QDoubleSpinBox *glial_pop3_processes_icvf_SpinBox;
    QDoubleSpinBox *blood_vessels_icvf_SpinBox;
    QDoubleSpinBox *blood_vessels_processes_icvf_SpinBox;
    QDoubleSpinBox *blood_vessel_voxel_size_SpinBox;
    QDoubleSpinBox *blood_vessel_mean_radius_SpinBox;
    QDoubleSpinBox *blood_vessel_std_radius_SpinBox;
    QDoubleSpinBox *blood_vessel_capillary_radius_SpinBox;
    QDoubleSpinBox *blood_vessel_max_generations_SpinBox;
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
    QDoubleSpinBox *glial_pop3_mean_process_length_SpinBox;
    QDoubleSpinBox *glial_pop3_std_process_length_SpinBox;
    QDoubleSpinBox *beading_amplitude_SpinBox;
    QDoubleSpinBox *beading_std_SpinBox;
    QDoubleSpinBox *alpha_SpinBox;
    QDoubleSpinBox *beta_SpinBox;
    QDoubleSpinBox *glial_pop1_radius_mean_SpinBox;
    QDoubleSpinBox *glial_pop1_radius_std_SpinBox;
    QDoubleSpinBox *glial_pop2_radius_mean_SpinBox;
    QDoubleSpinBox *glial_pop2_radius_std_SpinBox;
    QDoubleSpinBox *glial_pop3_radius_mean_SpinBox;
    QDoubleSpinBox *glial_pop3_radius_std_SpinBox;
    QDoubleSpinBox *glial_pop1_minimum_process_radius_SpinBox;
    QDoubleSpinBox *glial_pop2_minimum_process_radius_SpinBox;
    QDoubleSpinBox *glial_pop3_minimum_process_radius_SpinBox;
    QDoubleSpinBox *glial_pop1_nbr_primary_processes_SpinBox;
    QDoubleSpinBox *glial_pop2_nbr_primary_processes_SpinBox;
    QDoubleSpinBox *glial_pop3_nbr_primary_processes_SpinBox;
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
    


    // Additional member variables to store values
    int nbr_repetitions;
    bool visualise_voxel;
    Parameters parameters;

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

    std::vector<std::vector<double>> X_glial_pop3;
    std::vector<std::vector<double>> Y_glial_pop3;
    std::vector<std::vector<double>> Z_glial_pop3;
    std::vector<std::vector<double>> R_glial_pop3;
    std::vector<std::vector<int>> Branch_glial_pop3;

    std::vector<std::vector<double>> X_blood_vessels;
    std::vector<std::vector<double>> Y_blood_vessels;
    std::vector<std::vector<double>> Z_blood_vessels;
    std::vector<std::vector<double>> R_blood_vessels;


    // Member for the 3D scatter plot and modifier
    QtDataVisualization::Q3DScatter *graph;
    ScatterDataModifier *modifier;

    
};

#endif // MAINWINDOW_H