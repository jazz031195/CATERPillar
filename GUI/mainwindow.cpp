#include <iostream>
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "slidergroup.h"
#include "../src/axongammadistribution.h"
#include "../src/parameters.h"
#include "ScatterDataModifier.h"
#include <fstream>
#include <QFile>
#include <QTextStream>
#include <QMessageBox>
#include <QDebug>
#include <QMessageBox>
#include <QFileDialog>
#include <chrono>
#include <variant>
#include <QComboBox>
#include <QFontComboBox>
#include "qcustomplot-source/qcustomplot.h"


Window::Window(QWidget *parent)
    : QWidget(parent)
{
    createControls(tr("Parameters"));

    // Create a vertical layout for the entire window
    QVBoxLayout *mainLayout = new QVBoxLayout(this);
    setLayout(mainLayout);

    mainLayout->addWidget(controlsGroup);

    okButton = new QPushButton("OK", this);
    selectDirectoryButton = new QPushButton("Select Directory", this); // Create button for selecting directory

    QHBoxLayout *buttonLayout = new QHBoxLayout;
    buttonLayout->addWidget(okButton);
    buttonLayout->addWidget(selectDirectoryButton); // Add the new button to the layout
    mainLayout->addLayout(buttonLayout);

    connect(okButton, &QPushButton::clicked, this, &Window::onSaveButtonClicked);
    connect(selectDirectoryButton, &QPushButton::clicked, this, &Window::onSelectDirectoryButtonClicked); // Connect button to slot
}
void Window::PlotCells(){
    // Plotting cells
    //QtDataVisualization::Q3DScatter *graph = new QtDataVisualization::Q3DScatter();
    OpenGLWindow *openglWindow = new OpenGLWindow();

    // Ensure that the window title and size are set before displaying
    openglWindow->setTitle("3D Spheres Visualization");
    openglWindow->resize(800, 600);

    // Pass the spheres data to the OpenGL window and trigger an update
    std::vector<std::vector<double>> X = X_axons;
    X.insert(X.end(), X_astrocytes.begin(), X_astrocytes.end());
    X.insert(X.end(), X_oligodendrocytes.begin(), X_oligodendrocytes.end());
    std::vector<std::vector<double>> Y = Y_axons;
    Y.insert(Y.end(), Y_astrocytes.begin(), Y_astrocytes.end());
    Y.insert(Y.end(), Y_oligodendrocytes.begin(), Y_oligodendrocytes.end());
    std::vector<std::vector<double>> Z = Z_axons;
    Z.insert(Z.end(), Z_astrocytes.begin(), Z_astrocytes.end());
    Z.insert(Z.end(), Z_oligodendrocytes.begin(), Z_oligodendrocytes.end());
    std::vector<std::vector<double>> R = R_axons;
    R.insert(R.end(), R_astrocytes.begin(), R_astrocytes.end());
    R.insert(R.end(), R_oligodendrocytes.begin(), R_oligodendrocytes.end());
    openglWindow->setSpheres(X, Y, Z, R);

    // This will ensure the window is updated with new data
    //openglWindow->update();

    QWidget *container = QWidget::createWindowContainer(openglWindow);
    QWidget *widget = new QWidget;
    QHBoxLayout *hLayout = new QHBoxLayout(widget);
    QVBoxLayout *vLayout = new QVBoxLayout();
    hLayout->addWidget(container, 1);
    hLayout->addLayout(vLayout);
    ScatterDataModifier *modifier = new ScatterDataModifier(openglWindow);

    QPushButton *plotRadiiButton = new QPushButton(widget);  // New button for radii distribution
    plotRadiiButton->setText(QStringLiteral("Plot Radii Distribution"));

    QPushButton *plotTortuosityButton = new QPushButton(widget);  // New button for tortuosity distribution
    plotTortuosityButton->setText(QStringLiteral("Plot Tortuosity Distribution"));

    QPushButton *plotShollAnalysisButton = new QPushButton(widget);  // New button for Sholl analysis
    plotShollAnalysisButton->setText(QStringLiteral("Plot Sholl Analysis"));

    QPushButton *resetCameraButton = new QPushButton(widget);  // New button for Sholl analysis
    resetCameraButton->setText(QStringLiteral("Reset Camera"));


    // Add widgets to the vertical layout

    vLayout->addWidget(plotRadiiButton, 0, Qt::AlignTop);  // Add plot radii button to layout
    vLayout->addWidget(plotTortuosityButton, 0, Qt::AlignTop);  // Add plot tortuosity button to layout
    vLayout->addWidget(plotShollAnalysisButton, 0, Qt::AlignTop);  // Add plot Sholl analysis button to layout
    vLayout->addWidget(resetCameraButton, 0, Qt::AlignTop);  // Add reset camera button to layout

    // Connect the "Plot Radii Distribution" button to the plotRadiusDistribution function
    QObject::connect(plotRadiiButton, &QPushButton::clicked, this, &Window::plotRadiusDistribution);
    QObject::connect(plotTortuosityButton, &QPushButton::clicked, this, &Window::plotTortuosityDistribution);
    QObject::connect(plotShollAnalysisButton, &QPushButton::clicked, this, &Window::ShollAnalysis);
    QObject::connect(resetCameraButton, &QPushButton::clicked, openglWindow, &OpenGLWindow::resetCamera);

    widget->show();
}

void Window::resetCamera(){
    openglWindow->resetCamera();
}

void Window::createControls(const QString &title)
{
    controlsGroup = new QGroupBox(title);

     // --- Create the GroupBoxes ---
    QGroupBox *generalGroup = new QGroupBox("General Parameters");
    QGroupBox *axonsGroup = new QGroupBox("Axon Parameters");
    QGroupBox *glialGroup = new QGroupBox("Glial Cell Parameters");

    // add ticked box
    visualise_voxel_qlabel = new QLabel(tr("Visualise Voxel:"));
    axons_icvf_qlabel = new QLabel(tr("Axons ICVF (%):"));
    axons_w_myelin_icvf_qlabel = new QLabel(tr("Axons with myelin ICVF (%):"));
    astrocyte_soma_icvf_qlabel = new QLabel(tr("Astrocyte somas ICVF (%):"));
    astrocyte_processes_icvf_qlabel = new QLabel(tr("Astrocyte processes ICVF (%):"));
    oligodendrocyte_soma_icvf_qlabel = new QLabel(tr("Oligodendrocyte somas ICVF (%):"));
    oligodendrocyte_processes_icvf_qlabel = new QLabel(tr("Oligodendrocyte processes ICVF (%):"));
    voxel_size_qlabel = new QLabel(tr("Voxel Edge Length (μm):"));
    minimum_radius_qlabel = new QLabel(tr("Minimum Sphere Radius (μm):"));
    nbr_threads_qlabel = new QLabel(tr("Number of Threads:"));
    overlapping_factor_qlabel = new QLabel(tr("Overlapping Factor (R/f):"));
    c2_qlabel = new QLabel(tr("c2 (fODF):"));
    nbr_axons_populations_qlabel = new QLabel(tr("Number of axons populations:"));
    epsilon_qlabel = new QLabel(tr("ε (tortuousity):"));
    beading_amplitude_qlabel = new QLabel(tr("Beading Amplitude :"));
    mean_process_length_qlabel = new QLabel(tr("Mean Process Length (μm):"));
    std_process_length_qlabel = new QLabel(tr("Standard Deviation Process Length (μm):"));
    alpha_qlabel = new QLabel(tr("Gamma Distribution for radii α:"));
    beta_qlabel = new QLabel(tr("Gamma Distribution for radii β:"));
    astrocyte_radius_mean_qlabel = new QLabel(tr("Astrocyte Radius Mean:"));
    astrocyte_radius_std_qlabel = new QLabel(tr("Astrocyte Radius Standard Deviation:"));
    oligodendrocyte_radius_mean_qlabel = new QLabel(tr("Oligodendrocyte Radius Mean:"));
    oligodendrocyte_radius_std_qlabel = new QLabel(tr("Oligodendrocyte Radius Standard Deviation:"));

    visualise_voxel_checkbox = new QCheckBox;
    visualise_voxel_checkbox->setChecked(true);

    beading_amplitude_SpinBox = new QDoubleSpinBox;
    beading_amplitude_SpinBox->setRange(0, 1);
    beading_amplitude_SpinBox->setSingleStep(0.1);
    beading_amplitude_SpinBox->setValue(0.3);

    alpha_SpinBox = new QDoubleSpinBox;
    alpha_SpinBox->setRange(0, 10);
    alpha_SpinBox->setSingleStep(0.1);
    alpha_SpinBox->setValue(4);

    beta_SpinBox = new QDoubleSpinBox;
    beta_SpinBox->setRange(0, 10);
    beta_SpinBox->setSingleStep(0.1);
    beta_SpinBox->setValue(0.25);

    epsilon_SpinBox = new QDoubleSpinBox;
    epsilon_SpinBox->setRange(0, 2);
    epsilon_SpinBox->setSingleStep(0.1);
    epsilon_SpinBox->setValue(0.4);

    mean_process_length_SpinBox = new QDoubleSpinBox;
    mean_process_length_SpinBox->setRange(0, 100);
    mean_process_length_SpinBox->setSingleStep(1);
    mean_process_length_SpinBox->setValue(20);

    std_process_length_SpinBox = new QDoubleSpinBox;
    std_process_length_SpinBox->setRange(0, 100);
    std_process_length_SpinBox->setSingleStep(1);
    std_process_length_SpinBox->setValue(10);

    axons_icvf_SpinBox = new QDoubleSpinBox;
    axons_icvf_SpinBox->setRange(0, 100);
    axons_icvf_SpinBox->setSingleStep(1);

    axons_w_myelin_icvf_SpinBox = new QDoubleSpinBox;
    axons_w_myelin_icvf_SpinBox->setRange(0, 100);
    axons_w_myelin_icvf_SpinBox->setSingleStep(1);

    astrocyte_soma_icvf_SpinBox = new QDoubleSpinBox;
    astrocyte_soma_icvf_SpinBox->setRange(0, 100);
    astrocyte_soma_icvf_SpinBox->setSingleStep(1);

    astrocyte_processes_icvf_SpinBox = new QDoubleSpinBox;
    astrocyte_processes_icvf_SpinBox->setRange(0, 100);
    astrocyte_processes_icvf_SpinBox->setSingleStep(1);

    oligodendrocyte_soma_icvf_SpinBox = new QDoubleSpinBox;
    oligodendrocyte_soma_icvf_SpinBox->setRange(0, 100);
    oligodendrocyte_soma_icvf_SpinBox->setSingleStep(1);

    oligodendrocyte_processes_icvf_SpinBox = new QDoubleSpinBox;
    oligodendrocyte_processes_icvf_SpinBox->setRange(0, 100);
    oligodendrocyte_processes_icvf_SpinBox->setSingleStep(1);

    nbr_threads_SpinBox = new QDoubleSpinBox;
    nbr_threads_SpinBox->setRange(1, 1000);
    nbr_threads_SpinBox->setSingleStep(1);

    voxel_size_SpinBox = new QDoubleSpinBox;
    voxel_size_SpinBox->setRange(10, 1000);
    voxel_size_SpinBox->setSingleStep(1);
    voxel_size_SpinBox->setValue(30);

    minimum_radius_SpinBox = new QDoubleSpinBox;
    minimum_radius_SpinBox->setRange(0.05, 10);
    minimum_radius_SpinBox->setSingleStep(0.05);
    minimum_radius_SpinBox->setValue(0.15);

    overlapping_factor_SpinBox = new QDoubleSpinBox;
    overlapping_factor_SpinBox->setRange(1, 64);
    overlapping_factor_SpinBox->setSingleStep(1);
    overlapping_factor_SpinBox->setValue(4);

    astrocyte_radius_mean_SpinBox = new QDoubleSpinBox;
    astrocyte_radius_mean_SpinBox->setRange(0, 10);
    astrocyte_radius_mean_SpinBox->setSingleStep(0.1);
    astrocyte_radius_mean_SpinBox->setValue(5);

    astrocyte_radius_std_SpinBox = new QDoubleSpinBox;
    astrocyte_radius_std_SpinBox->setRange(0, 10);
    astrocyte_radius_std_SpinBox->setSingleStep(0.1);
    astrocyte_radius_std_SpinBox->setValue(0.5);

    oligodendrocyte_radius_mean_SpinBox = new QDoubleSpinBox;
    oligodendrocyte_radius_mean_SpinBox->setRange(0, 10);
    oligodendrocyte_radius_mean_SpinBox->setSingleStep(0.1);
    oligodendrocyte_radius_mean_SpinBox->setValue(5);

    oligodendrocyte_radius_std_SpinBox = new QDoubleSpinBox;
    oligodendrocyte_radius_std_SpinBox->setRange(0, 10);
    oligodendrocyte_radius_std_SpinBox->setSingleStep(0.1);
    oligodendrocyte_radius_std_SpinBox->setValue(0.5);




    // --- Configuration ComboBox (Initially Hidden) ---
    configurationComboBox = new QComboBox;
    configurationComboBox->addItem("Sheet Configuration");
    configurationComboBox->addItem("Interwoven Configuration");
    configurationComboBox->setVisible(false); // Initially hidden

    nbr_axons_populations_SpinBox = new QDoubleSpinBox;
    nbr_axons_populations_SpinBox->setRange(1, 3);
    nbr_axons_populations_SpinBox->setSingleStep(1);

    // --- Connect the SpinBox Signal to a Slot Function ---
    connect(nbr_axons_populations_SpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            this, &Window::updateConfigurationSelectionVisibility);



    c2_SpinBox = new QDoubleSpinBox;
    c2_SpinBox->setRange(0, 1);
    c2_SpinBox->setSingleStep(0.05);
    c2_SpinBox->setValue(1);

    QGridLayout *controlsLayout = new QGridLayout;

    std::vector<QLabel*> general_labels = {voxel_size_qlabel, overlapping_factor_qlabel, minimum_radius_qlabel};
    std::vector <QDoubleSpinBox*> general_spinBoxes = { voxel_size_SpinBox, overlapping_factor_SpinBox, minimum_radius_SpinBox};
    
    std::vector<QLabel*> axons_labels = {axons_w_myelin_icvf_qlabel, axons_icvf_qlabel ,nbr_threads_qlabel, epsilon_qlabel, c2_qlabel, nbr_axons_populations_qlabel, beading_amplitude_qlabel, alpha_qlabel, beta_qlabel};
    std::vector <QDoubleSpinBox*> axons_spinBoxes = {axons_w_myelin_icvf_SpinBox, axons_icvf_SpinBox, nbr_threads_SpinBox, epsilon_SpinBox, c2_SpinBox, nbr_axons_populations_SpinBox, beading_amplitude_SpinBox, alpha_SpinBox, beta_SpinBox};
    
    std::vector<QLabel*> glials_labels = {astrocyte_soma_icvf_qlabel, astrocyte_processes_icvf_qlabel, oligodendrocyte_soma_icvf_qlabel, oligodendrocyte_processes_icvf_qlabel, astrocyte_radius_mean_qlabel, astrocyte_radius_std_qlabel, oligodendrocyte_radius_mean_qlabel, oligodendrocyte_radius_std_qlabel, mean_process_length_qlabel, std_process_length_qlabel};
    std::vector <QDoubleSpinBox*> glials_spinBoxes = {astrocyte_soma_icvf_SpinBox, astrocyte_processes_icvf_SpinBox, oligodendrocyte_soma_icvf_SpinBox, oligodendrocyte_processes_icvf_SpinBox, astrocyte_radius_mean_SpinBox, astrocyte_radius_std_SpinBox, oligodendrocyte_radius_mean_SpinBox, oligodendrocyte_radius_std_SpinBox,  mean_process_length_SpinBox, std_process_length_SpinBox};
    
    
    QGridLayout *generalLayout = new QGridLayout;

    for (int i = 0; i < general_labels.size(); i++){
        generalLayout->addWidget(general_labels[i], i, 0);
        generalLayout->addWidget(general_spinBoxes[i], i, 1);
    }
    generalLayout->addWidget(visualise_voxel_qlabel, general_labels.size(), 0);
    generalLayout->addWidget(visualise_voxel_checkbox, general_labels.size(), 1);

    generalGroup->setLayout(generalLayout);

    QGridLayout *axonsLayout = new QGridLayout;
    for (int i = 0; i < axons_labels.size(); i++){
        axonsLayout->addWidget(axons_labels[i], i, 0);
        axonsLayout->addWidget(axons_spinBoxes[i], i, 1);
    }
    axonsGroup->setLayout(axonsLayout);

    QGridLayout *glialLayout = new QGridLayout;
    for (int i = 0; i < glials_labels.size(); i++){
        glialLayout->addWidget(glials_labels[i], i, 0);
        glialLayout->addWidget(glials_spinBoxes[i], i, 1);
    }
    glialGroup->setLayout(glialLayout);

    // --- Main Layout ---
    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addWidget(generalGroup);
    mainLayout->addWidget(axonsGroup);
    mainLayout->addWidget(glialGroup);

    controlsGroup->setLayout(mainLayout);

}
void Window::updateConfigurationSelectionVisibility(double value)
{
    configurationComboBox->setVisible(static_cast<int>(value) == 2);
}

void Window::resizeEvent(QResizeEvent *)
{
    if (width() == 0 || height() == 0)
        return;

}

void Window::onSaveButtonClicked()
{
    // Retrieve values from spin boxes and checkboxes
    axons_icvf = axons_icvf_SpinBox->value();
    axons_w_myelin_icvf = axons_w_myelin_icvf_SpinBox->value();
    astrocyte_soma_icvf = astrocyte_soma_icvf_SpinBox->value();
    astrocyte_processes_icvf = astrocyte_processes_icvf_SpinBox->value();
    oligodendrocyte_soma_icvf = oligodendrocyte_soma_icvf_SpinBox->value();
    oligodendrocyte_processes_icvf = oligodendrocyte_processes_icvf_SpinBox->value();
    nbr_threads = nbr_threads_SpinBox->value();
    overlapping_factor = overlapping_factor_SpinBox->value();
    voxel_size = voxel_size_SpinBox->value();
    minimum_radius = minimum_radius_SpinBox->value();
    c2 = c2_SpinBox->value();
    nbr_axons_populations = nbr_axons_populations_SpinBox->value();
    beading_amplitude = beading_amplitude_SpinBox->value();
    mean_process_length = mean_process_length_SpinBox->value();
    std_process_length = std_process_length_SpinBox->value();
    epsilon = epsilon_SpinBox->value();
    minimum_radius = minimum_radius_SpinBox->value();
    alpha = alpha_SpinBox->value();
    beta = beta_SpinBox->value();
    visualise_voxel = visualise_voxel_checkbox->isChecked();
    astrocyte_radius_mean = astrocyte_radius_mean_SpinBox->value();
    astrocyte_radius_std = astrocyte_radius_std_SpinBox->value();
    oligodendrocyte_radius_mean = oligodendrocyte_radius_mean_SpinBox->value();
    oligodendrocyte_radius_std = oligodendrocyte_radius_std_SpinBox->value();
    
    // Close the parameter input dialog
    this->close();

    StartSimulation();


}

void Window::StartSimulation(){

    Parameters parameters;

    int repetitions = parameters.repetitions;
    std::vector<double> vox_sizes = {voxel_size/1.0};
    double astrocyte_icvf_soma = astrocyte_soma_icvf/100.0;
    double astrocyte_icvf_branches = astrocyte_processes_icvf/100.0;
    double oligodendrocyte_icvf_soma = oligodendrocyte_soma_icvf/100.0;
    double oligodendrocyte_icvf_branches = oligodendrocyte_processes_icvf/100.0;
    double axons_wo_myelin_icvf = axons_icvf/100.0;
    double axons_with_myelin_icvf = axons_w_myelin_icvf/100.0;
    int spheres_overlap_factor = overlapping_factor;
    std::string directory = selectedDirectory.toStdString();
    double cosphisquared = c2;
    double std_dev = epsilon;
    double mean_glial_process_length = mean_process_length;
    double std_glial_process_length = std_process_length;
    double min_radius = minimum_radius;

    int regrow_thr = parameters.regrow_thr;
    int ondulation_factor = parameters.ondulation_factor;
    bool can_shrink = parameters.can_shrink;
    int crossing_fibers_type = 0;

    if (nbr_axons_populations == 2) {
        QString PopConfiguration = configurationComboBox->currentText();
        if (PopConfiguration == "Sheet Configuration") {
            crossing_fibers_type = 0;
        } else if (PopConfiguration == "Interwoven Configuration") {
            crossing_fibers_type = 1;
        }
    }

    for (int rep = 0; rep < repetitions; rep++){
        for (unsigned long i = 0; i < vox_sizes.size(); i++){

            double vox_size= vox_sizes[i];

            // min and max limits of voxel
            Eigen::Vector3d min_l = {0, 0, 0};
            Eigen::Vector3d max_l = {vox_size, vox_size, vox_size}; // um

            auto startTime = std::chrono::high_resolution_clock::now();
            cout << "STARTING SIMULATION " << endl;
            // create distribution of axons
            AxonGammaDistribution AxonDistribution = AxonGammaDistribution(axons_wo_myelin_icvf, axons_with_myelin_icvf, astrocyte_icvf_soma, astrocyte_icvf_branches, oligodendrocyte_icvf_soma, oligodendrocyte_icvf_branches, alpha, beta, min_l, max_l, min_radius, regrow_thr, beading_amplitude, std_dev, ondulation_factor, spheres_overlap_factor, can_shrink, cosphisquared, nbr_threads, nbr_axons_populations, 
            crossing_fibers_type, mean_glial_process_length, std_glial_process_length, astrocyte_radius_mean, astrocyte_radius_std, oligodendrocyte_radius_mean, oligodendrocyte_radius_std);

            AxonDistribution.createSubstrate();
            // saving spheres
            for (unsigned i=0; i< AxonDistribution.axons.size(); ++i){
                std::vector<double> x_;
                std::vector<double> y_;
                std::vector<double> z_;
                std::vector<double> r_;
                for (unsigned j=0; j< AxonDistribution.axons[i].outer_spheres.size(); ++j){
                    x_.push_back(AxonDistribution.axons[i].outer_spheres[j].center[0]);
                    y_.push_back(AxonDistribution.axons[i].outer_spheres[j].center[1]);
                    z_.push_back(AxonDistribution.axons[i].outer_spheres[j].center[2]);
                    r_.push_back(AxonDistribution.axons[i].outer_spheres[j].radius);
                }
                X_axons.push_back(x_);
                Y_axons.push_back(y_);
                Z_axons.push_back(z_);
                R_axons.push_back(r_);
                x_.clear();
                y_.clear();
                z_.clear();
                r_.clear();
            }

            for (unsigned i=0; i< AxonDistribution.astrocytes.size(); ++i){
                std::vector<double> x_;
                std::vector<double> y_;
                std::vector<double> z_;
                std::vector<double> r_;
                std::vector<int> b_;

                x_.push_back(AxonDistribution.astrocytes[i].soma.center[0]);
                y_.push_back(AxonDistribution.astrocytes[i].soma.center[1]);
                z_.push_back(AxonDistribution.astrocytes[i].soma.center[2]);
                r_.push_back(AxonDistribution.astrocytes[i].soma.radius);
                b_.push_back(0);

                for (unsigned j=0; j< AxonDistribution.astrocytes[i].ramification_spheres.size(); ++j){

                    for (unsigned k=0; k< AxonDistribution.astrocytes[i].ramification_spheres[j].size(); ++k){
                        x_.push_back(AxonDistribution.astrocytes[i].ramification_spheres[j][k].center[0]);
                        y_.push_back(AxonDistribution.astrocytes[i].ramification_spheres[j][k].center[1]);
                        z_.push_back(AxonDistribution.astrocytes[i].ramification_spheres[j][k].center[2]);
                        r_.push_back(AxonDistribution.astrocytes[i].ramification_spheres[j][k].radius);
                        b_.push_back(j);
                    }
                }

                X_astrocytes.push_back(x_);
                Y_astrocytes.push_back(y_);
                Z_astrocytes.push_back(z_);
                R_astrocytes.push_back(r_);
                Branch_astrocytes.push_back(b_);
                x_.clear();
                y_.clear();
                z_.clear();
                r_.clear();
                b_.clear();
            }

            for (unsigned i=0; i< AxonDistribution.oligodendrocytes.size(); ++i){
                std::vector<double> x_;
                std::vector<double> y_;
                std::vector<double> z_;
                std::vector<double> r_;
                std::vector<int> b_;

                x_.push_back(AxonDistribution.oligodendrocytes[i].soma.center[0]);
                y_.push_back(AxonDistribution.oligodendrocytes[i].soma.center[1]);
                z_.push_back(AxonDistribution.oligodendrocytes[i].soma.center[2]);
                r_.push_back(AxonDistribution.oligodendrocytes[i].soma.radius);
                b_.push_back(0);

                for (unsigned j=0; j< AxonDistribution.oligodendrocytes[i].ramification_spheres.size(); ++j){

                    for (unsigned k=0; k< AxonDistribution.oligodendrocytes[i].ramification_spheres[j].size(); ++k){
                        x_.push_back(AxonDistribution.oligodendrocytes[i].ramification_spheres[j][k].center[0]);
                        y_.push_back(AxonDistribution.oligodendrocytes[i].ramification_spheres[j][k].center[1]);
                        z_.push_back(AxonDistribution.oligodendrocytes[i].ramification_spheres[j][k].center[2]);
                        r_.push_back(AxonDistribution.oligodendrocytes[i].ramification_spheres[j][k].radius);
                        b_.push_back(j);
                    }
                }

                X_oligodendrocytes.push_back(x_);
                Y_oligodendrocytes.push_back(y_);
                Z_oligodendrocytes.push_back(z_);
                R_oligodendrocytes.push_back(r_);
                Branch_oligodendrocytes.push_back(b_);
                x_.clear();
                y_.clear();
                z_.clear();
                r_.clear();
                b_.clear();
            }

            std::string simulation_file_name;
            std::string swc_file_name;

            std::string filename;
            std::ifstream file(filename);

            simulation_file_name = (directory + "/growth_info.txt");
            swc_file_name = (directory + "/Voxel.swc");
            std::ofstream swc_file(swc_file_name);
            std::ofstream simulation_file(simulation_file_name);

            // Check if files opened successfully
            if (!swc_file)
            {
                std::cerr << "Error opening output file : "<< swc_file_name << std::endl;

            }
            cout << "Creating file: " << swc_file_name << endl;
            // write to file
            AxonDistribution.create_SWC_file(swc_file);
            swc_file.close();

            // Check if files opened successfully

            if (!simulation_file)
            {
                std::cerr << "Error opening output file : " << simulation_file_name <<std::endl;
            }

            auto endTime = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime);
            AxonDistribution.simulation_file(simulation_file, duration);
            simulation_file.close();

            cout << "End of simulation!" << endl;

            
        }
    }

    if (visualise_voxel) {
        // After simulation completes, call PlotCells to display the data
        PlotCells();
    }
    else{
        // Display a message box to inform the user that the simulation is complete
        QMessageBox::information(this, "Simulation Complete", "Simulation complete! Please check the output directory for the results.");
    }

}

void Window::onSelectDirectoryButtonClicked() {
    QString dirPath = QFileDialog::getExistingDirectory(this, tr("Select Directory"), "", QFileDialog::ShowDirsOnly);
    if (!dirPath.isEmpty()) {
        selectedDirectory = dirPath;
        qDebug() << "Selected directory:" << selectedDirectory;
    }
}


void Window::createStatisticsMenu()
{
    statisticsButton = new QPushButton("Statistics", this);
    plotRadiusDistributionButton = new QPushButton("Plot Radius Distribution", this);
    plotTortuosityDistributionButton = new QPushButton("Plot Tortuosity Distribution", this);  // New Button
    plotShollAnalysisButton = new QPushButton("Plot Sholl Analysis", this);
    resetCameraButton = new QPushButton("Reset Camera", this);

    QVBoxLayout *statisticsLayout = new QVBoxLayout;
    statisticsLayout->addWidget(statisticsButton);
    statisticsLayout->addWidget(plotRadiusDistributionButton);
    statisticsLayout->addWidget(plotTortuosityDistributionButton);  // Add the new button
    statisticsLayout->addWidget(plotShollAnalysisButton); 
    statisticsLayout->addWidget(resetCameraButton);

    controlsGroup->setLayout(statisticsLayout);

    // Connect the button clicks to the appropriate functions
    connect(plotRadiusDistributionButton, &QPushButton::clicked, this, &Window::plotRadiusDistribution);
    connect(plotTortuosityDistributionButton, &QPushButton::clicked, this, &Window::plotTortuosityDistribution);  // Connect the new button
    connect(plotShollAnalysisButton, &QPushButton::clicked, this, &Window::ShollAnalysis);
    connect(resetCameraButton, &QPushButton::clicked, openglWindow, &OpenGLWindow::resetCamera);
}


void Window::plotRadiusDistribution()
{

    if (X_axons.size() == 0){
        QMessageBox::warning(this, "Error", "No Axons to plot!");
        return;
    }
    // Map each axon ID to a vector of radii
    std::map<int, std::vector<double>> axonRadiiMap;

    // Iterate over all the spheres and collect radii for each axon
    for (size_t i = 0; i < X_axons.size(); ++i) {
        for (size_t j = 0; j < X_axons[i].size(); ++j) {
            int axonID = i; // Change this to the correct axon ID reference
            axonRadiiMap[axonID].push_back(R_axons[i][j]);  // Add radius to the corresponding axon
        }
    }

    // Calculate mean radius for each axon
    std::vector<double> meanRadii;
    for (const auto& axon : axonRadiiMap) {
        double sum = std::accumulate(axon.second.begin(), axon.second.end(), 0.0);
        double mean = sum / axon.second.size();
        meanRadii.push_back(mean);
    }

    // Sort the meanRadii for binning
    std::sort(meanRadii.begin(), meanRadii.end());

    // Calculate the histogram (binning)
    int binCount = 10;  // Number of bins, adjust this as needed
    double minRadius = *std::min_element(meanRadii.begin(), meanRadii.end());
    double maxRadius = *std::max_element(meanRadii.begin(), meanRadii.end());

    if (maxRadius == minRadius) {
        maxRadius += 1.0;  // Avoid division by zero in case all radii are the same
    }

    double binWidth = (maxRadius - minRadius) / binCount;

    QVector<double> bins(binCount, 0);  // Initialize bin counts to zero
    QVector<double> tickPositions(binCount);  // Positions on the x-axis

    // Generate tick positions (center of each bin)
    for (int i = 0; i < binCount; ++i) {
        tickPositions[i] = minRadius + binWidth * (i + 0.5);  // Center of each bin
    }

    // Assign mean radii to bins
    for (double radius : meanRadii) {
        int binIndex = static_cast<int>((radius - minRadius) / binWidth);
        // Clamp the bin index to make sure it is within bounds
        binIndex = std::min(std::max(binIndex, 0), binCount - 1);
        bins[binIndex]++;
    }

    // Create QCustomPlot and QCPBars for the histogram
    QCustomPlot *customPlot = new QCustomPlot;

    // Prepare the bars for the histogram
    QCPBars *histogram = new QCPBars(customPlot->xAxis, customPlot->yAxis);

    // Set data for the histogram
    histogram->setData(tickPositions, bins);

    // Set the width of each bar to match the bin width
    histogram->setWidth(binWidth);  // Use the binWidth for the bar width

    // Configure axis labels and ranges
    customPlot->xAxis->setLabel("Mean Radius");
    customPlot->yAxis->setLabel("Count");

    customPlot->xAxis->setRange(minRadius, maxRadius);
    customPlot->yAxis->setRange(0, *std::max_element(bins.begin(), bins.end()));

    // Display the plot in a dialog window
    QDialog *dialog = new QDialog(this);
    QVBoxLayout *layout = new QVBoxLayout;
    layout->addWidget(customPlot);
    dialog->setLayout(layout);
    dialog->setWindowTitle("Mean Radius Distribution");
    dialog->exec();
}

void Window::plotTortuosityDistribution()
{

    if (X_axons.size() == 0){
        QMessageBox::warning(this, "Error", "No Axons to plot!");
        return;
    }
    // Map each axon ID to a vector of sphere positions
    std::map<int, std::vector<Eigen::Vector3d>> axonPositionMap;

    // Iterate over all the spheres and collect positions for each axon
    if (X_axons.size() != Y_axons.size() || X_axons.size() != Z_axons.size()) {
        qDebug() << "Error: X, Y, Z vectors have different sizes!";
        return;
    }

    for (size_t i = 0; i < X_axons.size(); ++i) {
        if (X_axons[i].size() != Y_axons[i].size() || X_axons[i].size() != Z_axons[i].size()) {
            qDebug() << "Error: Mismatch in sphere sizes in axon " << i;
            continue;  // Skip this axon if sizes don't match
        }

        for (size_t j = 0; j < X_axons[i].size(); ++j) {
            Eigen::Vector3d position(X_axons[i][j], Y_axons[i][j], Z_axons[i][j]);
            int axonID = i;
            axonPositionMap[axonID].push_back(position);  // Add position to the corresponding axon
        }
    }

    // Calculate tortuosity for each axon
    std::vector<double> tortuosities;
    for (const auto& axon : axonPositionMap) {
        const std::vector<Eigen::Vector3d>& positions = axon.second;

        if (positions.size() < 2) {
            continue;  // Skip if there are less than 2 spheres
        }

        // Calculate the total length of the axon (sum of distances between consecutive spheres)
        double totalLength = 0.0;
        for (size_t i = 1; i < positions.size(); ++i) {
            totalLength += (positions[i] - positions[i - 1]).norm();
        }

        // Calculate the direct distance between the first and last sphere
        double directDistance = (positions.back() - positions.front()).norm();

        // Avoid division by zero (in case direct distance is 0)
        if (directDistance > 0) {
            // Calculate tortuosity: total length / direct distance
            double tortuosity = totalLength / directDistance;
            tortuosities.push_back(tortuosity);
        } else {
            qDebug() << "Warning: Direct distance is zero for axon " << axon.first;
        }
    }

    // Check if we have tortuosity values
    if (tortuosities.empty()) {
        qDebug() << "No valid tortuosity values calculated!";
        return;
    }

    // Sort the tortuosities for binning
    std::sort(tortuosities.begin(), tortuosities.end());

    // Calculate the histogram (binning)
    int binCount = 10;  // Number of bins, adjust this as needed
    double minTortuosity = *std::min_element(tortuosities.begin(), tortuosities.end());
    double maxTortuosity = *std::max_element(tortuosities.begin(), tortuosities.end());

    if (minTortuosity == maxTortuosity) {
        qDebug() << "Tortuosity range is zero. Cannot create a meaningful histogram.";
        return;
    }

    double binWidth = (maxTortuosity - minTortuosity) / binCount;
    QVector<double> bins(binCount, 0);  // Initialize bin counts to zero
    QVector<double> tickPositions(binCount);  // Positions on the x-axis

    // Assign tortuosities to bins
    for (double tortuosity : tortuosities) {
        int binIndex = std::min(static_cast<int>((tortuosity - minTortuosity) / binWidth), binCount - 1);
        bins[binIndex]++;
    }

    // Generate tick positions (the center of each bin)
    for (int i = 0; i < binCount; ++i) {
        tickPositions[i] = minTortuosity + binWidth * (i + 0.5);  // Center of each bin
    }

    // Create QCustomPlot and QCPBars for the histogram
    QCustomPlot *customPlot = new QCustomPlot;

    // Prepare the bars for the histogram
    QCPBars *histogram = new QCPBars(customPlot->xAxis, customPlot->yAxis);

    // Set data for the histogram
    histogram->setData(tickPositions, bins);

    // Set the width of each bar to match the bin width
    histogram->setWidth(binWidth);  // Use the binWidth for the bar width

    // Configure axis labels and ranges
    customPlot->xAxis->setLabel("Tortuosity");
    customPlot->yAxis->setLabel("Count");

    customPlot->xAxis->setRange(minTortuosity, maxTortuosity);
    customPlot->yAxis->setRange(0, *std::max_element(bins.begin(), bins.end()));

    // Display the plot in a dialog window
    QDialog *dialog = new QDialog(this);
    QVBoxLayout *layout = new QVBoxLayout;
    layout->addWidget(customPlot);
    dialog->setLayout(layout);
    dialog->setWindowTitle("Tortuosity Distribution");
    dialog->exec();
}


void Window::ShollAnalysis() {

    if (X_astrocytes.size() == 0) {
        QMessageBox::warning(this, "Error", "No Astrocytes to plot!");
        return;
    }

    for (unsigned long i = 0; i < X_astrocytes.size(); ++i) {

        // Soma position of the current astrocyte
        Eigen::Vector3d soma_position = {X_astrocytes[i][0], Y_astrocytes[i][0], Z_astrocytes[i][0]};

        // Radii for Sholl analysis
        std::vector<double> sphere_around_soma_radii = {5, 7, 10, 15, 20 , 25, 30, 40, 50, 60, 80};
        
        // Initialize intersections list with zeros
        std::vector<double> intersections_list(sphere_around_soma_radii.size(), 0);

        std::vector<int> branches_list;

        // Iterate through all spheres (excluding the soma) to compute intersections
        for (unsigned long r = 0; r < sphere_around_soma_radii.size(); ++r) {
            for (unsigned long j = 1; j < X_astrocytes[i].size(); ++j) {
                Eigen::Vector3d position = {X_astrocytes[i][j], Y_astrocytes[i][j], Z_astrocytes[i][j]};
                double distance = (position - soma_position).norm();

                if (distance < sphere_around_soma_radii[r] + R_astrocytes[i][j] && distance > sphere_around_soma_radii[r] - R_astrocytes[i][j]) {
                    // check if branch is already in list
                    if (std::find(branches_list.begin(), branches_list.end(), Branch_astrocytes[i][j]) == branches_list.end()) {
                        intersections_list[r] += 1;
                        branches_list.push_back(Branch_astrocytes[i][j]);
                    }
                }
            }
            branches_list.clear();
        }

        // Create QCustomPlot for Sholl analysis
        QCustomPlot *customPlot = new QCustomPlot;

        // Convert the data to QVector for QCustomPlot
        QVector<double> x = QVector<double>::fromStdVector(sphere_around_soma_radii);
        QVector<double> y = QVector<double>::fromStdVector(intersections_list);

        // Create a graph and set the data
        customPlot->addGraph();
        customPlot->graph(0)->setData(x, y);

        // Set axis labels
        customPlot->xAxis->setLabel("Distance to Soma (μm)");
        customPlot->yAxis->setLabel("Number of Intersections");

        // Set axis ranges
        customPlot->xAxis->setRange(0, *std::max_element(sphere_around_soma_radii.begin(), sphere_around_soma_radii.end()));
        customPlot->yAxis->setRange(0, *std::max_element(intersections_list.begin(), intersections_list.end()));

        // Add a title
        customPlot->plotLayout()->insertRow(0);
        //customPlot->plotLayout()->addElement(0, 0, new QCPPlotTitle(customPlot, "Sholl Analysis"));

        // Display the plot in a dialog window
        QDialog *dialog = new QDialog(this);
        //adjust size
        dialog->resize(800, 600);
        QVBoxLayout *layout = new QVBoxLayout;
        layout->addWidget(customPlot);
        dialog->setLayout(layout);
        dialog->setWindowTitle("Sholl Analysis for Astrocyte");
        dialog->exec();
    }
}

