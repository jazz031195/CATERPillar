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
#include <QJsonDocument>
#include <QJsonObject>


Window::Window(QWidget *parent)
    : QWidget(parent)
{

        // Inside your Window constructor or init function
    QVBoxLayout *startupLayout = new QVBoxLayout;

    // 1. Welcome Image
    QLabel *imageLabel = new QLabel(this);
    QPixmap image("../logo_catrepillar.png"); // replace with the actual image path in your resource
    imageLabel->setPixmap(image.scaled(800, 800, Qt::KeepAspectRatio));
    imageLabel->setAlignment(Qt::AlignCenter);
    startupLayout->addWidget(imageLabel);

    // 2. Welcome Buttons
    growButton = new QPushButton("Grow Substrate", this);
    visualiseButton = new QPushButton("Visualise Substrate", this);

    QHBoxLayout *buttonLayout = new QHBoxLayout;
    buttonLayout->addWidget(growButton);
    buttonLayout->addWidget(visualiseButton);

    startupLayout->addLayout(buttonLayout);

    // 3. Central widget stack (to switch between views)
    QWidget *startupWidget = new QWidget;
    startupWidget->setLayout(startupLayout);


    QWidget *mainWidget = new QWidget;
    QVBoxLayout *mainLayout = new QVBoxLayout(mainWidget);

    // Your parameter group
    controlsGroup = createControls(tr("Parameters"));

    mainLayout->addWidget(controlsGroup);

    // OK and Select Directory buttons
    okButton = new QPushButton("OK", this);
    selectDirectoryButton = new QPushButton("Select Directory", this);
    loadConfigButton = new QPushButton("Load Config", this);
    QHBoxLayout *mainButtonLayout = new QHBoxLayout;
    mainButtonLayout->addWidget(loadConfigButton);
    mainButtonLayout->addWidget(selectDirectoryButton);
    mainButtonLayout->addWidget(okButton);
    mainLayout->addLayout(mainButtonLayout);

    // 4. Stacked layout to switch views
    QStackedLayout *stack = new QStackedLayout(this);
    stack->addWidget(startupWidget); // index 0
    stack->addWidget(mainWidget);    // index 1
    this->setLayout(stack);

    // 5. Connections
    connect(growButton, &QPushButton::clicked, [stack]() {
        stack->setCurrentIndex(1); // Show parameters layout
    });


    connect(visualiseButton, &QPushButton::clicked, this, &Window::SelectSWCFileButton);

    connect(okButton, &QPushButton::clicked, this, &Window::onSaveButtonClicked);
    connect(selectDirectoryButton, &QPushButton::clicked, this, &Window::onSelectDirectoryButtonClicked);
    connect(loadConfigButton, &QPushButton::clicked,
            this, &Window::loadConfigFromFile);

}

void Window::SelectSWCFileButton() {
    QString filePath = QFileDialog::getOpenFileName(
        this,
        tr("Select SWC or CSV File"),
        "",
        tr("Data Files (*.csv *.swc);;CSV Files (*.csv);;SWC Files (*.swc)")
    );
    if (!filePath.isEmpty()) {
        SWCFile = filePath;  // Assuming you have a QString member named SWCFile
        qDebug() << "Selected CSV or SWC file:" << SWCFile;
        ReadAxonsFromFile(SWCFile);       // Load axons
        ReadGlialCellsFromFile(SWCFile);  // Load glial cells
        ReadBloodVesselsFromFile(SWCFile); // Load blood vessels
        PlotCells(true, true, true, true);              // Visualize
    }
}
void Window::PlotCells(const bool& axons_plot,
                       const bool& glial_pop1_plot,
                       const bool& glial_pop2_plot,
                       const bool& blood_vessels_plot)
{
    OpenGLWindow *openglWindow = new OpenGLWindow();
    openglWindow->setTitle("3D Spheres Visualization");
    openglWindow->resize(800, 600);

    std::vector<std::vector<double>> X, Y, Z, R;
    std::vector<int> groupIds;  // one entry per bundle (i in setSpheres)

    if (axons_plot){
        X.insert(X.end(), X_axons.begin(), X_axons.end());
        Y.insert(Y.end(), Y_axons.begin(), Y_axons.end());
        Z.insert(Z.end(), Z_axons.begin(), Z_axons.end());
        R.insert(R.end(), R_axons.begin(), R_axons.end());
        groupIds.insert(groupIds.end(), X_axons.size(), static_cast<int>(OpenGLWindow::SphereGroup::Axon));
    }
    if (glial_pop1_plot){
        X.insert(X.end(), X_glial_pop1.begin(), X_glial_pop1.end());
        Y.insert(Y.end(), Y_glial_pop1.begin(), Y_glial_pop1.end());
        Z.insert(Z.end(), Z_glial_pop1.begin(), Z_glial_pop1.end());
        R.insert(R.end(), R_glial_pop1.begin(), R_glial_pop1.end());
        groupIds.insert(groupIds.end(), X_glial_pop1.size(), static_cast<int>(OpenGLWindow::SphereGroup::Glial1));
    }
    if (glial_pop2_plot){
        X.insert(X.end(), X_glial_pop2.begin(), X_glial_pop2.end());
        Y.insert(Y.end(), Y_glial_pop2.begin(), Y_glial_pop2.end());
        Z.insert(Z.end(), Z_glial_pop2.begin(), Z_glial_pop2.end());
        R.insert(R.end(), R_glial_pop2.begin(), R_glial_pop2.end());
        groupIds.insert(groupIds.end(), X_glial_pop2.size(), static_cast<int>(OpenGLWindow::SphereGroup::Glial2));
    }
    if (blood_vessels_plot){
        X.insert(X.end(), X_blood_vessels.begin(), X_blood_vessels.end());
        Y.insert(Y.end(), Y_blood_vessels.begin(), Y_blood_vessels.end());
        Z.insert(Z.end(), Z_blood_vessels.begin(), Z_blood_vessels.end());
        R.insert(R.end(), R_blood_vessels.begin(), R_blood_vessels.end());
        groupIds.insert(groupIds.end(), X_blood_vessels.size(), static_cast<int>(OpenGLWindow::SphereGroup::Blood));
    }

    openglWindow->setSpheres(X, Y, Z, R, groupIds);

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

    QPushButton *hideAxonsButton = new QPushButton(widget);  // New button for Sholl analysis
    hideAxonsButton->setText(QStringLiteral("Hide Axons"));

    QPushButton *hideGlialCellsButton = new QPushButton(widget);  // New button for Sholl analysis
    hideGlialCellsButton->setText(QStringLiteral("Hide Glial Cells"));

    QPushButton *showall = new QPushButton(widget);  // New button for Sholl analysis
    showall->setText(QStringLiteral("Show All Cells"));


    // Add widgets to the vertical layout

    vLayout->addWidget(plotRadiiButton, 0, Qt::AlignTop);  // Add plot radii button to layout
    vLayout->addWidget(plotTortuosityButton, 0, Qt::AlignTop);  // Add plot tortuosity button to layout
    vLayout->addWidget(plotShollAnalysisButton, 0, Qt::AlignTop);  // Add plot Sholl analysis button to layout
    vLayout->addWidget(resetCameraButton, 0, Qt::AlignTop);  // Add reset camera button to layout
    vLayout->addWidget(hideAxonsButton, 0, Qt::AlignTop);  // Add hide axons button to layout
    vLayout->addWidget(hideGlialCellsButton, 0, Qt::AlignTop);  // Add hide glial cells button to layout
    vLayout->addWidget(showall, 0, Qt::AlignTop);  // Add show all button to layout

    // Connect the "Plot Radii Distribution" button to the plotRadiusDistribution function
    QObject::connect(plotRadiiButton, &QPushButton::clicked, this, &Window::plotRadiusDistribution);
    QObject::connect(plotTortuosityButton, &QPushButton::clicked, this, &Window::plotTortuosityDistribution);
    QObject::connect(plotShollAnalysisButton, &QPushButton::clicked, this, &Window::ShollAnalysis);
    QObject::connect(resetCameraButton, &QPushButton::clicked, openglWindow, &OpenGLWindow::resetCamera);
    QObject::connect(hideAxonsButton, &QPushButton::clicked, this, &Window::HideAxons);
    QObject::connect(hideGlialCellsButton, &QPushButton::clicked, this, &Window::HideGlialCells);
    QObject::connect(showall, &QPushButton::clicked, this, &Window::ShowAllCells);


    widget->show();
}

void Window::resetCamera(){
    openglWindow->resetCamera();
}


void Window::ShowAllCells(){

    std::vector<int> groupIds;  // one entry per bundle (i in setSpheres)

    std::vector<std::vector<double>> X = X_glial_pop1;
    std::vector<std::vector<double>> Y = Y_glial_pop1;
    std::vector<std::vector<double>> Z = Z_glial_pop1;
    std::vector<std::vector<double>> R = R_glial_pop1;
    groupIds.insert(groupIds.end(), X_glial_pop1.size(), static_cast<int>(OpenGLWindow::SphereGroup::Glial1));

    X.insert(X.end(), X_glial_pop2.begin(), X_glial_pop2.end());
    Y.insert(Y.end(), Y_glial_pop2.begin(), Y_glial_pop2.end());
    Z.insert(Z.end(), Z_glial_pop2.begin(), Z_glial_pop2.end());
    R.insert(R.end(), R_glial_pop2.begin(), R_glial_pop2.end());
    groupIds.insert(groupIds.end(), X_glial_pop2.size(), static_cast<int>(OpenGLWindow::SphereGroup::Glial2));

    X.insert(X.end(), X_axons.begin(), X_axons.end());
    Y.insert(Y.end(), Y_axons.begin(), Y_axons.end());
    Z.insert(Z.end(), Z_axons.begin(), Z_axons.end());
    R.insert(R.end(), R_axons.begin(), R_axons.end());
    groupIds.insert(groupIds.end(), X_axons.size(), static_cast<int>(OpenGLWindow::SphereGroup::Axon));

    X.insert(X.end(), X_blood_vessels.begin(), X_blood_vessels.end());
    Y.insert(Y.end(), Y_blood_vessels.begin(), Y_blood_vessels.end());
    Z.insert(Z.end(), Z_blood_vessels.begin(), Z_blood_vessels.end());
    R.insert(R.end(), R_blood_vessels.begin(), R_blood_vessels.end());
    groupIds.insert(groupIds.end(), X_blood_vessels.size(), static_cast<int>(OpenGLWindow::SphereGroup::Blood));

    openglWindow->setSpheres(X, Y, Z, R, groupIds);

    openglWindow->update();
}


void Window::HideAxons(){

    std::vector<int> groupIds = {};

    std::vector<std::vector<double>> X = X_glial_pop1;
    std::vector<std::vector<double>> Y = Y_glial_pop1;
    std::vector<std::vector<double>> Z = Z_glial_pop1;
    std::vector<std::vector<double>> R = R_glial_pop1;
    groupIds.insert(groupIds.end(), X_glial_pop1.size(), static_cast<int>(OpenGLWindow::SphereGroup::Glial1));

    X.insert(X.end(), X_glial_pop2.begin(), X_glial_pop2.end());
    Y.insert(Y.end(), Y_glial_pop2.begin(), Y_glial_pop2.end());
    Z.insert(Z.end(), Z_glial_pop2.begin(), Z_glial_pop2.end());
    R.insert(R.end(), R_glial_pop2.begin(), R_glial_pop2.end());
    groupIds.insert(groupIds.end(), X_glial_pop2.size(), static_cast<int>(OpenGLWindow::SphereGroup::Glial2));

    openglWindow->setSpheres(X, Y, Z, R, groupIds);

    openglWindow->update();
}

void Window::HideGlialCells(){

    std::vector<std::vector<double>> X = X_axons;
    std::vector<std::vector<double>> Y = Y_axons;
    std::vector<std::vector<double>> Z = Z_axons;
    std::vector<std::vector<double>> R = R_axons;
    std::vector<int> groupIds;
    groupIds.insert(groupIds.end(), X_axons.size(), static_cast<int>(OpenGLWindow::SphereGroup::Axon));

    openglWindow->setSpheres(X, Y, Z, R, groupIds);

    openglWindow->update();
}

QGroupBox* Window::createControls(const QString &title)
{
    controlsGroup = new QGroupBox(title);

     // --- Create the GroupBoxes ---
    QGroupBox *generalGroup = new QGroupBox("General Parameters");
    QGroupBox *axonsGroup = new QGroupBox("Axon Parameters");
    QGroupBox *glialGroup1 = new QGroupBox("Glial Cell Population 1 Parameters");
    QGroupBox *glialGroup2 = new QGroupBox("Glial Cell Population 2 Parameters");

    // add ticked box    
    nbr_repetitions_qlabel = new QLabel(tr("Number of Repetitions:"));
    visualise_voxel_qlabel = new QLabel(tr("Visualise Voxel:"));
    axons_icvf_qlabel = new QLabel(tr("Axons ICVF (%):"));
    axons_w_myelin_icvf_qlabel = new QLabel(tr("Axons with myelin ICVF (%):"));
    k1_qlabel = new QLabel(tr("K1 :"));
    k2_qlabel = new QLabel(tr("K2 :"));
    k3_qlabel = new QLabel(tr("K3 :"));
    glial_pop1_soma_icvf_qlabel = new QLabel(tr("Glial Cell somas ICVF (%):"));
    glial_pop1_processes_icvf_qlabel = new QLabel(tr("Glial Cell processes ICVF (%):"));
    glial_pop2_soma_icvf_qlabel = new QLabel(tr("Glial Cell somas ICVF (%):"));
    glial_pop2_processes_icvf_qlabel = new QLabel(tr("Glial Cell processes ICVF (%):"));
    blood_vessels_icvf_qlabel = new QLabel(tr("Blood Vessels ICVF (%):"));
    voxel_size_qlabel = new QLabel(tr("Voxel Edge Length (μm):"));
    minimum_radius_qlabel = new QLabel(tr("Minimum Sphere Radius (μm):"));
    nbr_threads_qlabel = new QLabel(tr("Number of Threads:"));
    overlapping_factor_qlabel = new QLabel(tr("Overlapping Factor (R/f):"));
    c2_qlabel = new QLabel(tr("c2 (fODF):"));
    nbr_axons_populations_qlabel = new QLabel(tr("Number of populations:"));
    epsilon_qlabel = new QLabel(tr("ε (tortuousity):"));
    beading_amplitude_qlabel = new QLabel(tr("Beading Amplitude :"));
    beading_std_qlabel = new QLabel(tr("Beading Standard Deviation :"));
    glial_pop1_mean_process_length_qlabel = new QLabel(tr("Mean Process Length (μm):"));
    glial_pop1_std_process_length_qlabel = new QLabel(tr("Standard Deviation Process Length (μm):"));
    glial_pop2_mean_process_length_qlabel = new QLabel(tr("Mean Process Length (μm):"));
    glial_pop2_std_process_length_qlabel = new QLabel(tr("Standard Deviation Process Length (μm):"));
    alpha_qlabel = new QLabel(tr("α:"));
    beta_qlabel = new QLabel(tr("β:"));
    glial_pop1_radius_mean_qlabel = new QLabel(tr("Glial Cell Soma Radius Mean:"));
    glial_pop1_radius_std_qlabel = new QLabel(tr("Glial Cell Soma Radius Standard Deviation:"));
    glial_pop2_radius_mean_qlabel = new QLabel(tr("Glial Cell Soma Radius Mean:"));
    glial_pop2_radius_std_qlabel = new QLabel(tr("Glial Cell Soma Radius Standard Deviation:"));
    glial_pop1_nbr_primary_processes_qlabel = new QLabel(tr("Number of Primary Processes:"));
    glial_pop2_nbr_primary_processes_qlabel = new QLabel(tr("Number of Primary Processes:"));
    glial_pop1_branching_qlabel = new QLabel(tr("Can Glial Cell Population have branching ? "));
    glial_pop2_branching_qlabel = new QLabel(tr("Can Glial Cell Population have branching ? "));

    nbr_repetitions_SpinBox = new QDoubleSpinBox;
    nbr_repetitions_SpinBox->setRange(1, 100);
    nbr_repetitions_SpinBox->setSingleStep(1);
    nbr_repetitions_SpinBox->setValue(1);

    visualise_voxel_checkbox = new QCheckBox;
    visualise_voxel_checkbox->setChecked(true);

    glial_pop1_branching_checkbox = new QCheckBox;
    glial_pop1_branching_checkbox->setChecked(true);

    glial_pop2_branching_checkbox = new QCheckBox;
    glial_pop2_branching_checkbox->setChecked(true);

    beading_amplitude_SpinBox = new QDoubleSpinBox;
    beading_amplitude_SpinBox->setRange(0, 1);
    beading_amplitude_SpinBox->setSingleStep(0.1);
    beading_amplitude_SpinBox->setValue(0.3);

    beading_std_SpinBox = new QDoubleSpinBox;
    beading_std_SpinBox->setRange(0, 1);
    beading_std_SpinBox->setSingleStep(0.1);
    beading_std_SpinBox->setValue(0.1);

    alpha_SpinBox = new QDoubleSpinBox;
    alpha_SpinBox->setRange(0, 10);
    alpha_SpinBox->setSingleStep(0.1);
    alpha_SpinBox->setValue(4);

    beta_SpinBox = new QDoubleSpinBox;
    beta_SpinBox->setRange(0, 10);
    beta_SpinBox->setSingleStep(0.001);
    beta_SpinBox->setValue(0.25);

    epsilon_SpinBox = new QDoubleSpinBox;
    epsilon_SpinBox->setRange(0, 2);
    epsilon_SpinBox->setSingleStep(0.1);
    epsilon_SpinBox->setValue(0.4);

    glial_pop1_mean_process_length_SpinBox = new QDoubleSpinBox;
    glial_pop1_mean_process_length_SpinBox->setRange(0, 100);
    glial_pop1_mean_process_length_SpinBox->setSingleStep(1);
    glial_pop1_mean_process_length_SpinBox->setValue(0);

    glial_pop1_std_process_length_SpinBox = new QDoubleSpinBox;
    glial_pop1_std_process_length_SpinBox->setRange(0, 100);
    glial_pop1_std_process_length_SpinBox->setSingleStep(1);
    glial_pop1_std_process_length_SpinBox->setValue(15);

    glial_pop2_mean_process_length_SpinBox = new QDoubleSpinBox;
    glial_pop2_mean_process_length_SpinBox->setRange(0, 100);
    glial_pop2_mean_process_length_SpinBox->setSingleStep(1);
    glial_pop2_mean_process_length_SpinBox->setValue(0);

    glial_pop2_std_process_length_SpinBox = new QDoubleSpinBox;
    glial_pop2_std_process_length_SpinBox->setRange(0, 100);
    glial_pop2_std_process_length_SpinBox->setSingleStep(1);
    glial_pop2_std_process_length_SpinBox->setValue(15);


    glial_pop1_nbr_primary_processes_SpinBox = new QDoubleSpinBox;
    glial_pop1_nbr_primary_processes_SpinBox->setRange(1, 20);
    glial_pop1_nbr_primary_processes_SpinBox->setSingleStep(1);
    glial_pop1_nbr_primary_processes_SpinBox->setValue(10);

    glial_pop2_nbr_primary_processes_SpinBox = new QDoubleSpinBox;
    glial_pop2_nbr_primary_processes_SpinBox->setRange(1, 20);
    glial_pop2_nbr_primary_processes_SpinBox->setSingleStep(1);
    glial_pop2_nbr_primary_processes_SpinBox->setValue(10);

    axons_icvf_SpinBox = new QDoubleSpinBox;
    axons_icvf_SpinBox->setRange(0, 100);
    axons_icvf_SpinBox->setSingleStep(1);

    axons_w_myelin_icvf_SpinBox = new QDoubleSpinBox;
    axons_w_myelin_icvf_SpinBox->setRange(0, 100);
    axons_w_myelin_icvf_SpinBox->setSingleStep(1);

    blood_vessels_icvf_SpinBox = new QDoubleSpinBox;
    blood_vessels_icvf_SpinBox->setRange(0, 100);
    blood_vessels_icvf_SpinBox->setSingleStep(1);

    glial_pop1_soma_icvf_SpinBox = new QDoubleSpinBox;
    glial_pop1_soma_icvf_SpinBox->setRange(0, 100);
    glial_pop1_soma_icvf_SpinBox->setSingleStep(1);

    glial_pop1_processes_icvf_SpinBox = new QDoubleSpinBox;
    glial_pop1_processes_icvf_SpinBox->setRange(0, 100);
    glial_pop1_processes_icvf_SpinBox->setSingleStep(1);

    glial_pop2_soma_icvf_SpinBox = new QDoubleSpinBox;
    glial_pop2_soma_icvf_SpinBox->setRange(0, 100);
    glial_pop2_soma_icvf_SpinBox->setSingleStep(1);

    glial_pop2_processes_icvf_SpinBox = new QDoubleSpinBox;
    glial_pop2_processes_icvf_SpinBox->setRange(0, 100);
    glial_pop2_processes_icvf_SpinBox->setSingleStep(1);

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

    glial_pop1_radius_mean_SpinBox = new QDoubleSpinBox;
    glial_pop1_radius_mean_SpinBox->setRange(0, 10);
    glial_pop1_radius_mean_SpinBox->setSingleStep(0.1);
    glial_pop1_radius_mean_SpinBox->setValue(5);

    glial_pop1_radius_std_SpinBox = new QDoubleSpinBox;
    glial_pop1_radius_std_SpinBox->setRange(0, 10);
    glial_pop1_radius_std_SpinBox->setSingleStep(0.1);
    glial_pop1_radius_std_SpinBox->setValue(0.5);

    glial_pop2_radius_mean_SpinBox = new QDoubleSpinBox;
    glial_pop2_radius_mean_SpinBox->setRange(0, 10);
    glial_pop2_radius_mean_SpinBox->setSingleStep(0.1);
    glial_pop2_radius_mean_SpinBox->setValue(5);

    glial_pop2_radius_std_SpinBox = new QDoubleSpinBox;
    glial_pop2_radius_std_SpinBox->setRange(0, 10);
    glial_pop2_radius_std_SpinBox->setSingleStep(0.1);
    glial_pop2_radius_std_SpinBox->setValue(0.5);

    k1_SpinBox = new QDoubleSpinBox;
    k1_SpinBox->setRange(0, 10);
    k1_SpinBox->setSingleStep(0.05);
    k1_SpinBox->setValue(0.35);
    k1_SpinBox->setDecimals(3); // Set at least 3 decimals

    k2_SpinBox = new QDoubleSpinBox;
    k2_SpinBox->setRange(0, 10);
    k2_SpinBox->setSingleStep(0.001);
    k2_SpinBox->setValue(0.006);
    k2_SpinBox->setDecimals(4); // Needed for values like 0.006

    k3_SpinBox = new QDoubleSpinBox;
    k3_SpinBox->setRange(0, 10);
    k3_SpinBox->setSingleStep(0.001);
    k3_SpinBox->setValue(0.024);
    k3_SpinBox->setDecimals(4); // Shows 0.024 cleanly

    
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

    std::vector<QLabel*> general_labels = { nbr_repetitions_qlabel, voxel_size_qlabel, overlapping_factor_qlabel, minimum_radius_qlabel, blood_vessels_icvf_qlabel};
    std::vector <QDoubleSpinBox*> general_spinBoxes = { nbr_repetitions_SpinBox, voxel_size_SpinBox, overlapping_factor_SpinBox, minimum_radius_SpinBox, blood_vessels_icvf_SpinBox};
    
    std::vector<QLabel*> axons_labels = {axons_w_myelin_icvf_qlabel, k1_qlabel, k2_qlabel, k3_qlabel, axons_icvf_qlabel ,nbr_threads_qlabel, epsilon_qlabel, c2_qlabel, nbr_axons_populations_qlabel, beading_amplitude_qlabel, beading_std_qlabel, alpha_qlabel, beta_qlabel};
    std::vector <QDoubleSpinBox*> axons_spinBoxes = {axons_w_myelin_icvf_SpinBox, k1_SpinBox, k2_SpinBox, k3_SpinBox, axons_icvf_SpinBox, nbr_threads_SpinBox, epsilon_SpinBox, c2_SpinBox, nbr_axons_populations_SpinBox, beading_amplitude_SpinBox, beading_std_SpinBox, alpha_SpinBox, beta_SpinBox};
    
    std::vector<QLabel*> glials_labels1 = {glial_pop1_soma_icvf_qlabel, glial_pop1_processes_icvf_qlabel, glial_pop1_radius_mean_qlabel, glial_pop1_radius_std_qlabel, glial_pop1_mean_process_length_qlabel, glial_pop1_std_process_length_qlabel, glial_pop1_nbr_primary_processes_qlabel};
    std::vector <QDoubleSpinBox*> glials_spinBoxes1 = {glial_pop1_soma_icvf_SpinBox, glial_pop1_processes_icvf_SpinBox, glial_pop1_radius_mean_SpinBox, glial_pop1_radius_std_SpinBox, glial_pop1_mean_process_length_SpinBox, glial_pop1_std_process_length_SpinBox, glial_pop1_nbr_primary_processes_SpinBox};
    
    std::vector<QLabel*> glials_labels2 = {glial_pop2_soma_icvf_qlabel, glial_pop2_processes_icvf_qlabel, glial_pop2_radius_mean_qlabel, glial_pop2_radius_std_qlabel, glial_pop2_mean_process_length_qlabel, glial_pop2_std_process_length_qlabel, glial_pop2_nbr_primary_processes_qlabel};
    std::vector <QDoubleSpinBox*> glials_spinBoxes2 = {glial_pop2_soma_icvf_SpinBox, glial_pop2_processes_icvf_SpinBox, glial_pop2_radius_mean_SpinBox, glial_pop2_radius_std_SpinBox, glial_pop2_mean_process_length_SpinBox, glial_pop2_std_process_length_SpinBox, glial_pop2_nbr_primary_processes_SpinBox};

    
    QGridLayout *generalLayout = new QGridLayout;

    for (int i = 0; i < general_labels.size(); i++){
        generalLayout->addWidget(general_labels[i], i, 0);
        generalLayout->addWidget(general_spinBoxes[i], i, 1);
    }
    generalLayout->addWidget(visualise_voxel_qlabel, general_labels.size(), 0);
    generalLayout->addWidget(visualise_voxel_checkbox, general_labels.size(), 1);

    generalGroup->setLayout(generalLayout);

    QGridLayout *axonsLayout = new QGridLayout;
    int row = 0;
    
    for (int i = 0; i < axons_labels.size(); i++) {
        if (i == 1) {
            // K1, K2, K3 all in the same row
            axonsLayout->addWidget(axons_labels[1], row, 0); // K1
            axonsLayout->addWidget(axons_spinBoxes[1], row, 1);

            axonsLayout->addWidget(axons_labels[2], row, 2); // K2
            axonsLayout->addWidget(axons_spinBoxes[2], row, 3);

            axonsLayout->addWidget(axons_labels[3], row, 4); // K3
            axonsLayout->addWidget(axons_spinBoxes[3], row, 5);

            row++;

            // Add formula below
            QLabel *formulaLabel = new QLabel("Myelin thickness = K1 + K2 × Inner diameter + K3 × log(Inner diameter)");
            QFont formulaFont = formulaLabel->font();
            formulaFont.setItalic(true);
            formulaLabel->setFont(formulaFont);
            formulaLabel->setAlignment(Qt::AlignCenter);
            axonsLayout->addWidget(formulaLabel, row, 0, 1, 6); // Span all 6 columns
            row++;
            i = 3; // Skip already handled K1, K2, K3
        } 
        else if (i == axons_labels.size()-2) {
            QLabel *formulaLabel = new QLabel("Gamma Distribution parameters for inner radii : ");
            QFont formulaFont = formulaLabel->font();
            formulaLabel->setFont(formulaFont);
            axonsLayout->addWidget(formulaLabel, row, 0, 1, 6); // Span all 6 columns
            row++;
            axonsLayout->addWidget(axons_labels[i], row, 0);
            axonsLayout->addWidget(axons_spinBoxes[i], row, 1);
            row++;
        }
        
        else {
            axonsLayout->addWidget(axons_labels[i], row, 0);
            axonsLayout->addWidget(axons_spinBoxes[i], row, 1);
            row++;
        }
    }
    axonsGroup->setLayout(axonsLayout);

    QGridLayout *glialLayout1 = new QGridLayout;

    for (int i = 0; i < glials_labels1.size(); i++){

        glialLayout1->addWidget(glials_labels1[i], i, 0);
        glialLayout1->addWidget(glials_spinBoxes1[i], i, 1);
    }


    glialLayout1->addWidget(glial_pop1_branching_qlabel, glials_labels1.size(), 0);
    glialLayout1->addWidget(glial_pop1_branching_checkbox, glials_labels1.size(), 1);

    glialGroup1->setLayout(glialLayout1);

    QGridLayout *glialLayout2 = new QGridLayout;
    for (int i = 0; i < glials_labels2.size(); i++){
        glialLayout2->addWidget(glials_labels2[i], i, 0);
        glialLayout2->addWidget(glials_spinBoxes2[i], i, 1);
    }
    glialLayout2->addWidget(glial_pop2_branching_qlabel, glials_labels2.size(), 0);
    glialLayout2->addWidget(glial_pop2_branching_checkbox, glials_labels2.size(), 1);

    glialGroup2->setLayout(glialLayout2);

    // Arrange groups in a 2x2 grid within `controlsGroup`
    QGridLayout *mainControlsLayout = new QGridLayout;
    mainControlsLayout->addWidget(generalGroup, 0, 0);
    mainControlsLayout->addWidget(axonsGroup, 0, 1);
    mainControlsLayout->addWidget(glialGroup1, 1, 0);
    mainControlsLayout->addWidget(glialGroup2, 1, 1);

    controlsGroup->setLayout(mainControlsLayout);

    return controlsGroup;


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
    updateMembersFromGui();  
    
    // Close the parameter input dialog
    this->close();

    StartSimulation();

}

void Window::updateMembersFromGui()
{
    // Retrieve values from spin boxes and checkboxes
    nbr_repetitions = nbr_repetitions_SpinBox->value();
    axons_icvf = axons_icvf_SpinBox->value();
    axons_w_myelin_icvf = axons_w_myelin_icvf_SpinBox->value();
    blood_vessels_icvf = blood_vessels_icvf_SpinBox->value();
    k1 = k1_SpinBox->value();
    k2 = k2_SpinBox->value();
    k3 = k3_SpinBox->value();
    glial_pop1_soma_icvf = glial_pop1_soma_icvf_SpinBox->value();
    glial_pop1_processes_icvf = glial_pop1_processes_icvf_SpinBox->value();
    glial_pop2_soma_icvf = glial_pop2_soma_icvf_SpinBox->value();
    glial_pop2_processes_icvf = glial_pop2_processes_icvf_SpinBox->value();
    nbr_threads = nbr_threads_SpinBox->value();
    overlapping_factor = overlapping_factor_SpinBox->value();
    voxel_size = voxel_size_SpinBox->value();
    minimum_radius = minimum_radius_SpinBox->value();
    c2 = c2_SpinBox->value();
    nbr_axons_populations = nbr_axons_populations_SpinBox->value();
    beading_amplitude = beading_amplitude_SpinBox->value();
    beading_std = beading_std_SpinBox->value();
    glial_pop1_mean_process_length = glial_pop1_mean_process_length_SpinBox->value();
    glial_pop1_std_process_length = glial_pop1_std_process_length_SpinBox->value();
    glial_pop2_mean_process_length = glial_pop2_mean_process_length_SpinBox->value();
    glial_pop2_std_process_length = glial_pop2_std_process_length_SpinBox->value();

    epsilon = epsilon_SpinBox->value();
    minimum_radius = minimum_radius_SpinBox->value();
    alpha = alpha_SpinBox->value();
    beta = beta_SpinBox->value();
    visualise_voxel = visualise_voxel_checkbox->isChecked();
    glial_pop1_radius_mean = glial_pop1_radius_mean_SpinBox->value();
    glial_pop1_radius_std = glial_pop1_radius_std_SpinBox->value();
    glial_pop2_radius_mean = glial_pop2_radius_mean_SpinBox->value();
    glial_pop2_radius_std = glial_pop2_radius_std_SpinBox->value();
    glial_pop1_nbr_primary_processes = glial_pop1_nbr_primary_processes_SpinBox->value();
    glial_pop2_nbr_primary_processes = glial_pop2_nbr_primary_processes_SpinBox->value();
    glial_pop1_branching = glial_pop1_branching_checkbox->isChecked();
    glial_pop2_branching = glial_pop2_branching_checkbox->isChecked();

}


void Window::ReadAxonsFromFile(const QString& fileName){

    if (fileName.endsWith(".csv")) {
        ReadAxonsFromCSV(fileName);
    } else if (fileName.endsWith(".swc")) {
        ReadAxonsFromSWC(fileName);
    } else {
        QMessageBox::warning(this, tr("Error"), tr("Unsupported file format. Please select a CSV or SWC file."));
    }
}
void Window::ReadAxonsFromCSV(const QString& fileName){
    X_axons.clear();
    Y_axons.clear();
    Z_axons.clear();
    R_axons.clear();

    if (fileName.isEmpty()) {
        return; // User canceled the file dialog
    }


    std::ifstream swcFile(fileName.toStdString());
    if (!swcFile.is_open()) {
        QMessageBox::warning(this, tr("Error"), tr("Could not open the CSV file."));
        return;
    }

    std::vector<double> x_ = {};
    std::vector<double> y_ = {};
    std::vector<double> z_ = {};
    std::vector<double> r_ = {};

    int old_cell_id = -1;

    std::string line;
    while (std::getline(swcFile, line)) {
        std::istringstream iss(line);
        double id_cell, component_id;
        std::string type, component;
        double x, y, z, radius_in, parent, radius_out;

        //skip first line
        if (line[0] == 'c') {
            continue;
        }


        if (!(iss >> type >> id_cell >> component >> component_id >> x >> y >> z >> radius_in >> radius_out)) {
            QMessageBox::warning(this, tr("Error"), tr("Invalid CSV file format for Axons."));
            return;
        }
        
        if (type == "axon") {
            
            if (old_cell_id != id_cell) {
                if (old_cell_id != -1){
                    if (x_.size() > 0) {
                        X_axons.push_back(x_);
                        Y_axons.push_back(y_);
                        Z_axons.push_back(z_);
                        R_axons.push_back(r_);
                    }
                }
                x_.clear();
                y_.clear();
                z_.clear();
                r_.clear();
            }

            x_.push_back(x);
            y_.push_back(y);
            z_.push_back(z);
            r_.push_back(radius_out);
            old_cell_id = id_cell;
        }
    }
    // Add the last axon
    if (x_.size() > 0) {
        X_axons.push_back(x_);
        Y_axons.push_back(y_);
        Z_axons.push_back(z_);
        R_axons.push_back(r_);
    }

    swcFile.close();
}

void Window::ReadAxonsFromSWC(const QString& fileName){

    X_axons.clear();
    Y_axons.clear();
    Z_axons.clear();
    R_axons.clear();

    if (fileName.isEmpty()) {
        return; // User canceled the file dialog
    }


    std::ifstream swcFile(fileName.toStdString());
    if (!swcFile.is_open()) {
        QMessageBox::warning(this, tr("Error"), tr("Could not open the SWC file."));
        return;
    }

    std::vector<double> x_ = {};
    std::vector<double> y_ = {};
    std::vector<double> z_ = {};
    std::vector<double> r_ = {};

    std::string line;
    while (std::getline(swcFile, line)) {
        std::istringstream iss(line);
        int id_branch;
        double id_cell, id_sphere;
        std::string type;
        double x, y, z, radius_in, parent, radius_out;

        //skip first line
        if (line[0] == 'i'|| line[0] == 'a') {
            continue;
        }

        if (!(iss >> id_cell >> id_sphere >> id_branch >> type >> x >> y >> z >> radius_in >> radius_out >> parent)) {
            QMessageBox::warning(this, tr("Error"), tr("Invalid SWC file format for Axons."));
            return;
        }
        
        if (type == "axon") {
            
            if (id_sphere == 0) {
                if (x_.size() > 0) {
                    X_axons.push_back(x_);
                    Y_axons.push_back(y_);
                    Z_axons.push_back(z_);
                    R_axons.push_back(r_);
                }
                x_.clear();
                y_.clear();
                z_.clear();
                r_.clear();
                x_.push_back(x);
                y_.push_back(y);
                z_.push_back(z);
                r_.push_back(radius_out);
            }
            else {
                x_.push_back(x);
                y_.push_back(y);
                z_.push_back(z);
                r_.push_back(radius_out);
            }
        }
    }
    // Add the last axon
    if (x_.size() > 0) {
        X_axons.push_back(x_);
        Y_axons.push_back(y_);
        Z_axons.push_back(z_);
        R_axons.push_back(r_);
    }

    swcFile.close();
}
void Window::ReadGlialCellsFromFile(const QString& fileName){
    if (fileName.endsWith(".csv")) {
        ReadGlialCellsFromCSV(fileName);
    } else if (fileName.endsWith(".swc")) {
        ReadGlialCellsFromSWC(fileName);
    } else {
        QMessageBox::warning(this, tr("Error"), tr("Unsupported file format. Please select a CSV or SWC file."));
    }
}
void Window::ReadGlialCellsFromCSV(const QString& fileName){

    X_glial_pop1.clear();
    Y_glial_pop1.clear();
    Z_glial_pop1.clear();
    R_glial_pop1.clear();
    Branch_glial_pop1.clear();

    if (fileName.isEmpty()) {
        return; // User canceled the file dialog
    }

    std::ifstream swcFile(fileName.toStdString());
    if (!swcFile.is_open()) {
        QMessageBox::warning(this, tr("Error"), tr("Could not open the CSV file."));
        assert(0);
        return;
    }

    int previous_branch_id = -2;
    int id_previous_cell = -1;

    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> z_;
    std::vector<double> r_;
    std::vector<int> b_;

    std::string line;
    while (std::getline(swcFile, line)) {

        //skip first line
        if (line[0] == 'c') {
            continue;
        }

        std::istringstream iss(line);
        double component_id, id_cell;
        std::string type, component;

        double x, y, z, radius_in, parent, radius_out;

        if (!(iss >> type >> id_cell >> component >> component_id >> x >> y >> z >> radius_in >> radius_out)) {
            QMessageBox::warning(this, tr("Error"), tr("Invalid SWC file format for Glial Cells."));
            return;
        }

        if (type == "glial_cell") {

            if (id_cell != id_previous_cell) {
                if (id_previous_cell != -1) {
                    X_glial_pop1.push_back(x_);
                    Y_glial_pop1.push_back(y_);
                    Z_glial_pop1.push_back(z_);
                    R_glial_pop1.push_back(r_);
                    Branch_glial_pop1.push_back(b_);

                }
                x_.clear();
                y_.clear();
                z_.clear();
                r_.clear();
                b_.clear();
            }
            x_.push_back(x);
            y_.push_back(y);
            z_.push_back(z);
            r_.push_back(radius_out);
            b_.push_back(component_id);
            id_previous_cell = id_cell;
        }
        
    }

    if (!x_.empty()) {
        X_glial_pop1.push_back(x_);
        Y_glial_pop1.push_back(y_);
        Z_glial_pop1.push_back(z_);
        R_glial_pop1.push_back(r_);
        Branch_glial_pop1.push_back(b_);
    }

    swcFile.close();
}

void Window::ReadGlialCellsFromSWC(const QString& fileName){

    X_glial_pop1.clear();
    Y_glial_pop1.clear();
    Z_glial_pop1.clear();
    R_glial_pop1.clear();
    Branch_glial_pop1.clear();

    if (fileName.isEmpty()) {
        return; // User canceled the file dialog
    }

    std::ifstream swcFile(fileName.toStdString());
    if (!swcFile.is_open()) {
        QMessageBox::warning(this, tr("Error"), tr("Could not open the SWC file."));
        assert(0);
        return;
    }

    int previous_branch_id = -2;
    int id_previous_cell = -1;

    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> z_;
    std::vector<double> r_;
    std::vector<int> b_;

    std::string line;
    while (std::getline(swcFile, line)) {

        //skip first line
        if (line[0] == 'i' || line[0] == 'a') {
            continue;
        }

        std::istringstream iss(line);
        double id_cell, id_sphere;
        int id_branch;
        std::string type;

        double x, y, z, radius_in, parent, radius_out;

        if (!(iss >> id_cell >> id_sphere >> id_branch >> type >> x >> y >> z >> radius_in >> radius_out >> parent)) {
            QMessageBox::warning(this, tr("Error"), tr("Invalid SWC file format for Glial Cells."));
            return;
        }

        if (type == "Process" || type == "CellSoma" || type == "glial_cell") {

            if (id_cell != id_previous_cell) {
                if (id_previous_cell != -1) {
                    X_glial_pop1.push_back(x_);
                    Y_glial_pop1.push_back(y_);
                    Z_glial_pop1.push_back(z_);
                    R_glial_pop1.push_back(r_);
                    Branch_glial_pop1.push_back(b_);

                }
                x_.clear();
                y_.clear();
                z_.clear();
                r_.clear();
                b_.clear();
            }
            x_.push_back(x);
            y_.push_back(y);
            z_.push_back(z);
            r_.push_back(radius_out);
            b_.push_back(id_branch);
            id_previous_cell = id_cell;
        }
        
    }

    if (!x_.empty()) {
        X_glial_pop1.push_back(x_);
        Y_glial_pop1.push_back(y_);
        Z_glial_pop1.push_back(z_);
        R_glial_pop1.push_back(r_);
        Branch_glial_pop1.push_back(b_);
    }

    swcFile.close();
}


void Window::ReadBloodVesselsFromFile(const QString& fileName){

    X_blood_vessels.clear();
    Y_blood_vessels.clear();
    Z_blood_vessels.clear();
    R_blood_vessels.clear();

    if (fileName.isEmpty()) {
        return; // User canceled the file dialog
    }


    std::ifstream swcFile(fileName.toStdString());
    if (!swcFile.is_open()) {
        QMessageBox::warning(this, tr("Error"), tr("Could not open the SWC file."));
        return;
    }

    std::vector<double> x_ = {};
    std::vector<double> y_ = {};
    std::vector<double> z_ = {};
    std::vector<double> r_ = {};

    std::string line;
    while (std::getline(swcFile, line)) {
        std::istringstream iss(line);
        int id_branch;
        double id_cell, id_sphere;
        std::string type;
        double x, y, z, radius_in, parent, radius_out;

        //skip first line
        if (line[0] == 'i'|| line[0] == 'a') {
            continue;
        }

        if (!(iss >> id_cell >> id_sphere >> id_branch >> type >> x >> y >> z >> radius_in >> radius_out >> parent)) {
            QMessageBox::warning(this, tr("Error"), tr("Invalid SWC file format for Axons."));
            return;
        }
        
        if (type == "blood_vessel") {
            
            if (id_sphere == 0) {
                if (x_.size() > 0) {
                    X_blood_vessels.push_back(x_);
                    Y_blood_vessels.push_back(y_);
                    Z_blood_vessels.push_back(z_);
                    R_blood_vessels.push_back(r_);
                }
                x_.clear();
                y_.clear();
                z_.clear();
                r_.clear();
                x_.push_back(x);
                y_.push_back(y);
                z_.push_back(z);
                r_.push_back(radius_out);
            }
            else {
                x_.push_back(x);
                y_.push_back(y);
                z_.push_back(z);
                r_.push_back(radius_out);
            }
        }
    }
    // Add the last axon
    if (x_.size() > 0) {
        X_blood_vessels.push_back(x_);
        Y_blood_vessels.push_back(y_);
        Z_blood_vessels.push_back(z_);
        R_blood_vessels.push_back(r_);
    }

    swcFile.close();
}

// Function to check if a point is inside a dilated box
bool Window::check_borders(const Eigen::Vector3d&  min_l, const Eigen::Vector3d&  max_l, const Eigen::Vector3d& pos, const double& distance_to_border) {

    
    // Check if the point is inside the dilated box
    for (int i = 0; i < 3; ++i) {
        double min_bound = min_l[i] - distance_to_border;
        double max_bound = max_l[i] + distance_to_border;
        if (pos[i] < min_bound || pos[i] > max_bound) {
            return false; // Point is outside the dilated box
        }
    }
    
    return true; // Point is inside the dilated box
}


void Window::StartSimulation(){

    Parameters parameters;

    int repetitions = nbr_repetitions;
    std::vector<double> vox_sizes = {voxel_size/1.0};
    double glial_pop1_icvf_soma = glial_pop1_soma_icvf/100.0;
    double glial_pop1_icvf_processes = glial_pop1_processes_icvf/100.0;
    double glial_pop2_icvf_soma = glial_pop2_soma_icvf/100.0;
    double glial_pop2_icvf_processes = glial_pop2_processes_icvf/100.0;
    double axons_wo_myelin_icvf = axons_icvf/100.0;
    double axons_with_myelin_icvf = axons_w_myelin_icvf/100.0;
    double blood_vessels_icvf_ = double(blood_vessels_icvf/100.0);
    int spheres_overlap_factor = overlapping_factor;
    std::string directory = selectedDirectory.toStdString();
    double cosPhiSquared = c2;
    double std_dev = epsilon;
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
            // create distribution of axons

            AxonGammaDistribution AxonDistribution = AxonGammaDistribution(axons_wo_myelin_icvf, axons_with_myelin_icvf, glial_pop1_icvf_soma, glial_pop1_icvf_processes, glial_pop2_icvf_soma, glial_pop2_icvf_processes, blood_vessels_icvf_, alpha, beta,
                                             min_l, max_l, min_radius, regrow_thr, beading_amplitude, beading_std, std_dev, ondulation_factor, spheres_overlap_factor, can_shrink, cosPhiSquared, nbr_threads, nbr_axons_populations, crossing_fibers_type, 
                                              glial_pop1_mean_process_length, glial_pop1_std_process_length, glial_pop2_mean_process_length, glial_pop2_std_process_length,
                                              glial_pop1_radius_mean, glial_pop1_radius_std, glial_pop2_radius_mean, glial_pop2_radius_std, glial_pop1_branching, glial_pop2_branching, glial_pop1_nbr_primary_processes, glial_pop2_nbr_primary_processes, k1, k2, k3);
            AxonDistribution.createSubstrate();
            // saving spheres
            if (rep == 0){
                for (unsigned i=0; i< AxonDistribution.axons.size(); ++i){
                    std::vector<double> x_;
                    std::vector<double> y_;
                    std::vector<double> z_;
                    std::vector<double> r_;
                    for (unsigned j=0; j< AxonDistribution.axons[i].outer_spheres.size(); ++j){

                        double _x_ = AxonDistribution.axons[i].outer_spheres[j].center[0];
                        double _y_ = AxonDistribution.axons[i].outer_spheres[j].center[1];
                        double _z_ = AxonDistribution.axons[i].outer_spheres[j].center[2];
                        double _r_ = AxonDistribution.axons[i].outer_spheres[j].radius;
                        Eigen::Vector3d pos = {_x_, _y_, _z_};

                        if (check_borders(min_l, max_l, pos, 0.0)) {
                            x_.push_back(_x_);
                            y_.push_back(_y_);
                            z_.push_back(_z_);
                            r_.push_back(_r_);
                        }
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

                for (unsigned i=0; i< AxonDistribution.glial_pop1.size(); ++i){
                    std::vector<double> x_;
                    std::vector<double> y_;
                    std::vector<double> z_;
                    std::vector<double> r_;
                    std::vector<int> b_;

                    double _x_ = AxonDistribution.glial_pop1[i].soma.center[0];
                    double _y_ = AxonDistribution.glial_pop1[i].soma.center[1];
                    double _z_ = AxonDistribution.glial_pop1[i].soma.center[2];
                    double _r_ = AxonDistribution.glial_pop1[i].soma.radius;
                    Eigen::Vector3d pos = {_x_, _y_, _z_};

                    if (check_borders(min_l, max_l, pos, 0.0)) {
                        x_.push_back(_x_);
                        y_.push_back(_y_);
                        z_.push_back(_z_);
                        r_.push_back(_r_);
                        b_.push_back(0);
                    }


                    for (unsigned j=0; j< AxonDistribution.glial_pop1[i].ramification_spheres.size(); ++j){

                        for (unsigned k=0; k< AxonDistribution.glial_pop1[i].ramification_spheres[j].size(); ++k){

                            double _x_ = AxonDistribution.glial_pop1[i].ramification_spheres[j][k].center[0];
                            double _y_ = AxonDistribution.glial_pop1[i].ramification_spheres[j][k].center[1];
                            double _z_ = AxonDistribution.glial_pop1[i].ramification_spheres[j][k].center[2];
                            double _r_ = AxonDistribution.glial_pop1[i].ramification_spheres[j][k].radius;
                            Eigen::Vector3d pos = {_x_, _y_, _z_};
                            if (!check_borders(min_l, max_l, pos, 0.0)) {
                                continue;
                            }
                            x_.push_back(_x_);
                            y_.push_back(_y_);
                            z_.push_back(_z_);
                            r_.push_back(_r_);
                            b_.push_back(j);
                        }
                    }

                    X_glial_pop1.push_back(x_);
                    Y_glial_pop1.push_back(y_);
                    Z_glial_pop1.push_back(z_);
                    R_glial_pop1.push_back(r_);
                    Branch_glial_pop1.push_back(b_);
                    x_.clear();
                    y_.clear();
                    z_.clear();
                    r_.clear();
                    b_.clear();
                }

                for (unsigned i=0; i< AxonDistribution.glial_pop2.size(); ++i){
                    std::vector<double> x_;
                    std::vector<double> y_;
                    std::vector<double> z_;
                    std::vector<double> r_;
                    std::vector<int> b_;

                    double _x_ = AxonDistribution.glial_pop2[i].soma.center[0];
                    double _y_ = AxonDistribution.glial_pop2[i].soma.center[1];
                    double _z_ = AxonDistribution.glial_pop2[i].soma.center[2];
                    double _r_ = AxonDistribution.glial_pop2[i].soma.radius;

                    if (check_borders(min_l, max_l, {_x_, _y_, _z_}, 0.0)) {
                        x_.push_back(_x_);
                        y_.push_back(_y_);
                        z_.push_back(_z_);
                        r_.push_back(_r_);
                        b_.push_back(0);
                    }


                    for (unsigned j=0; j< AxonDistribution.glial_pop2[i].ramification_spheres.size(); ++j){

                        for (unsigned k=0; k< AxonDistribution.glial_pop2[i].ramification_spheres[j].size(); ++k){

                            double _x_ = AxonDistribution.glial_pop2[i].ramification_spheres[j][k].center[0];
                            double _y_ = AxonDistribution.glial_pop2[i].ramification_spheres[j][k].center[1];
                            double _z_ = AxonDistribution.glial_pop2[i].ramification_spheres[j][k].center[2];
                            double _r_ = AxonDistribution.glial_pop2[i].ramification_spheres[j][k].radius;

                            if (!check_borders(min_l, max_l, {_x_, _y_, _z_}, 0.0)) {
                                continue;
                            }
                            x_.push_back(_x_);
                            y_.push_back(_y_);
                            z_.push_back(_z_);
                            r_.push_back(_r_);
                            b_.push_back(j);
                        }
                    }

                    X_glial_pop2.push_back(x_);
                    Y_glial_pop2.push_back(y_);
                    Z_glial_pop2.push_back(z_);
                    R_glial_pop2.push_back(r_);
                    Branch_glial_pop2.push_back(b_);
                    x_.clear();
                    y_.clear();
                    z_.clear();
                    r_.clear();
                    b_.clear();
                }

                for (unsigned i=0; i< AxonDistribution.blood_vessels.size(); ++i){
                    std::vector<double> x_;
                    std::vector<double> y_;
                    std::vector<double> z_;
                    std::vector<double> r_;
                    for (unsigned j=0; j< AxonDistribution.blood_vessels[i].spheres.size(); ++j){

                        double _x_ = AxonDistribution.blood_vessels[i].spheres[j].center[0];
                        double _y_ = AxonDistribution.blood_vessels[i].spheres[j].center[1];
                        double _z_ = AxonDistribution.blood_vessels[i].spheres[j].center[2];
                        double _r_ = AxonDistribution.blood_vessels[i].spheres[j].radius;
                        Eigen::Vector3d pos = {_x_, _y_, _z_};

                        if (check_borders(min_l, max_l, pos, 0.0)) {
                            x_.push_back(_x_);
                            y_.push_back(_y_);
                            z_.push_back(_z_);
                            r_.push_back(_r_);
                        }
                    }
                    X_blood_vessels.push_back(x_);
                    Y_blood_vessels.push_back(y_);
                    Z_blood_vessels.push_back(z_);
                    R_blood_vessels.push_back(r_);
                    x_.clear();
                    y_.clear();
                    z_.clear();
                    r_.clear();
                }
            }

            std::string simulation_file_name;
            std::string swc_file_name;

            std::string filename;
            std::ifstream file(filename);

         if (rep ==0){
            simulation_file_name = (directory + "/growth_info.txt");
            swc_file_name = (directory + "/Voxel.csv");
            }
        else{
            simulation_file_name = (directory + "/growth_info_" + std::to_string(rep) + ".txt");
            swc_file_name = (directory + "/Voxel_" + std::to_string(rep) + ".csv");
        }
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
        PlotCells(true, true, true, true);
    }
    else{
        // Display a message box to inform the user that the simulation is complete
        QMessageBox::information(this, "Simulation Complete", "Simulation complete! Please check the output directory for the results.");
    }

    // Save config file with used parameters for generatign voxel
    QString baseDir = selectedDirectory.isEmpty()
                  ? QDir::currentPath()
                  : selectedDirectory;

    QString cfgPath = QDir(baseDir).filePath("Voxel_config.json");
    saveConfigToPath(cfgPath);


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
    int binCount = 50;  // Number of bins, adjust this as needed
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

    if (X_glial_pop1.size() == 0) {
        QMessageBox::warning(this, "Error", "No Glial cells to plot!");
        return;
    }

    // Radii for Sholl analysis
    std::vector<double> sphere_around_soma_radii = {5, 7, 10, 15, 20, 25, 30, 40, 50, 60, 80};
    std::vector<double> mean_intersections(sphere_around_soma_radii.size(), 0);

    for (unsigned long i = 0; i < X_glial_pop1.size(); ++i) {
        // Soma position of the current glial_pop1
        Eigen::Vector3d soma_position = {X_glial_pop1[i][0], Y_glial_pop1[i][0], Z_glial_pop1[i][0]};
        std::vector<double> intersections_list(sphere_around_soma_radii.size(), 0);
        std::vector<int> branches_list;

        // Iterate through all spheres (excluding the soma) to compute intersections
        for (unsigned long r = 0; r < sphere_around_soma_radii.size(); ++r) {
            for (unsigned long j = 1; j < X_glial_pop1[i].size(); ++j) {
                Eigen::Vector3d position = {X_glial_pop1[i][j], Y_glial_pop1[i][j], Z_glial_pop1[i][j]};
                double distance = (position - soma_position).norm();
                
                if (distance < sphere_around_soma_radii[r] + R_glial_pop1[i][j] && distance > sphere_around_soma_radii[r] - R_glial_pop1[i][j]) {
                    if (std::find(branches_list.begin(), branches_list.end(), Branch_glial_pop1[i][j]) == branches_list.end()) {
                        intersections_list[r] += 1;
                        branches_list.push_back(Branch_glial_pop1[i][j]);
                    }
                }
            }
            branches_list.clear();
        }

        // Accumulate values for mean calculation
        for (size_t r = 0; r < sphere_around_soma_radii.size(); ++r) {
            mean_intersections[r] += intersections_list[r];
        }
    }

    // Compute mean intersections
    for (size_t r = 0; r < mean_intersections.size(); ++r) {
        mean_intersections[r] /= X_glial_pop1.size();
    }

    // Create QCustomPlot for mean Sholl analysis
    QCustomPlot *customPlot = new QCustomPlot;

    // Convert the data to QVector for QCustomPlot
    QVector<double> x = QVector<double>::fromStdVector(sphere_around_soma_radii);
    QVector<double> y = QVector<double>::fromStdVector(mean_intersections);

    // Create a graph and set the data
    customPlot->addGraph();
    customPlot->graph(0)->setData(x, y);
    customPlot->graph(0)->setLineStyle(QCPGraph::lsLine);
    customPlot->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 5));

    // Set axis labels
    customPlot->xAxis->setLabel("Distance to Soma (μm)");
    customPlot->yAxis->setLabel("Mean Number of Intersections");

    // Set axis ranges
    customPlot->xAxis->setRange(0, *std::max_element(sphere_around_soma_radii.begin(), sphere_around_soma_radii.end()));
    customPlot->yAxis->setRange(0, *std::max_element(mean_intersections.begin(), mean_intersections.end()));

    // Display the plot in a dialog window
    QDialog *dialog = new QDialog(this);
    dialog->resize(800, 600);
    QVBoxLayout *layout = new QVBoxLayout;
    layout->addWidget(customPlot);
    dialog->setLayout(layout);
    dialog->setWindowTitle("Mean Sholl Analysis for glial_pop1");
    dialog->exec();
}

void Window::loadConfigFromPath(const QString &fileName)
{
    QFile f(fileName);
    if (!f.open(QIODevice::ReadOnly)) {
        qWarning() << "Cannot open config file" << fileName;
        return;
    }

    QJsonParseError err;
    QJsonDocument doc = QJsonDocument::fromJson(f.readAll(), &err);
    if (err.error != QJsonParseError::NoError || !doc.isObject()) {
        qWarning() << "JSON parse error in" << fileName << ":" << err.errorString();
        return;
    }

    QJsonObject o = doc.object();

    // ---- Put JSON values to the spinboxes / checkboxes ----
    voxel_size_SpinBox->setValue(o.value("voxel_edge_length").toDouble(voxel_size_SpinBox->value()));
    overlapping_factor_SpinBox->setValue(o.value("overlapping_factor").toDouble(overlapping_factor_SpinBox->value()));
    minimum_radius_SpinBox->setValue(o.value("minimum_sphere_radius").toDouble(minimum_radius_SpinBox->value()));
    blood_vessels_icvf_SpinBox->setValue(o.value("blood_vessels_icvf").toDouble(blood_vessels_icvf_SpinBox->value()));

    axons_w_myelin_icvf_SpinBox->setValue(o.value("axons_w_myelin_icvf").toDouble(axons_w_myelin_icvf_SpinBox->value()));
    k1_SpinBox->setValue(o.value("k1").toDouble(k1_SpinBox->value()));
    k2_SpinBox->setValue(o.value("k2").toDouble(k2_SpinBox->value()));
    k3_SpinBox->setValue(o.value("k3").toDouble(k3_SpinBox->value()));
    axons_icvf_SpinBox->setValue(o.value("axons_icvf").toDouble(axons_icvf_SpinBox->value()));
    nbr_threads_SpinBox->setValue(o.value("nbr_threads").toDouble(nbr_threads_SpinBox->value()));
    epsilon_SpinBox->setValue(o.value("epsilon").toDouble(epsilon_SpinBox->value()));
    c2_SpinBox->setValue(o.value("c2").toDouble(c2_SpinBox->value()));
    nbr_axons_populations_SpinBox->setValue(o.value("nbr_axons_populations").toDouble(nbr_axons_populations_SpinBox->value()));
    beading_amplitude_SpinBox->setValue(o.value("beading_amplitude").toDouble(beading_amplitude_SpinBox->value()));
    beading_std_SpinBox->setValue(o.value("beading_std").toDouble(beading_std_SpinBox->value()));
    alpha_SpinBox->setValue(o.value("radius_dist_alpha").toDouble(alpha_SpinBox->value()));
    beta_SpinBox->setValue(o.value("radius_dist_beta").toDouble(beta_SpinBox->value()));

    glial_pop1_soma_icvf_SpinBox->setValue(o.value("glial_pop1_soma_icvf").toDouble(glial_pop1_soma_icvf_SpinBox->value()));
    glial_pop1_processes_icvf_SpinBox->setValue(o.value("glial_pop1_processes_icvf").toDouble(glial_pop1_processes_icvf_SpinBox->value()));
    glial_pop1_radius_mean_SpinBox->setValue(o.value("glial_pop1_soma_radius_mean").toDouble(glial_pop1_radius_mean_SpinBox->value()));
    glial_pop1_radius_std_SpinBox->setValue(o.value("glial_pop1_soma_radius_std").toDouble(glial_pop1_radius_std_SpinBox->value()));
    glial_pop1_mean_process_length_SpinBox->setValue(o.value("glial_pop1_process_length_mean").toDouble(glial_pop1_mean_process_length_SpinBox->value()));
    glial_pop1_std_process_length_SpinBox->setValue(o.value("glial_pop1_process_length_std").toDouble(glial_pop1_std_process_length_SpinBox->value()));
    glial_pop1_nbr_primary_processes_SpinBox->setValue(o.value("glial_pop1_nbr_primary_processes").toDouble(glial_pop1_nbr_primary_processes_SpinBox->value()));
    glial_pop1_branching_checkbox->setChecked(o.value("glial_pop1_branching").toBool(glial_pop1_branching_checkbox->isChecked()));

    
    glial_pop2_soma_icvf_SpinBox->setValue(o.value("glial_pop2_soma_icvf").toDouble(glial_pop2_soma_icvf_SpinBox->value()));
    glial_pop2_processes_icvf_SpinBox->setValue(o.value("glial_pop2_processes_icvf").toDouble(glial_pop2_processes_icvf_SpinBox->value()));
    glial_pop2_radius_mean_SpinBox->setValue(o.value("glial_pop2_soma_radius_mean").toDouble(glial_pop2_radius_mean_SpinBox->value()));
    glial_pop2_radius_std_SpinBox->setValue(o.value("glial_pop2_soma_radius_std").toDouble(glial_pop2_radius_std_SpinBox->value()));
    glial_pop2_mean_process_length_SpinBox->setValue(o.value("glial_pop2_process_length_mean").toDouble(glial_pop2_mean_process_length_SpinBox->value()));
    glial_pop2_std_process_length_SpinBox->setValue(o.value("glial_pop2_process_length_std").toDouble(glial_pop2_std_process_length_SpinBox->value()));
    glial_pop2_nbr_primary_processes_SpinBox->setValue(o.value("glial_pop2_nbr_primary_processes").toDouble(glial_pop2_nbr_primary_processes_SpinBox->value()));
    glial_pop2_branching_checkbox->setChecked(o.value("glial_pop2_branching").toBool(glial_pop2_branching_checkbox->isChecked()));

    nbr_repetitions_SpinBox->setValue(o.value("nbr_repetitions").toDouble(nbr_repetitions_SpinBox->value()));

    selectedDirectory = o.value("directory").toString();
    qInfo() << "Loaded configuration from" << fileName;
}


void Window::loadConfigFromFile()
{

    QString fileName = QFileDialog::getOpenFileName(
        this,
        tr("Open configuration file"),
        QString(),
        tr("Config files (*.json);;All files (*.*)")
    );
    if (fileName.isEmpty())
        return;

    loadConfigFromPath(fileName);
}

void Window::saveConfigToPath(const QString &fileName)
{
   
    QJsonObject o;

    // ---- Map widgets into JSON ----
    o["voxel_edge_length"]            = voxel_size_SpinBox->value();
    o["overlapping_factor"]           = overlapping_factor_SpinBox->value();
    o["minimum_sphere_radius"]        = minimum_radius_SpinBox->value();
    o["blood_vessels_icvf"]           = blood_vessels_icvf_SpinBox->value();

    o["axons_w_myelin_icvf"]          = axons_w_myelin_icvf_SpinBox->value();
    o["k1"]                           = k1_SpinBox->value();
    o["k2"]                           = k2_SpinBox->value();
    o["k3"]                           = k3_SpinBox->value();
    o["axons_icvf"]                   = axons_icvf_SpinBox->value();
    o["nbr_threads"]                  = nbr_threads_SpinBox->value();
    o["epsilon"]                      = epsilon_SpinBox->value();
    o["c2"]                           = c2_SpinBox->value();
    o["nbr_axons_populations"]        = nbr_axons_populations_SpinBox->value();
    o["beading_amplitude"]            = beading_amplitude_SpinBox->value();
    o["beading_std"]                  = beading_std_SpinBox->value();
    o["radius_dist_alpha"]            = alpha_SpinBox->value();
    o["radius_dist_beta"]             = beta_SpinBox->value();

    o["glial_pop1_soma_icvf"]             = glial_pop1_soma_icvf_SpinBox->value();
    o["glial_pop1_processes_icvf"]        = glial_pop1_processes_icvf_SpinBox->value();
    o["glial_pop1_soma_radius_mean"]      = glial_pop1_radius_mean_SpinBox->value();
    o["glial_pop1_soma_radius_std"]       = glial_pop1_radius_std_SpinBox->value();
    o["glial_pop1_process_length_mean"]   = glial_pop1_mean_process_length_SpinBox->value();
    o["glial_pop1_process_length_std"]    = glial_pop1_std_process_length_SpinBox->value();
    o["glial_pop1_nbr_primary_processes"] = glial_pop1_nbr_primary_processes_SpinBox->value();
    o["glial_pop1_branching"]             = glial_pop1_branching_checkbox->isChecked();

    o["glial_pop2_soma_icvf"]             = glial_pop2_soma_icvf_SpinBox->value();
    o["glial_pop2_processes_icvf"]        = glial_pop2_processes_icvf_SpinBox->value();
    o["glial_pop2_soma_radius_mean"]      = glial_pop2_radius_mean_SpinBox->value();
    o["glial_pop2_soma_radius_std"]       = glial_pop2_radius_std_SpinBox->value();
    o["glial_pop2_process_length_mean"]   = glial_pop2_mean_process_length_SpinBox->value();
    o["glial_pop2_process_length_std"]    = glial_pop2_std_process_length_SpinBox->value();
    o["glial_pop2_nbr_primary_processes"] = glial_pop2_nbr_primary_processes_SpinBox->value();
    o["glial_pop2_branching"]             = glial_pop2_branching_checkbox->isChecked();

    o["nbr_repetitions"]                  = nbr_repetitions_SpinBox->value();

    o["directory"] = selectedDirectory;

    QJsonDocument doc(o);
    QFile f(fileName);

    if (!f.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
        qWarning() << "Cannot write config file:" << fileName;
        return;
    }

    f.write(doc.toJson(QJsonDocument::Indented));
    f.close();

    qInfo() << "Auto-saved voxel configuration to:" << fileName;
}


