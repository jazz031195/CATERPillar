#include <iostream>
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "slidergroup.h"
#include "../src/core_logic.h"
#include "../src/Axon.h"
#include "../src/Blood_Vessel.h"
#include "../src/Glial.h"
#include "ScatterDataModifier.h"
#include <QDir>
#include <QFileInfo>
#include <fstream>
#include <QFile>
#include <QTextStream>
#include <QMessageBox>
#include <QDebug>
#include <QFileDialog>
#include <chrono>
#include <variant>
#include <QComboBox>
#include <QFontComboBox>
#include <QThread>
#include <QtCharts/QChartView>
#include <QtCharts/QBarSeries>
#include <QtCharts/QBarSet>
#include <QtCharts/QBarCategoryAxis>
#include <QtCharts/QValueAxis>
#include <QtCharts/QLineSeries>

QT_CHARTS_USE_NAMESPACE

Window::Window(QWidget *parent)
    : QWidget(parent)
{
    // Initialize your windows/processes
    this->openglWindow = new OpenGLWindow(); // Make sure this is instantiated
    this->visualizationWidget = nullptr;
    this->simulatorProcess = new QProcess(this);
    connect(simulatorProcess, static_cast<void(QProcess::*)(int, QProcess::ExitStatus)>(&QProcess::finished),
            this, [this](int exitCode, QProcess::ExitStatus exitStatus) {
        if (exitStatus == QProcess::NormalExit && exitCode == 0) {
            QMessageBox::information(this, "Success", "Simulation finished successfully.");
        } else {
            QMessageBox::warning(this, "Error", "Simulation crashed or exited with an error.");
        }
    });

    // Main layout for the entire window
    QVBoxLayout *mainLayout = new QVBoxLayout(this);
    QTabWidget *mainTabs = new QTabWidget(this);
    mainLayout->addWidget(mainTabs);

    // =========================================================
    // STEP 1: WHITE MATTER CELL GENERATION
    // =========================================================
    QWidget *tabWhiteMatter = new QWidget();
    QVBoxLayout *wmLayout = new QVBoxLayout(tabWhiteMatter);

    // 1A. Clickable Images using QToolButton
    QHBoxLayout *imagesLayout = new QHBoxLayout();
    
    auto createPictureButton = [](const QString& text, const QString& iconPath) {
        QToolButton *btn = new QToolButton();
        btn->setText(text);
        btn->setIcon(QIcon(iconPath)); // Make sure to add these to your Qt Resource file (.qrc)
        btn->setIconSize(QSize(100, 100)); // Adjust size as needed
        btn->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);
        return btn;
    };

    QToolButton *btnGeneral = createPictureButton("General", ":/images/general.png");
    QToolButton *btnAxons = createPictureButton("Axons", ":/images/axons.png");
    QToolButton *btnMyelin = createPictureButton("Myelinated", ":/images/myelin.png");
    QToolButton *btnGlial1 = createPictureButton("Glial Pop 1", ":/images/glial_cells.jpeg");
    QToolButton *btnGlial2 = createPictureButton("Glial Pop 2", ":/images/glial_cells.jpeg");
    QToolButton *btnGlial3 = createPictureButton("Glial Pop 3", ":/images/glial_cells.jpeg");
    QToolButton *btnBlood = createPictureButton("Blood Vessels", ":/images/blood_vessels.png");

    imagesLayout->addWidget(btnGeneral);
    imagesLayout->addWidget(btnAxons);
    imagesLayout->addWidget(btnMyelin);
    imagesLayout->addWidget(btnGlial1);
    imagesLayout->addWidget(btnGlial2);
    imagesLayout->addWidget(btnGlial3);
    imagesLayout->addWidget(btnBlood);

    wmLayout->addLayout(imagesLayout);

    initParameters();
    buildParameterStack(wmLayout);

    // 1C. Connect Buttons to Stacked Widget
    connect(btnGeneral, &QToolButton::clicked, [this]() { cellParamsStack->setCurrentIndex(0); });
    connect(btnAxons,   &QToolButton::clicked, [this]() { cellParamsStack->setCurrentIndex(1); });
    connect(btnMyelin,  &QToolButton::clicked, [this]() { cellParamsStack->setCurrentIndex(2); });
    connect(btnGlial1,  &QToolButton::clicked, [this]() { cellParamsStack->setCurrentIndex(3); });
    connect(btnGlial2,  &QToolButton::clicked, [this]() { cellParamsStack->setCurrentIndex(4); });
    connect(btnGlial3,  &QToolButton::clicked, [this]() { cellParamsStack->setCurrentIndex(5); });
    connect(btnBlood,   &QToolButton::clicked, [this]() { cellParamsStack->setCurrentIndex(6); });

    cellParamsStack->setCurrentIndex(0);

    QPushButton *btnGrow = new QPushButton("Grow White Matter Substrate", this);
    btnGrow->setStyleSheet("font-weight: bold; padding: 10px; margin-top: 10px;");
    wmLayout->addWidget(btnGrow);

    layerProgressBar = new QProgressBar(this);
    layerProgressBar->setRange(0, 100);
    layerProgressBar->setValue(0);
    layerProgressBar->setFormat("Growing axons: 0%");
    layerProgressBar->hide();
    wmLayout->addWidget(layerProgressBar);

    swellingProgressBar = new QProgressBar(this);
    swellingProgressBar->setRange(0, 100);
    swellingProgressBar->setValue(0);
    swellingProgressBar->setFormat("Swelling: not started yet");
    swellingProgressBar->hide();
    wmLayout->addWidget(swellingProgressBar);

    mainTabs->addTab(tabWhiteMatter, "1. White Matter Cell Generation");

    connect(btnGrow, &QPushButton::clicked, this, &Window::onSaveButtonClicked);

    // =========================================================
    // STEP 2: VISUALISATION
    // =========================================================
    QWidget *tabVisualisation = new QWidget();
    QVBoxLayout *visLayout = new QVBoxLayout(tabVisualisation);

    // Crucial: QOpenGLWindow inherits from QWindow, not QWidget. 
    // We must wrap it in a container to embed it inside our QTabWidget.
    QWidget *glContainer = QWidget::createWindowContainer(openglWindow, this);
    glContainer->setMinimumSize(600, 400); // Set an appropriate minimum size
    visLayout->addWidget(glContainer, 1);

    // Create the button and label it for your SWC/CSV functionality
    QPushButton *btnLoadVisFile = new QPushButton("Load SWC/CSV for Visualisation", this);
    btnLoadVisFile->setStyleSheet("padding: 8px; font-weight: bold;");

    QPushButton *btnPlotSholl = new QPushButton("Plot Mean Sholl Curve", this);
    btnPlotSholl->setStyleSheet("padding: 8px; font-weight: bold;");

    QHBoxLayout *visButtonsLayout = new QHBoxLayout();
    visButtonsLayout->addWidget(btnLoadVisFile);
    visButtonsLayout->addWidget(btnPlotSholl);
    visLayout->addLayout(visButtonsLayout);

    // CONNECT DIRECTLY TO YOUR EXISTING SLOT
    connect(btnLoadVisFile, &QPushButton::clicked, this, &Window::SelectSWCFileButton);
    connect(btnPlotSholl, &QPushButton::clicked, this, &Window::ShollAnalysis);

    mainTabs->addTab(tabVisualisation, "2. Visualisation");

    // =========================================================
    // STEP 3: MONTE CARLO SIMULATIONS
    // =========================================================
    QWidget *tabMonteCarlo = new QWidget();
    QVBoxLayout *mcLayout = new QVBoxLayout(tabMonteCarlo);

    // Path to the MC-DC simulator executable (machine-specific, so no default value)
    QFormLayout *executableForm = new QFormLayout();
    inputExecutablePath = new QLineEdit();
    inputExecutablePath->setPlaceholderText("/path/to/MC-DC_Simulator");
    QHBoxLayout *executableLayout = new QHBoxLayout();
    executableLayout->addWidget(inputExecutablePath);
    QPushButton *btnBrowseExecutable = new QPushButton("Browse...");
    executableLayout->addWidget(btnBrowseExecutable);
    executableForm->addRow("Simulator Executable:", executableLayout);
    mcLayout->addLayout(executableForm);

    connect(btnBrowseExecutable, &QPushButton::clicked, [this]() {
        QString file = QFileDialog::getOpenFileName(this, "Select MC-DC Simulator Executable", "", "All Files (*)");
        if (!file.isEmpty()) inputExecutablePath->setText(file);
    });

    // Load an existing config file to run as-is, bypassing the parameter form below
    QFormLayout *loadConfigForm = new QFormLayout();
    inputLoadConfigPath = new QLineEdit();
    inputLoadConfigPath->setPlaceholderText("Leave empty to build a config file from the parameters below");
    QHBoxLayout *loadConfigLayout = new QHBoxLayout();
    loadConfigLayout->addWidget(inputLoadConfigPath);
    QPushButton *btnBrowseLoadConfig = new QPushButton("Browse...");
    loadConfigLayout->addWidget(btnBrowseLoadConfig);
    loadConfigForm->addRow("Load Existing Config File:", loadConfigLayout);
    mcLayout->addLayout(loadConfigForm);

    connect(btnBrowseLoadConfig, &QPushButton::clicked, [this]() {
        QString file = QFileDialog::getOpenFileName(this, "Select Existing Config File", "", "Config Files (*.conf);;All Files (*)");
        if (!file.isEmpty()) inputLoadConfigPath->setText(file);
    });

    QFormLayout *mcForm = new QFormLayout();

    // Default values
    inputN = new QLineEdit("6044");
    inputT = new QLineEdit("55200");
    inputDuration = new QLineEdit("0.092");
    inputDiffIntra = new QLineEdit("2.5e-9");
    inputDiffExtra = new QLineEdit("1.5e-9");
    inputSchemeFile = new QLineEdit();
    inputSchemeFile->setPlaceholderText("/path/to/protocol.scheme (diffusion acquisition scheme)");
    inputCsvPath = new QLineEdit();
    inputCsvPath->setPlaceholderText("/path/to/obstacles.csv (cell geometry file)");

    inputVoxelSizeMC = new QDoubleSpinBox();
    inputVoxelSizeMC->setRange(0.001, 10.0);
    inputVoxelSizeMC->setDecimals(3);
    inputVoxelSizeMC->setSingleStep(0.01);
    inputVoxelSizeMC->setValue(0.1);
    inputVoxelSizeMC->setSuffix(" mm");

    inputNumThreadsMC = new QSpinBox();
    inputNumThreadsMC->setRange(1, 1000);
    int idealThreads = QThread::idealThreadCount();
    inputNumThreadsMC->setValue(idealThreads > 0 ? idealThreads : 1);

    mcForm->addRow("N (Walkers):", inputN);
    mcForm->addRow("T (Time steps):", inputT);
    mcForm->addRow("Duration:", inputDuration);
    mcForm->addRow("Diffusivity Intra:", inputDiffIntra);
    mcForm->addRow("Diffusivity Extra:", inputDiffExtra);

    // Initial walker compartment: both ticked (default) or neither means no
    // restriction is written to the config (voxel-wide sampling); exactly one
    // ticked writes "ini_walkers_pos intra"/"extra" -- see the config-writing block.
    checkIniWalkersIntra = new QCheckBox("Intra");
    checkIniWalkersIntra->setChecked(true);
    checkIniWalkersExtra = new QCheckBox("Extra");
    checkIniWalkersExtra->setChecked(true);

    QHBoxLayout *iniWalkersLayout = new QHBoxLayout();
    iniWalkersLayout->addWidget(checkIniWalkersIntra);
    iniWalkersLayout->addWidget(checkIniWalkersExtra);
    mcForm->addRow("Initial Walker Compartment:", iniWalkersLayout);

    mcForm->addRow("Voxel Size (edge length):", inputVoxelSizeMC);
    mcForm->addRow("Number of Threads:", inputNumThreadsMC);

    QHBoxLayout *SchemeLayout = new QHBoxLayout();
    SchemeLayout->addWidget(inputSchemeFile);
    QPushButton *btnBrowseScheme = new QPushButton("Browse...");
    SchemeLayout->addWidget(btnBrowseScheme);
    mcForm->addRow("Scheme File Path:", SchemeLayout);

    connect(btnBrowseScheme, &QPushButton::clicked, [this]() {
        QString file = QFileDialog::getOpenFileName(this, "Select Scheme File", "", "Scheme Files (*.scheme);;All Files (*)");
        if (!file.isEmpty()) inputSchemeFile->setText(file);
    });

    mcLayout->addLayout(mcForm);

    // Obstacles to include: the user must pick at least one, checked at Run time.
    // A single CSV path is shared across whichever obstacle types are checked.
    QGroupBox *obstacleGroup = new QGroupBox("Obstacles to Include");
    QVBoxLayout *obstacleLayout = new QVBoxLayout();

    checkIncludeAxons = new QCheckBox("Axons");
    checkIncludeAxons->setChecked(true);
    checkIncludeGlial = new QCheckBox("Glial Cells");
    checkIncludeBloodVessels = new QCheckBox("Blood Vessels");

    QLabel *bloodVesselsMCWarningLabel = new QLabel("Blood Vessels: work in progress, not validated yet");
    bloodVesselsMCWarningLabel->setStyleSheet("color: #b36b00;");

    QHBoxLayout *obstacleCheckboxLayout = new QHBoxLayout();
    obstacleCheckboxLayout->addWidget(checkIncludeAxons);
    obstacleCheckboxLayout->addWidget(checkIncludeGlial);
    obstacleCheckboxLayout->addWidget(checkIncludeBloodVessels);
    obstacleCheckboxLayout->addWidget(bloodVesselsMCWarningLabel);
    obstacleLayout->addLayout(obstacleCheckboxLayout);

    QFormLayout *obstacleCsvForm = new QFormLayout();
    QHBoxLayout *csvLayout = new QHBoxLayout();
    csvLayout->addWidget(inputCsvPath);
    QPushButton *btnBrowseCsv = new QPushButton("Browse...");
    csvLayout->addWidget(btnBrowseCsv);
    obstacleCsvForm->addRow("Obstacle CSV Path:", csvLayout);
    obstacleLayout->addLayout(obstacleCsvForm);

    obstacleGroup->setLayout(obstacleLayout);
    mcLayout->addWidget(obstacleGroup);

    connect(btnBrowseCsv, &QPushButton::clicked, [this]() {
        QString file = QFileDialog::getOpenFileName(this, "Select Obstacle CSV", "", "CSV Files (*.csv)");
        if (!file.isEmpty()) inputCsvPath->setText(file);
    });

    QPushButton *btnRun = new QPushButton("Run", this);
    btnRun->setStyleSheet("font-weight: bold; padding: 10px;"); // Make it stand out
    mcLayout->addWidget(btnRun);

    connect(btnRun, &QPushButton::clicked, this, &Window::runMCSimulation);

    mainTabs->addTab(tabMonteCarlo, "3. Monte Carlo Simulation");
}

Window::~Window()
{
    // A run kicked off by StartSimulation() may still be growing on
    // growthThread when the window closes -- std::thread's destructor calls
    // std::terminate() on a still-joinable thread, so wait for it here
    // instead (blocks briefly rather than crashing).
    if (growthThread.joinable()) {
        growthThread.join();
    }
}

void Window::buildParameterStack(QVBoxLayout *wmLayout)
{
    // 1. DEFINE VECTORS
    std::vector<QLabel*> general_labels = { nbr_repetitions_qlabel, voxel_size_qlabel, overlapping_factor_qlabel, minimum_radius_qlabel, nbr_threads_qlabel };
    std::vector<QDoubleSpinBox*> general_spinBoxes = { nbr_repetitions_SpinBox, voxel_size_SpinBox, overlapping_factor_SpinBox, minimum_radius_SpinBox, nbr_threads_SpinBox };
    
    std::vector<QLabel*> axons_labels = { axons_icvf_qlabel, nbr_axons_populations_qlabel, epsilon_qlabel, c2_qlabel, beading_amplitude_qlabel, beading_std_qlabel, alpha_qlabel, beta_qlabel };
    std::vector<QDoubleSpinBox*> axons_spinBoxes = { axons_icvf_SpinBox, nbr_axons_populations_SpinBox, epsilon_SpinBox, c2_SpinBox, beading_amplitude_SpinBox, beading_std_SpinBox, alpha_SpinBox, beta_SpinBox };
    
    std::vector<QLabel*> myelin_labels = { axons_w_myelin_icvf_qlabel, k1_qlabel, k2_qlabel, k3_qlabel };
    std::vector<QDoubleSpinBox*> myelin_spinBoxes = { axons_w_myelin_icvf_SpinBox, k1_SpinBox, k2_SpinBox, k3_SpinBox };

    std::vector<QLabel*> glials_labels1 = { glial_pop1_soma_icvf_qlabel, glial_pop1_processes_icvf_qlabel, glial_pop1_radius_mean_qlabel, glial_pop1_radius_std_qlabel, glial_pop1_mean_process_length_qlabel, glial_pop1_std_process_length_qlabel, glial_pop1_minimum_process_radius_qlabel, glial_pop1_nbr_primary_processes_qlabel };
    std::vector<QDoubleSpinBox*> glials_spinBoxes1 = { glial_pop1_soma_icvf_SpinBox, glial_pop1_processes_icvf_SpinBox, glial_pop1_radius_mean_SpinBox, glial_pop1_radius_std_SpinBox, glial_pop1_mean_process_length_SpinBox, glial_pop1_std_process_length_SpinBox, glial_pop1_minimum_process_radius_SpinBox, glial_pop1_nbr_primary_processes_SpinBox };

    std::vector<QLabel*> glials_labels2 = { glial_pop2_soma_icvf_qlabel, glial_pop2_processes_icvf_qlabel, glial_pop2_radius_mean_qlabel, glial_pop2_radius_std_qlabel, glial_pop2_mean_process_length_qlabel, glial_pop2_std_process_length_qlabel, glial_pop2_minimum_process_radius_qlabel, glial_pop2_nbr_primary_processes_qlabel };
    std::vector<QDoubleSpinBox*> glials_spinBoxes2 = { glial_pop2_soma_icvf_SpinBox, glial_pop2_processes_icvf_SpinBox, glial_pop2_radius_mean_SpinBox, glial_pop2_radius_std_SpinBox, glial_pop2_mean_process_length_SpinBox, glial_pop2_std_process_length_SpinBox, glial_pop2_minimum_process_radius_SpinBox, glial_pop2_nbr_primary_processes_SpinBox };

    std::vector<QLabel*> glials_labels3 = { glial_pop3_soma_icvf_qlabel, glial_pop3_processes_icvf_qlabel, glial_pop3_radius_mean_qlabel, glial_pop3_radius_std_qlabel, glial_pop3_mean_process_length_qlabel, glial_pop3_std_process_length_qlabel, glial_pop3_minimum_process_radius_qlabel, glial_pop3_nbr_primary_processes_qlabel };
    std::vector<QDoubleSpinBox*> glials_spinBoxes3 = { glial_pop3_soma_icvf_SpinBox, glial_pop3_processes_icvf_SpinBox, glial_pop3_radius_mean_SpinBox, glial_pop3_radius_std_SpinBox, glial_pop3_mean_process_length_SpinBox, glial_pop3_std_process_length_SpinBox, glial_pop3_minimum_process_radius_SpinBox, glial_pop3_nbr_primary_processes_SpinBox };

    // 2. INITIALIZE STACK
    cellParamsStack = new QStackedWidget();

    // --- PAGE 0: GENERAL ---
    QGroupBox *generalGroup = new QGroupBox("General Parameters");
    QFormLayout *generalLayout = new QFormLayout; 
    for (size_t i = 0; i < general_labels.size(); i++) {
        generalLayout->addRow(general_labels[i], general_spinBoxes[i]);
    }
    generalLayout->addRow(visualise_voxel_qlabel, visualise_voxel_checkbox);
    generalGroup->setLayout(generalLayout);
    cellParamsStack->addWidget(generalGroup);

    // --- PAGE 1: AXONS ---
    QGroupBox *axonsGroup = new QGroupBox("Axon Parameters");
    QFormLayout *axonsLayout = new QFormLayout;
    for (size_t i = 0; i < axons_labels.size(); i++) {
        if (i == axons_labels.size() - 2) {
            QLabel *gammaLabel = new QLabel("<b>Gamma Distribution parameters for non-myelinated inner radii:</b>");
            axonsLayout->addRow(gammaLabel);
        }
        axonsLayout->addRow(axons_labels[i], axons_spinBoxes[i]);
    }
    axonsGroup->setLayout(axonsLayout);
    cellParamsStack->addWidget(axonsGroup);

    // --- PAGE 2: MYELINATED AXONS ---
    QGroupBox *myelinGroup = new QGroupBox("Myelinated Axon Parameters");
    QGridLayout *myelinLayout = new QGridLayout;
    int mRow = 0;
    QLabel *inheritLabel = new QLabel("<i>Note: Beading is inherited from the Axons tab. Myelinated axons draw "
                                       "their radii from their own Gamma distribution below, independent of the "
                                       "non-myelinated one on the Axons tab.</i>");
    inheritLabel->setWordWrap(true);
    myelinLayout->addWidget(inheritLabel, mRow++, 0, 1, 6);

    QLabel *myelinGammaLabel = new QLabel("<b>Gamma Distribution parameters for myelinated inner radii:</b>");
    myelinLayout->addWidget(myelinGammaLabel, mRow++, 0, 1, 6);
    myelinLayout->addWidget(alpha_myelin_qlabel, mRow, 0);
    myelinLayout->addWidget(alpha_myelin_SpinBox, mRow, 1);
    myelinLayout->addWidget(beta_myelin_qlabel, mRow, 2);
    myelinLayout->addWidget(beta_myelin_SpinBox, mRow++, 3);

    myelinLayout->addWidget(myelin_labels[0], mRow, 0);
    myelinLayout->addWidget(myelin_spinBoxes[0], mRow++, 1);

    myelinLayout->addWidget(myelin_labels[1], mRow, 0); // K1
    myelinLayout->addWidget(myelin_spinBoxes[1], mRow, 1);
    myelinLayout->addWidget(myelin_labels[2], mRow, 2); // K2
    myelinLayout->addWidget(myelin_spinBoxes[2], mRow, 3);
    myelinLayout->addWidget(myelin_labels[3], mRow, 4); // K3
    myelinLayout->addWidget(myelin_spinBoxes[3], mRow++, 5);

    QLabel *formulaLabel = new QLabel("<b>Myelin thickness = K1 + K2 × Inner diameter + K3 × log(Inner diameter)</b>");
    formulaLabel->setAlignment(Qt::AlignCenter);
    myelinLayout->addWidget(formulaLabel, mRow++, 0, 1, 6);
    myelinGroup->setLayout(myelinLayout);
    cellParamsStack->addWidget(myelinGroup);

    // --- PAGE 3: GLIAL POPULATION 1 ---
    QGroupBox *glialGroup1 = new QGroupBox("Glial Cell Population 1 Parameters");
    QFormLayout *glialLayout1 = new QFormLayout;
    for (size_t i = 0; i < glials_labels1.size(); i++) {
        glialLayout1->addRow(glials_labels1[i], glials_spinBoxes1[i]);
    }
    glialLayout1->addRow(glial_pop1_branching_qlabel, glial_pop1_branching_checkbox);
    glialGroup1->setLayout(glialLayout1);
    cellParamsStack->addWidget(glialGroup1);

    // --- PAGE 4: GLIAL POPULATION 2 ---
    QGroupBox *glialGroup2 = new QGroupBox("Glial Cell Population 2 Parameters");
    QFormLayout *glialLayout2 = new QFormLayout;
    for (size_t i = 0; i < glials_labels2.size(); i++) {
        glialLayout2->addRow(glials_labels2[i], glials_spinBoxes2[i]);
    }
    glialLayout2->addRow(glial_pop2_branching_qlabel, glial_pop2_branching_checkbox);
    glialGroup2->setLayout(glialLayout2);
    cellParamsStack->addWidget(glialGroup2);

    // --- PAGE 5: GLIAL POPULATION 3 ---
    QGroupBox *glialGroup3 = new QGroupBox("Glial Cell Population 3 Parameters");
    QFormLayout *glialLayout3 = new QFormLayout;
    for (size_t i = 0; i < glials_labels3.size(); i++) {
        glialLayout3->addRow(glials_labels3[i], glials_spinBoxes3[i]);
    }
    glialLayout3->addRow(glial_pop3_branching_qlabel, glial_pop3_branching_checkbox);
    glialGroup3->setLayout(glialLayout3);
    cellParamsStack->addWidget(glialGroup3);

    // --- PAGE 6: BLOOD VESSELS ---
    QGroupBox *bloodVesselGroup = new QGroupBox("Blood Vessel Parameters");
    QFormLayout *bloodLayout = new QFormLayout;
    QLabel *bloodVesselsWarningLabel = new QLabel("<b>Blood Vessels: work in progress, not validated yet.</b>");
    bloodVesselsWarningLabel->setStyleSheet("color: #b36b00;");
    bloodLayout->addRow(bloodVesselsWarningLabel);
    bloodLayout->addRow(blood_vessels_icvf_qlabel, blood_vessels_icvf_SpinBox);
    bloodLayout->addRow(blood_vessels_processes_icvf_qlabel, blood_vessels_processes_icvf_SpinBox);
    bloodLayout->addRow(blood_vessel_voxel_size_qlabel, blood_vessel_voxel_size_SpinBox);
    bloodLayout->addRow(blood_vessel_mean_radius_qlabel, blood_vessel_mean_radius_SpinBox);
    bloodLayout->addRow(blood_vessel_std_radius_qlabel, blood_vessel_std_radius_SpinBox);
    bloodLayout->addRow(blood_vessel_gamma_qlabel, blood_vessel_gamma_SpinBox);
    bloodLayout->addRow(blood_vessel_max_generations_qlabel, blood_vessel_max_generations_SpinBox);
    bloodVesselGroup->setLayout(bloodLayout);
    cellParamsStack->addWidget(bloodVesselGroup);

    // 3. ADD TO THE MAIN TAB LAYOUT
    wmLayout->addWidget(cellParamsStack);
}

void Window::runMCSimulation()
{
    if (simulatorProcess->state() != QProcess::NotRunning) {
        QMessageBox::warning(this, "Error", "A Monte Carlo simulation is already running.");
        return;
    }

    QString confFilePath = inputLoadConfigPath->text().trimmed();

    if (!confFilePath.isEmpty()) {
        // An existing config file was supplied: read it directly, skipping the form below entirely.
        if (!QFile::exists(confFilePath)) {
            QMessageBox::critical(this, "Error", "Could not find the config file at:\n" + confFilePath);
            return;
        }
    } else {
        // No existing config supplied: the user must pick at least one obstacle type to build the config from.
        if (!checkIncludeAxons->isChecked() && !checkIncludeGlial->isChecked() && !checkIncludeBloodVessels->isChecked()) {
            QMessageBox::critical(this, "Error", "Please select at least one obstacle type to include (Axons, Glial Cells, or Blood Vessels).");
            return;
        }
        if (inputCsvPath->text().trimmed().isEmpty()) {
            QMessageBox::critical(this, "Error", "Please specify the Obstacle CSV path.");
            return;
        }

        // No existing config supplied: build one from the form, into a user-chosen output directory.
        QString outputDir = QFileDialog::getExistingDirectory(this,
                                                               tr("Select Output Directory for Monte Carlo Simulation"),
                                                               QDir::homePath());
        if (outputDir.isEmpty()) {
            return; // User canceled the dialog
        }

        confFilePath = outputDir + "/config.conf";
        if (QFile::exists(confFilePath)) {
            int rep = 0;
            QString candidate;
            do {
                candidate = outputDir + QString("/config_rep_%1.conf").arg(rep, 2, 10, QChar('0'));
                rep++;
            } while (QFile::exists(candidate));
            confFilePath = candidate;
        }
        QString expPrefix = outputDir + "/run";

        QFile file(confFilePath);
        if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
            QMessageBox::critical(this, "Error", "Could not create configuration file.");
            return;
        }

        QTextStream out(&file);
        out << "N " << inputN->text() << "\n";
        out << "T " << inputT->text() << "\n";
        out << "duration " << inputDuration->text() << "\n";
        out << "diffusivity_intra " << inputDiffIntra->text() << "\n";
        out << "diffusivity_extra " << inputDiffExtra->text() << "\n";

        // Initial walker compartment: only write ini_walkers_pos when exactly one of
        // Intra/Extra is ticked. Both ticked (or neither) means no restriction, so the
        // simulator's default voxel-wide sampling is used and nothing is written.
        bool iniIntra = checkIniWalkersIntra->isChecked();
        bool iniExtra = checkIniWalkersExtra->isChecked();
        if (iniIntra && !iniExtra) {
            out << "ini_walkers_pos intra\n";
        } else if (iniExtra && !iniIntra) {
            out << "ini_walkers_pos extra\n";
        }

        out << "num_process " << inputNumThreadsMC->value() << "\n\n";

        out << "scheme_file " << inputSchemeFile->text() << "\n\n";

        // exp_prefix tells the simulator where to save the config's output files (traj, hit, signals, ...)
        out << "exp_prefix " << expPrefix << "\n\n";

        // Inject the same CSV path into a list block for each checked obstacle type
        out << "<obstacle>\n";
        if (checkIncludeAxons->isChecked()) {
            out << "<axons_list>\n";
            out << inputCsvPath->text() << "\n";
            out << "permeability global 0\n";
            out << "</axons_list>\n";
        }
        if (checkIncludeGlial->isChecked()) {
            out << "<glials_list>\n";
            out << inputCsvPath->text() << "\n";
            out << "permeability global 0\n";
            out << "</glials_list>\n";
        }
        if (checkIncludeBloodVessels->isChecked()) {
            out << "<blood_vessels_list>\n";
            out << inputCsvPath->text() << "\n";
            out << "permeability global 0\n";
            out << "</blood_vessels_list>\n";
        }
        out << "</obstacle>\n\n";

        // Add the voxel block, sized from the user-defined voxel edge length.
        double voxelSize = inputVoxelSizeMC->value();
        QString voxelSizeStr = QString::number(voxelSize);
        out << "<voxel>\n0.0 0.0 0.0\n" << voxelSizeStr << " " << voxelSizeStr << " " << voxelSizeStr << "\n</voxel>\n\n";
        out << "<END>\n";
        file.close();
    }

    // Execute the simulation
    QString executablePath = inputExecutablePath->text();
    if (executablePath.isEmpty()) {
        QMessageBox::critical(this, "Error", "Please specify the path to the MC-DC simulator executable.");
        return;
    }
    if (!QFile::exists(executablePath)) {
        QMessageBox::critical(this, "Error", "Could not find the Monte Carlo simulator executable at:\n" + executablePath);
        return;
    }

    QStringList arguments;
    arguments << confFilePath;

    simulatorProcess->start(executablePath, arguments);

    if (!simulatorProcess->waitForStarted()) {
        QMessageBox::critical(this, "Error", "Could not start the simulator executable at:\n" + executablePath);
    }
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

        // The loaded file has no header describing the real voxel's actual
        // placement (it can be anywhere inside the padded blood-vessel voxel,
        // see PlaceSmallVoxel) -- try the accompanying growth_info.txt this
        // simulation would have written alongside it (same basename, per
        // CoreLogic::runSimulation's own naming convention) for the "Small
        // voxel min/max limits" lines it records; fall back to an
        // origin-anchored box of edge parameters.voxel_size, the only
        // reasonable guess with no growth_info.txt to read.
        voxelBoundsMin = QVector3D(0.0f, 0.0f, 0.0f);
        voxelBoundsMax = QVector3D(parameters.voxel_size, parameters.voxel_size, parameters.voxel_size);
        QFileInfo fi(filePath);
        QString growthInfoPath = fi.absolutePath() + "/" + fi.completeBaseName() + "_growth_info.txt";
        QFile growthInfoFile(growthInfoPath);
        if (growthInfoFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
            QTextStream in(&growthInfoFile);
            bool haveMin = false, haveMax = false;
            while (!in.atEnd()) {
                QString line = in.readLine();
                if (line.startsWith("Small voxel min limits ")) {
                    QStringList parts = line.split(' ', Qt::SkipEmptyParts);
                    if (parts.size() >= 6) {
                        voxelBoundsMin = QVector3D(parts[3].toFloat(), parts[4].toFloat(), parts[5].toFloat());
                        haveMin = true;
                    }
                } else if (line.startsWith("Small voxel max limits ")) {
                    QStringList parts = line.split(' ', Qt::SkipEmptyParts);
                    if (parts.size() >= 6) {
                        voxelBoundsMax = QVector3D(parts[3].toFloat(), parts[4].toFloat(), parts[5].toFloat());
                        haveMax = true;
                    }
                }
                if (haveMin && haveMax) break;
            }
            growthInfoFile.close();
        }

        PlotCells(true, true, true, true, true);              // Visualize
    }
}
void Window::PlotCells(const bool& axons_plot,
                       const bool& glial_pop1_plot,
                       const bool& glial_pop2_plot,
                       const bool& glial_pop3_plot,
                       const bool& blood_vessels_plot)
{
    // 1. DATA PREPARATION
    std::vector<std::vector<double>> X, Y, Z, R;
    std::vector<int> groupIds;  

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
    if (glial_pop3_plot){
        X.insert(X.end(), X_glial_pop3.begin(), X_glial_pop3.end());
        Y.insert(Y.end(), Y_glial_pop3.begin(), Y_glial_pop3.end());
        Z.insert(Z.end(), Z_glial_pop3.begin(), Z_glial_pop3.end());
        R.insert(R.end(), R_glial_pop3.begin(), R_glial_pop3.end());
        groupIds.insert(groupIds.end(), X_glial_pop3.size(), static_cast<int>(OpenGLWindow::SphereGroup::Glial3));
    }
    if (blood_vessels_plot){
        X.insert(X.end(), X_blood_vessels.begin(), X_blood_vessels.end());
        Y.insert(Y.end(), Y_blood_vessels.begin(), Y_blood_vessels.end());
        Z.insert(Z.end(), Z_blood_vessels.begin(), Z_blood_vessels.end());
        R.insert(R.end(), R_blood_vessels.begin(), R_blood_vessels.end());
        groupIds.insert(groupIds.end(), X_blood_vessels.size(), static_cast<int>(OpenGLWindow::SphereGroup::Blood));
    }

   // 2. WINDOW CREATION (Only if it doesn't exist)
    if (this->openglWindow == nullptr) {
        this->openglWindow = new OpenGLWindow();

        QSurfaceFormat format;
        format.setDepthBufferSize(24);
        this->openglWindow->setFormat(format);
        
        this->openglWindow->setTitle("3D Spheres Visualization");
        this->openglWindow->resize(800, 600);

        visualizationWidget = new QWidget; 
        QHBoxLayout *hLayout = new QHBoxLayout(visualizationWidget);
        QVBoxLayout *vLayout = new QVBoxLayout();
        
        QWidget *container = QWidget::createWindowContainer(this->openglWindow);
        hLayout->addWidget(container, 1);
        hLayout->addLayout(vLayout);

        // ---> CRITICAL RESTORE: Add your modifier back here <---
        ScatterDataModifier *modifier = new ScatterDataModifier(this->openglWindow);

        // Buttons
        QPushButton *plotRadiiButton = new QPushButton("Plot Radii Distribution", visualizationWidget);
        QPushButton *plotTortuosityButton = new QPushButton("Plot Tortuosity Distribution", visualizationWidget);
        QPushButton *plotShollAnalysisButton = new QPushButton("Plot Sholl Analysis", visualizationWidget);
        QPushButton *resetCameraButton = new QPushButton("Reset Camera", visualizationWidget);
        QPushButton *hideAxonsButton = new QPushButton("Hide Axons", visualizationWidget);
        QPushButton *hideGlialCellsButton = new QPushButton("Hide Glial Cells", visualizationWidget);
        QPushButton *showall = new QPushButton("Show All Cells", visualizationWidget);

        vLayout->addWidget(plotRadiiButton, 0, Qt::AlignTop);
        vLayout->addWidget(plotTortuosityButton, 0, Qt::AlignTop);
        vLayout->addWidget(plotShollAnalysisButton, 0, Qt::AlignTop);
        vLayout->addWidget(resetCameraButton, 0, Qt::AlignTop);
        vLayout->addWidget(hideAxonsButton, 0, Qt::AlignTop);
        vLayout->addWidget(hideGlialCellsButton, 0, Qt::AlignTop);
        vLayout->addWidget(showall, 0, Qt::AlignTop);

        QObject::connect(plotRadiiButton, &QPushButton::clicked, this, &Window::plotRadiusDistribution);
        QObject::connect(plotTortuosityButton, &QPushButton::clicked, this, &Window::plotTortuosityDistribution);
        QObject::connect(plotShollAnalysisButton, &QPushButton::clicked, this, &Window::ShollAnalysis);
        QObject::connect(resetCameraButton, &QPushButton::clicked, this->openglWindow, &OpenGLWindow::resetCamera);
        QObject::connect(hideAxonsButton, &QPushButton::clicked, this, &Window::HideAxons);
        QObject::connect(hideGlialCellsButton, &QPushButton::clicked, this, &Window::HideGlialCells);
        QObject::connect(showall, &QPushButton::clicked, this, &Window::ShowAllCells);

        visualizationWidget->show();
    }

    if (X.empty()) {
        QMessageBox::warning(this, "Empty Simulation", "No cells were generated in the selected voxel. Adjust parameters and try again.");
        return; // Stop here so we don't crash the renderer
    }

    // 3. DATA UPDATE
    if (this->openglWindow) {
        this->openglWindow->setSpheres(X, Y, Z, R, groupIds);
        this->openglWindow->setVoxelBounds(voxelBoundsMin, voxelBoundsMax);
        this->openglWindow->update();
    }

    if (this->visualizationWidget) {
        this->visualizationWidget->show();
        this->visualizationWidget->raise();
        this->visualizationWidget->activateWindow();
    }
}
void Window::resetCamera(){
    if (openglWindow) {
        openglWindow->resetCamera();
    }
}

void Window::ShowAllCells(){

    if (!openglWindow) return;

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

    X.insert(X.end(), X_glial_pop3.begin(), X_glial_pop3.end());
    Y.insert(Y.end(), Y_glial_pop3.begin(), Y_glial_pop3.end());
    Z.insert(Z.end(), Z_glial_pop3.begin(), Z_glial_pop3.end());
    R.insert(R.end(), R_glial_pop3.begin(), R_glial_pop3.end());
    groupIds.insert(groupIds.end(), X_glial_pop3.size(), static_cast<int>(OpenGLWindow::SphereGroup::Glial3));

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

    X.insert(X.end(), X_glial_pop3.begin(), X_glial_pop3.end());
    Y.insert(Y.end(), Y_glial_pop3.begin(), Y_glial_pop3.end());
    Z.insert(Z.end(), Z_glial_pop3.begin(), Z_glial_pop3.end());
    R.insert(R.end(), R_glial_pop3.begin(), R_glial_pop3.end());
    groupIds.insert(groupIds.end(), X_glial_pop3.size(), static_cast<int>(OpenGLWindow::SphereGroup::Glial3));

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

void Window::initParameters()
{
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
    glial_pop3_soma_icvf_qlabel = new QLabel(tr("Glial Cell somas ICVF (%):"));
    glial_pop2_processes_icvf_qlabel = new QLabel(tr("Glial Cell processes ICVF (%):"));
    glial_pop3_processes_icvf_qlabel = new QLabel(tr("Glial Cell processes ICVF (%):"));
    blood_vessels_icvf_qlabel = new QLabel(tr("Arteriole ICVF (%):"));
    blood_vessels_processes_icvf_qlabel = new QLabel(tr("Capillaries ICVF (%):"));
    blood_vessel_voxel_size_qlabel = new QLabel(tr("Blood Vessel Voxel Edge Length (μm):"));
    blood_vessel_mean_radius_qlabel = new QLabel(tr("Arteriole Radius Mean (μm):"));
    blood_vessel_std_radius_qlabel = new QLabel(tr("Arteriole Radius Standard Deviation (μm):"));
    blood_vessel_gamma_qlabel = new QLabel(tr("Capillary Branching Exponent (γ):"));
    blood_vessel_max_generations_qlabel = new QLabel(tr("Max Capillary Generations:"));
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
    glial_pop3_mean_process_length_qlabel = new QLabel(tr("Mean Process Length (μm):"));
    glial_pop2_std_process_length_qlabel = new QLabel(tr("Standard Deviation Process Length (μm):"));
    glial_pop3_std_process_length_qlabel = new QLabel(tr("Standard Deviation Process Length (μm):"));
    alpha_qlabel = new QLabel(tr("α:"));
    beta_qlabel = new QLabel(tr("β:"));
    alpha_myelin_qlabel = new QLabel(tr("α (myelinated):"));
    beta_myelin_qlabel = new QLabel(tr("β (myelinated):"));
    glial_pop1_radius_mean_qlabel = new QLabel(tr("Glial Cell Soma Radius Mean:"));
    glial_pop1_radius_std_qlabel = new QLabel(tr("Glial Cell Soma Radius Standard Deviation:"));
    glial_pop2_radius_mean_qlabel = new QLabel(tr("Glial Cell Soma Radius Mean:"));
    glial_pop3_radius_mean_qlabel = new QLabel(tr("Glial Cell Soma Radius Mean:"));
    glial_pop2_radius_std_qlabel = new QLabel(tr("Glial Cell Soma Radius Standard Deviation:"));
    glial_pop3_radius_std_qlabel = new QLabel(tr("Glial Cell Soma Radius Standard Deviation:"));
    glial_pop1_minimum_process_radius_qlabel = new QLabel(tr("Minimum Process Radius (μm):"));
    glial_pop2_minimum_process_radius_qlabel = new QLabel(tr("Minimum Process Radius (μm):"));
    glial_pop3_minimum_process_radius_qlabel = new QLabel(tr("Minimum Process Radius (μm):"));
    glial_pop1_nbr_primary_processes_qlabel = new QLabel(tr("Number of Primary Processes:"));
    glial_pop2_nbr_primary_processes_qlabel = new QLabel(tr("Number of Primary Processes:"));
    glial_pop3_nbr_primary_processes_qlabel = new QLabel(tr("Number of Primary Processes:"));
    glial_pop1_branching_qlabel = new QLabel(tr("Can Glial Cell Population have branching ? "));
    glial_pop2_branching_qlabel = new QLabel(tr("Can Glial Cell Population have branching ? "));
    glial_pop3_branching_qlabel = new QLabel(tr("Can Glial Cell Population have branching ? "));

    nbr_repetitions_SpinBox = new QDoubleSpinBox;
    nbr_repetitions_SpinBox->setRange(1, 100);
    nbr_repetitions_SpinBox->setSingleStep(1);
    nbr_repetitions_SpinBox->setValue(1);

    visualise_voxel_checkbox = new QCheckBox;
    visualise_voxel_checkbox->setChecked(true);

    glial_pop1_branching_checkbox = new QCheckBox;
    glial_pop1_branching_checkbox->setChecked(true);

    glial_pop2_branching_checkbox = new QCheckBox;
    glial_pop3_branching_checkbox = new QCheckBox;
    glial_pop2_branching_checkbox->setChecked(true);
    glial_pop3_branching_checkbox->setChecked(true);

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

    alpha_myelin_SpinBox = new QDoubleSpinBox;
    alpha_myelin_SpinBox->setRange(0, 10);
    alpha_myelin_SpinBox->setSingleStep(0.1);
    alpha_myelin_SpinBox->setValue(2);

    beta_myelin_SpinBox = new QDoubleSpinBox;
    beta_myelin_SpinBox->setRange(0, 10);
    beta_myelin_SpinBox->setSingleStep(0.001);
    beta_myelin_SpinBox->setValue(0.25);

    epsilon_SpinBox = new QDoubleSpinBox;
    epsilon_SpinBox->setRange(0, 2);
    epsilon_SpinBox->setSingleStep(0.1);
    epsilon_SpinBox->setValue(0.4);

    glial_pop1_mean_process_length_SpinBox = new QDoubleSpinBox;
    glial_pop1_mean_process_length_SpinBox->setRange(0, 100);
    glial_pop1_mean_process_length_SpinBox->setSingleStep(1);
    glial_pop1_mean_process_length_SpinBox->setValue(30);

    glial_pop1_std_process_length_SpinBox = new QDoubleSpinBox;
    glial_pop1_std_process_length_SpinBox->setRange(0, 100);
    glial_pop1_std_process_length_SpinBox->setSingleStep(1);
    glial_pop1_std_process_length_SpinBox->setValue(15);

    glial_pop2_mean_process_length_SpinBox = new QDoubleSpinBox;
    glial_pop3_mean_process_length_SpinBox = new QDoubleSpinBox;
    glial_pop2_mean_process_length_SpinBox->setRange(0, 100);
    glial_pop3_mean_process_length_SpinBox->setRange(0, 100);
    glial_pop2_mean_process_length_SpinBox->setSingleStep(1);
    glial_pop3_mean_process_length_SpinBox->setSingleStep(1);
    glial_pop2_mean_process_length_SpinBox->setValue(30);
    glial_pop3_mean_process_length_SpinBox->setValue(30);

    glial_pop2_std_process_length_SpinBox = new QDoubleSpinBox;
    glial_pop3_std_process_length_SpinBox = new QDoubleSpinBox;
    glial_pop2_std_process_length_SpinBox->setRange(0, 100);
    glial_pop3_std_process_length_SpinBox->setRange(0, 100);
    glial_pop2_std_process_length_SpinBox->setSingleStep(1);
    glial_pop3_std_process_length_SpinBox->setSingleStep(1);
    glial_pop2_std_process_length_SpinBox->setValue(15);
    glial_pop3_std_process_length_SpinBox->setValue(15);


    glial_pop1_nbr_primary_processes_SpinBox = new QDoubleSpinBox;
    glial_pop1_nbr_primary_processes_SpinBox->setRange(1, 20);
    glial_pop1_nbr_primary_processes_SpinBox->setSingleStep(1);
    glial_pop1_nbr_primary_processes_SpinBox->setValue(10);

    glial_pop2_nbr_primary_processes_SpinBox = new QDoubleSpinBox;
    glial_pop3_nbr_primary_processes_SpinBox = new QDoubleSpinBox;
    glial_pop2_nbr_primary_processes_SpinBox->setRange(1, 20);
    glial_pop3_nbr_primary_processes_SpinBox->setRange(1, 20);
    glial_pop2_nbr_primary_processes_SpinBox->setSingleStep(1);
    glial_pop3_nbr_primary_processes_SpinBox->setSingleStep(1);
    glial_pop2_nbr_primary_processes_SpinBox->setValue(10);
    glial_pop3_nbr_primary_processes_SpinBox->setValue(10);

    axons_icvf_SpinBox = new QDoubleSpinBox;
    axons_icvf_SpinBox->setRange(0, 100);
    axons_icvf_SpinBox->setSingleStep(1);

    axons_w_myelin_icvf_SpinBox = new QDoubleSpinBox;
    axons_w_myelin_icvf_SpinBox->setRange(0, 100);
    axons_w_myelin_icvf_SpinBox->setSingleStep(1);

    blood_vessels_icvf_SpinBox = new QDoubleSpinBox;
    blood_vessels_icvf_SpinBox->setRange(0, 100);
    blood_vessels_icvf_SpinBox->setSingleStep(1);

    blood_vessels_processes_icvf_SpinBox = new QDoubleSpinBox;
    blood_vessels_processes_icvf_SpinBox->setRange(0, 100);
    blood_vessels_processes_icvf_SpinBox->setSingleStep(1);
    blood_vessels_processes_icvf_SpinBox->setValue(0);

    // 0 means "no padding, same box as Voxel Edge Length"; matches
    // voxel_size_SpinBox's own range/step/default so leaving it untouched
    // reproduces today's behavior exactly.
    blood_vessel_voxel_size_SpinBox = new QDoubleSpinBox;
    blood_vessel_voxel_size_SpinBox->setRange(0, 1000);
    blood_vessel_voxel_size_SpinBox->setSingleStep(1);
    blood_vessel_voxel_size_SpinBox->setValue(30);

    blood_vessel_mean_radius_SpinBox = new QDoubleSpinBox;
    blood_vessel_mean_radius_SpinBox->setRange(0, 20);
    blood_vessel_mean_radius_SpinBox->setSingleStep(0.1);
    blood_vessel_mean_radius_SpinBox->setValue(6);

    blood_vessel_std_radius_SpinBox = new QDoubleSpinBox;
    blood_vessel_std_radius_SpinBox->setRange(0, 20);
    blood_vessel_std_radius_SpinBox->setSingleStep(0.1);
    blood_vessel_std_radius_SpinBox->setValue(1);

    blood_vessel_gamma_SpinBox = new QDoubleSpinBox;
    blood_vessel_gamma_SpinBox->setRange(0.5, 10);
    blood_vessel_gamma_SpinBox->setSingleStep(0.1);
    blood_vessel_gamma_SpinBox->setValue(3);

    blood_vessel_max_generations_SpinBox = new QDoubleSpinBox;
    blood_vessel_max_generations_SpinBox->setRange(1, 20);
    blood_vessel_max_generations_SpinBox->setSingleStep(1);
    blood_vessel_max_generations_SpinBox->setDecimals(0);
    blood_vessel_max_generations_SpinBox->setValue(7);

    glial_pop1_soma_icvf_SpinBox = new QDoubleSpinBox;
    glial_pop1_soma_icvf_SpinBox->setRange(0, 100);
    glial_pop1_soma_icvf_SpinBox->setSingleStep(1);

    glial_pop1_processes_icvf_SpinBox = new QDoubleSpinBox;
    glial_pop1_processes_icvf_SpinBox->setRange(0, 100);
    glial_pop1_processes_icvf_SpinBox->setSingleStep(1);

    glial_pop2_soma_icvf_SpinBox = new QDoubleSpinBox;
    glial_pop3_soma_icvf_SpinBox = new QDoubleSpinBox;
    glial_pop2_soma_icvf_SpinBox->setRange(0, 100);
    glial_pop3_soma_icvf_SpinBox->setRange(0, 100);
    glial_pop2_soma_icvf_SpinBox->setSingleStep(1);
    glial_pop3_soma_icvf_SpinBox->setSingleStep(1);

    glial_pop2_processes_icvf_SpinBox = new QDoubleSpinBox;
    glial_pop3_processes_icvf_SpinBox = new QDoubleSpinBox;
    glial_pop2_processes_icvf_SpinBox->setRange(0, 100);
    glial_pop3_processes_icvf_SpinBox->setRange(0, 100);
    glial_pop2_processes_icvf_SpinBox->setSingleStep(1);
    glial_pop3_processes_icvf_SpinBox->setSingleStep(1);

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
    glial_pop1_radius_mean_SpinBox->setValue(3);

    glial_pop1_radius_std_SpinBox = new QDoubleSpinBox;
    glial_pop1_radius_std_SpinBox->setRange(0, 10);
    glial_pop1_radius_std_SpinBox->setSingleStep(0.1);
    glial_pop1_radius_std_SpinBox->setValue(0.5);

    glial_pop2_radius_mean_SpinBox = new QDoubleSpinBox;
    glial_pop3_radius_mean_SpinBox = new QDoubleSpinBox;
    glial_pop2_radius_mean_SpinBox->setRange(0, 10);
    glial_pop3_radius_mean_SpinBox->setRange(0, 10);
    glial_pop2_radius_mean_SpinBox->setSingleStep(0.1);
    glial_pop3_radius_mean_SpinBox->setSingleStep(0.1);
    glial_pop2_radius_mean_SpinBox->setValue(3);
    glial_pop3_radius_mean_SpinBox->setValue(3);

    glial_pop2_radius_std_SpinBox = new QDoubleSpinBox;
    glial_pop3_radius_std_SpinBox = new QDoubleSpinBox;
    glial_pop2_radius_std_SpinBox->setRange(0, 10);
    glial_pop3_radius_std_SpinBox->setRange(0, 10);
    glial_pop2_radius_std_SpinBox->setSingleStep(0.1);
    glial_pop3_radius_std_SpinBox->setSingleStep(0.1);
    glial_pop2_radius_std_SpinBox->setValue(0.5);
    glial_pop3_radius_std_SpinBox->setValue(0.5);

    glial_pop1_minimum_process_radius_SpinBox = new QDoubleSpinBox;
    glial_pop1_minimum_process_radius_SpinBox->setRange(0, 10);
    glial_pop1_minimum_process_radius_SpinBox->setSingleStep(0.01);
    glial_pop1_minimum_process_radius_SpinBox->setValue(0.15);

    glial_pop2_minimum_process_radius_SpinBox = new QDoubleSpinBox;
    glial_pop3_minimum_process_radius_SpinBox = new QDoubleSpinBox;
    glial_pop2_minimum_process_radius_SpinBox->setRange(0, 10);
    glial_pop3_minimum_process_radius_SpinBox->setRange(0, 10);
    glial_pop2_minimum_process_radius_SpinBox->setSingleStep(0.01);
    glial_pop3_minimum_process_radius_SpinBox->setSingleStep(0.01);
    glial_pop2_minimum_process_radius_SpinBox->setValue(0.15);
    glial_pop3_minimum_process_radius_SpinBox->setValue(0.15);

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
    blood_vessels_icvf_qlabel = new QLabel(tr("Arteriole ICVF (%):"));
    blood_vessels_processes_icvf_qlabel = new QLabel(tr("Capillaries ICVF (%):"));
    blood_vessel_voxel_size_qlabel = new QLabel(tr("Blood Vessel Voxel Edge Length (μm):"));
    blood_vessel_mean_radius_qlabel = new QLabel(tr("Arteriole Radius Mean (μm):"));
    blood_vessel_std_radius_qlabel = new QLabel(tr("Arteriole Radius Standard Deviation (μm):"));
    blood_vessel_gamma_qlabel = new QLabel(tr("Capillary Branching Exponent (γ):"));
    blood_vessel_max_generations_qlabel = new QLabel(tr("Max Capillary Generations:"));
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
    alpha_myelin_qlabel = new QLabel(tr("α (myelinated):"));
    beta_myelin_qlabel = new QLabel(tr("β (myelinated):"));
    glial_pop1_radius_mean_qlabel = new QLabel(tr("Glial Cell Soma Radius Mean:"));
    glial_pop1_radius_std_qlabel = new QLabel(tr("Glial Cell Soma Radius Standard Deviation:"));
    glial_pop2_radius_mean_qlabel = new QLabel(tr("Glial Cell Soma Radius Mean:"));
    glial_pop2_radius_std_qlabel = new QLabel(tr("Glial Cell Soma Radius Standard Deviation:"));
    glial_pop1_minimum_process_radius_qlabel = new QLabel(tr("Minimum Process Radius (μm):"));
    glial_pop2_minimum_process_radius_qlabel = new QLabel(tr("Minimum Process Radius (μm):"));
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

    alpha_myelin_SpinBox = new QDoubleSpinBox;
    alpha_myelin_SpinBox->setRange(0, 10);
    alpha_myelin_SpinBox->setSingleStep(0.1);
    alpha_myelin_SpinBox->setValue(2);

    beta_myelin_SpinBox = new QDoubleSpinBox;
    beta_myelin_SpinBox->setRange(0, 10);
    beta_myelin_SpinBox->setSingleStep(0.001);
    beta_myelin_SpinBox->setValue(0.25);

    epsilon_SpinBox = new QDoubleSpinBox;
    epsilon_SpinBox->setRange(0, 2);
    epsilon_SpinBox->setSingleStep(0.1);
    epsilon_SpinBox->setValue(0.4);

    glial_pop1_mean_process_length_SpinBox = new QDoubleSpinBox;
    glial_pop1_mean_process_length_SpinBox->setRange(0, 100);
    glial_pop1_mean_process_length_SpinBox->setSingleStep(1);
    glial_pop1_mean_process_length_SpinBox->setValue(30);

    glial_pop1_std_process_length_SpinBox = new QDoubleSpinBox;
    glial_pop1_std_process_length_SpinBox->setRange(0, 100);
    glial_pop1_std_process_length_SpinBox->setSingleStep(1);
    glial_pop1_std_process_length_SpinBox->setValue(15);

    glial_pop2_mean_process_length_SpinBox = new QDoubleSpinBox;
    glial_pop2_mean_process_length_SpinBox->setRange(0, 100);
    glial_pop2_mean_process_length_SpinBox->setSingleStep(1);
    glial_pop2_mean_process_length_SpinBox->setValue(30);

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

    blood_vessels_processes_icvf_SpinBox = new QDoubleSpinBox;
    blood_vessels_processes_icvf_SpinBox->setRange(0, 100);
    blood_vessels_processes_icvf_SpinBox->setSingleStep(1);
    blood_vessels_processes_icvf_SpinBox->setValue(0);

    // 0 means "no padding, same box as Voxel Edge Length"; matches
    // voxel_size_SpinBox's own range/step/default so leaving it untouched
    // reproduces today's behavior exactly.
    blood_vessel_voxel_size_SpinBox = new QDoubleSpinBox;
    blood_vessel_voxel_size_SpinBox->setRange(0, 1000);
    blood_vessel_voxel_size_SpinBox->setSingleStep(1);
    blood_vessel_voxel_size_SpinBox->setValue(30);

    blood_vessel_mean_radius_SpinBox = new QDoubleSpinBox;
    blood_vessel_mean_radius_SpinBox->setRange(0, 20);
    blood_vessel_mean_radius_SpinBox->setSingleStep(0.1);
    blood_vessel_mean_radius_SpinBox->setValue(6);

    blood_vessel_std_radius_SpinBox = new QDoubleSpinBox;
    blood_vessel_std_radius_SpinBox->setRange(0, 20);
    blood_vessel_std_radius_SpinBox->setSingleStep(0.1);
    blood_vessel_std_radius_SpinBox->setValue(1);

    blood_vessel_gamma_SpinBox = new QDoubleSpinBox;
    blood_vessel_gamma_SpinBox->setRange(0.5, 10);
    blood_vessel_gamma_SpinBox->setSingleStep(0.1);
    blood_vessel_gamma_SpinBox->setValue(3);

    blood_vessel_max_generations_SpinBox = new QDoubleSpinBox;
    blood_vessel_max_generations_SpinBox->setRange(1, 20);
    blood_vessel_max_generations_SpinBox->setSingleStep(1);
    blood_vessel_max_generations_SpinBox->setDecimals(0);
    blood_vessel_max_generations_SpinBox->setValue(7);

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
    glial_pop1_radius_mean_SpinBox->setValue(3);

    glial_pop1_radius_std_SpinBox = new QDoubleSpinBox;
    glial_pop1_radius_std_SpinBox->setRange(0, 10);
    glial_pop1_radius_std_SpinBox->setSingleStep(0.1);
    glial_pop1_radius_std_SpinBox->setValue(0.5);

    glial_pop2_radius_mean_SpinBox = new QDoubleSpinBox;
    glial_pop2_radius_mean_SpinBox->setRange(0, 10);
    glial_pop2_radius_mean_SpinBox->setSingleStep(0.1);
    glial_pop2_radius_mean_SpinBox->setValue(3);

    glial_pop2_radius_std_SpinBox = new QDoubleSpinBox;
    glial_pop2_radius_std_SpinBox->setRange(0, 10);
    glial_pop2_radius_std_SpinBox->setSingleStep(0.1);
    glial_pop2_radius_std_SpinBox->setValue(0.5);

    glial_pop1_minimum_process_radius_SpinBox = new QDoubleSpinBox;
    glial_pop1_minimum_process_radius_SpinBox->setRange(0, 10);
    glial_pop1_minimum_process_radius_SpinBox->setSingleStep(0.01);
    glial_pop1_minimum_process_radius_SpinBox->setValue(0.15);

    glial_pop2_minimum_process_radius_SpinBox = new QDoubleSpinBox;
    glial_pop2_minimum_process_radius_SpinBox->setRange(0, 10);
    glial_pop2_minimum_process_radius_SpinBox->setSingleStep(0.01);
    glial_pop2_minimum_process_radius_SpinBox->setValue(0.15);

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
    
    std::vector<QLabel*> glials_labels1 = {glial_pop1_soma_icvf_qlabel, glial_pop1_processes_icvf_qlabel, glial_pop1_radius_mean_qlabel, glial_pop1_radius_std_qlabel, glial_pop1_mean_process_length_qlabel, glial_pop1_std_process_length_qlabel, glial_pop1_minimum_process_radius_qlabel, glial_pop1_nbr_primary_processes_qlabel};
    std::vector <QDoubleSpinBox*> glials_spinBoxes1 = {glial_pop1_soma_icvf_SpinBox, glial_pop1_processes_icvf_SpinBox, glial_pop1_radius_mean_SpinBox, glial_pop1_radius_std_SpinBox, glial_pop1_mean_process_length_SpinBox, glial_pop1_std_process_length_SpinBox, glial_pop1_minimum_process_radius_SpinBox, glial_pop1_nbr_primary_processes_SpinBox};
    
    std::vector<QLabel*> glials_labels2 = {glial_pop2_soma_icvf_qlabel, glial_pop2_processes_icvf_qlabel, glial_pop2_radius_mean_qlabel, glial_pop2_radius_std_qlabel, glial_pop2_mean_process_length_qlabel, glial_pop2_std_process_length_qlabel, glial_pop2_minimum_process_radius_qlabel, glial_pop2_nbr_primary_processes_qlabel};
    std::vector <QDoubleSpinBox*> glials_spinBoxes2 = {glial_pop2_soma_icvf_SpinBox, glial_pop2_processes_icvf_SpinBox, glial_pop2_radius_mean_SpinBox, glial_pop2_radius_std_SpinBox, glial_pop2_mean_process_length_SpinBox, glial_pop2_std_process_length_SpinBox, glial_pop2_minimum_process_radius_SpinBox, glial_pop2_nbr_primary_processes_SpinBox};

    
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

    onSelectDirectoryButtonClicked();


    // Retrieve values from spin boxes and checkboxes
    parameters.repetitions = nbr_repetitions_SpinBox->value();
    visualise_voxel = visualise_voxel_checkbox->isChecked();

    parameters.data_directory = selectedDirectory.toStdString();
    parameters.filename = "Voxel";

    parameters.axons_wo_myelin_icvf = axons_icvf_SpinBox->value()/100.0;
    parameters.axons_w_myelin_icvf = axons_w_myelin_icvf_SpinBox->value()/100.0;
    parameters.blood_vessels_icvf = blood_vessels_icvf_SpinBox->value()/100.0;
    parameters.blood_vessels_processes_icvf = blood_vessels_processes_icvf_SpinBox->value()/100.0;
    parameters.blood_vessels_voxel_size = blood_vessel_voxel_size_SpinBox->value();
    parameters.mean_vessel_rad = blood_vessel_mean_radius_SpinBox->value();
    parameters.std_vessel_rad = blood_vessel_std_radius_SpinBox->value();
    parameters.blood_vessel_gamma = blood_vessel_gamma_SpinBox->value();
    parameters.max_generations = blood_vessel_max_generations_SpinBox->value();
    parameters.glial_pop1_soma_icvf = glial_pop1_soma_icvf_SpinBox->value()/100.0;
    parameters.glial_pop1_processes_icvf = glial_pop1_processes_icvf_SpinBox->value()/100.0;
    parameters.glial_pop2_soma_icvf = glial_pop2_soma_icvf_SpinBox->value()/100.0;
    parameters.glial_pop2_processes_icvf = glial_pop2_processes_icvf_SpinBox->value()/100.0;
    parameters.glial_pop3_soma_icvf = glial_pop3_soma_icvf_SpinBox->value()/100.0;
    parameters.glial_pop3_processes_icvf = glial_pop3_processes_icvf_SpinBox->value()/100.0;

    parameters.c1 = k1_SpinBox->value();
    parameters.c2 = k2_SpinBox->value();
    parameters.c3 = k3_SpinBox->value();

    parameters.nbr_threads = nbr_threads_SpinBox->value();
    parameters.spheres_overlap_factor = overlapping_factor_SpinBox->value();
    parameters.voxel_size = voxel_size_SpinBox->value();
    parameters.min_rad = minimum_radius_SpinBox->value();
    parameters.cosPhiSquared = c2_SpinBox->value();
    parameters.nbr_axons_populations = nbr_axons_populations_SpinBox->value();
    parameters.beading_amplitude = beading_amplitude_SpinBox->value();
    parameters.beading_std = beading_std_SpinBox->value();

    parameters.mean_glial_pop1_process_length = glial_pop1_mean_process_length_SpinBox->value();
    parameters.std_glial_pop1_process_length = glial_pop1_std_process_length_SpinBox->value();
    parameters.mean_glial_pop2_process_length = glial_pop2_mean_process_length_SpinBox->value();
    parameters.std_glial_pop2_process_length = glial_pop2_std_process_length_SpinBox->value();
    parameters.mean_glial_pop3_process_length = glial_pop3_mean_process_length_SpinBox->value();
    parameters.std_glial_pop3_process_length = glial_pop3_std_process_length_SpinBox->value();

    parameters.epsilon = epsilon_SpinBox->value();
    parameters.alpha = alpha_SpinBox->value();
    parameters.beta = beta_SpinBox->value();
    parameters.alpha_myelin = alpha_myelin_SpinBox->value();
    parameters.beta_myelin = beta_myelin_SpinBox->value();

    parameters.glial_pop1_radius_mean = glial_pop1_radius_mean_SpinBox->value();
    parameters.glial_pop1_radius_std = glial_pop1_radius_std_SpinBox->value();
    parameters.glial_pop2_radius_mean = glial_pop2_radius_mean_SpinBox->value();
    parameters.glial_pop2_radius_std = glial_pop2_radius_std_SpinBox->value();
    parameters.glial_pop3_radius_mean = glial_pop3_radius_mean_SpinBox->value();
    parameters.glial_pop3_radius_std = glial_pop3_radius_std_SpinBox->value();
    parameters.glial_pop1_minimum_process_radius = glial_pop1_minimum_process_radius_SpinBox->value();
    parameters.glial_pop2_minimum_process_radius = glial_pop2_minimum_process_radius_SpinBox->value();
    parameters.glial_pop3_minimum_process_radius = glial_pop3_minimum_process_radius_SpinBox->value();
    parameters.glial_pop1_nbr_primary_processes = glial_pop1_nbr_primary_processes_SpinBox->value();
    parameters.glial_pop2_nbr_primary_processes = glial_pop2_nbr_primary_processes_SpinBox->value();
    parameters.glial_pop3_nbr_primary_processes = glial_pop3_nbr_primary_processes_SpinBox->value();
    parameters.glial_pop1_branching = glial_pop1_branching_checkbox->isChecked();
    parameters.glial_pop2_branching = glial_pop2_branching_checkbox->isChecked();
    parameters.glial_pop3_branching = glial_pop3_branching_checkbox->isChecked();

    parameters.crossing_fibers_type = 0;

    if (parameters.nbr_axons_populations == 2) {
        QString PopConfiguration = configurationComboBox->currentText();
        if (PopConfiguration == "Sheet Configuration") {
            parameters.crossing_fibers_type = 0;
        } else if (PopConfiguration == "Interwoven Configuration") {
            parameters.crossing_fibers_type = 1;
        }
    }

    // Close the parameter input dialog
    this->setEnabled(false);

    StartSimulation();

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
bool Window::parseSubstrateCsvLine(const std::string &line, std::string &type, double &id_cell,
                                    std::string &component, double &component_id, double &parent,
                                    double &x, double &y, double &z, double &radius_in, double &radius_out) {
    std::istringstream iss(line);
    std::vector<std::string> tokens;
    std::string tok;
    while (iss >> tok) {
        tokens.push_back(tok);
    }

    size_t n = tokens.size();
    if (n != 9 && n != 10) {
        return false;
    }

    size_t idx = 0;
    type = tokens[idx++];
    id_cell = std::stod(tokens[idx++]);
    component = tokens[idx++];
    component_id = std::stod(tokens[idx++]);
    if (n == 10) {
        parent = std::stod(tokens[idx++]);
    } else {
        // Pre-parent_component_id format: treat as having no distinct
        // parent, matching how the writer fills that column today for
        // objects like axons (parent_component_id = component_id).
        parent = component_id;
    }
    x = std::stod(tokens[idx++]);
    y = std::stod(tokens[idx++]);
    z = std::stod(tokens[idx++]);
    radius_in = std::stod(tokens[idx++]);
    radius_out = std::stod(tokens[idx++]);
    return true;
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
        double id_cell, component_id;
        std::string type, component;
        double x, y, z, radius_in, parent, radius_out;

        //skip first line
        if (line[0] == 'c') {
            continue;
        }


        if (!parseSubstrateCsvLine(line, type, id_cell, component, component_id, parent, x, y, z, radius_in, radius_out)) {
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

        double component_id, id_cell;
        std::string type, component;

        double x, y, z, radius_in, parent, radius_out;

        if (!parseSubstrateCsvLine(line, type, id_cell, component, component_id, parent, x, y, z, radius_in, radius_out)) {
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
        QMessageBox::warning(this, tr("Error"), tr("Could not open the CSV file."));
        return;
    }

    std::vector<double> x_ = {};
    std::vector<double> y_ = {};
    std::vector<double> z_ = {};
    std::vector<double> r_ = {};

    double old_cell_id = -1;

    std::string line;
    while (std::getline(swcFile, line)) {
        //skip first line
        if (line[0] == 'c') {
            continue;
        }

        double id_cell, component_id, parent;
        std::string type, component;
        double x, y, z, radius_in, radius_out;

        if (!parseSubstrateCsvLine(line, type, id_cell, component, component_id, parent, x, y, z, radius_in, radius_out)) {
            QMessageBox::warning(this, tr("Error"), tr("Invalid CSV file format for Blood Vessels."));
            return;
        }

        if (type == "blood_vessel") {
            if (id_cell != old_cell_id) {
                if (old_cell_id != -1 && x_.size() > 0) {
                    X_blood_vessels.push_back(x_);
                    Y_blood_vessels.push_back(y_);
                    Z_blood_vessels.push_back(z_);
                    R_blood_vessels.push_back(r_);
                }
                x_.clear();
                y_.clear();
                z_.clear();
                r_.clear();
                old_cell_id = id_cell;
            }
            x_.push_back(x);
            y_.push_back(y);
            z_.push_back(z);
            r_.push_back(radius_out);
        }
    }
    // Add the last vessel
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


void Window::onGrowthProgress(double completed_depth, double total_depth)
{
    int pct = total_depth > 0.0 ? int(100.0 * std::min(1.0, completed_depth / total_depth)) : 0;
    layerProgressBar->setValue(pct);
    layerProgressBar->setFormat(QString("Growing axons: %1 / %2 um (%3%)")
                                    .arg(completed_depth, 0, 'f', 1)
                                    .arg(total_depth, 0, 'f', 1)
                                    .arg(pct));
}

void Window::onSwellingProgress(double current_icvf, double target_icvf)
{
    int pct = target_icvf > 0.0 ? int(100.0 * std::min(1.0, current_icvf / target_icvf)) : 100;
    swellingProgressBar->setValue(pct);
    swellingProgressBar->setFormat(QString("Swelling: ICVF %1 / %2 (%3%)")
                                       .arg(current_icvf, 0, 'f', 4)
                                       .arg(target_icvf, 0, 'f', 4)
                                       .arg(pct));
}

namespace {
// RAII guard redirecting std::cout to another stream for its lifetime --
// used below to capture a growth run's console output into its own .log
// file, mirroring the .json/.log pair already saved for CLI-driven runs.
// Safe here specifically because this app only ever runs one growth job at
// a time (growthThread is joined before a new one starts) and the main GUI
// thread doesn't write to std::cout while a growth is in flight.
class CoutRedirect {
public:
    explicit CoutRedirect(std::ostream& new_stream)
        : old_buf(std::cout.rdbuf(new_stream.rdbuf())) {}
    ~CoutRedirect() { std::cout.rdbuf(old_buf); }
private:
    std::streambuf* old_buf;
};
}

void Window::StartSimulation(){

    X_axons.clear(); Y_axons.clear(); Z_axons.clear(); R_axons.clear();
    X_glial_pop1.clear(); Y_glial_pop1.clear(); Z_glial_pop1.clear(); R_glial_pop1.clear(); Branch_glial_pop1.clear();
    X_glial_pop2.clear(); Y_glial_pop2.clear(); Z_glial_pop2.clear(); R_glial_pop2.clear(); Branch_glial_pop2.clear();
    X_glial_pop3.clear(); Y_glial_pop3.clear(); Z_glial_pop3.clear(); R_glial_pop3.clear(); Branch_glial_pop3.clear();
    X_blood_vessels.clear(); Y_blood_vessels.clear(); Z_blood_vessels.clear(); R_blood_vessels.clear();

    layerProgressBar->setValue(0);
    layerProgressBar->setFormat("Growing axons: 0%");
    layerProgressBar->show();
    swellingProgressBar->setValue(0);
    swellingProgressBar->setFormat("Swelling: not started yet");
    swellingProgressBar->show();

    Parameters localParams = parameters; // by-value copy: safe to read from the worker thread

    // Save a record of exactly what this run was configured with, alongside
    // its console output (.log, captured below) -- mirroring the .json/.log
    // pair already saved for CLI-driven runs, so a GUI run's setup and full
    // growth log are both on disk even if nothing else was recorded.
    // Written synchronously here (not in the worker thread) so it exists the
    // moment the run starts, even if growth is later interrupted.
    std::string basePath = localParams.data_directory + "/" + localParams.filename;
    CoreLogic::writeParametersToJson(localParams, basePath + ".json");

    if (growthThread.joinable()) {
        growthThread.join(); // previous run's thread, if any, has already finished by now
    }

    growthThread = std::thread([this, localParams, basePath]() {
        auto growthCb = [this](double completed, double total) {
            QMetaObject::invokeMethod(this, "onGrowthProgress", Qt::QueuedConnection,
                                       Q_ARG(double, completed), Q_ARG(double, total));
        };
        auto swellCb = [this](double current_icvf, double target_icvf) {
            QMetaObject::invokeMethod(this, "onSwellingProgress", Qt::QueuedConnection,
                                       Q_ARG(double, current_icvf), Q_ARG(double, target_icvf));
        };

        std::ofstream log_file(basePath + ".log");
        CoutRedirect redirect(log_file);

        auto [axons, blood_vessels, glial_pop1, glial_pop2, glial_pop3, voxel_min, voxel_max] =
            CoreLogic::runSimulation(localParams, growthCb, swellCb);

        this->pendingAxons = std::move(axons);
        this->pendingBloodVessels = std::move(blood_vessels);
        this->pendingGlialPop1 = std::move(glial_pop1);
        this->pendingGlialPop2 = std::move(glial_pop2);
        this->pendingGlialPop3 = std::move(glial_pop3);
        this->pendingVoxelMin = voxel_min;
        this->pendingVoxelMax = voxel_max;

        QMetaObject::invokeMethod(this, "onGrowthFinished", Qt::QueuedConnection);
    });
}

void Window::onGrowthFinished(){

    if (growthThread.joinable()) {
        growthThread.join();
    }

    layerProgressBar->hide();
    swellingProgressBar->hide();

    auto &axons = pendingAxons;
    auto &blood_vessels = pendingBloodVessels;
    auto &glial_pop1 = pendingGlialPop1;
    auto &glial_pop2 = pendingGlialPop2;
    auto &glial_pop3 = pendingGlialPop3;

    // Real (small) voxel's actual grown placement (see PlaceSmallVoxel) --
    // not necessarily an origin-anchored box, e.g. once blood vessels are
    // involved (its own face-seeded plane can put it anywhere inside the
    // padded blood-vessel voxel). Used below for the 3D view's wireframe box.
    voxelBoundsMin = QVector3D(pendingVoxelMin[0], pendingVoxelMin[1], pendingVoxelMin[2]);
    voxelBoundsMax = QVector3D(pendingVoxelMax[0], pendingVoxelMax[1], pendingVoxelMax[2]);

    Eigen::Vector3d min_l = {0,0,0};
    Eigen::Vector3d max_l ={parameters.voxel_size, parameters.voxel_size, parameters.voxel_size};

    for (unsigned i=0; i< axons.size(); ++i){
        std::vector<double> x_;
        std::vector<double> y_;
        std::vector<double> z_;
        std::vector<double> r_;
        for (unsigned j=0; j< axons[i].outer_spheres.size(); ++j){

            double _x_ = axons[i].outer_spheres[j].center[0];
            double _y_ = axons[i].outer_spheres[j].center[1];
            double _z_ = axons[i].outer_spheres[j].center[2];
            double _r_ = axons[i].outer_spheres[j].radius;
            Eigen::Vector3d pos = {_x_, _y_, _z_};

            //if (check_borders(min_l, max_l, pos, 0.0)) {
            x_.push_back(_x_);
            y_.push_back(_y_);
            z_.push_back(_z_);
            r_.push_back(_r_);
            //}
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

    for (unsigned i=0; i< glial_pop1.size(); ++i){
        std::vector<double> x_;
        std::vector<double> y_;
        std::vector<double> z_;
        std::vector<double> r_;
        std::vector<int> b_;

        double _x_ = glial_pop1[i].soma.center[0];
        double _y_ = glial_pop1[i].soma.center[1];
        double _z_ = glial_pop1[i].soma.center[2];
        double _r_ = glial_pop1[i].soma.radius;
        Eigen::Vector3d pos = {_x_, _y_, _z_};

        //if (check_borders(min_l, max_l, pos, 0.0)) {
        x_.push_back(_x_);
        y_.push_back(_y_);
        z_.push_back(_z_);
        r_.push_back(_r_);
        b_.push_back(0);
        //}


        for (unsigned j=0; j< glial_pop1[i].ramification_spheres.size(); ++j){

            for (unsigned k=0; k< glial_pop1[i].ramification_spheres[j].size(); ++k){

                double _x_ = glial_pop1[i].ramification_spheres[j][k].center[0];
                double _y_ = glial_pop1[i].ramification_spheres[j][k].center[1];
                double _z_ = glial_pop1[i].ramification_spheres[j][k].center[2];
                double _r_ = glial_pop1[i].ramification_spheres[j][k].radius;
                Eigen::Vector3d pos = {_x_, _y_, _z_};
                //if (!check_borders(min_l, max_l, pos, 0.0)) {
                //    continue;
                //}
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

    for (unsigned i=0; i< glial_pop2.size(); ++i){
        std::vector<double> x_;
        std::vector<double> y_;
        std::vector<double> z_;
        std::vector<double> r_;
        std::vector<int> b_;

        double _x_ = glial_pop2[i].soma.center[0];
        double _y_ = glial_pop2[i].soma.center[1];
        double _z_ = glial_pop2[i].soma.center[2];
        double _r_ = glial_pop2[i].soma.radius;

        //if (check_borders(min_l, max_l, {_x_, _y_, _z_}, 0.0)) {
        x_.push_back(_x_);
        y_.push_back(_y_);
        z_.push_back(_z_);
        r_.push_back(_r_);
        b_.push_back(0);
        //}


        for (unsigned j=0; j< glial_pop2[i].ramification_spheres.size(); ++j){

            for (unsigned k=0; k< glial_pop2[i].ramification_spheres[j].size(); ++k){

                double _x_ = glial_pop2[i].ramification_spheres[j][k].center[0];
                double _y_ = glial_pop2[i].ramification_spheres[j][k].center[1];
                double _z_ = glial_pop2[i].ramification_spheres[j][k].center[2];
                double _r_ = glial_pop2[i].ramification_spheres[j][k].radius;

                //if (!check_borders(min_l, max_l, {_x_, _y_, _z_}, 0.0)) {
                //    continue;
                //}
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

    for (unsigned i=0; i< glial_pop3.size(); ++i){
        std::vector<double> x_;
        std::vector<double> y_;
        std::vector<double> z_;
        std::vector<double> r_;
        std::vector<int> b_;

        double _x_ = glial_pop3[i].soma.center[0];
        double _y_ = glial_pop3[i].soma.center[1];
        double _z_ = glial_pop3[i].soma.center[2];
        double _r_ = glial_pop3[i].soma.radius;

        //if (check_borders(min_l, max_l, {_x_, _y_, _z_}, 0.0)) {
        x_.push_back(_x_);
        y_.push_back(_y_);
        z_.push_back(_z_);
        r_.push_back(_r_);
        b_.push_back(0);
        //}


        for (unsigned j=0; j< glial_pop3[i].ramification_spheres.size(); ++j){

            for (unsigned k=0; k< glial_pop3[i].ramification_spheres[j].size(); ++k){

                double _x_ = glial_pop3[i].ramification_spheres[j][k].center[0];
                double _y_ = glial_pop3[i].ramification_spheres[j][k].center[1];
                double _z_ = glial_pop3[i].ramification_spheres[j][k].center[2];
                double _r_ = glial_pop3[i].ramification_spheres[j][k].radius;

                //if (!check_borders(min_l, max_l, {_x_, _y_, _z_}, 0.0)) {
                //    continue;
                //}
                x_.push_back(_x_);
                y_.push_back(_y_);
                z_.push_back(_z_);
                r_.push_back(_r_);
                b_.push_back(j);
            }
        }

        X_glial_pop3.push_back(x_);
        Y_glial_pop3.push_back(y_);
        Z_glial_pop3.push_back(z_);
        R_glial_pop3.push_back(r_);
        Branch_glial_pop3.push_back(b_);
        x_.clear();
        y_.clear();
        z_.clear();
        r_.clear();
        b_.clear();
    }

    for (unsigned i=0; i< blood_vessels.size(); ++i){
        std::vector<double> x_;
        std::vector<double> y_;
        std::vector<double> z_;
        std::vector<double> r_;
        for (const auto &branch : blood_vessels[i].ramification_spheres){
        for (unsigned j=0; j< branch.size(); ++j){

            double _x_ = branch[j].center[0];
            double _y_ = branch[j].center[1];
            double _z_ = branch[j].center[2];
            double _r_ = branch[j].radius;
            Eigen::Vector3d pos = {_x_, _y_, _z_};

            //if (check_borders(min_l, max_l, pos, 0.0)) {
            x_.push_back(_x_);
            y_.push_back(_y_);
            z_.push_back(_z_);
            r_.push_back(_r_);
            //}
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


    if (visualise_voxel) {
        // After simulation completes, call PlotCells to display the data
        PlotCells(true, true, true, true, true);
    }
    else{
        // Display a message box to inform the user that the simulation is complete
        QMessageBox::information(this, "Simulation Complete", "Simulation complete! Please check the output directory for the results.");
    }
    // Add this to the very end of StartSimulation:
    this->setEnabled(true);

}

void Window::onSelectDirectoryButtonClicked() {
    QString dirPath = QFileDialog::getExistingDirectory(
        this, 
        tr("Select Output Directory for Substrate"), 
        "", 
        QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks
    );

    // 2. Check if the user clicked "Cancel" or closed the window
    if (dirPath.isEmpty()) {
        qDebug() << "Growth aborted: No directory selected.";
        return; // Halt execution cleanly
    }

    // 3. Save the path for your backend to use
    selectedDirectory = dirPath;
    
    qDebug() << "Proceeding with growth. Output directory:" << selectedDirectory;
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

    // Create the histogram as a bar chart: one category per bin, labelled
    // with that bin's center value.
    QBarSet *barSet = new QBarSet("Count");
    QStringList categories;
    for (int i = 0; i < binCount; ++i) {
        *barSet << bins[i];
        categories << QString::number(tickPositions[i], 'f', 2);
    }

    QBarSeries *series = new QBarSeries();
    series->append(barSet);

    QChart *chart = new QChart();
    chart->addSeries(series);
    chart->setTitle("Mean Radius Distribution");
    chart->legend()->hide();

    QBarCategoryAxis *axisX = new QBarCategoryAxis();
    axisX->append(categories);
    axisX->setTitleText("Mean Radius");
    axisX->setLabelsAngle(-90);
    chart->addAxis(axisX, Qt::AlignBottom);
    series->attachAxis(axisX);

    QValueAxis *axisY = new QValueAxis();
    axisY->setTitleText("Count");
    axisY->setRange(0, *std::max_element(bins.begin(), bins.end()));
    chart->addAxis(axisY, Qt::AlignLeft);
    series->attachAxis(axisY);

    QChartView *chartView = new QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);

    // Display the plot in a dialog window
    QDialog *dialog = new QDialog(this);
    QVBoxLayout *layout = new QVBoxLayout;
    layout->addWidget(chartView);
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

    // Create the histogram as a bar chart: one category per bin, labelled
    // with that bin's center value.
    QBarSet *barSet = new QBarSet("Count");
    QStringList categories;
    for (int i = 0; i < binCount; ++i) {
        *barSet << bins[i];
        categories << QString::number(tickPositions[i], 'f', 2);
    }

    QBarSeries *series = new QBarSeries();
    series->append(barSet);

    QChart *chart = new QChart();
    chart->addSeries(series);
    chart->setTitle("Tortuosity Distribution");
    chart->legend()->hide();

    QBarCategoryAxis *axisX = new QBarCategoryAxis();
    axisX->append(categories);
    axisX->setTitleText("Tortuosity");
    chart->addAxis(axisX, Qt::AlignBottom);
    series->attachAxis(axisX);

    QValueAxis *axisY = new QValueAxis();
    axisY->setTitleText("Count");
    axisY->setRange(0, *std::max_element(bins.begin(), bins.end()));
    chart->addAxis(axisY, Qt::AlignLeft);
    series->attachAxis(axisY);

    QChartView *chartView = new QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);

    // Display the plot in a dialog window
    QDialog *dialog = new QDialog(this);
    QVBoxLayout *layout = new QVBoxLayout;
    layout->addWidget(chartView);
    dialog->setLayout(layout);
    dialog->setWindowTitle("Tortuosity Distribution");
    dialog->exec();
}


void Window::ShollAnalysis() {

    size_t total_cells = X_glial_pop1.size() + X_glial_pop2.size() + X_glial_pop3.size();
    if (total_cells == 0) {
        QMessageBox::warning(this, "Error", "No Glial cells to plot!");
        return;
    }

    // Radii for Sholl analysis
    std::vector<double> sphere_around_soma_radii = {5, 7, 10, 15, 20, 25, 30, 40, 50, 60, 80};

    // Computes one population's mean Sholl curve and opens its own dialog --
    // called once per population below, and skipped entirely for any
    // population with no cells, so (e.g.) growing only pop1 doesn't pop up
    // empty/meaningless windows for pop2 and pop3.
    auto plot_population = [&](const std::vector<std::vector<double>> &X,
                                const std::vector<std::vector<double>> &Y,
                                const std::vector<std::vector<double>> &Z,
                                const std::vector<std::vector<double>> &R,
                                const std::vector<std::vector<int>> &Branch,
                                const QString &title) {
        if (X.empty()) return;

        std::vector<double> mean_intersections(sphere_around_soma_radii.size(), 0);
        unsigned long nbr_cells_included = 0;

        for (unsigned long i = 0; i < X.size(); ++i) {
            Eigen::Vector3d soma_position = {X[i][0], Y[i][0], Z[i][0]};

            // Only cells whose soma actually lies within the (small, real)
            // voxel currently displayed -- somas seeded/swollen outside it
            // (see PlaceGlialCells/SwellGlialSomas) would otherwise skew the
            // curve with cells that aren't really part of this substrate.
            if (soma_position[0] < voxelBoundsMin.x() || soma_position[0] > voxelBoundsMax.x() ||
                soma_position[1] < voxelBoundsMin.y() || soma_position[1] > voxelBoundsMax.y() ||
                soma_position[2] < voxelBoundsMin.z() || soma_position[2] > voxelBoundsMax.z()) {
                continue;
            }
            ++nbr_cells_included;

            std::vector<double> intersections_list(sphere_around_soma_radii.size(), 0);
            std::vector<int> branches_list;

            // Iterate through all spheres (excluding the soma) to compute intersections
            for (unsigned long r = 0; r < sphere_around_soma_radii.size(); ++r) {
                for (unsigned long j = 1; j < X[i].size(); ++j) {
                    Eigen::Vector3d position = {X[i][j], Y[i][j], Z[i][j]};
                    double distance = (position - soma_position).norm();

                    if (distance < sphere_around_soma_radii[r] + R[i][j] && distance > sphere_around_soma_radii[r] - R[i][j]) {
                        if (std::find(branches_list.begin(), branches_list.end(), Branch[i][j]) == branches_list.end()) {
                            intersections_list[r] += 1;
                            branches_list.push_back(Branch[i][j]);
                        }
                    }
                }
                branches_list.clear();
            }

            for (size_t r = 0; r < sphere_around_soma_radii.size(); ++r) {
                mean_intersections[r] += intersections_list[r];
            }
        }

        if (nbr_cells_included == 0) return;

        for (size_t r = 0; r < mean_intersections.size(); ++r) {
            mean_intersections[r] /= nbr_cells_included;
        }

        // Create a line series (with markers at each point) for the mean Sholl curve
        QLineSeries *series = new QLineSeries();
        for (size_t r = 0; r < sphere_around_soma_radii.size(); ++r) {
            series->append(sphere_around_soma_radii[r], mean_intersections[r]);
        }
        series->setPointsVisible(true);

        QChart *chart = new QChart();
        chart->addSeries(series);
        chart->setTitle(title);
        chart->legend()->hide();

        QValueAxis *axisX = new QValueAxis();
        axisX->setTitleText("Distance to Soma (μm)");
        axisX->setRange(0, *std::max_element(sphere_around_soma_radii.begin(), sphere_around_soma_radii.end()));
        chart->addAxis(axisX, Qt::AlignBottom);
        series->attachAxis(axisX);

        QValueAxis *axisY = new QValueAxis();
        axisY->setTitleText("Mean Number of Intersections");
        axisY->setRange(0, *std::max_element(mean_intersections.begin(), mean_intersections.end()));
        chart->addAxis(axisY, Qt::AlignLeft);
        series->attachAxis(axisY);

        QChartView *chartView = new QChartView(chart);
        chartView->setRenderHint(QPainter::Antialiasing);

        // Display the plot in a dialog window
        QDialog *dialog = new QDialog(this);
        dialog->resize(800, 600);
        QVBoxLayout *layout = new QVBoxLayout;
        layout->addWidget(chartView);
        dialog->setLayout(layout);
        dialog->setWindowTitle(title);
        dialog->exec();
    };

    plot_population(X_glial_pop1, Y_glial_pop1, Z_glial_pop1, R_glial_pop1, Branch_glial_pop1, "Mean Sholl Analysis - Glial Population 1");
    plot_population(X_glial_pop2, Y_glial_pop2, Z_glial_pop2, R_glial_pop2, Branch_glial_pop2, "Mean Sholl Analysis - Glial Population 2");
    plot_population(X_glial_pop3, Y_glial_pop3, Z_glial_pop3, R_glial_pop3, Branch_glial_pop3, "Mean Sholl Analysis - Glial Population 3");
}
