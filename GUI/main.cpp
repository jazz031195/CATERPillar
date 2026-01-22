#include <QApplication>
#include <QHBoxLayout>
#include <QPushButton>
#include <QCheckBox>
#include <QLabel>
#include "mainwindow.h"

int main(int argc, char *argv[]) {

    // If a config file is passed, we know we are in batch/CLI mode
    if (argc > 1) {
        qputenv("QT_QPA_PLATFORM", "offscreen");  
        qputenv("SESSION_MANAGER", QByteArray()); 
    }

    QApplication app(argc, argv);
    Window window;

    // If config file is given in command line arguments, load it and start simulation directly
    if (argc > 1) {
        QString cfgPath = QString::fromLocal8Bit(argv[1]);
        window.loadConfigFromPath(cfgPath); // Load configuration file from the provided path
        window.updateMembersFromGui(); // Update variables for simulations from GUI elements
        window.StartSimulation(); // Start the simulation
        return 0;
    }

    // Otherwise, show the GUI
    window.show();

    return app.exec();
}
