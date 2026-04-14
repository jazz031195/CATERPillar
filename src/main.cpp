#include <QApplication>
#include <iostream>
#include <string>
#include "mainwindow.h" 
#include "core_logic.h" 

int main(int argc, char *argv[]) {
    // 1. Check if the user passed arguments in the terminal for JSON mode
    if (argc >= 3 && std::string(argv[1]) == "--config") {
        std::string jsonFilePath = argv[2];
        std::cout << "Starting in headless mode using: " << jsonFilePath << std::endl;
        
        // Pass the file path to your CoreLogic class
        CoreLogic::runSimulationFromJson(jsonFilePath); 
        
        return 0; // Exit without ever loading Qt graphical libraries
    }

    // 2. If no arguments were passed, launch the GUI normally
    QApplication app(argc, argv);
    
    // Make sure this matches exactly what is in your mainwindow.h file!
    // If your class is named Window instead of MainWindow, change it here.
    Window w; 
    w.show();
    
    return app.exec();
}