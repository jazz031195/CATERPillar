QT       += core gui datavisualization
QT += printsupport
QT += opengl


greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++11

INCLUDEPATH += ../src

INCLUDEPATH += /usr/include/GL

SOURCES += \
    main.cpp \
    mainwindow.cpp \
    slidergroup.cpp \
    openglwindow.cpp \
    ../src/Axon.cpp \
    ../src/axongammadistribution.cpp \
    ../src/Glial.cpp \
    ../src/grow_axons.cpp \
    ../src/grow_cells.cpp \
    ../src/grow_glial_cells.cpp \
    ../src/obstacle.cpp \
    ../src/parameters.cpp \
    ../src/sphere.cpp \
    ../src/threads.cpp \
    
SOURCES += qcustomplot-source/qcustomplot.cpp

HEADERS += \
    mainwindow.h \
    slidergroup.h \
    openglwindow.h \
    ../src/Axon.h \
    ../src/axongammadistribution.h \
    ../src/Glial.h \
    ../src/grow_axons.h \
    ../src/grow_cells.h \
    ../src/grow_glial_cells.h \
    ../src/obstacle.h \
    ../src/parameters.h \
    ../src/sphere.h \
    ../src/threads.h \
    ../src/constants.h \
    qcustomplot-source/qcustomplot.h


FORMS += \
    mainwindow.ui

LIBS += -lGL -lGLU -lglut

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

