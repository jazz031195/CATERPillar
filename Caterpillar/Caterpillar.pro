QT       += core gui datavisualization
QT += printsupport
QT += opengl


greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++11

INCLUDEPATH += ../src/Eigen

SOURCES += \
    ScatterDataModifier.cpp \
    main.cpp \
    mainwindow.cpp \
    slidergroup.cpp \
    openglwindow.cpp \
    ../src/Axon.cpp \
    ../src/axongammadistribution.cpp \
    ../src/Glial.cpp \
    ../src/grow_axons.cpp \
    ../src/obstacle.cpp \
    ../src/parameters.cpp \
    ../src/sphere.cpp \
    ../src/threads.cpp \
    qcustomplot-source/qcustomplot.cpp

HEADERS += \
    ScatterDataModifier.h \
    mainwindow.h \
    slidergroup.h \
    openglwindow.h \
    ../src/Axon.h \
    ../src/axongammadistribution.h \
    ../src/Glial.h \
    ../src/grow_axons.h \
    ../src/obstacle.h \
    ../src/parameters.h \
    ../src/sphere.h \
    ../src/threads.h \
    ../src/constants.h \
    qcustomplot-source/qcustomplot.h


FORMS += \
    mainwindow.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

