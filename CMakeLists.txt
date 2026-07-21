cmake_minimum_required(VERSION 3.10)
project(CATERPillar VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BUILD_GUI "Build the Qt Widgets GUI application (CATERPillar)" ON)

# Eigen is vendored at src/Eigen -- already on the include path below, no system
# Eigen3 install needed. The headless CLI target below needs no external
# libraries at all: just a C++17 compiler (std::filesystem) and CMake.

file(GLOB CORE_SOURCES
    "src/*.cpp"
    "src/*.h"
)
list(REMOVE_ITEM CORE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/main_cli.cpp
)

# --- Headless CLI target: just the growth simulation, no external libraries required ---
add_executable(CATERPillar-cli ${CORE_SOURCES} src/main_cli.cpp)

target_include_directories(CATERPillar-cli PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# --- GUI target: full CATERPillar application ---
if(BUILD_GUI)
    find_package(Qt5 COMPONENTS Widgets PrintSupport DataVisualization REQUIRED)
    find_package(OpenGL REQUIRED)

    file(GLOB_RECURSE GUI_SOURCES
        "GUI/*.cpp"
        "GUI/*.h"
        "GUI/*.ui"
        "GUI/*.qrc"
    )

    add_executable(CATERPillar ${CORE_SOURCES} src/main.cpp ${GUI_SOURCES})

    set_target_properties(CATERPillar PROPERTIES
        AUTOMOC ON
        AUTOUIC ON
        AUTORCC ON
    )

    target_include_directories(CATERPillar PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/GUI
    )

    target_link_libraries(CATERPillar PRIVATE Qt5::Widgets Qt5::PrintSupport Qt5::DataVisualization OpenGL::GL)

    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/logo_catrepillar.png DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()
