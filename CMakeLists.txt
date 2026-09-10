cmake_minimum_required(VERSION 3.10)
project(CATERPillar VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Default to an optimised build. Without this, a plain `cmake ..` leaves
# CMAKE_BUILD_TYPE empty, which means NO -O flag at all (i.e. -O0) -- and
# nothing here or in the toolchain file supplies one. That is not a marginal
# difference for this workload: the same 40um/60%-ICVF substrate takes 19
# seconds at -O3 and had not finished after 48 minutes unoptimised, roughly a
# 150x gap. A whole dictionary was generated (and its runtime repeatedly
# analysed) on an unoptimised binary before this was noticed, so make the
# optimised build the thing you get by default and require an explicit
# -DCMAKE_BUILD_TYPE=Debug to opt out.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type (Release/Debug/RelWithDebInfo)" FORCE)
    message(STATUS "CMAKE_BUILD_TYPE not set -- defaulting to Release (-O3)")
endif()

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
    find_package(Qt5 COMPONENTS Widgets Charts DataVisualization REQUIRED)
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

    target_link_libraries(CATERPillar PRIVATE Qt5::Widgets Qt5::Charts Qt5::DataVisualization OpenGL::GL)

    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/logo_catrepillar.png DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()
