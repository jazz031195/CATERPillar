cmake_minimum_required(VERSION 3.10)
project(CATERPillar VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

find_package(Qt5 COMPONENTS Widgets PrintSupport DataVisualization REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(OpenGL REQUIRED)

file(GLOB_RECURSE ALL_SOURCES 
    "src/*.cpp" 
    "src/*.h" 
    "GUI/*.cpp" 
    "GUI/*.h" 
    "GUI/*.ui"
    "GUI/*.qrc"
)

add_executable(CATERPillar ${ALL_SOURCES})

target_include_directories(CATERPillar PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR} 
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/GUI
)

target_link_libraries(CATERPillar PRIVATE Qt5::Widgets Qt5::PrintSupport Qt5::DataVisualization Eigen3::Eigen OpenGL::GL)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/logo_catrepillar.png DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
