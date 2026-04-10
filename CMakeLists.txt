# 1. Minimum version and Project Name
cmake_minimum_required(VERSION 3.10)
project(CATERPillar VERSION 1.0 LANGUAGES CXX)

# 2. C++ Standard (Using 14 or 17 based on how you fixed Eigen)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 3. Qt-Specific Magic
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

# 4. Find Qt5, Eigen3, AND OpenGL (NEW)
find_package(Qt5 COMPONENTS Widgets PrintSupport DataVisualization REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(OpenGL REQUIRED)

# 5. SAFE file searching! 
# This only looks in your src and GUI folders, ignoring the build folder and backups
file(GLOB_RECURSE ALL_SOURCES 
    "src/*.cpp" 
    "src/*.h" 
    "GUI/*.cpp" 
    "GUI/*.h" 
    "GUI/*.ui"
)

# 6. Build executable
add_executable(CATERPillar ${ALL_SOURCES})

# 7. Include directories (tell the compiler where to find headers)
target_include_directories(CATERPillar PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR} 
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/GUI
)

# 8. Link Qt5, Eigen3, AND OpenGL to your executable (NEW)
target_link_libraries(CATERPillar PRIVATE Qt5::Widgets Qt5::PrintSupport Qt5::DataVisualization Eigen3::Eigen OpenGL::GL)