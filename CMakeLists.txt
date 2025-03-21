cmake_minimum_required(VERSION 3.10)
project(LatticeTriangulation)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(Eigen3 REQUIRED)
find_package(CGAL REQUIRED)
find_package(Boost REQUIRED)

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${EIGEN3_INCLUDE_DIR}
    ${CGAL_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

# Create a list of library source files (excluding main.cpp)
file(GLOB_RECURSE LIB_SOURCES 
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)
list(FILTER LIB_SOURCES EXCLUDE REGEX ".*main\\.cpp$")

# Create a library
add_library(lattice_lib ${LIB_SOURCES})

# Create main executable with main.cpp
add_executable(lattice_triangulation src/main.cpp)
target_link_libraries(lattice_triangulation 
    lattice_lib
    ${CGAL_LIBRARIES}
    ${Boost_LIBRARIES}
)

# Create test executable with test files only
file(GLOB TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp")
if(TEST_SOURCES)
    add_executable(run_tests ${TEST_SOURCES})
    target_link_libraries(run_tests 
        lattice_lib
        ${CGAL_LIBRARIES}
        ${Boost_LIBRARIES}
    )
endif()
