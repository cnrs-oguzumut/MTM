cmake_minimum_required(VERSION 3.10)
project(LatticeTriangulation)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(Eigen3 REQUIRED)
find_package(CGAL REQUIRED)
find_package(Boost REQUIRED)
find_package(OpenMP)

# Set ALGLIB directory
set(ALGLIB_DIR "/Users/usalman/programming/alglib-cpp")
include_directories(${ALGLIB_DIR})

# Add ALGLIB source files
set(ALGLIB_SOURCES
    ${ALGLIB_DIR}/src/alglibinternal.cpp
    ${ALGLIB_DIR}/src/alglibmisc.cpp
    ${ALGLIB_DIR}/src/ap.cpp
    ${ALGLIB_DIR}/src/dataanalysis.cpp
    ${ALGLIB_DIR}/src/diffequations.cpp
    ${ALGLIB_DIR}/src/fasttransforms.cpp
    ${ALGLIB_DIR}/src/integration.cpp
    ${ALGLIB_DIR}/src/interpolation.cpp
    ${ALGLIB_DIR}/src/optimization.cpp
    ${ALGLIB_DIR}/src/linalg.cpp
    ${ALGLIB_DIR}/src/solvers.cpp
    ${ALGLIB_DIR}/src/specialfunctions.cpp
    ${ALGLIB_DIR}/src/statistics.cpp
)

# Create ALGLIB library
add_library(alglib STATIC ${ALGLIB_SOURCES})

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${EIGEN3_INCLUDE_DIR}
    ${CGAL_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

# Add OpenMP flags if available
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Create a library from source files (excluding main.cpp)
file(GLOB_RECURSE LIB_SOURCES 
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)
list(FILTER LIB_SOURCES EXCLUDE REGEX ".*main\\.cpp$")

add_library(lattice_lib STATIC ${LIB_SOURCES})

# Create main executable
add_executable(lattice_triangulation "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")
target_link_libraries(lattice_triangulation 
    lattice_lib
    ${CGAL_LIBRARIES}
    ${Boost_LIBRARIES}
    alglib
)

# Add OpenMP if available
if(OpenMP_CXX_FOUND)
    target_link_libraries(lattice_triangulation OpenMP::OpenMP_CXX)
endif()

# Add test executable if tests directory exists
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests")
    file(GLOB TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp")
    if(TEST_SOURCES)
        add_executable(run_tests ${TEST_SOURCES})
        target_link_libraries(run_tests 
            lattice_lib
            ${CGAL_LIBRARIES}
            ${Boost_LIBRARIES}
            alglib
        )
        if(OpenMP_CXX_FOUND)
            target_link_libraries(run_tests OpenMP::OpenMP_CXX)
        endif()
    endif()
endif()