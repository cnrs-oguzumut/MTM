cmake_minimum_required(VERSION 3.10)
project(LatticeTriangulation)

#-----------------------------------------------------------------------------
# Global Settings
#-----------------------------------------------------------------------------

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

#-----------------------------------------------------------------------------
# Compiler-Specific Optimizations
#-----------------------------------------------------------------------------

function(set_target_optimizations TARGET)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
        target_compile_options(${TARGET} PRIVATE 
            -O3
            -march=native
            -flto
            -ffast-math
            -funroll-loops
            -fstrict-aliasing
            -fvisibility=hidden
        )
        if(APPLE)
            target_link_options(${TARGET} PRIVATE -dead_strip)
        endif()
    elseif(MSVC)
        target_compile_options(${TARGET} PRIVATE /O2 /Ob2 /GL /fp:fast)
        target_link_options(${TARGET} PRIVATE /LTCG)
    endif()
endfunction()

# Global release flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-flto")
elseif(MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /Ob2 /DNDEBUG /GL")
    set(CMAKE_C_FLAGS_RELEASE "/O2 /Ob2 /DNDEBUG /GL")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "/LTCG")
endif()

#-----------------------------------------------------------------------------
# Dependencies
#-----------------------------------------------------------------------------

find_package(Eigen3 REQUIRED)
find_package(CGAL REQUIRED)
find_package(Boost REQUIRED)
find_package(OpenMP)

#-----------------------------------------------------------------------------
# ALGLIB Library
#-----------------------------------------------------------------------------

set(ALGLIB_DIR "/Users/usalman/programming/alglib-cpp")
include_directories(${ALGLIB_DIR})

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

add_library(alglib STATIC ${ALGLIB_SOURCES})
set_target_optimizations(alglib)
target_include_directories(alglib PUBLIC ${ALGLIB_DIR})

#-----------------------------------------------------------------------------
# Main Library
#-----------------------------------------------------------------------------

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${EIGEN3_INCLUDE_DIR}
    ${CGAL_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

file(GLOB_RECURSE LIB_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
list(FILTER LIB_SOURCES EXCLUDE REGEX ".*main\\.cpp$")

add_library(lattice_lib STATIC ${LIB_SOURCES})
set_target_optimizations(lattice_lib)
target_link_libraries(lattice_lib PUBLIC ${CGAL_LIBRARIES} ${Boost_LIBRARIES})

if(OpenMP_CXX_FOUND)
    target_link_libraries(lattice_lib PUBLIC OpenMP::OpenMP_CXX)
endif()

#-----------------------------------------------------------------------------
# Main Executable
#-----------------------------------------------------------------------------

add_executable(lattice_triangulation "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")
set_target_optimizations(lattice_triangulation)
target_link_libraries(lattice_triangulation PRIVATE lattice_lib alglib)

#-----------------------------------------------------------------------------
# Tests (if present)
#-----------------------------------------------------------------------------

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests")
    file(GLOB TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp")
    if(TEST_SOURCES)
        add_executable(run_tests ${TEST_SOURCES})
        set_target_optimizations(run_tests)
        target_link_libraries(run_tests PRIVATE lattice_lib alglib)
    endif()
endif()

#-----------------------------------------------------------------------------
# Build Information
#-----------------------------------------------------------------------------

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "Release flags: ${CMAKE_CXX_FLAGS_RELEASE}")