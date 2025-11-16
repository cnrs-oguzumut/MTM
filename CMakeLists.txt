cmake_minimum_required(VERSION 3.10)
project(LatticeTriangulation)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set aggressive optimization flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    # For GCC, Clang, and AppleClang
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -flto=auto -ffast-math -funroll-loops -fstrict-aliasing -w")
    set(CMAKE_C_FLAGS_RELEASE   "-O3 -DNDEBUG -march=native -flto=auto -ffast-math -funroll-loops -fstrict-aliasing -w")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-flto=auto ")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    # For Microsoft Visual C++
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /Ob2 /DNDEBUG /GL")
    set(CMAKE_C_FLAGS_RELEASE "/O2 /Ob2 /DNDEBUG /GL")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /LTCG")
endif()

# Make sure to build in Release mode by default
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

# Find required packages
find_package(Eigen3 REQUIRED)
find_package(CGAL REQUIRED)
find_package(Boost REQUIRED)
find_package(OpenMP)

# Option to enable/disable ITensor support
option(USE_ITENSOR "Enable ITensor support" ON)

# ITensor detection (optional)
set(ITENSOR_FOUND FALSE)
set(ITENSOR_LIBRARIES "")

if(USE_ITENSOR)
    # Set ITensor directory - adjust path as needed
    set(ITENSOR_DIR "/Users/usalman/programming/itensor")
    set(ITENSOR_INCLUDE_DIR "${ITENSOR_DIR}")
    set(ITENSOR_LIB_DIR "${ITENSOR_DIR}/lib")
    
    # Check if ITensor directory exists
    if(EXISTS "${ITENSOR_DIR}" AND EXISTS "${ITENSOR_INCLUDE_DIR}/itensor")
        # Check if ITensor library exists
        if(EXISTS "${ITENSOR_LIB_DIR}/libitensor-g.a")
            set(ITENSOR_LIBRARY "${ITENSOR_LIB_DIR}/libitensor-g.a")
            set(ITENSOR_FOUND TRUE)
            message(STATUS "Found ITensor debug library: ${ITENSOR_LIBRARY}")
        elseif(EXISTS "${ITENSOR_LIB_DIR}/libitensor.a")
            set(ITENSOR_LIBRARY "${ITENSOR_LIB_DIR}/libitensor.a")
            set(ITENSOR_FOUND TRUE)
            message(STATUS "Found ITensor release library: ${ITENSOR_LIBRARY}")
        else()
            message(WARNING "ITensor library not found in ${ITENSOR_LIB_DIR}")
            message(WARNING "ITensor support will be disabled")
        endif()
    else()
        message(WARNING "ITensor directory not found at ${ITENSOR_DIR}")
        message(WARNING "ITensor support will be disabled")
    endif()
    
    # Try to find BLAS/LAPACK for ITensor if ITensor was found
    if(ITENSOR_FOUND)
        find_package(BLAS)
        find_package(LAPACK)
        
        # Set up ITensor linking libraries
        set(ITENSOR_LIBRARIES ${ITENSOR_LIBRARY})
        
        # Add BLAS/LAPACK libraries if found (ITensor dependency)
        if(BLAS_FOUND AND LAPACK_FOUND)
            list(APPEND ITENSOR_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
            message(STATUS "Using BLAS/LAPACK for ITensor")
        else()
            # Try to find common BLAS implementations
            find_library(OPENBLAS_LIB openblas)
            find_library(CBLAS_LIB cblas)
            find_library(BLAS_LIB blas)
            
            if(OPENBLAS_LIB)
                list(APPEND ITENSOR_LIBRARIES ${OPENBLAS_LIB})
                message(STATUS "Using OpenBLAS for ITensor: ${OPENBLAS_LIB}")
            elseif(CBLAS_LIB AND BLAS_LIB)
                list(APPEND ITENSOR_LIBRARIES ${CBLAS_LIB} ${BLAS_LIB})
                message(STATUS "Using CBLAS/BLAS for ITensor")
            else()
                message(WARNING "No BLAS/LAPACK found - ITensor may not link properly")
                set(ITENSOR_FOUND FALSE)
                message(WARNING "ITensor support disabled due to missing BLAS/LAPACK")
            endif()
        endif()
    endif()
else()
    message(STATUS "ITensor support manually disabled")
endif()

# Option to enable/disable Spectra support
option(USE_SPECTRA "Enable Spectra eigenvalue solver support" ON)

# Spectra detection (optional, header-only library)
set(SPECTRA_FOUND FALSE)

if(USE_SPECTRA)
    # Set Spectra directory - adjust path as needed
    set(SPECTRA_DIR "/Users/usalman/programming/FEM_2D/factorized/lattice_triangulation/spectra")
    set(SPECTRA_INCLUDE_DIR "${SPECTRA_DIR}/include")
    
    # Check if Spectra directory exists
    if(EXISTS "${SPECTRA_INCLUDE_DIR}" AND EXISTS "${SPECTRA_INCLUDE_DIR}/Spectra")
        set(SPECTRA_FOUND TRUE)
        message(STATUS "Found Spectra at: ${SPECTRA_INCLUDE_DIR}")
    elseif(EXISTS "${SPECTRA_DIR}/include/Spectra")
        set(SPECTRA_INCLUDE_DIR "${SPECTRA_DIR}/include")
        set(SPECTRA_FOUND TRUE)
        message(STATUS "Found Spectra at: ${SPECTRA_INCLUDE_DIR}")
    elseif(EXISTS "${SPECTRA_DIR}/Spectra")
        set(SPECTRA_INCLUDE_DIR "${SPECTRA_DIR}")
        set(SPECTRA_FOUND TRUE)
        message(STATUS "Found Spectra at: ${SPECTRA_INCLUDE_DIR}")
    else()
        message(WARNING "Spectra directory not found at ${SPECTRA_DIR}")
        message(WARNING "Spectra support will be disabled")
    endif()
else()
    message(STATUS "Spectra support manually disabled")
endif()

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

# Create ALGLIB library with optimization flags
add_library(alglib STATIC ${ALGLIB_SOURCES})
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(alglib PRIVATE -O3 -march=native -flto -ffast-math)
endif()

# Include directories
set(INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${EIGEN3_INCLUDE_DIR}
    ${CGAL_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

# Add ITensor include directory if found
if(ITENSOR_FOUND)
    list(APPEND INCLUDE_DIRS ${ITENSOR_INCLUDE_DIR})
endif()

# Add Spectra include directory if found
if(SPECTRA_FOUND)
    list(APPEND INCLUDE_DIRS ${SPECTRA_INCLUDE_DIR})
endif()

include_directories(${INCLUDE_DIRS})

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

# Add compiler flags and definitions
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(lattice_lib PRIVATE -O3 -march=native -flto -ffast-math)
endif()

# Add ITensor preprocessor definition if found
if(ITENSOR_FOUND)
    target_compile_definitions(lattice_lib PRIVATE USE_ITENSOR=1)
    message(STATUS "ITensor support ENABLED - USE_ITENSOR=1 defined")
else()
    message(STATUS "ITensor support DISABLED")
endif()

# Add Spectra preprocessor definition if found
if(SPECTRA_FOUND)
    target_compile_definitions(lattice_lib PRIVATE USE_SPECTRA=1)
    message(STATUS "Spectra support ENABLED - USE_SPECTRA=1 defined")
else()
    message(STATUS "Spectra support DISABLED")
endif()

# Create main executable
add_executable(lattice_triangulation "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(lattice_triangulation PRIVATE -O3 -march=native -flto -ffast-math)
endif()

# Add ITensor definition to main executable too
if(ITENSOR_FOUND)
    target_compile_definitions(lattice_triangulation PRIVATE USE_ITENSOR=1)
endif()

# Add Spectra definition to main executable too
if(SPECTRA_FOUND)
    target_compile_definitions(lattice_triangulation PRIVATE USE_SPECTRA=1)
endif()

# Set up base libraries for linking
set(BASE_LIBRARIES
    lattice_lib
    ${CGAL_LIBRARIES}
    ${Boost_LIBRARIES}
    alglib
)

# Add ITensor libraries if found
if(ITENSOR_FOUND)
    list(APPEND BASE_LIBRARIES ${ITENSOR_LIBRARIES})
endif()

# Note: Spectra is header-only, so no libraries to link

# Link libraries to main executable
target_link_libraries(lattice_triangulation ${BASE_LIBRARIES})

# Add OpenMP if available
if(OpenMP_CXX_FOUND)
    target_link_libraries(lattice_triangulation OpenMP::OpenMP_CXX)
endif()

# Add test executable if tests directory exists
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests")
    file(GLOB TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp")
    if(TEST_SOURCES)
        add_executable(run_tests ${TEST_SOURCES})
        if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(run_tests PRIVATE -O3 -march=native -flto -ffast-math)
        endif()
        
        # Add definitions to tests
        if(ITENSOR_FOUND)
            target_compile_definitions(run_tests PRIVATE USE_ITENSOR=1)
        endif()
        if(SPECTRA_FOUND)
            target_compile_definitions(run_tests PRIVATE USE_SPECTRA=1)
        endif()
        
        target_link_libraries(run_tests
            lattice_lib
            ${ITENSOR_LIBRARIES}
            ${CGAL_LIBRARIES}
            ${Boost_LIBRARIES}
            alglib
        )
        if(OpenMP_CXX_FOUND)
            target_link_libraries(run_tests OpenMP::OpenMP_CXX)
        endif()
    endif()
endif()

# Print build information
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "Release flags: ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "ITensor directory: ${ITENSOR_DIR}")
message(STATUS "ITensor library: ${ITENSOR_LIBRARY}")
message(STATUS "Spectra directory: ${SPECTRA_DIR}")
if(SPECTRA_FOUND)
    message(STATUS "Spectra include: ${SPECTRA_INCLUDE_DIR}")
endif()