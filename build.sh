#!/bin/bash

# Create build directory
mkdir -p build
cd build

# Run CMake and build
cmake ..
make -j4

echo "Build complete. Run with:"
echo "./build/lattice_triangulation"
