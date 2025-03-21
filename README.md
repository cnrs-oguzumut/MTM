# Lattice Triangulation Project

This project provides tools for generating 2D lattices, creating periodic boundary conditions, 
and generating triangular meshes using Delaunay triangulation.

## Features

- Generate square and triangular lattices in 2D
- Create periodic copies of lattice points
- Generate triangular meshes using CGAL Delaunay triangulation
- Create finite element triangles with periodic boundary conditions
- Optimize lattice configurations with energy minimization (using ALGLIB)
- Parallel computation support with OpenMP

## Dependencies

- Eigen3: For linear algebra operations
- CGAL: For computational geometry algorithms
- Boost: For various utilities
- ALGLIB: For optimization routines
- OpenMP: For parallel computation (optional)## Building

### Prerequisites

Make sure you have all dependencies installed on your system.

For ALGLIB, you need to:
1. Download ALGLIB from https://www.alglib.net/
2. Place it in a known location
3. Update the `ALGLIB_DIR` in CMakeLists.txt to point to your ALGLIB directory

```cmake
# In CMakeLists.txt
set(ALGLIB_DIR "/path/to/your/alglib-cpp")


```bash
mkdir build
cd build
cmake ..
make
```

## Running

```bash
./lattice_triangulation
```

## Running Tests

```bash
./run_tests
```
