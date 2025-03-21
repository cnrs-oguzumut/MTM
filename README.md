# Lattice Triangulation Project

This project provides tools for generating 2D lattices, creating periodic boundary conditions, 
and generating triangular meshes using Delaunay triangulation.

## Features

- Generate square and triangular lattices in 2D
- Create periodic copies of lattice points
- Generate triangular meshes using CGAL Delaunay triangulation
- Create finite element triangles with periodic boundary conditions

## Dependencies

- Eigen3
- CGAL
- Boost

## Building

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
