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
## Example: 2D Square Crystal under simple shear load (for the strain-energy density, see Conti&Zanzotto 2024 )

- Stress-field of the  initial dislocated configuration

<img width="527" alt="Screenshot 2025-03-25 at 00 04 19" src="https://github.com/user-attachments/assets/e16c2de0-439b-498a-933d-20fa60365ce2" />

- Stress-field after 100% of deformation
  
<img width="615" alt="Screenshot 2025-03-25 at 00 00 33" src="https://github.com/user-attachments/assets/48c4a1a5-08f4-4198-b379-b7f132814216" />

- Stress-field on the mesh  after 100% of deformation

<img width="612" alt="Screenshot 2025-03-25 at 00 00 41" src="https://github.com/user-attachments/assets/54654d83-4e9f-492c-bc60-df7bc3c5cb7d" />
