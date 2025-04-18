# Lattice Triangulation Project

This project offers comprehensive tools for 2D lattice generation and manipulation. It enables the creation of lattices with or without periodic boundary conditions and produces optimized triangular meshes through Delaunay triangulation. These high-quality meshes serve as the foundation for minimizing strain-energy functionals that possess GL(2,Z) invariance.

## Features

- Generate square and triangular lattices in 2D
- Create periodic copies of lattice points
- Generate triangular meshes using CGAL Delaunay triangulation
- Create finite element triangles with periodic boundary conditions
- Remeshing based on triangle angles
- Minimize strain-energy functionals using a highly efficient L-BFGS algorithm
- Includes a coarse-graining procedure to calculate strain-energy functionals from pair interatomic potentials
- Nano indentation  can now be applied and crystal orientation can be chosen 

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
## Example: 2D Square Crystal under simple shear load 
(for the strain-energy density, see Conti&Zanzotto,  Archive for Rational Mechanics and Analysis, Volume 173, pages 69–88, 2004, also Baggio et al., 
Phys. Rev. Lett. 123, 205501, 2019)

- Stress field of the  initial dislocated configuration

<img width="527" alt="Screenshot 2025-03-25 at 00 04 19" src="https://github.com/user-attachments/assets/e16c2de0-439b-498a-933d-20fa60365ce2" />

- Stress-field after 100% of deformation
  
<img width="615" alt="Screenshot 2025-03-25 at 00 00 33" src="https://github.com/user-attachments/assets/48c4a1a5-08f4-4198-b379-b7f132814216" />

- Stress-field on the mesh  after 100% of deformation

<img width="612" alt="Screenshot 2025-03-25 at 00 00 41" src="https://github.com/user-attachments/assets/54654d83-4e9f-492c-bc60-df7bc3c5cb7d" />

- Nano-indentation of a Lennard-Jones crystal

<img width="779" alt="Screenshot 2025-04-14 at 01 13 26" src="https://github.com/user-attachments/assets/7ca8c635-8c63-451d-8176-92c2ed6694d9" />
