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
- ITensor (for tensor computations in acoustic analysis)
- SuiteSparse/CHOLMOD (optional; faster sparse Cholesky for `--precond=stiffness`)

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

### Preconditioned L-BFGS

The energy relaxation can use the analytic FEM stiffness as an L-BFGS preconditioner
(sparse Cholesky + change of variables, `src/optimization/PreconditionedLBFGS.cpp`):

```bash
./lattice_triangulation 150 150 positive 42 --precond=stiffness --precond-from-step=1
```

- `--precond=stiffness|laplacian|diag|none` (default `none`: the original plain L-BFGS)
- `--precond-tol=1e-6` stop when max|dE/dx| < tol
- `--precond-refresh=N` reuse the stiffness factorization for N load steps (default 1)
- `--precond-from-step=N` keep the plain solver for load steps < N; `1` keeps the initial
  relaxation of the noisy lattice identical to plain-solver runs

After an instability the preconditioned solver can settle in a different metastable
state than plain L-BFGS, so compare statistics, not individual avalanches.

`build/benchmark_preconditioner` compares the solvers from identical start states
(`--validate` also checks the analytic stiffness against the ITensor assembler).

### Stability monitor (lowest eigenvalues of the stiffness)

The lowest eigenvalues of the Hessian K at the relaxed states can be computed during the
loading (`src/optimization/StiffnessSpectrum.cpp`: analytic K, Cholesky shift-invert
Lanczos, rigid translations projected out; ~0.15 s per state at 150x150):

```bash
./lattice_triangulation 150 150 positive 42 1 --precond=stiffness --precond-from-step=1 \
    --eig-every=5 --eig-at-avalanche=1
```

- `--eig-every=N` every N-th relaxed state (default 0 = off)
- `--eig-at-avalanche=1` also the last stable state before and the state after each avalanche
- `--eig-modes=5` eigenvalues per computation, `--eig-vectors=2` soft modes per avalanche
- `--eig-refine=0.2` every step while lambda_min < 0.2 x its value after the last avalanche
- `--eig-refine-ahead=2` every step while lambda_min^2, extrapolated linearly, reaches zero
  within 2 x N steps
- `--eig-retro=4` at each instability (stress jump or avalanche) also the last 4 relaxed
  states before it, kept in memory, so the approach to every instability is resolved

Output: `eigen_log.csv` and `eigen_modes/soft_modes_XXXXX.vtk` (XXXXX = id of the
PRE-avalanche configuration; translations are never written). Post-processing of saved
configurations (`analyze_data_from_folder`) uses the same solver, or the old ITensor path
with `--eig-solver=legacy`.

`python3 plot_saddle_node.py <run_dir>` finds every instability of a run (stress jumps,
saved or not), fits lambda_min^2 = s (alpha_c - alpha) on the last points before it and
plots lambda_min with the fitted square roots, the collapse on slope 1/2 and where the
predicted alpha_c falls (`saddle_node.png`, `saddle_node_fits.csv`).

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
