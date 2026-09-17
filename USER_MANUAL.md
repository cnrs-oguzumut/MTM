# Lattice Triangulation 2D - User Manual & Simulation Guide

This manual covers the installation, compilation, command-line usage, high-performance solver options, bitwise-exact checkpoint restart system, output formats, and post-processing tools for the 2D Lattice Triangulation FEM framework.

---

## Table of Contents

1. [Architecture & Overview](#1-architecture--overview)
2. [Compilation & Building](#2-compilation--building)
   - [Local Build (CMake)](#local-build-cmake)
   - [Cluster Build (Makefile.cluster with METIS)](#cluster-build-makefilecluster-with-metis)
3. [Command-Line Interface (CLI) Reference](#3-command-line-interface-cli-reference)
   - [Positional Arguments](#positional-arguments)
   - [Solver & Preconditioning Options](#solver--preconditioning-options)
   - [Stability Monitor & Soft-Mode Options](#stability-monitor--soft-mode-options)
   - [Checkpoint & Restart Options](#checkpoint--restart-options)
4. [Checkpoint & Restart System](#4-checkpoint--restart-system)
   - [How Checkpointing Works](#how-checkpointing-works)
   - [Bitwise IEEE-754 Precision (`std::hexfloat`)](#bitwise-ieee-754-precision-stdhexfloat)
   - [Resuming Interrupted Runs](#resuming-interrupted-runs)
   - [Disk Space Savings (Elimination of `triangle_data/`)](#disk-space-savings-elimination-of-triangle_data)
5. [Preconditioned L-BFGS & METIS Acceleration](#5-preconditioned-l-bfgs--metis-acceleration)
6. [Simulation Outputs & File Layout](#6-simulation-outputs--file-layout)
7. [Post-Processing & Visualization](#7-post-processing--visualization)
8. [HPC & Cluster Workflow (SLURM)](#8-hpc--cluster-workflow-slurm)
9. [Dislocation Studies & Multi-Shift Core Analysis](#9-dislocation-studies--multi-shift-core-analysis)
   - [Overview & Physical Context](#overview--physical-context)
   - [Running & Re-running Simulations](#running--re-running-simulations)
   - [Core Center Tracking: Peak vs. Centroid](#core-center-tracking-peak-vs-centroid)
   - [Core Broadening, FWHM & Adaptive Core Radius](#core-broadening-fwhm--adaptive-core-radius)
   - [Aligned Core Energy Profiles $E(x - x_c)$](#aligned-core-energy-profiles-ex---xc)
   - [Automated Analysis & Plotting Script](#automated-analysis--plotting-script)

---

## 1. Architecture & Overview

The code simulates non-convex crystal lattices under simple shear loading using continuum nonlinear elasticity with GL(2,Z) modular invariance (Conti & Zanzotto 2004, Baggio et al. 2019). The simulation captures:
- Elastic shear deformation and lattice rotations.
- Nonlinear bifurcation, dislocation nucleation, and plastic slip avalanches.
- Dynamic adaptive remeshing upon topological connectivity changes.
- Preconditioned energy minimization (Stiffness / Laplacian Cholesky preconditioning).
- Real-time stability tracking via shift-invert Lanczos spectral analysis.

---

## 2. Compilation & Building

### Local Build (CMake)
Recommended on macOS or personal workstations:
```bash
mkdir -p build && cd build
cmake ..
make -j4
```
Dependencies:
- C++17 compiler (`clang++` or `g++`)
- `Eigen3`
- `CGAL`
- `Boost`
- `SuiteSparse` (optional, for CHOLMOD)
- `ITensor`

### Cluster Build (Makefile.cluster with METIS)
Recommended for high-performance computing clusters (e.g. `Magi`, `Tchenla`, SLURM clusters):
```bash
# Clean previous objects and build with system SuiteSparse + METIS
make -f Makefile.cluster clean
make -f Makefile.cluster -j4
```

> **METIS Support**: `Makefile.cluster` is configured to detect system `libsuitesparse-dev` at `/usr/include/suitesparse` and `/usr/lib/x86_64-linux-gnu`. When compiled this way, CHOLMOD automatically uses **METIS nested dissection** ordering, reducing Cholesky factorization non-zeros by 25–35% and accelerating L-BFGS back-solves by 25–35% on large systems ($600 \times 600$ to $1000 \times 1000$).

---

## 3. Command-Line Interface (CLI) Reference

### Basic Execution
```bash
./lattice_triangulation <nx> <ny> [mode] [seed] [enable_remeshing] [options...]
```

### Positional Arguments
| Argument | Type | Default | Description |
| :--- | :---: | :---: | :--- |
| `nx` | Integer | `150` | Number of lattice unit cells along the $x$-axis. |
| `ny` | Integer | `150` | Number of lattice unit cells along the $y$-axis. |
| `mode` | String | `positive` | Loading direction: `positive` ($\alpha > 0$) or `negative` ($\alpha < 0$). |
| `seed` | Integer | `50` | Deterministic random seed for initial perturbation noise. |
| `enable_remeshing` | Boolean/Int | `1` | `1` / `true` enables adaptive topological remeshing; `0` / `false` disables it. |

---

### Solver & Preconditioning Options

| Option | Values | Default | Description |
| :--- | :---: | :---: | :--- |
| `--precond=<type>` | `stiffness`, `laplacian`, `diag`, `none` | `none` | Type of L-BFGS preconditioner. Use `stiffness` for 3x–10x faster convergence on large meshes. |
| `--precond-tol=<tol>` | Float | `1e-6` | Convergence gradient norm tolerance ($\max \|\nabla E\|_\infty$). |
| `--precond-refresh=<N>` | Integer | `1` | Recompute and factorize the preconditioner stiffness matrix every $N$ load steps. |
| `--precond-from-step=<N>` | Integer | `0` | Keep plain unconditioned L-BFGS for load steps $< N$. Useful to keep initial relaxation identical across runs. |

---

### Stability Monitor & Soft-Mode Options

Tracks the lowest eigenvalues ($\lambda_{\min}$) of the analytic stiffness matrix $K$ at relaxed states via shift-invert Lanczos:

| Option | Type | Default | Description |
| :--- | :---: | :---: | :--- |
| `--eig-every=<N>` | Integer | `0` (off) | Compute lowest eigenvalues every $N$ load steps. |
| `--eig-at-avalanche=1` | Flag | `0` (off) | Compute lowest eigenvalues before and after each avalanche event. |
| `--eig-modes=<M>` | Integer | `5` | Number of lowest eigenvalues to compute. |
| `--eig-vectors=<V>` | Integer | `2` | Number of eigenvector soft modes saved to VTK. |
| `--eig-refine=<frac>` | Float | `0.0` | Refine step rate when $\lambda_{\min} < \text{frac} \times \lambda_{\text{post-avalanche}}$. |
| `--eig-refine-ahead=<N>`| Integer | `0` | Refine step rate when linearly extrapolated $\lambda_{\min}^2 \to 0$ within $N \times \text{every}$ steps. |
| `--eig-retro=<K>` | Integer | `0` | Keep the last $K$ relaxed states in memory before an instability and evaluate eigenvalues retroactively. |

---

### Checkpoint & Restart Options

| Option | Syntax | Description |
| :--- | :--- | :--- |
| `--restart=<file>` | `--restart=checkpoints/latest.chk`<br>`--restart checkpoints/checkpoint_00015.chk` | Directly resumes an interrupted simulation from the specified single-file checkpoint. |

---

## 4. Checkpoint & Restart System

### How Checkpointing Works
During continuous loading, the engine automatically saves checkpoints inside the `checkpoints/` directory:
1. **Avalanche Checkpoints**:
   - Every time a stress drop or avalanche is detected, the engine saves the exact post-avalanche state to `checkpoints/checkpoint_XXXXX.chk` (where `XXXXX` matches the VTK file ID).
2. **Periodic Checkpoints**:
   - During long elastic loading segments, the engine automatically saves a periodic checkpoint every 50 steps: `checkpoints/checkpoint_step_XXXXX.chk`.
3. **Latest Pointer**:
   - `checkpoints/latest.chk` is automatically updated on every checkpoint save.

### Bitwise IEEE-754 Precision (`std::hexfloat`)
In chaotic nonlinear systems with dislocation nucleation, decimal truncation introduces slight roundoff errors ($\sim 10^{-10}$) that grow exponentially near bifurcation points and cause restarted runs to diverge over time.

Our checkpoint system uses C++17 `std::hexfloat` to write the **exact 64-bit IEEE-754 binary representation** of:
- Atomic / nodal coordinates $(x, y)$
- Boundary periodic translations $(t_x, t_y)$
- Reference triangle areas $A_0$
- External deformation gradient $F_{\text{ext}}$

> **Zero Trajectory Drift**: Verification across 100+ steps confirms that restarting from a `.chk` file produces an energy and stress trajectory matching the continuous run down to **machine precision** ($\Delta E < 2.75 \times 10^{-14}$, $\Delta \sigma \sim 10^{-16}$).

### Resuming Interrupted Runs
To resume a simulation that was stopped or preempted by a cluster walltime limit:
```bash
./lattice_triangulation --restart=checkpoints/latest.chk --precond=stiffness
```
All system parameters (`nx`, `ny`, `mode`, `seed`, `current_alpha`, `step_size`, `alpha_end`, `enable_remeshing`, mesh elements, and nodal coordinates) are restored automatically from the checkpoint header.

### Disk Space Savings (Elimination of `triangle_data/`)
Previously, simulations wrote 30-column ASCII files (`triangle_data/triangles_XXXXX.dat`) consuming ~360 MB per step at $600 \times 600$. 
- These files are now **completely eliminated**.
- Deformation gradients and metric tensors are reconstructed in memory on-the-fly.
- A single-file checkpoint is only **~150 KB** at $30 \times 30$ and **~35 MB** at $600 \times 600$, reducing overall disk usage by **> 85%**.

---

## 5. Preconditioned L-BFGS & METIS Acceleration

For systems larger than $100 \times 100$, standard unconditioned L-BFGS suffers from severe ill-conditioning due to acoustic phonon modes (eigenvalues scale as $(h/L)^2 \to 0$).

### Enabling Stiffness Preconditioning
```bash
./lattice_triangulation 600 600 positive 50 1 --precond=stiffness
```
- **Stiffness Preconditioner**: Solves $P z = r$ using the sparse analytic tangent stiffness matrix assembled from the 2D finite element triangles.
- **Sparse Cholesky**: Handled by SuiteSparse `CHOLMOD`.
- **METIS Ordering**: With `Makefile.cluster`, CHOLMOD evaluates METIS nested dissection, minimizing matrix fill-in during Cholesky decomposition.

### Performance Gains
| Grid Size | Unconditioned L-BFGS Iterations | Preconditioned L-BFGS Iterations | Total Speedup |
| :---: | :---: | :---: | :---: |
| $150 \times 150$ | 800 – 1,200 | 40 – 80 | **~4x – 6x** |
| $600 \times 600$ | 5,000 – 15,000 | 80 – 150 | **~10x – 25x** |

---

## 6. Simulation Outputs & File Layout

When running a simulation, the output directory contains:

```
├── energy_stress_log.csv          # High-precision CSV of load, energy, stress, and area
├── fractional_drops.csv           # Detailed energy and stress drop diagnostic metrics
├── simulation.log                 # Console output redirect
├── checkpoints/
│   ├── latest.chk                 # Symlink/copy of the most recent valid checkpoint
│   ├── checkpoint_step_00050.chk  # Periodic elastic checkpoints
│   └── checkpoint_00001.chk       # Pre- and post-avalanche state checkpoints
├── vtk_output/
│   ├── configuration_00000.vtk    # Initial configuration (Binary legacy VTK, 32-bit floats)
│   ├── configuration_00001.vtk    # Pre-avalanche state (Binary legacy VTK, 32-bit floats)
│   └── configuration_00002.vtk    # Post-avalanche state (Binary legacy VTK, 32-bit floats)
├── eigen_log.csv                  # (Optional) Lowest stiffness eigenvalues
└── eigen_modes/                   # (Optional) Soft-mode displacement VTK fields
```

> **Binary VTK Storage**: `.vtk` files are stored in binary big-endian legacy format using 32-bit floats for visualization fields (nodal energies, Cauchy stresses, valence, and reference coordination). This reduces individual file size by **~65–70%** and speeds up I/O by **5x–10x**, while remaining 100% compatible with ParaView and `plot_vtk_output.py`.

### `energy_stress_log.csv` Columns
1. `Iteration`: Global step index.
2. `Alpha`: Applied shear strain $\alpha$.
3. `PreEnergy`: Total system energy before L-BFGS minimization.
4. `PreStress`: Average shear stress $\sigma_{xy}$ before minimization.
5. `PostEnergy`: Total system energy after minimization.
6. `PostStress`: Average shear stress $\sigma_{xy}$ after minimization.
7. `EnergyChange`: Energy change during relaxation ($-\Delta E$).
8. `StressChange`: Stress change during relaxation ($-\Delta \sigma$).
9. `PreArea`: Total domain area before remeshing.
10. `PostArea`: Total domain area after remeshing.
11. `shouldRemesh`: Boolean indicator whether topological remeshing was triggered.

---

## 7. Post-Processing & Visualization

### Plotting Stress-Strain & Energy Curves
```bash
python3 plot_vtk_output.py
```
Or parse `energy_stress_log.csv` using Python:
```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("energy_stress_log.csv")
plt.plot(df['Alpha'], df['PostStress'], label=r'$\sigma_{xy}$')
plt.xlabel(r'Shear Strain $\alpha$')
plt.ylabel(r'Shear Stress $\sigma_{xy}$')
plt.title('Stress-Strain Response')
plt.grid(True)
plt.show()
```

### Visualizing Dislocation Microstructures in ParaView
1. Open ParaView.
2. Load `vtk_output/configuration_*.vtk`.
3. Apply the **Warp By Vector** filter using nodal displacements.
4. Color by `coordination` (highlights dislocation cores) or `stress_tensor_xy`.

---

## 8. HPC & Cluster Workflow (SLURM)

### Example SLURM Batch Script (`run_cluster.sbatch`)
```bash
#!/bin/bash
#SBATCH --job-name=fem_zanzotto_600
#SBATCH --output=simulation_%j.log
#SBATCH --error=simulation_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=32G

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# 1. Check if a valid checkpoint exists to resume
if [ -f "checkpoints/latest.chk" ]; then
    echo "Resuming from existing checkpoint: checkpoints/latest.chk"
    ./lattice_triangulation --restart=checkpoints/latest.chk --precond=stiffness
else
    echo "Starting fresh simulation..."
    ./lattice_triangulation 600 600 positive 50 1 --precond=stiffness
fi
```

### Resubmitting Chained Jobs
If your cluster enforces an 8-hour or 24-hour walltime limit:
```bash
sbatch run_cluster.sbatch
```
When the time limit is reached, simply resubmit the same script. It will automatically detect `checkpoints/latest.chk` and resume seamless, bitwise-exact integration from the exact point of interruption.

---

## 9. Dislocation Studies & Multi-Shift Core Analysis

### Overview & Physical Context
The framework includes dedicated modules to investigate isolated edge dislocation cores, their non-linear core energies, and their response to applied simple shear loading:
- **Lattice Setup**: Square $100 \times 100$ crystal ($10,000$ atoms).
- **Deformation Protocol**:
  1. Upper half ($y > y_{\rm mid}$) shifted horizontally in $x$ by $s \times h$ ($s \in \{0, 1, 2, 3, 4, 5\}$).
  2. Analytical Volterra edge dislocation field installed at $(x_0, y_0) = (49.5, 49.5)$ with Burgers vector $\mathbf{b} = (1.0, 0)$.
  3. Fixed outer frame boundary conditions.
  4. Energy relaxed to equilibrium via stiffness-preconditioned L-BFGS.
  5. Connectivity remains fixed on the pristine crystal (no remeshing).

### Running & Re-running Simulations
Simulations can be executed or re-run directly with a single command:
```bash
# Run multi-shift dislocation study across all shifts (s = 0, 1, 2, 3, 4, 5)
./build/lattice_triangulation 100 100 shifted_dislocation
```
This automatically:
- Creates the study folder `shifted_dislocation_study/`.
- Computes states for shifts $s = 0, 1, 2, 3, 4, 5$ into individual subdirectories `shift_0/` through `shift_5/`.
- Exports VTK meshes (`configuration_00000.vtk` unrelaxed, `configuration_00001.vtk` relaxed, and `defects_00001.vtk`).
- Outputs the master data table `shifted_dislocation_summary.csv`.

---

### Core Center Tracking: Peak vs. Centroid
Under applied upper crystal shear, the dislocation core moves (glides) along the slip plane. Therefore, measuring core energies centered at the fixed initial position $(49.5, 49.5)$ produces inaccurate off-center results.

Two primary metrics are used to locate the core center $x_c$:
1. **Discrete Energy Peak on Slip Plane ($x_{\rm peak}$)**:
   $$x_{\rm peak} = \arg\max_x E(x, y_{\rm slip})$$
   Identifies the single node with the highest strain energy density.
2. **Core Energy Centroid / Center of Mass ($x_{\rm cm}$)**:
   $$x_{\rm cm} = \frac{\sum_{i \in \text{core}} x_i \cdot E(x_i)}{\sum_{i \in \text{core}} E(x_i)} \quad \text{for nodes with } E(x_i) \ge \frac{1}{2} E_{\max}$$

> **Why the Centroid is Superior**: At higher shifts ($s \ge 3$), the top half is shifted by multiple lattice spacings relative to the bottom half. The energy density flattens into a plateau with distinct shoulders on the top and bottom atomic layers. The centroid $x_{\rm cm}$ cleanly balances these contributions, providing a stable, symmetric reference for the true core position.

| Shift $s$ | Upper Shift $u_x$ | Peak Slip $x_{\rm peak}$ | Centroid $x_{\rm cm}$ | Full Mesh Peak $(x, y)$ | Core Width $\text{FWHM}$ | Adaptive $R_{\rm core}(s)$ |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **$0$** | $0.0h$ | $49.76h$ | **$49.3h$** | $(49.76, 48.77)$ | **$5.3h$** | $3.64h$ |
| **$1$** | $1.0h$ | $49.73h$ | **$50.2h$** | $(49.73, 48.78)$ | **$5.3h$** | $3.64h$ |
| **$2$** | $2.0h$ | $50.76h$ | **$50.4h$** | $(49.69, 48.79)$ | **$5.3h$** | $3.66h$ |
| **$3$** | $3.0h$ | $51.80h$ | **$51.3h$** | $(53.13, 49.77)$ | **$7.4h$** | $4.70h$ |
| **$4$** | $4.0h$ | $53.89h$ | **$51.3h$** | $(54.08, 49.76)$ | **$9.4h$** | $5.72h$ |
| **$5$** | $5.0h$ | $55.89h$ | **$52.8h$** | $(56.10, 49.76)$ | **$10.5h$** | $6.23h$ |

---

### Core Broadening, FWHM & Adaptive Core Radius
As shear is applied, the dislocation core delocalizes across the slip plane:
- **Core Width Broadening**: The Full Width at Half Maximum ($\text{FWHM}$) expands from **$5.3h$** at $s = 0$ to **$10.5h$** at $s = 5$ (nearly doubling in width).
- **Peak Softening**: The maximum nodal energy density drops from $0.0636$ to $0.0399$ as mismatch strain distributes over a larger plastic zone.
- **Adaptive Core Radius**: To prevent truncating the broadening core, an adaptive integration cutoff radius is used:
  $$R_{\rm core}(s) = \max\left(3.0,\; \frac{\text{FWHM}(s)}{2} + 1.0\right)h$$
  This expands from $3.64h$ to $6.23h$, correctly bounding the nonlinear core across all deformation levels.

---

### Aligned Core Energy Profiles $E(x - x_c)$
To directly compare core structures across shifts, profiles along the slip plane are shifted by their core center:
$$x_{\rm rel} = x - x_c$$
This superimposes all 6 dislocation profiles at $x = 0$:
1. **Absolute Profile $E(x - x_c)$**: Clearly demonstrates peak drop and symmetrical tail broadening.
2. **Normalized Profile $E(x - x_c) / E_{\max}$**: Normalizes all peaks to $1.0$, giving a direct geometric visualization of the $\text{FWHM}$ core delocalization under shear.

---

### Automated Analysis & Plotting Script
Run the automated post-processing script from the project root:
```bash
python3 shifted_dislocation_study/plot_shifted_dislocation_study.py
```
This script automatically computes the centroid positions, widths, and adaptive radii, generating:
- `shifted_dislocation_adaptive_core_summary.csv`: Tabulated values of $E_{\rm tot}$, $x_{\rm peak}$, $x_{\rm cm}$, $\text{FWHM}$, $R_{\rm core}$, and core energies.
- `shifted_dislocation_energy_profile_aligned.png`: Side-by-side aligned $E(x - x_c)$ and normalized shape comparisons.
- `shifted_dislocation_energy_profile_x.png`: Slip-plane energy profiles showing horizontal core gliding.
- `shifted_dislocation_core_energy_vs_R.png`: Cumulative radial energy $E(R)$ centered at each core's true position.
- `shifted_dislocation_energy_vs_shift.png`: 3-panel figure tracking core glide $x_c(s)$, core width $\text{FWHM}(s)$, and total vs. adaptive core energy.

