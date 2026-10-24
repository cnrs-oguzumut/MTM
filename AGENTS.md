# Antigravity Rules & Engineering Conventions for Triangular Lattices

These rules document essential conventions and known pitfalls in this repository (`lattice_triangulation`) to prevent repeating bugs across simulation experiments.

---

## 1. Reference Element Leg Length ($h$) & Lagrange Reduction
- **Principle:** Lagrange reduction ($\mathrm{GL}(2, \mathbb{Z})$) operates strictly on the integer reference grid $\mathbb{Z}^2$. We ALWAYS evaluate atomistic pair potentials using `SquareLatticeCalculator` over $\mathbb{Z}^2$.
- **Triangular Mapping:** For triangular lattices, the optimal lattice parameter $a_0$ is mapped to $\mathbb{Z}^2$ via:
  $$\text{calc\_scale} = \frac{a_0}{(4/3)^{1/4}}.$$
- **Reference Triangle:** The reference triangle is ALWAYS a right isosceles triangle:
  $$p_1 = (0, 0), \quad p_2 = (h, 0), \quad p_3 = (0, h).$$
- **CRITICAL RULE (Do not use $h = 1.0$):**
  $h$ MUST match `calc_scale` whenever atomistic potentials (modes 1, 2, 4) are used:
  ```cpp
  double calc_scale = (potential_mode == 3)
      ? 1.0
      : ((lattice_type == "triangular") ? (optimal_lattice_parameter / std::pow(4.0 / 3.0, 0.25)) : optimal_lattice_parameter);
  double h = (potential_mode == 3) ? 1.0 : calc_scale;
  ```
  If $h = 1.0$ is used when nodal coordinates in physical space are scaled by $a_0$, the deformation gradient $\mathbf{F} = \mathbf{x}_e \cdot \nabla_{\mathbf{X}}\mathbf{N}$ retains an extra factor of $a_0$, compressing all interatomic bonds by $a_0$ ($8\%$ artificial hydrostatic compression for LJ, generating massive compressive stress $\sim -7.0$).

---

## 2. Unit-Cell Area Normalization
- **Rule:** NEVER normalize strain energy or stress derivatives by `pow(ideal_lattice_parameter, 2.0)` or hardcoded $a_0^2$.
- **Correct Formulation:** ALWAYS use `userData->calculator.getUnitCellArea()`.
  - For square lattice: `getUnitCellArea()` returns $a_0^2$.
  - For triangular lattice: `getUnitCellArea()` returns $\text{calc\_scale}^2 = \frac{a_0^2}{\sqrt{4/3}} = \frac{\sqrt{3}}{2} a_0^2$.
  This ensures energy density and Cauchy stress match between solver totals and nodal averages.

---

## 3. Acoustic Tensor Spatial Stiffness Ground State ($\mathbf{F}_0$)
- **Principle:** `SquareLatticeCalculator` operates relative to $\mathbb{Z}^2$. At $\mathbf{C} = \mathbf{I}$, a square lattice under triangular potentials exhibits negative shear modulus (shear instability).
- **Rule:** The ground-state deformation gradient mapping $\mathbb{Z}^2$ to equilateral triangles is:
  $$\mathbf{F}_0 = (4/3)^{1/4} \begin{pmatrix} 1 & 0.5 \\ 0 & \frac{\sqrt{3}}{2} \end{pmatrix}.$$
  When evaluating spatial tangent stiffness $\mathbb{C}^{\text{spatial}}(\bar{\mathbf{F}})$ via `AcousticTensor`, pass the total deformation relative to $\mathbb{Z}^2$:
  $$\mathbf{F}_{\text{tot}} = \bar{\mathbf{F}} \, \mathbf{F}_0, \quad \mathbf{C}_{\text{tot}} = \mathbf{F}_{\text{tot}}^T \mathbf{F}_{\text{tot}}.$$
  This guarantees $\mathbf{C}_{\text{tot}} = \mathbf{C}_{\text{ref}}$ at $\bar{\mathbf{F}} = \mathbf{I}$ and yields a strictly positive-definite tangent.

---

## 4. Adaptive Delaunay Mesh Edge Filtering (`max_edge`)
- **Principle:** In non-periodic meshes, Delaunay circumcircles can span across rows to connect second-nearest neighbors ($a_0 \sqrt{3} \approx 1.732 a_0$), creating artificial boundary slivers that collapse under displacement cuts.
- **Rule:** In `AdaptiveMesher::generateMesh` ([`src/Remesher.cpp`](file:///Users/usalman/programming/FEM_2D/factorized/lattice_triangulation/src/Remesher.cpp)), separate square and triangular thresholds:
  ```cpp
  bool is_triangular = (offsets[0] > 1e-6 && std::abs(offsets[1] / offsets[0] - std::sqrt(3.0)) < 0.25);
  double max_edge = is_triangular ? (median_nn * 1.35) : (median_nn * (std::sqrt(2.0) + 0.35));
  ```
  $1.35 a_0$ accepts equilateral edges ($1.0 a_0$) while cleanly rejecting boundary slivers ($1.732 a_0$).

---

## 5. Periodic Boundary Condition Rows ($n_y$ Parity)
- **Rule:** For periodic triangular lattices, $n_y$ MUST be an even integer so that row $n_y - 1$ (which has $0.5 a_0$ stagger) alternates seamlessly with row $0$ across the vertical boundary. If $n_y$ is odd, increment $n_y \gets n_y + 1$.

---

## 6. Checkpoint Restart Scaling
- **Rule:** When calling `ConfigurationSaver::loadCheckpoint`, calculate `dndx` using $h = \text{calc\_scale} = a_0 / (4/3)^{0.25}$ for triangular lattices, rather than hardcoded square lattice parameters.
