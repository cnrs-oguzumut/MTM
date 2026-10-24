# Engineering Guide: Conventions & Pitfalls for Triangular Lattice Simulations

This document details the mathematical formulation, kinematic conventions, and common pitfalls encountered when implementing triangular lattice simulations with atomistic potentials and Lagrange reduction in this finite-element framework.

---

## 1. The Reference Configuration & Lagrange Reduction over $\mathbb{Z}^2$

### Mathematical Background
In our nonlinear finite-element framework:
1. **Integer Reference Grid ($\mathbb{Z}^2$):** Lagrange reduction operates via integer unimodular transformations $\mathbf{m} \in \mathrm{GL}(2, \mathbb{Z})$ on the metric tensor $\mathbf{C} = \mathbf{F}^T \mathbf{F}$ relative to the standard integer lattice $\mathbb{Z}^2$.
2. **Equilateral Ground State:** An equilateral triangular lattice has primitive basis vectors:
   $$\mathbf{a}_1 = a_0 \begin{pmatrix} 1 \\ 0 \end{pmatrix}, \quad \mathbf{a}_2 = a_0 \begin{pmatrix} 1/2 \\ \sqrt{3}/2 \end{pmatrix}.$$
   Under the mapping $\mathbf{F}_0$ from $\mathbb{Z}^2$:
   $$\mathbf{F}_0 = (4/3)^{1/4} \begin{pmatrix} 1 & 1/2 \\ 0 & \sqrt{3}/2 \end{pmatrix}.$$
   The reference metric tensor in $\mathbb{Z}^2$ is:
   $$\mathbf{C}_{\text{ref}} = \mathbf{F}_0^T \mathbf{F}_0 = \begin{pmatrix} \sqrt{4/3} & \frac{1}{2}\sqrt{4/3} \\ \frac{1}{2}\sqrt{4/3} & \sqrt{4/3} \end{pmatrix}.$$

### Potential Scaling Factor ($\text{calc\_scale}$)
In `SquareLatticeCalculator`, interatomic distances are evaluated as:
$$r = \text{scal} \cdot \sqrt{\mathbf{p}^T \mathbf{C} \mathbf{p}}, \quad \mathbf{p} \in \mathbb{Z}^2.$$
For the nearest neighbor $\mathbf{p} = (1, 0)^T$, $\mathbf{p}^T \mathbf{C}_{\text{ref}} \mathbf{p} = \sqrt{4/3} = (4/3)^{1/2}$.
To ensure that $r$ equals the true physical equilibrium lattice parameter $a_0$:
$$r = \text{scal} \cdot (4/3)^{1/4} = a_0 \implies \text{calc\_scale} = \frac{a_0}{(4/3)^{1/4}}.$$

---

## 2. Pitfall 1: Double-Scaling of Reference Element Legs $h$

### The Mechanism
The local deformation gradient inside element $e$ is:
$$\mathbf{F}^{(e)} = \mathbf{x}_e \cdot \nabla_{\mathbf{X}} \mathbf{N}^{(e)}.$$
The reference triangle is ALWAYS a right isosceles triangle with legs $h$:
$$p_1 = (0, 0), \quad p_2 = (h, 0), \quad p_3 = (0, h).$$
Its shape function derivatives are:
$$\nabla_{\mathbf{X}} \mathbf{N} = \begin{pmatrix} -1/h & -1/h \\ 1/h & 0 \\ 0 & 1/h \end{pmatrix}.$$
In physical space, the atomic coordinates $\mathbf{x}_e$ for an undeformed equilateral triangle have edge length $a_0$:
$$\mathbf{x}_e = \begin{pmatrix} 0 & a_0 & a_0/2 \\ 0 & 0 & a_0 \sqrt{3}/2 \end{pmatrix}.$$
Thus:
$$\mathbf{F}^{(e)} = \mathbf{x}_e \cdot \nabla_{\mathbf{X}} \mathbf{N} = \frac{a_0}{h} \begin{pmatrix} 1 & 1/2 \\ 0 & \sqrt{3}/2 \end{pmatrix}.$$

### The Bug
If $h = 1.0$ is used instead of $h = \text{calc\_scale}$:
$$\frac{a_0}{h} = a_0 \implies \mathbf{F}^{(e)} = a_0 \begin{pmatrix} 1 & 1/2 \\ 0 & \sqrt{3}/2 \end{pmatrix} = a_0 (4/3)^{-1/4} \mathbf{F}_0.$$
Then:
$$\mathbf{C}^{(e)} = a_0^2 (4/3)^{-1/2} \mathbf{C}_{\text{ref}}.$$
When passed into `SquareLatticeCalculator(calc_scale)`:
$$r = \text{calc\_scale} \cdot \sqrt{\mathbf{p}^T \mathbf{C} \mathbf{p}} = \frac{a_0}{(4/3)^{1/4}} \cdot \left[ a_0 (4/3)^{-1/4} (4/3)^{1/4} \right] = a_0^2 \approx 0.9964^2 \approx 0.924 < a_0.$$
This introduces an artificial $8\%$ hydrostatic compression on all bonds, generating:
- Huge initial compressive stress ($\sigma_{xx} = -6.9, \sigma_{yy} = -6.9$).
- Severe boundary deformation and shear.

### The Correct Implementation
```cpp
double calc_scale = (potential_mode == 3)
    ? 1.0
    : ((lattice_type == "triangular") ? (optimal_lattice_parameter / std::pow(4.0 / 3.0, 0.25)) : optimal_lattice_parameter);

double h = (potential_mode == 3) ? 1.0 : calc_scale;

Eigen::Vector2d p1(0.0, 0.0);
Eigen::Vector2d p2(h, 0.0);
Eigen::Vector2d p3(0.0, h);
Eigen::Matrix<double, 3, 2> dndx = calculateShapeDerivatives(p1, p2, p3);
```

With $h = \text{calc\_scale}$:
$$\frac{a_0}{h} = \frac{a_0}{a_0 / (4/3)^{1/4}} = (4/3)^{1/4} \implies \mathbf{F}^{(e)} = \mathbf{F}_0, \quad \mathbf{C}^{(e)} = \mathbf{C}_{\text{ref}}, \quad r = a_0.$$
Initial stress is identically zero.

---

## 3. Pitfall 2: Unit-Cell Area Normalization

### The Bug
Using `pow(ideal_lattice_parameter, 2.0)` ($a_0^2$) in nodal fields or stress evaluations:
- For square lattice: $\Omega_{\text{cell}} = a_0^2$.
- For triangular lattice: $\Omega_{\text{cell}} = \frac{\sqrt{3}}{2} a_0^2$.
Using $a_0^2$ for triangular lattice causes an error of $\frac{\sqrt{3}}{2} \approx 0.866025$ in total energy and fitted core parameters.

### The Rule
Always call `calculator.getUnitCellArea()`:
```cpp
double norm = userData->calculator.getUnitCellArea();
if (norm <= 1e-12) norm = 1.0;
```
For `SquareLatticeCalculator(calc_scale)`:
$$\text{getUnitCellArea}() = \text{calc\_scale}^2 = \frac{a_0^2}{\sqrt{4/3}} = \frac{\sqrt{3}}{2} a_0^2.$$
This automatically produces the exact unit cell area for any lattice type.

---

## 4. Pitfall 3: Acoustic Tensor Reference State Mapping $\mathbf{F}_0$

### The Bug
In macroscopic stress relaxation, `AcousticTensor::getSpatialTangent` evaluates spatial tangent stiffness $\mathbb{C}^{\text{spatial}}(\bar{\mathbf{F}})$.
If $\mathbf{C}_{\text{bar}} = \bar{\mathbf{F}}^T \bar{\mathbf{F}}$ is evaluated at $\bar{\mathbf{F}} = \mathbf{I}$ relative to $\mathbb{Z}^2$, the acoustic tensor evaluates the stiffness of a square lattice under the hexagonal potential, which has negative shear modulus ($C_{66} < 0$). This causes Newton-Raphson to diverge or oscillate wildly.

### The Rule
Pass the total deformation gradient relative to $\mathbb{Z}^2$:
```cpp
Eigen::Matrix2d F0 = Eigen::Matrix2d::Identity();
if (lattice_type == "triangular") {
    double s = std::pow(4.0 / 3.0, 0.25);
    F0 << s, 0.5 * s,
          0.0, s * std::sqrt(3.0) / 2.0;
}
Eigen::Matrix2d F_tot = F_bar * F0;
lagrange::Result red_res = lagrange::reduce(F_tot.transpose() * F_tot);
AcousticTensor acoustic_tensor(F_tot, red_res.C_reduced, red_res.m_matrix);
```
At $\bar{\mathbf{F}} = \mathbf{I}$, $\mathbf{F}_{\text{tot}} = \mathbf{F}_0 \implies \mathbf{C}_{\text{tot}} = \mathbf{C}_{\text{ref}}$, evaluating the true positive-definite stiffness of the triangular ground state.

---

## 5. Pitfall 4: Adaptive Delaunay Boundary Edge Filtering (`max_edge`)

### The Bug
In non-periodic Delaunay triangulation, outer boundary points can form long, non-physical sliver triangles connecting second-nearest neighbors ($a_0 \sqrt{3} \approx 1.732 a_0$).
The previous threshold:
$$\text{max\_edge} = \text{median\_nn} \cdot (\sqrt{2} + 0.35) = 1.764 a_0 > 1.732 a_0$$
accepted these boundary slivers. Under Volterra displacement step cuts, these slivers collapsed to $\det \mathbf{C} \to 0$, causing energy to explode to $10^8$ and reporting 75 spurious defects.

### The Rule
In `AdaptiveMesher::generateMesh`:
```cpp
bool is_triangular = (offsets[0] > 1e-6 && std::abs(offsets[1] / offsets[0] - std::sqrt(3.0)) < 0.25);
double max_edge = is_triangular ? (median_nn * 1.35) : (median_nn * (std::sqrt(2.0) + 0.35));
```
$1.35 a_0$ cleanly accepts equilateral edges ($1.0 a_0$) and cleanly rejects boundary slivers ($1.732 a_0$).

---

## 6. Summary Checklist for New Simulation Functions

When creating or modifying an experiment function:
- [ ] Is `calc_scale` defined as `optimal_lattice_parameter / std::pow(4.0 / 3.0, 0.25)` for triangular?
- [ ] Is $h$ set to `calc_scale` (or `1.0` only for continuum mode 3)?
- [ ] Is `SquareLatticeCalculator` initialized with `calc_scale`?
- [ ] Is `getUnitCellArea()` used for normalization instead of `pow(optimal_lattice_parameter, 2.0)`?
- [ ] If evaluating `AcousticTensor`, is $\mathbf{F}_{\text{tot}} = \bar{\mathbf{F}} \mathbf{F}_0$ used?
- [ ] If periodic boundaries are used, is $n_y$ guaranteed to be an even integer?
