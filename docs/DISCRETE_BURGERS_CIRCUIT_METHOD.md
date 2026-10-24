# Theoretical Formulation of Discrete Burgers Circuit Dislocation Analysis (Method 2)

## 1. Motivation: Why Topological Triangle Valence Fails on Square Lattices

In a 2D crystal modeled via a Delaunay triangulation $\mathcal{T}$ of atomic positions $\{\mathbf{y}_i\}_{i=1}^N \subset \mathbb{R}^2$, a classic heuristic for identifying dislocations is **topological triangle coordination** (`NodeValence` / `ReferenceCoordination`), defined for each interior vertex $i$ as the number of triangles incident to $i$:
$$z_i^{\text{tri}} = \#\{T \in \mathcal{T} : i \in T\}.$$

### 1.1 Six-Fold Triangular Lattices ($p = 6$)
In an ideal 2D triangular (hexagonal) Bravais lattice with spacing $a_0$, each atom has $6$ equidistant nearest neighbors at distance $a_0$, and every elementary Delaunay triangle is equilateral with edge lengths $(a_0, a_0, a_0)$. Consequently, in a defect-free crystal, $z_i^{\text{tri}} = 6$ everywhere in the bulk, and an isolated edge dislocation corresponds to a $5/7$ disclination dipole ($z_i^{\text{tri}} = 5$ adjacent to $z_j^{\text{tri}} = 7$).

### 1.2 Four-Fold Square Lattices ($p = 4$) and the Diagonal Ambiguity
In a 2D square Bravais lattice $\mathbb{Z}^2 a_0$, each atom has only **$4$ nearest neighbors** at distance $a_0$ and **$4$ second-nearest neighbors** at distance $\sqrt{2}\,a_0$. Every elementary $a_0 \times a_0$ square plaquette is cyclic (its four vertices lie on a circle of diameter $\sqrt{2}\,a_0$). To triangulate the square lattice, each square plaquette must be split into two right-angled isosceles triangles of edge lengths $(a_0, a_0, \sqrt{2}\,a_0)$ by choosing one of the two diagonals:
- the **positive diagonal** `/` connecting $(m, n) \leftrightarrow (m+1, n+1)$, or
- the **negative diagonal** `\` connecting $(m+1, n) \leftrightarrow (m, n+1)$.

For any interior vertex $i$, the four surrounding squares automatically contribute $4$ nearest-neighbor edges ($2$ horizontal and $2$ vertical). In addition, anywhere from $d_i \in \{0, 1, 2, 3, 4\}$ of the four surrounding square diagonals may be incident to vertex $i$. Because the number of incident triangles equals the number of incident mesh edges in the bulk, the triangle valence at vertex $i$ is:
$$z_i^{\text{tri}} = 4 + d_i \in \{4, 5, 6, 7, 8\}.$$

| Surrounding Diagonals Touching Node $i$ ($d_i$) | Triangle Valence $z_i^{\text{tri}} = 4 + d_i$ | Physical State of the Square Crystal |
| :---: | :---: | :--- |
| $d_i = 2$ (uniform `/` or uniform `\` grid) | $6$ | Perfect square lattice |
| $d_i = 1$ (boundary between `/` and `\` domains) | $5$ | **Perfect square lattice (false $5$-defect)** |
| $d_i = 3$ (boundary between `/` and `\` domains) | $7$ | **Perfect square lattice (false $7$-defect)** |
| $d_i = 0$ or $d_i = 4$ (alternating checkerboard diagonals) | $4$ or $8$ | **Perfect square lattice (false $4/8$-defect)** |

Any microscopic numerical rounding or tiny elastic strain flips the Delaunay diagonal of a square plaquette between `/` and `\` without moving the atoms out of the square ground state. Whenever adjacent squares have different diagonal orientations, $z_i^{\text{tri}}$ becomes $5$ and $7$, creating **spurious $5/7$ dislocation pairs in a 100% defect-free crystal**.

---

## 2. Continuum vs. Discrete Multiplicative Plasticity ($\mathbf{F} = \mathbf{F}_e \mathbf{F}_p$)

### 2.1 Continuum Kröner–Bilby Decomposition and Nye's Tensor
In continuum crystal plasticity, the total deformation gradient $\mathbf{F} = \nabla \mathbf{y}$ on a fixed reference mesh is factored via the Kröner–Lee multiplicative decomposition:
$$\mathbf{F} = \mathbf{F}_e \mathbf{F}_p,$$
where:
- $\mathbf{F}_p \in \mathrm{GL}(2, \mathbb{Z})$ (via Lagrange reduction $\mathbf{C}_{\text{red}} = \mathbf{M}^{\mathsf{T}} \mathbf{C} \mathbf{M}$ with $\mathbf{M} = \mathbf{F}_p^{-1} \in \mathrm{SL}(2, \mathbb{Z})$) is the lattice-invariant plastic slip mapping the reference Bravais lattice back to itself, and
- $\mathbf{F}_e = \mathbf{F} \mathbf{F}_p^{-1} = \mathbf{F} \mathbf{M}$ is the **elastic distortion tensor** measuring distortion relative to the nearest stress-free Bravais frame.

In a defect-free elastically deformed crystal, $\mathbf{F}_e^{-1} d\mathbf{y} = d\mathbf{X}$ is an exact differential ($\operatorname{Curl}(\mathbf{F}_e^{-1}) = \mathbf{0}$). Around a dislocation core, the line integral of the elastic lattice frame around a closed contour $\mathcal{C}$ defines the continuum **Burgers vector**:
$$\mathbf{b}(\mathcal{C}) = -\oint_{\mathcal{C}} \mathbf{F}_e^{-1}\,d\mathbf{y} = -\int_{\Omega} \operatorname{Curl}(\mathbf{F}_e^{-1})\,dA,$$
where $\bm{\alpha} = -\operatorname{Curl}(\mathbf{F}_e^{-1})$ is **Nye's dislocation density tensor**.

### 2.2 Why Remeshing Absorbs $\mathbf{F}_p$ into the Mesh Topology
In the MTM code, when `remeshing` is enabled, topological edge flips reconnect atoms to their **new physical neighbors** after each plastic slip event. Thus, on every triangle away from an active dislocation core, the remeshed element deformation gradient $\mathbf{F}_T$ is already close to the identity ($\mathbf{M}_T = \mathbf{I}$). The plastic history $\mathbf{F}_p$ is no longer stored as a large shear inside $\mathbf{F}_T$; instead, it is stored in the **discrete graph connectivity** of the remeshed triangulation.

**Method 2 (Discrete Burgers Circuit)** is the exact discrete counterpart of $-\oint_{\mathcal{C}} \mathbf{F}_e^{-1}\,d\mathbf{y}$ evaluated directly on the atomic bonds of the current triangulation, requiring neither a fixed reference mesh nor an explicit element-wise $\mathbf{F}_e / \mathbf{F}_p$ split.

---

## 3. Mathematical Formulation of Method 2 (Discrete Exterior Calculus)

Let $\mathcal{T} = (\mathcal{V}, \mathcal{E}, \mathcal{F})$ be the oriented 2D triangulation of the crystal in the externally unloaded frame:
$$\mathbf{x}_a = \mathbf{F}_{\text{ext}}^{-1}\,\mathbf{y}_{n_a} + \mathbf{t}_a, \qquad a \in \mathcal{V},$$
where $\mathbf{F}_{\text{ext}}$ is the macroscopic affine loading tensor and $\mathbf{t}_a$ is the periodic translation vector of vertex $a$.

### 3.1 Ideal Bravais Basis Matrix $\mathbf{A}_0$
Let $\mathbf{A}_0 = [\mathbf{e}_1 \;|\; \mathbf{e}_2] \in \mathbb{R}^{2 \times 2}$ be the matrix whose columns are the primitive Bravais basis vectors of the undeformed crystal with lattice spacing $a_0$ and global orientation angle $\theta_0$:
$$\mathbf{A}_0 = a_0\,\mathbf{R}(\theta_0)\,\mathbf{B}_{\text{lat}}, \qquad
\mathbf{R}(\theta_0) = \begin{pmatrix} \cos\theta_0 & -\sin\theta_0 \\ \sin\theta_0 & \cos\theta_0 \end{pmatrix},$$
where the dimensionless primitive basis $\mathbf{B}_{\text{lat}}$ is:
$$\mathbf{B}_{\text{lat}}^{\text{square}} = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} \quad (p = 4), \qquad
\mathbf{B}_{\text{lat}}^{\text{triangular}} = \begin{pmatrix} 1 & 1/2 \\ 0 & \sqrt{3}/2 \end{pmatrix} \quad (p = 6).$$

Both symmetry class $p \in \{4, 6\}$, spacing $a_0$, and global rotation angle $\theta_0$ are identified automatically from the active triangulation:
1. For each triangle $T \in \mathcal{F}$ with sorted edge lengths $L_1(T) \le L_2(T) \le L_3(T)$, the median aspect ratio $\operatorname{median}_T(L_3 / L_1)$ is $\approx 1.0$ for a triangular lattice and $\approx \sqrt{2} \approx 1.414$ for a square lattice (threshold $1.20$).
2. The global orientation $\theta_0 \in (-\pi/p, \pi/p]$ is extracted via the $p$-atic bond-orientational order parameter over all nearest-neighbor bonds ($L_{ab} \in [0.75\,a_0, 1.22\,a_0]$):
   $$Z_p = \sum_{(a,b) \in \mathcal{E}_{\text{NN}}} e^{i p \phi_{ab}}, \qquad \theta_0 = \frac{1}{p}\operatorname{atan2}(\operatorname{Im} Z_p, \operatorname{Re} Z_p), \qquad \phi_{ab} = \operatorname{atan2}(y_b - y_a, x_b - x_a).$$

### 3.2 Discrete Lattice 1-Cochain $\Delta\mathbf{n}_{ab} \in \mathbb{Z}^2$
For every directed mesh edge $e_{ab} = (a \to b) \in \mathcal{E}$, the continuous relative bond vector in the Bravais frame is:
$$\mathbf{s}_{ab} = \mathbf{A}_0^{-1}(\mathbf{x}_b - \mathbf{x}_a) \in \mathbb{R}^2.$$
We project $\mathbf{s}_{ab}$ onto the integer Bravais lattice $\mathbb{Z}^2$ to define the **discrete lattice step 1-cochain** $\Delta\mathbf{n}_{ab} \in \mathbb{Z}^2$. To guarantee strict anti-symmetry ($\Delta\mathbf{n}_{ba} \equiv -\Delta\mathbf{n}_{ab}$) even at half-integer rounding ties, we impose a canonical total ordering $\prec$ on vertices:
$$\Delta\mathbf{n}_{ab} = \begin{cases}
\operatorname{round}\!\left(\mathbf{A}_0^{-1}(\mathbf{x}_b - \mathbf{x}_a)\right) & \text{if } a \prec b, \\[4pt]
-\operatorname{round}\!\left(\mathbf{A}_0^{-1}(\mathbf{x}_a - \mathbf{x}_b)\right) & \text{if } b \prec a.
\end{cases}$$

### 3.3 Elementary Triangle Burgers Vector (2-Cochain $\mathbf{B}_T = (d\Delta\mathbf{n})_T$)
For every counter-clockwise oriented triangle $T = (0, 1, 2) \in \mathcal{F}$, the discrete exterior derivative (coboundary) of $\Delta\mathbf{n}$ is the closed 3-step Burgers circuit around $\partial T$:
$$\mathbf{B}_T = (d\Delta\mathbf{n})(T) = \Delta\mathbf{n}_{01} + \Delta\mathbf{n}_{12} + \Delta\mathbf{n}_{20} \in \mathbb{Z}^2.$$

The corresponding physical Burgers vector (in dimensionless units of $a_0$) is:
$$\mathbf{b}_T = \frac{1}{a_0}\mathbf{A}_0 \mathbf{B}_T \in \mathbb{R}^2, \qquad b_T = \|\mathbf{b}_T\|_2.$$

### 3.4 Exact Invariance Under `/` vs. `\` Diagonal Flips and Surface Steps

**Proposition 1 (Zero False Positives in Any Defect-Free Region).**
*Let $\Omega \subset \mathcal{T}$ be any defect-free region of the crystal (including regions with `/` vs. `\` diagonal flips or integer surface steps), meaning there exists a single-valued integer lattice assignment $\mathbf{m}_a \in \mathbb{Z}^2$ for each vertex $a \in \Omega$ such that $\mathbf{x}_a = \mathbf{A}_0 \mathbf{m}_a + \mathbf{u}_a$ with differential elastic strain $\|\mathbf{A}_0^{-1}(\mathbf{u}_b - \mathbf{u}_a)\|_\infty < \frac{1}{2}$ on every edge $(a,b)$ of $T \subset \Omega$. Then:*
$$\mathbf{B}_T \equiv \begin{pmatrix} 0 \\ 0 \end{pmatrix} \qquad \forall T \subset \Omega.$$

*Proof.* On each directed edge $(a \to b)$ of $T$,
$$\Delta\mathbf{n}_{ab} = \operatorname{round}\!\left(\mathbf{m}_b - \mathbf{m}_a + \mathbf{A}_0^{-1}(\mathbf{u}_b - \mathbf{u}_a)\right) = \mathbf{m}_b - \mathbf{m}_a = (d\mathbf{m})_{ab}.$$
Because $\Delta\mathbf{n} = d\mathbf{m}$ is an exact 1-coboundary and $d^2 \equiv 0$, its curl around every triangle $T = (0, 1, 2)$ telescopes identically to zero:
$$\mathbf{B}_T = (\mathbf{m}_1 - \mathbf{m}_0) + (\mathbf{m}_2 - \mathbf{m}_1) + (\mathbf{m}_0 - \mathbf{m}_2) \equiv \mathbf{0}.\quad\square$$

In particular, in a square lattice:
- On a `/` diagonal, $\Delta\mathbf{n}_{\text{diag}} = (1, 1) = (1, 0) + (0, 1) \implies \mathbf{B}_T = \mathbf{0}$.
- On a `\` diagonal, $\Delta\mathbf{n}_{\text{diag}} = (-1, 1) = (-1, 0) + (0, 1) \implies \mathbf{B}_T = \mathbf{0}$.

Thus `/` and `\` diagonals produce **zero** fake dislocations.

---

## 4. Discrete Stokes' Theorem, Multi-Resolution Circuits, and Junctions

Because $\Delta\mathbf{n}_{ab} + \Delta\mathbf{n}_{ba} \equiv \mathbf{0}$ in exact integer arithmetic on every interior edge, summing $\mathbf{B}_T$ over any patch of triangles $\Omega = \bigcup_{k} T_k$ cancels all shared interior edges identically:

**Proposition 2 (Discrete Stokes' Theorem & Frank's Rule at Junctions).**
*For any collection of triangles $\Omega \subset \mathcal{F}$ with counter-clockwise outer boundary loop $\partial\Omega$:*
$$\mathbf{B}(\partial\Omega) \equiv \oint_{\partial\Omega} \Delta\mathbf{n}_{ab} = \sum_{T \subset \Omega} \mathbf{B}_T \in \mathbb{Z}^2.$$

### 4.1 Unit-Cell Plaquette Cancellation Across Shared Square Diagonals
In a square lattice, the diagonal edge of a unit cell has length $\sqrt{2}\,a_0 \approx 1.414\,a_0$, whereas the four outer perimeter edges of the unit cell have length $1.0\,a_0$. Under severe local elastic shear near a defect, the longer internal diagonal edge $e_{\text{diag}}$ shared by two triangles $T_1, T_2$ can cross a rounding threshold even when the 4-atom square unit cell $Q = T_1 \cup T_2$ contains no dislocation. In that event, the outer 4-edge circuit has zero Burgers vector:
$$\mathbf{B}_{\text{quad}} = \mathbf{B}_{T_1} + \mathbf{B}_{T_2} = \mathbf{0}, \qquad \mathbf{B}_{T_1} = -\mathbf{B}_{T_2} \neq \mathbf{0}.$$
Whenever two triangles sharing a diagonal edge ($L_{\text{shared}} > 1.15\,a_0$) satisfy $\mathbf{B}_{T_1} + \mathbf{B}_{T_2} = \mathbf{0}$, both are cancelled ($\mathbf{B}_{T_1} = \mathbf{B}_{T_2} = \mathbf{0}$), eliminating single-diagonal elastic rounding artifacts without altering $\sum_{T \subset \Omega} \mathbf{B}_T$ on any larger circuit.

### 4.2 Nodal 1-Ring Circuits and Dislocation Junctions
For each vertex $i \in \mathcal{V}$, the **1-ring Burgers circuit** $\partial\Omega_i$ is the closed polygon formed by the neighbor vertices of $i$:
$$\mathbf{B}_i^{(1\text{-ring})} = \sum_{T \ni i} \mathbf{B}_T \in \mathbb{Z}^2, \qquad \mathbf{b}_i^{(1\text{-ring})} = \frac{1}{a_0}\mathbf{A}_0 \mathbf{B}_i^{(1\text{-ring})}.$$
- **Isolated dislocations:** A single elementary dislocation core resides inside a triangle $T^*$ with $\|\mathbf{b}_{T^*}\|_2 = 1$, and its three vertices have $\|\mathbf{b}_i^{(1\text{-ring})}\|_2 = 1$.
- **Dislocation reactions & junctions:** If two dislocations with Burgers vectors $\mathbf{B}_1$ and $\mathbf{B}_2$ react or form a junction (e.g., $\mathbf{B}_1 = (1,0)$ from the horizontal shift and $\mathbf{B}_2 = (0,1)$ from the vertical shift in `Shift / Staircase`), Discrete Stokes' Theorem automatically enforces **Frank's nodal conservation law**:
  $$\mathbf{B}_{\text{junction}} = \mathbf{B}_1 + \mathbf{B}_2 = (1, 1).$$

### 4.3 Summary of the 3 Nested Loop Sizes Used in the Code

| Level | Loop Size | Formula | Why We Use It |
| :--- | :--- | :--- | :--- |
| **Level 1: Single Triangle** | **3 atoms** ($1$ triangle $T$) | $\mathbf{B}_T = \Delta\mathbf{n}_{01} + \Delta\mathbf{n}_{12} + \Delta\mathbf{n}_{20}$ | Elementary building block (`ElementBurgersMagnitude`). Pinpoints the exact triangle containing the core. |
| **Level 2: Square Unit Cell** | **4 atoms** ($2$ triangles $T_1 \cup T_2$ sharing a diagonal) | $\mathbf{B}_{\text{square}} = \mathbf{B}_{T_1} + \mathbf{B}_{T_2}$ | **Anti-noise filter:** If elastic shear causes *only* a shared square diagonal ($\sqrt{2}\,a_0$) to round by $\pm 1$, then $\mathbf{B}_{T_1} = -\mathbf{B}_{T_2}$ and the 4-atom square loop has $\mathbf{B}_{\text{square}} = (0,0)$, so the code cancels both to zero. |
| **Level 3: Nodal 1-Ring Circuit** | **5–7 neighbor atoms** surrounding atom $i$ ($\bigcup_{T \ni i} T$) | $\mathbf{B}_i^{(1\text{-ring})} = \sum_{T \ni i} \mathbf{B}_T$ | **True multi-triangle Burgers circuit around atom $i$** (`BurgersMagnitude` & `ReferenceCoordination`). All spoke edges touching atom $i$ cancel out, leaving the closed polygon loop of neighbors walking around atom $i$. |

---

## 5. Implementation & VTK Output Reference

The algorithm is implemented in `DefectAnalysis::computeBurgersCircuits` (`src/DefectAnalysis.cpp`) and exported to VTK by `ConfigurationSaver::writeToVTKFile` and `ConfigurationSaver::writeToVTK_DefectAnalysis` (`src/configuration_saver.cpp`).

| VTK Field Name | Location | Type | Description |
| :--- | :--- | :--- | :--- |
| `ElementBurgersMagnitude` | `CELL_DATA` | `float` (scalar) | $\|\mathbf{b}_T\|_2 = \|\mathbf{A}_0 \mathbf{B}_T\|_2 / a_0$ on each triangle $T$ ($0$ in defect-free crystal, $1$ at an elementary dislocation core, $\sqrt{2}$ at a $(1,1)$ junction) |
| `ElementBurgersVector` | `CELL_DATA` | `float` (3-vector) | $(b_{T,x}, b_{T,y}, 0)$ physical Burgers vector of triangle $T$ in units of $a_0$ |
| `BurgersMagnitude` | `POINT_DATA` | `float` (scalar) | $\|\mathbf{b}_i^{(1\text{-ring})}\|_2$ around the 1-ring of node $i$ in units of $a_0$ |
| `BurgersVector` | `POINT_DATA` | `float` (3-vector) | $(b_{i,x}^{(1\text{-ring})}, b_{i,y}^{(1\text{-ring})}, 0)$ 1-ring Burgers vector at node $i$ |
| `ReferenceCoordination` | `POINT_DATA` | `float` (scalar) | Burgers-filtered coordination: `6` at all nodes with $\mathbf{B}_i^{(1\text{-ring})} = \mathbf{0}$, and $5/7$ only at genuine dislocation core nodes ($\mathbf{B}_i^{(1\text{-ring})} \neq \mathbf{0}$) |
| `NodeValence` | `POINT_DATA` | `float` (scalar) | Raw unfiltered mesh triangle count per node (preserved for mesh connectivity inspection) |

---

## 6. Historical Background and References

1. **J. M. Burgers (1939)**, *"Some considerations on the fields of stress connected with dislocations in a regular crystal lattice,"* Proc. Kon. Ned. Akad. Wetensch. **42**, 293–325 & 378–399.
2. **F. C. Frank (1951)**, *"Crystal dislocations — Elementary concepts and definitions,"* Philosophical Magazine **42**(331), 809–819.
3. **J. F. Nye (1953)**, *"Some geometrical relations in dislocated crystals,"* Acta Metallurgica **1**(2), 153–162.
4. **B. A. Bilby, R. Bullough, and E. Smith (1955)**, *"Continuous distributions of dislocations: a new application of the methods of non-Riemannian geometry,"* Proc. R. Soc. Lond. A **231**(1185), 263–273.
5. **E. Kröner (1960)**, *"Allgemeine Kontinuumstheorie der Versetzungen und Eigenspannungen,"* Archive for Rational Mechanics and Analysis **4**, 273–334.
6. **S. Conti and G. Zanzotto (2004)**, *"A variational model for reconstructive phase transformations in crystals, and their relation to dislocations and plasticity,"* Archive for Rational Mechanics and Analysis **173**, 69–88.
7. **C. S. Hartley and Y. Mishin (2005)**, *"Characterization and visualization of the lattice misfit associated with dislocation cores,"* Acta Materialia **53**(5), 1313–1321.
8. **M. P. Ariza and M. Ortiz (2005)**, *"Discrete crystal elasticity,"* Archive for Rational Mechanics and Analysis **178**(2), 179–226.
9. **A. Stukowski and K. Albe (2010)**, *"Extracting dislocations and non-dislocation crystal defects from atomistic simulation data,"* Modelling Simul. Mater. Sci. Eng. **18**(8), 085001.
10. **A. Stukowski, V. V. Bulatov, and A. Arsenlis (2012)**, *"Automated identification and indexing of dislocations in crystal interfaces,"* Modelling Simul. Mater. Sci. Eng. **20**(8), 085007.

