# Localized Topological Energy Changes ($\Delta E_{\text{topo}}$)

This document summarizes the mathematical formulation and methods for extracting **localized energy changes** during Delaunay reconnection / remeshing in 2D crystal lattices.

---

## 1. Why Direct Element-by-Element Subtraction Fails

In finite element calculations, energy is typically stored on elements $e \in \mathcal{T}$:
$$E_{\text{total}} = \sum_{e \in \mathcal{T}} E_e$$

When Delaunay remeshing occurs:
1. The old mesh $\mathcal{T}_{\text{old}}$ is replaced by a new mesh $\mathcal{T}_{\text{new}}$.
2. Even though the nodes $\{ \mathbf{y}_a \}$ do not change positions at the instant of reconnection, the element topology changes (edge flips, rewiring around defect cores).
3. Triangles do not have a 1-to-1 correspondence between $\mathcal{T}_{\text{old}}$ and $\mathcal{T}_{\text{new}}$. Hence, $\Delta E_e = E_e^{\text{new}} - E_e^{\text{old}}$ is ill-defined on the element indices.

---

## 2. Does it Give Multiple / Many Energy Changes?

**Yes, absolutely.** Instead of a single lumped global scalar $\Delta E_{\text{topo}}$, local analysis produces:

1. **A Spatial Distribution / Spectrum of Jumps Across Defect Sites:**
   - When an avalanche occurs across a large specimen, reconnections happen simultaneously at multiple disconnected dislocation cores, slip bands, or dipole pairs.
   - Each local defect cavity $p$ undergoes its own distinct topological jump: $\Delta E_{\text{topo}}^{(p)}$.
   - One obtains a **multiset / histogram / spatial field** of discrete energy increments, revealing the individual energetic cost/barrier of each local reconnection event.

2. **Nodal Energy Change Field ($\Delta E_a$ for each node $a$):**
   - For all nodes $a = 1, \dots, N$, $\Delta E_a$ provides an $N$-dimensional vector / field.
   - Nodes far from the reconnections experience $\Delta E_a = 0$.
   - Nodes directly involved in edge flips and core reconfigurations exhibit sharp localized jumps.

3. **Multi-Pass Avalanche Sequence (Temporal Jumps):**
   - During complex plastic avalanches, remeshing runs in consecutive passes ($k = 1, 2, \dots, K$).
   - Each pass $k$ produces its own topological jump $\Delta E_{\text{topo}}^{(k)}$, followed by a continuous relaxation drop $\Delta E_{\text{relax}}^{(k)}$.

---

## 3. Localization Formulations

### Method A: Cavity / Patch Decomposition (Exact Discrete Groups)
1. **Identify Unchanged vs. Changed Elements:**
   - An element $T = (i, j, k)$ is **untouched** if $\{i, j, k\} \in \mathcal{T}_{\text{old}} \cap \mathcal{T}_{\text{new}}$. For these, $\Delta E_T = 0$ identically at fixed nodal positions.
   - The set of removed elements $\mathcal{T}_{\text{removed}} = \mathcal{T}_{\text{old}} \setminus \mathcal{T}_{\text{new}}$ and newly formed elements $\mathcal{T}_{\text{added}} = \mathcal{T}_{\text{new}} \setminus \mathcal{T}_{\text{old}}$ partition into spatially disjoint connected components (cavities/patches $p = 1, \dots, P$).
2. **Cavity Energy Jump:**
   $$\Delta E_{\text{topo}}^{(p)} = \sum_{e \in \mathcal{T}_{\text{added}}^{(p)}} E_e^{\text{new}} - \sum_{e \in \mathcal{T}_{\text{removed}}^{(p)}} E_e^{\text{old}}$$
   $$\Delta E_{\text{topo}} = \sum_{p=1}^P \Delta E_{\text{topo}}^{(p)}$$
   Each patch $p$ has an exact physical location (its centroid or bounding box) and its own localized energy jump.

### Method B: Nodal Lumping (Scalar Field on Deformed Mesh)
1. **Lump Element Energies to Vertices:**
   $$E_a = \sum_{e \ni a} \frac{1}{3} E_e \quad (\text{or area-weighted})$$
2. **Nodal Difference:**
   $$\Delta E_a = E_a(\mathcal{T}_{\text{new}}) - E_a(\mathcal{T}_{\text{old}})$$
3. **Properties:**
   - Identical node ordering $a = 1, \dots, N$ is preserved before and after reconnection.
   - $\sum_{a=1}^N \Delta E_a = \Delta E_{\text{topo}}$.
   - Can be exported directly as a point field (`vtkPointData`) for contour plots, heatmaps, and spatial correlation with plastic slip planes.

### Method C: Continuous Energy Density Projection
- Energy density $w(\mathbf{x})$ is piecewise-constant on triangles.
- Integrating over sub-regions or dual Voronoi cells gives local volumetric energy increments $\Delta w(\mathbf{x})$.
