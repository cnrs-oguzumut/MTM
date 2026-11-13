# How the lattice rotation ω (and orientation φ) is measured

Code: `ConfigurationSaver::measureLatticeOrientation` and `ConfigurationSaver::logLatticeRotation`
in `src/configuration_saver.cpp`. Outputs: `lattice_rotation_log.csv` (one row per load step) and the
VTK fields `ElementRotationAngleDeg`, `ElementLatticeOrientationDeg`, `ElementLatticeDistortion`,
`NodalLatticeRotationDeg`.

ω and φ are a new, independent measure: they do not use F_aligned or anything in
`docs/aligned_deformation_gradient.tex`. Section 3 only sets the two side by side, because the
old `ElementRotationAngleDeg` field was derived from F_aligned.

---

## 1. Per element

1. **Current edges.** The element's deformation gradient F is built from the current positions of
   its three atoms (`F = x_e · dN/dX`). It is turned back into the real bond vectors:

       E = [x1 − x0, x2 − x0] = F · E_ref,   E_ref⁻¹ = rows (∇N1, ∇N2) of dN/dX

   E is the current geometry only: it does not depend on the reference frame (SquareZ2 or
   Triangular), on the parent-triangle shape, or on the starting positions X0. Periodic images are
   already handled inside F. If the node order is clockwise, the two columns are swapped.

2. **Best-fitting lattice cell.** The two element edges are not always the best pair of lattice
   vectors (e.g. after remeshing in a sheared region). Every unimodular basis change M with
   entries in {−1, 0, 1} and det M = 1 (about 20 candidates, identity first) is tried:

       G = E · M · B⁻¹,   B = [(1, 0), (1/2, √3/2)] (triangular)  or  B = I (square)

   Scale is removed (Gn = G / √det G) and the candidate closest to a pure rotation is kept:

       score = ‖Gnᵀ Gn − I‖²,   distortion = √(score / 2)  (≈ equivalent shear)

   `distortion` is written as `ElementLatticeDistortion`.

3. **Rotation of the fit.** Closed-form 2D polar decomposition G = R U:

       φ = atan2(G21 − G12, G11 + G22)

   φ is the current lattice orientation, defined only modulo the lattice symmetry (60° triangular,
   90° square).

4. **Rotation since the start.**

       ω = wrap(φ − θ)   into (−30°, 30°] (triangular) or (−45°, 45°] (square)

   θ = the grain's starting angle (majority grain label of the 3 nodes; for a rotated single
   crystal the reference rotation). Written as `ElementRotationAngleDeg`.

5. **Reported orientation.** `ElementLatticeOrientationDeg` = θ + ω, i.e. φ folded around the
   grain's own starting angle rather than around 0°. A grain whose starting angle is near the fold
   (e.g. θ = −26° triangular) therefore stays one continuous colour instead of splitting into
   ≈ −30° and ≈ +30° cells; values can lie outside ±30° / ±45°. φ only wraps if a cell rotates
   more than half the symmetry angle from its grain's start.

**Excluded from statistics (still written to VTK):** inverted elements; grain-boundary elements
(nodes from more than one grain); elements with distortion > 0.35 (dislocation cores,
phase-transformed regions) — counted in `NumUnrecognizedElements`.

## 2. Per grain and per step (`lattice_rotation_log.csv`)

- **Grain mean**: area-weighted circular mean with the lattice symmetry (k = 6 triangular, 4 square),
  w = current element area:

      ω̄ = atan2(Σ w sin kω, Σ w cos kω) / k

  so values near ±30° do not cancel.
- **Unwrapping**: each step's mean is joined to the previous one
  (ω̄_n = ω̄_{n−1} + wrap(ω̄ − ω̄_{n−1})), so a grain rotating past 30° gives a continuous curve →
  `Grain{k}RotationDeg`.
- **Orientation**: `Grain{k}OrientationDeg` = θ_k + `Grain{k}RotationDeg` (not wrapped).
- **Applied rotation**: polar angle of the macroscopic F of the loading (`F_ext`, or
  `computeMacroscopicF(α)` for custom geometries); −atan(γ/2) for simple shear.
  `MeanRotationNetDeg` = mean − applied.
- Also: `StdRotationDeg`, `MaxAbsRotationDeg`, `FracAbsRotationAbove2Deg`, element counts.
- `NodalLatticeRotationDeg`: per-atom area-weighted circular mean of ω over its non-boundary,
  recognised elements.

**Physical meaning.** Atoms sliding along a slip plane do not change bond directions, so dislocation
glide alone gives ω = 0. ω measures the turning of the lattice (elastic rotation, grain rotation,
bending, kinking): the lattice rotation an EBSD map would show, not the material rotation.

---

## 3. Comparison with `aligned_deformation_gradient.tex` (F_aligned)

| | **F_aligned** (stress post-processing) | **ω / φ** (this measure) |
|---|---|---|
| Purpose | Remove the arbitrary triangle labelling so that P and the projected stress do not cancel between elements | Measure how the lattice is oriented / has rotated |
| Input | Raw F relative to the parent triangle (frame-dependent) | Current edge vectors E = F·E_ref (frame-independent) |
| Candidates | SquareZ2: 12 = 3 cyclic vertex shifts M_p × 4 rotations R_r. Triangular reference: 6 rotations R(k·60°) | ~20 unimodular M with entries in {−1,0,1}, applied to the edges, then compared to the ideal cell B |
| Selection rule | Closest to a **fixed target** R(θ0)·F_target: argmin ‖F M_p R_rᵀ − F_target‖ | Closest to **any rotation** (shape only, scale removed): argmin ‖Gnᵀ Gn − I‖ |
| Depends on θ | Yes: the target contains R(θ0) | No: θ only enters afterwards, ω = φ − θ |
| Result for a pristine triangular grain in SquareZ2 | F_aligned = R(θ0)·F_hex (not a rotation) | G = R(θ0) exactly |
| Lattice type in the formula | Unified: no branching (square and hexagonal wells handled by one target) | Branches on the lattice: B = equilateral or square, symmetry 60° or 90° (lattice type recorded by `make_reference_setup`) |

**Why the old `ElementRotationAngleDeg` was wrong for triangular lattices.** It took the polar
rotation of R(θ)ᵀ·F_aligned. In SquareZ2 mode a pristine triangular element has
F_aligned = R(θ0)·F_hex, and F_hex = s0 [[1, ½], [0, √3/2]] is not a pure stretch: its own polar
rotation is atan2(−½, 1 + √3/2) = **−15°**. On top of that, the snapping can pick a candidate on the
other side of the fundamental domain, so every grain carried a different spurious offset (−13°, +9°,
… in the periodic Voronoi test). The new measure divides by the ideal cell B, so a pristine grain
gives exactly G = R(θ), i.e. ω = 0. F_aligned itself is still the right object for stress: the
Cauchy stress, energy and reduced metric are invariant under its choice (Sec. 6 of the PDF), and it
is left unchanged.

**Discrepancy found between the PDF and the code.** The PDF (Sec. 5, "Unified Formula") places the
SquareZ2 target at the midpoint F12 = 1/4, so that square and triangular phases are both at
distance ≈ 1/4. The code in `computeAlignmentMatrixInternal` uses
`B_target = [[1, 1e-6], [0, 1]]`, i.e. F12 ≈ 0: the square well, with a tiny tie-breaker. One of
the two is out of date; the stress outputs follow the code.

## 4. Validation (this branch)

| Case | Expected | Measured |
|---|---|---|
| Single crystal rotated 17°, triangular LJ, SquareZ2 and Triangular reference; Zanzotto square | ω = 0 at α = 0 | < 1e-9° |
| Same, sheared to γ = 0.01 | ω ≈ −atan(γ/2) = −0.287° | −0.286°, net < 0.001° |
| Periodic triangular box under shear | ω = applied | equal to 1e-15° |
| 6-grain periodic Voronoi after relaxation | ω ≈ 0 per grain | within ±1.4° (old field: −13°, +9°, …) |
| Zanzotto triangular, β = −0.25 (C66 < 0, transforms) | lattice not recognised | 48 of 1678 elements still triangular |

Known limits: per-element ω is shown only within ±30° (±45°) (the per-grain columns are
unwrapped); atoms right at a grain boundary can read large ω because they carry one grain's
label but sit in the neighbour's lattice (use `ElementLatticeOrientationDeg` there).
