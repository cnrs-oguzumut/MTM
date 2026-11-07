# Changelog

## Unreleased (feature/gui-onboarding-redesign)

### Added — stiffness preconditioner for PBC stress control
- With `--precond=stiffness` (the GUI default), stress-controlled PBC relaxations (nominal and Cauchy) are now
  preconditioned too. Before, they silently ran plain L-BFGS whatever the preconditioner setting. The
  preconditioner is the Gauss-Newton stiffness in the solver's unknowns (w, √A₀ F̄): the atom block DᵀKD
  (x = F̄w) bordered by the box rows K_xF and K_FF, assembled from the same PSD-projected element tangents.
  It is factored as a bordered Cholesky: a sparse CHOLMOD factor of the atom block, ≤ 3 triangular solves for
  the border, and a dense Cholesky of the ≤ 3×3 Schur complement (the preconditioner's relaxed box stiffness).
- The relaxation stops at max|∂H/∂z| < `--precond-tol` (1e-6), so the box stress matches its target to
  1e-6/√A₀. It is rebuilt every relaxation (each Cauchy pass), and at the current point when L-BFGS stalls
  during an avalanche (≤ 3 times). The 0.5 a₀ step cap of the plain solver is not used in the preconditioned
  variables. If a factorization fails it falls back to plain L-BFGS. `MTM_STRESS_PLAIN=1` forces the plain
  solver. Other preconditioner types still run plain L-BFGS under stress control. The strain-control code is
  unchanged.
- 100×100 defected periodic cell (relaxed at F₁₂ = 0.14 with noise, Cauchy shear in steps of 0.005, remeshing
  on): 35 s vs 119 s (seed 42) and 51 s vs 134 s (seed 7). Elastic steps are 3–5× faster, avalanche steps 2–7×.
  Each relaxation is now ≤ ~40 iterations and ≤ ~0.1 s, so an elastic step's time is mostly the remesh trial
  and output. Stresses match as before: σ_xy to ~1e-8 of the target, normal components ~1e-9.
- The two solvers give the same states until the first avalanche and can then settle into different metastable
  states (seed 42: limit load at 0.047 plain, none up to 0.057 preconditioned). That is avalanche path
  dependence, not an error.
- Unchanged code paths are byte-identical to HEAD: stress control with `--precond=none`, and stress control with
  `MTM_STRESS_PLAIN=1`. Strain-controlled runs with the stiffness preconditioner change by ~1e-14. The strain
  code is untouched; the build uses `-flto -ffast-math`, so adding code changes whole-program inlining and
  rounding.
- docs/stress_and_traction_control.tex explains the construction (Section "Variant C").

### Added — force control (dead-load tractions, plan Phase A)
- New boundary kind `traction` (`traction_mode`: `cartesian` with `tx`/`ty`, `radial_outward`, `radial_inward`,
  `circumferential_shear`). α is the nominal traction. Loaded atoms stay free, and the solver minimizes
  Π = E_el − α Σ f0·(x − X). The force on each outermost contour atom is half the traction times the
  reference length of each boundary edge it touches. Edge rules load only faces that look that way.
- The run is refused (exit code 1) when a dead load is unbalanced and nothing supports the specimen. A balanced
  load with no support gets a 3-DOF isostatic anchor.
- Remesh acceptance uses Π. Under force control an avalanche is a strain burst at fixed α: an enthalpy drop in
  the cascade plus a conjugate-strain jump beyond the elastic secant. A displacement guard stops the run
  cleanly at the limit load, and L-BFGS is capped at 20 000 iterations per relaxation.
- `energy_stress_log.csv` gains `EffectiveStrain`, `ElasticEnergy` and `StrainBurstFraction` columns on runs
  with traction rules only. Strain-controlled outputs are byte-identical to `90d9727`.
- GUI: Traction kind in the rule editor, magenta force arrows on the canvas, Dead-Load Tension/Shear and
  Pressurized Hole presets, a balance warning that disables Run, and a "Dead-Load Strain ε" chart series.
  Set `MTM_TRACTION_GRADCHECK=1` to print a finite-difference check of the enthalpy gradient.

### Added — PBC nominal-stress control (plan Phase C)
- `--stress-control=shear|uniaxial_y|uniaxial_x|relax_normal` on periodic cells (plain PBC and periodic
  dislocations). The box F̄ is upper triangular (F₂₁ = 0). Its stress-driven components are solved for with the
  atoms on H = E − A₀ Pᵗ : (F̄ − I), using internal coordinates w (x = F̄w). The gradients are exact:
  ∂E/∂w = F̄ᵀ∂E/∂x and ∂E/∂F̄ = (Σ A_e P_e F_eᵀ) F̄⁻ᵀ, from a virial accumulated in the energy loop.
  The variables are scaled to lengths, and the L-BFGS step is capped at 0.5 a₀: the homogeneous energy is
  periodic in shear, so H has no lower bound along the box components. The extra variables live only
  inside `relax_configuration`, and remesh rejection restores F̄.
- `--relax-stress=1` on a periodic loading run relaxes σ_xx = σ_yy = σ_xy = 0 at α = 0 before loading.
  Volterra mode keeps its own outer loop.
- Shares the force-control avalanche test, the limit-load stop (a stress-driven F̄ component changing by more
  than 0.5 in one step) and the CSV columns. Adds `Fbar11`, `Fbar12` and `Fbar22` columns. Prints the bulk
  Voigt stiffness and the uniaxial Poisson ratio it predicts. Set `MTM_STRESS_GRADCHECK=1` /
  `MTM_STRESS_TRACE=1` for finite-difference and convergence diagnostics.
- **Cauchy (true-stress) presets** `--stress-control=cauchy_shear|cauchy_uniaxial_y|cauchy_uniaxial_x` set
  σ_xy, σ_yy or σ_xx = α, with the other in-plane Cauchy stresses at zero. The nominal presets are kept;
  under nominal shear, σ_xx = P₁₂F₁₂/J is not zero. A fixed Cauchy stress is a follower load, so each
  relaxation is a fixed-point loop on the nominal target Pᵗ = J σᵗ F̄⁻ᵀ, re-evaluated at the relaxed box
  (about 6 passes). The controlled σ holds to about 1e-8. `MTM_STRESS_TRACE=1` prints the passes and the
  box stress from the solver's virial. In the GUI, the periodic box control picker lists the Cauchy and
  nominal presets.
- **Relaxation at a nonzero shear strain** (`--relax-strain=<F12>`) for plain periodic cells with a shear
  stress preset (`shear`, `cauchy_shear`). The relaxation step holds the box at F₁₂ under strain control
  with the initial noise and a remesh cascade, so dislocations nucleate on the unstable homogeneous branch
  (e.g. 0.14 in a large box). Stress control then starts from the relaxed state, with targets σ₀ + k·|dα|
  up to α_end (σ₀ is the relaxed stress). The conjugate strain is the absolute F₁₂. In the GUI, the
  "1. Relaxation Step" panel has a "Relax at shear strain F₁₂" field.
- `--legacy-stress-servo` runs the old `example_3_stress_controlled_final_clean` Cauchy servo as a baseline.
  It starts from a hard-coded F₁₂ = 0.14 pre-strain on the softening branch, so it can't be compared
  point by point with the minimizer.
- GUI: periodic box control picker (strain / shear stress / uniaxial stress / free normals), a "Relax box
  stress to zero at α = 0" toggle, and nominal-stress labels.

### Added — rigid force grip (plan Phase B)
- New boundary kind `force_grip`. The gripped atoms move together, x = X + U, and U is solved for together
  with the atoms: Π = E_el − α F0·U, with F0 = (tx, ty)·L_face. U lives only inside `relax_configuration`
  (plain L-BFGS; the stiffness preconditioner covers node DOFs only), so the rest of the pipeline is
  unchanged. A loaded component of U is always free. An unloaded one is held (`lock_x` / `lock_y`) or slides
  at zero force.
- A grip needs another kinematic support (the run is refused otherwise). It shares the force-control
  avalanche, limit-load and CSV machinery with dead-load tractions. Set `MTM_GRIP_GRADCHECK=1` for a
  finite-difference check of ∂Π/∂U.
- GUI: Force Grip kind (tx/ty, hold toggles for the unloaded axis), grip bar and force arrow on the canvas, a
  Force-Grip Tension preset, and the schedule reads "nominal stress α" under force control.

### Changed — physics
- **Affine origin.** Centre-symmetric affine modes (biaxial, radial in/out, pure shear, torsion) pivot
  about the specimen centre for `outer_perimeter` rules, and about each hole's own centre for
  `hole_surfaces` rules (circle/ellipse/rectangle centre, polygon vertex average, pre-crack midpoint).
  Edge and custom-box affine rules use the bottom-left corner again, so a one-sided grip (e.g. top
  `biaxial` + bottom fixed) is bit-identical to `ec7d83c`.
- **Pure shear reports ½(σxx − σyy).** This is the stress conjugate to F = diag(1+α, 1−α); σxy is ≈ 0
  under pure shear. Earlier builds reported σxy, so pure-shear curves from older runs are not comparable.
- **Torsion is an exact rotation.** Rim nodes driven by `circumferential_shear` / `torsion` follow
  x = X_c + R(α)(X − X_c) with α in radians, instead of the linearised X + α·v (which added a spurious
  √(1+α²) stretch). The free-node predictor uses the matching increment R(α) − R(α − Δα).
- **Torsion reports the rim torque** T = Σ (x − X_c) × f over the driven rim nodes
  (`primary_stress_component: "torque"`, also chosen by `auto`). The solver warns when torsion runs
  without any fixed or roller node, since the specimen then only rotates rigidly.

### Fixed
- **Limit load keeps the last equilibrium.** Past the limit load a dead load or stress target has no
  equilibrium, and the run used to log and save the runaway configuration as its final point. In a sheared
  periodic cell that was a box collapsed to a line in ParaView, and a vertical drop on the chart. Now the
  run restores the state at the start of the step and saves it as the final `configuration_*.vtk` and
  `checkpoints/latest.chk` (labelled with its own α). Nothing is logged for the load without equilibrium.
- **GUI chart under force / stress control** plots against the conjugate strain ε by default, so strain
  bursts at fixed load show as jumps along x. An `x: ε | x: α` toggle appears for such runs, the readout
  shows both, and the series is renamed "Conjugate Strain ε (force / stress control)".
- **GUI energy under force / stress control.** The logged energy there is the enthalpy H = E_el − W_ext,
  dominated by the load term (≈ −A₀σε). The Energy tab now plots the elastic energy E_el (`ElasticEnergy`
  column) by default, with an `E_el | H` toggle, and the readout shows both. The status pill shows the
  conjugate strain ε instead of α (which is the stress there).
- **GUI avalanche count under force / stress control.** The count was based on drops of the logged energy.
  Under force or stress control that energy is the enthalpy, which falls at every load step, so every step
  counted. Such runs now count the solver's own flags (`fractional_drops.csv`).
- **Element reference areas on pre-displaced meshes.** Every triangle now carries the reference element's
  area h²/2, the area its F is measured against (`ElementTriangle2D::referenceElementArea(dndx)`).
  - Before, the mesher used each triangle's current area at meshing time, so dislocation-dipole and
    grain-boundary meshes weighted the energy with distorted areas (0.39–0.63 instead of 0.5).
  - Each remesh trial then stamped all elements with `elements[0]`'s area, so remesh acceptance compared
    energies weighted differently.

  Single-crystal outputs are byte-identical, apart from a 1e-5 change in a triangular Lennard-Jones
  dead-load case. Polycrystal runs change: on the dogbone the earlier "plastic
  relaxation" remeshes at α = 0 and 0.01 are now rejected. Periodic dipoles now log σ = α under stress
  control, with no spurious σ_xx or σ_yy.
- **PBC stress control A0** is the periodic cell's reference area, not the summed element areas. These
  differ when a dipole inserts extra half-planes (dipole: 576 vs 586).
- **Checkpoint restart of custom-geometry / fixed-BC runs.** `--restart` hard-coded a periodic cell, so
  any other run restarted as a different problem. Now it takes the original run's `--geometry-json`, `--bc`,
  `--lattice` and rotation flags. It refuses a checkpoint whose node count or geometry doesn't match.
  Several further restart faults are fixed:
  - Driven and fixed atoms were read at their α = 0 positions inside the solver. The specimen unloaded
    completely.
  - Checkpoints now save the re-anchored Stage-1 reference positions in an optional `REF_NODES` section,
    which older builds ignore. Before, the reference was rebuilt from the bare lattice.
  - The loader's placeholders for inactive elements made remeshed elements inherit a zero reference area.
  - The last load step could be dropped through floating-point truncation.
  - The Stage-1 → Stage-2 handover no longer runs again on restart.

  Restarts now reproduce the uninterrupted run bit for bit on the tension, polycrystal, disk,
  dead-load and remeshing-avalanche tests. A force grip matches to the solver tolerance.
- Checkpoint restart re-applied the full α as the first increment (curve restarted from zero).
- Stage-2 re-anchoring replaced the specimen frame with the carved bounding box, shrinking disk radius,
  polygon and cutouts for non-rectangular shapes.
- Voronoi seeds inside a central hole were reported as relocated but never moved.

### GUI
- Stage-1 (relaxation) BCs stay **linked** to the loading BCs while they are empty or equal (ignoring
  rule ids); once edited they are **Custom** and are kept. The Protocol panel shows the state with a
  Re-link button. Engineering presets still reset Stage 1; bicrystal runs are never auto-linked.
- Picking Disk in the shape picker swaps plain tension/shear grips to a radial rim load, with an Undo
  notice. Loading a config or preset never changes rules, and compression (diametral) grips are kept.
- The Torsion chip needs an area cutout with a fixed `hole_surfaces` rule; "Make Annulus + Torsion"
  sets this up in one click (centred hole r = 0.35 R).
