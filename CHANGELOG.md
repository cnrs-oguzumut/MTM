# Changelog

## Unreleased (feature/gui-onboarding-redesign)

### Fixed — cylinder load stopped at the first step with --potential-mode=2
- The ring's deformation is found by Newton iteration to a stress residual of 1e-11. The polynomial
  Lennard-Jones of mode 2 has a rounding-noise floor of about 1e-10, so the solve never converged and the
  load stopped at step 1 ("Could not find the ring deformation"). After four iterations a residual below
  1e-9 is now accepted.
- A mode 2 load cycle now completes: triangular, R_free = 35, first jump at τ ∈ (0.062, 0.063]. Other
  potentials change at the relaxation-tolerance level (~1e-6).

### Added — Peierls protocol for the dislocation cylinder: stress against dislocation position
- `--cyl-cycle=N` loads until the dislocation has glided N Burgers vectors, then reverses the load until it is
  back at its start (the protocol of Fig. 8 in the MTM cylinder paper); `--alpha-end` then only caps |α|.
- `--cyl-sxx`, `--cyl-syy` set the normal stresses of the ring's Cauchy–Born stress, in both load modes,
  applied at the first load step (e.g. σ_yy = −0.1 for the hard-core case).
- The load log adds a `Burgers` column (b), the unit of the dislocation position x_d / b.
- GUI:
  - The Volterra load chart has an x-axis choice "x_d / σ / α": dislocation position x_d / b, the ring's
    applied stress, or the load.
  - It opens on the shear stress inside the free region (Cauchy σ_xy, averaged over the free disc) against
    x_d / b, as in Fig. 8a.
  - The ring's Cauchy–Born stress is available as "Applied σ_xy (frozen ring)".
  - The Volterra card has the cycle (default 2 b) and the normal stresses.
- Test: triangular, β = 4, R_free = 30, no remeshing, cycle to 2b. The core hops one b at a time, forward and
  back, and the free-region stress traces the hysteresis loop. Within a valley it rises with x_d and drops at
  each hop.

### Added — GUI: core position against applied stress
- In a Volterra cylinder load run, the "Core glide x_c − x₀" curve can be plotted against the far-field shear
  stress σ_xy applied by the frozen ring (`AppliedSigmaXY`; "x: σ", the default), also under strain control,
  or against the load α ("x: α").
- A dashed orange curve shows the core position after the fixed-mesh relaxation, before that step's
  remeshing (`CoreShiftRelaxed`). A jump only on the solid curve came from the reconnection; a jump on the
  dashed curve is an instability of the lattice energy.

### Added — cylinder load diagnostic: jump in the relaxation, or after remeshing?
- `--cyl-step-remesh=0` keeps the initial mesh (the Delaunay mesh of the dislocated crystal) for the whole
  load, so the core moves only through the lattice energy at fixed connectivity. The GUI has a "Remesh after
  each load step" switch.
- The load log adds the core position after the fixed-connectivity relaxation, before any remeshing of the
  step (`CoreShiftRelaxed`, `Core2ShiftRelaxed`), and `JumpBeforeRemesh`. The console marks each jump "in the
  relaxation" or "after remeshing".
- Triangular, β = 4, R_free = 30, under stress control:
  - With remeshing, the jump at τ = 0.01124 happens only after the remesh: before it the core was at 0.534a,
    still in its valley.
  - Keeping the initial mesh, the core shifts smoothly past 0.53a (1.59a at τ = 0.0148). It becomes unstable
    in the relaxation at τ ∈ (0.015, 0.0152], with the core mode's eigenvalue falling from 0.021 to 0.010
    just before, and glides about 18a.
  - So remeshing pre-empts the elastic instability of the fixed mesh.

### Added — load loop and dipole option for the dislocation cylinder (Peierls stress)
- `--volterra` (single dislocation in a cylinder, analytic field frozen outside R_free) gains a load loop:
  - `--cyl-load=strain` shears the frozen ring by F₁₂ = γ (F₁₁, F₂₂ such that its Cauchy–Born normal stress
    vanishes);
  - `--cyl-load=stress` sets the ring's homogeneous F̄ so that its Cauchy–Born Cauchy stress is
    (σ_xx, σ_xy, σ_yy) = (0, τ, 0);
  - the schedule comes from `--alpha-start` (default 0 here), `--alpha-end` and `--step-size`, and a
    negative end loads in reverse.
- Each step relaxes the core with the ring fixed (and remeshes when enabled). The core position is where the
  slip across its glide plane equals b/2, interpolated below one lattice spacing, in the unsheared frame. A
  step that moves a core more than a/2 is flagged as a jump.
- energy_stress_log.csv gains `AppliedSigmaXY, Fbar11, Fbar12, Fbar22, Core1X, Core1Y, [Core2X, Core2Y],
  CoreShift, [Core2Shift], CoreJump`.
- The stability monitor (`--eig-every`, …) runs at every step. The run stops when a core reaches 3a from the
  ring.
- `--dipole-dx=N --dipole-dy=N` inserts a +b / −b pair separated by N lattice spacings and N rows.
- A loaded run is refused when the frozen ring is thinner than 3a, since free surfaces would bias the result.
  A dipole that annihilates during the initial relaxation is reported and nothing is loaded.
- Without `--cyl-load` the cylinder output is unchanged: byte-identical on the square lattice, within the
  relaxation tolerance (~1e-6) on the triangular one.
- GUI: load picker (None / Strain γ / Stress τ), End and Step, a dipole toggle with dx and dy (default
  dy = 4: a bound dipole), and a "Load Curve" view of the core glide against α. Optional log columns are now
  read by header name.
- Findings (Zanzotto model):
  - Triangular, β = 4: the core shifts elastically within its valley, then jumps. (τ₊ − τ₋)/2 = 0.0116,
    0.01095 and 0.0101 for R_free = 20, 30 and 60 (≈ 0.0095 extrapolated).
  - The ring's pull toward the centre raises |τ| in both directions, so it does not cancel in
    (τ₊ − τ₋)/2: extrapolate in R_free.
  - With remeshing on, the jump coincides with the first accepted mesh reconnection, and the core mode's
    eigenvalue stays positive. Without remeshing, the core glides from τ ≈ 0.0003.
  - Square, β = −0.25: smooth glide, no threshold above 5e-5.

### Fixed — "Relax box stress to zero" now works together with a relaxation at prescribed strain
- With `--relax-strain` set, `--relax-stress=1` had no effect: there was no α = 0 step, so the stress
  schedule started at the relaxed stress σ₀. σ₀ can be negative, e.g. −0.0048 in a 150×150 cell at
  F₁₂ = 0.14, and the run then began with reverse loading and reverse avalanches.
- Now step 0 stays the relaxation at F₁₂, step 1 relaxes the cell to σ = 0 in all three components, and the
  schedule rises from 0 (dα, 2dα, …). The zero-stress step is not a load step: it cannot be flagged as an
  avalanche, it gives no stress predictor history, and it sets the strain origin of the first strain jump.
  Unchecked, the schedule still starts from σ₀ as before.
- The GUI labels the option "Relax box stress to zero after the relaxation step" and shows "Start stress: 0"
  in that case. The limit-load message after the first stress step reports σ₀ instead of the strain as
  the last equilibrium load.
- Checked on 40×40 cells, seeds 42 and 7 (σ₀ = +0.012, +0.040): step 1 has σ ≈ 1e-9 in all components, and
  loading rises from 0.

### Added — GUI: x-axis zoom on the results chart
- Drag across the stress / energy / λ_min chart to zoom the x axis to that range. The y axis refits to the
  points inside it. A selection that reaches the newest point keeps following the live data. Double-click
  the chart, click the "Zoom ✕" chip, or pick a scale mode to show all points again. Useful when the
  relaxation step (e.g. F₁₂ = 0.14 before stress loading) squeezes the loading curve into a corner.

### Added — avalanche-size columns under force / stress control
- energy_stress_log.csv of force- and stress-controlled runs gains three columns at the end: `StrainJump`
  (Δε of the step), `CascadeEnthalpyChange` (ΔΠ of the reconnection cascade at fixed load,
  Cauchy-corrected; the value the console prints for each avalanche), and `Avalanche` (1 when the step is
  flagged). They are 0 on steps that are not force-controlled. Existing columns keep their positions.
- 40×40 Cauchy-shear check: the 5 flagged avalanches match fractional_drops.csv and the console, StrainJump
  equals the difference of EffectiveStrain on every row, and ΔΠ < 0 on each avalanche.

### Changed — α is a stress throughout a relaxation-at-strain stress run
- With `--relax-strain`, step 0 holds the cell at F₁₂ = relax_strain, and the loading after it is a stress
  schedule. The logs wrote the strain (e.g. 0.14) as α at step 0 and the target stress (e.g. 0.0069) from
  step 1 on. Step 0 now logs its relaxed stress σ₀, the start of the schedule. This applies to
  energy_stress_log.csv, eigen_log.csv, fractional_drops.csv, the dislocation log and configuration_00000.vtk.
  The strain is in EffectiveStrain / Fbar12 as before.

### Fixed — eigen_log.csv cell-stability columns at the relaxation step
- The Strain / lambda_box / box_tangent cells of the fixed-strain relaxation step held garbage instead of
  being empty. The values used NaN as "not computed", which `-ffast-math` (assumes no NaNs) does not keep.
  They now carry explicit flags. Refinement could also have been triggered by that garbage.

### Changed — GUI: one place for the boundary regime and loading
- The "Boundary Regime & Mechanical Loading" card is shown only in step 2 (Shape & BCs); it was repeated in
  steps 1 and 3. The canvas toolbar's "Surface Grips / Custom BCs" switch is removed, and the canvas
  regime/loading pill is a status label (hover: edit in step 2).
- The "Carve Macroscopic Shape & Custom BCs" switch is removed: the PBC / Surface choice decides it. In PBC,
  step 2 says there is no specimen to carve. "Show Specimen on Canvas" no longer switches the regime to
  surface BCs as a side effect.

### Docs — derivations in stress_and_traction_control.tex
- Step-by-step derivations, each in index and matrix form, of:
  - the cell derivatives (7)–(8) and the rotational balance;
  - the scaled gradient;
  - the Hessian blocks K^FF, K^Fw, K^ww, with the geometric term shown to be the atomic force;
  - the relaxed tangent K_eff and why the homogeneous lattice has no relaxation term;
  - the Cauchy targets via cof F̄ and the load-stiffness matrices 𝕃 (and why only the hydrostatic target
    has a potential);
  - the x-form of K_eff and the load-path tangent.
- Fixed Eq. (P21), the rotational balance in the upper-triangular gauge. It now reads
  P̄₂₁F̄₁₁ + P̄₂₂F̄₁₂ = P̄₁₂F̄₂₂; the old form, P̄₂₁F̄₁₁ = P̄₁₁F̄₁₂ + P̄₁₂F̄₂₂, was wrong. Its only use, the nominal-shear
  example with P̄₂₂ = 0, is unaffected, and the code does not use it.

### Added — free-cell stability under PBC stress control
- Under stress control the stability monitor (`--eig-every`, `--eig-at-avalanche`) now also evaluates the
  stability of the free cell. That is λ_box, the lowest eigenvalue (real part) of the relaxed cell Jacobian
  J = K_eff − ∂Pᵗ/∂F̄ on the free components. K_eff = (K_FF − K_xFᵀ K⁺ K_xF)/A₀ is the relaxed (Schur
  complement) cell stiffness, from the unprojected element tangents. ∂Pᵗ/∂F̄ is the Cauchy follower-load
  stiffness (zero for nominal presets). λ_box reaches 0 at the limit load, where the fixed-cell λ_min(K) can
  still be positive.
- eigen_log.csv gains three columns after the lambda_k columns under stress control: `Strain` (as
  EffectiveStrain), `lambda_box`, and `box_tangent` = dα/dε along the load path, the slope of the
  load–strain curve. Their cells are empty at the strain-controlled relaxation step. Refinement (`--eig-refine`,
  `--eig-refine-ahead`) also follows λ_box. Runs without stress control write the same file as before
  (byte-identical in a strain-control check).
- Checks: a perfect 20×20 cell under nominal shear has λ_box 2.18 → 0.31 from α = 0.25 to 0.315, with the limit
  load in (0.315, 0.32]. Extrapolating λ_box² gives α_c ≈ 0.317. In a defected 40×40 cell under Cauchy shear,
  box_tangent matches the secant slope of the logged curve to < 1 % on the elastic steps.
- GUI: under stress control the λ_min chart has a Cell / Atoms picker (λ_box or λ_min(K)) and the ε / α axis
  picker. The x-axis label no longer calls a stress load "Shear Strain". The hover readout shows ε, λ_box
  and dα/dε.
- docs/stress_and_traction_control.tex: new subsection "Stability of the relaxed states".

### Fixed — review of the stress-control and traction work
- Cauchy control: a rejected reconnection pass restored the cell but not the nominal target, so the rest of
  the step logged H and E_el with the trial target. It now restores both.
- Cauchy control: reconnection acceptance, and the avalanche test "cascade lowered the enthalpy", compared H
  values taken with different nominal targets. They now use ΔΠ = ΔH + A₀(Pᵗ₁ − Pᵗ₀):(F̄_mid − I), the
  trapezoid rule for the work of the target change; nominal control is unchanged. In a 40×40 Cauchy-shear
  test the avalanches fall on the same steps, a few cascade passes are now rejected, and the final F̄₁₂
  moves by 1e-4.
- Cauchy control: `PreEnergy` (and `EnergyChange`) was the elastic energy, because the target was set only
  inside the relaxation. It is now the enthalpy at the predicted cell, like `PostEnergy`.
- `--relax-stress=1` was ignored when α_start ≠ 0 on a plain periodic cell. The zero-stress relaxation
  now runs at α = 0 first.
- Plain periodic stress runs starting at α = 0 applied the initial noise at that zero-load step, where it
  relaxed away. It is now applied at the first loaded step, as after a Stage 1. Strain-control runs are
  unchanged.
- Restart:
  - checkpoints carry the Stage-1 cell (`F_stage1`), so a resumed stress run keeps its strain origin;
  - the first resumed step no longer sees a bogus strain increment (a mid-run restart now reproduces the
    original EffectiveStrain exactly);
  - checkpoints written at the relaxation step before loading (Stage 1, relaxation at a prescribed shear
    strain, zero-stress relaxation) are marked, and a restart from them is refused with a message, since
    the loading setup (re-anchoring, first-step noise) would be skipped.
  - Older checkpoints load as before.
- Stress preconditioner:
  - with `--precond-refresh > 1` the factor is reused only on the same mesh and cell components;
  - `--precond-refresh-iters` is honoured;
  - `MTM_STRESS_TRACE` prints the real iteration count on the preconditioned path.
- Dead loads:
  - supports that leave a rigid-body mode free (rollers along one direction, a single row of rollers, a
    single pinned atom) are now checked: a load with a net force along a free translation, or a net moment
    about a free rotation, is refused instead of running away to a false "limit load".
  - The GUI preview weights edge tractions by edge length, checks the torque, and handles roller-only
    supports.
- GUI: the launch button names the PBC stress preset ("PBC Cauchy Shear (σ_xy = α)", …) instead of
  "PBC Simple Shear".
- Not changed: force grips still relax with plain L-BFGS, a performance limit rather than a bug.
- Checks: strain-control output is byte-identical; nominal stress control matches the previous build to the
  solver tolerance (σ ~1e-8).

### Fixed — eigenvalue runs crashed at startup
- `--eig-at-avalanche=` and `--eig-refine-ahead=` were read from the wrong character offset. The GUI always passes
  `--eig-at-avalanche`, so every run with "Compute Lowest Stiffness Eigenvalues" on stopped at startup with
  `stoi: no conversion`.

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
