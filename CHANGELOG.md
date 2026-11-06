# Changelog

## Unreleased (feature/gui-onboarding-redesign)

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
