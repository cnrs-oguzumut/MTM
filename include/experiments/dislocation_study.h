#pragma once

// Stripped-down Conti-Zanzotto square-crystal setup that installs a single
// edge dislocation, relaxes once with L-BFGS (no loading, no noise, no
void single_dislocation_study(int caller_id, int nx, int ny);

// Installs an analytical Volterra edge dislocation solution on an nx x ny lattice,
// freezes all atoms outside radius R_free as Dirichlet boundary conditions,
// and relaxes the atoms within R_free around the dislocation core.
void single_dislocation_cylinder_relaxation(int caller_id, int nx, int ny,
                                           double R_free = 20.0,
                                           bool enable_remeshing = false,
                                           bool export_full_mesh = true);
