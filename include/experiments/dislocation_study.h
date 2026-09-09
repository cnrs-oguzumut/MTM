#pragma once

// Stripped-down Conti-Zanzotto square-crystal setup that installs a single
// edge dislocation, relaxes once with L-BFGS (no loading, no noise, no
// remeshing), and writes the initial and relaxed states to VTK.
void single_dislocation_study(int caller_id, int nx, int ny);
