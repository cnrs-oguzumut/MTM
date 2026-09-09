#pragma once

// Stripped-down Conti-Zanzotto square-crystal setup that rigidly shifts the
// upper half of the crystal by one lattice spacing in x, relaxes once with
// L-BFGS (no loading, no noise, remeshing off by default), and writes the
// initial and relaxed states to VTK.
void shifted_upper_crystal_study(int caller_id, int nx, int ny);
