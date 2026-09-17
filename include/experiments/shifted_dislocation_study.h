#pragma once

#include <vector>

/**
 * Runs a multi-shift dislocation study without remeshing:
 * For each shift count s in {0, 1, 2, 3, 4, 5}:
 * 1. Generates an nx x ny square lattice.
 * 2. Meshes once on the pristine reference lattice (fixed topology, no remeshing).
 * 3. Shifts the upper half of the crystal (y > y_mid) horizontally by s * h.
 * 4. Applies analytical Volterra edge dislocation field centered at (x0, y0).
 * 5. Minimizes energy using preconditioned L-BFGS without remeshing.
 * 6. Writes configuration VTK, defect analysis, and logs.
 */
void run_shifted_dislocation_study(int caller_id, int nx, int ny,
                                   const std::vector<int> &shift_counts = {0, 1, 2, 3, 4, 5},
                                   bool use_cylinder = false,
                                   double R_free = 20.0);
