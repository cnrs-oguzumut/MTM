#pragma once

#include <string>

// Stripped-down Conti-Zanzotto crystal setup that installs a single
// edge dislocation (Stroh for square lattice, Volterra for triangular lattice),
// relaxes once with L-BFGS (no loading, no noise)
void single_dislocation_study(int caller_id, int nx, int ny,
                              const std::string &lattice_type = "square",
                              int potential_mode = 3);

// Installs an analytical edge dislocation solution (Stroh for square lattice,
// Volterra for triangular lattice) on an nx x ny lattice, freezes all atoms
// outside radius R_free as Dirichlet boundary conditions, and relaxes the atoms
// within R_free around the dislocation core.
void single_dislocation_cylinder_relaxation(int caller_id, int nx, int ny,
                                            double R_free = 20.0,
                                            bool enable_remeshing = false,
                                            bool export_full_mesh = true,
                                            const std::string &lattice_type = "square",
                                            int potential_mode = 3,
                                            double beta = -0.25,
                                            double bulk_K = 5.0);

