#pragma once

#include <string>
#include "experiments/stability_monitor.h"

// Stripped-down Conti-Zanzotto crystal setup that installs a single
// edge dislocation (Stroh for square lattice, Volterra for triangular lattice),
// relaxes once with L-BFGS (no loading, no noise)
void single_dislocation_study(int caller_id, int nx, int ny,
                              const std::string &lattice_type = "square",
                              int potential_mode = 3);

// Load and dipole options of the cylinder study (defaults: single dislocation, relaxation only)
struct CylinderLoadOptions {
  // "none": relax once (as before); "strain": the frozen ring is sheared by F12 = gamma per step, with
  // F11, F22 such that its Cauchy-Born normal stresses vanish; "stress": the ring's homogeneous F_bar is
  // chosen per step so that its Cauchy-Born Cauchy stress is (sigma_xx, sigma_xy, sigma_yy) = (0, tau, 0)
  std::string mode = "none";
  double alpha_start = 0.0; // schedule of gamma or tau (alpha_end < alpha_start loads in reverse)
  double alpha_end = 0.0;
  double step = 0.0;
  int save_vtk_every = 50;
  // Remesh after each load step's relaxation (when remeshing is on). false keeps the initial mesh (the
  // Delaunay mesh of the dislocated crystal) for the whole load: the core can then move only through the
  // lattice energy at fixed connectivity, which isolates an elastic instability from a mesh reconnection.
  bool step_remesh = true;
  // > 0: load cycle (the Peierls protocol of Fig. 8 in the MTM cylinder paper): load until the core has glided
  // cycle_b Burgers vectors, then reverse until it is back at its start; |alpha_end| caps the load
  double cycle_b = 0.0;
  // Normal stresses (sigma_xx, sigma_yy) of the ring's Cauchy-Born stress, in both load modes (applied at the
  // first load step)
  double sigma_xx = 0.0;
  double sigma_yy = 0.0;
  // Dipole: +b and -b cores separated by dipole_dx lattice spacings along x and dipole_dy rows along y
  // (both 0: a single dislocation)
  int dipole_dx = 0;
  int dipole_dy = 0;
  bool loading() const { return mode == "strain" || mode == "stress"; }
  bool dipole() const { return dipole_dx != 0 || dipole_dy != 0; }
};

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
                                            double bulk_K = 5.0,
                                            bool relax_macro_stress = false,
                                            double stress_tol = 1e-5,
                                            int max_macro_iters = 30,
                                            const StabilityMonitorOptions &stability_options = {},
                                            const CylinderLoadOptions &load = {});
