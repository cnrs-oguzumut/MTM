#pragma once
// Stability monitor for the quasi-static loading: lowest eigenvalues of the stiffness
// (Hessian) K at relaxed POST states, computed during the run with lowest_stiffness_modes
// (analytic K + Cholesky shift-invert Lanczos, translations projected out).
//
//   * every N-th load step (--eig-every=N), and every step while lambda_min is below
//     refine * (lambda_min right after the last avalanche) (--eig-refine),
//   * at each avalanche (--eig-at-avalanche=1): at POST(i-1), the last stable state before
//     the avalanche (kept in memory), and at POST(i), the state after it.
//
// Output (in the run directory):
//   eigen_log.csv                       step, alpha, trigger, lowest eigenvalues
//   eigen_modes/soft_modes_XXXXX.vtk    lowest nontrivial modes at POST(i-1) of each
//                                       avalanche; XXXXX = file id of the PRE-avalanche
//                                       configuration in vtk_output/. Rigid translations are
//                                       never written.

#include <Eigen/Dense>
#include <cstddef>
#include <fstream>
#include <vector>

#include "../geometry/Point2D.h"
#include "../mesh/ElementTriangle2D.h"
#include "../optimization/LatticeOptimizer.h"
#include "../optimization/PreconditionedLBFGS.h"
#include "../optimization/StiffnessSpectrum.h"

struct StabilityMonitorOptions {
  int every = 0;            // eigenvalues at every N-th POST state (0 = off)
  bool at_avalanche = false; // POST(i-1) and POST(i) of every avalanche
  int modes = 5;            // eigenvalues per computation
  int vectors = 2;          // soft modes written at POST(i-1) of each avalanche
  double refine = 0.2;      // every step while lambda_min < refine * lambda_ref (0 = off)
  bool enabled() const { return every > 0 || at_avalanche; }
};

// Options used by every StabilityMonitor created afterwards (default: off).
void configure_stability_monitor(const StabilityMonitorOptions &options);
const StabilityMonitorOptions &stability_monitor_options();

class StabilityMonitor {
public:
  StabilityMonitor(); // takes the configured options

  // Call once per load step, after relaxation, remeshing and the avalanche decision.
  //   post        UserData of the relaxed state (points, mesh, F_ext of this step)
  //   prev_mesh   mesh at the start of this step (= mesh of POST(step-1)), used at avalanches
  //   pre_file_id file id of the PRE-avalanche configuration written for this avalanche
  void end_of_step(int step, double alpha, UserData &post, bool avalanche,
                   const std::vector<ElementTriangle2D> &prev_mesh,
                   const std::vector<size_t> &prev_active, int pre_file_id);

private:
  StiffnessSpectrum compute(UserData &state, bool with_vectors);
  void log(int step, double alpha, const char *trigger, const StiffnessSpectrum &s,
           double seconds);
  void write_soft_modes(int file_id, int step, double alpha, const StiffnessSpectrum &s,
                        const std::vector<Point2D> &points,
                        const std::vector<ElementTriangle2D> &mesh,
                        const std::vector<size_t> &active,
                        const std::vector<std::pair<int, int>> &full_mapping);
  void update_refinement(const StiffnessSpectrum &s);

  StabilityMonitorOptions options_;
  FastStiffnessAssembler assembler_;
  std::ofstream csv_;

  // POST(step-1), kept for the avalanche case
  int prev_step_ = -1;
  double prev_alpha_ = 0.0;
  Eigen::Matrix2d prev_F_ = Eigen::Matrix2d::Identity();
  std::vector<Point2D> prev_points_;
  bool prev_computed_ = false;
  StiffnessSpectrum prev_spectrum_; // valid if prev_computed_

  double lambda_ref_ = -1.0; // lambda_min of the first computation after an avalanche
  bool refining_ = false;
};
