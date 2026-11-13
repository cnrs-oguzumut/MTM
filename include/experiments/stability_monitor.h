#pragma once
// Stability monitor for the quasi-static loading: lowest eigenvalues of the stiffness
// (Hessian) K at relaxed POST states, computed during the run with lowest_stiffness_modes
// (analytic K + Cholesky shift-invert Lanczos, translations projected out).
//
//   * every N-th load step (--eig-every=N), and every step while lambda_min is below
//     refine * (lambda_min right after the last avalanche) (--eig-refine), or while the
//     linear extrapolation of lambda_min^2 through the last two points reaches zero within
//     refine_ahead grid intervals (--eig-refine-ahead; near a saddle-node
//     lambda_min^2 ~ s (alpha_c - alpha), so the approach to an instability is resolved),
//   * at each instability (a stress jump or a saved avalanche): the last `retro` relaxed
//     states before it that were not computed yet (--eig-retro=K). They are kept in memory
//     (positions; on an elastic branch the mesh does not change), so the approach to every
//     instability is resolved step by step without computing at every step,
//   * at each saved avalanche (--eig-at-avalanche=1): also the state after it, POST(i), and
//     the soft modes of POST(i-1), the last stable state.
//
// Under PBC stress control the cell is free, and lambda_min(K) (fixed cell) is the strain-control
// criterion only. set_box_stability() adds, at the same states, the lowest eigenvalue of the relaxed
// cell Jacobian (pbc_box_stability), which vanishes at the limit load; refinement then also follows it.
//
// Output (in the run directory):
//   eigen_log.csv                       step, alpha, trigger, lowest eigenvalues (rows of
//                                       retroactive computations come after later steps:
//                                       sort by Iteration); with set_box_stability also
//                                       Strain, lambda_box, box_tangent
//   eigen_modes/soft_modes_XXXXX.vtk    lowest nontrivial modes at POST(i-1) of each
//                                       avalanche; XXXXX = file id of the PRE-avalanche
//                                       configuration in vtk_output/. Rigid translations are
//                                       never written.

#include <Eigen/Dense>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <fstream>
#include <functional>
#include <limits>
#include <utility>
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
  double refine_ahead = 2.0; // every step while lambda^2 extrapolates to 0 within
                             // refine_ahead * every steps (0 = off)
  int retro = 4;            // at each instability, the last `retro` states before it
  double jump_tol = 2.0;    // instability: stress decrease > jump_tol x median elastic
                            // stress increment per step (as in plot_saddle_node.py)
  bool enabled() const { return every > 0 || at_avalanche; }
};

// PBC stress control: the free-cell stability of a relaxed state (see pbc_box_stability). Each value is
// written only when its flag is set (no NaN sentinels: the build uses -ffast-math).
struct BoxStabilityPoint {
  double strain = 0.0;     // as EffectiveStrain in the stress log
  double lambda_min = 0.0; // lowest eigenvalue of the cell Jacobian
  double tangent = 0.0;    // d alpha / d strain along the load path
  bool has_strain = false, has_lambda = false, has_tangent = false;
};
using BoxStabilityFunction = std::function<BoxStabilityPoint(UserData &state, int step, double alpha)>;

// Options used by every StabilityMonitor created afterwards (default: off).
void configure_stability_monitor(const StabilityMonitorOptions &options);
const StabilityMonitorOptions &stability_monitor_options();

class StabilityMonitor {
public:
  StabilityMonitor(); // takes the configured options

  // Adds the Strain, lambda_box, box_tangent columns (call before the first step)
  void set_box_stability(BoxStabilityFunction fn) { box_fn_ = std::move(fn); }

  // Call once per load step, after relaxation, remeshing and the avalanche decision.
  //   post_stress   shear stress of the relaxed state (instability detection)
  //   post          UserData of the relaxed state (points, mesh, F_ext of this step)
  //   avalanche     this step is a saved avalanche
  //   prev_mesh     mesh at the start of this step (= mesh of POST(step-1))
  //   pre_file_id   file id of the PRE-avalanche configuration written for this avalanche
  //   is_final_step true on the final step of the schedule (or Relax Only completion)
  void end_of_step(int step, double alpha, double post_stress, UserData &post, bool avalanche,
                   const std::vector<ElementTriangle2D> &prev_mesh,
                   const std::vector<size_t> &prev_active, int pre_file_id,
                   bool is_final_step = false);

  // Drop cached states when boundary DOF mappings change (e.g. Stage 1 -> Stage 2 BC handover)
  void clear_kept_states();

private:
  StiffnessSpectrum compute(UserData &state, bool with_vectors);
  BoxStabilityPoint compute_box(UserData &state, int step, double alpha);
  void log(int step, double alpha, const char *trigger, const StiffnessSpectrum &s,
           const BoxStabilityPoint &box, double seconds);
  void write_soft_modes(const std::string &file_stem, const char *trigger, int step, double alpha,
                        const StiffnessSpectrum &s,
                        const std::vector<Point2D> &points,
                        const std::vector<ElementTriangle2D> &mesh,
                        const std::vector<size_t> &active,
                        const std::vector<std::pair<int, int>> &full_mapping);
  void update_refinement(int step, const StiffnessSpectrum &s, const BoxStabilityPoint &box);
  void start_branch();

  StabilityMonitorOptions options_;
  FastStiffnessAssembler assembler_;
  std::ofstream csv_;
  bool header_written_ = false;
  BoxStabilityFunction box_fn_;

  bool stress_jump(double alpha, double stress);

  // Last relaxed states, for the computations after an instability
  struct KeptState {
    int step;
    double alpha;
    Eigen::Matrix2d F;
    std::vector<Point2D> points;
    std::uint64_t mesh; // connectivity hash
    bool computed;
  };
  std::deque<KeptState> kept_;
  int branch_start_ = 0;             // first step after the last instability
  StiffnessSpectrum last_spectrum_;  // most recent computation of the current state
  int last_spectrum_step_ = -1;

  // Instability detection from the stress
  bool have_prev_ = false;
  double prev_alpha_ = 0.0, prev_stress_ = 0.0;
  std::vector<double> increments_; // recent elastic stress increments (loading direction)
  size_t increment_pos_ = 0;

  double lambda_ref_ = -1.0; // lambda_min of the first computation after an avalanche
  bool refining_ = false;
  int last_step_ = -1;       // last computation on the current branch (for extrapolation)
  double last_lambda_ = 0.0;
  // The same for lambda_box
  double box_ref_ = -1.0;
  int box_last_step_ = -1;
  double box_last_ = 0.0;
};
