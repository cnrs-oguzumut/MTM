#include "../../include/experiments/stability_monitor.h"
#include "../../include/geometry/MacroscopicGeometry.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

namespace {

StabilityMonitorOptions g_stability_options;

using Clock = std::chrono::high_resolution_clock;

std::uint64_t mesh_hash(const std::vector<ElementTriangle2D> &elements,
                        const std::vector<size_t> &active) {
  std::uint64_t h = 1469598103934665603ULL ^ active.size();
  for (size_t idx : active)
    for (int a = 0; a < 3; a++)
      h = (h ^ static_cast<std::uint64_t>(elements[idx].getNodeIndex(a) + 1)) * 1099511628211ULL;
  return h;
}

// Refinement test on one stability measure lambda (see StabilityMonitor::update_refinement)
bool approaching_instability(const StabilityMonitorOptions &o, double &ref, int &last_step,
                             double &last, int step, double lambda) {
  if (ref < 0.0) ref = lambda;
  const bool below = o.refine > 0.0 && ref > 0.0 && lambda < o.refine * ref;

  // Saddle-node: lambda^2 ~ s (alpha_c - alpha). Extrapolate lambda^2 linearly through the
  // previous point of this branch to estimate the steps left until it vanishes.
  bool approaching = false;
  if (o.refine_ahead > 0.0 && last_step >= 0 && step > last_step && lambda < last) {
    const double l2 = lambda * lambda;
    const double drop_per_step = (last * last - l2) / static_cast<double>(step - last_step);
    approaching = l2 / drop_per_step < o.refine_ahead * o.every;
  }
  last_step = step;
  last = lambda;
  return below || approaching || lambda <= 0.0;
}

// Fraction of |v|^2 carried by the two uniform translations (solver layout [u..., v...]).
double translation_fraction(const Eigen::VectorXd &v) {
  const Eigen::Index m = v.size() / 2;
  const double mu = v.head(m).mean();
  const double mv = v.tail(m).mean();
  return (mu * mu + mv * mv) * static_cast<double>(m) / std::max(v.squaredNorm(), 1e-300);
}

} // namespace

void configure_stability_monitor(const StabilityMonitorOptions &options) {
  g_stability_options = options;
}

const StabilityMonitorOptions &stability_monitor_options() { return g_stability_options; }

StabilityMonitor::StabilityMonitor() : options_(g_stability_options) {
  if (!options_.enabled()) return;
  csv_.open("eigen_log.csv"); // header at the first row, once set_box_stability had its chance
  csv_ << std::scientific << std::setprecision(12);
}

StiffnessSpectrum StabilityMonitor::compute(UserData &state, bool with_vectors) {
  const int n_vars = static_cast<int>(state.interior_mapping.size());
  alglib::real_1d_array x;
  x.setlength(2 * n_vars);
  map_points_to_solver_array(x, state.points, state.interior_mapping, n_vars);
  for (auto &element : state.elements) element.setExternalDeformation(state.F_external);

  const Eigen::SparseMatrix<double> &K = assembler_.assembleStiffness(x, &state, false);
  StiffnessSpectrumOptions opt;
  opt.n_modes = std::max(options_.modes, with_vectors ? options_.vectors : 0);
  opt.compute_vectors = with_vectors;

  // When Roller (partial DOF) constraints are active, FastStiffnessAssembler places 1.0 on
  // locked scalar diagonals. Lift those isolated diagonals above max_free_diag so stiff
  // potentials (where physical lambda_min > 1.0) never return locked boundary DOFs as soft modes.
  if (state.macro_geom_ptr != nullptr && state.macro_geom_ptr->hasPartialDOFConstraints()) {
    Eigen::SparseMatrix<double> K_adj = K;
    double max_free_diag = 1.0;
    for (int i = 0; i < n_vars; ++i) {
      const int node_idx = state.interior_mapping[i].first;
      const int dof_idx = state.interior_mapping[i].second;
      if (dof_idx < 0 || dof_idx >= n_vars) continue;
      if (!state.macro_geom_ptr->isNodeXConstrained(node_idx))
        max_free_diag = std::max(max_free_diag, std::abs(K_adj.coeff(dof_idx, dof_idx)));
      if (!state.macro_geom_ptr->isNodeYConstrained(node_idx))
        max_free_diag = std::max(max_free_diag, std::abs(K_adj.coeff(n_vars + dof_idx, n_vars + dof_idx)));
    }
    const double locked_diag = std::max(10.0, 2.0 * max_free_diag);
    for (int i = 0; i < n_vars; ++i) {
      const int node_idx = state.interior_mapping[i].first;
      const int dof_idx = state.interior_mapping[i].second;
      if (dof_idx < 0 || dof_idx >= n_vars) continue;
      if (state.macro_geom_ptr->isNodeXConstrained(node_idx))
        K_adj.coeffRef(dof_idx, dof_idx) = locked_diag;
      if (state.macro_geom_ptr->isNodeYConstrained(node_idx))
        K_adj.coeffRef(n_vars + dof_idx, n_vars + dof_idx) = locked_diag;
    }
    return lowest_stiffness_modes(K_adj, opt);
  }

  return lowest_stiffness_modes(K, opt);
}

BoxStabilityPoint StabilityMonitor::compute_box(UserData &state, int step, double alpha) {
  return box_fn_ ? box_fn_(state, step, alpha) : BoxStabilityPoint{};
}

void StabilityMonitor::log(int step, double alpha, const char *trigger,
                           const StiffnessSpectrum &s, const BoxStabilityPoint &box,
                           double seconds) {
  if (!header_written_) {
    csv_ << "Iteration,Alpha,Trigger,FailedShifts,Seconds";
    for (int k = 1; k <= options_.modes; k++) csv_ << ",lambda_" << k;
    if (box_fn_) csv_ << ",Strain,lambda_box,box_tangent";
    csv_ << "\n";
    header_written_ = true;
  }
  csv_ << step << "," << alpha << "," << trigger << "," << s.failed_shifts << "," << seconds;
  for (int k = 0; k < options_.modes; k++) {
    csv_ << ",";
    if (k < s.num_computed) csv_ << s.eigenvalues(k);
  }
  if (box_fn_) {
    csv_ << ",";
    if (box.has_strain) csv_ << box.strain;
    csv_ << ",";
    if (box.has_lambda) csv_ << box.lambda_min;
    csv_ << ",";
    if (box.has_tangent) csv_ << box.tangent;
  }
  csv_ << "\n";
  csv_.flush();
  std::cout << "Stability monitor [" << trigger << "] step " << step << ": lambda_min = "
            << (s.num_computed ? s.eigenvalues(0) : std::numeric_limits<double>::quiet_NaN())
            << (s.failed_shifts ? " (UNSTABLE: eigenvalue below the first trial shift)" : "");
  if (box.has_lambda)
    std::cout << ", lambda_box = " << box.lambda_min << (box.lambda_min <= 0.0 ? " (UNSTABLE cell)" : "");
  std::cout << ", " << seconds << " s" << std::endl;
}

void StabilityMonitor::update_refinement(int step, const StiffnessSpectrum &s,
                                         const BoxStabilityPoint &box) {
  const bool have_box = box.has_lambda;
  if (s.num_computed == 0 && !have_box) return;
  bool refine = false;
  if (s.num_computed > 0)
    refine = approaching_instability(options_, lambda_ref_, last_step_, last_lambda_, step,
                                     s.eigenvalues(0));
  if (have_box)
    refine = approaching_instability(options_, box_ref_, box_last_step_, box_last_, step,
                                     box.lambda_min) || refine;
  refining_ = options_.every > 0 && refine;
}

void StabilityMonitor::start_branch() {
  lambda_ref_ = -1.0;
  refining_ = false;
  last_step_ = -1;
  box_ref_ = -1.0;
  box_last_step_ = -1;
}

void StabilityMonitor::clear_kept_states() {
  kept_.clear();
  last_spectrum_step_ = -1;
  start_branch();
}

void StabilityMonitor::write_soft_modes(const std::string &file_stem, const char *trigger,
                                        int step, double alpha,
                                        const StiffnessSpectrum &s,
                                        const std::vector<Point2D> &points,
                                        const std::vector<ElementTriangle2D> &mesh,
                                        const std::vector<size_t> &active,
                                        const std::vector<std::pair<int, int>> &full_mapping) {
  // Lowest nontrivial modes; true rigid-body translations have both |lambda| ~ 0 and
  // near-unit translation fraction (> 0.95). Do not use a loose 0.5 cutoff, which falsely
  // rejects bicrystal grain-sliding soft modes where one grain (~50% of atoms) translates!
  std::vector<int> modes;
  for (int k = 0; k < s.num_computed && (int)modes.size() < options_.vectors; k++) {
    if (s.eigenvectors.cols() <= k) continue;
    const bool is_rigid_translation =
        (std::abs(s.eigenvalues(k)) < 1e-8 && translation_fraction(s.eigenvectors.col(k)) > 0.95);
    if (!is_rigid_translation) modes.push_back(k);
  }
  if (modes.empty()) return;

  std::filesystem::create_directories("eigen_modes");
  std::ostringstream name;
  name << "eigen_modes/" << file_stem << ".vtk";
  std::ofstream file(name.str());

  const int n_nodes = static_cast<int>(points.size());
  const int n_free = static_cast<int>(s.eigenvectors.rows() / 2);

  std::vector<size_t> cells; // triangles of the fundamental domain, as in exportSingleModeToVTK
  for (size_t idx : active) {
    const ElementTriangle2D &e = mesh[idx];
    bool fundamental = e.isInitialized();
    for (int a = 0; a < 3 && fundamental; a++)
      if (e.getTranslation(a).norm() > 1e-10) fundamental = false;
    if (fundamental) cells.push_back(idx);
  }

  file << "# vtk DataFile Version 3.0\n";
  file << "Soft modes (" << (trigger ? trigger : "state") << "): step " << step << " alpha " << alpha << " lambda";
  for (int k : modes) file << " " << s.eigenvalues(k);
  file << "\nASCII\nDATASET UNSTRUCTURED_GRID\n";
  file << std::setprecision(10);
  file << "POINTS " << n_nodes << " double\n";
  for (const auto &p : points) file << p.coord.x() << " " << p.coord.y() << " 0\n";
  file << "\nCELLS " << cells.size() << " " << 4 * cells.size() << "\n";
  for (size_t idx : cells)
    file << "3 " << mesh[idx].getNodeIndex(0) << " " << mesh[idx].getNodeIndex(1) << " "
         << mesh[idx].getNodeIndex(2) << "\n";
  file << "\nCELL_TYPES " << cells.size() << "\n";
  for (size_t c = 0; c < cells.size(); c++) file << "5\n";

  // Each mode scaled to max nodal amplitude 1 (the sign of an eigenvector is arbitrary).
  file << "\nPOINT_DATA " << n_nodes << "\n";
  for (size_t r = 0; r < modes.size(); r++) {
    const Eigen::VectorXd v = s.eigenvectors.col(modes[r]);
    std::vector<Eigen::Vector2d> d(n_nodes, Eigen::Vector2d::Zero());
    double amax = 0.0;
    for (int j = 0; j < n_nodes; j++) {
      const int dof = full_mapping[j].second;
      if (dof < 0 || dof >= n_free) continue;
      d[j] = Eigen::Vector2d(v(dof), v(dof + n_free));
      amax = std::max(amax, d[j].norm());
    }
    const double scale = (amax > 1e-30) ? amax : 1.0;
    file << "VECTORS mode_" << r + 1 << " double\n";
    for (const auto &dj : d) file << dj.x() / scale << " " << dj.y() / scale << " 0\n";
    file << "SCALARS mode_" << r + 1 << "_amplitude double 1\nLOOKUP_TABLE default\n";
    for (const auto &dj : d) file << dj.norm() / scale << "\n";
  }
  std::cout << "Stability monitor: soft modes written to " << name.str() << std::endl;
}

bool StabilityMonitor::stress_jump(double alpha, double stress) {
  bool jump = false;
  if (have_prev_ && alpha != prev_alpha_) {
    const double direction = alpha > prev_alpha_ ? 1.0 : -1.0;
    const double increment = direction * (stress - prev_stress_); // > 0 when elastic
    double median = 0.0;
    if (!increments_.empty()) {
      std::vector<double> sorted = increments_;
      std::nth_element(sorted.begin(), sorted.begin() + sorted.size() / 2, sorted.end());
      median = sorted[sorted.size() / 2];
    }
    if (median > 0.0 && increment < -options_.jump_tol * median) {
      jump = true;
    } else if (increment > 0.0) {
      constexpr size_t kWindow = 256;
      if (increments_.size() < kWindow) increments_.push_back(increment);
      else increments_[increment_pos_++ % kWindow] = increment;
    }
  }
  have_prev_ = true;
  prev_alpha_ = alpha;
  prev_stress_ = stress;
  return jump;
}

void StabilityMonitor::end_of_step(int step, double alpha, double post_stress, UserData &post,
                                   bool avalanche,
                                   const std::vector<ElementTriangle2D> &prev_mesh,
                                   const std::vector<size_t> &prev_active, int pre_file_id,
                                   bool is_final_step) {
  if (!options_.enabled()) return;
  using seconds = std::chrono::duration<double>;

  const bool jump = stress_jump(alpha, post_stress);
  const bool instability = step > 0 && (avalanche || jump);
  const bool avalanche_modes = avalanche && options_.at_avalanche && options_.vectors > 0;

  // 1) Before an instability: the kept states of the branch that ends here. Their mesh is
  //    the one this step started with unless a remeshing happened in between.
  if (instability) {
    const int depth = std::max(options_.retro, avalanche && options_.at_avalanche ? 1 : 0);
    const std::uint64_t mesh = mesh_hash(prev_mesh, prev_active);
    for (KeptState &st : kept_) { // oldest first, so the log stays in load order
      if (st.step < step - depth || st.step < branch_start_) continue;
      const bool last = (st.step == step - 1);
      const bool vectors = avalanche_modes && last;
      StiffnessSpectrum s;
      if (st.computed) {
        if (!vectors) continue;
        if (last_spectrum_step_ == st.step &&
            last_spectrum_.eigenvectors.cols() >= options_.vectors)
          s = last_spectrum_;
      }
      if (s.num_computed == 0) {
        if (st.mesh != mesh) continue;
        std::vector<ElementTriangle2D> elements = prev_mesh;
        std::vector<size_t> active = prev_active;
        UserData kept(st.points, elements, post.calculator, post.energy_function,
                      post.derivative_function, post.zero_energy, post.ideal_lattice_parameter,
                      st.F, post.interior_mapping, post.full_mapping, active, false,
                      post.second_derivative_function);
        kept.macro_geom_ptr = post.macro_geom_ptr;
        kept.current_alpha = st.alpha;
        kept.stress_ctrl = post.stress_ctrl;
        const auto t0 = Clock::now();
        s = compute(kept, vectors);
        if (!st.computed) {
          const BoxStabilityPoint box = compute_box(kept, st.step, st.alpha);
          log(st.step, st.alpha, avalanche && last ? "before_avalanche" : "before_instability", s,
              box, seconds(Clock::now() - t0).count());
        }
        st.computed = true;
      }
      if (vectors) {
        std::ostringstream stem;
        stem << "soft_modes_step_" << std::setw(5) << std::setfill('0') << st.step
             << "_pre_avalanche_" << std::setw(5) << std::setfill('0') << std::max(0, pre_file_id);
        write_soft_modes(stem.str(), "before_avalanche", st.step, st.alpha, s, st.points,
                         prev_mesh, prev_active, post.full_mapping);
      }
    }
    start_branch();
    branch_start_ = step;
  }

  // 2) The current state: after a saved avalanche, on the grid, while refining, or at the final step
  const char *trigger = nullptr;
  if (avalanche && options_.at_avalanche) trigger = "after_avalanche";
  else if (options_.every > 0 && step % options_.every == 0) trigger = "grid";
  else if (step == 0) trigger = "initial";
  else if (options_.every > 0 && refining_) trigger = "refine";
  else if (is_final_step && options_.every > 0) trigger = "final";

  bool computed = false;
  if (trigger) {
    const auto t0 = Clock::now();
    const bool write_here = options_.vectors > 0 && std::string(trigger) != "refine";
    const bool want_vectors = options_.vectors > 0 && (write_here || options_.at_avalanche);
    last_spectrum_ = compute(post, want_vectors);
    last_spectrum_step_ = step;
    const BoxStabilityPoint box = compute_box(post, step, alpha);
    log(step, alpha, trigger, last_spectrum_, box, seconds(Clock::now() - t0).count());
    update_refinement(step, last_spectrum_, box);
    computed = true;

    if (write_here && last_spectrum_.num_computed > 0) {
      std::ostringstream stem;
      stem << "soft_modes_step_" << std::setw(5) << std::setfill('0') << step;
      write_soft_modes(stem.str(), trigger, step, alpha, last_spectrum_, post.points,
                       post.elements, post.active_elements, post.full_mapping);
    }
  }

  // 3) Keep POST(step)
  const int keep = std::max(options_.retro, options_.at_avalanche ? 1 : 0);
  if (keep > 0) {
    kept_.push_back({step, alpha, post.F_external, post.points,
                     mesh_hash(post.elements, post.active_elements), computed});
    while (static_cast<int>(kept_.size()) > keep) kept_.pop_front();
  }
}
