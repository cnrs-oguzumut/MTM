#include "../../include/experiments/stability_monitor.h"

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
  csv_.open("eigen_log.csv");
  csv_ << "Iteration,Alpha,Trigger,FailedShifts,Seconds";
  for (int k = 1; k <= options_.modes; k++) csv_ << ",lambda_" << k;
  csv_ << "\n" << std::scientific << std::setprecision(12);
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
  return lowest_stiffness_modes(K, opt);
}

void StabilityMonitor::log(int step, double alpha, const char *trigger,
                           const StiffnessSpectrum &s, double seconds) {
  csv_ << step << "," << alpha << "," << trigger << "," << s.failed_shifts << "," << seconds;
  for (int k = 0; k < options_.modes; k++) {
    csv_ << ",";
    if (k < s.num_computed) csv_ << s.eigenvalues(k);
  }
  csv_ << "\n";
  csv_.flush();
  std::cout << "Stability monitor [" << trigger << "] step " << step << ": lambda_min = "
            << (s.num_computed ? s.eigenvalues(0) : std::numeric_limits<double>::quiet_NaN())
            << (s.failed_shifts ? " (UNSTABLE: eigenvalue below the first trial shift)" : "")
            << ", " << seconds << " s" << std::endl;
}

void StabilityMonitor::update_refinement(int step, const StiffnessSpectrum &s) {
  if (s.num_computed == 0) return;
  const double lambda_min = s.eigenvalues(0);
  if (lambda_ref_ < 0.0) lambda_ref_ = lambda_min;
  const bool below = options_.refine > 0.0 && lambda_ref_ > 0.0 &&
                     lambda_min < options_.refine * lambda_ref_;

  // Saddle-node: lambda^2 ~ s (alpha_c - alpha). Extrapolate lambda^2 linearly through the
  // previous point of this branch to estimate the steps left until it vanishes.
  bool approaching = false;
  if (options_.refine_ahead > 0.0 && last_step_ >= 0 && step > last_step_ &&
      lambda_min < last_lambda_) {
    const double l2 = lambda_min * lambda_min;
    const double drop_per_step =
        (last_lambda_ * last_lambda_ - l2) / static_cast<double>(step - last_step_);
    approaching = l2 / drop_per_step < options_.refine_ahead * options_.every;
  }
  refining_ = options_.every > 0 && (below || approaching || lambda_min <= 0.0);
  last_step_ = step;
  last_lambda_ = lambda_min;
}

void StabilityMonitor::start_branch() {
  lambda_ref_ = -1.0;
  refining_ = false;
  last_step_ = -1;
}

void StabilityMonitor::write_soft_modes(int file_id, int step, double alpha,
                                        const StiffnessSpectrum &s,
                                        const std::vector<Point2D> &points,
                                        const std::vector<ElementTriangle2D> &mesh,
                                        const std::vector<size_t> &active,
                                        const std::vector<std::pair<int, int>> &full_mapping) {
  // Lowest nontrivial modes; translations are projected out of the solver under PBC, the
  // check below also excludes them if they were ever returned.
  std::vector<int> modes;
  for (int k = 0; k < s.num_computed && (int)modes.size() < options_.vectors; k++) {
    if (s.eigenvectors.cols() > k && translation_fraction(s.eigenvectors.col(k)) < 0.5)
      modes.push_back(k);
  }
  if (modes.empty()) return;

  std::filesystem::create_directories("eigen_modes");
  std::ostringstream name;
  name << "eigen_modes/soft_modes_" << std::setw(5) << std::setfill('0') << file_id << ".vtk";
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
  file << "Soft modes before avalanche: step " << step << " alpha " << alpha << " lambda";
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
      if (dof < 0) continue;
      d[j] = Eigen::Vector2d(v(dof), v(dof + n_free));
      amax = std::max(amax, d[j].norm());
    }
    file << "VECTORS mode_" << r + 1 << " double\n";
    for (const auto &dj : d) file << dj.x() / amax << " " << dj.y() / amax << " 0\n";
    file << "SCALARS mode_" << r + 1 << "_amplitude double 1\nLOOKUP_TABLE default\n";
    for (const auto &dj : d) file << dj.norm() / amax << "\n";
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
                                   const std::vector<size_t> &prev_active, int pre_file_id) {
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
                      st.F, post.interior_mapping, post.full_mapping, active, false);
        const auto t0 = Clock::now();
        s = compute(kept, vectors);
        if (!st.computed)
          log(st.step, st.alpha, avalanche && last ? "before_avalanche" : "before_instability", s,
              seconds(Clock::now() - t0).count());
        st.computed = true;
      }
      if (vectors)
        write_soft_modes(pre_file_id, st.step, st.alpha, s, st.points, prev_mesh, prev_active,
                         post.full_mapping);
    }
    start_branch();
    branch_start_ = step;
  }

  // 2) The current state: after a saved avalanche, on the grid, or while refining
  const char *trigger = nullptr;
  if (avalanche && options_.at_avalanche) trigger = "after_avalanche";
  else if (options_.every > 0 && step % options_.every == 0) trigger = "grid";
  else if (options_.every > 0 && refining_) trigger = "refine";
  bool computed = false;
  if (trigger) {
    const auto t0 = Clock::now();
    // Vectors are kept so that they can be written if the next step is an avalanche.
    last_spectrum_ = compute(post, options_.at_avalanche && options_.vectors > 0);
    last_spectrum_step_ = step;
    log(step, alpha, trigger, last_spectrum_, seconds(Clock::now() - t0).count());
    update_refinement(step, last_spectrum_);
    computed = true;
  }

  // 3) Keep POST(step)
  const int keep = std::max(options_.retro, options_.at_avalanche ? 1 : 0);
  if (keep > 0) {
    kept_.push_back({step, alpha, post.F_external, post.points,
                     mesh_hash(post.elements, post.active_elements), computed});
    while (static_cast<int>(kept_.size()) > keep) kept_.pop_front();
  }
}
