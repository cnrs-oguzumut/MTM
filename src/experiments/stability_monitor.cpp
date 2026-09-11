#include "../../include/experiments/stability_monitor.h"

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

void StabilityMonitor::update_refinement(const StiffnessSpectrum &s) {
  if (s.num_computed == 0) return;
  const double lambda_min = s.eigenvalues(0);
  if (lambda_ref_ < 0.0) lambda_ref_ = lambda_min;
  refining_ = options_.every > 0 && options_.refine > 0.0 && lambda_ref_ > 0.0 &&
              lambda_min < options_.refine * lambda_ref_;
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

void StabilityMonitor::end_of_step(int step, double alpha, UserData &post, bool avalanche,
                                   const std::vector<ElementTriangle2D> &prev_mesh,
                                   const std::vector<size_t> &prev_active, int pre_file_id) {
  if (!options_.enabled()) return;

  bool computed = false;
  StiffnessSpectrum current;

  if (avalanche && options_.at_avalanche && prev_step_ == step - 1) {
    // POST(step-1): last stable state, reconstructed from the kept positions and the mesh
    // this step started with.
    StiffnessSpectrum before;
    if (prev_computed_ && prev_spectrum_.eigenvectors.cols() >= options_.vectors) {
      before = prev_spectrum_;
    } else {
      std::vector<ElementTriangle2D> mesh = prev_mesh;
      std::vector<size_t> active = prev_active;
      UserData prev(prev_points_, mesh, post.calculator, post.energy_function,
                    post.derivative_function, post.zero_energy, post.ideal_lattice_parameter,
                    prev_F_, post.interior_mapping, post.full_mapping, active, false);
      const auto t0 = Clock::now();
      before = compute(prev, options_.vectors > 0);
      if (!prev_computed_)
        log(prev_step_, prev_alpha_, "before_avalanche", before,
            std::chrono::duration<double>(Clock::now() - t0).count());
    }
    if (options_.vectors > 0)
      write_soft_modes(pre_file_id, prev_step_, prev_alpha_, before, prev_points_, prev_mesh,
                       prev_active, post.full_mapping);

    const auto t0 = Clock::now();
    current = compute(post, true);
    log(step, alpha, "after_avalanche", current,
        std::chrono::duration<double>(Clock::now() - t0).count());
    computed = true;
    lambda_ref_ = -1.0; // new elastic branch: reference = lambda_min right after it
    refining_ = false;
    update_refinement(current);
  } else if (options_.every > 0 && (step % options_.every == 0 || refining_)) {
    const bool grid = (step % options_.every == 0);
    const auto t0 = Clock::now();
    // Vectors are kept so that they can be written if the next step is an avalanche.
    current = compute(post, options_.at_avalanche && options_.vectors > 0);
    log(step, alpha, grid ? "grid" : "refine", current,
        std::chrono::duration<double>(Clock::now() - t0).count());
    computed = true;
    if (avalanche) {
      lambda_ref_ = -1.0;
      refining_ = false;
    }
    update_refinement(current);
  } else if (avalanche) {
    lambda_ref_ = -1.0;
    refining_ = false;
  }

  // Keep POST(step) for the next step
  if (options_.at_avalanche) {
    prev_step_ = step;
    prev_alpha_ = alpha;
    prev_F_ = post.F_external;
    prev_points_ = post.points;
    prev_computed_ = computed;
    if (computed) prev_spectrum_ = std::move(current);
  }
}
