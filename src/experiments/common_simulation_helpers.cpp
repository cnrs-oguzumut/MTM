#include "../../include/geometry/MacroscopicGeometry.h"
#include "../../include/experiments/pbc_stress_control.h"
#include "../../include/experiments/experiment_includes.h"
#include "../../include/optimization/SparseCholesky.h"
#include <chrono>

bool hasConnectivityChanged(const std::vector<ElementTriangle2D> &old_elements,
                            const std::vector<ElementTriangle2D> &new_elements,
                            const std::vector<size_t> &old_active,
                            const std::vector<size_t> &new_active) {
  // Check if number of active elements changed
  if (old_active.size() != new_active.size()) {
    return true;
  }

  // Compare connectivity of each active element
  for (size_t i = 0; i < old_active.size(); i++) {
    const auto &old_elem = old_elements[old_active[i]];
    const auto &new_elem = new_elements[new_active[i]];

    // Compare the three node indices
    for (int j = 0; j < 3; j++) {
      if (old_elem.getNodeIndex(j) != new_elem.getNodeIndex(j)) {
        return true;
      }
    }
  }

  return false;
}

// Order-independent topological comparison: compares the SET of active
// triangles (as sorted node-index triplets). Robust to element re-ordering
// (std::sort by first node is not stable for ties), so it only reports a change
// when the Delaunay connectivity of the current configuration truly differs.
static std::vector<std::array<int, 3>> getTopologyKeySet(
    const std::vector<ElementTriangle2D> &elems,
    const std::vector<size_t> &active) {
  std::vector<std::array<int, 3>> keys;
  keys.reserve(active.size());
  for (size_t idx : active) {
    std::array<int, 3> k = {static_cast<int>(elems[idx].getNodeIndex(0)),
                            static_cast<int>(elems[idx].getNodeIndex(1)),
                            static_cast<int>(elems[idx].getNodeIndex(2))};
    std::sort(k.begin(), k.end());
    keys.push_back(k);
  }
  std::sort(keys.begin(), keys.end());
  return keys;
}

static int countTopologyChanges(const std::vector<ElementTriangle2D> &old_elements,
                                const std::vector<ElementTriangle2D> &new_elements,
                                const std::vector<size_t> &old_active,
                                const std::vector<size_t> &new_active) {
  auto ko = getTopologyKeySet(old_elements, old_active);
  auto kn = getTopologyKeySet(new_elements, new_active);
  std::vector<std::array<int, 3>> diff;
  std::set_symmetric_difference(ko.begin(), ko.end(), kn.begin(), kn.end(),
                                std::back_inserter(diff));
  if (std::getenv("MTM_VERBOSE_LOG") != nullptr) {
    std::cout << "[RemeshTopo] old_active=" << old_active.size()
              << " new_active=" << new_active.size()
              << " differing_triangles=" << diff.size() << std::endl;
  }
  if (diff.empty() && old_active.size() != new_active.size()) {
    return static_cast<int>(std::abs(static_cast<long>(new_active.size()) -
                                     static_cast<long>(old_active.size())));
  }
  return static_cast<int>((diff.size() + 1) / 2);
}

std::tuple<double, Eigen::Matrix2d, int> perform_remeshing_loop_reduction(
    alglib::real_1d_array &x, const UserData *userData,
    const std::vector<int> &contact_atoms,
    const std::vector<int> &boundary_fixed_nodes, const Eigen::Matrix2d &F_ext,
    const Eigen::Matrix<double, 3, 2> &dndx,
    const std::array<double, 2> &offsets,
    const std::vector<int> &original_domain_map,
    const std::vector<std::tuple<double, double>> &translation_map,
    const Point2D &domain_dims_point, int &has_changes,
    int max_iterations = 100, double reference_area = 0.5, bool pbc = false,
    bool optimize_interior = true) {

  bool should_remesh = true;
  int mesh_iteration = 0;
  int accepted_remesh_count = 0;
  double final_energy = 0.0;
  double final_stress = 0.0;
  Eigen::Matrix2d stress_tensor;

  std::vector<Point2D> &square_points = userData->points;
  std::vector<ElementTriangle2D> &elements = userData->elements;
  std::vector<size_t> &active_elements = userData->active_elements;
  const auto &interior_mapping = userData->interior_mapping;
  const auto &full_mapping = userData->full_mapping;
  const int n_vars = interior_mapping.size();
  const int n_points = userData->points.size();
  std::function<double(double)> &potential_func = userData->energy_function;
  std::function<double(double)> &potential_func_der =
      userData->derivative_function;
  double zero = userData->zero_energy;
  double ideal_lattice_parameter = userData->ideal_lattice_parameter;
  double plasticity;

  std::set<std::vector<std::array<int, 3>>> visited_topologies;
  visited_topologies.insert(getTopologyKeySet(elements, active_elements));

  while (should_remesh && mesh_iteration < max_iterations) {
    // === SAVE OLD MESH STATE ===
    std::vector<ElementTriangle2D> old_elements = elements;
    std::vector<size_t> old_active_elements = active_elements;
    std::vector<Point2D> old_points = square_points;
    alglib::real_1d_array old_x;
    old_x.setcontent(x.length(), x.getcontent());
    // PBC stress control: the relaxation also moves F_bar, which a rejected pass must restore
    const Eigen::Matrix2d old_F_bar = (userData->stress_ctrl != nullptr && userData->stress_ctrl->F_bar != nullptr)
        ? *userData->stress_ctrl->F_bar : Eigen::Matrix2d::Identity();
    // ...and, under Cauchy control, the nominal target it was in equilibrium with
    const Eigen::Matrix2d old_P_target =
        (userData->stress_ctrl != nullptr) ? userData->stress_ctrl->P_target : Eigen::Matrix2d::Zero();

    // ============================================
    // ENERGY OF OLD STATE (before any remeshing)
    // ============================================
    UserData oldUserData(square_points, elements, userData->calculator,
                         potential_func, potential_func_der, zero,
                         ideal_lattice_parameter, F_ext, interior_mapping,
                         full_mapping, active_elements, plasticity,
                         userData->second_derivative_function);
    oldUserData.macro_geom_ptr = userData->macro_geom_ptr;
    oldUserData.current_alpha = userData->current_alpha;
    oldUserData.stress_ctrl = userData->stress_ctrl;

    double energy_old = 0.0;
    Eigen::Matrix2d stress_old;
    ConfigurationSaver::calculateEnergyAndStress(&oldUserData, energy_old,
                                                 stress_old, true);

    std::string current_remesh_phase = "REMESH_PASS_" + std::to_string(mesh_iteration + 1);
    AvalancheRecorder::instance().setPhase(current_remesh_phase);
    AvalancheRecorder::instance().recordStep(
        AvalancheEventType::BEFORE_REMESH, -1, energy_old, stress_old(0, 1), 0.0,
        elements.size(), 0, "before topological reconnection");

    if (AvalancheRecorder::instance().isSurgeryVTKEnabled()) {
      int load_step = AvalancheRecorder::instance().currentLoadStep();
      const std::string &dir = AvalancheRecorder::instance().getSurgeryVTKDir();
      std::stringstream ss_before;
      ss_before << dir << "/step_" << std::setw(5) << std::setfill('0') << load_step
                << "_pass_" << std::setw(2) << std::setfill('0') << (mesh_iteration + 1)
                << "_1_before_remesh.vtk";
      ConfigurationSaver::writeToVTKFile(ss_before.str(), old_points, old_elements, &oldUserData, true, {}, 0.0);
    }

    // === PERFORM REMESHING ===
    AdaptiveMesher mesher(domain_dims_point, offsets, original_domain_map,
                          translation_map, full_mapping,
                          1e-6, // Tolerance
                          pbc   // Use periodic copies
    );
    mesher.setUsePeriodicCopies(pbc);

    // Save original positions before remeshing
    alglib::real_1d_array original_x_remesh = mesher.saveOriginalPositions(x);

    std::tie(elements, active_elements) =
        mesher.createMesh(square_points, x, F_ext, &dndx);

    for (auto &element : elements) {
      element.setReferenceArea(reference_area);
      element.set_dof_mapping(full_mapping);
    }

    // === TOPOLOGICAL TRIGGER ===
    // If the Delaunay triangulation of the current configuration is identical
    // to the existing mesh, there is nothing to reconnect: stop immediately
    // (avoids accepting round-off ΔE ~ 1e-13 as a "successful" remesh pass).
    const int changed_tris =
        countTopologyChanges(old_elements, elements, old_active_elements, active_elements);
    if (changed_tris == 0) {
      elements = old_elements;
      active_elements = old_active_elements;
      final_energy = energy_old;
      stress_tensor = stress_old;
      has_changes = (accepted_remesh_count > 0) ? 1 : 0;
      should_remesh = false;
      if (accepted_remesh_count > 1) {
        std::cout << "  ✓ [Remesh Cascade] Stabilized after " << accepted_remesh_count
                  << " passes (" << active_elements.size() << " active triangles)" << std::endl;
      }
      break;
    }

    auto current_topo_key = getTopologyKeySet(elements, active_elements);
    if (visited_topologies.count(current_topo_key)) {
      std::cout << "  ↻ [Remesh Pass " << (mesh_iteration + 1)
                << "] Cycle detected (topology already visited) -> stopping cascade" << std::endl;
      elements = old_elements;
      active_elements = old_active_elements;
      final_energy = energy_old;
      stress_tensor = stress_old;
      has_changes = (accepted_remesh_count > 0) ? 1 : 0;
      should_remesh = false;
      break;
    }
    visited_topologies.insert(current_topo_key);

    UserData newUserData(square_points, elements, userData->calculator,
                         potential_func, potential_func_der, zero,
                         ideal_lattice_parameter, F_ext, interior_mapping,
                         full_mapping, active_elements, plasticity,
                         userData->second_derivative_function);
    newUserData.macro_geom_ptr = userData->macro_geom_ptr;
    newUserData.current_alpha = userData->current_alpha;
    newUserData.stress_ctrl = userData->stress_ctrl;

    // Evaluate immediate post-reconnection energy (coordinates x frozen)
    double energy_reconnected = 0.0;
    Eigen::Matrix2d stress_reconnected;
    ConfigurationSaver::calculateEnergyAndStress(&newUserData, energy_reconnected,
                                                 stress_reconnected, true);
    double delta_E_topo = energy_reconnected - energy_old;

    AvalancheRecorder::instance().recordStep(
        AvalancheEventType::AFTER_REMESH, -1, energy_reconnected, stress_reconnected(0, 1), 0.0,
        elements.size(), (hasConnectivityChanged(old_elements, elements, old_active_elements, active_elements) ? 1 : 0),
        "after topological reconnection (coordinates frozen)");

    if (AvalancheRecorder::instance().isSurgeryVTKEnabled()) {
      int load_step = AvalancheRecorder::instance().currentLoadStep();
      const std::string &dir = AvalancheRecorder::instance().getSurgeryVTKDir();
      std::stringstream ss_after;
      ss_after << dir << "/step_" << std::setw(5) << std::setfill('0') << load_step
               << "_pass_" << std::setw(2) << std::setfill('0') << (mesh_iteration + 1)
               << "_2_after_remesh.vtk";
      ConfigurationSaver::writeToVTKFile(ss_after.str(), square_points, elements, &newUserData, true, {}, 0.0);
    }

    // === OPTIMIZE ON NEW MESH ===
    if (optimize_interior)
      relax_configuration(x, &newUserData, 13);

    map_solver_array_to_points(x, square_points, interior_mapping, n_vars);

    // ============================================
    // ENERGY OF NEW STATE (after remeshing + optimization)
    // ============================================
    UserData newUserData_after(square_points, elements, userData->calculator,
                               potential_func, potential_func_der, zero,
                               ideal_lattice_parameter, F_ext, interior_mapping,
                               full_mapping, active_elements, plasticity,
                               userData->second_derivative_function);
    newUserData_after.macro_geom_ptr = userData->macro_geom_ptr;
    newUserData_after.current_alpha = userData->current_alpha;
    newUserData_after.stress_ctrl = userData->stress_ctrl;

    double energy_new = 0.0;
    Eigen::Matrix2d stress_new;
    ConfigurationSaver::calculateEnergyAndStress(&newUserData_after, energy_new,
                                                 stress_new, true);

    // ============================================
    // COMPARE OLD vs NEW
    // ============================================
    double energy_change = energy_new - energy_old;
    // Cauchy control: H_old and H_new use different nominal targets; add the work of the target change
    if (userData->stress_ctrl != nullptr && userData->stress_ctrl->F_bar != nullptr) {
      energy_change += userData->stress_ctrl->cauchyWorkCorrection(old_F_bar, old_P_target,
                                                                   *userData->stress_ctrl->F_bar,
                                                                   userData->stress_ctrl->P_target);
    }

    if (AvalancheRecorder::instance().isSurgeryVTKEnabled()) {
      int load_step = AvalancheRecorder::instance().currentLoadStep();
      const std::string &dir = AvalancheRecorder::instance().getSurgeryVTKDir();
      std::stringstream ss_rel;
      ss_rel << dir << "/step_" << std::setw(5) << std::setfill('0') << load_step
             << "_pass_" << std::setw(2) << std::setfill('0') << (mesh_iteration + 1)
             << (energy_change >= 0 ? "_3_relaxed_rejected.vtk" : "_3_relaxed_accepted.vtk");
      ConfigurationSaver::writeToVTKFile(ss_rel.str(), square_points, elements, &newUserData_after, true, {}, 0.0);
    }

    // Energy acceptance & cascade continuation:
    // 1) Genuine energy increase (> remesh_energy_tol) -> REJECT and restore previous state.
    // 2) Energy neutral reconnection (within [-cascade_threshold, +remesh_energy_tol]) -> ACCEPT, but
    //    TERMINATE the cascade loop (should_remesh = false) to prevent infinite chatter on
    //    degenerate/co-circular lattice diagonals (e.g. square lattice).
    // 3) Significant energy drop (< -cascade_threshold) -> ACCEPT and CONTINUE cascade loop
    //    to track multi-slip dislocation avalanche relaxation.
    const double cascade_energy_threshold = 1e-6 * std::max(1.0, std::abs(energy_old));
    const double remesh_energy_tol = 1e-8 * std::max(1.0, std::abs(energy_old));
    std::cout << "  ↻ [Remesh Pass " << (mesh_iteration + 1) << "] Delaunay reconnected "
              << changed_tris << " tris | E: " << std::scientific << std::setprecision(5)
              << energy_old << " -> " << energy_new
              << " (ΔE_topo=" << std::showpos << std::setprecision(2) << delta_E_topo
              << ", ΔE_rel=" << energy_change << std::noshowpos << ") -> ";
    if (energy_change > remesh_energy_tol) {
      std::cout << "REJECTED (restored previous mesh)" << std::endl;

      AvalancheRecorder::instance().recordStep(
          AvalancheEventType::REMESH_REJECTED, -1, energy_old, stress_old(0, 1), 0.0,
          old_elements.size(), 0, "remesh rejected: energy increased");

      // RESTORE OLD STATE
      x = old_x;
      if (userData->stress_ctrl != nullptr && userData->stress_ctrl->F_bar != nullptr) {
        *userData->stress_ctrl->F_bar = old_F_bar;
        userData->stress_ctrl->P_target = old_P_target;
      }
      square_points = old_points;
      elements = old_elements;
      active_elements = old_active_elements;

      // Final values from restored state
      final_energy = energy_old;
      stress_tensor = stress_old;

      // Stop remeshing - we couldn't improve further
      should_remesh = false;
      has_changes = (accepted_remesh_count > 0) ? 1 : 0;

    } else if (energy_change >= -cascade_energy_threshold) {
      std::cout << "ACCEPTED (neutral)" << std::endl;

      AvalancheRecorder::instance().recordStep(
          AvalancheEventType::REMESH_ACCEPTED, -1, energy_new, stress_new(0, 1), 0.0,
          elements.size(), 0, "remesh accepted: energy neutral");

      final_energy = energy_new;
      stress_tensor = stress_new;
      has_changes = 1;
      accepted_remesh_count++;
      should_remesh = false;

    } else {
      std::cout << "ACCEPTED (plastic relaxation)" << std::endl;

      AvalancheRecorder::instance().recordStep(
          AvalancheEventType::REMESH_ACCEPTED, -1, energy_new, stress_new(0, 1), 0.0,
          elements.size(), 0, "remesh accepted: energy decreased");

      // Keep new state (already in place)
      final_energy = energy_new;
      stress_tensor = stress_new;
      has_changes = 1;
      accepted_remesh_count++;

      // Continue cascade loop to see if further plastic relaxation occurs
      should_remesh = true;

      // Force control past the limit load: the dead load has no equilibrium, stop the cascade
      if ((userData->macro_geom_ptr != nullptr &&
           userData->macro_geom_ptr->exceedsRunawayLimit(square_points)) ||
          (userData->stress_ctrl != nullptr && userData->stress_ctrl->F_bar != nullptr &&
           userData->stress_ctrl->exceedsLimit(*userData->stress_ctrl->F_bar))) {
        std::cout << "  ↻ [Remesh Cascade] Stopped: nodal displacement since the step start exceeds the "
                     "limit-load guard (dead-load runaway)" << std::endl;
        should_remesh = false;
      }
    }

    mesh_iteration++;
  }

  return {final_energy, stress_tensor, accepted_remesh_count > 0 ? accepted_remesh_count : -mesh_iteration};
}

// Final calculation if loop completes without breaking
// if (mesh_iteration >= max_iterations) {
//   UserData finalUserData(square_points, elements, userData->calculator,
//                        potential_func, potential_func_der, zero,
//                        ideal_lattice_parameter, F_ext, interior_mapping,
//                        full_mapping, active_elements, plasticity);
//   ConfigurationSaver::calculateEnergyAndStress(&finalUserData, final_energy,
//                                                 stress_tensor, true);
//   final_stress = stress_tensor(0, 1);
// }

namespace {
std::unique_ptr<PreconditionedLBFGS> g_relaxation_solver;
int g_relaxation_first_step = 0;
int g_relaxation_step = 0;
int g_relaxation_maxits = 0; // 0 = unlimited (default); set only for force-controlled runs
}

void set_relaxation_max_iterations(int maxits) {
  g_relaxation_maxits = std::max(0, maxits);
  if (g_relaxation_solver) g_relaxation_solver->setMaxIterations(g_relaxation_maxits);
}

void configure_relaxation_solver(const PreconditionedLBFGSOptions &options,
                                 int first_step) {
  g_relaxation_first_step = first_step;
  if (options.type == LBFGSPreconditioner::None) {
    g_relaxation_solver.reset();
    return;
  }
  PreconditionedLBFGSOptions opts = options;
  if (!opts.progress_callback) {
    opts.progress_callback = [](int nfev, double energy, double grad_max) {
      if (AvalancheRecorder::instance().isEnabled()) {
        AvalancheRecorder::instance().recordStep(
            AvalancheEventType::LBFGS_ITER, nfev, energy, 0.0, grad_max,
            0, 0, "L-BFGS evaluation");
      }
    };
  }
  if (g_relaxation_maxits > 0) opts.maxits = g_relaxation_maxits;
  g_relaxation_solver = std::make_unique<PreconditionedLBFGS>(opts);
}

void relaxation_begin_step(int step) { g_relaxation_step = step; }

namespace {
// Rigid force grips: the solver vector is z = [x (2 n_vars), U components]. Gripped atoms sit at X + U;
// Pi = E_el(x, U) - alpha * F0 . U, so dPi/dU = sum of the element forces on the grip atoms - alpha * F0.
struct GripRelaxContext {
  UserData *ud = nullptr;
  const MacroscopicGeometryConfig *mg = nullptr;
  int n2 = 0;
  alglib::real_1d_array x_nodes, g_nodes;
  std::vector<int> probe_nodes;
  std::vector<Eigen::Vector2d> probe_out;
  std::vector<std::pair<size_t, size_t>> probe_range; // per grip [begin, end) in probe_nodes
  std::vector<std::array<int, 2>> var_index;           // per grip index of U_x / U_y in z, or -1 (held)
  std::vector<Eigen::Vector2d> U_held;
};

Eigen::Vector2d gripUFromVector(const GripRelaxContext &ctx, size_t g, const alglib::real_1d_array &z) {
  Eigen::Vector2d U = ctx.U_held[g];
  if (ctx.var_index[g][0] >= 0) U.x() = z[ctx.var_index[g][0]];
  if (ctx.var_index[g][1] >= 0) U.y() = z[ctx.var_index[g][1]];
  return U;
}

void grip_enthalpy(const alglib::real_1d_array &z, double &func, alglib::real_1d_array &grad, void *ptr) {
  auto *ctx = static_cast<GripRelaxContext *>(ptr);
  const size_t n_grips = ctx->mg->cached_force_grips.size();
  for (size_t g = 0; g < n_grips; ++g) {
    ctx->mg->setGripPositions(static_cast<int>(g), gripUFromVector(*ctx, g, z), ctx->ud->points);
  }
  for (int q = 0; q < ctx->n2; ++q) ctx->x_nodes[q] = z[q];
  ctx->ud->force_probe_nodes = &ctx->probe_nodes;
  ctx->ud->force_probe_out = &ctx->probe_out;
  minimize_energy_with_triangles(ctx->x_nodes, func, ctx->g_nodes, ctx->ud);
  ctx->ud->force_probe_nodes = nullptr;
  ctx->ud->force_probe_out = nullptr;
  for (int q = 0; q < ctx->n2; ++q) grad[q] = ctx->g_nodes[q];
  const double alpha = ctx->ud->current_alpha;
  for (size_t g = 0; g < n_grips; ++g) {
    const Eigen::Vector2d F = alpha * ctx->mg->cached_force_grips[g].F0;
    const Eigen::Vector2d U = gripUFromVector(*ctx, g, z);
    func -= F.dot(U);
    Eigen::Vector2d f_sum = Eigen::Vector2d::Zero();
    for (size_t k = ctx->probe_range[g].first; k < ctx->probe_range[g].second; ++k) f_sum += ctx->probe_out[k];
    if (ctx->var_index[g][0] >= 0) grad[ctx->var_index[g][0]] = f_sum.x() - F.x();
    if (ctx->var_index[g][1] >= 0) grad[ctx->var_index[g][1]] = f_sum.y() - F.y();
  }
}

void relax_with_force_grips(alglib::real_1d_array &x, UserData *userData, int corrections, int maxits) {
  GripRelaxContext ctx;
  ctx.ud = userData;
  ctx.mg = userData->macro_geom_ptr;
  ctx.n2 = static_cast<int>(x.length());
  const auto &grips = ctx.mg->cached_force_grips;
  int n_extra = 0;
  for (const auto &grip : grips) {
    const size_t begin = ctx.probe_nodes.size();
    ctx.probe_nodes.insert(ctx.probe_nodes.end(), grip.nodes.begin(), grip.nodes.end());
    ctx.probe_range.emplace_back(begin, ctx.probe_nodes.size());
    std::array<int, 2> idx = {-1, -1};
    if (grip.free_x) idx[0] = ctx.n2 + n_extra++;
    if (grip.free_y) idx[1] = ctx.n2 + n_extra++;
    ctx.var_index.push_back(idx);
  }
  alglib::real_1d_array z;
  z.setlength(ctx.n2 + n_extra);
  for (int q = 0; q < ctx.n2; ++q) z[q] = x[q];
  for (size_t g = 0; g < grips.size(); ++g) {
    const Eigen::Vector2d U0 = ctx.mg->gripDisplacement(static_cast<int>(g), userData->points);
    ctx.U_held.push_back(U0);
    if (ctx.var_index[g][0] >= 0) z[ctx.var_index[g][0]] = U0.x();
    if (ctx.var_index[g][1] >= 0) z[ctx.var_index[g][1]] = U0.y();
  }
  ctx.x_nodes.setlength(ctx.n2);
  ctx.g_nodes.setlength(ctx.n2);

  // Optional finite-difference check of dPi/dU (MTM_GRIP_GRADCHECK=1)
  if (std::getenv("MTM_GRIP_GRADCHECK") != nullptr && n_extra > 0) {
    alglib::real_1d_array g0, gt, zt;
    g0.setlength(z.length()); gt.setlength(z.length()); zt.setlength(z.length());
    double f0 = 0.0, fp = 0.0, fm = 0.0;
    grip_enthalpy(z, f0, g0, &ctx);
    const double h = 1e-6;
    for (int q = ctx.n2; q < z.length(); ++q) {
      for (int r = 0; r < z.length(); ++r) zt[r] = z[r];
      zt[q] = z[q] + h; grip_enthalpy(zt, fp, gt, &ctx);
      zt[q] = z[q] - h; grip_enthalpy(zt, fm, gt, &ctx);
      const double fd = (fp - fm) / (2.0 * h);
      std::cout << "[Grip GradCheck] alpha = " << userData->current_alpha << " | dPi/dU[" << (q - ctx.n2)
                << "] = " << std::scientific << g0[q] << " fd = " << fd << " rel = "
                << std::abs(fd - g0[q]) / std::max(1e-12, std::abs(fd) + std::abs(g0[q]))
                << std::defaultfloat << std::endl;
    }
  }

  // The stiffness preconditioner assembles node DOFs only, so grips use the plain L-BFGS
  LBFGSOptimizer optimizer(corrections, 0.0, 0.0, 0.0, maxits);
  optimizer.optimize(z, grip_enthalpy, &ctx);

  for (int q = 0; q < ctx.n2; ++q) x[q] = z[q];
  for (size_t g = 0; g < grips.size(); ++g) {
    ctx.mg->setGripPositions(static_cast<int>(g), gripUFromVector(ctx, g, z), userData->points);
  }
}
// PBC nominal-stress control: z = [w (2 n_vars), free components of F_bar], x = F_bar w.
// dE/dw = F_bar^T dE/dx;  dE/dF_bar = (sum_e A_e P_e F_e^T) F_bar^{-T};  H = E - A0 Pt : (F_bar - I).
struct StressRelaxContext {
  UserData *ud = nullptr;
  PBCStressControl *sc = nullptr;
  int n = 0; // n_vars
  Eigen::Matrix2d F_base = Eigen::Matrix2d::Identity();
  std::vector<int> comp_index; // per component k (11, 12, 22): index in z, or -1
  double f_scale = 1.0;         // z_F = f_scale * F (f_scale = sqrt(A0)): every unknown is a length in a0
  alglib::real_1d_array x_nodes, g_nodes;
};

Eigen::Matrix2d stressFbarFromVector(const StressRelaxContext &ctx, const alglib::real_1d_array &z) {
  Eigen::Matrix2d F = ctx.F_base;
  for (int k = 0; k < 3; ++k) {
    if (ctx.comp_index[k] >= 0)
      F(PBCStressControl::compRow(k), PBCStressControl::compCol(k)) = z[ctx.comp_index[k]] / ctx.f_scale;
  }
  F(1, 0) = 0.0; // gauge: upper triangular
  return F;
}

void applyFbar(StressRelaxContext &ctx, const Eigen::Matrix2d &F) {
  *ctx.sc->F_bar = F;
  for (auto &el : ctx.ud->elements) el.setExternalDeformation(F);
}

void stress_enthalpy(const alglib::real_1d_array &z, double &func, alglib::real_1d_array &grad, void *ptr) {
  auto *ctx = static_cast<StressRelaxContext *>(ptr);
  const int n = ctx->n;
  const Eigen::Matrix2d F = stressFbarFromVector(*ctx, z);
  applyFbar(*ctx, F);
  for (int i = 0; i < n; ++i) {
    const Eigen::Vector2d xi = F * Eigen::Vector2d(z[i], z[n + i]);
    ctx->x_nodes[i] = xi.x();
    ctx->x_nodes[n + i] = xi.y();
  }
  Eigen::Matrix2d virial = Eigen::Matrix2d::Zero();
  ctx->ud->virial_out = &virial;
  minimize_energy_with_triangles(ctx->x_nodes, func, ctx->g_nodes, ctx->ud);
  ctx->ud->virial_out = nullptr;
  for (int i = 0; i < n; ++i) {
    const Eigen::Vector2d gw = F.transpose() * Eigen::Vector2d(ctx->g_nodes[i], ctx->g_nodes[n + i]);
    grad[i] = gw.x();
    grad[n + i] = gw.y();
  }
  const Eigen::Matrix2d dEdF = virial * F.inverse().transpose();
  for (int k = 0; k < 3; ++k) {
    if (ctx->comp_index[k] < 0) continue;
    const int r = PBCStressControl::compRow(k), c = PBCStressControl::compCol(k);
    grad[ctx->comp_index[k]] = (dEdF(r, c) - ctx->sc->A0 * ctx->sc->P_target(r, c)) / ctx->f_scale;
  }
  func += ctx->sc->workTerm(F);
  if (std::getenv("MTM_STRESS_TRACE") != nullptr) {
    static long n_eval = 0;
    if (++n_eval % 200 == 0) {
      double gmax = 0.0;
      for (int q = 0; q < grad.length(); ++q) gmax = std::max(gmax, std::abs(grad[q]));
      std::cout << "[Stress Trace] eval " << n_eval << " H = " << std::scientific << func << " F = (" << F(0, 0)
                << ", " << F(0, 1) << ", " << F(1, 1) << ") |g|max = " << gmax << std::defaultfloat << std::endl;
    }
  }
}

// Stiffness preconditioner of the stress-controlled solve, bordered by the free cell components.
// In z = [w, f_scale F_bar_k] (x = F_bar w, D = F_bar (x) I_n in the solver layout):
//   M_ww = D^T K D,   M_wz = D^T K_xF / f_scale,   M_zz = K_FF / f_scale^2,
// with K, K_xF, K_FF assembled from the PSD-projected element tangents (Gauss-Newton). Factored as
//   P = R R^T,  R = [[R_w, 0], [Bt^T, L_S]],  R_w R_w^T = M_ww + shift,  Bt = R_w^{-1} M_wz,
//   L_S L_S^T = M_zz - Bt^T Bt (Schur complement, + shift).
// Border of the stiffness by the free cell components (Gauss-Newton, PSD-projected element tangents,
// the same as FastStiffnessAssembler's): with F_e = F_bar G_e, dF_e/dF_bar_(r,c) = e_r (x) G_e^T e_c, and
//   K_xF[d][k] = sum_e area_e (dF_e/dx_d) : A_e : (dF_e/dF_bar_k),  K_FF[k][l] = sum_e area_e (dF_e/dF_bar_k) : A_e : (dF_e/dF_bar_l)
// in the x layout [u_0..u_{n-1}, v_0..v_{n-1}]. Uses the F_e cached by the stiffness assembly at the same x.
void assemble_box_border(UserData *ud, const Eigen::Matrix2d &F_bar, const std::vector<std::pair<int, int>> &comps,
                         bool project_psd, int n, Eigen::MatrixXd &K_xF, Eigen::MatrixXd &K_FF) {
  const int m = static_cast<int>(comps.size());
  const auto &active = ud->active_elements;
  const double normalisation = ud->calculator.getUnitCellArea();
  std::function<double(double)> fd_d2pot;
  const std::function<double(double)> *d2pot = ud->second_derivative_function;
  if (d2pot == nullptr) {
    if (dynamic_cast<Strain_Energy_LatticeCalculator *>(&ud->calculator) == nullptr) {
      const auto &df = ud->derivative_function;
      fd_d2pot = [&df](double r) {
        const double eps = 1e-6 * std::max(1.0, std::abs(r));
        return (df(r + eps) - df(r - eps)) / (2.0 * eps);
      };
      d2pot = &fd_d2pot;
    } else {
      d2pot = &ud->derivative_function; // ignored by the analytic strain energy
    }
  }
  const Eigen::Matrix2d F_bar_inv = F_bar.inverse();
  K_xF.setZero(2 * n, m);
  K_FF.setZero(m, m);
  std::vector<Eigen::MatrixXd> K_FF_thread(omp_get_max_threads(), Eigen::MatrixXd::Zero(m, m));
  double *xF = K_xF.data();
#pragma omp parallel for schedule(guided, 256)
  for (size_t e = 0; e < active.size(); e++) {
    const ElementTriangle2D &element = ud->elements[active[e]];
    const Eigen::Matrix2d &F = element.getDeformationGradient();
    const auto red = reduce_metric(element.getMetricTensor());
    const Eigen::Matrix2d dE_dC =
        ud->calculator.calculate_derivative(red.C_reduced, ud->derivative_function) / normalisation;
    const HessianComponents d2E =
        ud->calculator.calculate_dseconderivative_components(red.C_reduced, ud->derivative_function, *d2pot);
    Eigen::Matrix4d A = element_lagrangian_tangent(F, red.m_matrix, dE_dC, d2E, normalisation);
    if (project_psd) {
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> eig(A);
      const Eigen::Vector4d lambda = eig.eigenvalues().cwiseMax(0.0);
      A = eig.eigenvectors() * lambda.asDiagonal() * eig.eigenvectors().transpose();
    }
    const Eigen::Matrix2d G = F_bar_inv * F;
    Eigen::Matrix<double, 4, Eigen::Dynamic> gF(4, m);
    gF.setZero();
    for (int k = 0; k < m; k++) {
      const int r = comps[k].first, c = comps[k].second;
      for (int L = 0; L < 2; L++) gF(2 * r + L, k) = G(c, L); // row 2i+K = F_iK
    }
    const Eigen::Matrix<double, 3, 2> &dN = element.getDNdX();
    Eigen::Matrix<double, 4, 6> Gm = Eigen::Matrix<double, 4, 6>::Zero();
    for (int a = 0; a < 3; a++)
      for (int i = 0; i < 2; i++)
        for (int K = 0; K < 2; K++) Gm(2 * i + K, 2 * a + i) = dN(a, K);
    const double area = element.getReferenceArea();
    const Eigen::Matrix<double, 4, Eigen::Dynamic> AgF = A * gF;
    const Eigen::MatrixXd Kxk = area * Gm.transpose() * AgF; // 6 x m
    K_FF_thread[omp_get_thread_num()] += area * gF.transpose() * AgF;
    for (int a = 0; a < 3; a++) {
      const int s = ud->full_mapping[element.getNodeIndex(a)].second;
      if (s < 0) continue;
      for (int i = 0; i < 2; i++) {
        const int d = s + i * n;
        for (int k = 0; k < m; k++) {
#pragma omp atomic
          xF[static_cast<size_t>(k) * 2 * n + d] += Kxk(2 * a + i, k);
        }
      }
    }
  }
  for (const auto &Kt : K_FF_thread) K_FF += Kt;
}

std::uint64_t sparse_pattern_hash(const Eigen::SparseMatrix<double> &M) {
  std::uint64_t h = 1469598103934665603ULL;
  auto mix = [&h](std::uint64_t v) { h = (h ^ v) * 1099511628211ULL; };
  mix(static_cast<std::uint64_t>(M.rows()));
  for (int c = 0; c <= M.outerSize(); ++c) mix(static_cast<std::uint64_t>(M.outerIndexPtr()[c]));
  for (long q = 0; q < M.nonZeros(); ++q) mix(static_cast<std::uint64_t>(M.innerIndexPtr()[q]));
  return h;
}

// Connectivity of the active mesh in solver dofs (a remesh changes it at fixed node count)
std::uint64_t stress_mesh_signature(const UserData *ud) {
  std::uint64_t h = 1469598103934665603ULL ^ ud->active_elements.size();
  for (size_t idx : ud->active_elements)
    for (int a = 0; a < 3; a++)
      h = (h ^ static_cast<std::uint64_t>(ud->full_mapping[ud->elements[idx].getNodeIndex(a)].second + 2)) *
          1099511628211ULL;
  return h;
}

struct StressPreconditioner {
  FastStiffnessAssembler assembler;
  std::unique_ptr<SparseCholeskyBackend> chol;
  std::uint64_t analyzed_signature = 0;
  int n2 = 0, m = 0;
  Eigen::MatrixXd Bt;  // n2 x m
  Eigen::MatrixXd L_S; // m x m, lower
  std::vector<double> tmp;
  int calls_since_refresh = 0;
  int last_iterations = 0;                    // iterations of the last minimization (refresh_iterations)
  std::uint64_t mesh_sig = 0;                 // mesh the factor was built on
  std::vector<std::pair<int, int>> comps_built; // free cell components it was built for
  double assemble_seconds = 0.0, factorize_seconds = 0.0;

  bool build(const alglib::real_1d_array &x_nodes, UserData *ud, const Eigen::Matrix2d &F,
             const std::vector<std::pair<int, int>> &comps, double f_scale,
             const PreconditionedLBFGSOptions &opt) {
    using SpMat = Eigen::SparseMatrix<double>;
    const auto t0 = std::chrono::high_resolution_clock::now();
    const int n = static_cast<int>(ud->interior_mapping.size());
    n2 = 2 * n;
    m = static_cast<int>(comps.size());
    mesh_sig = stress_mesh_signature(ud);
    comps_built = comps;
    const SpMat &K = assembler.assembleStiffness(x_nodes, ud, opt.project_psd); // also caches F_e at x
    Eigen::MatrixXd K_xF, K_FF;
    assemble_box_border(ud, F, comps, opt.project_psd, n, K_xF, K_FF);
    std::vector<Eigen::Triplet<double>> trip;
    trip.reserve(4 * static_cast<size_t>(n));
    for (int r = 0; r < 2; ++r)
      for (int c = 0; c < 2; ++c)
        for (int i = 0; i < n; ++i) trip.emplace_back(r * n + i, c * n + i, F(r, c)); // explicit zeros kept
    SpMat D(n2, n2);
    D.setFromTriplets(trip.begin(), trip.end());
    const SpMat Dt = D.transpose();
    SpMat Mww = Dt * K * D;
    Mww.makeCompressed();
    const Eigen::MatrixXd Mwz = (Dt * K_xF) / f_scale;
    const Eigen::MatrixXd Mzz = K_FF / (f_scale * f_scale);
    const auto t1 = std::chrono::high_resolution_clock::now();
    assemble_seconds = std::chrono::duration<double>(t1 - t0).count();

    const std::uint64_t pattern = sparse_pattern_hash(Mww);
    if (!chol || analyzed_signature != pattern) {
      chol = make_sparse_cholesky();
      chol->analyze(Mww);
      analyzed_signature = pattern;
    }
    double shift = opt.shift_rel * Mww.diagonal().mean();
    bool ok = false;
    for (int attempt = 0; attempt < 6 && !ok; ++attempt, shift *= 100.0) ok = chol->factorize(Mww, shift);
    if (!ok) {
      chol.reset();
      return false;
    }
    Bt.resize(n2, m);
    for (int k = 0; k < m; ++k) chol->apply_inverse(Mwz.col(k).data(), Bt.col(k).data());
    Eigen::MatrixXd S = Mzz - Bt.transpose() * Bt;
    S = 0.5 * (S + S.transpose());
    double s_shift = opt.shift_rel * std::max(1e-300, Mzz.diagonal().cwiseAbs().mean());
    ok = false;
    for (int attempt = 0; attempt < 8 && !ok; ++attempt, s_shift *= 100.0) {
      Eigen::LLT<Eigen::MatrixXd> llt(S + s_shift * Eigen::MatrixXd::Identity(m, m));
      if (llt.info() == Eigen::Success) {
        L_S = llt.matrixL();
        ok = true;
      }
    }
    if (!ok) {
      chol.reset();
      return false;
    }
    tmp.assign(n2, 0.0);
    factorize_seconds = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t1).count();
    return true;
  }

  // gz = R^{-1} g
  void applyInverse(const double *g, double *gz) {
    chol->apply_inverse(g, gz);
    Eigen::Map<const Eigen::VectorXd> gw(gz, n2);
    const Eigen::VectorXd rhs = Eigen::Map<const Eigen::VectorXd>(g + n2, m) - Bt.transpose() * gw;
    Eigen::Map<Eigen::VectorXd>(gz + n2, m) = L_S.triangularView<Eigen::Lower>().solve(rhs);
  }
  // dz = R^{-T} t
  void applyInverseTranspose(const double *t, double *dz) {
    const Eigen::VectorXd dF =
        L_S.transpose().triangularView<Eigen::Upper>().solve(Eigen::Map<const Eigen::VectorXd>(t + n2, m));
    Eigen::Map<Eigen::VectorXd> tw(tmp.data(), n2);
    tw = Eigen::Map<const Eigen::VectorXd>(t, n2) - Bt * dF;
    chol->apply_inverse_transpose(tmp.data(), dz);
    Eigen::Map<Eigen::VectorXd>(dz + n2, m) = dF;
  }
};

StressPreconditioner g_stress_precond;

struct PrecondStressContext {
  StressRelaxContext *inner = nullptr;
  StressPreconditioner *P = nullptr;
  std::vector<double> z0;
  alglib::real_1d_array z_work, g_work;
  double grad_tol = 0.0;
  bool hit = false;
  std::vector<double> z_hit;
  alglib::minlbfgsstate *state = nullptr;
};

void precond_stress_enthalpy(const alglib::real_1d_array &t, double &func, alglib::real_1d_array &grad, void *ptr) {
  auto *pc = static_cast<PrecondStressContext *>(ptr);
  const int N = static_cast<int>(t.length());
  double *z = pc->z_work.getcontent();
  pc->P->applyInverseTranspose(t.getcontent(), z);
  for (int q = 0; q < N; ++q) z[q] += pc->z0[q];
  stress_enthalpy(pc->z_work, func, pc->g_work, pc->inner);
  pc->P->applyInverse(pc->g_work.getcontent(), grad.getcontent());
  if (pc->grad_tol > 0.0 && !pc->hit) {
    double gmax = 0.0;
    for (int q = 0; q < N; ++q) gmax = std::max(gmax, std::abs(pc->g_work[q]));
    if (gmax < pc->grad_tol) {
      pc->hit = true;
      pc->z_hit.assign(z, z + N);
      alglib::minlbfgsrequesttermination(*pc->state);
    }
  }
}

// Stiffness-preconditioned L-BFGS on z (same change of variables as PreconditionedLBFGS); the
// preconditioner is rebuilt at the current point when the run stalls before grad_tol (avalanche).
// Returns false (z untouched) when the preconditioner cannot be factorized.
bool minimize_stress_enthalpy_preconditioned(alglib::real_1d_array &z, StressRelaxContext &ctx, int corrections,
                                             int maxits, const PreconditionedLBFGSOptions &opt,
                                             alglib::minlbfgsreport &out_rep) {
  const int N = static_cast<int>(z.length());
  const int n = ctx.n;
  std::vector<std::pair<int, int>> comps;
  for (int k = 0; k < 3; ++k)
    if (ctx.comp_index[k] >= 0) comps.emplace_back(PBCStressControl::compRow(k), PBCStressControl::compCol(k));

  PrecondStressContext pc;
  pc.inner = &ctx;
  pc.P = &g_stress_precond;
  pc.z_work.setlength(N);
  pc.g_work.setlength(N);
  pc.grad_tol = opt.grad_tol;
  alglib::real_1d_array t;
  t.setlength(N);
  alglib::minlbfgsstate state;
  alglib::minlbfgsreport rep;
  pc.state = &state;
  const auto start = std::chrono::high_resolution_clock::now();
  double setup = 0.0;
  int iterations = 0, nfev = 0, restarts = 0, termination = 0;

  auto build_at = [&](const alglib::real_1d_array &zz, bool force) {
    // Reuse (--precond-refresh > 1) only on the same mesh and cell components, and, with
    // --precond-refresh-iters, while the last minimization was short
    if (!force && g_stress_precond.chol && g_stress_precond.n2 == 2 * n &&
        g_stress_precond.comps_built == comps && g_stress_precond.mesh_sig == stress_mesh_signature(ctx.ud) &&
        g_stress_precond.calls_since_refresh < std::max(1, opt.refresh_every) &&
        (opt.refresh_iterations <= 0 || g_stress_precond.last_iterations <= opt.refresh_iterations))
      return true;
    const Eigen::Matrix2d F = stressFbarFromVector(ctx, zz);
    applyFbar(ctx, F);
    for (int i = 0; i < n; ++i) {
      const Eigen::Vector2d xi = F * Eigen::Vector2d(zz[i], zz[n + i]);
      ctx.x_nodes[i] = xi.x();
      ctx.x_nodes[n + i] = xi.y();
    }
    const bool ok = g_stress_precond.build(ctx.x_nodes, ctx.ud, F, comps, ctx.f_scale, opt);
    g_stress_precond.calls_since_refresh = 0;
    setup += g_stress_precond.assemble_seconds + g_stress_precond.factorize_seconds;
    return ok;
  };

  if (!build_at(z, false)) return false;
  g_stress_precond.calls_since_refresh++;
  for (int run = 0;; ++run) {
    pc.z0.assign(z.getcontent(), z.getcontent() + N);
    pc.hit = false;
    for (int q = 0; q < N; ++q) t[q] = 0.0;
    alglib::minlbfgscreate(corrections, t, state);
    alglib::minlbfgssetcond(state, 0.0, 0.0, 0.0, maxits > 0 ? maxits : (opt.grad_tol > 0.0 ? 200000 : 0));
    alglib::minlbfgsoptimize(state, precond_stress_enthalpy, nullptr, &pc);
    alglib::minlbfgsresults(state, t, rep);
    if (pc.hit) {
      for (int q = 0; q < N; ++q) z[q] = pc.z_hit[q];
    } else {
      g_stress_precond.applyInverseTranspose(t.getcontent(), z.getcontent());
      for (int q = 0; q < N; ++q) z[q] += pc.z0[q];
    }
    iterations += static_cast<int>(rep.iterationscount);
    nfev += static_cast<int>(rep.nfev);
    termination = static_cast<int>(rep.terminationtype);
    const bool stalled = opt.grad_tol > 0.0 && !pc.hit && rep.terminationtype != 5;
    if (!stalled || run >= opt.max_restarts) break;
    if (!build_at(z, true)) break;
    g_stress_precond.calls_since_refresh = 1;
    ++restarts;
  }
  g_stress_precond.last_iterations = iterations;
  out_rep.iterationscount = iterations;
  out_rep.nfev = nfev;
  out_rep.terminationtype = termination;
  if (opt.verbose) {
    double gmax = 0.0, h = 0.0;
    stress_enthalpy(z, h, pc.g_work, &ctx);
    for (int q = 0; q < N; ++q) gmax = std::max(gmax, std::abs(pc.g_work[q]));
    const double total = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
    std::cout << "Optimization completed [stiffness+cell]: iterations=" << iterations << ", nfev=" << nfev
              << ", termination type=" << termination << ", max|g|=" << gmax << ", setup=" << setup
              << "s, solve=" << (total - setup) << "s" << (restarts ? ", restarts=" + std::to_string(restarts) : "")
              << std::endl;
  }
  return true;
}

void relax_with_stress_control(alglib::real_1d_array &x, UserData *userData, int corrections, int maxits) {
  StressRelaxContext ctx;
  ctx.ud = userData;
  ctx.sc = userData->stress_ctrl;
  ctx.n = static_cast<int>(x.length() / 2);
  ctx.F_base = *ctx.sc->F_bar;
  ctx.F_base(1, 0) = 0.0;
  ctx.f_scale = std::sqrt(std::max(1e-12, ctx.sc->A0));
  const Eigen::Matrix2d F_inv = ctx.F_base.inverse();
  int n_extra = 0;
  ctx.comp_index.assign(3, -1);
  for (int k = 0; k < 3; ++k) {
    if (ctx.sc->comp_free[k]) ctx.comp_index[k] = 2 * ctx.n + n_extra++;
  }
  alglib::real_1d_array z;
  z.setlength(2 * ctx.n + n_extra);
  for (int i = 0; i < ctx.n; ++i) {
    const Eigen::Vector2d w = F_inv * Eigen::Vector2d(x[i], x[ctx.n + i]);
    z[i] = w.x();
    z[ctx.n + i] = w.y();
  }
  for (int k = 0; k < 3; ++k) {
    if (ctx.comp_index[k] >= 0)
      z[ctx.comp_index[k]] = ctx.f_scale * ctx.F_base(PBCStressControl::compRow(k), PBCStressControl::compCol(k));
  }
  ctx.x_nodes.setlength(2 * ctx.n);
  ctx.g_nodes.setlength(2 * ctx.n);

  // Optional finite-difference check of dH/dF_bar and a few dH/dw (MTM_STRESS_GRADCHECK=1)
  if (std::getenv("MTM_STRESS_GRADCHECK") != nullptr) {
    alglib::real_1d_array g0, gt, zt;
    g0.setlength(z.length()); gt.setlength(z.length()); zt.setlength(z.length());
    double f0 = 0.0, fp = 0.0, fm = 0.0;
    stress_enthalpy(z, f0, g0, &ctx);
    std::vector<int> probe;
    for (int k = 0; k < 3; ++k) if (ctx.comp_index[k] >= 0) probe.push_back(ctx.comp_index[k]);
    for (int i = 0; i < ctx.n && probe.size() < static_cast<size_t>(n_extra) + 6; i += std::max(1, ctx.n / 5)) {
      probe.push_back(i);
      probe.push_back(ctx.n + i);
    }
    const double h = 1e-6;
    for (int q : probe) {
      for (int r = 0; r < z.length(); ++r) zt[r] = z[r];
      zt[q] = z[q] + h; stress_enthalpy(zt, fp, gt, &ctx);
      zt[q] = z[q] - h; stress_enthalpy(zt, fm, gt, &ctx);
      const double fd = (fp - fm) / (2.0 * h);
      std::cout << "[Stress GradCheck] alpha = " << userData->current_alpha << " | "
                << (q >= 2 * ctx.n ? "dH/dF" : "dH/dw") << "[" << q << "] = " << std::scientific << g0[q]
                << " fd = " << fd << " rel = " << std::abs(fd - g0[q]) / std::max(1e-12, std::abs(fd) + std::abs(g0[q]))
                << std::defaultfloat << std::endl;
    }
    applyFbar(ctx, ctx.F_base);
  }

  // Stiffness preconditioner (--precond=stiffness, from --precond-from-step): bordered by the free
  // cell components, see StressPreconditioner. MTM_STRESS_PLAIN=1 forces the plain solver.
  alglib::minlbfgsstate state;
  alglib::minlbfgsreport rep;
  const bool use_precond = g_relaxation_solver && g_relaxation_step >= g_relaxation_first_step &&
                           g_relaxation_solver->options().type == LBFGSPreconditioner::Stiffness &&
                           std::getenv("MTM_STRESS_PLAIN") == nullptr;
  if (!use_precond ||
      !minimize_stress_enthalpy_preconditioned(z, ctx, corrections, maxits, g_relaxation_solver->options(), rep)) {
    // Plain L-BFGS. The homogeneous energy is periodic in shear under Lagrange reduction, so
    // H = E - A0 Pt : F is unbounded along the box components: cap the step (0.5 a0) so the line
    // search cannot extrapolate past the shear barrier.
    if (use_precond) std::cerr << "[Stress Control] preconditioner factorization failed; plain L-BFGS" << std::endl;
    alglib::minlbfgscreate(corrections, z, state);
    alglib::minlbfgssetcond(state, 0.0, 0.0, 0.0, maxits);
    alglib::minlbfgssetstpmax(state, 0.5);
    alglib::minlbfgsoptimize(state, stress_enthalpy, nullptr, &ctx);
    alglib::minlbfgsresults(state, z, rep);
  }
  if (std::getenv("MTM_STRESS_TRACE") != nullptr) {
    alglib::real_1d_array gz;
    gz.setlength(z.length());
    double hz = 0.0;
    stress_enthalpy(z, hz, gz, &ctx);
    double gF = 0.0, gW = 0.0;
    for (int q = 0; q < z.length(); ++q) (q >= 2 * ctx.n ? gF : gW) = std::max(q >= 2 * ctx.n ? gF : gW, std::abs(gz[q]));
    std::cout << "[Stress Relax] alpha = " << userData->current_alpha << " iterations = " << rep.iterationscount
              << " termination = " << rep.terminationtype << " H = " << std::scientific << hz
              << " max|dH/dzF| = " << gF << " max|dH/dw| = " << gW << std::defaultfloat << std::endl;
  }

  const Eigen::Matrix2d F = stressFbarFromVector(ctx, z);
  applyFbar(ctx, F);
  for (int i = 0; i < ctx.n; ++i) {
    const Eigen::Vector2d xi = F * Eigen::Vector2d(z[i], z[ctx.n + i]);
    x[i] = xi.x();
    x[ctx.n + i] = xi.y();
  }
}

// Cauchy (true-stress) target: the nominal target conjugate to the free box components depends on the
// box, Pt = J sigma_t F_bar^{-T}. Relax at fixed Pt, re-evaluate it at the relaxed box, repeat (each
// pass warm-starts from the last). The map contracts with ratio ~ sigma / modulus and stops converging
// only at the limit load, where the limit guard takes over.
void relax_with_cauchy_stress_control(alglib::real_1d_array &x, UserData *userData, int corrections, int maxits) {
  PBCStressControl *sc = userData->stress_ctrl;
  const double sigma_scale = sc->sigma_target.cwiseAbs().maxCoeff();
  const double tol = 1e-11 + 1e-9 * sigma_scale;
  constexpr int kMaxPasses = 60;
  const bool trace = std::getenv("MTM_STRESS_TRACE") != nullptr;
  double change = 0.0;
  int pass = 0;
  for (; pass < kMaxPasses; ++pass) {
    sc->setNominalFromCauchy(*sc->F_bar);
    const Eigen::Matrix2d P_used = sc->P_target;
    relax_with_stress_control(x, userData, corrections, maxits);
    sc->setNominalFromCauchy(*sc->F_bar);
    change = (sc->P_target - P_used).cwiseAbs().maxCoeff();
    if (trace) {
      std::cout << "[Cauchy Control] alpha = " << userData->current_alpha << " pass " << pass + 1
                << " |dPt| = " << std::scientific << change << std::defaultfloat << std::endl;
    }
    if (change <= tol) break;
    if (sc->exceedsLimit(*sc->F_bar)) return; // past the limit load: the loading loop stops the run
    sc->P_target = P_used; // keep the target the configuration is in equilibrium with
  }
  if (trace) {
    // Box Cauchy stress from the solver's own virial: sigma = (sum_e A_e P_e F_e^T) / (A0 det F_bar)
    alglib::real_1d_array g;
    g.setlength(x.length());
    double e = 0.0;
    Eigen::Matrix2d virial = Eigen::Matrix2d::Zero();
    userData->virial_out = &virial;
    minimize_energy_with_triangles(x, e, g, userData);
    userData->virial_out = nullptr;
    const Eigen::Matrix2d sig = virial / (sc->A0 * sc->F_bar->determinant());
    std::cout << "[Cauchy Control] alpha = " << userData->current_alpha << " box sigma (virial) xx = "
              << std::scientific << sig(0, 0) << " xy = " << sig(0, 1) << " yx = " << sig(1, 0)
              << " yy = " << sig(1, 1) << std::defaultfloat << std::endl;
  }
  if (pass == kMaxPasses) {
    std::cout << "[Cauchy Control] WARNING: alpha = " << userData->current_alpha << ": the nominal target did not "
              << "settle in " << kMaxPasses << " passes (|dPt| = " << std::scientific << change
              << std::defaultfloat << "); near the limit load." << std::endl;
  }
}

PBCStressControl g_pbc_stress_options;
} // namespace

void set_pbc_stress_control_options(const std::string &preset, bool relax_stage1, double relax_strain) {
  g_pbc_stress_options.preset = preset;
  g_pbc_stress_options.relax_stage1 = relax_stage1;
  g_pbc_stress_options.relax_strain = relax_strain;
}

const PBCStressControl &pbc_stress_control_options() { return g_pbc_stress_options; }

PBCBoxStability pbc_box_stability(UserData &state, const PBCStressControl &sc, double alpha, bool loaded) {
  using SpMat = Eigen::SparseMatrix<double>;
  PBCBoxStability out;
  std::vector<std::pair<int, int>> comps;
  for (int k = 0; k < 3; ++k)
    if (sc.comp_free[k]) comps.emplace_back(PBCStressControl::compRow(k), PBCStressControl::compCol(k));
  const int m = static_cast<int>(comps.size());
  const int n = static_cast<int>(state.interior_mapping.size());
  if (m == 0 || n < 2 || sc.A0 <= 0.0) return out;
  const Eigen::Matrix2d F = state.F_external;

  alglib::real_1d_array x;
  x.setlength(2 * n);
  map_points_to_solver_array(x, state.points, state.interior_mapping, n);
  for (auto &element : state.elements) element.setExternalDeformation(F);
  static FastStiffnessAssembler assembler;
  const SpMat &K = assembler.assembleStiffness(x, &state, false); // also caches F_e at x
  Eigen::MatrixXd K_xF, K_FF;
  assemble_box_border(&state, F, comps, false, n, K_xF, K_FF);

  // K is singular along the two translations and K_xF is orthogonal to them (the forces sum to zero),
  // so K^+ K_xF is exact with one node pinned (dofs 0 and n removed)
  auto reduced = [n](int i) { return i - (i > 0 ? 1 : 0) - (i > n ? 1 : 0); };
  std::vector<Eigen::Triplet<double>> trip;
  trip.reserve(static_cast<size_t>(K.nonZeros()));
  for (int c = 0; c < K.outerSize(); ++c)
    for (SpMat::InnerIterator it(K, c); it; ++it) {
      const int r = static_cast<int>(it.row()), cc = static_cast<int>(it.col());
      if (r == 0 || r == n || cc == 0 || cc == n) continue;
      trip.emplace_back(reduced(r), reduced(cc), it.value());
    }
  SpMat K_r(2 * n - 2, 2 * n - 2);
  K_r.setFromTriplets(trip.begin(), trip.end());
  Eigen::MatrixXd B_r(2 * n - 2, m);
  for (int i = 0; i < 2 * n; ++i)
    if (i != 0 && i != n) B_r.row(reduced(i)) = K_xF.row(i);
  Eigen::SimplicialLDLT<SpMat> ldlt(K_r);
  if (ldlt.info() != Eigen::Success) return out;
  const Eigen::MatrixXd Y = ldlt.solve(B_r);
  if (ldlt.info() != Eigen::Success) return out;
  Eigen::MatrixXd S = (K_FF - B_r.transpose() * Y) / sc.A0;
  S = 0.5 * (S + S.transpose());

  // Load stiffness and load direction dPt/dalpha
  const bool force = loaded && PBCStressControl::isForcePreset(sc.preset);
  const bool cauchy = force && PBCStressControl::isCauchyPreset(sc.preset);
  const int tr = PBCStressControl::targetRow(sc.preset), tc = PBCStressControl::targetCol(sc.preset);
  Eigen::Matrix2d pattern = Eigen::Matrix2d::Zero(); // d sigma_t / d alpha (Cauchy) or d Pt / d alpha
  pattern(tr, tc) = 1.0;
  pattern(tc, tr) = 1.0;
  auto Pt = [](const Eigen::Matrix2d &Fb, const Eigen::Matrix2d &sigma) {
    return Fb.determinant() * sigma * Fb.inverse().transpose();
  };
  Eigen::MatrixXd J = S;
  Eigen::VectorXd b = Eigen::VectorXd::Zero(m);
  if (cauchy) {
    const Eigen::Matrix2d sigma = alpha * pattern;
    const double h = 1e-6;
    for (int q = 0; q < m; ++q) {
      Eigen::Matrix2d Fp = F, Fm = F;
      Fp(comps[q].first, comps[q].second) += h;
      Fm(comps[q].first, comps[q].second) -= h;
      const Eigen::Matrix2d dP = (Pt(Fp, sigma) - Pt(Fm, sigma)) / (2.0 * h);
      for (int a = 0; a < m; ++a) J(a, q) -= dP(comps[a].first, comps[a].second);
    }
    const Eigen::Matrix2d dP_dalpha = Pt(F, pattern);
    for (int a = 0; a < m; ++a) b(a) = dP_dalpha(comps[a].first, comps[a].second);
  } else if (force) {
    for (int a = 0; a < m; ++a)
      if (comps[a].first == tr && comps[a].second == tc) b(a) = 1.0;
  }

  if (cauchy) {
    Eigen::EigenSolver<Eigen::MatrixXd> eig(J, false);
    if (eig.info() == Eigen::Success) {
      out.lambda_min = eig.eigenvalues().real().minCoeff();
      out.has_lambda = true;
    }
  } else {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(J, Eigen::EigenvaluesOnly);
    if (eig.info() == Eigen::Success) {
      out.lambda_min = eig.eigenvalues().minCoeff();
      out.has_lambda = true;
    }
  }
  if (force) {
    int t = -1;
    for (int a = 0; a < m; ++a)
      if (comps[a].first == tr && comps[a].second == tc) t = a;
    if (t >= 0) {
      const Eigen::VectorXd dF_dalpha = J.fullPivLu().solve(b);
      if (std::abs(dF_dalpha(t)) > 0.0) {
        out.tangent = 1.0 / dF_dalpha(t);
        out.has_tangent = true;
      }
    }
  }
  return out;
}

void relax_configuration(alglib::real_1d_array &x, UserData *userData, int corrections) {
  if (userData->macro_geom_ptr != nullptr && userData->macro_geom_ptr->hasForceGrips()) {
    relax_with_force_grips(x, userData, corrections, g_relaxation_maxits);
    return;
  }
  if (userData->stress_ctrl != nullptr && userData->stress_ctrl->active &&
      userData->stress_ctrl->F_bar != nullptr) {
    if (userData->stress_ctrl->cauchy) {
      relax_with_cauchy_stress_control(x, userData, corrections, g_relaxation_maxits);
    } else {
      relax_with_stress_control(x, userData, corrections, g_relaxation_maxits);
    }
    return;
  }
  if (g_relaxation_solver && g_relaxation_step >= g_relaxation_first_step) {
    g_relaxation_solver->optimize(x, userData);
    return;
  }
  LBFGSOptimizer optimizer(corrections, 0.0, 0.0, 0.0, g_relaxation_maxits);
  optimizer.optimize(x, minimize_energy_with_triangles, userData);
}

void writeSizesToFile(int Nx, int Ny) {
  std::ofstream file("sizes.dat");

  if (file.is_open()) {
    file << Nx << std::endl;
    file << Ny << std::endl;
    file.close();
    std::cout << "Successfully wrote sizes to sizes.dat" << std::endl;
  } else {
    std::cerr << "Error: Unable to open sizes.dat for writing" << std::endl;
  }
}


Eigen::Matrix<double, 3, 2>
calculateShapeDerivatives(const Eigen::Vector2d &p1, const Eigen::Vector2d &p2,
                          const Eigen::Vector2d &p3) {

  // Create matrix of nodal coordinates as columns
  Eigen::Matrix<double, 2, 3> X;
  X.col(0) = p1;
  X.col(1) = p2;
  X.col(2) = p3;

  // Natural derivatives (in reference element)
  Eigen::Matrix<double, 3, 2> dNdxi;
  dNdxi << -1.0, -1.0, // dN1/dξ, dN1/dη
      1.0, 0.0,        // dN2/dξ, dN2/dη
      0.0, 1.0;        // dN3/dξ, dN3/dη

  // Calculate Jacobian matrix
  Eigen::Matrix2d J = X * dNdxi;

  // Verify Jacobian is not singular
  double det = J.determinant();
  if (std::abs(det) < 1e-10) {
    throw std::runtime_error(
        "Near-singular Jacobian detected in shape function calculation");
  }

  // Calculate physical derivatives
  Eigen::Matrix<double, 3, 2> dNdX = dNdxi * J.inverse();

  return dNdX;
}
// Function to scale a lattice by a factor
std::vector<Point2D> scaleLattice(const std::vector<Point2D> &original_points,
                                  double scale_factor) {
  std::vector<Point2D> scaled_points = original_points;

  // Scale each point by the factor
  for (auto &point : scaled_points) {
    point.coord *= scale_factor;
  }

  return scaled_points;
}

std::vector<Point2D>
scaleLatticeAroundPoint(const std::vector<Point2D> &original_points,
                        double scale_factor,
                        const Eigen::Vector2d &reference_point) {

  std::vector<Point2D> scaled_points = original_points;

  for (auto &point : scaled_points) {
    // Shift to origin
    point.coord -= reference_point;
    // Scale
    point.coord *= scale_factor;
    // Shift back
    point.coord += reference_point;
  }

  return scaled_points;
}

std::tuple<double, Eigen::Matrix2d, int> perform_remeshing_loop(
    alglib::real_1d_array &x, UserData *userData,
    const std::vector<int> &contact_atoms,
    const std::vector<int> &boundary_fixed_nodes, const Eigen::Matrix2d &F_ext,
    const Eigen::Matrix<double, 3, 2> &dndx,
    const std::array<double, 2> &offsets,
    const std::vector<int> &original_domain_map,
    const std::vector<std::tuple<double, double>> &translation_map,
    const Point2D &domain_dims_point, int max_iterations = 100,
    double reference_area = 0.5) {
  bool should_remesh = true;
  int mesh_iteration = 0;
  double final_energy = 0.0;
  double final_stress = 0.0;
  Eigen::Matrix2d stress_tensor;
  // const int n_vars = x.length();

  std::vector<Point2D> &square_points = userData->points;
  std::vector<ElementTriangle2D> &elements = userData->elements;
  std::vector<size_t> &active_elements = userData->active_elements;
  const auto &interior_mapping = userData->interior_mapping;
  const auto &full_mapping = userData->full_mapping;
  const int n_vars = interior_mapping.size();
  const int n_points = userData->points.size();
  // const double normalisation = pow(userData->ideal_lattice_parameter, 2.0);
  std::function<double(double)> &potential_func = userData->energy_function;
  std::function<double(double)> &potential_func_der =
      userData->derivative_function;
  double zero = userData->zero_energy;
  double ideal_lattice_parameter = userData->ideal_lattice_parameter;
  double plasticity;
  TriangularLatticeCalculator calculator(ideal_lattice_parameter);

  std::cout << "REMESHING STARTED " << std::endl;

  std::set<std::vector<std::array<int, 3>>> visited_topologies;
  visited_topologies.insert(getTopologyKeySet(elements, active_elements));

  while (should_remesh && mesh_iteration < max_iterations) {
    std::cout << "REMESHING iteration: " << mesh_iteration << std::endl;

    // 1. Save original state
    alglib::real_1d_array original_x = x;
    // auto m3_before = analyzeElementReduction(elements, square_points,
    // &userData);

    // 2. Generate new mesh (using existing mesher)
    // 1. Create the AdaptiveMesher instance
    AdaptiveMesher mesher(domain_dims_point, offsets, original_domain_map,
                          translation_map, full_mapping,
                          1e-6 // Tolerance
    );
    mesher.setUsePeriodicCopies(true); // Switch to using original domain only

    alglib::real_1d_array original_x_remesh = mesher.saveOriginalPositions(x);
    std::tie(elements, active_elements) =
        mesher.createMesh(square_points, x, F_ext, &dndx);
    for (auto &element : elements) {
      element.setReferenceArea(
          reference_area); // or interior_mapping depending on needs
    }

    // It is called to find the number of nodes inside elements that touch the
    // boundary This is requited in mesh filtering auto [interior_mapping_dummy,
    // full_mapping_dummy] = create_dof_mapping_with_boundaries(
    //     square_points, elements,contact_atoms,boundary_fixed_nodes);

    // // 3. Re-optimize with new mesh
    UserData newUserData(square_points, elements, calculator, potential_func,
                         potential_func_der, zero, ideal_lattice_parameter,
                         F_ext, interior_mapping, full_mapping, active_elements,
                         plasticity,
                         userData->second_derivative_function);

    relax_configuration(x, &newUserData, 13);
    map_solver_array_to_points(x, square_points, interior_mapping, n_vars);

    // // 4. Check convergence
    auto change_result = computeChangeMeasures(
        x, original_x_remesh, userData->ideal_lattice_parameter, elements,
        &newUserData, square_points, true, &F_ext);
    should_remesh = change_result.has_distorted_triangles;

    if (!should_remesh) {

      ConfigurationSaver::calculateEnergyAndStress(&newUserData, final_energy,
                                                   stress_tensor, true);
      final_stress = stress_tensor(0, 1);
      break;
    }

    mesh_iteration++;
  }

  return {final_energy, stress_tensor, mesh_iteration};
}
