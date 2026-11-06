#include "../../include/geometry/MacroscopicGeometry.h"
#include "../../include/experiments/pbc_stress_control.h"
#include "../../include/experiments/experiment_includes.h"

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

  // Plain L-BFGS (the stiffness preconditioner assembles node DOFs only). The homogeneous energy is
  // periodic in shear under Lagrange reduction, so H = E - A0 Pt : F is unbounded along the box
  // components: cap the step (0.5 a0) so the line search cannot extrapolate past the shear barrier.
  alglib::minlbfgsstate state;
  alglib::minlbfgsreport rep;
  alglib::minlbfgscreate(corrections, z, state);
  alglib::minlbfgssetcond(state, 0.0, 0.0, 0.0, maxits);
  alglib::minlbfgssetstpmax(state, 0.5);
  alglib::minlbfgsoptimize(state, stress_enthalpy, nullptr, &ctx);
  alglib::minlbfgsresults(state, z, rep);
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

void set_pbc_stress_control_options(const std::string &preset, bool relax_stage1) {
  g_pbc_stress_options.preset = preset;
  g_pbc_stress_options.relax_stage1 = relax_stage1;
}

const PBCStressControl &pbc_stress_control_options() { return g_pbc_stress_options; }

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
