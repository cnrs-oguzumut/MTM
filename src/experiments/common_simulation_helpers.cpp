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
  g_relaxation_solver = std::make_unique<PreconditionedLBFGS>(opts);
}

void relaxation_begin_step(int step) { g_relaxation_step = step; }

void relax_configuration(alglib::real_1d_array &x, UserData *userData, int corrections) {
  if (g_relaxation_solver && g_relaxation_step >= g_relaxation_first_step) {
    g_relaxation_solver->optimize(x, userData);
    return;
  }
  LBFGSOptimizer optimizer(corrections, 0.0, 0.0, 0.0, 0);
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
