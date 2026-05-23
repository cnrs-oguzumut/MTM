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

  std::cout << "REMESHING STARTED " << std::endl;

  while (should_remesh && mesh_iteration < max_iterations) {
    std::cout << "REMESHING iteration: " << mesh_iteration << std::endl;

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
                         full_mapping, active_elements, plasticity);

    double energy_old = 0.0;
    Eigen::Matrix2d stress_old;
    ConfigurationSaver::calculateEnergyAndStress(&oldUserData, energy_old,
                                                 stress_old, true);

    std::cout << "Energy of OLD state: " << std::scientific
              << std::setprecision(8) << energy_old << std::endl;

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

    UserData newUserData(square_points, elements, userData->calculator,
                         potential_func, potential_func_der, zero,
                         ideal_lattice_parameter, F_ext, interior_mapping,
                         full_mapping, active_elements, plasticity);

    // === OPTIMIZE ON NEW MESH ===
    std::cout << "Optimization in REMESHING loop" << std::endl;
    LBFGSOptimizer optimizer(13, 0, 0, 0, 0);
    if (optimize_interior)
      optimizer.optimize(x, minimize_energy_with_triangles, &newUserData);

    map_solver_array_to_points(x, square_points, interior_mapping, n_vars);

    // ============================================
    // ENERGY OF NEW STATE (after remeshing + optimization)
    // ============================================
    UserData newUserData_after(square_points, elements, userData->calculator,
                               potential_func, potential_func_der, zero,
                               ideal_lattice_parameter, F_ext, interior_mapping,
                               full_mapping, active_elements, plasticity);

    double energy_new = 0.0;
    Eigen::Matrix2d stress_new;
    ConfigurationSaver::calculateEnergyAndStress(&newUserData_after, energy_new,
                                                 stress_new, true);

    // ============================================
    // COMPARE OLD vs NEW
    // ============================================
    double energy_change = energy_new - energy_old;
    double energy_change_percent =
        (energy_old != 0.0) ? (energy_change / std::abs(energy_old)) * 100.0
                            : 0.0;

    std::cout << "\n┌─────────────────────────────────────────────────────┐"
              << std::endl;
    std::cout << "│  REMESH ENERGY COMPARISON (Iteration " << std::setw(3)
              << mesh_iteration << ")       │" << std::endl;
    std::cout << "├─────────────────────────────────────────────────────┤"
              << std::endl;
    std::cout << "│  Energy OLD:     " << std::setw(24) << std::scientific
              << std::setprecision(8) << energy_old << "  │" << std::endl;
    std::cout << "│  Energy NEW:     " << std::setw(24) << std::scientific
              << std::setprecision(8) << energy_new << "  │" << std::endl;
    std::cout << "│  Change (ΔE):    " << std::setw(24) << std::scientific
              << std::setprecision(8) << energy_change << "  │" << std::endl;
    std::cout << "│  Change (%):     " << std::setw(23) << std::fixed
              << std::setprecision(6) << energy_change_percent << "%  │"
              << std::endl;
    std::cout << "└─────────────────────────────────────────────────────┘\n"
              << std::endl;

    if (energy_change >= 0) {
      std::cout << "⚠️  REJECTING remesh - energy increased!" << std::endl;
      std::cout << "    Restoring previous state..." << std::endl;

      // RESTORE OLD STATE
      x = old_x;
      square_points = old_points;
      elements = old_elements;
      active_elements = old_active_elements;

      // Final values from restored state
      final_energy = energy_old;
      stress_tensor = stress_old;

      // Stop remeshing - we couldn't improve
      should_remesh = false;
      has_changes = 0;

    } else {
      std::cout << "✓ ACCEPTING remesh - energy decreased" << std::endl;

      // Keep new state (already in place)
      final_energy = energy_new;
      stress_tensor = stress_new;
      has_changes = 1;

      // Continue to see if we can improve further
      // (or set should_remesh = false if you only want one successful remesh)
    }

    mesh_iteration++;
  }

  return {final_energy, stress_tensor, mesh_iteration};
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
                         plasticity);

    LBFGSOptimizer optimizer(13, 0, 0, 0, 0);
    optimizer.optimize(x, minimize_energy_with_triangles, &newUserData);
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
