#include "../../include/experiments/experiment_includes.h"

// ============================================================================
// single_dislocation_study
//
// Stripped-down version of the Conti-Zanzotto square-crystal experiment.
// It keeps the original setup pipeline (reference geometry, energy potential,
// lattice generation, periodic domain maps, DOF mapping, meshing and the
// strain-energy calculator) but REMOVES the external loading schedule, the
// initial random noise and the avalanche/stress-drop bookkeeping.
//
// Instead it installs a single edge dislocation (Volterra field) at the
// central atom, performs an initial L-BFGS relaxation followed by the
// remesh/minimize loop (remesh -> minimize -> remesh -> ...) from the
// Conti-Zanzotto experiment, and writes the resulting state to VTK together
// with a defect analysis.
// ============================================================================
void single_dislocation_study(int caller_id, int nx, int ny) {

  // Toggle for the remesh/minimize loop after the initial relaxation.
  const bool remeshing = true;

  // Parameters for lattice
  if (nx <= 0 || ny <= 0) {
    std::cerr << "Error: nx and ny must be positive integers." << std::endl;
    exit(EXIT_FAILURE);
  }

  // ==================== SETUP REFERENCE GEOMETRY ====================
  writeSizesToFile(nx, ny);

  std::string lattice_type = "square"; // Options: "square" or "triangular"
  double h = 1.0;                      // Reference element size

  // Reference triangle vertices (for shape function derivatives)
  Eigen::Vector2d p1(0.0, 0.0);
  Eigen::Vector2d p2(h, 0.0);
  Eigen::Vector2d p3(0.0, h);

  Eigen::Matrix<double, 3, 2> dndx = calculateShapeDerivatives(p1, p2, p3);
  std::cout << "Reference element shape derivatives (dN/dx):\n"
            << dndx << std::endl;

  // ==================== SETUP ENERGY POTENTIAL ====================
  std::function<double(double)> potential_func = square_energy;
  std::function<double(double)> potential_func_der = square_energy_der;
  std::function<double(double)> potential_func_sder = square_energy_der;

  // ==================== DETERMINE OPTIMAL LATTICE PARAMETER ============
  double symmetry_constantx =
      (lattice_type == "triangular") ? pow(4.0 / 3.0, 1.0 / 4.0) : 1.0;
  double optimal_lattice_parameter = symmetry_constantx * 1.0;
  double lattice_constant = optimal_lattice_parameter;
  std::cout << "Optimal lattice parameter: " << lattice_constant << std::endl;

  // ==================== GENERATE INITIAL LATTICE ====================
  std::vector<Point2D> square_points = LatticeGenerator::generate_2d_lattice(
      nx, ny, lattice_constant, lattice_type);

  std::vector<Point2D> square_points_ref =
      LatticeGenerator::generate_2d_lattice(nx, ny, lattice_constant,
                                            lattice_type);

  int original_domain_size = square_points.size();
  std::cout << "Generated lattice with " << original_domain_size << " points"
            << std::endl;

  // Calculate domain properties
  DomainInfo domain_size = compute_domain_size(square_points);

  const std::array<double, 2> offsets =
      (lattice_type == "square")
          ? std::array<double, 2>{lattice_constant, lattice_constant}
          : std::array<double, 2>{lattice_constant / 2.0,
                                  (sqrt(3.0) / 2.0) * lattice_constant};
  std::cout << "PBC offsets: [" << offsets[0] << ", " << offsets[1] << "]"
            << std::endl;

  DomainDimensions domain_dims(domain_size.get_width(),
                               domain_size.get_height());

  bool pbc = false;

  // Create domain maps
  auto [original_domain_map, translation_map] =
      MeshGenerator::create_domain_maps(original_domain_size, domain_dims,
                                        offsets);

  // Boundary conditions / DOF mapping
  auto [interior_mapping, full_mapping] =
      create_dof_mapping_original(square_points, 0.5 * lattice_constant, pbc);
  std::cout << "interior_mapping.size(): " << interior_mapping.size()
            << std::endl;
  std::cout << "full_mapping.size(): " << full_mapping.size() << std::endl;

  Point2D domain_dims_point(domain_dims.size_x, domain_dims.size_y);

  // ==================== BUILD MESH (on the pristine lattice) ============
  // Mesh is generated BEFORE introducing the dislocation, so the connectivity
  // and reference configuration come from the perfect crystal.
  AdaptiveMesher mesher(domain_dims_point, offsets, original_domain_map,
                        translation_map, full_mapping,
                        1e-6, // Tolerance
                        pbc   // Use periodic copies
  );
  mesher.setUsePeriodicCopies(pbc);

  alglib::real_1d_array free_dofs;
  int n_free_nodes = interior_mapping.size();
  free_dofs.setlength(2 * interior_mapping.size());
  map_points_to_solver_array(free_dofs, square_points, interior_mapping,
                             n_free_nodes);

  alglib::real_1d_array original_x_remesh =
      mesher.saveOriginalPositions(free_dofs);

  auto [elements, active_elements] = mesher.createMesh(
      square_points, free_dofs, Eigen::Matrix2d::Identity(), &dndx);
  double element_area = elements[0].getReferenceArea();

  for (auto &element : elements) {
    element.set_dof_mapping(full_mapping);
  }
  std::cout << "Created " << elements.size() << " element triangles"
            << std::endl;

  // ==================== INSTALL A SINGLE DISLOCATION ====================
  // Applied after meshing: displaces the (already meshed) atoms by the
  // Volterra edge-dislocation field.
  Eigen::Vector2d burgers_vector(lattice_constant, 0.0); // unit Burgers, x-dir
  double core_radius = 6.0;                              // core regularization
  double poisson_ratio = 0.3;                            // material Poisson
  size_t middle_atom_index = findMiddleAtom(square_points, true);

  std::vector<Point2D> dislocated_points =
      createSingleDislocation(square_points_ref, burgers_vector,
                              middle_atom_index, core_radius, poisson_ratio);

  // Use the dislocated configuration as the current state.
  square_points = dislocated_points;
  std::cout << "Installed single edge dislocation at atom " << middle_atom_index
            << " (|b| = " << burgers_vector.norm() << ")" << std::endl;

  // ==================== SETUP ENERGY CALCULATION ====================
  Strain_Energy_LatticeCalculator calculator(1.0);

  Eigen::Matrix2d F_I = Eigen::Matrix2d::Identity();
  Eigen::Matrix2d C_I = F_I.transpose() * F_I; // C = F^T * F
  double zero = calculator.calculate_energy(C_I, potential_func, 0);
  std::cout << "zero energy value: " << zero << std::endl;

  // ==================== CREATE USER DATA (no external load) ============
  Eigen::Matrix2d F_ext = Eigen::Matrix2d::Identity();
  bool plasticity = false;
  UserData userData(square_points, elements, calculator, potential_func,
                    potential_func_der, zero, optimal_lattice_parameter, F_ext,
                    interior_mapping, full_mapping, active_elements,
                    plasticity);

  // ==================== PREPARE OPTIMIZATION ====================
  alglib::real_1d_array x;
  int n_vars = interior_mapping.size();
  x.setlength(2 * n_vars);
  map_points_to_solver_array(x, square_points, interior_mapping, n_vars);

  // Pre-relaxation energy / stress
  double pre_energy = 0.0;
  Eigen::Matrix2d stress_tensor = Eigen::Matrix2d::Zero();
  ConfigurationSaver::calculateEnergyAndStress(&userData, pre_energy,
                                               stress_tensor, true);
  double pre_stress = stress_tensor(0, 1);
  std::cout << "Pre-relaxation - Energy: " << pre_energy
            << ", Stress: " << pre_stress << std::endl;

  // ==================== SAVE INITIAL DISLOCATED STATE ====================
  int file_id = caller_id;
  double saving_value = 0.0;

  ConfigurationSaver::saveConfigurationWithStressAndEnergy2D(
      &userData, file_id, pre_energy, pre_stress, true);
  ConfigurationSaver::saveTriangleData(&userData, file_id, domain_dims, offsets,
                                       full_mapping);
  ConfigurationSaver::saveElements(elements, active_elements, file_id);

  auto [num_dislocations_pre, coordination_pre] =
      DefectAnalysis::analyzeDefectsInReferenceConfig(
          &userData, file_id, dndx, offsets, original_domain_map,
          translation_map, domain_dims_point, element_area, pbc, true);

  ConfigurationSaver::writeToVTK(userData.points, userData.elements, &userData,
                                 file_id, true, coordination_pre, saving_value);
  ConfigurationSaver::logDislocationData(0.0, num_dislocations_pre);
  std::cout << "Saved initial dislocated config " << file_id << std::endl;

  double pre_area = ConfigurationSaver::calculateTotalArea2D(&userData);

  // ==================== RELAX (initial L-BFGS minimization) ============
  userData.third_condition_flag = false;
  LBFGSOptimizer optimizer(13, 0.00001, 0, 0, 0);
  optimizer.optimize(x, minimize_energy_with_triangles, &userData);
  map_solver_array_to_points(x, square_points, interior_mapping, n_vars);

  // ==================== POST-RELAXATION ENERGY / STRESS ============
  double post_energy = 0.0;
  stress_tensor.setZero();
  ConfigurationSaver::calculateEnergyAndStress(&userData, post_energy,
                                               stress_tensor, true);
  double post_stress = stress_tensor(0, 1);
  std::cout << "Post-relaxation - Energy: " << post_energy
            << ", Stress: " << post_stress << std::endl;

  // ==================== REMESHING DECISION ====================
  // Energy-based decision, as in the Conti-Zanzotto loop. The baseline is the
  // dislocated configuration energy before relaxation; this single-shot study
  // plays the role of the "first iteration" (i == 0), which always permits an
  // initial remeshing attempt.
  double post_energy_previous = pre_energy;
  bool first_pass = true;
  bool shouldRemesh = (post_energy < post_energy_previous) || first_pass;

  std::cout << "=== REMESH DECISION ===" << std::endl;
  std::cout << "  pre_energy = " << std::scientific << pre_energy << std::endl;
  std::cout << "  post_energy = " << std::scientific << post_energy
            << std::endl;
  std::cout << "  post_energy_previous = " << std::scientific
            << post_energy_previous << std::endl;
  std::cout << "  difference = " << (post_energy - post_energy_previous)
            << std::endl;
  std::cout << "  post_energy < post_energy_previous? "
            << (post_energy < post_energy_previous) << std::endl;
  std::cout << "  shouldRemesh = " << shouldRemesh << std::endl;

  if (remeshing && shouldRemesh) {
    std::cout << "REMESHING STARTS" << std::endl;

    std::vector<int> contact_atoms;
    std::vector<int> boundary_fixed_nodes;
    int hasChanges = 0;
    int max_iterations = 1000;

    // Iterative remesh/minimize loop (remesh -> minimize -> remesh -> ...).
    auto [post_energy_re, stress_tensor_re, iterations] =
        perform_remeshing_loop_reduction(
            x, &userData, contact_atoms, boundary_fixed_nodes, F_ext, dndx,
            offsets, original_domain_map, translation_map, domain_dims_point,
            hasChanges, max_iterations, element_area, pbc, true);

    post_energy = post_energy_re;
    post_stress = stress_tensor_re(0, 1);

    // Sync local containers with the (possibly) remeshed UserData state.
    square_points = userData.points;
    elements = userData.elements;
    active_elements = userData.active_elements;
    if (!elements.empty()) {
      element_area = elements[0].getReferenceArea();
    }
    for (auto &element : elements) {
      element.set_dof_mapping(full_mapping);
    }

    if (hasChanges) {
      std::cout << "Remeshing accepted - energy decreased (" << iterations
                << " iterations)" << std::endl;
    } else {
      std::cout << "Remeshing rejected - original mesh kept" << std::endl;
    }
    std::cout << "Final energy: " << post_energy << ", stress: " << post_stress
              << std::endl;
  }

  // ==================== REPORT / LOG ENERGY DROP ====================
  double post_area = ConfigurationSaver::calculateTotalArea2D(&userData);
  std::cout << "Energy dropped from " << post_energy_previous << " to "
            << post_energy << std::endl;
  ConfigurationSaver::logEnergyAndStress_v2(0, 0.0, pre_energy, pre_stress,
                                            post_energy, post_stress, pre_area,
                                            post_area, shouldRemesh);

  // ==================== SAVE RELAXED STATE ====================
  int post_file_id = caller_id + 1;
  UserData postOptUserData(square_points, elements, calculator, potential_func,
                           potential_func_der, zero, optimal_lattice_parameter,
                           F_ext, interior_mapping, full_mapping,
                           active_elements, plasticity);

  ConfigurationSaver::saveTriangleData(&postOptUserData, post_file_id,
                                       domain_dims, offsets, full_mapping);
  ConfigurationSaver::saveElements(elements, active_elements, post_file_id);

  auto [num_dislocations_post, coordination_post] =
      DefectAnalysis::analyzeDefectsInReferenceConfig(
          &postOptUserData, post_file_id, dndx, offsets, original_domain_map,
          translation_map, domain_dims_point, element_area, pbc, true);

  ConfigurationSaver::writeToVTK(
      postOptUserData.points, postOptUserData.elements, &postOptUserData,
      post_file_id, true, coordination_post, saving_value);
  ConfigurationSaver::logDislocationData(0.0, num_dislocations_post);
  std::cout << "Saved relaxed dislocated config " << post_file_id << std::endl;
}
