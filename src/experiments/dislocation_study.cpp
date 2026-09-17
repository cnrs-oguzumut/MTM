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
  double poisson_ratio = 1.0 / 3.0;                      // material Poisson
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

// ============================================================================
// single_dislocation_cylinder_relaxation
//
// Installs an analytical Volterra edge dislocation on an nx x ny lattice:
// - All atoms are initially given the analytical Volterra displacement field.
// - Atoms outside radius R_free from the core are FROZEN (Dirichlet boundary).
// - Atoms inside radius R_free are FREE to relax by energy minimization.
// - Supports L-BFGS preconditioning (--precond=stiffness/laplacian/diag/none).
// ============================================================================
void single_dislocation_cylinder_relaxation(int caller_id, int nx, int ny,
                                           double R_free,
                                           bool enable_remeshing,
                                           bool export_full_mesh) {
  if (nx <= 0 || ny <= 0) {
    std::cerr << "Error: nx and ny must be positive integers." << std::endl;
    exit(EXIT_FAILURE);
  }

  writeSizesToFile(nx, ny);

  std::string lattice_type = "square";
  double h = 1.0;

  // Reference triangle shape derivatives
  Eigen::Vector2d p1(0.0, 0.0);
  Eigen::Vector2d p2(h, 0.0);
  Eigen::Vector2d p3(0.0, h);
  Eigen::Matrix<double, 3, 2> dndx = calculateShapeDerivatives(p1, p2, p3);

  // Energy potential
  std::function<double(double)> potential_func = square_energy;
  std::function<double(double)> potential_func_der = square_energy_der;

  double optimal_lattice_parameter = 1.0;
  double lattice_constant = optimal_lattice_parameter;

  // Generate pristine reference lattice
  std::vector<Point2D> square_points_ref =
      LatticeGenerator::generate_2d_lattice(nx, ny, lattice_constant, lattice_type);
  std::vector<Point2D> square_points = square_points_ref;
  int n_points = square_points.size();

  DomainInfo domain_size = compute_domain_size(square_points_ref);
  const std::array<double, 2> offsets = {lattice_constant, lattice_constant};
  DomainDimensions domain_dims(domain_size.get_width(), domain_size.get_height());
  Point2D domain_dims_point(domain_dims.size_x, domain_dims.size_y);
  bool pbc = false;

  auto [original_domain_map, translation_map] =
      MeshGenerator::create_domain_maps(n_points, domain_dims, offsets);

  // Calculate bounding box and place dislocation core at center on slip plane
  double x_min = std::numeric_limits<double>::max();
  double x_max = std::numeric_limits<double>::lowest();
  double y_min = std::numeric_limits<double>::max();
  double y_max = std::numeric_limits<double>::lowest();
  for (const auto &p : square_points_ref) {
    x_min = std::min(x_min, p.coord.x());
    x_max = std::max(x_max, p.coord.x());
    y_min = std::min(y_min, p.coord.y());
    y_max = std::max(y_max, p.coord.y());
  }
  double core_x = 0.5 * (x_min + x_max);
  double core_y = 0.5 * (y_min + y_max);
  Eigen::Vector2d dislocation_core(core_x, core_y);

  std::cout << "\n=== VOLTERRA DISLOCATION CYLINDER RELAXATION ===" << std::endl;
  std::cout << "Lattice: " << nx << " x " << ny << " (" << n_points << " atoms)" << std::endl;
  std::cout << "Dislocation core at: (" << core_x << ", " << core_y << ")" << std::endl;
  std::cout << "Relaxation radius R_free: " << R_free << " h" << std::endl;

  // Apply analytical anisotropic Stroh / Eshelby-Read-Shockley displacement to all atoms
  Eigen::Vector2d burgers_vector(lattice_constant, 0.0);

  for (size_t i = 0; i < square_points.size(); ++i) {
    Eigen::Vector2d rel_pos = square_points_ref[i].coord - dislocation_core;
    Eigen::Vector2d u_dislocation = VolterraDisplacement::calculateAnisotropicEdgeDisplacement(
        rel_pos, burgers_vector);
    square_points[i].coord = square_points_ref[i].coord + u_dislocation;
  }

  // Partition into free nodes (r <= R_free) and fixed nodes (r > R_free)
  std::vector<std::pair<int, int>> interior_mapping;
  std::vector<std::pair<int, int>> full_mapping;
  interior_mapping.reserve(n_points);
  full_mapping.reserve(n_points);

  std::vector<int> fixed_nodes;
  int solver_idx = 0;
  for (size_t i = 0; i < square_points.size(); ++i) {
    double dist = (square_points_ref[i].coord - dislocation_core).norm();
    if (dist <= R_free) {
      interior_mapping.push_back({static_cast<int>(i), solver_idx});
      full_mapping.push_back({static_cast<int>(i), solver_idx});
      solver_idx++;
    } else {
      full_mapping.push_back({static_cast<int>(i), -1});
      fixed_nodes.push_back(static_cast<int>(i));
    }
  }

  int n_free_nodes = interior_mapping.size();
  int n_fixed_nodes = fixed_nodes.size();
  std::cout << "Partition: " << n_free_nodes << " free (relaxed) nodes, "
            << n_fixed_nodes << " frozen (Volterra Dirichlet) nodes." << std::endl;

  // Build mesh with full DOF mapping
  AdaptiveMesher mesher(domain_dims_point, offsets, original_domain_map,
                        translation_map, full_mapping, 1e-6, pbc);
  mesher.setUsePeriodicCopies(pbc);

  alglib::real_1d_array x;
  int n_vars = interior_mapping.size();
  x.setlength(2 * n_vars);
  map_points_to_solver_array(x, square_points, interior_mapping, n_vars);

  alglib::real_1d_array x_ref;
  x_ref.setlength(2 * n_vars);
  map_points_to_solver_array(x_ref, square_points_ref, interior_mapping, n_vars);

  // Both simulations start from the exact same pristine crystal mesh,
  // exactly like the avalanche / shifted crystal simulations.
  // When enable_remeshing is true, perform_remeshing_loop_reduction will perform
  // successive remesh-relaxation cycles on the relaxed configuration.
  auto [elements, active_elements] = mesher.createMesh(
      square_points_ref, x_ref,
      Eigen::Matrix2d::Identity(), &dndx);
  double element_area = elements.empty() ? 0.0 : elements[0].getReferenceArea();

  for (auto &element : elements) {
    element.set_reference_mesh(square_points);
    element.set_dof_mapping(full_mapping);
  }
  std::cout << "Mesh created: " << elements.size() << " triangular elements ("
            << (enable_remeshing ? "reconnected Delaunay" : "fixed pristine lattice") << ")." << std::endl;

  // Setup energy calculation
  Strain_Energy_LatticeCalculator calculator(1.0);
  Eigen::Matrix2d F_I = Eigen::Matrix2d::Identity();
  Eigen::Matrix2d C_I = F_I.transpose() * F_I;
  double zero = calculator.calculate_energy(C_I, potential_func, 0);

  Eigen::Matrix2d F_ext = Eigen::Matrix2d::Identity();
  bool plasticity = false;
  UserData userData(square_points, elements, calculator, potential_func,
                    potential_func_der, zero, optimal_lattice_parameter, F_ext,
                    interior_mapping, full_mapping, active_elements, plasticity);

  // Compute pre-relaxation energy and stress
  double pre_energy = 0.0;
  Eigen::Matrix2d stress_tensor = Eigen::Matrix2d::Zero();
  ConfigurationSaver::calculateEnergyAndStress(&userData, pre_energy, stress_tensor, true);
  double pre_stress = stress_tensor(0, 1);
  double pre_area = ConfigurationSaver::calculateTotalArea2D(&userData);

  std::cout << "State 0 (Unrelaxed Volterra): Energy = " << pre_energy
            << ", Stress_12 = " << pre_stress << std::endl;

  // Prepare mesh indices for visualization based on export_full_mesh:
  // export_full_mesh = true  -> exports entire crystal (all elements)
  // export_full_mesh = false -> exports only the active elements (cylinder core)
  std::vector<size_t> vis_elements;
  std::vector<std::pair<int, int>> vis_full_mapping;
  if (export_full_mesh) {
    vis_elements.resize(elements.size());
    std::iota(vis_elements.begin(), vis_elements.end(), 0);
    vis_full_mapping.resize(n_points);
    for (int i = 0; i < n_points; ++i) {
      vis_full_mapping[i] = {i, i};
    }
  } else {
    vis_elements = active_elements;
    vis_full_mapping = full_mapping;
  }

  // Save State 0 (Unrelaxed Volterra state)
  int initial_file_id = caller_id;
  UserData visUserDataPre(square_points, elements, calculator, potential_func,
                          potential_func_der, zero, optimal_lattice_parameter, F_ext,
                          interior_mapping, vis_full_mapping, vis_elements, plasticity);
  ConfigurationSaver::saveConfigurationWithStressAndEnergy2D(
      &visUserDataPre, initial_file_id, pre_energy, pre_stress, true);
  ConfigurationSaver::saveTriangleData(&visUserDataPre, initial_file_id, domain_dims, offsets, full_mapping);
  ConfigurationSaver::saveElements(elements, vis_elements, initial_file_id);

  auto [num_dislocations_pre, coordination_pre] =
      DefectAnalysis::analyzeDefectsInReferenceConfig(
          &visUserDataPre, initial_file_id, dndx, offsets, original_domain_map,
          translation_map, domain_dims_point, element_area, pbc, true);
  ConfigurationSaver::writeToVTK(visUserDataPre.points, visUserDataPre.elements, &visUserDataPre,
                                 initial_file_id, true, coordination_pre, 0.0);
  ConfigurationSaver::logDislocationData(0.0, num_dislocations_pre);
  std::cout << "Saved initial state as configuration_"
            << std::setw(5) << std::setfill('0') << initial_file_id << ".vtk ("
            << (export_full_mesh ? "full mesh" : "active elements only") << ")" << std::endl;

  // Energy minimization within R_free
  std::cout << "Minimizing energy of free nodes within R_free..." << std::endl;
  userData.third_condition_flag = false;
  relaxation_begin_step(0);
  relax_configuration(x, &userData, 13);

  map_solver_array_to_points(x, square_points, interior_mapping, n_vars);
  userData.points = square_points;

  // Compute post-relaxation energy and stress
  double post_energy = 0.0;
  stress_tensor.setZero();
  ConfigurationSaver::calculateEnergyAndStress(&userData, post_energy, stress_tensor, true);
  double post_stress = stress_tensor(0, 1);

  std::cout << "State 1 (Relaxed core): Energy = " << post_energy
            << " (drop: " << (pre_energy - post_energy) << ")"
            << ", Stress_12 = " << post_stress << std::endl;

  // Optional remeshing/reconnection loop around core
  if (enable_remeshing) {
    std::cout << "Running remeshing/reconnection around the relaxed core..." << std::endl;
    std::vector<int> contact_atoms;
    int max_iterations = 100;
    int hasChanges = 0;
    auto [post_energy_re, stress_tensor_re, iterations] =
        perform_remeshing_loop_reduction(
            x, &userData, contact_atoms, fixed_nodes, F_ext, dndx, offsets,
            original_domain_map, translation_map, domain_dims_point,
            hasChanges, max_iterations, element_area, pbc, true);

    post_energy = post_energy_re;
    post_stress = stress_tensor_re(0, 1);
    square_points = userData.points;
    elements = userData.elements;
    active_elements = userData.active_elements;
    if (!elements.empty()) {
      element_area = elements[0].getReferenceArea();
    }
    for (auto &element : elements) {
      element.set_dof_mapping(full_mapping);
    }
    std::cout << "Remeshing finished (" << iterations << " iterations, changes: "
              << hasChanges << ", final energy: " << post_energy << ")" << std::endl;
  }

  // Save State 1 (Relaxed state)
  int relaxed_file_id = caller_id + 1;
  std::vector<size_t> vis_elements_post;
  if (export_full_mesh) {
    vis_elements_post.resize(elements.size());
    std::iota(vis_elements_post.begin(), vis_elements_post.end(), 0);
  } else {
    vis_elements_post = active_elements;
  }
  UserData postOptUserData(square_points, elements, calculator, potential_func,
                           potential_func_der, zero, optimal_lattice_parameter,
                           F_ext, interior_mapping, vis_full_mapping, vis_elements_post, plasticity);

  ConfigurationSaver::saveConfigurationWithStressAndEnergy2D(
      &postOptUserData, relaxed_file_id, post_energy, post_stress, true);
  ConfigurationSaver::saveTriangleData(&postOptUserData, relaxed_file_id, domain_dims, offsets, full_mapping);
  ConfigurationSaver::saveElements(elements, vis_elements_post, relaxed_file_id);

  auto [num_dislocations_post, coordination_post] =
      DefectAnalysis::analyzeDefectsInReferenceConfig(
          &postOptUserData, relaxed_file_id, dndx, offsets, original_domain_map,
          translation_map, domain_dims_point, element_area, pbc, true);

  ConfigurationSaver::writeToVTK(postOptUserData.points, postOptUserData.elements,
                                 &postOptUserData, relaxed_file_id, true, coordination_post, 1.0);
  ConfigurationSaver::logDislocationData(1.0, num_dislocations_post);

  double post_area = ConfigurationSaver::calculateTotalArea2D(&postOptUserData);
  ConfigurationSaver::logEnergyAndStress_v2(1, 1.0, pre_energy, pre_stress,
                                            post_energy, post_stress, pre_area,
                                            post_area, enable_remeshing);

  std::cout << "Saved relaxed state as configuration_"
            << std::setw(5) << std::setfill('0') << relaxed_file_id << ".vtk ("
            << (export_full_mesh ? "full mesh" : "active elements only") << ")" << std::endl;
  std::cout << "=== VOLTERRA DISLOCATION RELAXATION COMPLETE ===\n" << std::endl;
}

