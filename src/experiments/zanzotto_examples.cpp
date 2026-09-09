#include "../../include/experiments/experiment_includes.h"

void memory(int caller_id, int nx, int ny, int restart_iteration) {
  // Restart simulation from a saved iteration

  if (nx <= 0 || ny <= 0) {
    std::cerr << "Error: nx and ny must be positive integers." << std::endl;
    exit(EXIT_FAILURE);
  }

  if (restart_iteration < 0) {
    std::cerr << "Error: restart_iteration must be non-negative." << std::endl;
    exit(EXIT_FAILURE);
  }

  std::cout << "\n" << std::string(60, '=') << std::endl;
  std::cout << "RESTARTING SIMULATION FROM ITERATION " << restart_iteration
            << std::endl;
  std::cout << std::string(60, '=') << std::endl;

  // ==================== SETUP REFERENCE GEOMETRY ====================
  writeSizesToFile(nx, ny);

  std::string lattice_type = "square";
  double h = 1.0;

  Eigen::Vector2d p1(0.0, 0.0);
  Eigen::Vector2d p2(h, 0.0);
  Eigen::Vector2d p3(0.0, h);

  Eigen::Matrix<double, 3, 2> dndx = calculateShapeDerivatives(p1, p2, p3);

  std::cout << "Reference element shape derivatives (dN/dx):" << std::endl;
  std::cout << dndx << std::endl;

  // ==================== SETUP ENERGY POTENTIAL ====================
  std::function<double(double)> potential_func = square_energy;
  std::function<double(double)> potential_func_der = square_energy_der;
  std::function<double(double)> potential_func_sder = square_energy_der;

  // ==================== LATTICE PARAMETERS ====================
  double symmetry_constantx =
      (lattice_type == "triangular") ? pow(4.0 / 3.0, 1.0 / 4.0) : 1.0;
  double optimal_lattice_parameter = symmetry_constantx * 1.0;
  double lattice_constant = optimal_lattice_parameter;

  std::cout << "Optimal lattice parameter: " << lattice_constant << std::endl;

  // ==================== LOAD SAVED CONFIGURATION ====================
  std::cout << "\nLoading configuration from iteration " << restart_iteration
            << "..." << std::endl;

  ConfigurationSaver::IterationData data =
      ConfigurationSaver::readIterationData(restart_iteration);

  // Convert loaded Eigen::Vector2d positions to Point2D objects
  std::vector<Point2D> square_points;
  square_points.reserve(data.positions.size());
  for (const auto &pos : data.positions) {
    square_points.emplace_back(pos.x(), pos.y());
  }

  // Generate reference configuration for element setup
  std::vector<Point2D> square_points_ref =
      LatticeGenerator::generate_2d_lattice(nx, ny, lattice_constant,
                                            lattice_type);

  // Extract loaded deformation gradient
  Eigen::Matrix2d F_ext = data.F_ext;
  double current_alpha = F_ext(0, 1); // Current shear strain

  std::cout << "Loaded " << square_points.size() << " points" << std::endl;
  std::cout << "F_ext matrix:\n" << F_ext << std::endl;
  std::cout << "Current shear strain α = " << current_alpha << std::endl;

  // ==================== SETUP DOMAIN ====================
  int original_domain_size = square_points.size();
  DomainInfo domain_size = compute_domain_size(square_points_ref);

  const std::array<double, 2> offsets =
      (lattice_type == "square")
          ? std::array<double, 2>{lattice_constant, lattice_constant}
          : std::array<double, 2>{lattice_constant / 2.0,
                                  (sqrt(3.0) / 2.0) * lattice_constant};

  std::cout << "PBC offsets: [" << offsets[0] << ", " << offsets[1] << "]"
            << std::endl;

  DomainDimensions domain_dims(domain_size.get_width(),
                               domain_size.get_height());
  Point2D domain_dims_point(domain_dims.size_x, domain_dims.size_y);

  bool pbc = true;

  auto [original_domain_map, translation_map] =
      MeshGenerator::create_domain_maps(original_domain_size, domain_dims,
                                        offsets);

  auto [interior_mapping, full_mapping] =
      create_dof_mapping_original(square_points, 0.5 * lattice_constant, pbc);

  std::cout << "Interior nodes: " << interior_mapping.size() << std::endl;
  std::cout << "Total nodes: " << full_mapping.size() << std::endl;

  // ==================== LOAD MESH ELEMENTS ====================
  std::cout << "\nLoading mesh elements from iteration " << restart_iteration
            << "..." << std::endl;

  auto [elements, active_elements] = ConfigurationSaver::loadElements(
      restart_iteration, dndx, full_mapping, square_points);

  std::cout << "Loaded " << elements.size() << " elements, "
            << active_elements.size() << " active" << std::endl;

  double element_area = elements[0].getReferenceArea();

  // Update element properties with loaded configuration
  for (auto &element : elements) {
    element.set_reference_mesh(square_points_ref);
    element.set_dof_mapping(full_mapping);
    element.setExternalDeformation(F_ext);
    element.calculate_deformation_gradient(square_points);
    element_area = element.getReferenceArea();
  }

  // element_area = element.getArea();
  // ref_element_area = element.getReferenceArea();

  // ==================== SETUP ENERGY CALCULATOR ====================
  Strain_Energy_LatticeCalculator calculator(1.0);
  Eigen::Matrix2d F_I = Eigen::Matrix2d::Identity();
  Eigen::Matrix2d C_I = F_I.transpose() * F_I;
  double zero = calculator.calculate_energy(C_I, potential_func, 0);

  std::cout << "Zero energy reference: " << zero << std::endl;

  // ==================== VERIFY LOADED STATE ====================
  bool plasticity = false;
  UserData checkUserData(square_points, elements, calculator, potential_func,
                         potential_func_der, zero, optimal_lattice_parameter,
                         F_ext, interior_mapping, full_mapping, active_elements,
                         plasticity);

  double loaded_energy = 0.0;
  Eigen::Matrix2d loaded_stress_tensor = Eigen::Matrix2d::Zero();
  ConfigurationSaver::calculateEnergyAndStress(&checkUserData, loaded_energy,
                                               loaded_stress_tensor, true);

  std::cout << "\nVerifying loaded state:" << std::endl;
  std::cout << "  Energy: " << loaded_energy << std::endl;
  std::cout << "  Shear stress σ_xy: " << loaded_stress_tensor(0, 1)
            << std::endl;

  // ==================== SETUP LOADING SCHEDULE FROM RESTART POINT
  // ====================
  double alpha_min = current_alpha; // Start from next increment
  double alpha_max = 0.145051;
  double step_size = 8e-7;
  if (restart_iteration % 2 != 0)
    alpha_min += step_size;

  int num_alpha_points =
      static_cast<int>((alpha_max - alpha_min) / step_size) + 1;
  std::cout << "\nResuming loading schedule:" << std::endl;
  std::cout << "  From α = " << alpha_min << " to α = " << alpha_max
            << std::endl;
  std::cout << "  Steps: " << num_alpha_points << " (step size: " << step_size
            << ")" << std::endl;

  // Generate loading sequence
  std::vector<double> alpha_values;
  alpha_values.reserve(num_alpha_points);
  for (int i = 0; i < num_alpha_points; i++) {
    alpha_values.push_back(alpha_min + i * step_size);
  }

  // ==================== INITIALIZE FILE COUNTER ====================
  static int file_counter = restart_iteration + 1;
  static int previous_file_id = restart_iteration;
  static double post_energy_previous = loaded_energy;

  std::cout << "Starting file counter at: " << file_counter << std::endl;
  std::cout << std::string(60, '=') << "\n" << std::endl;

  // ==================== RESUME SIMULATION LOOP ====================
  for (size_t i = 0; i < alpha_values.size(); i++) {
    double alpha = alpha_values[i];
    std::cout << "\n=== Processing α = " << alpha << " (restart step "
              << (restart_iteration + i + 1) << ") ===" << std::endl;

    double pre_area = 1.0;
    double post_area = 1.0;

    // ==================== APPLY DEFORMATION INCREMENT ====================
    Eigen::Matrix2d dF_ext;
    if (restart_iteration % 2 == 0 && i == 0)
      dF_ext << 1.0, 0, 0.0, 1.0;
    if (restart_iteration % 2 != 0 && i == 0)
      dF_ext << 1.0, step_size, 0.0, 1.0;
    else
      dF_ext << 1.0, step_size, 0.0, 1.0;

    F_ext << 1.0, alpha, 0.0, 1.0;

    // Apply incremental deformation
    for (size_t j = 0; j < square_points.size(); j++) {
      square_points[j].coord = dF_ext * square_points[j].coord;
    }

    // ==================== CREATE USER DATA ====================
    UserData userData(square_points, elements, calculator, potential_func,
                      potential_func_der, zero, optimal_lattice_parameter,
                      F_ext, interior_mapping, full_mapping, active_elements,
                      plasticity);

    // ==================== PREPARE OPTIMIZATION ====================
    alglib::real_1d_array x;
    int n_vars = interior_mapping.size();
    x.setlength(2 * n_vars);
    map_points_to_solver_array(x, square_points, interior_mapping, n_vars);

    // Calculate pre-optimization energy and stress
    double pre_energy = 0.0;
    double pre_stress = 0.0;
    Eigen::Matrix2d stress_tensor = Eigen::Matrix2d::Zero();

    ConfigurationSaver::calculateEnergyAndStress(&userData, pre_energy,
                                                 stress_tensor, true);
    pre_stress = stress_tensor(0, 1);
    pre_area = ConfigurationSaver::calculateTotalArea2D(&userData);

    std::cout << "Pre-optimization - Energy: " << pre_energy
              << ", Stress: " << pre_stress << std::endl;

    // Store original positions
    alglib::real_1d_array original_x;
    original_x.setlength(x.length());
    for (int j = 0; j < x.length(); j++) {
      original_x[j] = x[j];
    }

    // ==================== SAVE BEFORE OPTIMIZATION ====================
    int file_id = caller_id + file_counter;
    double saving_value = alpha;

    UserData preOptUserData(square_points, elements, calculator, potential_func,
                            potential_func_der, zero, optimal_lattice_parameter,
                            F_ext, interior_mapping, full_mapping,
                            active_elements, plasticity);

    ConfigurationSaver::saveConfigurationWithStressAndEnergy2D(
        &preOptUserData, file_id, pre_energy, pre_stress, true);
    pre_area = ConfigurationSaver::calculateTotalArea2D(&preOptUserData);

    ConfigurationSaver::saveTriangleData(&preOptUserData, file_id, domain_dims,
                                         offsets, full_mapping);
    ConfigurationSaver::saveElements(elements, active_elements, file_id);

    auto [num_dislocations_pre, coordination_pre] =
        DefectAnalysis::analyzeDefectsInReferenceConfig(
            &preOptUserData, file_id, dndx, offsets, original_domain_map,
            translation_map, domain_dims_point, element_area, pbc, true);

    ConfigurationSaver::writeToVTK(
        preOptUserData.points, preOptUserData.elements, &preOptUserData,
        file_id, true, coordination_pre, saving_value);

    ConfigurationSaver::logDislocationData(alpha, num_dislocations_pre);

    std::cout << "Saved PRE-optimization config " << file_id
              << " at load=" << saving_value << std::endl;

    NeighborAnalyzer analyzer(NeighborAnalyzer::SearchType::K_NEAREST);
    std::vector<Point2D> points_before_copy = square_points;
    analyzer.setKNearestSearch(4); // Just stores k=6, doesn't need points yet
    analyzer.setDebugMode(false);  // ← Active les prints

    // USAGE: Give it points when you want to build the neighbor list
    auto neighbors_before = analyzer.buildNeighbors(points_before_copy);

    // ==================== RUN OPTIMIZATION ====================
    auto wall_start = std::chrono::high_resolution_clock::now();
    clock_t cpu_start = clock();

    userData.third_condition_flag = false;
    LBFGSOptimizer optimizer(12, 0, pow(10., -13.), 0, 0);
    optimizer.optimize(x, minimize_energy_with_triangles, &userData);

    auto wall_end = std::chrono::high_resolution_clock::now();
    clock_t cpu_end = clock();

    double wall_time =
        std::chrono::duration<double>(wall_end - wall_start).count();
    double cpu_time = (double)(cpu_end - cpu_start) / CLOCKS_PER_SEC;
    std::cout << "Optimization wall-clock time: " << wall_time << " seconds\n";
    std::cout << "Optimization CPU time: " << cpu_time << " seconds\n";
    std::cout << "Optimization Ratio: " << cpu_time / wall_time << "\n";

    map_solver_array_to_points(x, square_points, interior_mapping, n_vars);

    auto neighbors_after = analyzer.buildNeighbors(square_points); // ← AND HERE

    auto change_info = NeighborAnalyzer::compareNeighborsWithTolerance(
        points_before_copy, // Points AVANT
        square_points,      // Points APRÈS
        neighbors_before,   // Voisins AVANT
        neighbors_after,    // Voisins APRÈS
        0.1                 // 10% de tolérance (optionnel, défaut = 0.1)
    );

    if (change_info.has_changed) {
      std::cout << "✗ Neighbor connectivity CHANGED" << std::endl;
      std::cout << "  → Nodes affected: " << change_info.total_nodes_changed
                << std::endl;
      std::cout << "  → Connections added: "
                << change_info.total_connections_added << std::endl;
      std::cout << "  → Connections removed: "
                << change_info.total_connections_removed << std::endl;
    } else {
      std::cout << "✓ Neighbor connectivity UNCHANGED" << std::endl;
    }

    // ==================== POST-OPTIMIZATION ENERGY ====================
    double post_energy = 0.0;
    double post_stress = 0.0;
    stress_tensor.setZero();
    ConfigurationSaver::calculateEnergyAndStress(&userData, post_energy,
                                                 stress_tensor, true);
    post_stress = stress_tensor(0, 1);
    post_area = ConfigurationSaver::calculateTotalArea2D(&userData);

    std::cout << "Post-optimization - Energy: " << post_energy
              << ", Stress: " << post_stress << std::endl;
    std::cout << "Energy change: " << (post_energy - pre_energy)
              << ", Stress change: " << (post_stress - pre_stress) << std::endl;

    // ==================== REMESHING (if needed) ====================
    // ChangeMeasures result = computeChangeMeasures(
    //     x, original_x, lattice_constant, elements, &userData, square_points,
    //     true, &F_ext);
    // bool shouldRemesh = result.has_distorted_triangles;

    bool shouldRemesh =
        change_info
            .has_changed; // disable remeshing for continous Zanzotto test

    if (shouldRemesh) {
      std::cout << "REMESHING STARTS" << std::endl;

      std::vector<int> contact_atoms;
      std::vector<int> boundary_fixed_nodes;
      int max_iterations = 1000;
      int hasChanges = 0;

      auto [post_energy_re, stress_tensor_re, iterations] =
          perform_remeshing_loop_reduction(
              x, &userData, contact_atoms, boundary_fixed_nodes, F_ext, dndx,
              offsets, original_domain_map, translation_map, domain_dims_point,
              hasChanges, max_iterations, element_area, pbc, true);

      for (auto &element : elements) {
        element.set_dof_mapping(full_mapping);
      }

      post_energy = post_energy_re;
      post_stress = stress_tensor_re(0, 1);

      std::cout << "Post-remeshing - Energy: " << post_energy
                << ", Stress: " << post_stress << std::endl;
    }

    // ==================== CHECK FOR STRESS DROP ====================
    bool stress_drop_detected = shouldRemesh;

    UserData postOptUserData(square_points, elements, calculator,
                             potential_func, potential_func_der, zero,
                             optimal_lattice_parameter, F_ext, interior_mapping,
                             full_mapping, active_elements, plasticity);
    post_area = ConfigurationSaver::calculateTotalArea2D(&postOptUserData);

    if (stress_drop_detected || i >= 0) {
      std::cout << "=== STRESS DROP DETECTED ===" << std::endl;
      std::cout << "Energy dropped from " << post_energy_previous << " to "
                << post_energy << std::endl;
      std::cout << "PRE-avalanche LOCKED as file " << file_id
                << " at load=" << saving_value << std::endl;

      // Save POST-avalanche state
      file_counter++;
      int post_file_id = caller_id + file_counter;

      ConfigurationSaver::saveTriangleData(&postOptUserData, post_file_id,
                                           domain_dims, offsets, full_mapping);
      ConfigurationSaver::saveElements(elements, active_elements, post_file_id);

      auto [num_dislocations_post, coordination_post] =
          DefectAnalysis::analyzeDefectsInReferenceConfig(
              &postOptUserData, post_file_id, dndx, offsets,
              original_domain_map, translation_map, domain_dims_point,
              element_area, pbc, true);

      ConfigurationSaver::writeToVTK(
          postOptUserData.points, postOptUserData.elements, &postOptUserData,
          post_file_id, true, coordination_post, saving_value);
      ConfigurationSaver::logDislocationData(alpha, num_dislocations_post);

      std::cout << "POST-avalanche saved as file " << post_file_id
                << " at load=" << saving_value << std::endl;

      file_counter++;
      previous_file_id = -1;

    } else {
      // No stress drop - delete previous file if it exists
      if (previous_file_id >= 0 || i >= 0) {
        std::cout << "Deleting previous file " << previous_file_id
                  << " (no avalanche)" << std::endl;

        std::stringstream vtk_file;
        vtk_file << "vtk_output/configuration_" << std::setw(5)
                 << std::setfill('0') << previous_file_id << ".vtk";
        std::filesystem::remove(vtk_file.str());
      }

      previous_file_id = file_id;
    }

    // ==================== LOG DATA ====================
    ConfigurationSaver::logEnergyAndStress_v2(
        restart_iteration + i, alpha, pre_energy, pre_stress, post_energy,
        post_stress, pre_area, post_area, shouldRemesh);

    post_energy_previous = post_energy;

    std::cout << "Iteration " << (restart_iteration + i + 1)
              << " completed successfully" << std::endl;
  }

  std::cout << "\n" << std::string(60, '=') << std::endl;
  std::cout << "RESTART SIMULATION COMPLETED" << std::endl;
  std::cout << std::string(60, '=') << std::endl;
}

// NOTE: example_1_conti_zanzotto now lives in
// src/experiments/dislocation_study.cpp (stripped, dislocation-only variant).
// This original loading-driven version is renamed to avoid a duplicate symbol.
void example_1_conti_zanzotto_loading(
    int caller_id, int nx, int ny,
    double alpha_min,
    double alpha_max,
    double step_size,
    double triangulation_perturbation,
    unsigned int seed) {

  //     auto compute_even_ny = [](int nx) {
  //     int ny = std::round(2.0 * nx / std::sqrt(3));
  //     return (ny % 2 == 0) ? ny : ny + 1; // Ensure ny is even
  // };

  // ny = compute_even_ny(nx);

  // Parameters for lattice
  if (nx <= 0 || ny <= 0) {
    std::cerr << "Error: nx and ny must be positive integers." << std::endl;
    exit(EXIT_FAILURE);
  }

  // ==================== SETUP REFERENCE GEOMETRY ====================
  // Save domain dimensions
  writeSizesToFile(nx, ny);

  // Define lattice type and reference element
  std::string lattice_type = "square"; // Options: "square" or "triangular"
  double h = 1.0;                      // Reference element size

  // Reference triangle vertices (for shape function derivatives)
  Eigen::Vector2d p1(0.0, 0.0);
  Eigen::Vector2d p2(h, 0.0);
  Eigen::Vector2d p3(0.0, h);

  // Calculate shape function derivatives for reference element
  Eigen::Matrix<double, 3, 2> dndx = calculateShapeDerivatives(p1, p2, p3);

  std::cout << "Reference element shape derivatives (dN/dx):" << std::endl;
  std::cout << dndx << std::endl;

  // ==================== SETUP ENERGY POTENTIAL ====================
  // Define DUMMY ATOMISTIC energy functions (currently using square potential)
  std::function<double(double)> potential_func = square_energy;
  std::function<double(double)> potential_func_der = square_energy_der;
  std::function<double(double)> potential_func_sder = square_energy_der;

  // ==================== DETERMINE OPTIMAL LATTICE PARAMETER
  // ====================
  std::cout << "STEP 1: Finding optimal lattice parameter..." << std::endl;

  // Lattice symmetry correction factor
  double symmetry_constantx =
      (lattice_type == "triangular") ? pow(4.0 / 3.0, 1.0 / 4.0) : 1.0;
  std::cout << "Symmetry constant for " << lattice_type
            << " lattice: " << symmetry_constantx << std::endl;

  // Set optimal lattice parameter
  double optimal_lattice_parameter = symmetry_constantx * 1.0;
  double lattice_constant = optimal_lattice_parameter;
  std::cout << "Optimal lattice parameter: " << lattice_constant << std::endl;

  // ==================== GENERATE INITIAL LATTICE ====================
  // Generate current and reference lattice configurations
  std::vector<Point2D> square_points = LatticeGenerator::generate_2d_lattice(
      nx, ny, lattice_constant, lattice_type);

  std::vector<Point2D> square_points_ref =
      LatticeGenerator::generate_2d_lattice(nx, ny, lattice_constant,
                                            lattice_type);

  // auto [square_points_ref, gap_offsets]
  // =LatticeGenerator::generate_periodic_rotated_lattice_v2( 1, 1,
  // nx,ny,lattice_constant, lattice_type) ;

  // auto [square_points, gap_offsets2]
  // =LatticeGenerator::generate_periodic_rotated_lattice_v2( 1, 1,
  // nx,ny,lattice_constant, lattice_type) ;

  // std::cout << "User-defined gap offset vector: ("
  //           << gap_offsets.coord.x() << ", "
  //           << gap_offsets.coord.y() << ")" << std::endl;

  int original_domain_size = square_points.size();
  std::cout << "Generated lattice with " << original_domain_size << " points"
            << std::endl;

  // Calculate domain properties
  DomainInfo domain_size = compute_domain_size(square_points);

  // Set periodic boundary offsets based on lattice type

  const std::array<double, 2> offsets =
      (lattice_type == "square")
          ? std::array<double, 2>{lattice_constant, lattice_constant}
          : std::array<double, 2>{lattice_constant / 2.0,
                                  (sqrt(3.0) / 2.0) * lattice_constant};

  // const std::array<double, 2> offsets = (lattice_type == "square") ?
  //     std::array<double, 2>{sqrt(2)*lattice_constant/2.,
  //     sqrt(2)*lattice_constant/2.} : std::array<double, 2>{lattice_constant
  //     / 2.0, (sqrt(3.0) / 2.0) * lattice_constant};

  // const std::array<double, 2> offsets =
  // {gap_offsets.coord.x(),gap_offsets.coord.y()};
  std::cout << "PBC offsets: [" << offsets[0] << ", " << offsets[1] << "]"
            << std::endl;

  DomainDimensions domain_dims(domain_size.get_width(),
                               domain_size.get_height());
  std::cout << "domain_size.get_width(): " << domain_size.get_width()
            << std::endl;
  std::cout << "domain_size.get_height(): " << domain_size.get_height()
            << std::endl;

  // Setup triangulation variables
  bool pbc = true;

  // Create domain maps
  auto [original_domain_map, translation_map] =
      MeshGenerator::create_domain_maps(original_domain_size, domain_dims,
                                        offsets);

  // Boundary conditions
  auto [interior_mapping, full_mapping] =
      create_dof_mapping_original(square_points, 0.5 * lattice_constant, pbc);
  // auto [interior_mapping, full_mapping] =
  // create_dof_mapping_original(square_points, lattice_constant, pbc);

  // auto [interior_mapping, full_mapping] =
  // create_dof_mapping_with_radius(square_points, 90, pbc);

  std::cout << "interior_mapping.size(): " << interior_mapping.size()
            << std::endl;
  std::cout << "full_mapping.size(): " << full_mapping.size() << std::endl;

  Point2D domain_dims_point(domain_dims.size_x, domain_dims.size_y);

  // // Call the function
  // Eigen::Vector2d burgers_vector(lattice_constant, 0.0); // Example: unit
  // Burgers vector in x-direction double core_radius = 6.0; // Example core
  // radius value double poisson_ratio = 0.3; // Typical value for many
  // materials size_t middle_atom_index = findMiddleAtom(square_points, true);
  // // true enables verbose output

  // std::vector<Point2D> dislocated_points = createSingleDislocation(
  //     square_points_ref,
  //     burgers_vector,
  //     middle_atom_index,
  //     core_radius,
  //     poisson_ratio
  // );

  // auto dipole_points = createDislocationDipole(
  //     square_points_ref,
  //     Eigen::Vector2d(lattice_constant, 0.0),
  //     Eigen::Vector2d(square_points[middle_atom_index].coord.x(),
  //     square_points[middle_atom_index].coord.y()), 180.0,  // separation
  //     distance 1.0,  // core radius 0.33,  // Poisson's ratio
  //     Eigen::Vector2d(1.0, 0.0) // dipole direction (45 degrees)
  // );

  // square_points = dipole_points;

  // 1. Create the AdaptiveMesher instance
  AdaptiveMesher mesher(domain_dims_point, offsets, original_domain_map,
                        translation_map, full_mapping,
                        1e-6, // Tolerance
                        pbc   // Use periodic copies
  );
  mesher.setUsePeriodicCopies(
      pbc); // Switch to using original domain only not necessary

  // Apply initial mesh orientation perturbation if specified.
  // For a square lattice, Delaunay triangulation has degenerate 4-node circumcircles,
  // making the diagonal selection ambiguous. A slight shear perturbation (e.g. -1e-7),
  // as done in the shifting study (shift_vertical_horizontal), breaks the tie to
  // select the opposite diagonal orientation across the initial mesh.
  if (triangulation_perturbation != 0.0) {
    mesher.setTriangulationPerturbation(triangulation_perturbation);
    std::cout << "Applied initial mesh orientation perturbation: "
              << triangulation_perturbation << std::endl;
  }
  alglib::real_1d_array free_dofs;
  int n_free_nodes = interior_mapping.size();
  free_dofs.setlength(2 *
                      interior_mapping.size()); // [u0, u1, ..., v0, v1, ...]
  map_points_to_solver_array(free_dofs, square_points, interior_mapping,
                             n_free_nodes);

  alglib::real_1d_array original_x_remesh =
      mesher.saveOriginalPositions(free_dofs);

  auto [elements, active_elements] = mesher.createMesh(
      square_points, free_dofs, Eigen::Matrix2d::Identity(), &dndx);
  double element_area = elements[0].getReferenceArea();

  for (auto &element : elements) {
    // element.set_reference_mesh(square_points);
    element.set_dof_mapping(
        full_mapping); // or interior_mapping depending on needs
    // double jac = element.calculate_shape_derivatives(x);  // current
    // positions
    const Eigen::Matrix<double, 3, 2> &dndx = element.getDNdX();
    // std::cout<< "dndx: " << dndx << std::endl;
  }

  // square_points = dipole_points;

  std::cout << "Created " << elements.size() << " element triangles"
            << std::endl;

  // Setup energy calculation
  Strain_Energy_LatticeCalculator calculator(1.0);

  Eigen::Matrix2d F_I = Eigen::Matrix2d::Identity();
  // F_I *= symmetry_constantx;
  Eigen::Matrix2d C_I = F_I.transpose() * F_I; // C = F^T * F

  double zero = calculator.calculate_energy(C_I, potential_func, 0);

  // Eigen::Matrix2d C_I = Eigen::Matrix2d::Identity();
  // double zero = calculator.calculate_energy(C_I, potential_func, 0);
  std::cout << "debugging simple shear test" << std::endl;
  std::cout << "zero energy value: " << zero << std::endl;
  // debug_deformation_tests();

  // ==================== SETUP LOADING SCHEDULE ====================
  // Calculate number of loading steps supporting both positive and negative loading paths
  int num_alpha_points =
      static_cast<int>(std::abs(alpha_max - alpha_min) / std::abs(step_size)) + 1;
  std::cout << "Loading schedule: " << num_alpha_points << " steps from "
            << alpha_min << " to " << alpha_max << " (step size: " << step_size
            << ")" << std::endl;

  // Generate loading sequence
  std::vector<double> alpha_values;
  alpha_values.reserve(num_alpha_points);
  for (int i = 0; i < num_alpha_points; i++) {
    alpha_values.push_back(alpha_min + i * step_size);
  }

  // Process each alpha value
  for (size_t i = 0; i < alpha_values.size(); i++) {
    double alpha = alpha_values[i];
    std::cout << "\n=== Processing alpha = " << alpha << " ===" << std::endl;
    double pre_area = 1.0;
    double post_area = 1.0;

    // ==================== SETUP DEFORMATION ====================
    Eigen::Matrix2d F_ext;
    F_ext << 1.0, alpha, 0.0, 1.0;

    Eigen::Matrix2d dF_ext;
    dF_ext << 1.0, step_size, 0.0, 1.0;

    // Apply initial noise (only for first iteration)
    if (i == 0) {
      // Use deterministic seed so both positive and negative runs share the exact same initial noise
      std::mt19937 gen(seed);
      double noise_level = 0.04;
      std::normal_distribution<double> noise_dist(0.0, noise_level);
      std::cout << "Seeded initial noise generator with seed: " << seed << std::endl;

      for (size_t j = 0; j < square_points.size(); j++) {
        Eigen::Vector2d noise(noise_dist(gen), noise_dist(gen));
        square_points[j].coord = F_ext * square_points[j].coord + noise;
      }
    } else {
      for (size_t j = 0; j < square_points.size(); j++) {
        square_points[j].coord = dF_ext * square_points[j].coord;
      }
    }

    // ==================== CREATE USER DATA ====================
    bool plasticity = false;
    UserData userData(square_points, elements, calculator, potential_func,
                      potential_func_der, zero, optimal_lattice_parameter,
                      F_ext, interior_mapping, full_mapping, active_elements,
                      plasticity);

    // ==================== PREPARE OPTIMIZATION ====================
    alglib::real_1d_array x;
    int n_vars = interior_mapping.size();
    x.setlength(2 * n_vars);
    map_points_to_solver_array(x, square_points, interior_mapping, n_vars);

    // Calculate pre-optimization energy and stress
    double pre_energy = 0.0;
    double pre_stress = 0.0;
    Eigen::Matrix2d stress_tensor = Eigen::Matrix2d::Zero();

    ConfigurationSaver::calculateEnergyAndStress(&userData, pre_energy,
                                                 stress_tensor, true);
    pre_stress = stress_tensor(0, 1);
    pre_area = ConfigurationSaver::calculateTotalArea2D(&userData);

    std::cout << "Pre-optimization - Energy: " << pre_energy
              << ", Stress: " << pre_stress << std::endl;

    // Store original positions
    alglib::real_1d_array original_x;
    original_x.setlength(x.length());
    for (int j = 0; j < x.length(); j++) {
      original_x[j] = x[j];
    }

    // ==================== SAVE BEFORE OPTIMIZATION ====================
    static int file_counter = 0;
    static int previous_file_id = -1;
    static double post_energy_previous = 0.0;

    int file_id = caller_id + file_counter;
    double saving_value = alpha_values[i];

    // Save pre-optimization state (potential pre-avalanche)
    UserData preOptUserData(square_points, elements, calculator, potential_func,
                            potential_func_der, zero, optimal_lattice_parameter,
                            F_ext, interior_mapping, full_mapping,
                            active_elements, plasticity);

    ConfigurationSaver::saveConfigurationWithStressAndEnergy2D(
        &preOptUserData, file_id, pre_energy, pre_stress, true);
    pre_area = ConfigurationSaver::calculateTotalArea2D(&preOptUserData);

    ConfigurationSaver::saveTriangleData(&preOptUserData, file_id, domain_dims,
                                         offsets, full_mapping);
    ConfigurationSaver::saveElements(elements, active_elements, file_id);

    std::cout << "Pre-optimization2 - Energy: " << pre_energy
              << ", Stress: " << pre_stress << std::endl;

    auto [num_dislocations_pre, coordination_pre] =
        DefectAnalysis::analyzeDefectsInReferenceConfig(
            &preOptUserData, file_id, dndx, offsets, original_domain_map,
            translation_map, domain_dims_point, element_area, pbc, true);

    ConfigurationSaver::writeToVTK(
        preOptUserData.points, preOptUserData.elements, &preOptUserData,
        file_id, true, coordination_pre, saving_value);

    ConfigurationSaver::logDislocationData(alpha, num_dislocations_pre);

    std::cout << "Saved PRE-optimization config " << file_id
              << " at load=" << saving_value << std::endl;

    std::cout << "================================\n" << std::endl;

    // ==================== RUN OPTIMIZATION ====================
    NeighborAnalyzer analyzer(NeighborAnalyzer::SearchType::K_NEAREST);
    std::vector<Point2D> points_before_copy = square_points;
    analyzer.setKNearestSearch(4);
    analyzer.setDebugMode(false);

    auto neighbors_before = analyzer.buildNeighbors(points_before_copy);

    auto wall_start = std::chrono::high_resolution_clock::now();
    clock_t cpu_start = clock();

    userData.third_condition_flag = false;
    LBFGSOptimizer optimizer(13, 0.00001, 0, 0, 0);
    optimizer.optimize(x, minimize_energy_with_triangles, &userData);

    auto wall_end = std::chrono::high_resolution_clock::now();
    clock_t cpu_end = clock();

    double wall_time =
        std::chrono::duration<double>(wall_end - wall_start).count();
    double cpu_time = (double)(cpu_end - cpu_start) / CLOCKS_PER_SEC;
    std::cout << "Optimization wall-clock time: " << wall_time << " seconds\n";
    std::cout << "Optimization CPU time: " << cpu_time << " seconds\n";
    std::cout << "Optimization Ratio: " << cpu_time / wall_time << "\n";

    map_solver_array_to_points(x, square_points, interior_mapping, n_vars);

    auto neighbors_after = analyzer.buildNeighbors(square_points);

    auto change_info = NeighborAnalyzer::compareNeighborsWithTolerance(
        points_before_copy, square_points, neighbors_before, neighbors_after,
        0.135);

    if (change_info.has_changed) {
      std::cout << "✗ Neighbor connectivity CHANGED" << std::endl;
      std::cout << "  → Nodes affected: " << change_info.total_nodes_changed
                << std::endl;
      std::cout << "  → Connections added: "
                << change_info.total_connections_added << std::endl;
      std::cout << "  → Connections removed: "
                << change_info.total_connections_removed << std::endl;
    } else {
      std::cout << "✓ Neighbor connectivity UNCHANGED" << std::endl;
    }

    int hasChanges = 0;

    // ==================== POST-OPTIMIZATION ENERGY ====================
    double post_energy = 0.0;
    double post_stress = 0.0;
    stress_tensor.setZero();
    ConfigurationSaver::calculateEnergyAndStress(&userData, post_energy,
                                                 stress_tensor, true);
    post_stress = stress_tensor(0, 1);
    post_area = ConfigurationSaver::calculateTotalArea2D(&userData);

    std::cout << "Post-optimization - Energy: " << post_energy
              << ", Stress: " << post_stress << std::endl;
    std::cout << "Energy change: " << (post_energy - pre_energy)
              << ", Stress change: " << (post_stress - pre_stress) << std::endl;

    // ==================== REMESHING DECISION ====================
    bool shouldRemesh2 =
        change_info.has_changed && post_energy < post_energy_previous;

    bool shouldRemesh = post_energy < post_energy_previous || i == 0;

    // Debug output
    std::cout << "=== REMESH DECISION ===" << std::endl;
    std::cout << "  i = " << i << std::endl;
    std::cout << "  pre_energy = " << std::scientific << pre_energy
              << std::endl;

    std::cout << "  post_energy = " << std::scientific << post_energy
              << std::endl;
    std::cout << "  post_energy_previous = " << std::scientific
              << post_energy_previous << std::endl;
    std::cout << "  difference = " << (post_energy - post_energy_previous)
              << std::endl;
    std::cout << "  post_energy < post_energy_previous? "
              << (post_energy < post_energy_previous) << std::endl;
    std::cout << "  shouldRemesh = " << shouldRemesh << std::endl;

    if (shouldRemesh) {
      std::cout << "REMESHING STARTS" << std::endl;

      std::vector<int> contact_atoms;
      std::vector<int> boundary_fixed_nodes;
      int max_iterations = 1000;

      auto [post_energy_re, stress_tensor_re, iterations] =
          perform_remeshing_loop_reduction(
              x, &userData, contact_atoms, boundary_fixed_nodes, F_ext, dndx,
              offsets, original_domain_map, translation_map, domain_dims_point,
              hasChanges, max_iterations, element_area, pbc, true);

      post_energy = post_energy_re;
      post_stress = stress_tensor_re(0, 1);

      // ============================================
      // SYNC LOCAL VARIABLES WITH USERDATA
      // (Not strictly necessary since UserData uses references,
      //  but explicit for clarity)
      // ============================================
      square_points = userData.points;
      elements = userData.elements;
      active_elements = userData.active_elements;

      // Update element_area after remeshing
      if (!elements.empty()) {
        element_area = elements[0].getReferenceArea();
      }

      if (hasChanges) {
        std::cout << "✓ Remeshing accepted - energy decreased" << std::endl;
      } else {
        std::cout << "⚠️ Remeshing rejected - original mesh kept" << std::endl;
      }

      std::cout << "Final energy: " << post_energy
                << ", stress: " << post_stress << std::endl;

      // Update DOF mapping
      for (auto &element : elements) {
        element.set_dof_mapping(full_mapping);
      }
    }

    // ==================== CHECK FOR STRESS DROP ====================
    bool stress_drop_detected = shouldRemesh2;

    UserData postOptUserData(square_points, elements, calculator,
                             potential_func, potential_func_der, zero,
                             optimal_lattice_parameter, F_ext, interior_mapping,
                             full_mapping, active_elements, plasticity);
    post_area = ConfigurationSaver::calculateTotalArea2D(&postOptUserData);

    std::cout << "Energy dropped from " << post_energy_previous << " to "
              << post_energy << std::endl;

    if (stress_drop_detected || i == 0) {
      std::cout << "=== STRESS DROP DETECTED ===" << std::endl;
      std::cout << "Energy dropped from " << post_energy_previous << " to "
                << post_energy << std::endl;
      std::cout << "PRE-avalanche LOCKED as file " << file_id
                << " at load=" << saving_value << std::endl;

      // Save POST-avalanche state
      file_counter++;
      int post_file_id = caller_id + file_counter;

      ConfigurationSaver::saveTriangleData(&postOptUserData, post_file_id,
                                           domain_dims, offsets, full_mapping);
      ConfigurationSaver::saveElements(elements, active_elements, post_file_id);

      auto [num_dislocations_post, coordination_post] =
          DefectAnalysis::analyzeDefectsInReferenceConfig(
              &postOptUserData, post_file_id, dndx, offsets,
              original_domain_map, translation_map, domain_dims_point,
              element_area, pbc, true);

      ConfigurationSaver::writeToVTK(
          postOptUserData.points, postOptUserData.elements, &postOptUserData,
          post_file_id, true, coordination_post, saving_value);
      ConfigurationSaver::logDislocationData(alpha, num_dislocations_post);

      std::cout << "POST-avalanche saved as file " << post_file_id
                << " at load=" << saving_value << std::endl;

      file_counter++;
      previous_file_id = -1;

    } else {
      // No stress drop - delete previous file if it exists
      if (previous_file_id >= 0 || i == 0) {
        std::cout << "Deleting previous file " << previous_file_id
                  << " (no avalanche)" << std::endl;

        std::stringstream vtk_file;
        vtk_file << "vtk_output/configuration_" << std::setw(5)
                 << std::setfill('0') << previous_file_id << ".vtk";
        std::filesystem::remove(vtk_file.str());
      }

      previous_file_id = file_id;
    }

    // ==================== LOG DATA ====================
    ConfigurationSaver::logEnergyAndStress_v2(
        i, alpha, pre_energy, pre_stress, post_energy, post_stress, pre_area,
        post_area, shouldRemesh);

    post_energy_previous = post_energy;

    std::cout << "Iteration " << i
              << " ended: setting post_energy_previous = " << post_energy
              << std::endl;
  }
}

void example_1_conti_zanzotto_negative_loading(int caller_id, int nx, int ny, unsigned int seed) {
  // Negative continuous shear loading:
  // - Starts at load alpha = -0.14 and increments negatively with step_size = -3e-5 down to -0.85
  // - Flips initial mesh orientation by applying a -1e-7 shear perturbation to the Delaunay mesher,
  //   matching the orientation change technique established in the shifting experiments.
  // - Uses deterministic seed for identical initial noise generation.
  example_1_conti_zanzotto_loading(caller_id, nx, ny,
                                   /*alpha_min=*/-0.14,
                                   /*alpha_max=*/-0.85,
                                   /*step_size=*/-3e-5,
                                   /*triangulation_perturbation=*/-1e-7,
                                   /*seed=*/seed);
}
