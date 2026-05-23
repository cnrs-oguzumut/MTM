#include "../../include/experiments/experiment_includes.h"

void example_1_shifting(int caller_id, int nx, int ny) {
  // Parameters for lattice
  if (nx <= 0 || ny <= 0) {
    std::cerr << "Error: nx and ny must be positive integers." << std::endl;
    exit(EXIT_FAILURE);
  }
  writeSizesToFile(nx, ny);
  std::string lattice_type = "triangular"; // "square" or "triangular"
  double symmetry_constantx = 1.;
  if (lattice_type == "triangular")
    symmetry_constantx = pow(4. / 3., 1. / 4.);
  double h = 1.;
  Eigen::Vector2d p1(0, 0);
  Eigen::Vector2d p2(h, 0);
  Eigen::Vector2d p3(0, h);
  Eigen::Matrix<double, 3, 2> dndx = calculateShapeDerivatives(p1, p2, p3);
  std::cout << dndx << std::endl;

  // Energy functions
  std::function<double(double)> potential_func = square_energy;
  std::function<double(double)> potential_func_der = square_energy_der;

  std::cout << "STEP 1: Finding optimal lattice parameter...\n";
  // double optimal_lattice_parameter = 1.;
  // double lattice_constant = optimal_lattice_parameter;

  double optimal_lattice_parameter = symmetry_constantx * 1.0;
  double lattice_constant = optimal_lattice_parameter;

  // Generate initial lattice
  std::vector<Point2D> square_points = LatticeGenerator::generate_2d_lattice(
      nx, ny, lattice_constant, lattice_type);
  std::vector<Point2D> square_points_ref =
      LatticeGenerator::generate_2d_lattice(nx, ny, lattice_constant,
                                            lattice_type);

  int original_domain_size = square_points.size();
  DomainInfo domain_size = compute_domain_size(square_points);

  const std::array<double, 2> offsets = {lattice_constant, lattice_constant};
  // const std::array<double, 2> offsets = {lattice_constant,
  // sqrt(3.)/2*lattice_constant};

  std::cout << "offsets: " << offsets[0] << " " << offsets[1] << std::endl;
  DomainDimensions domain_dims(domain_size.get_width(),
                               domain_size.get_height());
  std::cout << "domain_size.get_width(): " << domain_size.get_width()
            << std::endl;
  std::cout << "domain_size.get_height(): " << domain_size.get_height()
            << std::endl;

  // Setup triangulation variables
  bool pbc = false;

  // Create domain maps
  auto [original_domain_map, translation_map] =
      MeshGenerator::create_domain_maps(original_domain_size, domain_dims,
                                        offsets);

  // Boundary conditions
  auto [interior_mapping, full_mapping] =
      create_dof_mapping_original(square_points, 0.5 * lattice_constant, pbc);
  // auto [interior_mapping, full_mapping] =
  // create_dof_mapping_with_radius(square_points, 90, pbc);

  std::cout << "interior_mapping.size(): " << interior_mapping.size()
            << std::endl;
  std::cout << "full_mapping.size(): " << full_mapping.size() << std::endl;

  Point2D domain_dims_point(domain_dims.size_x, domain_dims.size_y);

  // 1. Create the AdaptiveMesher instance
  AdaptiveMesher mesher(domain_dims_point, offsets, original_domain_map,
                        translation_map, full_mapping,
                        1e-6, // Tolerance
                        pbc   // Use periodic copies
  );
  mesher.setUsePeriodicCopies(pbc); // Switch to using original domain only
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

  Eigen::Matrix2d F_ext1;
  F_ext1 << 1.0, 0.0, 0.0, 1.0;

  for (size_t elem_idx : active_elements) {
    if (elem_idx >= elements.size()) {
      std::cerr << "Warning: Element index " << elem_idx << " out of range."
                << std::endl;
      continue;
    }

    auto &element = elements[elem_idx];

    // Skip if the element isn't initialized
    if (!element.isInitialized()) {
      continue;
    }

    // Set external deformation and recalculate deformation gradient
    element.setExternalDeformation(F_ext1);
    element.calculate_deformation_gradient(free_dofs);

    // Get the deformation gradient
    const Eigen::Matrix2d &F = element.getDeformationGradient();
    std::cout << "Element " << elem_idx << " Deformation Gradient F:\n"
              << F << std::endl;

    int i1 = element.getNodeIndex(0);
    int i2 = element.getNodeIndex(1);
    int i3 = element.getNodeIndex(2);
    std::cout << "Node indices: " << i1 << ", " << i2 << ", " << i3
              << std::endl;
    std::cout << "Node positions:\n";
    std::cout << "  Node " << i1 << ": (" << square_points[i1].coord.x() << ", "
              << square_points[i1].coord.y() << ")\n";
    std::cout << "  Node " << i2 << ": (" << square_points[i2].coord.x() << ", "
              << square_points[i2].coord.y() << ")\n";
    std::cout << "  Node " << i3 << ": (" << square_points[i3].coord.x() << ", "
              << square_points[i3].coord.y() << ")\n";
  }

  // exit(0);

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

  // Set up alpha values for deformation steps
  double alpha_min = 0;
  double alpha_max = +2.0;
  double step_size = 5e-2;
  int num_alpha_points =
      static_cast<int>((alpha_max - alpha_min) / step_size) + 1;
  std::cout << "num_alpha_points: " << num_alpha_points << std::endl;
  // num_alpha_points = 1;
  //  Generate evenly spaced alpha values

  std::vector<double> alpha_values;
  alpha_values.reserve(num_alpha_points);
  for (int i = 0; i < num_alpha_points; i++) {
    alpha_values.push_back(alpha_min + i * step_size);
  }

  // First, find the maximum y-coordinate DELANUAY triangulation
  double max_y = std::numeric_limits<double>::lowest();
  for (size_t i1 = 0; i1 < square_points.size(); i1++) {
    max_y = std::max(max_y, square_points[i1].coord.y());
  }

  // Define the threshold as half of the maximum y-coordinate
  double threshold = max_y / 2.0;

  for (size_t i1 = 0; i1 < square_points.size(); i1++) {
    // Only deform points in the upper half of the crystal
    if (square_points[i1].coord.y() > threshold) {
      // Apply deformation: x_deformed = dF·x
      // square_points[i].coord = dF_ext * square_points[i].coord;
      square_points[i1].coord.x() =
          square_points[i1].coord.x() - h * symmetry_constantx;
    }
    // Points in the lower half remain unchanged
  }

  // Process each alpha value
  for (size_t i = 0; i < alpha_values.size(); i++) {
    double alpha = alpha_values[i];
    std::cout << "\n=== Processing alpha = " << alpha << " ===" << std::endl;

    // Create deformation gradient
    Eigen::Vector2d n1(0., 1.0);
    Eigen::Vector2d a1 = Eigen::Vector2d(1.0, 0.0);
    Eigen::Matrix2d F_ext;
    F_ext << 1.0, alpha, 0.0, 1.0;

    Eigen::Matrix2d dF_ext;
    dF_ext << 1.0, step_size, 0.0, 1.0;

    // Apply deformation only to points above the threshold
    for (size_t i1 = 0; i1 < square_points.size(); i1++) {
      // Only deform points in the upper half of the crystal
      if (square_points[i1].coord.y() > threshold) {
        // Apply deformation: x_deformed = dF·x
        // square_points[i].coord = dF_ext * square_points[i].coord;
        square_points[i1].coord.x() = square_points[i1].coord.x() + step_size;
      }
      // Points in the lower half remain unchanged
    }

    // Create user data
    bool plasticity = false;
    UserData userData(square_points, elements, calculator, potential_func,
                      potential_func_der, zero, optimal_lattice_parameter,
                      F_ext, interior_mapping, full_mapping, active_elements,
                      plasticity);

    // Prepare for optimization
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

    std::cout << "Pre-optimization - Energy: " << pre_energy
              << ", Stress: " << pre_stress << std::endl;

    if (i == 0) {
      ConfigurationSaver::saveConfigurationWithStressAndEnergy2D(
          &userData, i + 999, pre_energy, pre_stress, true);
      ConfigurationSaver::writeToVTK(userData.points, userData.elements,
                                     &userData, i + 999, true);
      ConfigurationSaver::saveTriangleData(&userData, i + 999, domain_dims,
                                           offsets, full_mapping);
    }

    // Store original positions
    alglib::real_1d_array original_x;
    original_x.setlength(x.length());
    for (int j = 0; j < x.length(); j++) {
      original_x[j] = x[j];
    }

    std::vector<int> m3_before =
        analyzeElementReduction(elements, square_points, &userData);
    // Run optimization
    // --- Start timing ---
    auto wall_start = std::chrono::high_resolution_clock::now();
    clock_t cpu_start = clock();

    // Run optimization
    LBFGSOptimizer optimizer(12, 0, pow(10, -13), 0, 0);
    // optimizer.optimize(x, minimize_energy_with_triangles, &userData);

    // --- Stop timing ---
    auto wall_end = std::chrono::high_resolution_clock::now();
    clock_t cpu_end = clock();

    // Compute durations
    double wall_time =
        std::chrono::duration<double>(wall_end - wall_start).count();
    double cpu_time = (double)(cpu_end - cpu_start) / CLOCKS_PER_SEC;

    // Print results
    std::cout << "Optimization wall-clock time: " << wall_time << " seconds\n";
    std::cout << "Optimization CPU time: " << cpu_time << " seconds\n";
    std::cout << "Optimization Ratio: " << cpu_time / wall_time << " seconds\n";

    map_solver_array_to_points(x, square_points, interior_mapping, n_vars);

    std::vector<int> m3_after =
        analyzeElementReduction(elements, square_points, &userData);
    int hasChanges = compareM3Activation(m3_before, m3_after);
    if (hasChanges > 0) {
      std::cout << "Changes in m3 detected! " << hasChanges << std::endl;
    }

    // Calculate post-optimization energy and stress
    double post_energy = 0.0;
    double post_stress = 0.0;
    double post_stress_special = 0.0;

    double post_energy_special = 0.0;
    static double post_stress_previous = 0.;

    static double post_energy_previous = 0.;

    stress_tensor.setZero();

    ConfigurationSaver::calculateEnergyAndStress(&userData, post_energy,
                                                 stress_tensor, true);

    post_stress = stress_tensor(0, 1);
    post_energy_special = post_energy;
    post_stress_special = stress_tensor(0, 1);
    std::cout << "Post-optimization - Energy: " << post_energy
              << ", Stress: " << post_stress << std::endl;
    std::cout << "Energy change: " << (post_energy - pre_energy)
              << ", Stress change: " << (post_stress - pre_stress) << std::endl;

    // Calculate change measures and check for remeshing need
    ChangeMeasures result =
        computeChangeMeasures(x, original_x, lattice_constant, elements,
                              &userData, square_points, true, &F_ext);

    std::cout << "max_abs_change: " << result.max_abs_change << std::endl;
    if (result.has_distorted_triangles) {
      std::cout << "Distorted triangles detected!" << std::endl;
    }

    // Determine if remeshing is needed
    bool shouldRemesh = false; // result.has_distorted_triangles ;
    std::cout << "shouldRemesh: " << shouldRemesh << std::endl;
    int mesh_iteration = 0;
    shouldRemesh = checkSquareDomainViolation(elements);
    // shouldRemesh= true;
    if (!shouldRemesh) {
      std::cout << "No remeshing needed based on domain check." << std::endl;
    }

    // Remeshing if needed
    if (shouldRemesh || i == 0) {
      std::cout << "REMESHING STARTS iteration: " << mesh_iteration++
                << std::endl;

      for (int j = 0; j < x.length(); j++) {
        original_x_remesh[j] = x[j];
      }

      post_energy = 0.0;
      post_stress = 0.0;
      std::vector<int> contact_atoms;
      contact_atoms.resize(0);
      std::vector<int> boundary_fixed_nodes;
      boundary_fixed_nodes.resize(0);
      int max_iterations = 1000;

      auto [post_energy_re, stress_tensor_re, iterations] =
          perform_remeshing_loop_reduction(
              x,
              &userData, // Note: removed & if userData is already a pointer
              contact_atoms, boundary_fixed_nodes, F_ext, dndx, offsets,
              original_domain_map, translation_map, domain_dims_point,
              hasChanges, max_iterations, element_area, pbc, false);

      // this is crucial; to be resolved
      for (auto &element : elements) {
        // element.set_reference_mesh(square_points);
        element.set_dof_mapping(
            full_mapping); // or interior_mapping depending on needs
        // double jac = element.calculate_shape_derivatives(x);  // current
        // positions
      }

      post_energy = post_energy_re;
      post_stress = stress_tensor_re(0, 1);
      post_stress_special = stress_tensor(0, 1);
      post_energy_special = post_energy;

      std::cout << "Post-remeshing - Energy: " << post_energy
                << ", Stress: " << post_stress << std::endl;
    }

    // Update plasticity flag for logging
    bool updated_plasticity = plasticity;

    // Log data to file
    ConfigurationSaver::logEnergyAndStress(
        i, alpha, pre_energy, pre_stress, post_energy, post_stress, hasChanges);

    // Save configuration periodically
    // if(hasChanges > 10 || i % 1000 == 0) {
    UserData finalUserData(square_points, elements, calculator, potential_func,
                           potential_func_der, zero, optimal_lattice_parameter,
                           F_ext, interior_mapping, full_mapping,
                           active_elements, plasticity);

    static int file_counter = 0;
    static bool stress_drop_detected_last_iteration = false;

    // Always save current configuration
    int file_id = i;
    ConfigurationSaver::saveConfigurationWithStressAndEnergy2D(
        &finalUserData, file_id, post_energy, post_stress, true);
    ConfigurationSaver::writeToVTK(finalUserData.points, finalUserData.elements,
                                   &finalUserData, file_id, true);
    ConfigurationSaver::saveTriangleData(&finalUserData, file_id, domain_dims,
                                         offsets, full_mapping);
    std::cout << "Saved configuration " << file_id << " (i=" << i
              << ", stress=" << post_stress << ")" << std::endl;

    post_stress_previous = post_stress_special;
    post_energy_previous = post_energy_special;

    post_energy = 0;
    post_stress = 0;

    std::cout << "caller_id " << caller_id << " completed successfully"
              << std::endl;

    std::cout << "Iteration " << i << " completed successfully" << std::endl;
  }
}
