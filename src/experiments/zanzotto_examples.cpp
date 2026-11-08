#include "../../include/experiments/experiment_includes.h"
#include "../../include/experiments/stability_monitor.h"
#include "../../include/geometry/MacroscopicGeometry.h"
#include "../../include/experiments/pbc_stress_control.h"

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
      create_dof_mapping_original(square_points, 0.25 * lattice_constant, pbc);

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
    relaxation_begin_step(restart_iteration + 1 + static_cast<int>(i));
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
    relax_configuration(x, &userData, 12);

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

    int hasChanges = 0;
    if (shouldRemesh) {
      std::cout << "REMESHING STARTS" << std::endl;

      std::vector<int> contact_atoms;
      std::vector<int> boundary_fixed_nodes;
      int max_iterations = 50;

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
    postOptUserData.second_derivative_function = &potential_func_sder;
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
        post_stress, pre_area, post_area, (hasChanges > 0), &postOptUserData);

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
    unsigned int seed,
    bool enable_remeshing,
    const ConfigurationSaver::CheckpointData* restart_chk,
    int checkpoint_interval,
    double stress_drop_threshold,
    bool save_triangle_data,
    int max_avalanches,
    const std::string& lattice_type,
    double noise_level,
    bool use_rotation,
    int csl_m,
    int csl_n,
    int csl_kx,
    int csl_ky,
    double rotation_angle,
    const std::string& bc,
    int potential_mode,
    const std::string& geometry_json,
    int save_vtk_every) {

  MacroscopicGeometryConfig macro_geom;
  bool has_custom_geom = false;
  if (!geometry_json.empty()) {
    if (MacroscopicGeometryConfig::loadFromJSON(geometry_json, macro_geom) && macro_geom.enabled) {
      has_custom_geom = true;
    }
  }

  bool is_periodic_dislocs = has_custom_geom &&
      macro_geom.microstructure.mode == MicrostructureMode::PeriodicDislocations &&
      macro_geom.microstructure.enforce_periodic_bc;

  bool pbc = is_periodic_dislocs ||
      (!has_custom_geom && (bc == "periodic" || bc == "Periodic" || bc == "pbc" || bc == "PBC"));

  // Parameters for lattice
  if (nx <= 0 || ny <= 0) {
    std::cerr << "Error: nx and ny must be positive integers." << std::endl;
    exit(EXIT_FAILURE);
  }

  // Periodic triangular lattice requires an even number of rows (ny) so that
  // row (ny - 1) (which has a 0.5*a stagger) alternates seamlessly with row 0 (0.0 stagger)
  // across the vertical periodic boundary, preventing artificial boundary dislocations.
  if (pbc && lattice_type == "triangular" && (ny % 2 != 0)) {
    int old_ny = ny;
    ny++;
    std::cout << "[INFO] Enforcing even ny for periodic triangular lattice: adjusted ny from "
              << old_ny << " to " << ny << std::endl;
  }

  // Print potential and calculated elastic constants
  if (potential_mode == 3) {
    auto ec = energy_functions::calculate_elastic_constants(
        energy_functions::current_beta, energy_functions::current_K, lattice_type);
    std::cout << "Zanzotto Potential: beta = " << energy_functions::current_beta
              << ", K = " << energy_functions::current_K << " | Lattice: " << lattice_type << std::endl;
    std::cout << "Ground-state Elastic Constants: C11=" << ec.C11 << ", C12=" << ec.C12 
              << ", C66=" << ec.C66 << ", Zener A=" << ec.Zener_A << std::endl;
  }

  // ==================== SETUP ENERGY POTENTIAL & OPTIMAL LATTICE PARAMETER ====================
  std::function<double(double)> potential_func;
  std::function<double(double)> potential_func_der;
  std::function<double(double)> potential_func_sder;
  double optimal_lattice_parameter = 1.0;
  double r_cutoff = 2.5;

  if (potential_mode == 3) {
    energy_functions::set_parameters(energy_functions::current_beta, energy_functions::current_K);
    optimal_lattice_parameter = (lattice_type == "triangular") ? std::pow(4.0 / 3.0, 0.25) : 1.0;
    potential_func = square_energy;
    potential_func_der = square_energy_der;
    potential_func_sder = square_energy_der;
    std::cout << "[POTENTIAL] Mode 3: Zanzotto Continuum Model (beta = " << energy_functions::current_beta
              << ", K = " << energy_functions::current_K << ")\n";
  } else if (potential_mode == 1) {
    potential_func = lennard_jones_energy_v2;
    potential_func_der = lennard_jones_energy_der_v2;
    potential_func_sder = lennard_jones_energy_sder_v2;
    r_cutoff = 2.5;
    optimal_lattice_parameter = find_optimal_lattice_parameter(potential_func, lattice_type, r_cutoff, potential_func_der, potential_func_sder);
    std::cout << "[POTENTIAL] Mode 1: Lennard-Jones v2 (cutoff = " << r_cutoff << ")\n";
  } else if (potential_mode == 2) {
    potential_func = lennard_jones_energy_v3;
    potential_func_der = lennard_jones_energy_der_v3;
    potential_func_sder = lennard_jones_energy_sder_v3;
    r_cutoff = 1.86602540378444;
    optimal_lattice_parameter = find_optimal_lattice_parameter(potential_func, lattice_type, r_cutoff, potential_func_der, potential_func_sder);
    std::cout << "[POTENTIAL] Mode 2: Lennard-Jones v3 (cutoff = " << r_cutoff << ")\n";
  } else if (potential_mode == 4) {
    potential_func = square_energy;
    potential_func_der = square_energy_der;
    potential_func_sder = square_energy_sder;
    r_cutoff = 2.5;
    optimal_lattice_parameter = find_optimal_lattice_parameter(potential_func, lattice_type, r_cutoff, potential_func_der, potential_func_sder);
    std::cout << "[POTENTIAL] Mode 4: Square 2-Well Gaussian (cutoff = " << r_cutoff << ")\n";
  } else {
    std::cerr << "Unknown potential mode: " << potential_mode << std::endl;
    exit(EXIT_FAILURE);
  }

  double lattice_constant = optimal_lattice_parameter;
  std::cout << "Optimal lattice parameter a_0: " << std::fixed << std::setprecision(8) << lattice_constant << std::endl;

  ReferenceSetup ref_setup = make_reference_setup(
      lattice_type, potential_mode, optimal_lattice_parameter, r_cutoff, potential_func);

  // ==================== SETUP REFERENCE GEOMETRY ====================
  // Save domain dimensions
  writeSizesToFile(nx, ny);

  std::cout << "Remeshing mode: "
            << (enable_remeshing ? "ENABLED (adaptive remeshing)" : "DISABLED (fixed mesh, no remeshing)")
            << std::endl;

  const Eigen::Matrix<double, 3, 2> &dndx = ref_setup.dndx;
  if (std::getenv("MTM_VERBOSE_LOG") != nullptr) {
    std::cout << "Reference element shape derivatives (dN/dx):\n" << dndx << std::endl;
  }

  // ==================== GENERATE INITIAL LATTICE ====================
  std::vector<Point2D> square_points;
  std::vector<Point2D> square_points_ref;
  DomainDimensions domain_dims(0.0, 0.0);
  std::array<double, 2> offsets = {0.0, 0.0};
  double theta_rot = 0.0;

  if (has_custom_geom && macro_geom.microstructure.mode != MicrostructureMode::SingleCrystal) {
    auto [pts, dims] = macro_geom.generateMicrostructureLattice(
        lattice_type, nx, ny, lattice_constant, use_rotation ? rotation_angle : 0.0, potential_mode);
    square_points = pts;
    square_points_ref = pts;
    domain_dims = dims;
    offsets = (lattice_type == "square")
        ? std::array<double, 2>{lattice_constant, lattice_constant}
        : std::array<double, 2>{lattice_constant / 2.0,
                                (sqrt(3.0) / 2.0) * lattice_constant};
    theta_rot = 0.0;
  } else if (use_rotation) {
    if (pbc) {
      int m = (csl_m == 0 && csl_n == 0) ? 1 : csl_m;
      int n = (csl_m == 0 && csl_n == 0) ? 0 : csl_n;
      if (m == 1 && n == 0 && std::abs(rotation_angle) > 0.01) {
        auto best_csl = LatticeGenerator::find_closest_csl_angle(lattice_type, rotation_angle, 80);
        m = best_csl.m;
        n = best_csl.n;
        std::cout << "[CSL] Snapped target angle " << rotation_angle << "° to closest CSL (m="
                  << m << ", n=" << n << ", θ=" << best_csl.theta_deg << "°, Σ=" << best_csl.sigma << ")\n";
      }
      auto cell_1x1 = LatticeGenerator::generate_csl_lattice(lattice_type, m, n, 1, 1, lattice_constant);
      int kx = csl_kx > 0 ? csl_kx : 1;
      int ky = csl_ky > 0 ? csl_ky : 1;
      if (kx <= 2 && ky <= 2 && nx >= 15 && ny >= 15) {
        int auto_kx = std::max(1, (int)std::round((nx * lattice_constant) / cell_1x1.Lx));
        int auto_ky = std::max(1, (int)std::round((ny * lattice_constant) / cell_1x1.Ly));
        if (auto_kx > kx || auto_ky > ky) {
          std::cout << "[CSL] Auto-scaling supercell from default (" << kx << "x" << ky
                    << ") to honor requested system size nx=" << nx << ", ny=" << ny
                    << " -> Kx=" << auto_kx << ", Ky=" << auto_ky << "\n";
          kx = auto_kx;
          ky = auto_ky;
        }
      }
      auto csl_res = LatticeGenerator::generate_csl_lattice(lattice_type, m, n, kx, ky, lattice_constant);
      square_points = csl_res.points;
      square_points_ref = csl_res.points;
      domain_dims = DomainDimensions(csl_res.Lx, csl_res.Ly);
      offsets = {0.0, 0.0};
      theta_rot = csl_res.theta_deg;
    } else {
      // Fixed Boundary Condition: ANY arbitrary rotation angle
      auto [pts, dims] = LatticeGenerator::generate_rotated_fixed_lattice(
          lattice_type, nx, ny, lattice_constant, rotation_angle);
      square_points = pts;
      square_points_ref = pts;
      domain_dims = dims;
      offsets = (lattice_type == "square")
          ? std::array<double, 2>{lattice_constant, lattice_constant}
          : std::array<double, 2>{lattice_constant / 2.0,
                                  (sqrt(3.0) / 2.0) * lattice_constant};
      theta_rot = rotation_angle;
      std::cout << "[Fixed BC] Generated rotated " << lattice_type << " lattice with arbitrary angle θ="
                << rotation_angle << "°: " << pts.size() << " points ("
                << dims.size_x << " x " << dims.size_y << ")\n";
    }
    for (auto& pt : square_points) {
      pt.grain_id = 0;
      pt.grain_theta_deg = theta_rot;
    }
    square_points_ref = square_points;
    // Note: Reference element is always kept unrotated square (p1=(0,0), p2=(h,0), p3=(0,h))
  } else {
    square_points = LatticeGenerator::generate_2d_lattice(
        nx, ny, lattice_constant, lattice_type);
    for (auto& pt : square_points) {
      pt.grain_id = 0;
      pt.grain_theta_deg = 0.0;
    }
    square_points_ref = square_points;
    DomainInfo domain_size = compute_domain_size(square_points);
    domain_dims = DomainDimensions(domain_size.get_width(),
                                   domain_size.get_height());
    offsets = (lattice_type == "square")
        ? std::array<double, 2>{lattice_constant, lattice_constant}
        : std::array<double, 2>{lattice_constant / 2.0,
                                (sqrt(3.0) / 2.0) * lattice_constant};
  }
  ConfigurationSaver::setReferenceRotationDeg(theta_rot);

  if (has_custom_geom && !is_periodic_dislocs) {
    square_points = macro_geom.carveLatticePoints(square_points, lattice_constant);
    square_points_ref = square_points;
  }

  // Stage 2 re-anchoring: X_0 becomes the relaxed Stage-1 state. The specimen frame (x_min, y_min,
  // Lx, Ly) is kept for non-rectangular shapes, whose radius / polygon / cutouts are normalized to it.
  const auto reanchor_frame_bounds = [&](const std::vector<Point2D> &ref) {
    const double frame_x_min = macro_geom.x_min, frame_y_min = macro_geom.y_min;
    const double frame_Lx = macro_geom.Lx, frame_Ly = macro_geom.Ly;
    macro_geom.setPhysicalBoundsFromPoints(ref, optimal_lattice_parameter);
    if (macro_geom.outer_shape != OuterDomainShapeType::Rectangle) {
      macro_geom.x_min = frame_x_min; macro_geom.y_min = frame_y_min;
      macro_geom.Lx = frame_Lx; macro_geom.Ly = frame_Ly;
      if (macro_geom.isIndenterActive()) {
        macro_geom.setupIndenterBounds(optimal_lattice_parameter, &ref);
      }
    }
  };
  // Restart: the lattice is regenerated from the command line, so it must describe the same specimen
  if (restart_chk != nullptr) {
    if (restart_chk->relaxation_step) {
      std::cerr << "[RESTART] Fatal: the checkpoint was written at the relaxation step before the loading (Stage 1, "
                   "relaxation at a prescribed shear strain, or zero-stress relaxation), which is not a load point to "
                   "resume from. Restart from a later checkpoint, or rerun from the start." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (!restart_chk->ref_points.empty() && !(has_custom_geom && !is_periodic_dislocs)) {
      std::cerr << "[RESTART] Fatal: the checkpoint was written by a custom-geometry run. Restart with the "
                   "original run's --geometry-json, --bc, --lattice and rotation flags." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (restart_chk->points.size() != square_points.size()) {
      std::cerr << "[RESTART] Fatal: the checkpoint has " << restart_chk->points.size()
                << " nodes but the regenerated specimen has " << square_points.size()
                << ". Restart with the original run's nx ny, --lattice, --bc, --geometry-json and rotation flags."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  // Restart: the checkpoint carries X_0 (re-anchored after Stage 1); a regenerated lattice does not
  if (restart_chk != nullptr && has_custom_geom && !is_periodic_dislocs) {
    if (restart_chk->ref_points.size() == square_points_ref.size()) {
      for (size_t j = 0; j < square_points_ref.size(); ++j) {
        square_points_ref[j].coord = restart_chk->ref_points[j];
      }
      square_points = square_points_ref;
      if ((macro_geom.microstructure.relax_only || macro_geom.microstructure.pre_relax_gb) &&
          macro_geom.microstructure.reanchor_ref_after_relax) {
        reanchor_frame_bounds(square_points_ref);
      }
    } else if (macro_geom.microstructure.pre_relax_gb && macro_geom.microstructure.reanchor_ref_after_relax) {
      std::cerr << "[RESTART] WARNING: checkpoint has no reference positions (written by an older build); "
                   "X_0 is the unrelaxed lattice, not the Stage-1 state, so driven boundaries may jump."
                << std::endl;
    }
  }

  int original_domain_size = square_points.size();

  // Setup triangulation variables and domain maps
  auto [original_domain_map, translation_map] =
      MeshGenerator::create_domain_maps(original_domain_size, domain_dims,
                                        offsets);

  // Boundary conditions & Two-Stage Protocol (Stage 1 Relaxation -> Stage 2 Loading)
  std::vector<std::pair<int, int>> interior_mapping;
  std::vector<std::pair<int, int>> full_mapping;
  std::vector<int> node_bc_role;
  std::vector<int> baseline_node_bc_role;
  std::vector<Eigen::Vector2d> node_bc_velocity;

  const std::vector<BoundaryZoneRule> original_studio_rules =
      has_custom_geom ? macro_geom.boundary_rules : std::vector<BoundaryZoneRule>{};
  const bool original_indenter_enabled =
      has_custom_geom ? macro_geom.indenter.enabled : false;
  const bool stage1_pre_relax_active =
      has_custom_geom && (restart_chk == nullptr) &&
      (macro_geom.microstructure.relax_only ||
       macro_geom.microstructure.pre_relax_gb);

  if (is_periodic_dislocs) {
    std::tie(interior_mapping, full_mapping) =
        create_dof_mapping_original(square_points, 0.25 * lattice_constant, /*pbc=*/true);
    node_bc_role.assign(square_points.size(), 0);
    baseline_node_bc_role = node_bc_role;
    node_bc_velocity.assign(square_points.size(), Eigen::Vector2d::Zero());
  } else if (has_custom_geom) {
    if (stage1_pre_relax_active) {
      macro_geom.boundary_rules = macro_geom.buildStage1RelaxationRules();
      macro_geom.indenter.enabled = false;
      std::cout << "[Stage 1/2] Equilibrating microstructure at alpha = 0 (Relaxation Step BCs: "
                << macro_geom.boundary_rules.size() << " rule(s), Indenter disengaged, Noise OFF)..." << std::endl;
    }
    macro_geom.buildBoundaryDOFs(square_points, interior_mapping, full_mapping,
                                 node_bc_role, node_bc_velocity);
    baseline_node_bc_role = node_bc_role;
  } else {
    std::tie(interior_mapping, full_mapping) =
        create_dof_mapping_original(square_points, 0.25 * lattice_constant, pbc);
  }

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

  std::vector<ElementTriangle2D> elements;
  std::vector<size_t> active_elements;

  if (has_custom_geom && !is_periodic_dislocs) {
    AdaptiveMesher::setGlobalGeometryFilter(&macro_geom, &square_points_ref);
  } else {
    AdaptiveMesher::setGlobalGeometryFilter(nullptr, nullptr);
  }

  if (restart_chk != nullptr) {
    square_points = restart_chk->points;
    // A checkpoint stores the active elements only; drop the placeholders the loader leaves in the
    // gaps (zero reference area, nodes {0,0,0}) so elements[0] is a real element
    elements.clear();
    active_elements.clear();
    for (size_t e : restart_chk->active_elements) {
      active_elements.push_back(elements.size());
      elements.push_back(restart_chk->elements[e]);
    }
    for (auto &element : elements) {
      element.set_shape_derivatives(dndx);
      // Older checkpoints carry the distorted current areas of the initial mesh
      element.setReferenceArea(ref_setup.element_area);
      element.set_reference_mesh(square_points); // fixed / driven nodes are read from here
      element.set_dof_mapping(full_mapping);
      element.setExternalDeformation(restart_chk->F_ext);
      element.calculate_deformation_gradient(square_points);
    }
    std::cout << "[RESTART] Loaded " << square_points.size() << " nodes, "
              << elements.size() << " elements (" << active_elements.size() << " active)" << std::endl;
  } else {
    auto [mesh_elems, mesh_act_elems] = mesher.createMesh(
        square_points, free_dofs, Eigen::Matrix2d::Identity(), &dndx);
    elements = std::move(mesh_elems);
    active_elements = std::move(mesh_act_elems);
    for (auto &element : elements) {
      element.set_dof_mapping(full_mapping);
    }
  }
  double element_area = elements.empty() ? ref_setup.element_area : elements[0].getReferenceArea();

  // Dead-load tractions need the reference boundary edges, so they are assembled once the mesh exists
  if (has_custom_geom && !is_periodic_dislocs &&
      !macro_geom.finalizeBoundaryTractions(square_points_ref, elements)) {
    std::cerr << "[Traction] Run refused: unbalanced dead load on an unsupported specimen." << std::endl;
    AdaptiveMesher::setGlobalGeometryFilter(nullptr, nullptr);
    std::exit(EXIT_FAILURE);
  }
  const auto rules_have_traction = [](const std::vector<BoundaryZoneRule> &rules) {
    for (const auto &r : rules) {
      if (r.enabled && (r.kind == BoundaryConditionKind::PrescribedTraction ||
                        r.kind == BoundaryConditionKind::RigidForceGrip)) return true;
    }
    return false;
  };
  // PBC nominal-stress control (--stress-control=<preset>, --relax-stress=1 on a periodic cell)
  PBCStressControl stress_ctrl = pbc_stress_control_options();
  const bool stress_ctrl_enabled = pbc && stress_ctrl.enabled();
  if (stress_ctrl.enabled() && !pbc) {
    std::cerr << "[Stress Control] WARNING: --stress-control / --relax-stress apply to periodic cells only; ignored."
              << std::endl;
  }
  const bool stress_preset_is_force = PBCStressControl::isForcePreset(stress_ctrl.preset);
  const bool stress_preset_is_cauchy = PBCStressControl::isCauchyPreset(stress_ctrl.preset);
  // Relaxation step at a nonzero shear strain (shear stress presets): step 0 holds F12 = relax_strain under
  // strain control with noise, then the stress schedule rises from the relaxed stress sigma_0
  const bool stress_relax_strain = stress_ctrl_enabled && stress_ctrl.relax_strain != 0.0 && restart_chk == nullptr &&
                                   !is_periodic_dislocs &&
                                PBCStressControl::basePreset(stress_ctrl.preset) == "shear";
  if (stress_ctrl_enabled && stress_ctrl.relax_strain != 0.0 && !stress_relax_strain) {
    std::cerr << "[Stress Control] WARNING: --relax-strain applies to fresh plain periodic runs with a shear stress "
                 "preset (shear, cauchy_shear); ignored." << std::endl;
  }
  double relax_strain_sigma0 = 0.0; // stress after the relaxation step (start of the stress schedule)
  Eigen::Matrix2d F_bar_state = (restart_chk != nullptr) ? restart_chk->F_ext : Eigen::Matrix2d::Identity();
  Eigen::Matrix2d F_stage1 = (restart_chk != nullptr && restart_chk->has_F_stage1)
      ? restart_chk->F_stage1 : Eigen::Matrix2d::Identity(); // box after the alpha = 0 relaxation (schedule origin)
  Eigen::Matrix2d F_bar_n = F_bar_state, F_bar_nm1 = F_bar_state;
  double F_bar_dalpha_n = 0.0;
  bool F_bar_hist_smooth = false;
  if (stress_ctrl_enabled) {
    // Reference area of the periodic cell (the work A0 Pt : dF_bar is per cell): the triangles tile the
    // current cell, area A0 det F_bar. Not the summed material areas, which exceed it when a dislocation
    // dipole inserts extra half-planes.
    stress_ctrl.A0 = 0.0;
    for (size_t e : active_elements) {
      elements[e].setExternalDeformation(F_bar_state);
      stress_ctrl.A0 += elements[e].calculateCurrentArea(square_points);
    }
    stress_ctrl.A0 /= F_bar_state.determinant();
    std::cout << "[Stress Control] PBC " << (stress_preset_is_cauchy ? "Cauchy (true)" : "nominal (Piola)")
              << "-stress control: preset = '"
              << (stress_ctrl.preset.empty() ? "none" : stress_ctrl.preset) << "', zero-stress Stage 1 = "
              << (stress_ctrl.relax_stage1 ? "yes" : "no") << ", A0 = " << stress_ctrl.A0
              << " | F_bar is upper triangular (F21 = 0)" << std::endl;
  }
  // Primary stress for periodic cells: sigma_xy, or the controlled normal component for uniaxial presets
  const auto pbc_primary_stress = [&](const Eigen::Matrix2d &st) {
    if (stress_ctrl_enabled && stress_preset_is_force) {
      return st(PBCStressControl::targetRow(stress_ctrl.preset), PBCStressControl::targetCol(stress_ctrl.preset));
    }
    return st(0, 1);
  };
  const bool run_has_tractions = stress_ctrl_enabled || (has_custom_geom && !is_periodic_dislocs &&
      (rules_have_traction(macro_geom.boundary_rules) || rules_have_traction(original_studio_rules)));

  std::cout << "[Triangulation] Initial Delaunay mesh: " << original_domain_size
            << " atoms (" << std::fixed << std::setprecision(1) << domain_dims.size_x
            << " x " << domain_dims.size_y << std::defaultfloat << std::setprecision(6)
            << ") -> " << elements.size()
            << " triangles (" << active_elements.size() << " active) | Free nodes: "
            << interior_mapping.size() << "/" << full_mapping.size()
            << " | Remeshing: " << (enable_remeshing ? "ENABLED" : "DISABLED")
            << std::endl;

  // Setup energy calculation
  BaseLatticeCalculator& calculator = *ref_setup.calc;
  double zero = ref_setup.zero;

  std::cout << "Reference state energy (zero): " << zero << std::endl;

  if (stress_ctrl_enabled) {
    // Bulk (Cauchy-Born) Voigt stiffness at the reference state, and the uniaxial-stress Poisson ratio
    // it predicts (small-alpha check for --stress-control=uniaxial_x/_y). F0 maps Z^2 to the lattice.
    const Eigen::Matrix2d &F0 = ref_setup.F0;
    const Eigen::Matrix2d F_tot = F_bar_state * F0;
    const lagrange::Result red_res = reduce_metric(Eigen::Matrix2d(F_tot.transpose() * F_tot));
    AcousticTensor acoustic_tensor(F_tot, red_res.C_reduced, red_res.m_matrix);
    acoustic_tensor.computeEnergyDerivatives(calculator, potential_func_der, potential_func_sder,
                                             calculator.getUnitCellArea());
    itensor::ITensor C_spatial = acoustic_tensor.getSpatialTangent();
    Eigen::Matrix3d C_voigt = Eigen::Matrix3d::Identity();
    auto inds = itensor::inds(C_spatial);
    if (inds.size() == 4) {
      itensor::Index ii = inds[0], jj = inds[1], kk = inds[2], ll = inds[3];
      const int vi[3][2] = {{1, 1}, {2, 2}, {1, 2}};
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
          C_voigt(a, b) = itensor::elt(C_spatial, ii = vi[a][0], jj = vi[a][1], kk = vi[b][0], ll = vi[b][1]);
    }
    const Eigen::Matrix3d S_voigt = C_voigt.inverse();
    std::cout << "[Stress Control] Bulk Voigt stiffness (xx, yy, xy): C11 = " << C_voigt(0, 0) << ", C22 = "
              << C_voigt(1, 1) << ", C12 = " << C_voigt(0, 1) << ", C66 = " << C_voigt(2, 2)
              << " | predicted Poisson ratio under sigma_yy: " << -S_voigt(0, 1) / S_voigt(1, 1)
              << ", under sigma_xx: " << -S_voigt(1, 0) / S_voigt(0, 0) << std::endl;
  }

  // ==================== SETUP LOADING SCHEDULE ====================
  // A restart starts at current_alpha + step_size, so the ratio can land just below an integer and
  // drop the last step; fresh runs keep their historical truncation
  const double schedule_tol = (restart_chk != nullptr) ? 1e-9 : 0.0;
  int num_alpha_points =
      static_cast<int>(std::abs(alpha_max - alpha_min) / std::abs(step_size) + schedule_tol) + 1;
  std::cout << "Loading schedule: " << num_alpha_points << " steps from "
            << alpha_min << " to " << alpha_max << " (step size: " << step_size
            << ")" << std::endl;

  // Generate loading sequence
  std::vector<double> alpha_values;
  if (has_custom_geom && macro_geom.microstructure.relax_only) {
    alpha_values.push_back(0.0);
    std::cout << "[Protocol] Relax-Only mode active: executing single Stage 1 equilibration at alpha = 0.0" << std::endl;
  } else {
    alpha_values.reserve(num_alpha_points + 1);
    if (stage1_pre_relax_active && std::abs(alpha_min) > 1e-12) {
      alpha_values.push_back(0.0);
    }
    for (int i = 0; i < num_alpha_points; i++) {
      alpha_values.push_back(alpha_min + i * step_size);
    }
  }
  if (stress_relax_strain) {
    alpha_values.assign(1, stress_ctrl.relax_strain); // the stress targets are appended after step 0
    std::cout << "[Stress Control] Relaxation step at F12 = " << stress_ctrl.relax_strain
              << " under strain control (noise " << noise_level << "), then "
              << (stress_ctrl.relax_stage1 ? "a zero-stress relaxation and " : "") << stress_ctrl.preset
              << " stress loading from " << (stress_ctrl.relax_stage1 ? "zero" : "the relaxed stress") << " up to "
              << alpha_max << " (step " << std::abs(step_size)
              << "); alpha_start is not used." << std::endl;
  }
  // PBC stress control with --relax-stress: the zero-stress relaxation runs at alpha = 0 also when the
  // schedule starts elsewhere (custom geometries prepend alpha = 0 for their Stage 1 already)
  if (stress_ctrl_enabled && stress_ctrl.relax_stage1 && !stress_relax_strain && !stage1_pre_relax_active &&
      restart_chk == nullptr && !alpha_values.empty() && std::abs(alpha_values.front()) > 1e-15) {
    alpha_values.insert(alpha_values.begin(), 0.0);
  }
  // PBC stress control starting at alpha = 0: that step is a zero-load relaxation, which would relax the
  // initial noise away; the noise is applied at the first loaded step instead (as after a Stage 1)
  const bool stress_zero_step0 = stress_ctrl_enabled && !stress_relax_strain && !stage1_pre_relax_active &&
                                 restart_chk == nullptr && alpha_values.size() > 1 &&
                                 std::abs(alpha_values.front()) < 1e-15;
  // Step 0 relaxes before the loading: its checkpoints are marked, a restart does not resume from them
  const bool step0_is_relaxation =
      restart_chk == nullptr && (stage1_pre_relax_active || stress_relax_strain || stress_zero_step0);

  int file_counter = (restart_chk != nullptr ? restart_chk->iteration + 1 : 0);
  double post_energy_previous = 0.0;
  double post_stress_previous = 0.0;

  if (restart_chk != nullptr) {
    bool plasticity = false;
    UserData checkUserData(square_points, elements, calculator, potential_func,
                           potential_func_der, zero, optimal_lattice_parameter,
                           restart_chk->F_ext, interior_mapping, full_mapping,
                           active_elements, plasticity);
    checkUserData.second_derivative_function = &potential_func_sder;
    const bool custom_bcs = has_custom_geom && !is_periodic_dislocs;
    if (custom_bcs) {
      checkUserData.macro_geom_ptr = &macro_geom;
      checkUserData.current_alpha = restart_chk->current_alpha;
    }
    Eigen::Matrix2d stress_tensor = Eigen::Matrix2d::Zero();
    ConfigurationSaver::calculateEnergyAndStress(&checkUserData, post_energy_previous,
                                                 stress_tensor, true);
    std::vector<int> no_contacts;
    post_stress_previous = custom_bcs
        ? macro_geom.extractPrimaryStress(stress_tensor, &square_points, restart_chk->current_alpha,
                                          &checkUserData, &no_contacts)
        : pbc_primary_stress(stress_tensor);
    std::cout << "[RESTART] Resumed state verified: Energy = " << post_energy_previous
              << ", Stress = " << post_stress_previous << std::endl;
  }

  // Force-control (dead-load) bookkeeping: conjugate strain, elastic secant, predictor, runaway
  double traction_strain_previous =
      (restart_chk != nullptr && has_custom_geom && !is_periodic_dislocs)
          ? macro_geom.computeTractionConjugateStrain(square_points) : 0.0;
  if (restart_chk != nullptr && stress_ctrl_enabled) {
    // PBC stress control: conjugate strain of the resumed cell, from the Stage-1 cell
    const int r = PBCStressControl::targetRow(stress_ctrl.preset), c = PBCStressControl::targetCol(stress_ctrl.preset);
    traction_strain_previous = restart_chk->F_ext(r, c) - F_stage1(r, c);
  }
  double elastic_dstrain_ref = 0.0, elastic_dalpha_ref = 0.0;
  bool have_elastic_ref = false;
  std::vector<Point2D> traction_points_n, traction_points_nm1;
  double traction_dalpha_n = 0.0;
  bool traction_prev_step_smooth = false;
  constexpr int kForceControlMaxIterations = 20000;

  // Lowest stiffness eigenvalues during the run (off unless --eig-every / --eig-at-avalanche)
  StabilityMonitor stability_monitor;
  if (stress_ctrl_enabled) {
    // Free-cell stability at the same states; strain from the Stage-1 cell, as EffectiveStrain
    stability_monitor.set_box_stability([&](UserData &state, int step, double alpha) {
      BoxStabilityPoint p;
      if (step == 0 && stress_relax_strain) return p; // relaxation step at fixed strain: the cell is held
      const int r = PBCStressControl::targetRow(stress_ctrl.preset);
      const int c = PBCStressControl::targetCol(stress_ctrl.preset);
      p.strain = state.F_external(r, c) - F_stage1(r, c);
      p.has_strain = true;
      const bool loaded = step > 0 && stress_preset_is_force;
      const PBCBoxStability b = pbc_box_stability(state, stress_ctrl, loaded ? alpha : 0.0, loaded);
      p.lambda_min = b.lambda_min;
      p.tangent = b.tangent;
      p.has_lambda = b.has_lambda;
      p.has_tangent = b.has_tangent;
      return p;
    });
  }
  int avalanche_counter = 0;

  // Process each alpha value
  int total_steps = static_cast<int>(alpha_values.size()) - 1; // grows after a nonzero-strain relaxation step
  const bool verbose_log = (std::getenv("MTM_VERBOSE_LOG") != nullptr);
  for (size_t i = 0; i < alpha_values.size(); i++) {
    double alpha = alpha_values[i];
    if (verbose_log) {
      std::cout << "\n=== Processing alpha = " << alpha << " ===" << std::endl;
    }
    relaxation_begin_step(static_cast<int>(i));
    double pre_area = 1.0;
    double post_area = 1.0;

    // ==================== SETUP DEFORMATION ====================
    // Stage 1 (i == 0 when stage1_pre_relax_active) is pure zero-strain relaxation (alpha = 0, NO noise).
    // Stage 2 starts at first_loading_step (i == 1 when Stage 1 ran, or i == 0 otherwise):
    // at first_loading_step, if alpha != 0, the entire crystal is affinely deformed by F(alpha)
    // and thermal/positional noise is activated only here at the start of loading.
    const bool is_stage1_relax_step =
        stage1_pre_relax_active && (i == 0) && (restart_chk == nullptr);
    const bool is_stress_zero_relax_step = stress_zero_step0 && (i == 0);
    const size_t first_loading_idx =
        ((stage1_pre_relax_active || stress_zero_step0) && restart_chk == nullptr) ? 1 : 0;
    const bool is_first_loading_step = !is_stage1_relax_step && !is_stress_zero_relax_step &&
                                       (i == first_loading_idx) && (restart_chk == nullptr);
    const double actual_dalpha =
        is_first_loading_step ? alpha
        : (i == 0 ? (restart_chk != nullptr ? step_size : alpha)
                  : (alpha - ((stress_relax_strain && i == 1) ? relax_strain_sigma0 : alpha_values[i - 1])));
    const double effective_noise_level = (is_stage1_relax_step || is_stress_zero_relax_step) ? 0.0 : noise_level;
    // PBC stress control: the alpha = 0 step is a zero-stress box relaxation; stress presets then drive
    const bool stress_relax_step = stress_ctrl_enabled && (i == 0) && std::abs(alpha) < 1e-15;
    const bool stress_relax_strain_step = stress_relax_strain && (i == 0);
    // --relax-stress with --relax-strain: step 1 relaxes the cell to sigma = 0 (all three components) after the
    // relaxation at the prescribed strain, and the stress schedule then rises from 0
    const bool stress_zero_after_strain_step = stress_relax_strain && stress_ctrl.relax_stage1 && (i == 1);
    const bool stress_force_step = stress_ctrl_enabled && !stress_relax_step && !stress_relax_strain_step &&
                                   !stress_zero_after_strain_step && stress_preset_is_force;
    // Force control: dead-load tractions or a nominal-stress target drive this step (never Stage 1)
    const bool force_control_step = stress_force_step || (has_custom_geom && !is_periodic_dislocs &&
                                    !is_stage1_relax_step && macro_geom.hasForceControl());
    set_relaxation_max_iterations(force_control_step ? kForceControlMaxIterations : 0);
    if (has_custom_geom && !is_periodic_dislocs) {
      macro_geom.setRunawayReference(square_points, force_control_step
          ? 0.35 * std::max(macro_geom.Lx, macro_geom.Ly) : 0.0);
    }
    // Last equilibrium state (start of a force-control step): past the limit load this step has no
    // equilibrium, and the run keeps this state as its result instead of the runaway configuration
    std::vector<Point2D> step_start_points;
    std::vector<ElementTriangle2D> step_start_elements;
    std::vector<size_t> step_start_active_elements;
    const Eigen::Matrix2d step_start_F_bar = F_bar_state;
    if (force_control_step) {
      step_start_points = square_points;
      step_start_elements = elements;
      step_start_active_elements = active_elements;
    }

    Eigen::Matrix2d F_ext;
    F_ext << 1.0, alpha, 0.0, 1.0;

    Eigen::Matrix2d dF_ext;
    dF_ext << 1.0, actual_dalpha, 0.0, 1.0;

    if (stress_ctrl_enabled) {
      // Which F_bar components are unknowns, and their nominal targets
      stress_ctrl.P_target.setZero();
      stress_ctrl.sigma_target.setZero();
      stress_ctrl.cauchy = false;
      for (bool &f : stress_ctrl.comp_free) f = false;
      const std::string &preset = stress_ctrl.preset;
      if ((stress_relax_step && stress_ctrl.relax_stage1) || stress_zero_after_strain_step) {
        for (bool &f : stress_ctrl.comp_free) f = true; // sigma_xx = sigma_yy = sigma_xy = 0
      } else if (stress_relax_strain_step) {
        // relaxation step: box held at F12 = relax_strain (strain control), nucleation on the unstable branch
      } else if (stress_preset_is_force) {
        for (bool &f : stress_ctrl.comp_free) f = true;
        const int r = PBCStressControl::targetRow(preset), c = PBCStressControl::targetCol(preset);
        if (stress_preset_is_cauchy) {
          stress_ctrl.cauchy = true; // P_target is re-evaluated from sigma_target in relax_configuration
          stress_ctrl.sigma_target(r, c) = alpha;
          stress_ctrl.sigma_target(c, r) = alpha;
        } else {
          stress_ctrl.P_target(r, c) = alpha;
        }
      } else if (preset == "relax_normal") {
        stress_ctrl.comp_free[0] = true; // F12 = alpha strain-driven, P11 = P22 = 0
        stress_ctrl.comp_free[2] = true;
      }
      stress_ctrl.active = stress_ctrl.comp_free[0] || stress_ctrl.comp_free[1] || stress_ctrl.comp_free[2];

      // Predictor: strain-driven components follow the schedule from the Stage-1 box, stress-driven ones
      // extrapolate the last smooth increment
      Eigen::Matrix2d F_pred = F_bar_state;
      if (!stress_ctrl.active) {
        F_pred = F_stage1;
        F_pred(0, 1) += alpha;
      } else {
        if (!stress_ctrl.comp_free[1]) F_pred(0, 1) = F_stage1(0, 1) + alpha;
        if (!stress_relax_step && !stress_zero_after_strain_step && F_bar_hist_smooth && F_bar_dalpha_n != 0.0) {
          const double ratio = actual_dalpha / F_bar_dalpha_n;
          for (int k = 0; k < 3; ++k) {
            if (!stress_ctrl.comp_free[k]) continue;
            const int r = PBCStressControl::compRow(k), c = PBCStressControl::compCol(k);
            F_pred(r, c) += ratio * (F_bar_n(r, c) - F_bar_nm1(r, c));
          }
        }
      }
      F_pred(1, 0) = 0.0;
      F_ext = F_pred;
      dF_ext = F_pred * F_bar_state.inverse();
      stress_ctrl.F_bar = &F_ext;
      // Cauchy: the nominal target at the predicted cell, so PreEnergy is an enthalpy like PostEnergy
      if (stress_ctrl.cauchy) stress_ctrl.setNominalFromCauchy(F_ext);
      stress_ctrl.F_step_start = F_bar_state;
      stress_ctrl.F_jump_limit = stress_force_step ? 0.5 : 0.0;
    }

    std::vector<int> indenter_contacts;
    if (has_custom_geom && !is_periodic_dislocs) {
      if (macro_geom.isIndenterActive()) {
        node_bc_role = baseline_node_bc_role;
      }
      F_ext = macro_geom.computeMacroscopicF(alpha);
      macro_geom.applyLoadStep(
          square_points, square_points_ref, baseline_node_bc_role, node_bc_velocity,
          alpha, actual_dalpha, is_first_loading_step, seed, effective_noise_level);
      // Pure force control has no affine predictor (F_macro = I): extrapolate the last smooth
      // (no topology change) increment instead, x0 = x_n + (dalpha / dalpha_n) (x_n - x_{n-1})
      if (force_control_step && traction_prev_step_smooth && traction_dalpha_n != 0.0 &&
          traction_points_n.size() == square_points.size() &&
          traction_points_nm1.size() == square_points.size() &&
          macro_geom.computeMacroscopicF(alpha).isIdentity(0.0)) {
        const double ratio = actual_dalpha / traction_dalpha_n;
        for (size_t j = 0; j < square_points.size(); ++j) {
          const int role = baseline_node_bc_role[j];
          const Eigen::Vector2d d = traction_points_n[j].coord - traction_points_nm1[j].coord;
          if (role == 0 || (role == 2 && macro_geom.isGripNode(j))) square_points[j].coord += ratio * d;
          else if (role == 3) square_points[j].coord.y() += ratio * d.y();
          else if (role == 4) square_points[j].coord.x() += ratio * d.x();
        }
      }
      if (macro_geom.isIndenterActive() &&
          macro_geom.indenter.contact_mode == IndenterContactMode::KinematicProjection) {
        macro_geom.updateIndenterContactKinematics(
            square_points, square_points_ref, baseline_node_bc_role, node_bc_role, indenter_contacts, alpha);
        macro_geom.rebuildBoundaryDOFs(node_bc_role, interior_mapping, full_mapping);
        for (auto &elem : elements) {
          elem.set_reference_mesh(square_points);
          elem.set_dof_mapping(full_mapping);
        }
        active_elements = initialize_active_elements(elements, full_mapping, square_points.size());
      }
    } else if (is_stage1_relax_step) {
      // Stage 1 relaxation at alpha = 0 (e.g. periodic dislocations): no affine deformation and no noise
      for (size_t j = 0; j < square_points.size(); j++) {
        square_points[j].coord = F_ext * square_points[j].coord;
      }
    } else if (is_first_loading_step) {
      if (effective_noise_level > 0.0 && !is_periodic_dislocs) {
        // Use deterministic seed so both positive and negative runs share the exact same initial noise
        std::mt19937 gen(seed);
        std::normal_distribution<double> noise_dist(0.0, effective_noise_level);
        std::cout << "[Setup] Seeded loading-stage thermal/positional noise (seed = " << seed
                  << ", σ = " << effective_noise_level << ")" << std::endl;

        for (size_t j = 0; j < square_points.size(); j++) {
          Eigen::Vector2d noise(noise_dist(gen), noise_dist(gen));
          square_points[j].coord = dF_ext * square_points[j].coord + noise;
        }
      } else {
        for (size_t j = 0; j < square_points.size(); j++) {
          square_points[j].coord = dF_ext * square_points[j].coord;
        }
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
    userData.second_derivative_function = &potential_func_sder;
    if (has_custom_geom && !is_periodic_dislocs) {
      userData.macro_geom_ptr = &macro_geom;
      userData.current_alpha = alpha;
    }
    if (stress_ctrl_enabled) {
      userData.stress_ctrl = &stress_ctrl;
      userData.current_alpha = alpha;
    }

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
    pre_stress = (has_custom_geom && !is_periodic_dislocs) ? macro_geom.extractPrimaryStress(stress_tensor, &square_points, alpha, &userData, &indenter_contacts) : pbc_primary_stress(stress_tensor);
    pre_area = ConfigurationSaver::calculateTotalArea2D(&userData);

    if (verbose_log) {
      std::cout << "Pre-optimization - Energy: " << pre_energy
                << ", Stress: " << pre_stress << std::endl;
    }

    // Start in-memory trace for this load step
    AvalancheRecorder::instance().startLoadStep(static_cast<int>(i));
    AvalancheRecorder::instance().recordStep(
        AvalancheEventType::INITIAL_GUESS, 0, pre_energy, pre_stress, 0.0,
        elements.size(), 0, "initial guess after affine deformation increment");

    // Store original positions
    alglib::real_1d_array original_x;
    original_x.setlength(x.length());
    for (int j = 0; j < x.length(); j++) {
      original_x[j] = x[j];
    }

    // ==================== CACHE PRE-OPTIMIZATION STATE IN RAM ====================
    std::vector<Point2D> pre_points = square_points;
    std::vector<ElementTriangle2D> pre_elements = elements;
    std::vector<size_t> pre_active_elements = active_elements;
    Eigen::Matrix2d pre_F_ext = F_ext;
    double pre_saving_value = alpha_values[i];
    double pre_energy_val = pre_energy;
    double pre_stress_val = pre_stress;

    // ==================== RUN OPTIMIZATION ====================
    auto wall_start = std::chrono::high_resolution_clock::now();

    userData.third_condition_flag = false;
    // Optional finite-difference check of the enthalpy gradient (MTM_TRACTION_GRADCHECK=1)
    if (force_control_step && std::getenv("MTM_TRACTION_GRADCHECK") != nullptr && n_vars > 0) {
      alglib::real_1d_array g0, gp, gm, xp;
      g0.setlength(x.length()); gp.setlength(x.length()); gm.setlength(x.length());
      xp.setlength(x.length());
      double f0 = 0.0, fp = 0.0, fm = 0.0;
      minimize_energy_with_triangles(x, f0, g0, &userData);
      std::vector<int> probe;
      for (int k = 0; k < n_vars && probe.size() < 5; ++k) { // loaded contour atoms first
        const int p_idx = interior_mapping[k].first;
        if (baseline_node_bc_role[p_idx] == 0 &&
            macro_geom.cached_node_ext_force_unit[p_idx].squaredNorm() > 0.0) probe.push_back(k);
      }
      for (int k = n_vars / 2; k < n_vars && probe.size() < 10; k += std::max(1, n_vars / 17)) {
        if (baseline_node_bc_role[interior_mapping[k].first] == 0) probe.push_back(k);
      }
      double worst = 0.0;
      const double h = 1e-6;
      double g_scale = 0.0; // relative-error floor: FD round-off is ~1e-10 where the gradient vanishes
      for (int q = 0; q < x.length(); ++q) g_scale = std::max(g_scale, std::abs(g0[q]));
      for (int k : probe) {
        for (int comp = 0; comp < 2; ++comp) {
          const int idx = comp * n_vars + interior_mapping[k].second;
          for (int q = 0; q < x.length(); ++q) xp[q] = x[q];
          xp[idx] = x[idx] + h; minimize_energy_with_triangles(xp, fp, gp, &userData);
          xp[idx] = x[idx] - h; minimize_energy_with_triangles(xp, fm, gm, &userData);
          const double fd = (fp - fm) / (2.0 * h);
          const double rel = std::abs(fd - g0[idx]) / std::max(1e-4 * g_scale, std::abs(g0[idx]) + std::abs(fd));
          if (std::getenv("MTM_TRACTION_GRADCHECK_VERBOSE") != nullptr) {
            std::cout << "    dof " << idx << " (node " << interior_mapping[k].first << ", comp " << comp
                      << ") grad = " << std::scientific << g0[idx] << " fd = " << fd << " rel = " << rel
                      << std::defaultfloat << std::endl;
          }
          worst = std::max(worst, rel);
        }
      }
      std::cout << "[Traction GradCheck] alpha = " << alpha << " | " << 2 * probe.size()
                << " DOFs | max relative error = " << std::scientific << worst << std::defaultfloat << std::endl;
    }
    relax_configuration(x, &userData, 13);

    auto wall_end = std::chrono::high_resolution_clock::now();
    double wall_time =
        std::chrono::duration<double>(wall_end - wall_start).count();

    map_solver_array_to_points(x, square_points, interior_mapping, n_vars);

    int hasChanges = 0;

    // ==================== POST-OPTIMIZATION ENERGY ====================
    double post_energy = 0.0;
    double post_stress = 0.0;
    stress_tensor.setZero();
    ConfigurationSaver::calculateEnergyAndStress(&userData, post_energy,
                                                 stress_tensor, true);
    post_stress = (has_custom_geom && !is_periodic_dislocs) ? macro_geom.extractPrimaryStress(stress_tensor, &square_points, alpha, &userData, &indenter_contacts) : pbc_primary_stress(stress_tensor);
    post_area = ConfigurationSaver::calculateTotalArea2D(&userData);

    AvalancheRecorder::instance().recordStep(
        AvalancheEventType::LBFGS_CONVERGED, -1, post_energy, post_stress, 0.0,
        elements.size(), 0, "initial continuous relaxation converged");

    // Enthalpy after the continuous relaxation, before any reconnection: the fixed-alpha reference
    // for the force-controlled avalanche test (Pi drops on every step under a dead load otherwise)
    const double enthalpy_before_cascade = post_energy;
    // Cauchy control: H after the cascade uses another nominal target; correct the comparison by the work
    // of the target change (PBCStressControl::cauchyWorkCorrection)
    const Eigen::Matrix2d F_before_cascade = F_ext;
    const Eigen::Matrix2d P_before_cascade = stress_ctrl.P_target;
    const auto cascade_enthalpy_change = [&]() {
      return post_energy - enthalpy_before_cascade +
             (stress_ctrl_enabled ? stress_ctrl.cauchyWorkCorrection(F_before_cascade, P_before_cascade, F_ext,
                                                                    stress_ctrl.P_target)
                                  : 0.0);
    };

    // Past the limit load the continuous relaxation itself runs away: skip the remesh cascade
    const bool runaway_before_cascade = force_control_step &&
        (macro_geom.exceedsRunawayLimit(square_points) || (stress_force_step && stress_ctrl.exceedsLimit(F_ext)));

    // ==================== REMESHING DECISION ====================
    bool shouldRemesh = enable_remeshing && !runaway_before_cascade;
    if (const char *gate = std::getenv("MTM_REMESH_GATE")) {
      if (std::string(gate) == "energy") // legacy gate, for A/B comparison only
        shouldRemesh = enable_remeshing && (post_energy < post_energy_previous || i == 0);
    }
    int remesh_iterations = 0;

    if (shouldRemesh) {
      auto remesh_start = std::chrono::high_resolution_clock::now();
      std::vector<int> contact_atoms = indenter_contacts;
      std::vector<int> boundary_fixed_nodes;
      int max_iterations = 50;

      auto [post_energy_re, stress_tensor_re, iterations] =
          perform_remeshing_loop_reduction(
              x, &userData, contact_atoms, boundary_fixed_nodes, F_ext, dndx,
              offsets, original_domain_map, translation_map, domain_dims_point,
              hasChanges, max_iterations, element_area, pbc, true);

      auto remesh_end = std::chrono::high_resolution_clock::now();
      wall_time += std::chrono::duration<double>(remesh_end - remesh_start).count();

      square_points = userData.points;
      elements = userData.elements;
      active_elements = userData.active_elements;

      // Update element_area after remeshing
      if (!elements.empty()) {
        element_area = elements[0].getReferenceArea();
      }

      // If KinematicProjection is active and remeshing occurred, update the active contact set
      // in case plastic pile-up brought new surface atoms into contact during the avalanche
      if (hasChanges > 0 && has_custom_geom && macro_geom.isIndenterActive() &&
          macro_geom.indenter.contact_mode == IndenterContactMode::KinematicProjection) {
        size_t prev_contact_count = indenter_contacts.size();
        macro_geom.updateIndenterContactKinematics(
            square_points, square_points_ref, baseline_node_bc_role, node_bc_role, indenter_contacts, alpha);
        if (indenter_contacts.size() != prev_contact_count) {
          macro_geom.rebuildBoundaryDOFs(node_bc_role, interior_mapping, full_mapping);
          for (auto &element : elements) {
            element.set_reference_mesh(square_points);
            element.set_dof_mapping(full_mapping);
          }
          active_elements = initialize_active_elements(elements, full_mapping, square_points.size());
          UserData activeSetUserData(square_points, elements, calculator, potential_func,
                                     potential_func_der, zero, optimal_lattice_parameter,
                                     F_ext, interior_mapping, full_mapping, active_elements,
                                     plasticity);
          activeSetUserData.second_derivative_function = &potential_func_sder;
          activeSetUserData.macro_geom_ptr = &macro_geom;
          activeSetUserData.current_alpha = alpha;

          n_vars = interior_mapping.size();
          x.setlength(2 * n_vars);
          map_points_to_solver_array(x, square_points, interior_mapping, n_vars);
          relax_configuration(x, &activeSetUserData, 13);
          map_solver_array_to_points(x, square_points, interior_mapping, n_vars);

          ConfigurationSaver::calculateEnergyAndStress(&activeSetUserData, post_energy_re,
                                                       stress_tensor_re, true);
        }
      }

      remesh_iterations = std::abs(iterations);
      post_energy = post_energy_re;
      post_stress = (has_custom_geom && !is_periodic_dislocs) ? macro_geom.extractPrimaryStress(stress_tensor_re, &square_points, alpha, &userData, &indenter_contacts) : pbc_primary_stress(stress_tensor_re);

      // Update DOF mapping
      for (auto &element : elements) {
        element.set_reference_mesh(square_points);
        element.set_dof_mapping(full_mapping);
      }
    }

    // ==================== CONCISE PER-STEP SUMMARY ====================
    {
      std::stringstream mesh_status;
      if (!enable_remeshing) {
        mesh_status << "fixed (" << active_elements.size() << " tris)";
      } else if (hasChanges > 0) {
        mesh_status << "remeshed (" << remesh_iterations
                    << (remesh_iterations == 1 ? " pass, " : " passes, ")
                    << active_elements.size() << " tris)";
      } else if (remesh_iterations > 0) {
        mesh_status << "trial rejected (kept " << active_elements.size() << " tris)";
      } else {
        mesh_status << "unchanged (" << active_elements.size() << " tris)";
      }

      std::cout << "Step " << std::setw(4) << i << "/" << total_steps
                << " | alpha = " << std::fixed << std::setprecision(6) << alpha
                << " | E = " << std::scientific << std::setprecision(5) << post_energy
                << " | stress = " << std::scientific << std::setprecision(5) << post_stress
                << " | mesh: " << mesh_status.str()
                << " | " << std::fixed << std::setprecision(2) << wall_time << "s"
                << std::endl;
    }

    // ==================== CHECK FOR STRESS DROP ====================
    bool energy_dropped = (post_energy < post_energy_previous);
    bool stress_dropped =
        (std::abs(post_stress) < std::abs(post_stress_previous));

    // Option 2: Fractional drops (for output / investigation)
    double fractional_energy_drop = 0.0;
    double fractional_stress_drop = 0.0;
    if (std::abs(post_energy_previous) > 1e-4) {
      fractional_energy_drop =
          (post_energy_previous - post_energy) / std::abs(post_energy_previous);
    }
    if (std::abs(post_stress_previous) > 1e-4) {
      fractional_stress_drop =
          (std::abs(post_stress_previous) - std::abs(post_stress)) /
          std::abs(post_stress_previous);
    }

    // Avalanche acceptance conditions:
    // 1) Remeshing accepted topological changes (hasChanges > 0) [if remeshing enabled]
    // 2) Energy decreased (post_energy < post_energy_previous)
    // 3) Stress magnitude dropped by at least threshold (fractional_stress_drop >= stress_drop_threshold)
    // (not at the zero-stress relaxation after a relaxation at prescribed strain: unloading, not an avalanche)
    const bool can_detect_avalanche = (i > 0) && !is_first_loading_step && !stress_zero_after_strain_step;
    bool remesh_accepted = (hasChanges > 0);
    bool stress_drop_criterion = (fractional_stress_drop >= stress_drop_threshold);
    bool stress_drop_detected = enable_remeshing
        ? (can_detect_avalanche && remesh_accepted && energy_dropped && stress_drop_criterion)
        : (can_detect_avalanche && energy_dropped && stress_drop_criterion);

    // Force control: an avalanche is a strain burst at fixed load, i.e. a reconnection that lowers the
    // enthalpy at this alpha plus a conjugate-strain increment beyond the elastic (secant) estimate.
    double traction_strain = 0.0, traction_dstrain = 0.0, fractional_strain_burst = 0.0;
    bool runaway_detected = false;
    if (force_control_step) {
      if (stress_force_step) {
        // Conjugate strain of the controlled component, measured from the Stage-1 box
        const int r = PBCStressControl::targetRow(stress_ctrl.preset);
        const int c = PBCStressControl::targetCol(stress_ctrl.preset);
        traction_strain = F_ext(r, c) - F_stage1(r, c);
      } else {
        traction_strain = macro_geom.computeTractionConjugateStrain(square_points);
      }
      traction_dstrain = traction_strain - traction_strain_previous;
      const double expected_dstrain = (have_elastic_ref && elastic_dalpha_ref != 0.0)
          ? elastic_dstrain_ref * (actual_dalpha / elastic_dalpha_ref) : 0.0;
      const double excess = std::abs(traction_dstrain) - std::abs(expected_dstrain);
      if (have_elastic_ref && std::abs(expected_dstrain) > 1e-15) {
        fractional_strain_burst = excess / std::abs(expected_dstrain);
      }
      const bool strain_burst = have_elastic_ref && excess > 1e-12 &&
                                fractional_strain_burst >= stress_drop_threshold;
      const bool enthalpy_dropped = cascade_enthalpy_change() < 0.0;
      stress_drop_detected = can_detect_avalanche && strain_burst &&
                             (enable_remeshing ? (remesh_accepted && enthalpy_dropped) : true);

      // Limit load: past the yield / Peierls stress the dead load has no equilibrium
      const double max_jump = macro_geom.maxDisplacementSinceStepStart(square_points);
      const double runaway_limit = macro_geom.runaway_limit;
      const bool box_runaway = stress_force_step && stress_ctrl.exceedsLimit(F_ext);
      if (runaway_before_cascade || max_jump > runaway_limit || box_runaway) {
        runaway_detected = true;
        stress_drop_detected = true; // not an elastic increment
        std::cout << std::defaultfloat << std::setprecision(6) << "  [LIMIT LOAD / YIELD RUNAWAY] alpha = " << alpha;
        if (stress_force_step) {
          std::cout << ": a stress-driven F_bar component changed by " << stress_ctrl.maxStressDrivenJump(F_ext)
                    << " in one step (limit " << stress_ctrl.F_jump_limit << ").";
        } else {
          std::cout << ": max nodal displacement " << max_jump << " in one step exceeds " << runaway_limit
                    << " (0.35 max(Lx, Ly)).";
        }
        std::cout << " No equilibrium at this load: keeping the last equilibrium state." << std::endl;
      }
      if (!stress_drop_detected && actual_dalpha != 0.0) {
        elastic_dstrain_ref = traction_dstrain;
        elastic_dalpha_ref = actual_dalpha;
        have_elastic_ref = true;
      }
      traction_strain_previous = traction_strain;
      traction_points_nm1 = std::move(traction_points_n);
      traction_points_n = square_points;
      traction_dalpha_n = actual_dalpha;
      traction_prev_step_smooth = (hasChanges == 0) && !stress_drop_detected;
    }
    if (runaway_detected) {
      // Limit load: no equilibrium at alpha. Restore the last equilibrium (start of this step), write it
      // as the final configuration and checkpoint, and log nothing for alpha (the runaway state is not
      // a point of the curve).
      const double alpha_eq = (i > 0) ? ((stress_relax_strain && i == 1) ? relax_strain_sigma0 : alpha_values[i - 1])
                                      : 0.0;
      square_points = std::move(step_start_points);
      elements = std::move(step_start_elements);
      active_elements = std::move(step_start_active_elements);
      F_ext = stress_ctrl_enabled ? step_start_F_bar
          : (has_custom_geom && !is_periodic_dislocs ? macro_geom.computeMacroscopicF(alpha_eq) : F_ext);
      for (auto &element : elements) {
        element.set_reference_mesh(square_points);
        element.set_dof_mapping(full_mapping);
        element.setExternalDeformation(F_ext);
      }
      UserData eqUserData(square_points, elements, calculator, potential_func, potential_func_der, zero,
                          optimal_lattice_parameter, F_ext, interior_mapping, full_mapping, active_elements,
                          plasticity);
      eqUserData.second_derivative_function = &potential_func_sder;
      if (has_custom_geom && !is_periodic_dislocs) {
        eqUserData.macro_geom_ptr = &macro_geom;
        eqUserData.current_alpha = alpha_eq;
      }
      const int eq_file_id = caller_id + file_counter;
      file_counter++;
      double eq_energy = post_energy_previous, eq_stress = post_stress_previous;
      ConfigurationSaver::saveConfigurationWithStressAndEnergy2D(&eqUserData, eq_file_id, eq_energy, eq_stress,
                                                                 true);
      auto [eq_dislocations, eq_coordination] = DefectAnalysis::analyzeDefectsInReferenceConfig(
          &eqUserData, eq_file_id, dndx, offsets, original_domain_map, translation_map, domain_dims_point,
          element_area, pbc, true);
      ConfigurationSaver::writeToVTK(square_points, elements, &eqUserData, eq_file_id, true, eq_coordination,
                                     alpha_eq);
      ConfigurationSaver::CheckpointData chk_eq;
      chk_eq.nx = nx; chk_eq.ny = ny; chk_eq.iteration = eq_file_id;
      chk_eq.current_alpha = alpha_eq; chk_eq.step_size = step_size; chk_eq.alpha_end = alpha_max;
      chk_eq.mode = (step_size > 0 ? "positive" : "negative");
      chk_eq.seed = seed; chk_eq.enable_remeshing = enable_remeshing;
      chk_eq.F_ext = F_ext; chk_eq.points = square_points;
      chk_eq.F_stage1 = F_stage1; chk_eq.has_F_stage1 = stress_ctrl_enabled;
      chk_eq.relaxation_step = (i == 0 && step0_is_relaxation);
      if (has_custom_geom && !is_periodic_dislocs) {
        for (const auto &p : square_points_ref) chk_eq.ref_points.push_back(p.coord);
      }
      chk_eq.elements = elements; chk_eq.active_elements = active_elements;
      std::stringstream chk_eq_ss;
      chk_eq_ss << "checkpoints/checkpoint_" << std::setw(5) << std::setfill('0') << eq_file_id << ".chk";
      ConfigurationSaver::saveCheckpoint(chk_eq_ss.str(), chk_eq);
      ConfigurationSaver::saveCheckpoint("checkpoints/latest.chk", chk_eq);
      std::cout << std::defaultfloat << std::setprecision(6)
                << "[Force Control] Limit load between alpha = " << alpha_eq << " and " << alpha
                << ". Last equilibrium state (alpha = " << alpha_eq << ") saved -> configuration_"
                << std::setw(5) << std::setfill('0') << eq_file_id << std::setfill(' ')
                << ".vtk (dislocations: " << eq_dislocations << "), checkpoints/latest.chk. Schedule stopped cleanly." << std::endl;
      break;
    }

    if (stress_ctrl_enabled) {
      F_bar_state = F_ext;
      if (stress_relax_step) F_stage1 = F_ext;
      F_bar_nm1 = F_bar_n;
      F_bar_n = F_ext;
      F_bar_dalpha_n = actual_dalpha;
      F_bar_hist_smooth = !stress_relax_step && !stress_zero_after_strain_step && (hasChanges == 0) &&
                          !stress_drop_detected;
    }

    UserData postOptUserData(square_points, elements, calculator,
                             potential_func, potential_func_der, zero,
                             optimal_lattice_parameter, F_ext, interior_mapping,
                             full_mapping, active_elements, plasticity);
    postOptUserData.second_derivative_function = &potential_func_sder;
    if (has_custom_geom && !is_periodic_dislocs) {
      postOptUserData.macro_geom_ptr = &macro_geom;
      postOptUserData.current_alpha = alpha;
    }
    if (stress_ctrl_enabled) postOptUserData.stress_ctrl = &stress_ctrl;
    post_area = ConfigurationSaver::calculateTotalArea2D(&postOptUserData);

    // ==================== GRAIN BOUNDARY ENERGY SUMMARY (alpha = 0) ====================
    if (i == 0 && has_custom_geom &&
        macro_geom.microstructure.mode == MicrostructureMode::Bicrystal &&
        !square_points.empty()) {
      const double a0_gb = (optimal_lattice_parameter > 1e-8) ? optimal_lattice_parameter : 1.0;
      const double phi_rad = macro_geom.microstructure.gb_angle_deg * M_PI / 180.0;
      const Eigen::Vector2d t_gb(-std::sin(phi_rad), std::cos(phi_rad));
      const Eigen::Vector2d p_gb(
          macro_geom.x_min + macro_geom.microstructure.gb_cx * macro_geom.Lx,
          macro_geom.y_min + macro_geom.microstructure.gb_cy * macro_geom.Ly);

      double s_min = 1e30, s_max = -1e30;
      for (const auto &pt : square_points) {
        const double s = (pt.coord - p_gb).dot(t_gb);
        s_min = std::min(s_min, s);
        s_max = std::max(s_max, s);
      }
      const double L_gb_phys = std::max(a0_gb, s_max - s_min);
      const double L_gb_a0 = L_gb_phys / a0_gb;
      const double collar_phys = std::min(4.0 * a0_gb, 0.15 * L_gb_phys);
      const double L_int_a0 = std::max(a0_gb, L_gb_phys - 2.0 * collar_phys) / a0_gb;

      double interior_energy = 0.0;
      const double norm_area = calculator.getUnitCellArea();
      for (size_t elem_idx : active_elements) {
        if (elem_idx >= elements.size()) continue;
        auto &el = elements[elem_idx];
        if (!el.isInitialized()) continue;
        const Eigen::Vector2d c_el =
            (square_points[el.getNodeIndex(0)].coord +
             square_points[el.getNodeIndex(1)].coord +
             square_points[el.getNodeIndex(2)].coord) / 3.0;
        const double s_c = (c_el - p_gb).dot(t_gb);
        if (s_c >= s_min + collar_phys && s_c <= s_max - collar_phys) {
          el.setExternalDeformation(F_ext);
          el.calculate_deformation_gradient(square_points);
          const Eigen::Matrix2d F_el = el.getDeformationGradient();
          const Eigen::Matrix2d C_el = F_el.transpose() * F_el;
          const auto red = reduce_metric(C_el);
          const double e_el = calculator.calculate_energy(
                                  red.C_reduced, potential_func, zero) /
                              norm_area;
          interior_energy += e_el * el.getReferenceArea();
        }
      }

      const double gamma_gb_total = post_energy / L_gb_a0;
      const double gamma_gb_interior = interior_energy / L_int_a0;
      const double th1_deg = macro_geom.microstructure.theta1_deg;
      const double th2_deg = macro_geom.microstructure.theta2_deg;
      const double misori_deg = std::abs(th2_deg - th1_deg);
      const double half_tilt_deg = 0.5 * misori_deg;

      std::cout << "\n======================================================================\n"
                << "[GRAIN BOUNDARY ENERGY SUMMARY (alpha = 0)]\n"
                << "  Misorientation Δθ       = " << std::fixed << std::setprecision(4) << misori_deg
                << "° (θ1 = " << th1_deg << "°, θ2 = " << th2_deg
                << "°, φ_GB = " << macro_geom.microstructure.gb_angle_deg << "°)\n"
                << "  GB Length L_GB          = " << std::fixed << std::setprecision(3) << L_gb_a0
                << " a0 (Interior Gauge L_int = " << L_int_a0 << " a0)\n"
                << "  Total Excess Energy E   = " << std::scientific << std::setprecision(6) << post_energy << "\n"
                << "  γ_GB (Total E / L_GB)   = " << std::fixed << std::setprecision(6) << gamma_gb_total << "\n"
                << "  γ_GB (Interior Gauge)   = " << std::fixed << std::setprecision(6) << gamma_gb_interior << "\n"
                << "======================================================================\n"
                << std::endl;

      // 1) Write per-run gb_energy_summary.csv
      {
        std::ofstream gb_out("gb_energy_summary.csv");
        gb_out << "MisorientationDeg,HalfAngleDeg,Theta1Deg,Theta2Deg,GBAngleDeg,Nx,Ny,L_GB_a0,L_Int_a0,TotalEnergy,InteriorEnergy,GBEnergyTotal,GBEnergyInterior,NormalGapA0,SeamDMinA0\n";
        gb_out << std::fixed << std::setprecision(8)
               << misori_deg << "," << half_tilt_deg << "," << th1_deg << "," << th2_deg << ","
               << macro_geom.microstructure.gb_angle_deg << "," << nx << "," << ny << ","
               << L_gb_a0 << "," << L_int_a0 << ","
               << post_energy << "," << interior_energy << ","
               << gamma_gb_total << "," << gamma_gb_interior << ","
               << macro_geom.microstructure.gb_normal_gap_a0 << ","
               << macro_geom.microstructure.min_seam_distance_a0 << "\n";
      }

      // 2) Update cumulative ../gb_energy_vs_misorientation.csv sorted by MisorientationDeg
      {
        const std::string cum_path = "../gb_energy_vs_misorientation.csv";
        struct GBRow {
          double misori;
          std::string line;
        };
        std::vector<GBRow> rows;
        std::string header = "MisorientationDeg,HalfAngleDeg,Theta1Deg,Theta2Deg,GBAngleDeg,Nx,Ny,L_GB_a0,L_Int_a0,TotalEnergy,InteriorEnergy,GBEnergyTotal,GBEnergyInterior,NormalGapA0,SeamDMinA0";
        if (std::filesystem::exists(cum_path)) {
          std::ifstream in_cum(cum_path);
          std::string ln;
          while (std::getline(in_cum, ln)) {
            if (ln.empty() || ln.rfind("MisorientationDeg", 0) == 0) continue;
            std::stringstream ss(ln);
            std::string first_tok;
            if (std::getline(ss, first_tok, ',')) {
              try {
                double m_val = std::stod(first_tok);
                if (std::abs(m_val - misori_deg) > 0.02) {
                  rows.push_back({m_val, ln});
                }
              } catch (...) {}
            }
          }
        }
        std::ostringstream new_ln;
        new_ln << std::fixed << std::setprecision(8)
               << misori_deg << "," << half_tilt_deg << "," << th1_deg << "," << th2_deg << ","
               << macro_geom.microstructure.gb_angle_deg << "," << nx << "," << ny << ","
               << L_gb_a0 << "," << L_int_a0 << ","
               << post_energy << "," << interior_energy << ","
               << gamma_gb_total << "," << gamma_gb_interior << ","
               << macro_geom.microstructure.gb_normal_gap_a0 << ","
               << macro_geom.microstructure.min_seam_distance_a0;
        rows.push_back({misori_deg, new_ln.str()});
        std::sort(rows.begin(), rows.end(), [](const GBRow &a, const GBRow &b) {
          return a.misori < b.misori;
        });
        std::ofstream out_cum(cum_path);
        out_cum << header << "\n";
        for (const auto &r : rows) {
          out_cum << r.line << "\n";
        }
      }
    }

    // Load parameter written to the logs. The relaxation step at a prescribed shear strain runs at
    // alpha = F12, but the loading that follows is a stress schedule: log its relaxed stress sigma_0 there,
    // so alpha is a stress throughout (the strain is in EffectiveStrain / Fbar12)
    const double alpha_log = stress_relax_strain_step ? post_stress : alpha;

    // Log fractional drops to separate CSV for investigation
    static std::ofstream frac_file;
    static std::filesystem::path frac_dir;
    if (!frac_file.is_open() || frac_dir != std::filesystem::current_path()) {
      if (frac_file.is_open()) {
        frac_file.close();
      }
      frac_file.open("fractional_drops.csv");
      frac_file << "Iteration,Alpha,PostEnergy,EnergyDrop,FractionalEnergyDrop,"
                << "PostStress,StressMagnitudeDrop,FractionalStressDrop,"
                << "RemeshIterations,hasChanges,StressDropDetected\n";
      frac_file << std::scientific
                << std::setprecision(std::numeric_limits<double>::max_digits10);
      frac_dir = std::filesystem::current_path();
    }
    frac_file << i << "," << alpha_log << "," << post_energy << ","
              << (post_energy_previous - post_energy) << ","
              << fractional_energy_drop << "," << post_stress << ","
              << (std::abs(post_stress_previous) - std::abs(post_stress)) << ","
              << fractional_stress_drop << "," << remesh_iterations << ","
              << hasChanges << "," << stress_drop_detected << "\n";
    frac_file.flush();

    const bool is_initial_loaded_snapshot = is_first_loading_step && (std::abs(alpha) > 1e-12);
    if (stress_drop_detected || i == 0 || is_initial_loaded_snapshot) {
      int eq_file_id = caller_id + file_counter;
      file_counter++;
      int pre_file_id = eq_file_id;
      int post_file_id = eq_file_id;
      int num_dislocations_pre = 0;

      // 1) Write PRE-avalanche state from cached RAM copy (with distinct pre_avalanche_XXXXX prefix)
      if (stress_drop_detected) {
        UserData preOptUserData(pre_points, pre_elements, calculator,
                                potential_func, potential_func_der, zero,
                                optimal_lattice_parameter, pre_F_ext,
                                interior_mapping, full_mapping,
                                pre_active_elements, plasticity);
        preOptUserData.second_derivative_function = &potential_func_sder;
        if (has_custom_geom) {
          preOptUserData.macro_geom_ptr = &macro_geom;
          preOptUserData.current_alpha = pre_saving_value;
        }
        if (stress_ctrl_enabled) preOptUserData.stress_ctrl = &stress_ctrl;

        double dummy_pre_e = pre_energy_val, dummy_pre_s = pre_stress_val;
        ConfigurationSaver::saveConfigurationWithStressAndEnergy2D(
            &preOptUserData, pre_file_id, dummy_pre_e, dummy_pre_s, true);

        // Save pre-avalanche single-file checkpoint under pre_avalanche_XXXXX.chk
        {
          ConfigurationSaver::CheckpointData chk_pre;
          chk_pre.nx = nx; chk_pre.ny = ny; chk_pre.iteration = pre_file_id;
          chk_pre.current_alpha = pre_saving_value; chk_pre.step_size = step_size; chk_pre.alpha_end = alpha_max;
          chk_pre.mode = (step_size > 0 ? "positive" : "negative");
          chk_pre.seed = seed; chk_pre.enable_remeshing = enable_remeshing;
          chk_pre.F_ext = pre_F_ext; chk_pre.points = pre_points;
          chk_pre.F_stage1 = F_stage1; chk_pre.has_F_stage1 = stress_ctrl_enabled;
          chk_pre.relaxation_step = (i == 0 && step0_is_relaxation);
          if (has_custom_geom && !is_periodic_dislocs) {
            for (const auto &p : square_points_ref) chk_pre.ref_points.push_back(p.coord);
          }
          chk_pre.elements = pre_elements; chk_pre.active_elements = pre_active_elements;

          std::stringstream chk_ss;
          chk_ss << "checkpoints/pre_avalanche_" << std::setw(5) << std::setfill('0') << pre_file_id << ".chk";
          ConfigurationSaver::saveCheckpoint(chk_ss.str(), chk_pre);
        }

        std::vector<int> coordination_pre;
        std::tie(num_dislocations_pre, coordination_pre) =
            DefectAnalysis::analyzeDefectsInReferenceConfig(
                &preOptUserData, pre_file_id, dndx, offsets, original_domain_map,
                translation_map, domain_dims_point, element_area, pbc, true);

        std::stringstream pre_vtk_ss;
        pre_vtk_ss << "vtk_output/pre_avalanche_" << std::setw(5) << std::setfill('0') << pre_file_id << ".vtk";
        ConfigurationSaver::writeToVTKFile(
            pre_vtk_ss.str(), preOptUserData.points, preOptUserData.elements,
            &preOptUserData, true, coordination_pre, pre_saving_value);
        ConfigurationSaver::logDislocationData(alpha, num_dislocations_pre);
      }

      // 2) Write POST-avalanche (equilibrium) state to configuration_XXXXX.vtk
      double dummy_post_e = post_energy, dummy_post_s = post_stress;
      ConfigurationSaver::saveConfigurationWithStressAndEnergy2D(
          &postOptUserData, post_file_id, dummy_post_e, dummy_post_s, true);
      if (save_triangle_data) {
        ConfigurationSaver::saveTriangleData(&postOptUserData, post_file_id,
                                             domain_dims, offsets, full_mapping);
        ConfigurationSaver::saveElements(elements, active_elements, post_file_id);
      }

      // Save post-avalanche single-file checkpoint & latest.chk
      {
        ConfigurationSaver::CheckpointData chk_post;
        chk_post.nx = nx; chk_post.ny = ny; chk_post.iteration = post_file_id;
        chk_post.current_alpha = alpha; chk_post.step_size = step_size; chk_post.alpha_end = alpha_max;
        chk_post.mode = (step_size > 0 ? "positive" : "negative");
        chk_post.seed = seed; chk_post.enable_remeshing = enable_remeshing;
        chk_post.F_ext = F_ext; chk_post.points = square_points;
        chk_post.F_stage1 = F_stage1; chk_post.has_F_stage1 = stress_ctrl_enabled;
        chk_post.relaxation_step = (i == 0 && step0_is_relaxation);
        if (has_custom_geom && !is_periodic_dislocs) {
          for (const auto &p : square_points_ref) chk_post.ref_points.push_back(p.coord);
        }
        chk_post.elements = elements; chk_post.active_elements = active_elements;

        std::stringstream chk_ss;
        chk_ss << "checkpoints/checkpoint_" << std::setw(5) << std::setfill('0') << post_file_id << ".chk";
        ConfigurationSaver::saveCheckpoint(chk_ss.str(), chk_post);
        ConfigurationSaver::saveCheckpoint("checkpoints/latest.chk", chk_post);
      }

      auto [num_dislocations_post, coordination_post] =
          DefectAnalysis::analyzeDefectsInReferenceConfig(
              &postOptUserData, post_file_id, dndx, offsets,
              original_domain_map, translation_map, domain_dims_point,
              element_area, pbc, true);

      ConfigurationSaver::writeToVTK(
          postOptUserData.points, postOptUserData.elements, &postOptUserData,
          post_file_id, true, coordination_post, alpha_log);
      ConfigurationSaver::logDislocationData(alpha_log, num_dislocations_post);

      if ((i == 0 || is_initial_loaded_snapshot) && !stress_drop_detected) {
        std::cout << "  💾 [VTK Snapshot] "
                  << (is_initial_loaded_snapshot ? "Initial Stage 2 loaded state (α = " + std::to_string(alpha) + ")" : "Initial equilibrium state")
                  << " saved -> configuration_"
                  << std::setw(5) << std::setfill('0') << post_file_id << ".vtk" << std::setfill(' ')
                  << " (dislocations: " << num_dislocations_post << ")" << std::endl;
      } else if (force_control_step) {
        std::cout << "  ⚡ [AVALANCHE #" << avalanche_counter << "] Strain burst at fixed load: Δε = "
                  << std::scientific << std::setprecision(3) << traction_dstrain << " ("
                  << std::fixed << std::setprecision(1) << (fractional_strain_burst * 100.0)
                  << "% above elastic), ΔΠ(cascade) = " << std::scientific << std::showpos
                  << cascade_enthalpy_change() << std::noshowpos
                  << " | Dislocations: " << num_dislocations_pre << " -> " << num_dislocations_post
                  << " | Saved pre_avalanche_" << std::setw(5) << std::setfill('0') << pre_file_id
                  << ".vtk & configuration_" << std::setw(5) << std::setfill('0') << post_file_id
                  << ".vtk" << std::setfill(' ') << std::endl;
      } else {
        std::cout << "  ⚡ [AVALANCHE #" << avalanche_counter << "] Stress dropped "
                  << std::fixed << std::setprecision(1) << (fractional_stress_drop * 100.0) << "% ("
                  << std::scientific << std::setprecision(3) << post_stress_previous
                  << " -> " << post_stress << "), ΔE = "
                  << std::showpos << (post_energy - post_energy_previous) << std::noshowpos
                  << " | Dislocations: " << num_dislocations_pre << " -> " << num_dislocations_post
                  << " | Saved pre_avalanche_" << std::setw(5) << std::setfill('0') << pre_file_id
                  << ".vtk & configuration_" << std::setw(5) << std::setfill('0') << post_file_id
                  << ".vtk" << std::setfill(' ') << std::endl;
      }

      AvalancheRecorder::instance().recordStep(
          AvalancheEventType::AVALANCHE_COMPLETE, -1, post_energy, post_stress, 0.0,
          elements.size(), 0, "avalanche cascade fully completed");
      AvalancheRecorder::instance().commitToFile("avalanche_trace", pre_file_id, post_file_id);

    } else {
      if (hasChanges > 0) {
        // Remeshing occurred and was accepted even if fractional stress drop didn't meet 10% threshold
        AvalancheRecorder::instance().recordStep(
            AvalancheEventType::AVALANCHE_COMPLETE, -1, post_energy, post_stress, 0.0,
            elements.size(), 0, "remeshing cascade completed (sub-avalanche threshold)");
        AvalancheRecorder::instance().commitToFile("avalanche_trace");
      } else {
        // Pure elastic step: discard in-memory buffer without disk I/O
        AvalancheRecorder::instance().discard();
      }

      bool is_periodic_vtk = (save_vtk_every > 0 && (i + 1) % save_vtk_every == 0);
      if (is_periodic_vtk) {
        int periodic_file_id = caller_id + file_counter;
        file_counter++;

        double dummy_per_e = post_energy, dummy_per_s = post_stress;
        ConfigurationSaver::saveConfigurationWithStressAndEnergy2D(
            &postOptUserData, periodic_file_id, dummy_per_e, dummy_per_s, true);
        if (save_triangle_data) {
          ConfigurationSaver::saveTriangleData(&postOptUserData, periodic_file_id,
                                               domain_dims, offsets, full_mapping);
          ConfigurationSaver::saveElements(elements, active_elements, periodic_file_id);
        }

        auto [num_dislocations_post, coordination_post] =
            DefectAnalysis::analyzeDefectsInReferenceConfig(
                &postOptUserData, periodic_file_id, dndx, offsets,
                original_domain_map, translation_map, domain_dims_point,
                element_area, pbc, true);

        ConfigurationSaver::writeToVTK(
            postOptUserData.points, postOptUserData.elements, &postOptUserData,
            periodic_file_id, true, coordination_post, alpha);
        ConfigurationSaver::logDislocationData(alpha, num_dislocations_post);

        std::cout << "  💾 [VTK Snapshot] Step " << (i + 1)
                  << " (alpha = " << std::fixed << std::setprecision(6) << alpha
                  << ") -> VTK #" << std::setw(5) << std::setfill('0') << periodic_file_id
                  << std::setfill(' ') << " (dislocations: " << num_dislocations_post << ")"
                  << std::endl;
      }

      // Elastic step: periodically save single-file checkpoint
      if (checkpoint_interval > 0 && (i + 1) % checkpoint_interval == 0) {
        ConfigurationSaver::CheckpointData chk_periodic;
        chk_periodic.nx = nx; chk_periodic.ny = ny;
        chk_periodic.iteration = (restart_chk != nullptr ? restart_chk->iteration + static_cast<int>(i) + 1 : static_cast<int>(i) + 1);
        chk_periodic.current_alpha = alpha; chk_periodic.step_size = step_size; chk_periodic.alpha_end = alpha_max;
        chk_periodic.mode = (step_size > 0 ? "positive" : "negative");
        chk_periodic.seed = seed; chk_periodic.enable_remeshing = enable_remeshing;
        chk_periodic.F_ext = F_ext; chk_periodic.points = square_points;
        chk_periodic.F_stage1 = F_stage1; chk_periodic.has_F_stage1 = stress_ctrl_enabled;
        chk_periodic.relaxation_step = (i == 0 && step0_is_relaxation);
        if (has_custom_geom && !is_periodic_dislocs) {
          for (const auto &p : square_points_ref) chk_periodic.ref_points.push_back(p.coord);
        }
        chk_periodic.elements = elements; chk_periodic.active_elements = active_elements;

        std::stringstream chk_ss;
        chk_ss << "checkpoints/checkpoint_step_" << std::setw(5) << std::setfill('0') << chk_periodic.iteration << ".chk";
        ConfigurationSaver::saveCheckpoint(chk_ss.str(), chk_periodic);
        ConfigurationSaver::saveCheckpoint("checkpoints/latest.chk", chk_periodic);
      }
    }

    // ==================== STABILITY MONITOR (optional) ====================
    stability_monitor.end_of_step(static_cast<int>(i), alpha_log, post_stress, postOptUserData,
                                  stress_drop_detected, pre_elements,
                                  pre_active_elements, caller_id + file_counter - 1);

    // ==================== LOG DATA ====================
    if (stress_relax_strain_step || stress_zero_after_strain_step) traction_strain = F_ext(0, 1) - F_stage1(0, 1);
    std::vector<std::pair<std::string, double>> traction_columns;
    if (run_has_tractions) {
      const double external_potential = stress_ctrl_enabled ? stress_ctrl.workTerm(F_ext)
          : (force_control_step ? macro_geom.computeExternalForcePotential(square_points, alpha) : 0.0);
      traction_columns = {{"EffectiveStrain", traction_strain},
                          {"ElasticEnergy", post_energy - external_potential},
                          {"StrainBurstFraction", fractional_strain_burst}};
      if (stress_ctrl_enabled) {
        traction_columns.push_back({"Fbar11", F_ext(0, 0)});
        traction_columns.push_back({"Fbar12", F_ext(0, 1)});
        traction_columns.push_back({"Fbar22", F_ext(1, 1)});
      }
      // Avalanche statistics under force control: the strain jump of the step, the enthalpy change of the
      // reconnection cascade at fixed load (Cauchy-corrected), and the avalanche flag; 0 on other steps
      traction_columns.push_back({"StrainJump", force_control_step ? traction_dstrain : 0.0});
      traction_columns.push_back({"CascadeEnthalpyChange", force_control_step ? cascade_enthalpy_change() : 0.0});
      traction_columns.push_back({"Avalanche", stress_drop_detected ? 1.0 : 0.0});
    }
    ConfigurationSaver::logEnergyAndStress_v2(
        i, alpha_log, pre_energy, pre_stress, post_energy, post_stress, pre_area,
        post_area, (hasChanges > 0), &postOptUserData, &indenter_contacts,
        run_has_tractions ? &traction_columns : nullptr);

    post_energy_previous = post_energy;
    post_stress_previous = post_stress;

    if (stress_relax_strain_step) {
      // Stress schedule from the relaxed state: sigma_0 + k |dalpha| up to alpha_end. With --relax-stress, a
      // zero-stress relaxation comes first and the schedule rises from 0 instead.
      relax_strain_sigma0 = post_stress;
      traction_strain_previous = traction_strain;
      const double d = std::abs(step_size);
      const double sigma_start = stress_ctrl.relax_stage1 ? 0.0 : relax_strain_sigma0;
      if (stress_ctrl.relax_stage1) alpha_values.push_back(0.0);
      for (int k = 1; sigma_start + k * d <= alpha_max + 1e-9 * d; ++k) {
        alpha_values.push_back(sigma_start + k * d);
      }
      total_steps = static_cast<int>(alpha_values.size()) - 1;
      std::cout << std::defaultfloat << std::setprecision(6) << "[Stress Control] Relaxation step at F12 = " << alpha
                << " relaxed: sigma = " << relax_strain_sigma0 << "."
                << (stress_ctrl.relax_stage1 ? " Relaxing the cell to zero stress, then switching to "
                                             : " Switching to ")
                << stress_ctrl.preset << " stress control from sigma = " << sigma_start << ": "
                << (alpha_values.size() - 1) << " steps up to " << alpha_max << "." << std::endl;
    }
    if (stress_zero_after_strain_step) traction_strain_previous = traction_strain; // strain origin of the loading


    // ==================== STAGE 1 -> STAGE 2 PROTOCOL HANDOVER ====================
    if (has_custom_geom && i == 0 && restart_chk == nullptr) {
      if (macro_geom.microstructure.relax_only || alpha_values.size() <= 1) {
        std::cout << "[Stage 1/1] Equilibration complete at alpha = 0 (Relax Only mode). Stopping cleanly." << std::endl;
        break;
      }
      if (is_periodic_dislocs) {
        if (macro_geom.microstructure.reanchor_ref_after_relax) {
          square_points_ref = square_points;
        }
        std::cout << "[Stage 2/2] Handover complete: Periodic dislocations equilibrated in full 2D PBCs at alpha = 0; starting homogeneous periodic simple shear F_ext(alpha)..." << std::endl;
      } else if (stage1_pre_relax_active) {
        macro_geom.indenter.enabled = original_indenter_enabled;
        if (macro_geom.microstructure.reanchor_ref_after_relax) {
          square_points_ref = square_points;
          reanchor_frame_bounds(square_points_ref);
        } else if (macro_geom.isIndenterActive()) {
          macro_geom.setupIndenterBounds(optimal_lattice_parameter, &square_points);
        }
        macro_geom.boundary_rules = original_studio_rules;
        if (macro_geom.isIndenterActive()) {
          macro_geom.primary_stress_component = "indentation_load";
        }
        macro_geom.buildBoundaryDOFs(square_points_ref, interior_mapping, full_mapping,
                                     node_bc_role, node_bc_velocity);
        baseline_node_bc_role = node_bc_role;
        for (auto &elem : elements) {
          elem.set_dof_mapping(full_mapping);
        }
        active_elements = initialize_active_elements(elements, full_mapping, square_points.size());
        if (!macro_geom.finalizeBoundaryTractions(square_points_ref, elements)) {
          std::cerr << "[Traction] Run refused at the Stage 2 handover: unbalanced dead load on an "
                       "unsupported specimen." << std::endl;
          AdaptiveMesher::setGlobalGeometryFilter(nullptr, nullptr);
          std::exit(EXIT_FAILURE);
        }
        traction_strain_previous = macro_geom.computeTractionConjugateStrain(square_points);
        traction_points_n.clear();
        traction_points_nm1.clear();
        traction_prev_step_smooth = false;
        have_elastic_ref = false;
        std::cout << "[Stage 2/2] Handover complete: X_0 re-anchored to relaxed equilibrium state, switched to Stage 2 BCs & Indenter ("
                  << macro_geom.boundary_rules.size() << " boundary rules"
                  << (macro_geom.isIndenterActive()
                      ? (", Nanoindenter ACTIVE at y_top = " + std::to_string(macro_geom.indenter.physical_y_top))
                      : "")
                  << ")." << std::endl;
      }
    }

    if (stress_drop_detected || i == 0) {
      avalanche_counter++;
      if (max_avalanches > 0 && avalanche_counter >= max_avalanches) {
        std::cout << "\n======================================================\n"
                  << ">>> Reached target of " << avalanche_counter << " / " << max_avalanches
                  << " avalanches. Stopping simulation. <<<\n"
                  << "======================================================\n" << std::endl;
        break;
      }
    }
  }
  ConfigurationSaver::setReferenceRotationDeg(theta_rot);
  set_relaxation_max_iterations(0);

  if (has_custom_geom) {
    AdaptiveMesher::setGlobalGeometryFilter(nullptr, nullptr);
  }
}

void example_1_conti_zanzotto_negative_loading(int caller_id, int nx, int ny,
                                               double alpha_min,
                                               double alpha_max,
                                               double step_size,
                                               unsigned int seed,
                                               bool enable_remeshing,
                                               const ConfigurationSaver::CheckpointData* restart_chk,
                                               int checkpoint_interval,
                                               double stress_drop_threshold,
                                               bool save_triangle_data,
                                               int max_avalanches,
                                               const std::string& lattice_type,
                                               double noise_level,
                                               bool use_rotation,
                                               int csl_m,
                                               int csl_n,
                                               int csl_kx,
                                               int csl_ky,
                                               double rotation_angle,
                                               const std::string& bc,
                                               int potential_mode,
                                               const std::string& geometry_json,
                                               int save_vtk_every) {
  example_1_conti_zanzotto_loading(caller_id, nx, ny,
                                   /*alpha_min=*/alpha_min,
                                   /*alpha_max=*/alpha_max,
                                   /*step_size=*/step_size,
                                   /*triangulation_perturbation=*/-1e-7,
                                   /*seed=*/seed,
                                   /*enable_remeshing=*/enable_remeshing,
                                   restart_chk,
                                   checkpoint_interval,
                                   stress_drop_threshold,
                                   save_triangle_data,
                                   max_avalanches,
                                   lattice_type,
                                   noise_level,
                                   use_rotation,
                                   csl_m,
                                   csl_n,
                                   csl_kx,
                                   csl_ky,
                                   rotation_angle,
                                   bc,
                                   potential_mode,
                                   geometry_json,
                                   save_vtk_every);
}

void restart_zanzotto_simulation(const std::string& checkpoint_path,
                                 int checkpoint_interval,
                                 double stress_drop_threshold,
                                 bool save_triangle_data,
                                 int max_avalanches,
                                 const std::string& lattice_type,
                                 int potential_mode,
                                 double noise_level,
                                 bool use_rotation,
                                 int csl_m,
                                 int csl_n,
                                 int csl_kx,
                                 int csl_ky,
                                 double rotation_angle,
                                 const std::string& bc,
                                 const std::string& geometry_json,
                                 int save_vtk_every) {
  std::string path = checkpoint_path;
  if (path == "latest" || path == "checkpoints/latest") {
    path = "checkpoints/latest.chk";
  }

  double optimal_lattice_parameter = 1.0;
  double r_cutoff = 2.5;
  std::function<double(double)> pot_fn = square_energy;
  if (potential_mode == 3) {
    optimal_lattice_parameter = (lattice_type == "triangular") ? std::pow(4.0 / 3.0, 0.25) : 1.0;
    pot_fn = square_energy;
  } else if (potential_mode == 1) {
    r_cutoff = 2.5;
    pot_fn = lennard_jones_energy_v2;
    optimal_lattice_parameter = find_optimal_lattice_parameter(
        lennard_jones_energy_v2, lattice_type, r_cutoff, lennard_jones_energy_der_v2, lennard_jones_energy_sder_v2);
  } else if (potential_mode == 2) {
    r_cutoff = 1.86602540378444;
    pot_fn = lennard_jones_energy_v3;
    optimal_lattice_parameter = find_optimal_lattice_parameter(
        lennard_jones_energy_v3, lattice_type, r_cutoff, lennard_jones_energy_der_v3, lennard_jones_energy_sder_v3);
  } else if (potential_mode == 4) {
    r_cutoff = 2.5;
    pot_fn = square_energy;
    optimal_lattice_parameter = find_optimal_lattice_parameter(
        square_energy, lattice_type, r_cutoff, square_energy_der, square_energy_sder);
  }

  ReferenceSetup ref_setup = make_reference_setup(
      lattice_type, potential_mode, optimal_lattice_parameter, r_cutoff, pot_fn);
  const Eigen::Matrix<double, 3, 2> &dndx = ref_setup.dndx;

  ConfigurationSaver::CheckpointData chk;
  if (!ConfigurationSaver::loadCheckpoint(path, chk, dndx)) {
    std::cerr << "Fatal Error: Could not load checkpoint from " << path << std::endl;
    exit(EXIT_FAILURE);
  }

  double next_alpha = chk.current_alpha + chk.step_size;
  std::cout << "\n" << std::string(60, '=') << std::endl;
  std::cout << "RESTARTING SIMULATION FROM: " << path << std::endl;
  std::cout << "  - System: " << chk.nx << "x" << chk.ny << " (" << chk.mode << " loading)" << std::endl;
  std::cout << "  - Resuming from iteration " << chk.iteration << " at alpha = " << next_alpha << std::endl;
  std::cout << "  - Target alpha_end = " << chk.alpha_end << " (step size: " << chk.step_size << ")" << std::endl;
  std::cout << std::string(60, '=') << "\n" << std::endl;

  if (chk.mode == "positive") {
    example_1_conti_zanzotto_loading(
        0, chk.nx, chk.ny,
        next_alpha, chk.alpha_end, chk.step_size,
        0.0, chk.seed, chk.enable_remeshing, &chk,
        checkpoint_interval, stress_drop_threshold, save_triangle_data,
        max_avalanches,
        lattice_type,
        noise_level, use_rotation, csl_m, csl_n, csl_kx, csl_ky, rotation_angle, bc,
        potential_mode, geometry_json, save_vtk_every);
  } else {
    example_1_conti_zanzotto_negative_loading(
        0, chk.nx, chk.ny,
        next_alpha, chk.alpha_end, chk.step_size,
        chk.seed, chk.enable_remeshing, &chk,
        checkpoint_interval, stress_drop_threshold, save_triangle_data,
        max_avalanches,
        lattice_type,
        noise_level, use_rotation, csl_m, csl_n, csl_kx, csl_ky, rotation_angle, bc,
        potential_mode, geometry_json, save_vtk_every);
  }
}
