#include "../../include/experiments/experiment_includes.h"

void example_3_stress_controlled_final_clean(int caller_id, int nx, int ny) {

  // =================================================================
  // APPLIED STRESS CONTROL EXAMPLE (WITH PRE-STEP REMESHING)
  // =================================================================

  if (nx <= 0 || ny <= 0) {
    std::cerr << "Error: nx and ny must be positive integers." << std::endl;
    exit(EXIT_FAILURE);
  }

  // --- SETUP (Assuming all types/functions are defined elsewhere) ---
  // ... (Lines 1-50: Initial setup, dndx, potential setup, etc. are implicitly
  // here) ...

  // ==================== SETUP REFERENCE GEOMETRY ====================
  writeSizesToFile(nx, ny);

  std::string lattice_type = "square";
  double h = 1.0;

  Eigen::Vector2d p1(0.0, 0.0);
  Eigen::Vector2d p2(h, 0.0);
  Eigen::Vector2d p3(0.0, h);

  Eigen::Matrix<double, 3, 2> dndx = calculateShapeDerivatives(p1, p2, p3);

  // ==================== SETUP ENERGY POTENTIAL ====================
  std::function<double(double)> potential_func = square_energy;
  std::function<double(double)> potential_func_der = square_energy_der;
  std::function<double(double)> potential_func_sder = square_energy_der;

  // ==================== LATTICE PARAMETER ====================
  double symmetry_constantx =
      (lattice_type == "triangular") ? pow(4.0 / 3.0, 1.0 / 4.0) : 1.0;
  double optimal_lattice_parameter = symmetry_constantx * 1.0;
  double lattice_constant = optimal_lattice_parameter;

  // ==================== GENERATE INITIAL LATTICE ====================
  std::vector<Point2D> square_points = LatticeGenerator::generate_2d_lattice(
      nx, ny, lattice_constant, lattice_type);
  std::vector<Point2D> square_points_ref =
      LatticeGenerator::generate_2d_lattice(nx, ny, lattice_constant,
                                            lattice_type);

  int original_domain_size = square_points.size();

  DomainInfo domain_size = compute_domain_size(square_points);

  const std::array<double, 2> offsets =
      (lattice_type == "square")
          ? std::array<double, 2>{lattice_constant, lattice_constant}
          : std::array<double, 2>{lattice_constant / 2.0,
                                  (sqrt(3.0) / 2.0) * lattice_constant};

  DomainDimensions domain_dims(domain_size.get_width(),
                               domain_size.get_height());
  bool pbc = true;

  auto [original_domain_map, translation_map] =
      MeshGenerator::create_domain_maps(original_domain_size, domain_dims,
                                        offsets);

  auto [interior_mapping, full_mapping] =
      create_dof_mapping_original(square_points, 0.5 * lattice_constant, pbc);

  Point2D domain_dims_point(domain_dims.size_x, domain_dims.size_y);

  // ==================== CREATE MESH ====================
  AdaptiveMesher mesher(domain_dims_point, offsets, original_domain_map,
                        translation_map, full_mapping, 1e-6, pbc);
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

  // ==================== SETUP ENERGY CALCULATION ====================
  Strain_Energy_LatticeCalculator calculator(1.0);
  double normalisation = calculator.getUnitCellArea();

  Eigen::Matrix2d F_I = Eigen::Matrix2d::Identity();
  Eigen::Matrix2d C_I = F_I.transpose() * F_I;
  double zero = calculator.calculate_energy(C_I, potential_func, 0);

  // Static logging variables
  static int file_counter = 0;
  static double post_energy_previous = 0.0;

  // ================================================================
  // PHASE 1: STRAIN-CONTROLLED EQUILIBRATION
  // ================================================================
  std::cout
      << "\n================================================================"
      << std::endl;
  std::cout << "PHASE 1: STRAIN-CONTROLLED EQUILIBRATION" << std::endl;
  std::cout
      << "================================================================"
      << std::endl;

  // Set initial deformation gradient
  Eigen::Matrix2d F_bar = Eigen::Matrix2d::Identity();
  double initial_shear = 0.14;
  F_bar(0, 1) = initial_shear;

  // --- ESSENTIAL DECLARATION ---
  // F_bar_old stores the converged F_bar from the *previous* load step.
  Eigen::Matrix2d F_bar_old = F_bar;
  // -----------------------------

  // Apply F_bar to reference positions + noise
  {
    std::random_device rd;
    std::mt19937 gen(rd());
    double noise_level = 0.04;
    std::normal_distribution<double> noise_dist(0.0, noise_level);

    for (size_t j = 0; j < square_points.size(); j++) {
      Eigen::Vector2d noise(noise_dist(gen), noise_dist(gen));
      square_points[j].coord = F_bar * square_points_ref[j].coord + noise;
    }

    // Minimize energy for the initial shear step
    bool plasticity = false;
    UserData userData(square_points, elements, calculator, potential_func,
                      potential_func_der, zero, optimal_lattice_parameter,
                      F_bar, interior_mapping, full_mapping, active_elements,
                      plasticity);

    alglib::real_1d_array x;
    int n_vars = interior_mapping.size();
    x.setlength(2 * n_vars);
    map_points_to_solver_array(x, square_points, interior_mapping, n_vars);

    userData.third_condition_flag = false;
    LBFGSOptimizer optimizer(13, 0, 0, 0, 0);
    optimizer.optimize(x, minimize_energy_with_triangles, &userData);

    map_solver_array_to_points(x, square_points, interior_mapping, n_vars);

    // Capture initial state for logging
    double equil_energy = 0.0;
    Eigen::Matrix2d equil_stress = Eigen::Matrix2d::Zero();
    ConfigurationSaver::calculateEnergyAndStress(&userData, equil_energy,
                                                 equil_stress, true);

    post_energy_previous = equil_energy;
    // Initialize F_bar_old with the converged state from Phase 1
    F_bar_old = F_bar;

    std::cout << "Equilibration complete! Initial F_12 = " << F_bar(0, 1)
              << std::endl;
  }

  // ================================================================
  // PHASE 2: STRESS-CONTROLLED LOADING
  // ================================================================
  std::cout
      << "\n================================================================"
      << std::endl;
  std::cout << "PHASE 2: STRESS-CONTROLLED LOADING" << std::endl;
  std::cout
      << "================================================================"
      << std::endl;

  // ==================== STRESS LOADING SCHEDULE ====================
  double sigma_app_min = 0.0;
  double sigma_app_max = 1.5;
  double sigma_app_step = 0.001;

  int num_stress_points =
      static_cast<int>((sigma_app_max - sigma_app_min) / sigma_app_step) + 1;

  std::vector<double> sigma_app_values;
  for (int i = 0; i < num_stress_points; i++) {
    sigma_app_values.push_back(sigma_app_min + i * sigma_app_step);
  }

  // ==================== STRESS CONTROL PARAMETERS ====================
  int max_stress_iterations = 100;
  double stress_tolerance = 1e-5;
  double alpha_relax = 0.1;

  // ------------------------------------------------------------
  // Configuration Flag (Set outside the iteration loop)
  // ------------------------------------------------------------
  const bool use_multiplicative_update =
      false; // Set to 'false' for standard additive update
  // ------------------------------------------------------------

  // ==================== STRESS CONTROL LOOP ====================
  for (size_t stress_idx = 0; stress_idx < sigma_app_values.size();
       stress_idx++) {
    double sigma_app_xy = sigma_app_values[stress_idx];

    Eigen::Matrix2d sigma_app = Eigen::Matrix2d::Zero();
    sigma_app(0, 1) = sigma_app_xy;
    sigma_app(1, 0) = sigma_app_xy;

    std::cout << "\n========================================" << std::endl;
    std::cout << "TARGET STRESS σ_xy = " << sigma_app_xy << std::endl;
    std::cout << "========================================" << std::endl;

    // --- LOGGING: CAPTURE PRE-STATE ---
    double pre_energy = 0.0;
    double pre_F_12 = F_bar(0, 1);
    double pre_area = 0.0;

    {
      bool plasticity = false;
      UserData preLogData(square_points, elements, calculator, potential_func,
                          potential_func_der, zero, optimal_lattice_parameter,
                          F_bar, interior_mapping, full_mapping,
                          active_elements, plasticity);

      Eigen::Matrix2d pre_stress_tensor = Eigen::Matrix2d::Zero();
      ConfigurationSaver::calculateEnergyAndStress(&preLogData, pre_energy,
                                                   pre_stress_tensor, true);
      pre_area = ConfigurationSaver::calculateTotalArea2D(&preLogData);
    }
    // ------------------------------------

    // ============================================================
    // ★★★ REMESHING BEFORE STRESS CONVERGENCE ★★★
    // ============================================================
    std::cout << "\n--- PRE-STEP REMESHING ---" << std::endl;

    std::vector<int> contact_atoms;
    std::vector<int> boundary_fixed_nodes;
    int hasChanges = 0;
    int max_remesh_iterations = 10;

    // Need current state in solver array format for remesh function
    alglib::real_1d_array x_temp;
    int n_vars = interior_mapping.size();
    x_temp.setlength(2 * n_vars);
    map_points_to_solver_array(x_temp, square_points, interior_mapping, n_vars);

    // Create temporary UserData for remeshing (F_bar is the current affine
    // state)
    bool plasticity_remesh = false;
    UserData remeshUserData(square_points, elements, calculator, potential_func,
                            potential_func_der, zero, optimal_lattice_parameter,
                            F_bar, interior_mapping, full_mapping,
                            active_elements, plasticity_remesh);

    auto [post_energy_re, stress_tensor_re, iterations] =
        perform_remeshing_loop_reduction(
            x_temp, // Optimized node positions
            &remeshUserData, contact_atoms, boundary_fixed_nodes, F_bar, dndx,
            offsets, original_domain_map, translation_map, domain_dims_point,
            hasChanges, max_remesh_iterations, element_area, pbc, true);

    // Update square_points, elements, and active_elements (done by reference
    // inside remeshing function)

    // Restore point positions from the solver array after remeshing
    map_solver_array_to_points(x_temp, square_points, interior_mapping, n_vars);

    // Set the energy of the newly optimized mesh/lattice as the starting energy
    // for the inner loop
    pre_energy = post_energy_re;

    std::cout << "--- Remeshing Complete. Starting Stress Convergence. ---"
              << std::endl;

    // Variables for inner loop results
    Eigen::Matrix2d final_stress = Eigen::Matrix2d::Zero();
    double post_energy = 0.0;
    double post_area = 0.0;
    bool converged = false;

    // Final F_12 must be visible after the iteration loop
    double final_F_12 = F_bar(0, 1);

    for (int iter = 0; iter < max_stress_iterations; iter++) {

      // ============================================================
      // STEP 1: Solve micro-equilibrium at current F_bar
      // ============================================================
      bool plasticity = false;
      UserData userData(square_points, elements, calculator, potential_func,
                        potential_func_der, zero, optimal_lattice_parameter,
                        F_bar, interior_mapping, full_mapping, active_elements,
                        plasticity);

      alglib::real_1d_array x;
      int n_vars_inner = interior_mapping.size();
      x.setlength(2 * n_vars_inner);
      map_points_to_solver_array(x, square_points, interior_mapping,
                                 n_vars_inner);

      userData.third_condition_flag = false;
      LBFGSOptimizer optimizer(13, 0, 0, 0, 0);
      optimizer.optimize(x, minimize_energy_with_triangles, &userData);

      map_solver_array_to_points(x, square_points, interior_mapping,
                                 n_vars_inner);

      // ============================================================
      // STEP 2: Compute current average stress (for convergence)
      // ============================================================
      Eigen::Matrix2d current_stress = Eigen::Matrix2d::Zero();
      double current_energy_iter = 0.0;
      ConfigurationSaver::calculateEnergyAndStress(
          &userData, current_energy_iter, current_stress, true);

      double stress_error = std::abs(current_stress(0, 1) - sigma_app_xy);

      // Update final state variables
      final_stress = current_stress;
      post_energy = current_energy_iter;
      post_area = ConfigurationSaver::calculateTotalArea2D(&userData);

      if (iter % 10 == 0 || stress_error < stress_tolerance) {
        std::cout << "  Iter " << iter << ": σ_xy = " << current_stress(0, 1)
                  << ", error = " << stress_error << std::endl;
      }

      if (stress_error < stress_tolerance) {
        converged = true;
        break;
      }

      // ============================================================
      // STEP 3: Compute spatial tangent (C_macro)
      // ============================================================
      Eigen::Matrix2d C_metric = Eigen::Matrix2d::Identity();
      Eigen::Matrix2d Z = Eigen::Matrix2d::Identity();

      AcousticTensor acoustic_tensor(F_bar, C_metric, Z);
      acoustic_tensor.computeEnergyDerivatives(
          calculator, potential_func_der, potential_func_sder, normalisation);

      itensor::ITensor C_spatial = acoustic_tensor.getSpatialTangent();

      // ============================================================
      // STEP 4: Convert spatial tangent to Voigt notation
      // ============================================================
      Eigen::Matrix3d C_voigt = Eigen::Matrix3d::Zero();

      auto inds = itensor::inds(C_spatial);
      if (inds.size() == 4) {
        itensor::Index i_idx = inds[0];
        itensor::Index j_idx = inds[1];
        itensor::Index k_idx = inds[2];
        itensor::Index l_idx = inds[3];

        // Direct extraction (Cauchy tangent)
        C_voigt(0, 0) =
            itensor::elt(C_spatial, i_idx = 1, j_idx = 1, k_idx = 1, l_idx = 1);
        C_voigt(0, 1) =
            itensor::elt(C_spatial, i_idx = 1, j_idx = 1, k_idx = 2, l_idx = 2);
        C_voigt(0, 2) =
            itensor::elt(C_spatial, i_idx = 1, j_idx = 1, k_idx = 1, l_idx = 2);

        C_voigt(1, 0) =
            itensor::elt(C_spatial, i_idx = 2, j_idx = 2, k_idx = 1, l_idx = 1);
        C_voigt(1, 1) =
            itensor::elt(C_spatial, i_idx = 2, j_idx = 2, k_idx = 2, l_idx = 2);
        C_voigt(1, 2) =
            itensor::elt(C_spatial, i_idx = 2, j_idx = 2, k_idx = 1, l_idx = 2);

        C_voigt(2, 0) =
            itensor::elt(C_spatial, i_idx = 1, j_idx = 2, k_idx = 1, l_idx = 1);
        C_voigt(2, 1) =
            itensor::elt(C_spatial, i_idx = 1, j_idx = 2, k_idx = 2, l_idx = 2);
        C_voigt(2, 2) =
            itensor::elt(C_spatial, i_idx = 1, j_idx = 2, k_idx = 1, l_idx = 2);

      } else {
        C_voigt = Eigen::Matrix3d::Identity();
      }

      // ============================================================
      // STEP 5: Solve C_voigt * Δε_voigt = -(σ̄ - σ_app)
      // ============================================================
      Eigen::Vector3d sigma_current_voigt;
      sigma_current_voigt(0) = current_stress(0, 0);
      sigma_current_voigt(1) = current_stress(1, 1);
      sigma_current_voigt(2) = current_stress(0, 1);
      Eigen::Vector3d sigma_app_voigt;
      sigma_app_voigt(0) = sigma_app(0, 0);
      sigma_app_voigt(1) = sigma_app(1, 1);
      sigma_app_voigt(2) = sigma_app(0, 1);

      // Negative residual: R = sigma_app - sigma_current
      Eigen::Vector3d delta_sigma_voigt = sigma_current_voigt - sigma_app_voigt;
      Eigen::Vector3d delta_eps_voigt;

      double det_C = C_voigt.determinant();
      if (std::abs(det_C) > 1e-12) {
        // Use robust LU solver for full matrix
        delta_eps_voigt = -C_voigt.fullPivLu().solve(delta_sigma_voigt);
      } else {
        // Simplified shear-only solve for near-singular matrix
        double C_1212 = C_voigt(2, 2);
        if (std::abs(C_1212) < 1e-12)
          C_1212 = 1.0;
        delta_eps_voigt.setZero();
        delta_eps_voigt(2) = -delta_sigma_voigt(2) / C_1212;
      }

      // ============================================================
      // STEP 6: Convert delta_eps_voigt to Update F_bar
      // ============================================================
      {
        // Apply relaxation to the strain correction
        Eigen::Vector3d E_correction = delta_eps_voigt;
        E_correction *= alpha_relax;

        // Limit step size for stability (applied to the strain correction)
        double max_dE = 0.01;
        for (int ii = 0; ii < 3; ii++) {
          if (std::abs(E_correction(ii)) > max_dE) {
            E_correction(ii) = (E_correction(ii) > 0) ? max_dE : -max_dE;
          }
        }

        Eigen::Matrix2d F_inc;

        if (use_multiplicative_update) {
          // --- OPTION B: TOTAL MULTIPLICATIVE STRETCH UPDATE (F = sqrt(2E +
          // I)) ---

          // 1. Calculate TOTAL Green-Lagrange Strain E_new from F_bar_old
          // C_old = F_old^T * F_old
          Eigen::Matrix2d C_old = F_bar_old.transpose() * F_bar_old;
          // E_old = 0.5 * (C_old - I)
          Eigen::Matrix2d E_old = 0.5 * (C_old - Eigen::Matrix2d::Identity());

          // Convert E_old to Voigt vector for correction update: E_new = E_old
          // + Correction
          Eigen::Vector3d E_total_voigt;
          E_total_voigt << E_old(0, 0), E_old(1, 1), 2.0 * E_old(0, 1);
          Eigen::Vector3d E_new_voigt = E_total_voigt + E_correction;

          // 2. Reconstruct the TOTAL Right Cauchy-Green Tensor (C = 2E + I)
          double E11 = E_new_voigt(0);
          double E22 = E_new_voigt(1);
          double E12 = E_new_voigt(2) / 2.0; // Extract E12 from 2*E12

          Eigen::Matrix2d C_total;
          C_total << 1.0 + 2.0 * E11, 2.0 * E12, 2.0 * E12, 1.0 + 2.0 * E22;

          // 3. Compute the Total Right Stretch Tensor (U = sqrt(C))
          Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> solver(C_total);
          if (solver.info() != Eigen::Success ||
              solver.eigenvalues().minCoeff() < -1e-10) {
            std::cerr
                << "Eigen decomposition failed or yielded unphysical strain."
                << std::endl;
            converged = false;
            break;
          }

          // U is set as the new F_bar (U = F if R=I)
          F_bar = solver.eigenvectors() *
                  solver.eigenvalues().cwiseSqrt().asDiagonal() *
                  solver.eigenvectors().transpose();

          // 4. Calculate Incremental F for position update: F_inc = F_new *
          // F_old^{-1}
          F_inc = F_bar * F_bar_old.inverse();
        } else {
          // --- OPTION A: STANDARD ADDITIVE INCREMENTAL UPDATE (F_new = F_old +
          // Delta_F) ---

          // 1. Convert E_correction to Delta_F_symmetric
          Eigen::Matrix2d delta_F = Eigen::Matrix2d::Zero();
          delta_F(0, 0) = E_correction(0);       // Δε_11
          delta_F(1, 1) = E_correction(1);       // Δε_22
          delta_F(0, 1) = E_correction(2) / 2.0; // Extract Δε_12
          delta_F(1, 0) = delta_F(0, 1);         // Imposes symmetry

          // 2. Update Total F_bar: F_new = F_old + Delta_F
          F_bar = F_bar + delta_F;

          // 3. Calculate Incremental F for position update: F_inc = I + Delta_F
          F_inc = Eigen::Matrix2d::Identity() + delta_F;
        }

        // ============================================================
        // STEP 7: Update point positions
        // ============================================================

        // F_inc is calculated in both modes and updates the micro-structure.
        for (size_t j = 0; j < square_points.size(); j++) {
          square_points[j].coord = F_inc * square_points[j].coord;
        }
      }

      // Update final F_12 for current iteration's logging
      final_F_12 = F_bar(0, 1);

    } // End of inner stress iteration loop

    // ============================================================
    // POST-CONVERGENCE UPDATE AND LOGGING
    // ============================================================

    if (!converged) {
      std::cout << "  ⚠ Did not converge after " << max_stress_iterations
                << " iterations" << std::endl;
    } else {
      // F_bar_old stores the converged state for the next outer step
      F_bar_old = F_bar;
    }

    // --- LOGGING: SAVE FINAL CONVERGED STATE ---
    int log_iteration = caller_id * 1000 + stress_idx;
    double log_applied_stress = sigma_app_xy;
    bool shouldRemesh = false;

    // Log: (iteration, applied_stress, pre_energy, pre_F_12, post_energy,
    // final_F_12, ...)
    ConfigurationSaver::logEnergyAndStress_v2(
        log_iteration, log_applied_stress, pre_energy, pre_F_12, post_energy,
        final_F_12, pre_area, post_area, shouldRemesh);

    // Set post_energy_previous for the next load step's pre-state
    post_energy_previous = post_energy;

    std::cout << "Logged data for iteration " << log_iteration
              << " (Applied σ_xy=" << log_applied_stress
              << ", Measured F_12=" << final_F_12 << ")" << std::endl;
    // ---------------------------------------------------------------------------------------

    // ==================== SAVE CONFIGURATION ====================
    int file_id = caller_id + file_counter;

    UserData finalUserData(square_points, elements, calculator, potential_func,
                           potential_func_der, zero, optimal_lattice_parameter,
                           F_bar, interior_mapping, full_mapping,
                           active_elements, false);

    double final_energy_out = 0.0;
    Eigen::Matrix2d final_stress_out = Eigen::Matrix2d::Zero();
    ConfigurationSaver::calculateEnergyAndStress(
        &finalUserData, final_energy_out, final_stress_out, true);

    auto [num_dislocations, coordination] =
        DefectAnalysis::analyzeDefectsInReferenceConfig(
            &finalUserData, file_id, dndx, offsets, original_domain_map,
            translation_map, domain_dims_point, element_area, pbc, true);

    ConfigurationSaver::writeToVTK(finalUserData.points, finalUserData.elements,
                                   &finalUserData, file_id, true, coordination,
                                   sigma_app_xy);

    std::cout << "Saved configuration " << file_id
              << " at target σ_xy = " << sigma_app_xy
              << ", measured F_12 = " << final_F_12 << std::endl;

    file_counter++;
  }

  std::cout << "\n=== STRESS-CONTROLLED SIMULATION COMPLETE ===" << std::endl;
}
