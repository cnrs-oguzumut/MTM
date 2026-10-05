#include "../../include/experiments/experiment_includes.h"

void parametricAcousticStudy() {
  std::cout
      << "Starting parametric acoustic tensor study with Lagrange reduction..."
      << std::endl;

  // Create output file
  std::ofstream file("acoustic_study_results.dat");
  file << std::scientific << std::setprecision(8);
  file << "# t p c11 c22 c12 c11_red c22_red c12_red min_detAc angle_deg "
          "third_condition"
       << std::endl;

  // Create strain energy calculator once (outside the loop)
  // double scale = 1.0;
  // Strain_Energy_LatticeCalculator strain_calculator(scale);
  // double normalisation = strain_calculator.getUnitCellArea();
  double gamma = pow(4. / 3., 1. / 4.);
  Eigen::Matrix2d H;
  H << 1.0, 0.5, 0.0, std::sqrt(3.0) / 2.0;

  //  H << 0.5,            -0.5,
  //              std::sqrt(3.0)/2.0,  std::sqrt(3.0)/2.0;

  // H = gamma*H;
  H.setIdentity();
  // H=H.transpose();

  // H << gamma * std::sqrt(2.0 + std::sqrt(3.0)) / 2.0,  gamma * std::sqrt(2.0
  // - std::sqrt(3.0)) / 2.0,
  //      gamma * std::sqrt(2.0 - std::sqrt(3.0)) / 2.0,  gamma * std::sqrt(2.0
  //      + std::sqrt(3.0)) / 2.0; H=H.transpose();
  //     //  H.setIdentity();

  // double scale = 0.687204444204349/gamma;

  int mode = 3;

  double scale;
  double normalisation;
  double r_cutoff;
  std::function<double(double)> potential_func;
  std::function<double(double)> potential_func_der;
  std::function<double(double)> potential_func_sder;

  gamma = 1.0;

  if (mode == 1) {

    scale = 0.6872044091828517; // 0.687204444204349 ;
    r_cutoff = 2.5;

    potential_func = lennard_jones_energy_v2;
    potential_func_der = lennard_jones_energy_der_v2;
    potential_func_sder = lennard_jones_energy_sder_v2;

  }

  else if (mode == 2) {

    scale = 0.996407146941421;
    r_cutoff = 1.86602540378444;

    potential_func = lennard_jones_energy_v3;
    potential_func_der = lennard_jones_energy_der_v3;
    potential_func_sder = lennard_jones_energy_sder_v3;

  }

  else if (mode == 3) {
    // dummies
    scale = 1.;
    r_cutoff = 1.86602540378444;

    potential_func = [](double r) -> double { return 1.0; };
    potential_func_der = [](double r) -> double { return 1.0; };
    potential_func_sder = [](double r) -> double { return 1.0; };
  }

  // SquareLatticeCalculator strain_calculator(scale,r_cutoff);
  // //TriangularLatticeCalculator strain_calculator(scale,r_cutoff);
  // normalisation = strain_calculator.getUnitCellArea();
  // //since I use triangular lattice vectors
  // normalisation *= sqrt(3.)/2.;
  // std::cout << "normalisation: " << normalisation << std::endl;

  Strain_Energy_LatticeCalculator strain_calculator(scale);
  normalisation = strain_calculator.getUnitCellArea();

  // Define dummy potential functions once (ignored by strain energy calculator)
  // auto dummy_dpot = [](double r) -> double { return 1.0; };
  // auto dummy_d2pot = [](double r) -> double { return 0.1; };

  int count = 0;

  Eigen::Matrix2d C_ref;
  Eigen::Matrix2d F_ref;
  Eigen::Matrix2d Z_ref;
  C_ref.setIdentity();
  F_ref.setIdentity();
  Z_ref.setIdentity();

  F_ref.setIdentity();
  C_ref = F_ref.transpose() * F_ref;

  AcousticTensor acoustic_tensor(F_ref, C_ref, Z_ref);
  acoustic_tensor.computeEnergyDerivatives(strain_calculator,
                                           potential_func_der,
                                           potential_func_sder, normalisation);
  // Cxxxx = Cyyyy = 37.2120020466688
  // Cxxyy = Cxxyy = 12.4040009178261
  acoustic_tensor.printHessianComponents();

  // Main parametric loop - same structure as your working code but in a loop
  for (double t = 0.0; t <= 1.; t += 0.01) {             // radius
    for (double p = -M_PI; p <= M_PI; p += M_PI / 128) { // angle

      try {
        // Compute C matrix components exactly as specified
        double c11 = cosh(t) + sinh(t) * sin(p);
        double c22 = cosh(t) - sinh(t) * sin(p);
        double c12 = sinh(t) * cos(p);
        if (c12 < 0 && mode != 3)
          continue; // only need half the space due to symmetry

        // Construct original C matrix
        Eigen::Matrix2d C_original;
        C_original << c11, c12, c12, c22;
        // Eigen::Matrix2d C_original = H.transpose() * C_original * H;

        // Check if C is positive definite
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver(C_original);
        if (eigensolver.eigenvalues().minCoeff() <= 1e-10) {
          continue; // Skip non-positive definite matrices
        }

        // Apply Lagrange reduction to get Z matrix and reduced C
        lagrange::Result reduction_result = lagrange::reduce(C_original);

        Eigen::Matrix2d C, Z;

        if (mode == 1 || mode == 2) {
          // use original C and identity Z
          // C = C_original;
          // Z = Eigen::Matrix2d::Identity();
          C = reduction_result.C_reduced;
          Z = reduction_result.m_matrix;

        } else {
          // For mode 3 (or any other mode), Use reduced C and Z from reduction
          // result
          C = reduction_result.C_reduced;
          Z = reduction_result.m_matrix;
        }

        bool third_condition = reduction_result.third_condition_satisfied;

        // Polar decomposition: find F such that C = F^T * F
        // Using your working approach
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> C_eigensolver(
            C_original);
        if (C_eigensolver.eigenvalues().minCoeff() <= 1e-10) {
          continue;
        }

        Eigen::Matrix2d S_sqrt =
            C_eigensolver.eigenvalues().cwiseSqrt().asDiagonal();
        // This F is indeed U such that F=QU
        Eigen::Matrix2d F = C_eigensolver.eigenvectors() * S_sqrt *
                            C_eigensolver.eigenvectors().transpose();

        // Verify: C should equal F^T * F
        Eigen::Matrix2d C_check = F.transpose() * F;
        if ((C_original - C_check).norm() > 1e-10) {
          std::cout << "Warning: Polar decomposition check failed at t=" << t
                    << ", p=" << p << std::endl;
        }

        // std::cout << "\n=== DEBUG: Point " << count+1 << " at t=" << t << ",
        // p=" << p << " ===" << std::endl; std::cout << "Original C matrix:\n"
        // << C_original << std::endl; std::cout << "Reduced C matrix:\n" << C
        // << std::endl; std::cout << "F matrix:\n" << F << std::endl; std::cout
        // << "Z matrix (from reduction):\n" << Z << std::endl; std::cout <<
        // "F^T * F check:\n" << C_check << std::endl; std::cout << "C - F^T*F
        // norm: " << (C - C_check).norm() << std::endl; std::cout << "C
        // eigenvalues: " << C_eigensolver.eigenvalues().transpose() <<
        // std::endl; std::cout << "C determinant: " << C.determinant() <<
        // std::endl; std::cout << "Z determinant: " << Z.determinant() <<
        // std::endl; std::cout << "About to create AcousticTensor..." <<
        // std::endl;

        // Create acoustic tensor object (same as your working code)
        AcousticTensor acoustic_tensor(F, C, Z);

        // Compute energy derivatives using the uniform interface (same as
        // working code)
        acoustic_tensor.computeEnergyDerivatives(
            strain_calculator, potential_func_der, potential_func_sder,
            normalisation);

        // Perform acoustic tensor analysis (same as working code)
        bool lagrangian = false;
        AcousticAnalysis result =
            acoustic_tensor.analyzeAcousticTensor(lagrangian);
        Eigen::Matrix2d C_original_new = H.transpose() * C_original * H;

        // Write results to file
        file << t << " " << p << " " << C_original_new(0, 0)
             << " " // Original C components
             << C_original_new(1, 1) << " " << C_original_new(0, 1) << " "
             << C(0, 0) << " " // Reduced C components
             << C(1, 1) << " " << C(0, 1) << " " << result.detAc << " "
             << result.xsi << " " << (third_condition ? 1 : 0) << std::endl;

        count++;

        // Progress update
        if (count % 100 == 0) {
          std::cout << "Processed " << count
                    << " points successfully. Current: t=" << std::fixed
                    << std::setprecision(3) << t << ", p=" << p << std::endl;
        }

      } catch (const std::exception &e) {
        std::cout << "Error at t=" << t << ", p=" << p << ": " << e.what()
                  << std::endl;
        continue;
      }
    }
  }

  file.close();
  std::cout << "Parametric study completed successfully!" << std::endl;
  std::cout << "Total points processed: " << count << std::endl;
  std::cout << "Results saved to 'acoustic_study_results.dat'" << std::endl;
}

// Compute F = κI + γR(θ)e₁⊗e₂
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

// Compute F = κI + γR(θ)[e₁ ⊗ e₂]
Eigen::Matrix2d compute_F(double kappa, double gamma, double theta) {
  Eigen::Matrix2d I = Eigen::Matrix2d::Identity();

  double cos_theta = std::cos(theta);
  double sin_theta = std::sin(theta);

  // Shear direction
  Eigen::Vector2d e(cos_theta, sin_theta);

  // Normal direction (perpendicular)
  Eigen::Vector2d n(-sin_theta, cos_theta);

  // Outer product e ⊗ n
  Eigen::Matrix2d e_outer_n = e * n.transpose();

  // F = I + γ(e ⊗ n)
  Eigen::Matrix2d F = I + gamma * e_outer_n;

  return F;
}

void parametricAcousticStudy_v2() {
  std::cout
      << "Starting parametric acoustic tensor study with Lagrange reduction..."
      << std::endl;

  // Create output file
  std::ofstream file("acoustic_study_results.dat");
  file << std::scientific << std::setprecision(8);
  file << "# t p c11 c22 c12 c11_red c22_red c12_red min_detAc angle_deg "
          "third_condition"
       << std::endl;

  const double kappa = 1.0;

  // Define ranges
  const int n_theta = 400;
  const int n_gamma = 400;

  std::vector<double> theta_values;
  std::vector<double> gamma_values;

  // Theta from 0 to 2π
  for (int i = 0; i < n_theta; ++i) {
    theta_values.push_back(2.0 * M_PI * i / n_theta);
  }

  // Gamma from 0 to a
  double a = 1.0;
  for (int i = 0; i < n_gamma; ++i) {
    gamma_values.push_back(a * i / (n_gamma - 1));
  }

  double gamma = pow(4. / 3., 1. / 4.);
  Eigen::Matrix2d H;
  H.setIdentity();
  int mode = 1;

  double scale;
  double normalisation;
  double r_cutoff;
  std::function<double(double)> potential_func;
  std::function<double(double)> potential_func_der;
  std::function<double(double)> potential_func_sder;

  gamma = 1.0;

  if (mode == 1) {

    scale = 0.6872044091828517; // 0.687204444204349 ;
    r_cutoff = 2.5;

    potential_func = lennard_jones_energy_v2;
    potential_func_der = lennard_jones_energy_der_v2;
    potential_func_sder = lennard_jones_energy_sder_v2;

  }

  else if (mode == 2) {

    scale = 0.996407146941421;
    r_cutoff = 1.86602540378444;

    potential_func = lennard_jones_energy_v3;
    potential_func_der = lennard_jones_energy_der_v3;
    potential_func_sder = lennard_jones_energy_sder_v3;

  }

  else if (mode == 3) {
    // dummies
    scale = 1.;
    r_cutoff = 1.86602540378444;

    potential_func = [](double r) -> double { return 1.0; };
    potential_func_der = [](double r) -> double { return 1.0; };
    potential_func_sder = [](double r) -> double { return 1.0; };
  }

  SquareLatticeCalculator strain_calculator(scale, r_cutoff);
  // TriangularLatticeCalculator strain_calculator(scale,r_cutoff);
  normalisation = strain_calculator.getUnitCellArea();
  // since I use triangular lattice vectors
  normalisation *= sqrt(3.) / 2.;
  std::cout << "normalisation: " << normalisation << std::endl;

  // Strain_Energy_LatticeCalculator strain_calculator(scale);
  // normalisation = strain_calculator.getUnitCellArea();

  // Define dummy potential functions once (ignored by strain energy calculator)
  // auto dummy_dpot = [](double r) -> double { return 1.0; };
  // auto dummy_d2pot = [](double r) -> double { return 0.1; };

  int count = 0;

  Eigen::Matrix2d C_ref;
  Eigen::Matrix2d F_ref;
  Eigen::Matrix2d Z_ref;
  C_ref.setIdentity();
  F_ref.setIdentity();
  Z_ref.setIdentity();

  F_ref.setIdentity();
  C_ref = F_ref.transpose() * F_ref;

  AcousticTensor acoustic_tensor(F_ref, C_ref, Z_ref);
  acoustic_tensor.computeEnergyDerivatives(strain_calculator,
                                           potential_func_der,
                                           potential_func_sder, normalisation);
  // Cxxxx = Cyyyy = 37.2120020466688
  // Cxxyy = Cxxyy = 12.4040009178261
  acoustic_tensor.printHessianComponents();

  // Main parametric loop - same structure as your working code but in a loop
  // Loop over all combinations
  std::cout << "Computing F matrices for κ = " << kappa << std::endl;
  std::cout << "Theta range: [0, 2π] with " << n_theta << " points"
            << std::endl;
  std::cout << "Gamma range: [0, 2] with " << n_gamma << " points" << std::endl;
  std::cout << "Total: " << n_theta * n_gamma << " combinations\n" << std::endl;

  for (double t : gamma_values) {
    for (double p : theta_values) {
      // Compute F
      Eigen::Matrix2d F = compute_F(1, t, p);

      try {
        // Compute C matrix components exactly as specified
        Eigen::Matrix2d C_original = F.transpose() * F;

        // Apply Lagrange reduction to get Z matrix and reduced C
        lagrange::Result reduction_result = lagrange::reduce(C_original);

        Eigen::Matrix2d C, Z;

        if (mode == 1 || mode == 2) {
          // use original C and identity Z
          // C = C_original;
          // Z = Eigen::Matrix2d::Identity();
          C = reduction_result.C_reduced;
          Z = reduction_result.m_matrix;

        } else {
          // For mode 3 (or any other mode), Use reduced C and Z from reduction
          // result
          C = reduction_result.C_reduced;
          Z = reduction_result.m_matrix;
        }

        bool third_condition = reduction_result.third_condition_satisfied;

        // Create acoustic tensor object (same as your working code)
        AcousticTensor acoustic_tensor(F, C, Z);

        // Compute energy derivatives using the uniform interface (same as
        // working code)
        acoustic_tensor.computeEnergyDerivatives(
            strain_calculator, potential_func_der, potential_func_sder,
            normalisation);

        // Perform acoustic tensor analysis (same as working code)
        bool lagrangian = false;
        AcousticAnalysis result =
            acoustic_tensor.analyzeAcousticTensor(lagrangian);
        Eigen::Matrix2d C_original_new = H.transpose() * C_original * H;

        // Write results to file
        file << t << " " << p << " " << C_original_new(0, 0)
             << " " // Original C components
             << C_original_new(1, 1) << " " << C_original_new(0, 1) << " "
             << C(0, 0) << " " // Reduced C components
             << C(1, 1) << " " << C(0, 1) << " " << result.detAc << " "
             << result.xsi << " " << (third_condition ? 1 : 0) << std::endl;

        count++;

        // Progress update
        if (count % 100 == 0) {
          std::cout << "Processed " << count
                    << " points successfully. Current: t=" << std::fixed
                    << std::setprecision(3) << t << ", p=" << p
                    << ", det F=" << F.determinant() << std::endl;
        }

      } catch (const std::exception &e) {
        std::cout << "Error at t=" << t << ", p=" << p << ": " << e.what()
                  << std::endl;
        continue;
      }
    }
  }

  file.close();
  std::cout << "Parametric study completed successfully!" << std::endl;
  std::cout << "Total points processed: " << count << std::endl;
  std::cout << "Results saved to 'acoustic_study_results.dat'" << std::endl;
}

// void example_3_stress_controlled_with_tangent(int caller_id, int nx, int ny)
// {

//   // Stress-controlled loading using proper spatial tangent from
//   AcousticTensor
//   //
//   // PHASE 1: Strain-controlled equilibration at F_bar = [1, 0.14; 0, 1]
//   // PHASE 2: Stress-controlled loading starting from equilibrated state

//   if (nx <= 0 || ny <= 0) {
//     std::cerr << "Error: nx and ny must be positive integers." << std::endl;
//     exit(EXIT_FAILURE);
//   }

//   // ==================== SETUP REFERENCE GEOMETRY ====================
//   writeSizesToFile(nx, ny);

//   std::string lattice_type = "square";
//   double h = 1.0;

//   Eigen::Vector2d p1(0.0, 0.0);
//   Eigen::Vector2d p2(h, 0.0);
//   Eigen::Vector2d p3(0.0, h);

//   Eigen::Matrix<double, 3, 2> dndx = calculateShapeDerivatives(p1, p2, p3);

//   // ==================== SETUP ENERGY POTENTIAL ====================
//   std::function<double(double)> potential_func = square_energy;
//   std::function<double(double)> potential_func_der = square_energy_der;
//   std::function<double(double)> potential_func_sder = square_energy_der;  //
//   Second derivative

//   // ==================== LATTICE PARAMETER ====================
//   double symmetry_constantx =
//       (lattice_type == "triangular") ? pow(4.0 / 3.0, 1.0 / 4.0) : 1.0;
//   double optimal_lattice_parameter = symmetry_constantx * 1.0;
//   double lattice_constant = optimal_lattice_parameter;

//   // ==================== GENERATE INITIAL LATTICE ====================
//   std::vector<Point2D> square_points = LatticeGenerator::generate_2d_lattice(
//       nx, ny, lattice_constant, lattice_type);
//   std::vector<Point2D> square_points_ref =
//       LatticeGenerator::generate_2d_lattice(nx, ny, lattice_constant,
//                                             lattice_type);

//   int original_domain_size = square_points.size();

//   DomainInfo domain_size = compute_domain_size(square_points);

//   const std::array<double, 2> offsets =
//       (lattice_type == "square")
//           ? std::array<double, 2>{lattice_constant, lattice_constant}
//           : std::array<double, 2>{lattice_constant / 2.0,
//                                   (sqrt(3.0) / 2.0) * lattice_constant};

//   DomainDimensions domain_dims(domain_size.get_width(),
//                                domain_size.get_height());
//   bool pbc = true;

//   auto [original_domain_map, translation_map] =
//       MeshGenerator::create_domain_maps(original_domain_size, domain_dims,
//                                         offsets);

//   auto [interior_mapping, full_mapping] =
//       create_dof_mapping_original(square_points, 0.5 * lattice_constant,
//       pbc);

//   Point2D domain_dims_point(domain_dims.size_x, domain_dims.size_y);

//   // ==================== CREATE MESH ====================
//   AdaptiveMesher mesher(domain_dims_point, offsets, original_domain_map,
//                         translation_map, full_mapping, 1e-6, pbc);
//   mesher.setUsePeriodicCopies(pbc);

//   alglib::real_1d_array free_dofs;
//   int n_free_nodes = interior_mapping.size();
//   free_dofs.setlength(2 * interior_mapping.size());
//   map_points_to_solver_array(free_dofs, square_points, interior_mapping,
//                              n_free_nodes);

//   alglib::real_1d_array original_x_remesh =
//       mesher.saveOriginalPositions(free_dofs);

//   auto [elements, active_elements] = mesher.createMesh(
//       square_points, free_dofs, Eigen::Matrix2d::Identity(), &dndx);
//   double element_area = elements[0].getReferenceArea();

//   for (auto &element : elements) {
//     element.set_dof_mapping(full_mapping);
//   }

//   // ==================== SETUP ENERGY CALCULATION ====================
//   Strain_Energy_LatticeCalculator calculator(1.0);
//   double normalisation = calculator.getUnitCellArea();

//   Eigen::Matrix2d F_I = Eigen::Matrix2d::Identity();
//   Eigen::Matrix2d C_I = F_I.transpose() * F_I;
//   double zero = calculator.calculate_energy(C_I, potential_func, 0);

//   // ================================================================
//   // PHASE 1: STRAIN-CONTROLLED EQUILIBRATION
//   // Start at F_bar = [1, 0.14; 0, 1] and minimize energy
//   // ================================================================
//   std::cout <<
//   "\n================================================================" <<
//   std::endl; std::cout << "PHASE 1: STRAIN-CONTROLLED EQUILIBRATION" <<
//   std::endl; std::cout <<
//   "================================================================" <<
//   std::endl;

//   // Set initial deformation gradient
//   Eigen::Matrix2d F_bar = Eigen::Matrix2d::Identity();
//   double initial_shear = 0.14;
//   F_bar(0, 1) = initial_shear;

//   std::cout << "Initial F_bar:\n" << F_bar << std::endl;

//   // Apply F_bar to reference positions + noise
//   {
//     std::random_device rd;
//     std::mt19937 gen(rd());
//     double noise_level = 0.04;
//     std::normal_distribution<double> noise_dist(0.0, noise_level);

//     for (size_t j = 0; j < square_points.size(); j++) {
//       Eigen::Vector2d noise(noise_dist(gen), noise_dist(gen));
//       // Apply F_bar to reference positions, then add noise
//       square_points[j].coord = F_bar * square_points_ref[j].coord + noise;
//     }
//   }

//   // Perform strain-controlled minimization
//   double initial_stress_xy = 0.0;
//   {
//     bool plasticity = false;
//     UserData userData(square_points, elements, calculator, potential_func,
//                       potential_func_der, zero, optimal_lattice_parameter,
//                       F_bar, interior_mapping, full_mapping, active_elements,
//                       plasticity);

//     alglib::real_1d_array x;
//     int n_vars = interior_mapping.size();
//     x.setlength(2 * n_vars);
//     map_points_to_solver_array(x, square_points, interior_mapping, n_vars);

//     std::cout << "Running strain-controlled minimization at F_12 = " <<
//     initial_shear << "..." << std::endl;

//     userData.third_condition_flag = false;
//     LBFGSOptimizer optimizer(13, 0, 0, 0, 0);
//     optimizer.optimize(x, minimize_energy_with_triangles, &userData);

//     map_solver_array_to_points(x, square_points, interior_mapping, n_vars);

//     // Compute equilibrated stress
//     double equil_energy = 0.0;
//     Eigen::Matrix2d equil_stress = Eigen::Matrix2d::Zero();
//     ConfigurationSaver::calculateEnergyAndStress(&userData, equil_energy,
//     equil_stress, true);

//     initial_stress_xy = equil_stress(0, 1);

//     std::cout << "Equilibration complete!" << std::endl;
//     std::cout << "  Energy = " << equil_energy << std::endl;
//     std::cout << "  Stress σ_xy = " << initial_stress_xy << std::endl;
//     std::cout << "  Full stress tensor:\n" << equil_stress << std::endl;
//     std::cout << "  F_bar:\n" << F_bar << std::endl;
//   }

//   // ================================================================
//   // PHASE 2: STRESS-CONTROLLED LOADING
//   // Use equilibrated state as initial condition
//   // ================================================================
//   std::cout <<
//   "\n================================================================" <<
//   std::endl; std::cout << "PHASE 2: STRESS-CONTROLLED LOADING" << std::endl;
//   std::cout <<
//   "================================================================" <<
//   std::endl;

//   // ==================== STRESS LOADING SCHEDULE ====================
//   double sigma_app_min = 0.0;
//   double sigma_app_max = 1.5;
//   double sigma_app_step = 0.01;

//   int num_stress_points =
//       static_cast<int>((sigma_app_max - sigma_app_min) / sigma_app_step) + 1;

//   std::vector<double> sigma_app_values;
//   sigma_app_values.reserve(num_stress_points);
//   for (int i = 0; i < num_stress_points; i++) {
//     sigma_app_values.push_back(sigma_app_min + i * sigma_app_step);
//   }

//   std::cout << "Stress control: " << num_stress_points << " steps from "
//             << sigma_app_min << " to " << sigma_app_max << std::endl;
//   std::cout << "Starting from equilibrated state with σ_xy = " <<
//   initial_stress_xy
//             << ", F_12 = " << F_bar(0,1) << std::endl;

//   // ==================== STRESS CONTROL PARAMETERS ====================
//   int max_stress_iterations = 100;
//   double stress_tolerance = 1e-5;
//   double alpha_relax = 0.1;

//   static int file_counter = 0;
//   static double post_energy_previous = 0.0;

//   // ==================== STRESS CONTROL LOOP ====================
//   for (size_t stress_idx = 0; stress_idx < sigma_app_values.size();
//   stress_idx++) {
//     double sigma_app_xy = sigma_app_values[stress_idx];

//     Eigen::Matrix2d sigma_app = Eigen::Matrix2d::Zero();
//     sigma_app(0, 1) = sigma_app_xy;
//     sigma_app(1, 0) = sigma_app_xy;

//     std::cout << "\n========================================" << std::endl;
//     std::cout << "TARGET STRESS σ_xy = " << sigma_app_xy << std::endl;
//     std::cout << "========================================" << std::endl;

//     std::vector<Point2D> points_at_load_start = square_points;
//     Eigen::Matrix2d F_bar_at_load_start = F_bar;

//     bool converged = false;

//     for (int iter = 0; iter < max_stress_iterations; iter++) {

//       // ============================================================
//       // STEP 1: Solve micro-equilibrium at current F_bar
//       // ============================================================
//       bool plasticity = false;
//       UserData userData(square_points, elements, calculator, potential_func,
//                         potential_func_der, zero, optimal_lattice_parameter,
//                         F_bar, interior_mapping, full_mapping,
//                         active_elements, plasticity);

//       alglib::real_1d_array x;
//       int n_vars = interior_mapping.size();
//       x.setlength(2 * n_vars);
//       map_points_to_solver_array(x, square_points, interior_mapping, n_vars);

//       userData.third_condition_flag = false;
//       LBFGSOptimizer optimizer(13, 0, 0, 0, 0);
//       optimizer.optimize(x, minimize_energy_with_triangles, &userData);

//       map_solver_array_to_points(x, square_points, interior_mapping, n_vars);

//       // ============================================================
//       // STEP 2: Compute current average stress
//       // ============================================================
//       double current_energy = 0.0;
//       Eigen::Matrix2d current_stress = Eigen::Matrix2d::Zero();
//       ConfigurationSaver::calculateEnergyAndStress(&userData, current_energy,
//                                                    current_stress, true);

//       double stress_error = std::abs(current_stress(0, 1) - sigma_app_xy);

//       if (iter % 10 == 0 || stress_error < stress_tolerance) {
//         std::cout << "  Iter " << iter
//                   << ": σ_xy = " << current_stress(0, 1)
//                   << ", error = " << stress_error
//                   << std::endl;
//         std::cout << "    F_bar = [" << F_bar(0,0) << ", " << F_bar(0,1) <<
//         "; "
//                   << F_bar(1,0) << ", " << F_bar(1,1) << "]" << std::endl;
//       }

//       if (stress_error < stress_tolerance) {
//         std::cout << "  ✓ Converged in " << iter << " iterations" <<
//         std::endl; converged = true; break;
//       }

//       // ============================================================
//       // STEP 3: Compute spatial tangent using AcousticTensor
//       // ============================================================
//       Eigen::Matrix2d C_metric = Eigen::Matrix2d::Identity();
//       Eigen::Matrix2d Z = Eigen::Matrix2d::Identity();

//       AcousticTensor acoustic_tensor(F_bar, C_metric, Z);
//       acoustic_tensor.computeEnergyDerivatives(
//           calculator,
//           potential_func_der,
//           potential_func_sder,
//           normalisation);

//       itensor::ITensor C_spatial = acoustic_tensor.getSpatialTangent();

//       // ============================================================
//       // STEP 4: Convert spatial tangent to Voigt notation (3x3 in 2D)
//       // NO SYMMETRIZATION - direct extraction only
//       // ============================================================
//       Eigen::Matrix3d C_voigt = Eigen::Matrix3d::Zero();

//       auto inds = itensor::inds(C_spatial);
//       if (inds.size() == 4) {
//         itensor::Index i_idx = inds[0];
//         itensor::Index j_idx = inds[1];
//         itensor::Index k_idx = inds[2];
//         itensor::Index l_idx = inds[3];

//         // Direct extraction - NO symmetrization factors!
//         C_voigt(0, 0) = itensor::elt(C_spatial, i_idx=1, j_idx=1, k_idx=1,
//         l_idx=1); C_voigt(0, 1) = itensor::elt(C_spatial, i_idx=1, j_idx=1,
//         k_idx=2, l_idx=2); C_voigt(0, 2) = itensor::elt(C_spatial, i_idx=1,
//         j_idx=1, k_idx=1, l_idx=2);

//         C_voigt(1, 0) = itensor::elt(C_spatial, i_idx=2, j_idx=2, k_idx=1,
//         l_idx=1); C_voigt(1, 1) = itensor::elt(C_spatial, i_idx=2, j_idx=2,
//         k_idx=2, l_idx=2); C_voigt(1, 2) = itensor::elt(C_spatial, i_idx=2,
//         j_idx=2, k_idx=1, l_idx=2);

//         C_voigt(2, 0) = itensor::elt(C_spatial, i_idx=1, j_idx=2, k_idx=1,
//         l_idx=1); C_voigt(2, 1) = itensor::elt(C_spatial, i_idx=1, j_idx=2,
//         k_idx=2, l_idx=2); C_voigt(2, 2) = itensor::elt(C_spatial, i_idx=1,
//         j_idx=2, k_idx=1, l_idx=2);

//       } else {
//         std::cerr << "  Error: Expected 4th order tensor, got " <<
//         inds.size() << " indices" << std::endl; C_voigt =
//         Eigen::Matrix3d::Identity();
//       }

//       if (iter % 20 == 0) {
//         std::cout << "    C_voigt:\n" << C_voigt << std::endl;
//       }

//       // ============================================================
//       // STEP 5: Solve C_voigt * Δε_voigt = -(σ̄ - σ_app)
//       // ============================================================
//       Eigen::Vector3d sigma_current_voigt;
//       sigma_current_voigt(0) = current_stress(0, 0);
//       sigma_current_voigt(1) = current_stress(1, 1);
//       sigma_current_voigt(2) = current_stress(0, 1);

//       Eigen::Vector3d sigma_app_voigt;
//       sigma_app_voigt(0) = sigma_app(0, 0);
//       sigma_app_voigt(1) = sigma_app(1, 1);
//       sigma_app_voigt(2) = sigma_app(0, 1);

//       Eigen::Vector3d delta_sigma_voigt = sigma_current_voigt -
//       sigma_app_voigt;

//       Eigen::Vector3d delta_eps_voigt;

//       double det_C = C_voigt.determinant();
//       if (std::abs(det_C) > 1e-12) {
//         delta_eps_voigt = -C_voigt.fullPivLu().solve(delta_sigma_voigt);
//       } else {
//         std::cout << "  Warning: Singular C_voigt (det = " << det_C << "),
//         using shear-only" << std::endl; double C_1212 = C_voigt(2, 2); if
//         (std::abs(C_1212) < 1e-12) C_1212 = 1.0; delta_eps_voigt.setZero();
//         delta_eps_voigt(2) = -delta_sigma_voigt(2) / C_1212;
//       }

//       // ============================================================
//       // STEP 6: Convert Δε_voigt to ΔF
//       // ============================================================
//       Eigen::Matrix2d delta_F = Eigen::Matrix2d::Zero();
//       delta_F(0, 0) = delta_eps_voigt(0);
//       delta_F(1, 1) = delta_eps_voigt(1);
//       delta_F(0, 1) = delta_eps_voigt(2);
//       delta_F(1, 0) = 0.0;

//       // Apply relaxation
//       delta_F *= alpha_relax;

//       // Limit step size for stability
//       double max_dF = 0.01;
//       for (int ii = 0; ii < 2; ii++) {
//         for (int jj = 0; jj < 2; jj++) {
//           if (std::abs(delta_F(ii, jj)) > max_dF) {
//             delta_F(ii, jj) = (delta_F(ii, jj) > 0) ? max_dF : -max_dF;
//           }
//         }
//       }

//       // ============================================================
//       // STEP 7: Update F_bar and point positions
//       // ============================================================
//       F_bar = F_bar + delta_F;

//       Eigen::Matrix2d F_inc = Eigen::Matrix2d::Identity();
//       F_inc(0, 0) = 1.0 + delta_F(0, 0);
//       F_inc(1, 1) = 1.0 + delta_F(1, 1);
//       F_inc(0, 1) = delta_F(0, 1);
//       F_inc(1, 0) = delta_F(1, 0);

//       for (size_t j = 0; j < square_points.size(); j++) {
//         square_points[j].coord = F_inc * square_points[j].coord;
//       }

//       if (iter % 10 == 0) {
//         std::cout << "    C_1212 = " << C_voigt(2, 2)
//                   << ", ΔF_12 = " << delta_F(0, 1)
//                   << ", det(C_voigt) = " << det_C << std::endl;
//       }
//     }

//     if (!converged) {
//       std::cout << "  ⚠ Did not converge after " << max_stress_iterations
//                 << " iterations" << std::endl;
//     }

//     // ==================== PRINT FINAL F_bar ====================
//     std::cout << "  Final F_bar:\n" << F_bar << std::endl;

//     // ==================== NEIGHBOR CONNECTIVITY CHECK ====================
//     NeighborAnalyzer analyzer(NeighborAnalyzer::SearchType::K_NEAREST);
//     analyzer.setKNearestSearch(4);

//     auto neighbors_before = analyzer.buildNeighbors(points_at_load_start);
//     auto neighbors_after = analyzer.buildNeighbors(square_points);
//     auto change_info = NeighborAnalyzer::compareNeighborsWithTolerance(
//         points_at_load_start, square_points, neighbors_before,
//         neighbors_after, 0.135);

//     if (change_info.has_changed) {
//       std::cout << "✗ Neighbor connectivity CHANGED" << std::endl;
//     } else {
//       std::cout << "✓ Neighbor connectivity UNCHANGED" << std::endl;
//     }

//     // ==================== FINAL ENERGY/STRESS ====================
//     UserData finalUserData(square_points, elements, calculator,
//     potential_func,
//                            potential_func_der, zero,
//                            optimal_lattice_parameter, F_bar,
//                            interior_mapping, full_mapping, active_elements,
//                            false);

//     double post_energy = 0.0;
//     Eigen::Matrix2d final_stress = Eigen::Matrix2d::Zero();
//     ConfigurationSaver::calculateEnergyAndStress(&finalUserData, post_energy,
//                                                  final_stress, true);

//     std::cout << "Final state: F_12 = " << F_bar(0, 1)
//               << ", σ_xy = " << final_stress(0, 1)
//               << ", Energy = " << post_energy << std::endl;

//     // ==================== SAVE CONFIGURATION ====================
//     int file_id = caller_id + file_counter;

//     auto [num_dislocations, coordination] =
//     DefectAnalysis::analyzeDefectsInReferenceConfig(
//         &finalUserData, file_id, dndx, offsets, original_domain_map,
//         translation_map, domain_dims_point, element_area, pbc, true);

//     ConfigurationSaver::writeToVTK(
//         finalUserData.points, finalUserData.elements, &finalUserData,
//         file_id, true, coordination, sigma_app_xy);

//     std::cout << "Saved configuration " << file_id
//               << " at target σ_xy = " << sigma_app_xy << std::endl;

//     file_counter++;
//     post_energy_previous = post_energy;
//   }

//   std::cout << "\n=== STRESS-CONTROLLED SIMULATION COMPLETE ===" <<
//   std::endl;
// }

AcousticElasticResults computeAcousticElasticConstants(const AcousticStudyConfig& config) {
    AcousticElasticResults res;
    
    if (config.mode == 3) {
        // Zanzotto Potential Mode 3
        energy_functions::set_parameters(config.beta, config.K);
        auto ec = energy_functions::calculate_elastic_constants(config.beta, config.K, config.lattice_type);
        res.C11 = ec.C11;
        res.C22 = ec.C22;
        res.C12 = ec.C12;
        res.C66 = ec.C66;
        res.Zener_A = ec.Zener_A;
    } else {
        // Mode 1 (LJ v2) or Mode 2 (LJ v3)
        double scale = 1.0;
        double r_cutoff = 2.5;
        std::function<double(double)> dpot;
        std::function<double(double)> d2pot;
        
        if (config.mode == 1) {
            scale = 0.6872044091828517;
            r_cutoff = 2.5;
            dpot = lennard_jones_energy_der_v2;
            d2pot = lennard_jones_energy_sder_v2;
        } else {
            scale = 0.996407146941421;
            r_cutoff = 1.86602540378444;
            dpot = lennard_jones_energy_der_v3;
            d2pot = lennard_jones_energy_sder_v3;
        }
        
        Eigen::Matrix2d C_ref = Eigen::Matrix2d::Identity();
        if (config.lattice_type == "triangular") {
            TriangularLatticeCalculator calc(scale, r_cutoff);
            double norm = calc.getUnitCellArea();
            HessianComponents hess = calc.calculate_dseconderivative_components(C_ref, dpot, d2pot);
            res.C11 = 4.0 * hess.c11_c11 / norm;
            res.C22 = 4.0 * hess.c22_c22 / norm;
            res.C12 = 4.0 * hess.c11_c22 / norm;
            res.C66 = 1.0 * hess.c12_c12 / norm;
        } else {
            SquareLatticeCalculator calc(scale, r_cutoff);
            double norm = calc.getUnitCellArea();
            HessianComponents hess = calc.calculate_dseconderivative_components(C_ref, dpot, d2pot);
            res.C11 = 4.0 * hess.c11_c11 / norm;
            res.C22 = 4.0 * hess.c22_c22 / norm;
            res.C12 = 4.0 * hess.c11_c22 / norm;
            res.C66 = 1.0 * hess.c12_c12 / norm;
        }
        
        double diff = res.C11 - res.C12;
        res.Zener_A = (std::abs(diff) > 1e-12) ? (2.0 * res.C66 / diff) : 1.0;
    }
    
    // 2D Bulk modulus K_2D = (C11 + C22 + 2*C12) / 4
    res.bulk_K = 0.25 * (res.C11 + res.C22 + 2.0 * res.C12);
    res.shear_mu = res.C66;
    
    // Born stability check in 2D: C11 > 0, C22 > 0, C11*C22 - C12^2 > 0, C66 > 0
    res.born_stable = (res.C11 > 0.0) && (res.C22 > 0.0) && 
                      (res.C11 * res.C22 - res.C12 * res.C12 > 0.0) && 
                      (res.C66 > 0.0);
                      
    res.is_isotropic = (std::abs(res.Zener_A - 1.0) < 0.01) && (std::abs(res.C11 - res.C22) < 1e-4);
    
    // Build informative commentary
    std::ostringstream comment;
    if (res.is_isotropic) {
        comment << "Elastically isotropic in the linear regime (Zener ratio A = " 
                << std::fixed << std::setprecision(3) << res.Zener_A << " ≈ 1.000). "
                << "Plane wave speeds and shear stiffness are isotropic in the plane. ";
        if (config.lattice_type == "triangular") {
            comment << "Triangular 2D lattices with central pair forces inherently satisfy hexagonal isotropy (C11 - C12 = 2*C66) and the Cauchy relation (C12 = C66).";
        }
    } else if (res.Zener_A > 1.0) {
        comment << "Elastically anisotropic with Zener ratio A = " 
                << std::fixed << std::setprecision(3) << res.Zener_A << " > 1. "
                << "Shear resistance is softer along principal <10> axes and stiffer along diagonal <11> directions (C66 > (C11-C12)/2).";
    } else {
        comment << "Elastically anisotropic with Zener ratio A = " 
                << std::fixed << std::setprecision(3) << res.Zener_A << " < 1. "
                << "Shear resistance is softer along diagonal <11> directions than principal axes, typical for square lattices susceptible to diagonal shear instability.";
    }
    
    if (res.born_stable) {
        comment << " Born stability criteria satisfied: locally stable ground state.";
    } else {
        comment << " WARNING: Born stability criteria violated: reference lattice is unstable to affine perturbations.";
    }
    
    res.anisotropy_comment = comment.str();
    return res;
}

void printAcousticElasticConstantsJSON(const AcousticStudyConfig& config) {
    AcousticElasticResults ec = computeAcousticElasticConstants(config);
    std::cout << "{\n";
    std::cout << "  \"mode\": " << config.mode << ",\n";
    std::cout << "  \"lattice_type\": \"" << config.lattice_type << "\",\n";
    std::cout << "  \"beta\": " << config.beta << ",\n";
    std::cout << "  \"K\": " << config.K << ",\n";
    std::cout << "  \"C11\": " << ec.C11 << ",\n";
    std::cout << "  \"C22\": " << ec.C22 << ",\n";
    std::cout << "  \"C12\": " << ec.C12 << ",\n";
    std::cout << "  \"C66\": " << ec.C66 << ",\n";
    std::cout << "  \"Zener_A\": " << ec.Zener_A << ",\n";
    std::cout << "  \"bulk_modulus\": " << ec.bulk_K << ",\n";
    std::cout << "  \"shear_modulus\": " << ec.shear_mu << ",\n";
    std::cout << "  \"born_stable\": " << (ec.born_stable ? "true" : "false") << ",\n";
    std::cout << "  \"is_isotropic\": " << (ec.is_isotropic ? "true" : "false") << ",\n";
    std::cout << "  \"anisotropy_comment\": \"" << ec.anisotropy_comment << "\"\n";
    std::cout << "}\n";
}

ShearPathResult computeSimpleShearPath(const AcousticStudyConfig& config, double alpha_min, double alpha_max, double dalpha) {
    ShearPathResult res;
    std::unique_ptr<BaseLatticeCalculator> calculator_ptr;
    std::function<double(double)> potential_func;
    std::function<double(double)> potential_func_der;
    std::function<double(double)> potential_func_sder;
    double normalisation = 1.0;

    if (config.mode == 3) {
        energy_functions::set_parameters(config.beta, config.K);
        auto calc = std::make_unique<Strain_Energy_LatticeCalculator>(1.0);
        normalisation = calc->getUnitCellArea();
        calculator_ptr = std::move(calc);
        potential_func = [](double r) { return 1.0; };
        potential_func_der = [](double r) { return 1.0; };
        potential_func_sder = [](double r) { return 1.0; };
    } else {
        double scale = (config.mode == 1) ? 0.6872044091828517 : 0.996407146941421;
        double r_cutoff = (config.mode == 1) ? 2.5 : 1.86602540378444;
        potential_func = (config.mode == 1) ? lennard_jones_energy_v2 : lennard_jones_energy_v3;
        potential_func_der = (config.mode == 1) ? lennard_jones_energy_der_v2 : lennard_jones_energy_der_v3;
        potential_func_sder = (config.mode == 1) ? lennard_jones_energy_sder_v2 : lennard_jones_energy_sder_v3;

        auto calc = std::make_unique<SquareLatticeCalculator>(scale, r_cutoff);
        normalisation = calc->getUnitCellArea();
        if (config.lattice_type == "triangular") {
            normalisation *= std::sqrt(3.0) / 2.0;
        }
        calculator_ptr = std::move(calc);
    }

    Eigen::Matrix2d C_ref_energy = Eigen::Matrix2d::Identity();
    if (config.lattice_type == "triangular" && config.mode != 3) {
        double a2 = std::sqrt(4.0 / 3.0);
        C_ref_energy << a2, 0.5 * a2, 0.5 * a2, a2;
    }
    double energy_ref = calculator_ptr->calculate_energy(C_ref_energy, potential_func, 0.0) / normalisation;

    // Starting deformation gradient F0
    Eigen::Matrix2d F0 = Eigen::Matrix2d::Identity();
    if (config.lattice_type == "triangular") {
        double a = std::pow(4.0 / 3.0, 0.25);
        F0 << a, a * 0.5,
              0.0, a * std::sqrt(3.0) / 2.0;
    }

    double prev_det = 0.0;
    double prev_alpha = 0.0;
    ShearPathPoint prev_pt;
    bool has_prev = false;

    if (dalpha <= 0.0) dalpha = 0.005;

    for (double alpha = alpha_min; alpha <= alpha_max + 1e-9; alpha += dalpha) {
        Eigen::Matrix2d F_shear;
        F_shear << 1.0, alpha,
                   0.0, 1.0;

        Eigen::Matrix2d F_tot = F_shear * F0;
        Eigen::Matrix2d C_orig = F_tot.transpose() * F_tot;

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver(C_orig);
        if (eigensolver.eigenvalues().minCoeff() <= 1e-10) continue;

        lagrange::Result reduction = lagrange::reduce(C_orig);
        Eigen::Matrix2d C_red = reduction.C_reduced;
        Eigen::Matrix2d Z = reduction.m_matrix;

        Eigen::Matrix2d S_sqrt = eigensolver.eigenvalues().cwiseSqrt().asDiagonal();
        Eigen::Matrix2d F_eval = eigensolver.eigenvectors() * S_sqrt * eigensolver.eigenvectors().transpose();

        AcousticTensor acoustic_tensor(F_eval, C_red, Z);
        acoustic_tensor.computeEnergyDerivatives(*calculator_ptr, potential_func_der, potential_func_sder, normalisation);
        AcousticAnalysis res_ac = acoustic_tensor.analyzeAcousticTensor(false);

        double t_proj = 2.0 / (2.0 + C_orig(0, 0) + C_orig(1, 1));
        double x = t_proj * (C_orig(0, 0) - C_orig(1, 1)) * 0.5;
        double y = t_proj * C_orig(0, 1);

        double t_proj_red = 2.0 / (2.0 + C_red(0, 0) + C_red(1, 1));
        double x_red = t_proj_red * (C_red(0, 0) - C_red(1, 1)) * 0.5;
        double y_red = t_proj_red * C_red(0, 1);

        double raw_e = calculator_ptr->calculate_energy(C_red, potential_func, 0.0) / normalisation;
        double strain_e = raw_e - energy_ref;

        ShearPathPoint pt;
        pt.alpha = alpha;
        pt.c11 = C_orig(0, 0); pt.c22 = C_orig(1, 1); pt.c12 = C_orig(0, 1);
        pt.c11_red = C_red(0, 0); pt.c22_red = C_red(1, 1); pt.c12_red = C_red(0, 1);
        pt.x = x; pt.y = y;
        pt.x_red = x_red; pt.y_red = y_red;
        pt.min_detAc = res_ac.detAc;
        pt.energy = strain_e;
        pt.is_stable = (res_ac.detAc > 0.0);

        if (has_prev) {
            if ((prev_det > 0.0 && res_ac.detAc <= 0.0) || (prev_det < 0.0 && res_ac.detAc >= 0.0)) {
                double frac = (0.0 - prev_det) / (res_ac.detAc - prev_det);
                double alpha_c = prev_alpha + frac * (alpha - prev_alpha);
                res.critical_alphas.push_back(alpha_c);

                ShearPathPoint cp;
                cp.alpha = alpha_c;
                cp.x = prev_pt.x + frac * (pt.x - prev_pt.x);
                cp.y = prev_pt.y + frac * (pt.y - prev_pt.y);
                cp.x_red = prev_pt.x_red + frac * (pt.x_red - prev_pt.x_red);
                cp.y_red = prev_pt.y_red + frac * (pt.y_red - prev_pt.y_red);
                cp.c11 = prev_pt.c11 + frac * (pt.c11 - prev_pt.c11);
                cp.c22 = prev_pt.c22 + frac * (pt.c22 - prev_pt.c22);
                cp.c12 = prev_pt.c12 + frac * (pt.c12 - prev_pt.c12);
                cp.min_detAc = 0.0;
                cp.energy = prev_pt.energy + frac * (pt.energy - prev_pt.energy);
                cp.is_stable = false;
                res.critical_points.push_back(cp);
            }
        }

        prev_det = res_ac.detAc;
        prev_alpha = alpha;
        prev_pt = pt;
        has_prev = true;

        res.path.push_back(pt);
    }

    return res;
}

void runParametricAcousticStudy(const AcousticStudyConfig& config) {
    std::cout << "=========================================================\n";
    std::cout << "Starting Parametric Acoustic Study\n";
    std::cout << "  Mode: " << config.mode;
    if (config.mode == 3) std::cout << " (Zanzotto Continuum, beta=" << config.beta << ", K=" << config.K << ")\n";
    else if (config.mode == 1) std::cout << " (Lennard-Jones v2 V(r))\n";
    else if (config.mode == 2) std::cout << " (Lennard-Jones v3 V(r))\n";
    std::cout << "  Lattice Symmetry: " << config.lattice_type << "\n";
    std::cout << "  Radius t_max: " << config.t_max << ", dt: " << config.dt << "\n";
    std::cout << "  Angle dp: " << config.dp << " rad\n";
    std::cout << "  Output files: " << config.dat_output << " | " << config.json_output << "\n";
    std::cout << "=========================================================\n";

    AcousticElasticResults ec = computeAcousticElasticConstants(config);
    std::cout << "\n--- Ground State Elastic Constants (C = I) ---\n";
    std::cout << "  C11: " << ec.C11 << " | C22: " << ec.C22 << "\n";
    std::cout << "  C12: " << ec.C12 << " | C66: " << ec.C66 << "\n";
    std::cout << "  Zener Anisotropy Ratio A = 2*C66/(C11-C12): " << ec.Zener_A << "\n";
    std::cout << "  Bulk Modulus K_2D: " << ec.bulk_K << " | Shear Modulus: " << ec.shear_mu << "\n";
    std::cout << "  Born Stability: " << (ec.born_stable ? "STABLE" : "UNSTABLE") << "\n";
    std::cout << "  Analysis: " << ec.anisotropy_comment << "\n\n";

    if (config.only_elastic_constants) {
        printAcousticElasticConstantsJSON(config);
        return;
    }

    std::unique_ptr<BaseLatticeCalculator> calculator_ptr;
    std::function<double(double)> potential_func;
    std::function<double(double)> potential_func_der;
    std::function<double(double)> potential_func_sder;
    double normalisation = 1.0;

    if (config.mode == 3) {
        energy_functions::set_parameters(config.beta, config.K);
        auto calc = std::make_unique<Strain_Energy_LatticeCalculator>(1.0);
        normalisation = calc->getUnitCellArea();
        calculator_ptr = std::move(calc);
        potential_func = [](double r) -> double { return 1.0; };
        potential_func_der = [](double r) -> double { return 1.0; };
        potential_func_sder = [](double r) -> double { return 1.0; };
    } else {
        double scale = (config.mode == 1) ? 0.6872044091828517 : 0.996407146941421;
        double r_cutoff = (config.mode == 1) ? 2.5 : 1.86602540378444;
        potential_func = (config.mode == 1) ? lennard_jones_energy_v2 : lennard_jones_energy_v3;
        potential_func_der = (config.mode == 1) ? lennard_jones_energy_der_v2 : lennard_jones_energy_der_v3;
        potential_func_sder = (config.mode == 1) ? lennard_jones_energy_sder_v2 : lennard_jones_energy_sder_v3;

        auto calc = std::make_unique<SquareLatticeCalculator>(scale, r_cutoff);
        normalisation = calc->getUnitCellArea();
        if (config.lattice_type == "triangular") {
            normalisation *= std::sqrt(3.0) / 2.0;
        }
        calculator_ptr = std::move(calc);
    }

    Eigen::Matrix2d C_ref_energy = Eigen::Matrix2d::Identity();
    if (config.lattice_type == "triangular" && config.mode != 3) {
        double a2 = std::sqrt(4.0 / 3.0);
        C_ref_energy << a2, 0.5 * a2, 0.5 * a2, a2;
    }
    double energy_ref = calculator_ptr->calculate_energy(C_ref_energy, potential_func, 0.0) / normalisation;

    // Open .dat file
    std::ofstream dat_file(config.dat_output);
    dat_file << std::scientific << std::setprecision(8);
    dat_file << "# t p c11 c22 c12 c11_red c22_red c12_red min_detAc angle_deg third_condition x y x_red y_red energy\n";

    // Estimate total points
    int t_steps = static_cast<int>(std::ceil(config.t_max / config.dt)) + 1;
    int p_steps = static_cast<int>(std::ceil((2.0 * M_PI) / config.dp)) + 1;
    int total_estimated = t_steps * p_steps;
    std::cout << "Scanning grid: ~" << total_estimated << " points...\n";

    struct OutputPoint {
        double t, p;
        double c11, c22, c12;
        double c11_red, c22_red, c12_red;
        double min_detAc, angle_deg;
        bool third_condition;
        double x, y;
        double x_red, y_red;
        bool is_stable;
        double energy;
    };

    std::vector<OutputPoint> points;
    points.reserve(total_estimated);

    int count = 0;
    int stable_count = 0;
    int unstable_count = 0;
    double min_det_min = 1e20;
    double min_det_max = -1e20;
    double energy_min = 1e20;
    double energy_max = -1e20;

    for (double t = 0.0; t <= config.t_max + 1e-9; t += config.dt) {
        for (double p = -M_PI; p <= M_PI + 1e-9; p += config.dp) {
            try {
                double c11 = std::cosh(t) + std::sinh(t) * std::sin(p);
                double c22 = std::cosh(t) - std::sinh(t) * std::sin(p);
                double c12 = std::sinh(t) * std::cos(p);

                Eigen::Matrix2d C_original;
                C_original << c11, c12, c12, c22;

                Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver(C_original);
                if (eigensolver.eigenvalues().minCoeff() <= 1e-10) {
                    continue;
                }

                lagrange::Result reduction_result = lagrange::reduce(C_original);
                Eigen::Matrix2d C = reduction_result.C_reduced;
                Eigen::Matrix2d Z = reduction_result.m_matrix;
                bool third_condition = reduction_result.third_condition_satisfied;

                Eigen::Matrix2d S_sqrt = eigensolver.eigenvalues().cwiseSqrt().asDiagonal();
                Eigen::Matrix2d F = eigensolver.eigenvectors() * S_sqrt * eigensolver.eigenvectors().transpose();

                AcousticTensor acoustic_tensor(F, C, Z);
                acoustic_tensor.computeEnergyDerivatives(*calculator_ptr, potential_func_der, potential_func_sder, normalisation);

                AcousticAnalysis result = acoustic_tensor.analyzeAcousticTensor(false);

                // User's Poincaré stereographic projection formula:
                // t_proj = 2. / (2. + c11 + c22)
                // x = t_proj * (c11 - c22) / 2.
                // y = t_proj * c12
                double t_proj = 2.0 / (2.0 + c11 + c22);
                double x = t_proj * (c11 - c22) * 0.5;
                double y = t_proj * c12;

                double t_proj_red = 2.0 / (2.0 + C(0, 0) + C(1, 1));
                double x_red = t_proj_red * (C(0, 0) - C(1, 1)) * 0.5;
                double y_red = t_proj_red * C(0, 1);

                double raw_energy = calculator_ptr->calculate_energy(C, potential_func, 0.0) / normalisation;
                double strain_energy = raw_energy - energy_ref;

                bool is_stable = (result.detAc > 0.0);
                if (is_stable) stable_count++;
                else unstable_count++;

                min_det_min = std::min(min_det_min, result.detAc);
                min_det_max = std::max(min_det_max, result.detAc);
                energy_min = std::min(energy_min, strain_energy);
                energy_max = std::max(energy_max, strain_energy);

                // Write to dat file
                dat_file << t << " " << p << " "
                         << c11 << " " << c22 << " " << c12 << " "
                         << C(0, 0) << " " << C(1, 1) << " " << C(0, 1) << " "
                         << result.detAc << " " << result.xsi << " " << (third_condition ? 1 : 0) << " "
                         << x << " " << y << " " << x_red << " " << y_red << " "
                         << strain_energy << "\n";

                OutputPoint pt;
                pt.t = t; pt.p = p;
                pt.c11 = c11; pt.c22 = c22; pt.c12 = c12;
                pt.c11_red = C(0, 0); pt.c22_red = C(1, 1); pt.c12_red = C(0, 1);
                pt.min_detAc = result.detAc;
                pt.angle_deg = result.xsi;
                pt.third_condition = third_condition;
                pt.x = x; pt.y = y;
                pt.x_red = x_red; pt.y_red = y_red;
                pt.is_stable = is_stable;
                pt.energy = strain_energy;
                points.push_back(pt);

                count++;
                if (count % 250 == 0) {
                    std::cout << "Computed " << count << " points (t=" 
                              << std::fixed << std::setprecision(2) << t 
                              << ", p=" << std::setprecision(2) << p << ")..." << std::endl;
                }
            } catch (const std::exception& e) {
                continue;
            }
        }
    }

    dat_file.close();

    // Write JSON file
    std::ofstream json_file(config.json_output);
    json_file << "{\n";
    json_file << "  \"mode\": " << config.mode << ",\n";
    json_file << "  \"mode_name\": \"" << (config.mode == 3 ? "Zanzotto Continuum" : (config.mode == 1 ? "Lennard-Jones v2 V(r)" : "Lennard-Jones v3 V(r)")) << "\",\n";
    json_file << "  \"lattice_type\": \"" << config.lattice_type << "\",\n";
    json_file << "  \"beta\": " << config.beta << ",\n";
    json_file << "  \"K\": " << config.K << ",\n";
    json_file << "  \"t_max\": " << config.t_max << ",\n";
    json_file << "  \"dt\": " << config.dt << ",\n";
    json_file << "  \"dp\": " << config.dp << ",\n";
    json_file << "  \"elastic_constants\": {\n";
    json_file << "    \"C11\": " << ec.C11 << ",\n";
    json_file << "    \"C22\": " << ec.C22 << ",\n";
    json_file << "    \"C12\": " << ec.C12 << ",\n";
    json_file << "    \"C66\": " << ec.C66 << ",\n";
    json_file << "    \"Zener_A\": " << ec.Zener_A << ",\n";
    json_file << "    \"bulk_modulus\": " << ec.bulk_K << ",\n";
    json_file << "    \"shear_modulus\": " << ec.shear_mu << ",\n";
    json_file << "    \"born_stable\": " << (ec.born_stable ? "true" : "false") << ",\n";
    json_file << "    \"is_isotropic\": " << (ec.is_isotropic ? "true" : "false") << ",\n";
    json_file << "    \"anisotropy_comment\": \"" << ec.anisotropy_comment << "\"\n";
    json_file << "  },\n";
    json_file << "  \"summary\": {\n";
    json_file << "    \"total_points\": " << count << ",\n";
    json_file << "    \"stable_points\": " << stable_count << ",\n";
    json_file << "    \"unstable_points\": " << unstable_count << ",\n";
    json_file << "    \"min_detAc_min\": " << (count > 0 ? min_det_min : 0.0) << ",\n";
    json_file << "    \"min_detAc_max\": " << (count > 0 ? min_det_max : 0.0) << ",\n";
    json_file << "    \"energy_min\": " << (count > 0 ? energy_min : 0.0) << ",\n";
    json_file << "    \"energy_max\": " << (count > 0 ? energy_max : 0.0) << "\n";
    json_file << "  },\n";
    json_file << "  \"points\": [\n";

    for (size_t i = 0; i < points.size(); ++i) {
        const auto& pt = points[i];
        json_file << "    {"
                  << "\"t\":" << pt.t << ","
                  << "\"p\":" << pt.p << ","
                  << "\"c11\":" << pt.c11 << ","
                  << "\"c22\":" << pt.c22 << ","
                  << "\"c12\":" << pt.c12 << ","
                  << "\"c11_red\":" << pt.c11_red << ","
                  << "\"c22_red\":" << pt.c22_red << ","
                  << "\"c12_red\":" << pt.c12_red << ","
                  << "\"min_detAc\":" << pt.min_detAc << ","
                  << "\"energy\":" << pt.energy << ","
                  << "\"angle_deg\":" << pt.angle_deg << ","
                  << "\"third_condition\":" << (pt.third_condition ? "true" : "false") << ","
                  << "\"x\":" << pt.x << ","
                  << "\"y\":" << pt.y << ","
                  << "\"x_red\":" << pt.x_red << ","
                  << "\"y_red\":" << pt.y_red << ","
                  << "\"is_stable\":" << (pt.is_stable ? "true" : "false")
                  << "}" << (i + 1 < points.size() ? ",\n" : "\n");
    }

    json_file << "  ],\n";

    // Compute simple shear path on Poincaré space and find critical alpha values
    ShearPathResult shear_res = computeSimpleShearPath(config, config.alpha_min, config.alpha_max, config.dalpha);
    std::cout << "\n--- Simple Shear Path Analysis (" << config.lattice_type << " reference) ---\n";
    std::cout << "  Alpha range: [" << config.alpha_min << ", " << config.alpha_max << "] step: " << config.dalpha << "\n";
    std::cout << "  Total trajectory points: " << shear_res.path.size() << "\n";
    if (shear_res.critical_alphas.empty()) {
        std::cout << "  No acoustic tensor zero-crossings detected in range.\n";
    } else {
        std::cout << "  Found " << shear_res.critical_alphas.size() << " critical alpha transition(s) (Loss of Ellipticity):\n";
        for (size_t i = 0; i < shear_res.critical_alphas.size(); ++i) {
            std::cout << "    • α_c[" << (i + 1) << "] = " << std::fixed << std::setprecision(5) << shear_res.critical_alphas[i] << "\n";
        }
    }

    json_file << "  \"shear_path\": {\n";
    json_file << "    \"alpha_min\": " << config.alpha_min << ",\n";
    json_file << "    \"alpha_max\": " << config.alpha_max << ",\n";
    json_file << "    \"dalpha\": " << config.dalpha << ",\n";
    json_file << "    \"critical_alphas\": [";
    for (size_t i = 0; i < shear_res.critical_alphas.size(); ++i) {
        json_file << shear_res.critical_alphas[i] << (i + 1 < shear_res.critical_alphas.size() ? ", " : "");
    }
    json_file << "],\n";
    json_file << "    \"critical_points\": [\n";
    for (size_t i = 0; i < shear_res.critical_points.size(); ++i) {
        const auto& cp = shear_res.critical_points[i];
        json_file << "      {\"alpha\":" << cp.alpha << ",\"x\":" << cp.x << ",\"y\":" << cp.y
                  << ",\"x_red\":" << cp.x_red << ",\"y_red\":" << cp.y_red
                  << ",\"c11\":" << cp.c11 << ",\"c22\":" << cp.c22 << ",\"c12\":" << cp.c12
                  << ",\"energy\":" << cp.energy << "}" << (i + 1 < shear_res.critical_points.size() ? ",\n" : "\n");
    }
    json_file << "    ],\n";
    json_file << "    \"points\": [\n";
    for (size_t i = 0; i < shear_res.path.size(); ++i) {
        const auto& pt = shear_res.path[i];
        json_file << "      {\"alpha\":" << pt.alpha << ",\"x\":" << pt.x << ",\"y\":" << pt.y
                  << ",\"x_red\":" << pt.x_red << ",\"y_red\":" << pt.y_red
                  << ",\"c11\":" << pt.c11 << ",\"c22\":" << pt.c22 << ",\"c12\":" << pt.c12
                  << ",\"min_detAc\":" << pt.min_detAc << ",\"energy\":" << pt.energy
                  << ",\"is_stable\":" << (pt.is_stable ? "true" : "false") << "}"
                  << (i + 1 < shear_res.path.size() ? ",\n" : "\n");
    }
    json_file << "    ]\n";
    json_file << "  }\n";
    json_file << "}\n";
    json_file.close();

    std::cout << "\n✓ Parametric Acoustic Study Complete!\n";
    std::cout << "  Processed: " << count << " points (Stable: " << stable_count 
              << ", Unstable: " << unstable_count << ")\n";
    std::cout << "  Data saved to: " << config.dat_output << "\n";
    std::cout << "  JSON saved to: " << config.json_output << "\n";
}
