#include "../../include/experiments/experiment_includes.h"
#include "../../include/optimization/StiffnessSpectrum.h"

StiffnessEigenSolver parse_stiffness_eigen_solver(const std::string &name) {
  if (name == "fast" || name == "cholesky") return StiffnessEigenSolver::Fast;
  if (name == "legacy" || name == "old" || name == "itensor")
    return StiffnessEigenSolver::Legacy;
  throw std::invalid_argument("Unknown eigen solver: " + name + " (fast | legacy)");
}

std::string to_string(StiffnessEigenSolver solver) {
  return solver == StiffnessEigenSolver::Fast ? "fast" : "legacy";
}

namespace {

// n_eig lowest eigenpairs of the stiffness at x with lowest_stiffness_modes, returned in
// the layout of the legacy path: the two uniform translations, which the fast solver
// projects out, are added back explicitly (eigenvalue = Rayleigh quotient, ~1e-15), so
// the rigid-mode detection and the output columns below stay the same.
EigenResults fast_lowest_eigenpairs(const alglib::real_1d_array &x, UserData *userData,
                                    int n_eig) {
  FastStiffnessAssembler fast_assembler;
  const Eigen::SparseMatrix<double> &K =
      fast_assembler.assembleStiffness(x, userData, /*project_psd=*/false);
  const int n = static_cast<int>(K.rows());
  const int n_translations = stiffness_has_translation_modes(K) ? 2 : 0;

  StiffnessSpectrumOptions opt;
  opt.n_modes = n_eig - n_translations;
  opt.deflate_translations = n_translations ? 1 : 0;
  opt.verbose = true;
  const StiffnessSpectrum spectrum = lowest_stiffness_modes(K, opt);

  EigenResults results;
  results.num_computed = n_translations + spectrum.num_computed;
  results.eigenvalues.resize(results.num_computed);
  results.eigenvectors.resize(n, results.num_computed);
  const int m = n / 2;
  for (int block = 0; block < n_translations; block++) {
    Eigen::VectorXd t = Eigen::VectorXd::Zero(n);
    t.segment(block * m, m).setConstant(1.0 / std::sqrt(static_cast<double>(m)));
    results.eigenvalues(block) = t.dot(K * t);
    results.eigenvectors.col(block) = t;
  }
  for (int k = 0; k < spectrum.num_computed; k++) {
    results.eigenvalues(n_translations + k) = spectrum.eigenvalues(k);
    results.eigenvectors.col(n_translations + k) = spectrum.eigenvectors.col(k);
  }
  return results;
}

} // namespace

void analyze_data_from_folder(int caller_id, int nx, int ny, int iter_start,
                              int iter_end, int n_eig,
                              StiffnessEigenSolver eig_solver) {
  std::cout << "Stiffness eigen solver: " << to_string(eig_solver) << std::endl;

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

  std::cout << "interior_mapping.size(): " << interior_mapping.size()
            << std::endl;
  std::cout << "full_mapping.size(): " << full_mapping.size() << std::endl;

  Point2D domain_dims_point(domain_dims.size_x, domain_dims.size_y);

  // ==================== SPECIFY ITERATION RANGE ====================
  int iteration_start = iter_start;
  int iteration_end = iter_end; // Change these values as needed
  int iteration_step =
      2; // Step size (e.g., 1 for every iteration, 2 for every other)

  std::cout << "Processing iterations from " << iteration_start << " to "
            << iteration_end << " (step: " << iteration_step << ")"
            << std::endl;

  // ==================== OPEN OUTPUT FILE ====================
  std::string output_dir = "./second_order_analysis";
  std::filesystem::create_directories(output_dir);
  output_dir = "./second_order_analysis";
  std::filesystem::create_directories(output_dir);

  std::string eigenvalue_file = output_dir + "/eigenvalues_vs_shear.dat";
  std::ofstream eigen_out(eigenvalue_file);
  if (!eigen_out.is_open()) {
    std::cerr << "Error: Could not open file " << eigenvalue_file << std::endl;
    exit(EXIT_FAILURE);
  }

  // Write header
  eigen_out << "# F_ext[0,1] (shear), min_eigenvalue, "
               "first_non_rigid_eigenvalue, num_rigid_modes\n";
  eigen_out << std::scientific << std::setprecision(12);

  std::cout << "Will save eigenvalue data to: " << eigenvalue_file << std::endl;

  // ==================== LOOP OVER ITERATIONS ====================
  for (int iteration = iteration_start; iteration <= iteration_end;
       iteration += iteration_step) {

    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "PROCESSING ITERATION " << iteration << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    ConfigurationSaver::IterationData data =
        ConfigurationSaver::readIterationData(iteration);

    // Convert Eigen::Vector2d positions to Point2D objects
    square_points.resize(0);
    square_points.reserve(data.positions.size());

    square_points_ref.resize(0);
    square_points_ref.reserve(data.positions.size());

    for (const auto &pos : data.positions) {
      square_points.emplace_back(pos.x(), pos.y());
      square_points_ref.emplace_back(pos.x(), pos.y());
    }

    // F_ext is already Eigen::Matrix2d
    Eigen::Matrix2d F_ext = data.F_ext;
    double shear_component = F_ext(0, 1); // Extract shear component

    std::cout << "Loaded " << square_points.size() << " points" << std::endl;
    std::cout << "F_ext matrix:\n" << F_ext << std::endl;
    std::cout << "Shear component F_ext[0,1] = " << shear_component
              << std::endl;

    // Example: Access first point
    std::cout << "First point: (" << square_points[0].coord.x() << ", "
              << square_points[0].coord.y() << ")" << std::endl;

    alglib::real_1d_array free_dofs;
    int n_free_nodes = interior_mapping.size();
    free_dofs.setlength(2 *
                        interior_mapping.size()); // [u0, u1, ..., v0, v1, ...]
    map_points_to_solver_array(free_dofs, square_points, interior_mapping,
                               n_free_nodes);

    // Load saved elements instead of remeshing
    auto [elements, active_elements] = ConfigurationSaver::loadElements(
        iteration, dndx, full_mapping, square_points_ref);

    double element_area, ref_element_area;
    for (auto &element : elements) {
      element.set_reference_mesh(square_points);
      element.set_dof_mapping(
          full_mapping); // or interior_mapping depending on needs
      // double jac = element.calculate_shape_derivatives(x);  // current
      // positions
      const Eigen::Matrix<double, 3, 2> &dndx = element.getDNdX();
      element.setExternalDeformation(F_ext);
      element.calculate_deformation_gradient(square_points);

      element.calculateCurrentArea(square_points);
      element_area = element.getArea();
      ref_element_area = element.getReferenceArea();

      // std::cout<<"ref_element_area: "<<ref_element_area<<std::endl;
    }

    std::cout << "Loaded " << elements.size() << " elements, "
              << active_elements.size() << " active" << std::endl;

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

    // ==================== CREATE USER DATA ====================
    bool plasticity = false;
    UserData userData(square_points, elements, calculator, potential_func,
                      potential_func_der, zero, optimal_lattice_parameter,
                      F_ext, interior_mapping, full_mapping, active_elements,
                      plasticity);

    // ConfigurationSaver::saveTriangleData(&userData, 100 + iteration,
    //                                      domain_dims, offsets, full_mapping);

    auto [num_dislocations_post, coordination_post] =
        DefectAnalysis::analyzeDefectsInReferenceConfig(
            &userData, 100 + iteration, dndx, offsets, original_domain_map,
            translation_map, domain_dims_point, element_area, pbc, true);

    // ConfigurationSaver::writeToVTK(userData.points, userData.elements,
    //                                &userData, 100 + iteration, true,
    //                                coordination_post, shear_component);

    // ==================== PREPARE OPTIMIZATION ====================
    alglib::real_1d_array x;
    int n_vars = interior_mapping.size();
    x.setlength(2 * n_vars);
    map_points_to_solver_array(x, square_points, interior_mapping, n_vars);

    // 1. Create FEM Hessian assembler
    FEMHessianAssembler assembler;

    // 2. Set energy parameters (REQUIRED before calling
    // assembleGlobalStiffness)
    double normalisation = calculator.getUnitCellArea();
    assembler.setEnergyParameters(&calculator, potential_func_der,
                                  potential_func_sder,
                                  calculator.getUnitCellArea());

    // Compute eigenvalues (both paths return the translations among the n_eig values)
    const auto eig_start = std::chrono::high_resolution_clock::now();
    EigenResults results;
    if (eig_solver == StiffnessEigenSolver::Fast) {
      results = fast_lowest_eigenpairs(x, &userData, n_eig);
    } else {
      Eigen::SparseMatrix<double> global_stiffness =
          assembler.assembleGlobalStiffness(
              elements,
              square_points, // <- ADDED: Current node positions
              2 * n_free_nodes, full_mapping);
      std::cout << "--- Sparse Matrix Output (Row, Col, Value) ---" << std::endl;

      // Iterate over the sparse matrix efficiently
      // outerSize() is usually the number of columns (for column-major matrices)
      // for (int k = 0; k < global_stiffness.outerSize(); ++k) {
      //     // InnerIterator iterates over non-zero entries of the k-th column
      //     for (Eigen::SparseMatrix<double>::InnerIterator it(global_stiffness,
      //     k); it; ++it) {

      //         // it.row()   = row index
      //         // it.col()   = column index
      //         // it.value() = the value of the element
      //         std::cout << it.row() << " " << it.col() << " " << it.value() <<
      //         "\n";
      //     }
      // }
      std::cout << "----------------------------------------------" << std::endl;

      // EigenResults results =
      //     assembler.computeSmallestEigenvaluesIterative_armadillo(global_stiffness,
      //                                                           100,0.0);

      results =
          assembler.computeSmallestEigenvaluesIterative_spectra(global_stiffness,
                                                                n_eig, -1);
    }
    std::cout << "Eigenvalues (" << to_string(eig_solver) << "): "
              << std::chrono::duration<double>(
                     std::chrono::high_resolution_clock::now() - eig_start)
                     .count()
              << " s" << std::endl;

    // STEP 1: Sort by ABSOLUTE VALUE to identify rigid body modes
    std::vector<std::pair<double, int>> eigen_pairs;
    for (int i = 0; i < results.num_computed; i++) {
      eigen_pairs.push_back({results.eigenvalues(i), i});
    }
    std::sort(eigen_pairs.begin(), eigen_pairs.end(),
              [](const auto &a, const auto &b) {
                return std::abs(a.first) < std::abs(b.first);
              });

    // STEP 2: Build temporarily sorted arrays for rigid body mode detection
    Eigen::VectorXd sorted_eigenvalues(results.num_computed);
    Eigen::MatrixXd sorted_eigenvectors(results.eigenvectors.rows(),
                                        results.num_computed);
    for (int i = 0; i < results.num_computed; i++) {
      sorted_eigenvalues(i) = results.eigenvalues(eigen_pairs[i].second);
      sorted_eigenvectors.col(i) =
          results.eigenvectors.col(eigen_pairs[i].second);
    }
    results.eigenvalues = sorted_eigenvalues;
    results.eigenvectors = sorted_eigenvectors;

    // STEP 3: Detect rigid body modes
    int num_rigid = assembler.detectRigidBodyModes(results);
    std::cout << "Detected " << num_rigid << " rigid body modes" << std::endl;

    // STEP 4: Create NEW pairs for non-rigid modes and sort them
    std::vector<std::pair<double, int>> non_rigid_pairs;
    for (int i = num_rigid; i < results.num_computed; i++) {
      non_rigid_pairs.push_back({results.eigenvalues(i), i});
    }

    // Sort non-rigid modes from negative to positive
    std::sort(non_rigid_pairs.begin(), non_rigid_pairs.end(),
              [](const auto &a, const auto &b) { return a.first < b.first; });

    // STEP 5: Build FINAL sorted arrays
    // Keep rigid modes at the beginning, then sorted non-rigid modes
    for (int i = 0; i < num_rigid; i++) {
      // Rigid modes stay in their positions (already sorted by abs value)
      sorted_eigenvalues(i) = results.eigenvalues(i);
      sorted_eigenvectors.col(i) = results.eigenvectors.col(i);
    }

    for (int i = 0; i < non_rigid_pairs.size(); i++) {
      sorted_eigenvalues(num_rigid + i) =
          results.eigenvalues(non_rigid_pairs[i].second);
      sorted_eigenvectors.col(num_rigid + i) =
          results.eigenvectors.col(non_rigid_pairs[i].second);
    }

    results.eigenvalues = sorted_eigenvalues;
    results.eigenvectors = sorted_eigenvectors;

    // Print eigenvalue ordering for verification
    std::cout << "\nEigenvalue ordering:" << std::endl;
    std::cout << "Rigid body modes:" << std::endl;
    for (int i = 0; i < num_rigid; i++) {
      std::cout << "  Mode " << i << ": λ = " << sorted_eigenvalues(i)
                << std::endl;
    }
    std::cout << "Non-rigid modes (first 15):" << std::endl;
    for (int i = num_rigid;
         i < std::min(num_rigid + 15, (int)results.num_computed); i++) {
      std::cout << "  Mode " << i << ": λ = " << sorted_eigenvalues(i)
                << std::endl;
    }

    // Extract minimum eigenvalues
    double min_eigenvalue = sorted_eigenvalues(0); // Smallest rigid mode
    double first_non_rigid_eigenvalue =
        (num_rigid < results.num_computed)
            ? sorted_eigenvalues(num_rigid)
            : std::numeric_limits<double>::quiet_NaN();

    std::cout << "\nEigenvalue summary:" << std::endl;
    std::cout << "  Num rigid modes: " << num_rigid << std::endl;
    std::cout << "  Min rigid eigenvalue: " << min_eigenvalue << std::endl;
    std::cout << "  First non-rigid eigenvalue: " << first_non_rigid_eigenvalue
              << std::endl;

    // Write to file
    eigen_out << shear_component << " " << min_eigenvalue << " "
              << first_non_rigid_eigenvalue << " " << num_rigid << "\n";
    eigen_out.flush();

    std::cout << "Saved: shear=" << shear_component
              << ", min_eig=" << min_eigenvalue
              << ", first_non_rigid=" << first_non_rigid_eigenvalue
              << std::endl;
    continue;

    // ONE LINE TO DO EVERYTHING!
    assembler.exportCompleteEigenmodeAnalysis(
        output_dir + "/analysis_" +
            std::to_string(10000 + iteration), // Base filename
        results,                               // Sorted eigenvalue results
        square_points,                         // Mesh points
        full_mapping,                          // DOF mapping
        elements,                              // Elements
        num_rigid,                             // Number of rigid body modes
        0.3,                                   // Participation threshold
        1.                                     // VTK scale factor
    );

    // Check stability
    if (first_non_rigid_eigenvalue < 0) {
      std::cout << "\n*** ITERATION " << iteration
                << ": SYSTEM IS UNSTABLE! ***" << std::endl;
    } else {
      std::cout << "\nIteration " << iteration << ": System is stable."
                << std::endl;
    }

    std::vector<int> localized_modes = assembler.getLocalizedModes(
        results, full_mapping, square_points.size(), num_rigid, 0.3);

    std::cout << "Found " << localized_modes.size() << " localized modes"
              << std::endl;

    // 3. Compute energy landscapes for each localized mode
    std::string energy_dir =
        output_dir + "/energy_iter_" + std::to_string(iteration);
    std::filesystem::create_directories(energy_dir);

    for (int i = 0; i < std::min(20, (int)localized_modes.size()); i++) {
      int mode = i + 2; // localized_modes[i];

      // Convert Eigen to alglib
      // alglib::real_1d_array eigenvec;
      // eigenvec.setlength(results.eigenvectors.rows());
      // for (int j = 0; j < results.eigenvectors.rows(); j++) {
      //     eigenvec[j] = results.eigenvectors(j, mode);
      // }

      Eigen::VectorXd eigenvec = results.eigenvectors.col(mode);

      // Compute and save
      assembler.compute_energy_landscape(free_dofs, eigenvec, &userData, mode,
                                         results.eigenvalues[mode], 1, 160,
                                         energy_dir);
    }

  } // End of iteration loop

  // Close output file
  eigen_out.close();
  std::cout << "\nEigenvalue data saved to: " << eigenvalue_file << std::endl;

  std::cout << "\n" << std::string(60, '=') << std::endl;
  std::cout << "COMPLETED ALL ITERATIONS" << std::endl;
  std::cout << std::string(60, '=') << std::endl;
  // exit(0);
}
