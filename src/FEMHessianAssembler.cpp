#include "FEMHessianAssembler.h"
#include "../include/reductions/LagrangeReduction.h"

#include <algorithm>
#include <iostream>
// SPECTRA INCLUDES
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>

#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>





double FEMHessianAssembler::calculate_total_energy(
    const alglib::real_1d_array &x,
    UserData* userData)
{
    std::vector<ElementTriangle2D>& elements = userData->elements;
    const double normalisation = pow(userData->ideal_lattice_parameter, 2.0);

    double total_energy = 0.0;
    
    #pragma omp parallel reduction(+:total_energy)
    {
        #pragma omp for schedule(guided)
        for (size_t idx = 0; idx < userData->active_elements.size(); idx++) {
            ElementTriangle2D& element = elements[userData->active_elements[idx]];
            
            element.calculate_deformation_gradient(x);
            const Eigen::Matrix2d F = element.getDeformationGradient();
            Eigen::Matrix2d C = F.transpose() * F;
            
            const auto result = lagrange::reduce(C);
            
            // ✅ CORRECT: Use energy_function (V(r)), not derivative_function (V'(r))
            const double element_energy = userData->calculator.calculate_energy(
                result.C_reduced, 
                userData->energy_function,    // ← CORRECT!
                userData->zero_energy
            ) / normalisation;
            
            total_energy += element_energy * element.getReferenceArea();
        }
    }
    
    return total_energy;
}
void FEMHessianAssembler::compute_energy_landscape(
    const alglib::real_1d_array &x0,
    const Eigen::VectorXd &eigenvector,
    UserData* userData,
    int mode_number,
    double eigenvalue,
    double alpha_max ,
    int num_points ,
    const std::string& output_folder) // Changed from int itereration
{
    const int n_vars = x0.length();
    
    // Verify size
    if (eigenvector.size() != n_vars) {
        std::cerr << "ERROR: Size mismatch!" << std::endl;
        return;
    }

    // Reference energy at α=0
    double E0 = calculate_total_energy(x0, userData);

    // ---------------------------------------------------------
    // Construct the Full Filename
    // Use the passed output_folder directly.
    // ---------------------------------------------------------
    std::stringstream filename_ss;
    // Assumes output_folder does not have a trailing slash (based on your snippet)
    filename_ss << output_folder << "/energy_mode_" << std::setfill('0') << std::setw(4) << mode_number << ".dat";
    std::string full_filepath = filename_ss.str();

    std::ofstream file(full_filepath);
    
    if (!file.is_open()) {
        std::cerr << "ERROR: Could not open file for writing: " << full_filepath << std::endl;
        // This might happen if the directory wasn't created externally before calling this
        return;
    }
    
    // ---------------------------------------------------------
    // Write Data
    // ---------------------------------------------------------
    file << "# Mode " << mode_number << ", eigenvalue = " << eigenvalue << std::endl;
    file << "# alpha  E(alpha)  deltaE(alpha)" << std::endl;
    file << std::setprecision(12) << std::scientific;

    std::cout << "\nMode " << mode_number << " (λ = " << eigenvalue << ")" << std::endl;
    std::cout << "E₀ = " << E0 << std::endl;

    // Perturbed configuration array
    alglib::real_1d_array x_perturbed;
    x_perturbed.setlength(n_vars);

    // Loop over α
    for (int i = 0; i < num_points; i++) {
        double alpha = -alpha_max + (2.0 * alpha_max) * i / (num_points - 1);

        // Perturb: x = x₀ + α·u
        for (int j = 0; j < n_vars; j++) {
            x_perturbed[j] = x0[j] + alpha * eigenvector(j);
        }

        // Calculate energy
        double E_alpha = calculate_total_energy(x_perturbed, userData);
        double deltaE = E_alpha - E0;

        // Write to file
        file << alpha << "  " << E_alpha << "  " << deltaE << std::endl;

        // Progress
        if (i % 10 == 0) {
            std::cout << "  α = " << std::setw(7) << std::fixed << std::setprecision(3) 
                      << alpha << ",  ΔE = " << std::scientific << deltaE << std::endl;
        }
    }

    file.close();
    std::cout << "Saved: " << full_filepath << std::endl;
}

void FEMHessianAssembler::extractAcousticTensor(
    const itensor::ITensor &A_tensor, double A[2][2][2][2]) {
  // Extract A_{iKjL} from ITensor with indices (r_, i_, s_, j_)
  // where r_=K (material coord 1), s_=L (material coord 2)
  for (int i = 1; i <= 2; i++) {
    for (int K = 1; K <= 2; K++) {
      for (int j = 1; j <= 2; j++) {
        for (int L = 1; L <= 2; L++) {
          // Access pattern: (r_=K, i_=i, s_=L, j_=j)
          // ITensor uses 1-based indexing, array uses 0-based
          A[i - 1][K - 1][j - 1][L - 1] = A_tensor.elt(K, i, L, j);
        }
      }
    }
  }
}

Eigen::MatrixXd
FEMHessianAssembler::computeElementStiffness(AcousticTensor &acoustic_tensor,
                                             const ElementTriangle2D &element,
                                             double area_weight) {
  // Safety check: ensure energy parameters were set
  if (strain_calculator == nullptr) {
    std::cerr
        << "Error: Energy parameters not set! Call setEnergyParameters() first."
        << std::endl;
    throw std::runtime_error(
        "FEMHessianAssembler: energy parameters not initialized");
  }

  const int num_nodes = 3;                        // Triangular element
  const int spatial_dim = 2;                      // 2D
  const int total_dofs = num_nodes * spatial_dim; // 6 DOFs

  Eigen::MatrixXd K_elem = Eigen::MatrixXd::Zero(total_dofs, total_dofs);

  // Compute energy derivatives using member variables
  // CRITICAL: Dereference pointer with * to pass as reference
  acoustic_tensor.computeEnergyDerivatives(*strain_calculator,
                                           potential_func_der,
                                           potential_func_sder, normalisation);

  // Get shape function derivatives: dN^a/dX_K (3 nodes × 2 coords)
  const Eigen::Matrix<double, 3, 2> &dN_dX = element.getDNdX();

  // Get acoustic tensor A_{iKjL} (Lagrangian frame)
  // Returns ITensor with indices (r_, i_, s_, j_) where r_=K, s_=L
  itensor::ITensor A_tensor = acoustic_tensor.getAcousticTensor(true);

  // Extract to 4D array for easier manipulation
  double A[2][2][2][2];
  extractAcousticTensor(A_tensor, A);

  // Assemble element stiffness: K^{ab}_{ij} = A_{iKjL} * dN^a_K * dN^b_L
  for (int a = 0; a < num_nodes; a++) {         // Node a
    for (int b = 0; b < num_nodes; b++) {       // Node b
      for (int i = 0; i < spatial_dim; i++) {   // DOF direction i
        for (int j = 0; j < spatial_dim; j++) { // DOF direction j

          double K_ab_ij = 0.0;

          // Contract over material coordinates K, L
          for (int K = 0; K < spatial_dim; K++) {
            for (int L = 0; L < spatial_dim; L++) {
              K_ab_ij += A[i][K][j][L] * dN_dX(a, K) * dN_dX(b, L);
            }
          }

          // Multiply by integration weight (includes det(J) for reference
          // element)
          K_ab_ij *= area_weight;

          // Place in element matrix
          // Row: node a, direction i → global index: a*spatial_dim + i
          // Col: node b, direction j → global index: b*spatial_dim + j
          int row = a * spatial_dim + i;
          int col = b * spatial_dim + j;
          K_elem(row, col) = K_ab_ij;
        }
      }
    }
  }

  return K_elem;
}

Eigen::SparseMatrix<double> FEMHessianAssembler::assembleGlobalStiffness(
    std::vector<ElementTriangle2D> &elements,
    const std::vector<Point2D> &current_points, int num_total_dofs,
    const std::vector<std::pair<int, int>> &dof_mapping) {
  // Safety check: ensure energy parameters were set
  if (strain_calculator == nullptr) {
    std::cerr
        << "Error: Energy parameters not set! Call setEnergyParameters() first."
        << std::endl;
    throw std::runtime_error(
        "FEMHessianAssembler: energy parameters not initialized");
  }

  // Create triplet list for sparse matrix assembly
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(elements.size() * 36); // 6x6 entries per element

  for (size_t elem_idx = 0; elem_idx < elements.size(); elem_idx++) {
    auto &element = elements[elem_idx];

    // Calculate deformation gradient for this element
    // FIXED: Now passing current_points as required
    element.calculate_deformation_gradient(current_points);
    Eigen::Matrix2d F = element.getDeformationGradient();

    // Calculate right Cauchy-Green tensor
    Eigen::Matrix2d C = F.transpose() * F;

    // Perform Lagrange reduction
    lagrange::Result reduction_result = lagrange::reduce(C);
    Eigen::Matrix2d Z = reduction_result.m_matrix;

    // Create acoustic tensor for this element
    AcousticTensor acoustic_tensor(F, C, Z);

    // Get element area for integration weight
    double area = element.getReferenceArea();
    // std::cout<<"area: "<<area<<std::endl;

    // Compute element stiffness matrix
    // (computeEnergyDerivatives is called inside computeElementStiffness)
    Eigen::MatrixXd K_elem =
        computeElementStiffness(acoustic_tensor, element, area);

    // Get global DOF indices for this element
    std::vector<int> global_dof_indices;
    global_dof_indices.reserve(6);

    for (int local_node = 0; local_node < 3; local_node++) {
      int global_node = element.getNodeIndex(local_node);

      // Get DOF mapping for this node
      auto [orig_idx, solver_idx] = dof_mapping[global_node];

      if (solver_idx != -1) { // Free DOF
        // Add x and y DOF indices
        // Convention: [u0, u1, ..., u_n, v0, v1, ..., v_n]
        int num_free_nodes = num_total_dofs / 2;
        global_dof_indices.push_back(solver_idx);                  // x DOF
        global_dof_indices.push_back(solver_idx + num_free_nodes); // y DOF
      } else {
        // Fixed DOF - use -1 as sentinel
        global_dof_indices.push_back(-1);
        global_dof_indices.push_back(-1);
      }
    }

    // Assemble into global matrix using triplets
    for (int i = 0; i < 6; i++) {
      int global_i = global_dof_indices[i];
      if (global_i < 0)
        continue; // Skip fixed DOFs

      for (int j = 0; j < 6; j++) {
        int global_j = global_dof_indices[j];
        if (global_j < 0)
          continue; // Skip fixed DOFs

        if (std::abs(K_elem(i, j)) > 1e-14) { // Only add non-zero entries
          triplets.push_back(
              Eigen::Triplet<double>(global_i, global_j, K_elem(i, j)));
        }
      }
    }
  }

  // Build sparse matrix from triplets
  Eigen::SparseMatrix<double> K_global(num_total_dofs, num_total_dofs);
  K_global.setFromTriplets(triplets.begin(), triplets.end());

  // Make sure matrix is symmetric (numerical check)
  Eigen::SparseMatrix<double> K_transpose = K_global.transpose();
  if ((K_global - K_transpose).norm() > 1e-10) {
    std::cerr << "Warning: Global stiffness matrix is not symmetric!"
              << std::endl;
    std::cerr << "Asymmetry norm: " << (K_global - K_transpose).norm()
              << std::endl;
    exit(0);
  }

  return K_global;
}


Eigen::VectorXd FEMHessianAssembler::assembleGlobalResidual(
    std::vector<ElementTriangle2D> &elements,
    const std::vector<Point2D> &current_points, 
    int num_total_dofs,
    const std::vector<std::pair<int, int>> &dof_mapping) {
  // Safety check: ensure energy parameters were set
  if (strain_calculator == nullptr) {
    std::cerr
        << "Error: Energy parameters not set! Call setEnergyParameters() first."
        << std::endl;
    throw std::runtime_error(
        "FEMHessianAssembler: energy parameters not initialized");
  }

  // Initialize global residual vector
  Eigen::VectorXd R_global = Eigen::VectorXd::Zero(num_total_dofs);

  for (size_t elem_idx = 0; elem_idx < elements.size(); elem_idx++) {
    auto &element = elements[elem_idx];

    // Compute element residual vector (6x1)
    // (calculates deformation gradient inside using current_points)
    Eigen::VectorXd R_elem = computeElementResidual(element, current_points);

    // Get global DOF indices for this element
    std::vector<int> global_dof_indices;
    global_dof_indices.reserve(6);

    for (int local_node = 0; local_node < 3; local_node++) {
      int global_node = element.getNodeIndex(local_node);

      // Get DOF mapping for this node
      auto [orig_idx, solver_idx] = dof_mapping[global_node];

      if (solver_idx != -1) { // Free DOF
        // Add x and y DOF indices
        // Convention: [u0, u1, ..., u_n, v0, v1, ..., v_n]
        int num_free_nodes = num_total_dofs / 2;
        global_dof_indices.push_back(solver_idx);                  // x DOF
        global_dof_indices.push_back(solver_idx + num_free_nodes); // y DOF
      } else {
        // Fixed DOF - use -1 as sentinel
        global_dof_indices.push_back(-1);
        global_dof_indices.push_back(-1);
      }
    }

    // Assemble into global residual vector
    for (int i = 0; i < 6; i++) {
      int global_i = global_dof_indices[i];
      if (global_i < 0)
        continue; // Skip fixed DOFs

      R_global(global_i) += R_elem(i);
    }
  }

  return R_global;
}

Eigen::VectorXd
FEMHessianAssembler::computeElementResidual(const ElementTriangle2D &element,
                                            const std::vector<Point2D>& current_points) {
  // Safety check: ensure energy parameters were set
  if (strain_calculator == nullptr) {
    std::cerr
        << "Error: Energy parameters not set! Call setEnergyParameters() first."
        << std::endl;
    throw std::runtime_error(
        "FEMHessianAssembler: energy parameters not initialized");
  }

  const int num_nodes = 3;                        // Triangular element
  const int spatial_dim = 2;                      // 2D
  const int total_dofs = num_nodes * spatial_dim; // 6 DOFs

  Eigen::VectorXd f_elem = Eigen::VectorXd::Zero(total_dofs);

  // Calculate deformation gradient from current positions
  const_cast<ElementTriangle2D&>(element).calculate_deformation_gradient(current_points);
  
  // Get deformation gradient
  const Eigen::Matrix2d F = element.getDeformationGradient();
  
  // Compute right Cauchy-Green tensor
  Eigen::Matrix2d C = F.transpose() * F;
  
  // Reduce to fundamental domain
  const auto result = lagrange::reduce(C);
  
  // Calculate energy derivative using member variables
  const Eigen::Matrix2d dE_dC = strain_calculator->calculate_derivative(
      result.C_reduced, potential_func_der) / normalisation;
  
  // Compute first Piola-Kirchhoff stress tensor
  // (includes reference area)
  const Eigen::Matrix2d P = 2.0 * F * result.m_matrix * dE_dC * 
                           result.m_matrix.transpose() * element.getReferenceArea();

  // Get shape function derivatives: dN^a/dX_K (3 nodes × 2 coords)
  const Eigen::Matrix<double, 3, 2> &dN_dX = element.getDNdX();

  // Assemble element residual: r^a_i = P_{iK} * dN^a_K
  for (int a = 0; a < num_nodes; a++) {         // Node a
    for (int i = 0; i < spatial_dim; i++) {     // DOF direction i

      double r_a_i = 0.0;

      // Contract over material coordinates K
      for (int K = 0; K < spatial_dim; K++) {
        r_a_i += P(i, K) * dN_dX(a, K);
      }

      // Place in element vector
      int row = a * spatial_dim + i;
      f_elem(row) = r_a_i;
    }
  }

  return f_elem;
}


std::vector<Point2D> FEMHessianAssembler::solveNewtonRaphson(
    std::vector<ElementTriangle2D> &elements,
    const std::vector<Point2D> &initial_points,
    const std::vector<std::pair<int, int>> &dof_mapping,
    int num_total_dofs,
    double tolerance,
    int max_iterations,
    bool verbose) {
  
  // Safety check
  if (strain_calculator == nullptr) {
    std::cerr << "Error: Energy parameters not set! Call setEnergyParameters() first." << std::endl;
    throw std::runtime_error("FEMHessianAssembler: energy parameters not initialized");
  }

  // Initialize current positions
  std::vector<Point2D> current_points = initial_points;
  
  // Number of free DOFs per direction
  int num_free_nodes = num_total_dofs / 2;

  if (verbose) {
    std::cout << "\n=== Newton-Raphson Solver ===" << std::endl;
    std::cout << "Number of free DOFs: " << num_total_dofs << std::endl;
    std::cout << "Tolerance: " << tolerance << std::endl;
    std::cout << "Max iterations: " << max_iterations << std::endl;
  }

  // Newton-Raphson iteration
  for (int iter = 0; iter < max_iterations; iter++) {
    
    // 1. Assemble global residual vector R
    Eigen::VectorXd R = assembleGlobalResidual(elements, current_points, 
                                               num_total_dofs, dof_mapping);
    
    // 2. Check convergence
    double residual_norm = R.norm();
    
    if (verbose) {
      std::cout << "Iteration " << iter << ": ||R|| = " << residual_norm << std::endl;
    }
    
    if (residual_norm < tolerance) {
      if (verbose) {
        std::cout << "Converged in " << iter << " iterations!" << std::endl;
      }
      return current_points;
    }
    
    // 3. Assemble global stiffness matrix K (Jacobian)
    Eigen::SparseMatrix<double> K = assembleGlobalStiffness(elements, current_points,
                                                            num_total_dofs, dof_mapping);
    
    // 4. Solve linear system: K * delta_u = -R
    Eigen::VectorXd delta_u;
    
    // ============================================================
    // CHOOSE YOUR SOLVER (uncomment one)
    // ============================================================
    
    // OPTION 1: SparseLU (general purpose, slower)
    // {
    //   Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    //   solver.analyzePattern(K);
    //   solver.factorize(K);
    //   if (solver.info() != Eigen::Success) {
    //     std::cerr << "Error: SparseLU factorization failed at iteration " << iter << std::endl;
    //     return current_points;
    //   }
    //   delta_u = solver.solve(-R);
    //   if (solver.info() != Eigen::Success) {
    //     std::cerr << "Error: SparseLU solve failed at iteration " << iter << std::endl;
    //     return current_points;
    //   }
    // }
    
    // OPTION 2: SimplicialLDLT (fast for symmetric matrices) - RECOMMENDED
    {
      Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
      solver.analyzePattern(K);
      solver.factorize(K);
      if (solver.info() != Eigen::Success) {
        std::cerr << "Error: SimplicialLDLT factorization failed at iteration " << iter << std::endl;
        return current_points;
      }
      delta_u = solver.solve(-R);
      if (solver.info() != Eigen::Success) {
        std::cerr << "Error: SimplicialLDLT solve failed at iteration " << iter << std::endl;
        return current_points;
      }
    }
    
    // OPTION 3: ConjugateGradient (iterative, good for large systems)
    // {
    //   Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
    //   solver.setMaxIterations(3000);
    //   solver.setTolerance(1e-6);
    //   solver.compute(K);
    //   if (solver.info() != Eigen::Success) {
    //     std::cerr << "Error: CG initialization failed at iteration " << iter << std::endl;
    //     return current_points;
    //   }
    //   delta_u = solver.solve(-R);
    //   if (solver.info() != Eigen::Success) {
    //     std::cerr << "Error: CG solve failed at iteration " << iter << std::endl;
    //     return current_points;
    //   }
    //   if (verbose) {
    //     std::cout << "  CG iterations: " << solver.iterations() 
    //               << ", error: " << solver.error() << std::endl;
    //   }
    // }
    
    // ============================================================
    
    // 5. Update positions: u_new = u_old + delta_u
    for (size_t node_idx = 0; node_idx < current_points.size(); node_idx++) {
      auto [orig_idx, solver_idx] = dof_mapping[node_idx];
      
      if (solver_idx != -1) { // Free DOF
        // Extract displacement from solution vector
        double delta_x = delta_u(solver_idx);
        double delta_y = delta_u(solver_idx + num_free_nodes);
        
        // Update position using Eigen::Vector2d interface
        current_points[node_idx].coord.x() += 0.01*delta_x;
        current_points[node_idx].coord.y() += 0.01*delta_y;
      }
      // Fixed DOFs: do nothing
    }
    
    if (verbose && iter % 1 == 0) {
      std::cout << "  |delta_u| = " << delta_u.norm() << std::endl;
    }
  }
  
  // Did not converge
  std::cerr << "Warning: Newton-Raphson did not converge in " << max_iterations << " iterations" << std::endl;
  Eigen::VectorXd R_final = assembleGlobalResidual(elements, current_points,
                                                   num_total_dofs, dof_mapping);
  std::cerr << "Final residual norm: " << R_final.norm() << std::endl;
  
  return current_points;
}
















/// @brief ////////////////////////////////////////////////////////

EigenResults FEMHessianAssembler::computeSmallestEigenvalues(
    const Eigen::SparseMatrix<double> &K_global, int N) {
  EigenResults results;

  int n = K_global.rows();

  // Check if N is valid
  if (N <= 0 || N > n) {
    std::cerr << "Error: N must be between 1 and " << n << std::endl;
    return results;
  }

  std::cout << "Computing " << N << " smallest eigenvalues..." << std::endl;
  std::cout << "Matrix size: " << n << " x " << n << std::endl;

  // For moderate size matrices (< 5000 DOFs), convert to dense
  if (n <= 5000) {
    std::cout << "Using dense eigenvalue solver..." << std::endl;

    // Convert sparse to dense
    Eigen::MatrixXd K_dense = Eigen::MatrixXd(K_global);

    // Compute all eigenvalues and eigenvectors
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(K_dense);

    if (eigensolver.info() != Eigen::Success) {
      std::cerr << "Error: Eigenvalue computation failed!" << std::endl;
      return results;
    }

    // Extract N smallest eigenvalues and corresponding eigenvectors
    // Eigenvalues are already sorted in increasing order
    results.eigenvalues = eigensolver.eigenvalues().head(N);
    results.eigenvectors = eigensolver.eigenvectors().leftCols(N);
    results.num_computed = N;

    std::cout << "Successfully computed " << N << " eigenvalues" << std::endl;
    std::cout << "Smallest eigenvalue: " << results.eigenvalues(0) << std::endl;
    if (N > 1) {
      std::cout << "Largest of the " << N
                << " smallest: " << results.eigenvalues(N - 1) << std::endl;
    }

  } else {
    std::cerr << "Warning: Matrix too large (" << n
              << " DOFs) for dense solver." << std::endl;
    std::cerr << "Use computeSmallestEigenvaluesIterative() instead."
              << std::endl;
  }

  return results;
}

EigenResults FEMHessianAssembler::computeSmallestEigenvaluesIterative(
    const Eigen::SparseMatrix<double> &K_global, int N, double shift) {
  EigenResults results;

  int n = K_global.rows();

  // Check if N is valid
  if (N <= 0 || N > n) {
    std::cerr << "Error: N must be between 1 and " << n << std::endl;
    return results;
  }

  std::cout << "Computing " << N
            << " smallest eigenvalues using iterative method..." << std::endl;
  std::cout << "Matrix size: " << n << " x " << n << std::endl;
  std::cout << "Shift value: " << shift << std::endl;

  // Create shifted matrix: K_shifted = K - shift * I
  Eigen::SparseMatrix<double> K_shifted = K_global;
  for (int i = 0; i < n; i++) {
    K_shifted.coeffRef(i, i) -= shift;
  }

  // Set up sparse linear solver for (K - shift*I)
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(K_shifted);

  if (solver.info() != Eigen::Success) {
    std::cerr << "Error: Failed to factorize shifted matrix!" << std::endl;
    std::cerr
        << "Try a different shift value or use computeSmallestEigenvalues()"
        << std::endl;
    return results;
  }

  std::cout << "Note: Full iterative eigenvalue solver (Arnoldi/Lanczos) "
               "requires additional libraries."
            << std::endl;
  std::cout
      << "For now, falling back to dense solver for moderate-sized problems."
      << std::endl;

  // For a complete implementation, you would use libraries like ARPACK or
  // Spectra here For now, provide a basic implementation that works for
  // moderate sizes

  if (n <= 5000) {
    Eigen::MatrixXd K_dense = Eigen::MatrixXd(K_global);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(K_dense);

    if (eigensolver.info() == Eigen::Success) {
      results.eigenvalues = eigensolver.eigenvalues().head(N);
      results.eigenvectors = eigensolver.eigenvectors().leftCols(N);
      results.num_computed = N;

      std::cout << "Successfully computed " << N << " eigenvalues" << std::endl;
    }
  } else {
    std::cerr << "Matrix too large for current implementation." << std::endl;
    std::cerr << "Consider using external libraries like ARPACK or Spectra."
              << std::endl;
  }

  return results;
}

#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <Spectra/SymEigsShiftSolver.h>

EigenResults FEMHessianAssembler::computeSmallestEigenvaluesIterative_spectra(
    const Eigen::SparseMatrix<double> &K_global, int N, double shift) {
  EigenResults results;
  int n = K_global.rows();

  // Check if N is valid
  if (N <= 0 || N > n) {
    std::cerr << "Error: N must be between 1 and " << n << std::endl;
    return results;
  }

  std::cout << "Computing " << N
            << " smallest eigenvalues using Spectra with shift-invert..."
            << std::endl;
  std::cout << "Matrix size: " << n << " x " << n << std::endl;

  // Auto-select shift if not provided
  if (shift == -1) {
    // Estimate shift: use small negative value
    shift = -0.01 * K_global.diagonal().cwiseAbs().maxCoeff();
    std::cout << "Auto-selected shift: " << shift << std::endl;
  } else {
    std::cout << "Using shift: " << shift << std::endl;
  }

  // Number of Lanczos vectors (ncv must be > N)
  //int ncv = std::min(2 * N + 1, n - 1);
  // int ncv = std::min(4 * N, n - 1);  // Give it plenty of room
  int ncv = std::min(std::max(2 * N + 1, 20), n - 1);

  std::cout << "Spectra: N=" << N << ", ncv=" << ncv << ", shift=" << shift
            << std::endl;

  try {
    std::cout << "Performing sparse LU factorization of (K - σI)..."
              << std::endl;

    // Construct shift-invert matrix operation object
    // This performs LU factorization of (K - shift*I)
    Spectra::SparseSymShiftSolve<double> op(K_global);

    // Construct symmetric shift-invert eigen solver
    Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<double>> eigs(
        op, N, ncv, shift);

    std::cout << "LU factorization complete. Computing eigenvalues..."
              << std::endl;

    // Initialize and compute
    eigs.init();

    // Use LargestMagn because shift-invert transforms smallest to largest
    int nconv = eigs.compute(Spectra::SortRule::LargestMagn, 1000, 1e-6);

    std::cout << "Spectra: Converged " << nconv << " eigenvalues" << std::endl;

    // Check convergence
    if (eigs.info() != Spectra::CompInfo::Successful) {
      std::cerr << "Error: Spectra failed to converge!" << std::endl;
      std::cerr << "Info code: " << static_cast<int>(eigs.info()) << std::endl;
      return results;
    }

    // Retrieve results
    // Eigenvalues are returned in order relative to the shift
    results.eigenvalues = eigs.eigenvalues();
    results.eigenvectors = eigs.eigenvectors();
    results.num_computed = nconv;

    // Sort eigenvalues in ascending order (they come sorted by distance from
    // shift)
    std::vector<std::pair<double, int>> sorted_indices;
    for (int i = 0; i < nconv; i++) {
      sorted_indices.push_back({results.eigenvalues(i), i});
    }
    std::sort(sorted_indices.begin(), sorted_indices.end());

    // Reorder eigenvalues and eigenvectors
    Eigen::VectorXd sorted_eigenvalues(nconv);
    Eigen::MatrixXd sorted_eigenvectors(n, nconv);

    for (int i = 0; i < nconv; i++) {
      sorted_eigenvalues(i) = sorted_indices[i].first;
      sorted_eigenvectors.col(i) =
          results.eigenvectors.col(sorted_indices[i].second);
    }

    results.eigenvalues = sorted_eigenvalues;
    results.eigenvectors = sorted_eigenvectors;

    // Print results
    std::cout << "Successfully computed " << results.num_computed
              << " eigenvalues." << std::endl;
    if (results.num_computed > 0) {
      std::cout << "Smallest eigenvalue: " << results.eigenvalues(0)
                << std::endl;
      if (results.num_computed > 1) {
        std::cout << "Largest of computed eigenvalues: "
                  << results.eigenvalues(results.num_computed - 1) << std::endl;
      }
    }

  } catch (const std::exception &e) {
    std::cerr << "Exception in Spectra shift-invert computation: " << e.what()
              << std::endl;
  }

  return results;
}

#include <Accelerate/Accelerate.h>

EigenResults FEMHessianAssembler::computeSmallestEigenvalues_Accelerate(
    const Eigen::SparseMatrix<double> &K_global, int N) {
  EigenResults results;
  int n = K_global.rows();

  // Check if N is valid
  if (N <= 0 || N > n) {
    std::cerr << "Error: N must be between 1 and " << n << std::endl;
    return results;
  }

  std::cout << "Computing " << N
            << " smallest eigenvalues using Apple Accelerate..." << std::endl;
  std::cout << "Matrix size: " << n << " x " << n << std::endl;

  try {
    // Convert sparse to dense (column-major for LAPACK)
    std::vector<double> K_dense(n * n, 0.0);
    for (int k = 0; k < K_global.outerSize(); ++k) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(K_global, k); it;
           ++it) {
        K_dense[it.col() * n + it.row()] = it.value(); // Column-major
      }
    }

    std::cout << "Converted to dense matrix" << std::endl;

    // Prepare output arrays
    std::vector<double> eigenvalues(n);

    // LAPACK parameters for dsyevd (divide-and-conquer, faster than dsyev)
    char jobz = 'V'; // Compute eigenvalues and eigenvectors
    char uplo = 'L'; // Lower triangle of matrix
    __CLPK_integer n_lapack = n;
    __CLPK_integer lda = n;
    __CLPK_integer info = 0;

    // Query optimal workspace size
    __CLPK_integer lwork = -1;
    __CLPK_integer liwork = -1;
    double work_query;
    __CLPK_integer iwork_query;

    std::cout << "Querying optimal workspace..." << std::endl;
    dsyevd_(&jobz, &uplo, &n_lapack, K_dense.data(), &lda, eigenvalues.data(),
            &work_query, &lwork, &iwork_query, &liwork, &info);

    if (info != 0) {
      std::cerr << "Error in workspace query: info = " << info << std::endl;
      return results;
    }

    // Allocate optimal workspace
    lwork = static_cast<__CLPK_integer>(work_query);
    liwork = iwork_query;
    std::vector<double> work(lwork);
    std::vector<__CLPK_integer> iwork(liwork);

    std::cout << "Workspace allocated: lwork=" << lwork << ", liwork=" << liwork
              << std::endl;
    std::cout << "Computing eigenvalues..." << std::endl;

    // Compute eigenvalues and eigenvectors
    dsyevd_(&jobz, &uplo, &n_lapack, K_dense.data(), &lda, eigenvalues.data(),
            work.data(), &lwork, iwork.data(), &liwork, &info);

    if (info != 0) {
      if (info < 0) {
        std::cerr << "Error: Invalid argument at position " << -info
                  << std::endl;
      } else {
        std::cerr << "Error: Algorithm failed to converge. info = " << info
                  << std::endl;
      }
      return results;
    }

    std::cout << "Successfully computed all eigenvalues" << std::endl;

    // Extract smallest N eigenvalues and eigenvectors
    // Eigenvalues are already sorted in ascending order
    results.eigenvalues.resize(N);
    results.eigenvectors.resize(n, N);

    for (int i = 0; i < N; ++i) {
      results.eigenvalues(i) = eigenvalues[i];

      // Copy eigenvector (stored column-major in K_dense)
      for (int j = 0; j < n; ++j) {
        results.eigenvectors(j, i) = K_dense[i * n + j];
      }
    }

    results.num_computed = N;

    // Print results
    std::cout << "Successfully extracted " << N << " smallest eigenvalues."
              << std::endl;
    std::cout << "Smallest eigenvalue: " << results.eigenvalues(0) << std::endl;
    if (N > 1) {
      std::cout << "Largest of extracted eigenvalues: "
                << results.eigenvalues(N - 1) << std::endl;
    }

  } catch (const std::exception &e) {
    std::cerr << "Exception in Accelerate computation: " << e.what()
              << std::endl;
  }

  return results;
}



#include <armadillo>

EigenResults FEMHessianAssembler::computeSmallestEigenvaluesIterative_armadillo(
    const Eigen::SparseMatrix<double> &K_global, int N, double shift) {
  EigenResults results;
  int n = K_global.rows();

  // Check if N is valid
  if (N <= 0 || N > n) {
    std::cerr << "Error: N must be between 1 and " << n << std::endl;
    return results;
  }

  std::cout << "Computing " << N
            << " smallest eigenvalues using Armadillo..."
            << std::endl;
  std::cout << "Matrix size: " << n << " x " << n << std::endl;

  // Auto-select shift if not provided
  if (shift == -1) {
    shift = 0.0; // Better default for smallest eigenvalues
    std::cout << "Auto-selected shift: " << shift << std::endl;
  } else {
    std::cout << "Using shift: " << shift << std::endl;
  }

  try {
    // Convert Eigen sparse matrix to Armadillo sparse matrix
    std::cout << "Converting Eigen sparse matrix to Armadillo..." << std::endl;
    
    // Create vectors for Armadillo batch constructor
    arma::umat locations(2, K_global.nonZeros());
    arma::vec values(K_global.nonZeros());
    
    int idx = 0;
    for (int k = 0; k < K_global.outerSize(); ++k) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(K_global, k); it; ++it) {
        locations(0, idx) = it.row();
        locations(1, idx) = it.col();
        values(idx) = it.value();
        idx++;
      }
    }
    
    arma::sp_mat K_arma(locations, values, n, n);
    
    std::cout << "Conversion complete. Non-zeros: " << K_arma.n_nonzero 
              << std::endl;

    // Prepare for eigenvalue computation
    arma::vec eigval;
    arma::mat eigvec;
    
    std::cout << "Computing eigenvalues with Armadillo eigs_sym..." << std::endl;
    
    // Set up eigs_opts for tolerance and max iterations
    arma::eigs_opts opts;
    opts.tol = 1e-10;
    opts.maxiter = 10000;
    
    // Compute smallest eigenvalues
    bool success;
    if (shift == 0.0) {
      // Direct computation of smallest magnitude eigenvalues
      // "sm" = smallest magnitude
      success = arma::eigs_sym(eigval, eigvec, K_arma, N, "sm", opts);
    } else {
      // Use shift-invert mode
      // Pass shift as sigma parameter
      success = arma::eigs_sym(eigval, eigvec, K_arma, N, shift, opts);
    }

    if (!success) {
      std::cerr << "Error: Armadillo eigs_sym failed to converge!" << std::endl;
      return results;
    }

    std::cout << "Armadillo: Successfully computed " << eigval.n_elem 
              << " eigenvalues" << std::endl;

    // Convert results back to Eigen format
    results.num_computed = eigval.n_elem;
    results.eigenvalues = Eigen::VectorXd(results.num_computed);
    results.eigenvectors = Eigen::MatrixXd(n, results.num_computed);

    for (size_t i = 0; i < eigval.n_elem; ++i) {
      results.eigenvalues(i) = eigval(i);
      for (int j = 0; j < n; ++j) {
        results.eigenvectors(j, i) = eigvec(j, i);
      }
    }

    // Sort eigenvalues in ascending order
    std::vector<std::pair<double, int>> sorted_indices;
    for (int i = 0; i < results.num_computed; i++) {
      sorted_indices.push_back({results.eigenvalues(i), i});
    }
    std::sort(sorted_indices.begin(), sorted_indices.end());

    // Reorder
    Eigen::VectorXd sorted_eigenvalues(results.num_computed);
    Eigen::MatrixXd sorted_eigenvectors(n, results.num_computed);

    for (int i = 0; i < results.num_computed; i++) {
      sorted_eigenvalues(i) = sorted_indices[i].first;
      sorted_eigenvectors.col(i) = results.eigenvectors.col(sorted_indices[i].second);
    }

    results.eigenvalues = sorted_eigenvalues;
    results.eigenvectors = sorted_eigenvectors;

    // Print results
    std::cout << "Successfully computed " << results.num_computed
              << " eigenvalues." << std::endl;
    if (results.num_computed > 0) {
      std::cout << "Smallest eigenvalue: " << results.eigenvalues(0) << std::endl;
      if (results.num_computed > 1) {
        std::cout << "Largest of computed eigenvalues: "
                  << results.eigenvalues(results.num_computed - 1) << std::endl;
      }
    }

  } catch (const std::exception &e) {
    std::cerr << "Exception in Armadillo computation: " << e.what() << std::endl;
  }

  return results;
}


#include <fstream>
#include <iomanip>

#include <fstream>
#include <iomanip>
#include <sstream>

void FEMHessianAssembler::exportSingleModeToVTK(
    const std::string &filename, const std::vector<Point2D> &current_points,
    const Eigen::VectorXd &eigenvector, double eigenvalue, int mode_number,
    double scale_factor, const std::vector<std::pair<int, int>> &dof_mapping,
    const std::vector<ElementTriangle2D> &elements) {
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error: Could not open file " << filename << std::endl;
    return;
  }

  int num_nodes = current_points.size();
  int num_free_nodes = eigenvector.size() / 2;

  // Reconstruct full displacement field from eigenvector
  // eigenvector format: [u0, u1, ..., u_n, v0, v1, ..., v_n]
  std::vector<Eigen::Vector2d> displacements(num_nodes,
                                             Eigen::Vector2d::Zero());

  for (int node = 0; node < num_nodes; node++) {
    auto [orig_idx, solver_idx] = dof_mapping[node];

    if (solver_idx != -1) {                                // Free DOF
      double u = eigenvector(solver_idx);                  // x displacement
      double v = eigenvector(solver_idx + num_free_nodes); // y displacement
      displacements[node] = Eigen::Vector2d(u, v);
    }
    // Fixed DOFs remain zero
  }

  // Normalize displacement field for better visualization
  double max_displacement = 0.0;
  for (const auto &disp : displacements) {
    max_displacement = std::max(max_displacement, disp.norm());
  }

  if (max_displacement > 1e-12) {
    for (auto &disp : displacements) {
      disp /= max_displacement;
      disp *= scale_factor;
    }
  }

  // Filter elements: only keep those with ZERO translation (fundamental domain
  // only)
  std::vector<int> fundamental_domain_elements;
  for (size_t i = 0; i < elements.size(); i++) {
    const auto &elem = elements[i];
    if (!elem.isInitialized())
      continue;

    // Check if all nodes have zero translation
    bool is_fundamental = true;
    for (int j = 0; j < 3; j++) {
      Eigen::Vector2d translation = elem.getTranslation(j);
      if (translation.norm() > 1e-10) { // Non-zero translation
        is_fundamental = false;
        break;
      }
    }

    if (is_fundamental) {
      fundamental_domain_elements.push_back(i);
    }
  }

  // Write VTK header
  file << "# vtk DataFile Version 3.0\n";
  file << "Eigenmode " << mode_number << " (lambda = " << std::scientific
       << eigenvalue << ")\n";
  file << "ASCII\n";
  file << "DATASET UNSTRUCTURED_GRID\n";

  // Write points (current configuration)
  file << "POINTS " << num_nodes << " double\n";
  for (int i = 0; i < num_nodes; i++) {
    file << std::setprecision(10) << current_points[i].coord.x() << " "
         << current_points[i].coord.y() << " " << 0.0 << "\n";
  }

  // Write cells (only fundamental domain - no periodic images)
  file << "\nCELLS " << fundamental_domain_elements.size() << " "
       << (fundamental_domain_elements.size() * 4) << "\n";
  for (int elem_idx : fundamental_domain_elements) {
    const auto &elem = elements[elem_idx];
    file << "3 " << elem.getNodeIndex(0) << " " << elem.getNodeIndex(1) << " "
         << elem.getNodeIndex(2) << "\n";
  }

  file << "\nCELL_TYPES " << fundamental_domain_elements.size() << "\n";
  for (size_t i = 0; i < fundamental_domain_elements.size(); i++) {
    file << "5\n"; // VTK_TRIANGLE
  }

  // Write point data
  file << "\nPOINT_DATA " << num_nodes << "\n";

  // Write eigenvector as displacement vectors (for Glyph/Arrow visualization)
  file << "VECTORS eigenvector double\n";
  for (int i = 0; i < num_nodes; i++) {
    file << std::setprecision(10) << displacements[i](0) << " "
         << displacements[i](1) << " " << 0.0 << "\n";
  }

  // Write displacement magnitude as scalar (for coloring)
  file << "\nSCALARS eigenvector_magnitude double 1\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_nodes; i++) {
    file << std::setprecision(10) << displacements[i].norm() << "\n";
  }

  // Write x and y components separately (useful for analysis)
  file << "\nSCALARS eigenvector_x double 1\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_nodes; i++) {
    file << std::setprecision(10) << displacements[i](0) << "\n";
  }

  file << "\nSCALARS eigenvector_y double 1\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_nodes; i++) {
    file << std::setprecision(10) << displacements[i](1) << "\n";
  }

  // Write deformed positions as vectors (for Warp By Vector)
  file << "\nVECTORS deformed_position double\n";
  for (int i = 0; i < num_nodes; i++) {
    file << std::setprecision(10)
         << (current_points[i].coord.x() + displacements[i](0)) << " "
         << (current_points[i].coord.y() + displacements[i](1)) << " " << 0.0
         << "\n";
  }

  file.close();
  std::cout << "  Exported mode " << mode_number << " to " << filename
            << std::endl;
}

void FEMHessianAssembler::exportEigenvectorsToVTK(
    const std::string &filename, const std::vector<Point2D> &current_points,
    const EigenResults &eigen_results, const std::vector<int> &mode_indices,
    double scale_factor, const std::vector<std::pair<int, int>> &dof_mapping,
    const std::vector<ElementTriangle2D> &elements) {
  if (eigen_results.num_computed == 0) {
    std::cerr << "Error: No eigenvalues computed!" << std::endl;
    return;
  }

  if (dof_mapping.empty()) {
    std::cerr << "Error: DOF mapping required for VTK export!" << std::endl;
    return;
  }

  // Create vtk_eigen directory
  std::filesystem::path filepath(filename);
  std::filesystem::path parent_dir = filepath.parent_path();
  std::filesystem::path vtk_eigen_dir = parent_dir / "vtk_eigen";

  // Create directory if it doesn't exist
  if (!std::filesystem::exists(vtk_eigen_dir)) {
    std::filesystem::create_directories(vtk_eigen_dir);
    std::cout << "Created directory: " << vtk_eigen_dir << std::endl;
  }

  // Get base filename without directory
  std::string base_filename = filepath.filename().string();

  std::cout << "\n=== Exporting Eigenmodes to VTK ===" << std::endl;
  std::cout << "Output directory: " << vtk_eigen_dir << std::endl;
  std::cout << "Base filename: " << base_filename << std::endl;
  std::cout << "Scale factor: " << scale_factor << std::endl;
  std::cout << "Number of modes to export: " << mode_indices.size()
            << std::endl;

  // Count fundamental domain elements (for reporting)
  int fundamental_count = 0;
  for (size_t i = 0; i < elements.size(); i++) {
    const auto &elem = elements[i];
    if (!elem.isInitialized())
      continue;

    bool is_fundamental = true;
    for (int j = 0; j < 3; j++) {
      Eigen::Vector2d translation = elem.getTranslation(j);
      if (translation.norm() > 1e-10) {
        is_fundamental = false;
        break;
      }
    }
    if (is_fundamental)
      fundamental_count++;
  }

  std::cout << "Total elements: " << elements.size()
            << ", Fundamental domain elements: " << fundamental_count
            << std::endl;

  for (int mode_idx : mode_indices) {
    if (mode_idx < 0 || mode_idx >= eigen_results.num_computed) {
      std::cerr << "Warning: Mode index " << mode_idx << " out of range [0, "
                << eigen_results.num_computed << ")" << std::endl;
      continue;
    }

    // Create filename for this mode inside vtk_eigen directory
    std::stringstream ss;
    ss << base_filename << "_mode_" << std::setfill('0') << std::setw(3)
       << mode_idx << ".vtk";
    std::filesystem::path mode_filepath = vtk_eigen_dir / ss.str();
    std::string mode_filename = mode_filepath.string();

    // Export this mode (automatically handles PBC filtering)
    exportSingleModeToVTK(mode_filename, current_points,
                          eigen_results.eigenvectors.col(mode_idx),
                          eigen_results.eigenvalues(mode_idx), mode_idx,
                          scale_factor, dof_mapping, elements);
  }

  std::cout << "\nSuccessfully exported " << mode_indices.size()
            << " eigenmodes" << std::endl;
  std::cout << "All modes show ONLY fundamental domain (no periodic images)"
            << std::endl;
  std::cout << "\n=== Paraview Visualization Guide ===" << std::endl;
  std::cout << "1. Open the .vtk files in Paraview" << std::endl;
  std::cout << "2. For ARROWS: Apply 'Glyph' filter" << std::endl;
  std::cout << "   - Glyph Type: Arrow (or 2D Glyph)" << std::endl;
  std::cout << "   - Orientation Array: eigenvector" << std::endl;
  std::cout << "   - Scale Array: eigenvector_magnitude" << std::endl;
  std::cout << "   - Adjust Scale Factor to see arrows clearly" << std::endl;
  std::cout << "3. For DEFORMATION: Apply 'Warp By Vector' filter" << std::endl;
  std::cout << "   - Vectors: eigenvector" << std::endl;
  std::cout << "   - Adjust Scale Factor to see deformation" << std::endl;
  std::cout << "4. Color by: eigenvector_magnitude" << std::endl;
}

double FEMHessianAssembler::computeParticipationRatio(
    const Eigen::VectorXd &eigenvector,
    const std::vector<std::pair<int, int>> &dof_mapping, int num_nodes) {
  if (eigenvector.size() == 0) {
    std::cerr << "Error: Empty eigenvector!" << std::endl;
    return 0.0;
  }

  int num_free_nodes = eigenvector.size() / 2;

  // Reconstruct displacement magnitudes at each node
  std::vector<double> displacement_magnitude(num_nodes, 0.0);

  for (int node = 0; node < num_nodes; node++) {
    auto [orig_idx, solver_idx] = dof_mapping[node];

    if (solver_idx != -1) {                                // Free DOF
      double u = eigenvector(solver_idx);                  // x displacement
      double v = eigenvector(solver_idx + num_free_nodes); // y displacement
      displacement_magnitude[node] = std::sqrt(u * u + v * v);
    }
    // Fixed DOFs have zero displacement
  }

  // Compute participation ratio: P = (Σu²)² / (N * Σu⁴)
  double sum_u2 = 0.0;
  double sum_u4 = 0.0;

  for (int i = 0; i < num_nodes; i++) {
    double u2 = displacement_magnitude[i] * displacement_magnitude[i];
    sum_u2 += u2;
    sum_u4 += u2 * u2;
  }

  // Avoid division by zero
  if (sum_u4 < 1e-20) {
    return 0.0;
  }

  double P = (sum_u2 * sum_u2) / (num_nodes * sum_u4);

  return P;
}

ParticipationAnalysis FEMHessianAssembler::analyzeParticipationRatios(
    const EigenResults &eigen_results,
    const std::vector<std::pair<int, int>> &dof_mapping, int num_nodes,
    double threshold) {
  ParticipationAnalysis analysis;
  analysis.threshold = threshold;

  if (eigen_results.num_computed == 0) {
    std::cerr << "Error: No eigenmodes to analyze!" << std::endl;
    return analysis;
  }

  std::cout << "\n=== Analyzing Participation Ratios ===" << std::endl;
  std::cout << "Number of modes: " << eigen_results.num_computed << std::endl;
  std::cout << "Number of nodes: " << num_nodes << std::endl;
  std::cout << "Classification threshold: " << threshold << std::endl;

  analysis.participation_ratios.reserve(eigen_results.num_computed);
  analysis.is_extended.reserve(eigen_results.num_computed);
  analysis.localization_lengths.reserve(eigen_results.num_computed);

  int num_extended = 0;
  int num_localized = 0;

  for (int mode = 0; mode < eigen_results.num_computed; mode++) {
    // Compute participation ratio
    double P = computeParticipationRatio(eigen_results.eigenvectors.col(mode),
                                         dof_mapping, num_nodes);

    analysis.participation_ratios.push_back(P);

    // Classify as extended or localized
    bool extended = (P >= threshold);
    analysis.is_extended.push_back(extended);

    // Compute localization length (effective number of participating nodes)
    double L_eff = P * num_nodes;
    analysis.localization_lengths.push_back(L_eff);

    if (extended) {
      num_extended++;
    } else {
      num_localized++;
    }
  }

  std::cout << "\nResults:" << std::endl;
  std::cout << "  Extended modes (phonons): " << num_extended << std::endl;
  std::cout << "  Localized modes: " << num_localized << std::endl;
  std::cout << "  Fraction extended: "
            << static_cast<double>(num_extended) / eigen_results.num_computed
            << std::endl;

  return analysis;
}

void FEMHessianAssembler::exportParticipationAnalysis(
    const std::string &filename, const ParticipationAnalysis &analysis,
    const EigenResults &eigen_results) {
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error: Could not open file " << filename << std::endl;
    return;
  }

  // Write header
  file << "# Participation Ratio Analysis\n";
  file << "# Threshold for extended modes: " << analysis.threshold << "\n";
  file << "# Columns: Mode_Index, Eigenvalue, Participation_Ratio, "
          "Localization_Length, Is_Extended\n";
  file << "# P ≈ 1: Extended mode (phonon), P ≈ 0: Localized mode\n";
  file << "#\n";
  file << std::setw(12) << "Mode" << std::setw(20) << "Eigenvalue"
       << std::setw(20) << "P_ratio" << std::setw(20) << "L_eff"
       << std::setw(15) << "Type"
       << "\n";

  // Write data
  for (size_t i = 0; i < analysis.participation_ratios.size(); i++) {
    file << std::setw(12) << i << std::setw(20) << std::setprecision(10)
         << eigen_results.eigenvalues(i) << std::setw(20)
         << std::setprecision(6) << analysis.participation_ratios[i]
         << std::setw(20) << std::setprecision(2)
         << analysis.localization_lengths[i] << std::setw(15)
         << (analysis.is_extended[i] ? "Extended" : "Localized") << "\n";
  }

  file.close();
  std::cout << "Participation analysis exported to: " << filename << std::endl;
}

void FEMHessianAssembler::printParticipationSummary(
    const ParticipationAnalysis &analysis, const EigenResults &eigen_results) {
  std::cout << "\n=== Participation Ratio Summary ===" << std::endl;
  std::cout << std::setw(8) << "Mode" << std::setw(16) << "Eigenvalue"
            << std::setw(12) << "P_ratio" << std::setw(12) << "L_eff"
            << std::setw(12) << "Type" << std::endl;
  std::cout << std::string(60, '-') << std::endl;

  for (size_t i = 0; i < analysis.participation_ratios.size(); i++) {
    std::cout << std::setw(8) << i << std::setw(16) << std::scientific
              << std::setprecision(4) << eigen_results.eigenvalues(i)
              << std::setw(12) << std::fixed << std::setprecision(4)
              << analysis.participation_ratios[i] << std::setw(12) << std::fixed
              << std::setprecision(1) << analysis.localization_lengths[i]
              << std::setw(12)
              << (analysis.is_extended[i] ? "Extended" : "Localized")
              << std::endl;
  }

  // Statistics
  std::cout << std::string(60, '-') << std::endl;

  int num_extended = std::count(analysis.is_extended.begin(),
                                analysis.is_extended.end(), true);
  int num_localized = analysis.is_extended.size() - num_extended;

  std::cout << "\nStatistics:" << std::endl;
  std::cout << "  Total modes analyzed: "
            << analysis.participation_ratios.size() << std::endl;
  std::cout << "  Extended modes (phonons): " << num_extended << std::endl;
  std::cout << "  Localized modes: " << num_localized << std::endl;
  std::cout << "  Fraction extended: "
            << static_cast<double>(num_extended) /
                   analysis.participation_ratios.size()
            << std::endl;

  // Average participation ratios
  double avg_P_extended = 0.0;
  double avg_P_localized = 0.0;

  for (size_t i = 0; i < analysis.participation_ratios.size(); i++) {
    if (analysis.is_extended[i]) {
      avg_P_extended += analysis.participation_ratios[i];
    } else {
      avg_P_localized += analysis.participation_ratios[i];
    }
  }

  if (num_extended > 0) {
    avg_P_extended /= num_extended;
    std::cout << "  Average P (extended): " << avg_P_extended << std::endl;
  }
  if (num_localized > 0) {
    avg_P_localized /= num_localized;
    std::cout << "  Average P (localized): " << avg_P_localized << std::endl;
  }

  std::cout << std::endl;
}

#include <fstream>
#include <iomanip>
#include <sstream>

void FEMHessianAssembler::exportNonRigidEigenvalues(
    const std::string &filename, const EigenResults &eigen_results,
    int num_rigid_body) {
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error: Could not open file " << filename << std::endl;
    return;
  }

  file << "# Non-rigid eigenvalues analysis\n";
  file << "# Number of rigid body modes skipped: " << num_rigid_body << "\n";
  file << "# Columns: Mode_Index, Eigenvalue(λ), Frequency(ω=√λ), "
          "Lambda_Squared(λ²)\n";
  file << "#\n";
  file << std::setw(12) << "Mode" << std::setw(20) << "Lambda" << std::setw(20)
       << "Omega" << std::setw(20) << "Lambda_Squared"
       << "\n";

  for (int i = num_rigid_body; i < eigen_results.num_computed; i++) {
    double lambda = eigen_results.eigenvalues(i);
    double omega = std::sqrt(std::abs(lambda));
    if (lambda < 0)
      omega = -omega; // Keep sign for negative eigenvalues
    double lambda_sq = lambda * lambda;

    file << std::setw(12) << i << std::setw(20) << std::scientific
         << std::setprecision(10) << lambda << std::setw(20) << std::scientific
         << std::setprecision(10) << omega << std::setw(20) << std::scientific
         << std::setprecision(10) << lambda_sq << "\n";
  }

  file.close();
  std::cout << "Exported non-rigid eigenvalues to: " << filename << std::endl;
}

std::vector<int> FEMHessianAssembler::exportPhononModes(
    const std::string &filename, const EigenResults &eigen_results,
    const ParticipationAnalysis &participation, int num_rigid_body) {
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error: Could not open file " << filename << std::endl;
    return {};
  }

  file << "# Phonon modes (extended modes with P >= " << participation.threshold
       << ")\n";
  file << "# Columns: Mode_Index, Eigenvalue(λ), Frequency(ω), "
          "Lambda_Squared(λ²), P_ratio, L_eff\n";
  file << "#\n";
  file << std::setw(12) << "Mode" << std::setw(20) << "Lambda" << std::setw(20)
       << "Omega" << std::setw(20) << "Lambda_Squared" << std::setw(15)
       << "P_ratio" << std::setw(15) << "L_eff"
       << "\n";

  std::vector<int> phonon_modes;

  for (int i = num_rigid_body; i < eigen_results.num_computed; i++) {
    if (participation.is_extended[i]) {
      double lambda = eigen_results.eigenvalues(i);
      double omega = std::sqrt(std::abs(lambda));
      if (lambda < 0)
        omega = -omega;
      double lambda_sq = lambda * lambda;

      file << std::setw(12) << i << std::setw(20) << std::scientific
           << std::setprecision(10) << lambda << std::setw(20)
           << std::scientific << std::setprecision(10) << omega << std::setw(20)
           << std::scientific << std::setprecision(10) << lambda_sq
           << std::setw(15) << std::fixed << std::setprecision(6)
           << participation.participation_ratios[i] << std::setw(15)
           << std::fixed << std::setprecision(2)
           << participation.localization_lengths[i] << "\n";

      phonon_modes.push_back(i);
    }
  }

  file.close();
  std::cout << "Exported " << phonon_modes.size()
            << " phonon modes to: " << filename << std::endl;

  return phonon_modes;
}

std::vector<int> FEMHessianAssembler::exportLocalizedModes(
    const std::string &filename, const EigenResults &eigen_results,
    const ParticipationAnalysis &participation, int num_rigid_body) {
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error: Could not open file " << filename << std::endl;
    return {};
  }

  file << "# Localized modes (non-phonon modes with P < "
       << participation.threshold << ")\n";
  file << "# Columns: Mode_Index, Eigenvalue(λ), Frequency(ω), "
          "Lambda_Squared(λ²), P_ratio, L_eff\n";
  file << "#\n";
  file << std::setw(12) << "Mode" << std::setw(20) << "Lambda" << std::setw(20)
       << "Omega" << std::setw(20) << "Lambda_Squared" << std::setw(15)
       << "P_ratio" << std::setw(15) << "L_eff"
       << "\n";

  std::vector<int> localized_modes;

  for (int i = num_rigid_body; i < eigen_results.num_computed; i++) {
    if (!participation.is_extended[i]) {
      double lambda = eigen_results.eigenvalues(i);
      double omega = std::sqrt(std::abs(lambda));
      if (lambda < 0)
        omega = -omega;
      double lambda_sq = lambda * lambda;

      file << std::setw(12) << i << std::setw(20) << std::scientific
           << std::setprecision(10) << lambda << std::setw(20)
           << std::scientific << std::setprecision(10) << omega << std::setw(20)
           << std::scientific << std::setprecision(10) << lambda_sq
           << std::setw(15) << std::fixed << std::setprecision(6)
           << participation.participation_ratios[i] << std::setw(15)
           << std::fixed << std::setprecision(2)
           << participation.localization_lengths[i] << "\n";

      localized_modes.push_back(i);
    }
  }

  file.close();
  std::cout << "Exported " << localized_modes.size()
            << " localized modes to: " << filename << std::endl;

  return localized_modes;
}

int FEMHessianAssembler::detectRigidBodyModes(const EigenResults &eigen_results,
                                              double threshold) {
  int num_rigid_body = 0;
  for (int i = 0; i < eigen_results.num_computed; i++) {
    if (std::abs(eigen_results.eigenvalues(i)) < threshold) {
      num_rigid_body++;
    } else {
      break; // Assumes sorted
    }
  }
  return num_rigid_body;
}

void FEMHessianAssembler::exportCompleteEigenmodeAnalysis(
    const std::string &base_filename, const EigenResults &eigen_results,
    const std::vector<Point2D> &current_points,
    const std::vector<std::pair<int, int>> &dof_mapping,
    const std::vector<ElementTriangle2D> &elements, int num_rigid_body,
    double participation_threshold, double vtk_scale_factor) {
  std::cout << "\n" << std::string(70, '=') << std::endl;
  std::cout << "COMPLETE EIGENMODE ANALYSIS AND EXPORT" << std::endl;
  std::cout << std::string(70, '=') << std::endl;

  // Export all non-rigid eigenvalues
  std::string eigen_filename = base_filename + "_eigenvalues_nonrigid.dat";
  exportNonRigidEigenvalues(eigen_filename, eigen_results, num_rigid_body);

  // Perform participation ratio analysis
  std::cout << "\nPerforming participation ratio analysis..." << std::endl;
  ParticipationAnalysis participation = analyzeParticipationRatios(
      eigen_results, dof_mapping, current_points.size(),
      participation_threshold);

  // Print summary
  printParticipationSummary(participation, eigen_results);

  // Export participation analysis
  std::string participation_filename =
      base_filename + "_participation_analysis.dat";
  exportParticipationAnalysis(participation_filename, participation,
                              eigen_results);

  // Export phonon modes
  std::string phonon_filename = base_filename + "_phonon_modes.dat";
  std::vector<int> phonon_modes = exportPhononModes(
      phonon_filename, eigen_results, participation, num_rigid_body);

  // Export localized modes
  std::string localized_filename = base_filename + "_localized_modes.dat";
  std::vector<int> localized_modes = exportLocalizedModes(
      localized_filename, eigen_results, participation, num_rigid_body);

  // NEW: Export density of states
  exportAllDensityOfStates(
      base_filename, eigen_results, participation, num_rigid_body,
      50,  // Number of bins
      true // Gaussian broadening (set to true for smoother curves)
  );

  // Export VTK files for phonon modes
  if (!phonon_modes.empty()) {
    std::cout << "\nExporting phonon mode VTK files..." << std::endl;
    exportEigenvectorsToVTK(base_filename + "_phonon", current_points,
                            eigen_results, phonon_modes, vtk_scale_factor,
                            dof_mapping, elements);
  }

  // Export VTK files for localized modes
  if (!localized_modes.empty()) {
    std::cout << "\nExporting localized mode VTK files..." << std::endl;
    exportEigenvectorsToVTK(base_filename + "_localized", current_points,
                            eigen_results, localized_modes, vtk_scale_factor,
                            dof_mapping, elements);
  }

  // Summary statistics
  std::cout << "\n" << std::string(70, '=') << std::endl;
  std::cout << "EXPORT SUMMARY" << std::endl;
  std::cout << std::string(70, '=') << std::endl;
  std::cout << "Total modes computed: " << eigen_results.num_computed
            << std::endl;
  std::cout << "Rigid body modes: " << num_rigid_body << std::endl;
  std::cout << "Non-rigid modes: "
            << (eigen_results.num_computed - num_rigid_body) << std::endl;
  std::cout << "  - Phonon modes (extended): " << phonon_modes.size()
            << std::endl;
  std::cout << "  - Localized modes: " << localized_modes.size() << std::endl;
  std::cout << "  - Fraction phonon: "
            << static_cast<double>(phonon_modes.size()) /
                   (eigen_results.num_computed - num_rigid_body)
            << std::endl;
  std::cout << "\nFiles created:" << std::endl;
  std::cout << "  " << eigen_filename << std::endl;
  std::cout << "  " << participation_filename << std::endl;
  std::cout << "  " << phonon_filename << std::endl;
  std::cout << "  " << localized_filename << std::endl;
  std::cout << "  " << base_filename << "_dos_all_modes.dat" << std::endl;
  std::cout << "  " << base_filename << "_dos_phonon_modes.dat" << std::endl;
  std::cout << "  " << base_filename << "_dos_localized_modes.dat" << std::endl;
  std::cout << "  " << base_filename << "_phonon_mode_*.vtk ("
            << phonon_modes.size() << " files)" << std::endl;
  std::cout << "  " << base_filename << "_localized_mode_*.vtk ("
            << localized_modes.size() << " files)" << std::endl;
  std::cout << std::string(70, '=') << std::endl;
}

DOSResults FEMHessianAssembler::computeDensityOfStates(
    const EigenResults &eigen_results, const std::vector<int> &mode_indices,
    int num_bins, bool use_gaussian_broadening, double sigma) {
  DOSResults dos_results;

  if (mode_indices.empty()) {
    std::cerr << "Warning: No modes to compute DOS" << std::endl;
    return dos_results;
  }

  // Extract frequencies (ω = √λ) for selected modes
  std::vector<double> frequencies;
  frequencies.reserve(mode_indices.size());

  for (int idx : mode_indices) {
    if (idx >= 0 && idx < eigen_results.num_computed) {
      double lambda = eigen_results.eigenvalues(idx);
      double omega = std::sqrt(std::abs(lambda));
      if (lambda < 0)
        omega = -omega; // Keep sign for unstable modes
      frequencies.push_back(omega);
    }
  }

  if (frequencies.empty()) {
    std::cerr << "Warning: No valid frequencies for DOS" << std::endl;
    return dos_results;
  }

  // Find frequency range
  dos_results.omega_min =
      *std::min_element(frequencies.begin(), frequencies.end());
  dos_results.omega_max =
      *std::max_element(frequencies.begin(), frequencies.end());

  // Add small padding to range
  double range = dos_results.omega_max - dos_results.omega_min;
  // dos_results.omega_min -= 0.05 * range;
  // dos_results.omega_max += 0.05 * range;

  dos_results.num_bins = num_bins;
  double bin_width = (dos_results.omega_max - dos_results.omega_min) / num_bins;

  // Initialize bins
  dos_results.omega_bins.resize(num_bins);
  dos_results.dos.resize(num_bins, 0.0);

  for (int i = 0; i < num_bins; i++) {
    dos_results.omega_bins[i] = dos_results.omega_min + (i + 0.5) * bin_width;
  }

  if (use_gaussian_broadening) {
    // Gaussian broadening for smoother DOS
    if (sigma < 0) {
      sigma = bin_width * 2.0; // Auto sigma = 2 * bin width
    }

    std::cout << "  Using Gaussian broadening with σ = " << sigma << std::endl;

    for (double freq : frequencies) {
      for (int i = 0; i < num_bins; i++) {
        double omega = dos_results.omega_bins[i];
        double gauss = std::exp(-0.5 * std::pow((omega - freq) / sigma, 2)) /
                       (sigma * std::sqrt(2.0 * M_PI));
        dos_results.dos[i] += gauss;
      }
    }
  } else {
    // Simple histogram binning
    for (double freq : frequencies) {
      int bin_idx =
          static_cast<int>((freq - dos_results.omega_min) / bin_width);
      if (bin_idx >= 0 && bin_idx < num_bins) {
        dos_results.dos[bin_idx] += 1.0;
      }
    }

    // Normalize by bin width to get density
    for (int i = 0; i < num_bins; i++) {
      dos_results.dos[i] /= bin_width;
    }
  }

  // Compute cumulative DOS
  dos_results.dos_cumulative.resize(num_bins);
  dos_results.dos_cumulative[0] = dos_results.dos[0] * bin_width;
  for (int i = 1; i < num_bins; i++) {
    dos_results.dos_cumulative[i] =
        dos_results.dos_cumulative[i - 1] + dos_results.dos[i] * bin_width;
  }

  // Compute g(ω)/ω
  dos_results.dos_over_omega.resize(num_bins);
  for (int i = 0; i < num_bins; i++) {
    double omega = dos_results.omega_bins[i];

    if (std::abs(omega) > 1e-8) {
      dos_results.dos_over_omega[i] = dos_results.dos[i] / omega;
    } else {
      // Fit g(ω) = a*ω + b*ω² near zero
      // Then g(ω)/ω = a + b*ω
      if (i + 2 < num_bins) {
        double w1 = dos_results.omega_bins[i + 1];
        double w2 = dos_results.omega_bins[i + 2];
        double g1 = dos_results.dos[i + 1];
        double g2 = dos_results.dos[i + 2];

        // Linear fit: a ≈ g1/w1
        double a = g1 / w1;
        dos_results.dos_over_omega[i] = a; // Limit as ω→0
      } else {
        dos_results.dos_over_omega[i] = 0.0;
      }
    }
  }

  return dos_results;
}

void FEMHessianAssembler::exportDensityOfStates(const std::string &filename,
                                                const DOSResults &dos_results,
                                                const std::string &label) {
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error: Could not open file " << filename << std::endl;
    return;
  }

  file << "# Density of States: " << label << "\n";
  file << "# Number of bins: " << dos_results.num_bins << "\n";
  file << "# Frequency range: [" << dos_results.omega_min << ", "
       << dos_results.omega_max << "]\n";
  file << "# Columns: Omega(ω), g(ω), g(ω)/ω, Cumulative_DOS, log10 Omega(ω), "
          "log10 g(ω), log10 Cumulative_DOS, log10 g(ω)/ω \n";
  file << "#\n";
  file << std::setw(20) << "Omega" << std::setw(20) << "g(omega)"
       << std::setw(20) << "g/omega" << std::setw(20) << "Cumulative"
       << std::setw(20) << "log10_Omega" << std::setw(20) << "log10_g(omega)"
       << std::setw(20) << "log10_Cumulative" << std::setw(20)
       << "log10_g/omega"
       << "\n";

  for (int i = 0; i < dos_results.num_bins; i++) {
    double omega = dos_results.omega_bins[i];
    double g = dos_results.dos[i];
    double g_over_omega = dos_results.dos_over_omega[i];
    double cumulative = dos_results.dos_cumulative[i];

    // Safe log10 (avoid log of zero or negative)
    auto safe_log10 = [](double x) {
      return (x > 1e-100) ? std::log10(x) : -100.0;
    };

    file << std::setw(20) << std::scientific << std::setprecision(10) << omega
         << std::setw(20) << std::scientific << std::setprecision(10) << g
         << std::setw(20) << std::scientific << std::setprecision(10)
         << g_over_omega << std::setw(20) << std::fixed << std::setprecision(10)
         << cumulative << std::setw(20) << std::fixed << std::setprecision(10)
         << safe_log10(std::abs(omega)) << std::setw(20) << std::fixed
         << std::setprecision(10) << safe_log10(g) << std::setw(20)
         << std::fixed << std::setprecision(10) << safe_log10(cumulative)
         << std::setw(20) << std::fixed << std::setprecision(10)
         << safe_log10(g_over_omega) << "\n";
  }

  file.close();
  std::cout << "Exported DOS to: " << filename << std::endl;
}

void FEMHessianAssembler::exportAllDensityOfStates(
    const std::string &base_filename, const EigenResults &eigen_results,
    const ParticipationAnalysis &participation, int num_rigid_body,
    int num_bins, bool use_gaussian_broadening) {
  std::cout << "\n" << std::string(70, '=') << std::endl;
  std::cout << "COMPUTING DENSITY OF STATES" << std::endl;
  std::cout << std::string(70, '=') << std::endl;

  // Prepare mode index lists
  std::vector<int> all_modes;
  std::vector<int> phonon_modes;
  std::vector<int> localized_modes;

  for (int i = num_rigid_body; i < eigen_results.num_computed; i++) {
    all_modes.push_back(i);
    if (participation.is_extended[i]) {
      phonon_modes.push_back(i);
    } else {
      localized_modes.push_back(i);
    }
  }

  std::cout << "Total non-rigid modes: " << all_modes.size() << std::endl;
  std::cout << "Phonon modes: " << phonon_modes.size() << std::endl;
  std::cout << "Localized modes: " << localized_modes.size() << std::endl;
  std::cout << "Number of bins: " << num_bins << std::endl;
  std::cout << "Gaussian broadening: "
            << (use_gaussian_broadening ? "Yes" : "No") << std::endl;

  // Compute DOS for all modes
  std::cout << "\nComputing DOS for all modes..." << std::endl;
  DOSResults dos_all = computeDensityOfStates(
      eigen_results, all_modes, num_bins, use_gaussian_broadening);
  exportDensityOfStates(base_filename + "_dos_all_modes.dat", dos_all,
                        "All non-rigid modes");

  // Compute DOS for phonon modes
  if (!phonon_modes.empty()) {
    std::cout << "\nComputing DOS for phonon modes..." << std::endl;
    DOSResults dos_phonon = computeDensityOfStates(
        eigen_results, phonon_modes, num_bins, use_gaussian_broadening);
    exportDensityOfStates(base_filename + "_dos_phonon_modes.dat", dos_phonon,
                          "Phonon (extended) modes");
  }

  // Compute DOS for localized modes
  if (!localized_modes.empty()) {
    std::cout << "\nComputing DOS for localized modes..." << std::endl;
    DOSResults dos_localized = computeDensityOfStates(
        eigen_results, localized_modes, num_bins, use_gaussian_broadening);
    exportDensityOfStates(base_filename + "_dos_localized_modes.dat",
                          dos_localized, "Localized modes");
  }

  std::cout << std::string(70, '=') << std::endl;
}

// Helper function for eigenvalue analysis - returns only min non-rigid
// eigenvalue
double FEMHessianAssembler::computeMinNonRigidEigenvalue(
    std::vector<ElementTriangle2D> &elements,
    const std::vector<Point2D> &square_points, int n_dofs,
    const std::vector<std::pair<int, int>> &full_mapping,
    Strain_Energy_LatticeCalculator *calculator,
    std::function<double(double)> potential_func_der,
    std::function<double(double)> potential_func_sder, int num_eigenvalues,
    double alpha, const std::string &output_filename) {
  // 1. Create FEM Hessian assembler
  FEMHessianAssembler assembler;

  // 2. Set energy parameters
  assembler.setEnergyParameters(calculator, potential_func_der,
                                potential_func_sder,
                                calculator->getUnitCellArea());

  // 3. Assemble global stiffness matrix
  Eigen::SparseMatrix<double> global_stiffness =
      assembler.assembleGlobalStiffness(elements, square_points, n_dofs,
                                        full_mapping);

  std::cout << std::scientific << std::setprecision(17);

  // std::cout << "--- Sparse Matrix Output (Row, Col, Value) ---" << std::endl;

  // // Iterate over the sparse matrix efficiently
  // // outerSize() is usually the number of columns (for column-major matrices)
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

  // 4. Compute eigenvalues
  EigenResults results = assembler.computeSmallestEigenvaluesIterative_spectra(
      global_stiffness, num_eigenvalues, 0);

  // 5. Sort by ABSOLUTE VALUE to identify rigid body modes
  std::vector<std::pair<double, int>> eigen_pairs;
  for (int i = 0; i < results.num_computed; i++) {
    eigen_pairs.push_back({results.eigenvalues(i), i});
  }
  std::sort(eigen_pairs.begin(), eigen_pairs.end(),
            [](const auto &a, const auto &b) {
              return std::abs(a.first) < std::abs(b.first);
            });

  // 6. Build temporarily sorted arrays for rigid body mode detection
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

  // 7. Detect rigid body modes
  int num_rigid = assembler.detectRigidBodyModes(results);

  // 8. Create NEW pairs for non-rigid modes and sort them
  std::vector<std::pair<double, int>> non_rigid_pairs;
  for (int i = num_rigid; i < results.num_computed; i++) {
    non_rigid_pairs.push_back({results.eigenvalues(i), i});
  }

  // Sort non-rigid modes from negative to positive
  std::sort(non_rigid_pairs.begin(), non_rigid_pairs.end(),
            [](const auto &a, const auto &b) { return a.first < b.first; });

  // 9. Build FINAL sorted arrays
  for (int i = 0; i < num_rigid; i++) {
    sorted_eigenvalues(i) = results.eigenvalues(i);
    sorted_eigenvectors.col(i) = results.eigenvectors.col(i);
  }

  for (int i = 0; i < non_rigid_pairs.size(); i++) {
    sorted_eigenvalues(num_rigid + i) =
        results.eigenvalues(non_rigid_pairs[i].second);
    sorted_eigenvectors.col(num_rigid + i) =
        results.eigenvectors.col(non_rigid_pairs[i].second);
  }

  // 10. Extract minimum, second, and third non-rigid eigenvalues
  double first_non_rigid_eigenvalue =
      (num_rigid < results.num_computed)
          ? sorted_eigenvalues(num_rigid)
          : std::numeric_limits<double>::quiet_NaN();

  double second_non_rigid_eigenvalue =
      (num_rigid + 1 < results.num_computed)
          ? sorted_eigenvalues(num_rigid + 1)
          : std::numeric_limits<double>::quiet_NaN();

  double third_non_rigid_eigenvalue =
      (num_rigid + 2 < results.num_computed)
          ? sorted_eigenvalues(num_rigid + 2)
          : std::numeric_limits<double>::quiet_NaN();

  std::cout << "Num rigid modes: " << num_rigid << std::endl;
  std::cout << "  1st non-rigid eigenvalue: " << first_non_rigid_eigenvalue
            << std::endl;
  std::cout << "  2nd non-rigid eigenvalue: " << second_non_rigid_eigenvalue
            << std::endl;
  std::cout << "  3rd non-rigid eigenvalue: " << third_non_rigid_eigenvalue
            << std::endl;

  // 11. Write to file
  std::ofstream outfile;

  // Check if file exists to decide whether to write header
  bool file_exists = std::ifstream(output_filename).good();

  outfile.open(output_filename, std::ios::app); // Append mode
  if (!outfile.is_open()) {
    std::cerr << "Error: Could not open file " << output_filename
              << " for writing" << std::endl;
    return first_non_rigid_eigenvalue;
  }

  // Write header if new file
  if (!file_exists) {
    outfile << "# Alpha, Eigenvalue_1st, Eigenvalue_2nd, Eigenvalue_3rd, "
               "Num_Rigid_Modes\n";
  }

  // Write data
  outfile << std::scientific << std::setprecision(12);
  outfile << alpha << " " << first_non_rigid_eigenvalue << " "
          << second_non_rigid_eigenvalue << " " << third_non_rigid_eigenvalue
          << " " << num_rigid << "\n";
  outfile.close();

  std::cout << "Saved to " << output_filename << ": alpha=" << alpha
            << ", eig1=" << first_non_rigid_eigenvalue
            << ", eig2=" << second_non_rigid_eigenvalue
            << ", eig3=" << third_non_rigid_eigenvalue << std::endl;

  return first_non_rigid_eigenvalue;
}



std::vector<int> FEMHessianAssembler::getLocalizedModes(
    const EigenResults& eigen_results,
    const std::vector<std::pair<int, int>>& dof_mapping,
    int num_nodes,
    int num_rigid,
    double threshold)
{
    // Compute participation ratios
    ParticipationAnalysis participation = analyzeParticipationRatios(
        eigen_results, dof_mapping, num_nodes, threshold);
    
    // Extract localized mode indices
    std::vector<int> localized_indices;
    for (int i = num_rigid; i < eigen_results.num_computed; i++) {
        if (!participation.is_extended[i]) {
            localized_indices.push_back(i);
        }
    }
    
    return localized_indices;
}

Eigen::VectorXd FEMHessianAssembler::computeDiagonalPreconditioner(
    std::vector<ElementTriangle2D> &elements,
    const std::vector<Point2D> &points,
    const std::vector<std::pair<int, int>> &dof_mapping,
    int num_total_dofs) {
  
  std::cout << "Computing diagonal preconditioner (fast version)..." << std::endl;
  
  // Assemble ONLY diagonal of Hessian - much faster!
  Eigen::VectorXd H_diag = assembleGlobalStiffnessDiagonal(
      elements, points, num_total_dofs, dof_mapping);
  
  // Compute preconditioner from diagonal
  Eigen::VectorXd diag_precond(num_total_dofs);
  
  for (int i = 0; i < num_total_dofs; i++) {
    double diag_val = H_diag(i);
    
    // Prevent division by zero and ensure positive
    if (std::abs(diag_val) < 1e-12) {
      diag_precond(i) = 1.0;
    } else {
      // Use inverse of square root of diagonal for preconditioning
      diag_precond(i) = 1.0 / std::sqrt(std::abs(diag_val));
    }
  }
  
  std::cout << "Preconditioner computed. Min: " << diag_precond.minCoeff() 
            << ", Max: " << diag_precond.maxCoeff() << std::endl;
  
  return diag_precond;
}

alglib::real_1d_array FEMHessianAssembler::eigenToAlglibArray(
    const Eigen::VectorXd& eigen_vec) {
  
  alglib::real_1d_array alglib_array;
  alglib_array.setlength(eigen_vec.size());
  
  for (int i = 0; i < eigen_vec.size(); i++) {
    alglib_array[i] = eigen_vec(i);
  }
  
  return alglib_array;
}



// Add to FEMHessianAssembler.cpp

Eigen::VectorXd FEMHessianAssembler::assembleGlobalStiffnessDiagonal(
    std::vector<ElementTriangle2D> &elements,
    const std::vector<Point2D> &current_points,
    int num_total_dofs,
    const std::vector<std::pair<int, int>> &dof_mapping) {
  
  // Safety check: ensure energy parameters were set
  if (strain_calculator == nullptr) {
    std::cerr << "Error: Energy parameters not set! Call setEnergyParameters() first." 
              << std::endl;
    throw std::runtime_error("FEMHessianAssembler: energy parameters not initialized");
  }

  // Initialize diagonal vector
  Eigen::VectorXd diagonal = Eigen::VectorXd::Zero(num_total_dofs);
  
  int num_free_nodes = num_total_dofs / 2;

  // Loop over all elements
  for (size_t elem_idx = 0; elem_idx < elements.size(); elem_idx++) {
    auto &element = elements[elem_idx];

    // Calculate deformation gradient for this element
    element.calculate_deformation_gradient(current_points);
    Eigen::Matrix2d F = element.getDeformationGradient();

    // Calculate right Cauchy-Green tensor
    Eigen::Matrix2d C = F.transpose() * F;

    // Perform Lagrange reduction
    lagrange::Result reduction_result = lagrange::reduce(C);
    Eigen::Matrix2d Z = reduction_result.m_matrix;

    // Create acoustic tensor for this element
    AcousticTensor acoustic_tensor(F, C, Z);

    // Get element area for integration weight
    double area = element.getReferenceArea();

    // Compute FULL element stiffness matrix (still need it for diagonal)
    Eigen::MatrixXd K_elem = computeElementStiffness(acoustic_tensor, element, area);

    // Get global DOF indices for this element
    std::vector<int> global_dof_indices;
    global_dof_indices.reserve(6);

    for (int local_node = 0; local_node < 3; local_node++) {
      int global_node = element.getNodeIndex(local_node);

      // Get DOF mapping for this node
      auto [orig_idx, solver_idx] = dof_mapping[global_node];

      if (solver_idx != -1) { // Free DOF
        global_dof_indices.push_back(solver_idx);                  // x DOF
        global_dof_indices.push_back(solver_idx + num_free_nodes); // y DOF
      } else {
        global_dof_indices.push_back(-1);
        global_dof_indices.push_back(-1);
      }
    }

    // Accumulate ONLY diagonal entries
    for (int i = 0; i < 6; i++) {
      int global_i = global_dof_indices[i];
      if (global_i < 0) continue; // Skip fixed DOFs

      // Add diagonal contribution
      diagonal(global_i) += K_elem(i, i);
    }
  }

  return diagonal;
}