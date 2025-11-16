#include "FEMHessianAssembler.h"
#include "../include/reductions/LagrangeReduction.h"

#include <iostream>
#include <algorithm>
// SPECTRA INCLUDES
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>


    #include <Spectra/SymEigsSolver.h>
    #include <Spectra/MatOp/SparseSymMatProd.h>



void FEMHessianAssembler::extractAcousticTensor(
    const itensor::ITensor& A_tensor,
    double A[2][2][2][2])
{
    // Extract A_{iKjL} from ITensor with indices (r_, i_, s_, j_)
    // where r_=K (material coord 1), s_=L (material coord 2)
    for (int i = 1; i <= 2; i++) {
        for (int K = 1; K <= 2; K++) {
            for (int j = 1; j <= 2; j++) {
                for (int L = 1; L <= 2; L++) {
                    // Access pattern: (r_=K, i_=i, s_=L, j_=j)
                    // ITensor uses 1-based indexing, array uses 0-based
                    A[i-1][K-1][j-1][L-1] = A_tensor.elt(K, i, L, j);
                }
            }
        }
    }
}

Eigen::MatrixXd FEMHessianAssembler::computeElementStiffness(
    AcousticTensor& acoustic_tensor,
    const ElementTriangle2D& element,
    double area_weight)
{
    // Safety check: ensure energy parameters were set
    if (strain_calculator == nullptr) {
        std::cerr << "Error: Energy parameters not set! Call setEnergyParameters() first." << std::endl;
        throw std::runtime_error("FEMHessianAssembler: energy parameters not initialized");
    }
    
    const int num_nodes = 3;      // Triangular element
    const int spatial_dim = 2;     // 2D
    const int total_dofs = num_nodes * spatial_dim;  // 6 DOFs
    
    Eigen::MatrixXd K_elem = Eigen::MatrixXd::Zero(total_dofs, total_dofs);
    
    
    // Compute energy derivatives using member variables
    // CRITICAL: Dereference pointer with * to pass as reference
    acoustic_tensor.computeEnergyDerivatives(
        *strain_calculator,
        potential_func_der,
        potential_func_sder,
        normalisation
    );
    
    // Get shape function derivatives: dN^a/dX_K (3 nodes × 2 coords)
    const Eigen::Matrix<double, 3, 2>& dN_dX = element.getDNdX();
    
    // Get acoustic tensor A_{iKjL} (Lagrangian frame)
    // Returns ITensor with indices (r_, i_, s_, j_) where r_=K, s_=L
    itensor::ITensor A_tensor = acoustic_tensor.getAcousticTensor(true);
    
    // Extract to 4D array for easier manipulation
    double A[2][2][2][2];
    extractAcousticTensor(A_tensor, A);
    
    // Assemble element stiffness: K^{ab}_{ij} = A_{iKjL} * dN^a_K * dN^b_L
    for (int a = 0; a < num_nodes; a++) {           // Node a
        for (int b = 0; b < num_nodes; b++) {       // Node b
            for (int i = 0; i < spatial_dim; i++) { // DOF direction i
                for (int j = 0; j < spatial_dim; j++) { // DOF direction j
                    
                    double K_ab_ij = 0.0;
                    
                    // Contract over material coordinates K, L
                    for (int K = 0; K < spatial_dim; K++) {
                        for (int L = 0; L < spatial_dim; L++) {
                            K_ab_ij += A[i][K][j][L] * dN_dX(a, K) * dN_dX(b, L);
                        }
                    }
                    
                    // Multiply by integration weight (includes det(J) for reference element)
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
    std::vector<ElementTriangle2D>& elements,
    const std::vector<Point2D>& current_points,
    int num_total_dofs,
    const std::vector<std::pair<int, int>>& dof_mapping)
{
    // Safety check: ensure energy parameters were set
    if (strain_calculator == nullptr) {
        std::cerr << "Error: Energy parameters not set! Call setEnergyParameters() first." << std::endl;
        throw std::runtime_error("FEMHessianAssembler: energy parameters not initialized");
    }
    
    // Create triplet list for sparse matrix assembly
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(elements.size() * 36); // 6x6 entries per element
    
    for (size_t elem_idx = 0; elem_idx < elements.size(); elem_idx++) {
        auto& element = elements[elem_idx];
        
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
        double area = element.getArea();
        
        // Compute element stiffness matrix
        // (computeEnergyDerivatives is called inside computeElementStiffness)
        Eigen::MatrixXd K_elem = computeElementStiffness(
            acoustic_tensor, element, area);
        
        // Get global DOF indices for this element
        std::vector<int> global_dof_indices;
        global_dof_indices.reserve(6);
        
        for (int local_node = 0; local_node < 3; local_node++) {
            int global_node = element.getNodeIndex(local_node);
            
            // Get DOF mapping for this node
            auto [orig_idx, solver_idx] = dof_mapping[global_node];
            
            if (solver_idx != -1) {  // Free DOF
                // Add x and y DOF indices
                // Convention: [u0, u1, ..., u_n, v0, v1, ..., v_n]
                int num_free_nodes = num_total_dofs / 2;
                global_dof_indices.push_back(solver_idx);              // x DOF
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
            if (global_i < 0) continue;  // Skip fixed DOFs
            
            for (int j = 0; j < 6; j++) {
                int global_j = global_dof_indices[j];
                if (global_j < 0) continue;  // Skip fixed DOFs
                
                if (std::abs(K_elem(i, j)) > 1e-14) {  // Only add non-zero entries
                    triplets.push_back(
                        Eigen::Triplet<double>(global_i, global_j, K_elem(i, j))
                    );
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
        std::cerr << "Warning: Global stiffness matrix is not symmetric!" << std::endl;
        std::cerr << "Asymmetry norm: " << (K_global - K_transpose).norm() << std::endl;
    }
    
    return K_global;
}

EigenResults FEMHessianAssembler::computeSmallestEigenvalues(
    const Eigen::SparseMatrix<double>& K_global,
    int N)
{
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
            std::cout << "Largest of the " << N << " smallest: " << results.eigenvalues(N-1) << std::endl;
        }
        
    } else {
        std::cerr << "Warning: Matrix too large (" << n << " DOFs) for dense solver." << std::endl;
        std::cerr << "Use computeSmallestEigenvaluesIterative() instead." << std::endl;
    }
    
    return results;
}

EigenResults FEMHessianAssembler::computeSmallestEigenvaluesIterative(
    const Eigen::SparseMatrix<double>& K_global,
    int N,
    double shift)
{
    EigenResults results;
    
    int n = K_global.rows();
    
    // Check if N is valid
    if (N <= 0 || N > n) {
        std::cerr << "Error: N must be between 1 and " << n << std::endl;
        return results;
    }
    
    std::cout << "Computing " << N << " smallest eigenvalues using iterative method..." << std::endl;
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
        std::cerr << "Try a different shift value or use computeSmallestEigenvalues()" << std::endl;
        return results;
    }
    
    std::cout << "Note: Full iterative eigenvalue solver (Arnoldi/Lanczos) requires additional libraries." << std::endl;
    std::cout << "For now, falling back to dense solver for moderate-sized problems." << std::endl;
    
    // For a complete implementation, you would use libraries like ARPACK or Spectra here
    // For now, provide a basic implementation that works for moderate sizes
    
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
        std::cerr << "Consider using external libraries like ARPACK or Spectra." << std::endl;
    }
    
    return results;
}


EigenResults FEMHessianAssembler::computeSmallestEigenvaluesIterative_spectra(
    const Eigen::SparseMatrix<double>& K_global,
    int N,
    double shift)
{
    EigenResults results;
    int n = K_global.rows();

    // Check if N is valid
    if (N <= 0 || N > n) {
        std::cerr << "Error: N must be between 1 and " << n << std::endl;
        return results;
    }

    std::cout << "Computing " << N << " smallest eigenvalues using Spectra..." << std::endl;
    std::cout << "Matrix size: " << n << " x " << n << std::endl;

    // Number of Lanczos vectors (ncv must be > N)
    // Rule of thumb: ncv = 2*N to 3*N, but ncv < n
    int ncv = std::max(2 * N + 1, 20);
    ncv = std::min(ncv, n - 1);
    
    std::cout << "Spectra: N=" << N << ", ncv=" << ncv << std::endl;

    try {
        // Construct matrix operation object using the wrapper class
        Spectra::SparseSymMatProd<double> op(K_global);
        
        // Construct symmetric eigen solver, requesting N smallest eigenvalues
        Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigs(op, N, ncv);
        
        // Initialize and compute
        eigs.init();
        std::cout << "Spectra: Computing eigenvalues..." << std::endl;
        int nconv = eigs.compute(Spectra::SortRule::SmallestAlge, 1000, 1e-10);
        
        std::cout << "Spectra: Converged " << nconv << " eigenvalues" << std::endl;
        
        // Check convergence
        if (eigs.info() != Spectra::CompInfo::Successful) {
            std::cerr << "Error: Spectra failed to converge!" << std::endl;
            std::cerr << "Info code: " << static_cast<int>(eigs.info()) << std::endl;
            return results;
        }
        
        // Retrieve results (already sorted by Spectra)
        results.eigenvalues = eigs.eigenvalues();
        results.eigenvectors = eigs.eigenvectors();
        results.num_computed = nconv;
        
        // Print results
        std::cout << "Successfully computed " << results.num_computed << " eigenvalues." << std::endl;
        if (results.num_computed > 0) {
            std::cout << "Smallest eigenvalue: " << results.eigenvalues(0) << std::endl;
            if (results.num_computed > 1) {
                std::cout << "Largest of computed eigenvalues: " 
                         << results.eigenvalues(results.num_computed - 1) << std::endl;
            }
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Exception in Spectra computation: " << e.what() << std::endl;
    }

    return results;
}


#include <fstream>
#include <iomanip>

#include <fstream>
#include <iomanip>
#include <sstream>

void FEMHessianAssembler::exportSingleModeToVTK(
    const std::string& filename,
    const std::vector<Point2D>& current_points,
    const Eigen::VectorXd& eigenvector,
    double eigenvalue,
    int mode_number,
    double scale_factor,
    const std::vector<std::pair<int, int>>& dof_mapping,
    const std::vector<ElementTriangle2D>& elements)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }
    
    int num_nodes = current_points.size();
    int num_free_nodes = eigenvector.size() / 2;
    
    // Reconstruct full displacement field from eigenvector
    // eigenvector format: [u0, u1, ..., u_n, v0, v1, ..., v_n]
    std::vector<Eigen::Vector2d> displacements(num_nodes, Eigen::Vector2d::Zero());
    
    for (int node = 0; node < num_nodes; node++) {
        auto [orig_idx, solver_idx] = dof_mapping[node];
        
        if (solver_idx != -1) {  // Free DOF
            double u = eigenvector(solver_idx);              // x displacement
            double v = eigenvector(solver_idx + num_free_nodes); // y displacement
            displacements[node] = Eigen::Vector2d(u, v);
        }
        // Fixed DOFs remain zero
    }
    
    // Normalize displacement field for better visualization
    double max_displacement = 0.0;
    for (const auto& disp : displacements) {
        max_displacement = std::max(max_displacement, disp.norm());
    }
    
    if (max_displacement > 1e-12) {
        for (auto& disp : displacements) {
            disp /= max_displacement;
            disp *= scale_factor;
        }
    }
    
    // Filter elements: only keep those with ZERO translation (fundamental domain only)
    std::vector<int> fundamental_domain_elements;
    for (size_t i = 0; i < elements.size(); i++) {
        const auto& elem = elements[i];
        if (!elem.isInitialized()) continue;
        
        // Check if all nodes have zero translation
        bool is_fundamental = true;
        for (int j = 0; j < 3; j++) {
            Eigen::Vector2d translation = elem.getTranslation(j);
            if (translation.norm() > 1e-10) {  // Non-zero translation
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
    file << "Eigenmode " << mode_number << " (lambda = " << std::scientific << eigenvalue << ")\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";
    
    // Write points (current configuration)
    file << "POINTS " << num_nodes << " double\n";
    for (int i = 0; i < num_nodes; i++) {
        file << std::setprecision(10) 
             << current_points[i].coord.x() << " " 
             << current_points[i].coord.y() << " " 
             << 0.0 << "\n";
    }
    
    // Write cells (only fundamental domain - no periodic images)
    file << "\nCELLS " << fundamental_domain_elements.size() 
         << " " << (fundamental_domain_elements.size() * 4) << "\n";
    for (int elem_idx : fundamental_domain_elements) {
        const auto& elem = elements[elem_idx];
        file << "3 " 
             << elem.getNodeIndex(0) << " "
             << elem.getNodeIndex(1) << " "
             << elem.getNodeIndex(2) << "\n";
    }
    
    file << "\nCELL_TYPES " << fundamental_domain_elements.size() << "\n";
    for (size_t i = 0; i < fundamental_domain_elements.size(); i++) {
        file << "5\n";  // VTK_TRIANGLE
    }
    
    // Write point data
    file << "\nPOINT_DATA " << num_nodes << "\n";
    
    // Write eigenvector as displacement vectors (for Glyph/Arrow visualization)
    file << "VECTORS eigenvector double\n";
    for (int i = 0; i < num_nodes; i++) {
        file << std::setprecision(10)
             << displacements[i](0) << " " 
             << displacements[i](1) << " " 
             << 0.0 << "\n";
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
             << (current_points[i].coord.y() + displacements[i](1)) << " "
             << 0.0 << "\n";
    }
    
    file.close();
    std::cout << "  Exported mode " << mode_number << " to " << filename << std::endl;
}

void FEMHessianAssembler::exportEigenvectorsToVTK(
    const std::string& filename,
    const std::vector<Point2D>& current_points,
    const EigenResults& eigen_results,
    const std::vector<int>& mode_indices,
    double scale_factor,
    const std::vector<std::pair<int, int>>& dof_mapping,
    const std::vector<ElementTriangle2D>& elements)
{
    if (eigen_results.num_computed == 0) {
        std::cerr << "Error: No eigenvalues computed!" << std::endl;
        return;
    }
    
    if (dof_mapping.empty()) {
        std::cerr << "Error: DOF mapping required for VTK export!" << std::endl;
        return;
    }
    
    std::cout << "\n=== Exporting Eigenmodes to VTK ===" << std::endl;
    std::cout << "Base filename: " << filename << std::endl;
    std::cout << "Scale factor: " << scale_factor << std::endl;
    std::cout << "Number of modes to export: " << mode_indices.size() << std::endl;
    
    // Count fundamental domain elements (for reporting)
    int fundamental_count = 0;
    for (size_t i = 0; i < elements.size(); i++) {
        const auto& elem = elements[i];
        if (!elem.isInitialized()) continue;
        
        bool is_fundamental = true;
        for (int j = 0; j < 3; j++) {
            Eigen::Vector2d translation = elem.getTranslation(j);
            if (translation.norm() > 1e-10) {
                is_fundamental = false;
                break;
            }
        }
        if (is_fundamental) fundamental_count++;
    }
    
    std::cout << "Total elements: " << elements.size() 
              << ", Fundamental domain elements: " << fundamental_count << std::endl;
    
    for (int mode_idx : mode_indices) {
        if (mode_idx < 0 || mode_idx >= eigen_results.num_computed) {
            std::cerr << "Warning: Mode index " << mode_idx 
                     << " out of range [0, " << eigen_results.num_computed << ")" << std::endl;
            continue;
        }
        
        // Create filename for this mode
        std::stringstream ss;
        ss << filename << "_mode_" << std::setfill('0') << std::setw(3) << mode_idx << ".vtk";
        std::string mode_filename = ss.str();
        
        // Export this mode (automatically handles PBC filtering)
        exportSingleModeToVTK(
            mode_filename,
            current_points,
            eigen_results.eigenvectors.col(mode_idx),
            eigen_results.eigenvalues(mode_idx),
            mode_idx,
            scale_factor,
            dof_mapping,
            elements
        );
    }
    
    std::cout << "\nSuccessfully exported " << mode_indices.size() << " eigenmodes" << std::endl;
    std::cout << "All modes show ONLY fundamental domain (no periodic images)" << std::endl;
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
    const Eigen::VectorXd& eigenvector,
    const std::vector<std::pair<int, int>>& dof_mapping,
    int num_nodes)
{
    if (eigenvector.size() == 0) {
        std::cerr << "Error: Empty eigenvector!" << std::endl;
        return 0.0;
    }
    
    int num_free_nodes = eigenvector.size() / 2;
    
    // Reconstruct displacement magnitudes at each node
    std::vector<double> displacement_magnitude(num_nodes, 0.0);
    
    for (int node = 0; node < num_nodes; node++) {
        auto [orig_idx, solver_idx] = dof_mapping[node];
        
        if (solver_idx != -1) {  // Free DOF
            double u = eigenvector(solver_idx);                    // x displacement
            double v = eigenvector(solver_idx + num_free_nodes);   // y displacement
            displacement_magnitude[node] = std::sqrt(u*u + v*v);
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
    const EigenResults& eigen_results,
    const std::vector<std::pair<int, int>>& dof_mapping,
    int num_nodes,
    double threshold)
{
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
        double P = computeParticipationRatio(
            eigen_results.eigenvectors.col(mode),
            dof_mapping,
            num_nodes
        );
        
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
    const std::string& filename,
    const ParticipationAnalysis& analysis,
    const EigenResults& eigen_results)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }
    
    // Write header
    file << "# Participation Ratio Analysis\n";
    file << "# Threshold for extended modes: " << analysis.threshold << "\n";
    file << "# Columns: Mode_Index, Eigenvalue, Participation_Ratio, Localization_Length, Is_Extended\n";
    file << "# P ≈ 1: Extended mode (phonon), P ≈ 0: Localized mode\n";
    file << "#\n";
    file << std::setw(12) << "Mode" 
         << std::setw(20) << "Eigenvalue"
         << std::setw(20) << "P_ratio"
         << std::setw(20) << "L_eff"
         << std::setw(15) << "Type"
         << "\n";
    
    // Write data
    for (size_t i = 0; i < analysis.participation_ratios.size(); i++) {
        file << std::setw(12) << i
             << std::setw(20) << std::setprecision(10) << eigen_results.eigenvalues(i)
             << std::setw(20) << std::setprecision(6) << analysis.participation_ratios[i]
             << std::setw(20) << std::setprecision(2) << analysis.localization_lengths[i]
             << std::setw(15) << (analysis.is_extended[i] ? "Extended" : "Localized")
             << "\n";
    }
    
    file.close();
    std::cout << "Participation analysis exported to: " << filename << std::endl;
}

void FEMHessianAssembler::printParticipationSummary(
    const ParticipationAnalysis& analysis,
    const EigenResults& eigen_results)
{
    std::cout << "\n=== Participation Ratio Summary ===" << std::endl;
    std::cout << std::setw(8) << "Mode" 
              << std::setw(16) << "Eigenvalue"
              << std::setw(12) << "P_ratio"
              << std::setw(12) << "L_eff"
              << std::setw(12) << "Type"
              << std::endl;
    std::cout << std::string(60, '-') << std::endl;
    
    for (size_t i = 0; i < analysis.participation_ratios.size(); i++) {
        std::cout << std::setw(8) << i
                  << std::setw(16) << std::scientific << std::setprecision(4) 
                  << eigen_results.eigenvalues(i)
                  << std::setw(12) << std::fixed << std::setprecision(4) 
                  << analysis.participation_ratios[i]
                  << std::setw(12) << std::fixed << std::setprecision(1) 
                  << analysis.localization_lengths[i]
                  << std::setw(12) << (analysis.is_extended[i] ? "Extended" : "Localized")
                  << std::endl;
    }
    
    // Statistics
    std::cout << std::string(60, '-') << std::endl;
    
    int num_extended = std::count(analysis.is_extended.begin(), 
                                   analysis.is_extended.end(), true);
    int num_localized = analysis.is_extended.size() - num_extended;
    
    std::cout << "\nStatistics:" << std::endl;
    std::cout << "  Total modes analyzed: " << analysis.participation_ratios.size() << std::endl;
    std::cout << "  Extended modes (phonons): " << num_extended << std::endl;
    std::cout << "  Localized modes: " << num_localized << std::endl;
    std::cout << "  Fraction extended: " 
              << static_cast<double>(num_extended) / analysis.participation_ratios.size() 
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
    const std::string& filename,
    const EigenResults& eigen_results,
    int num_rigid_body)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }
    
    file << "# Non-rigid eigenvalues analysis\n";
    file << "# Number of rigid body modes skipped: " << num_rigid_body << "\n";
    file << "# Columns: Mode_Index, Eigenvalue(λ), Frequency(ω=√λ), Lambda_Squared(λ²)\n";
    file << "#\n";
    file << std::setw(12) << "Mode"
         << std::setw(20) << "Lambda"
         << std::setw(20) << "Omega"
         << std::setw(20) << "Lambda_Squared"
         << "\n";
    
    for (int i = num_rigid_body; i < eigen_results.num_computed; i++) {
        double lambda = eigen_results.eigenvalues(i);
        double omega = std::sqrt(std::abs(lambda));
        if (lambda < 0) omega = -omega;  // Keep sign for negative eigenvalues
        double lambda_sq = lambda * lambda;
        
        file << std::setw(12) << i
             << std::setw(20) << std::scientific << std::setprecision(10) << lambda
             << std::setw(20) << std::scientific << std::setprecision(10) << omega
             << std::setw(20) << std::scientific << std::setprecision(10) << lambda_sq
             << "\n";
    }
    
    file.close();
    std::cout << "Exported non-rigid eigenvalues to: " << filename << std::endl;
}

std::vector<int> FEMHessianAssembler::exportPhononModes(
    const std::string& filename,
    const EigenResults& eigen_results,
    const ParticipationAnalysis& participation,
    int num_rigid_body)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return {};
    }
    
    file << "# Phonon modes (extended modes with P >= " << participation.threshold << ")\n";
    file << "# Columns: Mode_Index, Eigenvalue(λ), Frequency(ω), Lambda_Squared(λ²), P_ratio, L_eff\n";
    file << "#\n";
    file << std::setw(12) << "Mode"
         << std::setw(20) << "Lambda"
         << std::setw(20) << "Omega"
         << std::setw(20) << "Lambda_Squared"
         << std::setw(15) << "P_ratio"
         << std::setw(15) << "L_eff"
         << "\n";
    
    std::vector<int> phonon_modes;
    
    for (int i = num_rigid_body; i < eigen_results.num_computed; i++) {
        if (participation.is_extended[i]) {
            double lambda = eigen_results.eigenvalues(i);
            double omega = std::sqrt(std::abs(lambda));
            if (lambda < 0) omega = -omega;
            double lambda_sq = lambda * lambda;
            
            file << std::setw(12) << i
                 << std::setw(20) << std::scientific << std::setprecision(10) << lambda
                 << std::setw(20) << std::scientific << std::setprecision(10) << omega
                 << std::setw(20) << std::scientific << std::setprecision(10) << lambda_sq
                 << std::setw(15) << std::fixed << std::setprecision(6) << participation.participation_ratios[i]
                 << std::setw(15) << std::fixed << std::setprecision(2) << participation.localization_lengths[i]
                 << "\n";
            
            phonon_modes.push_back(i);
        }
    }
    
    file.close();
    std::cout << "Exported " << phonon_modes.size() << " phonon modes to: " << filename << std::endl;
    
    return phonon_modes;
}

std::vector<int> FEMHessianAssembler::exportLocalizedModes(
    const std::string& filename,
    const EigenResults& eigen_results,
    const ParticipationAnalysis& participation,
    int num_rigid_body)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return {};
    }
    
    file << "# Localized modes (non-phonon modes with P < " << participation.threshold << ")\n";
    file << "# Columns: Mode_Index, Eigenvalue(λ), Frequency(ω), Lambda_Squared(λ²), P_ratio, L_eff\n";
    file << "#\n";
    file << std::setw(12) << "Mode"
         << std::setw(20) << "Lambda"
         << std::setw(20) << "Omega"
         << std::setw(20) << "Lambda_Squared"
         << std::setw(15) << "P_ratio"
         << std::setw(15) << "L_eff"
         << "\n";
    
    std::vector<int> localized_modes;
    
    for (int i = num_rigid_body; i < eigen_results.num_computed; i++) {
        if (!participation.is_extended[i]) {
            double lambda = eigen_results.eigenvalues(i);
            double omega = std::sqrt(std::abs(lambda));
            if (lambda < 0) omega = -omega;
            double lambda_sq = lambda * lambda;
            
            file << std::setw(12) << i
                 << std::setw(20) << std::scientific << std::setprecision(10) << lambda
                 << std::setw(20) << std::scientific << std::setprecision(10) << omega
                 << std::setw(20) << std::scientific << std::setprecision(10) << lambda_sq
                 << std::setw(15) << std::fixed << std::setprecision(6) << participation.participation_ratios[i]
                 << std::setw(15) << std::fixed << std::setprecision(2) << participation.localization_lengths[i]
                 << "\n";
            
            localized_modes.push_back(i);
        }
    }
    
    file.close();
    std::cout << "Exported " << localized_modes.size() << " localized modes to: " << filename << std::endl;
    
    return localized_modes;
}

int FEMHessianAssembler::detectRigidBodyModes(
    const EigenResults& eigen_results,
    double threshold)
{
    int num_rigid_body = 0;
    for (int i = 0; i < eigen_results.num_computed; i++) {
        if (std::abs(eigen_results.eigenvalues(i)) < threshold) {
            num_rigid_body++;
        } else {
            break;  // Assumes sorted
        }
    }
    return num_rigid_body;
}

void FEMHessianAssembler::exportCompleteEigenmodeAnalysis(
    const std::string& base_filename,
    const EigenResults& eigen_results,
    const std::vector<Point2D>& current_points,
    const std::vector<std::pair<int, int>>& dof_mapping,
    const std::vector<ElementTriangle2D>& elements,
    int num_rigid_body,
    double participation_threshold,
    double vtk_scale_factor)
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "COMPLETE EIGENMODE ANALYSIS AND EXPORT" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Export all non-rigid eigenvalues
    std::string eigen_filename = base_filename + "_eigenvalues_nonrigid.dat";
    exportNonRigidEigenvalues(eigen_filename, eigen_results, num_rigid_body);
    
    // Perform participation ratio analysis
    std::cout << "\nPerforming participation ratio analysis..." << std::endl;
    ParticipationAnalysis participation = analyzeParticipationRatios(
        eigen_results,
        dof_mapping,
        current_points.size(),
        participation_threshold
    );
    
    // Print summary
    printParticipationSummary(participation, eigen_results);
    
    // Export participation analysis
    std::string participation_filename = base_filename + "_participation_analysis.dat";
    exportParticipationAnalysis(participation_filename, participation, eigen_results);
    
    // Export phonon modes
    std::string phonon_filename = base_filename + "_phonon_modes.dat";
    std::vector<int> phonon_modes = exportPhononModes(
        phonon_filename, eigen_results, participation, num_rigid_body
    );
    
    // Export localized modes
    std::string localized_filename = base_filename + "_localized_modes.dat";
    std::vector<int> localized_modes = exportLocalizedModes(
        localized_filename, eigen_results, participation, num_rigid_body
    );
    
    // NEW: Export density of states
    exportAllDensityOfStates(
        base_filename,
        eigen_results,
        participation,
        num_rigid_body,
        20,      // Number of bins
        true    // Gaussian broadening (set to true for smoother curves)
    );
    
    // Export VTK files for phonon modes
    if (!phonon_modes.empty()) {
        std::cout << "\nExporting phonon mode VTK files..." << std::endl;
        exportEigenvectorsToVTK(
            base_filename + "_phonon",
            current_points,
            eigen_results,
            phonon_modes,
            vtk_scale_factor,
            dof_mapping,
            elements
        );
    }
    
    // Export VTK files for localized modes
    if (!localized_modes.empty()) {
        std::cout << "\nExporting localized mode VTK files..." << std::endl;
        exportEigenvectorsToVTK(
            base_filename + "_localized",
            current_points,
            eigen_results,
            localized_modes,
            vtk_scale_factor,
            dof_mapping,
            elements
        );
    }
    
    // Summary statistics
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "EXPORT SUMMARY" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    std::cout << "Total modes computed: " << eigen_results.num_computed << std::endl;
    std::cout << "Rigid body modes: " << num_rigid_body << std::endl;
    std::cout << "Non-rigid modes: " << (eigen_results.num_computed - num_rigid_body) << std::endl;
    std::cout << "  - Phonon modes (extended): " << phonon_modes.size() << std::endl;
    std::cout << "  - Localized modes: " << localized_modes.size() << std::endl;
    std::cout << "  - Fraction phonon: " 
              << static_cast<double>(phonon_modes.size()) / (eigen_results.num_computed - num_rigid_body) 
              << std::endl;
    std::cout << "\nFiles created:" << std::endl;
    std::cout << "  " << eigen_filename << std::endl;
    std::cout << "  " << participation_filename << std::endl;
    std::cout << "  " << phonon_filename << std::endl;
    std::cout << "  " << localized_filename << std::endl;
    std::cout << "  " << base_filename << "_dos_all_modes.dat" << std::endl;
    std::cout << "  " << base_filename << "_dos_phonon_modes.dat" << std::endl;
    std::cout << "  " << base_filename << "_dos_localized_modes.dat" << std::endl;
    std::cout << "  " << base_filename << "_phonon_mode_*.vtk (" << phonon_modes.size() << " files)" << std::endl;
    std::cout << "  " << base_filename << "_localized_mode_*.vtk (" << localized_modes.size() << " files)" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
}

DOSResults FEMHessianAssembler::computeDensityOfStates(
    const EigenResults& eigen_results,
    const std::vector<int>& mode_indices,
    int num_bins,
    bool use_gaussian_broadening,
    double sigma)
{
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
            if (lambda < 0) omega = -omega;  // Keep sign for unstable modes
            frequencies.push_back(omega);
        }
    }
    
    if (frequencies.empty()) {
        std::cerr << "Warning: No valid frequencies for DOS" << std::endl;
        return dos_results;
    }
    
    // Find frequency range
    dos_results.omega_min = *std::min_element(frequencies.begin(), frequencies.end());
    dos_results.omega_max = *std::max_element(frequencies.begin(), frequencies.end());
    
    // Add small padding to range
    double range = dos_results.omega_max - dos_results.omega_min;
    dos_results.omega_min -= 0.05 * range;
    dos_results.omega_max += 0.05 * range;
    
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
            sigma = bin_width * 2.0;  // Auto sigma = 2 * bin width
        }
        
        std::cout << "  Using Gaussian broadening with σ = " << sigma << std::endl;
        
        for (double freq : frequencies) {
            for (int i = 0; i < num_bins; i++) {
                double omega = dos_results.omega_bins[i];
                double gauss = std::exp(-0.5 * std::pow((omega - freq) / sigma, 2)) 
                              / (sigma * std::sqrt(2.0 * M_PI));
                dos_results.dos[i] += gauss;
            }
        }
    } else {
        // Simple histogram binning
        for (double freq : frequencies) {
            int bin_idx = static_cast<int>((freq - dos_results.omega_min) / bin_width);
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
        dos_results.dos_cumulative[i] = dos_results.dos_cumulative[i-1] + dos_results.dos[i] * bin_width;
    }
    
    return dos_results;
}

void FEMHessianAssembler::exportDensityOfStates(
    const std::string& filename,
    const DOSResults& dos_results,
    const std::string& label)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }
    
    file << "# Density of States: " << label << "\n";
    file << "# Number of bins: " << dos_results.num_bins << "\n";
    file << "# Frequency range: [" << dos_results.omega_min << ", " << dos_results.omega_max << "]\n";
    file << "# Columns: Omega(ω), g(ω), Cumulative_DOS, log10 Omega(ω), log10 g(ω), log10 Cumulative_DOS\n";
    file << "#\n";
    file << std::setw(20) << "Omega"
         << std::setw(20) << "g(omega)"
         << std::setw(20) << "Cumulative"
         << std::setw(20) << "log10 Omega"
         << std::setw(20) << "log10 g(omega)"
         << std::setw(20) << "log10 Cumulative"

         << "\n";
    
    for (int i = 0; i < dos_results.num_bins; i++) {
        file << std::setw(20) << std::scientific << std::setprecision(10) << dos_results.omega_bins[i]
             << std::setw(20) << std::scientific << std::setprecision(10) << dos_results.dos[i]
             << std::setw(20) << std::fixed << std::setprecision(10) << dos_results.dos_cumulative[i]
             << std::setw(20) << std::fixed << std::setprecision(10) << std::log10(dos_results.omega_bins[i])
             << std::setw(20) << std::fixed << std::setprecision(10) << std::log10(dos_results.dos[i])
             << std::setw(20) << std::fixed << std::setprecision(10) << std::log10(dos_results.dos_cumulative[i])

             << "\n";
    }
    
    file.close();
    std::cout << "Exported DOS to: " << filename << std::endl;
}

void FEMHessianAssembler::exportAllDensityOfStates(
    const std::string& base_filename,
    const EigenResults& eigen_results,
    const ParticipationAnalysis& participation,
    int num_rigid_body,
    int num_bins,
    bool use_gaussian_broadening)
{
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
    std::cout << "Gaussian broadening: " << (use_gaussian_broadening ? "Yes" : "No") << std::endl;
    
    // Compute DOS for all modes
    std::cout << "\nComputing DOS for all modes..." << std::endl;
    DOSResults dos_all = computeDensityOfStates(
        eigen_results, all_modes, num_bins, use_gaussian_broadening
    );
    exportDensityOfStates(
        base_filename + "_dos_all_modes.dat", dos_all, "All non-rigid modes"
    );
    
    // Compute DOS for phonon modes
    if (!phonon_modes.empty()) {
        std::cout << "\nComputing DOS for phonon modes..." << std::endl;
        DOSResults dos_phonon = computeDensityOfStates(
            eigen_results, phonon_modes, num_bins, use_gaussian_broadening
        );
        exportDensityOfStates(
            base_filename + "_dos_phonon_modes.dat", dos_phonon, "Phonon (extended) modes"
        );
    }
    
    // Compute DOS for localized modes
    if (!localized_modes.empty()) {
        std::cout << "\nComputing DOS for localized modes..." << std::endl;
        DOSResults dos_localized = computeDensityOfStates(
            eigen_results, localized_modes, num_bins, use_gaussian_broadening
        );
        exportDensityOfStates(
            base_filename + "_dos_localized_modes.dat", dos_localized, "Localized modes"
        );
    }
    
    std::cout << std::string(70, '=') << std::endl;
}