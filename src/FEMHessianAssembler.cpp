#include "FEMHessianAssembler.h"
#include <iostream>

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
    const int num_nodes = 3;      // Triangular element
    const int spatial_dim = 2;     // 2D
    const int total_dofs = num_nodes * spatial_dim;  // 6 DOFs
    
    Eigen::MatrixXd K_elem = Eigen::MatrixXd::Zero(total_dofs, total_dofs);
    
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
    std::vector<AcousticTensor>& acoustic_tensors,
    int num_total_dofs,
    const std::vector<std::pair<int, int>>& dof_mapping)
{
    if (elements.size() != acoustic_tensors.size()) {
        std::cerr << "Error: Number of elements and acoustic tensors must match!" << std::endl;
        throw std::runtime_error("Size mismatch in global assembly");
    }
    
    // Create triplet list for sparse matrix assembly
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(elements.size() * 36); // 6x6 entries per element
    
    for (size_t elem_idx = 0; elem_idx < elements.size(); elem_idx++) {
        auto& element = elements[elem_idx];
        auto& acoustic_tensor = acoustic_tensors[elem_idx];
        
        // Get element area for integration weight
        double area = element.getArea();
        
        // Compute element stiffness matrix
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