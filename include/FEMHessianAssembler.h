#ifndef FEM_HESSIAN_ASSEMBLY_H
#define FEM_HESSIAN_ASSEMBLY_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include "acoustic_tensor.h"
#include "../include/mesh/ElementTriangle2D.h"


class FEMHessianAssembler {
public:
    /**
     * Compute element stiffness matrix K^{ab}_{ij} = A_{iKjL} * (dN^a/dX_K) * (dN^b/dX_L)
     * 
     * @param acoustic_tensor The acoustic tensor object containing A_{iKjL}
     * @param element The triangular element with shape function derivatives
     * @param area_weight Integration weight (typically element area for single Gauss point)
     * @return 6x6 element stiffness matrix
     */
    static Eigen::MatrixXd computeElementStiffness(
        AcousticTensor& acoustic_tensor,
        const ElementTriangle2D& element,
        double area_weight = 1.0);
    
    /**
     * Assemble global stiffness matrix from all elements
     * 
     * @param elements Vector of all triangular elements
     * @param acoustic_tensors Vector of acoustic tensors (one per element)
     * @param num_total_dofs Total number of free DOFs in the system
     * @param dof_mapping DOF mapping (original_idx, solver_idx) for each node
     * @return Sparse global stiffness matrix
     */
    static Eigen::SparseMatrix<double> assembleGlobalStiffness(
        std::vector<ElementTriangle2D>& elements,
        std::vector<AcousticTensor>& acoustic_tensors,
        int num_total_dofs,
        const std::vector<std::pair<int, int>>& dof_mapping);

private:
    /**
     * Extract acoustic tensor ITensor to 4D array for easier manipulation
     * 
     * @param A_tensor The ITensor acoustic tensor with indices (r_, i_, s_, j_)
     * @param A Output 4D array A[i][K][j][L]
     */
    static void extractAcousticTensor(
        const itensor::ITensor& A_tensor,
        double A[2][2][2][2]);
};

#endif // FEM_HESSIAN_ASSEMBLY_H