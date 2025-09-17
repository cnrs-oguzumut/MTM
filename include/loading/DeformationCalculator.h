#ifndef DEFORMATION_CALCULATOR_H
#define DEFORMATION_CALCULATOR_H

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

/**
 * @brief Class for calculating deformation gradients from stress states
 * 
 * This class implements the computation of deformation gradients using
 * linear elasticity theory and finite strain kinematics. It converts
 * Cauchy stress to strain via the compliance matrix and then computes
 * the right stretch tensor as an approximation to the deformation gradient.
 */
class DeformationCalculator {
private:
    Eigen::Matrix3d stiffness;  ///< Material stiffness matrix in Voigt notation

public:
    /**
     * @brief Default constructor with predefined stiffness matrix
     * 
     * Initializes the calculator with a default stiffness matrix for
     * an isotropic or orthotropic material.
     */
    DeformationCalculator();
    
    /**
     * @brief Constructor with custom stiffness matrix
     * 
     * @param customStiffness 3x3 stiffness matrix in Voigt notation
     *                        [C11 C12 C16]
     *                        [C12 C22 C26] 
     *                        [C16 C26 C66]
     */
    explicit DeformationCalculator(const Eigen::Matrix3d& customStiffness);
    
    /**
     * @brief Calculate deformation gradient from stress state
     * 
     * Computes the deformation gradient tensor by:
     * 1. Converting Cauchy stress to infinitesimal strain via compliance
     * 2. Approximating Green-Lagrange strain with infinitesimal strain
     * 3. Computing right Cauchy-Green tensor C = I + 2E
     * 4. Extracting right stretch tensor U from C via eigendecomposition
     * 
     * @param sigma Base stress state (2x2 Cauchy stress tensor)
     * @param cauchy_applied Applied stress increment (2x2 Cauchy stress tensor)
     * @return 2x2 deformation gradient tensor (approximated as right stretch tensor U)
     * 
     * @note This implementation assumes moderate strains where the approximation
     *       E ≈ ε is valid. For large deformations, a more sophisticated approach
     *       involving iterative stress-strain relationships would be needed.
     */
    Eigen::Matrix2d calculateDeformationGradient(
        const Eigen::Matrix2d& sigma, 
        const Eigen::Matrix2d& cauchy_applied) const;
    
    /**
     * @brief Get the current stiffness matrix
     * @return const reference to the 3x3 stiffness matrix
     */
    const Eigen::Matrix3d& getStiffnessMatrix() const { return stiffness; }
    
    /**
     * @brief Set a new stiffness matrix
     * @param newStiffness 3x3 stiffness matrix in Voigt notation
     */
    void setStiffnessMatrix(const Eigen::Matrix3d& newStiffness) { 
        stiffness = newStiffness; 
    }
};

#endif // DEFORMATION_CALCULATOR_H