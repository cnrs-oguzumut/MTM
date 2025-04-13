// DeformationCalculator.cpp
#include "../include/loading/DeformationCalculator.h"

DeformationCalculator::DeformationCalculator() {
    // Initialize stiffness matrix (for linear elasticity)
    stiffness << 110.9363766697, 36.9787922232, 0,
                 36.9787922232, 110.9363766697, 0,
                 0, 0, 36.9787922232;
}

DeformationCalculator::DeformationCalculator(const Eigen::Matrix3d& customStiffness) 
    : stiffness(customStiffness) {}

Eigen::Matrix2d DeformationCalculator::calculateDeformationGradient(
    const Eigen::Matrix2d& sigma, const Eigen::Matrix2d& cauchy_applied) const {
    
    // Convert Cauchy stress to Voigt notation
    Eigen::Vector3d sigma_voigt(sigma(0,0) + cauchy_applied(0,0), 
                                sigma(1,1) + cauchy_applied(1,1), 
                                sigma(0,1) + cauchy_applied(0,1));
    
    // Compute infinitesimal strain in Voigt notation
    Eigen::Vector3d eps_voigt = -stiffness.inverse() * sigma_voigt;
    
    // Convert strain to tensor form
    Eigen::Matrix2d eps;
    eps << eps_voigt(0), eps_voigt(2)/2.0,
           eps_voigt(2)/2.0, eps_voigt(1);
    
    // For large deformations, compute the right Cauchy-Green tensor C
    // C = I + 2E where E is the Green-Lagrange strain tensor
    // For moderate strains, we can approximate E ≈ ε (the infinitesimal strain)
    Eigen::Matrix2d C = Eigen::Matrix2d::Identity() + 2.0 * eps;
    
    // Eigendecomposition of C
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> solver(C);
    Eigen::Vector2d eigenvalues = solver.eigenvalues();
    Eigen::Matrix2d eigenvectors = solver.eigenvectors();
    
    // Compute principal stretches (square root of eigenvalues of C)
    Eigen::Vector2d stretches = eigenvalues.cwiseSqrt();
    
    // Construct the right stretch tensor U
    Eigen::Matrix2d U = eigenvectors * stretches.asDiagonal() * eigenvectors.transpose();
    
    // For pure stretch (no rotation), F = U
    return U;
}