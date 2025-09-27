#ifndef ACOUSTIC_TENSOR_H
#define ACOUSTIC_TENSOR_H

#include "itensor/all.h"
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <iomanip>
#include <functional>

// Forward declaration to avoid circular includes
// class Strain_Energy_LatticeCalculator;
#include "../include/lattice_energy/TriangularLatticeCalculator.h"
#include "../include/lattice_energy/SquareLatticeCalculator.h"
#include "../include/lattice_energy/Strain_Energy_LatticeCalculator.h"

/**
 * Structure to hold acoustic tensor analysis results
 */
struct AcousticAnalysis {
    double detAc;  // Determinant of acoustic tensor
    double xsi;    // Angle at minimum
};

/**
 * Class for comprehensive acoustic tensor calculations from physical quantities.
 * Replaces struct-based approach with Eigen matrices.
 * Supports both first derivatives (2x2 matrix) and second derivatives (4th order ITensor).
 */
class AcousticTensor {
private:
    Eigen::Matrix2d F_;           // Deformation gradient (2x2) - replaces base_stru v
    Eigen::Matrix2d C_;           // Cell/geometry matrix (2x2) - replaces cella_stru c  
    Eigen::Matrix2d Z_;           // Material properties matrix (2x2) - replaces matrix_stru z
    
    // Additional members for derivative calculation
    Eigen::Matrix2d dE_dC_;         // First derivative of energy w.r.t metric tensor
    itensor::ITensor dE2_dC2_;      // Second derivative (Hessian) of energy w.r.t metric tensor
    double normalisation_;          // Normalization factor
    bool has_first_derivative_;     // Flag to track if first derivative is available
    bool has_second_derivative_;    // Flag to track if second derivative is available
    
    // ITensor indices for various tensors
    itensor::Index i_, j_, k_, l_, m_, n_, r_, s_;
    itensor::Index a_, b_, w1_, w2_, w3_, w4_;

public:
    /**
     * Basic constructor that initializes with physical quantities only.
     * 
     * @param F Deformation gradient (2x2 Eigen matrix) - replaces base_stru v
     * @param C Cell/geometry matrix (2x2 Eigen matrix) - replaces cella_stru c
     * @param Z Material properties matrix (2x2 Eigen matrix) - replaces matrix_stru z
     */
    AcousticTensor(const Eigen::Matrix2d& F, const Eigen::Matrix2d& C, const Eigen::Matrix2d& Z);

    /**
     * Constructor with first derivative
     * 
     * @param F Deformation gradient
     * @param C Cell/geometry matrix  
     * @param Z Material properties matrix
     * @param dE_dC First derivative of energy w.r.t metric tensor
     * @param normalisation Normalization factor
     */
    AcousticTensor(const Eigen::Matrix2d& F, const Eigen::Matrix2d& C, const Eigen::Matrix2d& Z,
                   const Eigen::Matrix2d& dE_dC, double normalisation);

    /**
     * Constructor with both first and second derivatives
     * 
     * @param F Deformation gradient
     * @param C Cell/geometry matrix  
     * @param Z Material properties matrix
     * @param dE_dC First derivative of energy w.r.t metric tensor
     * @param dE2_dC2 Second derivative (Hessian) of energy w.r.t metric tensor (4th order ITensor)
     * @param normalisation Normalization factor
     */
    AcousticTensor(const Eigen::Matrix2d& F, const Eigen::Matrix2d& C, const Eigen::Matrix2d& Z,
                   const Eigen::Matrix2d& dE_dC, const itensor::ITensor& dE2_dC2, double normalisation);

    /**
     * Set the first energy derivative only
     */
    void setEnergyDerivative(const Eigen::Matrix2d& dE_dC, double normalisation);

    /**
     * Set both first and second derivatives
     */
    void setEnergyDerivatives(const Eigen::Matrix2d& dE_dC, const itensor::ITensor& dE2_dC2, double normalisation);

    /**
     * Compute first energy derivative using any calculator (uniform interface)
     * Works with all lattice calculators: Square, Triangular, and Strain_Energy
     * 
     * @param calculator Reference to lattice calculator
     * @param derivative_function Energy derivative function (scalar input)
     * @param normalisation Normalization factor
     */
    template<typename CalculatorType>
    void computeEnergyDerivative(CalculatorType& calculator,
                                const std::function<double(double)>& derivative_function,
                                double normalisation) {
        dE_dC_ = calculator.calculate_derivative(C_, derivative_function) / normalisation;
        normalisation_ = normalisation;
        has_first_derivative_ = true;
    }

    /**
     * Compute both first and second energy derivatives using any calculator (uniform interface)
     * Works with all lattice calculators: Square, Triangular, and Strain_Energy
     * 
     * @param calculator Reference to lattice calculator  
     * @param derivative_function First derivative function (scalar input)
     * @param second_derivative_function Second derivative function (scalar input)
     * @param normalisation Normalization factor
     */
    template<typename CalculatorType>
    void computeEnergyDerivatives(CalculatorType& calculator,
                                 const std::function<double(double)>& derivative_function,
                                 const std::function<double(double)>& second_derivative_function,
                                 double normalisation) {
        dE_dC_ = calculator.calculate_derivative(C_, derivative_function) / normalisation;
        dE2_dC2_ = calculator.calculate_dseconderivative(C_, derivative_function, second_derivative_function) / normalisation;
        normalisation_ = normalisation;
        has_first_derivative_ = true;
        has_second_derivative_ = true;
    }

    /**
     * Compute energy derivative using analytical functions (for Strain_Energy_LatticeCalculator)
     * @param calculator Reference to strain energy calculator
     * @param dphi_func Analytical derivative function (matrix input)
     * @param normalisation Normalization factor
     */
    void computeEnergyDerivativeAnalytical(Strain_Energy_LatticeCalculator& calculator,
                                         const std::function<Eigen::Matrix2d(const Eigen::Matrix2d&)>& dphi_func,
                                         double normalisation);

    /**
     * Create the T54 tensor.
     */
    void createT54(itensor::ITensor& T54);

    /**
     * Create the T53 tensor and F tensor from deformation gradient.
     */
    void createT53(itensor::ITensor& T53, itensor::ITensor& F_tensor);

    /**
     * Compute Jacobian tensor from cell matrix and deformation gradient.
     * Uses first derivatives if available, otherwise fallback calculation.
     */
    void computeJacobian(itensor::ITensor& dPhi);

    /**
     * Compute Hessian tensor from cell matrix and deformation gradient.
     * Uses second derivatives if available, enhanced first derivatives, or fallback calculation.
     */
    void computeHessian(itensor::ITensor& ddPhi);

    /**
     * Main acoustic tensor analysis - equivalent to acoustic_tensor() function.
     * 
     * @param alpha Parameter for analysis
     * @param theta Angle parameter 
     * @return AcousticAnalysis results
     */
    AcousticAnalysis analyzeAcousticTensor(double alpha, double theta);

    /**
     * Compute Aijkl component - equivalent to Aijkl() function.
     * 
     * @param dF Differential deformation gradient matrix
     * @param alpha Parameter for analysis
     * @return Scalar result
     */
    double computeAijkl(const Eigen::Matrix2d& dF, double alpha);

    /**
     * Find minimum determinant of acoustic tensor - equivalent to min_det_acoustic_tensor().
     * 
     * @param xsi Angle parameter (>= 10 for search, < 10 for specific angle)
     * @return AcousticAnalysis results
     */
    AcousticAnalysis findMinDetAcousticTensor(double xsi);

    /**
     * Update deformation gradient and related calculations.
     */
    void updateDeformationGradient(const Eigen::Matrix2d& F);

    /**
     * Update cell matrix and related calculations.
     */
    void updateCellMatrix(const Eigen::Matrix2d& C);

    /**
     * Update material properties matrix.
     */
    void updateMaterialMatrix(const Eigen::Matrix2d& Z);

    /**
     * Get current matrices.
     */
    const Eigen::Matrix2d& getDeformationGradient() const { return F_; }
    const Eigen::Matrix2d& getCellMatrix() const { return C_; }
    const Eigen::Matrix2d& getMaterialMatrix() const { return Z_; }

    /**
     * Get current energy derivative matrix (first derivative)
     */
    const Eigen::Matrix2d& getEnergyDerivative() const { return dE_dC_; }

    /**
     * Get current energy second derivative tensor (Hessian)
     */
    const itensor::ITensor& getEnergySecondDerivative() const { return dE2_dC2_; }

    /**
     * Get normalization factor
     */
    double getNormalization() const { return normalisation_; }

    /**
     * Check if first derivative is available
     */
    bool hasFirstDerivative() const { return has_first_derivative_; }

    /**
     * Check if second derivative is available
     */
    bool hasSecondDerivative() const { return has_second_derivative_; }

    /**
     * Print tensor information and material properties.
     * Shows available derivatives and their status.
     */
    void printInfo() const;

private:
    /**
     * Initialize ITensor indices.
     */
    void initializeIndices();

    /**
     * Convert Eigen matrix to ITensor with specified indices.
     */
    itensor::ITensor eigenToITensor(const Eigen::Matrix2d& mat, const itensor::Index& idx1, const itensor::Index& idx2) const;

    /**
     * Helper function to convert angle to degrees.
     */
    double toDegrees(double radians) const { return radians * 180.0 / M_PI; }

    /**
     * Helper function to find minimum value in vector.
     */
    double findMinValue(const std::vector<double>& vec) const;

    /**
     * Helper function to find index of minimum value in vector.
     */
    size_t findMinIndex(const std::vector<double>& vec) const;
};

#endif // ACOUSTIC_TENSOR_H