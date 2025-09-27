#ifndef SQUARE_LATTICE_CALCULATOR_H
#define SQUARE_LATTICE_CALCULATOR_H

#include "BaseLatticeCalculator.h"
#include "itensor/all.h"
#include <Eigen/Dense>
#include <functional>
#include <vector>

/**
 * Structure to hold second derivatives (Hessian) of energy w.r.t. metric tensor components
 * This replaces the old double_cella_stru structure
 */
// struct HessianComponents {
//     double c11_c11;  // ∂²E/∂c₁₁²
//     double c22_c22;  // ∂²E/∂c₂₂²
//     double c12_c12;  // ∂²E/∂c₁₂²
//     double c11_c22;  // ∂²E/∂c₁₁∂c₂₂
//     double c11_c12;  // ∂²E/∂c₁₁∂c₁₂
//     double c22_c12;  // ∂²E/∂c₂₂∂c₁₂
    
//     HessianComponents() : c11_c11(0), c22_c22(0), c12_c12(0), 
//                          c11_c22(0), c11_c12(0), c22_c12(0) {}
// };

class SquareLatticeCalculator : public BaseLatticeCalculator {
public:
    SquareLatticeCalculator(double scale);
    
    double calculate_energy(const Eigen::Matrix2d& C,
                           const std::function<double(double)>& pot,
                           double zero) override;
                           
    Eigen::Matrix2d calculate_derivative(const Eigen::Matrix2d& C,
                                       const std::function<double(double)>& dpot) override;

    /**
     * Calculate second derivatives as 4th order ITensor (for acoustic tensor calculations)
     * This is the main method you'll use: dE2_dC2 = calculator.calculate_dseconderivative(C, dpot, d2pot)
     * 
     * @param C Metric tensor (2x2)
     * @param dpot First derivative function (distance → dE/dr)
     * @param d2pot Second derivative function (distance → d²E/dr²)
     * @return 4th order ITensor representing ∂²E/∂C_ij∂C_kl
     */
    itensor::ITensor calculate_dseconderivative(const Eigen::Matrix2d& C,
                                              const std::function<double(double)>& dpot,
                                              const std::function<double(double)>& d2pot);

    /**
     * Calculate second derivatives (Hessian) of energy w.r.t. metric tensor components
     * This is the acoustic tensor calculation functionality from your old code
     * 
     * @param C Metric tensor (2x2)
     * @param pot Energy function (distance → energy)
     * @param dpot First derivative function (distance → dE/dr)
     * @param d2pot Second derivative function (distance → d²E/dr²)
     * @return HessianComponents containing all second derivatives
     */
    HessianComponents calculate_hessian_components(const Eigen::Matrix2d& C,
                                                 const std::function<double(double)>& pot,
                                                 const std::function<double(double)>& dpot,
                                                 const std::function<double(double)>& d2pot);

    /**
     * Convert HessianComponents to 4th order ITensor for acoustic tensor calculations
     * Indices: (i,j,k,l) where ∂²E/∂C_ij∂C_kl
     * This is the natural representation for tensor operations
     */
    itensor::ITensor hessianComponentsToITensor(const HessianComponents& hess) const;

    /**
     * Calculate acoustic tensor components directly (alternative interface)
     * Returns the complete 4th-order acoustic tensor flattened to 2x2 blocks
     */
    std::vector<Eigen::Matrix2d> calculateAcousticTensorBlocks(const Eigen::Matrix2d& C,
                                                              const std::function<double(double)>& pot,
                                                              const std::function<double(double)>& dpot,
                                                              const std::function<double(double)>& d2pot);
                                       
    double getNearestNeighborDistance() const override;
    
    double getUnitCellArea() const override;
    
private:
    const double rcut;
    const double scal;
    const double burgers;
    const int nb_atoms;
    const double normalisation;
    std::vector<std::vector<double>> square_basis;
    
    // Cutoff radius squared for efficiency (from your old code)
    const double r_cutoff_sq = 4.5 * 4.5;
};

#endif // SQUARE_LATTICE_CALCULATOR_H