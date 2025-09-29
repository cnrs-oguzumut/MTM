#ifndef TRIANGULAR_LATTICE_CALCULATOR_H
#define TRIANGULAR_LATTICE_CALCULATOR_H

#include "BaseLatticeCalculator.h"
#include "itensor/all.h"
#include <Eigen/Dense>
#include <functional>
#include <vector>
#include <array>

class TriangularLatticeCalculator : public BaseLatticeCalculator {
public:
    TriangularLatticeCalculator(double scale,double cutoff_radius=2.5);
    
    double calculate_energy(const Eigen::Matrix2d& C,
                           const std::function<double(double)>& pot,
                           double zero) override;
                           
    Eigen::Matrix2d calculate_derivative(const Eigen::Matrix2d& C,
                                       const std::function<double(double)>& dpot) override;

    /**
     * Calculate second derivatives as 4th order ITensor (for acoustic tensor calculations)
     * UNIFORM INTERFACE - same signature as other lattice calculators
     * 
     * @param C Metric tensor (2x2)
     * @param dpot First derivative function (distance → dE/dr)
     * @param d2pot Second derivative function (distance → d²E/dr²)
     * @return 4th order ITensor representing ∂²E/∂C_ij∂C_kl
     */
    itensor::ITensor calculate_dseconderivative(const Eigen::Matrix2d& C,
                                              const std::function<double(double)>& dpot,
                                              const std::function<double(double)>& d2pot) override;
                                       
    double getNearestNeighborDistance() const override;
    
    double getUnitCellArea() const override;
    
private:
    const double rcut;
    const double scal;
    const double burgers;
    const int nb_atoms;
    const double normalisation;
    const double r_cutoff_sq;  // Added for efficiency (matching SquareLatticeCalculator)
    std::vector<std::vector<double>> triangular_basis;  // Changed from array to vector for consistency

    // Private helper methods (same pattern as SquareLatticeCalculator)
    HessianComponents calculate_hessian_components(const Eigen::Matrix2d& C,
                                                 const std::function<double(double)>& pot,
                                                 const std::function<double(double)>& dpot,
                                                 const std::function<double(double)>& d2pot);

    itensor::ITensor hessianComponentsToITensor(const HessianComponents& hess) const;

    /**
     * Calculate acoustic tensor components directly (alternative interface)
     * Returns the complete 4th-order acoustic tensor flattened to 2x2 blocks
     */
    std::vector<Eigen::Matrix2d> calculateAcousticTensorBlocks(const Eigen::Matrix2d& C,
                                                              const std::function<double(double)>& pot,
                                                              const std::function<double(double)>& dpot,
                                                              const std::function<double(double)>& d2pot);
};

#endif // TRIANGULAR_LATTICE_CALCULATOR_H