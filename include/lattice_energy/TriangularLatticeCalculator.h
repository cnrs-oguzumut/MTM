#ifndef TRIANGULAR_LATTICE_CALCULATOR_H
#define TRIANGULAR_LATTICE_CALCULATOR_H

#include "BaseLatticeCalculator.h"
#include "itensor/all.h"
#include <Eigen/Dense>
#include <functional>
#include <vector>
#include <array>

/**
 * @class TriangularLatticeCalculator
 * 
 * ARCHITECTURAL USAGE NOTE:
 * -------------------------
 * TriangularLatticeCalculator evaluates pair potentials directly on the physical
 * triangular lattice:
 *   p = m * [1, 0] + n * [1/2, sqrt(3)/2],   r = scal * sqrt(p^T * C * p),
 * using an adaptive summation window over (m, n) that covers the cutoff ellipse
 * for any positive-definite metric C (making the lattice sum GL(2, Z)-invariant).
 *
 * VALID USAGE:
 * 1. Triangular-reference mode (--triangular-reference=1, ReferenceFrame::Triangular):
 *    - Reference parent triangle is equilateral: (0,0), (a_0, 0), (a_0/2, a_0*sqrt(3)/2),
 *      so an undeformed crystal has C = C_ref = I and F_0 = I.
 *    - Lagrange reduction is NOT applied (reduce_metric returns m = I).
 *    - Unit-cell area normalisation is getUnitCellArea() = (sqrt(3)/2) * a_0^2.
 * 2. Direct physical-coordinate evaluations at C = I (e.g. Volterra elastic constants,
 *    nanoindentation).
 *
 * DO NOT MIX WITH Z^2 METRICS OR LAGRANGE REDUCTION:
 * - Never pass the Z^2 reference metric C_ref = F_0^T * F_0 or use right-isosceles
 *   shape derivatives with TriangularLatticeCalculator (that applies the triangular
 *   distortion twice).
 * - Never apply Z^2 Lagrange reduction (lagrange::reduce) before passing C to
 *   TriangularLatticeCalculator. For Z^2-reference + reduction studies, use
 *   SquareLatticeCalculator(a_0 / (4/3)^(1/4), r_cutoff).
 */
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

    HessianComponents calculate_dseconderivative_components(const Eigen::Matrix2d& C,
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