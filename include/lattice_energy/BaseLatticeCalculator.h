// BaseLatticeCalculator.h
#ifndef BASE_LATTICE_CALCULATOR_H
#define BASE_LATTICE_CALCULATOR_H

#include <Eigen/Dense>
#include <functional>
#include "itensor/all.h"

struct HessianComponents;

class BaseLatticeCalculator {
public:
    virtual ~BaseLatticeCalculator() = default;
    
    virtual double calculate_energy(const Eigen::Matrix2d& C,
                                   const std::function<double(double)>& pot,
                                   double zero) = 0;
                                   
    virtual Eigen::Matrix2d calculate_derivative(const Eigen::Matrix2d& C,
                                               const std::function<double(double)>& dpot) = 0;

    virtual itensor::ITensor calculate_dseconderivative(const Eigen::Matrix2d& C,
                                                      const std::function<double(double)>& dpot,
                                                      const std::function<double(double)>& d2pot) = 0;

    // Same second derivatives as calculate_dseconderivative, returned as the six
    // independent raw partials (c11, c22, c12) without building an ITensor.
    // The default reads them back from the ITensor; calculators override it to
    // skip the ITensor entirely (used by the fast FEM Hessian assembly).
    virtual HessianComponents calculate_dseconderivative_components(
        const Eigen::Matrix2d& C,
        const std::function<double(double)>& dpot,
        const std::function<double(double)>& d2pot);

    virtual double getNearestNeighborDistance() const = 0;
    
    virtual double getUnitCellArea() const = 0;
};

struct HessianComponents {
    double c11_c11;  // ∂²E/∂c₁₁²
    double c22_c22;  // ∂²E/∂c₂₂²
    double c12_c12;  // ∂²E/∂c₁₂²
    double c11_c22;  // ∂²E/∂c₁₁∂c₂₂
    double c11_c12;  // ∂²E/∂c₁₁∂c₁₂
    double c22_c12;  // ∂²E/∂c₂₂∂c₁₂
    
    HessianComponents() : c11_c11(0), c22_c22(0), c12_c12(0),
                         c11_c22(0), c11_c12(0), c22_c12(0) {}
};

inline HessianComponents BaseLatticeCalculator::calculate_dseconderivative_components(
    const Eigen::Matrix2d& C,
    const std::function<double(double)>& dpot,
    const std::function<double(double)>& d2pot) {
    const itensor::ITensor t = calculate_dseconderivative(C, dpot, d2pot);
    HessianComponents hess;
    hess.c11_c11 = t.elt(1, 1, 1, 1);
    hess.c22_c22 = t.elt(2, 2, 2, 2);
    hess.c12_c12 = t.elt(1, 2, 1, 2);
    hess.c11_c22 = t.elt(1, 1, 2, 2);
    hess.c11_c12 = t.elt(1, 1, 1, 2);
    hess.c22_c12 = t.elt(2, 2, 1, 2);
    return hess;
}

#endif // BASE_LATTICE_CALCULATOR_H