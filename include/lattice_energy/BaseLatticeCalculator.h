// BaseLatticeCalculator.h
#ifndef BASE_LATTICE_CALCULATOR_H
#define BASE_LATTICE_CALCULATOR_H

#include <Eigen/Dense>
#include <functional>
#include "itensor/all.h"

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

#endif // BASE_LATTICE_CALCULATOR_H