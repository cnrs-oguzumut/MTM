// BaseLatticeCalculator.h
#ifndef BASE_LATTICE_CALCULATOR_H
#define BASE_LATTICE_CALCULATOR_H

#include <Eigen/Dense>
#include <functional>

class BaseLatticeCalculator {
public:
    virtual ~BaseLatticeCalculator() = default;
    
    virtual double calculate_energy(const Eigen::Matrix2d& C,
                                   const std::function<double(double)>& pot,
                                   double zero) = 0;
                                   
    virtual Eigen::Matrix2d calculate_derivative(const Eigen::Matrix2d& C,
                                               const std::function<double(double)>& dpot) = 0;
                                               
    virtual double getNearestNeighborDistance() const = 0;
    
    virtual double getUnitCellArea() const = 0;
};

#endif // BASE_LATTICE_CALCULATOR_H