#ifndef SQUARE_LATTICE_CALCULATOR_H
#define SQUARE_LATTICE_CALCULATOR_H

#include "BaseLatticeCalculator.h"
#include <Eigen/Dense>
#include <functional>
#include <vector>

class SquareLatticeCalculator : public BaseLatticeCalculator {
public:
    SquareLatticeCalculator(double scale);
    
    double calculate_energy(const Eigen::Matrix2d& C,
                           const std::function<double(double)>& pot,
                           double zero) override;
                           
    Eigen::Matrix2d calculate_derivative(const Eigen::Matrix2d& C,
                                       const std::function<double(double)>& dpot) override;
                                       
    double getNearestNeighborDistance() const override;
    
    double getUnitCellArea() const override;
    
private:
    const double rcut;
    const double scal;
    const double burgers;
    const int nb_atoms;
    const double normalisation;
    std::vector<std::vector<double>> square_basis;
};

#endif // SQUARE_LATTICE_CALCULATOR_H