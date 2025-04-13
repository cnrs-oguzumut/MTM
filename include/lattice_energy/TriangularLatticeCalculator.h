// TriangularLatticeCalculator.h
#ifndef TRIANGULAR_LATTICE_CALCULATOR_H
#define TRIANGULAR_LATTICE_CALCULATOR_H

#include "BaseLatticeCalculator.h"
#include <array>

class TriangularLatticeCalculator : public BaseLatticeCalculator {
private:
    const double rcut;
    const double scal;
    const double burgers;
    const int nb_atoms;
    const double normalisation;
    const std::array<std::array<double, 2>, 1> triangular_basis;
    
public:
    TriangularLatticeCalculator(double scale = 1.0);
    
    double calculate_energy(const Eigen::Matrix2d& C,
                           const std::function<double(double)>& pot,
                           double zero = 0.0) override;
                           
    Eigen::Matrix2d calculate_derivative(const Eigen::Matrix2d& C,
                                       const std::function<double(double)>& dpot) override;
                                       
    double getNearestNeighborDistance() const override;
    
    double getUnitCellArea() const override;
};

#endif // TRIANGULAR_LATTICE_CALCULATOR_H