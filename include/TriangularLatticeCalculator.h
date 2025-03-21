#pragma once
#ifndef TRIANGULAR_LATTICE_CALCULATOR_H
#define TRIANGULAR_LATTICE_CALCULATOR_H

#include <array>
#include <functional>
#include <cmath>
#include <Eigen/Dense>

class TriangularLatticeCalculator {
private:
    const double rcut;
    const double scal;
    const double burgers;
    const int nb_atoms;
    const double normalisation;
    
    // For a triangular lattice, we only need one basis point
    // The hexagonal pattern comes from the lattice vectors themselves
    const std::array<std::array<double, 2>, 1> triangular_basis;

public:
    TriangularLatticeCalculator(double scale = 1.0);
    
    // Calculate energy
    double calculate_energy(const Eigen::Matrix2d& C, 
                        const std::function<double(double)>& pot, 
                        double zero = 0.0);
    
    // Calculate energy derivative
    Eigen::Matrix2d calculate_derivative(const Eigen::Matrix2d& C,
                                      const std::function<double(double)>& dpot);
    
    // Helper method to get nearest neighbor distance in perfect lattice
    double getNearestNeighborDistance() const;
    
    // Helper method to get area of unit cell
    double getUnitCellArea() const;
};

#endif // TRIANGULAR_LATTICE_CALCULATOR_H
