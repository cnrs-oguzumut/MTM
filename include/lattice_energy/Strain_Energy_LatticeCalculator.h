#ifndef STRAIN_ENERGY_LATTICE_CALCULATOR_H
#define STRAIN_ENERGY_LATTICE_CALCULATOR_H

#include <Eigen/Dense>
#include <cmath>
#include <functional>

class Strain_Energy_LatticeCalculator {
private:
    double scal;           // Scale factor
    double burgers;        // Burgers vector magnitude
    double normalisation;  // Normalization factor (unit cell area)
    double rcut;           // Cut-off radius (added for compatibility)
    int nb_atoms;          // Number of atoms (added for compatibility)

public:
    // Constructor - matching SquareLatticeCalculator constructor
    Strain_Energy_LatticeCalculator(double scale);
    
    // Methods with exact same signatures as SquareLatticeCalculator
    double calculate_energy(const Eigen::Matrix2d& C,
                           const std::function<double(double)>& pot,
                           double zero);
    
    Eigen::Matrix2d calculate_derivative(const Eigen::Matrix2d& C,
                                        const std::function<double(double)>& dpot);
    
    // New methods that use direct analytical functions
    double calculate_energy_analytical(const Eigen::Matrix2d& C,
                                     const std::function<double(const Eigen::Matrix2d&)>& phi_func,
                                     double zero = 0.0);
    
    Eigen::Matrix2d calculate_derivative_analytical(const Eigen::Matrix2d& C,
                                                  const std::function<Eigen::Matrix2d(const Eigen::Matrix2d&)>& dphi_func);
    
    // Lagrange reduction method that works directly with Eigen::Matrix2d
    Eigen::Matrix2d lagrange_reduction(const Eigen::Matrix2d& C);
    
    // Helper method to get scale factor
    double getScale() const;
    
    // Helper method to get Burgers vector magnitude
    double getBurgersVector() const;
    
    // Helper method to get area of unit cell
    double getUnitCellArea() const;
    
    // Helper methods from SquareLatticeCalculator (for full compatibility)
    double getNearestNeighborDistance() const;
    
    // Set scale factor
    void setScale(double scale);
};

#endif // STRAIN_ENERGY_LATTICE_CALCULATOR_H