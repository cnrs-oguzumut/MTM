#include "../include/lattice_energy/Strain_Energy_LatticeCalculator.h"
#include "../include/lattice_energy/EnergyFunctions.h"

// Constructor - matching SquareLatticeCalculator constructor
Strain_Energy_LatticeCalculator::Strain_Energy_LatticeCalculator(double scale) :
    rcut(2.5),
    scal(scale),
    burgers(scale),
    nb_atoms(10),
    normalisation(std::pow(burgers, 2.0)) {} // Unit cell area is a²

// Calculate energy - same signature as SquareLatticeCalculator
double Strain_Energy_LatticeCalculator::calculate_energy(const Eigen::Matrix2d& C,
                       const std::function<double(double)>& pot,
                       double zero) {
    // Apply Lagrange reduction
    //Eigen::Matrix2d C_reduced = lagrange_reduction(C);
    
    // Use the analytical energy function
    double phi = energy_functions::phi_func(C);
    return phi - zero;
}

// Calculate energy derivative - same signature as SquareLatticeCalculator
Eigen::Matrix2d Strain_Energy_LatticeCalculator::calculate_derivative(const Eigen::Matrix2d& C,
                                    const std::function<double(double)>& dpot) {
    // Apply Lagrange reduction
    //Eigen::Matrix2d C_reduced = lagrange_reduction(C);
    
    // Use the analytical derivative function
    return energy_functions::dphi_func(C);
}

// New analytical energy calculation method
double Strain_Energy_LatticeCalculator::calculate_energy_analytical(const Eigen::Matrix2d& C,
                       const std::function<double(const Eigen::Matrix2d&)>& phi_func,
                       double zero) {
    // Apply Lagrange reduction
    //Eigen::Matrix2d C_reduced = lagrange_reduction(C);
    
    // Apply the analytical function directly to the metric tensor
    double phi = phi_func(C);
    return phi - zero;
}

// New analytical derivative calculation method
Eigen::Matrix2d Strain_Energy_LatticeCalculator::calculate_derivative_analytical(const Eigen::Matrix2d& C,
                                    const std::function<Eigen::Matrix2d(const Eigen::Matrix2d&)>& dphi_func) {
    // Apply Lagrange reduction
    //Eigen::Matrix2d C_reduced = lagrange_reduction(C);
    
    // Apply the analytical derivative function directly to the metric tensor
    return dphi_func(C);
}

// Lagrange reduction for metric tensor C
Eigen::Matrix2d Strain_Energy_LatticeCalculator::lagrange_reduction(const Eigen::Matrix2d& C) {
    double c11 = C(0, 0);
    double c22 = C(1, 1);
    double c12 = C(0, 1);
    
    // Iteratively apply Lagrange reduction
    while (c12 < 0 || c22 < c11 || 2 * c12 > c11) {
        if (c12 < 0) {
            c12 = -c12;
        }
        
        if (c22 < c11) {
            double temp = c11;
            c11 = c22;
            c22 = temp;
        }
        
        if (2 * c12 > c11) {
            double d11 = c11;
            double d12 = c12 - c11;
            double d22 = c22 + c11 - 2 * c12;
            
            c11 = d11;
            c12 = d12;
            c22 = d22;
        }
    }
    
    Eigen::Matrix2d C_reduced;
    C_reduced(0, 0) = c11;
    C_reduced(1, 1) = c22;
    C_reduced(0, 1) = c12;
    C_reduced(1, 0) = c12;  // Symmetric matrix
    
    return C_reduced;
}

// Helper method to get scale factor
double Strain_Energy_LatticeCalculator::getScale() const {
    return scal;
}

// Helper method to get Burgers vector magnitude
double Strain_Energy_LatticeCalculator::getBurgersVector() const {
    return burgers;
}

// Helper method to get area of unit cell
double Strain_Energy_LatticeCalculator::getUnitCellArea() const {
    return normalisation;
}

// Helper method for compatibility with SquareLatticeCalculator
double Strain_Energy_LatticeCalculator::getNearestNeighborDistance() const {
    return scal;
}

// Set scale factor
void Strain_Energy_LatticeCalculator::setScale(double scale) {
    scal = scale;
    burgers = scale;
    normalisation = std::pow(burgers, 2.0);
}