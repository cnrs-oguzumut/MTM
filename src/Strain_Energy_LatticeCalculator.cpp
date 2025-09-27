// Strain_Energy_LatticeCalculator.cpp
#include "../include/lattice_energy/Strain_Energy_LatticeCalculator.h"
#include "../include/lattice_energy/EnergyFunctions.h"
#include <cmath>

Strain_Energy_LatticeCalculator::Strain_Energy_LatticeCalculator(double scale) :
    scal(scale),
    burgers(scale),
    normalisation(std::pow(burgers, 2.0)),
    rcut(2.5),  // Added for compatibility
    nb_atoms(10) // Added for compatibility
{}

// Implementation using distance-based potential (for compatibility with base class)
double Strain_Energy_LatticeCalculator::calculate_energy(const Eigen::Matrix2d& C,
                                                        const std::function<double(double)>& pot,
                                                        double zero) {
    // For strain energy, we bypass the distance-based potential and use analytical formula
    return calculate_energy_analytical(C, energy_functions::phi_func, zero);
}

// Implementation using distance-based potential derivative (for compatibility with base class)
Eigen::Matrix2d Strain_Energy_LatticeCalculator::calculate_derivative(const Eigen::Matrix2d& C,
                                                                     const std::function<double(double)>& dpot) {
    // For strain energy, we bypass the distance-based derivative and use analytical formula
    return calculate_derivative_analytical(C, energy_functions::dphi_func);
}

// UNIFORM INTERFACE: Calculate second derivatives as 4th order ITensor
itensor::ITensor Strain_Energy_LatticeCalculator::calculate_dseconderivative(const Eigen::Matrix2d& C,
                                                                           const std::function<double(double)>& dpot,
                                                                           const std::function<double(double)>& d2pot) {
    // For strain energy, we ignore dpot and d2pot and use analytical second derivative function
    // This provides the uniform interface while using analytical functions internally
    return energy_functions::ddphi_func(C);
}

// Direct analytical energy calculation
double Strain_Energy_LatticeCalculator::calculate_energy_analytical(const Eigen::Matrix2d& C,
                                                                   const std::function<double(const Eigen::Matrix2d&)>& phi_func,
                                                                   double zero) {
    return phi_func(C) - zero;
}

// Direct analytical derivative calculation
Eigen::Matrix2d Strain_Energy_LatticeCalculator::calculate_derivative_analytical(const Eigen::Matrix2d& C,
                                                                                const std::function<Eigen::Matrix2d(const Eigen::Matrix2d&)>& dphi_func) {
    return dphi_func(C);
}

// Lagrange reduction method
Eigen::Matrix2d Strain_Energy_LatticeCalculator::lagrange_reduction(const Eigen::Matrix2d& C) {
    // Implementation of Lagrange reduction algorithm
    // This is a placeholder - you would implement the actual algorithm here
    return C;
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

// Helper method to get nearest neighbor distance
double Strain_Energy_LatticeCalculator::getNearestNeighborDistance() const {
    return scal;
}

// Set scale factor
void Strain_Energy_LatticeCalculator::setScale(double scale) {
    scal = scale;
    burgers = scale;
    normalisation = std::pow(burgers, 2.0);
}