// SquareLatticeCalculator.cpp
#include "../include/lattice_energy/SquareLatticeCalculator.h"

SquareLatticeCalculator::SquareLatticeCalculator(double scale) :
    rcut(2.5),
    scal(scale),
    burgers(scale),
    nb_atoms(10),
    normalisation(std::pow(burgers, 2.0)), // Unit cell area is a²
    square_basis({{
        {0.0, 0.0}
    }}) {}

// Calculate energy
double SquareLatticeCalculator::calculate_energy(const Eigen::Matrix2d& C,
                                                const std::function<double(double)>& pot,
                                                double zero) {
    double phi = 0.0;
    // In a square lattice, points are at integer positions
    for (int m = -nb_atoms; m <= nb_atoms; ++m) {
        for (int n = -nb_atoms; n <= nb_atoms; ++n) {
            // Convert from lattice coordinates to Cartesian coordinates
            double px = m;
            double py = n;
            
            // Skip the origin
            if (std::abs(px) < 1e-10 && std::abs(py) < 1e-10) continue;
            
            Eigen::Vector2d p(px, py);
            // Calculate distance with scaling
            double r_squared = p.transpose() * C * p;
            double r = scal * std::sqrt(r_squared);
            
            if (r <= rcut) {
                phi += 0.5 * pot(r);
            }
        }
    }
    return phi - zero;
}

// Calculate energy derivative
Eigen::Matrix2d SquareLatticeCalculator::calculate_derivative(const Eigen::Matrix2d& C,
                                                            const std::function<double(double)>& dpot) {
    Eigen::Matrix2d phi = Eigen::Matrix2d::Zero();
    
    for (int m = -nb_atoms; m <= nb_atoms; ++m) {
        for (int n = -nb_atoms; n <= nb_atoms; ++n) {
            // Convert from lattice coordinates to Cartesian coordinates
            double px = m;
            double py = n;
            
            // Skip the origin
            if (std::abs(px) < 1e-10 && std::abs(py) < 1e-10) continue;
            
            Eigen::Vector2d p(px, py);
            // Calculate distance with scaling
            double r_squared = p.transpose() * C * p;
            double r = scal * std::sqrt(r_squared);
            
            if (r <= rcut) {
                double tmp_d = dpot(r);
                double dphi = 0.5 * tmp_d / r;
                phi(0, 0) += 0.5 * dphi * px * px * scal * scal;
                phi(1, 1) += 0.5 * dphi * py * py * scal * scal;
                phi(0, 1) += dphi * px * py * scal * scal;
                phi(1, 0) += dphi * px * py * scal * scal;
            }
        }
    }
    
    // Apply symmetry factors
    phi(0, 1) *= 0.5;
    phi(1, 0) *= 0.5;
    
    return phi;
}

// Helper method to get nearest neighbor distance in perfect lattice
double SquareLatticeCalculator::getNearestNeighborDistance() const {
    return scal;
}

// Helper method to get area of unit cell
double SquareLatticeCalculator::getUnitCellArea() const {
    return normalisation;
}