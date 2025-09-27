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

// Calculate second derivatives (Hessian) of energy w.r.t. metric tensor components
HessianComponents SquareLatticeCalculator::calculate_hessian_components(const Eigen::Matrix2d& C,
                                                                       const std::function<double(double)>& pot,
                                                                       const std::function<double(double)>& dpot,
                                                                       const std::function<double(double)>& d2pot) {
    HessianComponents hessian;
    
    // Loop over lattice points (similar to the original C code)
    for (int m = -nb_atoms; m <= nb_atoms; ++m) {
        for (int n = -nb_atoms; n <= nb_atoms; ++n) {
            // Skip the origin
            if (m == 0 && n == 0) continue;
            
            // Convert from lattice coordinates to Cartesian coordinates
            double dx2 = static_cast<double>(m);
            double dy2 = static_cast<double>(n);
            
            // Calculate squared distance using metric tensor
            double r2 = dx2 * dx2 * C(0,0) + dy2 * dy2 * C(1,1) + 2.0 * dx2 * dy2 * C(0,1);
            
            // Apply cutoff
            if (r2 >= r_cutoff_sq) continue;
            
            double r = std::sqrt(r2) * scal;
            
            // Calculate potential derivatives
            double tempf = 0.5 * dpot(r) / r;        // First derivative factor
            double tempf2 = 0.5 * d2pot(r) / r;     // Second derivative factor
            
            // Calculate geometric factors (same as A, B, C, D, E in original code)
            double A = dx2 * dx2 * dx2 * dx2;        // dx⁴
            double B = dy2 * dy2 * dy2 * dy2;        // dy⁴  
            double C_geom = dx2 * dx2 * dy2 * dy2;   // dx²dy²
            double D = dx2 * dx2 * dx2 * dy2;        // dx³dy
            double E = dy2 * dy2 * dy2 * dx2;        // dy³dx
            
            // Calculate Hessian components (following the original formulas)
            hessian.c11_c11 += 0.25 * tempf2 * A / r - 0.25 * tempf * A / (r * r);
            hessian.c22_c22 += 0.25 * tempf2 * B / r - 0.25 * tempf * B / (r * r);
            hessian.c12_c12 += tempf2 * C_geom / r - tempf * C_geom / (r * r);
            
            hessian.c11_c12 += 0.5 * tempf2 * D / r - 0.5 * tempf * D / (r * r);
            hessian.c22_c12 += 0.5 * tempf2 * E / r - 0.5 * tempf * E / (r * r);
            hessian.c11_c22 += 0.25 * tempf2 * C_geom / r - 0.25 * tempf * C_geom / (r * r);
        }
    }
    
    return hessian;
}

// Calculate second derivatives as 4th order ITensor
itensor::ITensor SquareLatticeCalculator::calculate_dseconderivative(const Eigen::Matrix2d& C,
                                                                    const std::function<double(double)>& dpot,
                                                                    const std::function<double(double)>& d2pot) {
    // For this implementation, we'll use a dummy pot function since we only need dpot and d2pot
    auto dummy_pot = [](double r) { return 0.0; };
    
    HessianComponents hess = calculate_hessian_components(C, dummy_pot, dpot, d2pot);
    return hessianComponentsToITensor(hess);
}

// Convert HessianComponents to 4th order ITensor
itensor::ITensor SquareLatticeCalculator::hessianComponentsToITensor(const HessianComponents& hess) const {
    // Create indices for the 4th order tensor (i,j,k,l) where ∂²E/∂C_ij∂C_kl
    auto i = itensor::Index(2, "i");
    auto j = itensor::Index(2, "j");
    auto k = itensor::Index(2, "k");
    auto l = itensor::Index(2, "l");
    
    itensor::ITensor tensor(i, j, k, l);
    
    // Set tensor components based on Hessian components
    // Note: ITensor uses 1-based indexing
    tensor.set(i=1, j=1, k=1, l=1, hess.c11_c11);  // ∂²E/∂c₁₁²
    tensor.set(i=2, j=2, k=2, l=2, hess.c22_c22);  // ∂²E/∂c₂₂²
    tensor.set(i=1, j=2, k=1, l=2, hess.c12_c12);  // ∂²E/∂c₁₂²
    tensor.set(i=2, j=1, k=2, l=1, hess.c12_c12);  // ∂²E/∂c₂₁² (symmetry)
    
    // Mixed derivatives
    tensor.set(i=1, j=1, k=2, l=2, hess.c11_c22);  // ∂²E/∂c₁₁∂c₂₂
    tensor.set(i=2, j=2, k=1, l=1, hess.c11_c22);  // ∂²E/∂c₂₂∂c₁₁ (symmetry)
    
    tensor.set(i=1, j=1, k=1, l=2, hess.c11_c12);  // ∂²E/∂c₁₁∂c₁₂
    tensor.set(i=1, j=1, k=2, l=1, hess.c11_c12);  // ∂²E/∂c₁₁∂c₂₁ (symmetry)
    tensor.set(i=1, j=2, k=1, l=1, hess.c11_c12);  // ∂²E/∂c₁₂∂c₁₁ (symmetry)
    tensor.set(i=2, j=1, k=1, l=1, hess.c11_c12);  // ∂²E/∂c₂₁∂c₁₁ (symmetry)
    
    tensor.set(i=2, j=2, k=1, l=2, hess.c22_c12);  // ∂²E/∂c₂₂∂c₁₂
    tensor.set(i=2, j=2, k=2, l=1, hess.c22_c12);  // ∂²E/∂c₂₂∂c₂₁ (symmetry)
    tensor.set(i=1, j=2, k=2, l=2, hess.c22_c12);  // ∂²E/∂c₁₂∂c₂₂ (symmetry)
    tensor.set(i=2, j=1, k=2, l=2, hess.c22_c12);  // ∂²E/∂c₂₁∂c₂₂ (symmetry)
    
    return tensor;
}

// Calculate acoustic tensor components as 2x2 blocks
std::vector<Eigen::Matrix2d> SquareLatticeCalculator::calculateAcousticTensorBlocks(const Eigen::Matrix2d& C,
                                                                                   const std::function<double(double)>& pot,
                                                                                   const std::function<double(double)>& dpot,
                                                                                   const std::function<double(double)>& d2pot) {
    HessianComponents hess = calculate_hessian_components(C, pot, dpot, d2pot);
    
    std::vector<Eigen::Matrix2d> blocks(4);
    
    // Block (0,0): ∂²E/∂c₁ᵢ∂c₁ⱼ
    blocks[0] << hess.c11_c11, hess.c11_c12,
                 hess.c11_c12, hess.c12_c12;
    
    // Block (0,1): ∂²E/∂c₁ᵢ∂c₂ⱼ  
    blocks[1] << hess.c11_c12, hess.c11_c22,
                 hess.c12_c12, hess.c22_c12;
    
    // Block (1,0): ∂²E/∂c₂ᵢ∂c₁ⱼ (should be transpose of block (0,1))
    blocks[2] << hess.c11_c12, hess.c12_c12,
                 hess.c11_c22, hess.c22_c12;
    
    // Block (1,1): ∂²E/∂c₂ᵢ∂c₂ⱼ
    blocks[3] << hess.c12_c12, hess.c22_c12,
                 hess.c22_c12, hess.c22_c22;
    
    return blocks;
}

// Helper method to get nearest neighbor distance in perfect lattice
double SquareLatticeCalculator::getNearestNeighborDistance() const {
    return scal;
}

// Helper method to get area of unit cell
double SquareLatticeCalculator::getUnitCellArea() const {
    return normalisation;
}