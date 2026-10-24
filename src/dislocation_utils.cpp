#include "../include/utils/dislocation_utils.h"
#include "../include/lattice_energy/EnergyFunctions.h"
#include "../include/lattice_energy/TriangularLatticeCalculator.h"
#include "../include/lattice_energy/SquareLatticeCalculator.h"
#include "../include/interatomic/inter_atomic.h"
#include <limits>
#include <iostream>
#include <iomanip>

DislocationSolutionParams DislocationSolutionParams::fromElasticConstants(
    double c11, double c12, double c66,
    const std::string& lattice
) {
    DislocationSolutionParams p;
    p.lattice_type = lattice;
    p.C11 = c11;
    p.C22 = c11;
    p.C12 = c12;
    p.C66 = c66;

    double diff = c11 - c12;
    p.Zener_A = (std::abs(diff) > 1e-12) ? (2.0 * c66 / diff) : 1.0;
    p.born_stable = (c11 > 0.0) && (c11 * c11 - c12 * c12 > 0.0) && (c66 > 0.0);

    // 2D Poisson ratio (Baggio, Patinet, Salman, Truskinovsky, Sec. 4.2 & Eq. 15):
    // nu = C12 / C11 (equals 1/3 when Cauchy relation C11 = 3*C12 holds)
    p.poisson_ratio = (std::abs(c11) > 1e-12) ? (c12 / c11) : (1.0 / 3.0);

    // Stroh biquadratic parameter kappa (Eq. A5):
    // kappa = (C11^2 - C12^2 - 2*C12*C66) / (2*C11*C66)
    double denom = 2.0 * c11 * c66;
    if (std::abs(denom) > 1e-12) {
        p.kappa = (c11 * c11 - c12 * c12 - 2.0 * c12 * c66) / denom;
    } else {
        p.kappa = 1.0;
    }

    // Characteristic roots p_{1,2} (Eq. A8) and Stroh coefficients (Eq. A16)
    if (p.kappa > -1.0 && p.kappa < 1.0 - 1e-10) {
        // Regime A_Z > 1 (-1 < kappa < 1): complex conjugate roots in upper half-plane p_{1,2} = \mp xi + i*zeta
        p.xi   = std::sqrt(0.5 * (1.0 - p.kappa));
        p.zeta = std::sqrt(0.5 * (1.0 + p.kappa));
        p.g1   = (c11 * c11 - c12 * c12) / (16.0 * p.xi * p.zeta * c11 * c66);
        p.g2   = (c11 + c12) / (8.0 * p.xi * c11);
        p.g3   = (c11 - c12) / (8.0 * p.zeta * c11);
    } else if (p.kappa > 1.0 + 1e-10) {
        // Regime 0 < A_Z < 1 (kappa > 1): purely imaginary roots p_{1,2} = i*(zeta \mp xi_imag)
        p.xi   = std::sqrt(0.5 * (p.kappa - 1.0));
        p.zeta = std::sqrt(0.5 * (p.kappa + 1.0));
        p.g1   = (c11 * c11 - c12 * c12) / (16.0 * p.xi * p.zeta * c11 * c66);
        p.g2   = (c11 + c12) / (8.0 * p.xi * c11);
        p.g3   = (c11 - c12) / (8.0 * p.zeta * c11);
    } else {
        // Degenerate isotropic limit (A_Z -> 1, kappa -> 1)
        p.xi   = 0.0;
        p.zeta = 1.0;
        p.g1   = 0.0;
        p.g2   = 0.0;
        p.g3   = (std::abs(c11) > 1e-12) ? ((c11 - c12) / (8.0 * c11)) : 0.0;
    }

    // Foreman's prelogarithmic energy factor K (Eq. A17), which reduces in the
    // isotropic triangular case (2*C66 = C11 - C12) to mu*(1 + nu) in Eq. (18)
    double foreman_denom = c11 * (c11 + c12 + 2.0 * c66);
    double foreman_num   = c66 * (c11 - c12);
    if (foreman_denom > 0.0 && foreman_num >= 0.0) {
        p.K_foreman = (c11 + c12) * std::sqrt(foreman_num / foreman_denom);
    } else {
        p.K_foreman = 0.0;
    }

    return p;
}

DislocationSolutionParams DislocationSolutionParams::fromZanzottoPotential(
    double beta, double K,
    const std::string& lattice
) {
    auto ec = energy_functions::calculate_elastic_constants(beta, K, lattice);
    return fromElasticConstants(ec.C11, ec.C12, ec.C66, lattice);
}

DislocationSolutionParams DislocationSolutionParams::fromCurrentZanzottoPotential(
    const std::string& lattice
) {
    return fromZanzottoPotential(
        energy_functions::current_beta,
        energy_functions::current_K,
        lattice
    );
}

DislocationSolutionParams DislocationSolutionParams::fromAtomisticPotential(
    BaseLatticeCalculator& calculator,
    const std::function<double(double)>& dpot,
    const std::function<double(double)>& d2pot,
    const std::string& lattice
) {
    Eigen::Matrix2d C_ref = Eigen::Matrix2d::Identity();
    double norm = calculator.getUnitCellArea();
    HessianComponents hess = calculator.calculate_dseconderivative_components(C_ref, dpot, d2pot);
    double c11 = 4.0 * hess.c11_c11 / norm;
    double c12 = 4.0 * hess.c11_c22 / norm;
    double c66 = 1.0 * hess.c12_c12 / norm;
    return fromElasticConstants(c11, c12, c66, lattice);
}

DislocationSolutionParams DislocationSolutionParams::fromPotentialMode(
    int potential_mode,
    const std::string& lattice
) {
    if (potential_mode == 3) {
        return fromCurrentZanzottoPotential(lattice);
    }

    double r_cutoff = 2.5;
    std::function<double(double)> pot;
    std::function<double(double)> dpot;
    std::function<double(double)> d2pot;

    if (potential_mode == 1) {
        r_cutoff = 2.5;
        pot = lennard_jones_energy_v2;
        dpot = lennard_jones_energy_der_v2;
        d2pot = lennard_jones_energy_sder_v2;
    } else if (potential_mode == 2) {
        r_cutoff = 1.86602540378444;
        pot = lennard_jones_energy_v3;
        dpot = lennard_jones_energy_der_v3;
        d2pot = lennard_jones_energy_sder_v3;
    } else if (potential_mode == 4) {
        r_cutoff = 2.5;
        pot = square_energy;
        dpot = square_energy_der;
        d2pot = square_energy_sder;
    } else {
        std::cerr << "Unknown potential mode: " << potential_mode << std::endl;
        exit(EXIT_FAILURE);
    }

    double scale = find_optimal_lattice_parameter(pot, lattice, r_cutoff, dpot, d2pot);

    if (lattice == "triangular") {
        TriangularLatticeCalculator calc(scale, r_cutoff);
        return fromAtomisticPotential(calc, dpot, d2pot, lattice);
    } else {
        SquareLatticeCalculator calc(scale, r_cutoff);
        return fromAtomisticPotential(calc, dpot, d2pot, lattice);
    }
}

void DislocationSolutionParams::printSummary(double burgers_magnitude) const {
    double b2_over_4pi = (burgers_magnitude * burgers_magnitude) / (4.0 * M_PI);
    std::cout << "  --- Dislocation Solution Parameters (" << lattice_type << " lattice) ---\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  Elastic Moduli: C11 = " << C11
              << ", C12 = " << C12
              << ", C66 = " << C66
              << " | Zener A_Z = " << Zener_A << "\n";
    if (lattice_type == "triangular") {
        std::cout << "  Formulation:    Isotropic 2D Volterra (Eq. 15 of Baggio et al.)\n";
        std::cout << "  Parameters:     nu = C12/C11 = " << poisson_ratio
                  << " | (1+nu)/4 = " << (0.25 * (1.0 + poisson_ratio))
                  << " | (1-nu)/4 = " << (0.25 * (1.0 - poisson_ratio)) << "\n";
    } else {
        std::cout << "  Formulation:    Anisotropic Stroh / Eshelby-Read-Shockley (Eqs. A5-A16)\n";
        std::cout << "  Stroh Roots:    kappa = " << kappa
                  << ", xi = " << xi
                  << ", zeta = " << zeta << "\n";
        std::cout << "  Coefficients:   g1 = " << g1
                  << ", g2 = " << g2
                  << ", g3 = " << g3 << "\n";
    }
    std::cout << "  Energy Factor:  K = " << K_foreman
              << " | K*b^2/(4*pi) = " << (K_foreman * b2_over_4pi)
              << " (|b| = " << burgers_magnitude << ")\n";
}

// Nabarro smoothing factor implementation
double NabarroSmoothing::smoothingFactor(
    double distance,
    double core_radius
) {
    if (core_radius <= 0) {
        std::cerr << "Warning: Invalid core radius. Using default value." << std::endl;
        core_radius = 1.0;
    }
    if (distance < 0) {
        std::cerr << "Warning: Negative distance. Using absolute value." << std::endl;
        distance = std::abs(distance);
    }
    if (distance <= 1e-10) {
        return 1.0;
    }
    if (distance <= core_radius) {
        return 1.0;
    } else {
        return core_radius / distance;
    }
}

// Isotropic 2D Volterra displacement calculation (Eq. 15 of Baggio, Patinet, Salman, Truskinovsky):
//   u1 = (b / 2pi) * [ atan2(X2, X1) + (1 + nu)*X1*X2 / (2*(X1^2 + X2^2)) ]
//   u2 = -(b / 2pi) * [ (1 - nu)*ln(X1^2 + X2^2)/4 + (1 + nu)*(X1^2 - X2^2) / (4*(X1^2 + X2^2)) ]
// with 2D Poisson ratio nu = C12 / C11 (e.g. nu = 1/3 for LJ triangular lattice)
Eigen::Vector2d VolterraDisplacement::calculateEdgeDisplacement(
    const Eigen::Vector2d& position,
    const Eigen::Vector2d& burgers_vector,
    double poisson_ratio
) {
    double xx = position.x();
    double yy = position.y();
    double r2 = xx * xx + yy * yy;

    if (r2 < 1e-12) {
        return Eigen::Vector2d::Zero();
    }

    double b = burgers_vector.x();
    double nu = poisson_ratio;

    double ux = (b / (2.0 * M_PI)) * (
        std::atan2(yy, xx) + (1.0 + nu) * xx * yy / (2.0 * r2)
    );
    double uy = -(b / (2.0 * M_PI)) * (
        0.25 * (1.0 - nu) * std::log(r2) +
        (1.0 + nu) * (xx * xx - yy * yy) / (4.0 * r2)
    );

    return Eigen::Vector2d(ux, uy);
}

Eigen::Vector2d VolterraDisplacement::calculateEdgeDisplacement(
    const Eigen::Vector2d& position,
    const Eigen::Vector2d& burgers_vector,
    const DislocationSolutionParams& params
) {
    return calculateEdgeDisplacement(position, burgers_vector, params.poisson_ratio);
}

// Anisotropic Stroh displacement calculation (for square lattice, Eqs. A14-A16)
Eigen::Vector2d VolterraDisplacement::calculateAnisotropicEdgeDisplacement(
    const Eigen::Vector2d& position,
    const Eigen::Vector2d& burgers_vector,
    const DislocationSolutionParams& params
) {
    double x1 = position.x();
    double x2 = position.y();
    double b = burgers_vector.x();

    double r2 = x1 * x1 + x2 * x2;
    if (r2 < 1e-12) {
        return Eigen::Vector2d::Zero();
    }

    // If the crystal is at the isotropic limit (xi -> 0), use the smooth limit Eq. (15)
    if (params.xi < 1e-7) {
        return calculateEdgeDisplacement(position, burgers_vector, params.poisson_ratio);
    }

    // Regime 0 < A_Z < 1 (kappa > 1): purely imaginary Stroh roots p_{1,2} = i*beta_{1,2}
    if (params.kappa > 1.0 + 1e-10) {
        double beta1 = params.zeta - params.xi;
        double beta2 = params.zeta + params.xi;
        double r1_sq = x1 * x1 + (beta1 * x2) * (beta1 * x2);
        double r2_sq = x1 * x1 + (beta2 * x2) * (beta2 * x2);
        double th1   = std::atan2(beta1 * x2, x1);
        double th2   = std::atan2(beta2 * x2, x1);
        double gamma = (params.C11 * params.C11 - params.C12 * params.C12) /
                       (2.0 * params.C11 * params.C66);
        double nu    = params.poisson_ratio;

        double u1 = (b / (2.0 * M_PI * (beta2 - beta1))) * (
            (beta2 - gamma / beta2) * th1 - (beta1 - gamma / beta1) * th2
        );
        double u2 = -(b / (4.0 * M_PI * (beta2 - beta1))) * (
            ((beta1 * beta1 + nu) / beta1) * std::log(r1_sq) -
            ((beta2 * beta2 + nu) / beta2) * std::log(r2_sq)
        );
        return Eigen::Vector2d(u1, u2);
    }

    // Regime A_Z > 1 (-1 < kappa < 1): complex Stroh roots p_{1,2} = \mp xi + i*zeta (Eqs. A14-A15)
    double dx1 = x1 - params.xi * x2;
    double dx2 = x1 + params.xi * x2;
    double by  = params.zeta * x2;

    double r1_sq = dx1 * dx1 + by * by;
    double r2_sq = dx2 * dx2 + by * by;

    double th1 = std::atan2(by, dx1);
    double th2 = std::atan2(by, dx2);

    double u1 = (b / M_PI) * (0.25 * (th1 + th2) - 0.5 * params.g1 * std::log(r1_sq / r2_sq));
    double u2 = (b / M_PI) * (params.g2 * (th1 - th2) - 0.5 * params.g3 * std::log(r1_sq * r2_sq));

    return Eigen::Vector2d(u1, u2);
}

Eigen::Vector2d VolterraDisplacement::calculateAnisotropicEdgeDisplacement(
    const Eigen::Vector2d& position,
    const Eigen::Vector2d& burgers_vector,
    double C11,
    double C12,
    double C66
) {
    DislocationSolutionParams params = DislocationSolutionParams::fromElasticConstants(C11, C12, C66, "square");
    return calculateAnisotropicEdgeDisplacement(position, burgers_vector, params);
}

Eigen::Vector2d VolterraDisplacement::calculateAnisotropicEdgeDisplacement(
    const Eigen::Vector2d& position,
    const Eigen::Vector2d& burgers_vector
) {
    // Dynamically compute Stroh parameters from the active Zanzotto potential
    DislocationSolutionParams params = DislocationSolutionParams::fromCurrentZanzottoPotential("square");
    return calculateAnisotropicEdgeDisplacement(position, burgers_vector, params);
}

Eigen::Vector2d VolterraDisplacement::calculateDisplacement(
    const Eigen::Vector2d& position,
    const Eigen::Vector2d& burgers_vector,
    const DislocationSolutionParams& params
) {
    if (params.lattice_type == "triangular") {
        return calculateEdgeDisplacement(position, burgers_vector, params);
    }
    return calculateAnisotropicEdgeDisplacement(position, burgers_vector, params);
}

// Non-singular dislocation displacement calculation implementation (regularized Eq. 15)
Eigen::Vector2d NonSingularDisplacement::calculateEdgeDisplacement(
    const Eigen::Vector2d& position,
    const Eigen::Vector2d& burgers_vector,
    double core_radius,
    double poisson_ratio
) {
    double xx = position.x();
    double yy = position.y();
    double r2 = xx * xx + yy * yy;
    double r2_smoothed = r2 + core_radius * core_radius;

    double b = burgers_vector.x();
    double nu = poisson_ratio;

    double ux = (b / (2.0 * M_PI)) * (
        std::atan2(yy, xx) + (1.0 + nu) * xx * yy / (2.0 * r2_smoothed)
    );
    double uy = -(b / (2.0 * M_PI)) * (
        0.25 * (1.0 - nu) * std::log(r2_smoothed) +
        (1.0 + nu) * (xx * xx - yy * yy) / (4.0 * r2_smoothed)
    );

    return Eigen::Vector2d(ux, uy);
}

Eigen::Vector2d NonSingularDisplacement::calculateEdgeDisplacement(
    const Eigen::Vector2d& position,
    const Eigen::Vector2d& burgers_vector,
    double core_radius,
    const DislocationSolutionParams& params
) {
    if (params.lattice_type == "triangular" || params.xi < 1e-7) {
        return calculateEdgeDisplacement(position, burgers_vector, core_radius, params.poisson_ratio);
    }

    double x1 = position.x();
    double x2 = position.y();
    double b = burgers_vector.x();
    double r0_sq = core_radius * core_radius;

    if (params.kappa > 1.0 + 1e-10) {
        double beta1 = params.zeta - params.xi;
        double beta2 = params.zeta + params.xi;
        double r1_sq = x1 * x1 + (beta1 * x2) * (beta1 * x2) + r0_sq;
        double r2_sq = x1 * x1 + (beta2 * x2) * (beta2 * x2) + r0_sq;
        double th1   = std::atan2(beta1 * x2, x1);
        double th2   = std::atan2(beta2 * x2, x1);
        double gamma = (params.C11 * params.C11 - params.C12 * params.C12) /
                       (2.0 * params.C11 * params.C66);
        double nu    = params.poisson_ratio;

        double u1 = (b / (2.0 * M_PI * (beta2 - beta1))) * (
            (beta2 - gamma / beta2) * th1 - (beta1 - gamma / beta1) * th2
        );
        double u2 = -(b / (4.0 * M_PI * (beta2 - beta1))) * (
            ((beta1 * beta1 + nu) / beta1) * std::log(r1_sq) -
            ((beta2 * beta2 + nu) / beta2) * std::log(r2_sq)
        );
        return Eigen::Vector2d(u1, u2);
    }

    double dx1 = x1 - params.xi * x2;
    double dx2 = x1 + params.xi * x2;
    double by  = params.zeta * x2;

    double r1_sq = dx1 * dx1 + by * by + r0_sq;
    double r2_sq = dx2 * dx2 + by * by + r0_sq;

    double th1 = std::atan2(by, dx1);
    double th2 = std::atan2(by, dx2);

    double u1 = (b / M_PI) * (0.25 * (th1 + th2) - 0.5 * params.g1 * std::log(r1_sq / r2_sq));
    double u2 = (b / M_PI) * (params.g2 * (th1 - th2) - 0.5 * params.g3 * std::log(r1_sq * r2_sq));

    return Eigen::Vector2d(u1, u2);
}

// Create single dislocation implementation
std::vector<Point2D> createSingleDislocation(
    const std::vector<Point2D>& original_points,
    const Eigen::Vector2d& burgers_vector,
    size_t middle_atom_index,
    double core_radius,
    double poisson_ratio
) {
    std::vector<Point2D> modified_points = original_points;
    if (middle_atom_index >= original_points.size()) {
        std::cerr << "Error: Middle atom index out of bounds" << std::endl;
        return modified_points;
    }

    const Point2D& dislocation_core = original_points[middle_atom_index];
    for (size_t i = 0; i < modified_points.size(); ++i) {
        if (i == middle_atom_index) continue;
        Eigen::Vector2d relative_pos = modified_points[i].coord - dislocation_core.coord;
        Eigen::Vector2d volterra_disp = NonSingularDisplacement::calculateEdgeDisplacement(
            relative_pos,
            burgers_vector,
            core_radius,
            poisson_ratio
        );
        modified_points[i].coord += volterra_disp;
    }
    return modified_points;
}

std::vector<Point2D> createSingleDislocation(
    const std::vector<Point2D>& original_points,
    const Eigen::Vector2d& burgers_vector,
    size_t middle_atom_index,
    double core_radius,
    const DislocationSolutionParams& params
) {
    std::vector<Point2D> modified_points = original_points;
    if (middle_atom_index >= original_points.size()) {
        std::cerr << "Error: Middle atom index out of bounds" << std::endl;
        return modified_points;
    }

    const Point2D& dislocation_core = original_points[middle_atom_index];
    for (size_t i = 0; i < modified_points.size(); ++i) {
        if (i == middle_atom_index) continue;
        Eigen::Vector2d relative_pos = modified_points[i].coord - dislocation_core.coord;
        Eigen::Vector2d disp = NonSingularDisplacement::calculateEdgeDisplacement(
            relative_pos,
            burgers_vector,
            core_radius,
            params
        );
        modified_points[i].coord += disp;
    }
    return modified_points;
}

// Create dislocation dipole implementation centered around a point
std::vector<Point2D> createDislocationDipole(
    const std::vector<Point2D>& original_points,
    const Eigen::Vector2d& burgers_vector,
    const Eigen::Vector2d& center_position,
    double separation_distance,
    double core_radius,
    double poisson_ratio,
    const Eigen::Vector2d& dipole_direction
) {
    std::vector<Point2D> modified_points = original_points;
    Eigen::Vector2d dir_normalized = dipole_direction.normalized();
    Eigen::Vector2d dislocation_core1 = center_position - 0.5 * separation_distance * dir_normalized;
    Eigen::Vector2d dislocation_core2 = center_position + 0.5 * separation_distance * dir_normalized;

    for (size_t i = 0; i < modified_points.size(); ++i) {
        Eigen::Vector2d relative_pos1 = modified_points[i].coord - dislocation_core1;
        Eigen::Vector2d relative_pos2 = modified_points[i].coord - dislocation_core2;

        Eigen::Vector2d disp1 = NonSingularDisplacement::calculateEdgeDisplacement(
            relative_pos1,
            burgers_vector,
            core_radius,
            poisson_ratio
        );
        Eigen::Vector2d disp2 = NonSingularDisplacement::calculateEdgeDisplacement(
            relative_pos2,
            -burgers_vector,
            core_radius,
            poisson_ratio
        );
        modified_points[i].coord += disp1 + disp2;
    }
    return modified_points;
}