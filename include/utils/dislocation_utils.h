#ifndef DISLOCATION_UTILS_H
#define DISLOCATION_UTILS_H

#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <string>
#include <functional>
#include "../include/geometry/Point2D.h"
#include "../include/lattice_energy/BaseLatticeCalculator.h"

/**
 * @brief Elastic constants and derived dislocation solution parameters
 *
 * Computes and stores the exact dislocation displacement field parameters
 * from the ground-state elastic moduli (C11, C12, C66) of any Zanzotto continuum
 * potential or atomistic pair potential:
 * - Triangular lattice: isotropic 2D Volterra solution (Eq. 15 of Baggio et al.)
 *   with 2D Poisson ratio nu = C12 / C11 (e.g. nu = 1/3 for Cauchy central-force lattices).
 * - Square lattice: anisotropic Stroh / Eshelby-Read-Shockley solution with
 *   characteristic roots p_{1,2} = \mp xi + i*zeta and coefficients (g1, g2, g3)
 *   (Eqs. A5-A16).
 */
struct DislocationSolutionParams {
    std::string lattice_type = "square"; // "square" (Stroh) or "triangular" (Volterra)
    double C11 = 0.0;
    double C22 = 0.0;
    double C12 = 0.0;
    double C66 = 0.0;
    double Zener_A = 1.0;

    // 2D Isotropic Volterra parameter (triangular lattice, Eq. 15 of Baggio et al.)
    double poisson_ratio = 1.0 / 3.0; // nu = C12 / C11

    // Anisotropic Stroh parameters (square lattice, Eqs. A5-A16)
    double kappa = 1.0;
    double xi    = 0.0; // real part of Stroh root (alpha in legacy notation)
    double zeta  = 1.0; // imag part of Stroh root (beta in legacy notation)
    double g1    = 0.0; // B1 in legacy notation
    double g2    = 0.0; // C2 in legacy notation
    double g3    = 0.0; // D2 in legacy notation

    // Foreman prelogarithmic energy factor K (Eq. A17)
    double K_foreman = 0.0;
    bool born_stable = true;

    /**
     * @brief Construct dislocation solution parameters from elastic moduli (C11, C12, C66)
     */
    static DislocationSolutionParams fromElasticConstants(
        double c11, double c12, double c66,
        const std::string& lattice = "square"
    );

    /**
     * @brief Compute elastic constants and dislocation parameters from Zanzotto potential (beta, K)
     */
    static DislocationSolutionParams fromZanzottoPotential(
        double beta, double K,
        const std::string& lattice = "square"
    );

    /**
     * @brief Compute elastic constants and dislocation parameters from active Zanzotto parameters
     */
    static DislocationSolutionParams fromCurrentZanzottoPotential(
        const std::string& lattice = "square"
    );

    /**
     * @brief Compute elastic constants and dislocation parameters from an atomistic pair potential
     */
    static DislocationSolutionParams fromAtomisticPotential(
        BaseLatticeCalculator& calculator,
        const std::function<double(double)>& dpot,
        const std::function<double(double)>& d2pot,
        const std::string& lattice = "square"
    );

    /**
     * @brief Compute elastic constants and dislocation parameters from potential mode
     *        (mode 3 = Zanzotto continuum, mode 1 = Lennard-Jones v2, mode 2 = Lennard-Jones v3)
     */
    static DislocationSolutionParams fromPotentialMode(
        int potential_mode,
        const std::string& lattice = "square"
    );

    /**
     * @brief Print summary of elastic moduli and dislocation field parameters
     */
    void printSummary(double burgers_magnitude = 1.0) const;
};

/**
 * @brief Nabarro smoothing functions for dislocation calculations
 * 
 * Provides smoothing factors to regularize the singular behavior
 * of dislocation displacement fields near the core region.
 */
class NabarroSmoothing {
public:
    static double smoothingFactor(double distance, double core_radius);
};

/**
 * @brief Dislocation displacement field calculations:
 * - Isotropic Volterra field for triangular lattices
 * - Anisotropic Stroh field for square lattices
 */
class VolterraDisplacement {
public:
    /**
     * @brief Calculate edge dislocation displacement for triangular lattice (Volterra Eq. A23)
     * 
     * @param position Position vector relative to dislocation core
     * @param burgers_vector Burgers vector of the dislocation
     * @param poisson_ratio Plane-strain Poisson's ratio nu = C12 / (C11 + C12)
     * @return Displacement vector (ux, uy)
     */
    static Eigen::Vector2d calculateEdgeDisplacement(
        const Eigen::Vector2d& position,
        const Eigen::Vector2d& burgers_vector,
        double poisson_ratio
    );

    /**
     * @brief Calculate edge dislocation displacement for triangular lattice from elastic constants
     */
    static Eigen::Vector2d calculateEdgeDisplacement(
        const Eigen::Vector2d& position,
        const Eigen::Vector2d& burgers_vector,
        const DislocationSolutionParams& params
    );

    /**
     * @brief Calculate edge dislocation displacement in an anisotropic square crystal (Stroh Eqs. A14-A15)
     * 
     * Uses dynamically computed Stroh parameters from the active potential if none are passed.
     */
    static Eigen::Vector2d calculateAnisotropicEdgeDisplacement(
        const Eigen::Vector2d& position,
        const Eigen::Vector2d& burgers_vector
    );

    /**
     * @brief Calculate edge dislocation displacement in an anisotropic square crystal using given Stroh parameters
     */
    static Eigen::Vector2d calculateAnisotropicEdgeDisplacement(
        const Eigen::Vector2d& position,
        const Eigen::Vector2d& burgers_vector,
        const DislocationSolutionParams& params
    );

    /**
     * @brief Dispatch to Volterra (for triangular lattice) or Stroh (for square lattice)
     *        using dynamically computed parameters from the potential
     */
    static Eigen::Vector2d calculateDisplacement(
        const Eigen::Vector2d& position,
        const Eigen::Vector2d& burgers_vector,
        const DislocationSolutionParams& params
    );
};

/**
 * @brief Non-singular dislocation displacement field calculations
 * 
 * Implements regularized displacement fields that avoid singularities
 * by smoothing the core region while preserving far-field behavior.
 */
class NonSingularDisplacement {
public:
    static Eigen::Vector2d calculateEdgeDisplacement(
        const Eigen::Vector2d& position,
        const Eigen::Vector2d& burgers_vector,
        double core_radius,
        double poisson_ratio
    );

    static Eigen::Vector2d calculateEdgeDisplacement(
        const Eigen::Vector2d& position,
        const Eigen::Vector2d& burgers_vector,
        double core_radius,
        const DislocationSolutionParams& params
    );
};

std::vector<Point2D> createSingleDislocation(
    const std::vector<Point2D>& original_points,
    const Eigen::Vector2d& burgers_vector,
    size_t middle_atom_index,
    double core_radius,
    double poisson_ratio
);

std::vector<Point2D> createSingleDislocation(
    const std::vector<Point2D>& original_points,
    const Eigen::Vector2d& burgers_vector,
    size_t middle_atom_index,
    double core_radius,
    const DislocationSolutionParams& params
);

std::vector<Point2D> createDislocationDipole(
    const std::vector<Point2D>& original_points,
    const Eigen::Vector2d& burgers_vector,
    const Eigen::Vector2d& center_position,
    double separation_distance,
    double core_radius,
    double poisson_ratio,
    const Eigen::Vector2d& dipole_direction = Eigen::Vector2d(1.0, 0.0)
);

#endif // DISLOCATION_UTILS_H