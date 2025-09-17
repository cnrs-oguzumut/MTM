#ifndef DISLOCATION_UTILS_H
#define DISLOCATION_UTILS_H

#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include "../include/geometry/Point2D.h"



/**
 * @brief Nabarro smoothing functions for dislocation calculations
 * 
 * Provides smoothing factors to regularize the singular behavior
 * of dislocation displacement fields near the core region.
 */
class NabarroSmoothing {
public:
    /**
     * @brief Calculate smoothing factor for dislocation displacement
     * 
     * Applies smoothing to avoid singularities in the dislocation core
     * while preserving the far-field Volterra solution.
     * 
     * @param distance Distance from the dislocation core
     * @param core_radius Radius of the dislocation core region
     * @return Smoothing factor (1.0 in core, core_radius/distance outside)
     */
    static double smoothingFactor(double distance, double core_radius);
};

/**
 * @brief Classic Volterra dislocation displacement field calculations
 * 
 * Implements the singular elastic displacement field for edge dislocations
 * based on classical continuum mechanics theory.
 */
class VolterraDisplacement {
public:
    /**
     * @brief Calculate edge dislocation displacement using Volterra solution
     * 
     * Computes the displacement field for an edge dislocation using the
     * classical Volterra solution with logarithmic and arctangent terms.
     * 
     * @param position Position vector relative to dislocation core
     * @param burgers_vector Burgers vector of the dislocation
     * @param poisson_ratio Poisson's ratio of the material (default: 0.3)
     * @return Displacement vector (ux, uy)
     */
    static Eigen::Vector2d calculateEdgeDisplacement(
        const Eigen::Vector2d& position,
        const Eigen::Vector2d& burgers_vector,
        double poisson_ratio = 0.3
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
    /**
     * @brief Calculate regularized edge dislocation displacement
     * 
     * Computes displacement field with core radius smoothing to eliminate
     * the logarithmic singularity at the dislocation center.
     * 
     * @param position Position vector relative to dislocation core
     * @param burgers_vector Burgers vector of the dislocation
     * @param core_radius Core radius for regularization
     * @param poisson_ratio Poisson's ratio of the material (default: 0.3)
     * @return Smoothed displacement vector (ux, uy)
     */
    static Eigen::Vector2d calculateEdgeDisplacement(
        const Eigen::Vector2d& position,
        const Eigen::Vector2d& burgers_vector,
        double core_radius,
        double poisson_ratio = 0.3
    );
};

/**
 * @brief Create a single edge dislocation in an atomic configuration
 * 
 * Applies dislocation displacement field to all atoms in the system
 * based on their distance from a specified dislocation core atom.
 * 
 * @param original_points Vector of initial atomic positions
 * @param burgers_vector Burgers vector defining dislocation strength and direction
 * @param middle_atom_index Index of atom to serve as dislocation core
 * @param core_radius Core radius for displacement smoothing
 * @param poisson_ratio Material Poisson's ratio (default: 0.3)
 * @return Vector of modified atomic positions with dislocation
 */
std::vector<Point2D> createSingleDislocation(
    const std::vector<Point2D>& original_points,
    const Eigen::Vector2d& burgers_vector,
    size_t middle_atom_index,
    double core_radius,
    double poisson_ratio = 0.3
);

/**
 * @brief Create a dislocation dipole in an atomic configuration
 * 
 * Creates two dislocations with opposite Burgers vectors separated by
 * a specified distance, resulting in a dipole configuration with
 * localized strain fields.
 * 
 * @param original_points Vector of initial atomic positions
 * @param burgers_vector Burgers vector for the positive dislocation
 * @param center_position Center point between the two dislocations
 * @param separation_distance Total distance between dislocation cores
 * @param core_radius Core radius for displacement smoothing
 * @param poisson_ratio Material Poisson's ratio (default: 0.3)
 * @param dipole_direction Direction vector for dipole orientation (default: x-direction)
 * @return Vector of modified atomic positions with dislocation dipole
 */
std::vector<Point2D> createDislocationDipole(
    const std::vector<Point2D>& original_points,
    const Eigen::Vector2d& burgers_vector,
    const Eigen::Vector2d& center_position,
    double separation_distance,
    double core_radius,
    double poisson_ratio = 0.3,
    const Eigen::Vector2d& dipole_direction = Eigen::Vector2d(1.0, 0.0)
);

#endif // DISLOCATION_UTILS_H