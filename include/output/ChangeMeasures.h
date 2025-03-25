#ifndef CHANGE_MEASURES_H
#define CHANGE_MEASURES_H

#include <vector>
#include <Eigen/Dense>
#include <stdexcept>
#include <cmath>
#include "../include/geometry/Point2D.h"
#include "../include/reductions/LagrangeReduction.h"
#include "../include/mesh/Triangle.h"
#include "../include/mesh/ElementTriangle2D.h"
#include "../include/mesh/MeshGenerator.h"
#include "../include/lattice_energy/TriangularLatticeCalculator.h"
#include "../include/lattice_energy/SquareLatticeCalculator.h"
#include "../include/lattice_energy/Strain_Energy_LatticeCalculator.h"
#include "../include/optimization/LatticeOptimizer.h"
#include "src/optimization.h" // This path is relative to ALGLIB_DIR
#include "src/ap.h" // This path is relative to ALGLIB_DIR

// Structure to store various measures of change between point sets
struct ChangeMeasures {
    double euclidean_norm;
    double max_abs_change;
    double relative_change;
    double normalized_euclidean_norm;
    double mean_abs_change;
    bool displacement_exceeds_half_lattice;
    bool has_distorted_triangles;
};

/*
 * Computes various measures of change between two alglib arrays
 *
 * @param current_points The current (or new) array of points
 * @param reference_points The reference (or old) array of points to compare against
 * @param lattice_constant The constant representing the lattice spacing
 * @param elements Optional collection of finite elements to check for angle distortion
 * @param points Optional collection of points for angle calculations
 * @param check_angles Whether to perform triangle angle checks
 * @return A ChangeMeasures structure containing the computed metrics
 * @throws std::invalid_argument if the arrays have different sizes or are empty
 */
ChangeMeasures computeChangeMeasures(
    const alglib::real_1d_array& current_points,
    const alglib::real_1d_array& reference_points,
    double lattice_constant = 1.0,
    const std::vector<ElementTriangle2D>& elements = {},  // Keep as const
    const std::vector<Point2D>& points = {},
    bool check_angles = false,
    const Eigen::Matrix2d* F_ext = nullptr);
    
    
/*
 * Analyzes the reduction metrics for finite elements
 * 
 * @param elements The collection of finite elements to analyze
 * @param points The points that make up the vertices of the elements
 * @param userData User data containing active elements and other properties
 * @return The number of elements with reduction
 */
std::vector<int> analyzeElementReduction(
    std::vector<ElementTriangle2D>& elements,
    std::vector<Point2D>& points,
    const UserData* userData);


int compareM3Activation(const std::vector<int>& m3_before, const std::vector<int>& m3_after);


#endif // CHANGE_MEASURES_H