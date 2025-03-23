#ifndef CHANGE_MEASURES_H
#define CHANGE_MEASURES_H

#include <vector>
#include <Eigen/Dense>
#include <stdexcept>
#include <cmath>
#include "../include/geometry/Point2D.h"
#include "src/optimization.h" // This path is relative to ALGLIB_DIR
#include "src/ap.h" // This path is relative to ALGLIB_DIR


// Structure to store various measures of change between point sets
// Structure to store various measures of change between point sets
struct ChangeMeasures {
    double euclidean_norm;        // L2 norm of the differences
    double max_abs_change;        // Maximum absolute change
    double relative_change;       // Average squared relative change
    double normalized_euclidean_norm; // Euclidean norm divided by sqrt(n)
    double mean_abs_change;       // Mean absolute change
};

/**
 * Computes various measures of change between two alglib arrays
 * 
 * @param current_points The current (or new) array of points
 * @param reference_points The reference (or old) array of points to compare against
 * @return A ChangeMeasures structure containing the computed metrics
 * @throws std::invalid_argument if the arrays have different sizes or are empty
 */
ChangeMeasures computeChangeMeasures(const alglib::real_1d_array& current_points,
                                     const alglib::real_1d_array& reference_points);

#endif // CHANGE_MEASURES_H