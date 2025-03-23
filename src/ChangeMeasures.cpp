#include "../include/output/ChangeMeasures.h"

ChangeMeasures computeChangeMeasures(const alglib::real_1d_array& current_points,
                                     const alglib::real_1d_array& reference_points) {
    int n = current_points.length();
    
    if (n != reference_points.length() || n == 0) {
        throw std::invalid_argument("Arrays must be of the same nonzero length");
    }
    
    double sum_squared_diff = 0.0;
    double max_abs_diff = 0.0;
    double sum_relative_change = 0.0;
    double sum_abs_diff = 0.0;
    
    for (int i = 0; i < n; ++i) {
        double diff = current_points[i] - reference_points[i];
        sum_squared_diff += diff * diff;
        max_abs_diff = std::max(max_abs_diff, std::abs(diff));
        
        if (std::abs(reference_points[i]) > 1e-12) {  // Avoid division by zero
            sum_relative_change += (diff * diff) / (reference_points[i] * reference_points[i]);
        }
        
        sum_abs_diff += std::abs(diff);
    }
    
    ChangeMeasures measures;
    measures.euclidean_norm = std::sqrt(sum_squared_diff);
    measures.max_abs_change = max_abs_diff;
    measures.relative_change = sum_relative_change / n;
    measures.normalized_euclidean_norm = measures.euclidean_norm / std::sqrt(n);
    measures.mean_abs_change = sum_abs_diff / n;
    
    return measures;
}