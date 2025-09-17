#include "../include/utils/dislocation_utils.h"
#include <limits>
#include <iostream>

// Nabarro smoothing factor implementation
double NabarroSmoothing::smoothingFactor(
    double distance,
    double core_radius
) {
    // Handle edge cases as before
    if (core_radius <= 0) {
        std::cerr << "Warning: Invalid core radius. Using default value." << std::endl;
        core_radius = 1.0;
    }
    
    if (distance < 0) {
        std::cerr << "Warning: Negative distance. Using absolute value." << std::endl;
        distance = std::abs(distance);
    }
    
    // At the very core, return a limited maximum value to avoid singularity
    if (distance <= 1e-10) {
        return 1.0;
    }
    
    // For distances beyond the core, use a more gradually decreasing function
    // This ensures the dislocation structure is preserved
    if (distance <= core_radius) {
        // Inside core: gradual transition from 1.0 to something less
        return 1.0;
    } else {
        // Outside core: preserve more of the Volterra solution
        // Use core_radius/distance for slower decay
        return core_radius / distance;
    }
}

// Volterra displacement calculation implementation
Eigen::Vector2d VolterraDisplacement::calculateEdgeDisplacement(
    const Eigen::Vector2d& position,
    const Eigen::Vector2d& burgers_vector,
    double poisson_ratio
) {
    // Extract coordinates
    double xx = position.x();
    double yy = position.y();
    
    // Calculate r^2 to avoid repeated calculation
    double r2 = xx*xx + yy*yy;
    
    // Prevent division by zero
    if (r2 < 1e-10) {
        return Eigen::Vector2d::Zero();
    }
    
    // Get Burgers vector magnitude (assuming x-component is the one we care about)
    double b = burgers_vector.x();
    double nu = poisson_ratio;
    
    // Calculate displacement using the provided formula
    double term1_x = std::atan2(yy, xx);
    double term2_x = (1 + nu) * xx * yy / (2.0 * r2);
    double ux = (b / (2.0 * M_PI)) * (term1_x + term2_x);
    
    double term1_y = (1 - nu) * 0.25 * std::log(r2);
    double term2_y = (1 + nu) * (xx * xx - yy * yy) / (4.0 * r2);
    double uy = -(b / (2.0 * M_PI)) * (term1_y + term2_y);
    
    return Eigen::Vector2d(ux, uy);
}

// Non-singular dislocation displacement calculation implementation
Eigen::Vector2d NonSingularDisplacement::calculateEdgeDisplacement(
    const Eigen::Vector2d& position,
    const Eigen::Vector2d& burgers_vector,
    double core_radius,
    double poisson_ratio
) {
    // Extract coordinates
    double xx = position.x();
    double yy = position.y();
    
    // Calculate r^2 to avoid repeated calculation
    double r2 = xx*xx + yy*yy;
    
    // Add core radius squared to smooth the singularity
    double r2_smoothed = r2 + core_radius*core_radius;
    
    // Get Burgers vector magnitude
    double b = burgers_vector.x();
    double nu = poisson_ratio;
    
    // Calculate smoothed displacement using modified formula
    double term1_x = std::atan2(yy, xx); // Arc tangent is already well-behaved
    double term2_x = (1 + nu) * xx * yy / (2.0 * r2_smoothed);
    double ux = (b / (2.0 * M_PI)) * (term1_x + term2_x);
    
    double term1_y = (1 - nu) * 0.25 * std::log(r2_smoothed);
    double term2_y = (1 + nu) * (xx * xx - yy * yy) / (4.0 * r2_smoothed);
    double uy = -(b / (2.0 * M_PI)) * (term1_y + term2_y);
    
    return Eigen::Vector2d(ux, uy);
}

// Create single dislocation implementation
std::vector<Point2D> createSingleDislocation(
    const std::vector<Point2D>& original_points,
    const Eigen::Vector2d& burgers_vector,
    size_t middle_atom_index,
    double core_radius,
    double poisson_ratio 
) {
    // Create a copy of the original points to modify
    std::vector<Point2D> modified_points = original_points;
    
    // Ensure the middle atom index is valid
    if (middle_atom_index >= original_points.size()) {
        std::cerr << "Error: Middle atom index out of bounds" << std::endl;
        return modified_points;
    }
    
    // Get the coordinates of the middle atom (dislocation core)
    const Point2D& dislocation_core = original_points[middle_atom_index];
    
    // Iterate through all points
    for (size_t i = 0; i < modified_points.size(); ++i) {
        // Skip the middle atom (dislocation core)
        if (i == middle_atom_index) continue;
        
        // Calculate relative position from dislocation core
        Eigen::Vector2d relative_pos = 
            modified_points[i].coord - dislocation_core.coord;
        
        // Calculate distance from dislocation core
        double distance = relative_pos.norm();
        
        // Calculate Volterra displacement
        Eigen::Vector2d volterra_disp = 
            NonSingularDisplacement::calculateEdgeDisplacement(
                relative_pos, 
                burgers_vector, 
                core_radius,
                poisson_ratio
            );
        
        // Calculate Nabarro smoothing factor
        // double smooth_factor = NabarroSmoothing::smoothingFactor(
        //     distance, 
        //     core_radius
        // );
        
        // Apply smoothed Volterra displacement
        modified_points[i].coord += volterra_disp * 1;
    }
    
    return modified_points;
}

// Create dislocation dipole implementation centered around a point
std::vector<Point2D> createDislocationDipole(
    const std::vector<Point2D>& original_points,
    const Eigen::Vector2d& burgers_vector,
    const Eigen::Vector2d& center_position,  // Center point between dipoles
    double separation_distance,              // Total distance between dipoles
    double core_radius,
    double poisson_ratio,
    const Eigen::Vector2d& dipole_direction  // Default x-direction
) {
    // Create a copy of the original points to modify
    std::vector<Point2D> modified_points = original_points;
    
    // Normalize the dipole direction vector
    Eigen::Vector2d dir_normalized = dipole_direction.normalized();
    
    // Calculate the two dislocation core positions
    Eigen::Vector2d dislocation_core1 = center_position - 0.5 * separation_distance * dir_normalized;
    Eigen::Vector2d dislocation_core2 = center_position + 0.5 * separation_distance * dir_normalized;
    
    // Find nearest atoms to core positions (optional)
    // size_t nearest_idx1 = findNearestAtom(original_points, dislocation_core1);
    // size_t nearest_idx2 = findNearestAtom(original_points, dislocation_core2);
    
    // Iterate through all points
    for (size_t i = 0; i < modified_points.size(); ++i) {
        // Calculate relative positions from both dislocation cores
        Eigen::Vector2d relative_pos1 = modified_points[i].coord - dislocation_core1;
        Eigen::Vector2d relative_pos2 = modified_points[i].coord - dislocation_core2;
        
        // Calculate displacements from both dislocations
        Eigen::Vector2d disp1 = NonSingularDisplacement::calculateEdgeDisplacement(
            relative_pos1, 
            burgers_vector, 
            core_radius,
            poisson_ratio
        );
        
        Eigen::Vector2d disp2 = NonSingularDisplacement::calculateEdgeDisplacement(
            relative_pos2, 
            -burgers_vector,  // opposite Burgers vector
            core_radius,
            poisson_ratio
        );
        
        // Apply the combined displacement
        modified_points[i].coord += disp1 + disp2;
    }
    
    return modified_points;
}

/* Optional helper function to find nearest atom to a position
size_t findNearestAtom(const std::vector<Point2D>& points, const Eigen::Vector2d& position) {
    size_t nearest_idx = 0;
    double min_dist = std::numeric_limits<double>::max();
    
    for (size_t i = 0; i < points.size(); ++i) {
        double dist = (points[i].coord - position).norm();
        if (dist < min_dist) {
            min_dist = dist;
            nearest_idx = i;
        }
    }
    return nearest_idx;
}
*/