#include <cmath>
#include <algorithm>
#include <array>
#include <map>
#include <iostream>
#include <vector>
#include <stdexcept>
#include "../include/output/ChangeMeasures.h"



// bool checkSquareDomainViolation(const std::vector<ElementTriangle2D>& elements) {
//     for (const auto& element : elements) {
//         if (!element.isInitialized()) continue;
        
//         const Eigen::Matrix2d& C = element.getMetricTensor();
//         double c11 = C(0, 0);
//         double c22 = C(1, 1); 
//         double c12 = C(0, 1);
        
//         // Check both conditions: C12 ≥ 0 and 2|C12| ≤ min(C11, C22)
//         if (c12 < 0.0 || 2.0 * std::abs(c12) > std::min(c11, c22)) {
            
//             return true;
//         }
//     }
//     return false;
// }

//checkTriangularDomainViolation
bool checkSquareDomainViolation(const std::vector<ElementTriangle2D>& elements) {
    
    for (const auto& element : elements) {
        if (!element.isInitialized()) continue;
        
        const Eigen::Matrix2d& C = element.getMetricTensor();
        double detC = C.determinant();
        
        // Early exit for non-physical deformation
        if (detC <= 0.0) {
            std::cout << "Violation: Non-physical deformation (det(C) <= 0)" << std::endl;
            return true;
        }
        
        // Normalize metric tensor components
        double inv_sqrt_detC = 1.0;  // Set to 1.0/sqrt(detC) if normalization needed
        double c11 = C(0, 0) * inv_sqrt_detC;
        double c22 = C(1, 1) * inv_sqrt_detC;
        double c12 = C(0, 1) * inv_sqrt_detC;
        
        double min_val = std::min(c11, c22);
        
        bool condition1 = (c12 > 0 && c12 < min_val);
        bool condition2 = (c12 < 0 && 2 * std::abs(c12) < min_val);

        
        // Violation if either condition is NOT satisfied
        bool violation_found = (!condition2 && !condition1 );
        
        // if( !condition1 &&  condition2 ){
        //     violation_found = false;
        // }

        // if( !condition1 && !condition2 ){
        //     violation_found = true;
        // }

        // if( condition1 &&  condition2 ){
        //     violation_found = false;
        // }
        
        // if( condition1 &&  !condition2 ){
        //     violation_found = false;
        // }






        if (violation_found) {
            // std::cout << "Violation found:" << std::endl;
            // std::cout << "  c11: " << c11 << ", c22: " << c22 << ", c12: " << c12 << std::endl;
            // std::cout << "  min(c11,c22): " << min_val << std::endl;
            // std::cout << "  Condition 1 (0 < c12 < min): " << (condition1 ? "PASS" : "FAIL") << std::endl;
            // std::cout << "  Condition 2 (c122|c12| < min): " << (condition2 ? "PASS" : "FAIL") << std::endl;
            return true;
        }
    }
    
    return false;  // No violation found in any element
}

ChangeMeasures computeChangeMeasures(const alglib::real_1d_array& current_points,
    const alglib::real_1d_array& reference_points,
    double lattice_constant,
    const std::vector<ElementTriangle2D>& elements,
    const UserData* userData,
    const std::vector<Point2D>& points,
    bool check_angles,
    const Eigen::Matrix2d* F_ext) {
    int n = current_points.length();
    if (n != reference_points.length() || n == 0) {
    throw std::invalid_argument("Arrays must be of the same nonzero length");
    }

    double sum_squared_diff = 0.0;
    double max_abs_diff = 0.0;
    double sum_relative_change = 0.0;
    double sum_abs_diff = 0.0;
    bool has_distorted_triangles = false;

    for (int i = 0; i < n; ++i) {

        double diff = current_points[i] - reference_points[i];
        sum_squared_diff += diff * diff;
        max_abs_diff = std::max(max_abs_diff, std::abs(diff));

        if (std::abs(reference_points[i]) > 1e-12) { // Avoid division by zero
            sum_relative_change += (diff * diff) / (reference_points[i] * reference_points[i]);
        }

        sum_abs_diff += std::abs(diff);
    }   

// Check triangle angles if requested
std::vector<double> v_min_angle(0.);
std::vector<double> v_max_angle(0.);
if (check_angles && !elements.empty() && !points.empty()) {
    
//for (const auto& element : elements) {
    for (size_t idx = 0; idx < userData->active_elements.size(); idx++) {
    
        const ElementTriangle2D& element = elements[userData->active_elements[idx]];

    // Create a temporary copy of the element if we need to apply external deformation
     ElementTriangle2D temp_element   = element;
     const ElementTriangle2D* element_to_use = &element;

    // this is a security to avoid meshing due the Boundary triangles
    // if(temp_element.getBoundaryNodeNumber() < 3) {
    //     std::cout << "Skipping element with fewer than 3 nodes." << std::endl;

    //     continue; // Skip boundary touching elements 
    // }

    // Apply external deformation to the temporary element if provided
    if (F_ext != nullptr) {
    temp_element.setExternalDeformation(*F_ext);
    element_to_use = &temp_element;
    }

    // Get vertex indices using the appropriate element
    int idx1 = element_to_use->getNodeIndex(0);
    int idx2 = element_to_use->getNodeIndex(1);
    int idx3 = element_to_use->getNodeIndex(2);

    // Skip if any index is out of range
    if (idx1 >= points.size() || idx2 >= points.size() || idx3 >= points.size()) {
        continue;
    }

    // Get vertex coordinates with effective translations using the appropriate element
    Eigen::Vector2d v1 = points[idx1].coord + element_to_use->getEffectiveTranslation(0);
    Eigen::Vector2d v2 = points[idx2].coord + element_to_use->getEffectiveTranslation(1);
    Eigen::Vector2d v3 = points[idx3].coord + element_to_use->getEffectiveTranslation(2);

    // Calculate sides of the triangle
    Eigen::Vector2d side1 = v2 - v1;
    Eigen::Vector2d side2 = v3 - v2;
    Eigen::Vector2d side3 = v1 - v3;

    // Calculate lengths of sides
    double a = side1.norm();
    double b = side2.norm();
    double c = side3.norm();

    // Use law of cosines to find angles
    double angle1 = std::acos((-side3.dot(side1)) / (side3.norm() * side1.norm())) * 180.0 / M_PI;
    double angle2 = std::acos((-side1.dot(side2)) / (side1.norm() * side2.norm())) * 180.0 / M_PI;
    double angle3 = std::acos((-side2.dot(side3)) / (side2.norm() * side3.norm())) * 180.0 / M_PI;

    // For square lattice, ideal triangles have angles of 90-45-45
    // Check for deviations from ideal angles
    double max_angle = std::max({angle1, angle2, angle3});
    double min_angle = std::min({angle1, angle2, angle3});
    v_max_angle.push_back(max_angle);
    v_min_angle.push_back(min_angle);
    // If the largest angle is over 110 degrees or smallest angle is under 25 degrees,
    // consider the triangle distorted
    //if (max_angle > 119 || min_angle < 31.0) { inden
    //if (max_angle > 100.0 || min_angle < 40.0) {
        if (max_angle > 110.0 || min_angle < 25.0) {

        has_distorted_triangles = true;
        std::cout << "max_angle: " << max_angle << " min_angle: " << min_angle << std::endl;
        break;
    }

    // Also check if the remaining angle deviates too much from 45 degrees
    // Sort angles to find the "middle" angle
    std::array<double, 3> angles = {angle1, angle2, angle3};
    std::sort(angles.begin(), angles.end());
    double middle_angle = angles[1];

    // If the middle angle deviates too much from 45 degrees
    // if (std::abs(middle_angle - 45.0) > 20.0) {
    //     has_distorted_triangles = true;
    //     break;
    // }
    }
}

ChangeMeasures measures;
measures.euclidean_norm = std::sqrt(sum_squared_diff);
measures.max_abs_change = max_abs_diff;
measures.relative_change = sum_relative_change / n;
measures.normalized_euclidean_norm = measures.euclidean_norm / std::sqrt(n);
measures.mean_abs_change = sum_abs_diff / n;
measures.displacement_exceeds_half_lattice = max_abs_diff > (0.5 * lattice_constant);
measures.has_distorted_triangles = has_distorted_triangles;

if (!v_max_angle.empty()) {
    // Find the maximum value in max_angle
    auto max_it = std::max_element(v_max_angle.begin(), v_max_angle.end());
    std::cout << "Maximum of max_angle: " << *max_it << std::endl;
    
    // Find the minimum value in max_angle
    auto min_it = std::min_element(v_max_angle.begin(), v_max_angle.end());
    std::cout << "Minimum of max_angle: " << *min_it << std::endl;
} else {
    std::cout << "max_angle vector is empty" << std::endl;
}

// For min_angle vector
if (!v_min_angle.empty()) {
    // Find the maximum value in min_angle
    auto max_it = std::max_element(v_min_angle.begin(), v_min_angle.end());
    std::cout << "Maximum of min_angle: " << *max_it << std::endl;
    
    // Find the minimum value in min_angle
    auto min_it = std::min_element(v_min_angle.begin(), v_min_angle.end());
    std::cout << "Minimum of min_angle: " << *min_it << std::endl;
} else {
    std::cout << "min_angle vector is empty" << std::endl;
}
return measures;
}


std::vector<int> analyzeElementReduction(
    std::vector<ElementTriangle2D>& elements,
    std::vector<Point2D>& points,
    const UserData* userData)
{
    std::vector<int> m3_activation;
    m3_activation.reserve(userData->active_elements.size());
    
    // Check each element for third reduction operation activation
    for (size_t elem_idx : userData->active_elements) {
        if (elem_idx >= elements.size()) {
            m3_activation.push_back(0); // Invalid element index
            continue;
        }
        
        auto& element = elements[elem_idx];
        if (!element.isInitialized()) {
            m3_activation.push_back(0); // Uninitialized element
            continue;
        }
        
        // Set external deformation and calculate deformation gradient
        element.setExternalDeformation(userData->F_external);
        element.calculate_deformation_gradient(points);
        
        // Get the metric tensor
        Eigen::Matrix2d C = element.getMetricTensor();
        
        // Apply Lagrange reduction with counting
        lagrange::ResultWithCount result = lagrange::reduceAndCount(C);
        
        // Record if third condition is triggered (1) or not (0)
        m3_activation.push_back(result.third_condition_count > 0 ? 1 : 0);
    }
    
    return m3_activation;
}

int compareM3Activation(const std::vector<int>& m3_before, const std::vector<int>& m3_after) 
{
    
    // Make sure vectors have the same size
    if (m3_before.size() != m3_after.size()) {
        std::cerr << "Error: Vector size mismatch in compareM3Activation" << std::endl;
        return false;
    }
    
    // Count differences
    int diff_count = 0;
    for (size_t i = 0; i < m3_before.size(); i++) {
        if (m3_before[i] != m3_after[i]) {
            diff_count++;
        }
    }
    
    // Report if differences were found
    if (diff_count > 0) {
        std::cout << "Found " << diff_count << " differences in M3 activation" << std::endl;
        //return true;
    }
    
    return diff_count;  // No differences found
}

