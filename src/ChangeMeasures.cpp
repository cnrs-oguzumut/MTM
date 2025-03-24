#include <cmath>
#include <algorithm>
#include <array>
#include <map>
#include <iostream>
#include <vector>
#include <stdexcept>
#include "../include/output/ChangeMeasures.h"


ChangeMeasures computeChangeMeasures(const alglib::real_1d_array& current_points,
    const alglib::real_1d_array& reference_points,
    double lattice_constant,
    const std::vector<ElementTriangle2D>& elements,
    const std::vector<Point2D>& points,
    bool check_angles,
    const Eigen::Matrix2d* F_ext ) {  // Fixed extra parenthesis
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
if (check_angles && !elements.empty() && !points.empty()) {
for (const auto& element : elements) {
// Create a temporary copy of the element if we need to apply external deformation
ElementTriangle2D temp_element = element;
const ElementTriangle2D* element_to_use = &element;

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

// If the largest angle is over 110 degrees or smallest angle is under 25 degrees,
// consider the triangle distorted
if (max_angle > 110.0 || min_angle < 25.0) {
has_distorted_triangles = true;
std::cout<<"max_angle: "<<max_angle<<" min_angle: "<<min_angle<<std::endl;
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

return measures;
}
int analyzeElementReduction(
    std::vector<ElementTriangle2D>& elements, 
    std::vector<Point2D>& points,
    const UserData* userData) 
{
    std::cout << "Analyzing Lagrange reduction metrics for finite elements..." << std::endl;
    
    int total_elements = 0;
    int elements_with_reduction = 0;
    int total_third_condition_triggers = 0;
    int max_triggers_in_element = 0;
    
    // For histogram of third condition count occurrences
    std::map<int, int> count_histogram;
    
    // Loop through all elements
    for (size_t elem_idx : userData->active_elements) {
        if (elem_idx >= elements.size()) continue;
        
        auto& element = elements[elem_idx];
        if (!element.isInitialized()) continue;
        
        total_elements++;
        
        // Set external deformation and calculate deformation gradient
        element.setExternalDeformation(userData->F_external);
        element.calculate_deformation_gradient(points);
        
        // Get the metric tensor
        Eigen::Matrix2d C = element.getMetricTensor();
        std::cout << "C: " << C << std::endl;
        
        // Apply Lagrange reduction with counting
        lagrange::ResultWithCount result = lagrange::reduceAndCount(C);
        
        // Record statistics
        if (result.third_condition_count > 0) {
            elements_with_reduction++;
            total_third_condition_triggers += result.third_condition_count;
            max_triggers_in_element = std::max(max_triggers_in_element, result.third_condition_count);
        }
        
        // Add to histogram
        count_histogram[result.third_condition_count]++;
    }
    
    // Print statistics
    std::cout << "Total elements analyzed: " << total_elements << std::endl;
    std::cout << "Elements with third condition triggered: " << elements_with_reduction
              << " (" << (100.0 * elements_with_reduction / total_elements) << "%)" << std::endl;
    std::cout << "Total third condition triggers: " << total_third_condition_triggers << std::endl;
    
    // Return the number of elements with reduction
    return elements_with_reduction;
}