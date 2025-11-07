#include "../include/defects/DefectAnalysis.h"

#include <iostream>
#include <map>
#include <set>




std::pair<int, std::vector<int>> DefectAnalysis::analyzeDefectsInReferenceConfig(
    const UserData* userData,
    int file_id,
    const Eigen::Matrix<double, 3, 2>& dndx,
    const std::array<double, 2>& offsets,
    const std::vector<int>& original_domain_map,
    const std::vector<std::tuple<double, double>>& translation_map,
    const Point2D& domain_dims_point,
    double element_area,
    bool pbc,
    bool reduction
) {
    std::cout << "\n=== Analyzing defects in reference configuration ===" << std::endl;
    
    // 1. COPY the current points (don't modify originals!)
    std::vector<Point2D> reference_points = userData->points;
    
    // 2. Transform to reference configuration
    Eigen::Matrix2d F_ext_inv = userData->F_external.inverse();
    transformToReferenceConfig(reference_points, F_ext_inv);
    
    std::cout << "Transformed " << reference_points.size() 
              << " points to reference configuration" << std::endl;
    
    // 3. Map to solver array
    alglib::real_1d_array x_ref;
    int n_vars = userData->interior_mapping.size();
    x_ref.setlength(2 * n_vars);
    map_points_to_solver_array(x_ref, reference_points, 
                               userData->interior_mapping, n_vars);
    
    // 4. REMESH in reference configuration
    AdaptiveMesher mesher(
        domain_dims_point,
        offsets,
        original_domain_map,
        translation_map,
        userData->full_mapping,
        1e-6,
        pbc
    );
    
    Eigen::Matrix2d F_identity = Eigen::Matrix2d::Identity();
    auto [ref_elements, ref_active_elements] = mesher.createMesh(
        reference_points, 
        x_ref, 
        F_identity,  // Identity because we're in reference config
        &dndx
    );
    
    std::cout << "Created " << ref_elements.size() 
              << " elements in reference configuration" << std::endl;
    
    // 5. Create temporary UserData for reference config
    UserData refUserData(
        reference_points, 
        ref_elements, 
        userData->calculator,
        userData->energy_function, 
        userData->derivative_function,
        userData->zero_energy, 
        userData->ideal_lattice_parameter, 
        F_identity,  // Identity in reference config
        userData->interior_mapping,
        userData->full_mapping, 
        ref_active_elements, 
        false  // plasticity flag
    );
    
    // 6. Write VTK with proper coordination numbers (simplified output)
    ConfigurationSaver::writeToVTK_DefectAnalysis(
        reference_points, 
        ref_elements, 
        &refUserData, 
        file_id
    );
    
    std::cout << "Saved reference configuration VTK to file " << file_id << std::endl;
    
    // 7. Analyze 5/7 pairs for dislocation detection
    int num_57_pairs = count57Pairs(reference_points, ref_elements, ref_active_elements);
    
    // Calculate coordination numbers for all nodes
    std::vector<int> coordination = calculateCoordinationNumbers(
        reference_points, ref_elements, ref_active_elements
    );
    
    std::cout << "Found " << num_57_pairs << " dislocation (5/7) pairs" << std::endl;
    std::cout << "=== Reference configuration analysis complete ===" << std::endl;
    
    return std::make_pair(num_57_pairs, coordination);
}

void DefectAnalysis::transformToReferenceConfig(
    std::vector<Point2D>& points,
    const Eigen::Matrix2d& F_ext_inv
) {
    for (auto& point : points) {
        Eigen::Vector2d pos(point.coord.x(), point.coord.y());
        Eigen::Vector2d transformed_pos = F_ext_inv * pos;
        point.coord.x() = transformed_pos(0);
        point.coord.y() = transformed_pos(1);
    }
}

std::vector<int> DefectAnalysis::calculateCoordinationNumbers(
    const std::vector<Point2D>& points,
    const std::vector<ElementTriangle2D>& elements,
    const std::vector<size_t>& active_elements
) {
    std::vector<int> coordination(points.size(), 0);
    
    // Count how many elements each node belongs to
    for (size_t elem_idx : active_elements) {
        if (elem_idx >= elements.size()) continue;
        
        const auto& element = elements[elem_idx];
        if (!element.isInitialized()) continue;
        
        for (int i = 0; i < 3; i++) {
            int node_idx = element.getNodeIndex(i);
            if (node_idx < static_cast<int>(coordination.size())) {
                coordination[node_idx]++;
            }
        }
    }
    
    return coordination;
}

// int DefectAnalysis::count57Pairs(
//     const std::vector<Point2D>& points,
//     const std::vector<ElementTriangle2D>& elements,
//     const std::vector<size_t>& active_elements
// ) {
//     // Calculate coordination numbers
//     std::vector<int> coordination = calculateCoordinationNumbers(
//         points, elements, active_elements
//     );
    
//     // Build adjacency map (which nodes are connected through elements)
//     std::map<int, std::set<int>> adjacency;
    
//     for (size_t elem_idx : active_elements) {
//         if (elem_idx >= elements.size()) continue;
        
//         const auto& element = elements[elem_idx];
//         if (!element.isInitialized()) continue;
        
//         int nodes[3];
//         for (int i = 0; i < 3; i++) {
//             nodes[i] = element.getNodeIndex(i);
//         }
        
//         // Add edges between all pairs of nodes in this element
//         for (int i = 0; i < 3; i++) {
//             for (int j = i + 1; j < 3; j++) {
//                 adjacency[nodes[i]].insert(nodes[j]);
//                 adjacency[nodes[j]].insert(nodes[i]);
//             }
//         }
//     }
    
//     // Count 5/7 pairs (adjacent nodes with coordination 5 and 7)
//     int pair_count = 0;
//     std::set<std::pair<int, int>> counted_pairs;
    
//     for (size_t i = 0; i < points.size(); i++) {
//         if (coordination[i] == 5) {
//             // Check if any neighbor has coordination 7
//             for (int neighbor : adjacency[i]) {
//                 if (coordination[neighbor] == 7) {
//                     // Ensure we don't count the same pair twice
//                     int min_idx = std::min(static_cast<int>(i), neighbor);
//                     int max_idx = std::max(static_cast<int>(i), neighbor);
//                     std::pair<int, int> pair(min_idx, max_idx);
                    
//                     if (counted_pairs.find(pair) == counted_pairs.end()) {
//                         counted_pairs.insert(pair);
//                         pair_count++;
                        
//                         std::cout << "5/7 pair found: node " << i 
//                                   << " (coord=" << coordination[i] << ") - node " 
//                                   << neighbor << " (coord=" << coordination[neighbor] 
//                                   << ")" << std::endl;
//                     }
//                 }
//             }
//         }
//     }
    
//     return pair_count;
// }



// int DefectAnalysis::count57Pairs(
//     const std::vector<Point2D>& points,
//     const std::vector<ElementTriangle2D>& elements,
//     const std::vector<size_t>& active_elements
// ) {
//     // Calculate coordination numbers
//     std::vector<int> coordination = calculateCoordinationNumbers(
//         points, elements, active_elements
//     );
    
//     // Build adjacency map (which nodes are connected through elements)
//     std::map<int, std::set<int>> adjacency;
    
//     for (size_t elem_idx : active_elements) {
//         if (elem_idx >= elements.size()) continue;
        
//         const auto& element = elements[elem_idx];
//         if (!element.isInitialized()) continue;
        
//         int nodes[3];
//         for (int i = 0; i < 3; i++) {
//             nodes[i] = element.getNodeIndex(i);
//         }
        
//         // Add edges between all pairs of nodes in this element
//         for (int i = 0; i < 3; i++) {
//             for (int j = i + 1; j < 3; j++) {
//                 adjacency[nodes[i]].insert(nodes[j]);
//                 adjacency[nodes[j]].insert(nodes[i]);
//             }
//         }
//     }
    
//     // Count 5/7 pairs (adjacent nodes with coordination 5 and 7)
//     int pair_count = 0;
//     std::set<std::pair<int, int>> counted_pairs;
    
//     for (size_t i = 0; i < points.size(); i++) {
//         if (coordination[i] == 5) {
//             // Check if any neighbor has coordination 7
//             for (int neighbor : adjacency[i]) {
//                 if (coordination[neighbor] == 7) {
//                     // Ensure we don't count the same pair twice
//                     int min_idx = std::min(static_cast<int>(i), neighbor);
//                     int max_idx = std::max(static_cast<int>(i), neighbor);
//                     std::pair<int, int> pair(min_idx, max_idx);
                    
//                     if (counted_pairs.find(pair) == counted_pairs.end()) {
//                         counted_pairs.insert(pair);
//                         pair_count++;
                        
//                         std::cout << "5/7 pair found: node " << i 
//                                   << " (coord=" << coordination[i] << ") - node " 
//                                   << neighbor << " (coord=" << coordination[neighbor] 
//                                   << ")" << std::endl;
//                     }
//                 }
//             }
//         }
//     }
    
//     return pair_count;
// }

int DefectAnalysis::count57Pairs(
    const std::vector<Point2D>& points,
    const std::vector<ElementTriangle2D>& elements,
    const std::vector<size_t>& active_elements
) {
    // Calculate coordination numbers
    std::vector<int> coordination = calculateCoordinationNumbers(
        points, elements, active_elements
    );
    
    // --- Adjacency map logic removed ---

    // Count all nodes with coordination 5 and 7
    int count_5 = 0;
    int count_7 = 0;
    
    // Iterate through all coordination numbers
    // (Assuming coordination.size() == points.size())
    for (int coord : coordination) {
        if (coord == 5) {
            count_5++;
        } else if (coord == 7) {
            count_7++;
        }
    }
    
    // Sum the counts as requested
    int total_defects = count_5 + count_7;
    
    // Return the total count divided by 2
    return total_defects / 2;
}