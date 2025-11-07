#include "../include/mesh/MeshGenerator.h"
#include "../include/mesh/MeshGenerator.h"

// std::vector<Triangle> MeshGenerator::createTrianglesFromPoints_original(const std::vector<Point2D>& points) {
//     // Convert to CGAL points
//     std::vector<CGALPoint> cgal_points;
//     cgal_points.reserve(points.size());
//     for (const auto& p : points) {
//         cgal_points.push_back(p.to_cgal_point());
//     }
    
//     // Create triangulation with vertex info
//     Delaunay triangulation;
    
//     // Insert points and assign indices
//     std::vector<std::pair<CGALPoint, int>> points_with_info;
//     points_with_info.reserve(cgal_points.size());
    
//     for (size_t i = 0; i < cgal_points.size(); ++i) {
//         points_with_info.emplace_back(cgal_points[i], i);
//     }
    
//     triangulation.insert(points_with_info.begin(), points_with_info.end());
    
//     // Extract triangles
//     std::vector<Triangle> triangles;
//     triangles.reserve(triangulation.number_of_faces());
    
//     for (auto face_it = triangulation.finite_faces_begin();
//          face_it != triangulation.finite_faces_end(); ++face_it) {
//         int v0 = face_it->vertex(0)->info();
//         int v1 = face_it->vertex(1)->info();
//         int v2 = face_it->vertex(2)->info();
//         // Ensure counter-clockwise ordering using cross product
//         const Point2D& p0 = points[v0];
//         const Point2D& p1 = points[v1];
//         const Point2D& p2 = points[v2];

//         // Cross product: (p1-p0) × (p2-p0)
//         double cross_product = (p1.coord.x() - p0.coord.x()) * (p2.coord.y() - p0.coord.y()) - 
//                             (p2.coord.x() - p0.coord.x()) * (p1.coord.y() - p0.coord.y());

//         if (cross_product > 0.0) {
//             // Counter-clockwise - correct order
//             triangles.emplace_back(v0, v1, v2);
//         } else {
//             // Clockwise - swap v1 and v2 to make counter-clockwise
//             triangles.emplace_back(v0, v2, v1);
//         }        
//                 //triangles.emplace_back(v0, v1, v2);
//     }
    
//     return triangles;
// }


//chooses the shortest edges to be from the origin
// std::vector<Triangle> MeshGenerator::createTrianglesFromPoints(const std::vector<Point2D>& points) {
//     // Convert to CGAL points
//     std::vector<CGALPoint> cgal_points;
//     cgal_points.reserve(points.size());
//     for (const auto& p : points) {
//         cgal_points.push_back(p.to_cgal_point());
//     }
    
//     // Create triangulation with vertex info
//     Delaunay triangulation;
    
//     // Insert points and assign indices
//     std::vector<std::pair<CGALPoint, int>> points_with_info;
//     points_with_info.reserve(cgal_points.size());
//     for (size_t i = 0; i < cgal_points.size(); ++i) {
//         points_with_info.emplace_back(cgal_points[i], i);
//     }
//     triangulation.insert(points_with_info.begin(), points_with_info.end());
    
//     // Helper function to reorder triangle with shortest edges from origin
//     auto reorderToRightTriangle = [&](int v0, int v1, int v2) -> std::array<int, 3> {
//         std::array<int, 3> vertices = {v0, v1, v2};
        
//         // Calculate all three edge lengths
//         struct EdgeInfo {
//             int from, to;
//             double length;
//             EdgeInfo(int f, int t, double l) : from(f), to(t), length(l) {}
//         };
        
//         std::vector<EdgeInfo> edges;
//         for (int i = 0; i < 3; ++i) {
//             int v_from = vertices[i];
//             int v_to = vertices[(i+1)%3];
            
//             double dx = points[v_to].coord.x() - points[v_from].coord.x();
//             double dy = points[v_to].coord.y() - points[v_from].coord.y();
//             double length = std::sqrt(dx*dx + dy*dy);
            
//             edges.emplace_back(v_from, v_to, length);
//         }
        
//         // Sort edges by length (shortest first)
//         std::sort(edges.begin(), edges.end(), 
//                   [](const EdgeInfo& a, const EdgeInfo& b) { return a.length < b.length; });
        
//         // Find the vertex that is connected to the two shortest edges
//         // This vertex will be our origin
//         int origin = -1;
        
//         // Check if the two shortest edges share a vertex
//         if (edges[0].from == edges[1].from || edges[0].from == edges[1].to) {
//             origin = edges[0].from;
//         } else if (edges[0].to == edges[1].from || edges[0].to == edges[1].to) {
//             origin = edges[0].to;
//         }
        
//         if (origin == -1) {
//             // Fallback: shouldn't happen in a triangle, but use first vertex
//             origin = vertices[0];
//         }
        
//         // Find the other two vertices
//         std::vector<int> others;
//         for (int v : vertices) {
//             if (v != origin) {
//                 others.push_back(v);
//             }
//         }
        
//         // Order the other two vertices by distance from origin (shortest first)
//         Point2D origin_pos = points[origin];
        
//         auto distance = [&](int vertex) {
//             double dx = points[vertex].coord.x() - origin_pos.coord.x();
//             double dy = points[vertex].coord.y() - origin_pos.coord.y();
//             return std::sqrt(dx*dx + dy*dy);
//         };
        
//         if (distance(others[0]) <= distance(others[1])) {
//             return {origin, others[0], others[1]}; // origin -> closest -> farthest
//         } else {
//             return {origin, others[1], others[0]}; // origin -> closest -> farthest
//         }
//     };
    
//     // Extract triangles with right triangle ordering
//     std::vector<Triangle> triangles;
//     triangles.reserve(triangulation.number_of_faces());
    
//     for (auto face_it = triangulation.finite_faces_begin();
//          face_it != triangulation.finite_faces_end(); ++face_it) {
        
//         int v0 = face_it->vertex(0)->info();
//         int v1 = face_it->vertex(1)->info();
//         int v2 = face_it->vertex(2)->info();
        
//         // Reorder to match right triangle pattern
//         auto ordered = reorderToRightTriangle(v0, v1, v2);
        
//         // Ensure counter-clockwise orientation
//         const Point2D& p0 = points[ordered[0]];
//         const Point2D& p1 = points[ordered[1]];
//         const Point2D& p2 = points[ordered[2]];
        
//         double cross_product = (p1.coord.x() - p0.coord.x()) * (p2.coord.y() - p0.coord.y()) -
//                               (p2.coord.x() - p0.coord.x()) * (p1.coord.y() - p0.coord.y());
        
//         if (cross_product > 0.0) {
//             triangles.emplace_back(ordered[0], ordered[1], ordered[2]);
//         } else {
//             // Flip to maintain CCW while keeping origin first
//             triangles.emplace_back(ordered[0], ordered[2], ordered[1]);
//         }
//     }
    
//     return triangles;
// }



//chooses the ordering such that vectors from origin have the largest angle between them
std::vector<Triangle> MeshGenerator::createTrianglesFromPoints(const std::vector<Point2D>& points) {
    // Convert to CGAL points
    std::vector<CGALPoint> cgal_points;
    cgal_points.reserve(points.size());
    for (const auto& p : points) {
        cgal_points.push_back(p.to_cgal_point());
    }
    
    // Create triangulation with vertex info
    Delaunay triangulation;
    
    // Insert points and assign indices
    std::vector<std::pair<CGALPoint, int>> points_with_info;
    points_with_info.reserve(cgal_points.size());
    for (size_t i = 0; i < cgal_points.size(); ++i) {
        points_with_info.emplace_back(cgal_points[i], i);
    }
    triangulation.insert(points_with_info.begin(), points_with_info.end());
    
    // Helper function to reorder triangle for maximum angle between vectors from origin
    auto reorderForMaxAngle = [&](int v0, int v1, int v2) -> std::array<int, 3> {
        std::array<int, 3> vertices = {v0, v1, v2};
        
        // Try each vertex as origin and find the one that gives maximum angle
        double max_angle = -1.0;
        std::array<int, 3> best_ordering = {v0, v1, v2};
        
        for (int origin_idx = 0; origin_idx < 3; ++origin_idx) {
            int origin = vertices[origin_idx];
            int other1 = vertices[(origin_idx + 1) % 3];
            int other2 = vertices[(origin_idx + 2) % 3];
            
            // Calculate vectors from origin to the other two vertices
            Point2D origin_pos = points[origin];
            Point2D vec1(points[other1].coord.x() - origin_pos.coord.x(),
                        points[other1].coord.y() - origin_pos.coord.y());
            Point2D vec2(points[other2].coord.x() - origin_pos.coord.x(),
                        points[other2].coord.y() - origin_pos.coord.y());
            
            // Calculate angle between vectors using dot product
            double dot_product = vec1.coord.x() * vec2.coord.x() + vec1.coord.y() * vec2.coord.y();
            double mag1 = std::sqrt(vec1.coord.x() * vec1.coord.x() + vec1.coord.y() * vec1.coord.y());
            double mag2 = std::sqrt(vec2.coord.x() * vec2.coord.x() + vec2.coord.y() * vec2.coord.y());
            
            if (mag1 > 1e-12 && mag2 > 1e-12) { // Avoid division by zero
                double cos_angle = dot_product / (mag1 * mag2);
                // Clamp to [-1, 1] to handle numerical errors
                cos_angle = std::max(-1.0, std::min(1.0, cos_angle));
                double angle = std::acos(cos_angle); // Angle in radians [0, π]
                
                if (angle > max_angle) {
                    max_angle = angle;
                    best_ordering = {origin, other1, other2};
                }
            }
        }
        
        return best_ordering;
    };
    
    // Extract triangles with maximum angle ordering
    std::vector<Triangle> triangles;
    triangles.reserve(triangulation.number_of_faces());
    
    for (auto face_it = triangulation.finite_faces_begin();
         face_it != triangulation.finite_faces_end(); ++face_it) {
        
        int v0 = face_it->vertex(0)->info();
        int v1 = face_it->vertex(1)->info();
        int v2 = face_it->vertex(2)->info();
        
        // Reorder for maximum angle between vectors from origin
        auto ordered = reorderForMaxAngle(v0, v1, v2);
        
        // Ensure counter-clockwise orientation
        const Point2D& p0 = points[ordered[0]];
        const Point2D& p1 = points[ordered[1]];
        const Point2D& p2 = points[ordered[2]];
        
        double cross_product = (p1.coord.x() - p0.coord.x()) * (p2.coord.y() - p0.coord.y()) -
                              (p2.coord.x() - p0.coord.x()) * (p1.coord.y() - p0.coord.y());
        
        if (cross_product > 0.0) {
            triangles.emplace_back(ordered[0], ordered[1], ordered[2]);
        } else {
            // Flip to maintain CCW while keeping origin first
            triangles.emplace_back(ordered[0], ordered[2], ordered[1]);
        }
    }
    
    return triangles;
}

//chooses the ordering such that vectors from origin have the largest angle between them <90

// std::vector<Triangle> MeshGenerator::createTrianglesFromPoints(const std::vector<Point2D>& points) {
//     // Convert to CGAL points
//     std::vector<CGALPoint> cgal_points;
//     cgal_points.reserve(points.size());
//     for (const auto& p : points) {
//         cgal_points.push_back(p.to_cgal_point());
//     }
    
//     // Create triangulation with vertex info
//     Delaunay triangulation;
    
//     // Insert points and assign indices
//     std::vector<std::pair<CGALPoint, int>> points_with_info;
//     points_with_info.reserve(cgal_points.size());
//     for (size_t i = 0; i < cgal_points.size(); ++i) {
//         points_with_info.emplace_back(cgal_points[i], i);
//     }
//     triangulation.insert(points_with_info.begin(), points_with_info.end());
    
//     // Helper function to reorder triangle for maximum angle ≤ 90° between vectors from origin
//     auto reorderForMaxAngle = [&](int v0, int v1, int v2) -> std::array<int, 3> {
//         std::array<int, 3> vertices = {v0, v1, v2};
        
//         // Try each vertex as origin and find the one that gives maximum angle ≤ 90°
//         double max_valid_angle = -1.0;  // Maximum angle ≤ π/2
//         std::array<int, 3> best_ordering = {v0, v1, v2};
        
//         for (int origin_idx = 0; origin_idx < 3; ++origin_idx) {
//             int origin = vertices[origin_idx];
//             int other1 = vertices[(origin_idx + 1) % 3];
//             int other2 = vertices[(origin_idx + 2) % 3];
            
//             // Calculate vectors from origin to the other two vertices
//             Point2D origin_pos = points[origin];
//             Point2D vec1(points[other1].coord.x() - origin_pos.coord.x(),
//                         points[other1].coord.y() - origin_pos.coord.y());
//             Point2D vec2(points[other2].coord.x() - origin_pos.coord.x(),
//                         points[other2].coord.y() - origin_pos.coord.y());
            
//             // Calculate angle between vectors using dot product
//             double dot_product = vec1.coord.x() * vec2.coord.x() + vec1.coord.y() * vec2.coord.y();
//             double mag1 = std::sqrt(vec1.coord.x() * vec1.coord.x() + vec1.coord.y() * vec1.coord.y());
//             double mag2 = std::sqrt(vec2.coord.x() * vec2.coord.x() + vec2.coord.y() * vec2.coord.y());
            
//             if (mag1 > 1e-12 && mag2 > 1e-12) { // Avoid division by zero
//                 double cos_angle = dot_product / (mag1 * mag2);
//                 // Clamp to [-1, 1] to handle numerical errors
//                 cos_angle = std::max(-1.0, std::min(1.0, cos_angle));
//                 double angle = std::acos(cos_angle); // Angle in radians [0, π]
                
//                 // Only consider angles ≤ π/2 (90 degrees) to ensure C12 ≥ 0
//                 if (angle <= M_PI/2 + 1e-9) { // Small tolerance for numerical errors
//                     if (angle > max_valid_angle) {
//                         max_valid_angle = angle;
//                         best_ordering = {origin, other1, other2};
//                     }
//                 }
//             }
//         }
        
//         // If no valid angle ≤ 90° was found, fall back to smallest angle > 90°
//         // (this handles degenerate cases where all angles are obtuse)
//         if (max_valid_angle < 0) {
//             double min_obtuse_angle = M_PI + 1.0; // Initialize to > π
            
//             for (int origin_idx = 0; origin_idx < 3; ++origin_idx) {
//                 int origin = vertices[origin_idx];
//                 int other1 = vertices[(origin_idx + 1) % 3];
//                 int other2 = vertices[(origin_idx + 2) % 3];
                
//                 Point2D origin_pos = points[origin];
//                 Point2D vec1(points[other1].coord.x() - origin_pos.coord.x(),
//                             points[other1].coord.y() - origin_pos.coord.y());
//                 Point2D vec2(points[other2].coord.x() - origin_pos.coord.x(),
//                             points[other2].coord.y() - origin_pos.coord.y());
                
//                 double dot_product = vec1.coord.x() * vec2.coord.x() + vec1.coord.y() * vec2.coord.y();
//                 double mag1 = std::sqrt(vec1.coord.x() * vec1.coord.x() + vec1.coord.y() * vec1.coord.y());
//                 double mag2 = std::sqrt(vec2.coord.x() * vec2.coord.x() + vec2.coord.y() * vec2.coord.y());
                
//                 if (mag1 > 1e-12 && mag2 > 1e-12) {
//                     double cos_angle = dot_product / (mag1 * mag2);
//                     cos_angle = std::max(-1.0, std::min(1.0, cos_angle));
//                     double angle = std::acos(cos_angle);
                    
//                     if (angle > M_PI/2 && angle < min_obtuse_angle) {
//                         min_obtuse_angle = angle;
//                         best_ordering = {origin, other1, other2};
//                     }
//                 }
//             }
//         }
        
//         return best_ordering;
//     };
    
//     // Extract triangles with maximum angle ordering
//     std::vector<Triangle> triangles;
//     triangles.reserve(triangulation.number_of_faces());
    
//     for (auto face_it = triangulation.finite_faces_begin();
//          face_it != triangulation.finite_faces_end(); ++face_it) {
        
//         int v0 = face_it->vertex(0)->info();
//         int v1 = face_it->vertex(1)->info();
//         int v2 = face_it->vertex(2)->info();
        
//         // Reorder for maximum angle ≤ 90° between vectors from origin
//         auto ordered = reorderForMaxAngle(v0, v1, v2);
        
//         // Ensure counter-clockwise orientation
//         const Point2D& p0 = points[ordered[0]];
//         const Point2D& p1 = points[ordered[1]];
//         const Point2D& p2 = points[ordered[2]];
        
//         double cross_product = (p1.coord.x() - p0.coord.x()) * (p2.coord.y() - p0.coord.y()) -
//                               (p2.coord.x() - p0.coord.x()) * (p1.coord.y() - p0.coord.y());
        
//         if (cross_product > 0.0) {
//             triangles.emplace_back(ordered[0], ordered[1], ordered[2]);
//         } else {
//             // Flip to maintain CCW while keeping origin first
//             triangles.emplace_back(ordered[0], ordered[2], ordered[1]);
//         }
//     }
    
//     return triangles;
// }

std::pair<std::vector<int>, std::vector<std::tuple<double, double>>>
MeshGenerator::create_domain_maps(int original_domain_size, const DomainDimensions& domain_dims, 
                  const std::array<double, 2>& offsets) {
    // Pre-compute maps for all unique indices - 9 domains in 2D
    std::vector<int> original_domain_map(9 * original_domain_size, 0);
    std::vector<std::tuple<double, double>> translation_map(9 * original_domain_size, {0.0, 0.0});
    
    // Define all 9 possible domain translations (including original)
    const std::vector<std::tuple<int, int>> translations = {
        {0, 0},   // Original domain
        // 4 edge neighbors
        {1, 0}, {-1, 0}, {0, 1}, {0, -1},
        // 4 corner neighbors
        {1, 1}, {1, -1}, {-1, 1}, {-1, -1}
    };
    
    // Process each translation domain
    for (int i = 0; i < translations.size(); i++) {
        for (int j = 0; j < original_domain_size; j++) {
            // Create a unique index for each element across all domains
            int idx = i * original_domain_size + j;
            
            // Map each element back to its corresponding element in original domain
            original_domain_map[idx] = j;
            
            // Calculate translation vector using full domain dimensions + offsets
            int tx = std::get<0>(translations[i]);
            int ty = std::get<1>(translations[i]);
            
            translation_map[idx] = {
                tx * (domain_dims.size_x + offsets[0]),
                ty * (domain_dims.size_y + offsets[1])
            };
            //std::cout<<"translation_map[idx]: "<<std::get<0>(translation_map[idx])<<" "<<std::get<1>(translation_map[idx])<<std::endl;
        }
    }
    
    return {original_domain_map, translation_map};
}


// std::vector<Triangle> MeshGenerator::select_unique_connected_triangles(
//     const std::vector<Point2D>& all_points,
//     const std::vector<Triangle>& duplicated_triangles,
//     const std::vector<int>& original_domain_map,
//     int original_domain_size,
//     double min_jacobian_threshold,
//     double max_edge_length ) {  // Default to negative value to indicate no length filtering
    
//     // First, collect all triangles that meet our connection criteria
//     std::vector<Triangle> connected_triangles;
//     std::cout << "Original domain size: " << original_domain_size << std::endl;
//     std::cout << "duplicated_triangles size: " << duplicated_triangles.size() << std::endl;

//     connected_triangles.clear();
    
//     // Process triangles fully within the original domain
//     for (const auto& tri : duplicated_triangles) {
//         bool fully_inside = true;
//         for (int node_idx = 0; node_idx < 3; node_idx++) {
//             int vertex_value = tri.vertex_indices[node_idx];
            
//             if (vertex_value >= original_domain_size) {
//                 fully_inside = false;
//                 break;
//             }
//         }
        
//         if (fully_inside) {
//             connected_triangles.push_back(tri);
//         }
//     }
    
//     std::cout << "Connected triangles before cross domain boundaries filtering: " 
//               << connected_triangles.size() << std::endl;

//     // Process triangles that cross domain boundaries
//     for (const auto& tri : duplicated_triangles) {
//         // Skip if already processed (fully inside)
//         bool fully_inside = true;
//         for (int node_idx = 0; node_idx < 3; node_idx++) {
//             int vertex_value = tri.vertex_indices[node_idx];
//             if (vertex_value >= original_domain_size) {
//                 fully_inside = false;
//                 break;
//             }
//         }
        
//         if (fully_inside) {
//             continue;
//         }
        
//         // Check if at least one node is in the main domain
//         bool connected_to_main = false;
//         for (int node_idx = 0; node_idx < 3; node_idx++) {
//             int vertex_value = tri.vertex_indices[node_idx];
//             if (vertex_value < original_domain_size) {
//                 connected_to_main = true;
//                 break;
//             }
//         }
        
//         if (connected_to_main) {
//             connected_triangles.push_back(tri);
//         }
//     }
    
//     std::cout << "Connected triangles before uniqueness and quality filtering: " << connected_triangles.size() << std::endl;
    
//     // Now filter for uniqueness and quality
//     std::vector<Triangle> selected_triangles;
//     std::set<std::array<int, 3>> unique_node_combinations;
//     int rejected_by_jacobian = 0;
//     int rejected_by_edge_length = 0;  // Counter for edge length rejections
    
//     // Check if edge length filtering is enabled
//     bool filter_by_edge_length = (max_edge_length > 0.0);
//     if (filter_by_edge_length) {
//         std::cout << "Edge length filtering enabled with max length: " << max_edge_length << std::endl;
//     }
    
//     for (const auto& tri : connected_triangles) {



//         // Create normalized representation for uniqueness checking
//         std::array<int, 3> normalized_nodes;
//         for (int j = 0; j < 3; j++) {
//             normalized_nodes[j] = original_domain_map[tri.vertex_indices[j]];
//         }
//         // Sort to ensure consistent representation
//         std::sort(normalized_nodes.begin(), normalized_nodes.end());
        
//         // Add only if this is a new combination
//         if (unique_node_combinations.find(normalized_nodes) == unique_node_combinations.end()) {
//             // Get vertices
//             const Point2D& v0 = all_points[tri.vertex_indices[0]];
//             const Point2D& v1 = all_points[tri.vertex_indices[1]];
//             const Point2D& v2 = all_points[tri.vertex_indices[2]];
            
//             // Check edge lengths only if filtering is enabled
//             if (filter_by_edge_length) {
//                 double edge1 = (v1.coord - v0.coord).norm();
//                 double edge2 = (v2.coord - v1.coord).norm();
//                 double edge3 = (v0.coord - v2.coord).norm();
                
//                 // Check if any edge exceeds the maximum length
//                 if ((edge1 > max_edge_length || edge2 > max_edge_length || edge3 > max_edge_length) ) {
//                     rejected_by_edge_length++;
//                     continue;
//                 }
//             }
            
//             // Calculate Jacobian for triangle
//             Eigen::Vector2d e1 = Eigen::Vector2d(v1.coord.x() - v0.coord.x(), v1.coord.y() - v0.coord.y());
//             Eigen::Vector2d e2 = Eigen::Vector2d(v2.coord.x() - v0.coord.x(), v2.coord.y() - v0.coord.y());
            
//             Eigen::Matrix2d jacobianMatrix;
//             jacobianMatrix.col(0) = e1;
//             jacobianMatrix.col(1) = e2;
            
//             double detJ = jacobianMatrix.determinant();
            
//             if (detJ > min_jacobian_threshold) {
//                 // Add to selected list if the quality is good
//                 unique_node_combinations.insert(normalized_nodes);
//                 selected_triangles.push_back(tri);
//             } else {
//                 rejected_by_jacobian++;
//             }
//         }
//     }
    
//     std::cout << "Triangles after uniqueness filtering: " << selected_triangles.size() << std::endl;
//     std::cout << "Rejected triangles with small/negative Jacobians: " << rejected_by_jacobian << std::endl;
    
//     if (filter_by_edge_length) {
//         std::cout << "Rejected triangles with edges exceeding max length: " << rejected_by_edge_length << std::endl;
//     }
    
//     return selected_triangles;
// }

// std::vector<Triangle> MeshGenerator::select_unique_connected_triangles_newversion(
//     const std::vector<Point2D>& all_points,
//     const std::vector<Triangle>& duplicated_triangles,
//     const std::vector<int>& original_domain_map,
//     int original_domain_size,
//     double min_jacobian_threshold,
//     double max_edge_length ) {
    
//     // First, collect all triangles that meet our connection criteria
//     std::vector<Triangle> connected_triangles;
//     std::cout << "Original domain size: " << original_domain_size << std::endl;
//     std::cout << "duplicated_triangles size: " << duplicated_triangles.size() << std::endl;

//     connected_triangles.clear();
    
//     // Helper function to determine which copy domain a vertex belongs to
//     auto get_copy_number = [&](int vertex_idx) -> int {
//         return vertex_idx / original_domain_size;
//     };
    
//     // Define allowed copy domains: 0=original, 1=right, 3=top, 5=top-right
//     std::set<int> allowed_copies = {0, 1, 3, 5};
    
//     // Process triangles fully within the original domain
//     for (const auto& tri : duplicated_triangles) {
//         bool fully_inside = true;
//         for (int node_idx = 0; node_idx < 3; node_idx++) {
//             int vertex_value = tri.vertex_indices[node_idx];
            
//             if (vertex_value >= original_domain_size) {
//                 fully_inside = false;
//                 break;
//             }
//         }
        
//         if (fully_inside) {
//             connected_triangles.push_back(tri);
//         }
//     }
    
//     std::cout << "Connected triangles before cross domain boundaries filtering: " 
//               << connected_triangles.size() << std::endl;

//     // Process triangles that cross domain boundaries
//     int rejected_forbidden_copies = 0;
//     int accepted_boundary = 0;
//     int accepted_corner = 0;
    
//     for (const auto& tri : duplicated_triangles) {
//         // Skip if already processed (fully inside)
//         bool fully_inside = true;
//         for (int node_idx = 0; node_idx < 3; node_idx++) {
//             int vertex_value = tri.vertex_indices[node_idx];
//             if (vertex_value >= original_domain_size) {
//                 fully_inside = false;
//                 break;
//             }
//         }
        
//         if (fully_inside) {
//             continue;
//         }
        
//         // Get copy numbers for all three vertices
//         int copy0 = get_copy_number(tri.vertex_indices[0]);
//         int copy1 = get_copy_number(tri.vertex_indices[1]);
//         int copy2 = get_copy_number(tri.vertex_indices[2]);
        
//         std::set<int> vertex_copies = {copy0, copy1, copy2};
        
//         // Special case: Corner triangle with nodes from right(1), top(3), and top-right(5)
//         bool is_corner_triangle = (vertex_copies == std::set<int>{1, 3, 5});
        
//         if (is_corner_triangle) {
//             connected_triangles.push_back(tri);
//             accepted_corner++;
//             continue;
//         }
        
//         // Regular case: at least one node in original, all nodes in allowed copies
//         bool has_original = (copy0 == 0 || copy1 == 0 || copy2 == 0);
//         bool all_in_allowed = (allowed_copies.count(copy0) && 
//                                allowed_copies.count(copy1) && 
//                                allowed_copies.count(copy2));
        
//         if (has_original && all_in_allowed) {
//             connected_triangles.push_back(tri);
//             accepted_boundary++;
//         } else {
//             rejected_forbidden_copies++;
//         }
//     }
    
//     std::cout << "Accepted boundary triangles: " << accepted_boundary << std::endl;
//     std::cout << "Accepted corner triangle (1-3-5): " << accepted_corner << std::endl;
//     std::cout << "Rejected triangles: " << rejected_forbidden_copies << std::endl;
//     std::cout << "Connected triangles before uniqueness and quality filtering: " << connected_triangles.size() << std::endl;
    
//     // Now filter for uniqueness and quality
//     std::vector<Triangle> selected_triangles;
//     std::set<std::array<int, 3>> unique_node_combinations;
//     int rejected_by_jacobian = 0;
//     int rejected_by_edge_length = 0;
    
//     // Check if edge length filtering is enabled
//     bool filter_by_edge_length = (max_edge_length > 0.0);
//     if (filter_by_edge_length) {
//         std::cout << "Edge length filtering enabled with max length: " << max_edge_length << std::endl;
//     }
    
//     for (const auto& tri : connected_triangles) {
//         // Create normalized representation for uniqueness checking
//         std::array<int, 3> normalized_nodes;
//         for (int j = 0; j < 3; j++) {
//             normalized_nodes[j] = original_domain_map[tri.vertex_indices[j]];
//         }
//         // Sort to ensure consistent representation
//         std::sort(normalized_nodes.begin(), normalized_nodes.end());
        
//         // Add only if this is a new combination
//         if (unique_node_combinations.find(normalized_nodes) == unique_node_combinations.end()) {
//             // Get vertices
//             const Point2D& v0 = all_points[tri.vertex_indices[0]];
//             const Point2D& v1 = all_points[tri.vertex_indices[1]];
//             const Point2D& v2 = all_points[tri.vertex_indices[2]];
            
//             // Check edge lengths only if filtering is enabled
//             if (filter_by_edge_length) {
//                 double edge1 = (v1.coord - v0.coord).norm();
//                 double edge2 = (v2.coord - v1.coord).norm();
//                 double edge3 = (v0.coord - v2.coord).norm();
                
//                 // Check if any edge exceeds the maximum length
//                 if ((edge1 > max_edge_length || edge2 > max_edge_length || edge3 > max_edge_length) ) {
//                     rejected_by_edge_length++;
//                     continue;
//                 }
//             }
            
//             // Calculate Jacobian for triangle
//             Eigen::Vector2d e1 = Eigen::Vector2d(v1.coord.x() - v0.coord.x(), v1.coord.y() - v0.coord.y());
//             Eigen::Vector2d e2 = Eigen::Vector2d(v2.coord.x() - v0.coord.x(), v2.coord.y() - v0.coord.y());
            
//             Eigen::Matrix2d jacobianMatrix;
//             jacobianMatrix.col(0) = e1;
//             jacobianMatrix.col(1) = e2;
            
//             double detJ = jacobianMatrix.determinant();
            
//             if (detJ > min_jacobian_threshold) {
//                 // Add to selected list if the quality is good
//                 unique_node_combinations.insert(normalized_nodes);
//                 selected_triangles.push_back(tri);
//             } else {
//                 rejected_by_jacobian++;
//             }
//         }
//     }
    
//     std::cout << "Triangles after uniqueness filtering: " << selected_triangles.size() << std::endl;
//     std::cout << "Rejected triangles with small/negative Jacobians: " << rejected_by_jacobian << std::endl;
    
//     if (filter_by_edge_length) {
//         std::cout << "Rejected triangles with edges exceeding max length: " << rejected_by_edge_length << std::endl;
//     }
    
//     return selected_triangles;

//}

std::vector<Triangle> MeshGenerator::select_unique_connected_triangles(
    const std::vector<Point2D>& all_points,
    const std::vector<Triangle>& duplicated_triangles,
    const std::vector<int>& original_domain_map,
    int original_domain_size,
    double min_jacobian_threshold,
    double max_edge_length ) {  // Default to negative value to indicate no length filtering
    
    // First, collect all triangles that meet our connection criteria
    std::vector<Triangle> connected_triangles;
    std::cout << "Original domain size: " << original_domain_size << std::endl;
    std::cout << "duplicated_triangles size: " << duplicated_triangles.size() << std::endl;

    connected_triangles.clear();
    
    // Process triangles fully within the original domain
    for (const auto& tri : duplicated_triangles) {
        bool fully_inside = true;
        for (int node_idx = 0; node_idx < 3; node_idx++) {
            int vertex_value = tri.vertex_indices[node_idx];
            
            if (vertex_value >= original_domain_size) {
                fully_inside = false;
                break;
            }
        }
        
        if (fully_inside) {
            connected_triangles.push_back(tri);
        }
    }
    
    std::cout << "Connected triangles before cross domain boundaries filtering: " 
              << connected_triangles.size() << std::endl;

    // Process triangles that cross domain boundaries
    for (const auto& tri : duplicated_triangles) {
        // Skip if already processed (fully inside)
        bool fully_inside = true;
        for (int node_idx = 0; node_idx < 3; node_idx++) {
            int vertex_value = tri.vertex_indices[node_idx];
            if (vertex_value >= original_domain_size) {
                fully_inside = false;
                break;
            }
        }
        
        if (fully_inside) {
            continue;
        }
        
        // Check if at least one node is in the main domain
        bool connected_to_main = false;
        for (int node_idx = 0; node_idx < 3; node_idx++) {
            int vertex_value = tri.vertex_indices[node_idx];
            if (vertex_value < original_domain_size) {
                connected_to_main = true;
                break;
            }
        }
        
        if (connected_to_main) {
            connected_triangles.push_back(tri);
        }
    }
    
    std::cout << "Connected triangles before uniqueness and quality filtering: " << connected_triangles.size() << std::endl;
    
    // Helper function to calculate domain priority (higher is better)
    auto calculate_priority = [&](const Triangle& tri) -> int {
        int copy0 = tri.vertex_indices[0] / original_domain_size;
        int copy1 = tri.vertex_indices[1] / original_domain_size;
        int copy2 = tri.vertex_indices[2] / original_domain_size;
        
        std::set<int> preferred_copies = {0, 1, 3, 5};
        bool all_preferred = (preferred_copies.count(copy0) && 
                             preferred_copies.count(copy1) && 
                             preferred_copies.count(copy2));
        bool has_original = (copy0 == 0 || copy1 == 0 || copy2 == 0);
        
        // Priority scheme:
        // 3: All in domain 0 (fully inside original)
        // 2: Has original vertex, all in preferred domains
        // 1: All in preferred domains, no original vertex
        // 0: Has vertices outside preferred domains
        
        if (copy0 == 0 && copy1 == 0 && copy2 == 0) {
            return 3;  // Fully inside original
        } else if (has_original && all_preferred) {
            return 2;  // Boundary triangle in preferred domains
        } else if (all_preferred) {
            return 1;  // Corner triangle in preferred domains
        } else {
            return 0;  // Has vertices in non-preferred domains
        }
    };
    
    // Now filter for uniqueness and quality with domain prioritization
    std::vector<Triangle> selected_triangles;
    std::map<std::array<int, 3>, std::pair<Triangle, int>> unique_triangles_with_priority;
    int rejected_by_jacobian = 0;
    int rejected_by_edge_length = 0;  // Counter for edge length rejections
    
    // Check if edge length filtering is enabled
    bool filter_by_edge_length = (max_edge_length > 0.0);
    if (filter_by_edge_length) {
        std::cout << "Edge length filtering enabled with max length: " << max_edge_length << std::endl;
    }
    
    for (const auto& tri : connected_triangles) {
        // Create normalized representation for uniqueness checking
        std::array<int, 3> normalized_nodes;
        for (int j = 0; j < 3; j++) {
            normalized_nodes[j] = original_domain_map[tri.vertex_indices[j]];
        }
        // Sort to ensure consistent representation
        std::sort(normalized_nodes.begin(), normalized_nodes.end());
        
        // Get vertices
        const Point2D& v0 = all_points[tri.vertex_indices[0]];
        const Point2D& v1 = all_points[tri.vertex_indices[1]];
        const Point2D& v2 = all_points[tri.vertex_indices[2]];
        
        // Check edge lengths only if filtering is enabled
        if (filter_by_edge_length) {
            double edge1 = (v1.coord - v0.coord).norm();
            double edge2 = (v2.coord - v1.coord).norm();
            double edge3 = (v0.coord - v2.coord).norm();
            
            // Check if any edge exceeds the maximum length
            if ((edge1 > max_edge_length || edge2 > max_edge_length || edge3 > max_edge_length) ) {
                rejected_by_edge_length++;
                continue;
            }
        }
        
        // Calculate Jacobian for triangle
        Eigen::Vector2d e1 = Eigen::Vector2d(v1.coord.x() - v0.coord.x(), v1.coord.y() - v0.coord.y());
        Eigen::Vector2d e2 = Eigen::Vector2d(v2.coord.x() - v0.coord.x(), v2.coord.y() - v0.coord.y());
        
        Eigen::Matrix2d jacobianMatrix;
        jacobianMatrix.col(0) = e1;
        jacobianMatrix.col(1) = e2;
        
        double detJ = jacobianMatrix.determinant();
        
        if (detJ > min_jacobian_threshold) {
            // Calculate priority for this triangle
            int priority = calculate_priority(tri);
            
            // Check if this combination exists
            auto it = unique_triangles_with_priority.find(normalized_nodes);
            if (it == unique_triangles_with_priority.end()) {
                // New combination - add it using emplace
                unique_triangles_with_priority.emplace(normalized_nodes, std::make_pair(tri, priority));
            } else {
                // Duplicate found - keep the one with higher priority
                if (priority > it->second.second) {
                    it->second = std::make_pair(tri, priority);
                }
            }
        } else {
            rejected_by_jacobian++;
        }
    }
    
    // Extract the selected triangles from the map
    for (const auto& pair : unique_triangles_with_priority) {
        selected_triangles.push_back(pair.second.first);
    }
    
    std::cout << "Triangles after uniqueness filtering: " << selected_triangles.size() << std::endl;
    std::cout << "Rejected triangles with small/negative Jacobians: " << rejected_by_jacobian << std::endl;
    
    if (filter_by_edge_length) {
        std::cout << "Rejected triangles with edges exceeding max length: " << rejected_by_edge_length << std::endl;
    }
    
    return selected_triangles;
}

std::vector<ElementTriangle2D> MeshGenerator::createElementTri2D(
    const std::vector<Triangle>& unique_triangles,
    const std::vector<Point2D>& points,
    const std::vector<int>& original_domain_map,
    const std::vector<std::tuple<double, double>>& translation_map) {
    
    std::vector<CGALPoint> cgal_points_o;
    cgal_points_o.reserve(points.size());
    for (const auto& p : points) {
        cgal_points_o.push_back(p.to_cgal_point());
    }

    // Process each triangle
    std::vector<ElementTriangle2D> triangles;
    triangles.clear();
    
    for (const auto& tri : unique_triangles) {
        // Get the coordinates of the triangle vertices
        std::array<Eigen::Vector2d, 3> vertex_coords;
        for (int v = 0; v < 3; v++) {
            // Get vertex coordinates from your data structure
            int vertex_idx = tri.vertex_indices[v];
            const auto& point = cgal_points_o[vertex_idx];
            
            // Convert to Eigen::Vector2d
            vertex_coords[v] = Eigen::Vector2d(point.x(), point.y());
        }
        
        // Calculate the Jacobian determinant
        // Using the edge vectors from vertex 0 to vertices 1 and 2
        Eigen::Vector2d v01 = vertex_coords[1] - vertex_coords[0];
        Eigen::Vector2d v02 = vertex_coords[2] - vertex_coords[0];
        
        // Form the Jacobian matrix (2x2 in 2D)
        Eigen::Matrix2d jacobian;
        jacobian.col(0) = v01;
        jacobian.col(1) = v02;
        
        // Calculate determinant (area of the parallelogram, proportional to triangle area)
        double jacobian_det = jacobian.determinant();
        
        // Skip triangles with small or negative Jacobian determinant
        const double min_determinant = 1e-6; // Adjust this threshold as needed
        // if (jacobian_det < min_determinant) {
        //     std::cout << "Skipping degenerate triangle with Jacobian determinant: " 
        //               << jacobian_det << std::endl;
        //     //continue;
        // }
        
        // If the Jacobian is acceptable, create the triangle element
        ElementTriangle2D new_tri;
        // std::cout<<"-------"<<std::endl;
        for (int v = 0; v < 3; v++) {
            new_tri.setNodeIndex(v, original_domain_map[tri.vertex_indices[v]]);
            // std::cout<<"original_domain_map[tri.vertex_indices[v]]: "<<original_domain_map[tri.vertex_indices[v]]<<std::endl;
            
            // Get the translation from translation_map
            auto& trans = translation_map[tri.vertex_indices[v]];
            
            // Convert to Eigen::Vector2d
            Eigen::Vector2d trans_vector;
            trans_vector(0) = std::get<0>(trans);
            trans_vector(1) = std::get<1>(trans);
            // std::cout<<"trans_vector x: "<<trans_vector<<std::endl;
            
            new_tri.setTranslation(v, trans_vector);
        }
        
        triangles.push_back(new_tri);
    }
    
    return triangles;
}
