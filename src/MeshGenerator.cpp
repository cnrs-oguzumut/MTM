#include "../include/mesh/MeshGenerator.h"
#include "../include/mesh/MeshGenerator.h"

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
    
    // Extract triangles
    std::vector<Triangle> triangles;
    triangles.reserve(triangulation.number_of_faces());
    
    for (auto face_it = triangulation.finite_faces_begin();
         face_it != triangulation.finite_faces_end(); ++face_it) {
        int v0 = face_it->vertex(0)->info();
        int v1 = face_it->vertex(1)->info();
        int v2 = face_it->vertex(2)->info();
        
        triangles.emplace_back(v0, v1, v2);
    }
    
    return triangles;
}

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

std::vector<Triangle> MeshGenerator::select_unique_connected_triangles(
    const std::vector<Point2D>& all_points,
    const std::vector<Triangle>& duplicated_triangles,
    const std::vector<int>& original_domain_map,
    int original_domain_size,
    double min_jacobian_threshold) {
    
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
    
    // Now filter for uniqueness and quality
    std::vector<Triangle> selected_triangles;
    std::set<std::array<int, 3>> unique_node_combinations;
    int rejected_by_jacobian = 0;
    
    for (const auto& tri : connected_triangles) {
        // Create normalized representation for uniqueness checking
        std::array<int, 3> normalized_nodes;
        for (int j = 0; j < 3; j++) {
            normalized_nodes[j] = original_domain_map[tri.vertex_indices[j]];
        }
        // Sort to ensure consistent representation
        std::sort(normalized_nodes.begin(), normalized_nodes.end());
        
        // Add only if this is a new combination
        if (unique_node_combinations.find(normalized_nodes) == unique_node_combinations.end()) {
            // Calculate Jacobian
            const Point2D& v0 = all_points[tri.vertex_indices[0]];
            const Point2D& v1 = all_points[tri.vertex_indices[1]];
            const Point2D& v2 = all_points[tri.vertex_indices[2]];
            
            // Calculate Jacobian for triangle
            Eigen::Vector2d e1 = Eigen::Vector2d(v1.coord.x() - v0.coord.x(), v1.coord.y() - v0.coord.y());
            Eigen::Vector2d e2 = Eigen::Vector2d(v2.coord.x() - v0.coord.x(), v2.coord.y() - v0.coord.y());
            
            Eigen::Matrix2d jacobianMatrix;
            jacobianMatrix.col(0) = e1;
            jacobianMatrix.col(1) = e2;
            
            double detJ = jacobianMatrix.determinant();
            
            if (detJ > min_jacobian_threshold) {
                // Add to selected list if the quality is good
                unique_node_combinations.insert(normalized_nodes);
                selected_triangles.push_back(tri);
            } else {
                rejected_by_jacobian++;
            }
        }
    }
    
    std::cout << "Triangles after uniqueness filtering: " << selected_triangles.size() << std::endl;
    std::cout << "Rejected triangles with small/negative Jacobians: " << rejected_by_jacobian << std::endl;
    
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
