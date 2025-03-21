#include "../include/geometry/LatticeGenerator.h"

// Generate 2D lattice points (square or triangular)
std::vector<Point2D> LatticeGenerator::generate_2d_lattice(
    int nx, int ny, double lattice_constant, std::string lattice_type) {
    
    std::vector<Point2D> points;
    
    // Calculate domain sizes
    double size_x = nx * lattice_constant;
    double size_y = ny * lattice_constant;
    
    std::cout << "Generating " << lattice_type << " lattice with " << nx << "x" << ny << " cells" << std::endl;
    std::cout << "Lattice constant: " << lattice_constant << std::endl;
    std::cout << "Domain size: " << size_x << " x " << size_y << std::endl;
    
    if (lattice_type == "square") {
        // Generate square lattice points
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                double x = i * lattice_constant;
                double y = j * lattice_constant;
                
                // Add point
                points.push_back(Point2D(x, y));
            }
        }
    }
    else if (lattice_type == "triangular") {
        // Generate triangular lattice points
        // Triangular lattice has the second row offset by a/2
        // where a is the lattice constant
        double hex_height = lattice_constant * std::sqrt(3.0) / 2.0;
        
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                double x = i * lattice_constant + (j % 2) * (lattice_constant / 2.0);
                double y = j * hex_height;
                
                // Add point
                points.push_back(Point2D(x, y));
            }
        }
    }
    else {
        std::cerr << "Error: Unknown lattice type '" << lattice_type << "'. "
                << "Supported types are 'square' and 'triangular'." << std::endl;
    }
    
    return points;
}

// Create periodic copies of the lattice points
std::vector<Point2D> LatticeGenerator::create_periodic_copies(
    const std::vector<Point2D>& original_points,
    const DomainDimensions& domain_dims,
    const std::array<double, 2>& offsets) {
    
    // Define all 9 possible domain translations (including original)
    const std::vector<std::tuple<int, int>> translations = {
        {0, 0},   // Original domain
        // 4 edge neighbors
        {1, 0}, {-1, 0}, {0, 1}, {0, -1},
        // 4 corner neighbors
        {1, 1}, {1, -1}, {-1, 1}, {-1, -1}
    };
    
    // Calculate total number of points
    const size_t total_points = translations.size() * original_points.size();
    
    // Pre-allocate memory for efficiency
    std::vector<Point2D> all_points;
    all_points.reserve(total_points);
    
    // Process each translation
    for (const auto& translation : translations) {
        // Calculate the translation vector
        double tx = std::get<0>(translation) * (domain_dims.size_x + offsets[0]);
        double ty = std::get<1>(translation) * (domain_dims.size_y + offsets[1]);
        Eigen::Vector2d trans_vector(tx, ty);
        
        // Add translated copies of all original points
        for (const auto& point : original_points) {
            // Create translated point
            Eigen::Vector2d translated_coord = point.coord + trans_vector;
            all_points.push_back(Point2D(translated_coord.x(), translated_coord.y()));
        }
    }
    
    std::cout << "Created " << all_points.size() << " points (including "
              << original_points.size() << " original points and "
              << all_points.size() - original_points.size() << " periodic copies)" << std::endl;
    
    return all_points;
}
