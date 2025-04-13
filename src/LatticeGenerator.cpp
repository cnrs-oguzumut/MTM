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
        for (int i = 0; i < nx * ny; i++) {
            int iy = i / nx;
            int ix = i % nx;
            double x = ix * lattice_constant;
            double y = iy * lattice_constant;
            
            // Add point
            points.push_back(Point2D(x, y));
        }
    }
    else if (lattice_type == "triangular") {
        // Generate triangular lattice points
        double hex_height = lattice_constant * std::sqrt(3.0) / 2.0;
        
        for (int i = 0; i < nx * ny; i++) {
            int iy = i / nx;
            int ix = i % nx;
            double x = ix * lattice_constant + (iy % 2) * (lattice_constant / 2.0);
            double y = iy * hex_height;
            
            // Add point
            points.push_back(Point2D(x, y));
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
    const std::array<double, 2>& offsets,
    const Eigen::Matrix2d& F_ext) {
    
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
        Eigen::Vector2d trans_vector_def =F_ext*trans_vector;
        
        // Add translated copies of all original points
        for (const auto& point : original_points) {
            // Create translated point
            Eigen::Vector2d translated_coord = point.coord + trans_vector_def;
            all_points.push_back(Point2D(translated_coord.x(), translated_coord.y()));
        }
    }
    
    std::cout << "Created " << all_points.size() << " points (including "
              << original_points.size() << " original points and "
              << all_points.size() - original_points.size() << " periodic copies)" << std::endl;
    
    return all_points;
}

std::vector<Point2D> LatticeGenerator::generate_2d_lattice_rotated(
    int nx, int ny, double lattice_constant, std::string lattice_type, double rotation_angle_degrees) {
    
    std::vector<Point2D> points;
    
    // Calculate domain sizes
    double size_x = nx * lattice_constant;
    double size_y = ny * lattice_constant;
    
    std::cout << "Generating " << lattice_type << " lattice with " << nx << "x" << ny << " cells" << std::endl;
    std::cout << "Lattice constant: " << lattice_constant << std::endl;
    std::cout << "Domain size: " << size_x << " x " << size_y << std::endl;
    std::cout << "Rotation angle: " << rotation_angle_degrees << " degrees" << std::endl;
    
    // Convert rotation angle to radians
    double rotation_angle_radians = rotation_angle_degrees * M_PI / 180.0;
    
    // Calculate the center of the domain for rotation
    double center_x = size_x / 2.0;
    double center_y = size_y / 2.0;
    
    // Create rotation matrix using Eigen
    Eigen::Matrix2d rotation_matrix;
    rotation_matrix << std::cos(rotation_angle_radians), -std::sin(rotation_angle_radians),
                       std::sin(rotation_angle_radians), std::cos(rotation_angle_radians);
    
    // Calculate expanded grid size to ensure coverage after rotation
    // Use a more precise calculation for the expanded grid
    double diagonal = std::sqrt(size_x * size_x + size_y * size_y);
    double padding_factor = 1.5; // Increase this for better edge coverage
    
    // Calculate how many lattice constants fit across the diagonal (with padding)
    int expanded_nx = std::ceil((diagonal * padding_factor) / lattice_constant);
    int expanded_ny = std::ceil((diagonal * padding_factor) / lattice_constant);
    
    // Generate grid of unrotated points (centered at the origin)
    std::vector<Point2D> unrotated_points;
    
    if (lattice_type == "square") {
        // Calculate a finer lattice constant for better coverage
        double adjusted_lattice_constant = lattice_constant;
        
        // Generate square lattice points with finer spacing
        for (int iy = -expanded_ny; iy <= expanded_ny; iy++) {
            for (int ix = -expanded_nx; ix <= expanded_nx; ix++) {
                double x = center_x + ix * adjusted_lattice_constant;
                double y = center_y + iy * adjusted_lattice_constant;
                unrotated_points.push_back(Point2D(x, y));
            }
        }
    }
    else if (lattice_type == "triangular") {
        // Calculate a finer lattice constant for better coverage
        double adjusted_lattice_constant = lattice_constant;
        double hex_height = adjusted_lattice_constant * std::sqrt(3.0) / 2.0;
        
        // Generate triangular lattice points with finer spacing
        for (int iy = -expanded_ny; iy <= expanded_ny; iy++) {
            for (int ix = -expanded_nx; ix <= expanded_nx; ix++) {
                double x = center_x + ix * adjusted_lattice_constant + 
                          (iy % 2 == 0 ? 0 : adjusted_lattice_constant / 2.0);
                double y = center_y + iy * hex_height;
                unrotated_points.push_back(Point2D(x, y));
            }
        }
    }
    else {
        std::cerr << "Error: Unknown lattice type '" << lattice_type << "'. "
                 << "Supported types are 'square' and 'triangular'." << std::endl;
        return points;
    }
    
    // Apply rotation to all points
    for (const auto& point : unrotated_points) {
        // Translate point to center for rotation
        Eigen::Vector2d translated = point.coord;
        translated.x() -= center_x;
        translated.y() -= center_y;
        
        // Rotate point
        Eigen::Vector2d rotated = rotation_matrix * translated;
        
        // Translate back
        rotated.x() += center_x;
        rotated.y() += center_y;
        
        // Use a small epsilon to ensure points exactly on the boundary are included
        const double epsilon = 1e-10;
        
        // Check if the rotated point is inside the domain (with epsilon)
        if (rotated.x() >= -epsilon && rotated.x() <= size_x + epsilon && 
            rotated.y() >= -epsilon && rotated.y() <= size_y + epsilon) {
            points.push_back(Point2D(rotated.x(), rotated.y()));
        }
    }
    
    std::cout << "Generated " << points.size() << " points after rotation and boundary check" << std::endl;
    
    return points;
}