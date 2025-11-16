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
    
    std::cout << "Generating " << lattice_type << " lattice" << std::endl;
    std::cout << "Requested: nx=" << nx << ", ny=" << ny << std::endl;
    std::cout << "Lattice constant: " << lattice_constant << std::endl;
    std::cout << "Rotation angle: " << rotation_angle_degrees << " degrees" << std::endl;
    
    // Convert rotation angle to radians
    double rotation_angle = rotation_angle_degrees * M_PI / 180.0;
    
    // Adjust nx, ny for CSL commensurate periodic boundaries
    int nx_adjusted = nx;
    int ny_adjusted = ny;
    double offset_x = 0.0;
    double offset_y = 0.0;
    
    if (lattice_type == "square" && std::abs(rotation_angle_degrees) > 0.01) {
        
        double angle_tol = 0.5; // tolerance in degrees
        
        if (std::abs(rotation_angle_degrees - 45.0) < angle_tol) {
            std::cout << "\nΣ=2 rotation (45°, m=1, n=1)" << std::endl;
            
            double sqrt2 = std::sqrt(2.0);
            offset_x = lattice_constant / sqrt2;
            offset_y = lattice_constant / sqrt2;
            
            // For 45°, commensurate condition: (nx + 1/√2)/√2 ≈ integer
            // This gives valid sequence: 1, 2, 5, 7, 9, 12, 14, 16, 19, 21, 23, 26, ...
            
            int best_nx = nx;
            int best_ny = ny;
            double best_error_x = 1.0;
            double best_error_y = 1.0;
            
            // Search in range [nx-2, nx+2] for best commensurate value
            for (int test_nx = std::max(1, nx - 2); test_nx <= nx + 2; test_nx++) {
                double k_x = (test_nx + 1.0/sqrt2) / sqrt2;
                double error_x = std::abs(k_x - std::round(k_x));
                if (error_x < best_error_x) {
                    best_error_x = error_x;
                    best_nx = test_nx;
                }
            }
            
            for (int test_ny = std::max(1, ny - 2); test_ny <= ny + 2; test_ny++) {
                double k_y = (test_ny + 1.0/sqrt2) / sqrt2;
                double error_y = std::abs(k_y - std::round(k_y));
                if (error_y < best_error_y) {
                    best_error_y = error_y;
                    best_ny = test_ny;
                }
            }
            
            nx_adjusted = best_nx;
            ny_adjusted = best_ny;
            
            int k_x = static_cast<int>(std::round((nx_adjusted + 1.0/sqrt2) / sqrt2));
            int k_y = static_cast<int>(std::round((ny_adjusted + 1.0/sqrt2) / sqrt2));
            
            std::cout << "CSL parameters: Σ=2, √Σ=√2" << std::endl;
            std::cout << "Valid nx sequence: 1,2,5,7,9,12,14,16,19,21,23,26,..." << std::endl;
            std::cout << "Periodicity check X: (nx+1/√2)/√2 = " << (nx_adjusted + 1.0/sqrt2)/sqrt2 
                      << " ≈ " << k_x << " (error=" << best_error_x << ")" << std::endl;
            std::cout << "Periodicity check Y: (ny+1/√2)/√2 = " << (ny_adjusted + 1.0/sqrt2)/sqrt2 
                      << " ≈ " << k_y << " (error=" << best_error_y << ")" << std::endl;
            std::cout << "Commensurate multipliers: k_x=" << k_x << ", k_y=" << k_y << std::endl;
            
            if (nx_adjusted != nx || ny_adjusted != ny) {
                std::cout << "*** ADJUSTED: nx=" << nx_adjusted << " (was " << nx 
                          << "), ny=" << ny_adjusted << " (was " << ny << ") ***" << std::endl;
            } else {
                std::cout << "Values are commensurate: nx=" << nx_adjusted << ", ny=" << ny_adjusted << std::endl;
            }
        }
        else if (std::abs(rotation_angle_degrees - 26.565) < angle_tol) {
            std::cout << "\nΣ=5 rotation (26.565°, m=1, n=2)" << std::endl;
            
            double sqrt5 = std::sqrt(5.0);
            double cos_theta = std::cos(rotation_angle);
            double sin_theta = std::sin(rotation_angle);
            offset_x = lattice_constant * cos_theta;
            offset_y = lattice_constant * sin_theta;
            
            // For Σ=5: (nx*a + offset_x) / (a*√5) ≈ integer
            int best_nx = nx;
            int best_ny = ny;
            double best_error_x = 1.0;
            double best_error_y = 1.0;
            
            for (int test_nx = std::max(1, nx - 2); test_nx <= nx + 2; test_nx++) {
                double k_x = (test_nx + offset_x / lattice_constant) / sqrt5;
                double error_x = std::abs(k_x - std::round(k_x));
                if (error_x < best_error_x) {
                    best_error_x = error_x;
                    best_nx = test_nx;
                }
            }
            
            for (int test_ny = std::max(1, ny - 2); test_ny <= ny + 2; test_ny++) {
                double k_y = (test_ny + offset_y / lattice_constant) / sqrt5;
                double error_y = std::abs(k_y - std::round(k_y));
                if (error_y < best_error_y) {
                    best_error_y = error_y;
                    best_ny = test_ny;
                }
            }
            
            nx_adjusted = best_nx;
            ny_adjusted = best_ny;
            
            int k_x = static_cast<int>(std::round((nx_adjusted + offset_x/lattice_constant) / sqrt5));
            int k_y = static_cast<int>(std::round((ny_adjusted + offset_y/lattice_constant) / sqrt5));
            
            std::cout << "CSL parameters: Σ=5, √Σ=√5" << std::endl;
            std::cout << "Commensurate multipliers: k_x=" << k_x << ", k_y=" << k_y << std::endl;
            
            if (nx_adjusted != nx || ny_adjusted != ny) {
                std::cout << "*** ADJUSTED: nx=" << nx_adjusted << " (was " << nx 
                          << "), ny=" << ny_adjusted << " (was " << ny << ") ***" << std::endl;
            } else {
                std::cout << "Values are commensurate: nx=" << nx_adjusted << ", ny=" << ny_adjusted << std::endl;
            }
        }
        else if (std::abs(rotation_angle_degrees - 63.435) < angle_tol) {
            std::cout << "\nΣ=5 rotation (63.435°, m=2, n=1)" << std::endl;
            
            double sqrt5 = std::sqrt(5.0);
            double cos_theta = std::cos(rotation_angle);
            double sin_theta = std::sin(rotation_angle);
            offset_x = lattice_constant * cos_theta;
            offset_y = lattice_constant * sin_theta;
            
            int best_nx = nx;
            int best_ny = ny;
            double best_error_x = 1.0;
            double best_error_y = 1.0;
            
            for (int test_nx = std::max(1, nx - 2); test_nx <= nx + 2; test_nx++) {
                double k_x = (test_nx + offset_x / lattice_constant) / sqrt5;
                double error_x = std::abs(k_x - std::round(k_x));
                if (error_x < best_error_x) {
                    best_error_x = error_x;
                    best_nx = test_nx;
                }
            }
            
            for (int test_ny = std::max(1, ny - 2); test_ny <= ny + 2; test_ny++) {
                double k_y = (test_ny + offset_y / lattice_constant) / sqrt5;
                double error_y = std::abs(k_y - std::round(k_y));
                if (error_y < best_error_y) {
                    best_error_y = error_y;
                    best_ny = test_ny;
                }
            }
            
            nx_adjusted = best_nx;
            ny_adjusted = best_ny;
            
            int k_x = static_cast<int>(std::round((nx_adjusted + offset_x/lattice_constant) / sqrt5));
            int k_y = static_cast<int>(std::round((ny_adjusted + offset_y/lattice_constant) / sqrt5));
            
            std::cout << "CSL parameters: Σ=5, √Σ=√5" << std::endl;
            std::cout << "Commensurate multipliers: k_x=" << k_x << ", k_y=" << k_y << std::endl;
            
            if (nx_adjusted != nx || ny_adjusted != ny) {
                std::cout << "*** ADJUSTED: nx=" << nx_adjusted << " (was " << nx 
                          << "), ny=" << ny_adjusted << " (was " << ny << ") ***" << std::endl;
            } else {
                std::cout << "Values are commensurate: nx=" << nx_adjusted << ", ny=" << ny_adjusted << std::endl;
            }
        }
        else if (std::abs(rotation_angle_degrees - 36.87) < angle_tol) {
            std::cout << "\nΣ=25 rotation (36.87°, m=3, n=4)" << std::endl;
            
            double sqrt25 = 5.0;
            double cos_theta = std::cos(rotation_angle);
            double sin_theta = std::sin(rotation_angle);
            offset_x = lattice_constant * cos_theta;
            offset_y = lattice_constant * sin_theta;
            
            int best_nx = nx;
            int best_ny = ny;
            double best_error_x = 1.0;
            double best_error_y = 1.0;
            
            for (int test_nx = std::max(1, nx - 2); test_nx <= nx + 2; test_nx++) {
                double k_x = (test_nx + offset_x / lattice_constant) / sqrt25;
                double error_x = std::abs(k_x - std::round(k_x));
                if (error_x < best_error_x) {
                    best_error_x = error_x;
                    best_nx = test_nx;
                }
            }
            
            for (int test_ny = std::max(1, ny - 2); test_ny <= ny + 2; test_ny++) {
                double k_y = (test_ny + offset_y / lattice_constant) / sqrt25;
                double error_y = std::abs(k_y - std::round(k_y));
                if (error_y < best_error_y) {
                    best_error_y = error_y;
                    best_ny = test_ny;
                }
            }
            
            nx_adjusted = best_nx;
            ny_adjusted = best_ny;
            
            int k_x = static_cast<int>(std::round((nx_adjusted + offset_x/lattice_constant) / sqrt25));
            int k_y = static_cast<int>(std::round((ny_adjusted + offset_y/lattice_constant) / sqrt25));
            
            std::cout << "CSL parameters: Σ=25, √Σ=5" << std::endl;
            std::cout << "Commensurate multipliers: k_x=" << k_x << ", k_y=" << k_y << std::endl;
            
            if (nx_adjusted != nx || ny_adjusted != ny) {
                std::cout << "*** ADJUSTED: nx=" << nx_adjusted << " (was " << nx 
                          << "), ny=" << ny_adjusted << " (was " << ny << ") ***" << std::endl;
            } else {
                std::cout << "Values are commensurate: nx=" << nx_adjusted << ", ny=" << ny_adjusted << std::endl;
            }
        }
        else if (std::abs(rotation_angle_degrees - 53.13) < angle_tol) {
            std::cout << "\nΣ=25 rotation (53.13°, m=4, n=3)" << std::endl;
            
            double sqrt25 = 5.0;
            double cos_theta = std::cos(rotation_angle);
            double sin_theta = std::sin(rotation_angle);
            offset_x = lattice_constant * cos_theta;
            offset_y = lattice_constant * sin_theta;
            
            int best_nx = nx;
            int best_ny = ny;
            double best_error_x = 1.0;
            double best_error_y = 1.0;
            
            for (int test_nx = std::max(1, nx - 2); test_nx <= nx + 2; test_nx++) {
                double k_x = (test_nx + offset_x / lattice_constant) / sqrt25;
                double error_x = std::abs(k_x - std::round(k_x));
                if (error_x < best_error_x) {
                    best_error_x = error_x;
                    best_nx = test_nx;
                }
            }
            
            for (int test_ny = std::max(1, ny - 2); test_ny <= ny + 2; test_ny++) {
                double k_y = (test_ny + offset_y / lattice_constant) / sqrt25;
                double error_y = std::abs(k_y - std::round(k_y));
                if (error_y < best_error_y) {
                    best_error_y = error_y;
                    best_ny = test_ny;
                }
            }
            
            nx_adjusted = best_nx;
            ny_adjusted = best_ny;
            
            int k_x = static_cast<int>(std::round((nx_adjusted + offset_x/lattice_constant) / sqrt25));
            int k_y = static_cast<int>(std::round((ny_adjusted + offset_y/lattice_constant) / sqrt25));
            
            std::cout << "CSL parameters: Σ=25, √Σ=5" << std::endl;
            std::cout << "Commensurate multipliers: k_x=" << k_x << ", k_y=" << k_y << std::endl;
            
            if (nx_adjusted != nx || ny_adjusted != ny) {
                std::cout << "*** ADJUSTED: nx=" << nx_adjusted << " (was " << nx 
                          << "), ny=" << ny_adjusted << " (was " << ny << ") ***" << std::endl;
            } else {
                std::cout << "Values are commensurate: nx=" << nx_adjusted << ", ny=" << ny_adjusted << std::endl;
            }
        }
        else {
            std::cout << "\nWarning: Non-standard CSL angle. No adjustment applied." << std::endl;
            std::cout << "For commensurate boundaries, use standard CSL angles:" << std::endl;
            std::cout << "  45° (Σ=2), 26.565° or 63.435° (Σ=5)" << std::endl;
            std::cout << "  36.87° or 53.13° (Σ=25), 18.435° or 71.565° (Σ=10)" << std::endl;
        }
        
        std::cout << "Offsets for periodic copies: (" << offset_x << ", " << offset_y << ")" << std::endl;
    }
    
    // Calculate domain sizes with adjusted values
    double size_x = (nx_adjusted) * lattice_constant + 2*lattice_constant;
    double size_y = (ny_adjusted) * lattice_constant + lattice_constant;
    
    std::cout << "Final domain size: " << size_x << " x " << size_y << std::endl;
    
    // Domain center for rotation
    double cx = size_x / 2.0;
    double cy = size_y / 2.0;
    
    // Rotation matrix components
    double cos_theta = std::cos(rotation_angle);
    double sin_theta = std::sin(rotation_angle);
    
    // Generate extended lattice to ensure coverage after rotation
    int margin = static_cast<int>(std::ceil(std::max(nx_adjusted, ny_adjusted) * 0.8));
    int nx_extended = nx_adjusted + 2 * margin;
    int ny_extended = ny_adjusted + 2 * margin;
    
    if (lattice_type == "square") {
        for (int iy = 0; iy < ny_extended; iy++) {
            for (int ix = 0; ix < nx_extended; ix++) {
                double x = (ix - margin) * lattice_constant;
                double y = (iy - margin) * lattice_constant;
                
                double x_shifted = x - cx;
                double y_shifted = y - cy;
                double x_rot = cos_theta * x_shifted - sin_theta * y_shifted + cx;
                double y_rot = sin_theta * x_shifted + cos_theta * y_shifted + cy;
                
                if (x_rot >= 0.0 && x_rot <= size_x && 
                    y_rot >= 0.0 && y_rot <= size_y) {
                    points.push_back(Point2D(x_rot, y_rot));
                }
            }
        }
    }
    else if (lattice_type == "triangular") {
        double hex_height = lattice_constant * std::sqrt(3.0) / 2.0;
        
        for (int iy = 0; iy < ny_extended; iy++) {
            for (int ix = 0; ix < nx_extended; ix++) {
                double x = (ix - margin) * lattice_constant + (iy % 2) * (lattice_constant / 2.0);
                double y = (iy - margin) * hex_height;
                
                double x_shifted = x - cx;
                double y_shifted = y - cy;
                double x_rot = cos_theta * x_shifted - sin_theta * y_shifted + cx;
                double y_rot = sin_theta * x_shifted + cos_theta * y_shifted + cy;
                
                if (x_rot >= 0.0 && x_rot <= size_x && 
                    y_rot >= 0.0 && y_rot <= size_y) {
                    points.push_back(Point2D(x_rot, y_rot));
                }
            }
        }
    }
    else {
        std::cerr << "Error: Unknown lattice type '" << lattice_type << "'. "
                  << "Supported types are 'square' and 'triangular'." << std::endl;
    }
    
    std::cout << "Generated " << points.size() << " points inside target domain" << std::endl;
    std::cout << "\n========================================\n" << std::endl;
    
    return points;
}






 
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <set>
#include <algorithm>
#include <stdexcept>

// Use M_PI if available, otherwise define it
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- CHANGE 1: Add a custom comparator for std::set ---
// This tells std::set how to sort Point2D objects
// since we cannot add operator< to the class itself.
struct ComparePoint2D {
    bool operator()(const Point2D& a, const Point2D& b) const {
        double tol = 1e-9; // Tolerance for floating point
        
        // Compare x-coordinates
        if (std::abs(a.coord.x() - b.coord.x()) > tol) {
            return a.coord.x() < b.coord.x();
        }
        
        // If x is "equal", compare y-coordinates
        if (std::abs(a.coord.y() - b.coord.y()) > tol) {
            return a.coord.y() < b.coord.y();
        }
        
        // If both are "equal", 'a' is not less than 'b'
        return false;
    }
};
// --- END CHANGE 1 ---


/**
 * @brief Generates a periodic rotated lattice in a rectangular domain
 * based on Coincident Site Lattice (CSL) constraints.
 */
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <set>          // We still use set
#include <utility>      // For std::pair
#include <algorithm>
#include <stdexcept>



std::pair<std::vector<Point2D>, Point2D> 
LatticeGenerator::generate_periodic_rotated_lattice_v2(
    int m, int n, int kx, int ky,
    double lattice_constant, std::string lattice_type) {

    if (lattice_type != "square") {
        throw std::runtime_error("Error: Only 'square' lattice type is supported.");
    }
    if ((m == 0 && n == 0) || kx <= 0 || ky <= 0) {
        throw std::runtime_error("Error: CSL params invalid.");
    }

    // --- 1. Calculate CSL/Periodicity Parameters ---
    int sigma_int = m * m + n * n;
    double sigma = static_cast<double>(sigma_int);
    double sqrt_sigma = std::sqrt(sigma);
    double rotation_angle = std::atan2(static_cast<double>(n), static_cast<double>(m));
    double cos_theta = static_cast<double>(m) / sqrt_sigma;
    double sin_theta = static_cast<double>(n) / sqrt_sigma;
    
    // --- 2. Define Rotated Lattice Vectors (v1, v2) ---
    Point2D v1(lattice_constant * cos_theta, lattice_constant * sin_theta);
    Point2D v2(lattice_constant * -sin_theta, lattice_constant * cos_theta);

    // --- 3. Define the Periodic Box Dimensions ---
    double csl_cell_size = lattice_constant * sqrt_sigma;
    double size_x = static_cast<double>(kx) * csl_cell_size;
    double size_y = static_cast<double>(ky) * csl_cell_size;

    std::cout << "--- Periodic Rotated SQUARE Lattice Generation ---" << std::endl;
    std::cout << "Unit Cell Size (shortestLx): " << csl_cell_size << std::endl;
    std::cout << "Final Domain Size (Total PBC Offset): " << size_x << " x " << size_y << std::endl;

    // ======================================================================
    // --- STEP 4: Find Basis Atoms and Max Coordinates ---
    // ======================================================================
    
    std::set<std::pair<long long, long long>> basis_atoms_int;
    const double precision = 1e9;
    long long L_int = static_cast<long long>(std::round(csl_cell_size * precision));
    int search_range = std::abs(m) + std::abs(n) + 1;

    for (int i = -search_range; i <= search_range; ++i) {
        for (int j = -search_range; j <= search_range; ++j) {
            
            double x = i * v1.coord.x() + j * v2.coord.x();
            double y = i * v1.coord.y() + j * v2.coord.y();

            double x_pbc = std::fmod(x, csl_cell_size);
            double y_pbc = std::fmod(y, csl_cell_size);
            if (x_pbc < 0) x_pbc += csl_cell_size;
            if (y_pbc < 0) y_pbc += csl_cell_size;

            long long x_int = static_cast<long long>(std::round(x_pbc * precision));
            long long y_int = static_cast<long long>(std::round(y_pbc * precision));

            if (x_int == L_int) x_int = 0;
            if (y_int == L_int) y_int = 0;
            
            basis_atoms_int.insert({x_int, y_int});
        }
    }
    std::cout << "Found " << basis_atoms_int.size() << " basis atoms (Sigma)." << std::endl;

    // --- 4b. Find max(P.x), max(P.y) from basis atoms ---
    double max_px = 0.0;
    double max_py = 0.0;
    std::vector<Point2D> basis_atoms; // Temporary vector to store basis

    for (const auto& p_int : basis_atoms_int) {
        double x_basis = static_cast<double>(p_int.first) / precision;
        double y_basis = static_cast<double>(p_int.second) / precision;
        
        basis_atoms.push_back(Point2D(x_basis, y_basis));
        
        if (x_basis > max_px) max_px = x_basis;
        if (y_basis > max_py) max_py = y_basis;
    }

    // --- 4c. Calculate your "user_defined_offset" ---
    double user_offset_x = csl_cell_size - max_px;
    double user_offset_y = csl_cell_size - max_py;
    Point2D user_defined_offset(user_offset_x, user_offset_y);
    
    std::cout << "Max basis atom coord (max_Px, max_Py): (" << max_px << ", " << max_py << ")" << std::endl;
    std::cout << "User Offset (L_unit - max_P): (" << user_offset_x << ", " << user_offset_y << ")" << std::endl;


    // ======================================================================
    // --- STEP 5: Tile the Basis Atoms ---
    // ======================================================================
    
    std::vector<Point2D> final_points;
    for (int ix = 0; ix < kx; ++ix) {
        for (int iy = 0; iy < ky; ++iy) {
            double offset_x = ix * csl_cell_size;
            double offset_y = iy * csl_cell_size;

            for (const auto& p_basis : basis_atoms) {
                final_points.push_back(
                    Point2D(p_basis.coord.x() + offset_x, p_basis.coord.y() + offset_y)
                );
            }
        }
    }

    std::cout << "Generated " << final_points.size() << " total points." << std::endl;
    std::cout << "========================================\n" << std::endl;

    // --- 6. Return the pair ---
    return std::make_pair(final_points, user_defined_offset);
}


// #include <vector>
// #include <string>
// #include <cmath>
// #include <iostream>
// #include <set>
// #include <utility>      // For std::pair
// #include <stdexcept>    // For std::runtime_error
// #include <algorithm>    // For std::abs


// // Use M_PI if available
// #ifndef M_PI
// #define M_PI 3.14159265358979323846
// #endif

// /**
//  * @brief Generates a periodic ROTATED TRIANGULAR lattice in a 
//  * rectangular domain based on Coincident Site Lattice (CSL) theory.
//  *
//  * NOTE: The CSL theory for triangular lattices is different from square lattices.
//  * m and n define the rotation, Sigma, and the rectangular unit cell.
//  * - Angle: tan(theta) = (n*sqrt(3)) / (2m + n)
//  * - Sigma: Sigma = m^2 + n^2 + mn
//  * - Rectangular Unit Cell: (a*sqrt(Sigma)) x (a*sqrt(3*Sigma))
//  * - Atoms in Unit Cell: 2*Sigma
//  *
//  * @param m CSL integer parameter
//  * @param n CSL integer parameter
//  * @param kx Integer multiplier for the CSL supercell in X
//  * @param ky Integer multiplier for the CSL supercell in Y
//  * @param lattice_constant The 'a' spacing (atom-atom distance)
//  * @return std::vector<Point2D> List of atom positions
//  */
// std::vector<Point2D> LatticeGenerator::generate_periodic_rotated_triangular_lattice(
//     int m, int n, int kx, int ky,
//     double lattice_constant) {

//     std::vector<Point2D> final_points;
//     const double a = lattice_constant;
//     const double sqrt3 = std::sqrt(3.0);

//     if ((m == 0 && n == 0) || kx <= 0 || ky <= 0) {
//         throw std::runtime_error("Error: 'm' and 'n' cannot both be zero. 'kx' and 'ky' must be > 0.");
//     }

//     // --- 1. Calculate Triangular CSL Parameters ---
//     int sigma_int = m * m + n * n + m * n;
//     double sigma = static_cast<double>(sigma_int);
    
//     double rotation_angle = std::atan2(static_cast<double>(n) * sqrt3, static_cast<double>(2 * m + n));
//     double rotation_angle_deg = rotation_angle * 180.0 / M_PI;

//     double cos_theta = std::cos(rotation_angle);
//     double sin_theta = std::sin(rotation_angle);
    
//     // --- 2. Define UNROTATED Triangular Basis Vectors ---
//     // v1 = a(1, 0)
//     // v2 = a(1/2, sqrt(3)/2)
//     Point2D v1_unrot(a, 0.0);
//     Point2D v2_unrot(a * 0.5, a * sqrt3 * 0.5);

//     // --- 3. Define ROTATED Triangular Basis Vectors ---
//     Point2D v1_rot(
//         v1_unrot.coord.x() * cos_theta - v1_unrot.coord.y() * sin_theta,
//         v1_unrot.coord.x() * sin_theta + v1_unrot.coord.y() * cos_theta
//     );
//     Point2D v2_rot(
//         v2_unrot.coord.x() * cos_theta - v2_unrot.coord.y() * sin_theta,
//         v2_unrot.coord.x() * sin_theta + v2_unrot.coord.y() * cos_theta
//     );

//     // --- 4. Define the Periodic Box Dimensions ---
//     // The rectangular CSL unit cell for a triangular lattice has
//     // dimensions L_x = a*sqrt(Sigma) and L_y = a*sqrt(3*Sigma)
//     double L_x_unit = a * std::sqrt(sigma);
//     double L_y_unit = a * std::sqrt(3.0 * sigma);
    
//     double size_x = static_cast<double>(kx) * L_x_unit; // Final Lx
//     double size_y = static_cast<double>(ky) * L_y_unit; // Final Ly

//     std::cout << "--- Periodic Rotated TRIANGULAR Lattice Generation ---" << std::endl;
//     std::cout << "Parameters: m=" << m << ", n=" << n << " -> Sigma = " << sigma_int << std::endl;
//     std::cout << "Rotation Angle: " << rotation_angle_deg << " degrees" << std::endl;
//     std::cout << "Unit Cell (Lx, Ly): " << L_x_unit << " x " << L_y_unit << std::endl;
//     std::cout << "Final Domain Size: " << size_x << " x " << size_y << std::endl;

//     // ======================================================================
//     // --- STEP 5: Find Basis Atoms using an Integer Grid ---
//     // ======================================================================
    
//     std::set<std::pair<long long, long long>> basis_atoms_int;
//     const double precision = 1e9;
    
//     long long L_int_x = static_cast<long long>(std::round(L_x_unit * precision));
//     long long L_int_y = static_cast<long long>(std::round(L_y_unit * precision));

//     // Search range needs to be larger for triangular lattices
//     int search_range = std::abs(m) + std::abs(n) + 4; 
//     // Heuristic search range, may need adjustment for complex CSLs
//     if (sigma_int > 10) search_range *= 2; 

//     for (int i = -search_range; i <= search_range; ++i) {
//         for (int j = -search_range; j <= search_range; ++j) {
            
//             // Atom position in the ROTATED lattice
//             double x = i * v1_rot.coord.x() + j * v2_rot.coord.x();
//             double y = i * v1_rot.coord.y() + j * v2_rot.coord.y();

//             // Fold the atom back into the [0, L_x) x [0, L_y) CSL cell
//             double x_pbc = std::fmod(x, L_x_unit);
//             double y_pbc = std::fmod(y, L_y_unit);

//             if (x_pbc < 0) x_pbc += L_x_unit;
//             if (y_pbc < 0) y_pbc += L_y_unit;

//             long long x_int = static_cast<long long>(std::round(x_pbc * precision));
//             long long y_int = static_cast<long long>(std::round(y_pbc * precision));

//             if (x_int == L_int_x) x_int = 0;
//             if (y_int == L_int_y) y_int = 0;
            
//             basis_atoms_int.insert({x_int, y_int});
//         }
//     }

//     // This rectangular cell should contain 2*Sigma atoms
//     int expected_atoms = 2 * sigma_int;
//     if (basis_atoms_int.size() != expected_atoms) {
//          std::cout << "\n*** Warning: Found " << basis_atoms_int.size() << " basis atoms, "
//                    << "but expected 2*Sigma=" << expected_atoms << ". " 
//                    << "(This is normal if m,n are not coprime)." << std::endl;
//     } else {
//          std::cout << "Found " << basis_atoms_int.size() << " basis atoms (2*Sigma)." << std::endl;
//     }

//     // ======================================================================
//     // --- STEP 6: Tile the Basis Atoms ---
//     // ======================================================================
    
//     for (int ix = 0; ix < kx; ++ix) {
//         for (int iy = 0; iy < ky; ++iy) {
//             double offset_x = ix * L_x_unit;
//             double offset_y = iy * L_y_unit;

//             for (const auto& p_int : basis_atoms_int) {
//                 double x_basis = static_cast<double>(p_int.first) / precision;
//                 double y_basis = static_cast<double>(p_int.second) / precision;
                
//                 final_points.push_back(
//                     Point2D(x_basis + offset_x, y_basis + offset_y)
//                 );
//             }
//         }
//     }

//     std::cout << "Generated " << final_points.size() << " total points." << std::endl;
//     std::cout << "========================================\n" << std::endl;

//     return final_points;
// }



// --- Example main() function to demonstrate usage ---
// (No changes needed here)
