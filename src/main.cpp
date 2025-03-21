#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <tuple>

#include "../include/geometry/Point2D.h"
#include "../include/geometry/DomainDimensions.h"
#include "../include/geometry/LatticeGenerator.h"
#include "../include/mesh/Triangle.h"
#include "../include/mesh/ElementTriangle2D.h"
#include "../include/mesh/MeshGenerator.h"
#include "../include/TriangularLatticeCalculator.h"

int main() {
    // Parameters for lattice
    int nx = 3;
    int ny = 3;
    double lattice_constant = 0.9545707000000000;
    
    std::string lattice_type = "square"; // "square" or "triangular"
    std::vector<Point2D> square_points = LatticeGenerator::generate_2d_lattice(
        nx, ny, lattice_constant, lattice_type);
    
    int original_domain_size = square_points.size();
    
    const std::array<double, 2> offsets = {0.0, 0.0};
    DomainDimensions domain_dims(nx, ny);

    // Generate periodic copies
    std::vector<Point2D> square_points_periodic = LatticeGenerator::create_periodic_copies(
        square_points,
        domain_dims,
        offsets);
    
    // Create Delaunay triangulation
    std::vector<Triangle> triangulation = MeshGenerator::createTrianglesFromPoints(
        square_points_periodic);

    // Create domain maps
    auto [original_domain_map, translation_map] = MeshGenerator::create_domain_maps(
        original_domain_size, domain_dims, offsets);

    // Select unique connected triangles
    std::vector<Triangle> unique_triangles = MeshGenerator::select_unique_connected_triangles(
        square_points_periodic,
        triangulation,
        original_domain_map,
        square_points.size(), // Original domain size
        1e-6 // Minimum Jacobian threshold
    );

    // Create element triangles
    std::vector<ElementTriangle2D> elements = MeshGenerator::createElementTri2D(
        unique_triangles,
        square_points_periodic,
        original_domain_map,
        translation_map
    );
    
    std::cout << "Created " << elements.size() << " element triangles" << std::endl;

    return 0;
}
