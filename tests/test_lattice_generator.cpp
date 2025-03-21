#include <iostream>
#include <cassert>
#include "../include/geometry/LatticeGenerator.h"
#include "../include/geometry/DomainDimensions.h"

int main() {
    // Test square lattice generation
    int nx = 2;
    int ny = 2;
    double lattice_constant = 1.0;
    std::string lattice_type = "square";
    
    std::vector<Point2D> points = LatticeGenerator::generate_2d_lattice(
        nx, ny, lattice_constant, lattice_type);
    
    // Check number of points
    assert(points.size() == nx * ny);
    
    // Check a few point coordinates
    assert(points[0].coord.x() == 0.0);
    assert(points[0].coord.y() == 0.0);
    
    assert(points[1].coord.x() == 1.0);
    assert(points[1].coord.y() == 0.0);
    
    assert(points[2].coord.x() == 0.0);
    assert(points[2].coord.y() == 1.0);
    
    assert(points[3].coord.x() == 1.0);
    assert(points[3].coord.y() == 1.0);
    
    std::cout << "All tests passed!" << std::endl;
    return 0;
}
