#pragma once
#ifndef LATTICE_GENERATOR_H
#define LATTICE_GENERATOR_H

#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include "../geometry/Point2D.h"
#include "../geometry/DomainDimensions.h"

class LatticeGenerator {
public:
    // Generate 2D lattice points (square or triangular)
    static std::vector<Point2D> generate_2d_lattice(
        int nx, int ny, double lattice_constant, std::string lattice_type);
    
    // Create periodic copies of the lattice points
    static std::vector<Point2D> create_periodic_copies(
        const std::vector<Point2D>& original_points,
        const DomainDimensions& domain_dims,
        const std::array<double, 2>& offsets,
        const Eigen::Matrix2d& F_ext = Eigen::Matrix2d::Identity());
};

#endif // LATTICE_GENERATOR_H
