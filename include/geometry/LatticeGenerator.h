#pragma once
#ifndef LATTICE_GENERATOR_H
#define LATTICE_GENERATOR_H

#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <utility>
#include <array>
#include <Eigen/Dense>
#include "../geometry/Point2D.h"
#include "../geometry/DomainDimensions.h"

struct CSLAngleInfo {
    int m = 1;
    int n = 0;
    int sigma = 1;
    double theta_rad = 0.0;
    double theta_deg = 0.0;
    double Lx_over_a = 1.0;
    double Ly_over_a = 1.0;
    int atoms_per_cell = 1;
    bool is_reduced = false;
};

struct CSLLatticeResult {
    std::vector<Point2D> points;
    double Lx = 0.0;
    double Ly = 0.0;
    int m = 1;
    int n = 0;
    int kx = 1;
    int ky = 1;
    int sigma = 1;
    double theta_deg = 0.0;
    int atoms_cell = 0;
    int total_atoms = 0;
    bool is_reduced = false;
};

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
    
    static std::vector<Point2D> generate_2d_lattice_rotated(
        int nx, int ny, double lattice_constant, std::string lattice_type, double rotation_angle_degrees);

    static std::pair<std::vector<Point2D>, Point2D> generate_periodic_rotated_lattice_v2(
        int m, int n, int kx, int ky,
        double lattice_constant, std::string lattice_type);

    // Modern CSL (Coincident Site Lattice) Recipe (Paper Oct 19, 2026)
    static std::vector<CSLAngleInfo> get_allowed_csl_angles(
        const std::string& lattice_type, int max_sigma = 80);

    static CSLAngleInfo find_closest_csl_angle(
        const std::string& lattice_type, double target_angle_deg, int max_sigma = 80);

    static CSLLatticeResult generate_csl_lattice(
        const std::string& lattice_type, int m, int n, int kx, int ky, double lattice_constant = 1.0);

    // Fixed boundary condition rotated lattice (any arbitrary angle theta in degrees)
    static std::pair<std::vector<Point2D>, DomainDimensions> generate_rotated_fixed_lattice(
        const std::string& lattice_type, int nx, int ny, double lattice_constant, double rotation_angle_degrees);
};

#endif // LATTICE_GENERATOR_H
