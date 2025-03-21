#pragma once
#ifndef LATTICE_OPTIMIZER_H
#define LATTICE_OPTIMIZER_H

#include <vector>
#include <functional>
#include <Eigen/Dense>
#include <omp.h>

#include "../geometry/Point2D.h"
#include "../mesh/ElementTriangle2D.h"
#include "../TriangularLatticeCalculator.h"

#include "src/optimization.h" // This path is relative to ALGLIB_DIR
// UserData structure for optimization
struct UserData {
    std::vector<Point2D>& points;
    std::vector<ElementTriangle2D>& elements;
    TriangularLatticeCalculator& calculator;
    std::function<double(double)>& energy_function;
    std::function<double(double)>& derivative_function;
    double zero_energy;
    double ideal_lattice_parameter;
    Eigen::Matrix2d F_external;
    const std::vector<std::pair<int, int>>& interior_mapping;
    const std::vector<std::pair<int, int>>& full_mapping;
    
    UserData(std::vector<Point2D>& pts, 
            std::vector<ElementTriangle2D>& elems,
            TriangularLatticeCalculator& calc,
            std::function<double(double)>& energy_func,
            std::function<double(double)>& deriv_func,
            double zero_e,
            double ideal_lat_param,
            const Eigen::Matrix2d& F_ext,
            const std::vector<std::pair<int, int>>& int_map,
            const std::vector<std::pair<int, int>>& full_map)
        : points(pts), elements(elems), calculator(calc),
          energy_function(energy_func), derivative_function(deriv_func),
          zero_energy(zero_e), ideal_lattice_parameter(ideal_lat_param),
          F_external(F_ext), interior_mapping(int_map), full_mapping(full_map) {}
};

// DOF mapping function
std::pair<std::vector<std::pair<int, int>>, std::vector<std::pair<int, int>>> 
create_dof_mapping_original(
    const std::vector<Point2D>& points,
    double boundary_tolerance, 
    int pbc);

// Helper functions
void map_solver_array_to_points(
    const alglib::real_1d_array &x, 
    std::vector<Point2D>& points, 
    const std::vector<std::pair<int, int>>& mapping, 
    int n_vars);

void map_points_to_solver_array(
    alglib::real_1d_array &grad, 
    const std::vector<Eigen::Vector2d>& forces,
    const std::vector<std::pair<int, int>>& mapping, 
    int n_vars);

// Save configuration
void saveConfigurationToXY(
    const std::vector<Point2D>& points, 
    int iteration, 
    double energy);

// Main optimization function
void minimize_energy_with_triangles(
    const alglib::real_1d_array &x, 
    double &func, 
    alglib::real_1d_array &grad, 
    void *ptr);

#endif // LATTICE_OPTIMIZER_H
