#pragma once
#ifndef LATTICE_OPTIMIZER_H
#define LATTICE_OPTIMIZER_H

#include <vector>
#include <functional>
#include <type_traits> // For std::is_same_v
#include <Eigen/Dense>
#include <omp.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <string>

#include "../geometry/Point2D.h"
#include "../mesh/ElementTriangle2D.h"
#include "../lattice_energy/TriangularLatticeCalculator.h"
#include "../lattice_energy/SquareLatticeCalculator.h"
#include "../lattice_energy/Strain_Energy_LatticeCalculator.h"
#include "../include/reductions/LagrangeReduction.h"


#include "src/optimization.h" // This path is relative to ALGLIB_DIR
// UserData structure for optimization
struct UserData {
    std::vector<Point2D>& points;
    std::vector<ElementTriangle2D>& elements;
    Strain_Energy_LatticeCalculator& calculator;
    std::function<double(double)>& energy_function;
    std::function<double(double)>& derivative_function;
    double zero_energy;
    double ideal_lattice_parameter;
    const Eigen::Matrix2d& F_external;
    const std::vector<std::pair<int, int>>& interior_mapping;
    const std::vector<std::pair<int, int>>& full_mapping;
    const std::vector<size_t>& active_elements;  // Added active_elements

    UserData(std::vector<Point2D>& pts,
             std::vector<ElementTriangle2D>& elems,
             Strain_Energy_LatticeCalculator& calc,
             std::function<double(double)>& energy_func,
             std::function<double(double)>& deriv_func,
             double zero_e,
             double ideal_lat_param,
             const Eigen::Matrix2d& F_ext,
             const std::vector<std::pair<int, int>>& int_map,
             const std::vector<std::pair<int, int>>& full_map,
             const std::vector<size_t>& active_elems)  // Added parameter
        : points(pts), elements(elems), calculator(calc),
          energy_function(energy_func), derivative_function(deriv_func),
          zero_energy(zero_e), ideal_lattice_parameter(ideal_lat_param),
          F_external(F_ext), interior_mapping(int_map), full_mapping(full_map),
          active_elements(active_elems) {}  // Initialize active_elements
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

// Map points to solver array in 2D
template <typename PointType>
void map_points_to_solver_array(
    alglib::real_1d_array& x,
    const std::vector<PointType>& all_points,
    const std::vector<std::pair<int, int>>& interior_mapping,
    int n_vars)
{
    // Safety checks
    if (all_points.empty()) {
        std::cerr << "[ERROR] Points vector is empty in map_points_to_solver_array" << std::endl;
        return;
    }
    
    if (interior_mapping.empty()) {
        std::cerr << "[ERROR] Mapping vector is empty in map_points_to_solver_array" << std::endl;
        return;
    }
    
    if (x.length() < 2 * n_vars) {
        std::cerr << "[ERROR] Array too small: x.length()=" << x.length() 
                 << ", needed=" << 2 * n_vars << std::endl;
        return;
    }
    
    // std::cout << "[DEBUG] Starting map_points_to_solver_array with " << interior_mapping.size() 
    //           << " mappings, n_vars=" << n_vars << ", x.length()=" << x.length() << std::endl;
    
    // Map coordinates to the solver array
    for (const auto& pair : interior_mapping) {
        int original_idx = pair.first;
        int solver_idx   = pair.second;
        
        // Check point index validity
        if (original_idx < 0 || original_idx >= all_points.size()) {
            std::cerr << "[ERROR] Point index out of bounds: original_idx=" << original_idx 
                     << ", all_points.size()=" << all_points.size() << std::endl;
            exit(0);
            
        }
        
        // Check array index validity for both coordinates
        if (solver_idx < 0 || solver_idx >= n_vars || solver_idx + n_vars >= x.length()) {
            std::cerr << "[ERROR] Array index out of bounds: solver_idx=" << solver_idx 
                     << ", solver_idx + n_vars=" << (solver_idx + n_vars)
                     << ", x.length()=" << x.length() << std::endl;
            exit(0);
        }
        
        // Handle different point types through template specialization
        if constexpr (std::is_same_v<PointType, Point2D>) {
            // For Point2D with Eigen::Vector2d coord
            x[solver_idx]          = all_points[original_idx].coord(0);
            x[solver_idx + n_vars] = all_points[original_idx].coord(1);
        }
        else if constexpr (std::is_same_v<PointType, Eigen::Vector2d>) {
            // For Eigen::Vector2d
            x[solver_idx]          = all_points[original_idx](0);
            x[solver_idx + n_vars] = all_points[original_idx](1);
        }
        // Add more point types as needed
    }
    
    //std::cout << "[DEBUG] Completed map_points_to_solver_array successfully" << std::endl;
}

    bool createFolder(std::string path);

    // Save configuration
    void saveConfigurationToXY(
        const std::vector<Point2D>& points, 
        int iteration);

    // Function to initialize and get active elements
    const std::vector<size_t>& initialize_active_elements(
        const std::vector<ElementTriangle2D>& elements,
        const std::vector<std::pair<int, int>>& full_mapping,
        int num_points);

    
    // Main optimization function
    void minimize_energy_with_triangles(
        const alglib::real_1d_array &x, 
        double &func, 
        alglib::real_1d_array &grad, 
        void *ptr);

    // Add this declaration to LatticeOptimizer.h
    void deform_boundary_nodes(
        std::vector<Point2D>& points,
        const std::vector<std::pair<int, int>>& dof_mapping,
        const Eigen::Matrix2d& F_ext);

 double find_optimal_lattice_parameter(const std::function<double(double)>& potential, std::string& lattice_type);


#endif // LATTICE_OPTIMIZER_H
