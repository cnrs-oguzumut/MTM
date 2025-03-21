#include "../include/optimization/LatticeOptimizer.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <limits>

std::pair<std::vector<std::pair<int, int>>, std::vector<std::pair<int, int>>> 
create_dof_mapping_original(
    const std::vector<Point2D>& points,
    double boundary_tolerance, int pbc)
{
    // Interior points mapping (original_idx, solver_idx)
    std::vector<std::pair<int, int>> interior_mapping;
    // Full mapping including boundary points (original_idx, solver_idx or -1)
    std::vector<std::pair<int, int>> full_mapping;
    
    if (points.empty()) {
        return {interior_mapping, full_mapping};
    }
    
    // 1. Find system bounds (min and max coordinates)
    double min_x = std::numeric_limits<double>::max();
    double min_y = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::lowest();
    double max_y = std::numeric_limits<double>::lowest();
    
    for (const auto& p : points) {
        min_x = std::min(min_x, p.coord.x());
        min_y = std::min(min_y, p.coord.y());
        max_x = std::max(max_x, p.coord.x());
        max_y = std::max(max_y, p.coord.y());
    }
    
    // Pre-allocate space
    interior_mapping.reserve(points.size() * 0.75); // Estimate 75% interior points
    full_mapping.reserve(points.size());
    
    // 2. Create both mappings
    int solver_index = 0;
    for (size_t i = 0; i < points.size(); i++) {
        const auto& p = points[i];
        
        // Check if point is on or near the boundary
        bool is_boundary =
            std::abs(p.coord.x() - min_x) <= boundary_tolerance ||
            std::abs(p.coord.x() - max_x) <= boundary_tolerance ||
            std::abs(p.coord.y() - min_y) <= boundary_tolerance ||
            std::abs(p.coord.y() - max_y) <= boundary_tolerance;
            
        if(pbc==1)
            is_boundary = false;
            
        if (is_boundary) {
            // Boundary node - only add to full_mapping with -1
            full_mapping.push_back(std::make_pair(i, -1));
        } else {
            // Interior node - add to both mappings
            interior_mapping.push_back(std::make_pair(i, solver_index));
            full_mapping.push_back(std::make_pair(i, solver_index));
            solver_index++;
        }
    }
    
    // Optimize memory usage
    interior_mapping.shrink_to_fit();
    full_mapping.shrink_to_fit();
    
    return {interior_mapping, full_mapping};
}

// Map solver array to points in 2D
void map_solver_array_to_points(
    const alglib::real_1d_array &x, 
    std::vector<Point2D>& points, 
    const std::vector<std::pair<int, int>>& mapping, 
    int n_vars) 
{
    for (const auto& map_pair : mapping) {
        int point_idx = map_pair.first;
        int var_idx = map_pair.second;
        
        if (var_idx >= 0) {
            points[point_idx].coord(0) = x[var_idx];
            points[point_idx].coord(1) = x[var_idx + n_vars];
        }
    }
}

// Map points to solver array in 2D
void map_points_to_solver_array(
    alglib::real_1d_array &grad, 
    const std::vector<Eigen::Vector2d>& forces,
    const std::vector<std::pair<int, int>>& mapping, 
    int n_vars) 
{
    // Set gradient to zero
    for (int i = 0; i < grad.length(); i++) {
        grad[i] = 0.0;
    }
    
    // Map forces to gradient
    for (const auto& map_pair : mapping) {
        int point_idx = map_pair.first;
        int var_idx = map_pair.second;
        
        if (var_idx >= 0) {
            grad[var_idx] = -forces[point_idx](0);
            grad[var_idx + n_vars] = -forces[point_idx](1);
        }
    }
}

// Save configuration to XY file (2D version)
void saveConfigurationToXY(
    const std::vector<Point2D>& points, 
    int iteration, 
    double energy) 
{
    std::ostringstream filename;
    filename << "config_" << std::setfill('0') << std::setw(4) << iteration << ".xy";
    
    std::ofstream file(filename.str());
    file << points.size() << std::endl;
    file << "# Energy: " << energy << std::endl;
    
    for (const auto& point : points) {
        file << "P " << point.coord.x() << " " << point.coord.y() << std::endl;
    }
    
    file.close();
}

// Main optimization function
void minimize_energy_with_triangles(
    const alglib::real_1d_array &x, 
    double &func, 
    alglib::real_1d_array &grad, 
    void *ptr) 
{
    // Get user data from pointer
    auto* userData = reinterpret_cast<UserData*>(ptr);

    std::vector<Point2D>& points = userData->points;
    std::vector<ElementTriangle2D>& elements = userData->elements;
    TriangularLatticeCalculator& calculator = userData->calculator;
    std::function<double(double)>& energy_function = userData->energy_function;
    std::function<double(double)>& derivative_function = userData->derivative_function;
    double zero_energy = userData->zero_energy;
    double ideal_lattice_parameter = userData->ideal_lattice_parameter;
    const Eigen::Matrix2d F_external = userData->F_external;
    // Access the mapping
    const std::vector<std::pair<int, int>>& interior_mapping = userData->interior_mapping;
    const std::vector<std::pair<int, int>>& full_mapping = userData->full_mapping;

    // For triangular lattice in 2D, normalization is different
    static double normalisation = pow(ideal_lattice_parameter, 2.0) * sqrt(3.0) / 4.0;
    
    int n_points = points.size();
    int n_vars = x.length() / 2; // 2D has only x and y coordinates
    static int iteration = 0;
    
    // Reset function value
    func = 0.0;
    
    // Clear forces before recalculation
    std::vector<Eigen::Vector2d> global_forces(n_points, Eigen::Vector2d::Zero());
    
    // Set positions from optimization variables
    map_solver_array_to_points(x, points, interior_mapping, n_vars);
    
    // Reset energy
    double total_energy = 0.0;
    
    // Reset global forces
    for (auto& force : global_forces) {
        force.setZero();
    }

    // Thread-local storage for forces to avoid race conditions
    std::vector<std::vector<Eigen::Vector2d>> thread_local_forces;

    // Static variables that persist between function calls
    static std::vector<bool> is_dof;
    static std::vector<size_t> active_elements;
    static bool is_initialized = false;
    
    // Initialize only on first call
    if (!is_initialized) {
        // Create DOF lookup array
        is_dof.resize(points.size(), false);
        for (const auto& mapping : full_mapping) {
            if (mapping.second != -1) {
                is_dof[mapping.first] = true;
            }
        }
        
        // Pre-filter elements that have at least one DOF node
        active_elements.clear();
        active_elements.reserve(elements.size());
        
        for (size_t i = 0; i < elements.size(); i++) {
            const auto& element = elements[i];
            bool hasActiveDof = false;
            
            for (int j = 0; j < 3; j++) { // 3 nodes in a triangle
                int nodeIndex = element.getNodeIndex(j);
                if (is_dof[nodeIndex]) {
                    hasActiveDof = true;
                    break;
                }
            }
            
            if (hasActiveDof) {
                active_elements.push_back(i);
            }
        }
        
        is_initialized = true;
        std::cout << "Active elements: " << active_elements.size() << " out of " 
                  << elements.size() << " (" 
                  << (100.0 * active_elements.size() / elements.size()) << "%)" << std::endl;
    }  

    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        int num_threads = omp_get_num_threads();
        
        #pragma omp single
        {
            thread_local_forces.resize(num_threads, 
                                    std::vector<Eigen::Vector2d>(global_forces.size(), 
                                                                Eigen::Vector2d::Zero()));
        }
        
        #pragma omp for reduction(+:total_energy)
        for (size_t idx = 0; idx < active_elements.size(); idx++) {
            size_t i = active_elements[idx];
            auto& element = elements[i];
            
            element.setExternalDeformation(F_external);
            element.calculate_deformation_gradient(points);
            Eigen::Matrix2d F = element.getDeformationGradient();
            
            Eigen::Matrix2d C = element.getMetricTensor();
            double element_energy = calculator.calculate_energy(C, energy_function, zero_energy)/normalisation;
            total_energy += element_energy * element.getArea(); // Use area instead of volume
            
            Eigen::Matrix2d dE_dC = calculator.calculate_derivative(C, derivative_function)/normalisation;
            Eigen::Matrix2d P = 2.0 * F * dE_dC;
            
            // Add forces to thread-local storage
            element.assemble_forces(P, thread_local_forces[thread_id]);
        }
    }

    // Combine thread-local forces
    for (size_t i = 0; i < global_forces.size(); i++) {
        for (size_t t = 0; t < thread_local_forces.size(); t++) {
            global_forces[i] += thread_local_forces[t][i];
        }
    }
    
    // Set function value to total energy
    func = total_energy;
    if (iteration % 2 == 0) {
        std::cout << "\rEnergy: " << func << " at iteration: " << iteration << "    " << std::flush;
    }
    
    // Map forces to gradient array
    map_points_to_solver_array(grad, global_forces, interior_mapping, n_vars);
    
    if (iteration % 2 == 0) {
        saveConfigurationToXY(points, iteration, func);
        //saveMetricTensors(elements, iteration);
    }
    
    iteration++;
}
