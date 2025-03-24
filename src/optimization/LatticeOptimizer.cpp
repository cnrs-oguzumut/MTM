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

bool createFolder(std::string path)
{
    // Check if folder already exists
    struct stat info;
    if (stat(path.c_str(), &info) == 0 && S_ISDIR(info.st_mode))
    {
        std::cout << "Folder already exists" << std::endl;
        return true;
    }

    // Create the folder
    int status = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (status != 0)
    {
        std::cerr << "Error creating folder" << std::endl;
        return false;
    }

    std::cout << "Folder created successfully" << std::endl;
    return true;
}


// Save configuration to XY file (2D version)
void saveConfigurationToXY(
    const std::vector<Point2D>& points,
    int iteration)
{
    
    createFolder("./configurations");
    
    std::ostringstream filename;
    filename << "./configurations/config_" << std::setfill('0') << std::setw(4) << 10000+iteration << ".xy";
    std::ofstream file(filename.str());
    
    // XY format requires a header line with the number of atoms
    file << points.size() << "\n";
    // Add a comment line with information (Ovito can handle this)
    file << "# Configuration at iteration " << iteration << "\n";
    
    // Write each point with a particle type identifier
    for (const auto& point : points) {
        // Format: "Type X Y"
        // Use "1" as the type instead of "P" for better compatibility
        file << "1 " << point.coord.x() << " " << point.coord.y() << "\n";
    }
    
    file.close();
}

// Function to initialize and get active elements - returns a reference to a static vector
std::vector<size_t> initialize_active_elements(
    const std::vector<ElementTriangle2D>& elements,
    const std::vector<std::pair<int, int>>& full_mapping,
    int num_points) 
{
    // Static variables that persist between function calls
    std::vector<bool> is_dof;
    std::vector<size_t> active_elements;
    bool is_initialized = false;
    
    // === DEBUGGING ===
    std::cout << "[DEBUG] is_initialized: " << (is_initialized ? "true" : "false") << std::endl;
    
    // Initialize only on first call
    if (!is_initialized) {
        // === DEBUGGING ===
        std::cout << "[DEBUG] Initializing static data structures..." << std::endl;
        
        // Create DOF lookup array
        is_dof.resize(num_points, false);
        
        // === DEBUGGING ===
        std::cout << "[DEBUG] Resized is_dof vector to " << is_dof.size() << std::endl;
        
        for (const auto& mapping : full_mapping) {
            if (mapping.second != -1) {
                // === DEBUGGING ===
                if (mapping.first >= is_dof.size()) {
                    std::cerr << "[ERROR] Index out of bounds in is_dof! mapping.first="
                              << mapping.first << ", is_dof.size()=" << is_dof.size() << std::endl;
                    continue;
                }
                
                is_dof[mapping.first] = true;
            }
        }
        
        // === DEBUGGING ===
        std::cout << "[DEBUG] Filtering active elements..." << std::endl;
        
        // Pre-filter elements that have at least one DOF node
        active_elements.clear();
        active_elements.reserve(elements.size());
        
        for (size_t i = 0; i < elements.size(); i++) {
            const auto& element = elements[i];
            bool hasActiveDof = false;
            
            for (int j = 0; j < 3; j++) { // 3 nodes in a triangle
                int nodeIndex = element.getNodeIndex(j);
                
                // === DEBUGGING ===
                if (nodeIndex < 0 || nodeIndex >= is_dof.size()) {
                    std::cerr << "[ERROR] Node index out of bounds! nodeIndex="
                              << nodeIndex << ", is_dof.size()=" << is_dof.size() << std::endl;
                    continue;
                }
                
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
    
    // === DEBUGGING ===
    std::cout << "[DEBUG] Starting computation with " << active_elements.size()
              << " active elements" << std::endl;
              
    // Return a reference to the static vector
    return active_elements;
}
// Main optimization function
// Main optimization function
void minimize_energy_with_triangles(
    const alglib::real_1d_array &x, 
    double &func, 
    alglib::real_1d_array &grad, 
    void *ptr) 
{
    static int iteration = 0;
    static std::ofstream energy_log;
    
    // Open log file on first iteration
    if (iteration == 0) {
        energy_log.open("energy_evolution.csv");
        if (energy_log.is_open()) {
            energy_log << "Iteration,Energy" << std::endl;
        } else {
            std::cerr << "[ERROR] Could not open energy_evolution.csv for writing!" << std::endl;
        }
    }

    // Check if ptr is null
    if (ptr == nullptr) {
        std::cerr << "[ERROR] User data pointer is NULL!" << std::endl;
        return;
    }
    
    // Get user data from pointer
    auto* userData = reinterpret_cast<UserData*>(ptr);
    
    std::vector<Point2D>& points = userData->points;
    std::vector<ElementTriangle2D>& elements = userData->elements;
    Strain_Energy_LatticeCalculator& calculator = userData->calculator;
    std::function<double(double)>& energy_function = userData->energy_function;
    std::function<double(double)>& derivative_function = userData->derivative_function;
    double zero_energy = userData->zero_energy;
    double ideal_lattice_parameter = userData->ideal_lattice_parameter;
    const Eigen::Matrix2d& F_external = userData->F_external;
    const std::vector<std::pair<int, int>>& interior_mapping = userData->interior_mapping;
    const std::vector<std::pair<int, int>>& full_mapping = userData->full_mapping;
    const std::vector<size_t>& active_elements = userData->active_elements;
    bool third_condition_flag = userData->third_condition_flag;
    
    // For square lattice in 2D
    static double normalisation = pow(ideal_lattice_parameter, 2.0);
    
    int n_points = points.size();
    int n_vars = x.length() / 2; // 2D has only x and y coordinates
    
    // Reset function value
    func = 0.0;
    
    // Clear forces before recalculation
    std::vector<Eigen::Vector2d> global_forces(n_points, Eigen::Vector2d::Zero());
    
    // Set positions from optimization variables
    map_solver_array_to_points(x, points, interior_mapping, n_vars);
    
    // Reset energy
    double total_energy = 0.0;
    
    // Thread-local storage for forces to avoid race conditions
    std::vector<std::vector<Eigen::Vector2d>> thread_local_forces;

    // PARALLEL VERSION
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
            
            if (i >= elements.size()) {
                continue;
                std::cout<<"Error: Element index out of bounds"<<std::endl;
                exit(0);
            }
            
            auto& element = elements[i];
            
            element.setExternalDeformation(F_external);
            element.calculate_deformation_gradient(points);
            Eigen::Matrix2d F = element.getDeformationGradient();
            
            Eigen::Matrix2d C = element.getMetricTensor();
            lagrange::Result result = lagrange::reduce(C);
            C = result.C_reduced; 

            
            double element_energy = calculator.calculate_energy(C, energy_function, zero_energy)/normalisation;
            total_energy += element_energy * element.getArea(); // Use area instead of volume
            
            Eigen::Matrix2d dE_dC = calculator.calculate_derivative(C, derivative_function)/normalisation;
            Eigen::Matrix2d P = 2.0 * F * result.m_matrix * dE_dC * result.m_matrix.transpose();

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
    
    // Save energy to log file
    if (energy_log.is_open()) {
        Eigen::Matrix2d C_ext = F_external.transpose() * F_external;
        lagrange::Result result_ext = lagrange::reduce(C_ext);
        C_ext = result_ext.C_reduced; 
        double homogeneous_energy = calculator.calculate_energy(C_ext, energy_function, zero_energy)/normalisation;
        energy_log << iteration << "," << func-0.5*homogeneous_energy*active_elements.size() << std::endl;
    }
    
    // if (iteration % 10 == 0) {
    //     std::cout << "Energy: " << func << " at iteration: " << iteration << std::endl;
    // }
    
    // Map forces to gradient array
    map_points_to_solver_array(grad, global_forces, interior_mapping, n_vars);
    
    if (iteration % 10 == 0) {
        //saveConfigurationToXY(points, iteration);
    }
    
    iteration++;
}

// Add this implementation to LatticeOptimizer.cpp
void deform_boundary_nodes(
    std::vector<Point2D>& points,
    const std::vector<std::pair<int, int>>& dof_mapping,
    const Eigen::Matrix2d& F_ext)
{
    for (const auto& pair : dof_mapping) {
        int original_idx = pair.first;
        int solver_idx = pair.second;
        
        // Only apply deformation to boundary nodes
        if (solver_idx == -1) {
            // Apply deformation: x_deformed = F·x
            points[original_idx].coord = F_ext * points[original_idx].coord;
        }
    }
}


double find_optimal_lattice_parameter(const std::function<double(double)>& potential, std::string& lattice_type) {
    std::vector<double> scales;
    std::vector<double> energies;
    const int num_points = 100000;
    const double scale_min = 1.05;
    const double scale_max = 1.07;
    double min_energy = std::numeric_limits<double>::max();
    double optimal_scale = 0.0;
    Eigen::Matrix2d C;
    C.setIdentity();
    
    // Calculate energy for different lattice parameters
    for (int i = 0; i <= num_points; ++i) {
        double scale = scale_min + (scale_max - scale_min) * i / num_points;
        scales.push_back(scale);
        
        double energy = 0.0;
        
        // Calculate energy using 2D lattice calculator
        if (lattice_type == "triangular") {
            TriangularLatticeCalculator calculator(scale);
            energy = calculator.calculate_energy(C, potential, 0);
        } else if (lattice_type == "square") {
            SquareLatticeCalculator calculator(scale);
            energy = calculator.calculate_energy(C, potential, 0);
        }
        
        energies.push_back(energy);
        
        if (energy < min_energy) {
            min_energy = energy;
            optimal_scale = scale;
        }
        
        if (i % 1000 == 0) {
            std::cout << "Progress: " << (i * 100.0 / num_points) << "%, Current scale: "
                      << scale << ", Energy: " << energy << std::endl;
        }
    }
    
    std::cout << "\n*** OPTIMAL LATTICE PARAMETER FOUND ***" << std::endl;
    std::cout << "Optimal lattice parameter: " << optimal_scale << std::endl;
    std::cout << "Minimum energy: " << min_energy << std::endl;
    
    // Write results to file
    std::ofstream outfile("scale_energy_2d.dat");
    outfile << "# Scale Energy\n";
    for (size_t i = 0; i < scales.size(); ++i) {
        outfile << scales[i] << " " << energies[i] << "\n";
    }
    outfile.close();
    
    // Mark optimal point
    std::ofstream outfile_optimal("optimal_scale_2d.dat");
    outfile_optimal << "# Scale Energy\n";
    outfile_optimal << optimal_scale << " " << min_energy * 1.05 << "\n";
    outfile_optimal << optimal_scale << " " << min_energy * 0.95 << "\n";
    outfile_optimal.close();
    
    std::cout << "Lattice parameter analysis written to scale_energy_2d.dat\n";
    return optimal_scale;
}