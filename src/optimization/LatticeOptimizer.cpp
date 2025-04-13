#include "../include/optimization/LatticeOptimizer.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <limits>
//Fixed boundary conditions on the bo
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


size_t findMiddleAtom(const std::vector<Point2D>& , bool  ) ;

std::pair<std::vector<std::pair<int, int>>, std::vector<std::pair<int, int>>>
create_dof_mapping_with_radius(
    const std::vector<Point2D>& points,
    double radius,
    int pbc
) {
    // Interior points mapping (original_idx, solver_idx)
    std::vector<std::pair<int, int>> interior_mapping;
    // Full mapping including boundary points (original_idx, solver_idx or -1)
    std::vector<std::pair<int, int>> full_mapping;
    
    if (points.empty()) {
        return {interior_mapping, full_mapping};
    }
    
    // Find the middle atom index
    size_t middle_index = findMiddleAtom(points, true);
    
    // Use middle atom's coordinates as the center
    double center_x = points[middle_index].coord.x();
    double center_y = points[middle_index].coord.y();
    
    std::cout << "System center (middle atom): (" << center_x << ", " << center_y << ")" << std::endl;
    std::cout << "Using radius: " << radius << " for fixed boundary" << std::endl;
    
    // Pre-allocate space
    interior_mapping.reserve(points.size() * 0.75); // Estimate 75% interior points
    full_mapping.reserve(points.size());
    
    // Create both mappings
    int solver_index = 0;
    for (size_t i = 0; i < points.size(); i++) {
        const auto& p = points[i];
        
        // Calculate distance from middle atom
        double dx = p.coord.x() - center_x;
        double dy = p.coord.y() - center_y;
        double distance_from_center = std::sqrt(dx*dx + dy*dy);
        
        // Check if point is outside the radius
        bool is_outside_radius = distance_from_center > radius;
        
        if (pbc == 1) {
            is_outside_radius = false; // Override for periodic boundary conditions
        }
        
        if (is_outside_radius) {
            // Outside radius node - only add to full_mapping with -1
            full_mapping.push_back(std::make_pair(i, -1));
        } else {
            // Inside radius node - add to both mappings
            interior_mapping.push_back(std::make_pair(i, solver_index));
            full_mapping.push_back(std::make_pair(i, solver_index));
            solver_index++;
        }
    }
    
    // Optimize memory usage
    interior_mapping.shrink_to_fit();
    full_mapping.shrink_to_fit();
    
    std::cout << "Total nodes: " << points.size() << std::endl;
    std::cout << "Interior nodes: " << interior_mapping.size() << " ("
              << (100.0 * interior_mapping.size() / points.size()) << "%)" << std::endl;
    std::cout << "Fixed boundary nodes: " << (points.size() - interior_mapping.size()) << " ("
              << (100.0 * (points.size() - interior_mapping.size()) / points.size()) << "%)" << std::endl;
    
    return {interior_mapping, full_mapping};
}




std::pair<std::vector<std::pair<int, int>>, std::vector<std::pair<int, int>>> 
create_dof_mapping_with_boundaries(
    const std::vector<Point2D>& points,
    std::vector<ElementTriangle2D>& element_triangles,
    const std::vector<int>& fixed_atoms_from_indenter,
    const std::vector<int>& boundary_fixed_nodes  // New parameter
) {
    // Interior points mapping (original_idx, solver_idx)
    std::vector<std::pair<int, int>> interior_mapping;
    
    // Full mapping including boundary points (original_idx, solver_idx or -1)
    std::vector<std::pair<int, int>> full_mapping;
    
    if (points.empty()) {
        return {interior_mapping, full_mapping};
    }
    

    
    // Pre-allocate space
    interior_mapping.reserve(points.size());
    full_mapping.reserve(points.size());
    
    // Combine fixed atoms from boundary and indenter
    std::set<int> fixed_atom_set;
    fixed_atom_set.insert(boundary_fixed_nodes.begin(), boundary_fixed_nodes.end());
    fixed_atom_set.insert(fixed_atoms_from_indenter.begin(), fixed_atoms_from_indenter.end());
    
    // Create both mappings
    int solver_index = 0;
    int fixed_boundary_count = 0;
    int fixed_indenter_count = 0;
    
    for (size_t i = 0; i < points.size(); i++) {
        if (fixed_atom_set.count(i) > 0) {
            // Fixed node (either boundary or indenter)
            full_mapping.push_back(std::make_pair(i, -1));
            
            // Distinguish between boundary and indenter fixed nodes
            if (std::find(boundary_fixed_nodes.begin(), boundary_fixed_nodes.end(), i) 
                != boundary_fixed_nodes.end()) {
                fixed_boundary_count++;
            } else {
                fixed_indenter_count++;
            }
        } else {
            // Interior node - free
            interior_mapping.push_back(std::make_pair(i, solver_index));
            full_mapping.push_back(std::make_pair(i, solver_index));
            solver_index++;
        }
    }
    
    // Optimize memory usage
    interior_mapping.shrink_to_fit();
    full_mapping.shrink_to_fit();
    
    std::cout << "Total nodes: " << points.size() << std::endl;
    std::cout << "Free nodes: " << interior_mapping.size() << " ("
              << (100.0 * interior_mapping.size() / points.size()) << "%)" << std::endl;
    std::cout << "Fixed boundary nodes: " << fixed_boundary_count << " ("
              << (100.0 * fixed_boundary_count / points.size()) << "%)" << std::endl;
    std::cout << "Fixed indenter nodes: " << fixed_indenter_count << " ("
              << (100.0 * fixed_indenter_count / points.size()) << "%)" << std::endl;
   
    std::unordered_map<int, int> boundary_node_distribution;

    for (size_t tri_idx = 0; tri_idx < element_triangles.size(); tri_idx++) {
        int boundary_count = 0;
        
        for (int node_idx = 0; node_idx < 3; node_idx++) {
            int global_node_idx = element_triangles[tri_idx].getNodeIndex(node_idx);
            
            // Use full_mapping to check if the node is fixed (boundary or indenter)
            if (full_mapping[global_node_idx].second == -1) {
                boundary_count++;
            }
        }
        
        // Set the number of boundary nodes for this triangle
        element_triangles[tri_idx].setBoundaryNodeNumber(boundary_count);
        
        // Track distribution
        boundary_node_distribution[boundary_count]++;
    }

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
    size_t num_points) 
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

void minimize_energy_with_triangles(
    const alglib::real_1d_array &x, 
    double &func, 
    alglib::real_1d_array &grad, 
    void *ptr) 
{
    static int iteration = 0;
    static std::vector<std::vector<Eigen::Vector2d>> thread_local_forces;

    // if (ptr == nullptr) {
    //     std::cerr << "[ERROR] User data pointer is NULL!" << std::endl;
    //     return;
    // }
    
    auto* userData = reinterpret_cast<UserData*>(ptr);
    std::vector<ElementTriangle2D>& elements = userData->elements;
    const auto& interior_mapping = userData->interior_mapping;
    const int n_vars = interior_mapping.size();
    const int n_points = userData->points.size();
    //const double normalisation = pow(userData->ideal_lattice_parameter, 2.0);
    //const double normalisation = 1.;// (sqrt(3.)/2.)*pow(userData->ideal_lattice_parameter, 2.0);

    const double normalisation = userData->calculator.getUnitCellArea();
    // Initialize thread storage (forces for all points, including boundaries)
    #pragma omp parallel
    {
        const int thread_id = omp_get_thread_num();
        #pragma omp single
        {
            thread_local_forces.resize(omp_get_num_threads(), 
                                     std::vector<Eigen::Vector2d>(n_points, Eigen::Vector2d::Zero()));
        }
        std::fill(thread_local_forces[thread_id].begin(), 
                 thread_local_forces[thread_id].end(), 
                 Eigen::Vector2d::Zero());
    }

    // Initialize gradient
    for(int i = 0; i < grad.length(); i++) {
        grad[i] = 0.0;
    }

    // Main computation
    double total_energy = 0.0;
    #pragma omp parallel reduction(+:total_energy)
    {
        const int thread_id = omp_get_thread_num();
        auto& my_forces = thread_local_forces[thread_id];
        
        #pragma omp for schedule(guided)
        for (size_t idx = 0; idx < userData->active_elements.size(); idx++) {
            ElementTriangle2D& element = elements[userData->active_elements[idx]];
            
            //element.setExternalDeformation(userData->F_external);
            element.calculate_deformation_gradient(x);
            
            const Eigen::Matrix2d F = element.getDeformationGradient();
            const Eigen::Matrix2d C =  element.getMetricTensor(); //F.transpose() * F;
            
            const auto result = lagrange::reduce(C);
            const double element_energy = userData->calculator.calculate_energy(
                result.C_reduced, userData->energy_function, userData->zero_energy) / normalisation;
            
            total_energy += element_energy * element.getReferenceArea();
            
            const Eigen::Matrix2d dE_dC = userData->calculator.calculate_derivative(
                result.C_reduced, userData->derivative_function) / normalisation;
            const Eigen::Matrix2d P = 2.0 * F * result.m_matrix * dE_dC * result.m_matrix.transpose();
            
            element.assemble_forces(P, my_forces);
        }
    }

    // Directly accumulate forces into grad using mapping
    #pragma omp parallel for
    for (int i = 0; i < n_vars; i++) {
        const auto& [point_idx, dof_idx] = interior_mapping[i];
        
        // Sum forces from all threads for this point
        Eigen::Vector2d force_sum = Eigen::Vector2d::Zero();
        for (const auto& forces : thread_local_forces) {
            force_sum += forces[point_idx];
        }
        
        grad[dof_idx] = force_sum.x();
        grad[n_vars + dof_idx] = force_sum.y();
    }

    // Set energy result
    func = total_energy;


    
    iteration++;
}


void minimize_energy_with_triangles_noreduction(
    const alglib::real_1d_array &x, 
    double &func, 
    alglib::real_1d_array &grad, 
    void *ptr) 
{
    static int iteration = 0;
    static std::vector<std::vector<Eigen::Vector2d>> thread_local_forces;

    // if (ptr == nullptr) {
    //     std::cerr << "[ERROR] User data pointer is NULL!" << std::endl;
    //     return;
    // }
    
    auto* userData = reinterpret_cast<UserData*>(ptr);
    //std::vector<Point2D>& points = userData->points;
    std::vector<ElementTriangle2D>& elements = userData->elements;
    const auto& interior_mapping = userData->interior_mapping;
    const int n_vars = interior_mapping.size();
    const int n_points = userData->points.size();
    //const double normalisation = pow(userData->ideal_lattice_parameter, 2.0);
    //const double normalisation = 1.;// (sqrt(3.)/2.)*pow(userData->ideal_lattice_parameter, 2.0);

    const double normalisation =(sqrt(3.)/2.)*pow(userData->ideal_lattice_parameter, 2.0);
    // Set positions from optimization variables
    //map_solver_array_to_points(x, points, interior_mapping, n_vars);
    #pragma omp parallel
    {
        const int thread_id = omp_get_thread_num();
        #pragma omp single
        {
            thread_local_forces.resize(omp_get_num_threads(), 
                                     std::vector<Eigen::Vector2d>(n_points, Eigen::Vector2d::Zero()));
        }
        std::fill(thread_local_forces[thread_id].begin(), 
                 thread_local_forces[thread_id].end(), 
                 Eigen::Vector2d::Zero());
    }

    // Initialize gradient
    for(int i = 0; i < grad.length(); i++) {
        grad[i] = 0.0;
    }

    // Main computation
    double total_energy = 0.0;
    #pragma omp parallel reduction(+:total_energy)
    {
        const int thread_id = omp_get_thread_num();
        auto& my_forces = thread_local_forces[thread_id];
        
        #pragma omp for schedule(guided)
        for (size_t idx = 0; idx < userData->active_elements.size(); idx++) {
            ElementTriangle2D& element = elements[userData->active_elements[idx]];
            
            //element.setExternalDeformation(userData->F_external);
            element.calculate_deformation_gradient(x);
            //element.calculate_deformation_gradient(points);
            
            const Eigen::Matrix2d F = element.getDeformationGradient();
            Eigen::Matrix2d C = F.transpose() * F;
            // std::cout<<"F: "<<F<<std::endl;
            //std::cout<<"C: "<<C<<std::endl;
            
            const double element_energy = userData->calculator.calculate_energy(
                C, userData->energy_function, userData->zero_energy)/ normalisation ;
            
            total_energy += element_energy * element.getReferenceArea();
            
            const Eigen::Matrix2d dE_dC = userData->calculator.calculate_derivative(
                C, userData->derivative_function) / normalisation;
            const Eigen::Matrix2d P = 2.0 * F  * dE_dC ;
            
            element.assemble_forces(P, my_forces);
        }
    }


    // Directly accumulate forces into grad using mapping
    #pragma omp parallel for
    for (int i = 0; i < n_vars; i++) {
        const auto& [point_idx, dof_idx] = interior_mapping[i];
        
        // Sum forces from all threads for this point
        Eigen::Vector2d force_sum = Eigen::Vector2d::Zero();
        for (const auto& forces : thread_local_forces) {
            force_sum += forces[point_idx];
        }
        
        grad[dof_idx] = force_sum.x();
        grad[n_vars + dof_idx] = force_sum.y();
    }

    // Set energy result
    func = total_energy;


    
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


void deform_boundary_nodes_ref(
    std::vector<Point2D>& points,
    const std::vector<Point2D>& points_ref,
    const std::vector<std::pair<int, int>>& dof_mapping,
    const Eigen::Matrix2d& F_ext)
{
    for (const auto& pair : dof_mapping) {
        int original_idx = pair.first;
        int solver_idx = pair.second;
        
        // Only apply deformation to boundary nodes
        if (solver_idx == -1) {
            // Apply deformation: x_deformed = F·x
            points[original_idx].coord = F_ext * points_ref[original_idx].coord;
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