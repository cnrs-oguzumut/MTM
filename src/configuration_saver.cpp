#include "../include/output/configuration_saver.h"
#include "../include/reductions/LagrangeReduction.h"

// Function to save configuration with nodal stress and energy values to XYZ file for OVITO (2D version)
void saveConfigurationWithStressAndEnergy2D(
    UserData* userData,
    int iteration, 
    double total_energy) {
    
    // Perform null check
    if (!userData) {
        std::cerr << "Error: userData is null in saveConfigurationWithStressAndEnergy2D" << std::endl;
        return;
    }
    
    // Extract needed values from userData with non-const references
    std::vector<Point2D>& points = userData->points;
    std::vector<ElementTriangle2D>& elements = userData->elements;
    Strain_Energy_LatticeCalculator& calculator = userData->calculator;
    std::function<double(double)>& energy_function = userData->energy_function;
    std::function<double(double)>& derivative_function = userData->derivative_function;
    double zero_energy = userData->zero_energy;
    double ideal_lattice_parameter = userData->ideal_lattice_parameter;
    const Eigen::Matrix2d& F_ext = userData->F_external;
    const std::vector<std::pair<int, int>>& interior_mapping = userData->interior_mapping;
    const std::vector<std::pair<int, int>>& full_mapping = userData->full_mapping;
    const std::vector<size_t>& active_elements = userData->active_elements;
    
    // For square lattice in 2D - match the normalization from minimize_energy_with_triangles
    double normalisation = pow(ideal_lattice_parameter, 2.0);    
    
    // Create directory if it doesn't exist
    std::filesystem::create_directory("trajectory");
    
    // Create filename with iteration number
    std::stringstream filename;
    filename << "trajectory/configuration_full_2d_" << std::setw(5) << std::setfill('0') << iteration << ".xyz";
    
    // Open file for writing
    std::ofstream file(filename.str());
    if (!file) {
        std::cerr << "Error: Could not open file " << filename.str() << " for writing." << std::endl;
        return;
    }
    
    // Calculate nodal stress and energy values
    std::vector<double> nodal_stress(points.size(), 0.0);
    std::vector<double> nodal_energy(points.size(), 0.0);
    std::vector<int> node_count(points.size(), 0);
    
    // DEBUG: Print basic info
    std::cout << "DEBUG: Number of points: " << points.size() << std::endl;
    std::cout << "DEBUG: Number of elements: " << elements.size() << std::endl;
    std::cout << "DEBUG: Number of active elements: " << active_elements.size() << std::endl;
    std::cout << "DEBUG: Normalization value: " << normalisation << std::endl;
    std::cout << "DEBUG: F_external matrix:\n" << F_ext << std::endl;
    
    int processed_elements = 0;
    double recalculated_total_energy = 0.0;
    
    // Loop through active elements to accumulate stress and energy at nodes
    for (size_t elem_idx : active_elements) {
        if (elem_idx >= elements.size()) {
            std::cerr << "Warning: Element index " << elem_idx << " out of range." << std::endl;
            continue;
        }
        
        auto& element = elements[elem_idx]; // Remove const to allow setting external deformation
        
        // Skip if the element isn't initialized
        if (!element.isInitialized()) {
            std::cout << "DEBUG: Element " << elem_idx << " is not initialized." << std::endl;
            continue;
        }
        
        // Set external deformation and recalculate deformation gradient - important!
        element.setExternalDeformation(F_ext);
        element.calculate_deformation_gradient(points);
        
        // Get the deformation gradient and metric tensor
        const Eigen::Matrix2d& F = element.getDeformationGradient();
        Eigen::Matrix2d C = element.getMetricTensor();
        
        // Apply Lagrange reduction like in minimize_energy_with_triangles
        lagrange::Result result = lagrange::reduce(C);
        C = result.C_reduced;
        
        // Use reference area instead of current area
        double element_area = element.getReferenceArea();
        
        // Calculate element energy
        double element_energy = calculator.calculate_energy(C, energy_function, zero_energy)/normalisation;
        element_energy *= element_area;  // Use reference area
        recalculated_total_energy += element_energy;
        
        // Calculate first Piola-Kirchhoff stress tensor with Lagrange reduction
        Eigen::Matrix2d dE_dC = calculator.calculate_derivative(C, derivative_function)/normalisation;
        Eigen::Matrix2d P = 2.0 * F * result.m_matrix * dE_dC * result.m_matrix.transpose();
        
        // Calculate the stress = P:F_ext (double contraction of P and F_ext)
        double stress_value = 0.0;
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                stress_value += P(i, j) * F_ext(i, j);
            }
        }
        
        // Distribute to all nodes of this element
        for (int node_idx = 0; node_idx < 3; node_idx++) {
            int global_idx = element.getNodeIndex(node_idx);
            
            nodal_stress[global_idx] += stress_value;
            nodal_energy[global_idx] += element_energy / 3.0;  // Divide by 3 nodes for triangle
            node_count[global_idx]++;
        }
        
        processed_elements++;
    }
    
    // DEBUG: Print number of elements that were actually processed
    std::cout << "DEBUG: Processed " << processed_elements << " elements out of " 
              << active_elements.size() << " active elements." << std::endl;
    std::cout << "DEBUG: Recalculated total energy: " << recalculated_total_energy << std::endl;
    std::cout << "DEBUG: Passed total energy: " << total_energy << std::endl;
    
    // Count nodes that received data
    int nodes_with_data = 0;
    for (size_t i = 0; i < points.size(); i++) {
        if (node_count[i] > 0) {
            nodes_with_data++;
        }
    }
    std::cout << "DEBUG: " << nodes_with_data << " nodes out of " << points.size() 
              << " received stress/energy data." << std::endl;
    
    // Average stress values where nodes were counted multiple times
    for (size_t i = 0; i < points.size(); i++) {
        if (node_count[i] > 0) {
            nodal_stress[i] /= node_count[i];
            // Note: energy is already distributed per node, no need to average
        }
    }
    
    // Write XYZ format (number of atoms followed by comment line)
    file << points.size() << std::endl;
    file << "Iteration " << iteration << " Total_Energy " << recalculated_total_energy << std::endl;
    
    // Write atom positions with stress and energy as custom properties
    double value = 0;
    for (size_t i = 0; i < points.size(); i++) {
        const auto& point = points[i];
        const auto& pair = full_mapping[i];
        value = 0;
        if(pair.second == -1) 
            value = -1; 
        file << "A " << std::fixed << std::setprecision(8)
             << point.coord.x() << " "
             << point.coord.y() << " "
             << 0.0 << " "                // z-coordinate is zero for 2D
             << nodal_stress[i] << " "    // Add stress value as fourth column
             << nodal_energy[i] << " "
             << value << std::endl;       // Add energy value as fifth column
    }    
    file.close();
    
    // Find min and max values for reporting
    double min_stress = *std::min_element(nodal_stress.begin(), nodal_stress.end());
    double max_stress = *std::max_element(nodal_stress.begin(), nodal_stress.end());
    
    double min_energy = *std::min_element(nodal_energy.begin(), nodal_energy.end());
    double max_energy = *std::max_element(nodal_energy.begin(), nodal_energy.end());
    
    std::cout << "Saved 2D configuration with stress and energy data to " << filename.str() << std::endl;
    std::cout << "Stress range: [" << min_stress << ", " << max_stress << "]" << std::endl;
    std::cout << "Energy range: [" << min_energy << ", " << max_energy << "]" << std::endl;
}