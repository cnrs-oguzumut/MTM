#include "../include/output/configuration_saver.h"
#include "../include/reductions/LagrangeReduction.h"

void ConfigurationSaver::ensureDirectoryExists(const std::string& directory) {
    std::filesystem::create_directory(directory);
}

// Function to save configuration with nodal stress and energy values to XYZ file for OVITO (2D version)
void ConfigurationSaver::saveConfigurationWithStressAndEnergy2D(
    UserData* userData,
    int iteration, 
    double& total_energy,
    double& total_stress,bool reduction) {
        
    // Perform null check
    if (!userData) {
        std::cerr << "Error: userData is null in saveConfigurationWithStressAndEnergy2D" << std::endl;
        return;
    }
    
    // Extract needed values from userData with non-const references
    std::vector<Point2D>& points = userData->points;
    std::vector<ElementTriangle2D>& elements = userData->elements;
    //TriangularLatticeCalculator& calculator = userData->calculator;
    BaseLatticeCalculator& calculator = userData->calculator;
    std::function<double(double)>& energy_function = userData->energy_function;
    std::function<double(double)>& derivative_function = userData->derivative_function;
    double zero_energy = userData->zero_energy;
    double ideal_lattice_parameter = userData->ideal_lattice_parameter;
    const Eigen::Matrix2d& F_ext = userData->F_external;
    const std::vector<std::pair<int, int>>& interior_mapping = userData->interior_mapping;
    const std::vector<std::pair<int, int>>& full_mapping = userData->full_mapping;
    const std::vector<size_t>& active_elements = userData->active_elements;
    
    // For square lattice in 2D - match the normalization from minimize_energy_with_triangles
    //double normalisation = pow(ideal_lattice_parameter, 2.0); 
    double normalisation = calculator.getUnitCellArea();
   
    
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
    // Add vectors for Cauchy stress components
    std::vector<double> nodal_cauchy_xx(points.size(), 0.0);
    std::vector<double> nodal_cauchy_xy(points.size(), 0.0);
    std::vector<double> nodal_cauchy_yy(points.size(), 0.0);
    std::vector<int> node_count(points.size(), 0);
    
    // DEBUG: Print basic info
    // std::cout << "DEBUG: Number of points: " << points.size() << std::endl;
    // std::cout << "DEBUG: Number of elements: " << elements.size() << std::endl;
    // std::cout << "DEBUG: Number of active elements: " << active_elements.size() << std::endl;
    // std::cout << "DEBUG: Normalization value: " << normalisation << std::endl;
    // std::cout << "DEBUG: F_external matrix:\n" << F_ext << std::endl;
    
    int processed_elements = 0;
    double recalculated_total_energy = 0.0;
    double recalculated_total_stress = 0.0;
    
    // Loop through active elements to accumulate stress and energy at nodes
    for (size_t elem_idx : active_elements) {
        if (elem_idx >= elements.size()) {
            std::cerr << "Warning: Element index " << elem_idx << " out of range." << std::endl;
            continue;
        }
        
        auto& element = elements[elem_idx]; // Remove const to allow setting external deformation
        Eigen::Matrix2d P;
        lagrange::Result result;

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
        
        if(reduction==true) {
            // Apply Lagrange reduction
            result = lagrange::reduce(C);
            C = result.C_reduced;
        }
        

        //C = result.C_reduced;
        
        // Use reference area instead of current area
        double element_area = element.getReferenceArea();

        //std::cout<<elements[elem_idx].getArea()<<" "<<elements[elem_idx].getReferenceArea()<<std::endl;
        
        // Calculate element energy
        double element_energy = calculator.calculate_energy(C, energy_function, zero_energy)/normalisation;
        element_energy *= element_area;  // Use reference area
        recalculated_total_energy += element_energy;
        
        // Calculate first Piola-Kirchhoff stress tensor with Lagrange reduction
        Eigen::Matrix2d dE_dC = calculator.calculate_derivative(C, derivative_function)/normalisation;
        if(reduction==true) {
            P= 2.0 * F * result.m_matrix * dE_dC * result.m_matrix.transpose();
        }
        else {
            P = 2.0 * F * dE_dC;
        }

        // Calculate Cauchy stress tensor from PK1
        // σ = (1/det(F)) * P * F^T
        double detF = F.determinant();
        Eigen::Matrix2d cauchy_stress = (1.0 / detF) * P * F.transpose();
        
        // Calculate the stress = P:F_ext (double contraction of P and F_ext)
        Eigen::Matrix2d dF_d_alpha = Eigen::Matrix2d::Zero();
        dF_d_alpha(0, 1) = 1.0;
        double stress_value = P.cwiseProduct(dF_d_alpha).sum();
        
        // Accumulate total cauchy stress weighted by element area
        //recalculated_total_stress += stress_value * element.getArea();
        double current_area = element.calculateCurrentArea(points);
        recalculated_total_stress += cauchy_stress(0,1) * current_area;

        // Distribute to all nodes of this element
        for (int node_idx = 0; node_idx < 3; node_idx++) {
            int global_idx = element.getNodeIndex(node_idx);
            
            nodal_stress[global_idx] += stress_value;
            nodal_energy[global_idx] += element_energy / 3.0;  // Divide by 3 nodes for triangle
            // Add Cauchy stress components
            nodal_cauchy_xx[global_idx] += cauchy_stress(0, 0);
            nodal_cauchy_xy[global_idx] += cauchy_stress(0, 1);
            nodal_cauchy_yy[global_idx] += cauchy_stress(1, 1);
            node_count[global_idx]++;
        }
        
        processed_elements++;
    }
    
    // Set the output values
    total_energy = recalculated_total_energy;
    total_stress = recalculated_total_stress;
    
    // Normalize total stress by total area if needed
    double total_area = 0.0;
    for (size_t elem_idx : active_elements) {
        if (elem_idx < elements.size()) {
            total_area += elements[elem_idx].calculateCurrentArea(points);
        }
    }
    if (total_area > 0.0) {
        double avg_stress = recalculated_total_stress / total_area;
        //std::cout << "DEBUG: Average stress (total/area): " << avg_stress << std::endl;
    }
    
    // DEBUG: Print number of elements that were actually processed
    // std::cout << "DEBUG: Processed " << processed_elements << " elements out of " 
    //           << active_elements.size() << " active elements." << std::endl;
    // std::cout << "DEBUG: Recalculated total energy: " << recalculated_total_energy << std::endl;
    // std::cout << "DEBUG: Recalculated total stress: " << recalculated_total_stress << std::endl;
    
    // Count nodes that received data
    int nodes_with_data = 0;
    for (size_t i = 0; i < points.size(); i++) {
        if (node_count[i] > 0) {
            nodes_with_data++;
        }
    }
    // std::cout << "DEBUG: " << nodes_with_data << " nodes out of " << points.size() 
    //           << " received stress/energy data." << std::endl;
    
    // Average stress values where nodes were counted multiple times
    for (size_t i = 0; i < points.size(); i++) {
        if (node_count[i] > 0) {
            nodal_stress[i] /= node_count[i];
            // Average Cauchy stress components
            nodal_cauchy_xx[i] /= node_count[i];
            nodal_cauchy_xy[i] /= node_count[i];
            nodal_cauchy_yy[i] /= node_count[i];
            // Note: energy is already distributed per node, no need to average
        }
    }
    
    // Write XYZ format (number of atoms followed by comment line)
    file << points.size() << std::endl;
    file << "Iteration " << iteration << " Total_Energy " << recalculated_total_energy 
         << " Total_Stress " << recalculated_total_stress << std::endl;
    
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
             << nodal_stress[i] << " "    // Stress value (4th column)
             << nodal_energy[i] << " "    // Energy value (5th column)
             << nodal_cauchy_xx[i] << " " // Cauchy xx component (6th column)
             << nodal_cauchy_xy[i] << " " // Cauchy xy component (7th column)
             << nodal_cauchy_yy[i] << " " // Cauchy yy component (8th column)
             << value << std::endl;       // Boundary marker (9th column)
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

void ConfigurationSaver::calculateEnergyAndStress(UserData* userData, 
    double& total_energy, 
    Eigen::Matrix2d& total_stress, bool reduction) {
        
    // Check if userData is valid
    if (!userData) {
    std::cerr << "Error: userData is null in calculateEnergyAndStress" << std::endl;
    total_energy = 0.0;
    total_stress.setZero();
    return;
    }

    // Extract needed values from userData
    std::vector<Point2D>& points = userData->points;
    std::vector<ElementTriangle2D>& elements = userData->elements;
    //TriangularLatticeCalculator& calculator = userData->calculator;
    BaseLatticeCalculator& calculator = userData->calculator;
    std::function<double(double)>& energy_function = userData->energy_function;
    std::function<double(double)>& derivative_function = userData->derivative_function;
    double zero_energy = userData->zero_energy;
    double ideal_lattice_parameter = userData->ideal_lattice_parameter;
    const Eigen::Matrix2d& F_ext = userData->F_external;
    const std::vector<size_t>& active_elements = userData->active_elements;

    // For square lattice in 2D - match the normalization from minimize_energy_with_triangles
    double normalisation = calculator.getUnitCellArea();
    //normalisation = sqrt(3.)/2.*pow(ideal_lattice_parameter, 2.0);

    // Initialize totals
    total_energy = 0.0;
    total_stress.setZero();
    double total_area = 0.0;

    // Loop through active elements to calculate energy and stress
    for (size_t elem_idx : active_elements) {
        if (elem_idx >= elements.size()) {
            std::cerr << "Warning: Element index " << elem_idx << " out of range." << std::endl;
            continue;
        }
        auto& element = elements[elem_idx];
        Eigen::Matrix2d P;
        lagrange::Result result;
        

        // Skip if the element isn't initialized
        if (!element.isInitialized()) {
            continue;
        }

        // Set external deformation and recalculate deformation gradient
        element.setExternalDeformation(F_ext);
        element.calculate_deformation_gradient(points);

        // Get the deformation gradient and metric tensor
        const Eigen::Matrix2d& F = element.getDeformationGradient();
        Eigen::Matrix2d C = element.getMetricTensor();
        if(reduction==true) {
            // Apply Lagrange reduction
            result = lagrange::reduce(C);
            C = result.C_reduced;
        }


        // Use reference area instead of current area
        double element_area = element.getReferenceArea();

        // Calculate element energy
        double element_energy = calculator.calculate_energy(C, energy_function, zero_energy)/normalisation;
        element_energy *= element_area; // Use reference area
        total_energy += element_energy;

        // Calculate first Piola-Kirchhoff stress tensor
        Eigen::Matrix2d dE_dC = calculator.calculate_derivative(C, derivative_function)/normalisation;
        
        if(reduction==true) {
            P= 2.0 * F * result.m_matrix * dE_dC * result.m_matrix.transpose();
        }
        else {
            P = 2.0 * F * dE_dC;
        }


        // Calculate Cauchy stress tensor
        double detF = F.determinant();
        Eigen::Matrix2d cauchy_stress = (1.0 / detF) * F * P.transpose();

        // Accumulate total stress weighted by element area
        double current_area = element.calculateCurrentArea(points);
        total_stress += cauchy_stress*current_area;
        total_area += current_area;
        }

        // Normalize stress by total area
        if (total_area > 0.0) {
            total_stress /= total_area; // Normalize to get average stress
        }
    }

void ConfigurationSaver::logEnergyAndStress(
    int iteration, 
    double alpha, 
    double pre_energy, 
    double pre_stress,
    double post_energy, 
    double post_stress,
    int plasticity_flag)  // Added plasticity flag parameter
{
    static bool first_call = true;
    static std::ofstream log_file;
    
    if (first_call) {
        log_file.open("energy_stress_log.csv");
        log_file << "Iteration,Alpha,PreEnergy,PreStress,PostEnergy,PostStress,EnergyChange,StressChange,Plasticity\n";
        first_call = false;
    }
    
    log_file << iteration << "," << alpha << "," 
             << pre_energy << "," << pre_stress << "," 
             << post_energy << "," << post_stress << "," 
             << -(post_energy - pre_energy) << "," << -(post_stress - pre_stress) << ","
             << plasticity_flag << "\n";  // Added plasticity flag to CSV
    log_file.flush();
}



void ConfigurationSaver::writeToVTK(
    const std::vector<Point2D>& points,
    const std::vector<ElementTriangle2D>& elements,
    const UserData* userData,
    int iteration, bool reduction)
{
    // Create directory if it doesn't exist
    std::filesystem::create_directory("vtk_output");
    
    // Create filename with iteration number
    std::stringstream filename;
    filename << "vtk_output/configuration_" << std::setw(5) << std::setfill('0') << iteration << ".vtk";
    
    // Get original domaingit points count
    int original_points_count = points.size();
    
    // Define a custom key for map to handle Eigen::Vector2d
    struct PairKey {
        int index;
        double tx, ty; // Translation components
        
        PairKey(int idx, const Eigen::Vector2d& trans) : 
            index(idx), tx(trans.x()), ty(trans.y()) {}
            
        bool operator<(const PairKey& other) const {
            if (index != other.index) return index < other.index;
            if (tx != other.tx) return tx < other.tx;
            return ty < other.ty;
        }
    };
    
    // Build a map of unique points including periodically translated ones
    std::map<PairKey, int> node_map;
    int extended_point_count = 0;
    
    // First pass: identify all unique points (original + translated)
    for (const auto& element : elements) {
        if (!element.isInitialized()) continue;
        
        for (int i = 0; i < 3; i++) {
            int node_idx = element.getNodeIndex(i);
            Eigen::Vector2d translation = element.getTranslation(i);
            
            PairKey node_key(node_idx, translation);
            if (node_map.find(node_key) == node_map.end()) {
                node_map[node_key] = extended_point_count++;
            }
        }
    }
    
    std::cout << "Original points: " << original_points_count 
              << ", Extended points: " << extended_point_count << std::endl;
    
    // Prepare arrays for extended point coordinates
    std::vector<Eigen::Vector2d> extended_points(extended_point_count);
    
// Inside the writeToVTK function, when filling extended point coordinates:

// Fill extended point coordinates
for (const auto& entry : node_map) {
    int original_idx = entry.first.index;
    Eigen::Vector2d translation(entry.first.tx, entry.first.ty);
    int extended_idx = entry.second;
    
    // Apply translation to the original point
    if (original_idx < original_points_count) {
        // Apply external deformation to the translation vector
        Eigen::Vector2d deformed_translation = userData->F_external * translation;
        
        // Add the deformed translation to the original point
        extended_points[extended_idx] = points[original_idx].coord + deformed_translation;
    } else {
        std::cerr << "Warning: Invalid point index " << original_idx << std::endl;
        extended_points[extended_idx] = Eigen::Vector2d::Zero();
    }
}    
    // Prepare nodal data vectors
    std::vector<double> nodal_energy(extended_point_count, 0.0);
    std::vector<double> nodal_cauchy_xy(extended_point_count, 0.0);
    std::vector<double> nodal_cauchy_xx(extended_point_count, 0.0);
    std::vector<double> nodal_cauchy_yy(extended_point_count, 0.0);
    std::vector<double> nodal_stress(extended_point_count, 0.0);
    std::vector<int> triangle_count(extended_point_count, 0);
    std::vector<int> node_count(extended_point_count, 0);
    
    // Count triangles per node and calculate nodal values
    for (size_t elem_idx : userData->active_elements) {
        if (elem_idx >= elements.size()) continue;
        
        const auto& element = elements[elem_idx];
        if (!element.isInitialized()) continue;
        
        // Calculate element quantities
        const Eigen::Matrix2d& F = element.getDeformationGradient();
        Eigen::Matrix2d C = element.getMetricTensor();
        lagrange::Result result;
        Eigen::Matrix2d P;
    
        
        // Apply Lagrange reduction 
        if(reduction==true) {
            result = lagrange::reduce(C);
            C = result.C_reduced;
        }
        C = result.C_reduced;
        
        
        // Calculate element energy
        double element_area = element.getReferenceArea();
        double element_energy = userData->calculator.calculate_energy(
            C, userData->energy_function, userData->zero_energy
        ) / pow(userData->ideal_lattice_parameter, 2.0);
        element_energy *= element_area;
        
        // Calculate Cauchy stress
        Eigen::Matrix2d dE_dC = userData->calculator.calculate_derivative(
            C, userData->derivative_function
        ) / pow(userData->ideal_lattice_parameter, 2.0);
 //       Eigen::Matrix2d P = 2.0 * F * result.m_matrix * dE_dC * result.m_matrix.transpose();
        if(reduction==true) {
            P= 2.0 * F * result.m_matrix * dE_dC * result.m_matrix.transpose();
        }
        else {
            P = 2.0 * F * dE_dC;
        }

        double detF = F.determinant();
        Eigen::Matrix2d cauchy_stress = (1.0 / detF) * P * F.transpose();
        
        // Calculate projected stress
        Eigen::Matrix2d dF_d_alpha = Eigen::Matrix2d::Zero();
        dF_d_alpha(0, 1) = 1.0;
        double stress_value = P.cwiseProduct(dF_d_alpha).sum();
        
        // Distribute to nodes (using the extended node indices)
        for (int node_idx = 0; node_idx < 3; node_idx++) {
            int original_idx = element.getNodeIndex(node_idx);
            Eigen::Vector2d translation = element.getTranslation(node_idx);
            
            PairKey node_key(original_idx, translation);
            int extended_idx = node_map[node_key];
            
            triangle_count[extended_idx]++;
            nodal_energy[extended_idx] += element_energy / 3.0;
            nodal_cauchy_xy[extended_idx] += cauchy_stress(0, 1);
            nodal_cauchy_xx[extended_idx] += cauchy_stress(0, 0);
            nodal_cauchy_yy[extended_idx] += cauchy_stress(1, 1);

            nodal_stress[extended_idx] += stress_value;
            node_count[extended_idx]++;
        }
    }
    
    // Average nodal values
    for (int i = 0; i < extended_point_count; i++) {
        if (node_count[i] > 0) {
            nodal_energy[i] /= node_count[i];
            nodal_cauchy_xy[i] /= node_count[i];
            nodal_cauchy_xx[i] /= node_count[i];
            nodal_cauchy_yy[i] /= node_count[i];

            nodal_stress[i] /= node_count[i];
        }
    }
    
    // Open file for writing
    std::ofstream file(filename.str());
    if (!file) {
        std::cerr << "Error: Could not open file " << filename.str() << " for writing." << std::endl;
        return;
    }
    
    // Write VTK header
    file << "# vtk DataFile Version 1.0\n";
    file << "2D Unstructured Grid of Linear Triangles\n";
    file << "ASCII\n\n";
    file << "DATASET UNSTRUCTURED_GRID\n";
    
    // Write points (extended set)
    file << "POINTS " << extended_point_count << " float\n";
    for (const auto& point : extended_points) {
        file << point.x() << " " << point.y() << " " << 0.0 << "\n";
    }
    
    // Count valid elements
    int valid_elements = 0;
    for (size_t elem_idx : userData->active_elements) {
        if (elem_idx < elements.size() && elements[elem_idx].isInitialized()) {
            valid_elements++;
        }
    }
    
    // Write cells
    file << "CELLS " << valid_elements << " " << (4 * valid_elements) << "\n";
    for (size_t elem_idx : userData->active_elements) {
        if (elem_idx >= elements.size() || !elements[elem_idx].isInitialized()) continue;
        
        const auto& element = elements[elem_idx];
        file << "3";
        
        // Get the extended indices for each node of this element
        for (int i = 0; i < 3; i++) {
            int original_idx = element.getNodeIndex(i);
            Eigen::Vector2d translation = element.getTranslation(i);
            PairKey node_key(original_idx, translation);
            int extended_idx = node_map[node_key];
            
            file << " " << extended_idx;
        }
        file << "\n";
    }
    
    // Write cell types
    file << "CELL_TYPES " << valid_elements << "\n";
    for (int i = 0; i < valid_elements; i++) {
        file << "5\n"; // Triangle type
    }
    
    // Write point data
    file << "POINT_DATA " << extended_point_count << "\n";
    
    // Energy data
    file << "SCALARS Energy float\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < extended_point_count; i++) {
        file << nodal_energy[i] << "\n";
    }
    
    // Cauchy stress xy component
    file << "SCALARS cauchy_xy float\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < extended_point_count; i++) {
        file << nodal_cauchy_xy[i] << "\n";
    }
    
    // Cauchy stress xx component
    file << "SCALARS cauchy_xx float\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < extended_point_count; i++) {
        file << nodal_cauchy_xx[i] << "\n";
    }

    // Cauchy stress xx component
    file << "SCALARS cauchy_yy float\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < extended_point_count; i++) {
        file << nodal_cauchy_yy[i] << "\n";
    }

    // Projected stress
    file << "SCALARS projected_stress float\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < extended_point_count; i++) {
        file << nodal_stress[i] << "\n";
    }
    
    // // Triangle count for topological analysis
    // file << "SCALARS coordination int\n";
    // file << "LOOKUP_TABLE default\n";
    // for (int i = 0; i < extended_point_count; i++) {
    //     // Adjust this condition for your specific coordination criteria
    //     int coordination = triangle_count[i];
    //     if (coordination != 5 && coordination != 7) {
    //         coordination = 0;
    //     }
    //     file << coordination << "\n";
    // }
    
    // // Write cell data
    // file << "CELL_DATA " << valid_elements << "\n";
    // file << "TENSORS Stress float\n";
    
    // for (size_t elem_idx : userData->active_elements) {
    //     if (elem_idx >= elements.size() || !elements[elem_idx].isInitialized()) continue;
        
    //     const auto& element = elements[elem_idx];
        
    //     // Get stress tensors and other element data
    //     Eigen::Matrix2d C = element.getMetricTensor();
    //     lagrange::Result result = lagrange::reduce(C);
        
    //     // Write tensor data in VTK format (3x3 even for 2D)
    //     file << C(0, 1) << " " << C(0, 0) << " " << C(1, 1) << "\n";
    //     file << "0 " << result.third_condition_satisfied << " " << "0" << "\n";
    //     file << result.C_reduced(0, 1) << " " << C(0, 1) << " " << userData->F_external(0, 1) << "\n\n";
    // }
    
    file.close();
    std::cout << "Saved VTK file: " << filename.str() << std::endl;
}