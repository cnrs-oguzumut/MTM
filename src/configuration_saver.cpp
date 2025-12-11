#include "../include/output/configuration_saver.h"
#include "../include/reductions/LagrangeReduction.h"
#include "../include/geometry/DomainDimensions.h"


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
        total_stress = avg_stress;
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
file << "A " << std::fixed << std::setprecision(16)
     << point.coord.x() << " "
     << point.coord.y() << " "
     << 0.0 << " "                // z-coordinate is zero for 2D
     << nodal_stress[i] << " "    // Stress value (4th column)
     << nodal_energy[i] << " "    // Energy value (5th column)
     << nodal_cauchy_xx[i] << " " // Cauchy xx component (6th column)
     << nodal_cauchy_xy[i] << " " // Cauchy xy component (7th column)
     << nodal_cauchy_yy[i] << " " // Cauchy yy component (8th column)
     << value << " "              // Boundary marker (9th column)
     << F_ext(0, 0) << " "        // F_ext_xx component (10th column)
     << F_ext(0, 1) << " "        // F_ext_xy component (11th column)
     << F_ext(1, 0) << " "        // F_ext_yx component (12th column)
     << F_ext(1, 1) << std::endl; // F_ext_yy component (13th column)
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

// After saving the XYZ file, add this code to print row profiles

// Create directory for profiles if it doesn't exist
//std::filesystem::create_directory("profiles");
createFolder("./profiles");

// Group points into rows based on y-coordinates (within a tolerance)
const double row_tolerance = ideal_lattice_parameter * 0.6; // Adjust tolerance as needed
std::cout<<"Row tolerance: " << row_tolerance << std::endl;
std::map<double, std::vector<size_t>> rows;

// Group points by their y-coordinate
for (size_t i = 0; i < points.size(); i++) {
    bool added = false;
    const auto& point = points[i];
    
    // Skip boundary points if desired
    if (full_mapping[i].second == -1) continue;
    
    // Try to add to existing row
    for (auto& row : rows) {
        if (std::abs(point.coord.y() - row.first) < row_tolerance) {
            row.second.push_back(i);
            added = true;
            break;
        }
    }
    
    // If not added to existing row, create new row
    if (!added) {
        rows[point.coord.y()] = {i};
    }
}

// Convert map to vector for easier sorting
std::vector<std::pair<double, std::vector<size_t>>> sorted_rows;
for (const auto& row : rows) {
    sorted_rows.push_back(row);
}

// Sort rows by y-coordinate
std::sort(sorted_rows.begin(), sorted_rows.end(), 
    [](const auto& a, const auto& b) { return a.first < b.first; });

// Find middle row
size_t middle_row_idx = sorted_rows.size() / 2;
std::cout<<"middle_row_idx: " << middle_row_idx << std::endl;

int rows_to_show =3; // Number of rows to show above and below middle

// Output each row to a separate file
for (int r = -rows_to_show; r <= rows_to_show; r++) {
    int row_idx = middle_row_idx + r;
    if (row_idx >= 0 && row_idx < sorted_rows.size()) {
        const auto& row = sorted_rows[row_idx];
        
        // Create filename for this row
        std::stringstream row_filename;
        row_filename << "profiles/row_" 
                     << std::setw(2) << std::setfill('0') << (r + rows_to_show) 
                     << "_iter_" << std::setw(5) << std::setfill('0') << iteration 
                     << ".dat";
        
        std::ofstream row_file(row_filename.str());
        if (!row_file) {
            std::cerr << "Error: Could not open file " << row_filename.str() << " for writing." << std::endl;
            continue;
        }
        
        // Write header
        row_file << "# Row profile for y = " << row.first 
                << " (row " << (r + rows_to_show) << ", " 
                << (r < 0 ? std::abs(r) : 0) << " rows below middle, "
                << (r > 0 ? r : 0) << " rows above middle)\n";
        row_file << "# point_idx x_coord stress energy cauchy_xx cauchy_xy cauchy_yy\n";
        
        // Sort points in this row by x-coordinate
        std::vector<size_t> sorted_points = row.second;
        std::sort(sorted_points.begin(), sorted_points.end(), 
            [&points](size_t a, size_t b) { return points[a].coord.x() < points[b].coord.x(); });
        
        // Output each point in this row
        for (size_t point_idx : sorted_points) {
            row_file << std::fixed << std::setprecision(8)
                << point_idx << " "               // Point index
                << points[point_idx].coord.x() << " " // X-coordinate
                << nodal_stress[point_idx] << " "     // Stress
                << nodal_energy[point_idx] << " "     // Energy
                << nodal_cauchy_xx[point_idx] << " "  // Cauchy xx
                << nodal_cauchy_xy[point_idx] << " "  // Cauchy xy
                << nodal_cauchy_yy[point_idx] << " "  // Cauchy yy
                << "\n";
        }
        
        row_file.close();
        //std::cout << "Row " << (r + rows_to_show) << " profile saved to " << row_filename.str() << std::endl;
    }
}


}

void ConfigurationSaver::saveTriangleData(
    UserData* userData,
    int iteration,
    const DomainDimensions& domain_dims,
    const std::array<double, 2>& offsets,
    const std::vector<std::pair<int, int> >& full_mapping) {
        
    // Perform null check
    if (!userData) {
        std::cerr << "Error: userData is null in saveTriangleData" << std::endl;
        return;
    }
    
    // Extract needed values from userData
    std::vector<Point2D>& points = userData->points;
    std::vector<ElementTriangle2D>& elements = userData->elements;
    const Eigen::Matrix2d& F_ext = userData->F_external;
    const std::vector<size_t>& active_elements = userData->active_elements;
    
    // Create directory if it doesn't exist
    std::filesystem::create_directory("triangle_data");
    
    // Create filename with iteration number
    std::stringstream filename;
    filename << "triangle_data/triangles_" << std::setw(5) << std::setfill('0') << 10000+iteration << ".dat";
    
    // Open file for writing
    std::ofstream file(filename.str());
    if (!file) {
        std::cerr << "Error: Could not open file " << filename.str() << " for writing." << std::endl;
        return;
    }
    
    // Write header
    file << "# Triangle data for iteration " << iteration << std::endl;
    file << "# element_number F11 F12 F21 F22 node1_idx node2_idx node3_idx" << std::endl;
    
    // Set maximum precision
    file << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10);
    
    // Loop through all active elements
    for (size_t elem_idx : active_elements) {
        if (elem_idx >= elements.size()) {
            std::cerr << "Warning: Element index " << elem_idx << " out of range." << std::endl;
            continue;
        }
        
        auto& element = elements[elem_idx];
        
        // Skip if the element isn't initialized
        if (!element.isInitialized()) {
            continue;
        }
        
        // Set external deformation and recalculate deformation gradient
        element.setExternalDeformation(F_ext);
        element.calculate_deformation_gradient(points);

        // Get the deformation gradient
        const Eigen::Matrix2d& F = element.getDeformationGradient();
        
        // Write element data: element_number, F components, node indices
        Eigen::Vector2d translation0 = element.getTranslation(0);
        Eigen::Vector2d translation1 = element.getTranslation(1);
        Eigen::Vector2d translation2 = element.getTranslation(2);


        const Eigen::Matrix<double, 3, 2>& dndx = element.getDNdX();
        double detF = element.getDeformationGradient().determinant();
        if(detF <0 )  {
            std::cerr << "Warning: Negative determinant detected for element " << elem_idx << ", reflecting the triangle." << std::endl;
            std::cerr << "Warning: Negative determinant detected for element " << element.getDeformationGradient().determinant()  << ", reflecting the triangle." << std::endl;
            
        }

        if(F(0,0)*F(1,1) - F(0,1)*F(1,0) <0 )  {
            std::cerr << "#Warning: Negative determinant detected for element " << elem_idx << ", reflecting the triangle." << std::endl;
            std::cerr << "#Warning: Negative determinant detected for element " << element.getDeformationGradient().determinant()  << ", reflecting the triangle." << std::endl;
            
        }
    
        int i1 = element.getNodeIndex(0);
        int i2 = element.getNodeIndex(1);
        int i3 = element.getNodeIndex(2);
        if(full_mapping[i1].second == -1 || full_mapping[i2].second == -1|| full_mapping[i3].second  == -1) {
            // std::cerr << "Warning: boundary node index detected for element " << elem_idx << ", skipping." << std::endl;
            continue;
        }

        file << elem_idx << " "                                    // Element number
             << F(0,0) << " "                                      // F11
             << F(0,1) << " "                                      // F12
             << F(1,0) << " "                                      // F21
             << F(1,1) << " "                                      // F22
             << element.getNodeIndex(0) << " "                     // Global index of node 1
             << element.getNodeIndex(1) << " "                     // Global index of node 2
             << element.getNodeIndex(2) << " "                     // Global index of node 3
             << points[element.getNodeIndex(0)].coord.x() << " "   // x coordinate of node 1
             << points[element.getNodeIndex(0)].coord.y() << " "   // y coordinate of node 1
             << points[element.getNodeIndex(1)].coord.x() << " "   // x coordinate of node 2
             << points[element.getNodeIndex(1)].coord.y() << " "   // y coordinate of node 2
             << points[element.getNodeIndex(2)].coord.x() << " "   // x coordinate of node 3
             << points[element.getNodeIndex(2)].coord.y() << " "   // y coordinate of node 3
             << translation0.x() << " "                            // Effective translation x
             << translation0.y() << " "
             << translation1.x() << " "                            // Effective translation y
             << translation1.y() << " "                            // Effective translation z
             << translation2.x() << " "                            // Effective translation x
             << translation2.y() << " "                            // Effective translation y
             << domain_dims.size_x << " "                          // Domain size x
             << domain_dims.size_y << " "                          // Domain size y
             << offsets[0] << " "                                  // Offset x
             << offsets[1] << " "                                  // Offset y
             << F_ext(0,0) << " "                                  // F_ext(0,0)
             << F_ext(0,1) << " "                                  // F_ext(0,1)
             << F_ext(1,0) << " "                                  // F_ext(1,0)
             << F_ext(1,1) << " "                                  // F_ext(1,1)
             << dndx(0,0) << " "                                   // dN1/dx (shape function derivative)
             << dndx(0,1) << " "                                   // dN1/dy
             << dndx(1,0) << " "                                   // dN2/dx
             << dndx(1,1) << " "                                   // dN2/dy
             << dndx(2,0) << " "                                   // dN3/dx
             << dndx(2,1)                                          // dN3/dy
             << std::endl;                                         // End of line for this triangle
    }
   
    file.close();
    
    std::cout << "Saved triangle data to " << filename.str() << std::endl;
    std::cout << "Processed " << active_elements.size() << " triangles" << std::endl;
    
    // NOW SAVE POINT COORDINATES
    std::stringstream points_filename;
    points_filename << "triangle_data/points_" << std::setw(5) << std::setfill('0') << 10000+iteration << ".dat";
    
    std::ofstream points_file(points_filename.str());
    if (!points_file) {
        std::cerr << "Error: Could not open file " << points_filename.str() << " for writing." << std::endl;
        return;
    }
    
    // Write header with number of points
    points_file << "# Point coordinates for iteration " << iteration << std::endl;
    points_file << "# Number of points: " << points.size() << std::endl;
    points_file << "# Format: point_index x y" << std::endl;
    
    // Set maximum precision for points file
    points_file << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10);
    
    // Write each point
    for (size_t i = 0; i < points.size(); ++i) {
        points_file << i << " "
                    << points[i].coord.x() << " "
                    << points[i].coord.y() << " "
                    << full_mapping[i].second << " " // 1 = fixed, 0 = free
                    << F_ext(0,0) << " "                                  // F_ext(0,0)
                    << F_ext(0,1) << " "                                  // F_ext(0,1)
                    << F_ext(1,0) << " "                                  // F_ext(1,0)
                    << F_ext(1,1) << " "                                  // F_ext(1,1)
                    << std::endl;
    }
    points_file.close();
    
    std::cout << "Saved point coordinates to " << points_filename.str() << std::endl;
    std::cout << "Processed " << active_elements.size() << " triangles and " << points.size() << " points" << std::endl;
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
    //std::cout << "DEBUG: Normalization value: " << normalisation << std::endl;
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


void ConfigurationSaver::logEnergyAndStress_v2(
    int iteration, 
    double alpha, 
    double pre_energy, 
    double pre_stress,
    double post_energy, 
    double post_stress,
    double pre_area,
    double post_area,  
    bool shouldRemesh)  // Added plasticity flag parameter
{
    static bool first_call = true;
    static std::ofstream log_file;
    
    if (first_call) {
        log_file.open("energy_stress_log.csv");
        
        // Write header
        log_file << "Iteration,Alpha,PreEnergy,PreStress,PostEnergy,PostStress,"
                 << "EnergyChange,StressChange,PreArea,PostArea,shouldRemesh\n";
        
        // Set maximum precision for all subsequent writes
        log_file << std::scientific 
                 << std::setprecision(std::numeric_limits<double>::max_digits10);
        
        first_call = false;
    }
    
    // Write data (precision is already set)
    log_file << iteration << "," 
             << alpha << "," 
             << pre_energy << "," 
             << pre_stress << "," 
             << post_energy << "," 
             << post_stress << "," 
             << -(post_energy - pre_energy) << "," 
             << -(post_stress * post_area - pre_stress * pre_area) << ","
             << pre_area << ","
             << post_area << ","
             << shouldRemesh
             << "\n";
    
    log_file.flush();
}

// Start of the modified function body
void ConfigurationSaver::writeToVTK(
    const std::vector<Point2D>& points,
    const std::vector<ElementTriangle2D>& elements,
    const UserData* userData,
    int iteration, 
    bool reduction,
    const std::vector<int>& coordination,// NEW: Optional coordination vector
    double load_strength)  
{
    // Create directory if it doesn't exist
    std::filesystem::create_directory("vtk_output");
    
    // Create filename with iteration number
    std::stringstream filename;
    filename << "vtk_output/configuration_" << std::setw(5) << std::setfill('0') << iteration << ".vtk";
    
    // Get original domain points count
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
    
    // Fill extended point coordinates
    for (const auto& entry : node_map) {
        int original_idx = entry.first.index;
        Eigen::Vector2d translation(entry.first.tx, entry.first.ty);
        int extended_idx = entry.second;
        
        // Apply translation to the original point
        if (original_idx < original_points_count) {
            // Apply external deformation to the translation vector
            Eigen::Vector2d deformed_translation = userData->F_external * translation;
            extended_points[extended_idx] = points[original_idx].coord + deformed_translation;
        } else {
            std::cerr << "Warning: Invalid point index " << original_idx << std::endl;
            extended_points[extended_idx] = Eigen::Vector2d::Zero();
        }
    }    
    
    // NEW: Map coordination from original points to extended points (if provided)
    std::vector<int> extended_coordination;
    bool has_coordination = !coordination.empty() && 
                           coordination.size() == static_cast<size_t>(original_points_count);
    
    if (has_coordination) {
        extended_coordination.resize(extended_point_count, 6);  // Default to 6
        for (const auto& entry : node_map) {
            int original_idx = entry.first.index;
            int extended_idx = entry.second;
            if (original_idx < static_cast<int>(coordination.size())) {
                extended_coordination[extended_idx] = coordination[original_idx];
            }
        }
    }
    
    // Prepare vectors to store element-wise data
    std::vector<double> element_energy_values;
    std::vector<Eigen::Matrix2d> element_cauchy_stress_tensors;
    std::vector<double> element_projected_stress_values;
    
    // Initialize nodal value arrays
    std::vector<double> nodal_energy(extended_point_count, 0.0);
    std::vector<double> nodal_projected_stress(extended_point_count, 0.0);
    std::vector<double> nodal_cauchy_xx(extended_point_count, 0.0);
    std::vector<double> nodal_cauchy_xy(extended_point_count, 0.0);
    std::vector<double> nodal_cauchy_yy(extended_point_count, 0.0);
    std::vector<int> node_count(extended_point_count, 0);

    // Loop over active elements to calculate and store values
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
        
        // Calculate element energy
        double element_area = element.getReferenceArea();
        double element_energy = userData->calculator.calculate_energy(
            C, userData->energy_function, userData->zero_energy
        ) / pow(userData->ideal_lattice_parameter, 2.0);
        element_energy *= element_area;
        element_energy_values.push_back(element_energy);
        
        // Calculate Cauchy stress
        Eigen::Matrix2d dE_dC = userData->calculator.calculate_derivative(
            C, userData->derivative_function
        ) / pow(userData->ideal_lattice_parameter, 2.0);
        
        if(reduction==true) {
            P = 2.0 * F * result.m_matrix * dE_dC * result.m_matrix.transpose();
        } else {
            P = 2.0 * F * dE_dC;
        }

        double detF = F.determinant();
        Eigen::Matrix2d cauchy_stress = (1.0 / detF) * P * F.transpose();
        element_cauchy_stress_tensors.push_back(cauchy_stress);
        
        // Calculate projected stress
        Eigen::Matrix2d dF_d_alpha = Eigen::Matrix2d::Zero();
        dF_d_alpha(0, 1) = 1.0;
        double stress_value = P.cwiseProduct(dF_d_alpha).sum();
        element_projected_stress_values.push_back(stress_value);
        
        // Distribute to all 3 nodes of this triangle
        for (int i = 0; i < 3; i++) {
            int original_idx = element.getNodeIndex(i);
            Eigen::Vector2d translation = element.getTranslation(i);
            PairKey node_key(original_idx, translation);
            int extended_idx = node_map[node_key];
            
            nodal_energy[extended_idx] += element_energy / 3.0;  // Divide by 3 nodes
            nodal_projected_stress[extended_idx] += stress_value;
            nodal_cauchy_xx[extended_idx] += cauchy_stress(0, 0);
            nodal_cauchy_xy[extended_idx] += cauchy_stress(0, 1);
            nodal_cauchy_yy[extended_idx] += cauchy_stress(1, 1);
            node_count[extended_idx]++;
        }
    }
    
    // Average the stress values where nodes were counted multiple times
    for (int i = 0; i < extended_point_count; i++) {
        if (node_count[i] > 0) {
            nodal_projected_stress[i] /= node_count[i];
            nodal_cauchy_xx[i] /= node_count[i];
            nodal_cauchy_xy[i] /= node_count[i];
            nodal_cauchy_yy[i] /= node_count[i];
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
    

    // --- NEW: Write Field Data (Global parameters) ---
    file << "\nFIELD FieldData 1\n";
    file << "LoadParameter 1 1 float\n";
    file << load_strength << "\n";  // REPLACE with your actual load parameter

    // --- Write Cell Data (Per-Element Values) ---
    file << "\nCELL_DATA " << valid_elements << "\n";
    
    // Element Energy
    file << "SCALARS ElementEnergy float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (double val : element_energy_values) {
        file << val << "\n";
    }
    
    // Element Projected Stress
    file << "SCALARS ElementProjectedStress float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (double val : element_projected_stress_values) {
        file << val << "\n";
    }
    
    // Element Cauchy Stress Tensor
    file << "TENSORS CauchyStress float\n";
    for (const auto& stress : element_cauchy_stress_tensors) {
        file << stress(0, 0) << " " << stress(0, 1) << " " << 0.0 << "\n";
        file << stress(1, 0) << " " << stress(1, 1) << " " << 0.0 << "\n";
        file << 0.0 << " " << 0.0 << " " << 0.0 << "\n\n";
    }
    
    // --- Write Point Data (Per-Node Values) ---
    file << "\nPOINT_DATA " << extended_point_count << "\n";
    
    // Nodal Energy
    file << "SCALARS NodalEnergy float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (double val : nodal_energy) {
        file << val << "\n";
    }

    // Nodal Projected Stress
    file << "SCALARS NodalProjectedStress float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (double val : nodal_projected_stress) {
        file << val << "\n";
    }

    // Number of Elements per Node (Valence/Coordination) - from mesh
    file << "SCALARS NodeValence float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int val : node_count) {
        file << (val < 5 ? 6 : val) << "\n";
    }
    
    // NEW: Reference Coordination - from defect analysis (if provided)
    if (has_coordination) {
        file << "SCALARS ReferenceCoordination float 1\n";
        file << "LOOKUP_TABLE default\n";
        for (int val : extended_coordination) {
            file << (val < 5 ? 6 : val) << "\n";
        }
        std::cout << "  Included ReferenceCoordination field from defect analysis" << std::endl;
    }

    // Nodal Cauchy Stress Tensor
    file << "TENSORS NodalCauchyStress float\n";
    for (int i = 0; i < extended_point_count; i++) {
        file << nodal_cauchy_xx[i] << " " << nodal_cauchy_xy[i] << " " << 0.0 << "\n";
        file << nodal_cauchy_xy[i] << " " << nodal_cauchy_yy[i] << " " << 0.0 << "\n";
        file << 0.0 << " " << 0.0 << " " << 0.0 << "\n\n";
    }

    file.close();
    std::cout << "Saved VTK file: " << filename.str() << std::endl;
}



// Function to save configuration with nodal stress and energy values to XYZ file for OVITO (2D version)
// void ConfigurationSaver::saveConfigurationWithNeigborinformation(
//     UserData* userData,
//     int iteration, 
//     double& total_energy,
//     double& total_stress,bool reduction) {
        
//     // Perform null check
//     if (!userData) {
//         std::cerr << "Error: userData is null in saveConfigurationWithStressAndEnergy2D" << std::endl;
//         return;
//     }
    
//     // Extract needed values from userData with non-const references
//     std::vector<Point2D>& points = userData->points;
//     std::vector<ElementTriangle2D>& elements = userData->elements;
//     //TriangularLatticeCalculator& calculator = userData->calculator;
//     BaseLatticeCalculator& calculator = userData->calculator;
//     std::function<double(double)>& energy_function = userData->energy_function;
//     std::function<double(double)>& derivative_function = userData->derivative_function;
//     double zero_energy = userData->zero_energy;
//     double ideal_lattice_parameter = userData->ideal_lattice_parameter;
//     const Eigen::Matrix2d& F_ext = userData->F_external;
//     const std::vector<std::pair<int, int>>& interior_mapping = userData->interior_mapping;
//     const std::vector<std::pair<int, int>>& full_mapping = userData->full_mapping;
//     const std::vector<size_t>& active_elements = userData->active_elements;
    
//     // Transform nodal positions by inverse of F_ext to remove external load state
//     Eigen::Matrix2d F_ext_inv = F_ext.inverse();
//     for (auto& point : points) {
//         Eigen::Vector2d pos(point.x, point.y);
//         Eigen::Vector2d transformed_pos = F_ext_inv * pos;
//         point.x = transformed_pos(0);
//         point.y = transformed_pos(1);
//     }
    
//     // For square lattice in 2D - match the normalization from minimize_energy_with_triangles
//     //double normalisation = pow(ideal_lattice_parameter, 2.0); 
//     double normalisation = calculator.getUnitCellArea();
   
    
//     // Create directory if it doesn't exist
//     std::filesystem::create_directory("Inverse_F");
    
    
    
//     int processed_elements = 0;
//     double recalculated_total_energy = 0.0;
//     double recalculated_total_stress = 0.0;
    
//     // Loop through active elements to accumulate stress and energy at nodes
//     for (size_t elem_idx : active_elements) {
//         if (elem_idx >= elements.size()) {
//             std::cerr << "Warning: Element index " << elem_idx << " out of range." << std::endl;
//             continue;
//         }
//     }
// }



void ConfigurationSaver::writeToVTK_DefectAnalysis(
    const std::vector<Point2D>& points,
    const std::vector<ElementTriangle2D>& elements,
    const UserData* userData,
    int iteration)
{
    // Create directory if it doesn't exist
    std::filesystem::create_directory("vtk_defects");
    
    // Create filename with iteration number
    std::stringstream filename;
    filename << "vtk_defects/defects_" << std::setw(5) << std::setfill('0') << iteration << ".vtk";
    
    // Get original domain points count
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
    
    std::cout << "Defect analysis - Original points: " << original_points_count 
              << ", Extended points: " << extended_point_count << std::endl;
    
    // Prepare arrays for extended point coordinates
    std::vector<Eigen::Vector2d> extended_points(extended_point_count);
    
    // Fill extended point coordinates
    for (const auto& entry : node_map) {
        int original_idx = entry.first.index;
        Eigen::Vector2d translation(entry.first.tx, entry.first.ty);
        int extended_idx = entry.second;
        
        // Apply translation to the original point
        if (original_idx < original_points_count) {
            // Apply external deformation to the translation vector
            Eigen::Vector2d deformed_translation = userData->F_external * translation;
            extended_points[extended_idx] = points[original_idx].coord + deformed_translation;
        } else {
            std::cerr << "Warning: Invalid point index " << original_idx << std::endl;
            extended_points[extended_idx] = Eigen::Vector2d::Zero();
        }
    }    
    
    // Initialize nodal coordination number array
    std::vector<int> node_count(extended_point_count, 0);

    // Loop over active elements to count coordination numbers
    for (size_t elem_idx : userData->active_elements) {
        if (elem_idx >= elements.size()) continue;
        
        const auto& element = elements[elem_idx];
        if (!element.isInitialized()) continue;
        
        // Count coordination for each node of this triangle
        for (int i = 0; i < 3; i++) {
            int original_idx = element.getNodeIndex(i);
            Eigen::Vector2d translation = element.getTranslation(i);
            PairKey node_key(original_idx, translation);
            int extended_idx = node_map[node_key];
            
            node_count[extended_idx]++;
        }
    }
    
    // Fix valence for periodic copies - copy from original points
    for (const auto& entry : node_map) {
        int original_idx = entry.first.index;
        Eigen::Vector2d translation(entry.first.tx, entry.first.ty);
        int extended_idx = entry.second;
        
        // If this is a periodic copy (has non-zero translation)
        if (translation.norm() > 1e-10) {
            // Find the corresponding original point (same index, zero translation)
            PairKey original_key(original_idx, Eigen::Vector2d::Zero());
            auto it = node_map.find(original_key);
            if (it != node_map.end()) {
                int original_extended_idx = it->second;
                // Copy the valence from the original to the periodic copy
                node_count[extended_idx] = node_count[original_extended_idx];
            }
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
    file << "2D Defect Analysis - Coordination Numbers Only\n";
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
    
    // --- Write Point Data (Per-Node Values) - ONLY NodeValence ---
    file << "\nPOINT_DATA " << extended_point_count << "\n";
    
    // Number of Elements per Node (Valence/Coordination)
    file << "SCALARS NodeValence float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int val : node_count) {
        file << (val < 5 ? 6 : val) << "\n";  // Set < 5 to 6
    }
    file.close();
    std::cout << "Saved defect analysis VTK file: " << filename.str() << std::endl;
    
    // Print statistics
    std::map<int, int> valence_histogram;
    for (int val : node_count) {
        valence_histogram[val]++;
    }
    
    std::cout << "Coordination number statistics:" << std::endl;
    for (const auto& pair : valence_histogram) {
        std::cout << "  Coord " << pair.first << ": " << pair.second << " nodes";
        if (pair.first == 5 || pair.first == 7) {
            std::cout << " ← Defect!";
        }
        std::cout << std::endl;
    }
}

void ConfigurationSaver::logDislocationData(
    double alpha,
    int num_dislocations
) {
    static bool first_write = true;
    std::ofstream file;
    
    // Open file in append mode (or create if first time)
    if (first_write) {
        file.open("dislocation_log.txt", std::ios::out);
        if (!file.is_open()) {
            std::cerr << "Error: Could not create dislocation_log.txt" << std::endl;
            return;
        }
        // Write header
        file << "# Alpha\tDislocations\n";
        first_write = false;
    } else {
        file.open("dislocation_log.txt", std::ios::app);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open dislocation_log.txt" << std::endl;
            return;
        }
    }
    
    // Write data
    file << std::fixed << std::setprecision(8);
    file << alpha << "\t" << num_dislocations << "\n";
    
    file.close();
}


double ConfigurationSaver::calculateTotalArea2D(UserData* userData) {
    // Perform null check
    if (!userData) {
        std::cerr << "Error: userData is null in calculateTotalArea2D" << std::endl;
        return 0.0;
    }
    
    // Extract needed values from userData
    std::vector<Point2D>& points = userData->points;
    std::vector<ElementTriangle2D>& elements = userData->elements;
    const std::vector<size_t>& active_elements = userData->active_elements;
    
    // Calculate total area by summing current areas of all active elements
    double total_area = 0.0;
    for (size_t elem_idx : active_elements) {
        if (elem_idx < elements.size()) {
            total_area += elements[elem_idx].calculateCurrentArea(points);
        }
    }
    
    return total_area;
}

double ConfigurationSaver::calculateTotalReferenceArea2D(UserData* userData) {
    if (!userData) {
        std::cerr << "Error: userData is null in calculateTotalReferenceArea2D" << std::endl;
        return 0.0;
    }
    
    std::vector<ElementTriangle2D>& elements = userData->elements;
    const std::vector<size_t>& active_elements = userData->active_elements;
    
    double total_area = 0.0;
    for (size_t elem_idx : active_elements) {
        if (elem_idx < elements.size()) {
            total_area += elements[elem_idx].getReferenceArea();
        }
    }
    
    return total_area;
}


ConfigurationSaver::IterationData ConfigurationSaver::readIterationData(int iteration) {
    IterationData data;
    data.iteration = iteration;
    
    // Construct filename for points file
    std::stringstream points_filename;
    points_filename << "triangle_data/points_" << std::setw(5) << std::setfill('0') 
                    << 10000+iteration << ".dat";
    
    // Open points file
    std::ifstream points_file(points_filename.str());
    if (!points_file) {
        std::cerr << "Error: Could not open file " << points_filename.str() 
                  << " for reading." << std::endl;
        return data;
    }
    
    // ❌ REMOVE THIS - setprecision doesn't affect INPUT streams
    // points_file >> std::setprecision(16);
    
    // ✅ OPTIONAL: Set stream to handle scientific notation explicitly
    points_file >> std::scientific;
    
    std::string line;
    
    // Skip header lines (lines starting with #)
    while (std::getline(points_file, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        
        // Parse data line
        std::istringstream iss(line);
        int point_index, fixed_flag;
        double x, y;
        double F00, F01, F10, F11;
        
        if (iss >> point_index >> x >> y >> fixed_flag >> F00 >> F01 >> F10 >> F11) {
            // Store position
            data.positions.push_back(Eigen::Vector2d(x, y));
            
            // Store F_ext (same for all points, but we'll take the last one)
            data.F_ext << F00, F01,
                          F10, F11;
        } else {
            std::cerr << "Warning: Failed to parse line: " << line << std::endl;
        }
    }
    
    points_file.close();
    
    std::cout << "Read " << data.positions.size() << " points from iteration " 
              << iteration << std::endl;
    std::cout << "F_ext matrix:\n" << data.F_ext << std::endl;
    
    return data;
}
// Helper function to read multiple iterations
std::vector<ConfigurationSaver::IterationData> ConfigurationSaver::readIterationRange(int start_iter, int end_iter) {
    std::vector<IterationData> all_data;  // Inside the function, you can still use just IterationData
    
    for (int iter = start_iter; iter <= end_iter; ++iter) {
        all_data.push_back(readIterationData(iter));
    }
    
    return all_data;
}


void ConfigurationSaver::processIterationsInFolder(
    const std::string& folder_path,
    std::function<void(const IterationData&)> processFunction) {
    
    // Check if directory exists
    if (!std::filesystem::exists(folder_path)) {
        std::cerr << "Error: Directory " << folder_path << " does not exist." << std::endl;
        return;
    }
    
    // Map to store iteration number -> iteration for sorting
    std::map<int, int> iteration_list;
    
    // Scan directory to find all iterations
    for (const auto& entry : std::filesystem::directory_iterator(folder_path)) {
        if (entry.is_regular_file()) {
            std::string filename = entry.path().filename().string();
            
            // Check if it's a points file (points_XXXXX.dat)
            if (filename.find("points_") == 0 && filename.find(".dat") != std::string::npos) {
                // Extract iteration number
                std::string num_str = filename.substr(7, 5);
                int file_num = std::stoi(num_str);
                int iteration = file_num - 10000;
                
                iteration_list[iteration] = iteration;
            }
        }
    }
    
    // Process each iteration one at a time
    for (const auto& [iteration, _] : iteration_list) {
        std::cout << "Processing iteration " << iteration << std::endl;
        
        // Read this iteration's data
        IterationData data = readIterationData(iteration);
        
        // Process it with the provided function
        processFunction(data);
        
        // data goes out of scope here and memory is freed
    }
    
    std::cout << "Processed " << iteration_list.size() << " iterations total." << std::endl;
}



void ConfigurationSaver::saveElements(const std::vector<ElementTriangle2D>& elements, 
                                      const std::vector<size_t>& active_elements,
                                      int iteration) {
    std::stringstream filename;
    filename << "triangle_data/elements_" << std::setw(5) << std::setfill('0') 
             << 10000 + iteration << ".dat";
    
    std::ofstream elem_file(filename.str());
    if (!elem_file) {
        std::cerr << "Error: Could not open " << filename.str() << std::endl;
        return;
    }
    
    // Write header
    elem_file << "# Element data: elem_idx node0 node1 node2 trans0_x trans0_y trans1_x trans1_y trans2_x trans2_y ref_area\n";
    
    // Set maximum precision for double
    elem_file << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10);
    
    // Save number of active elements
    elem_file << "# Active elements: " << active_elements.size() << "\n";
    
    // Save only active elements
    for (int idx : active_elements) {
        const ElementTriangle2D& elem = elements[idx];
        
        // Node indices (integers don't need precision)
        elem_file << idx << " "
                  << elem.getNodeIndex(0) << " "
                  << elem.getNodeIndex(1) << " "
                  << elem.getNodeIndex(2) << " ";
        
        // Translations
        for (int i = 0; i < 3; i++) {
            Eigen::Vector2d trans = elem.getTranslation(i);
            elem_file << trans.x() << " " << trans.y() << " ";
        }
        
        // Reference area
        elem_file << elem.getReferenceArea() << "\n";
    }
    
    elem_file.close();
    std::cout << "Saved " << active_elements.size() << " elements to " << filename.str() << std::endl;
}
std::pair<std::vector<ElementTriangle2D>, std::vector<size_t>>
ConfigurationSaver::loadElements(int iteration, 
                                 const Eigen::Matrix<double, 3, 2>& dndx,
                                 const std::vector<std::pair<int, int>>& full_mapping,
                                 const std::vector<Point2D>& reference_points) {
    
    std::stringstream filename;
    filename << "triangle_data/elements_" << std::setw(5) << std::setfill('0') 
             << 10000 + iteration << ".dat";
    
    std::ifstream elem_file(filename.str());
    if (!elem_file) {
        std::cerr << "Error: Could not open " << filename.str() << std::endl;
        return {{}, {}};
    }
    
    // ✅ ADD: Handle scientific notation explicitly
    elem_file >> std::scientific;
    
    std::vector<ElementTriangle2D> elements;
    std::vector<size_t> active_elements;
    
    std::string line;
    int num_active = 0;
    int line_number = 0;
    
    // Read file
    while (std::getline(elem_file, line)) {
        line_number++;
        
        if (line.empty() || line[0] == '#') {
            if (line.find("Active elements:") != std::string::npos) {
                std::istringstream iss(line);
                std::string dummy1, dummy2;
                iss >> dummy1 >> dummy1 >> dummy2 >> num_active;
                std::cout << "Expecting " << num_active << " active elements" << std::endl;
            }
            continue;
        }
        
        std::istringstream iss(line);
        size_t elem_idx;
        int n0, n1, n2;
        double tx0, ty0, tx1, ty1, tx2, ty2, ref_area;
        
        if (iss >> elem_idx >> n0 >> n1 >> n2 
                >> tx0 >> ty0 >> tx1 >> ty1 >> tx2 >> ty2 >> ref_area) {
            
            // ✅ ADD: Verify node indices are valid
            if (n0 >= reference_points.size() || n1 >= reference_points.size() || 
                n2 >= reference_points.size()) {
                std::cerr << "Warning: Invalid node indices at line " << line_number 
                          << " (n0=" << n0 << ", n1=" << n1 << ", n2=" << n2 
                          << ", max=" << reference_points.size() << ")" << std::endl;
                continue;
            }
            
            ElementTriangle2D elem;
            elem.setNodeIndex(0, n0);
            elem.setNodeIndex(1, n1);
            elem.setNodeIndex(2, n2);
            
            elem.setTranslation(0, Eigen::Vector2d(tx0, ty0));
            elem.setTranslation(1, Eigen::Vector2d(tx1, ty1));
            elem.setTranslation(2, Eigen::Vector2d(tx2, ty2));
            
            elem.set_shape_derivatives(dndx);
            elem.set_reference_mesh(reference_points);
            elem.setReferenceArea(ref_area);
            
            // Ensure elements vector is large enough
            while (elements.size() <= elem_idx) {
                elements.push_back(ElementTriangle2D());
            }
            
            elements[elem_idx] = elem;
            active_elements.push_back(elem_idx);
            
        } else {
            std::cerr << "Warning: Failed to parse line " << line_number 
                      << ": " << line << std::endl;
        }
    }
    
    // Set DOF mapping for all elements
    for (auto &element : elements) {
        element.set_dof_mapping(full_mapping);
    }
    
    elem_file.close();
    
    std::cout << "Loaded " << active_elements.size() << " elements from " 
              << filename.str() << std::endl;
    
    // ✅ ADD: Verify we got expected number
    if (num_active > 0 && active_elements.size() != static_cast<size_t>(num_active)) {
        std::cerr << "Warning: Expected " << num_active << " elements but loaded " 
                  << active_elements.size() << std::endl;
    }
    
    // ✅ ADD: Print some statistics for verification
    if (!active_elements.empty()) {
        std::cout << "  First element index: " << active_elements[0] << std::endl;
        std::cout << "  Last element index: " << active_elements.back() << std::endl;
        std::cout << "  Reference area (first): " << std::scientific 
                  << std::setprecision(std::numeric_limits<double>::max_digits10)
                  << elements[active_elements[0]].getReferenceArea() << std::endl;
    }
    
    return {elements, active_elements};
}