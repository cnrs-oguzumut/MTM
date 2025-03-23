#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <tuple>
#include <random>


#include "../include/geometry/Point2D.h"
#include "../include/geometry/DomainDimensions.h"
#include "../include/geometry/LatticeGenerator.h"
#include "../include/geometry/DomainInfo.h"
#include "../include/mesh/MeshGenerator.h"
#include "../include/optimization/LatticeOptimizer.h"
#include "../include/optimization/LBFGSOptimizer.h"



#include "../include/interatomic/inter_atomic.h"

#include "../include/mesh/Triangle.h"
#include "../include/mesh/ElementTriangle2D.h"
#include "../include/mesh/MeshGenerator.h"
#include "../include/lattice_energy/TriangularLatticeCalculator.h"
#include "../include/lattice_energy/SquareLatticeCalculator.h"
#include "../include/lattice_energy/Strain_Energy_LatticeCalculator.h"

#include "../include/lattice_energy/EnergyFunctions.h"
#include "../include/lattice_energy/EnergyFunctions.h"
#include "../include/output/triangulation_to_vtk.h"
#include "../include/output/configuration_saver.h"
#include "../include/output/ChangeMeasures.h"




/*
void example_1_atomistic_square() {
    // Parameters for lattice
    int nx = 30;
    int ny = 30;
    std::function<double(double)> potential_func = square_energy;
    std::function<double(double)> potential_func_der = square_energy_der;
    std::string lattice_type = "square"; // "square" or "triangular"

    std::cout << "STEP 1: Finding optimal lattice parameter...\n";
    //double optimal_lattice_parameter = find_optimal_lattice_parameter(potential_func,lattice_type);
    double optimal_lattice_parameter = 1.06585;
    double lattice_constant = optimal_lattice_parameter;
    
    std::vector<Point2D> square_points = LatticeGenerator::generate_2d_lattice(
        nx, ny, lattice_constant, lattice_type);
    
    int original_domain_size = square_points.size();
    DomainInfo domain_size = compute_domain_size(square_points);
    
    const std::array<double, 2> offsets = {lattice_constant, lattice_constant};
    DomainDimensions domain_dims(domain_size.get_width() , domain_size.get_height());
    
    int pbc = 0;
    std::vector<Triangle> triangulation;
    std::vector<Point2D> points_used_in_triangulation;
    if(pbc==1){
        // Generate periodic copies
        std::vector<Point2D> square_points_periodic = LatticeGenerator::create_periodic_copies(
            square_points,
            domain_dims,
            offsets);
        // Create Delaunay triangulation
        std::vector<Triangle> triangulation = MeshGenerator::createTrianglesFromPoints(
        square_points_periodic);
        points_used_in_triangulation = square_points_periodic;
    }
    else{
        triangulation = MeshGenerator::createTrianglesFromPoints(square_points);
        points_used_in_triangulation = square_points;
    }

    // Create domain maps
    auto [original_domain_map, translation_map] = MeshGenerator::create_domain_maps(
        original_domain_size, domain_dims, offsets);

    // Select unique connected triangles
    std::vector<Triangle> unique_triangles = MeshGenerator::select_unique_connected_triangles(
        points_used_in_triangulation,
        triangulation,
        original_domain_map,
        square_points.size(), // Original domain size
        1e-6 // Minimum Jacobian threshold
    );

    // Create element triangles
    std::vector<ElementTriangle2D> elements = MeshGenerator::createElementTri2D(
        unique_triangles,
        points_used_in_triangulation,
        original_domain_map,
        translation_map
    );

    // Calculate shape derivatives
    for (auto& element : elements) {
        element.calculate_shape_derivatives(points_used_in_triangulation);
    }
    
    std::cout << "Created " << elements.size() << " element triangles" << std::endl;
    
    // Boundary conditions are implemented here 
    auto [interior_mapping, full_mapping] = create_dof_mapping_original(points_used_in_triangulation, 0.5*lattice_constant, 0);

    // Create calculator for energy calculations
    SquareLatticeCalculator calculator(optimal_lattice_parameter);
    Eigen::Matrix2d C_I = Eigen::Matrix2d::Identity();
    double zero = calculator.calculate_energy(C_I, potential_func, 0);

    // Set up alpha values for deformation steps
    double alpha_min = 0.0;   // Minimum alpha value
    double alpha_max = 1.0;   // Maximum alpha value
    int num_alpha_points = 1000; // Number of alpha values to test
    
    // Generate evenly spaced alpha values
    std::vector<double> alpha_values;
    alpha_values.reserve(num_alpha_points);
    double alpha_step = (alpha_max - alpha_min) / (num_alpha_points - 1);
    
    for (int i = 0; i < num_alpha_points; i++) {
        alpha_values.push_back(alpha_min + i * alpha_step);
    }

    // Precompute active elements
    const std::vector<size_t>& active_elements = 
        initialize_active_elements(elements, full_mapping, points_used_in_triangulation.size());

    // Process each alpha value
    for (size_t i = 0; i < alpha_values.size(); i++) {
        double alpha = alpha_values[i];
        std::cout << "\n=== Processing alpha = " << alpha << " ===" << std::endl;
        
        // Create deformation gradient
        Eigen::Vector2d n1(0., 1.0);
        Eigen::Vector2d a1 = Eigen::Vector2d(1.0, 0.0);
        Eigen::Matrix2d F_ext = Eigen::Matrix2d::Identity() + alpha * a1 * n1.transpose();
        
        // Create user data with current alpha's deformation gradient
        UserData userData(
            points_used_in_triangulation,
            elements,
            calculator,
            potential_func,
            potential_func_der,
            zero,
            optimal_lattice_parameter,
            F_ext, 
            interior_mapping, 
            full_mapping, 
            active_elements
        );
        
        // Prepare initial point
        alglib::real_1d_array x;
        int n_vars = interior_mapping.size();
        x.setlength(2*n_vars);
        
        // Apply boundary conditions
        deform_boundary_nodes(points_used_in_triangulation, full_mapping, F_ext);
        
        // Map points to solver array
        map_points_to_solver_array(x, points_used_in_triangulation, interior_mapping, n_vars);
        
        // Set up minimization parameters
        alglib::minlbfgsstate state;
        alglib::minlbfgsreport rep;
        double epsg = 0;
        double epsf = 0;
        double epsx = 0;
        alglib::ae_int_t maxits = 0;
        
        // Create and configure optimizer
        alglib::minlbfgscreate(10, x, state);
        alglib::minlbfgssetcond(state, epsg, epsf, epsx, maxits);
        
        // Run optimization
        alglib::minlbfgsoptimize(state, minimize_energy_with_triangles, nullptr, &userData);
        
        // Get results
        alglib::minlbfgsresults(state, x, rep);
        
        std::cout << "Optimization completed: iterations=" << rep.iterationscount 
                  << ", termination type=" << rep.terminationtype << std::endl;
        
        // Update points with final optimized positions
        map_solver_array_to_points(x, points_used_in_triangulation, interior_mapping, n_vars);
        
        // Save configuration
        saveConfigurationToXY(points_used_in_triangulation, i);
        
        std::cout << "Iteration " << i << " completed successfully" << std::endl;
    }
}
*/


#include <iostream>
#include <fstream>
#include <iomanip>

void debug_deformation_tests() {
    // Create calculator with scale = 1.0
    Strain_Energy_LatticeCalculator calculator(1.0);
    
    // Identity matrix as reference state
    Eigen::Matrix2d F_I = Eigen::Matrix2d::Identity();
    Eigen::Matrix2d C_I = F_I.transpose() * F_I; // C = F^T * F
    
    // Calculate reference energy
    double zero = calculator.calculate_energy(C_I, nullptr, 0);
    std::cout << "Reference energy at identity: " << zero << std::endl;
    
    // 1. Simple Shear Test (F = [1 alpha; 0 1])
    {
        std::ofstream outFile("simple_shear_F_debug.txt");
        
        if (!outFile.is_open()) {
            std::cerr << "Error: Could not open simple_shear_F_debug.txt for writing!" << std::endl;
        } else {
            // Write header
            outFile << "# Simple Shear Energy Debug (using F matrix)" << std::endl;
            outFile << "# alpha, C11, C12, C22, C11_reduced, C12_reduced, C22_reduced, raw_energy, energy_minus_zero" << std::endl;
            
            // Calculate energy for different alpha values
            for (double alpha = 0.0; alpha <= 1.0; alpha += 0.01) {
                // Create deformation gradient for simple shear: F = [1 alpha; 0 1]
                Eigen::Matrix2d F;
                F << 1.0, alpha,
                     0.0, 1.0;
                
                // Calculate metric tensor C = F^T * F
                Eigen::Matrix2d C = F.transpose() * F;
                
                // Apply Lagrange reduction
                Eigen::Matrix2d C_reduced = calculator.lagrange_reduction(C);
                
                // Calculate energy
                double raw_energy = energy_functions::phi_func(C_reduced);
                double energy = raw_energy - zero;
                
                // Write to file
                outFile << std::fixed << std::setprecision(6)
                        << alpha << ", " 
                        << C(0, 0) << ", " 
                        << C(0, 1) << ", " 
                        << C(1, 1) << ", "
                        << C_reduced(0, 0) << ", " 
                        << C_reduced(0, 1) << ", " 
                        << C_reduced(1, 1) << ", "
                        << raw_energy << ", " 
                        << energy << std::endl;
            }
            
            outFile.close();
            std::cout << "Debug data written to simple_shear_F_debug.txt" << std::endl;
        }
    }
    
    // 2. Volume Change Test (F = [lambda 0; 0 lambda])
    {
        std::ofstream outFile("volume_change_F_debug.txt");
        
        if (!outFile.is_open()) {
            std::cerr << "Error: Could not open volume_change_F_debug.txt for writing!" << std::endl;
        } else {
            // Write header
            outFile << "# Volume Change Energy Debug (using F matrix)" << std::endl;
            outFile << "# lambda, volume_ratio, C11, C22, C11_reduced, C12_reduced, C22_reduced, raw_energy, energy_minus_zero" << std::endl;
            
            // Calculate energy for different lambda values
            for (double lambda = 0.8; lambda <= 1.2; lambda += 0.01) {
                // Create deformation gradient for isotropic volume change: F = [lambda 0; 0 lambda]
                Eigen::Matrix2d F;
                F << lambda, 0.0,
                     0.0, lambda;
                
                // Calculate metric tensor C = F^T * F
                Eigen::Matrix2d C = F.transpose() * F;
                
                // Calculate volume ratio (det(F) in 2D)
                double volume_ratio = F.determinant();
                
                // Apply Lagrange reduction
                Eigen::Matrix2d C_reduced = calculator.lagrange_reduction(C);
                
                // Calculate energy
                double raw_energy = energy_functions::phi_func(C_reduced);
                double energy = raw_energy - zero;
                
                // Write to file
                outFile << std::fixed << std::setprecision(6)
                        << lambda << ", " 
                        << volume_ratio << ", "
                        << C(0, 0) << ", " 
                        << C(1, 1) << ", "
                        << C_reduced(0, 0) << ", " 
                        << C_reduced(0, 1) << ", " 
                        << C_reduced(1, 1) << ", "
                        << raw_energy << ", " 
                        << energy << std::endl;
            }
            
            outFile.close();
            std::cout << "Debug data written to volume_change_F_debug.txt" << std::endl;
        }
    }
    
    // 3. Uniaxial Strain Test (F = [lambda 0; 0 1])
    {
        std::ofstream outFile("uniaxial_strain_F_debug.txt");
        
        if (!outFile.is_open()) {
            std::cerr << "Error: Could not open uniaxial_strain_F_debug.txt for writing!" << std::endl;
        } else {
            // Write header
            outFile << "# Uniaxial Strain Energy Debug (using F matrix)" << std::endl;
            outFile << "# lambda, C11, C22, C11_reduced, C12_reduced, C22_reduced, raw_energy, energy_minus_zero" << std::endl;
            
            // Calculate energy for different lambda values
            for (double lambda = 0.8; lambda <= 1.2; lambda += 0.01) {
                // Create deformation gradient for uniaxial strain: F = [lambda 0; 0 1]
                Eigen::Matrix2d F;
                F << lambda, 0.0,
                     0.0, 1.0;
                
                // Calculate metric tensor C = F^T * F
                Eigen::Matrix2d C = F.transpose() * F;
                
                // Apply Lagrange reduction
                Eigen::Matrix2d C_reduced = calculator.lagrange_reduction(C);
                
                // Calculate energy
                double raw_energy = energy_functions::phi_func(C_reduced);
                double energy = raw_energy - zero;
                
                // Write to file
                outFile << std::fixed << std::setprecision(6)
                        << lambda << ", " 
                        << C(0, 0) << ", " 
                        << C(1, 1) << ", "
                        << C_reduced(0, 0) << ", " 
                        << C_reduced(0, 1) << ", " 
                        << C_reduced(1, 1) << ", "
                        << raw_energy << ", " 
                        << energy << std::endl;
            }
            
            outFile.close();
            std::cout << "Debug data written to uniaxial_strain_F_debug.txt" << std::endl;
        }
    }
}

void example_1_conti_zanzotto() {
    // Parameters for lattice
    int nx = 100;
    int ny = 100;
    std::string lattice_type = "square"; // "square" or "triangular"
    
    // Irrelevant in this example
    std::function<double(double)> potential_func = square_energy;
    std::function<double(double)> potential_func_der = square_energy_der;

    std::cout << "STEP 1: Finding optimal lattice parameter...\n";
    //double optimal_lattice_parameter = find_optimal_lattice_parameter(potential_func,lattice_type);
    double optimal_lattice_parameter = 1.;
    double lattice_constant = optimal_lattice_parameter;
    
    std::vector<Point2D> square_points = LatticeGenerator::generate_2d_lattice(
        nx, ny, lattice_constant, lattice_type);
    
    int original_domain_size = square_points.size();
    DomainInfo domain_size = compute_domain_size(square_points);
    
    const std::array<double, 2> offsets = {lattice_constant, lattice_constant};
    std::cout<<"offsets: "<<offsets[0]<<" "<<offsets[1]<<std::endl;
    DomainDimensions domain_dims(domain_size.get_width() , domain_size.get_height());
    std::cout<<"domain_size.get_width() : "<<domain_size.get_width() <<std::endl;
    std::cout<<"domain_size.get_height() : "<<domain_size.get_height() <<std::endl;
    
    int pbc = 1;
    std::vector<Triangle> triangulation;
    std::vector<Point2D> points_used_in_triangulation;
    std::vector<Point2D> original_points;


    if(pbc==1){

        Eigen::Matrix2d F_ext;
        F_ext << 1.0, 0.,
                 0.0, 1.0;        
        // Generate periodic copies
        std::vector<Point2D> square_points_periodic = LatticeGenerator::create_periodic_copies(
            square_points,
            domain_dims,
            offsets,F_ext);
        // Create Delaunay triangulation
        triangulation = MeshGenerator::createTrianglesFromPoints(
        square_points_periodic);
        //saveConfigurationToXY(original_points, 999);
        points_used_in_triangulation = square_points_periodic;
    }
    else{
        triangulation = MeshGenerator::createTrianglesFromPoints(square_points);
        points_used_in_triangulation = square_points;
    }

    // Create domain maps
    auto [original_domain_map, translation_map] = MeshGenerator::create_domain_maps(
        original_domain_size, domain_dims, offsets);

    // Select unique connected triangles:
    // Here,  we use points_used_in_triangulation instead of square_points
    // we can calculate correctly jacobian determinant (triangles are fully constructed in the
    // periodic domain) 
    std::vector<Triangle> unique_triangles = MeshGenerator::select_unique_connected_triangles(
        points_used_in_triangulation,
        triangulation,
        original_domain_map,
        square_points.size(), // Original domain size
        1e-6 // Minimum Jacobian threshold
    );
    
    //write_triangulation_to_vtk(unique_triangles, points_used_in_triangulation, 999);



    // Create finite element triangles
    // Here,  we use square_points instead of points_used_in_triangulation
    // The triangles are constructed in the original domain
    // Using original_domain_map and translation_map
    // No need to clean triangles here; it is already done in select_unique_connected_triangles

    std::vector<ElementTriangle2D> elements = MeshGenerator::createElementTri2D(
        unique_triangles,
        square_points,
        original_domain_map,
        translation_map
    );

    // Calculate shape derivatives
    for (auto& element : elements) {
        element.calculate_shape_derivatives(square_points);
    }
    
    std::cout << "Created " << elements.size() << " element triangles" << std::endl;
    
    // Boundary conditions are implemented here 
    auto [interior_mapping, full_mapping] = create_dof_mapping_original(square_points, 0.5*lattice_constant, pbc);
    std::cout<<"interior_mapping.size(): "<<interior_mapping.size()<<std::endl;
    std::cout<<"full_mapping.size(): "<<full_mapping.size()<<std::endl;
    
    // Create calculator for energy calculations
    Strain_Energy_LatticeCalculator calculator(1.0);
    Eigen::Matrix2d C_I = Eigen::Matrix2d::Identity();
    // Using the analytical version that takes a function of metric tensor
    double zero = calculator.calculate_energy(C_I, potential_func, 0);
    std::cout<<"debugging simple shear test"<<std::endl;
    //debug_deformation_tests();
    
    std::cout<<"zero enegy value: "<<zero<<std::endl;


    // Set up alpha values for deformation steps
 // Set range boundaries
    double alpha_min = 0.13;  // Minimum alpha value
    double alpha_max = 1.0;   // Maximum alpha value

    // Set desired step size
    double step_size = 0.0001;

    // Calculate required number of points based on step size
    int num_alpha_points = static_cast<int>((alpha_max - alpha_min) / step_size) + 1;

    // Generate evenly spaced alpha values
    std::vector<double> alpha_values;
    alpha_values.reserve(num_alpha_points);

    for (int i = 0; i < num_alpha_points; i++) {
        alpha_values.push_back(alpha_min + i * step_size);
    }

    // Precompute active elements
    std::vector<size_t> active_elements = 
        initialize_active_elements(elements, full_mapping, square_points.size());
    std::cout<<"active_elements.size(): "<<active_elements.size()<<std::endl;
    // Process each alpha value
    for (size_t i = 0; i < alpha_values.size(); i++) {
        double alpha = alpha_values[i];
        std::cout << "\n=== Processing alpha = " << alpha << " ===" << std::endl;
        
        // Create deformation gradient
        Eigen::Vector2d n1(0., 1.0);
        Eigen::Vector2d a1 = Eigen::Vector2d(1.0, 0.0);
        Eigen::Matrix2d F_ext;
        F_ext << 1.0, alpha,
                 0.0, 1.0;   
                 
        if(i==0){
            // Create random number generator
            std::random_device rd;
            std::mt19937 gen(rd());
            
            // Create Gaussian distribution with small standard deviation
            double noise_level = 0.05; // Adjust this value for more/less noise
            std::normal_distribution<double> noise_dist(0.0, noise_level);
            
            for (const auto& pair : full_mapping) {
                int original_idx = pair.first;
                
                // Add noise to the coordinates
                Eigen::Vector2d noise(noise_dist(gen), noise_dist(gen));

                
                // Apply deformation with noise: x_deformed = F·x + noise
                square_points[original_idx].coord = F_ext * square_points[original_idx].coord + noise;
            }
        }            

        
        // Create user data with current alpha's deformation gradient
        UserData userData(
            square_points,
            elements,
            calculator,
            potential_func,
            potential_func_der,
            zero,
            optimal_lattice_parameter,
            F_ext, 
            interior_mapping, 
            full_mapping, 
            active_elements
        );
        
        // Prepare initial point
        alglib::real_1d_array x;
        int n_vars = interior_mapping.size();
        x.setlength(2*n_vars);
        double total_energy2;

        // Apply boundary conditions
        //deform_boundary_nodes(square_points, full_mapping, F_ext);
        
        // Map points to solver array
        map_points_to_solver_array(x, square_points, interior_mapping, n_vars);

        double epsg = 0;
        double epsf = 0;
        double epsx = 0;
        alglib::ae_int_t maxits = 0;

        LBFGSOptimizer optimizer(
            10,     // corrections - history size for L-BFGS
            epsg,   // epsg - gradient norm stopping condition
            epsf,   // epsf - function value stopping condition
            epsx,   // epsx - argument changes stopping condition
            maxits    // maxits - maximum iterations
        );
       
        // Check if remeshing is necessary based on point changes
        alglib::real_1d_array original_x;
        original_x.setlength(x.length());
        // Store original positions before mapping back (if you need to keep them)
        for (int i = 0; i < x.length(); i++) {
            original_x[i] = x[i];
        }

        optimizer.optimize(x, minimize_energy_with_triangles, &userData);
        
        
        // Update points with final optimized positions
        map_solver_array_to_points(x, square_points, interior_mapping, n_vars);
        // Calculate change measures between original and optimized positions
        ChangeMeasures result = computeChangeMeasures(x, original_x);
        std::cout<<"relative_change: "<<result.relative_change<<std::endl;

        if (result.relative_change > 2.0){
            std::cout<<"REMESHING STARTS"<<std::endl;
            // Generate periodic copies
            std::vector<Point2D> new_square_points_periodic = LatticeGenerator::create_periodic_copies(
            square_points,
            domain_dims,
            offsets,F_ext);
            // Create Delaunay triangulation
            triangulation = MeshGenerator::createTrianglesFromPoints(
            new_square_points_periodic);
            //saveConfigurationToXY(original_points, 999);
            points_used_in_triangulation = new_square_points_periodic;

            unique_triangles = MeshGenerator::select_unique_connected_triangles(
            points_used_in_triangulation,
            triangulation,
            original_domain_map,
            square_points.size(), // Original domain size
            1e-6 // Minimum Jacobian threshold
            );

            elements = MeshGenerator::createElementTri2D(
            unique_triangles,
            square_points,
            original_domain_map,
            translation_map
            );

            const std::vector<size_t> new_active_elements = 
            initialize_active_elements(elements, full_mapping, square_points.size());
            active_elements = new_active_elements;  // Copy the new elements to the original vector
            // Calculate shape derivatives
            for (auto& element : elements) {
            element.calculate_shape_derivatives(square_points);
            }

            optimizer.optimize(x, minimize_energy_with_triangles, &userData);
            map_solver_array_to_points(x, square_points, interior_mapping, n_vars);

        }
        
        
        // Save configuration
        if(i%5==0){
            //saveConfigurationToXY(square_points, i);
            double total_energy = 0;
            saveConfigurationWithStressAndEnergy2D(&userData, i, total_energy);
            //write_triangulation_to_vtk(unique_triangles, points_used_in_triangulation, i);
        }
        
        std::cout << "Iteration " << i << " completed successfully" << std::endl;
    }
}

int main() {
    example_1_conti_zanzotto();
    return 0;
}