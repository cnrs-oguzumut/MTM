#ifndef CONFIGURATION_SAVER_H
#define CONFIGURATION_SAVER_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <filesystem>
#include <map>
#include <sstream>
#include <array>
#include <Eigen/Dense>
#include <functional>


// Include the header that defines UserData
#include "../include/optimization/LatticeOptimizer.h"




// Forward declarations for types used in function signatures
class DomainDimensions;

/**
 * ConfigurationSaver class for managing energy and stress calculation and saving
 * configuration data to files for visualization and analysis.
 */
class ConfigurationSaver {
public:
    

    struct IterationData {
        int iteration;
        std::vector<Eigen::Vector2d> positions;
        Eigen::Matrix2d F_ext;
    };

    /**
     * Calculate energy and stress from lattice configuration without saving to file
     * 
     * @param userData Pointer to user data containing lattice configuration
     * @param total_energy Reference to store calculated total energy
     * @param total_stress Reference to store calculated total stress tensor
     * @param reduction Flag to enable/disable Lagrange reduction
     */
    static void calculateEnergyAndStress(
        UserData* userData, 
        double& total_energy, 
        Eigen::Matrix2d& total_stress, 
        bool reduction = false);
    
    /**
     * Save configuration with nodal stress and energy values to XYZ file for OVITO
     * 
     * @param userData Pointer to user data containing lattice configuration
     * @param iteration Current iteration number for filename
     * @param total_energy Reference to store calculated total energy
     * @param total_stress Reference to store calculated total stress (scalar value)
     * @param reduction Flag to enable/disable Lagrange reduction
     */
    static void saveConfigurationWithStressAndEnergy2D(
        UserData* userData,
        int iteration,
        double& total_energy,
        double& total_stress, 
        bool reduction = false);
    
    /**
     * Save triangle data including deformation gradients and node information
     * 
     * @param userData Pointer to user data containing lattice configuration
     * @param iteration Current iteration number for filename
     * @param domain_dims Domain dimensions structure
     * @param offsets Array containing x and y offsets
     * @param full_mapping Vector of pairs mapping nodes to boundary conditions
     */
    static void saveTriangleData(
        UserData* userData,
        int iteration,
        const DomainDimensions& domain_dims,
        const std::array<double, 2>& offsets,
        const std::vector<std::pair<int, int>>& full_mapping);
    
    /**
     * Log energy and stress values to CSV file for later analysis
     * 
     * @param iteration Current iteration number
     * @param alpha Current deformation parameter
     * @param pre_energy Energy before optimization
     * @param pre_stress Stress before optimization
     * @param post_energy Energy after optimization
     * @param post_stress Stress after optimization
     * @param plasticity_flag Flag indicating plasticity state
     */
    static void logEnergyAndStress(
        int iteration, 
        double alpha, 
        double pre_energy, 
        double pre_stress,
        double post_energy, 
        double post_stress,
        int plasticity_flag);

    static void logEnergyAndStress_v2(
    int iteration, 
    double alpha, 
    double pre_energy, 
    double pre_stress,
    double post_energy, 
    double post_stress,
    double pre_area,
    double post_area,
    bool shouldRemesh);  // Added plasticity flag parameter

    /**
     * Write configuration data to VTK format for visualization
     * 
     * @param points Vector of 2D points
     * @param elements Vector of triangular elements
     * @param userData Pointer to user data (const)
     * @param iteration Current iteration number
     * @param reduction Flag to enable/disable Lagrange reduction
     */
    static void writeToVTK(
        const std::vector<Point2D>& points,
        const std::vector<ElementTriangle2D>& elements,
        const UserData* userData,
        const int iteration,
        bool reduction = false,
        const std::vector<int>& coordination= std::vector<int>(),
        double load_strength=0);


    static void writeToVTK_DefectAnalysis(
        const std::vector<Point2D>& points,
        const std::vector<ElementTriangle2D>& elements,
        const UserData* userData,
        int iteration);

    static void logDislocationData(
        double alpha,
        int num_dislocations);

    /**
     * Calculate total current area of all active elements
     * 
     * @param userData Pointer to user data containing lattice configuration
     * @return Total current area (deformed configuration)
     */
    static double calculateTotalArea2D(UserData* userData);
    
    /**
     * Calculate total reference area of all active elements
     * 
     * @param userData Pointer to user data containing lattice configuration
     * @return Total reference area (undeformed configuration)
     */
    static double calculateTotalReferenceArea2D(UserData* userData);

    static std::vector<IterationData> readIterationRange(int start_iter, int end_iter);
    static IterationData readIterationData(int iteration);


    static void processIterationsInFolder(
        const std::string& folder_path,
        std::function<void(const IterationData&)> processFunction);


    static void saveElements(const std::vector<ElementTriangle2D>& elements, 
                  const std::vector<size_t>& active_elements,
                  int iteration) ;


    static  std::pair<std::vector<ElementTriangle2D>, std::vector<size_t>> 
            loadElements(int iteration,
             const Eigen::Matrix<double, 3, 2>& dndx,
             const std::vector<std::pair<int, int>>& full_mapping,
             const std::vector<Point2D>& reference_points);
        
    
private:
    /**
     * Create output directory if it doesn't exist
     * 
     * @param directory Directory path to create
     */
    static void ensureDirectoryExists(const std::string& directory);
    

};

#endif // CONFIGURATION_SAVER_H