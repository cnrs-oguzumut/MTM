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
#include <Eigen/Dense>

// Include the header that defines UserData
#include "../include/optimization/LatticeOptimizer.h"

/**
 * ConfigurationSaver class for managing energy and stress calculation and saving
 * configuration data to files for visualization and analysis.
 */
class ConfigurationSaver {
public:
    /**
     * Calculate energy and stress from lattice configuration without saving to file
     * 
     * @param userData Pointer to user data containing lattice configuration
     * @param total_energy Reference to store calculated total energy
     * @param total_stress Reference to store calculated total stress
     */
    static void calculateEnergyAndStress(
        UserData* userData, 
        double& total_energy, 
        double& total_stress);
    
    /**
     * Save configuration with nodal stress and energy values to XYZ file for OVITO
     * 
     * @param userData Pointer to user data containing lattice configuration
     * @param iteration Current iteration number for filename
     * @param total_energy Reference to store calculated total energy
     * @param total_stress Reference to store calculated total stress
     */
    static void saveConfigurationWithStressAndEnergy2D(
        UserData* userData,
        int iteration,
        double& total_energy,
        double& total_stress);
    
    /**
     * Log energy and stress values to CSV file for later analysis
     * 
     * @param iteration Current iteration number
     * @param alpha Current deformation parameter
     * @param pre_energy Energy before optimization
     * @param pre_stress Stress before optimization
     * @param post_energy Energy after optimization
     * @param post_stress Stress after optimization
     */
    static void logEnergyAndStress(
        int iteration, 
        double alpha, 
        double pre_energy, 
        double pre_stress,
        double post_energy, 
        double post_stress,
        bool plasticity_flag);  // Added plasticity flag parameter

    static void writeToVTK(
        const std::vector<Point2D>& points,
        const std::vector<ElementTriangle2D>& elements,
        const UserData* userData,
        int iteration);
        
    
private:
    // Private method to create output directory if it doesn't exist
    static void ensureDirectoryExists(const std::string& directory);
};

#endif // CONFIGURATION_SAVER_H