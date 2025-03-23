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

// Function to save configuration with nodal stress and energy values to XYZ file for OVITO (2D version)
void saveConfigurationWithStressAndEnergy2D(
    UserData* userData,
    int iteration, 
    double total_energy
);

#endif // CONFIGURATION_SAVER_2D_H