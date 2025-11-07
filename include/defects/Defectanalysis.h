#ifndef DEFECT_ANALYSIS_H
#define DEFECT_ANALYSIS_H

#include <vector>
#include <array>
#include <Eigen/Dense>

#include "../include/geometry/Point2D.h"


#include "../include/mesh/Triangle.h"
#include "../include/mesh/ElementTriangle2D.h"
#include "../include/mesh/MeshGenerator.h"
#include "../include/mesh/Remesher.h"

#include "../include/output/configuration_saver.h"
#include "../include/optimization/LatticeOptimizer.h" //for userdata




#include "src/optimization.h" 
// This path is relative to ALGLIB_DIR
// Forward declarations
class Point2D;
class ElementTriangle2D;
class UserData;
class AdaptiveMesher;

/**
 * @class DefectAnalysis
 * @brief Handles defect analysis in reference configuration
 * 
 * This class provides utilities to analyze defects (dislocations, vacancies, etc.)
 * by transforming the current configuration to reference configuration and remeshing.
 */
class DefectAnalysis {
public:
    /**
     * @brief Analyze defects in the reference configuration
     * 
     * This function:
     * 1. Copies current points and transforms to reference config using F_ext^(-1)
     * 2. Remeshes in the reference configuration
     * 3. Calculates coordination numbers
     * 4. Saves VTK output for visualization
     * 5. Optionally counts 5/7 pairs for dislocation detection
     * 
     * @param userData Current simulation state (remains unchanged)
     * @param file_id Output file identifier
     * @param dndx Shape function derivatives
     * @param offsets Periodic offsets for domain
     * @param original_domain_map Mapping of original domain points
     * @param translation_map Translation vectors for periodic copies
     * @param domain_dims_point Domain dimensions
     * @param element_area Reference element area
     * @param pbc Whether to use periodic boundary conditions
     * @param reduction Whether to apply Lagrange reduction
     * @return Number of 5/7 dislocation pairs detected
     */
    static std::pair<int, std::vector<int>> analyzeDefectsInReferenceConfig(
        const UserData* userData,
        int file_id,
        const Eigen::Matrix<double, 3, 2>& dndx,
        const std::array<double, 2>& offsets,
        const std::vector<int>& original_domain_map,
        const std::vector<std::tuple<double, double>>& translation_map,
        const Point2D& domain_dims_point,
        double element_area,
        bool pbc,
        bool reduction = true
    );
    /**
     * @brief Count 5/7 coordination pairs in the mesh
     * 
     * Identifies dislocation dipoles by finding adjacent nodes with
     * coordination numbers 5 and 7.
     * 
     * @param points Point positions
     * @param elements Mesh elements
     * @param active_elements Active element indices
     * @return Number of 5/7 pairs found
     */
    static int count57Pairs(
        const std::vector<Point2D>& points,
        const std::vector<ElementTriangle2D>& elements,
        const std::vector<size_t>& active_elements
    );

    /**
     * @brief Calculate coordination number for each node
     * 
     * @param points Point positions
     * @param elements Mesh elements
     * @param active_elements Active element indices
     * @return Vector of coordination numbers indexed by node
     */
    static std::vector<int> calculateCoordinationNumbers(
        const std::vector<Point2D>& points,
        const std::vector<ElementTriangle2D>& elements,
        const std::vector<size_t>& active_elements
    );

private:
    /**
     * @brief Helper function to transform points to reference configuration
     * 
     * @param points Points to transform (modified in place)
     * @param F_ext_inv Inverse of external deformation gradient
     */
    static void transformToReferenceConfig(
        std::vector<Point2D>& points,
        const Eigen::Matrix2d& F_ext_inv
    );
};

#endif // DEFECT_ANALYSIS_H