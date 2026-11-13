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
    static void setReferenceRotationDeg(double deg);
    static double getReferenceRotationDeg();
    

    struct IterationData {
        int iteration;
        std::vector<Eigen::Vector2d> positions;
        Eigen::Matrix2d F_ext;
    };

    struct CheckpointData {
        int nx = 0;
        int ny = 0;
        int iteration = 0;
        double current_alpha = 0.0;
        double step_size = 6e-5;
        double alpha_end = 1.0;
        std::string mode = "positive";
        unsigned int seed = 50;
        bool enable_remeshing = true;
        Eigen::Matrix2d F_ext = Eigen::Matrix2d::Identity();
        std::vector<Point2D> points;
        std::vector<ElementTriangle2D> elements;
        std::vector<size_t> active_elements;
        // Reference positions X_0 (custom geometry: the Stage-1 re-anchored state). Optional
        // REF_NODES section; empty when absent (older checkpoints, periodic cells).
        std::vector<Eigen::Vector2d> ref_points;
        // PBC stress control: the cell after the alpha = 0 relaxation (origin of the conjugate strain).
        // Optional F_stage1 line; identity when absent.
        Eigen::Matrix2d F_stage1 = Eigen::Matrix2d::Identity();
        bool has_F_stage1 = false;
        // Written at a relaxation step that precedes the loading (Stage 1, relaxation at a prescribed shear
        // strain, zero-stress relaxation): not a load point to resume from. Optional relaxation_step line.
        bool relaxation_step = false;
        // Reference frame used when saving the checkpoint ("square" or "triangular"; missing line -> "square").
        std::string reference_frame_tag = "square";
    };

    static Eigen::Matrix2d alignDeformationGradient(const Eigen::Matrix2d& F, double theta_rad = 0.0);
    static Eigen::Matrix2d computeAlignmentMatrix(const Eigen::Matrix2d& F, double theta_rad = 0.0);
    static Eigen::Matrix2d computeElementGrainAlignedElasticF(
        const ElementTriangle2D& element,
        const std::vector<Point2D>& points,
        int* out_grain_id = nullptr,
        double* out_theta_rad = nullptr);

    // Lattice orientation of one element, measured from its CURRENT edge vectors (E = F * E_ref):
    // the least-distorted lattice cell spanned by the element is compared with the ideal cell
    // (equilateral or square), and the rotation of that fit is the lattice orientation. Independent of
    // the reference frame, the mesh history (remeshing) and the starting positions X_0.
    struct LatticeOrientationSample {
        bool valid = false;        // false for inverted / degenerate elements
        double phi_deg = 0.0;      // absolute lattice orientation = theta + omega (within sym/2 of the grain's start)
        double omega_deg = 0.0;    // rotation from the element's grain angle theta, wrapped to (-sym/2, sym/2]
        double theta_deg = 0.0;    // grain reference angle used
        int grain_id = 0;          // majority grain of the 3 nodes
        bool mixed_grain = false;  // nodes from 3 different grains (grain-boundary element)
        bool boundary_element = false; // nodes from more than one grain
        double area = 0.0;         // current element area
        double distortion = 0.0;   // equivalent shear of the best-fitting cell away from the ideal cell
        bool recognized = false;   // distortion <= kLatticeRecognitionShear (else core / transformed region)
    };
    // Cells sheared further than this from every ideal lattice cell (dislocation cores, phase-transformed
    // regions) have no meaningful lattice orientation and are left out of the rotation statistics.
    static constexpr double kLatticeRecognitionShear = 0.35;
    static double latticeSymmetryDeg();  // 60 (triangular) or 90 (square)
    static double wrapToSymmetry(double angle_deg, double symmetry_deg);
    static LatticeOrientationSample measureLatticeOrientation(const ElementTriangle2D& element,
                                                              const std::vector<Point2D>& points);
    // Per-step lattice rotation log (lattice_rotation_log.csv in the current directory): area-weighted
    // rotation statistics, per-grain mean rotations (unwrapped across steps) and the rigid rotation of
    // the applied macroscopic deformation. Called from logEnergyAndStress_v2.
    static void logLatticeRotation(int iteration, double alpha, UserData* userData);

    static void saveCheckpoint(const std::string& filename, const CheckpointData& data);
    static bool loadCheckpoint(const std::string& filename, CheckpointData& data,
                               const Eigen::Matrix<double, 3, 2>& dndx);

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
    bool shouldRemesh,
    UserData* postUserData = nullptr,
    const std::vector<int>* contact_atoms = nullptr,  // Added plasticity flag & extended stress metrics
    // Optional extra columns appended after ContactAtoms (force-controlled runs only). The header
    // is written when the file is created, so pass the same names on every call of a run.
    const std::vector<std::pair<std::string, double>>* extra_columns = nullptr);

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

    static void writeToVTKFile(
        const std::string& full_filepath,
        const std::vector<Point2D>& points,
        const std::vector<ElementTriangle2D>& elements,
        const UserData* userData,
        bool reduction = true,
        const std::vector<int>& coordination = std::vector<int>(),
        double load_strength = 0);


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