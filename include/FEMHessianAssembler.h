#ifndef FEM_HESSIAN_ASSEMBLY_H
#define FEM_HESSIAN_ASSEMBLY_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <vector>
#include <functional>
#include "acoustic_tensor.h"
#include "../include/mesh/ElementTriangle2D.h"
#include "../include/lattice_energy/BaseLatticeCalculator.h"

/**
 * Structure to hold eigenvalue analysis results
 */
struct EigenResults {
    Eigen::VectorXd eigenvalues;      // N smallest eigenvalues
    Eigen::MatrixXd eigenvectors;     // Corresponding eigenvectors (columns)
    int num_computed;                  // Number of eigenvalues successfully computed
    
    EigenResults() : num_computed(0) {}
};

/**
 * Structure to hold participation ratio analysis results
 */
struct ParticipationAnalysis {
    std::vector<double> participation_ratios;  // P for each mode
    std::vector<bool> is_extended;             // true if mode is extended (phonon-like)
    std::vector<double> localization_lengths;  // Effective number of participating nodes
    double threshold;                           // Threshold used for classification
    
    ParticipationAnalysis() : threshold(0.5) {}
};

/**
 * Structure to hold density of states results
 */
struct DOSResults {
    std::vector<double> omega_bins;      // Frequency bins (center values)
    std::vector<double> dos;             // g(ω) - density of states
    std::vector<double> dos_cumulative;  // Cumulative DOS
    int num_bins;                        // Number of bins
    double omega_min;                    // Minimum frequency
    double omega_max;                    // Maximum frequency
    
    DOSResults() : num_bins(0), omega_min(0.0), omega_max(0.0) {}
};

class FEMHessianAssembler {
private:
    // Energy derivative parameters - THESE ARE ESSENTIAL!
    BaseLatticeCalculator* strain_calculator;
    std::function<double(double)> potential_func_der;
    std::function<double(double)> potential_func_sder;
    double normalisation;
    
    /**
     * Extract acoustic tensor ITensor to 4D array for easier manipulation
     * 
     * @param A_tensor The ITensor acoustic tensor with indices (r_, i_, s_, j_)
     * @param A Output 4D array A[i][K][j][L]
     */
    void extractAcousticTensor(
        const itensor::ITensor& A_tensor,
        double A[2][2][2][2]);

public:
    /**
     * Constructor
     */
    FEMHessianAssembler() : strain_calculator(nullptr), normalisation(0.0) {}
    
    /**
     * Set energy derivative parameters (must be called before assembleGlobalStiffness)
     * 
     * @param calc Pointer to lattice calculator for strain energy
     * @param dpot First derivative of potential function
     * @param d2pot Second derivative of potential function
     * @param norm Normalization factor
     */
    void setEnergyParameters(
        BaseLatticeCalculator* calc,
        const std::function<double(double)>& dpot,
        const std::function<double(double)>& d2pot,
        double norm)
    {
        strain_calculator = calc;
        potential_func_der = dpot;
        potential_func_sder = d2pot;
        normalisation = norm;
    }
    
    /**
     * Compute element stiffness matrix K^{ab}_{ij} = A_{iKjL} * (dN^a/dX_K) * (dN^b/dX_L)
     * 
     * @param acoustic_tensor The acoustic tensor object containing A_{iKjL}
     * @param element The triangular element with shape function derivatives
     * @param area_weight Integration weight (typically element area for single Gauss point)
     * @return 6x6 element stiffness matrix
     */
    Eigen::MatrixXd computeElementStiffness(
        AcousticTensor& acoustic_tensor,
        const ElementTriangle2D& element,
        double area_weight = 1.0);
    
    /**
     * Assemble global stiffness matrix from all elements
     * 
     * @param elements Vector of all triangular elements
     * @param current_points Current node positions for deformation gradient
     * @param num_total_dofs Total number of free DOFs in the system
     * @param dof_mapping DOF mapping (original_idx, solver_idx) for each node
     * @return Sparse global stiffness matrix
     */
    Eigen::SparseMatrix<double> assembleGlobalStiffness(
        std::vector<ElementTriangle2D>& elements,
        const std::vector<Point2D>& current_points,
        int num_total_dofs,
        const std::vector<std::pair<int, int>>& dof_mapping);
    
    /**
     * Compute N smallest eigenvalues and eigenvectors of the stiffness matrix
     * 
     * For large sparse matrices, this converts to dense for eigenvalue computation.
     * For very large systems (>10000 DOFs), consider using iterative methods instead.
     * 
     * @param K_global The global stiffness matrix (typically from assembleGlobalStiffness)
     * @param N Number of smallest eigenvalues to compute
     * @return EigenResults structure containing eigenvalues and eigenvectors
     * 
     * Note: The stiffness matrix may have zero eigenvalues corresponding to rigid body modes.
     *       These will appear as the smallest eigenvalues.
     */
    EigenResults computeSmallestEigenvalues(
        const Eigen::SparseMatrix<double>& K_global,
        int N);
    
    /**
     * Compute N smallest eigenvalues and eigenvectors using shift-invert mode
     * More efficient for large sparse matrices than full conversion to dense.
     * 
     * @param K_global The global stiffness matrix
     * @param N Number of smallest eigenvalues to compute
     * @param shift Shift value (typically small positive number like 1e-6)
     * @return EigenResults structure containing eigenvalues and eigenvectors
     * 
     * This method uses (K - σI)^(-1) to find eigenvalues near σ,
     * which is more efficient than computing all eigenvalues.
     */
    EigenResults computeSmallestEigenvaluesIterative(
        const Eigen::SparseMatrix<double>& K_global,
        int N,
        double shift = 1e-6);

    /**
     * Compute N smallest eigenvalues using Spectra library
     * 
     * @param K_global The global stiffness matrix
     * @param N Number of smallest eigenvalues to compute
     * @param shift Shift value (not used in current implementation)
     * @return EigenResults structure containing eigenvalues and eigenvectors
     */
    EigenResults computeSmallestEigenvaluesIterative_spectra(
        const Eigen::SparseMatrix<double>& K_global,
        int N,
        double shift = 1e-6);

    /**
     * Export eigenvector modes to VTK for visualization in Paraview
     * 
     * @param filename Base filename (without extension)
     * @param current_points Current mesh node positions
     * @param eigen_results Results from eigenvalue computation
     * @param mode_indices Which eigenmodes to export (e.g., {0, 1, 2} for first 3 modes)
     * @param scale_factor Scaling factor for visualizing displacement (default: 1.0)
     * @param dof_mapping DOF mapping to convert from solver indices to node indices
     * @param elements Mesh elements (for connectivity)
     */
    void exportEigenvectorsToVTK(
        const std::string& filename,
        const std::vector<Point2D>& current_points,
        const EigenResults& eigen_results,
        const std::vector<int>& mode_indices,
        double scale_factor = 1.0,
        const std::vector<std::pair<int, int>>& dof_mapping = {},
        const std::vector<ElementTriangle2D>& elements = {});

    /**
     * Export a single eigenmode to VTK
     * 
     * @param filename Output filename
     * @param current_points Current mesh node positions
     * @param eigenvector The eigenvector to visualize
     * @param eigenvalue Corresponding eigenvalue
     * @param mode_number Mode index (for labeling)
     * @param scale_factor Scaling factor for visualization
     * @param dof_mapping DOF mapping
     * @param elements Mesh elements (for connectivity)
     */
    void exportSingleModeToVTK(
        const std::string& filename,
        const std::vector<Point2D>& current_points,
        const Eigen::VectorXd& eigenvector,
        double eigenvalue,
        int mode_number,
        double scale_factor,
        const std::vector<std::pair<int, int>>& dof_mapping,
        const std::vector<ElementTriangle2D>& elements);

    /**
     * Compute participation ratio for a single eigenvector
     * 
     * The participation ratio P is defined as:
     * P = (Σ|u_i|²)² / (N * Σ|u_i|⁴)
     * 
     * where u_i is the displacement amplitude at node i, and N is the number of nodes.
     * 
     * - P ≈ 1: Extended mode (phonon) - all nodes participate equally
     * - P ≈ 0: Localized mode - only few nodes participate
     * 
     * @param eigenvector The eigenvector to analyze
     * @param dof_mapping DOF mapping to reconstruct full displacement field
     * @param num_nodes Total number of nodes in the mesh
     * @return Participation ratio in range [0, 1]
     */
    double computeParticipationRatio(
        const Eigen::VectorXd& eigenvector,
        const std::vector<std::pair<int, int>>& dof_mapping,
        int num_nodes);

    /**
     * Compute participation ratios for all eigenmodes
     * 
     * @param eigen_results Results from eigenvalue computation
     * @param dof_mapping DOF mapping
     * @param num_nodes Total number of nodes
     * @param threshold Threshold for classifying extended vs localized modes (default: 0.5)
     * @return ParticipationAnalysis structure with results
     */
    ParticipationAnalysis analyzeParticipationRatios(
        const EigenResults& eigen_results,
        const std::vector<std::pair<int, int>>& dof_mapping,
        int num_nodes,
        double threshold = 0.5);

    /**
     * Export participation ratio analysis to file
     * 
     * @param filename Output filename
     * @param analysis Participation analysis results
     * @param eigen_results Eigenvalue results (for eigenvalues)
     */
    void exportParticipationAnalysis(
        const std::string& filename,
        const ParticipationAnalysis& analysis,
        const EigenResults& eigen_results);

    /**
     * Print participation ratio summary to console
     * 
     * @param analysis Participation analysis results
     * @param eigen_results Eigenvalue results
     */
    void printParticipationSummary(
        const ParticipationAnalysis& analysis,
        const EigenResults& eigen_results);

    /**
     * Export non-rigid eigenvalues with derived quantities (ω, λ²)
     * 
     * @param filename Output filename
     * @param eigen_results Eigenvalue results
     * @param num_rigid_body Number of rigid body modes to skip
     */
    void exportNonRigidEigenvalues(
        const std::string& filename,
        const EigenResults& eigen_results,
        int num_rigid_body);

    /**
     * Export phonon modes (extended modes) with participation ratios
     * 
     * @param filename Output filename
     * @param eigen_results Eigenvalue results
     * @param participation Participation analysis results
     * @param num_rigid_body Number of rigid body modes to skip
     * @return Vector of phonon mode indices
     */
    std::vector<int> exportPhononModes(
        const std::string& filename,
        const EigenResults& eigen_results,
        const ParticipationAnalysis& participation,
        int num_rigid_body);

    /**
     * Export localized modes (non-phonon modes) with participation ratios
     * 
     * @param filename Output filename
     * @param eigen_results Eigenvalue results
     * @param participation Participation analysis results
     * @param num_rigid_body Number of rigid body modes to skip
     * @return Vector of localized mode indices
     */
    std::vector<int> exportLocalizedModes(
        const std::string& filename,
        const EigenResults& eigen_results,
        const ParticipationAnalysis& participation,
        int num_rigid_body);

    /**
     * Complete eigenmode analysis and export
     * Performs participation analysis, exports all data files, and VTK files
     * 
     * @param base_filename Base name for output files
     * @param eigen_results Eigenvalue results (sorted)
     * @param current_points Current mesh positions
     * @param dof_mapping DOF mapping
     * @param elements Mesh elements
     * @param num_rigid_body Number of rigid body modes to skip
     * @param participation_threshold Threshold for phonon classification (default: 0.5)
     * @param vtk_scale_factor Scale factor for VTK visualization (default: 0.1)
     */
    void exportCompleteEigenmodeAnalysis(
        const std::string& base_filename,
        const EigenResults& eigen_results,
        const std::vector<Point2D>& current_points,
        const std::vector<std::pair<int, int>>& dof_mapping,
        const std::vector<ElementTriangle2D>& elements,
        int num_rigid_body,
        double participation_threshold = 0.5,
        double vtk_scale_factor = 0.1);

    /**
     * Detect number of rigid body modes automatically
     * 
     * @param eigen_results Eigenvalue results (should be sorted)
     * @param threshold Eigenvalue threshold for rigid body modes (default: 1e-8)
     * @return Number of rigid body modes
     */
    int detectRigidBodyModes(
        const EigenResults& eigen_results,
        double threshold = 1e-8);

    /**
     * Compute density of states g(ω) using histogram binning
     * 
     * @param eigen_results Eigenvalue results
     * @param mode_indices Indices of modes to include in DOS
     * @param num_bins Number of frequency bins (default: 50)
     * @param use_gaussian_broadening Apply Gaussian broadening (default: false)
     * @param sigma Gaussian width for broadening (default: auto)
     * @return DOSResults structure
     */
    DOSResults computeDensityOfStates(
        const EigenResults& eigen_results,
        const std::vector<int>& mode_indices,
        int num_bins = 50,
        bool use_gaussian_broadening = false,
        double sigma = -1.0);

    /**
     * Export density of states to file
     * 
     * @param filename Output filename
     * @param dos_results DOS results
     * @param label Label for the data (e.g., "All modes", "Phonon modes")
     */
    void exportDensityOfStates(
        const std::string& filename,
        const DOSResults& dos_results,
        const std::string& label = "");

    /**
     * Compute and export DOS for all mode types (all, phonon, localized)
     * 
     * @param base_filename Base filename for output
     * @param eigen_results Eigenvalue results
     * @param participation Participation analysis results
     * @param num_rigid_body Number of rigid body modes to skip
     * @param num_bins Number of frequency bins (default: 50)
     * @param use_gaussian_broadening Apply Gaussian broadening (default: false)
     */
    void exportAllDensityOfStates(
        const std::string& base_filename,
        const EigenResults& eigen_results,
        const ParticipationAnalysis& participation,
        int num_rigid_body,
        int num_bins = 50,
        bool use_gaussian_broadening = false);
};

#endif // FEM_HESSIAN_ASSEMBLY_H