#ifndef STIFFNESS_SHEAR_STUDY_H
#define STIFFNESS_SHEAR_STUDY_H

#include <string>
#include <vector>

struct StiffnessShearConfig {
    int nx = 10;
    int ny = 10;
    std::string lattice_type = "square"; // "square" or "triangular"
    int potential_mode = 3;              // 1: LJ v2, 2: LJ v3, 3: Zanzotto, 4: Square Gaussian
    std::string bc = "periodic";         // "fixed" or "periodic"
    
    double alpha_min = -0.3;
    double alpha_max = 0.3;
    double dalpha = 0.005;
    
    // Zanzotto potential parameters (Mode 3)
    double beta = -0.25;
    double K = 5.0;
    
    std::string csv_output = "stiffness_shear_study.csv";
    std::string json_output = "stiffness_shear_study.json";
    std::string fullscale_json_output = "scratch/stiffness_fullscale.json";
    bool single_alpha_only = false;
    double single_alpha = 0.0;
};

struct StiffnessShearPoint {
    double alpha = 0.0;
    double lambda_min = 0.0;
    double lambda_2 = 0.0;
    double lambda_3 = 0.0;
    bool is_stable = true;
    double energy = 0.0;           // Total integrated energy over all elements
    double energy_density = 0.0;   // Energy per unit area: total_energy / total_area
    double continuum_energy = 0.0; // Homogeneous continuum strain energy density W_cont(alpha)
    // Normalized nodal eigenvectors v_1, v_2, v_3 on the full N_nodes lattice (max |du_i| = 1)
    std::vector<double> mode1_ux;
    std::vector<double> mode1_uy;
    std::vector<double> mode2_ux;
    std::vector<double> mode2_uy;
    std::vector<double> mode3_ux;
    std::vector<double> mode3_uy;
};

struct StiffnessShearSummary {
    int total_nodes = 0;
    int free_nodes = 0;
    int total_dofs = 0;
    int total_elements = 0;
    double alpha_c_pos = 0.0;
    double alpha_c_neg = 0.0;
    double min_lambda = 0.0;
    double max_lambda = 0.0;
    double execution_seconds = 0.0;
};

struct StiffnessShearResults {
    StiffnessShearConfig config;
    StiffnessShearSummary summary;
    std::vector<StiffnessShearPoint> points;
};

StiffnessShearResults runStiffnessShearStudy(const StiffnessShearConfig &config);

#endif // STIFFNESS_SHEAR_STUDY_H
