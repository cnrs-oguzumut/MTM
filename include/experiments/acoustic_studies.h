#pragma once
#include <string>

struct AcousticStudyConfig {
    int mode = 3; // 1: LJ v2, 2: LJ v3, 3: Zanzotto, 4: Square Gaussian
    std::string lattice_type = "square"; // "square" or "triangular"
    double beta = -0.25;
    double K = 5.0;
    double t_max = 1.0;
    double dt = 0.02;
    double dp = 0.05; // radians (~2*pi / 125)
    std::string dat_output = "acoustic_study_results.dat";
    std::string json_output = "acoustic_study_results.json";
    std::string csv_output = "energy_path.csv";
    std::string path_type = "simple_shear"; // "simple_shear", "pure_shear", "diagonal_shear"
    double shear_theta_deg = 0.0; // Angle of slip direction theta in degrees
    double alpha_min = -1.5;
    double alpha_max = 1.5;
    double dalpha = 0.005;
    bool only_elastic_constants = false;
};

struct ShearPathPoint {
    double alpha = 0.0;
    double c11 = 1.0, c22 = 1.0, c12 = 0.0;
    double c11_red = 1.0, c22_red = 1.0, c12_red = 0.0;
    double x = 0.0, y = 0.0;
    double x_red = 0.0, y_red = 0.0;
    double min_detAc = 0.0;
    double energy = 0.0;
    bool is_stable = true;
    double angle_deg = 0.0;        // Lagrangian critical normal angle N* (deg)
    double angle_euler_deg = 0.0;  // Eulerian critical normal angle n* (deg)
    double lambda_min_Q = 0.0;     // Smallest eigenvalue of Q(n*)
    double lambda_max_Q = 0.0;     // Largest eigenvalue of Q(n*)
    double mode_mx = 1.0;          // Spatial polarization eigenvector m*_x
    double mode_my = 0.0;          // Spatial polarization eigenvector m*_y
    double f11 = 1.0, f12 = 0.0, f21 = 0.0, f22 = 1.0;
    double q11_xx = 0.0, q11_yy = 0.0, q11_xy = 0.0;
    double q22_xx = 0.0, q22_yy = 0.0, q22_xy = 0.0;
    double q12_xx = 0.0, q12_yy = 0.0, q12_xy = 0.0;
};

struct ShearPathResult {
    std::vector<ShearPathPoint> path;
    std::vector<double> critical_alphas;
    std::vector<ShearPathPoint> critical_points;
    std::string path_type = "simple_shear";
    double theta_deg = 0.0;
    double period_alpha = 1.0; // gamma^2 for triangular, 1.0 for square
};

struct AcousticElasticResults {
    double C11 = 0.0;
    double C22 = 0.0;
    double C12 = 0.0;
    double C66 = 0.0;
    double Zener_A = 1.0;
    double bulk_K = 0.0;
    double shear_mu = 0.0;
    bool born_stable = true;
    bool is_isotropic = true;
    std::string anisotropy_comment;
};

AcousticElasticResults computeAcousticElasticConstants(const AcousticStudyConfig& config);
void printAcousticElasticConstantsJSON(const AcousticStudyConfig& config);
ShearPathResult computeSimpleShearPath(const AcousticStudyConfig& config, double alpha_min, double alpha_max, double dalpha);
void runParametricAcousticStudy(const AcousticStudyConfig& config);

void parametricAcousticStudy();
void parametricAcousticStudy_v2();


