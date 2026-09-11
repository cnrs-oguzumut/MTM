#pragma once

#include <string>

// How analyze_data_from_folder computes the lowest eigenvalues of the stiffness:
//   Fast   - analytic K (FastStiffnessAssembler) + Cholesky shift-invert Lanczos,
//            lowest_stiffness_modes (StiffnessSpectrum.h)
//   Legacy - FEMHessianAssembler: ITensor K + Spectra shift-invert with SparseLU
enum class StiffnessEigenSolver { Fast, Legacy };

StiffnessEigenSolver parse_stiffness_eigen_solver(const std::string &name);
std::string to_string(StiffnessEigenSolver solver);

void analyze_data_from_folder(int caller_id, int nx, int ny, int iter_start,
                              int iter_end, int n_eig,
                              StiffnessEigenSolver eig_solver = StiffnessEigenSolver::Fast);
