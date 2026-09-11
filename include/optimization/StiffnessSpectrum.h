// StiffnessSpectrum.h
//
// Lowest eigenpairs of the analytic FEM stiffness (Hessian) K at a configuration x,
// for stability analysis. Same K as FEMHessianAssembler::assembleGlobalStiffness, but
//
//   * assembled by FastStiffnessAssembler (parallel, no ITensor, no PSD projection),
//   * shift-invert Lanczos (Spectra) with a sparse CHOLESKY of K - sigma I instead of a
//     sparse LU: sigma is placed just below the lowest eigenvalue, so K - sigma I is
//     positive definite; a failed Cholesky means "sigma is above lambda_min" and sigma
//     is lowered (this is also a free stability test, Sylvester's law of inertia),
//   * the two rigid translations of a periodic cell (exact null vectors of K) are
//     projected out, so no eigenpairs are spent on them and no detection is needed.
#ifndef STIFFNESS_SPECTRUM_H
#define STIFFNESS_SPECTRUM_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "../include/optimization/PreconditionedLBFGS.h"

struct StiffnessSpectrumOptions {
  int n_modes = 10;       // number of lowest eigenpairs
  int ncv = 0;            // Lanczos basis size; 0 = max(2 n_modes + 1, 20)
  double tol = 1e-8;      // Spectra relative tolerance
  int max_restarts = 1000;
  // Translations: 1 = always project out, 0 = never, -1 = if ||K t|| is at rounding level
  int deflate_translations = -1;
  // First trial shift sigma = -shift_rel * mean(diag K); lowered by shift_factor until
  // K - sigma I is positive definite.
  double shift_rel = 1e-4;
  double shift_factor = 10.0;
  int max_shift_attempts = 30;
  bool compute_vectors = true;
  bool verbose = false;
};

struct StiffnessSpectrum {
  Eigen::VectorXd eigenvalues;  // ascending; translations excluded if deflated
  Eigen::MatrixXd eigenvectors; // columns, solver dof layout [u_0..u_{n-1}, v_0..v_{n-1}]
  int num_computed = 0;
  bool converged = false;
  bool translations_deflated = false;
  double sigma = 0.0;         // shift of the successful factorization (sigma < lambda_min)
  int failed_shifts = 0;      // Cholesky failures while searching sigma
  int n_operations = 0;       // shift-invert solves
  double assemble_seconds = 0.0;
  double factorize_seconds = 0.0; // all attempts (incl. symbolic analysis)
  double lanczos_seconds = 0.0;
};

// Whether both uniform translations are null vectors of K (periodic cell, all nodes free).
bool stiffness_has_translation_modes(const Eigen::SparseMatrix<double> &K);

// Lowest eigenpairs of a given symmetric K (full storage, solver dof layout).
StiffnessSpectrum lowest_stiffness_modes(const Eigen::SparseMatrix<double> &K,
                                         const StiffnessSpectrumOptions &options = {});

// Assembles the analytic (unprojected) K at x with FastStiffnessAssembler, then as above.
// Leaves the elements' cached deformation gradients at x.
StiffnessSpectrum lowest_stiffness_modes(const alglib::real_1d_array &x, UserData *userData,
                                         const StiffnessSpectrumOptions &options = {});

#endif // STIFFNESS_SPECTRUM_H
