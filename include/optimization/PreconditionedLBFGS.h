// PreconditionedLBFGS.h
//
// L-BFGS with a sparse preconditioner built from the FEM stiffness.
//
// ALGLIB's minlbfgs only accepts diagonal or dense-Cholesky preconditioners. A sparse
// preconditioner P = R R^T is applied instead through a linear change of variables
//
//     x = x0 + R^{-T} z,        dE/dz = R^{-1} dE/dx,
//
// so ALGLIB runs unmodified on z. This is exactly L-BFGS with initial inverse Hessian
// H0 = P^{-1} (plus ALGLIB's usual gamma_k scaling), at the cost of two sparse
// triangular solves per energy evaluation.
//
// Preconditioners:
//   None       - plain L-BFGS on x (reference behaviour)
//   Diagonal   - ALGLIB minlbfgssetprecdiag with diag(K)
//   Laplacian  - reference-configuration P1 Laplacian (x and y decoupled); depends on the
//                mesh only, so it is factorized once per mesh
//   Stiffness  - analytic FEM stiffness K at the starting configuration, element tangents
//                projected to PSD, plus a small shift; see refresh_every / refresh_iterations
//
// The sparse Cholesky uses CHOLMOD (supernodal) when compiled with USE_CHOLMOD, otherwise
// Eigen's SimplicialLLT.
#ifndef PRECONDITIONED_LBFGS_H
#define PRECONDITIONED_LBFGS_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "src/optimization.h" // This path is relative to ALGLIB_DIR
#include "../include/optimization/LatticeOptimizer.h"

enum class LBFGSPreconditioner { None, Diagonal, Laplacian, Stiffness };

LBFGSPreconditioner parse_lbfgs_preconditioner(const std::string &name);
std::string to_string(LBFGSPreconditioner type);

// Lagrangian tangent A_iKjL = d2W/dF_iK dF_jL of one element (per unit reference area),
// stored as a 4x4 matrix with row 2*i+K and column 2*j+L. Same physics as
// AcousticTensor::getAcousticTensor(true), without ITensor:
//   A = B^T H B + delta_ij S_KL,   B = dC_red/dF,   S = 2 Z (dE/dC_red) Z^T,
// with H the tensor second derivative of the energy in the reduced metric C_red = Z^T C Z.
Eigen::Matrix4d element_lagrangian_tangent(const Eigen::Matrix2d &F,
                                           const Eigen::Matrix2d &Z,
                                           const Eigen::Matrix2d &dE_dC_reduced,
                                           const HessianComponents &d2E_dC2_raw,
                                           double normalisation);

// 6x6 element stiffness (local dof 2*a+i) from the element tangent.
Eigen::Matrix<double, 6, 6>
element_stiffness_from_tangent(const Eigen::Matrix4d &A,
                               const Eigen::Matrix<double, 3, 2> &dN_dX,
                               double area);

// Parallel assembly of the free-dof stiffness matrix in the solver layout
// [u_0..u_{n-1}, v_0..v_{n-1}], reusing the sparsity pattern while the mesh is unchanged.
class FastStiffnessAssembler {
public:
  // Full stiffness at configuration x. If project_psd is true, negative eigenvalues of
  // each element tangent are clamped to zero (K is then positive semi-definite).
  const Eigen::SparseMatrix<double> &assembleStiffness(const alglib::real_1d_array &x,
                                                       UserData *userData,
                                                       bool project_psd);

  // Scalar reference Laplacian sum_e area_e dN_a . dN_b on the free nodes (n x n).
  const Eigen::SparseMatrix<double> &assembleLaplacian(UserData *userData);

private:
  void ensurePattern(UserData *userData);

  std::uint64_t signature_ = 0;
  int n_nodes_ = -1;

  Eigen::SparseMatrix<double> K_;
  std::vector<std::array<int, 36>> k_slots_; // per active element: index into K_ values

  Eigen::SparseMatrix<double> L_;
  std::uint64_t laplacian_signature_ = 0;
};

struct PreconditionedLBFGSOptions {
  LBFGSPreconditioner type = LBFGSPreconditioner::Stiffness;
  int corrections = 13;
  double epsg = 0.0;
  double epsf = 0.0;
  double epsx = 0.0;
  alglib::ae_int_t maxits = 0;
  // If > 0: stop when max_i |dE/dx_i| < grad_tol (measured in x for every
  // preconditioner, so runs are comparable). ALGLIB's own criteria are then disabled.
  double grad_tol = 0.0;
  bool project_psd = true;
  // Shift added to the preconditioner: P = K + shift_rel * mean(diag K) * I.
  // Removes the two rigid translations under PBC and keeps the factorization stable.
  double shift_rel = 1e-6;
  // Reuse the Stiffness factorization for at most this many optimize() calls (same mesh)...
  int refresh_every = 1;
  // ...and only while the previous call needed at most this many iterations
  // (0 = no iteration criterion). With refresh_every > 1 this refactorizes lazily,
  // e.g. right after an avalanche, and keeps the factor through quiet elastic steps.
  int refresh_iterations = 0;
  // Stiffness only: rebuild the preconditioner at the current point and continue when
  // the run stalls before grad_tol (the start-point K went stale), at most this often.
  int max_restarts = 3;
  bool verbose = true;
};

struct PreconditionedLBFGSReport {
  int iterations = 0;
  int nfev = 0;
  int termination_type = 0;
  double energy = 0.0;
  double grad_max = 0.0; // max_i |dE/dx_i| at the returned x
  double setup_seconds = 0.0;     // assembly + factorization (+ Diagonal setup)
  double assemble_seconds = 0.0;  // part of setup: stiffness / Laplacian assembly
  double factorize_seconds = 0.0; // part of setup: sparse Cholesky
  double solve_seconds = 0.0;
  bool refactorized = false;
  int restarts = 0;
};

class PreconditionedLBFGS {
public:
  explicit PreconditionedLBFGS(PreconditionedLBFGSOptions options = {});
  ~PreconditionedLBFGS();

  PreconditionedLBFGSReport optimize(alglib::real_1d_array &x, UserData *userData);

  const PreconditionedLBFGSOptions &options() const { return options_; }

  // Forget the cached factorization (call after the energy or the mesh changes by hand).
  void invalidate();

  struct Factor; // implementation detail (sparse Cholesky + change of variables)

private:
  bool buildPreconditioner(const alglib::real_1d_array &x, UserData *userData);

  PreconditionedLBFGSOptions options_;
  FastStiffnessAssembler assembler_;
  std::unique_ptr<Factor> factor_;
  int calls_since_refresh_ = 0;
  int last_iterations_ = 0;
  double last_assemble_seconds_ = 0.0;
  double last_factorize_seconds_ = 0.0;
};

#endif // PRECONDITIONED_LBFGS_H
