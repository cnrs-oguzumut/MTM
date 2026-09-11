// SparseCholesky.h
//
// Sparse LL^T of a symmetric positive definite matrix, M + shift * I = R R^T, with the
// two half-solves used by the preconditioned L-BFGS change of variables and the full
// solve used by the shift-invert eigensolver (StiffnessSpectrum).
//
// Backend: CHOLMOD supernodal when compiled with USE_CHOLMOD, otherwise Eigen's
// SimplicialLLT (AMD ordering).
#ifndef SPARSE_CHOLESKY_H
#define SPARSE_CHOLESKY_H

#include <Eigen/Sparse>
#include <memory>

class SparseCholeskyBackend {
public:
  virtual ~SparseCholeskyBackend() = default;
  // Symbolic analysis (ordering, elimination tree) of the pattern of M; reusable for
  // every M with the same pattern. M is read through its lower triangle.
  virtual void analyze(const Eigen::SparseMatrix<double> &M) = 0;
  // Numeric factorization of M + shift * I; false if it is not positive definite.
  virtual bool factorize(const Eigen::SparseMatrix<double> &M, double shift) = 0;
  virtual void apply_inverse(const double *g, double *gz) = 0;           // R^{-1} g
  virtual void apply_inverse_transpose(const double *z, double *dx) = 0; // R^{-T} z
  // x = (M + shift I)^{-1} b = R^{-T} R^{-1} b  (x and b must not alias)
  virtual void solve(const double *b, double *x) = 0;
  virtual const char *name() const = 0;
};

std::unique_ptr<SparseCholeskyBackend> make_sparse_cholesky();

#endif // SPARSE_CHOLESKY_H
