// SparseCholesky.cpp
#include "../include/optimization/SparseCholesky.h"

#include <Eigen/SparseCholesky>
#include <algorithm>
#include <vector>

#ifdef USE_CHOLMOD
#include <cholmod.h>
#endif

namespace {

using SpMat = Eigen::SparseMatrix<double>;

// Eigen SimplicialLLT (portable, serial):  Perm M Perm^T = L L^T,  R = Perm^T L
class EigenSimplicialBackend final : public SparseCholeskyBackend {
public:
  void analyze(const SpMat &M) override { solver_.analyzePattern(M); }
  bool factorize(const SpMat &M, double shift) override {
    solver_.setShift(shift);
    solver_.factorize(M);
    return solver_.info() == Eigen::Success;
  }
  void apply_inverse(const double *g, double *gz) override {
    const int n = static_cast<int>(solver_.rows());
    work_ = solver_.permutationP() * Eigen::Map<const Eigen::VectorXd>(g, n);
    solver_.matrixL().solveInPlace(work_);
    Eigen::Map<Eigen::VectorXd>(gz, n) = work_;
  }
  void apply_inverse_transpose(const double *z, double *dx) override {
    const int n = static_cast<int>(solver_.rows());
    work_ = Eigen::Map<const Eigen::VectorXd>(z, n);
    solver_.matrixU().solveInPlace(work_);
    Eigen::Map<Eigen::VectorXd>(dx, n) = solver_.permutationPinv() * work_;
  }
  void solve(const double *b, double *x) override {
    const int n = static_cast<int>(solver_.rows());
    Eigen::Map<Eigen::VectorXd>(x, n) = solver_.solve(Eigen::Map<const Eigen::VectorXd>(b, n));
  }
  const char *name() const override { return "eigen-simplicial"; }

private:
  Eigen::SimplicialLLT<SpMat, Eigen::Lower, Eigen::AMDOrdering<int>> solver_;
  Eigen::VectorXd work_;
};

#ifdef USE_CHOLMOD
// CHOLMOD supernodal LL^T (BLAS-3, multithreaded through BLAS):
//   L L^T = Perm M Perm^T,  R = Perm^T L
class CholmodSupernodalBackend final : public SparseCholeskyBackend {
public:
  CholmodSupernodalBackend() {
    cholmod_start(&common_);
    common_.supernodal = CHOLMOD_SUPERNODAL;
    common_.final_ll = 1;
    common_.print = 0;
  }
  ~CholmodSupernodalBackend() override {
    free_work();
    if (factor_) cholmod_free_factor(&factor_, &common_);
    cholmod_finish(&common_);
  }
  void analyze(const SpMat &M) override {
    if (factor_) cholmod_free_factor(&factor_, &common_);
    cholmod_sparse A = view(M);
    factor_ = cholmod_analyze(&A, &common_);
    n_ = static_cast<int>(M.rows());
  }
  bool factorize(const SpMat &M, double shift) override {
    cholmod_sparse A = view(M);
    double beta[2] = {shift, 0.0};
    cholmod_factorize_p(&A, beta, nullptr, 0, factor_, &common_);
    return common_.status == CHOLMOD_OK && factor_->minor == factor_->n;
  }
  void apply_inverse(const double *g, double *gz) override {
    solve(CHOLMOD_P, g, buffer());
    solve(CHOLMOD_L, buffer(), gz);
  }
  void apply_inverse_transpose(const double *z, double *dx) override {
    solve(CHOLMOD_Lt, z, buffer());
    solve(CHOLMOD_Pt, buffer(), dx);
  }
  void solve(const double *b, double *x) override { solve(CHOLMOD_A, b, x); }
  const char *name() const override { return "cholmod-supernodal"; }

private:
  cholmod_sparse view(const SpMat &M) {
    cholmod_sparse A{};
    A.nrow = M.rows();
    A.ncol = M.cols();
    A.nzmax = M.nonZeros();
    A.p = const_cast<int *>(M.outerIndexPtr());
    A.i = const_cast<int *>(M.innerIndexPtr());
    A.x = const_cast<double *>(M.valuePtr());
    A.stype = -1; // symmetric, lower triangle used
    A.itype = CHOLMOD_INT;
    A.xtype = CHOLMOD_REAL;
    A.dtype = CHOLMOD_DOUBLE;
    A.sorted = 1;
    A.packed = 1;
    return A;
  }
  double *buffer() {
    tmp_.resize(n_);
    return tmp_.data();
  }
  void solve(int sys, const double *in, double *out) {
    cholmod_dense B{};
    B.nrow = n_;
    B.ncol = 1;
    B.nzmax = n_;
    B.d = n_;
    B.x = const_cast<double *>(in);
    B.xtype = CHOLMOD_REAL;
    B.dtype = CHOLMOD_DOUBLE;
    cholmod_solve2(sys, factor_, &B, nullptr, &X_, nullptr, &Y_, &E_, &common_);
    std::copy(static_cast<double *>(X_->x), static_cast<double *>(X_->x) + n_, out);
  }
  void free_work() {
    if (X_) cholmod_free_dense(&X_, &common_);
    if (Y_) cholmod_free_dense(&Y_, &common_);
    if (E_) cholmod_free_dense(&E_, &common_);
  }

  cholmod_common common_;
  cholmod_factor *factor_ = nullptr;
  cholmod_dense *X_ = nullptr, *Y_ = nullptr, *E_ = nullptr;
  std::vector<double> tmp_;
  int n_ = 0;
};
#endif

} // namespace

std::unique_ptr<SparseCholeskyBackend> make_sparse_cholesky() {
#ifdef USE_CHOLMOD
  return std::make_unique<CholmodSupernodalBackend>();
#else
  return std::make_unique<EigenSimplicialBackend>();
#endif
}
