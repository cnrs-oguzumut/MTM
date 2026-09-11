// StiffnessSpectrum.cpp
#include "../include/optimization/StiffnessSpectrum.h"
#include "../include/optimization/SparseCholesky.h"

#include <Spectra/SymEigsShiftSolver.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>

namespace {

using SpMat = Eigen::SparseMatrix<double>;
using Clock = std::chrono::high_resolution_clock;

double seconds_since(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}

// Removes the means of the u block and of the v block, i.e. applies
// Q = I - t_u t_u^T - t_v t_v^T with the unit translation vectors t_u, t_v.
void project_out_translations(double *v, int n_total) {
  const int m = n_total / 2;
  for (int block = 0; block < 2; block++) {
    double *b = v + block * m;
    const double mean = std::accumulate(b, b + m, 0.0) / m;
    for (int i = 0; i < m; i++) b[i] -= mean;
  }
}

// Spectra operator y = Q (K - sigma I)^{-1} Q x, using a Cholesky factorization of
// K - sigma I computed beforehand (set_shift only checks that it is the same sigma).
class CholeskyShiftInvertOp {
public:
  using Scalar = double;

  CholeskyShiftInvertOp(SparseCholeskyBackend &cholesky, int n, double sigma, bool deflate)
      : cholesky_(cholesky), n_(n), sigma_(sigma), deflate_(deflate), work_(n) {}

  Eigen::Index rows() const { return n_; }
  Eigen::Index cols() const { return n_; }

  void set_shift(const Scalar &sigma) {
    if (sigma != sigma_) {
      throw std::logic_error("CholeskyShiftInvertOp: factorization was done for another shift");
    }
  }

  void perform_op(const Scalar *x_in, Scalar *y_out) const {
    std::copy(x_in, x_in + n_, work_.data());
    if (deflate_) project_out_translations(work_.data(), n_);
    cholesky_.solve(work_.data(), y_out);
    if (deflate_) project_out_translations(y_out, n_);
    ++count_;
  }

  int count() const { return count_; }

private:
  SparseCholeskyBackend &cholesky_;
  int n_;
  double sigma_;
  bool deflate_;
  mutable Eigen::VectorXd work_;
  mutable int count_ = 0;
};

} // namespace

bool stiffness_has_translation_modes(const SpMat &K) {
  const int n = static_cast<int>(K.rows());
  if (n % 2 != 0) return false;
  const int m = n / 2;
  const double scale = K.coeffs().cwiseAbs().maxCoeff();
  for (int block = 0; block < 2; block++) {
    Eigen::VectorXd t = Eigen::VectorXd::Zero(n);
    t.segment(block * m, m).setOnes();
    if ((K * t).cwiseAbs().maxCoeff() > 1e-9 * scale) return false;
  }
  return true;
}

StiffnessSpectrum lowest_stiffness_modes(const SpMat &K, const StiffnessSpectrumOptions &opt) {
  StiffnessSpectrum out;
  const int n = static_cast<int>(K.rows());

  const bool deflate = opt.deflate_translations == 1 ||
                       (opt.deflate_translations == -1 && stiffness_has_translation_modes(K));
  out.translations_deflated = deflate;
  const int n_available = deflate ? n - 2 : n;
  const int nev = std::min(opt.n_modes, n_available - 1);
  const int ncv = std::min(opt.ncv > 0 ? opt.ncv : std::max(2 * nev + 1, 20), n_available);
  if (nev <= 0) return out;

  // Shift search: K - sigma I is positive definite  <=>  sigma < lambda_min(K).
  const auto factorize_start = Clock::now();
  std::unique_ptr<SparseCholeskyBackend> cholesky = make_sparse_cholesky();
  cholesky->analyze(K);
  const double scale = K.diagonal().cwiseAbs().mean();
  double sigma = -opt.shift_rel * scale;
  bool factorized = false;
  for (int attempt = 0; attempt < opt.max_shift_attempts; attempt++) {
    if (cholesky->factorize(K, -sigma)) {
      factorized = true;
      break;
    }
    out.failed_shifts++;
    sigma *= opt.shift_factor;
  }
  out.factorize_seconds = seconds_since(factorize_start);
  out.sigma = sigma;
  if (!factorized) {
    std::cerr << "lowest_stiffness_modes: no positive definite shift found down to sigma = "
              << sigma << std::endl;
    return out;
  }

  const auto lanczos_start = Clock::now();
  CholeskyShiftInvertOp op(*cholesky, n, sigma, deflate);
  Spectra::SymEigsShiftSolver<CholeskyShiftInvertOp> eigs(op, nev, ncv, sigma);
  eigs.init();
  const int nconv = static_cast<int>(
      eigs.compute(Spectra::SortRule::LargestMagn, opt.max_restarts, opt.tol,
                   Spectra::SortRule::SmallestAlge));
  out.lanczos_seconds = seconds_since(lanczos_start);
  out.n_operations = op.count();
  out.converged = eigs.info() == Spectra::CompInfo::Successful;
  out.num_computed = nconv;
  if (nconv > 0) {
    out.eigenvalues = eigs.eigenvalues();
    if (opt.compute_vectors) out.eigenvectors = eigs.eigenvectors();
    // SmallestAlge sorting already gives ascending order; enforce it anyway.
    std::vector<int> order(nconv);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](int a, int b) { return out.eigenvalues(a) < out.eigenvalues(b); });
    Eigen::VectorXd values(nconv);
    Eigen::MatrixXd vectors(opt.compute_vectors ? n : 0, opt.compute_vectors ? nconv : 0);
    for (int k = 0; k < nconv; k++) {
      values(k) = out.eigenvalues(order[k]);
      if (opt.compute_vectors) vectors.col(k) = out.eigenvectors.col(order[k]);
    }
    out.eigenvalues = values;
    out.eigenvectors = vectors;
  }

  if (opt.verbose) {
    std::cout << "lowest_stiffness_modes [" << cholesky->name() << "]: n=" << n
              << ", modes=" << nconv << "/" << nev << (deflate ? " (translations deflated)" : "")
              << ", sigma=" << sigma << " (" << out.failed_shifts << " failed shifts)"
              << ", lambda_min=" << (nconv ? out.eigenvalues(0) : NAN)
              << ", solves=" << out.n_operations << ", factorize=" << out.factorize_seconds
              << "s, lanczos=" << out.lanczos_seconds << "s" << std::endl;
  }
  return out;
}

StiffnessSpectrum lowest_stiffness_modes(const alglib::real_1d_array &x, UserData *userData,
                                         const StiffnessSpectrumOptions &opt) {
  const auto assemble_start = Clock::now();
  FastStiffnessAssembler assembler;
  const SpMat &K = assembler.assembleStiffness(x, userData, /*project_psd=*/false);
  const double assemble_seconds = seconds_since(assemble_start);
  StiffnessSpectrum out = lowest_stiffness_modes(K, opt);
  out.assemble_seconds = assemble_seconds;
  return out;
}
