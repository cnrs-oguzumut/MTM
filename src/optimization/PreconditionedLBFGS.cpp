// PreconditionedLBFGS.cpp
#include "../include/optimization/PreconditionedLBFGS.h"

#include <Eigen/Eigenvalues>
#include <Eigen/SparseCholesky>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <omp.h>
#include <stdexcept>

#ifdef USE_CHOLMOD
#include <cholmod.h>
#endif

namespace {

using SpMat = Eigen::SparseMatrix<double>;
using Clock = std::chrono::high_resolution_clock;

double seconds_since(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}

std::uint64_t hash_mix(std::uint64_t h, std::uint64_t v) {
  // FNV-1a style mixing
  h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
  return h;
}

// Solver dof of node `a` of an element, or -1 for a fixed node.
inline int node_dof(const ElementTriangle2D &element, int a, const UserData *userData) {
  return userData->full_mapping[element.getNodeIndex(a)].second;
}

std::uint64_t mesh_signature(const UserData *userData) {
  std::uint64_t h = hash_mix(0, userData->interior_mapping.size());
  h = hash_mix(h, userData->active_elements.size());
  for (size_t idx : userData->active_elements) {
    const ElementTriangle2D &element = userData->elements[idx];
    for (int a = 0; a < 3; a++) {
      h = hash_mix(h, static_cast<std::uint64_t>(node_dof(element, a, userData) + 1));
    }
  }
  return h;
}

// Position of entry (row, col) in the value array of a compressed column-major matrix.
int value_slot(const SpMat &M, int row, int col) {
  const int *inner = M.innerIndexPtr();
  const int begin = M.outerIndexPtr()[col];
  const int end = M.outerIndexPtr()[col + 1];
  const int *it = std::lower_bound(inner + begin, inner + end, row);
  if (it == inner + end || *it != row) {
    throw std::logic_error("FastStiffnessAssembler: entry missing from sparsity pattern");
  }
  return static_cast<int>(it - inner);
}

} // namespace

LBFGSPreconditioner parse_lbfgs_preconditioner(const std::string &name) {
  if (name == "none") return LBFGSPreconditioner::None;
  if (name == "diag" || name == "diagonal") return LBFGSPreconditioner::Diagonal;
  if (name == "laplacian" || name == "lap") return LBFGSPreconditioner::Laplacian;
  if (name == "stiffness" || name == "hessian") return LBFGSPreconditioner::Stiffness;
  throw std::invalid_argument("Unknown L-BFGS preconditioner: " + name);
}

std::string to_string(LBFGSPreconditioner type) {
  switch (type) {
  case LBFGSPreconditioner::None: return "none";
  case LBFGSPreconditioner::Diagonal: return "diag";
  case LBFGSPreconditioner::Laplacian: return "laplacian";
  case LBFGSPreconditioner::Stiffness: return "stiffness";
  }
  return "unknown";
}

// ============================================================================
// Element kernels
// ============================================================================

Eigen::Matrix4d element_lagrangian_tangent(const Eigen::Matrix2d &F,
                                           const Eigen::Matrix2d &Z,
                                           const Eigen::Matrix2d &dE_dC_reduced,
                                           const HessianComponents &h,
                                           double normalisation) {
  // Tensor second derivative d2E/dC_ab dC_cd (row/col 2a+b) from the raw scalar partials,
  // with the same factors as AcousticTensor::computeHessian: 1/2 when one index pair is
  // off-diagonal, 1/4 when both are.
  const double d11 = h.c11_c11;
  const double d22 = h.c22_c22;
  const double d1122 = h.c11_c22;
  const double d1112 = 0.5 * h.c11_c12;
  const double d2212 = 0.5 * h.c22_c12;
  const double d1212 = 0.25 * h.c12_c12;
  Eigen::Matrix4d H;
  H << d11,   d1112, d1112, d1122,
       d1112, d1212, d1212, d2212,
       d1112, d1212, d1212, d2212,
       d1122, d2212, d2212, d22;
  H /= normalisation;

  // B(2a+b, 2i+K) = dC_red_ab / dF_iK with C_red = Z^T F^T F Z and G = F Z:
  //   B = Z_Ka G_ib + Z_Kb G_ia
  const Eigen::Matrix2d G = F * Z;
  Eigen::Matrix4d B;
  for (int a = 0; a < 2; a++)
    for (int b = 0; b < 2; b++)
      for (int i = 0; i < 2; i++)
        for (int K = 0; K < 2; K++)
          B(2 * a + b, 2 * i + K) = Z(K, a) * G(i, b) + Z(K, b) * G(i, a);

  Eigen::Matrix4d A = B.transpose() * H * B;

  // Geometric term delta_ij S_KL with the second Piola-Kirchhoff stress S = 2 Z dE/dC Z^T
  const Eigen::Matrix2d S = 2.0 * Z * dE_dC_reduced * Z.transpose();
  for (int i = 0; i < 2; i++)
    for (int K = 0; K < 2; K++)
      for (int L = 0; L < 2; L++)
        A(2 * i + K, 2 * i + L) += S(K, L);

  return A;
}

Eigen::Matrix<double, 6, 6>
element_stiffness_from_tangent(const Eigen::Matrix4d &A,
                               const Eigen::Matrix<double, 3, 2> &dN_dX,
                               double area) {
  // dF_iK / du_aj = delta_ij dN_a/dX_K  ->  Gm(2i+K, 2a+j)
  Eigen::Matrix<double, 4, 6> Gm = Eigen::Matrix<double, 4, 6>::Zero();
  for (int a = 0; a < 3; a++)
    for (int i = 0; i < 2; i++)
      for (int K = 0; K < 2; K++)
        Gm(2 * i + K, 2 * a + i) = dN_dX(a, K);
  return area * Gm.transpose() * A * Gm;
}

// ============================================================================
// FastStiffnessAssembler
// ============================================================================

void FastStiffnessAssembler::ensurePattern(UserData *userData) {
  const std::uint64_t sig = mesh_signature(userData);
  if (sig == signature_ && K_.rows() > 0) return;

  signature_ = sig;
  n_nodes_ = static_cast<int>(userData->interior_mapping.size());
  const int n = n_nodes_;
  const auto &active = userData->active_elements;

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(active.size() * 36);
  for (size_t idx : active) {
    const ElementTriangle2D &element = userData->elements[idx];
    int dofs[6];
    for (int a = 0; a < 3; a++) {
      const int s = node_dof(element, a, userData);
      dofs[2 * a] = s;
      dofs[2 * a + 1] = (s >= 0) ? s + n : -1;
    }
    for (int r = 0; r < 6; r++)
      for (int c = 0; c < 6; c++)
        if (dofs[r] >= 0 && dofs[c] >= 0)
          triplets.emplace_back(dofs[r], dofs[c], 1.0);
  }
  K_.resize(2 * n, 2 * n);
  K_.setFromTriplets(triplets.begin(), triplets.end());
  K_.makeCompressed();

  k_slots_.assign(active.size(), {});
  for (size_t e = 0; e < active.size(); e++) {
    const ElementTriangle2D &element = userData->elements[active[e]];
    int dofs[6];
    for (int a = 0; a < 3; a++) {
      const int s = node_dof(element, a, userData);
      dofs[2 * a] = s;
      dofs[2 * a + 1] = (s >= 0) ? s + n : -1;
    }
    for (int r = 0; r < 6; r++)
      for (int c = 0; c < 6; c++)
        k_slots_[e][6 * r + c] =
            (dofs[r] >= 0 && dofs[c] >= 0) ? value_slot(K_, dofs[r], dofs[c]) : -1;
  }
}

const Eigen::SparseMatrix<double> &
FastStiffnessAssembler::assembleStiffness(const alglib::real_1d_array &x,
                                          UserData *userData, bool project_psd) {
  ensurePattern(userData);

  const auto &active = userData->active_elements;
  const double a0 = userData->ideal_lattice_parameter;
  const double normalisation = a0 * a0;

  const std::function<double(double)> *d2pot = userData->second_derivative_function;
  if (d2pot == nullptr) {
    if (dynamic_cast<Strain_Energy_LatticeCalculator *>(&userData->calculator) == nullptr) {
      throw std::runtime_error("FastStiffnessAssembler: set UserData::second_derivative_function "
                               "for lattice-sum calculators");
    }
    d2pot = &userData->derivative_function; // ignored by the analytic strain energy
  }

  double *values = K_.valuePtr();
  std::fill(values, values + K_.nonZeros(), 0.0);

#pragma omp parallel for schedule(static)
  for (size_t e = 0; e < active.size(); e++) {
    ElementTriangle2D &element = userData->elements[active[e]];
    element.calculate_deformation_gradient(x);
    const Eigen::Matrix2d &F = element.getDeformationGradient();
    const auto red = lagrange::reduce(element.getMetricTensor());

    const Eigen::Matrix2d dE_dC =
        userData->calculator.calculate_derivative(red.C_reduced,
                                                  userData->derivative_function) /
        normalisation;
    const HessianComponents d2E = userData->calculator.calculate_dseconderivative_components(
        red.C_reduced, userData->derivative_function, *d2pot);

    Eigen::Matrix4d A =
        element_lagrangian_tangent(F, red.m_matrix, dE_dC, d2E, normalisation);
    if (project_psd) {
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> eig(A);
      const Eigen::Vector4d lambda = eig.eigenvalues().cwiseMax(0.0);
      A = eig.eigenvectors() * lambda.asDiagonal() * eig.eigenvectors().transpose();
    }
    const Eigen::Matrix<double, 6, 6> Ke =
        element_stiffness_from_tangent(A, element.getDNdX(), element.getReferenceArea());

    const auto &slots = k_slots_[e];
    for (int r = 0; r < 6; r++)
      for (int c = 0; c < 6; c++) {
        const int slot = slots[6 * r + c];
        if (slot >= 0) {
#pragma omp atomic
          values[slot] += Ke(r, c);
        }
      }
  }
  return K_;
}

const Eigen::SparseMatrix<double> &
FastStiffnessAssembler::assembleLaplacian(UserData *userData) {
  const std::uint64_t sig = mesh_signature(userData);
  if (sig == laplacian_signature_ && L_.rows() > 0) return L_;
  laplacian_signature_ = sig;

  const int n = static_cast<int>(userData->interior_mapping.size());
  const auto &active = userData->active_elements;

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(active.size() * 9);
  for (size_t idx : active) {
    const ElementTriangle2D &element = userData->elements[idx];
    const Eigen::Matrix<double, 3, 2> &dN = element.getDNdX();
    const double area = element.getReferenceArea();
    for (int a = 0; a < 3; a++) {
      const int sa = node_dof(element, a, userData);
      if (sa < 0) continue;
      for (int b = 0; b < 3; b++) {
        const int sb = node_dof(element, b, userData);
        if (sb < 0) continue;
        triplets.emplace_back(sa, sb, area * dN.row(a).dot(dN.row(b)));
      }
    }
  }
  L_.resize(n, n);
  L_.setFromTriplets(triplets.begin(), triplets.end()); // duplicates are summed
  L_.makeCompressed();
  return L_;
}

// ============================================================================
// PreconditionedLBFGS
// ============================================================================

namespace {

// Sparse LL^T of a symmetric positive definite matrix M (+ shift * I), exposing the
// two half-solves needed for the change of variables:  M = R R^T  with
//   R^{-1} g   (gradient to z space)   and   R^{-T} z   (step back to x space).
class SparseCholeskyBackend {
public:
  virtual ~SparseCholeskyBackend() = default;
  virtual void analyze(const SpMat &M) = 0;
  virtual bool factorize(const SpMat &M, double shift) = 0;
  virtual void apply_inverse(const double *g, double *gz) = 0;           // R^{-1} g
  virtual void apply_inverse_transpose(const double *z, double *dx) = 0; // R^{-T} z
  virtual const char *name() const = 0;
};

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

std::unique_ptr<SparseCholeskyBackend> make_cholesky_backend() {
#ifdef USE_CHOLMOD
  return std::make_unique<CholmodSupernodalBackend>();
#else
  return std::make_unique<EigenSimplicialBackend>();
#endif
}

} // namespace

struct PreconditionedLBFGS::Factor {
  LBFGSPreconditioner type = LBFGSPreconditioner::None;
  std::unique_ptr<SparseCholeskyBackend> cholesky;
  bool analyzed = false;
  std::uint64_t signature = 0;
  bool scalar_blocks = false; // Laplacian: the same n x n factor acts on u and on v
  alglib::real_1d_array diagonal; // Diagonal preconditioner (x space)

  // dx = R^{-T} z  with  P = R R^T
  void apply_inverse_transpose(const double *z, double *dx, int n_total) {
    const int m = scalar_blocks ? n_total / 2 : n_total;
    for (int block = 0; block < (scalar_blocks ? 2 : 1); block++)
      cholesky->apply_inverse_transpose(z + block * m, dx + block * m);
  }

  // gz = R^{-1} g
  void apply_inverse(const double *g, double *gz, int n_total) {
    const int m = scalar_blocks ? n_total / 2 : n_total;
    for (int block = 0; block < (scalar_blocks ? 2 : 1); block++)
      cholesky->apply_inverse(g + block * m, gz + block * m);
  }

  bool uses_transform() const {
    return type == LBFGSPreconditioner::Laplacian || type == LBFGSPreconditioner::Stiffness;
  }
};

namespace {

struct TransformContext {
  UserData *userData = nullptr;
  PreconditionedLBFGS::Factor *factor = nullptr; // null when no change of variables
  std::vector<double> x0;
  alglib::real_1d_array x_work;
  alglib::real_1d_array gx_work;
  double grad_tol = 0.0;
  alglib::minlbfgsstate *state = nullptr;
  bool hit = false;
  std::vector<double> x_hit;
};

} // namespace

PreconditionedLBFGS::PreconditionedLBFGS(PreconditionedLBFGSOptions options)
    : options_(options) {}

PreconditionedLBFGS::~PreconditionedLBFGS() = default;

void PreconditionedLBFGS::invalidate() {
  factor_.reset();
  calls_since_refresh_ = 0;
}

bool PreconditionedLBFGS::buildPreconditioner(const alglib::real_1d_array &x,
                                              UserData *userData) {
  const LBFGSPreconditioner type = options_.type;
  if (type == LBFGSPreconditioner::None) return false;

  const std::uint64_t sig = mesh_signature(userData);
  const bool same_mesh = factor_ && factor_->type == type && factor_->signature == sig;

  if (type == LBFGSPreconditioner::Laplacian && same_mesh) return false;
  if (type == LBFGSPreconditioner::Stiffness && same_mesh &&
      calls_since_refresh_ < options_.refresh_every &&
      (options_.refresh_iterations <= 0 || last_iterations_ <= options_.refresh_iterations))
    return false;

  if (!same_mesh) {
    factor_ = std::make_unique<Factor>();
    factor_->type = type;
    factor_->signature = sig;
  }
  calls_since_refresh_ = 0;

  if (type == LBFGSPreconditioner::Diagonal) {
    const SpMat &K = assembler_.assembleStiffness(x, userData, options_.project_psd);
    const Eigen::VectorXd d = K.diagonal();
    const double floor = std::max(1e-12, options_.shift_rel * d.mean());
    factor_->diagonal.setlength(d.size());
    for (int i = 0; i < d.size(); i++) factor_->diagonal[i] = std::max(d(i), floor);
    return true;
  }

  const auto assemble_start = Clock::now();
  const SpMat &M = (type == LBFGSPreconditioner::Laplacian)
                       ? assembler_.assembleLaplacian(userData)
                       : assembler_.assembleStiffness(x, userData, options_.project_psd);
  last_assemble_seconds_ = seconds_since(assemble_start);
  factor_->scalar_blocks = (type == LBFGSPreconditioner::Laplacian);

  const auto factorize_start = Clock::now();
  if (!factor_->analyzed) {
    factor_->cholesky = make_cholesky_backend();
    factor_->cholesky->analyze(M);
    factor_->analyzed = true;
  }

  double shift = options_.shift_rel * M.diagonal().mean();
  for (int attempt = 0; attempt < 6; attempt++) {
    if (factor_->cholesky->factorize(M, shift)) {
      last_factorize_seconds_ = seconds_since(factorize_start);
      return true;
    }
    std::cerr << "PreconditionedLBFGS: Cholesky failed with shift " << shift
              << ", increasing shift" << std::endl;
    shift *= 100.0;
  }
  std::cerr << "PreconditionedLBFGS: factorization failed, running without preconditioner"
            << std::endl;
  factor_.reset();
  return false;
}

static void transformed_energy(const alglib::real_1d_array &z, double &func,
                               alglib::real_1d_array &grad, void *ptr) {
  auto *ctx = static_cast<TransformContext *>(ptr);
  PreconditionedLBFGS::Factor *factor = ctx->factor;
  const int n = static_cast<int>(z.length());
  double *x = ctx->x_work.getcontent();

  if (factor) {
    factor->apply_inverse_transpose(z.getcontent(), x, n);
    for (int i = 0; i < n; i++) x[i] += ctx->x0[i];
  } else {
    std::copy(z.getcontent(), z.getcontent() + n, x);
  }

  minimize_energy_with_triangles(ctx->x_work, func, ctx->gx_work, ctx->userData);

  const double *gx = ctx->gx_work.getcontent();
  if (factor) {
    factor->apply_inverse(gx, grad.getcontent(), n);
  } else {
    std::copy(gx, gx + n, grad.getcontent());
  }

  if (ctx->grad_tol > 0.0 && !ctx->hit) {
    double gmax = 0.0;
    for (int i = 0; i < n; i++) gmax = std::max(gmax, std::abs(gx[i]));
    if (gmax < ctx->grad_tol) {
      ctx->hit = true;
      ctx->x_hit.assign(x, x + n);
      alglib::minlbfgsrequesttermination(*ctx->state);
    }
  }
}

PreconditionedLBFGSReport PreconditionedLBFGS::optimize(alglib::real_1d_array &x,
                                                        UserData *userData) {
  PreconditionedLBFGSReport report;
  const int n = static_cast<int>(x.length());

  const auto setup_start = Clock::now();
  last_assemble_seconds_ = last_factorize_seconds_ = 0.0;
  report.refactorized = buildPreconditioner(x, userData);
  calls_since_refresh_++;
  report.setup_seconds = seconds_since(setup_start);
  report.assemble_seconds = last_assemble_seconds_;
  report.factorize_seconds = last_factorize_seconds_;

  const auto solve_start = Clock::now();
  double rebuild_seconds = 0.0;

  TransformContext ctx;
  ctx.userData = userData;
  ctx.x_work.setlength(n);
  ctx.gx_work.setlength(n);
  ctx.grad_tol = options_.grad_tol;

  alglib::real_1d_array z;
  z.setlength(n);
  alglib::minlbfgsstate state;
  alglib::minlbfgsreport rep;
  ctx.state = &state;

  for (int run = 0;; run++) {
    const bool transform = factor_ && factor_->uses_transform();
    ctx.factor = transform ? factor_.get() : nullptr;
    ctx.x0.assign(x.getcontent(), x.getcontent() + n);
    for (int i = 0; i < n; i++) z[i] = transform ? 0.0 : x[i];

    alglib::minlbfgscreate(options_.corrections, z, state);
    if (options_.grad_tol > 0.0) {
      // Only the x-space gradient test (and a safety cap) terminates the run.
      alglib::minlbfgssetcond(state, 0.0, 0.0, 0.0,
                              options_.maxits > 0 ? options_.maxits : 200000);
    } else {
      alglib::minlbfgssetcond(state, options_.epsg, options_.epsf, options_.epsx,
                              options_.maxits);
    }
    if (factor_ && options_.type == LBFGSPreconditioner::Diagonal) {
      alglib::minlbfgssetprecdiag(state, factor_->diagonal);
    }

    alglib::minlbfgsoptimize(state, transformed_energy, nullptr, &ctx);
    alglib::minlbfgsresults(state, z, rep);

    if (ctx.hit) {
      for (int i = 0; i < n; i++) x[i] = ctx.x_hit[i];
    } else if (transform) {
      factor_->apply_inverse_transpose(z.getcontent(), x.getcontent(), n);
      for (int i = 0; i < n; i++) x[i] += ctx.x0[i];
    } else {
      for (int i = 0; i < n; i++) x[i] = z[i];
    }
    report.iterations += static_cast<int>(rep.iterationscount);
    report.nfev += static_cast<int>(rep.nfev);
    report.termination_type = static_cast<int>(rep.terminationtype);

    // A Stiffness preconditioner built at the start point can go stale during a large
    // rearrangement (avalanche); ALGLIB then stalls before reaching grad_tol. Rebuild it
    // at the current point and continue.
    const bool stalled = options_.grad_tol > 0.0 && !ctx.hit && rep.terminationtype != 5;
    if (!stalled || run >= options_.max_restarts ||
        options_.type != LBFGSPreconditioner::Stiffness || !factor_)
      break;
    const auto rebuild_start = Clock::now();
    calls_since_refresh_ = options_.refresh_every; // force a refactorization
    buildPreconditioner(x, userData);
    report.assemble_seconds += last_assemble_seconds_;
    report.factorize_seconds += last_factorize_seconds_;
    calls_since_refresh_ = 1;
    const double rebuild = seconds_since(rebuild_start);
    report.setup_seconds += rebuild;
    rebuild_seconds += rebuild;
    report.restarts++;
  }

  // Final energy and x-space gradient at the returned point (also leaves the
  // elements' cached F consistent with x).
  minimize_energy_with_triangles(x, report.energy, ctx.gx_work, userData);
  for (int i = 0; i < n; i++)
    report.grad_max = std::max(report.grad_max, std::abs(ctx.gx_work[i]));

  report.solve_seconds = seconds_since(solve_start) - rebuild_seconds;
  last_iterations_ = report.iterations;

  if (options_.verbose) {
    std::cout << "Optimization completed [" << to_string(options_.type)
              << "]: iterations=" << report.iterations << ", nfev=" << report.nfev
              << ", termination type=" << report.termination_type
              << ", max|g|=" << report.grad_max << ", setup=" << report.setup_seconds
              << "s, solve=" << report.solve_seconds << "s"
              << (report.restarts ? ", restarts=" + std::to_string(report.restarts) : "")
              << std::endl;
  }
  return report;
}
