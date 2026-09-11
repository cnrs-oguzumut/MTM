// benchmark_preconditioner.cpp
//
// Compares L-BFGS preconditioners on the Conti-Zanzotto simple-shear loading
// (square lattice, PBC, fixed mesh). At every load step all methods start from the
// same configuration; the state is then advanced with the first method's solution,
// so the first method (default "prod") reproduces the production trajectory.
//
// Methods:
//   prod       plain L-BFGS with ALGLIB's automatic stopping (what production runs use)
//   none       plain L-BFGS, stopped at max|dE/dx| < grad_tol
//   diag       ALGLIB diagonal preconditioner (diag K), same stopping
//   laplacian  sparse reference Laplacian via change of variables, same stopping
//   stiffness  sparse analytic stiffness K via change of variables, same stopping
//   <m>@N      reuse the Stiffness factorization for N load steps (e.g. stiffness@10)
//
// Usage:
//   benchmark_preconditioner [nx ny] [--steps N] [--methods prod,none,diag,laplacian,stiffness]
//                            [--grad-tol 1e-6] [--alpha0 0.14] [--dalpha 6e-5] [--seed 42]
//                            [--corrections 13] [--refresh N] [--refresh-its N] [--validate]
//                            [--csv file] [--eig N]
//
//   --eig N   after the last load step, compute the N lowest eigenvalues of the stiffness
//             with FEMHessianAssembler (ITensor K + Spectra/SparseLU, as in data_analysis)
//             and with lowest_stiffness_modes (fast K + Cholesky shift-invert), and compare.

#include "../include/experiments/experiment_includes.h"
#include "../include/optimization/PreconditionedLBFGS.h"
#include "../include/optimization/StiffnessSpectrum.h"

#include <chrono>
#include <sstream>

namespace {

struct Args {
  int nx = 50, ny = 50;
  int steps = 20;
  std::vector<std::string> methods = {"prod", "none", "diag", "laplacian", "stiffness"};
  double grad_tol = 1e-6;
  double alpha0 = 0.14;
  double dalpha = 6e-5;
  unsigned int seed = 42;
  int corrections = 13;
  int refresh = 1;
  int refresh_its = 0;
  bool validate = false;
  int eig = 0;
  std::string csv = "precond_benchmark.csv";
};

Args parse_args(int argc, char **argv) {
  Args args;
  std::vector<std::string> positional;
  for (int i = 1; i < argc; i++) {
    std::string a = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc) throw std::invalid_argument("missing value for " + a);
      return argv[++i];
    };
    if (a == "--steps") args.steps = std::stoi(next());
    else if (a == "--grad-tol") args.grad_tol = std::stod(next());
    else if (a == "--alpha0") args.alpha0 = std::stod(next());
    else if (a == "--dalpha") args.dalpha = std::stod(next());
    else if (a == "--seed") args.seed = static_cast<unsigned int>(std::stoul(next()));
    else if (a == "--corrections") args.corrections = std::stoi(next());
    else if (a == "--refresh") args.refresh = std::stoi(next());
    else if (a == "--refresh-its") args.refresh_its = std::stoi(next());
    else if (a == "--csv") args.csv = next();
    else if (a == "--validate") args.validate = true;
    else if (a == "--eig") args.eig = std::stoi(next());
    else if (a == "--methods") {
      args.methods.clear();
      std::stringstream ss(next());
      std::string m;
      while (std::getline(ss, m, ',')) args.methods.push_back(m);
    } else positional.push_back(a);
  }
  if (positional.size() >= 2) {
    args.nx = std::stoi(positional[0]);
    args.ny = std::stoi(positional[1]);
  }
  return args;
}

// Verifies the fast analytic stiffness against (1) the ITensor-based
// FEMHessianAssembler::computeElementStiffness element by element and (2) central
// differences of the analytic gradient along random directions.
void validate_stiffness(const alglib::real_1d_array &x, UserData &userData,
                        BaseLatticeCalculator &calculator,
                        const std::function<double(double)> &dpot,
                        const std::function<double(double)> &d2pot) {
  std::cout << "\n=== VALIDATION: fast analytic stiffness ===" << std::endl;

  // (1) element-by-element against the ITensor assembler
  FEMHessianAssembler reference;
  reference.setEnergyParameters(&calculator, dpot, d2pot, calculator.getUnitCellArea());
  const double a0 = userData.ideal_lattice_parameter;
  double max_diff = 0.0, max_ref = 0.0;
  size_t n_checked = 0, n_reduced = 0;
  const auto t_ref = std::chrono::high_resolution_clock::now();
  for (size_t idx : userData.active_elements) {
    ElementTriangle2D &element = userData.elements[idx];
    element.calculate_deformation_gradient(x);
    const Eigen::Matrix2d F = element.getDeformationGradient();
    const auto red = lagrange::reduce(F.transpose() * F);
    if (!red.m_matrix.isIdentity(1e-12)) n_reduced++;
    AcousticTensor acoustic(F, red.C_reduced, red.m_matrix);
    const Eigen::MatrixXd K_ref =
        reference.computeElementStiffness(acoustic, element, element.getReferenceArea());
    n_checked++;
    const Eigen::Matrix2d dE_dC =
        calculator.calculate_derivative(red.C_reduced, dpot) / (a0 * a0);
    const HessianComponents d2 =
        calculator.calculate_dseconderivative_components(red.C_reduced, dpot, d2pot);
    const Eigen::Matrix4d A =
        element_lagrangian_tangent(F, red.m_matrix, dE_dC, d2, a0 * a0);
    const Eigen::Matrix<double, 6, 6> K_fast =
        element_stiffness_from_tangent(A, element.getDNdX(), element.getReferenceArea());
    max_diff = std::max(max_diff, (K_fast - K_ref).cwiseAbs().maxCoeff());
    max_ref = std::max(max_ref, K_ref.cwiseAbs().maxCoeff());
  }
  const double t_itensor =
      std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t_ref).count();
  std::cout << "Elements checked: " << n_checked << " (" << n_reduced
            << " with non-identity Lagrange reduction)" << std::endl;
  std::cout << "max |K_fast - K_itensor| = " << max_diff << "  (max |K_itensor| = " << max_ref
            << ", relative " << max_diff / max_ref << ")" << std::endl;
  std::cout << "Serial ITensor element loop (reference + fast): " << t_itensor << " s"
            << std::endl;

  // (2) global K v against central differences of the analytic gradient
  FastStiffnessAssembler assembler;
  const auto t_fast = std::chrono::high_resolution_clock::now();
  const Eigen::SparseMatrix<double> K = assembler.assembleStiffness(x, &userData, false);
  const double t_asm =
      std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t_fast).count();
  const auto t_fast2 = std::chrono::high_resolution_clock::now();
  assembler.assembleStiffness(x, &userData, false);
  const double t_asm2 =
      std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t_fast2).count();
  std::cout << "Fast parallel assembly: " << t_asm << " s (first, builds pattern), " << t_asm2
            << " s (pattern reused); nnz = " << K.nonZeros() << std::endl;
  std::cout << "Symmetry: max |K - K^T| = "
            << Eigen::SparseMatrix<double>(K - Eigen::SparseMatrix<double>(K.transpose()))
                   .coeffs()
                   .cwiseAbs()
                   .maxCoeff()
            << std::endl;

  const int n = static_cast<int>(x.length());
  std::mt19937 gen(7);
  std::normal_distribution<double> nd(0.0, 1.0);
  alglib::real_1d_array xp, xm, gp, gm;
  xp.setlength(n); xm.setlength(n); gp.setlength(n); gm.setlength(n);
  for (int trial = 0; trial < 3; trial++) {
    Eigen::VectorXd v(n);
    for (int i = 0; i < n; i++) v(i) = nd(gen);
    v.normalize();
    const double eps = 1e-6;
    for (int i = 0; i < n; i++) {
      xp[i] = x[i] + eps * v(i);
      xm[i] = x[i] - eps * v(i);
    }
    double fp, fm;
    minimize_energy_with_triangles(xp, fp, gp, &userData);
    minimize_energy_with_triangles(xm, fm, gm, &userData);
    Eigen::VectorXd fd(n);
    for (int i = 0; i < n; i++) fd(i) = (gp[i] - gm[i]) / (2.0 * eps);
    const Eigen::VectorXd Kv = K * v;

    // Energy/gradient consistency at a non-equilibrium point (at the relaxed x the
    // directional derivative is at rounding level): xb = x + 1e-2 w, w random.
    alglib::real_1d_array xb, gb;
    xb.setlength(n);
    gb.setlength(n);
    for (int i = 0; i < n; i++) {
      xb[i] = x[i] + 1e-2 * nd(gen);
    }
    double fb;
    minimize_energy_with_triangles(xb, fb, gb, &userData);
    for (int i = 0; i < n; i++) {
      xp[i] = xb[i] + eps * v(i);
      xm[i] = xb[i] - eps * v(i);
    }
    minimize_energy_with_triangles(xp, fp, gp, &userData);
    minimize_energy_with_triangles(xm, fm, gm, &userData);
    double gv = 0.0;
    for (int i = 0; i < n; i++) gv += gb[i] * v(i);
    const double dfdv = (fp - fm) / (2.0 * eps);
    std::cout << "Directional check " << trial << ": |Kv - dg/dv| / |Kv| = "
              << (Kv - fd).norm() / Kv.norm() << ",  (dg/dv . Kv)/|Kv|^2 = "
              << fd.dot(Kv) / Kv.squaredNorm() << ",  off-equilibrium dE/dv / (g.v) = " << dfdv / gv
              << std::endl;
  }
  // Restore element state at x
  double f;
  alglib::real_1d_array g;
  g.setlength(n);
  minimize_energy_with_triangles(x, f, g, &userData);
  std::cout << "===========================================\n" << std::endl;
}

// Lowest eigenvalues of K at the relaxed x: the data_analysis path (ITensor assembly +
// Spectra shift-invert with SparseLU, translations included) against lowest_stiffness_modes.
void compare_eigen_solvers(const alglib::real_1d_array &x, UserData &userData,
                           BaseLatticeCalculator &calculator,
                           const std::function<double(double)> &dpot,
                           const std::function<double(double)> &d2pot, int n_eig) {
  std::cout << "\n=== LOWEST EIGENVALUES OF K: " << n_eig << " modes ===" << std::endl;
  const int n = static_cast<int>(x.length());
  using clk = std::chrono::high_resolution_clock;
  auto secs = [](clk::time_point t0) {
    return std::chrono::duration<double>(clk::now() - t0).count();
  };

  // (1) data_analysis path
  std::vector<Point2D> points = userData.points;
  map_solver_array_to_points(x, points, userData.interior_mapping, n / 2);
  FEMHessianAssembler reference;
  reference.setEnergyParameters(&calculator, dpot, d2pot, calculator.getUnitCellArea());
  auto t0 = clk::now();
  const Eigen::SparseMatrix<double> K_ref =
      reference.assembleGlobalStiffness(userData.elements, points, n, userData.full_mapping);
  const double t_ref_asm = secs(t0);
  t0 = clk::now();
  const EigenResults ref = reference.computeSmallestEigenvaluesIterative_spectra(K_ref, n_eig, -1);
  const double t_ref_eig = secs(t0);

  // (2) fast path; the fast path skips the two translations, so ask for n_eig - 2
  StiffnessSpectrumOptions opt;
  opt.n_modes = n_eig - 2;
  opt.verbose = true;
  t0 = clk::now();
  const StiffnessSpectrum fast = lowest_stiffness_modes(x, &userData, opt);
  const double t_fast = secs(t0);

  // Same K?
  FastStiffnessAssembler assembler;
  const Eigen::SparseMatrix<double> &K_fast = assembler.assembleStiffness(x, &userData, false);
  const double k_diff = Eigen::SparseMatrix<double>(K_fast - K_ref).coeffs().cwiseAbs().maxCoeff();

  std::cout << "max |K_fast - K_itensor| = " << k_diff << std::endl;
  std::cout << "old: assembly " << t_ref_asm << " s + eigensolver " << t_ref_eig
            << " s = " << t_ref_asm + t_ref_eig << " s  (" << ref.num_computed
            << " values incl. translations)" << std::endl;
  std::cout << "new: total " << t_fast << " s = assembly " << fast.assemble_seconds
            << " + factorization(s) " << fast.factorize_seconds << " + Lanczos "
            << fast.lanczos_seconds << "  (" << fast.num_computed << " values, "
            << fast.n_operations << " solves, converged=" << fast.converged << ")" << std::endl;
  std::cout << "speedup: " << (t_ref_asm + t_ref_eig) / t_fast << "x" << std::endl;

  // Old list without the translations (|lambda| < 1e-8, as in detectRigidBodyModes).
  // Lanczos may return only one of the two (degenerate) translations when few values are
  // requested.
  std::vector<double> old_values;
  std::cout << "old translation eigenvalues: ";
  for (int i = 0; i < ref.num_computed; i++) {
    if (std::abs(ref.eigenvalues(i)) < 1e-8) std::cout << ref.eigenvalues(i) << " ";
    else old_values.push_back(ref.eigenvalues(i));
  }
  std::cout << std::endl;
  std::sort(old_values.begin(), old_values.end());

  double max_rel = 0.0;
  const int m = std::min<int>(old_values.size(), fast.num_computed);
  for (int k = 0; k < m; k++) {
    const double rel = std::abs(old_values[k] - fast.eigenvalues(k)) /
                       std::max(std::abs(old_values[k]), 1e-300);
    max_rel = std::max(max_rel, rel);
    if (k < 6 || k == m - 1)
      std::cout << "  mode " << std::setw(3) << k << ": old " << std::setw(14) << old_values[k]
                << "  new " << std::setw(14) << fast.eigenvalues(k) << std::endl;
  }
  std::cout << "max relative difference over " << m << " modes: " << max_rel << std::endl;

  // Residuals of the new eigenpairs on the unprojected K
  double max_res = 0.0;
  for (int k = 0; k < fast.num_computed; k++) {
    const Eigen::VectorXd v = fast.eigenvectors.col(k);
    max_res = std::max(max_res, (K_fast * v - fast.eigenvalues(k) * v).norm() /
                                    std::max(1e-300, std::abs(fast.eigenvalues(k))));
  }
  std::cout << "max |K v - lambda v| / |lambda| (new): " << max_res << std::endl;

  // Smaller requests, as used for a stability check
  for (int k : {1, 10}) {
    StiffnessSpectrumOptions o;
    o.n_modes = k;
    o.compute_vectors = true;
    t0 = clk::now();
    const StiffnessSpectrum s = lowest_stiffness_modes(x, &userData, o);
    std::cout << "new, " << k << " mode(s): " << secs(t0) << " s, lambda_min = "
              << (s.num_computed ? s.eigenvalues(0) : NAN) << ", solves = " << s.n_operations
              << std::endl;
  }

  // Unstable case: K - c I has negative eigenvalues lambda_k - c (translations -> -c,
  // still exact eigenvectors, so they are deflated explicitly).
  if (fast.num_computed >= 3) {
    const double c = 0.5 * (fast.eigenvalues(1) + fast.eigenvalues(2)) + 1e-3;
    Eigen::SparseMatrix<double> I(n, n);
    I.setIdentity();
    const Eigen::SparseMatrix<double> K_shifted = K_fast - c * I;
    StiffnessSpectrumOptions o;
    o.n_modes = 5;
    o.deflate_translations = 1;
    o.verbose = true;
    const StiffnessSpectrum s = lowest_stiffness_modes(K_shifted, o);
    std::cout << "unstable test (K - " << c << " I): expected";
    for (int k = 0; k < 5; k++) std::cout << " " << fast.eigenvalues(k) - c;
    std::cout << "\n                                  got     ";
    for (int k = 0; k < s.num_computed; k++) std::cout << " " << s.eigenvalues(k);
    std::cout << std::endl;
  }
  std::cout << "===========================================\n" << std::endl;
}

struct MethodStats {
  long iterations = 0, nfev = 0;
  double setup = 0.0, solve = 0.0;
  int steps = 0;
  double max_energy_diff = 0.0;
  int differing_minima = 0;
};

} // namespace

int main(int argc, char **argv) {
  const Args args = parse_args(argc, argv);
  std::cout << "Preconditioner benchmark: " << args.nx << "x" << args.ny
            << ", steps=" << args.steps << ", grad_tol=" << args.grad_tol
            << ", threads=" << omp_get_max_threads() << std::endl;

  // ==================== SETUP (same as example_1_conti_zanzotto_loading) ============
  const std::string lattice_type = "square";
  Eigen::Matrix<double, 3, 2> dndx = calculateShapeDerivatives(
      Eigen::Vector2d(0.0, 0.0), Eigen::Vector2d(1.0, 0.0), Eigen::Vector2d(0.0, 1.0));

  std::function<double(double)> potential_func = square_energy;
  std::function<double(double)> potential_func_der = square_energy_der;
  std::function<double(double)> potential_func_sder = square_energy_der;

  const double lattice_constant = 1.0;
  std::vector<Point2D> square_points =
      LatticeGenerator::generate_2d_lattice(args.nx, args.ny, lattice_constant, lattice_type);
  const int original_domain_size = square_points.size();
  DomainInfo domain_size = compute_domain_size(square_points);
  const std::array<double, 2> offsets = {lattice_constant, lattice_constant};
  DomainDimensions domain_dims(domain_size.get_width(), domain_size.get_height());
  const bool pbc = true;
  auto [original_domain_map, translation_map] =
      MeshGenerator::create_domain_maps(original_domain_size, domain_dims, offsets);
  auto [interior_mapping, full_mapping] =
      create_dof_mapping_original(square_points, 0.5 * lattice_constant, pbc);
  Point2D domain_dims_point(domain_dims.size_x, domain_dims.size_y);

  AdaptiveMesher mesher(domain_dims_point, offsets, original_domain_map, translation_map,
                        full_mapping, 1e-6, pbc);
  mesher.setUsePeriodicCopies(pbc);

  alglib::real_1d_array free_dofs;
  const int n_vars = interior_mapping.size();
  free_dofs.setlength(2 * n_vars);
  map_points_to_solver_array(free_dofs, square_points, interior_mapping, n_vars);
  auto [elements, active_elements] =
      mesher.createMesh(square_points, free_dofs, Eigen::Matrix2d::Identity(), &dndx);
  for (auto &element : elements) element.set_dof_mapping(full_mapping);
  std::cout << "Nodes: " << n_vars << ", DOFs: " << 2 * n_vars
            << ", active elements: " << active_elements.size() << std::endl;

  Strain_Energy_LatticeCalculator calculator(1.0);
  const double zero =
      calculator.calculate_energy(Eigen::Matrix2d::Identity(), potential_func, 0);

  // ==================== METHODS ====================
  std::vector<std::unique_ptr<PreconditionedLBFGS>> solvers;
  for (const std::string &m : args.methods) {
    PreconditionedLBFGSOptions opt;
    opt.corrections = args.corrections;
    opt.refresh_every = args.refresh;
    opt.refresh_iterations = args.refresh_its;
    opt.verbose = false;
    // Suffix "@N" reuses the Stiffness factorization for N load steps.
    std::string base = m;
    const size_t at = base.find('@');
    if (at != std::string::npos) {
      opt.refresh_every = std::stoi(base.substr(at + 1));
      base = base.substr(0, at);
    }
    if (base == "prod") {
      opt.type = LBFGSPreconditioner::None;
      opt.grad_tol = 0.0; // ALGLIB automatic stopping, as in production
    } else {
      opt.type = parse_lbfgs_preconditioner(base);
      opt.grad_tol = args.grad_tol;
    }
    solvers.push_back(std::make_unique<PreconditionedLBFGS>(opt));
  }
  std::vector<MethodStats> stats(args.methods.size());

  std::ofstream csv(args.csv);
  csv << "step,alpha,method,iterations,nfev,termination,setup_s,solve_s,energy,grad_max,"
         "energy_minus_ref,max_dx_vs_ref,assemble_s,factorize_s,restarts\n";
  csv << std::setprecision(12);

  // ==================== LOADING LOOP ====================
  for (int step = 0; step < args.steps; step++) {
    const double alpha = args.alpha0 + step * args.dalpha;
    Eigen::Matrix2d F_ext;
    F_ext << 1.0, alpha, 0.0, 1.0;
    Eigen::Matrix2d dF_ext;
    dF_ext << 1.0, args.dalpha, 0.0, 1.0;

    if (step == 0) {
      std::mt19937 gen(args.seed);
      std::normal_distribution<double> noise_dist(0.0, 0.04);
      for (auto &p : square_points) {
        Eigen::Vector2d noise(noise_dist(gen), noise_dist(gen));
        p.coord = F_ext * p.coord + noise;
      }
    } else {
      for (auto &p : square_points) p.coord = dF_ext * p.coord;
    }
    // Periodic images follow the macroscopic deformation. In the production loop this
    // happens as a side effect of ConfigurationSaver::calculateEnergyAndStress.
    for (auto &element : elements) element.setExternalDeformation(F_ext);

    UserData userData(square_points, elements, calculator, potential_func,
                      potential_func_der, zero, lattice_constant, F_ext, interior_mapping,
                      full_mapping, active_elements, false);
    userData.second_derivative_function = &potential_func_sder;

    alglib::real_1d_array x0;
    x0.setlength(2 * n_vars);
    map_points_to_solver_array(x0, square_points, interior_mapping, n_vars);

    std::cout << "\n--- step " << step << "  alpha=" << alpha << " ---" << std::endl;
    alglib::real_1d_array x_ref;
    double e_ref = 0.0;
    for (size_t k = 0; k < solvers.size(); k++) {
      alglib::real_1d_array x;
      x.setlength(x0.length());
      for (int i = 0; i < x0.length(); i++) x[i] = x0[i];

      const PreconditionedLBFGSReport r = solvers[k]->optimize(x, &userData);

      double dx = 0.0, de = 0.0;
      if (k == 0) {
        x_ref = x;
        e_ref = r.energy;
      } else {
        for (int i = 0; i < x.length(); i++) dx = std::max(dx, std::abs(x[i] - x_ref[i]));
        de = r.energy - e_ref;
      }
      MethodStats &s = stats[k];
      s.iterations += r.iterations;
      s.nfev += r.nfev;
      s.setup += r.setup_seconds;
      s.solve += r.solve_seconds;
      s.steps++;
      s.max_energy_diff = std::max(s.max_energy_diff, std::abs(de));
      if (dx > 1e-3) s.differing_minima++;

      std::cout << std::left << std::setw(10) << args.methods[k] << std::right
                << " its=" << std::setw(6) << r.iterations << " nfev=" << std::setw(6)
                << r.nfev << " term=" << std::setw(2) << r.termination_type
                << " setup=" << std::fixed << std::setprecision(4) << r.setup_seconds
                << " (asm " << r.assemble_seconds << " + chol " << r.factorize_seconds << ")"
                << "s solve=" << r.solve_seconds << "s" << std::scientific
                << std::setprecision(6) << " E=" << r.energy << " max|g|=" << r.grad_max
                << " dE=" << de << " max|dx|=" << dx
                << (r.restarts ? " restarts=" + std::to_string(r.restarts) : "") << std::endl;
      csv << step << "," << alpha << "," << args.methods[k] << "," << r.iterations << ","
          << r.nfev << "," << r.termination_type << "," << r.setup_seconds << ","
          << r.solve_seconds << "," << r.energy << "," << r.grad_max << "," << de << ","
          << dx << "," << r.assemble_seconds << "," << r.factorize_seconds << ","
          << r.restarts << "\n";

      if (args.validate && step == 0 && k == 0) {
        validate_stiffness(x, userData, calculator, potential_func_der, potential_func_sder);
      }
    }
    csv.flush();

    // Advance with the reference (first) method's solution
    map_solver_array_to_points(x_ref, square_points, interior_mapping, n_vars);

    if (args.eig > 0 && step == args.steps - 1) {
      compare_eigen_solvers(x_ref, userData, calculator, potential_func_der,
                            potential_func_sder, args.eig);
    }
  }

  // ==================== SUMMARY ====================
  std::cout << "\n==================== SUMMARY (" << args.steps << " steps, " << args.nx
            << "x" << args.ny << ") ====================" << std::endl;
  std::cout << std::left << std::setw(11) << "method" << std::right << std::setw(10)
            << "its" << std::setw(10) << "nfev" << std::setw(11) << "setup[s]"
            << std::setw(11) << "solve[s]" << std::setw(11) << "total[s]" << std::setw(10)
            << "speedup" << std::setw(14) << "max|dE|" << std::setw(10) << "#diff"
            << std::endl;
  const double base_total = stats[0].setup + stats[0].solve;
  for (size_t k = 0; k < stats.size(); k++) {
    const MethodStats &s = stats[k];
    const double total = s.setup + s.solve;
    std::cout << std::left << std::setw(11) << args.methods[k] << std::right
              << std::setw(10) << s.iterations << std::setw(10) << s.nfev << std::fixed
              << std::setprecision(3) << std::setw(11) << s.setup << std::setw(11)
              << s.solve << std::setw(11) << total << std::setw(10) << std::setprecision(2)
              << base_total / total << std::scientific << std::setprecision(2)
              << std::setw(14) << s.max_energy_diff << std::setw(10) << s.differing_minima
              << std::endl;
  }
  std::cout << "(speedup relative to '" << args.methods[0]
            << "'; #diff = steps whose solution differs from it by max|dx| > 1e-3)"
            << std::endl;
  std::cout << "Per-step results written to " << args.csv << std::endl;
  return 0;
}
