#include <cstdlib>
#include <iostream>

#include "../include/experiments/common_simulation_helpers.h"
#include "../include/experiments/dislocation_indentation.h"
#include "../include/experiments/dislocation_study.h"
#include "../include/experiments/shift_vertical_horizontal.h"
#include "../include/experiments/shifted_crystal_study.h"
#include "../include/experiments/zanzotto_examples.h"
#include "../include/experiments/stress_controlled_examples.h"
#include "../include/experiments/acoustic_studies.h"
#include "../include/experiments/data_analysis.h"

int main(int argc, char **argv) {
  // Default system size and parameters
  int nx = 150;
  int ny = 150;
  std::string mode = "negative"; // "positive" or "negative"
  unsigned int seed = 42;
  bool enable_remeshing = true;

  // Energy relaxation solver (default: plain L-BFGS, unchanged behaviour).
  //   --precond=stiffness|laplacian|diag|none   L-BFGS preconditioner
  //   --precond-tol=1e-6                        stop when max|dE/dx| < tol
  //   --precond-refresh=N                       reuse the stiffness factorization N steps
  //   --precond-from-step=N                     plain L-BFGS for load steps < N
  //                                             (1 keeps the plain initial relaxation)
  PreconditionedLBFGSOptions relax_options;
  int precond_from_step = 0;
  relax_options.type = LBFGSPreconditioner::None;
  relax_options.grad_tol = 1e-6;

  // Scan all arguments for flags
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--no-remesh" || arg == "--noremesh" || arg == "--without-remesh" || arg == "-nr") {
      enable_remeshing = false;
    } else if (arg == "--remesh" || arg == "-r") {
      enable_remeshing = true;
    } else if (arg.rfind("--precond=", 0) == 0) {
      relax_options.type = parse_lbfgs_preconditioner(arg.substr(10));
    } else if (arg.rfind("--precond-tol=", 0) == 0) {
      relax_options.grad_tol = std::stod(arg.substr(14));
    } else if (arg.rfind("--precond-refresh=", 0) == 0) {
      relax_options.refresh_every = std::stoi(arg.substr(18));
    } else if (arg.rfind("--precond-from-step=", 0) == 0) {
      precond_from_step = std::stoi(arg.substr(20));
    }
  }
  configure_relaxation_solver(relax_options, precond_from_step);

  if (argc >= 3) {
    nx = std::atoi(argv[1]);
    ny = std::atoi(argv[2]);
  }
  if (argc >= 4) {
    mode = argv[3];
  }
  if (argc >= 5) {
    seed = static_cast<unsigned int>(std::atoi(argv[4]));
  }
  if (argc >= 6) {
    std::string remesh_arg = argv[5];
    if (remesh_arg == "0" || remesh_arg == "false" || remesh_arg == "False" ||
        remesh_arg == "no" || remesh_arg == "no-remesh" || remesh_arg == "noremesh" ||
        remesh_arg == "without-remesh" || remesh_arg == "without-remeshing" ||
        remesh_arg == "off") {
      enable_remeshing = false;
    } else if (remesh_arg == "1" || remesh_arg == "true" || remesh_arg == "True" ||
               remesh_arg == "yes" || remesh_arg == "remesh" || remesh_arg == "on") {
      enable_remeshing = true;
    }
  }

  std::cout << "System size: nx=" << nx << ", ny=" << ny
            << " | mode=" << mode << " | seed=" << seed
            << " | remeshing=" << (enable_remeshing ? "enabled" : "disabled")
            << " | precond=" << to_string(relax_options.type);
  if (relax_options.type != LBFGSPreconditioner::None)
    std::cout << " (grad_tol=" << relax_options.grad_tol
              << ", refresh=" << relax_options.refresh_every
              << ", from step " << precond_from_step << ")";
  std::cout << std::endl;

  // =========================================================================
  // LIST OF EXPERIMENT EXAMPLES
  // Uncomment the desired example to run:
  // =========================================================================

  // 1. Shifting Examples:
  //    Simulates vertical and horizontal shifts of crystal blocks with relaxation.
  // run_final_shift_tests(nx, ny, 100, 100);

  // 2. Dislocation Studies:
  //    Simulates single dislocation nucleation and relaxation on square lattice.
  // single_dislocation_study(0, nx, ny);

  // 3. Shifted Upper Crystal Study:
  //    Studies relaxation under shifted upper crystal boundary conditions.
  // shifted_upper_crystal_study(0, nx, ny);

  // 4. Acoustic Studies:
  //    Calculates and logs acoustic tensor properties / stability across states.
  // parametricAcousticStudy();
  // parametricAcousticStudy_v2();

  // 5. Nano-Indentation:
  //    Simulates circular nano-indenter pushing into crystal lattice.
  // indentation();

  // 6. Stress-Controlled Loading:
  //    Applies stress-controlled shear/axial loading with adaptive remeshing.
  // example_3_stress_controlled_final_clean(0, nx, ny);

  // 7. Zanzotto Continuous Shear Loading:
  //    Configure loading schedule once; automatically negated for negative loading
  double alpha_start = 0.14;
  double alpha_end = 1.0;
  double step_size = 6e-5;

  if (mode == "positive") {
    example_1_conti_zanzotto_loading(0, nx, ny, alpha_start, alpha_end, step_size, 0.0, seed, enable_remeshing);
  } else {
    example_1_conti_zanzotto_negative_loading(0, nx, ny, -alpha_start, -alpha_end, -step_size, seed, enable_remeshing);
  }

  // 9. Zanzotto Continuous Loading (Triangular Lattice):
  //    Shear loading on a triangular lattice with adaptive remeshing.
  // example_2_conti_zanzotto_triangular();

  // 10. Restart / Memory Loading:
  //    Restarts a continuous loading simulation from a previous saved configuration.
  // memory(0, nx, ny, /*restart_iteration=*/3);

  // 11. Data Post-Processing / Analysis:
  //     Analyzes configurations and dislocation data from saved folder.
  // analyze_data_from_folder(0, nx, ny, 3301, 3651, 100);

  return 0;
}
