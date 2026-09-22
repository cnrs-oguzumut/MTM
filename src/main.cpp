#include <cstdlib>
#include <iostream>

#include "../include/experiments/common_simulation_helpers.h"
#include "../include/experiments/dislocation_indentation.h"
#include "../include/experiments/dislocation_study.h"
#include "../include/experiments/shift_vertical_horizontal.h"
#include "../include/experiments/shifted_crystal_study.h"
#include "../include/experiments/zanzotto_examples.h"
#include "../include/experiments/stress_controlled_examples.h"
#include "../include/experiments/data_analysis.h"
#include "../include/experiments/stability_monitor.h"
#include "../include/experiments/shifted_dislocation_study.h"
#include "../include/output/AvalancheRecorder.h"

int main(int argc, char **argv) {
  // Default system size and parameters
  int nx = 150;
  int ny = 150;
  std::string mode = "negative"; // "positive" or "negative"
  unsigned int seed = 42;
  bool enable_remeshing = true;
  double r_free = 20.0;
  bool use_cylinder = false;
  bool export_full_mesh = true;
  double alpha_start = 0.14;
  double alpha_end = 1.0;
  double step_size = 6e-5;

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

  // Eigenvalue solver of analyze_data_from_folder (example 11 below):
  //   --eig-solver=fast     analytic K + Cholesky shift-invert (default)
  //   --eig-solver=legacy   ITensor K + Spectra/SparseLU (FEMHessianAssembler)
  StiffnessEigenSolver eig_solver = StiffnessEigenSolver::Fast;

  // Stability monitor during the loading (off by default):
  //   --eig-every=N           lowest eigenvalues of K at every N-th relaxed state
  //   --eig-at-avalanche=1    also at the last stable state before and the state after
  //                           each avalanche, with the soft modes written to eigen_modes/
  //   --eig-modes=5           eigenvalues per computation (eigen_log.csv)
  //   --eig-vectors=2         soft modes written per avalanche (never translations)
  //   --eig-refine=0.2        every step while lambda_min < 0.2 x its value after the last
  //                           avalanche (0 = off)
  //   --eig-refine-ahead=2    every step while lambda_min^2, extrapolated linearly, reaches
  //                           zero within 2 x N steps (0 = off)
  //   --eig-retro=4           at each instability (stress jump or avalanche), also the last
  //                           4 relaxed states before it (kept in memory)
  StabilityMonitorOptions stability_options;
  std::string restart_checkpoint;

  // Data saving and checkpoint controls:
  //   --checkpoint-interval=N (or --chk-interval=N)  periodic elastic checkpoint every N steps (default: 500, 0 = disabled)
  //   --stress-drop-threshold=VAL (or --stress-drop=VAL) fractional stress drop to trigger avalanche save (default: 0.10)
  //   --triangle-data / --no-triangle-data           write legacy triangle_data/points_*.dat and elements_*.dat (default: disabled)
  int checkpoint_interval = 500;
  double stress_drop_threshold = 0.10;
  bool save_triangle_data = false;
  int max_avalanches = 0;
  double triangulation_perturbation = -1e-7; // default for negative loading

  // Scan all arguments for flags
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--no-remesh" || arg == "--noremesh" || arg == "--without-remesh" || arg == "-nr") {
      enable_remeshing = false;
    } else if (arg == "--remesh" || arg == "-r") {
      enable_remeshing = true;
    } else if (arg.rfind("--restart=", 0) == 0) {
      restart_checkpoint = arg.substr(10);
    } else if (arg == "--restart" && i + 1 < argc) {
      restart_checkpoint = argv[++i];
    } else if (arg.rfind("--precond=", 0) == 0) {
      relax_options.type = parse_lbfgs_preconditioner(arg.substr(10));
    } else if (arg.rfind("--precond-tol=", 0) == 0) {
      relax_options.grad_tol = std::stod(arg.substr(14));
    } else if (arg.rfind("--precond-refresh=", 0) == 0) {
      relax_options.refresh_every = std::stoi(arg.substr(18));
    } else if (arg.rfind("--precond-from-step=", 0) == 0) {
      precond_from_step = std::stoi(arg.substr(20));
    } else if (arg.rfind("--eig-solver=", 0) == 0) {
      eig_solver = parse_stiffness_eigen_solver(arg.substr(13));
    } else if (arg.rfind("--eig-every=", 0) == 0) {
      stability_options.every = std::stoi(arg.substr(12));
    } else if (arg.rfind("--eig-at-avalanche=", 0) == 0) {
      stability_options.at_avalanche = std::stoi(arg.substr(19)) != 0;
    } else if (arg.rfind("--eig-modes=", 0) == 0) {
      stability_options.modes = std::stoi(arg.substr(12));
    } else if (arg.rfind("--eig-vectors=", 0) == 0) {
      stability_options.vectors = std::stoi(arg.substr(14));
    } else if (arg.rfind("--eig-refine=", 0) == 0) {
      stability_options.refine = std::stod(arg.substr(13));
    } else if (arg.rfind("--eig-refine-ahead=", 0) == 0) {
      stability_options.refine_ahead = std::stod(arg.substr(19));
    } else if (arg.rfind("--eig-retro=", 0) == 0) {
      stability_options.retro = std::stoi(arg.substr(12));
    } else if (arg == "--shift" || arg == "--staircase") {
      mode = "shift";
    } else if (arg.rfind("--r-free=", 0) == 0) {
      r_free = std::stod(arg.substr(9));
      use_cylinder = true;
    } else if (arg == "--cylinder" || arg == "--cylinder=1") {
      use_cylinder = true;
    } else if (arg == "--volterra" || arg == "--dislocation") {
      mode = "volterra";
    } else if (arg == "--full-mesh" || arg == "--all-elements") {
      export_full_mesh = true;
    } else if (arg == "--circle" || arg == "--circle-only" || arg == "--active-only") {
      export_full_mesh = false;
    } else if (arg == "--trace-avalanche" || arg == "--trace-avalanche=1" || arg == "--trace=1") {
      AvalancheRecorder::instance().setEnabled(true);
    } else if (arg == "--trace-avalanche=0" || arg == "--trace=0") {
      AvalancheRecorder::instance().setEnabled(false);
    } else if (arg == "--save-surgery-vtk" || arg == "--save-surgery-vtk=1" ||
               arg == "--surgery-vtk" || arg == "--surgery-vtk=1" ||
               arg == "--save-remesh-vtk" || arg == "--save-remesh-vtk=1") {
      AvalancheRecorder::instance().setSaveSurgeryVTK(true);
    } else if (arg == "--save-surgery-vtk=0" || arg == "--surgery-vtk=0" || arg == "--no-surgery-vtk") {
      AvalancheRecorder::instance().setSaveSurgeryVTK(false);
    } else if (arg.rfind("--alpha-end=", 0) == 0) {
      alpha_end = std::stod(arg.substr(12));
    } else if (arg.rfind("--alpha-max=", 0) == 0) {
      alpha_end = std::stod(arg.substr(12));
    } else if (arg.rfind("--alpha-start=", 0) == 0) {
      alpha_start = std::stod(arg.substr(14));
    } else if (arg.rfind("--alpha-min=", 0) == 0) {
      alpha_start = std::stod(arg.substr(12));
    } else if (arg.rfind("--step-size=", 0) == 0) {
      step_size = std::stod(arg.substr(12));
    } else if (arg.rfind("--checkpoint-interval=", 0) == 0) {
      checkpoint_interval = std::stoi(arg.substr(22));
    } else if (arg.rfind("--chk-interval=", 0) == 0) {
      checkpoint_interval = std::stoi(arg.substr(15));
    } else if (arg.rfind("--stress-drop-threshold=", 0) == 0) {
      stress_drop_threshold = std::stod(arg.substr(24));
    } else if (arg.rfind("--stress-drop=", 0) == 0) {
      stress_drop_threshold = std::stod(arg.substr(14));
    } else if (arg == "--save-triangle-data" || arg == "--triangle-data" || arg == "--save-triangle-data=1" || arg == "--triangle-data=1") {
      save_triangle_data = true;
    } else if (arg == "--no-triangle-data" || arg == "--save-triangle-data=0" || arg == "--triangle-data=0") {
      save_triangle_data = false;
    } else if (arg.rfind("--max-avalanches=", 0) == 0) {
      max_avalanches = std::stoi(arg.substr(17));
    } else if (arg.rfind("--avalanches=", 0) == 0) {
      max_avalanches = std::stoi(arg.substr(13));
    } else if (arg == "--no-perturbation" || arg == "--no-perturb") {
      triangulation_perturbation = 0.0;
    } else if (arg.rfind("--perturbation=", 0) == 0) {
      triangulation_perturbation = std::stod(arg.substr(15));
    }
  }
  configure_relaxation_solver(relax_options, precond_from_step);
  configure_stability_monitor(stability_options);

  if (!restart_checkpoint.empty()) {
    restart_zanzotto_simulation(restart_checkpoint, checkpoint_interval, stress_drop_threshold, save_triangle_data, max_avalanches);
    return 0;
  }

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
            << " | alpha_start=" << alpha_start << " | alpha_end=" << alpha_end
            << " | step_size=" << step_size
            << " | trace=" << (AvalancheRecorder::instance().isEnabled() ? "enabled" : "disabled")
            << " | surgery_vtk=" << (AvalancheRecorder::instance().isSurgeryVTKEnabled() ? "enabled" : "disabled");
  if (max_avalanches > 0)
    std::cout << " | max_avalanches=" << max_avalanches;
  std::cout << " | precond=" << to_string(relax_options.type);
  if (relax_options.type != LBFGSPreconditioner::None)
    std::cout << " (grad_tol=" << relax_options.grad_tol
              << ", refresh=" << relax_options.refresh_every
              << ", from step " << precond_from_step << ")";
  std::cout << std::endl;
  std::cout << "Saving controls: chk_interval=";
  if (checkpoint_interval > 0) std::cout << checkpoint_interval << " steps";
  else std::cout << "disabled (avalanches only)";
  std::cout << " | stress_drop_threshold=" << (stress_drop_threshold * 100.0) << "%"
            << " | triangle_data=" << (save_triangle_data ? "enabled" : "disabled")
            << std::endl;
  if (stability_options.enabled()) {
    std::cout << "Stability monitor: every " << stability_options.every << " steps"
              << ", at avalanches " << (stability_options.at_avalanche ? "on" : "off")
              << ", modes=" << stability_options.modes
              << ", vectors=" << stability_options.vectors
              << ", refine=" << stability_options.refine
              << ", refine-ahead=" << stability_options.refine_ahead
              << ", retro=" << stability_options.retro << std::endl;
  }

  // =========================================================================
  // LIST OF EXPERIMENT EXAMPLES
  // Uncomment the desired example to run:
  // =========================================================================

  // 1. Shifting Examples:
  //    Simulates vertical and horizontal shifts of crystal blocks with relaxation.
  if (mode == "shift" || mode == "final_shift" || mode == "staircase") {
    std::cout << "\n>>> Running staircase shifting tests with nx=" << nx
              << ", ny=" << ny << " <<<\n" << std::endl;
    run_final_shift_tests(nx, ny, 100, 100);
    return 0;
  }

  // 2. Dislocation Studies:
  //    Simulates single dislocation in a cylinder / boundary-fixed crystal.
  if (mode == "volterra" || mode == "cylinder" || mode == "dislocation") {
    std::cout << "\n>>> Running single dislocation cylinder relaxation with nx=" << nx
              << ", ny=" << ny << ", R_free=" << r_free
              << ", mesh=" << (export_full_mesh ? "full" : "circle") << " <<<\n" << std::endl;
    single_dislocation_cylinder_relaxation(0, nx, ny, r_free, enable_remeshing, export_full_mesh);
    return 0;
  }

  // 3. Multi-Shift Dislocation Study (no remeshing):
  if (mode == "shifted_dislocation" || mode == "shift_dislocation" || mode == "shifts") {
    run_shifted_dislocation_study(0, nx, ny, {0, 1, 2, 3, 4, 5}, use_cylinder, r_free);
    return 0;
  }
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
  if (mode == "positive") {
    example_1_conti_zanzotto_loading(0, nx, ny, alpha_start, alpha_end, step_size, 0.0, seed, enable_remeshing, nullptr, checkpoint_interval, stress_drop_threshold, save_triangle_data, max_avalanches);
  } else {
    example_1_conti_zanzotto_loading(0, nx, ny, -alpha_start, -alpha_end, -step_size, triangulation_perturbation, seed, enable_remeshing, nullptr, checkpoint_interval, stress_drop_threshold, save_triangle_data, max_avalanches);
  }

  // 9. Zanzotto Continuous Loading (Triangular Lattice):
  //    Shear loading on a triangular lattice with adaptive remeshing.
  // example_2_conti_zanzotto_triangular();

  // 10. Restart / Memory Loading:
  //    Restarts a continuous loading simulation from a previous saved configuration.
  // memory(0, nx, ny, /*restart_iteration=*/3);

  // 11. Data Post-Processing / Analysis:
  //     Analyzes configurations and dislocation data from saved folder.
  // analyze_data_from_folder(0, nx, ny, 3301, 3651, 100, eig_solver);

  return 0;
}
