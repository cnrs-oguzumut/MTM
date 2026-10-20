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
#include "../include/lattice_energy/EnergyFunctions.h"
#include "../include/experiments/acoustic_studies.h"
#include "../include/experiments/stiffness_shear_study.h"

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
  double noise_level = 0.04;

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
  double beta = -0.25;
  double bulk_K = 5.0;
  std::string lattice_type = "square";

  // Acoustic study flags
  bool run_acoustic_study = false;
  int acoustic_potential_mode = 3; // 1: LJ v2, 2: LJ v3, 3: Zanzotto, 4: Square Gaussian
  double acoustic_t_max = 1.0;
  double acoustic_dt = 0.02;
  double acoustic_dp = 0.05;
  double acoustic_alpha_min = -1.5;
  double acoustic_alpha_max = 1.5;
  double acoustic_dalpha = 0.005;
  std::string acoustic_path_type = "simple_shear";
  double acoustic_shear_theta = 0.0;
  std::string acoustic_dat_output = "acoustic_study_results.dat";
  std::string acoustic_json_output = "acoustic_study_results.json";
  std::string acoustic_csv_output = "energy_path.csv";
  bool only_elastic_constants = false;

  // Stiffness study flags
  bool run_stiffness_study = false;
  std::string stiffness_bc = "periodic";
  std::string stiffness_csv_output = "stiffness_shear_study.csv";
  std::string stiffness_json_output = "stiffness_shear_study.json";
  std::string stiffness_fullscale_json = "scratch/stiffness_fullscale.json";
  bool stiffness_single_alpha_only = false;
  double stiffness_single_alpha = 0.0;
  bool stiffness_unique_levels = false;
  double stiffness_alpha_min = -0.3;
  double stiffness_alpha_max = 0.3;
  double stiffness_dalpha = 0.005;
  bool stiffness_use_rotation = false;
  double stiffness_rotation_angle = 0.0;
  int stiffness_csl_m = 1;
  int stiffness_csl_n = 0;
  int stiffness_csl_kx = 1;
  int stiffness_csl_ky = 1;
  bool stiffness_initial_mesh_only = false;
  std::string stiffness_initial_mesh_json = "scratch/initial_lattice.json";

  // Scan all arguments for flags
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--acoustic-study" || arg == "--acoustic" || arg == "--parametric-acoustic") {
      run_acoustic_study = true;
    } else if (arg == "--stiffness-study" || arg == "--stiffness-shear" || arg == "--stiffness") {
      run_stiffness_study = true;
    } else if (arg.rfind("--bc=", 0) == 0) {
      stiffness_bc = arg.substr(5);
    } else if (arg.rfind("--boundary=", 0) == 0) {
      stiffness_bc = arg.substr(11);
    } else if (arg.rfind("--stiffness-csv=", 0) == 0) {
      stiffness_csv_output = arg.substr(16);
    } else if (arg.rfind("--stiffness-json=", 0) == 0) {
      stiffness_json_output = arg.substr(17);
    } else if (arg.rfind("--stiffness-fullscale-json=", 0) == 0) {
      stiffness_fullscale_json = arg.substr(27);
    } else if (arg.rfind("--single-alpha=", 0) == 0) {
      stiffness_single_alpha_only = true;
      stiffness_single_alpha = std::stod(arg.substr(15));
      run_stiffness_study = true;
    } else if (arg == "--unique-levels" || arg == "--unique-levels-only" || arg == "--unique-modes" || arg == "--no-doublets") {
      stiffness_unique_levels = true;
    } else if (arg == "--initial-mesh-only" || arg == "--initial-mesh-preview" || arg == "--preview-initial-mesh") {
      stiffness_initial_mesh_only = true;
      run_stiffness_study = true;
    } else if (arg.rfind("--initial-mesh-json=", 0) == 0) {
      stiffness_initial_mesh_json = arg.substr(20);
    } else if (arg == "--use-rotation" || arg == "--rotated") {
      stiffness_use_rotation = true;
    } else if (arg.rfind("--rotation-angle=", 0) == 0) {
      stiffness_rotation_angle = std::stod(arg.substr(17));
      stiffness_use_rotation = true;
    } else if (arg.rfind("--rotation=", 0) == 0) {
      stiffness_rotation_angle = std::stod(arg.substr(11));
      stiffness_use_rotation = true;
    } else if (arg.rfind("--csl-m=", 0) == 0) {
      stiffness_csl_m = std::stoi(arg.substr(8));
      stiffness_use_rotation = true;
    } else if (arg.rfind("--csl-n=", 0) == 0) {
      stiffness_csl_n = std::stoi(arg.substr(8));
      stiffness_use_rotation = true;
    } else if (arg.rfind("--csl-kx=", 0) == 0) {
      stiffness_csl_kx = std::stoi(arg.substr(9));
      stiffness_use_rotation = true;
    } else if (arg.rfind("--csl-ky=", 0) == 0) {
      stiffness_csl_ky = std::stoi(arg.substr(9));
      stiffness_use_rotation = true;
    } else if (arg.rfind("--nx=", 0) == 0) {
      nx = std::stoi(arg.substr(5));
    } else if (arg.rfind("--ny=", 0) == 0) {
      ny = std::stoi(arg.substr(5));
    } else if (arg.rfind("--potential-mode=", 0) == 0) {
      acoustic_potential_mode = std::stoi(arg.substr(17));
    } else if (arg.rfind("--potential=", 0) == 0) {
      acoustic_potential_mode = std::stoi(arg.substr(12));
    } else if (arg.rfind("--acoustic-mode=", 0) == 0) {
      acoustic_potential_mode = std::stoi(arg.substr(16));
    } else if (arg.rfind("--t-max=", 0) == 0) {
      acoustic_t_max = std::stod(arg.substr(8));
    } else if (arg.rfind("--dt=", 0) == 0) {
      acoustic_dt = std::stod(arg.substr(5));
    } else if (arg.rfind("--dp=", 0) == 0) {
      acoustic_dp = std::stod(arg.substr(5));
    } else if (arg.rfind("--shear-alpha-min=", 0) == 0) {
      acoustic_alpha_min = std::stod(arg.substr(18));
      stiffness_alpha_min = acoustic_alpha_min;
    } else if (arg.rfind("--shear-alpha-max=", 0) == 0) {
      acoustic_alpha_max = std::stod(arg.substr(18));
      stiffness_alpha_max = acoustic_alpha_max;
    } else if (arg.rfind("--shear-dalpha=", 0) == 0) {
      acoustic_dalpha = std::stod(arg.substr(15));
      stiffness_dalpha = acoustic_dalpha;
    } else if (arg.rfind("--alpha-min=", 0) == 0) {
      acoustic_alpha_min = std::stod(arg.substr(12));
      stiffness_alpha_min = acoustic_alpha_min;
    } else if (arg.rfind("--alpha-max=", 0) == 0) {
      acoustic_alpha_max = std::stod(arg.substr(12));
      stiffness_alpha_max = acoustic_alpha_max;
    } else if (arg.rfind("--dalpha=", 0) == 0) {
      acoustic_dalpha = std::stod(arg.substr(9));
      stiffness_dalpha = acoustic_dalpha;
    } else if (arg.rfind("--shear-path-type=", 0) == 0) {
      acoustic_path_type = arg.substr(18);
    } else if (arg.rfind("--shear-path=", 0) == 0) {
      acoustic_path_type = arg.substr(13);
    } else if (arg.rfind("--path-type=", 0) == 0) {
      acoustic_path_type = arg.substr(12);
    } else if (arg.rfind("--shear-theta=", 0) == 0) {
      acoustic_shear_theta = std::stod(arg.substr(14));
    } else if (arg.rfind("--shear-angle=", 0) == 0) {
      acoustic_shear_theta = std::stod(arg.substr(14));
    } else if (arg.rfind("--csv-output=", 0) == 0) {
      acoustic_csv_output = arg.substr(13);
    } else if (arg.rfind("--dat-output=", 0) == 0) {
      acoustic_dat_output = arg.substr(13);
    } else if (arg.rfind("--json-output=", 0) == 0) {
      acoustic_json_output = arg.substr(14);
    } else if (arg == "--print-elastic-constants" || arg == "--elastic-constants") {
      run_acoustic_study = true;
      only_elastic_constants = true;
    } else if (arg == "--no-remesh" || arg == "--noremesh" || arg == "--without-remesh" || arg == "-nr") {
      enable_remeshing = false;
    } else if (arg == "--remesh" || arg == "-r") {
      enable_remeshing = true;
    } else if (arg.rfind("--restart=", 0) == 0) {
      restart_checkpoint = arg.substr(10);
    } else if (arg == "--restart" && i + 1 < argc) {
      restart_checkpoint = argv[++i];
    } else if (arg.rfind("--beta=", 0) == 0) {
      beta = std::stod(arg.substr(7));
    } else if (arg.rfind("--K=", 0) == 0) {
      bulk_K = std::stod(arg.substr(4));
    } else if (arg.rfind("--bulk-K=", 0) == 0) {
      bulk_K = std::stod(arg.substr(9));
    } else if (arg.rfind("--lattice=", 0) == 0) {
      lattice_type = arg.substr(10);
    } else if (arg.rfind("--symmetry=", 0) == 0) {
      lattice_type = arg.substr(11);
    } else if (arg.rfind("--lattice-symmetry=", 0) == 0) {
      lattice_type = arg.substr(19);
    } else if (arg == "--triangular") {
      lattice_type = "triangular";
    } else if (arg == "--square") {
      lattice_type = "square";
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
    } else if (arg.rfind("--noise-level=", 0) == 0) {
      noise_level = std::stod(arg.substr(14));
    } else if (arg.rfind("--noise=", 0) == 0) {
      noise_level = std::stod(arg.substr(8));
    }
  }

  // Set Zanzotto potential parameters
  energy_functions::set_parameters(beta, bulk_K);

  configure_relaxation_solver(relax_options, precond_from_step);
  configure_stability_monitor(stability_options);

  if (!restart_checkpoint.empty()) {
    restart_zanzotto_simulation(restart_checkpoint, checkpoint_interval, stress_drop_threshold, save_triangle_data, max_avalanches, lattice_type);
    return 0;
  }

  if (argc >= 3 && argv[1][0] != '-') {
    nx = std::atoi(argv[1]);
    ny = std::atoi(argv[2]);
  }
  if (argc >= 4 && argv[3][0] != '-') {
    mode = argv[3];
  }
  if (argc >= 5 && argv[4][0] != '-') {
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
            << " | lattice=" << lattice_type
            << " | beta=" << beta << " | K=" << bulk_K
            << " | mode=" << mode << " | seed=" << seed
            << " | noise=" << noise_level
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
              << ", ny=" << ny << ", lattice=" << lattice_type
              << ", potential_mode=" << acoustic_potential_mode
              << ", R_free=" << r_free
              << ", mesh=" << (export_full_mesh ? "full" : "circle") << " <<<\n" << std::endl;
    single_dislocation_cylinder_relaxation(0, nx, ny, r_free, enable_remeshing, export_full_mesh, lattice_type, acoustic_potential_mode);
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
  if (run_acoustic_study || mode == "acoustic" || mode == "acoustic_study" || mode == "parametric_acoustic") {
    AcousticStudyConfig acoustic_cfg;
    acoustic_cfg.mode = acoustic_potential_mode;
    acoustic_cfg.lattice_type = lattice_type;
    acoustic_cfg.beta = beta;
    acoustic_cfg.K = bulk_K;
    acoustic_cfg.t_max = acoustic_t_max;
    acoustic_cfg.dt = acoustic_dt;
    acoustic_cfg.dp = acoustic_dp;
    acoustic_cfg.dat_output = acoustic_dat_output;
    acoustic_cfg.json_output = acoustic_json_output;
    acoustic_cfg.only_elastic_constants = only_elastic_constants;
    acoustic_cfg.path_type = acoustic_path_type;
    acoustic_cfg.shear_theta_deg = acoustic_shear_theta;
    acoustic_cfg.csv_output = acoustic_csv_output;
    acoustic_cfg.alpha_min = acoustic_alpha_min;
    acoustic_cfg.alpha_max = acoustic_alpha_max;
    acoustic_cfg.dalpha = acoustic_dalpha;

    runParametricAcousticStudy(acoustic_cfg);
    return 0;
  }
  // parametricAcousticStudy();
  // parametricAcousticStudy_v2();

  // 4b. Stiffness Matrix K(α) Shear Study:
  //     Calculates minimum eigenvalues of global stiffness matrix K as a function of alpha.
  if (run_stiffness_study || mode == "stiffness" || mode == "stiffness_study" || mode == "stiffness_shear") {
    StiffnessShearConfig stiff_cfg;
    stiff_cfg.nx = nx;
    stiff_cfg.ny = ny;
    stiff_cfg.lattice_type = lattice_type;
    stiff_cfg.potential_mode = acoustic_potential_mode;
    stiff_cfg.bc = stiffness_bc;
    stiff_cfg.alpha_min = stiffness_alpha_min;
    stiff_cfg.alpha_max = stiffness_alpha_max;
    stiff_cfg.dalpha = stiffness_dalpha;
    stiff_cfg.beta = beta;
    stiff_cfg.K = bulk_K;
    stiff_cfg.csv_output = stiffness_csv_output;
    stiff_cfg.json_output = stiffness_json_output;
    stiff_cfg.fullscale_json_output = stiffness_fullscale_json;
    stiff_cfg.single_alpha_only = stiffness_single_alpha_only;
    stiff_cfg.single_alpha = stiffness_single_alpha;
    stiff_cfg.unique_levels_only = stiffness_unique_levels;
    stiff_cfg.use_rotation = stiffness_use_rotation;
    stiff_cfg.rotation_angle = stiffness_rotation_angle;
    stiff_cfg.csl_m = stiffness_csl_m;
    stiff_cfg.csl_n = stiffness_csl_n;
    stiff_cfg.csl_kx = stiffness_csl_kx;
    stiff_cfg.csl_ky = stiffness_csl_ky;
    stiff_cfg.initial_mesh_only = stiffness_initial_mesh_only;
    stiff_cfg.initial_mesh_json = stiffness_initial_mesh_json;

    runStiffnessShearStudy(stiff_cfg);
    return 0;
  }

  // 5. Nano-Indentation:
  //    Simulates circular nano-indenter pushing into crystal lattice.
  // indentation();

  // 6. Stress-Controlled Loading:
  //    Applies stress-controlled shear/axial loading with adaptive remeshing.
  // example_3_stress_controlled_final_clean(0, nx, ny);

  // 7. Zanzotto Continuous Shear Loading:
  //    Configure loading schedule once; automatically negated for negative loading
  if (mode == "positive") {
    example_1_conti_zanzotto_loading(0, nx, ny, alpha_start, alpha_end, step_size, 0.0, seed, enable_remeshing, nullptr, checkpoint_interval, stress_drop_threshold, save_triangle_data, max_avalanches, lattice_type, noise_level);
  } else {
    example_1_conti_zanzotto_loading(0, nx, ny, -alpha_start, -alpha_end, -step_size, triangulation_perturbation, seed, enable_remeshing, nullptr, checkpoint_interval, stress_drop_threshold, save_triangle_data, max_avalanches, lattice_type, noise_level);
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
