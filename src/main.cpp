#include <cstdlib>
#include <iostream>

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

  std::cout << "System size: nx=" << nx << ", ny=" << ny
            << " | mode=" << mode << " | seed=" << seed << std::endl;

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
  //    Runs either positive loading (alpha: +0.14 -> +1.0, default orientation)
  //    or negative loading (alpha: -0.14 -> -1.0, -1e-7 orientation perturbation)
  //    using the exact same deterministic random seed for initial noise.
  if (mode == "positive") {
    example_1_conti_zanzotto_loading(0, nx, ny, 0.14, 1.0, 6e-5, 0.0, seed);
  } else {
    example_1_conti_zanzotto_negative_loading(0, nx, ny, seed);
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
