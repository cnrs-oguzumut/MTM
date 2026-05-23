#include <cstdlib>

#include "../include/experiments/shift_vertical_horizontal.h"
#include "../include/experiments/dislocation_indentation.h"

int main(int argc, char **argv) {
  // TensorExample exple;
  // example.run();
  // exit(0);
  // for (int caller_id = 0; caller_id <= 0; ++caller_id) {
  // example_1_shifting(0,20,20);
  // single_dislo_LJ();

  // parametricAcousticStudy();./
  // parametricAcousticStudy();
  // memory(0,100,100,3);
  // example_3_stress_controlled_final_clean(0,100,100);
  // analyze_data_from_folder(0, 100,100,3301,3651,100);
  //  }
  //   indentation();
  int nx = 20;
  int ny = 20;
  int horizontal_steps = 100;
  int vertical_steps = 100;
  if (argc >= 3) {
    nx = std::atoi(argv[1]);
    ny = std::atoi(argv[2]);
  }
  if (argc >= 5) {
    horizontal_steps = std::atoi(argv[3]);
    vertical_steps = std::atoi(argv[4]);
  }

  run_final_shift_tests(nx, ny, horizontal_steps, vertical_steps);
  // parametricAcousticStudy();
  //  parametricAcousticStudy_v2();

  exit(0);

  indentation();
  return 0;
}
