#pragma once

#include <string>

void run_shift_staircase_simulation(
    int caller_id, int nx, int ny,
    int horizontal_steps = 2,
    int vertical_steps = 2,
    double amplitude_lattice_spacings = 2.0,
    bool remeshing_enabled = true,
    bool perturb_initial_triangulation = false,
    const std::string &lattice_type = "square",
    int potential_mode = 3,
    bool use_rotation = false,
    double rotation_angle = 0.0,
    bool save_triangle_data = false);

void run_final_shift_tests(
    int nx, int ny, int horizontal_steps = 2,
    int vertical_steps = 2,
    const std::string &lattice_type = "square",
    int potential_mode = 3,
    bool use_rotation = false,
    double rotation_angle = 0.0,
    bool save_triangle_data = false);


