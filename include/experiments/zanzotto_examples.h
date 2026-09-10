#pragma once

void memory(int caller_id, int nx, int ny, int restart_iteration);
void example_1_conti_zanzotto(int caller_id, int nx, int ny);
void example_1_conti_zanzotto_loading(
    int caller_id, int nx, int ny,
    double alpha_min = 0.14,
    double alpha_max = 1.0,
    double step_size = 6e-5,
    double triangulation_perturbation = 0.0,
    unsigned int seed = 42,
    bool enable_remeshing = true);
void example_1_conti_zanzotto_negative_loading(
    int caller_id, int nx, int ny,
    double alpha_min = -0.14,
    double alpha_max = -1.0,
    double step_size = -6e-5,
    unsigned int seed = 42,
    bool enable_remeshing = true);
void example_2_conti_zanzotto_triangular();

