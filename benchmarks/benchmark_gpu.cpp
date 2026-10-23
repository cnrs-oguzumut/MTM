// benchmark_gpu.cpp
//
// Tests and benchmarks CPU (minimize_energy_with_triangles) vs.
// GPU (minimize_energy_with_triangles_GPU) on Apple Silicon Metal.
//
// Validates:
//   1. Accuracy of total energy func (|func_CPU - func_GPU|).
//   2. Accuracy of gradient vector grad (|grad_CPU - grad_GPU|_inf).
//   3. Evaluation throughput (time per evaluation and speedup).
//   4. Optimization convergence with ALGLIB L-BFGS.

#include "../include/experiments/experiment_includes.h"
#include "../include/optimization/MetalEnergyMinimizer.h"
#include "../include/optimization/LatticeOptimizer.h"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

namespace {

using clock_type = std::chrono::high_resolution_clock;

double elapsed_ms(clock_type::time_point start, clock_type::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

} // namespace

int main(int argc, char** argv) {
    int nx = 60;
    int ny = 60;
    int n_evals = 100;
    bool run_lbfgs = true;

    if (argc >= 3) {
        nx = std::atoi(argv[1]);
        ny = std::atoi(argv[2]);
    }
    if (argc >= 4) {
        n_evals = std::atoi(argv[3]);
    }

    std::cout << "==========================================================" << std::endl;
    std::cout << "   Apple Silicon Metal GPU Acceleration Benchmark        " << std::endl;
    std::cout << "   Mesh Size: " << nx << " x " << ny << ", Timing Runs: " << n_evals << std::endl;
    std::cout << "==========================================================" << std::endl;

    if (!metal_gpu::is_available()) {
        std::cerr << "[ERROR] Metal GPU device is NOT available on this system!" << std::endl;
        return 1;
    }
    std::cout << "[OK] Apple Metal GPU backend is active." << std::endl;

    // 1. Setup Mesh & Initial Configuration
    const std::string lattice_type = "square";
    Eigen::Matrix<double, 3, 2> dndx = calculateShapeDerivatives(
        Eigen::Vector2d(0.0, 0.0), Eigen::Vector2d(1.0, 0.0), Eigen::Vector2d(0.0, 1.0));

    std::function<double(double)> potential_func = square_energy;
    std::function<double(double)> potential_func_der = square_energy_der;

    const double lattice_constant = 1.0;
    std::vector<Point2D> square_points =
        LatticeGenerator::generate_2d_lattice(nx, ny, lattice_constant, lattice_type);
    const int original_domain_size = static_cast<int>(square_points.size());
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

    const int n_vars = static_cast<int>(interior_mapping.size());
    alglib::real_1d_array free_dofs;
    free_dofs.setlength(2 * n_vars);
    map_points_to_solver_array(free_dofs, square_points, interior_mapping, n_vars);

    // Apply a periodic sinusoidal deformation to free DOFs to generate non-uniform strain & gradients
    const double L_y = domain_dims.size_y;
    for (int i = 0; i < n_vars; ++i) {
        double y_val = free_dofs[n_vars + i];
        free_dofs[i] += 0.05 * std::sin(2.0 * M_PI * y_val / L_y);
    }

    // Apply external deformation gradient (identity for pure internal perturbation)
    Eigen::Matrix2d F_ext = Eigen::Matrix2d::Identity();

    auto [elements, active_elements] =
        mesher.createMesh(square_points, free_dofs, F_ext, &dndx);
    for (auto &element : elements) element.set_dof_mapping(full_mapping);

    Strain_Energy_LatticeCalculator calculator(1.0);
    const double zero = calculator.calculate_energy(Eigen::Matrix2d::Identity(), potential_func, 0);

    UserData userData(square_points, elements, calculator, potential_func,
                      potential_func_der, zero, 1.0, F_ext,
                      interior_mapping, full_mapping, active_elements);

    std::cout << "Mesh Info:" << std::endl;
    std::cout << "  Nodes:           " << square_points.size() << std::endl;
    std::cout << "  Free DOFs:       " << 2 * n_vars << std::endl;
    std::cout << "  Active Elements: " << active_elements.size() << std::endl;
    std::cout << "  OpenMP Threads:  " << omp_get_max_threads() << std::endl;
    std::cout << "----------------------------------------------------------" << std::endl;

    // 2. Accuracy Comparison
    double f_cpu = 0.0;
    alglib::real_1d_array g_cpu;
    g_cpu.setlength(2 * n_vars);
    minimize_energy_with_triangles(free_dofs, f_cpu, g_cpu, &userData);

    double f_gpu = 0.0;
    alglib::real_1d_array g_gpu;
    g_gpu.setlength(2 * n_vars);
    minimize_energy_with_triangles_GPU(free_dofs, f_gpu, g_gpu, &userData);

    double energy_abs_diff = std::abs(f_cpu - f_gpu);
    double energy_scale = std::max(std::abs(f_cpu), 1.0);
    double energy_rel_diff = energy_abs_diff / energy_scale;

    double max_g_diff = 0.0;
    double max_g_val = 0.0;
    for (int i = 0; i < 2 * n_vars; ++i) {
        double diff = std::abs(g_cpu[i] - g_gpu[i]);
        if (diff > max_g_diff) max_g_diff = diff;
        if (std::abs(g_cpu[i]) > max_g_val) max_g_val = std::abs(g_cpu[i]);
    }
    double g_scale = std::max(max_g_val, 1.0);
    double g_rel_diff = max_g_diff / g_scale;

    std::cout << "\n--- Accuracy Validation ---" << std::endl;
    std::cout << std::setprecision(8);
    std::cout << "Energy CPU:         " << f_cpu << std::endl;
    std::cout << "Energy GPU:         " << f_gpu << std::endl;
    std::cout << "Energy Abs Diff:    " << energy_abs_diff << std::endl;
    std::cout << "Energy Rel Diff:    " << energy_rel_diff << std::endl;
    std::cout << "Grad Max Abs Diff:  " << max_g_diff << std::endl;
    std::cout << "Grad Max Rel Diff:  " << g_rel_diff << std::endl;

    if (energy_rel_diff < 1e-4 && g_rel_diff < 1e-3) {
        std::cout << "[PASS] GPU matches CPU within single precision numerical tolerance!\n" << std::endl;
    } else {
        std::cout << "[WARNING] Discrepancy above standard threshold.\n" << std::endl;
    }

    // 3. Performance Benchmark
    std::cout << "--- Performance Benchmark (" << n_evals << " evaluations) ---" << std::endl;

    // Warm-up
    minimize_energy_with_triangles(free_dofs, f_cpu, g_cpu, &userData);
    minimize_energy_with_triangles_GPU(free_dofs, f_gpu, g_gpu, &userData);

    // CPU Timing
    auto t_start = clock_type::now();
    for (int i = 0; i < n_evals; ++i) {
        minimize_energy_with_triangles(free_dofs, f_cpu, g_cpu, &userData);
    }
    auto t_end = clock_type::now();
    double cpu_total_ms = elapsed_ms(t_start, t_end);
    double cpu_avg_ms = cpu_total_ms / n_evals;

    // GPU Timing
    t_start = clock_type::now();
    for (int i = 0; i < n_evals; ++i) {
        minimize_energy_with_triangles_GPU(free_dofs, f_gpu, g_gpu, &userData);
    }
    t_end = clock_type::now();
    double gpu_total_ms = elapsed_ms(t_start, t_end);
    double gpu_avg_ms = gpu_total_ms / n_evals;

    double speedup = cpu_avg_ms / gpu_avg_ms;

    std::cout << "CPU Average: " << std::fixed << std::setprecision(3) << cpu_avg_ms << " ms / eval ("
              << cpu_total_ms << " ms total)" << std::endl;
    std::cout << "GPU Average: " << std::fixed << std::setprecision(3) << gpu_avg_ms << " ms / eval ("
              << gpu_total_ms << " ms total)" << std::endl;
    std::cout << ">>> GPU SPEEDUP: " << std::fixed << std::setprecision(2) << speedup << "x <<<\n" << std::endl;

    // 4. ALGLIB L-BFGS Optimization Comparison
    if (run_lbfgs) {
        std::cout << "--- ALGLIB L-BFGS Optimization Test ---" << std::endl;

        // CPU run
        alglib::real_1d_array x_cpu = free_dofs;
        alglib::minlbfgsstate state_cpu;
        alglib::minlbfgsreport rep_cpu;
        alglib::minlbfgscreate(2 * n_vars, 10, x_cpu, state_cpu);
        alglib::minlbfgssetcond(state_cpu, 1e-4, 0.0, 0.0, 50);

        auto opt_cpu_start = clock_type::now();
        alglib::minlbfgsoptimize(state_cpu, minimize_energy_with_triangles, nullptr, &userData);
        auto opt_cpu_end = clock_type::now();
        alglib::minlbfgsresults(state_cpu, x_cpu, rep_cpu);
        double opt_cpu_ms = elapsed_ms(opt_cpu_start, opt_cpu_end);

        // GPU run
        alglib::real_1d_array x_gpu = free_dofs;
        alglib::minlbfgsstate state_gpu;
        alglib::minlbfgsreport rep_gpu;
        alglib::minlbfgscreate(2 * n_vars, 10, x_gpu, state_gpu);
        alglib::minlbfgssetcond(state_gpu, 1e-4, 0.0, 0.0, 50);

        auto opt_gpu_start = clock_type::now();
        alglib::minlbfgsoptimize(state_gpu, minimize_energy_with_triangles_GPU, nullptr, &userData);
        auto opt_gpu_end = clock_type::now();
        alglib::minlbfgsresults(state_gpu, x_gpu, rep_gpu);
        double opt_gpu_ms = elapsed_ms(opt_gpu_start, opt_gpu_end);

        double f_min_cpu = 0.0, f_min_gpu = 0.0;
        minimize_energy_with_triangles(x_cpu, f_min_cpu, g_cpu, &userData);
        minimize_energy_with_triangles(x_gpu, f_min_gpu, g_gpu, &userData);

        std::cout << "CPU Optimization:" << std::endl;
        std::cout << "  Iterations:   " << rep_cpu.iterationscount << std::endl;
        std::cout << "  Function calls: " << rep_cpu.nfev << std::endl;
        std::cout << "  Final Energy: " << std::setprecision(8) << f_min_cpu << std::endl;
        std::cout << "  Time:         " << std::setprecision(2) << opt_cpu_ms << " ms" << std::endl;

        std::cout << "\nGPU Optimization:" << std::endl;
        std::cout << "  Iterations:   " << rep_gpu.iterationscount << std::endl;
        std::cout << "  Function calls: " << rep_gpu.nfev << std::endl;
        std::cout << "  Final Energy: " << std::setprecision(8) << f_min_gpu << std::endl;
        std::cout << "  Time:         " << std::setprecision(2) << opt_gpu_ms << " ms" << std::endl;

        std::cout << "\nOptimization Speedup: " << std::fixed << std::setprecision(2)
                  << (opt_cpu_ms / opt_gpu_ms) << "x" << std::endl;
        std::cout << "Energy Discrepancy:   " << std::abs(f_min_cpu - f_min_gpu) << std::endl;
    }

    std::cout << "\n==========================================================" << std::endl;
    std::cout << "   Benchmark Complete!                                    " << std::endl;
    std::cout << "==========================================================" << std::endl;

    return 0;
}
