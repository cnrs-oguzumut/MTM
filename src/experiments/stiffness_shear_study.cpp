#include "../../include/experiments/stiffness_shear_study.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <memory>
#include <algorithm>

#include "../../include/interatomic/inter_atomic.h"
#include "../../include/geometry/LatticeGenerator.h"
#include "../../include/mesh/MeshGenerator.h"
#include "../../include/mesh/Remesher.h"
#include "../../include/mesh/ElementTriangle2D.h"
#include "../../include/optimization/LatticeOptimizer.h"
#include "../../include/optimization/PreconditionedLBFGS.h"
#include "../../include/optimization/StiffnessSpectrum.h"
#include "../../include/lattice_energy/SquareLatticeCalculator.h"
#include "../../include/lattice_energy/TriangularLatticeCalculator.h"
#include "../../include/lattice_energy/Strain_Energy_LatticeCalculator.h"
#include "../../include/lattice_energy/EnergyFunctions.h"
#include "../../include/geometry/DomainInfo.h"
#include "../../include/experiments/common_simulation_helpers.h"

namespace {

using Clock = std::chrono::high_resolution_clock;

double seconds_since(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

} // namespace

StiffnessShearResults runStiffnessShearStudy(const StiffnessShearConfig &config) {
    const auto t_start = Clock::now();
    StiffnessShearResults results;
    results.config = config;

    std::cout << "\n=======================================================\n";
    std::cout << "       STIFFNESS MATRIX K(α) SHEAR STUDY               \n";
    std::cout << "=======================================================\n";
    std::cout << "  Grid size:          " << config.nx << " x " << config.ny << "\n";
    std::cout << "  Lattice symmetry:   " << config.lattice_type << "\n";
    std::cout << "  Potential mode:     " << config.potential_mode << "\n";
    std::cout << "  Boundary condition: " << config.bc << "\n";
    std::cout << "  Shear range:        [" << config.alpha_min << ", " << config.alpha_max
              << "] step " << config.dalpha << "\n";
    std::cout << "-------------------------------------------------------\n";

    // 1. Potential setup & Optimal lattice parameter
    std::function<double(double)> pot_func;
    std::function<double(double)> pot_der;
    std::function<double(double)> pot_sder;
    double a_0 = 1.0;
    double r_cutoff = 2.5;

    if (config.potential_mode == 3) {
        energy_functions::set_parameters(config.beta, config.K);
        a_0 = (config.lattice_type == "triangular") ? std::pow(4.0 / 3.0, 0.25) : 1.0;
        pot_func = square_energy;
        pot_der = square_energy_der;
        pot_sder = square_energy_der;
    } else if (config.potential_mode == 1) {
        pot_func = lennard_jones_energy_v2;
        pot_der = lennard_jones_energy_der_v2;
        pot_sder = lennard_jones_energy_sder_v2;
        r_cutoff = 2.5;
        a_0 = find_optimal_lattice_parameter(pot_func, config.lattice_type, r_cutoff, pot_der, pot_sder);
    } else if (config.potential_mode == 2) {
        pot_func = lennard_jones_energy_v3;
        pot_der = lennard_jones_energy_der_v3;
        pot_sder = lennard_jones_energy_sder_v3;
        r_cutoff = 1.86602540378444;
        a_0 = find_optimal_lattice_parameter(pot_func, config.lattice_type, r_cutoff, pot_der, pot_sder);
    } else if (config.potential_mode == 4) {
        pot_func = square_energy;
        pot_der = square_energy_der;
        pot_sder = square_energy_sder;
        r_cutoff = 2.5;
        a_0 = find_optimal_lattice_parameter(pot_func, config.lattice_type, r_cutoff, pot_der, pot_sder);
    } else {
        std::cerr << "Unknown potential mode: " << config.potential_mode << std::endl;
        return results;
    }

    std::cout << "  Optimal lattice parameter a_0 = " << std::fixed << std::setprecision(8) << a_0 << "\n";

    // 2. Calculator instantiation
    std::unique_ptr<BaseLatticeCalculator> calc;
    double normalisation = 1.0;
    Eigen::Matrix2d C_ref;

    double calc_scale = a_0;
    if (config.potential_mode == 3) {
        calc = std::make_unique<Strain_Energy_LatticeCalculator>(1.0);
        normalisation = calc->getUnitCellArea();
        if (config.lattice_type == "triangular") {
            double a2 = std::sqrt(4.0 / 3.0);
            C_ref << a2, 0.5 * a2,
                     0.5 * a2, a2;
            calc_scale = 1.0;
        } else {
            C_ref = Eigen::Matrix2d::Identity();
            calc_scale = 1.0;
        }
    } else {
        if (config.lattice_type == "triangular") {
            // NOTE: Do NOT use TriangularLatticeCalculator here.
            // Lagrange reduction (GL(2, Z)) operates on the integer reference grid Z^2.
            // Using SquareLatticeCalculator with reference metric C_ref = F0^T * F0 avoids
            // double-applying the triangular distortion. (TriangularLatticeCalculator is
            // reserved for direct physical-coordinate simulations such as nanoindentation).
            calc_scale = a_0 / std::pow(4.0 / 3.0, 0.25);
            calc = std::make_unique<SquareLatticeCalculator>(calc_scale, r_cutoff);
            normalisation = calc->getUnitCellArea() * (std::sqrt(3.0) / 2.0);
            double a2 = std::sqrt(4.0 / 3.0);
            C_ref << a2, 0.5 * a2,
                     0.5 * a2, a2;
        } else {
            calc_scale = a_0;
            calc = std::make_unique<SquareLatticeCalculator>(a_0, r_cutoff);
            normalisation = calc->getUnitCellArea();
            C_ref = Eigen::Matrix2d::Identity();
        }
    }
    double zero_energy = calc->calculate_energy(C_ref, pot_func, 0.0);

    int study_ny = config.ny;
    const bool pbc = (config.bc == "periodic" || config.bc == "Periodic" || config.bc == "pbc" || config.bc == "PBC");
    if (config.lattice_type == "triangular" && (study_ny % 2 != 0)) {
        study_ny++;
        std::cout << "  [INFO] Enforcing even ny for triangular lattice: adjusted from "
                  << config.ny << " to " << study_ny << "\n";
    }

    // 3. Pristine Reference Mesh Generation
    std::vector<Point2D> points_ref = LatticeGenerator::generate_2d_lattice(
        config.nx, study_ny, a_0, config.lattice_type);
    const int total_nodes = static_cast<int>(points_ref.size());

    DomainInfo domain_size = compute_domain_size(points_ref);
    DomainDimensions domain_dims(domain_size.get_width(), domain_size.get_height());
    std::array<double, 2> offsets = (config.lattice_type == "square")
        ? std::array<double, 2>{a_0, a_0}
        : std::array<double, 2>{a_0 / 2.0, (std::sqrt(3.0) / 2.0) * a_0};

    auto [original_domain_map, translation_map] =
        MeshGenerator::create_domain_maps(total_nodes, domain_dims, offsets);

    auto [interior_mapping, full_mapping] =
        create_dof_mapping_original(points_ref, 0.25 * a_0, pbc);

    const int free_nodes = static_cast<int>(interior_mapping.size());
    const int total_dofs = 2 * free_nodes;

    Point2D domain_dims_point(domain_dims.size_x, domain_dims.size_y);
    AdaptiveMesher mesher(domain_dims_point, offsets, original_domain_map,
                          translation_map, full_mapping, 1e-6, pbc);
    mesher.setUsePeriodicCopies(pbc);

    alglib::real_1d_array x0;
    x0.setlength(total_dofs);
    map_points_to_solver_array(x0, points_ref, interior_mapping, free_nodes);

    // Reference element shape function derivatives:
    // Scale h corresponds to calc_scale so that F_elem = F_shear * F0 at all scales
    double h = (config.potential_mode == 3) ? 1.0 : calc_scale;
    Eigen::Vector2d p1(0.0, 0.0);
    Eigen::Vector2d p2(h, 0.0);
    Eigen::Vector2d p3(0.0, h);
    Eigen::Matrix<double, 3, 2> dndx = calculateShapeDerivatives(p1, p2, p3);

    auto [elements, active_elements] =
        mesher.createMesh(points_ref, x0, Eigen::Matrix2d::Identity(), &dndx);

    for (auto &elem : elements) {
        elem.set_dof_mapping(full_mapping);
        elem.set_reference_mesh(points_ref);
    }

    const int total_elements = static_cast<int>(elements.size());
    results.summary.total_nodes = total_nodes;
    results.summary.free_nodes = free_nodes;
    results.summary.total_dofs = total_dofs;
    results.summary.total_elements = total_elements;

    std::cout << "  Mesh info:          " << total_nodes << " nodes, "
              << free_nodes << " free nodes (" << total_dofs << " DOFs), "
              << total_elements << " elements (" << active_elements.size() << " active)\n";
    std::cout << "-------------------------------------------------------\n";
    std::cout << std::setw(10) << "alpha"
              << std::setw(16) << "lambda_min"
              << std::setw(16) << "lambda_2"
              << std::setw(16) << "lambda_3"
              << std::setw(14) << "Energy"
              << std::setw(10) << "Status" << "\n";
    std::cout << "-------------------------------------------------------\n";

    // 4. Simple Shear Deformation Sweep
    if (config.dalpha <= 1e-9) {
        std::cerr << "  [ERROR] Step dalpha must be positive and non-zero (got " << config.dalpha << ")\n";
        return results;
    }
    FastStiffnessAssembler assembler;
    int num_steps = std::max(1, static_cast<int>(std::round((config.alpha_max - config.alpha_min) / config.dalpha))) + 1;

    double global_min_lambda = 1e30;
    double global_max_lambda = -1e30;
    double last_lambda = 0.0;
    bool has_last = false;
    double alpha_c_neg = 0.0;
    double alpha_c_pos = 0.0;
    bool found_neg_c = false;
    bool found_pos_c = false;

    for (int step = 0; step < num_steps; ++step) {
        double alpha = config.alpha_min + step * config.dalpha;
        if (std::abs(alpha) < 1e-12) alpha = 0.0;

        Eigen::Matrix2d F_shear;
        F_shear << 1.0, alpha,
                   0.0, 1.0;

        // Apply homogeneous affine shear to node positions
        std::vector<Point2D> current_points = points_ref;
        for (size_t i = 0; i < current_points.size(); ++i) {
            current_points[i].coord = F_shear * points_ref[i].coord;
        }

        alglib::real_1d_array x_alpha;
        x_alpha.setlength(total_dofs);
        map_points_to_solver_array(x_alpha, current_points, interior_mapping, free_nodes);

        // Update elements
        for (auto &elem : elements) {
            elem.setExternalDeformation(F_shear);
            elem.set_reference_mesh(current_points);
            elem.calculate_deformation_gradient(x_alpha);
        }

        // Assemble tangent stiffness matrix K
        UserData userData(current_points, elements, *calc, pot_func, pot_der, zero_energy,
                          a_0, F_shear, interior_mapping, full_mapping, active_elements, false);
        if (config.potential_mode != 3) {
            userData.second_derivative_function = &pot_sder;
        }

        const Eigen::SparseMatrix<double> &K = assembler.assembleStiffness(x_alpha, &userData, false);

        // Compute smallest eigenvalues of K
        double lambda_min = 0.0;
        double lambda_2 = 0.0;
        double lambda_3 = 0.0;

        if (total_dofs <= 1200) {
            Eigen::MatrixXd K_dense(K);
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(K_dense);
            const auto &evals = es.eigenvalues();
            const auto &evecs = es.eigenvectors();

            std::vector<double> non_rigid_evals;
            non_rigid_evals.reserve(evals.size());

            for (int k = 0; k < evals.size(); ++k) {
                if (pbc) {
                    const int m = static_cast<int>(evecs.rows() / 2);
                    const double mu = evecs.col(k).head(m).mean();
                    const double mv = evecs.col(k).tail(m).mean();
                    const double trans_frac = (mu * mu + mv * mv) * static_cast<double>(m) / std::max(evecs.col(k).squaredNorm(), 1e-300);
                    if (trans_frac > 0.5) {
                        continue; // Skip rigid translation mode
                    }
                }
                non_rigid_evals.push_back(evals(k));
            }

            if (!non_rigid_evals.empty()) {
                lambda_min = non_rigid_evals[0];
                lambda_2 = (non_rigid_evals.size() > 1) ? non_rigid_evals[1] : lambda_min;
                lambda_3 = (non_rigid_evals.size() > 2) ? non_rigid_evals[2] : lambda_2;
            }
        } else {
            StiffnessSpectrumOptions opt;
            opt.n_modes = 5;
            opt.deflate_translations = pbc ? 1 : 0;
            StiffnessSpectrum spec = lowest_stiffness_modes(K, opt);
            if (spec.num_computed > 0) {
                lambda_min = spec.eigenvalues(0);
                lambda_2 = (spec.num_computed > 1) ? spec.eigenvalues(1) : lambda_min;
                lambda_3 = (spec.num_computed > 2) ? spec.eigenvalues(2) : lambda_2;
            }
        }

        // Calculate total strain energy & normalize by finite element triangle number
        double total_energy = 0.0;
        for (size_t idx : active_elements) {
            const ElementTriangle2D &elem = elements[idx];
            const auto red = lagrange::reduce(elem.getMetricTensor());
            double elem_e = calc->calculate_energy(red.C_reduced, pot_func, zero_energy) / normalisation;
            total_energy += elem_e;
        }
        const double num_active = std::max(1.0, static_cast<double>(active_elements.size()));
        double energy_density = total_energy / num_active;

        // Continuum strain energy density
        Eigen::Matrix2d C_cont = F_shear.transpose() * F_shear;
        if (config.lattice_type == "triangular") {
            Eigen::Matrix2d F0;
            double a = std::pow(4.0 / 3.0, 0.25);
            F0 << a, a * 0.5,
                  0.0, a * std::sqrt(3.0) / 2.0;
            Eigen::Matrix2d F_tot = F_shear * F0;
            C_cont = F_tot.transpose() * F_tot;
        }
        const auto red_cont = lagrange::reduce(C_cont);
        double continuum_energy = calc->calculate_energy(red_cont.C_reduced, pot_func, zero_energy) / normalisation;

        bool is_stable = (lambda_min > 0.0);
        global_min_lambda = std::min(global_min_lambda, lambda_min);
        global_max_lambda = std::max(global_max_lambda, lambda_min);

        // Detect critical alpha (zero crossing)
        if (has_last) {
            if ((last_lambda > 0.0 && lambda_min <= 0.0) || (last_lambda < 0.0 && lambda_min >= 0.0)) {
                double frac = std::abs(last_lambda) / (std::abs(last_lambda) + std::abs(lambda_min) + 1e-15);
                double alpha_cross = (alpha - config.dalpha) + frac * config.dalpha;
                if (alpha_cross < 0.0 && !found_neg_c) {
                    alpha_c_neg = alpha_cross;
                    found_neg_c = true;
                } else if (alpha_cross > 0.0 && !found_pos_c) {
                    alpha_c_pos = alpha_cross;
                    found_pos_c = true;
                }
            }
        }
        last_lambda = lambda_min;
        has_last = true;

        StiffnessShearPoint pt;
        pt.alpha = alpha;
        pt.lambda_min = lambda_min;
        pt.lambda_2 = lambda_2;
        pt.lambda_3 = lambda_3;
        pt.is_stable = is_stable;
        pt.energy = total_energy;
        pt.energy_density = energy_density;
        pt.continuum_energy = continuum_energy;
        results.points.push_back(pt);

        // Print log every few steps or near transition
        if (step % 5 == 0 || !is_stable || (has_last && last_lambda * lambda_min <= 0.0)) {
            std::cout << std::fixed << std::setprecision(4) << std::setw(10) << alpha
                      << std::scientific << std::setprecision(5)
                      << std::setw(16) << lambda_min
                      << std::setw(16) << lambda_2
                      << std::setw(16) << lambda_3
                      << std::setw(14) << total_energy
                      << std::setw(10) << (is_stable ? "STABLE" : "UNSTABLE") << "\n";
        }
    }

    results.summary.alpha_c_neg = alpha_c_neg;
    results.summary.alpha_c_pos = alpha_c_pos;
    results.summary.min_lambda = global_min_lambda;
    results.summary.max_lambda = global_max_lambda;
    results.summary.execution_seconds = seconds_since(t_start);

    std::cout << "-------------------------------------------------------\n";
    std::cout << "  Summary:\n";
    if (found_neg_c) std::cout << "    Negative critical shear α_c(-) = " << std::fixed << std::setprecision(5) << alpha_c_neg << "\n";
    if (found_pos_c) std::cout << "    Positive critical shear α_c(+) = " << std::fixed << std::setprecision(5) << alpha_c_pos << "\n";
    std::cout << "    Global min λ_min(K) = " << std::scientific << std::setprecision(5) << global_min_lambda << "\n";
    std::cout << "    Execution time:       " << std::fixed << std::setprecision(3) << results.summary.execution_seconds << " s\n";
    std::cout << "=======================================================\n\n";

    // 5. Output CSV
    if (!config.csv_output.empty()) {
        std::ofstream csv(config.csv_output);
        if (csv.is_open()) {
            csv << "alpha,lambda_min,lambda_2,lambda_3,is_stable,energy,energy_density,continuum_energy\n";
            csv << std::scientific << std::setprecision(10);
            for (const auto &p : results.points) {
                csv << p.alpha << ","
                    << p.lambda_min << ","
                    << p.lambda_2 << ","
                    << p.lambda_3 << ","
                    << (p.is_stable ? 1 : 0) << ","
                    << p.energy << ","
                    << p.energy_density << ","
                    << p.continuum_energy << "\n";
            }
            csv.close();
            std::cout << "Wrote CSV results to " << config.csv_output << std::endl;
        }
    }

    // 6. Output JSON
    if (!config.json_output.empty()) {
        std::ofstream json(config.json_output);
        if (json.is_open()) {
            json << "{\n";
            json << "  \"config\": {\n";
            json << "    \"nx\": " << config.nx << ",\n";
            json << "    \"ny\": " << config.ny << ",\n";
            json << "    \"lattice_type\": \"" << config.lattice_type << "\",\n";
            json << "    \"potential_mode\": " << config.potential_mode << ",\n";
            json << "    \"bc\": \"" << config.bc << "\",\n";
            json << "    \"alpha_min\": " << config.alpha_min << ",\n";
            json << "    \"alpha_max\": " << config.alpha_max << ",\n";
            json << "    \"dalpha\": " << config.dalpha << ",\n";
            json << "    \"beta\": " << config.beta << ",\n";
            json << "    \"K\": " << config.K << "\n";
            json << "  },\n";
            json << "  \"summary\": {\n";
            json << "    \"total_nodes\": " << results.summary.total_nodes << ",\n";
            json << "    \"free_nodes\": " << results.summary.free_nodes << ",\n";
            json << "    \"total_dofs\": " << results.summary.total_dofs << ",\n";
            json << "    \"total_elements\": " << results.summary.total_elements << ",\n";
            json << "    \"alpha_c_pos\": " << results.summary.alpha_c_pos << ",\n";
            json << "    \"alpha_c_neg\": " << results.summary.alpha_c_neg << ",\n";
            json << "    \"min_lambda\": " << results.summary.min_lambda << ",\n";
            json << "    \"max_lambda\": " << results.summary.max_lambda << ",\n";
            json << "    \"execution_seconds\": " << results.summary.execution_seconds << "\n";
            json << "  },\n";
            json << "  \"points\": [\n";
            for (size_t i = 0; i < results.points.size(); ++i) {
                const auto &p = results.points[i];
                json << "    {\"alpha\": " << p.alpha
                     << ", \"lambda_min\": " << p.lambda_min
                     << ", \"lambda_2\": " << p.lambda_2
                     << ", \"lambda_3\": " << p.lambda_3
                     << ", \"is_stable\": " << (p.is_stable ? "true" : "false")
                     << ", \"energy\": " << p.energy
                     << ", \"energy_density\": " << p.energy_density
                     << ", \"continuum_energy\": " << p.continuum_energy << "}";
                if (i + 1 < results.points.size()) json << ",";
                json << "\n";
            }
            json << "  ]\n";
            json << "}\n";
            json.close();
            std::cout << "Wrote JSON results to " << config.json_output << std::endl;
        }
    }

    return results;
}
