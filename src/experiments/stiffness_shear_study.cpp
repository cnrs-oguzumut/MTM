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
    int num_steps = config.single_alpha_only
        ? 1
        : (std::max(1, static_cast<int>(std::round((config.alpha_max - config.alpha_min) / config.dalpha))) + 1);

    double global_min_lambda = 1e30;
    double global_max_lambda = -1e30;
    double last_lambda = 0.0;
    bool has_last = false;
    double alpha_c_neg = 0.0;
    double alpha_c_pos = 0.0;
    bool found_neg_c = false;
    bool found_pos_c = false;

    for (int step = 0; step < num_steps; ++step) {
        double alpha = config.single_alpha_only ? config.single_alpha : (config.alpha_min + step * config.dalpha);
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

        // Compute smallest eigenvalues & eigenvectors of K
        double lambda_min = 0.0;
        double lambda_2 = 0.0;
        double lambda_3 = 0.0;
        double lambda_4 = 0.0;
        double lambda_5 = 0.0;
        std::vector<Eigen::VectorXd> mode_vecs;

        std::vector<double> cand_evals;
        std::vector<Eigen::VectorXd> cand_vecs;

        if (total_dofs <= 1200) {
            Eigen::MatrixXd K_dense(K);
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(K_dense);
            const auto &evals = es.eigenvalues();
            const auto &evecs = es.eigenvectors();

            const int m = static_cast<int>(evecs.rows() / 2);
            for (int k = 0; k < evals.size(); ++k) {
                Eigen::VectorXd vec = evecs.col(k);
                if (pbc && m > 0) {
                    const double mu = vec.head(m).mean();
                    const double mv = vec.tail(m).mean();
                    const double trans_frac = (mu * mu + mv * mv) * static_cast<double>(m) / std::max(vec.squaredNorm(), 1e-300);
                    if (trans_frac > 0.5) {
                        continue; // Skip rigid translation mode
                    }
                    // Project out any residual uniform translation component
                    vec.head(m).array() -= mu;
                    vec.tail(m).array() -= mv;
                }
                cand_evals.push_back(evals(k));
                cand_vecs.push_back(vec);
            }
        } else {
            StiffnessSpectrumOptions opt;
            // When filtering unique physical levels (skipping sin/cos doublets), request more raw modes
            opt.n_modes = config.unique_levels_only ? 14 : 6;
            opt.deflate_translations = pbc ? 1 : 0;
            opt.compute_vectors = true;
            StiffnessSpectrum spec = lowest_stiffness_modes(K, opt);
            if (spec.num_computed > 0) {
                for (int k = 0; k < spec.num_computed; ++k) {
                    cand_evals.push_back(spec.eigenvalues(k));
                    if (k < spec.eigenvectors.cols()) {
                        cand_vecs.push_back(spec.eigenvectors.col(k));
                    }
                }
            }
        }

        // Filter and select up to 5 modes (optionally skipping degenerate sin/cos twins)
        std::vector<double> selected_evals;
        std::vector<Eigen::VectorXd> selected_vecs;
        for (size_t k = 0; k < cand_evals.size() && selected_evals.size() < 5; ++k) {
            double val = cand_evals[k];
            if (config.unique_levels_only && !selected_evals.empty()) {
                double last_val = selected_evals.back();
                double rel_diff = std::abs(val - last_val) / std::max(1e-12, std::abs(last_val));
                if (rel_diff < 1e-4) {
                    continue; // Skip sin/cos degenerate partner of previous level
                }
            }
            selected_evals.push_back(val);
            if (k < cand_vecs.size()) {
                selected_vecs.push_back(cand_vecs[k]);
            }
        }

        if (!selected_evals.empty()) {
            lambda_min = selected_evals[0];
            lambda_2 = (selected_evals.size() > 1) ? selected_evals[1] : lambda_min;
            lambda_3 = (selected_evals.size() > 2) ? selected_evals[2] : lambda_2;
            lambda_4 = (selected_evals.size() > 3) ? selected_evals[3] : lambda_3;
            lambda_5 = (selected_evals.size() > 4) ? selected_evals[4] : lambda_4;
        }
        mode_vecs = selected_vecs;

        // Map solver DOF eigenvectors to full N_nodes lattice and normalize so max_i |du_i| = 1
        auto extractFullNodeMode = [&](const Eigen::VectorXd &v,
                                       std::vector<double> &out_ux,
                                       std::vector<double> &out_uy) {
            out_ux.assign(total_nodes, 0.0);
            out_uy.assign(total_nodes, 0.0);
            if (v.size() < 2 * free_nodes) return;

            double max_amp = 0.0;
            int max_idx = 0;
            for (int i = 0; i < total_nodes; ++i) {
                int s_idx = (i < static_cast<int>(full_mapping.size())) ? full_mapping[i].second : -1;
                if (s_idx >= 0 && s_idx < free_nodes) {
                    double ux = v(s_idx);
                    double uy = v(s_idx + free_nodes);
                    out_ux[i] = ux;
                    out_uy[i] = uy;
                    double amp = std::hypot(ux, uy);
                    if (amp > max_amp) {
                        max_amp = amp;
                        max_idx = i;
                    }
                }
            }
            if (max_amp > 1e-15) {
                // Consistent sign orientation so the dominant node displacement points in +x/+y half-plane
                double sgn = (out_ux[max_idx] + 0.37 * out_uy[max_idx] >= 0.0) ? 1.0 : -1.0;
                double inv = sgn / max_amp;
                for (int i = 0; i < total_nodes; ++i) {
                    out_ux[i] *= inv;
                    out_uy[i] *= inv;
                }
            }
        };

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
        pt.lambda_4 = lambda_4;
        pt.lambda_5 = lambda_5;
        pt.is_stable = is_stable;
        pt.energy = total_energy;
        pt.energy_density = energy_density;
        pt.continuum_energy = continuum_energy;
        if (mode_vecs.size() > 0) extractFullNodeMode(mode_vecs[0], pt.mode1_ux, pt.mode1_uy);
        if (mode_vecs.size() > 1) extractFullNodeMode(mode_vecs[1], pt.mode2_ux, pt.mode2_uy);
        if (mode_vecs.size() > 2) extractFullNodeMode(mode_vecs[2], pt.mode3_ux, pt.mode3_uy);
        if (mode_vecs.size() > 3) extractFullNodeMode(mode_vecs[3], pt.mode4_ux, pt.mode4_uy);
        if (mode_vecs.size() > 4) extractFullNodeMode(mode_vecs[4], pt.mode5_ux, pt.mode5_uy);
        results.points.push_back(pt);

        // Print log every few steps or near transition, and flush progress for GUI
        if (step % 5 == 0 || !is_stable || (has_last && last_lambda * lambda_min <= 0.0)) {
            std::cout << std::fixed << std::setprecision(4) << std::setw(10) << alpha
                      << std::scientific << std::setprecision(5)
                      << std::setw(15) << lambda_min
                      << std::setw(15) << lambda_2
                      << std::setw(15) << lambda_3
                      << std::setw(15) << lambda_4
                      << std::setw(15) << lambda_5
                      << std::setw(14) << total_energy
                      << std::setw(10) << (is_stable ? "STABLE" : "UNSTABLE") << std::endl;
        }
        std::cout << "[PROGRESS] " << (step + 1) << "/" << num_steps
                  << " alpha=" << std::fixed << std::setprecision(4) << alpha
                  << " elapsed=" << std::fixed << std::setprecision(2) << seconds_since(t_start)
                  << " pt=" << std::scientific << std::setprecision(8)
                  << alpha << "," << lambda_min << "," << lambda_2 << "," << lambda_3 << ","
                  << lambda_4 << "," << lambda_5 << ","
                  << (is_stable ? 1 : 0) << "," << total_energy << "," << energy_density << "," << continuum_energy
                  << std::endl;
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

    // 5. Output CSV (skip if single-alpha full-scale evaluation)
    if (!config.single_alpha_only && !config.csv_output.empty()) {
        std::ofstream csv(config.csv_output);
        if (csv.is_open()) {
            csv << "alpha,lambda_min,lambda_2,lambda_3,lambda_4,lambda_5,is_stable,energy,energy_density,continuum_energy\n";
            csv << std::scientific << std::setprecision(10);
            for (const auto &p : results.points) {
                csv << p.alpha << ","
                    << p.lambda_min << ","
                    << p.lambda_2 << ","
                    << p.lambda_3 << ","
                    << p.lambda_4 << ","
                    << p.lambda_5 << ","
                    << (p.is_stable ? 1 : 0) << ","
                    << p.energy << ","
                    << p.energy_density << ","
                    << p.continuum_energy << "\n";
            }
            csv.close();
            std::cout << "Wrote CSV results to " << config.csv_output << std::endl;
        }
    }

    // 6. Output JSON (skip if single-alpha full-scale evaluation)
    if (!config.single_alpha_only && !config.json_output.empty()) {
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
            json << "    \"K\": " << config.K << ",\n";
            json << "    \"unique_levels_only\": " << (config.unique_levels_only ? "true" : "false") << "\n";
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

            // Build visualization node subset (full mesh if <= 625 nodes; stride-subsampled <= 24x24 grid if larger)
            std::vector<int> vis_nodes;
            std::vector<std::array<int, 3>> vis_tris;
            if (total_nodes <= 625) {
                vis_nodes.resize(total_nodes);
                for (int i = 0; i < total_nodes; ++i) vis_nodes[i] = i;
                for (size_t idx : active_elements) {
                    const ElementTriangle2D &elem = elements[idx];
                    double trans_norm = elem.getTranslation(0).norm() + elem.getTranslation(1).norm() + elem.getTranslation(2).norm();
                    if (trans_norm > 1e-5) continue;
                    int n0 = elem.getNodeIndex(0), n1 = elem.getNodeIndex(1), n2 = elem.getNodeIndex(2);
                    if (n0 >= 0 && n0 < total_nodes && n1 >= 0 && n1 < total_nodes && n2 >= 0 && n2 < total_nodes) {
                        vis_tris.push_back({n0, n1, n2});
                    }
                }
            } else {
                const int sx = std::max(1, config.nx / 22);
                const int sy = std::max(1, study_ny / 22);
                std::vector<int> grid_row_cols;
                int sub_nx = 0, sub_ny = 0;
                for (int j = 0; j < study_ny; j += sy) {
                    int row_cnt = 0;
                    for (int i = 0; i < config.nx; i += sx) {
                        int orig_idx = j * config.nx + i;
                        if (orig_idx < total_nodes) {
                            vis_nodes.push_back(orig_idx);
                            row_cnt++;
                        }
                    }
                    if (sub_nx == 0) sub_nx = row_cnt;
                    sub_ny++;
                }
                for (int j = 0; j + 1 < sub_ny; ++j) {
                    for (int i = 0; i + 1 < sub_nx; ++i) {
                        int n00 = j * sub_nx + i;
                        int n10 = j * sub_nx + (i + 1);
                        int n01 = (j + 1) * sub_nx + i;
                        int n11 = (j + 1) * sub_nx + (i + 1);
                        vis_tris.push_back({n00, n10, n11});
                        vis_tris.push_back({n00, n11, n01});
                    }
                }
            }
            const size_t n_vis = vis_nodes.size();

            // Write reference mesh geometry once for eigenvector mode visualization
            json << "  \"mesh\": {\n";
            json << "    \"a_0\": " << a_0 << ",\n";
            json << "    \"nodes_x\": [";
            for (size_t k = 0; k < n_vis; ++k) {
                json << std::fixed << std::setprecision(5) << points_ref[vis_nodes[k]].coord(0) << (k + 1 < n_vis ? "," : "");
            }
            json << "],\n";
            json << "    \"nodes_y\": [";
            for (size_t k = 0; k < n_vis; ++k) {
                json << std::fixed << std::setprecision(5) << points_ref[vis_nodes[k]].coord(1) << (k + 1 < n_vis ? "," : "");
            }
            json << "],\n";
            json << "    \"is_free\": [";
            for (size_t k = 0; k < n_vis; ++k) {
                int i = vis_nodes[k];
                int s_idx = (i < static_cast<int>(full_mapping.size())) ? full_mapping[i].second : -1;
                json << ((s_idx >= 0 && s_idx < free_nodes) ? 1 : 0) << (k + 1 < n_vis ? "," : "");
            }
            json << "],\n";
            json << "    \"triangles\": [";
            for (size_t k = 0; k < vis_tris.size(); ++k) {
                json << "[" << vis_tris[k][0] << "," << vis_tris[k][1] << "," << vis_tris[k][2] << "]"
                     << (k + 1 < vis_tris.size() ? "," : "");
            }
            json << "]\n";
            json << "  },\n";

            auto writeVecCompact = [&](const char *key, const std::vector<double> &vec) {
                if (vec.empty()) return;
                json << ", \"" << key << "\": [";
                json << std::fixed << std::setprecision(4);
                for (size_t k = 0; k < n_vis; ++k) {
                    int idx = vis_nodes[k];
                    double val = (idx >= 0 && idx < static_cast<int>(vec.size())) ? vec[idx] : 0.0;
                    json << val << (k + 1 < n_vis ? "," : "");
                }
                json << "]";
                json << std::defaultfloat << std::setprecision(6);
            };

            json << "  \"points\": [\n";
            for (size_t i = 0; i < results.points.size(); ++i) {
                const auto &p = results.points[i];
                json << std::defaultfloat << std::setprecision(6);
                json << "    {\"alpha\": " << p.alpha
                     << ", \"lambda_min\": " << p.lambda_min
                     << ", \"lambda_2\": " << p.lambda_2
                     << ", \"lambda_3\": " << p.lambda_3
                     << ", \"lambda_4\": " << p.lambda_4
                     << ", \"lambda_5\": " << p.lambda_5
                     << ", \"is_stable\": " << (p.is_stable ? "true" : "false")
                     << ", \"energy\": " << p.energy
                     << ", \"energy_density\": " << p.energy_density
                     << ", \"continuum_energy\": " << p.continuum_energy;
                writeVecCompact("mode1_ux", p.mode1_ux);
                writeVecCompact("mode1_uy", p.mode1_uy);
                writeVecCompact("mode2_ux", p.mode2_ux);
                writeVecCompact("mode2_uy", p.mode2_uy);
                writeVecCompact("mode3_ux", p.mode3_ux);
                writeVecCompact("mode3_uy", p.mode3_uy);
                writeVecCompact("mode4_ux", p.mode4_ux);
                writeVecCompact("mode4_uy", p.mode4_uy);
                writeVecCompact("mode5_ux", p.mode5_ux);
                writeVecCompact("mode5_uy", p.mode5_uy);
                json << "}";
                if (i + 1 < results.points.size()) json << ",";
                json << "\n";
            }
            json << "  ]\n";
            json << "}\n";
            json.close();
            std::cout << "Wrote JSON results to " << config.json_output << std::endl;
        }
    }

    // 7. Output Full-Scale Static Snapshot JSON (all N_x * N_y nodes at critical alphas or single_alpha)
    if (!config.fullscale_json_output.empty() && !results.points.empty()) {
        std::vector<size_t> snap_indices;
        if (config.single_alpha_only || results.points.size() <= 3) {
            for (size_t i = 0; i < results.points.size(); ++i) snap_indices.push_back(i);
        } else {
            // Find closest steps to alpha_c(-), alpha = 0, and alpha_c(+)
            auto addClosest = [&](double target_a) {
                size_t best_i = 0;
                double best_d = 1e30;
                for (size_t i = 0; i < results.points.size(); ++i) {
                    double d = std::abs(results.points[i].alpha - target_a);
                    if (d < best_d) { best_d = d; best_i = i; }
                }
                if (std::find(snap_indices.begin(), snap_indices.end(), best_i) == snap_indices.end()) {
                    snap_indices.push_back(best_i);
                }
            };
            if (found_neg_c) addClosest(alpha_c_neg);
            addClosest(0.0);
            if (found_pos_c) addClosest(alpha_c_pos);
            // Also add global min lambda step
            size_t min_l_idx = 0;
            for (size_t i = 1; i < results.points.size(); ++i) {
                if (results.points[i].lambda_min < results.points[min_l_idx].lambda_min) min_l_idx = i;
            }
            if (std::find(snap_indices.begin(), snap_indices.end(), min_l_idx) == snap_indices.end()) {
                snap_indices.push_back(min_l_idx);
            }
        }

        std::ofstream fs_json(config.fullscale_json_output);
        if (fs_json.is_open()) {
            fs_json << "{\n";
            fs_json << "  \"nx\": " << config.nx << ",\n";
            fs_json << "  \"ny\": " << study_ny << ",\n";
            fs_json << "  \"a_0\": " << a_0 << ",\n";
            fs_json << "  \"lattice_type\": \"" << config.lattice_type << "\",\n";
            fs_json << "  \"bc\": \"" << config.bc << "\",\n";
            fs_json << "  \"unique_levels_only\": " << (config.unique_levels_only ? "true" : "false") << ",\n";

            // Write exact reference mesh used in calculation (zero UI re-meshing)
            fs_json << "  \"mesh\": {\n";
            fs_json << "    \"a_0\": " << a_0 << ",\n";
            fs_json << "    \"nodes_x\": [";
            for (int k = 0; k < total_nodes; ++k) {
                fs_json << std::fixed << std::setprecision(5) << points_ref[k].coord(0) << (k + 1 < total_nodes ? "," : "");
            }
            fs_json << "],\n";
            fs_json << "    \"nodes_y\": [";
            for (int k = 0; k < total_nodes; ++k) {
                fs_json << std::fixed << std::setprecision(5) << points_ref[k].coord(1) << (k + 1 < total_nodes ? "," : "");
            }
            fs_json << "],\n";
            fs_json << "    \"is_free\": [";
            for (int k = 0; k < total_nodes; ++k) {
                int s_idx = (k < static_cast<int>(full_mapping.size())) ? full_mapping[k].second : -1;
                fs_json << ((s_idx >= 0 && s_idx < free_nodes) ? 1 : 0) << (k + 1 < total_nodes ? "," : "");
            }
            fs_json << "],\n";
            fs_json << "    \"triangles\": [";
            bool first_tri = true;
            for (size_t idx : active_elements) {
                const ElementTriangle2D &elem = elements[idx];
                double trans_norm = elem.getTranslation(0).norm() + elem.getTranslation(1).norm() + elem.getTranslation(2).norm();
                if (trans_norm > 1e-5) continue;
                int n0 = elem.getNodeIndex(0), n1 = elem.getNodeIndex(1), n2 = elem.getNodeIndex(2);
                if (n0 >= 0 && n0 < total_nodes && n1 >= 0 && n1 < total_nodes && n2 >= 0 && n2 < total_nodes) {
                    if (!first_tri) fs_json << ",";
                    fs_json << "[" << n0 << "," << n1 << "," << n2 << "]";
                    first_tri = false;
                }
            }
            fs_json << "]\n";
            fs_json << "  },\n";

            fs_json << "  \"snapshots\": [\n";

            auto writeFullVec = [&](const char *key, const std::vector<double> &vec) {
                if (vec.empty()) return;
                fs_json << ", \"" << key << "\": [";
                fs_json << std::fixed << std::setprecision(4);
                for (size_t k = 0; k < vec.size(); ++k) {
                    fs_json << vec[k] << (k + 1 < vec.size() ? "," : "");
                }
                fs_json << "]";
                fs_json << std::defaultfloat << std::setprecision(6);
            };

            for (size_t s = 0; s < snap_indices.size(); ++s) {
                const auto &p = results.points[snap_indices[s]];
                fs_json << std::defaultfloat << std::setprecision(6);
                fs_json << "    {\"alpha\": " << p.alpha
                        << ", \"lambda_min\": " << p.lambda_min
                        << ", \"lambda_2\": " << p.lambda_2
                        << ", \"lambda_3\": " << p.lambda_3
                        << ", \"lambda_4\": " << p.lambda_4
                        << ", \"lambda_5\": " << p.lambda_5
                        << ", \"is_stable\": " << (p.is_stable ? "true" : "false");
                writeFullVec("mode1_ux", p.mode1_ux);
                writeFullVec("mode1_uy", p.mode1_uy);
                writeFullVec("mode2_ux", p.mode2_ux);
                writeFullVec("mode2_uy", p.mode2_uy);
                writeFullVec("mode3_ux", p.mode3_ux);
                writeFullVec("mode3_uy", p.mode3_uy);
                writeFullVec("mode4_ux", p.mode4_ux);
                writeFullVec("mode4_uy", p.mode4_uy);
                writeFullVec("mode5_ux", p.mode5_ux);
                writeFullVec("mode5_uy", p.mode5_uy);
                fs_json << "}" << (s + 1 < snap_indices.size() ? ",\n" : "\n");
            }
            fs_json << "  ]\n";
            fs_json << "}\n";
            fs_json.close();
            std::cout << "Wrote Full-Scale (" << config.nx << "x" << study_ny << ") snapshots to "
                      << config.fullscale_json_output << std::endl;
        }
    }

    return results;
}
