#include "../../include/experiments/experiment_includes.h"
#include "../../include/experiments/shifted_dislocation_study.h"
#include "../../include/utils/dislocation_utils.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <vector>
#include <numeric>

void run_shifted_dislocation_study(int caller_id, int nx, int ny,
                                   const std::vector<int> &shift_counts,
                                   bool use_cylinder,
                                   double R_free) {
  if (nx <= 0 || ny <= 0) {
    std::cerr << "Error: nx and ny must be positive integers." << std::endl;
    exit(EXIT_FAILURE);
  }

  std::cout << "\n========================================================" << std::endl;
  std::cout << ">>> RUNNING MULTI-SHIFT DISLOCATION STUDY (NO REMESH) <<<" << std::endl;
  std::cout << "Lattice: " << nx << " x " << ny << " (" << (nx * ny) << " atoms)" << std::endl;
  std::cout << "Mode: Upper crystal shift s * h (s in {";
  for (size_t i = 0; i < shift_counts.size(); ++i) {
    std::cout << shift_counts[i] << (i + 1 < shift_counts.size() ? ", " : "");
  }
  std::cout << "}) + Central Volterra edge dislocation" << std::endl;
  std::cout << "Boundary: " << (use_cylinder ? "Cylinder R_free=" + std::to_string(R_free) : "Fixed outer perimeter") << std::endl;
  std::cout << "Remeshing: DISABLED (Pristine lattice connectivity fixed)" << std::endl;
  std::cout << "========================================================\n" << std::endl;

  writeSizesToFile(nx, ny);

  std::string lattice_type = "square";
  double h = 1.0;
  double lattice_constant = 1.0;

  Eigen::Vector2d p1(0.0, 0.0);
  Eigen::Vector2d p2(h, 0.0);
  Eigen::Vector2d p3(0.0, h);
  Eigen::Matrix<double, 3, 2> dndx = calculateShapeDerivatives(p1, p2, p3);

  std::function<double(double)> potential_func = square_energy;
  std::function<double(double)> potential_func_der = square_energy_der;

  // Generate pristine reference points
  std::vector<Point2D> pristine_points =
      LatticeGenerator::generate_2d_lattice(nx, ny, lattice_constant, lattice_type);
  const int n_points = pristine_points.size();

  DomainInfo domain_size = compute_domain_size(pristine_points);
  const std::array<double, 2> offsets = {lattice_constant, lattice_constant};
  DomainDimensions domain_dims(domain_size.get_width(), domain_size.get_height());
  Point2D domain_dims_point(domain_dims.size_x, domain_dims.size_y);
  bool pbc = false;

  auto [original_domain_map, translation_map] =
      MeshGenerator::create_domain_maps(n_points, domain_dims, offsets);

  // Bounding box & core coordinates
  double x_min = std::numeric_limits<double>::max();
  double x_max = std::numeric_limits<double>::lowest();
  double y_min = std::numeric_limits<double>::max();
  double y_max = std::numeric_limits<double>::lowest();
  for (const auto &p : pristine_points) {
    x_min = std::min(x_min, p.coord.x());
    x_max = std::max(x_max, p.coord.x());
    y_min = std::min(y_min, p.coord.y());
    y_max = std::max(y_max, p.coord.y());
  }
  double core_x = 0.5 * (x_min + x_max);
  double core_y = 0.5 * (y_min + y_max);
  double y_mid = core_y;
  Eigen::Vector2d dislocation_core(core_x, core_y);

  // Master summary CSV
  std::filesystem::path current_dir = std::filesystem::current_path();
  std::string study_name = (nx >= 300) ? "shifted_dislocation_study_300x300" : "shifted_dislocation_study";
  bool in_study_dir = (current_dir.filename() == study_name);
  std::string base_dir = in_study_dir ? "" : (study_name + "/");
  if (!in_study_dir) {
    std::filesystem::create_directory(study_name);
  }

  std::ofstream summary_csv(base_dir + "shifted_dislocation_summary.csv");
  summary_csv << "shift,shift_amount,unrelaxed_energy,unrelaxed_stress,relaxed_energy,relaxed_stress,energy_drop\n";

  // Loop over each shift count s
  for (int shift_s : shift_counts) {
    double shift_amount = shift_s * lattice_constant;
    std::cout << "\n--------------------------------------------------------" << std::endl;
    std::cout << ">>> Processing Shift s = " << shift_s << " (" << shift_amount << " h) <<<" << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;

    std::string shift_dir = base_dir + "shift_" + std::to_string(shift_s);
    std::filesystem::create_directory(shift_dir);

    // Copy pristine points as baseline
    std::vector<Point2D> current_points = pristine_points;

    // 1. Shift upper half of crystal (y > y_mid) in x-direction
    size_t shifted_count = 0;
    for (auto &p : current_points) {
      if (p.coord.y() > y_mid) {
        p.coord.x() += shift_amount;
        ++shifted_count;
      }
    }
    std::cout << "Shifted upper half (" << shifted_count << " atoms above y = "
              << y_mid << ") by " << shift_amount << " h" << std::endl;

    // 2. Install exact anisotropic Stroh edge dislocation field centered at dislocation_core
    Eigen::Vector2d burgers_vector(lattice_constant, 0.0);
    for (size_t i = 0; i < current_points.size(); ++i) {
      Eigen::Vector2d rel_pos = pristine_points[i].coord - dislocation_core;
      Eigen::Vector2d u_aniso = VolterraDisplacement::calculateAnisotropicEdgeDisplacement(
          rel_pos, burgers_vector);
      current_points[i].coord += u_aniso;
    }
    std::cout << "Installed Anisotropic Stroh edge dislocation at (" << core_x << ", " << core_y << ")" << std::endl;

    // 3. Setup boundary conditions (DOF mapping)
    std::vector<std::pair<int, int>> interior_mapping;
    std::vector<std::pair<int, int>> full_mapping;
    interior_mapping.reserve(n_points);
    full_mapping.reserve(n_points);

    int solver_idx = 0;
    for (size_t i = 0; i < current_points.size(); ++i) {
      bool is_free = false;
      if (use_cylinder) {
        double dist = (pristine_points[i].coord - dislocation_core).norm();
        is_free = (dist <= R_free);
      } else {
        // Fix boundary frame (within 0.5 * h of domain boundary in pristine coordinates)
        double px = pristine_points[i].coord.x();
        double py = pristine_points[i].coord.y();
        bool on_boundary = (px <= x_min + 0.5 * h) || (px >= x_max - 0.5 * h) ||
                           (py <= y_min + 0.5 * h) || (py >= y_max - 0.5 * h);
        is_free = !on_boundary;
      }

      if (is_free) {
        interior_mapping.push_back({static_cast<int>(i), solver_idx});
        full_mapping.push_back({static_cast<int>(i), solver_idx});
        solver_idx++;
      } else {
        full_mapping.push_back({static_cast<int>(i), -1});
      }
    }
    std::cout << "Partition: " << interior_mapping.size() << " free (relaxed) nodes, "
              << (n_points - interior_mapping.size()) << " frozen boundary nodes." << std::endl;

    // 4. Build mesh once on pristine crystal (no remeshing!)
    AdaptiveMesher mesher(domain_dims_point, offsets, original_domain_map,
                          translation_map, full_mapping, 1e-6, pbc);
    mesher.setUsePeriodicCopies(pbc);

    alglib::real_1d_array x_ref;
    x_ref.setlength(2 * interior_mapping.size());
    map_points_to_solver_array(x_ref, pristine_points, interior_mapping, interior_mapping.size());

    auto [elements, active_elements] = mesher.createMesh(
        pristine_points, x_ref, Eigen::Matrix2d::Identity(), &dndx);
    double element_area = elements.empty() ? 0.0 : elements[0].getReferenceArea();

    for (auto &element : elements) {
      element.set_reference_mesh(current_points);
      element.set_dof_mapping(full_mapping);
    }
    std::cout << "Mesh created: " << elements.size() << " elements (fixed crystal connectivity)." << std::endl;

    // 5. Setup energy calculator
    Strain_Energy_LatticeCalculator calculator(1.0);
    Eigen::Matrix2d F_I = Eigen::Matrix2d::Identity();
    Eigen::Matrix2d C_I = F_I.transpose() * F_I;
    double zero = calculator.calculate_energy(C_I, potential_func, 0);

    Eigen::Matrix2d F_ext = Eigen::Matrix2d::Identity();
    bool plasticity = false;
    UserData userData(current_points, elements, calculator, potential_func,
                      potential_func_der, zero, lattice_constant, F_ext,
                      interior_mapping, full_mapping, active_elements, plasticity);

    // Initial state energy & stress
    alglib::real_1d_array x;
    int n_vars = interior_mapping.size();
    x.setlength(2 * n_vars);
    map_points_to_solver_array(x, current_points, interior_mapping, n_vars);

    double pre_energy = 0.0;
    Eigen::Matrix2d stress_tensor = Eigen::Matrix2d::Zero();
    ConfigurationSaver::calculateEnergyAndStress(&userData, pre_energy, stress_tensor, true);
    double pre_stress = stress_tensor(0, 1);
    std::cout << "State 0 (Unrelaxed): Energy = " << pre_energy << ", Stress_12 = " << pre_stress << std::endl;

    // Enter subfolder to write State 0
    std::filesystem::path orig_dir = std::filesystem::current_path();
    std::filesystem::current_path(shift_dir);

    std::vector<size_t> all_elements_idx(elements.size());
    std::iota(all_elements_idx.begin(), all_elements_idx.end(), 0);
    std::vector<std::pair<int, int>> vis_full_mapping(n_points);
    for (int i = 0; i < n_points; ++i) vis_full_mapping[i] = {i, i};

    UserData visUserDataPre(current_points, elements, calculator, potential_func,
                            potential_func_der, zero, lattice_constant, F_ext,
                            interior_mapping, vis_full_mapping, all_elements_idx, plasticity);
    ConfigurationSaver::saveConfigurationWithStressAndEnergy2D(
        &visUserDataPre, 0, pre_energy, pre_stress, true);
    ConfigurationSaver::saveTriangleData(&visUserDataPre, 0, domain_dims, offsets, full_mapping);
    ConfigurationSaver::saveElements(elements, all_elements_idx, 0);

    auto [num_disloc_pre, coord_pre] = DefectAnalysis::analyzeDefectsInReferenceConfig(
        &visUserDataPre, 0, dndx, offsets, original_domain_map,
        translation_map, domain_dims_point, element_area, pbc, true);
    ConfigurationSaver::writeToVTK(visUserDataPre.points, visUserDataPre.elements, &visUserDataPre,
                                   0, true, coord_pre, 0.0);
    ConfigurationSaver::logDislocationData(0.0, num_disloc_pre);

    // 6. Relax using preconditioned L-BFGS
    std::cout << "Relaxing configuration with preconditioned L-BFGS..." << std::endl;
    userData.third_condition_flag = false;
    relaxation_begin_step(shift_s);
    relax_configuration(x, &userData, 13);

    map_solver_array_to_points(x, current_points, interior_mapping, n_vars);
    userData.points = current_points;

    double post_energy = 0.0;
    stress_tensor.setZero();
    ConfigurationSaver::calculateEnergyAndStress(&userData, post_energy, stress_tensor, true);
    double post_stress = stress_tensor(0, 1);
    std::cout << "State 1 (Relaxed): Energy = " << post_energy
              << " (Drop: " << (pre_energy - post_energy) << ")"
              << ", Stress_12 = " << post_stress << std::endl;

    // Save State 1 (Relaxed)
    UserData visUserDataPost(current_points, elements, calculator, potential_func,
                             potential_func_der, zero, lattice_constant, F_ext,
                             interior_mapping, vis_full_mapping, all_elements_idx, plasticity);
    ConfigurationSaver::saveConfigurationWithStressAndEnergy2D(
        &visUserDataPost, 1, post_energy, post_stress, true);
    ConfigurationSaver::saveTriangleData(&visUserDataPost, 1, domain_dims, offsets, full_mapping);
    ConfigurationSaver::saveElements(elements, all_elements_idx, 1);

    auto [num_disloc_post, coord_post] = DefectAnalysis::analyzeDefectsInReferenceConfig(
        &visUserDataPost, 1, dndx, offsets, original_domain_map,
        translation_map, domain_dims_point, element_area, pbc, true);
    ConfigurationSaver::writeToVTK(visUserDataPost.points, visUserDataPost.elements, &visUserDataPost,
                                   1, true, coord_post, 0.0);
    ConfigurationSaver::logDislocationData(0.0, num_disloc_post);
    ConfigurationSaver::logEnergyAndStress_v2(1, static_cast<double>(shift_s), pre_energy, pre_stress,
                                              post_energy, post_stress, 0.0, 0.0, false);

    // Return to main dir and write summary
    std::filesystem::current_path(orig_dir);
    summary_csv << shift_s << "," << shift_amount << ","
                << std::scientific << std::setprecision(8)
                << pre_energy << "," << pre_stress << ","
                << post_energy << "," << post_stress << ","
                << (pre_energy - post_energy) << "\n";
    summary_csv.flush();
  }

  summary_csv.close();
  std::cout << "\n========================================================" << std::endl;
  std::cout << ">>> MULTI-SHIFT DISLOCATION STUDY COMPLETE <<<" << std::endl;
  std::cout << "Summary saved to shifted_dislocation_summary.csv" << std::endl;
  std::cout << "========================================================\n" << std::endl;
}
