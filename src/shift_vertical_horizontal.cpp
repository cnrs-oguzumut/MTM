#include "../include/experiments/shift_vertical_horizontal.h"

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <iostream>
#include <limits>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

#include "../include/defects/Defectanalysis.h"
#include "../include/geometry/DomainDimensions.h"
#include "../include/geometry/DomainInfo.h"
#include "../include/geometry/LatticeGenerator.h"
#include "../include/geometry/NeighborAnalyzer.h"
#include "../include/interatomic/inter_atomic.h"
#include "../include/lattice_energy/Strain_Energy_LatticeCalculator.h"
#include "../include/mesh/ElementTriangle2D.h"
#include "../include/mesh/MeshGenerator.h"
#include "../include/mesh/Remesher.h"
#include "../include/optimization/LBFGSOptimizer.h"
#include "../include/optimization/LatticeOptimizer.h"
#include "../include/output/configuration_saver.h"

Eigen::Matrix<double, 3, 2>
calculateShapeDerivatives(const Eigen::Vector2d &p1, const Eigen::Vector2d &p2,
                          const Eigen::Vector2d &p3);

std::tuple<double, Eigen::Matrix2d, int> perform_remeshing_loop_reduction(
    alglib::real_1d_array &x, const UserData *userData,
    const std::vector<int> &contact_atoms,
    const std::vector<int> &boundary_fixed_nodes, const Eigen::Matrix2d &F_ext,
    const Eigen::Matrix<double, 3, 2> &dndx,
    const std::array<double, 2> &offsets,
    const std::vector<int> &original_domain_map,
    const std::vector<std::tuple<double, double>> &translation_map,
    const Point2D &domain_dims_point, int &has_changes, int max_iterations,
    double reference_area, bool pbc, bool optimize_interior);

void writeSizesToFile(int Nx, int Ny);

namespace {

enum class ShiftBoundaryTest {
  LeftBottomPositive,
  LeftBottomNegative
};

struct BoundaryBox {
  double min_x = std::numeric_limits<double>::max();
  double max_x = std::numeric_limits<double>::lowest();
  double min_y = std::numeric_limits<double>::max();
  double max_y = std::numeric_limits<double>::lowest();
};

BoundaryBox computeBoundaryBox(const std::vector<Point2D> &points) {
  BoundaryBox box;
  for (const auto &point : points) {
    box.min_x = std::min(box.min_x, point.coord.x());
    box.max_x = std::max(box.max_x, point.coord.x());
    box.min_y = std::min(box.min_y, point.coord.y());
    box.max_y = std::max(box.max_y, point.coord.y());
  }
  return box;
}

std::vector<int> collectBoundaryNodesInWindow(
    const std::vector<Point2D> &points, double boundary_value, bool use_y,
    double boundary_tolerance, double x_min, double x_max) {
  std::vector<int> selected_nodes;
  for (size_t i = 0; i < points.size(); ++i) {
    const auto &point = points[i];
    const double boundary_coordinate = use_y ? point.coord.y() : point.coord.x();
    if (std::abs(boundary_coordinate - boundary_value) <= boundary_tolerance &&
        point.coord.x() >= x_min && point.coord.x() <= x_max) {
      selected_nodes.push_back(static_cast<int>(i));
    }
  }
  return selected_nodes;
}

std::vector<int> collectBoundaryNodesInYWindow(
    const std::vector<Point2D> &points, double boundary_x,
    double boundary_tolerance, double y_min, double y_max) {
  std::vector<int> selected_nodes;
  for (size_t i = 0; i < points.size(); ++i) {
    const auto &point = points[i];
    if (std::abs(point.coord.x() - boundary_x) <= boundary_tolerance &&
        point.coord.y() >= y_min && point.coord.y() <= y_max) {
      selected_nodes.push_back(static_cast<int>(i));
    }
  }
  return selected_nodes;
}

std::pair<std::vector<std::pair<int, int>>, std::vector<std::pair<int, int>>>
create_dof_mapping_from_fixed_nodes(const std::vector<Point2D> &points,
                                    const std::vector<int> &fixed_nodes) {
  std::set<int> fixed_node_set(fixed_nodes.begin(), fixed_nodes.end());
  std::vector<std::pair<int, int>> interior_mapping;
  std::vector<std::pair<int, int>> full_mapping;
  interior_mapping.reserve(points.size());
  full_mapping.reserve(points.size());

  int solver_index = 0;
  for (size_t i = 0; i < points.size(); ++i) {
    if (fixed_node_set.count(static_cast<int>(i)) > 0) {
      full_mapping.push_back(std::make_pair(static_cast<int>(i), -1));
    } else {
      interior_mapping.push_back(
          std::make_pair(static_cast<int>(i), solver_index));
      full_mapping.push_back(std::make_pair(static_cast<int>(i), solver_index));
      solver_index++;
    }
  }

  return {interior_mapping, full_mapping};
}

std::vector<int> makeUniqueNodes(const std::vector<int> &nodes) {
  std::set<int> unique_nodes(nodes.begin(), nodes.end());
  return std::vector<int>(unique_nodes.begin(), unique_nodes.end());
}

template <typename Work>
void runInSubdirectory(const std::filesystem::path &directory, Work work) {
  const std::filesystem::path original_path = std::filesystem::current_path();
  std::filesystem::create_directories(directory);
  std::filesystem::current_path(directory);

  try {
    work();
  } catch (...) {
    std::filesystem::current_path(original_path);
    throw;
  }

  std::filesystem::current_path(original_path);
}

const char *shiftBoundaryTestName(ShiftBoundaryTest test) {
  switch (test) {
  case ShiftBoundaryTest::LeftBottomPositive:
    return "left_bottom_positive";
  case ShiftBoundaryTest::LeftBottomNegative:
    return "left_bottom_negative";
  }

  return "unknown";
}

void shift_vertical_horizontal(int caller_id, int nx, int ny,
                               int horizontal_steps = 100,
                               int vertical_steps = 100,
                               ShiftBoundaryTest boundary_test =
                                   ShiftBoundaryTest::LeftBottomPositive,
                               double amplitude_lattice_spacings = 1.0,
                               bool remeshing_enabled = true,
                               bool perturb_initial_triangulation = false) {
  if (nx <= 0 || ny <= 0) {
    std::cerr << "Error: nx and ny must be positive integers." << std::endl;
    exit(EXIT_FAILURE);
  }

  writeSizesToFile(nx, ny);

  std::string lattice_type = "square";
  double h = 1.0;
  Eigen::Vector2d p1(0.0, 0.0);
  Eigen::Vector2d p2(h, 0.0);
  Eigen::Vector2d p3(0.0, h);
  Eigen::Matrix<double, 3, 2> dndx = calculateShapeDerivatives(p1, p2, p3);

  std::function<double(double)> potential_func = square_energy;
  std::function<double(double)> potential_func_der = square_energy_der;
  double optimal_lattice_parameter = 1.0;
  double lattice_constant = optimal_lattice_parameter;

  std::vector<Point2D> square_points = LatticeGenerator::generate_2d_lattice(
      nx, ny, lattice_constant, lattice_type);
  const std::vector<Point2D> reference_points = square_points;

  int original_domain_size = square_points.size();
  DomainInfo domain_size = compute_domain_size(square_points);
  const std::array<double, 2> offsets = {lattice_constant, lattice_constant};
  DomainDimensions domain_dims(domain_size.get_width(),
                               domain_size.get_height());
  Point2D domain_dims_point(domain_dims.size_x, domain_dims.size_y);
  bool pbc = false;

  auto [original_domain_map, translation_map] =
      MeshGenerator::create_domain_maps(original_domain_size, domain_dims,
                                        offsets);

  BoundaryBox box = computeBoundaryBox(reference_points);
  double boundary_tolerance = 0.5 * lattice_constant;
  double mid_x = 0.5 * (box.min_x + box.max_x);
  double mid_y = 0.5 * (box.min_y + box.max_y);

  double horizontal_boundary_x = box.min_x;
  double vertical_boundary_y = box.min_y;
  double horizontal_window_min_y = mid_y;
  double horizontal_window_max_y = box.max_y;
  double vertical_window_min_x = mid_x;
  double vertical_window_max_x = box.max_x;
  double horizontal_sign = 1.0;
  double vertical_sign = 1.0;

  switch (boundary_test) {
  case ShiftBoundaryTest::LeftBottomPositive:
    horizontal_sign = 1.0;
    vertical_sign = 1.0;
    break;
  case ShiftBoundaryTest::LeftBottomNegative:
    horizontal_sign = -1.0;
    vertical_sign = -1.0;
    break;
  }

  std::vector<int> horizontal_loaded_nodes = collectBoundaryNodesInYWindow(
      reference_points, horizontal_boundary_x, boundary_tolerance,
      horizontal_window_min_y, horizontal_window_max_y);
  std::vector<int> vertical_loaded_nodes = collectBoundaryNodesInWindow(
      reference_points, vertical_boundary_y, true, boundary_tolerance,
      vertical_window_min_x, vertical_window_max_x);

  std::vector<int> fixed_nodes;
  fixed_nodes.reserve(reference_points.size());
  for (size_t i = 0; i < reference_points.size(); ++i) {
    const auto &point = reference_points[i];
    bool is_horizontal_support =
        std::abs(point.coord.x() - horizontal_boundary_x) <=
        boundary_tolerance;
    bool is_vertical_support =
        std::abs(point.coord.y() - vertical_boundary_y) <= boundary_tolerance;
    if (is_horizontal_support || is_vertical_support) {
      fixed_nodes.push_back(static_cast<int>(i));
    }
  }
  fixed_nodes.insert(fixed_nodes.end(), horizontal_loaded_nodes.begin(),
                     horizontal_loaded_nodes.end());
  fixed_nodes.insert(fixed_nodes.end(), vertical_loaded_nodes.begin(),
                     vertical_loaded_nodes.end());
  fixed_nodes = makeUniqueNodes(fixed_nodes);

  std::cout << "shift_vertical_horizontal "
            << shiftBoundaryTestName(boundary_test)
            << " boundary setup:" << std::endl;
  std::cout << "  fixed-or-prescribed nodes: " << fixed_nodes.size()
            << std::endl;
  std::cout << "  horizontal boundary loaded nodes: "
            << horizontal_loaded_nodes.size() << std::endl;
  std::cout << "  vertical boundary loaded nodes: "
            << vertical_loaded_nodes.size() << std::endl;
  std::cout << "  amplitude/lattice spacing: " << amplitude_lattice_spacings
            << std::endl;
  std::cout << "  remeshing: " << (remeshing_enabled ? "enabled" : "disabled")
            << std::endl;
  std::cout << "  initial triangulation perturbation: "
            << (perturb_initial_triangulation ? "enabled" : "disabled")
            << std::endl;

  if (horizontal_loaded_nodes.empty() || vertical_loaded_nodes.empty()) {
    std::cerr << "Error: loading node selection is empty. Increase the loading "
                 "window or boundary tolerance."
              << std::endl;
    exit(EXIT_FAILURE);
  }

  auto [interior_mapping, full_mapping] =
      create_dof_mapping_from_fixed_nodes(square_points, fixed_nodes);

  std::cout << "free nodes: " << interior_mapping.size()
            << ", fixed/prescribed nodes: "
            << (full_mapping.size() - interior_mapping.size()) << std::endl;

  AdaptiveMesher mesher(domain_dims_point, offsets, original_domain_map,
                        translation_map, full_mapping,
                        1e-6, // Tolerance
                        pbc);
  mesher.setUsePeriodicCopies(pbc);
  if (perturb_initial_triangulation) {
    mesher.setTriangulationPerturbation(-1e-7);
  }

  alglib::real_1d_array free_dofs;
  int n_free_nodes = interior_mapping.size();
  free_dofs.setlength(2 * interior_mapping.size());
  map_points_to_solver_array(free_dofs, square_points, interior_mapping,
                             n_free_nodes);

  auto [elements, active_elements] = mesher.createMesh(
      square_points, free_dofs, Eigen::Matrix2d::Identity(), &dndx);
  double element_area = elements.empty() ? 0.0 : elements[0].getReferenceArea();

  for (auto &element : elements) {
    element.set_dof_mapping(full_mapping);
  }

  Strain_Energy_LatticeCalculator calculator(1.0);
  Eigen::Matrix2d F_ext = Eigen::Matrix2d::Identity();
  Eigen::Matrix2d C_I = F_ext.transpose() * F_ext;
  double zero = calculator.calculate_energy(C_I, potential_func, 0);

  horizontal_steps = std::max(1, horizontal_steps);
  vertical_steps = std::max(1, vertical_steps);

  std::vector<std::pair<double, double>> load_path;
  double horizontal_level = 0.0;
  double vertical_level = 0.0;
  while (horizontal_level < amplitude_lattice_spacings ||
         vertical_level < amplitude_lattice_spacings) {
    if (horizontal_level < amplitude_lattice_spacings) {
      double next_horizontal_level =
          std::min(amplitude_lattice_spacings, horizontal_level + 1.0);
      for (int step = 1; step <= horizontal_steps; ++step) {
        double alpha = static_cast<double>(step) /
                       static_cast<double>(horizontal_steps);
        double h_level =
            horizontal_level +
            alpha * (next_horizontal_level - horizontal_level);
        load_path.push_back({h_level * lattice_constant,
                             vertical_level * lattice_constant});
      }
      horizontal_level = next_horizontal_level;
    }

    if (vertical_level < amplitude_lattice_spacings) {
      double next_vertical_level =
          std::min(amplitude_lattice_spacings, vertical_level + 1.0);
      for (int step = 1; step <= vertical_steps; ++step) {
        double alpha = static_cast<double>(step) /
                       static_cast<double>(vertical_steps);
        double v_level =
            vertical_level + alpha * (next_vertical_level - vertical_level);
        load_path.push_back({horizontal_level * lattice_constant,
                             v_level * lattice_constant});
      }
      vertical_level = next_vertical_level;
    }
  }

  int total_steps = static_cast<int>(load_path.size());

  int file_counter = 0;
  double post_energy_previous = 0.0;

  for (int step = 0; step < total_steps; ++step) {
    double horizontal_displacement = load_path[step].first;
    double vertical_displacement = load_path[step].second;

    std::vector<Eigen::Vector2d> prescribed_displacements(
        reference_points.size(), Eigen::Vector2d::Zero());
    for (int node_idx : horizontal_loaded_nodes) {
      prescribed_displacements[node_idx].x() +=
          horizontal_sign * horizontal_displacement;
    }
    for (int node_idx : vertical_loaded_nodes) {
      prescribed_displacements[node_idx].y() +=
          vertical_sign * vertical_displacement;
    }
    for (int node_idx : fixed_nodes) {
      square_points[node_idx].coord =
          reference_points[node_idx].coord + prescribed_displacements[node_idx];
    }

    bool plasticity = false;
    UserData userData(square_points, elements, calculator, potential_func,
                      potential_func_der, zero, optimal_lattice_parameter,
                      F_ext, interior_mapping, full_mapping, active_elements,
                      plasticity);

    alglib::real_1d_array x;
    int n_vars = interior_mapping.size();
    x.setlength(2 * n_vars);
    map_points_to_solver_array(x, square_points, interior_mapping, n_vars);

    double pre_energy = 0.0;
    double post_energy = 0.0;
    double pre_area = 1.0;
    double post_area = 1.0;
    Eigen::Matrix2d stress_tensor = Eigen::Matrix2d::Zero();
    ConfigurationSaver::calculateEnergyAndStress(&userData, pre_energy,
                                                 stress_tensor, true);
    double pre_stress = stress_tensor(0, 1);
    pre_area = ConfigurationSaver::calculateTotalArea2D(&userData);

    int file_id = caller_id + file_counter;
    double saving_value =
        (horizontal_displacement + vertical_displacement) / lattice_constant;

    UserData preOptUserData(square_points, elements, calculator, potential_func,
                            potential_func_der, zero, optimal_lattice_parameter,
                            F_ext, interior_mapping, full_mapping,
                            active_elements, plasticity);
    runInSubdirectory("pre_relaxation", [&]() {
      ConfigurationSaver::saveConfigurationWithStressAndEnergy2D(
          &preOptUserData, file_id, pre_energy, pre_stress, true);
      ConfigurationSaver::saveTriangleData(&preOptUserData, file_id,
                                           domain_dims, offsets, full_mapping);
      ConfigurationSaver::saveElements(elements, active_elements, file_id);

      auto [num_dislocations_pre, coordination_pre] =
          DefectAnalysis::analyzeDefectsInReferenceConfig(
              &preOptUserData, file_id, dndx, offsets, original_domain_map,
              translation_map, domain_dims_point, element_area, pbc, true);
      ConfigurationSaver::writeToVTK(
          preOptUserData.points, preOptUserData.elements, &preOptUserData,
          file_id, true, coordination_pre, saving_value);
      ConfigurationSaver::logDislocationData(saving_value,
                                             num_dislocations_pre);
    });

    NeighborAnalyzer analyzer(NeighborAnalyzer::SearchType::K_NEAREST);
    std::vector<Point2D> points_before_copy = square_points;
    analyzer.setKNearestSearch(4);
    analyzer.setDebugMode(false);
    auto neighbors_before = analyzer.buildNeighbors(points_before_copy);

    userData.third_condition_flag = false;
    LBFGSOptimizer optimizer(13, 0.00001, 0, 0, 0);
    optimizer.optimize(x, minimize_energy_with_triangles, &userData);

    map_solver_array_to_points(x, square_points, interior_mapping, n_vars);

    auto neighbors_after = analyzer.buildNeighbors(square_points);
    auto change_info = NeighborAnalyzer::compareNeighborsWithTolerance(
        points_before_copy, square_points, neighbors_before, neighbors_after,
        0.135);

    stress_tensor.setZero();
    ConfigurationSaver::calculateEnergyAndStress(&userData, post_energy,
                                                 stress_tensor, true);
    double post_stress = stress_tensor(0, 1);
    post_area = ConfigurationSaver::calculateTotalArea2D(&userData);

    bool shouldRemesh =
        remeshing_enabled && (post_energy < post_energy_previous || step == 0);

    if (shouldRemesh) {
      std::vector<int> contact_atoms;
      int max_iterations = 1000;
      int hasChanges = 0;
      auto [post_energy_re, stress_tensor_re, iterations] =
          perform_remeshing_loop_reduction(
              x, &userData, contact_atoms, fixed_nodes, F_ext, dndx, offsets,
              original_domain_map, translation_map, domain_dims_point,
              hasChanges, max_iterations, element_area, pbc, true);

      post_energy = post_energy_re;
      post_stress = stress_tensor_re(0, 1);
      square_points = userData.points;
      elements = userData.elements;
      active_elements = userData.active_elements;
      if (!elements.empty()) {
        element_area = elements[0].getReferenceArea();
      }
      for (auto &element : elements) {
        element.set_dof_mapping(full_mapping);
      }
    }

    UserData postOptUserData(square_points, elements, calculator,
                             potential_func, potential_func_der, zero,
                             optimal_lattice_parameter, F_ext, interior_mapping,
                             full_mapping, active_elements, plasticity);
    ConfigurationSaver::saveConfigurationWithStressAndEnergy2D(
        &postOptUserData, file_id, post_energy, post_stress, true);
    ConfigurationSaver::saveTriangleData(&postOptUserData, file_id,
                                         domain_dims, offsets, full_mapping);
    ConfigurationSaver::saveElements(elements, active_elements, file_id);

    auto [num_dislocations_post, coordination_post] =
        DefectAnalysis::analyzeDefectsInReferenceConfig(
            &postOptUserData, file_id, dndx, offsets, original_domain_map,
            translation_map, domain_dims_point, element_area, pbc, true);
    ConfigurationSaver::writeToVTK(
        postOptUserData.points, postOptUserData.elements, &postOptUserData,
        file_id, true, coordination_post, saving_value);
    ConfigurationSaver::logDislocationData(saving_value, num_dislocations_post);

    ConfigurationSaver::logEnergyAndStress_v2(
        step, saving_value, pre_energy, pre_stress, post_energy, post_stress,
        pre_area, post_area, shouldRemesh);

    std::cout << "Step " << step << "/" << (total_steps - 1)
              << " horizontal_shift=" << horizontal_displacement
              << " vertical_push=" << vertical_displacement
              << " energy=" << post_energy << " stress_xy=" << post_stress
              << std::endl;

    file_counter++;
    post_energy_previous = post_energy;
  }
}

} // namespace

void run_final_shift_tests(int nx, int ny, int horizontal_steps,
                           int vertical_steps) {
  struct FinalShiftTest {
    const char *folder_name;
    ShiftBoundaryTest boundary_test;
    bool remeshing_enabled;
    bool perturb_initial_triangulation;
  };

  const double amplitude_lattice_spacings = 2.0;
  const FinalShiftTest tests[] = {
      {"left_bottom_perturbed_no_remesh_amp2",
       ShiftBoundaryTest::LeftBottomPositive, false, true},
      {"left_bottom_perturbed_remesh_amp2",
       ShiftBoundaryTest::LeftBottomPositive, true, true},
  };

  runInSubdirectory("final_tests", [&]() {
    for (const auto &test : tests) {
      std::cout << "\n=== Running final test: " << test.folder_name
                << " ===" << std::endl;
      runInSubdirectory(test.folder_name, [&]() {
        shift_vertical_horizontal(0, nx, ny, horizontal_steps, vertical_steps,
                                  test.boundary_test,
                                  amplitude_lattice_spacings,
                                  test.remeshing_enabled,
                                  test.perturb_initial_triangulation);
      });
    }
  });
}

