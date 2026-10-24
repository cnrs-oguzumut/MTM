#include "../../include/experiments/experiment_includes.h"
#include "../../include/optimization/StiffnessSpectrum.h"
#include "../../include/experiments/stability_monitor.h"

// ============================================================================
// single_dislocation_study
//
// Stripped-down version of the Conti-Zanzotto square-crystal experiment.
// It keeps the original setup pipeline (reference geometry, energy potential,
// lattice generation, periodic domain maps, DOF mapping, meshing and the
// strain-energy calculator) but REMOVES the external loading schedule, the
// initial random noise and the avalanche/stress-drop bookkeeping.
//
// Instead it installs a single edge dislocation (Volterra field) at the
// central atom, performs an initial L-BFGS relaxation followed by the
// remesh/minimize loop (remesh -> minimize -> remesh -> ...) from the
// Conti-Zanzotto experiment, and writes the resulting state to VTK together
// with a defect analysis.
// ============================================================================
void single_dislocation_study(int caller_id, int nx, int ny,
                              const std::string &lattice_type,
                              int potential_mode) {

  // Toggle for the remesh/minimize loop after the initial relaxation.
  const bool remeshing = true;

  // Parameters for lattice
  if (nx <= 0 || ny <= 0) {
    std::cerr << "Error: nx and ny must be positive integers." << std::endl;
    exit(EXIT_FAILURE);
  }

  // ==================== SETUP REFERENCE GEOMETRY ====================
  writeSizesToFile(nx, ny);

  double h = 1.0; // Reference element size

  // Reference triangle vertices (for shape function derivatives)
  Eigen::Vector2d p1(0.0, 0.0);
  Eigen::Vector2d p2(h, 0.0);
  Eigen::Vector2d p3(0.0, h);

  Eigen::Matrix<double, 3, 2> dndx = calculateShapeDerivatives(p1, p2, p3);
  std::cout << "Reference element shape derivatives (dN/dx):\n"
            << dndx << std::endl;

  // ==================== SETUP ENERGY POTENTIAL ====================
  std::function<double(double)> potential_func = square_energy;
  std::function<double(double)> potential_func_der = square_energy_der;
  std::function<double(double)> potential_func_sder = square_energy_sder;

  // ==================== DETERMINE OPTIMAL LATTICE PARAMETER ============
  double symmetry_constantx =
      (lattice_type == "triangular") ? pow(4.0 / 3.0, 1.0 / 4.0) : 1.0;
  double optimal_lattice_parameter = symmetry_constantx * 1.0;
  double lattice_constant = optimal_lattice_parameter;
  std::cout << "Optimal lattice parameter: " << lattice_constant << std::endl;

  // ==================== GENERATE INITIAL LATTICE ====================
  std::vector<Point2D> square_points = LatticeGenerator::generate_2d_lattice(
      nx, ny, lattice_constant, lattice_type);

  std::vector<Point2D> square_points_ref =
      LatticeGenerator::generate_2d_lattice(nx, ny, lattice_constant,
                                            lattice_type);

  int original_domain_size = square_points.size();
  std::cout << "Generated lattice with " << original_domain_size << " points"
            << std::endl;

  // Calculate domain properties
  DomainInfo domain_size = compute_domain_size(square_points);

  const std::array<double, 2> offsets =
      (lattice_type == "square")
          ? std::array<double, 2>{lattice_constant, lattice_constant}
          : std::array<double, 2>{lattice_constant / 2.0,
                                  (sqrt(3.0) / 2.0) * lattice_constant};
  std::cout << "PBC offsets: [" << offsets[0] << ", " << offsets[1] << "]"
            << std::endl;

  DomainDimensions domain_dims(domain_size.get_width(),
                               domain_size.get_height());

  bool pbc = false;

  // Create domain maps
  auto [original_domain_map, translation_map] =
      MeshGenerator::create_domain_maps(original_domain_size, domain_dims,
                                        offsets);

  // Boundary conditions / DOF mapping
  auto [interior_mapping, full_mapping] =
      create_dof_mapping_original(square_points, 0.5 * lattice_constant, pbc);
  std::cout << "interior_mapping.size(): " << interior_mapping.size()
            << std::endl;
  std::cout << "full_mapping.size(): " << full_mapping.size() << std::endl;

  Point2D domain_dims_point(domain_dims.size_x, domain_dims.size_y);

  // ==================== BUILD MESH (on the pristine lattice) ============
  // Mesh is generated BEFORE introducing the dislocation, so the connectivity
  // and reference configuration come from the perfect crystal.
  AdaptiveMesher mesher(domain_dims_point, offsets, original_domain_map,
                        translation_map, full_mapping,
                        1e-6, // Tolerance
                        pbc   // Use periodic copies
  );
  mesher.setUsePeriodicCopies(pbc);

  alglib::real_1d_array free_dofs;
  int n_free_nodes = interior_mapping.size();
  free_dofs.setlength(2 * interior_mapping.size());
  map_points_to_solver_array(free_dofs, square_points, interior_mapping,
                             n_free_nodes);

  alglib::real_1d_array original_x_remesh =
      mesher.saveOriginalPositions(free_dofs);

  auto [elements, active_elements] = mesher.createMesh(
      square_points, free_dofs, Eigen::Matrix2d::Identity(), &dndx);
  double element_area = elements[0].getReferenceArea();

  for (auto &element : elements) {
    element.set_dof_mapping(full_mapping);
  }
  std::cout << "Created " << elements.size() << " element triangles"
            << std::endl;

  // ==================== INSTALL A SINGLE DISLOCATION ====================
  // Dynamically compute elastic constants and dislocation solution parameters
  // from the active potential (Volterra for triangular, Stroh for square).
  DislocationSolutionParams disloc_params =
      DislocationSolutionParams::fromPotentialMode(potential_mode, lattice_type);
  disloc_params.printSummary(lattice_constant);

  Eigen::Vector2d burgers_vector(lattice_constant, 0.0); // unit Burgers, x-dir
  double core_radius = 6.0;                              // core regularization
  size_t middle_atom_index = findMiddleAtom(square_points, true);

  std::vector<Point2D> dislocated_points =
      createSingleDislocation(square_points_ref, burgers_vector,
                              middle_atom_index, core_radius, disloc_params);

  // Use the dislocated configuration as the current state.
  square_points = dislocated_points;
  std::cout << "Installed single edge dislocation at atom " << middle_atom_index
            << " (|b| = " << burgers_vector.norm() << ")" << std::endl;

  // ==================== SETUP ENERGY CALCULATION ====================
  Strain_Energy_LatticeCalculator calculator(1.0);

  Eigen::Matrix2d C_ref = Eigen::Matrix2d::Identity();
  if (lattice_type == "triangular") {
    double a2 = std::sqrt(4.0 / 3.0);
    C_ref << a2, 0.5 * a2, 0.5 * a2, a2;
  }
  double zero = calculator.calculate_energy(C_ref, potential_func, 0);
  std::cout << "zero energy value: " << zero << std::endl;

  // ==================== CREATE USER DATA (no external load) ============
  Eigen::Matrix2d F_ext = Eigen::Matrix2d::Identity();
  bool plasticity = false;
  UserData userData(square_points, elements, calculator, potential_func,
                    potential_func_der, zero, optimal_lattice_parameter, F_ext,
                    interior_mapping, full_mapping, active_elements,
                    plasticity);

  // ==================== PREPARE OPTIMIZATION ====================
  alglib::real_1d_array x;
  int n_vars = interior_mapping.size();
  x.setlength(2 * n_vars);
  map_points_to_solver_array(x, square_points, interior_mapping, n_vars);

  // Pre-relaxation energy / stress
  double pre_energy = 0.0;
  Eigen::Matrix2d stress_tensor = Eigen::Matrix2d::Zero();
  ConfigurationSaver::calculateEnergyAndStress(&userData, pre_energy,
                                               stress_tensor, true);
  double pre_stress = stress_tensor(0, 1);
  std::cout << "Pre-relaxation - Energy: " << pre_energy
            << ", Stress: " << pre_stress << std::endl;

  // ==================== SAVE INITIAL DISLOCATED STATE ====================
  int file_id = caller_id;
  double saving_value = 0.0;

  ConfigurationSaver::saveConfigurationWithStressAndEnergy2D(
      &userData, file_id, pre_energy, pre_stress, true);
  ConfigurationSaver::saveTriangleData(&userData, file_id, domain_dims, offsets,
                                       full_mapping);
  ConfigurationSaver::saveElements(elements, active_elements, file_id);

  auto [num_dislocations_pre, coordination_pre] =
      DefectAnalysis::analyzeDefectsInReferenceConfig(
          &userData, file_id, dndx, offsets, original_domain_map,
          translation_map, domain_dims_point, element_area, pbc, true);

  ConfigurationSaver::writeToVTK(userData.points, userData.elements, &userData,
                                 file_id, true, coordination_pre, saving_value);
  ConfigurationSaver::logDislocationData(0.0, num_dislocations_pre);
  std::cout << "Saved initial dislocated config " << file_id << std::endl;

  double pre_area = ConfigurationSaver::calculateTotalArea2D(&userData);

  // ==================== RELAX (initial L-BFGS minimization) ============
  userData.third_condition_flag = false;
  LBFGSOptimizer optimizer(13, 0.00001, 0, 0, 0);
  optimizer.optimize(x, minimize_energy_with_triangles, &userData);
  map_solver_array_to_points(x, square_points, interior_mapping, n_vars);

  // ==================== POST-RELAXATION ENERGY / STRESS ============
  double post_energy = 0.0;
  stress_tensor.setZero();
  ConfigurationSaver::calculateEnergyAndStress(&userData, post_energy,
                                               stress_tensor, true);
  double post_stress = stress_tensor(0, 1);
  std::cout << "Post-relaxation - Energy: " << post_energy
            << ", Stress: " << post_stress << std::endl;

  // ==================== REMESHING DECISION ====================
  // Energy-based decision, as in the Conti-Zanzotto loop. The baseline is the
  // dislocated configuration energy before relaxation; this single-shot study
  // plays the role of the "first iteration" (i == 0), which always permits an
  // initial remeshing attempt.
  double post_energy_previous = pre_energy;
  bool first_pass = true;
  bool shouldRemesh = (post_energy < post_energy_previous) || first_pass;

  std::cout << "=== REMESH DECISION ===" << std::endl;
  std::cout << "  pre_energy = " << std::scientific << pre_energy << std::endl;
  std::cout << "  post_energy = " << std::scientific << post_energy
            << std::endl;
  std::cout << "  post_energy_previous = " << std::scientific
            << post_energy_previous << std::endl;
  std::cout << "  difference = " << (post_energy - post_energy_previous)
            << std::endl;
  std::cout << "  post_energy < post_energy_previous? "
            << (post_energy < post_energy_previous) << std::endl;
  std::cout << "  shouldRemesh = " << shouldRemesh << std::endl;

  if (remeshing && shouldRemesh) {
    std::cout << "REMESHING STARTS" << std::endl;

    std::vector<int> contact_atoms;
    std::vector<int> boundary_fixed_nodes;
    int hasChanges = 0;
    int max_iterations = 1000;

    // Iterative remesh/minimize loop (remesh -> minimize -> remesh -> ...).
    auto [post_energy_re, stress_tensor_re, iterations] =
        perform_remeshing_loop_reduction(
            x, &userData, contact_atoms, boundary_fixed_nodes, F_ext, dndx,
            offsets, original_domain_map, translation_map, domain_dims_point,
            hasChanges, max_iterations, element_area, pbc, true);

    post_energy = post_energy_re;
    post_stress = stress_tensor_re(0, 1);

    // Sync local containers with the (possibly) remeshed UserData state.
    square_points = userData.points;
    elements = userData.elements;
    active_elements = userData.active_elements;
    if (!elements.empty()) {
      element_area = elements[0].getReferenceArea();
    }
    for (auto &element : elements) {
      element.set_dof_mapping(full_mapping);
    }

    if (hasChanges) {
      std::cout << "Remeshing accepted - energy decreased (" << iterations
                << " iterations)" << std::endl;
    } else {
      std::cout << "Remeshing rejected - original mesh kept" << std::endl;
    }
    std::cout << "Final energy: " << post_energy << ", stress: " << post_stress
              << std::endl;
  }

  // ==================== REPORT / LOG ENERGY DROP ====================
  double post_area = ConfigurationSaver::calculateTotalArea2D(&userData);
  std::cout << "Energy dropped from " << post_energy_previous << " to "
            << post_energy << std::endl;
  ConfigurationSaver::logEnergyAndStress_v2(0, 0.0, pre_energy, pre_stress,
                                            post_energy, post_stress, pre_area,
                                            post_area, shouldRemesh);

  // ==================== SAVE RELAXED STATE ====================
  int post_file_id = caller_id + 1;
  UserData postOptUserData(square_points, elements, calculator, potential_func,
                           potential_func_der, zero, optimal_lattice_parameter,
                           F_ext, interior_mapping, full_mapping,
                           active_elements, plasticity);

  ConfigurationSaver::saveTriangleData(&postOptUserData, post_file_id,
                                       domain_dims, offsets, full_mapping);
  ConfigurationSaver::saveElements(elements, active_elements, post_file_id);

  auto [num_dislocations_post, coordination_post] =
      DefectAnalysis::analyzeDefectsInReferenceConfig(
          &postOptUserData, post_file_id, dndx, offsets, original_domain_map,
          translation_map, domain_dims_point, element_area, pbc, true);

  ConfigurationSaver::writeToVTK(
      postOptUserData.points, postOptUserData.elements, &postOptUserData,
      post_file_id, true, coordination_post, saving_value);
  ConfigurationSaver::logDislocationData(0.0, num_dislocations_post);
  std::cout << "Saved relaxed dislocated config " << post_file_id << std::endl;
}

// ============================================================================
// Structures and Helpers for Dislocation Core and Profile Analysis
// ============================================================================
struct DislocationNodalFields {
  std::vector<double> energy;
  std::vector<double> s_xx;
  std::vector<double> s_yy;
  std::vector<double> s_xy;
  std::vector<double> s_hydro;
  std::vector<double> s_vm;
};

struct RadialProfilePoint {
  double r = 0.0;
  double log_r = 0.0;
  double cumulative_energy = 0.0;
  double energy_density = 0.0;
  double avg_atom_energy = 0.0;
  double s_xx = 0.0;
  double s_yy = 0.0;
  double s_xy = 0.0;
  double s_hydro = 0.0;
  double s_vm = 0.0;
  int atom_count = 0;
};

struct SlipProfilePoint {
  double dx = 0.0; // x - xc
  double dy = 0.0; // y - yc
  double x = 0.0;
  double y = 0.0;
  double r = 0.0;
  double energy = 0.0;
  double s_xx = 0.0;
  double s_yy = 0.0;
  double s_xy = 0.0;
  double s_hydro = 0.0;
  double s_vm = 0.0;
  std::string row; // "upper", "lower", "core", "average"
};

struct CoreFitResult {
  double core_energy = 0.0;
  double elastic_K = 0.0;
  double R2 = 0.0;
};

struct StateProfileData {
  double core_x = 0.0;
  double core_y = 0.0;
  double total_energy = 0.0;
  double total_stress = 0.0;
  CoreFitResult core_fit;
  std::vector<RadialProfilePoint> radial_points;
  std::vector<SlipProfilePoint> slip_points;
  std::vector<SlipProfilePoint> slip_avg_points;
};

static DislocationNodalFields compute_dislocation_nodal_fields(UserData* userData, bool reduction = true) {
  DislocationNodalFields fields;
  size_t n_points = userData->points.size();
  fields.energy.assign(n_points, 0.0);
  fields.s_xx.assign(n_points, 0.0);
  fields.s_yy.assign(n_points, 0.0);
  fields.s_xy.assign(n_points, 0.0);
  fields.s_hydro.assign(n_points, 0.0);
  fields.s_vm.assign(n_points, 0.0);
  std::vector<int> node_count(n_points, 0);

  double norm = userData->calculator.getUnitCellArea();
  if (norm <= 1e-12) norm = 1.0;

  // Loop over ALL elements in the crystal mesh so nodal fields (energy, stresses)
  // are evaluated across the full crystal domain, not just the relaxed active cylinder
  for (size_t elem_idx = 0; elem_idx < userData->elements.size(); ++elem_idx) {
    auto& element = userData->elements[elem_idx];
    if (!element.isInitialized()) continue;

    element.setExternalDeformation(userData->F_external);
    element.calculate_deformation_gradient(userData->points);

    const Eigen::Matrix2d& F = element.getDeformationGradient();
    Eigen::Matrix2d C = element.getMetricTensor();
    lagrange::Result result;
    Eigen::Matrix2d P;

    if (reduction) {
      result = lagrange::reduce(C);
      C = result.C_reduced;
    }

    double element_area = element.getReferenceArea();
    double element_energy = (userData->calculator.calculate_energy(
        C, userData->energy_function, userData->zero_energy) / norm) * element_area;

    Eigen::Matrix2d dE_dC = userData->calculator.calculate_derivative(
        C, userData->derivative_function) / norm;

    if (reduction) {
      P = 2.0 * F * result.m_matrix * dE_dC * result.m_matrix.transpose();
    } else {
      P = 2.0 * F * dE_dC;
    }

    double detF = F.determinant();
    if (std::abs(detF) < 1e-12) detF = 1.0;
    Eigen::Matrix2d cauchy_stress = (1.0 / detF) * P * F.transpose();

    for (int i = 0; i < 3; ++i) {
      int node_idx = element.getNodeIndex(i);
      if (node_idx >= 0 && node_idx < static_cast<int>(n_points)) {
        fields.energy[node_idx] += element_energy / 3.0;
        fields.s_xx[node_idx] += cauchy_stress(0, 0);
        fields.s_xy[node_idx] += cauchy_stress(0, 1);
        fields.s_yy[node_idx] += cauchy_stress(1, 1);
        node_count[node_idx]++;
      }
    }
  }

  for (size_t i = 0; i < n_points; ++i) {
    if (node_count[i] > 0) {
      fields.s_xx[i] /= node_count[i];
      fields.s_xy[i] /= node_count[i];
      fields.s_yy[i] /= node_count[i];
      fields.s_hydro[i] = 0.5 * (fields.s_xx[i] + fields.s_yy[i]);
      fields.s_vm[i] = std::sqrt(std::max(0.0,
          fields.s_xx[i]*fields.s_xx[i] - fields.s_xx[i]*fields.s_yy[i] +
          fields.s_yy[i]*fields.s_yy[i] + 3.0*fields.s_xy[i]*fields.s_xy[i]));
    }
  }
  return fields;
}

static StateProfileData extract_dislocation_profiles(
    UserData* userData,
    double initial_core_x, double initial_core_y,
    double burgers_vector,
    double R_free,
    double lattice_constant,
    const std::string& lattice_type,
    double total_energy,
    double total_stress) {

  StateProfileData sp;
  sp.total_energy = total_energy;
  sp.total_stress = total_stress;

  DislocationNodalFields fields = compute_dislocation_nodal_fields(userData, true);
  const auto& pts = userData->points;
  size_t n_points = pts.size();

  double r_max = R_free * lattice_constant;

  // Locate the dislocation core (x_c, y_c) from the concentrated nodal energy peak
  // (matching plot_dislocation_profiles.py find_core algorithm: np.argmax(NodalEnergy))
  size_t peak_node = 0;
  double max_nodal_energy = -1e30;
  for (size_t i = 0; i < n_points; ++i) {
    double dist_from_nominal = std::hypot(pts[i].coord.x() - initial_core_x, pts[i].coord.y() - initial_core_y);
    if (dist_from_nominal <= 0.75 * r_max) {
      if (fields.energy[i] > max_nodal_energy) {
        max_nodal_energy = fields.energy[i];
        peak_node = i;
      }
    }
  }

  double core_x = initial_core_x;
  double core_y = initial_core_y;

  if (max_nodal_energy > -1e20) {
    double px = pts[peak_node].coord.x();
    double py = pts[peak_node].coord.y();
    double row_spacing = (lattice_type == "square") ? lattice_constant : (std::sqrt(3.0) / 2.0 * lattice_constant);

    // Look for paired core atom across the slip plane with highest energy
    size_t comp_node = peak_node;
    double max_comp_energy = -1e30;
    for (size_t i = 0; i < n_points; ++i) {
      if (i == peak_node) continue;
      double dy = std::abs(pts[i].coord.y() - py);
      double dx = std::abs(pts[i].coord.x() - px);
      if (dy >= 0.35 * row_spacing && dy <= 1.35 * row_spacing && dx <= 0.85 * lattice_constant) {
        if (fields.energy[i] > max_comp_energy) {
          max_comp_energy = fields.energy[i];
          comp_node = i;
        }
      }
    }

    if (comp_node != peak_node && max_comp_energy > 0.05 * max_nodal_energy) {
      core_x = 0.5 * (px + pts[comp_node].coord.x());
      core_y = 0.5 * (py + pts[comp_node].coord.y());
    } else {
      core_x = px;
      core_y = py;
    }
  }

  sp.core_x = core_x;
  sp.core_y = core_y;

  // 1. Gather all atoms with distance from the detected true core center
  struct AtomDist {
    size_t idx;
    double dx, dy, r;
  };
  std::vector<AtomDist> atoms;
  atoms.reserve(n_points);

  for (size_t i = 0; i < n_points; ++i) {
    double dx = pts[i].coord.x() - core_x;
    double dy = pts[i].coord.y() - core_y;
    double r = std::hypot(dx, dy);
    atoms.push_back({i, dx, dy, r});
  }

  // Sort atoms by radial distance from core
  std::sort(atoms.begin(), atoms.end(), [](const AtomDist& a, const AtomDist& b) {
    return a.r < b.r;
  });

  // 2. Compute Radial Profile: cumulative energy E(r), shell density, shell stresses
  // Use geometric (logarithmically equispaced) shells so that the dislocation core vicinity
  // (r <= 3b) is sampled densely with >50 points, yielding uniform resolution along ln(r/b).
  int num_bins = 100;
  double r_min_shell = std::max(0.12 * burgers_vector, 0.05);
  double running_cumulative_energy = 0.0;
  size_t current_atom_idx = 0;
  double prev_r_bin = 0.0;

  for (int b = 1; b <= num_bins; ++b) {
    double frac = static_cast<double>(b) / static_cast<double>(num_bins);
    double r_bin = r_min_shell * std::pow(r_max / r_min_shell, frac);
    double r_prev = prev_r_bin;
    prev_r_bin = r_bin;
    double shell_area = M_PI * (r_bin * r_bin - r_prev * r_prev);

    double shell_energy = 0.0;
    double shell_sxx = 0.0, shell_syy = 0.0, shell_sxy = 0.0;
    double shell_hydro = 0.0, shell_vm = 0.0;
    int shell_count = 0;

    while (current_atom_idx < atoms.size() && atoms[current_atom_idx].r <= r_bin) {
      size_t ai = atoms[current_atom_idx].idx;
      double en = fields.energy[ai];
      running_cumulative_energy += en;
      shell_energy += en;
      shell_sxx += fields.s_xx[ai];
      shell_syy += fields.s_yy[ai];
      shell_sxy += fields.s_xy[ai];
      shell_hydro += fields.s_hydro[ai];
      shell_vm += fields.s_vm[ai];
      shell_count++;
      current_atom_idx++;
    }

    RadialProfilePoint rp;
    rp.r = r_bin;
    rp.log_r = std::log(std::max(1e-4, r_bin / burgers_vector));
    rp.cumulative_energy = running_cumulative_energy;
    rp.atom_count = static_cast<int>(current_atom_idx);

    if (shell_count > 0) {
      rp.energy_density = (shell_area > 1e-8) ? (shell_energy / shell_area) : 0.0;
      rp.avg_atom_energy = shell_energy / shell_count;
      rp.s_xx = shell_sxx / shell_count;
      rp.s_yy = shell_syy / shell_count;
      rp.s_xy = shell_sxy / shell_count;
      rp.s_hydro = shell_hydro / shell_count;
      rp.s_vm = shell_vm / shell_count;
    } else {
      rp.energy_density = 0.0;
      rp.avg_atom_energy = 0.0;
      if (!sp.radial_points.empty()) {
        rp.s_xx = sp.radial_points.back().s_xx;
        rp.s_yy = sp.radial_points.back().s_yy;
        rp.s_xy = sp.radial_points.back().s_xy;
        rp.s_hydro = sp.radial_points.back().s_hydro;
        rp.s_vm = sp.radial_points.back().s_vm;
      }
    }
    sp.radial_points.push_back(rp);
  }

  // 3. Fit Core Energy: E(r) = K * ln(r/b) + E_core
  double r_fit_min = 1.5 * burgers_vector;
  double r_fit_max = 0.95 * r_max;
  std::vector<double> X_fit, Y_fit;
  for (const auto& rp : sp.radial_points) {
    if (rp.r >= r_fit_min && rp.r <= r_fit_max) {
      X_fit.push_back(rp.log_r);
      Y_fit.push_back(rp.cumulative_energy);
    }
  }

  if (X_fit.size() >= 3) {
    double sum_x = 0.0, sum_y = 0.0;
    for (size_t i = 0; i < X_fit.size(); ++i) {
      sum_x += X_fit[i];
      sum_y += Y_fit[i];
    }
    double mean_x = sum_x / X_fit.size();
    double mean_y = sum_y / X_fit.size();

    double num = 0.0, den = 0.0, ss_tot = 0.0;
    for (size_t i = 0; i < X_fit.size(); ++i) {
      double dx = X_fit[i] - mean_x;
      double dy = Y_fit[i] - mean_y;
      num += dx * dy;
      den += dx * dx;
      ss_tot += dy * dy;
    }

    if (den > 1e-12) {
      sp.core_fit.elastic_K = num / den;
      sp.core_fit.core_energy = mean_y - sp.core_fit.elastic_K * mean_x;
      double ss_res = 0.0;
      for (size_t i = 0; i < X_fit.size(); ++i) {
        double pred = sp.core_fit.elastic_K * X_fit[i] + sp.core_fit.core_energy;
        ss_res += (Y_fit[i] - pred) * (Y_fit[i] - pred);
      }
      sp.core_fit.R2 = (ss_tot > 1e-12) ? std::max(0.0, 1.0 - ss_res / ss_tot) : 1.0;
    }
  }

  // 4. Extract Slip Plane Profile: atoms within vertical band |y - core_y| <= row_spacing * 1.05
  double row_spacing = (lattice_type == "square") ? lattice_constant : (std::sqrt(3.0) / 2.0 * lattice_constant);
  double slip_y_tol = 1.05 * row_spacing;

  std::vector<SlipProfilePoint> upper_row, lower_row;
  for (size_t i = 0; i < n_points; ++i) {
    double dy = pts[i].coord.y() - core_y;
    double dx = pts[i].coord.x() - core_x;
    // Include all atoms along the slip plane across the full crystal domain
    if (std::abs(dy) <= slip_y_tol) {
      SlipProfilePoint p;
      p.dx = dx;
      p.dy = dy;
      p.x = pts[i].coord.x();
      p.y = pts[i].coord.y();
      p.r = std::hypot(dx, dy);
      p.energy = fields.energy[i];
      p.s_xx = fields.s_xx[i];
      p.s_yy = fields.s_yy[i];
      p.s_xy = fields.s_xy[i];
      p.s_hydro = fields.s_hydro[i];
      p.s_vm = fields.s_vm[i];

      if (dy > 1e-4) {
        p.row = "upper";
        upper_row.push_back(p);
      } else if (dy < -1e-4) {
        p.row = "lower";
        lower_row.push_back(p);
      } else {
        p.row = "core";
      }
      sp.slip_points.push_back(p);
    }
  }

  // Sort slip points by dx
  std::sort(sp.slip_points.begin(), sp.slip_points.end(), [](const SlipProfilePoint& a, const SlipProfilePoint& b) {
    return a.dx < b.dx;
  });
  std::sort(upper_row.begin(), upper_row.end(), [](const SlipProfilePoint& a, const SlipProfilePoint& b) {
    return a.dx < b.dx;
  });
  std::sort(lower_row.begin(), lower_row.end(), [](const SlipProfilePoint& a, const SlipProfilePoint& b) {
    return a.dx < b.dx;
  });

  // 5. Slip Plane Average: pair upper and lower atoms to get mean energy and stress fields along slip line
  for (const auto& up : upper_row) {
    auto it = std::min_element(lower_row.begin(), lower_row.end(),
        [&up](const SlipProfilePoint& a, const SlipProfilePoint& b) {
          return std::abs(a.dx - up.dx) < std::abs(b.dx - up.dx);
        });

    if (it != lower_row.end() && std::abs(it->dx - up.dx) <= 0.65 * lattice_constant) {
      SlipProfilePoint avg;
      avg.dx = 0.5 * (up.dx + it->dx);
      avg.dy = 0.0;
      avg.x = 0.5 * (up.x + it->x);
      avg.y = core_y;
      avg.r = std::abs(avg.dx);
      avg.energy = 0.5 * (up.energy + it->energy);
      avg.s_xx = 0.5 * (up.s_xx + it->s_xx);
      avg.s_yy = 0.5 * (up.s_yy + it->s_yy);
      avg.s_xy = 0.5 * (up.s_xy + it->s_xy);
      avg.s_hydro = 0.5 * (up.s_hydro + it->s_hydro);
      avg.s_vm = 0.5 * (up.s_vm + it->s_vm);
      avg.row = "average";
      sp.slip_avg_points.push_back(avg);
    }
  }

  std::sort(sp.slip_avg_points.begin(), sp.slip_avg_points.end(), [](const SlipProfilePoint& a, const SlipProfilePoint& b) {
    return a.dx < b.dx;
  });

  return sp;
}

static void write_dislocation_profiles_json_and_csv(
    const std::string& filename_json,
    const std::string& filename_csv,
    const std::string& lattice_type,
    int potential_mode,
    double beta,
    double bulk_K,
    double R_free,
    double lattice_constant,
    double burgers_vector,
    double core_x,
    double core_y,
    double C11, double C12, double C66, double poisson_ratio,
    const StateProfileData& unrelaxed,
    const StateProfileData* relaxed = nullptr,
    bool macro_stress_relaxed = false,
    const Eigen::Matrix2d* final_F = nullptr,
    const Eigen::Matrix2d* final_stress = nullptr,
    const std::vector<double>* lowest_eigenvalues = nullptr) {

  std::ofstream jf(filename_json);
  if (!jf.is_open()) {
    std::cerr << "Warning: Could not open " << filename_json << " for writing." << std::endl;
    return;
  }

  jf << "{\n";
  jf << "  \"status\": \"" << (relaxed ? "completed" : "in_progress") << "\",\n";
  jf << "  \"lattice_type\": \"" << lattice_type << "\",\n";
  jf << "  \"potential_mode\": " << potential_mode << ",\n";
  jf << "  \"beta\": " << beta << ",\n";
  jf << "  \"bulk_K\": " << bulk_K << ",\n";
  jf << "  \"R_free\": " << R_free << ",\n";
  jf << "  \"lattice_constant\": " << lattice_constant << ",\n";
  jf << "  \"burgers_vector\": " << burgers_vector << ",\n";
  double active_core_x = relaxed ? relaxed->core_x : unrelaxed.core_x;
  double active_core_y = relaxed ? relaxed->core_y : unrelaxed.core_y;
  jf << "  \"core_x\": " << active_core_x << ",\n";
  jf << "  \"core_y\": " << active_core_y << ",\n";
  jf << "  \"macro_stress_relaxed\": " << (macro_stress_relaxed ? "true" : "false") << ",\n";
  if (final_F) {
    jf << "  \"final_F\": [[" << (*final_F)(0, 0) << ", " << (*final_F)(0, 1) << "], ["
       << (*final_F)(1, 0) << ", " << (*final_F)(1, 1) << "]],\n";
  }
  if (final_stress) {
    jf << "  \"final_macro_stress\": {\"xx\": " << (*final_stress)(0, 0)
       << ", \"yy\": " << (*final_stress)(1, 1)
       << ", \"xy\": " << (*final_stress)(0, 1) << "},\n";
  }
  if (lowest_eigenvalues && !lowest_eigenvalues->empty()) {
    jf << "  \"lowest_eigenvalues\": [";
    for (size_t i = 0; i < lowest_eigenvalues->size(); ++i) {
      jf << (*lowest_eigenvalues)[i] << (i + 1 < lowest_eigenvalues->size() ? ", " : "");
    }
    jf << "],\n";
    jf << "  \"lambda_min\": " << (*lowest_eigenvalues)[0] << ",\n";
    jf << "  \"is_stable\": " << (((*lowest_eigenvalues)[0] > 0.0) ? "true" : "false") << ",\n";
  }
  jf << "  \"elastic_constants\": {\n";
  jf << "    \"C11\": " << C11 << ",\n";
  jf << "    \"C12\": " << C12 << ",\n";
  jf << "    \"C66\": " << C66 << ",\n";
  jf << "    \"poisson_ratio\": " << poisson_ratio << "\n";
  jf << "  },\n";

  auto write_state_json = [&jf](const std::string& key, const StateProfileData& sp, bool is_last) {
    jf << "  \"" << key << "\": {\n";
    jf << "    \"core_x\": " << sp.core_x << ",\n";
    jf << "    \"core_y\": " << sp.core_y << ",\n";
    jf << "    \"total_energy\": " << sp.total_energy << ",\n";
    jf << "    \"total_stress\": " << sp.total_stress << ",\n";
    jf << "    \"core_fit\": {\n";
    jf << "      \"core_energy\": " << sp.core_fit.core_energy << ",\n";
    jf << "      \"elastic_K\": " << sp.core_fit.elastic_K << ",\n";
    jf << "      \"R2\": " << sp.core_fit.R2 << "\n";
    jf << "    },\n";

    // Radial points
    jf << "    \"radial\": [\n";
    for (size_t i = 0; i < sp.radial_points.size(); ++i) {
      const auto& rp = sp.radial_points[i];
      jf << "      {\"r\":" << rp.r
         << ",\"log_r\":" << rp.log_r
         << ",\"cumulative_energy\":" << rp.cumulative_energy
         << ",\"energy_density\":" << rp.energy_density
         << ",\"avg_atom_energy\":" << rp.avg_atom_energy
         << ",\"s_xx\":" << rp.s_xx
         << ",\"s_yy\":" << rp.s_yy
         << ",\"s_xy\":" << rp.s_xy
         << ",\"s_hydro\":" << rp.s_hydro
         << ",\"s_vm\":" << rp.s_vm
         << ",\"atom_count\":" << rp.atom_count
         << "}" << (i + 1 < sp.radial_points.size() ? ",\n" : "\n");
    }
    jf << "    ],\n";

    // Slip plane points
    jf << "    \"slip_plane\": [\n";
    for (size_t i = 0; i < sp.slip_points.size(); ++i) {
      const auto& sp_pt = sp.slip_points[i];
      jf << "      {\"dx\":" << sp_pt.dx
         << ",\"dy\":" << sp_pt.dy
         << ",\"x\":" << sp_pt.x
         << ",\"y\":" << sp_pt.y
         << ",\"r\":" << sp_pt.r
         << ",\"energy\":" << sp_pt.energy
         << ",\"s_xx\":" << sp_pt.s_xx
         << ",\"s_yy\":" << sp_pt.s_yy
         << ",\"s_xy\":" << sp_pt.s_xy
         << ",\"s_hydro\":" << sp_pt.s_hydro
         << ",\"s_vm\":" << sp_pt.s_vm
         << ",\"row\":\"" << sp_pt.row << "\""
         << "}" << (i + 1 < sp.slip_points.size() ? ",\n" : "\n");
    }
    jf << "    ],\n";

    // Slip plane average points
    jf << "    \"slip_avg\": [\n";
    for (size_t i = 0; i < sp.slip_avg_points.size(); ++i) {
      const auto& avg = sp.slip_avg_points[i];
      jf << "      {\"dx\":" << avg.dx
         << ",\"energy\":" << avg.energy
         << ",\"s_xx\":" << avg.s_xx
         << ",\"s_yy\":" << avg.s_yy
         << ",\"s_xy\":" << avg.s_xy
         << ",\"s_hydro\":" << avg.s_hydro
         << ",\"s_vm\":" << avg.s_vm
         << "}" << (i + 1 < sp.slip_avg_points.size() ? ",\n" : "\n");
    }
    jf << "    ]\n";
    jf << "  }" << (is_last ? "\n" : ",\n");
  };

  write_state_json("unrelaxed", unrelaxed, relaxed == nullptr);
  if (relaxed) {
    write_state_json("relaxed", *relaxed, true);
  }
  jf << "}\n";
  jf.close();

  // Write CSV file
  std::ofstream cf(filename_csv);
  if (cf.is_open()) {
    cf << "state,profile_type,coordinate,dx,dy,energy,cumulative_energy,sigma_xx,sigma_yy,sigma_xy,sigma_hydro,sigma_vm,row\n";

    auto write_state_csv = [&cf](const std::string& state_name, const StateProfileData& sp) {
      for (const auto& rp : sp.radial_points) {
        cf << state_name << ",radial," << rp.r << ",0,0," << rp.avg_atom_energy << ","
           << rp.cumulative_energy << "," << rp.s_xx << "," << rp.s_yy << "," << rp.s_xy << ","
           << rp.s_hydro << "," << rp.s_vm << ",shell\n";
      }
      for (const auto& sp_pt : sp.slip_points) {
        cf << state_name << ",slip," << sp_pt.dx << "," << sp_pt.dx << "," << sp_pt.dy << ","
           << sp_pt.energy << ",0," << sp_pt.s_xx << "," << sp_pt.s_yy << "," << sp_pt.s_xy << ","
           << sp_pt.s_hydro << "," << sp_pt.s_vm << "," << sp_pt.row << "\n";
      }
      for (const auto& avg : sp.slip_avg_points) {
        cf << state_name << ",slip_avg," << avg.dx << "," << avg.dx << ",0,"
           << avg.energy << ",0," << avg.s_xx << "," << avg.s_yy << "," << avg.s_xy << ","
           << avg.s_hydro << "," << avg.s_vm << ",average\n";
      }
    };

    write_state_csv("unrelaxed", unrelaxed);
    if (relaxed) {
      write_state_csv("relaxed", *relaxed);
    }
    cf.close();
  }
}

// ============================================================================
// single_dislocation_cylinder_relaxation
//
// Installs an analytical edge dislocation on an nx x ny lattice:
// - Triangular lattice: isotropic Volterra solution (Eqs. A18, A23)
// - Square lattice: anisotropic Stroh / Eshelby-Read-Shockley solution (Eqs. A5-A16)
// Elastic moduli and dislocation parameters are computed dynamically from the
// chosen potential (Zanzotto continuum or atomistic pair potential).
// - Atoms outside radius R_free from the core are FROZEN (Dirichlet boundary).
// - Atoms inside radius R_free are FREE to relax by energy minimization.
// - Supports L-BFGS preconditioning (--precond=stiffness/laplacian/diag/none).
// ============================================================================
void single_dislocation_cylinder_relaxation(int caller_id, int nx, int ny,
                                            double R_free,
                                            bool enable_remeshing,
                                            bool export_full_mesh,
                                            const std::string &lattice_type,
                                            int potential_mode,
                                            double beta,
                                            double bulk_K,
                                            bool relax_macro_stress,
                                            double stress_tol,
                                            int max_macro_iters,
                                            const StabilityMonitorOptions &stability_options) {
  if (nx <= 0 || ny <= 0) {
    std::cerr << "Error: nx and ny must be positive integers." << std::endl;
    exit(EXIT_FAILURE);
  }

  writeSizesToFile(nx, ny);

  double h = 1.0;
  double lattice_constant = 1.0;
  double scale = 1.0;
  double r_cutoff = 1.86602540378444;

  // Determine lattice parameter and setup potential / calculator based on potential_mode
  double optimal_lattice_parameter = 1.0;
  std::function<double(double)> potential_func;
  std::function<double(double)> potential_func_der;
  std::function<double(double)> potential_func_sder;
  std::unique_ptr<BaseLatticeCalculator> calculator_ptr;

  if (potential_mode == 3) {
    // Zanzotto continuum potential
    energy_functions::current_beta = beta;
    energy_functions::current_K = bulk_K;
    double symmetry_constantx =
        (lattice_type == "triangular") ? std::pow(4.0 / 3.0, 0.25) : 1.0;
    optimal_lattice_parameter = symmetry_constantx;
    potential_func = square_energy;
    potential_func_der = square_energy_der;
    potential_func_sder = square_energy_sder;
    calculator_ptr = std::make_unique<Strain_Energy_LatticeCalculator>(1.0);
  } else {
    // Atomistic pair potential: mode 1 = LJ v2, mode 2 = LJ v3, mode 4 = Square (2-Well Gauss)
    if (potential_mode == 1) {
      r_cutoff = 2.5;
      potential_func = lennard_jones_energy_v2;
      potential_func_der = lennard_jones_energy_der_v2;
      potential_func_sder = lennard_jones_energy_sder_v2;
    } else if (potential_mode == 2) {
      r_cutoff = 1.86602540378444;
      potential_func = lennard_jones_energy_v3;
      potential_func_der = lennard_jones_energy_der_v3;
      potential_func_sder = lennard_jones_energy_sder_v3;
    } else if (potential_mode == 4) {
      r_cutoff = 2.5;
      potential_func = square_energy;
      potential_func_der = square_energy_der;
      potential_func_sder = square_energy_sder;
    } else {
      std::cerr << "Unknown potential mode: " << potential_mode << std::endl;
      exit(EXIT_FAILURE);
    }

    EquilibriumLatticeResult eq_lat = compute_equilibrium_lattice_parameter(
        potential_func, potential_func_der, potential_func_sder, lattice_type, r_cutoff);
    eq_lat.printSummary();
    optimal_lattice_parameter = eq_lat.lattice_parameter;
  }
  lattice_constant = optimal_lattice_parameter;

  // NOTE: Lagrange reduction (GL(2, Z)) in the minimization operates on the integer reference grid Z^2.
  // We ALWAYS use SquareLatticeCalculator over Z^2. For triangular lattice, calc_scale = a_0 / (4/3)^{1/4}
  // with reference metric C_ref = [[sqrt(4/3), 0.5*sqrt(4/3)], [0.5*sqrt(4/3), sqrt(4/3)]].
  // Reference element leg h matches calc_scale (or 1.0 for continuum model) so that undeformed physical elements evaluate to C_ref.
  double calc_scale = (potential_mode == 3)
      ? 1.0
      : ((lattice_type == "triangular") ? (optimal_lattice_parameter / std::pow(4.0 / 3.0, 0.25)) : optimal_lattice_parameter);
  h = (potential_mode == 3) ? 1.0 : calc_scale;

  if (potential_mode != 3) {
    calculator_ptr = std::make_unique<SquareLatticeCalculator>(calc_scale, r_cutoff);
  }
  BaseLatticeCalculator &calculator = *calculator_ptr;

  // Reference triangle shape derivatives: ALWAYS standard square reference triangle (p1=(0,0), p2=(h,0), p3=(0,h))
  // with leg length h matching reference element scale so that undeformed physical elements evaluate to C_ref
  Eigen::Vector2d p1(0.0, 0.0);
  Eigen::Vector2d p2(h, 0.0);
  Eigen::Vector2d p3(0.0, h);
  Eigen::Matrix<double, 3, 2> dndx = calculateShapeDerivatives(p1, p2, p3);

  // Generate pristine reference lattice
  std::vector<Point2D> square_points_ref =
      LatticeGenerator::generate_2d_lattice(nx, ny, lattice_constant, lattice_type);
  std::vector<Point2D> square_points = square_points_ref;
  int n_points = square_points.size();

  DomainInfo domain_size = compute_domain_size(square_points_ref);
  const std::array<double, 2> offsets =
      (lattice_type == "square")
          ? std::array<double, 2>{lattice_constant, lattice_constant}
          : std::array<double, 2>{lattice_constant / 2.0,
                                  (std::sqrt(3.0) / 2.0) * lattice_constant};
  DomainDimensions domain_dims(domain_size.get_width(), domain_size.get_height());
  Point2D domain_dims_point(domain_dims.size_x, domain_dims.size_y);
  bool pbc = false;

  auto [original_domain_map, translation_map] =
      MeshGenerator::create_domain_maps(n_points, domain_dims, offsets);

  // Calculate bounding box and place dislocation core at center on slip plane
  double x_min = std::numeric_limits<double>::max();
  double x_max = std::numeric_limits<double>::lowest();
  double y_min = std::numeric_limits<double>::max();
  double y_max = std::numeric_limits<double>::lowest();
  for (const auto &p : square_points_ref) {
    x_min = std::min(x_min, p.coord.x());
    x_max = std::max(x_max, p.coord.x());
    y_min = std::min(y_min, p.coord.y());
    y_max = std::max(y_max, p.coord.y());
  }
  double core_x = 0.5 * (x_min + x_max);
  double core_y = 0.5 * (y_min + y_max);
  Eigen::Vector2d dislocation_core(core_x, core_y);

  // Dynamically compute elastic moduli and dislocation field parameters from the potential
  DislocationSolutionParams disloc_params =
      (potential_mode == 3)
          ? DislocationSolutionParams::fromCurrentZanzottoPotential(lattice_type)
          : DislocationSolutionParams::fromPotentialMode(potential_mode, lattice_type);

  std::cout << "\n=== DISLOCATION CYLINDER RELAXATION ("
            << (lattice_type == "triangular" ? "VOLTERRA / TRIANGULAR" : "STROH / SQUARE")
            << ") ===" << std::endl;
  std::cout << "Lattice: " << lattice_type << " (" << nx << " x " << ny << ", " << n_points
            << " atoms, a = " << lattice_constant << ")" << std::endl;
  std::cout << "Potential Mode: " << potential_mode << std::endl;
  std::cout << "Dislocation core at: (" << core_x << ", " << core_y << ")" << std::endl;
  std::cout << "Relaxation radius R_free: " << R_free << " h" << std::endl;
  disloc_params.printSummary(lattice_constant);

  // Apply analytical displacement field to all atoms:
  // - Volterra (Eq. A23) for triangular lattice
  // - Stroh (Eqs. A14-A15) for square lattice
  Eigen::Vector2d burgers_vector(lattice_constant, 0.0);

  // Dynamically computed ground-state elastic constants for the selected potential
  double C11 = 0.0, C12 = 0.0, C66 = 0.0;
  double poisson_ratio = 0.0;

  if (potential_mode == 3) {
    auto ec = energy_functions::calculate_elastic_constants(beta, bulk_K, lattice_type);
    C11 = ec.C11;
    C12 = ec.C12;
    C66 = ec.C66;
    if (ec.C11 + ec.C12 > 1e-12) {
      poisson_ratio = ec.C12 / (ec.C11 + ec.C12);
    }
    std::cout << "Ground-state Elastic Constants: C11=" << C11 << ", C12=" << C12
              << ", C66=" << C66 << ", Zener A=" << ec.Zener_A
              << " | Poisson ratio nu=" << poisson_ratio << std::endl;
  } else {
    C11 = disloc_params.C11;
    C12 = disloc_params.C12;
    C66 = disloc_params.C66;
    poisson_ratio = disloc_params.poisson_ratio;
    std::cout << "Ground-state Elastic Constants: C11=" << C11 << ", C12=" << C12
              << ", C66=" << C66 << ", Zener A=" << disloc_params.Zener_A
              << " | Poisson ratio nu=" << poisson_ratio << std::endl;
  }

  for (size_t i = 0; i < square_points.size(); ++i) {
    Eigen::Vector2d rel_pos = square_points_ref[i].coord - dislocation_core;
    Eigen::Vector2d u_dislocation = VolterraDisplacement::calculateDisplacement(
        rel_pos, burgers_vector, disloc_params);
    square_points[i].coord = square_points_ref[i].coord + u_dislocation;
  }

  // Partition into free nodes (r <= R_free * h) and fixed nodes (r > R_free * h)
  std::vector<std::pair<int, int>> interior_mapping;
  std::vector<std::pair<int, int>> full_mapping;
  interior_mapping.reserve(n_points);
  full_mapping.reserve(n_points);

  std::vector<int> fixed_nodes;
  int solver_idx = 0;
  double r_free_dist = R_free * lattice_constant;
  for (size_t i = 0; i < square_points.size(); ++i) {
    double dist = (square_points_ref[i].coord - dislocation_core).norm();
    if (dist <= r_free_dist) {
      interior_mapping.push_back({static_cast<int>(i), solver_idx});
      full_mapping.push_back({static_cast<int>(i), solver_idx});
      solver_idx++;
    } else {
      full_mapping.push_back({static_cast<int>(i), -1});
      fixed_nodes.push_back(static_cast<int>(i));
    }
  }

  int n_free_nodes = interior_mapping.size();
  int n_fixed_nodes = fixed_nodes.size();
  std::cout << "Partition: " << n_free_nodes << " free (relaxed) nodes, "
            << n_fixed_nodes << " frozen (Dirichlet) nodes." << std::endl;

  // Build mesh with full DOF mapping
  AdaptiveMesher mesher(domain_dims_point, offsets, original_domain_map,
                        translation_map, full_mapping, 1e-6, pbc);
  mesher.setUsePeriodicCopies(pbc);

  alglib::real_1d_array x;
  int n_vars = interior_mapping.size();
  x.setlength(2 * n_vars);
  map_points_to_solver_array(x, square_points, interior_mapping, n_vars);

  alglib::real_1d_array x_ref;
  x_ref.setlength(2 * n_vars);
  map_points_to_solver_array(x_ref, square_points_ref, interior_mapping, n_vars);

  auto [elements, active_elements] = mesher.createMesh(
      square_points_ref, x_ref,
      Eigen::Matrix2d::Identity(), &dndx);
  double element_area = elements.empty() ? 0.0 : elements[0].getReferenceArea();

  for (auto &element : elements) {
    element.set_reference_mesh(square_points);
    element.set_dof_mapping(full_mapping);
  }
  std::cout << "Mesh created: " << elements.size() << " triangular elements ("
            << (enable_remeshing ? "reconnected Delaunay" : "fixed pristine lattice") << ")." << std::endl;

  // Setup energy calculation
  Eigen::Matrix2d C_ref = Eigen::Matrix2d::Identity();
  if (lattice_type == "triangular") {
    double a2 = std::sqrt(4.0 / 3.0);
    C_ref << a2, 0.5 * a2, 0.5 * a2, a2;
  }
  double zero = calculator.calculate_energy(C_ref, potential_func, 0);

  Eigen::Matrix2d F_ext = Eigen::Matrix2d::Identity();
  bool plasticity = false;
  UserData userData(square_points, elements, calculator, potential_func,
                    potential_func_der, zero, optimal_lattice_parameter, F_ext,
                    interior_mapping, full_mapping, active_elements, plasticity);
  userData.second_derivative_function = &potential_func_sder;

  // Compute pre-relaxation energy and stress
  double pre_energy = 0.0;
  Eigen::Matrix2d stress_tensor = Eigen::Matrix2d::Zero();
  ConfigurationSaver::calculateEnergyAndStress(&userData, pre_energy, stress_tensor, true);
  double pre_stress = stress_tensor(0, 1);
  double pre_area = ConfigurationSaver::calculateTotalArea2D(&userData);

  std::cout << "State 0 (Unrelaxed Volterra): Energy = " << pre_energy
            << ", Stress_12 = " << pre_stress << std::endl;

  // Prepare mesh indices for visualization based on export_full_mesh:
  // export_full_mesh = true  -> exports entire crystal (all elements)
  // export_full_mesh = false -> exports only the active elements (cylinder core)
  std::vector<size_t> vis_elements;
  std::vector<std::pair<int, int>> vis_full_mapping;
  if (export_full_mesh) {
    vis_elements.resize(elements.size());
    std::iota(vis_elements.begin(), vis_elements.end(), 0);
    vis_full_mapping.resize(n_points);
    for (int i = 0; i < n_points; ++i) {
      vis_full_mapping[i] = {i, i};
    }
  } else {
    vis_elements = active_elements;
    vis_full_mapping = full_mapping;
  }

  // Save State 0 (Unrelaxed Volterra state)
  int initial_file_id = caller_id;
  UserData visUserDataPre(square_points, elements, calculator, potential_func,
                          potential_func_der, zero, optimal_lattice_parameter, F_ext,
                          interior_mapping, vis_full_mapping, vis_elements, plasticity);
  visUserDataPre.second_derivative_function = &potential_func_sder;
  ConfigurationSaver::saveConfigurationWithStressAndEnergy2D(
      &visUserDataPre, initial_file_id, pre_energy, pre_stress, true);
  ConfigurationSaver::saveTriangleData(&visUserDataPre, initial_file_id, domain_dims, offsets, full_mapping);
  ConfigurationSaver::saveElements(elements, vis_elements, initial_file_id);

  auto [num_dislocations_pre, coordination_pre] =
      DefectAnalysis::analyzeDefectsInReferenceConfig(
          &visUserDataPre, initial_file_id, dndx, offsets, original_domain_map,
          translation_map, domain_dims_point, element_area, pbc, true);
  ConfigurationSaver::writeToVTK(visUserDataPre.points, visUserDataPre.elements, &visUserDataPre,
                                 initial_file_id, true, coordination_pre, 0.0);
  ConfigurationSaver::logDislocationData(0.0, num_dislocations_pre);
  std::cout << "Saved initial state as configuration_"
            << std::setw(5) << std::setfill('0') << initial_file_id << ".vtk ("
            << (export_full_mesh ? "full mesh" : "active elements only") << ")" << std::endl;

  // Extract and export State 0 (Unrelaxed Volterra) profile data
  StateProfileData unrelaxed_profile = extract_dislocation_profiles(
      &visUserDataPre, core_x, core_y, burgers_vector.x(), R_free,
      lattice_constant, lattice_type, pre_energy, pre_stress);
  write_dislocation_profiles_json_and_csv(
      "dislocation_profiles.json", "dislocation_profiles.csv",
      lattice_type, potential_mode, beta, bulk_K, R_free,
      lattice_constant, burgers_vector.x(), core_x, core_y,
      C11, C12, C66, poisson_ratio, unrelaxed_profile, nullptr);

  // Energy minimization: either standard fixed boundary or Full Cauchy Tensor macroscopic relaxation
  Eigen::Matrix2d F_bar = F_ext; // Starts at Identity
  double post_energy = 0.0;
  double post_stress = 0.0;

  if (relax_macro_stress) {
    std::cout << "\n========================================================" << std::endl;
    std::cout << "  STARTING FULL CAUCHY TENSOR MACRO-STRESS RELAXATION" << std::endl;
    std::cout << "  Target: sigma_xx -> 0, sigma_yy -> 0, sigma_xy -> 0" << std::endl;
    std::cout << "  Tolerance: " << stress_tol << " | Max Macro Iterations: " << max_macro_iters << std::endl;
    std::cout << "========================================================" << std::endl;

    bool macro_converged = false;
    for (int macro_iter = 0; macro_iter < max_macro_iters; ++macro_iter) {
      // Step A: Minimize internal free nodes with current boundary positions
      userData.third_condition_flag = false;
      relaxation_begin_step(macro_iter);
      relax_configuration(x, &userData, 13);

      map_solver_array_to_points(x, square_points, interior_mapping, n_vars);
      userData.points = square_points;

      // Step B: Calculate volume-averaged Cauchy stress tensor
      post_energy = 0.0;
      stress_tensor.setZero();
      ConfigurationSaver::calculateEnergyAndStress(&userData, post_energy, stress_tensor, true);
      post_stress = stress_tensor(0, 1);

      double max_stress_norm = std::max({
          std::abs(stress_tensor(0, 0)),
          std::abs(stress_tensor(1, 1)),
          std::abs(stress_tensor(0, 1))
      });

      std::cout << "  [Macro Step " << macro_iter << "] Energy = " << std::setprecision(10) << post_energy
                << " | Cauchy Stress: [xx=" << stress_tensor(0, 0)
                << ", yy=" << stress_tensor(1, 1)
                << ", xy=" << stress_tensor(0, 1) << "]"
                << " | Max|sigma| = " << max_stress_norm << std::endl;

      // Check convergence
      if (max_stress_norm < stress_tol) {
        macro_converged = true;
        std::cout << "  >>> Full Cauchy Tensor macro-stress relaxation CONVERGED in "
                  << (macro_iter + 1) << " iterations! Max|sigma| = " << max_stress_norm
                  << " < tol=" << stress_tol << std::endl;
        break;
      }

      if (macro_iter + 1 >= max_macro_iters) {
        std::cout << "  [Warning] Reached max macro iterations (" << max_macro_iters
                  << "). Stopping with Max|sigma| = " << max_stress_norm << std::endl;
        break;
      }

      // Step C: Compute instantaneous spatial tangent stiffness via AcousticTensor
      // For triangular lattice, the reference state mapping Z^2 to equilateral triangles
      // is F0 = (4/3)^{1/4} * [[1, 0.5], [0, sqrt(3)/2]]. The total deformation gradient
      // relative to Z^2 is F_tot = F_bar * F0, ensuring C_tot = C_ref at F_bar = I.
      Eigen::Matrix2d F0 = Eigen::Matrix2d::Identity();
      if (lattice_type == "triangular") {
        double s = std::pow(4.0 / 3.0, 0.25);
        F0 << s, 0.5 * s,
              0.0, s * std::sqrt(3.0) / 2.0;
      }
      Eigen::Matrix2d F_tot = F_bar * F0;
      Eigen::Matrix2d C_tot = F_tot.transpose() * F_tot;
      lagrange::Result red_res = lagrange::reduce(C_tot);
      Eigen::Matrix2d Z = red_res.m_matrix;
      double normalisation = calculator.getUnitCellArea();

      AcousticTensor acoustic_tensor(F_tot, red_res.C_reduced, Z);
      acoustic_tensor.computeEnergyDerivatives(
          calculator, potential_func_der, potential_func_sder, normalisation);

      itensor::ITensor C_spatial = acoustic_tensor.getSpatialTangent();

      // Step D: Extract 3x3 Voigt tangent stiffness
      Eigen::Matrix3d C_voigt = Eigen::Matrix3d::Zero();
      auto inds = itensor::inds(C_spatial);
      if (inds.size() == 4) {
        itensor::Index i_idx = inds[0];
        itensor::Index j_idx = inds[1];
        itensor::Index k_idx = inds[2];
        itensor::Index l_idx = inds[3];

        C_voigt(0, 0) = itensor::elt(C_spatial, i_idx = 1, j_idx = 1, k_idx = 1, l_idx = 1);
        C_voigt(0, 1) = itensor::elt(C_spatial, i_idx = 1, j_idx = 1, k_idx = 2, l_idx = 2);
        C_voigt(0, 2) = itensor::elt(C_spatial, i_idx = 1, j_idx = 1, k_idx = 1, l_idx = 2);

        C_voigt(1, 0) = itensor::elt(C_spatial, i_idx = 2, j_idx = 2, k_idx = 1, l_idx = 1);
        C_voigt(1, 1) = itensor::elt(C_spatial, i_idx = 2, j_idx = 2, k_idx = 2, l_idx = 2);
        C_voigt(1, 2) = itensor::elt(C_spatial, i_idx = 2, j_idx = 2, k_idx = 1, l_idx = 2);

        C_voigt(2, 0) = itensor::elt(C_spatial, i_idx = 1, j_idx = 2, k_idx = 1, l_idx = 1);
        C_voigt(2, 1) = itensor::elt(C_spatial, i_idx = 1, j_idx = 2, k_idx = 2, l_idx = 2);
        C_voigt(2, 2) = itensor::elt(C_spatial, i_idx = 1, j_idx = 2, k_idx = 1, l_idx = 2);
      } else {
        std::cerr << "Warning: C_spatial rank != 4. Using fallback identity tangent." << std::endl;
        C_voigt = Eigen::Matrix3d::Identity();
      }

      // Step E: Solve C_voigt * delta_eps = -sigma_residual
      Eigen::Vector3d delta_sigma;
      delta_sigma << stress_tensor(0, 0), stress_tensor(1, 1), stress_tensor(0, 1);

      Eigen::Vector3d delta_eps;
      if (std::abs(C_voigt.determinant()) > 1e-10) {
        delta_eps = -C_voigt.fullPivLu().solve(delta_sigma);
      } else {
        delta_eps = -C_voigt.bdcSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(delta_sigma);
      }

      // Step F: Apply relaxation factor and safeguard maximum strain increment per macro step
      double alpha_relax = 0.5;
      Eigen::Vector3d eps_corr = alpha_relax * delta_eps;
      double max_dE = 0.02;
      for (int ii = 0; ii < 3; ++ii) {
        if (std::abs(eps_corr(ii)) > max_dE) {
          eps_corr(ii) = (eps_corr(ii) > 0) ? max_dE : -max_dE;
        }
      }

      // Step G: Upper-triangular incremental deformation gradient
      // keeps the slip plane y = core_y strictly horizontal (F_21 = 0)
      Eigen::Matrix2d F_inc = Eigen::Matrix2d::Identity();
      F_inc(0, 0) += eps_corr(0);      // Delta eps_11
      F_inc(1, 1) += eps_corr(1);      // Delta eps_22
      F_inc(0, 1) += eps_corr(2);      // Delta gamma_12 = 2 * Delta eps_12

      Eigen::Matrix2d F_bar_new = F_inc * F_bar;

      std::cout << "  [Macro Step " << macro_iter << "] Delta eps: ["
                << eps_corr(0) << ", " << eps_corr(1) << ", " << eps_corr(2)
                << "] -> New F_bar: [" << F_bar_new(0, 0) << ", " << F_bar_new(0, 1)
                << "; " << F_bar_new(1, 0) << ", " << F_bar_new(1, 1) << "]" << std::endl;

      // Step H: Update node positions centered at the dislocation core
      for (size_t j = 0; j < square_points.size(); ++j) {
        Eigen::Vector2d rel = square_points[j].coord - dislocation_core;
        square_points[j].coord = dislocation_core + F_inc * rel;
      }

      // Re-map free nodes to solver array x
      map_points_to_solver_array(x, square_points, interior_mapping, n_vars);

      // Update state
      F_bar = F_bar_new;
      F_ext = F_bar;
      userData.points = square_points;
      for (auto &element : elements) {
        element.setExternalDeformation(F_ext);
        element.set_reference_mesh(square_points);
      }
    }
  } else {
    // Standard fixed-mantle relaxation (relax_macro_stress == false)
    std::cout << "Minimizing energy of free nodes within R_free..." << std::endl;
    userData.third_condition_flag = false;
    relaxation_begin_step(0);
    relax_configuration(x, &userData, 13);

    map_solver_array_to_points(x, square_points, interior_mapping, n_vars);
    userData.points = square_points;

    post_energy = 0.0;
    stress_tensor.setZero();
    ConfigurationSaver::calculateEnergyAndStress(&userData, post_energy, stress_tensor, true);
    post_stress = stress_tensor(0, 1);
  }

  std::cout << "State 1 (Relaxed core): Energy = " << post_energy
            << " (drop: " << (pre_energy - post_energy) << ")"
            << ", Stress_12 = " << post_stress << std::endl;

  // Optional remeshing/reconnection loop around core
  if (enable_remeshing) {
    std::cout << "Running remeshing/reconnection around the relaxed core..." << std::endl;
    std::vector<int> contact_atoms;
    int max_iterations = 100;
    int hasChanges = 0;
    auto [post_energy_re, stress_tensor_re, iterations] =
        perform_remeshing_loop_reduction(
            x, &userData, contact_atoms, fixed_nodes, F_ext, dndx, offsets,
            original_domain_map, translation_map, domain_dims_point,
            hasChanges, max_iterations, element_area, pbc, true);

    post_energy = post_energy_re;
    stress_tensor = stress_tensor_re;
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
    std::cout << "Remeshing finished (" << iterations << " iterations, changes: "
              << hasChanges << ", final energy: " << post_energy << ")" << std::endl;
  }

  // Save State 1 (Relaxed state)
  int relaxed_file_id = caller_id + 1;
  std::vector<size_t> vis_elements_post;
  if (export_full_mesh) {
    vis_elements_post.resize(elements.size());
    std::iota(vis_elements_post.begin(), vis_elements_post.end(), 0);
  } else {
    vis_elements_post = active_elements;
  }
  UserData postOptUserData(square_points, elements, calculator, potential_func,
                           potential_func_der, zero, optimal_lattice_parameter,
                           F_ext, interior_mapping, vis_full_mapping, vis_elements_post, plasticity);
  postOptUserData.second_derivative_function = &potential_func_sder;

  ConfigurationSaver::saveConfigurationWithStressAndEnergy2D(
      &postOptUserData, relaxed_file_id, post_energy, post_stress, true);
  ConfigurationSaver::saveTriangleData(&postOptUserData, relaxed_file_id, domain_dims, offsets, full_mapping);
  ConfigurationSaver::saveElements(elements, vis_elements_post, relaxed_file_id);

  auto [num_dislocations_post, coordination_post] =
      DefectAnalysis::analyzeDefectsInReferenceConfig(
          &postOptUserData, relaxed_file_id, dndx, offsets, original_domain_map,
          translation_map, domain_dims_point, element_area, pbc, true);

  ConfigurationSaver::writeToVTK(postOptUserData.points, postOptUserData.elements,
                                 &postOptUserData, relaxed_file_id, true, coordination_post, 1.0);
  ConfigurationSaver::logDislocationData(1.0, num_dislocations_post);

  double post_area = ConfigurationSaver::calculateTotalArea2D(&postOptUserData);
  ConfigurationSaver::logEnergyAndStress_v2(1, 1.0, pre_energy, pre_stress,
                                            post_energy, post_stress, pre_area,
                                            post_area, enable_remeshing);

  std::cout << "Saved relaxed state as configuration_"
            << std::setw(5) << std::setfill('0') << relaxed_file_id << ".vtk ("
            << (export_full_mesh ? "full mesh" : "active elements only") << ")" << std::endl;

  // Extract and export State 1 (Relaxed Core) profile data (and update combined JSON + CSV)
  StateProfileData relaxed_profile = extract_dislocation_profiles(
      &postOptUserData, core_x, core_y, burgers_vector.x(), R_free,
      lattice_constant, lattice_type, post_energy, post_stress);

  // ==================== STABILITY / EIGENVALUE ANALYSIS ====================
  std::vector<double> computed_eigenvalues;
  double computed_lambda_min = 0.0;
  bool has_eigenvalues = false;

  if (stability_options.enabled()) {
    std::cout << "\n========================================================" << std::endl;
    std::cout << "  STARTING DISLOCATION CORE STABILITY ANALYSIS (EIGENVALUES)" << std::endl;
    std::cout << "  Requested modes: " << stability_options.modes
              << " | Soft mode vectors: " << stability_options.vectors
              << " | Solver: Cholesky shift-invert Lanczos (Spectra)" << std::endl;
    std::cout << "========================================================" << std::endl;

    StiffnessSpectrumOptions eig_opt;
    eig_opt.n_modes = stability_options.modes;
    eig_opt.compute_vectors = (stability_options.vectors > 0);
    eig_opt.verbose = false;

    // Use fresh UserData with true solver DOF partition (Dirichlet fixed nodes have -1)
    UserData eigUserData(square_points, userData.elements, calculator, potential_func,
                         potential_func_der, zero, optimal_lattice_parameter,
                         F_ext, interior_mapping, full_mapping, active_elements, plasticity);
    eigUserData.second_derivative_function = &potential_func_sder;
    for (auto &elem : eigUserData.elements) {
      elem.setExternalDeformation(F_ext);
    }
    map_points_to_solver_array(x, square_points, interior_mapping, n_vars);

    StiffnessSpectrum spectrum = lowest_stiffness_modes(x, &eigUserData, eig_opt);

    if (spectrum.converged && spectrum.num_computed > 0) {
      has_eigenvalues = true;
      computed_lambda_min = spectrum.eigenvalues(0);
      computed_eigenvalues.resize(spectrum.num_computed);
      for (int k = 0; k < spectrum.num_computed; ++k) {
        computed_eigenvalues[k] = spectrum.eigenvalues(k);
      }

      std::cout << "Computed " << spectrum.num_computed << " lowest eigenvalues of relaxed core stiffness K:" << std::endl;
      for (int k = 0; k < spectrum.num_computed; ++k) {
        std::cout << "  λ_" << (k + 1) << " = " << std::scientific << std::setprecision(6)
                  << spectrum.eigenvalues(k)
                  << (k == 0 ? " (λ_min)" : "") << std::endl;
      }

      bool is_stable = (computed_lambda_min > 0.0);
      std::cout << "Core Stability Status:       " << (is_stable ? "STABLE (Positive Definite)" : "UNSTABLE / SADDLE (λ_min <= 0)")
                << " | λ_min = " << computed_lambda_min << std::endl;

      // Log to eigen_log.csv for Lattice Launcher live charting and live KPI monitoring
      std::ofstream eig_csv("eigen_log.csv");
      if (eig_csv.is_open()) {
        eig_csv << "Iteration,Alpha,Trigger,FailedShifts,Seconds";
        for (int k = 1; k <= spectrum.num_computed; ++k) {
          eig_csv << ",lambda_" << k;
        }
        eig_csv << "\n";
        eig_csv << "1,1.000000,dislocation_relaxed," << spectrum.failed_shifts << ","
                << std::scientific << std::setprecision(10)
                << (spectrum.assemble_seconds + spectrum.factorize_seconds + spectrum.lanczos_seconds);
        for (int k = 0; k < spectrum.num_computed; ++k) {
          eig_csv << "," << spectrum.eigenvalues(k);
        }
        eig_csv << "\n";
        eig_csv.close();
      }

      // Export soft modes to VTK in eigen_modes/
      if (stability_options.vectors > 0 && spectrum.eigenvectors.cols() > 0) {
        std::filesystem::create_directories("eigen_modes");
        std::string mode_vtk_path = "eigen_modes/soft_modes_00001.vtk";
        std::ofstream vtk_file(mode_vtk_path);
        if (vtk_file.is_open()) {
          const auto &points = square_points;
          const auto &mesh = userData.elements;
          const auto &active = active_elements;
          const int n_nodes = static_cast<int>(points.size());
          const int n_free = static_cast<int>(interior_mapping.size());
          const std::vector<size_t> &cells = active.empty() ? std::vector<size_t>() : active;

          vtk_file << "# vtk DataFile Version 3.0\n";
          vtk_file << "Soft modes of relaxed dislocation: lambda_min " << computed_lambda_min << "\n";
          vtk_file << "ASCII\nDATASET UNSTRUCTURED_GRID\n";
          vtk_file << std::setprecision(10);
          vtk_file << "POINTS " << n_nodes << " double\n";
          for (const auto &p : points) {
            vtk_file << p.coord.x() << " " << p.coord.y() << " 0\n";
          }
          vtk_file << "\nCELLS " << cells.size() << " " << 4 * cells.size() << "\n";
          for (size_t idx : cells) {
            vtk_file << "3 " << mesh[idx].getNodeIndex(0) << " " << mesh[idx].getNodeIndex(1) << " "
                     << mesh[idx].getNodeIndex(2) << "\n";
          }
          vtk_file << "\nCELL_TYPES " << cells.size() << "\n";
          for (size_t c = 0; c < cells.size(); ++c) {
            vtk_file << "5\n";
          }

          int n_export_modes = std::min(static_cast<int>(spectrum.eigenvalues.size()), stability_options.vectors);
          vtk_file << "\nPOINT_DATA " << n_nodes << "\n";
          for (int r = 0; r < n_export_modes; ++r) {
            const Eigen::VectorXd v = spectrum.eigenvectors.col(r);
            std::vector<Eigen::Vector2d> d(n_nodes, Eigen::Vector2d::Zero());
            double amax = 0.0;
            for (int j = 0; j < n_nodes; ++j) {
              const int dof = full_mapping[j].second;
              if (dof < 0 || dof >= n_free) continue;
              d[j] = Eigen::Vector2d(v(dof), v(dof + n_free));
              amax = std::max(amax, d[j].norm());
            }
            if (amax <= 1e-12) amax = 1.0;
            vtk_file << "VECTORS mode_" << (r + 1) << " double\n";
            for (const auto &dj : d) {
              vtk_file << (dj.x() / amax) << " " << (dj.y() / amax) << " 0\n";
            }
            vtk_file << "SCALARS mode_" << (r + 1) << "_amplitude double 1\nLOOKUP_TABLE default\n";
            for (const auto &dj : d) {
              vtk_file << (dj.norm() / amax) << "\n";
            }
          }
          vtk_file.close();
          std::cout << "  Exported " << n_export_modes << " soft modes to " << mode_vtk_path << std::endl;
        }
      }
    } else {
      std::cout << "Warning: Spectra failed to compute lowest stiffness eigenvalues." << std::endl;
    }
  }

  write_dislocation_profiles_json_and_csv(
      "dislocation_profiles.json", "dislocation_profiles.csv",
      lattice_type, potential_mode, beta, bulk_K, R_free,
      lattice_constant, burgers_vector.x(), core_x, core_y,
      C11, C12, C66, poisson_ratio, unrelaxed_profile, &relaxed_profile,
      relax_macro_stress, &F_ext, &stress_tensor,
      has_eigenvalues ? &computed_eigenvalues : nullptr);

  std::cout << "\n================ DISLOCATION ANALYSIS SUMMARY ================" << std::endl;
  std::cout << "Dislocation Core Position: (" << core_x << ", " << core_y << ")" << std::endl;
  std::cout << "State 0 (Unrelaxed): Total Energy = " << pre_energy
            << " | Core Fit E_core = " << unrelaxed_profile.core_fit.core_energy
            << ", Elastic K = " << unrelaxed_profile.core_fit.elastic_K
            << ", R² = " << unrelaxed_profile.core_fit.R2 << std::endl;
  std::cout << "State 1 (Relaxed):   Total Energy = " << post_energy
            << " | Core Fit E_core = " << relaxed_profile.core_fit.core_energy
            << ", Elastic K = " << relaxed_profile.core_fit.elastic_K
            << ", R² = " << relaxed_profile.core_fit.R2 << std::endl;
  std::cout << "Core Energy Relaxation Drop: ΔE_core = "
            << (unrelaxed_profile.core_fit.core_energy - relaxed_profile.core_fit.core_energy) << std::endl;
  if (relax_macro_stress) {
    std::cout << "Macro-Stress Relaxation:     ENABLED (Full Cauchy Tensor to zero stress)" << std::endl;
    std::cout << "Final F_bar:                 [" << F_ext(0, 0) << ", " << F_ext(0, 1)
              << "; " << F_ext(1, 0) << ", " << F_ext(1, 1) << "]" << std::endl;
    std::cout << "Final Cauchy Stress:         sigma_xx=" << stress_tensor(0, 0)
              << ", sigma_yy=" << stress_tensor(1, 1)
              << ", sigma_xy=" << stress_tensor(0, 1) << std::endl;
  }
  std::cout << "Dislocation profiles exported to dislocation_profiles.json and dislocation_profiles.csv" << std::endl;
  std::cout << "==============================================================\n" << std::endl;
  std::cout << "=== VOLTERRA DISLOCATION RELAXATION COMPLETE ===\n" << std::endl;
}

