#include "../../include/experiments/experiment_includes.h"

std::vector<Point2D> readPositionsFromFile(const std::string &filename) {
  // Maps in C++ are typically implemented
  // as self-balancing binary search trees (like red-black trees),
  // which dynamically allocate nodes as elements are inserted.
  std::map<int, Point2D> pointsMap; // Map to store points by atom ID
  std::ifstream file(filename);

  if (!file.is_open()) {
    std::cerr << "Error: Unable to open file " << filename << std::endl;
    return {}; // Return empty vector
  }

  std::string line;
  bool readingAtoms = false;

  while (std::getline(file, line)) {
    if (line.find("ITEM: ATOMS") != std::string::npos) {
      readingAtoms = true;
      continue;
    }

    if (readingAtoms) {
      int id;
      double x, y;
      std::istringstream iss(line);

      if (iss >> id >> x >> y) {
        pointsMap[id] = Point2D(x, y);
      }
    }
  }

  file.close();

  // Convert the map to a vector in order of atom ID
  std::vector<Point2D> points;
  points.reserve(pointsMap.size());

  // The std::map already keeps keys (atom IDs) in sorted order
  for (const auto &pair : pointsMap) {
    points.push_back(pair.second);
  }

  return points;
}

#include <Eigen/Dense>

// Eigen::Matrix2d calculateDeformationGradient(const Eigen::Matrix2d& P) {
//     // Elasticity matrix [A] (4x4 for non-symmetric P and F in 2D)
//     Eigen::Matrix4d A;
//     A << 27.734,  9.245,      0,      0,
//           9.245, 27.734,      0,      0,
//               0,      0, 36.979,      0,
//               0,      0,      0, 36.979;

//     // Convert input P (2x2) to Voigt notation (4x1)
//     Eigen::Vector4d P_voigt;
//     P_voigt << P(0,0), P(1,1), P(0,1), P(1,0);  // P11, P22, P12, P21

//     // Reference F0 = I (identity) in Voigt notation
//     Eigen::Vector4d I_voigt;
//     I_voigt << 1, 1, 0, 0;  // I11=1, I22=1, I12=0, I21=0

//     // Solve for F_voigt = A^{-1} * P_voigt + I_voigt
//     Eigen::Vector4d F_voigt = A.partialPivLu().solve(-P_voigt) + I_voigt;

//     // Convert F_voigt back to 2x2 matrix
//     Eigen::Matrix2d F;
//     F << F_voigt(0), F_voigt(2),  // F11, F12
//          F_voigt(3), F_voigt(1);  // F21, F22

//     return F;
// }

Eigen::Matrix2d
calculateDeformationGradientEigen(const Eigen::Matrix2d &sigma,
                                  const Eigen::Matrix2d &cauchy_applied) {
  // Stiffness matrix (for linear elasticity)
  Eigen::Matrix3d stiffness;
  stiffness << 110.9363766697, 36.9787922232, 0, 36.9787922232, 110.9363766697,
      0, 0, 0, 36.9787922232;

  // Convert Cauchy stress to Voigt notation
  Eigen::Vector3d sigma_voigt(sigma(0, 0) + cauchy_applied(0, 0),
                              sigma(1, 1) + cauchy_applied(1, 1),
                              sigma(0, 1) + cauchy_applied(0, 1));

  // Compute infinitesimal strain in Voigt notation
  Eigen::Vector3d eps_voigt = -stiffness.inverse() * sigma_voigt;

  // Convert strain to tensor form
  Eigen::Matrix2d eps;
  eps << eps_voigt(0), eps_voigt(2) / 2.0, eps_voigt(2) / 2.0, eps_voigt(1);

  // For large deformations, we compute the right Cauchy-Green tensor C directly
  // C = I + 2E where E is the Green-Lagrange strain tensor
  // For moderate strains, we can approximate E ≈ ε (the infinitesimal strain)
  Eigen::Matrix2d C = Eigen::Matrix2d::Identity() + 2.0 * eps;

  // Eigendecomposition of C
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> solver(C);
  Eigen::Vector2d eigenvalues = solver.eigenvalues();
  Eigen::Matrix2d eigenvectors = solver.eigenvectors();

  // Compute principal stretches (square root of eigenvalues of C)
  Eigen::Vector2d stretches = eigenvalues.cwiseSqrt();

  // Construct the right stretch tensor U
  Eigen::Matrix2d U =
      eigenvectors * stretches.asDiagonal() * eigenvectors.transpose();

  // For pure stretch (no rotation), F = U
  Eigen::Matrix2d F = U;

  return F;
}
Eigen::Matrix2d calculateDeformationGradient(const Eigen::Matrix2d &S) {
  // Stiffness matrix (Voigt notation)
  Eigen::Matrix3d stiffness;
  stiffness << 110.9363766697, 36.9787922232, 0, 36.9787922232, 110.9363766697,
      0, 0, 0, 36.9787922232;

  // Stress in Voigt notation [S11, S22, S12] (no factor of 2)
  Eigen::Vector3d S_voigt(S(0, 0), S(1, 1), S(0, 1));

  // Strain in Voigt notation [E11, E22, 2*E12]
  Eigen::Vector3d E_voigt = -stiffness.inverse() * S_voigt;

  // Convert to tensor form (E12 = 0.5 * E_voigt[2])
  Eigen::Matrix2d E;
  E << E_voigt(0), 0.5 * E_voigt(2), 0.5 * E_voigt(2), E_voigt(1);

  // Small strain: F ≈ I + E
  return Eigen::Matrix2d::Identity() + E;
}
//     // Eigen decomposition
//     Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> solver(FtF);
//     if (solver.info() != Eigen::Success) {
//         std::cerr << "Eigen decomposition failed!" << std::endl;
//         return Eigen::Matrix2d::Identity(); // Return identity in case of
//         failure
//     }

//     // Extract eigenvalues and eigenvectors
//     Eigen::Vector2d eigenvalues = solver.eigenvalues();
//     Eigen::Matrix2d eigenvectors = solver.eigenvectors();

//     // Compute the square root of the diagonal matrix of eigenvalues
//     Eigen::Matrix2d sqrtLambda = eigenvalues.cwiseSqrt().asDiagonal();

//     // Reconstruct the deformation gradient F
//     Eigen::Matrix2d F = eigenvectors * sqrtLambda * eigenvectors.transpose();

//     return F;
// }

#include <fstream>
#include <iomanip>

// void saveDisplacementsToFile(const std::vector<Point2D>& current_points,
//                            const std::vector<Point2D>& ref_points,
//                            const std::string& output_filename) {
//     std::ofstream outfile(output_filename);

//     // Write header
//     outfile << "# AtomID\tRef_x\tRef_y\tDx\tDy\n";

//     // Set precision for floating point output
//     outfile << std::fixed << std::setprecision(6);

//     // Write data for each point
//     size_t num_points = std::min(current_points.size(), ref_points.size());
//     for (size_t i = 0; i < num_points; ++i) {
//         double dx = current_points[i].coord.x() - ref_points[i].coord.x();
//         double dy = current_points[i].coord.y() - ref_points[i].coord.y();

//         outfile << i << "\t"                     // AtomID starting from 7069
//                 << ref_points[i].coord.x() << "\t"     // Ref_x
//                 << ref_points[i].coord.y() << "\t"     // Ref_y
//                 << dx << "\t"                         // Dx
//                 << dy << "\n";                        // Dy
//     }

//     outfile.close();
// }

std::vector<Point2D>
calculateDisplacements(const std::vector<Point2D> &current_points,
                       const std::vector<Point2D> &ref_points) {

  // Create vector to store displacements as Point2D objects
  std::vector<Point2D> displacements;
  displacements.reserve(std::min(current_points.size(), ref_points.size()));

  // Calculate displacement for each point
  size_t num_points = std::min(current_points.size(), ref_points.size());
  for (size_t i = 0; i < num_points; ++i) {
    double dx = current_points[i].coord.x() - ref_points[i].coord.x();
    double dy = current_points[i].coord.y() - ref_points[i].coord.y();
    displacements.push_back(Point2D(dx, dy));
  }

  return displacements;
}
// Sort Point2D objects by y-coordinate and then x-coordinate
void sortByCoordinates(std::vector<Point2D> &points, bool verbose = false) {
  std::sort(points.begin(), points.end(),
            [](const Point2D &a, const Point2D &b) {
              if (a.coord.y() != b.coord.y()) {
                return a.coord.y() < b.coord.y();
              }
              return a.coord.x() < b.coord.x();
            });

  if (verbose) {
    std::cout << "Sorted " << points.size()
              << " points by y-coordinate, then x-coordinate\n";
  }
}

std::vector<std::vector<Point2D>>
extractAtomRows(const std::vector<Point2D> &sorted_points,
                double y_tolerance = 0.1, bool verbose = false) {
  std::vector<std::vector<Point2D>> rows;
  if (sorted_points.empty()) {
    return rows;
  }

  std::vector<Point2D> current_row = {sorted_points[0]};
  double current_y = sorted_points[0].coord.y();

  for (size_t i = 1; i < sorted_points.size(); ++i) {
    if (std::abs(sorted_points[i].coord.y() - current_y) <= y_tolerance) {
      current_row.push_back(sorted_points[i]);
    } else {
      // Sort current row by x before saving
      std::sort(current_row.begin(), current_row.end(),
                [](const Point2D &a, const Point2D &b) {
                  return a.coord.x() < b.coord.x();
                });
      rows.push_back(std::move(current_row));
      current_row = {sorted_points[i]};
      current_y = sorted_points[i].coord.y();
    }
  }

  // Sort and add the last row
  if (!current_row.empty()) {
    std::sort(current_row.begin(), current_row.end(),
              [](const Point2D &a, const Point2D &b) {
                return a.coord.x() < b.coord.x();
              });
    rows.push_back(std::move(current_row));
  }

  // if (verbose) {
  // std::cout << "Extracted " << rows.size() << " rows (x-sorted within each
  // row)\n";
  // }

  return rows;
}

struct RowResults {
  std::vector<double> xmid;
  std::vector<double> deltax;
};

RowResults processPointRows(std::vector<Point2D> &square_points,
                            double lattice_constant) {
  RowResults results;

  // 2. Sort points by y-coordinate (ascending), then x-coordinate (ascending)
  sortByCoordinates(square_points, true); // `true` enables debug output

  // 3. Find the middle atom (optional, useful for symmetry analysis)
  size_t middle_idx = findMiddleAtom(square_points);

  // 4. Extract rows using a y-coordinate tolerance
  double y_tolerance = 0.42; // Adjust based on your system's row spacing
  std::vector<std::vector<Point2D>> rows =
      extractAtomRows(square_points, y_tolerance, false);

  // 5. Print the results (optional)
  for (size_t i = 0; i < rows.size(); ++i) {
    // std::cout << "Row " << i << " (" << rows[i].size() << " atoms): ";
    // std::cout << "\n";
    if (rows[i].size() != 99) {
      std::cout << "Row " << i << " (" << rows[i].size() << " atoms): ";
      exit(0);
    }
  }

  // 6. Find the middle row index
  // size_t middle_row_idx = rows.size() / 2-2 ;
  size_t middle_row_idx = rows.size() / 2 - 1;

  std::cout << "Middle row index: " << middle_row_idx << std::endl;

  const auto &mid_row = rows[middle_row_idx];
  const auto &upper_row =
      rows[middle_row_idx + 1]; // Assuming middle_row_idx is not the top row

  // Process and store x_mid and delta_x
  for (size_t i = 0; i < mid_row.size() - 20; ++i) {
    double x_mid = mid_row[i].coord.x();
    double y_mid = mid_row[i].coord.y();
    double y_upper = upper_row[i].coord.y();
    double x_upper =
        upper_row[i].coord.x(); // Assumes same number of atoms per row

    double delta_x =
        0 * lattice_constant + (x_upper - x_mid - lattice_constant / 2);

    // std::cout << "Middle row index x coo: " << mid_row[i].coord.x() <<
    // std::endl; std::cout << "Middle row index y coo: " <<
    // mid_row[i].coord.y() << std::endl;

    // std::cout << "Middle row up index x coo: " << upper_row[i].coord.x() <<
    // std::endl; std::cout << "Middle row up index y coo: " <<
    // upper_row[i].coord.y() << std::endl;
    // Store values in results
    results.xmid.push_back(x_mid);
    results.deltax.push_back(delta_x);

    // Optional: print for debugging
    // std::cout << x_mid << " " << delta_x << std::endl;
  }

  return results;
}

struct FitResult {
  double x_d;
  double eta;
};

struct NabarroModel {
  const std::vector<double> &x_data; // x-coordinates
  const std::vector<double> &y_data; // Displacements (simulation data)
  double b;                          // Fixed Burgers vector

  // Constructor
  NabarroModel(const std::vector<double> &x, const std::vector<double> &y,
               double fixed_b)
      : x_data(x), y_data(y), b(fixed_b) {}

  // Compute the cost function and its gradient
  void compute(const alglib::real_1d_array &params, double &func,
               alglib::real_1d_array &grad) const {
    double x_d = params[0];
    double eta = params[1];

    if (eta <= 0.0) {
      throw alglib::ap_error("Parameter eta must be positive.");
    }

    func = 0.0;
    grad[0] = 0.0;
    grad[1] = 0.0;

    for (size_t i = 0; i < x_data.size(); ++i) {
      double diff =
          y_data[i] -
          (-b / M_PI * (std::atan((x_data[i] - x_d) / eta) - M_PI / 2));
      func += diff * diff;

      // Partial derivatives
      double dx = x_data[i] - x_d;
      double denom = 1.0 + (dx / eta) * (dx / eta);

      grad[0] -= 2 * diff * b / (M_PI * denom) *
                 (1.0 / eta); // Partial derivative w.r.t x_d
      grad[1] -= 2 * diff * b / (M_PI * denom) *
                 (dx / (eta * eta)); // Partial derivative w.r.t eta
    }
  }
};

void callback(const alglib::real_1d_array &params, double &func,
              alglib::real_1d_array &grad, void *ptr) {
  auto *model = static_cast<NabarroModel *>(ptr);
  if (!model) {
    throw alglib::ap_error("Invalid model pointer.");
  }

  try {
    model->compute(params, func, grad);
  } catch (const std::exception &e) {
    throw alglib::ap_error(
        (std::string("Exception in callback: ") + e.what()).c_str());
  }
}

FitResult fitNabarroModel(const std::vector<double> &x_data,
                          const std::vector<double> &y_data, double b) {
  // Initialize the model
  NabarroModel model(x_data, y_data, b);

  // Print values in two columns
  //  std::cout<<"------"<<std::endl;
  //  for (size_t i = 0; i < x_data.size(); ++i) {
  //      std::cout << x_data[i] << "\t" << y_data[i] << std::endl;
  //  }

  // exit(0);

  // Initial parameter guesses: x_d and eta
  alglib::real_1d_array params;
  params.setlength(2);
  params[0] = 0.5;  // Initial guess for x_d
  params[1] = 10.0; // Initial guess for eta (must be positive)

  // Bounds: [x_d_min, eta_min], [x_d_max, eta_max]
  alglib::real_1d_array lower_bounds, upper_bounds;
  lower_bounds.setlength(2);
  upper_bounds.setlength(2);

  lower_bounds[0] = -1e6; // No restriction on x_d
  upper_bounds[0] = 1e6;  // No restriction on x_d

  lower_bounds[1] = 1e-6; // Minimum eta > 0
  upper_bounds[1] = 1e3;  // Large upper limit for eta

  // Set up ALGLIB optimizer
  alglib::minbleicstate state;
  alglib::minbleicreport report;

  try {
    alglib::minbleiccreate(params, state);
    alglib::minbleicsetbc(state, lower_bounds, upper_bounds); // Apply bounds
    alglib::minbleicsetcond(state, 1e-6, 0, 0, 10000); // Stopping criteria

    // Run optimization
    alglib::minbleicoptimize(state, callback, nullptr, &model);
    alglib::minbleicresults(state, params, report);
    std::cout << "BLEIC optimization terminated with code: "
              << report.terminationtype << "\n";
    std::cout << "Iterations: " << report.iterationscount << "\n";
    std::cout << "Function evaluations: " << report.nfev << "\n";

    // Output results
    double x_d_fit = params[0];
    double eta_fit = params[1];
    std::cout << "Optimization complete.\n";
    std::cout << "x_d = " << x_d_fit << "\n";
    std::cout << "eta = " << eta_fit << "\n";
    std::cout << "Burgers vector b = " << b << "\n";

    // Generate fitted data
    std::vector<double> y_fit(x_data.size());
    for (size_t i = 0; i < x_data.size(); ++i) {
      y_fit[i] =
          -b / M_PI * (std::atan((x_data[i] - x_d_fit) / eta_fit) - M_PI / 2);
    }

    // Write results to a file
    std::ofstream file("nabarro_fit_results.csv");
    file << "x_data,y_data,y_fit\n";
    for (size_t i = 0; i < x_data.size(); ++i) {
      file << x_data[i] << "," << y_data[i] << "," << y_fit[i] << "\n";
    }
    file.close();
    std::cout << "Results written to 'nabarro_fit_results.csv'.\n";

    // Return the fitted parameters
    return {x_d_fit, eta_fit};

  } catch (const alglib::ap_error &e) {
    std::cerr << "ALGLIB error: " << e.msg << "\n";
    return {std::nan(""), std::nan("")}; // Return NaN values on error
  }
}

void savePositionEtaStress(const std::string &filename, double x_d, double eta,
                           double stress_xy) {
  static bool first_call = true;

  try {
    // Open file in write mode for first call, append mode for subsequent calls
    std::ofstream file(filename, first_call ? std::ios::out : std::ios::app);

    if (!file.is_open()) {
      std::cerr << "Error: Cannot open file for writing: " << filename
                << std::endl;
      return;
    }

    // Set precision for floating-point output
    file << std::setprecision(10);

    // Write the three values in space-separated format
    file << x_d << " " << eta << " " << stress_xy << std::endl;
    file.close();

    // Reset first_call flag after first write
    first_call = false;

    // Also output to console
    std::cout << "Saved: x_d = " << x_d << ", eta = " << eta
              << ", stress_xy = " << stress_xy << std::endl;
  } catch (const std::exception &e) {
    std::cerr << "Error saving data: " << e.what() << std::endl;
  }
}

void single_dislo_LJ() {

  // Parameters for lattice
  std::string lattice_type = "triangular"; // "square" or "triangular"

  auto compute_even_ny = [](int nx) {
    int ny = std::round(2.0 * nx / std::sqrt(3));
    return (ny % 2 == 0) ? ny : ny + 1; // Ensure ny is even
  };

  int nx = 120;
  int ny = compute_even_ny(nx);

  // Energy functions (dumb)
  std::function<double(double)> potential_func = lennard_jones_energy_v3;
  std::function<double(double)> potential_func_der =
      lennard_jones_energy_der_v3;

  std::cout << "STEP 1: Finding optimal lattice parameter...\n";
  double h = 1.;
  double optimal_lattice_parameter = 0.996407146941421;
  double lattice_constant = 0.996407146941421;
  Eigen::Vector2d p1(0, 0);
  Eigen::Vector2d p2(h * lattice_constant, 0);
  Eigen::Vector2d p3(h * 0.5 * lattice_constant,
                     h * sqrt(3) * lattice_constant / 2);
  Eigen::Matrix<double, 3, 2> dndx = calculateShapeDerivatives(p1, p2, p3);
  std::cout << dndx << std::endl;

  // double symmetry_constantx = (lattice_type == "triangular") ?  pow(4.
  // / 3., 1. / 4.) : 1.0; double optimal_lattice_parameter = symmetry_constantx
  // * 1.0; double lattice_constant = optimal_lattice_parameter;

  // Generate initial lattice
  // std::vector<Point2D> square_points = LatticeGenerator::generate_2d_lattice(
  //     nx, ny, lattice_constant, lattice_type);
  // std::vector<Point2D> square_points_ref =
  // LatticeGenerator::generate_2d_lattice(
  //         nx, ny, lattice_constant, lattice_type);

  // Using the lambda function:

  std::string filename_ref =
      "/Users/usalman/draft_MTM_MS/Rfree_90/"
      "LocalState_Single_dislocation_cluster_0_Perfect_Crystal.dump";
  std::string filename =
      "/Users/usalman/draft_MTM_MS/Rfree_90/"
      "LocalState_Single_dislocation_cluster_1_With_Dislocation_Field.dump";

  std::vector<Point2D> square_points_ref = readPositionsFromFile(filename);
  std::vector<Point2D> square_points = readPositionsFromFile(filename_ref);
  double factor = h; // 90% of original size
  if (h != 1) {
    square_points =
        scaleLatticeAroundPoint(square_points, factor, Eigen::Vector2d(0, 0));
    std::cout << "STEP 2: No scaling...\n";
  }

  bool pbc = false;

  // RowResults registery = processPointRows(square_points, lattice_constant);
  // RowResults registery_diso = processPointRows(square_points_ref,
  // lattice_constant);

  // 2. Find the middle atom index
  size_t middle_atom_index =
      findMiddleAtom(square_points, true); // true enables verbose output
  // RowResults registery = processPointRows(square_points, lattice_constant);

  // std::vector<Point2D> square_points_ref =
  // readPositionsFromFile(filename_ref);

  // Now you can use the points vector
  std::cout << "Read " << square_points.size() << " points from file."
            << std::endl;

  std::cout << "Read " << square_points_ref.size() << " points from file."
            << std::endl;
  // Usage:
  // saveDisplacementsToFile(square_points, square_points_ref,
  // "displacement_field_volterra.txt");

  std::vector<Point2D> u =
      calculateDisplacements(square_points_ref, square_points);

  int original_domain_size = square_points.size();
  DomainInfo domain_size = compute_domain_size(square_points);

  const std::array<double, 2> offsets = {
      h * lattice_constant, h * (sqrt(3.) / 2.) * lattice_constant};
  std::cout << "offsets: " << offsets[0] << " " << offsets[1] << std::endl;
  DomainDimensions domain_dims(domain_size.get_width(),
                               domain_size.get_height());
  std::cout << "domain_size.get_width(): " << domain_size.get_width()
            << std::endl;
  std::cout << "domain_size.get_height(): " << domain_size.get_height()
            << std::endl;

  // Boundary conditions
  double radius = 90 * h;
  auto [interior_mapping, full_mapping] =
      create_dof_mapping_with_radius(square_points, radius, pbc);
  // Create domain maps
  auto [original_domain_map, translation_map] =
      MeshGenerator::create_domain_maps(original_domain_size, domain_dims,
                                        offsets);

  Eigen::Vector2d burgers_vector(
      lattice_constant, 0.0); // Example: unit Burgers vector in x-direction
  double core_radius = 1.0;   // Example core radius value
  double poisson_ratio = 0.3; // Typical value for many materials

  // Call the function
  // std::vector<Point2D> dislocated_points = createSingleDislocation(
  //     square_points_ref,
  //     burgers_vector,
  //     middle_atom_index,
  //     core_radius,
  //     poisson_ratio
  // );

  // square_points = dislocated_points;

  // Create ALGLIB array for free DOFs (displacements)
  alglib::real_1d_array free_dofs;
  int n_free_nodes = interior_mapping.size();
  free_dofs.setlength(2 *
                      interior_mapping.size()); // [u0, u1, ..., v0, v1, ...]
  map_points_to_solver_array(free_dofs, square_points, interior_mapping,
                             n_free_nodes);

  // Setup triangulation variables
  // Create a Point2D from your DomainDimensions locally
  Point2D domain_dims_point(domain_dims.size_x, domain_dims.size_y);

  // 1. Create the AdaptiveMesher instance
  AdaptiveMesher mesher(domain_dims_point, offsets, original_domain_map,
                        translation_map, full_mapping,
                        1e-6, // Tolerance
                        pbc);

  mesher.setUsePeriodicCopies(pbc); // Switch to using original domain only

  alglib::real_1d_array original_x_remesh =
      mesher.saveOriginalPositions(free_dofs);

  // Create mesh - this returns both elements and active_elements
  Eigen::Matrix2d F_temp;
  F_temp << 1.0, 0., // First row: identity
      0., 1;         // Second row: modified

  auto [elements, active_elements] = mesher.createMesh(
      square_points, free_dofs, Eigen::Matrix2d::Identity(), &dndx);

  std::cout << "Created " << elements.size() << " element triangles"
            << std::endl;

  // Add coordinates element-wise
  // Loop through each point in the square_points vector
  for (size_t i = 0; i < square_points.size(); ++i) {
    // Access the coord member of the current square_points element at index i
    // Access the coord member of the current u element at index i
    // Add the u[i] coordinate to the square_points[i] coordinate
    square_points[i].coord(0) += u[i].coord(0); // Add x-component
    square_points[i].coord(1) += u[i].coord(1); // Add y-component
  }
  // auto dipole_points = createDislocationDipole(
  //     square_points_ref,
  //     Eigen::Vector2d(lattice_constant, 0.0),
  //     Eigen::Vector2d(0.0, 0.0),
  //     40.0,  // separation distance
  //     2.0,  // core radius
  //     0.33,  // Poisson's ratio
  //     Eigen::Vector2d(1.0, 0.0) // dipole direction (45 degrees)
  // );
  // square_points =dipole_points;

  // Setup energy calculation
  TriangularLatticeCalculator calculator(lattice_constant);

  Eigen::Matrix2d F_I;
  F_I << 1.0, 0.,                              // First row: identity
      0., 1;                                   // Second row: modified
  Eigen::Matrix2d C_I = F_I.transpose() * F_I; // C = F^T * F
  double zero = calculator.calculate_energy(C_I, potential_func, 0);
  std::cout << "debugging simple shear test" << std::endl;
  std::cout << "zero energy value: " << zero << std::endl;

  // use this to find equlibrium atomistic distance
  // debug_deformation_tests_triangular() ;
  // std::cout << "debug_deformation_tests_triangular completed: "  <<
  // std::endl;
  //  exit(0);

  // Prepare for optimization
  alglib::real_1d_array x;
  int n_vars = interior_mapping.size();
  x.setlength(2 * n_vars);

  // Set up alpha values for deformation steps
  double alpha_min = 0.0;
  double alpha_max = 4.0;
  double step_size = 0.;

  int num_alpha_points;
  if (step_size > 0.00000000001)
    num_alpha_points =
        static_cast<int>((alpha_max - alpha_min) / step_size) + 1;
  else
    num_alpha_points = 80000;

  // Generate evenly spaced alpha values
  std::vector<double> alpha_values;
  double alpha_peirls = 0;
  double alpha_peirls_total = 0;

  alpha_values.reserve(num_alpha_points);
  for (int i = 0; i < num_alpha_points; i++) {
    alpha_values.push_back(alpha_min + i * step_size);
  }

  std::cout << "active_elements.size(): " << active_elements.size()
            << std::endl;
  std::vector<ElementTriangle2D> &elements_ref =
      elements; // Reference to elements
  std::vector<size_t> &active_elements_ref =
      active_elements; // Reference to active elements

  // Process each alpha value
  for (size_t i = 0; i < alpha_values.size(); i++) {
    double alpha = alpha_values[i];
    std::cout << "\n=== Processing alpha = " << alpha << " ===" << std::endl;

    // Create deformation gradient
    static Eigen::Matrix2d F_ext = Eigen::Matrix2d::Identity();
    // F_ext << 1.0, step_size,
    //          0.0, 1.0;
    std::cout << "F tensor [nounit]:" << std::endl;
    std::cout << "  Fxx = " << F_ext(0, 0) << std::endl;
    std::cout << "  Fxy = " << F_ext(0, 1) << std::endl;
    std::cout << "  Fyx = " << F_ext(1, 0) << std::endl;
    std::cout << "  Fyy = " << F_ext(1, 1) << std::endl;

    deform_boundary_nodes(square_points, full_mapping, F_ext);

    // for (size_t i = 0; i < square_points.size(); i++) {
    //     // Apply deformation: x_deformed = dF·x
    //     //square_points[i].coord = F_ext * square_points[i].coord;
    //     alpha*dislocated_points[i].coord
    // }

    // Prepare for optimization
    map_points_to_solver_array(x, square_points, interior_mapping, n_vars);
    for (auto &element : elements) {
      element.set_reference_mesh(square_points);
      element.set_dof_mapping(
          full_mapping); // or interior_mapping depending on needs
      // double jac = element.calculate_shape_derivatives(x);  // current
      // positions
    }
    const std::vector<size_t> new_active_elements = initialize_active_elements(
        elements, full_mapping, square_points.size());
    active_elements = new_active_elements;

    // Create user data
    bool plasticity = false;
    UserData userData(square_points, elements, calculator, potential_func,
                      potential_func_der, zero, optimal_lattice_parameter,
                      F_ext, interior_mapping, full_mapping, active_elements,
                      plasticity);
    double post_energyp = 0.0;
    double post_stressp = 0.0;
    if (i == 0) {
      ConfigurationSaver::saveConfigurationWithStressAndEnergy2D(
          &userData, -1, post_energyp, post_stressp);
      ConfigurationSaver::writeToVTK(userData.points, userData.elements,
                                     &userData, -1);
    }

    // Calculate pre-optimization energy and stress
    double pre_energy = 0.0;
    double pre_stress = 0.0;
    Eigen::Matrix2d pre_stress_tensor = Eigen::Matrix2d::Zero();
    ConfigurationSaver::calculateEnergyAndStress(&userData, pre_energy,
                                                 pre_stress_tensor);
    pre_stress = pre_stress_tensor(
        0, 1); // Extract xy component for backward compatibility std::cout <<
               // "Pre-optimization - Energy: " << pre_energy << ", Stress: " <<
               // pre_stress << std::endl;

    // Store original positions
    alglib::real_1d_array original_x;
    original_x.setlength(x.length());
    for (int j = 0; j < x.length(); j++) {
      original_x[j] = x[j];
    }

    std::vector<int> m3_before =
        analyzeElementReduction(elements, square_points, &userData);
    // Run optimization
    // --- Start timing ---
    auto wall_start = std::chrono::high_resolution_clock::now();
    clock_t cpu_start = clock();

    // Run optimization
    LBFGSOptimizer optimizer(10, 0, pow(12, -12), 0, 0);
    // CGOptimizer optimizer;
    optimizer.optimize(x, minimize_energy_with_triangles_noreduction,
                       &userData);

    // --- Stop timing ---
    auto wall_end = std::chrono::high_resolution_clock::now();
    clock_t cpu_end = clock();

    // Compute durations
    double wall_time =
        std::chrono::duration<double>(wall_end - wall_start).count();
    double cpu_time = (double)(cpu_end - cpu_start) / CLOCKS_PER_SEC;

    // Print results
    std::cout << "Optimization wall-clock time: " << wall_time << " seconds\n";
    std::cout << "Optimization CPU time: " << cpu_time << " seconds\n";
    std::cout << "Optimization Ratio: " << cpu_time / wall_time << " seconds\n";

    map_solver_array_to_points(x, square_points, interior_mapping, n_vars);

    std::vector<int> m3_after =
        analyzeElementReduction(elements, square_points, &userData);
    int hasChanges = compareM3Activation(m3_before, m3_after);
    if (hasChanges > 0) {
      std::cout << "Changes in m3 detected! " << hasChanges << std::endl;
    }

    // Calculate post-optimization energy and stress
    double post_energy = 0.0;
    Eigen::Matrix2d post_stress_tensor = Eigen::Matrix2d::Zero();
    double post_stress = 0.0;
    ConfigurationSaver::calculateEnergyAndStress(&userData, post_energy,
                                                 post_stress_tensor);
    post_stress = post_stress_tensor(
        0, 1); // Extract xy component for backward compatibility

    // Calculate change measures and check for remeshing need
    ChangeMeasures result =
        computeChangeMeasures(x, original_x, lattice_constant, elements,
                              &userData, square_points, true, &F_ext);
    std::cout << "max_abs_change: " << result.max_abs_change << std::endl;
    if (result.has_distorted_triangles) {
      std::cout << "Distorted triangles detected!" << std::endl;
    }

    // Determine if remeshing is needed
    bool shouldRemesh = result.has_distorted_triangles;
    std::cout << "shouldRemesh: " << shouldRemesh << std::endl;
    shouldRemesh = false;
    // Remeshing if needed
    int mesh_iteration = 0;
    shouldRemesh = true;
    shouldRemesh = checkSquareDomainViolation(elements);
    // shouldRemesh=false;
    const int MAX_ITERATIONS = 100;

    // 4. Remeshing loop
    while (shouldRemesh && mesh_iteration < MAX_ITERATIONS) {
      std::cout << "REMESHING STARTS iteration: " << mesh_iteration++
                << std::endl;

      // Save original positions
      alglib::real_1d_array original_x_remesh = mesher.saveOriginalPositions(x);

      // Create mesh - this returns both elements and active_elements
      auto [new_elements, new_active_elements] =
          mesher.createMesh(square_points, x, F_ext, &dndx);
      elements_ref = new_elements;               // Update through reference
      active_elements_ref = new_active_elements; // Update through reference

      std::vector<int> m3_before_remeshed =
          analyzeElementReduction(elements, square_points, &userData);

      // Re-optimize with new mesh
      plasticity = true;
      UserData newUserData(square_points, elements, calculator, potential_func,
                           potential_func_der, zero, optimal_lattice_parameter,
                           F_ext, interior_mapping, full_mapping,
                           active_elements, plasticity);

      LBFGSOptimizer optimizer_remesh(10, 0, pow(10, -13), 0, 0);
      optimizer_remesh.optimize(x, minimize_energy_with_triangles_noreduction,
                                &newUserData);
      map_solver_array_to_points(x, square_points, interior_mapping, n_vars);
      // recalculate
      std::vector<int> m3_after_remeshed =
          analyzeElementReduction(elements, square_points, &userData);
      hasChanges = compareM3Activation(m3_before_remeshed, m3_after_remeshed);

      ChangeMeasures result = computeChangeMeasures(
          x, original_x_remesh, lattice_constant, elements, &newUserData,
          square_points, true, &F_ext);
      std::cout << "max_abs_change: " << result.max_abs_change << std::endl;
      if (result.has_distorted_triangles) {
        std::cout << "Distorted triangles detected!" << std::endl;
      }

      // Determine if remeshing is needed
      // shouldRemesh = result.has_distorted_triangles ;
      shouldRemesh = checkSquareDomainViolation(elements);
      std::cout << "shouldRemesh: " << shouldRemesh << std::endl;

      mesh_iteration++;
      // Calculate post-remeshing energy and stress
      if (!shouldRemesh) {
        Eigen::Matrix2d post_stress_tensor = Eigen::Matrix2d::Zero();
        ConfigurationSaver::calculateEnergyAndStress(&userData, post_energy,
                                                     post_stress_tensor);
        post_stress = post_stress_tensor(
            0, 1); // Extract xy component for backward compatibility
        std::cout << "Post-remeshing - Energy: " << post_energy
                  << ", Stress: " << post_stress << std::endl;
      }
    }

    // Update plasticity flag for logging
    bool updated_plasticity = plasticity;

    // Log data to file
    ConfigurationSaver::logEnergyAndStress(
        i, alpha, pre_energy, pre_stress, post_energy, post_stress, hasChanges);
    hasChanges = 11;
    // Save configuration periodically
    if (hasChanges > 10 || i % 1000 == 0) {
      UserData finalUserData(square_points, elements, calculator,
                             potential_func, potential_func_der, zero,
                             optimal_lattice_parameter, F_ext, interior_mapping,
                             full_mapping, active_elements, plasticity);

      post_energy = 0;
      post_stress = 0;
      ConfigurationSaver::saveConfigurationWithStressAndEnergy2D(
          &finalUserData, i, post_energy, post_stress);
      ConfigurationSaver::writeToVTK(finalUserData.points,
                                     finalUserData.elements, &finalUserData, i);
    }

    exit(0);

    std::cout << "Stress tensor [MPa]:" << std::endl;
    std::cout << "  σxx = " << post_stress_tensor(0, 0) << std::endl;
    std::cout << "  σxy = " << post_stress_tensor(0, 1) << std::endl;
    std::cout << "  σyx = " << post_stress_tensor(1, 0) << std::endl;
    std::cout << "  σyy = " << post_stress_tensor(1, 1) << std::endl;

    // F_ext = calculateDeformationGradientEigen(post_stress_tensor);
    // Create calculator with default stiffness
    DeformationCalculator calculator;
    // DeformationCalculator calculator(customStiffness);

    // Calculate deformation gradient
    Eigen::Matrix2d cauchy_applied; /* your cauchy_applied matrix */
    ;

    // Parameters to control stress reduction
    double initialStress = 0.0; // Starting stress level

    int stressReductionInterval = 20;     // Reduce stress every X iterations
    double stressReductionFactor = 0.000; // How much to reduce stress by at
                                          // each step (0.5 = 50% reduction)

    // Calculate current stress level based on iteration count
    // if one wants stress loading
    int stressStage = i / stressReductionInterval;
    double stressLevel = initialStress;

    // Apply reduction for each completed stage
    for (int stage = 0; stage < stressStage; stage++) {
      // stressLevel *= (1.0 - stressReductionFactor);
      stressLevel -= stressReductionFactor;
    }

    // if(stressLevel >= 0 || stressStage<1 ){
    if (stressStage < 1) {

      std::cout << "Stress level: " << stressLevel << std::endl;
      std::cout << "Stress stage: " << stressStage << std::endl;

      // Apply deviatoric stress (trace = 0)
      cauchy_applied << .3, -0.00, -0.00, -.3;

      F_ext = calculator.calculateDeformationGradient(post_stress_tensor,
                                                      cauchy_applied);

    } else {
      std::vector<Point2D> square_points_back = square_points;

      RowResults registery =
          processPointRows(square_points_back, lattice_constant);
      FitResult result_fit =
          fitNabarroModel(registery.xmid, registery.deltax, lattice_constant);
      std::cout << "x_d_fit: " << result_fit.x_d
                << ", eta_fit: " << result_fit.eta << std::endl;
      std::string filename = "position_eta_stress.dat";

      savePositionEtaStress(filename, result_fit.x_d, result_fit.eta,
                            post_stress_tensor(0, 1));

      std::cout << "alpha level: " << alpha << std::endl;

      alpha_peirls = 0.00001;
      ;
      alpha_peirls_total += alpha_peirls;
      F_ext << 1, alpha_peirls, 0.0, 1;
      alpha_values[i + 1] = alpha_peirls_total;
    }

    std::cout << "Iteration " << i << " completed successfully" << std::endl;
  }
}
std::vector<int> find_boundary_fixed_nodes(const std::vector<Point2D> &points,
                                           double boundary_distance = 0.678,
                                           double tolerance = 1e-6) {
  std::vector<int> boundary_fixed_nodes;

  // Find the coordinate bounds
  double min_x = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();

  for (const auto &p : points) {
    min_x = std::min(min_x, p.coord.x());
    max_x = std::max(max_x, p.coord.x());
    min_y = std::min(min_y, p.coord.y());
    max_y = std::max(max_y, p.coord.y());
  }

  // Identify boundary nodes
  for (size_t i = 0; i < points.size(); i++) {
    const auto &p = points[i];

    bool is_left_boundary = std::abs(p.coord.x() - min_x) <= boundary_distance;
    bool is_right_boundary = std::abs(p.coord.x() - max_x) <= boundary_distance;
    bool is_bottom_boundary = std::abs(p.coord.y() - min_y) < boundary_distance;
    bool is_boundary =
        is_left_boundary || is_right_boundary || is_bottom_boundary;
    if (is_boundary) {
      boundary_fixed_nodes.push_back(i);
    }
  }

  return boundary_fixed_nodes;
}

void indentation() {
  // Parameters for lattice
  std::string lattice_type = "triangular"; // "square" or "triangular"

  // Energy functions (dumb)
  std::function<double(double)> potential_func = lennard_jones_energy_v2;
  std::function<double(double)> potential_func_der =
      lennard_jones_energy_der_v2;
  auto compute_even_ny = [](int nx) {
    int ny = std::round(2.0 * nx / std::sqrt(3));
    return (ny % 2 == 0) ? ny : ny + 1; // Ensure ny is even
  };

  int nx = 120;

  int ny = compute_even_ny(nx) / 2.;

  std::cout << "STEP 1: Finding optimal lattice parameter...\n";
  double h = 1.;
  double optimal_lattice_parameter = 0.6872044091828517;
  double lattice_constant = 0.6872044091828517;
  Eigen::Vector2d p1(0, 0);
  Eigen::Vector2d p2(h * lattice_constant, 0);
  Eigen::Vector2d p3(h * 0.5 * lattice_constant,
                     h * sqrt(3) * lattice_constant / 2);
  Eigen::Matrix<double, 3, 2> dndx = calculateShapeDerivatives(p1, p2, p3);
  std::cout << dndx << std::endl;

  // Using the lambda function:
  std::string filename_ref = "/Users/usalman/programming/FEM_2D/factorized/"
                             "lattice_triangulation/single/R30_ref.dat";
  std::string filename = "/Users/usalman/programming/FEM_2D/factorized/"
                         "lattice_triangulation/single/R30.dat";

  // std::vector<Point2D> square_points_ref = readPositionsFromFile(filename);
  // std::vector<Point2D> square_points = readPositionsFromFile(filename_ref);

  // Generate initial lattice
  std::vector<Point2D> square_points =
      LatticeGenerator::generate_2d_lattice_rotated(nx, ny, lattice_constant,
                                                    lattice_type, 90);

  // Generate initial lattice
  std::vector<Point2D> square_points_ref =
      LatticeGenerator::generate_2d_lattice_rotated(nx, ny, lattice_constant,
                                                    lattice_type, 90);

  // Create the indenter with a placeholder position (will be updated)
  double indenter_radius = 50.0;
  // NanoIndenter indenter_init(indenter_radius, Eigen::Vector2d(0.0, 0.0));

  // // Print initial state
  // std::cout << "Initial indenter position: ("
  //         << indenter.getPosition().x() << ", "
  //         << indenter.getPosition().y() << ")" << std::endl;

  // // Initialize indenter position based on the uppermost atoms
  // indenter_init.initializePosition(square_points);

  // // Print the initialized position
  // std::cout << "Initialized indenter position: ("
  //         << indenter.getPosition().x() << ", "
  //         << indenter.getPosition().y() << ")" << std::endl;
  // std::cout << "The indenter is positioned to touch the uppermost atom(s)" <<
  // std::endl; std::cout << "Bottom-most point of indenter: " <<
  // indenter.getBottomY() << std::endl;

  // // Set initial indentation depth to 0 to trigger contact detection
  // indenter_init.setIndentationDepth(-1, square_points);
  // const std::vector<int>& contacts = indenter_init.getFixedAtoms();
  // std::cout << "Number of atoms in contact (using class): " <<
  // contacts.size() << std::endl;

  // // Print positions of atoms in initial contact
  // std::cout << "Positions of atoms in initial contact:" << std::endl;
  // for (int idx : contacts) {
  //     std::cout << "Atom " << idx << ": ("
  //             << square_points[idx].coord.x() << ", "
  //             << square_points[idx].coord.y() << ")" << std::endl;

  //     // Also print the distance from this atom to the indenter center
  //     double distance = (square_points[idx].coord -
  //     indenter.getPosition()).norm(); std::cout << "  Distance to indenter
  //     center: " << distance
  //             << " (radius: " << indenter_radius << ")" << std::endl;
  // }

  // // Set the new indentation depth
  // indenter_init.setIndentationDepth(-1, square_points);

  // // Check for atoms in contact
  // const std::vector<int>& newcontacts = indenter_init.getFixedAtoms();
  // std::cout << "Number of atoms in contact: " << newcontacts.size() <<
  // std::endl;

  // // Print positions of atoms in contact after indentation
  // std::cout << "Positions of atoms in contact after indentation:" <<
  // std::endl;
  // for (int idx : newcontacts) {
  //     std::cout << "Atom " << idx << ": ("
  //             << square_points[idx].coord.x() << ", "
  //             << square_points[idx].coord.y() << ")" << std::endl;

  //     // Project the atom onto the indenter surface and show where it would
  //     be positioned Eigen::Vector2d projected_pos =
  //     indenter.projectOntoSurface(square_points[idx].coord); std::cout << "
  //     Projected position: ("
  //             << projected_pos.x() << ", "
  //             << projected_pos.y() << ")" << std::endl;

  //     // Calculate the movement required
  //     double movement = (projected_pos - square_points[idx].coord).norm();
  //     std::cout << "  Required movement: " << movement << std::endl;
  // }

  // // Print the indenter's new position
  // std::cout << "Current indenter position: ("
  //         << indenter.getPosition().x() << ", "
  //         << indenter.getPosition().y() << ")" << std::endl;
  // std::cout << "Current bottom-most point: " << indenter.getBottomY() <<
  // std::endl; exit(0);

  // RowResults registery = processPointRows(square_points, lattice_constant);
  // RowResults registery_diso = processPointRows(square_points_ref,
  // lattice_constant);

  // 2. Find the middle atom index
  size_t middle_atom_index =
      findMiddleAtom(square_points, true); // true enables verbose output
  // RowResults registery = processPointRows(square_points, lattice_constant);

  // std::vector<Point2D> square_points_ref =
  // readPositionsFromFile(filename_ref);

  // Now you can use the points vector
  std::cout << "Read " << square_points.size() << " points from file."
            << std::endl;

  std::cout << "Read " << square_points_ref.size() << " points from file."
            << std::endl;

  int original_domain_size = square_points.size();

  DomainInfo domain_size = compute_domain_size(square_points);

  const std::array<double, 2> offsets = {
      h * lattice_constant, h * (sqrt(3.) / 2.) * lattice_constant};
  std::cout << "offsets: " << offsets[0] << " " << offsets[1] << std::endl;
  DomainDimensions domain_dims(domain_size.get_width(),
                               domain_size.get_height());
  std::cout << "domain_size.get_width(): " << domain_size.get_width()
            << std::endl;
  std::cout << "domain_size.get_height(): " << domain_size.get_height()
            << std::endl;
  bool pbc = true;
  // Update DOF mapping for pbc
  auto [interior_mapping, full_mapping] = create_dof_mapping_original(
      square_points, 0.001, pbc); // Create domain maps

  auto [original_domain_map, translation_map] =
      MeshGenerator::create_domain_maps(original_domain_size, domain_dims,
                                        offsets);

  // Create ALGLIB array for free DOFs (displacements)

  // Setup triangulation variables
  // Create a Point2D from your DomainDimensions locally
  Point2D domain_dims_point(domain_dims.size_x, domain_dims.size_y);

  // 1. Create the AdaptiveMesher instance
  AdaptiveMesher mesher(domain_dims_point, offsets, original_domain_map,
                        translation_map, full_mapping,
                        1e-6, // Tolerance,
                        pbc);
  mesher.setUsePeriodicCopies(pbc); // Switch to using original domain only
  alglib::real_1d_array free_dofs;
  int n_free_nodes = interior_mapping.size();
  free_dofs.setlength(2 *
                      interior_mapping.size()); // [u0, u1, ..., v0, v1, ...]
  map_points_to_solver_array(free_dofs, square_points, interior_mapping,
                             n_free_nodes);

  alglib::real_1d_array original_x_remesh =
      mesher.saveOriginalPositions(free_dofs);

  std::vector<int> boundary_fixed_nodes =
      find_boundary_fixed_nodes(square_points, lattice_constant);

  auto [elements, active_elements] = mesher.createMesh(
      square_points, free_dofs, Eigen::Matrix2d::Identity(), &dndx);
  double element_area = elements[0].getReferenceArea();

  std::cout << "Created " << elements.size() << " element triangles"
            << std::endl;

  // Setup energy calculation
  TriangularLatticeCalculator calculator(lattice_constant);

  double zero = calculator.calculate_energy(Eigen::Matrix2d::Identity(),
                                            potential_func, 0);
  std::cout << "debugging simple shear test" << std::endl;
  std::cout << "zero energy value: " << zero << std::endl;

  // Prepare for optimization
  alglib::real_1d_array x;
  int n_vars = interior_mapping.size();
  x.setlength(2 * n_vars);

  std::cout << "active_elements.size(): " << active_elements.size()
            << std::endl;
  std::vector<ElementTriangle2D> &elements_ref =
      elements; // Reference to elements
  std::vector<size_t> &active_elements_ref =
      active_elements; // Reference to active elements

  // Define indenter parameters
  NanoIndenter indenter(indenter_radius, Eigen::Vector2d(0.0, 0.0));

  // Initialize indenter position based on the uppermost atoms
  indenter.initializePosition(square_points);

  // Define a series of increasing indentation depths
  std::vector<double> indentation_depths;
  for (double depth = 0.0; depth <= 5.0; depth += 0.05) {
    indentation_depths.push_back(depth);
  }
  // Process each indentation depth
  for (size_t i = 0; i < indentation_depths.size(); i++) {
    double depth = indentation_depths[i];
    std::cout << "\n=== Processing indentation depth = " << depth
              << " ===" << std::endl;

    // Set the indentation depth and identify atoms in contact
    indenter.setIndentationDepth(depth, square_points);

    const std::vector<int> &contact_atoms = indenter.getFixedAtoms();

    std::cout << "Total contact atoms: " << contact_atoms.size() << std::endl;
    std::cout << "Contact Atoms Details:" << std::endl;
    // for (size_t i = 0; i < contact_atoms.size(); ++i) {
    //     std::cout << "Index " << i << ": " << contact_atoms[i] << std::endl;
    // }
    std::cout << "Number of atoms in contact: " << contact_atoms.size()
              << std::endl;

    // Assuming `contact_atoms` is already populated with the indices of atoms
    // in contact with the indenter
    std::vector<int> neighboring_atoms = indenter.getNeighboringAtomIndices(
        square_points, contact_atoms, 1.1 * sqrt(3) * lattice_constant / 2,
        1.1 * lattice_constant);

    std::set<int> unique_atoms_set;

    // Insert contact atoms and neighboring atoms into the set
    unique_atoms_set.insert(contact_atoms.begin(), contact_atoms.end());
    unique_atoms_set.insert(neighboring_atoms.begin(), neighboring_atoms.end());

    // Convert the set back to a vector
    std::vector<int> combined_atoms(unique_atoms_set.begin(),
                                    unique_atoms_set.end());

    // Project atoms in contact onto the indenter surface
    for (int idx : contact_atoms) {
      square_points[idx].coord =
          indenter.projectOntoSurface(square_points[idx].coord);
    }
    for (int idx : neighboring_atoms) {
      square_points[idx].coord =
          indenter.projectOntoSurface_NN(square_points[idx].coord);
    }

    // Update domain mapping taking into account fixed atoms from indenter
    auto [interior_mapping, full_mapping] = create_dof_mapping_with_boundaries(
        square_points, elements, combined_atoms, boundary_fixed_nodes);
    // auto [interior_mapping, full_mapping] = create_dof_mapping_original(
    //     square_points,0.001,1);    // Create domain maps

    // Update n_vars based on new interior_mapping size
    n_vars = interior_mapping.size(); // 2D coordinates

    // Create deformation gradient (identity for indentation simulation)
    static Eigen::Matrix2d F_ext = Eigen::Matrix2d::Identity();
    x.setlength(2 * n_vars);

    // Prepare for optimization
    map_points_to_solver_array(x, square_points, interior_mapping, n_vars);

    for (auto &element : elements) {
      element.set_reference_mesh(square_points);
      element.set_dof_mapping(
          full_mapping); // or interior_mapping depending on needs
      // double jac = element.calculate_shape_derivatives(x);  // current
      // positions
    }
    const std::vector<size_t> new_active_elements = initialize_active_elements(
        elements, full_mapping, square_points.size());
    active_elements = new_active_elements;

    // Create user data
    bool plasticity = false;
    UserData userData(square_points, elements, calculator, potential_func,
                      potential_func_der, zero, optimal_lattice_parameter,
                      F_ext, interior_mapping, full_mapping, active_elements,
                      plasticity);
    double post_energyp = 0.0;
    double post_stressp = 0.0;
    if (i == 0) {
      ConfigurationSaver::saveConfigurationWithStressAndEnergy2D(
          &userData, 999, post_energyp, post_stressp);
      ConfigurationSaver::writeToVTK(userData.points, userData.elements,
                                     &userData, 999);
    }

    // Calculate pre-optimization energy and stress
    double pre_energy = 0.0;
    double pre_stress = 0.0;
    Eigen::Matrix2d pre_stress_tensor = Eigen::Matrix2d::Zero();
    ConfigurationSaver::calculateEnergyAndStress(&userData, pre_energy,
                                                 pre_stress_tensor);
    pre_stress = pre_stress_tensor(
        0, 1); // Extract xy component for backward compatibility std::cout <<
               // "Pre-optimization - Energy: " << pre_energy << ", Stress: " <<
               // pre_stress << std::endl;

    // Store original positions
    alglib::real_1d_array original_x;
    original_x.setlength(x.length());
    for (int j = 0; j < x.length(); j++) {
      original_x[j] = x[j];
    }

    // Run optimization
    // --- Start timing ---
    auto wall_start = std::chrono::high_resolution_clock::now();
    clock_t cpu_start = clock();

    // Run optimization
    LBFGSOptimizer optimizer(10, 0, 0, 0, 0);
    // CGOptimizer optimizer;
    optimizer.optimize(x, minimize_energy_with_triangles_noreduction,
                       &userData);

    // --- Stop timing ---
    auto wall_end = std::chrono::high_resolution_clock::now();
    clock_t cpu_end = clock();

    // Compute durations
    double wall_time =
        std::chrono::duration<double>(wall_end - wall_start).count();
    double cpu_time = (double)(cpu_end - cpu_start) / CLOCKS_PER_SEC;

    // Print results
    std::cout << "Optimization wall-clock time: " << wall_time << " seconds\n";
    std::cout << "Optimization CPU time: " << cpu_time << " seconds\n";
    std::cout << "Optimization Ratio: " << cpu_time / wall_time << " seconds\n";

    map_solver_array_to_points(x, square_points, interior_mapping, n_vars);

    if (i == 0) {
      ConfigurationSaver::saveConfigurationWithStressAndEnergy2D(
          &userData, 9990, post_energyp, post_stressp);
      ConfigurationSaver::writeToVTK(userData.points, userData.elements,
                                     &userData, 9009);
    }

    // Calculate change measures and check for remeshing need
    ChangeMeasures result =
        computeChangeMeasures(x, original_x, lattice_constant, elements,
                              &userData, square_points, true, &F_ext);
    std::cout << "max_abs_change: " << result.max_abs_change << std::endl;
    if (result.has_distorted_triangles) {
      std::cout << "Distorted triangles detected!" << std::endl;
    }

    // Determine if remeshing is needed
    bool shouldRemesh = result.has_distorted_triangles;
    std::cout << "shouldRemesh: " << shouldRemesh << std::endl;

    double post_energy = 0.0;
    double post_stress = 0.0;
    Eigen::Matrix2d post_stress_tensor = Eigen::Matrix2d::Zero();
    // 2. Call the function
    if (shouldRemesh) {

      // Calling the function
      int max_iterations = 10;
      auto [post_energy_re, post_stress_tensor_re, iterations] =
          perform_remeshing_loop(
              x,
              &userData, // Note: removed & if userData is already a pointer
              contact_atoms, boundary_fixed_nodes, F_ext, dndx, offsets,
              original_domain_map, translation_map, domain_dims_point,
              max_iterations, element_area);
      post_energy = post_energy_re;
      post_stress_tensor = post_stress_tensor_re; // Extract xy component for
                                                  // backward compatibility
      post_stress = post_stress_tensor(0, 1);
    }

    else {
      // Calculate post-optimization energy and stress

      ConfigurationSaver::calculateEnergyAndStress(&userData, post_energy,
                                                   post_stress_tensor);
      post_stress = post_stress_tensor(
          0, 1); // Extract xy component for backward compatibility
    }

    double hasChanges = 11;

    // Update plasticity flag for logging
    bool updated_plasticity = 1;

    // Log data to file
    ConfigurationSaver::logEnergyAndStress(
        i, depth, pre_energy, pre_stress, post_energy, post_stress, hasChanges);
    // Save configuration periodically
    if (hasChanges > 10 || i % 1000 == 0) {
      UserData finalUserData(square_points, elements, calculator,
                             potential_func, potential_func_der, zero,
                             optimal_lattice_parameter, F_ext, interior_mapping,
                             full_mapping, active_elements, plasticity);

      post_energy = 0;
      post_stress = 0;
      ConfigurationSaver::saveConfigurationWithStressAndEnergy2D(
          &finalUserData, i, post_energy, post_stress);
      ConfigurationSaver::writeToVTK(finalUserData.points,
                                     finalUserData.elements, &finalUserData, i);
    }

    std::cout << "Iteration " << i << " completed successfully" << std::endl;
  }
}
