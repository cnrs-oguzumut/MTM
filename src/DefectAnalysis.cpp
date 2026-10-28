#include "../include/defects/Defectanalysis.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <tuple>

namespace {

struct EdgeKey {
  int n_min;
  int n_max;
  long long dtx;
  long long dty;

  bool operator<(const EdgeKey &other) const {
    return std::tie(n_min, n_max, dtx, dty) <
           std::tie(other.n_min, other.n_max, other.dtx, other.dty);
  }
};

Eigen::Vector2i computeCanonicalLatticeStep(int n_a, const Eigen::Vector2d &x_a,
                                            int n_b, const Eigen::Vector2d &x_b,
                                            const Eigen::Matrix2d &A0_inv) {
  bool forward = false;
  if (n_a != n_b) {
    forward = (n_a < n_b);
  } else if (std::abs(x_a.x() - x_b.x()) > 1e-12) {
    forward = (x_a.x() < x_b.x());
  } else {
    forward = (x_a.y() < x_b.y());
  }

  Eigen::Vector2d d = forward ? (x_b - x_a) : (x_a - x_b);
  Eigen::Vector2d s = A0_inv * d;
  Eigen::Vector2i step(static_cast<int>(std::llround(s.x())),
                       static_cast<int>(std::llround(s.y())));
  return forward ? step : (-step);
}

} // namespace

DefectAnalysis::BurgersCircuitResult
DefectAnalysis::computeBurgersCircuits(const UserData *userData) {
  BurgersCircuitResult result;
  if (!userData) {
    return result;
  }

  const auto &points = userData->points;
  const auto &elements = userData->elements;
  const auto &active_elements = userData->active_elements;
  const size_t num_points = points.size();

  result.nodal_burgers_int.assign(num_points, Eigen::Vector2i::Zero());
  result.nodal_burgers_vec.assign(num_points, Eigen::Vector2d::Zero());
  result.nodal_burgers_mag.assign(num_points, 0.0);
  result.filtered_coordination.assign(num_points, 6);

  if (num_points == 0 || active_elements.empty()) {
    return result;
  }

  Eigen::Matrix2d F_ext_inv = Eigen::Matrix2d::Identity();
  if (std::abs(userData->F_external.determinant()) > 1e-12) {
    F_ext_inv = userData->F_external.inverse();
  }

  struct TriInfo {
    int nodes[3];
    Eigen::Vector2d trans[3];
    Eigen::Vector2d pos[3];
    double lengths[3]; // sorted L1 <= L2 <= L3
    double edge_len[3];
  };

  std::vector<TriInfo> valid_tris;
  valid_tris.reserve(active_elements.size());

  std::vector<double> shortest_edges;
  std::vector<double> aspect_ratios;
  shortest_edges.reserve(active_elements.size());
  aspect_ratios.reserve(active_elements.size());

  for (size_t elem_idx : active_elements) {
    if (elem_idx >= elements.size())
      continue;
    const auto &element = elements[elem_idx];
    if (!element.isInitialized())
      continue;

    TriInfo tri;
    bool valid_indices = true;
    for (int k = 0; k < 3; ++k) {
      tri.nodes[k] = element.getNodeIndex(k);
      if (tri.nodes[k] < 0 ||
          tri.nodes[k] >= static_cast<int>(num_points)) {
        valid_indices = false;
        break;
      }
      tri.trans[k] = element.getTranslation(k);
      tri.pos[k] = F_ext_inv * points[tri.nodes[k]].coord + tri.trans[k];
    }
    if (!valid_indices)
      continue;

    // Ensure counter-clockwise orientation in reference frame
    Eigen::Vector2d e01 = tri.pos[1] - tri.pos[0];
    Eigen::Vector2d e02 = tri.pos[2] - tri.pos[0];
    double cross = e01.x() * e02.y() - e01.y() * e02.x();
    if (cross < 0.0) {
      std::swap(tri.nodes[1], tri.nodes[2]);
      std::swap(tri.trans[1], tri.trans[2]);
      std::swap(tri.pos[1], tri.pos[2]);
    }

    tri.edge_len[0] = (tri.pos[1] - tri.pos[0]).norm();
    tri.edge_len[1] = (tri.pos[2] - tri.pos[1]).norm();
    tri.edge_len[2] = (tri.pos[0] - tri.pos[2]).norm();

    tri.lengths[0] = tri.edge_len[0];
    tri.lengths[1] = tri.edge_len[1];
    tri.lengths[2] = tri.edge_len[2];
    std::sort(tri.lengths, tri.lengths + 3);

    if (tri.lengths[0] > 1e-8) {
      shortest_edges.push_back(tri.lengths[0]);
      aspect_ratios.push_back(tri.lengths[2] / tri.lengths[0]);
    }

    valid_tris.push_back(tri);
  }

  const size_t num_valid_tris = valid_tris.size();
  result.element_burgers_int.assign(num_valid_tris, Eigen::Vector2i::Zero());
  result.element_burgers_vec.assign(num_valid_tris, Eigen::Vector2d::Zero());
  result.element_burgers_mag.assign(num_valid_tris, 0.0);

  if (num_valid_tris == 0 || shortest_edges.empty()) {
    return result;
  }

  std::sort(shortest_edges.begin(), shortest_edges.end());
  std::sort(aspect_ratios.begin(), aspect_ratios.end());
  double median_L1 = shortest_edges[shortest_edges.size() / 2];
  double median_ratio = aspect_ratios[aspect_ratios.size() / 2];

  double a0 = median_L1;
  if (userData->ideal_lattice_parameter > 1e-6 &&
      std::abs(userData->ideal_lattice_parameter - median_L1) / median_L1 <
          0.25) {
    a0 = userData->ideal_lattice_parameter;
  }

  bool is_square = (median_ratio > 1.20);
  int sym_p = is_square ? 4 : 6;

  // Estimate global crystal orientation angle theta0 from nearest-neighbor bonds
  double cos_sum = 0.0;
  double sin_sum = 0.0;
  for (const auto &tri : valid_tris) {
    for (int k = 0; k < 3; ++k) {
      double len = tri.edge_len[k];
      if (len >= 0.75 * a0 && len <= 1.22 * a0) {
        Eigen::Vector2d d = tri.pos[(k + 1) % 3] - tri.pos[k];
        double phi = std::atan2(d.y(), d.x());
        cos_sum += std::cos(sym_p * phi);
        sin_sum += std::sin(sym_p * phi);
      }
    }
  }
  double theta0 = (cos_sum * cos_sum + sin_sum * sin_sum > 1e-12)
                      ? (std::atan2(sin_sum, cos_sum) / sym_p)
                      : 0.0;
  if (std::abs(theta0) < 1e-4) {
    theta0 = 0.0;
  }

  Eigen::Matrix2d R0;
  R0 << std::cos(theta0), -std::sin(theta0),
        std::sin(theta0),  std::cos(theta0);

  Eigen::Matrix2d basis = Eigen::Matrix2d::Identity();
  if (!is_square) {
    basis << 1.0, 0.5,
             0.0, std::sqrt(3.0) / 2.0;
  }
  Eigen::Matrix2d A0 = a0 * R0 * basis;
  Eigen::Matrix2d A0_inv = A0.inverse();

  result.is_square_lattice = is_square;
  result.lattice_spacing = a0;
  result.lattice_angle_rad = theta0;
  result.A0 = A0;

  // Compute elementary triangle Burgers circuits B_T = dn_01 + dn_12 + dn_20
  std::map<EdgeKey, std::vector<std::pair<size_t, double>>> edge_to_tris;
  for (size_t t_idx = 0; t_idx < num_valid_tris; ++t_idx) {
    const auto &tri = valid_tris[t_idx];
    // Ignore non-physical outer boundary sliver triangles
    if (tri.lengths[2] > 2.15 * a0) {
      continue;
    }

    Eigen::Vector2i B_T = Eigen::Vector2i::Zero();
    for (int k = 0; k < 3; ++k) {
      int k_next = (k + 1) % 3;
      B_T += computeCanonicalLatticeStep(tri.nodes[k], tri.pos[k],
                                         tri.nodes[k_next], tri.pos[k_next],
                                         A0_inv);

      int n_a = tri.nodes[k];
      int n_b = tri.nodes[k_next];
      Eigen::Vector2d dt =
          (n_a <= n_b) ? (tri.trans[k_next] - tri.trans[k])
                       : (tri.trans[k] - tri.trans[k_next]);
      EdgeKey key{std::min(n_a, n_b), std::max(n_a, n_b),
                   static_cast<long long>(std::llround(dt.x() * 1e4)),
                  static_cast<long long>(std::llround(dt.y() * 1e4))};
      edge_to_tris[key].push_back({t_idx, tri.edge_len[k]});
    }
    result.element_burgers_int[t_idx] = B_T;
  }

  // Unit-cell plaquette cancellation across shared diagonal edges:
  // If two triangles sharing a diagonal edge (> 1.15 * a0) have B_T1 + B_T2 == 0,
  // their 4-atom unit cell has zero net Burgers vector (only the internal diagonal
  // crossed a rounding threshold).
  for (const auto &kv : edge_to_tris) {
    const auto &adj = kv.second;
    if (adj.size() == 2) {
      size_t t1 = adj[0].first;
      size_t t2 = adj[1].first;
      double shared_len = 0.5 * (adj[0].second + adj[1].second);
      if (shared_len > 1.15 * a0 &&
          result.element_burgers_int[t1] != Eigen::Vector2i::Zero() &&
          (result.element_burgers_int[t1] + result.element_burgers_int[t2]) ==
              Eigen::Vector2i::Zero()) {
        result.element_burgers_int[t1].setZero();
        result.element_burgers_int[t2].setZero();
      }
    }
  }

  // Accumulate element physical Burgers vectors and 1-ring nodal Burgers vectors
  int total_dislocations = 0;
  for (size_t t_idx = 0; t_idx < num_valid_tris; ++t_idx) {
    const Eigen::Vector2i &B_T = result.element_burgers_int[t_idx];
    if (B_T != Eigen::Vector2i::Zero()) {
      Eigen::Vector2d b_vec = (A0 * B_T.cast<double>()) / a0;
      double b_mag = b_vec.norm();
      result.element_burgers_vec[t_idx] = b_vec;
      result.element_burgers_mag[t_idx] = b_mag;
      total_dislocations += std::max(1, static_cast<int>(std::llround(b_mag)));

      const auto &tri = valid_tris[t_idx];
      for (int k = 0; k < 3; ++k) {
        int node_id = tri.nodes[k];
        result.nodal_burgers_int[node_id] += B_T;
      }
    }
  }
  result.num_dislocations = total_dislocations;

  std::vector<int> raw_coord =
      calculateCoordinationNumbers(points, elements, active_elements);

  for (size_t i = 0; i < num_points; ++i) {
    const Eigen::Vector2i &B_i = result.nodal_burgers_int[i];
    if (B_i != Eigen::Vector2i::Zero()) {
      Eigen::Vector2d b_vec = (A0 * B_i.cast<double>()) / a0;
      result.nodal_burgers_vec[i] = b_vec;
      result.nodal_burgers_mag[i] = b_vec.norm();
      if (i < raw_coord.size() && (raw_coord[i] == 5 || raw_coord[i] == 7)) {
        result.filtered_coordination[i] = raw_coord[i];
      } else {
        result.filtered_coordination[i] =
            (B_i.x() + B_i.y() >= 0) ? 7 : 5;
      }
    } else {
      result.filtered_coordination[i] = 6;
    }
  }

  return result;
}

std::pair<int, std::vector<int>>
DefectAnalysis::analyzeDefectsInReferenceConfig(
    const UserData *userData, int file_id,
    const Eigen::Matrix<double, 3, 2> &dndx,
    const std::array<double, 2> &offsets,
    const std::vector<int> &original_domain_map,
    const std::vector<std::tuple<double, double>> &translation_map,
    const Point2D &domain_dims_point, double element_area, bool pbc,
    bool reduction) {
  const bool verbose_log = (std::getenv("MTM_VERBOSE_LOG") != nullptr);
  if (verbose_log) {
    std::cout << "\n=== Analyzing defects in reference configuration (Method 2: Discrete Burgers Circuit) ==="
              << std::endl;
  }

  // 1. COPY the current points (don't modify originals!)
  std::vector<Point2D> reference_points = userData->points;

  // 2. Transform to reference configuration
  Eigen::Matrix2d F_ext_inv = userData->F_external.inverse();
  transformToReferenceConfig(reference_points, F_ext_inv);

  if (verbose_log) {
    std::cout << "Transformed " << reference_points.size()
              << " points to reference configuration" << std::endl;
  }

  // 3. Map to solver array
  alglib::real_1d_array x_ref;
  int n_vars = userData->interior_mapping.size();
  x_ref.setlength(2 * n_vars);
  map_points_to_solver_array(x_ref, reference_points,
                             userData->interior_mapping, n_vars);

  // 4. REMESH in reference configuration
  AdaptiveMesher mesher(domain_dims_point, offsets, original_domain_map,
                        translation_map, userData->full_mapping, 1e-6, pbc);

  Eigen::Matrix2d F_identity = Eigen::Matrix2d::Identity();
  auto [ref_elements, ref_active_elements] = mesher.createMesh(
      reference_points, x_ref,
      F_identity, // Identity because we're in reference config
      &dndx);

  if (verbose_log) {
    std::cout << "Created " << ref_elements.size()
              << " elements in reference configuration" << std::endl;
  }

  // 5. Create temporary UserData for reference config
  UserData refUserData(reference_points, ref_elements, userData->calculator,
                       userData->energy_function, userData->derivative_function,
                       userData->zero_energy, userData->ideal_lattice_parameter,
                       F_identity, // Identity in reference config
                       userData->interior_mapping, userData->full_mapping,
                       ref_active_elements,
                       false // plasticity flag
  );

  // 6. Compute Method 2 Discrete Burgers Circuits on reference mesh
  BurgersCircuitResult burgers = computeBurgersCircuits(&refUserData);

  // 7. Write VTK with Burgers circuit fields and filtered coordination numbers
  ConfigurationSaver::writeToVTK_DefectAnalysis(reference_points, ref_elements,
                                                &refUserData, file_id);

  if (verbose_log) {
    std::cout << "Saved reference configuration VTK to file " << file_id
              << std::endl;
    std::cout << "Found " << burgers.num_dislocations
              << " dislocations via Method 2 (Discrete Burgers Circuit, lattice="
              << (burgers.is_square_lattice ? "square" : "triangular")
              << ", a0=" << burgers.lattice_spacing << ")" << std::endl;
    std::cout << "=== Reference configuration analysis complete ===" << std::endl;
  }

  return std::make_pair(burgers.num_dislocations,
                        burgers.filtered_coordination);
}

void DefectAnalysis::transformToReferenceConfig(
    std::vector<Point2D> &points, const Eigen::Matrix2d &F_ext_inv) {
  for (auto &point : points) {
    Eigen::Vector2d pos(point.coord.x(), point.coord.y());
    Eigen::Vector2d transformed_pos = F_ext_inv * pos;
    point.coord.x() = transformed_pos(0);
    point.coord.y() = transformed_pos(1);
  }
}

std::vector<int> DefectAnalysis::calculateCoordinationNumbers(
    const std::vector<Point2D> &points,
    const std::vector<ElementTriangle2D> &elements,
    const std::vector<size_t> &active_elements) {
  std::vector<int> coordination(points.size(), 0);

  // Count how many elements each node belongs to
  for (size_t elem_idx : active_elements) {
    if (elem_idx >= elements.size())
      continue;

    const auto &element = elements[elem_idx];
    if (!element.isInitialized())
      continue;

    for (int i = 0; i < 3; i++) {
      int node_idx = element.getNodeIndex(i);
      if (node_idx < static_cast<int>(coordination.size())) {
        coordination[node_idx]++;
      }
    }
  }

  return coordination;
}

// int DefectAnalysis::count57Pairs(
//     const std::vector<Point2D>& points,
//     const std::vector<ElementTriangle2D>& elements,
//     const std::vector<size_t>& active_elements
// ) {
//     // Calculate coordination numbers
//     std::vector<int> coordination = calculateCoordinationNumbers(
//         points, elements, active_elements
//     );

//     // Build adjacency map (which nodes are connected through elements)
//     std::map<int, std::set<int>> adjacency;

//     for (size_t elem_idx : active_elements) {
//         if (elem_idx >= elements.size()) continue;

//         const auto& element = elements[elem_idx];
//         if (!element.isInitialized()) continue;

//         int nodes[3];
//         for (int i = 0; i < 3; i++) {
//             nodes[i] = element.getNodeIndex(i);
//         }

//         // Add edges between all pairs of nodes in this element
//         for (int i = 0; i < 3; i++) {
//             for (int j = i + 1; j < 3; j++) {
//                 adjacency[nodes[i]].insert(nodes[j]);
//                 adjacency[nodes[j]].insert(nodes[i]);
//             }
//         }
//     }

//     // Count 5/7 pairs (adjacent nodes with coordination 5 and 7)
//     int pair_count = 0;
//     std::set<std::pair<int, int>> counted_pairs;

//     for (size_t i = 0; i < points.size(); i++) {
//         if (coordination[i] == 5) {
//             // Check if any neighbor has coordination 7
//             for (int neighbor : adjacency[i]) {
//                 if (coordination[neighbor] == 7) {
//                     // Ensure we don't count the same pair twice
//                     int min_idx = std::min(static_cast<int>(i), neighbor);
//                     int max_idx = std::max(static_cast<int>(i), neighbor);
//                     std::pair<int, int> pair(min_idx, max_idx);

//                     if (counted_pairs.find(pair) == counted_pairs.end()) {
//                         counted_pairs.insert(pair);
//                         pair_count++;

//                         std::cout << "5/7 pair found: node " << i
//                                   << " (coord=" << coordination[i] << ") -
//                                   node "
//                                   << neighbor << " (coord=" <<
//                                   coordination[neighbor]
//                                   << ")" << std::endl;
//                     }
//                 }
//             }
//         }
//     }

//     return pair_count;
// }

// int DefectAnalysis::count57Pairs(
//     const std::vector<Point2D>& points,
//     const std::vector<ElementTriangle2D>& elements,
//     const std::vector<size_t>& active_elements
// ) {
//     // Calculate coordination numbers
//     std::vector<int> coordination = calculateCoordinationNumbers(
//         points, elements, active_elements
//     );

//     // Build adjacency map (which nodes are connected through elements)
//     std::map<int, std::set<int>> adjacency;

//     for (size_t elem_idx : active_elements) {
//         if (elem_idx >= elements.size()) continue;

//         const auto& element = elements[elem_idx];
//         if (!element.isInitialized()) continue;

//         int nodes[3];
//         for (int i = 0; i < 3; i++) {
//             nodes[i] = element.getNodeIndex(i);
//         }

//         // Add edges between all pairs of nodes in this element
//         for (int i = 0; i < 3; i++) {
//             for (int j = i + 1; j < 3; j++) {
//                 adjacency[nodes[i]].insert(nodes[j]);
//                 adjacency[nodes[j]].insert(nodes[i]);
//             }
//         }
//     }

//     // Count 5/7 pairs (adjacent nodes with coordination 5 and 7)
//     int pair_count = 0;
//     std::set<std::pair<int, int>> counted_pairs;

//     for (size_t i = 0; i < points.size(); i++) {
//         if (coordination[i] == 5) {
//             // Check if any neighbor has coordination 7
//             for (int neighbor : adjacency[i]) {
//                 if (coordination[neighbor] == 7) {
//                     // Ensure we don't count the same pair twice
//                     int min_idx = std::min(static_cast<int>(i), neighbor);
//                     int max_idx = std::max(static_cast<int>(i), neighbor);
//                     std::pair<int, int> pair(min_idx, max_idx);

//                     if (counted_pairs.find(pair) == counted_pairs.end()) {
//                         counted_pairs.insert(pair);
//                         pair_count++;

//                         std::cout << "5/7 pair found: node " << i
//                                   << " (coord=" << coordination[i] << ") -
//                                   node "
//                                   << neighbor << " (coord=" <<
//                                   coordination[neighbor]
//                                   << ")" << std::endl;
//                     }
//                 }
//             }
//         }
//     }

//     return pair_count;
// }

int DefectAnalysis::count57Pairs(const std::vector<Point2D> &points,
                                 const std::vector<ElementTriangle2D> &elements,
                                 const std::vector<size_t> &active_elements) {
  // Calculate coordination numbers
  std::vector<int> coordination =
      calculateCoordinationNumbers(points, elements, active_elements);

  // --- Adjacency map logic removed ---

  // Count all nodes with coordination 5 and 7
  int count_5 = 0;
  int count_7 = 0;

  // Iterate through all coordination numbers
  // (Assuming coordination.size() == points.size())
  for (int coord : coordination) {
    if (coord == 5) {
      count_5++;
    } else if (coord == 7) {
      count_7++;
    }
  }

  // Sum the counts as requested
  int total_defects = count_5 + count_7;

  // Return the total count divided by 2
  return total_defects / 2;
}