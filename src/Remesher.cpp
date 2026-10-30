#include "../include/mesh/Remesher.h"
#include "../include/geometry/DomainDimensions.h"
#include "../include/geometry/LatticeGenerator.h"
#include "../include/geometry/MacroscopicGeometry.h"
#include "../include/mesh/ElementTriangle2D.h"
#include "../include/mesh/MeshGenerator.h"
#include "../include/mesh/Triangle.h"

namespace {
const MacroscopicGeometryConfig* g_active_geometry_filter = nullptr;
const std::vector<Point2D>* g_active_geometry_ref_points = nullptr;
} // namespace

void AdaptiveMesher::setGlobalGeometryFilter(const MacroscopicGeometryConfig* geom,
                                             const std::vector<Point2D>* ref_points) {
  g_active_geometry_filter = geom;
  g_active_geometry_ref_points = ref_points;
}

/**
 * Constructor for AdaptiveMesher
 */
AdaptiveMesher::AdaptiveMesher(
    const Point2D &domain_dimensions,
    const std::array<double, 2> &periodic_offsets,
    const std::vector<int> &original_domain_mapping,
    const std::vector<std::tuple<double, double>> &translation_mapping,
    const std::vector<std::pair<int, int>> &dof_mapping, double tolerance,
    bool use_periodic_copies // Added parameter for periodic copies option
    )
    : domain_dims(domain_dimensions), offsets(periodic_offsets),
      original_domain_map(original_domain_mapping),
      translation_map(translation_mapping), full_mapping(dof_mapping),
      unique_triangle_tolerance(tolerance),
      use_periodic(use_periodic_copies), // Initialize the new member
      triangulation_perturbation_amplitude(0.0) {}

/**
 * Generate a new mesh based on the current configuration
 */
/**
 * Generate a new mesh based on the current configuration with optional manual
 * shape derivatives
 */
std::vector<size_t> AdaptiveMesher::generateMesh(
    const std::vector<Point2D> &square_points, const alglib::real_1d_array &x,
    const Eigen::Matrix2d &F_ext, std::vector<ElementTriangle2D> &elements,
    std::vector<Point2D> &points_used, std::vector<Triangle> &triangulation,
    std::vector<Triangle> &unique_triangles,
    const Eigen::Matrix<double, 3, 2> *manual_shape_derivative) {
  // Create a DomainDimensions object from the Point2D domain_dims
  DomainDimensions domainDimensions(domain_dims.coord.x(),
                                    domain_dims.coord.y());

  std::vector<Point2D> points_for_triangulation;

  if (use_periodic) {
    // Generate new periodic copies with current deformation
    points_for_triangulation = LatticeGenerator::create_periodic_copies(
        square_points, domainDimensions, offsets, F_ext);
  } else {
    // Use original domain without periodic copies
    points_for_triangulation = square_points;
  }

  if (triangulation_perturbation_amplitude != 0.0) {
    std::cout << "Applying triangulation-only shear perturbation: "
              << triangulation_perturbation_amplitude << std::endl;
    for (auto &point : points_for_triangulation) {
      point.coord.y() += triangulation_perturbation_amplitude * point.coord.x();
    }
  } else if (g_active_geometry_filter != nullptr &&
             g_active_geometry_filter->enabled &&
             g_active_geometry_filter->microstructure.mode == MicrostructureMode::Bicrystal) {
    // Mirror-symmetric Delaunay tie-breaker for bicrystals: breaks 4-node co-circular
    // unit-cell degeneracy identically in Grain 0 and Grain 1 while preserving exact
    // reflection symmetry across the grain boundary plane.
    const double cx_gb = g_active_geometry_filter->x_min +
                         g_active_geometry_filter->microstructure.gb_cx * g_active_geometry_filter->Lx;
    const double cy_gb = g_active_geometry_filter->y_min +
                         g_active_geometry_filter->microstructure.gb_cy * g_active_geometry_filter->Ly;
    const double phi_rad = g_active_geometry_filter->microstructure.gb_angle_deg * M_PI / 180.0;
    const double nx_gb = std::cos(phi_rad);
    const double ny_gb = std::sin(phi_rad);
    const double tx_gb = -ny_gb;
    const double ty_gb = nx_gb;
    const double sym_eps = -1e-6;
    for (auto &point : points_for_triangulation) {
      const double d_norm = std::abs((point.coord.x() - cx_gb) * nx_gb +
                                     (point.coord.y() - cy_gb) * ny_gb);
      point.coord.x() += sym_eps * d_norm * tx_gb;
      point.coord.y() += sym_eps * d_norm * ty_gb;
    }
  }

  // Create new triangulation
  triangulation =
      MeshGenerator::createTrianglesFromPoints(points_for_triangulation);
  points_used = points_for_triangulation;

  // Select triangles
  if (use_periodic) {
    // Select unique triangles if using periodic copies
    unique_triangles = MeshGenerator::select_unique_connected_triangles(
        points_used, triangulation, original_domain_map, square_points.size(),
        unique_triangle_tolerance);
  } else {
    // Reject triangles with overly long diagonals/edges in non-periodic
    // remeshing, scaled by the median nearest-neighbor distance in square_points
    // to be robust against localized core/cut atom pairs.
    int n_pts = static_cast<int>(square_points.size());
    int check_n = std::min(100, n_pts);
    int stride = std::max(1, n_pts / std::max(1, check_n));
    std::vector<double> nn_dists;
    nn_dists.reserve(check_n);
    for (int idx = 0; idx < n_pts && static_cast<int>(nn_dists.size()) < check_n; idx += stride) {
      double d_i = 1e9;
      for (int j = 0; j < n_pts; ++j) {
        if (idx == j) continue;
        double d = (square_points[idx].coord - square_points[j].coord).norm();
        if (d > 1e-4 && d < d_i) {
          d_i = d;
        }
      }
      if (d_i < 1e8) nn_dists.push_back(d_i);
    }
    std::sort(nn_dists.begin(), nn_dists.end());
    double median_nn = nn_dists.empty() ? 1.0 : nn_dists[nn_dists.size() / 2];
    if (g_active_geometry_filter != nullptr && g_active_geometry_filter->a0 > 1e-6) {
      median_nn = g_active_geometry_filter->a0;
    }
    bool is_triangular = (offsets[0] > 1e-6 && std::abs(offsets[1] / offsets[0] - std::sqrt(3.0)) < 0.25);
    // For triangular lattices, core dilated edges span up to ~1.58 a0, while boundary slivers
    // span across second-nearest neighbors (sqrt(3) a0 ≈ 1.732 a0). Using 1.65 * median_nn
    // preserves all core triangles while cleanly filtering out boundary slivers.
    double base_max_edge = is_triangular ? (median_nn * 1.65) : (median_nn * (std::sqrt(2.0) + 0.35));
    double max_edge = base_max_edge;
    if (g_active_geometry_filter != nullptr && g_active_geometry_filter->enabled) {
      double user_gb_max = median_nn * g_active_geometry_filter->microstructure.max_gb_edge_a0;
      if (g_active_geometry_filter->microstructure.mode != MicrostructureMode::SingleCrystal) {
        max_edge = user_gb_max;
      } else {
        max_edge = std::max(base_max_edge, user_gb_max);
      }
    }
    unique_triangles = MeshGenerator::select_unique_connected_triangles(
        points_used, triangulation, original_domain_map, square_points.size(),
        unique_triangle_tolerance, max_edge);

    if (g_active_geometry_filter != nullptr && g_active_geometry_filter->enabled) {
      const std::vector<Point2D>& check_pts =
          (g_active_geometry_ref_points != nullptr &&
           g_active_geometry_ref_points->size() == square_points.size())
              ? *g_active_geometry_ref_points
              : square_points;
      const double gb_max_edge = median_nn * g_active_geometry_filter->microstructure.max_gb_edge_a0;
      const double wall_tol = 0.85 * median_nn;
      const double cx_min = g_active_geometry_filter->carved_x_min;
      const double cx_max = g_active_geometry_filter->carved_x_max;
      const double cy_min = g_active_geometry_filter->carved_y_min;
      const double cy_max = g_active_geometry_filter->carved_y_max;

      auto onSameOuterWall = [&](const Eigen::Vector2d& a, const Eigen::Vector2d& b) -> bool {
        if (std::abs(a.x() - cx_min) < wall_tol && std::abs(b.x() - cx_min) < wall_tol) return true;
        if (std::abs(a.x() - cx_max) < wall_tol && std::abs(b.x() - cx_max) < wall_tol) return true;
        if (std::abs(a.y() - cy_min) < wall_tol && std::abs(b.y() - cy_min) < wall_tol) return true;
        if (std::abs(a.y() - cy_max) < wall_tol && std::abs(b.y() - cy_max) < wall_tol) return true;
        return false;
      };
      auto nearOuterWall = [&](const Eigen::Vector2d& a) -> bool {
        return std::abs(a.x() - cx_min) < wall_tol || std::abs(a.x() - cx_max) < wall_tol ||
               std::abs(a.y() - cy_min) < wall_tol || std::abs(a.y() - cy_max) < wall_tol;
      };

      std::vector<Triangle> filtered_tris;
      filtered_tris.reserve(unique_triangles.size());
      for (const auto& tri : unique_triangles) {
        int i0 = original_domain_map[tri.vertex_indices[0]];
        int i1 = original_domain_map[tri.vertex_indices[1]];
        int i2 = original_domain_map[tri.vertex_indices[2]];
        const Eigen::Vector2d& r0 = check_pts[i0].coord;
        const Eigen::Vector2d& r1 = check_pts[i1].coord;
        const Eigen::Vector2d& r2 = check_pts[i2].coord;
        if (g_active_geometry_filter->shouldRejectTriangle(r0, r1, r2)) {
          continue;
        }

        const Eigen::Vector2d& p0 = square_points[i0].coord;
        const Eigen::Vector2d& p1 = square_points[i1].coord;
        const Eigen::Vector2d& p2 = square_points[i2].coord;
        double l01 = (p1 - p0).norm();
        double l12 = (p2 - p1).norm();
        double l20 = (p0 - p2).norm();
        double l_max = std::max({l01, l12, l20});
        double two_area = std::abs((p1.x() - p0.x()) * (p2.y() - p0.y()) -
                                   (p1.y() - p0.y()) * (p2.x() - p0.x()));
        double h_min = two_area / std::max(1e-12, l_max);

        bool same_grain = (check_pts[i0].grain_id == check_pts[i1].grain_id &&
                           check_pts[i1].grain_id == check_pts[i2].grain_id);
        bool outer_wall_edge = onSameOuterWall(r0, r1) || onSameOuterWall(r1, r2) || onSameOuterWall(r2, r0);
        bool touches_outer_wall = nearOuterWall(r0) || nearOuterWall(r1) || nearOuterWall(r2);

        if (!same_grain && touches_outer_wall) {
          // Where a grain boundary meets an outer wall, strip V-notch bridge triangles
          // that span across the truncated kite (> 1.48 a0) while keeping genuine kite bonds.
          if (l_max > 1.48 * median_nn || h_min < 0.28 * median_nn) {
            continue;
          }
        } else if (outer_wall_edge) {
          bool is_rotated_grain = (std::abs(check_pts[i0].grain_theta_deg) > 1e-3);
          if (!is_rotated_grain) {
            if (l_max > base_max_edge) continue;
          } else {
            // Along a rotated outer wall, strip convex-hull staircase bridges (l_max > base_max_edge)
            if (l_max > base_max_edge || h_min < 0.28 * median_nn) {
              continue;
            }
          }
        } else {
          // Interior bulk and cross-grain-boundary elements: preserve up to gb_max_edge
          if (l_max > gb_max_edge || h_min < 0.08 * median_nn) {
            continue;
          }
        }

        filtered_tris.push_back(tri);
      }
      unique_triangles = std::move(filtered_tris);
    }
  }

  // Create new finite elements
  elements = MeshGenerator::createElementTri2D(
      unique_triangles, square_points, original_domain_map, translation_map);

  // Initialize elements with reference configuration
  for (auto &element : elements) {
    element.set_reference_mesh(square_points);
    element.set_dof_mapping(full_mapping);

    // If manual shape derivative is provided, use it for all elements
    if (manual_shape_derivative) {
      element.set_shape_derivatives(*manual_shape_derivative);
      element.setExternalDeformation(F_ext);

      // Calculate deformation gradient with current positions
      element.calculate_deformation_gradient(x);
      // Also update area calculations
      element.calculateReferenceArea(square_points);
      element.calculateCurrentArea(x);
    } else {
      // Otherwise calculate shape derivatives normally
      double jac = element.calculate_shape_derivatives(x);
      element.setExternalDeformation(F_ext);
    }
  }

  // Sort elements directly by their first node index
  std::sort(elements.begin(), elements.end(),
            [](const ElementTriangle2D &a, const ElementTriangle2D &b) {
              return a.getNodeIndex(0) < b.getNodeIndex(0);
            });

  // Update active elements
  return initialize_active_elements(elements, full_mapping,
                                    square_points.size());
}

/**
 * A simpler interface that only returns the elements and active elements
 */
std::pair<std::vector<ElementTriangle2D>, std::vector<size_t>>
AdaptiveMesher::createMesh(
    const std::vector<Point2D> &square_points, const alglib::real_1d_array &x,
    const Eigen::Matrix2d &F_ext,
    const Eigen::Matrix<double, 3, 2> *manual_shape_derivative) {
  std::vector<ElementTriangle2D> elements;
  std::vector<Point2D> points_used;
  std::vector<Triangle> triangulation;
  std::vector<Triangle> unique_triangles;

  std::vector<size_t> active_elements =
      generateMesh(square_points, x, F_ext, elements, points_used,
                   triangulation, unique_triangles, manual_shape_derivative);

  return {elements, active_elements};
}
/**
 * Save the original positions before remeshing
 */
alglib::real_1d_array
AdaptiveMesher::saveOriginalPositions(const alglib::real_1d_array &x) {
  alglib::real_1d_array original_x_remesh;
  original_x_remesh.setlength(x.length());
  for (int j = 0; j < x.length(); j++) {
    original_x_remesh[j] = x[j];
  }
  return original_x_remesh;
}

/**
 * Set whether to use periodic copies or just the original domain
 */
void AdaptiveMesher::setUsePeriodicCopies(bool use_periodic_copies) {
  use_periodic = use_periodic_copies;
}

/**
 * Get the current setting for using periodic copies
 */
bool AdaptiveMesher::getUsePeriodicCopies() const { return use_periodic; }

void AdaptiveMesher::setTriangulationPerturbation(
    double perturbation_amplitude) {
  triangulation_perturbation_amplitude = perturbation_amplitude;
}

double AdaptiveMesher::getTriangulationPerturbation() const {
  return triangulation_perturbation_amplitude;
}
