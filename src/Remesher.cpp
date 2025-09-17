#include "../include/mesh/ReMesher.h"
#include "../include/mesh/Triangle.h"
#include "../include/mesh/ElementTriangle2D.h"
#include "../include/geometry/LatticeGenerator.h"
#include "../include/mesh/MeshGenerator.h"
#include "../include/geometry/DomainDimensions.h"

/**
 * Constructor for AdaptiveMesher
 */
AdaptiveMesher::AdaptiveMesher(
    const Point2D& domain_dimensions,
    const std::array<double, 2>& periodic_offsets,
    const std::vector<int>& original_domain_mapping,
    const std::vector<std::tuple<double, double>>& translation_mapping,
    const std::vector<std::pair<int, int>>& dof_mapping,
    double tolerance,
    bool use_periodic_copies // Added parameter for periodic copies option
) : domain_dims(domain_dimensions),
    offsets(periodic_offsets),
    original_domain_map(original_domain_mapping),
    translation_map(translation_mapping),
    full_mapping(dof_mapping),
    unique_triangle_tolerance(tolerance),
    use_periodic(use_periodic_copies) // Initialize the new member
{
}

/**
 * Generate a new mesh based on the current configuration
 */
/**
 * Generate a new mesh based on the current configuration with optional manual shape derivatives
 */
std::vector<size_t> AdaptiveMesher::generateMesh(
    const std::vector<Point2D>& square_points,
    const alglib::real_1d_array& x,
    const Eigen::Matrix2d& F_ext,
    std::vector<ElementTriangle2D>& elements,
    std::vector<Point2D>& points_used,
    std::vector<Triangle>& triangulation,
    std::vector<Triangle>& unique_triangles,
    const Eigen::Matrix<double, 3, 2>* manual_shape_derivative 
) {
    // Create a DomainDimensions object from the Point2D domain_dims
    DomainDimensions domainDimensions(domain_dims.coord.x(), domain_dims.coord.y());
    
    std::vector<Point2D> points_for_triangulation;
    
    if (use_periodic) {
        // Generate new periodic copies with current deformation
        points_for_triangulation = LatticeGenerator::create_periodic_copies(
            square_points, domainDimensions, offsets, F_ext);
    } else {
        // Use original domain without periodic copies
        points_for_triangulation = square_points;
    }
    
    // Create new triangulation
    triangulation = MeshGenerator::createTrianglesFromPoints(points_for_triangulation);
    points_used = points_for_triangulation;
    
    // Select triangles
    if (use_periodic) {
        // Select unique triangles if using periodic copies
        unique_triangles = MeshGenerator::select_unique_connected_triangles(
            points_used, triangulation, original_domain_map,
            square_points.size(), unique_triangle_tolerance
        );
    } else {
        // Use all triangles if not using periodic copies
        unique_triangles = MeshGenerator::select_unique_connected_triangles(
            points_used, triangulation, original_domain_map,
            square_points.size(), unique_triangle_tolerance, false);

    }
    
    // Create new finite elements
    elements = MeshGenerator::createElementTri2D(
        unique_triangles, square_points, original_domain_map, translation_map
    );

    // Initialize elements with reference configuration
    for (auto& element : elements) {
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
    [](const ElementTriangle2D& a, const ElementTriangle2D& b) {
        return a.getNodeIndex(0) < b.getNodeIndex(0);
    });
    
    // Update active elements
    return initialize_active_elements(elements, full_mapping, square_points.size());
}

/**
 * A simpler interface that only returns the elements and active elements
 */
std::pair<std::vector<ElementTriangle2D>, std::vector<size_t>> AdaptiveMesher::createMesh(
    const std::vector<Point2D>& square_points,
    const alglib::real_1d_array& x,
    const Eigen::Matrix2d& F_ext,
    const Eigen::Matrix<double, 3, 2>* manual_shape_derivative 
) {
    std::vector<ElementTriangle2D> elements;
    std::vector<Point2D> points_used;
    std::vector<Triangle> triangulation;
    std::vector<Triangle> unique_triangles;
    
    std::vector<size_t> active_elements = generateMesh(
        square_points, x, F_ext, elements, points_used, triangulation, unique_triangles,
        manual_shape_derivative
    );
    
    return {elements, active_elements};
}
/**
 * Save the original positions before remeshing
 */
alglib::real_1d_array AdaptiveMesher::saveOriginalPositions(const alglib::real_1d_array& x) {
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
bool AdaptiveMesher::getUsePeriodicCopies() const {
    return use_periodic;
}