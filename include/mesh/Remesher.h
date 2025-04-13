#pragma once

#ifndef REMESHER_H
#define REMESHER_H

#include <vector>
#include <map>
#include <array>
#include <tuple>
#include <utility>
#include <Eigen/Dense>
#include "src/optimization.h"
#include "../include/geometry/Point2D.h"

// Forward declarations
class Triangle;
class ElementTriangle2D;

/**
 * AdaptiveMesher - A class for generating and updating finite element meshes
 * based on lattice structures with periodic boundary conditions.
 */
class AdaptiveMesher {
public:
    /**
     * Constructor for AdaptiveMesher
     * 
     * @param domain_dimensions The dimensions of the simulation domain
     * @param periodic_offsets Offset vectors for periodic copies
     * @param original_domain_mapping Mapping from original to periodic domain
     * @param translation_mapping Translation mapping for periodic copies
     * @param dof_mapping Degree of freedom mapping for elements
     * @param tolerance Tolerance for unique triangle selection
     * @param use_periodic_copies Whether to use periodic copies (default: true)
     */
    AdaptiveMesher(
        const Point2D& domain_dimensions,
        const std::array<double, 2>& periodic_offsets,
        const std::vector<int>& original_domain_mapping,
        const std::vector<std::tuple<double, double>>& translation_mapping,
        const std::vector<std::pair<int, int>>& dof_mapping,
        double tolerance = 1e-6,
        bool use_periodic_copies = true
    );

    /**
     * Generate a new mesh based on the current configuration
     * 
     * @param square_points Current point positions
     * @param x Current position vector for optimization
     * @param F_ext External deformation gradient
     * @param elements Output parameter for the generated elements
     * @param points_used Output parameter for points used in triangulation
     * @param triangulation Output parameter for the triangulation
     * @param unique_triangles Output parameter for unique triangles
     * @return Vector of active element indices
     */
    std::vector<size_t> generateMesh(
        const std::vector<Point2D>& square_points,
        const alglib::real_1d_array& x,
        const Eigen::Matrix2d& F_ext,
        std::vector<ElementTriangle2D>& elements,
        std::vector<Point2D>& points_used,
        std::vector<Triangle>& triangulation,
        std::vector<Triangle>& unique_triangles,
        const Eigen::Matrix<double, 3, 2>* manual_shape_derivative = nullptr
    );
    
    /**
     * A simpler interface that only returns the elements and active elements
     * 
     * @param square_points Current point positions
     * @param x Current position vector for optimization
     * @param F_ext External deformation gradient
     * @return Pair of element vector and active element indices
     */
    std::pair<std::vector<ElementTriangle2D>, std::vector<size_t>> createMesh(
        const std::vector<Point2D>& square_points,
        const alglib::real_1d_array& x,
        const Eigen::Matrix2d& F_ext,
        const Eigen::Matrix<double, 3, 2>* manual_shape_derivative = nullptr
    );
    /**
     * Save the original positions before remeshing
     * 
     * @param x Current position array
     * @return Copy of the original position array
     */
    alglib::real_1d_array saveOriginalPositions(const alglib::real_1d_array& x);

    /**
     * Set whether to use periodic copies
     * @param use_periodic_copies True to use periodic copies, false to use original domain only
     */
    void setUsePeriodicCopies(bool use_periodic_copies);

    /**
     * Get the current setting for using periodic copies
     * @return True if using periodic copies, false if using original domain only
     */
    bool getUsePeriodicCopies() const;

    /**
     * Get the current domain dimensions
     * @return Domain dimensions
     */
    Point2D getDomainDimensions() const { return domain_dims; }

    /**
     * Get the current periodic offsets
     * @return Array of periodic offsets
     */
    const std::array<double, 2>& getOffsets() const { return offsets; }

    /**
     * Get the original domain mapping
     * @return Vector of original domain mapping
     */
    const std::vector<int>& getOriginalDomainMap() const { return original_domain_map; }

    /**
     * Get the translation mapping
     * @return Translation mapping for periodic copies
     */
    const std::vector<std::tuple<double, double>>& getTranslationMap() const { return translation_map; }

    /**
     * Get the DOF mapping
     * @return Degree of freedom mapping
     */
    const std::vector<std::pair<int, int>>& getFullMapping() const { return full_mapping; }

    /**
     * Get the unique triangle tolerance
     * @return Tolerance value
     */
    double getUniqueTriangleTolerance() const { return unique_triangle_tolerance; }

    /**
     * Set the domain dimensions
     * @param dimensions New domain dimensions
     */
    void setDomainDimensions(const Point2D& dimensions) { domain_dims = dimensions; }

    /**
     * Set the periodic offsets
     * @param periodic_offsets New periodic offsets
     */
    void setOffsets(const std::array<double, 2>& periodic_offsets) { offsets = periodic_offsets; }

    /**
     * Set the original domain mapping
     * @param mapping New original domain mapping
     */
    void setOriginalDomainMap(const std::vector<int>& mapping) { original_domain_map = mapping; }

    /**
     * Set the translation mapping
     * @param mapping New translation mapping
     */
    void setTranslationMap(const std::vector<std::tuple<double, double>>& mapping) { translation_map = mapping; }

    /**
     * Set the DOF mapping
     * @param mapping New DOF mapping
     */
    void setFullMapping(const std::vector<std::pair<int, int>>& mapping) { full_mapping = mapping; }

    /**
     * Set the unique triangle tolerance
     * @param tolerance New tolerance value
     */
    void setUniqueTriangleTolerance(double tolerance) { unique_triangle_tolerance = tolerance; }

private:
    Point2D domain_dims;
    std::array<double, 2> offsets;
    std::vector<int> original_domain_map;
    std::vector<std::tuple<double, double>> translation_map;
    std::vector<std::pair<int, int>> full_mapping;
    double unique_triangle_tolerance;
    bool use_periodic;  // Flag to determine whether to use periodic copies
};

// Declaration for the initialize_active_elements function (defined elsewhere)
std::vector<size_t> initialize_active_elements(
    const std::vector<ElementTriangle2D>& elements, 
    const std::vector<std::pair<int, int>>& mapping,
    size_t num_points
);

#endif // REMESHER_H