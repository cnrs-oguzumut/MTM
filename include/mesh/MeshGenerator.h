#pragma once
#ifndef MESH_GENERATOR_H
#define MESH_GENERATOR_H

#include <vector>
#include <array>
#include <tuple>
#include <set>
#include <algorithm>
#include "../geometry/Point2D.h"
#include "../geometry/DomainDimensions.h"
#include "Triangle.h"
#include "CGALTypes.h"
#include "ElementTriangle2D.h"

class MeshGenerator {
public:
    // Create triangles from points using CGAL Delaunay triangulation
    static std::vector<Triangle> createTrianglesFromPoints(const std::vector<Point2D>& points);
    
    // Create domain maps for translating points
    static std::pair<std::vector<int>, std::vector<std::tuple<double, double>>>
    create_domain_maps(int original_domain_size, const DomainDimensions& domain_dims, 
                      const std::array<double, 2>& offsets);
    
    // Select unique connected triangles from the triangulation
    static std::vector<Triangle> select_unique_connected_triangles(
        const std::vector<Point2D>& all_points,
        const std::vector<Triangle>& duplicated_triangles,
        const std::vector<int>& original_domain_map,
        int original_domain_size,
        double min_jacobian_threshold = 1e-6);
    
    // Create ElementTriangle2D objects from triangles
    static std::vector<ElementTriangle2D> createElementTri2D(
        const std::vector<Triangle>& unique_triangles,
        const std::vector<Point2D>& points,
        const std::vector<int>& original_domain_map,
        const std::vector<std::tuple<double, double>>& translation_map);
};

#endif // MESH_GENERATOR_H
