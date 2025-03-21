#pragma once
#ifndef POINT2D_H
#define POINT2D_H

#include <Eigen/Dense>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_2 CGALPoint;

// Point coordinates class
class Point2D {
public:
    Point2D() : coord(0.0, 0.0) {}
    Point2D(double x, double y) : coord(x, y) {}
    Eigen::Vector2d coord;

    // Convert to CGAL Point_2
    CGALPoint to_cgal_point() const {
        return CGALPoint(coord.x(), coord.y());
    }
    
    // Optional: Constructor from CGAL Point_2
    Point2D(const CGALPoint& cgal_point) 
        : coord(cgal_point.x(), cgal_point.y()) {}
};

#endif // POINT2D_H
