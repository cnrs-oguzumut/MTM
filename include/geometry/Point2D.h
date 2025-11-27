#pragma once
#ifndef POINT2D_H
#define POINT2D_H
#include <Eigen/Dense>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include "src/ap.h" // Adjust this path relative to your ALGLIB setup

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_2 CGALPoint;

// Point coordinates class
class Point2D {
public:
    Eigen::Vector2d coord;
    
    Point2D() : coord(0.0, 0.0) {}
    Point2D(double x, double y) : coord(x, y) {}
    
    // Convenience accessors
    double x() const { return coord.x(); }
    double y() const { return coord.y(); }
    double& x() { return coord.x(); }
    double& y() { return coord.y(); }
    
    // Convert to CGAL Point_2
    CGALPoint to_cgal_point() const {
        return CGALPoint(coord.x(), coord.y());
    }
    
    // Constructor from CGAL Point_2
    Point2D(const CGALPoint& cgal_point) 
        : coord(cgal_point.x(), cgal_point.y()) {}
    
    // Convert to ALGLIB real_1d_array (for single point queries)
    alglib::real_1d_array to_alglib_point() const {
        alglib::real_1d_array point;
        point.setlength(2);
        point[0] = coord.x();
        point[1] = coord.y();
        return point;
    }
    
    // Constructor from ALGLIB real_1d_array
    Point2D(const alglib::real_1d_array& alglib_point) 
        : coord(alglib_point[0], alglib_point[1]) {}
    
    // Static method: Convert vector of Point2D to ALGLIB real_2d_array
    static alglib::real_2d_array to_alglib_array(const std::vector<Point2D>& points) {
        alglib::real_2d_array xy;
        const int n = points.size();
        xy.setlength(n, 2);
        
        for (int i = 0; i < n; ++i) {
            xy[i][0] = points[i].coord.x();
            xy[i][1] = points[i].coord.y();
        }
        
        return xy;
    }
    
    // Static method: Convert ALGLIB real_2d_array to vector of Point2D
    static std::vector<Point2D> from_alglib_array(const alglib::real_2d_array& xy) {
        std::vector<Point2D> points;
        const int n = xy.rows();
        points.reserve(n);
        
        for (int i = 0; i < n; ++i) {
            points.emplace_back(xy[i][0], xy[i][1]);
        }
        
        return points;
    }
};

#endif // POINT2D_H