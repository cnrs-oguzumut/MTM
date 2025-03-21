#pragma once
#ifndef DOMAIN_INFO_H
#define DOMAIN_INFO_H

#include <vector>
#include <Eigen/Dense>
#include "DomainDimensions.h"
#include "Point2D.h"

struct DomainInfo {
    double min_x, max_x;
    double min_y, max_y;
    
    Eigen::Vector2d translation_x;
    Eigen::Vector2d translation_y;
    
    // Optional: Add convenience methods
    double get_width() const { return max_x - min_x; }
    double get_height() const { return max_y - min_y; }
    
    DomainDimensions to_domain_dimensions() const {
        return DomainDimensions(get_width(), get_height());
    }
};

// Function declaration
DomainInfo compute_domain_size(const std::vector<Point2D>& points);

#endif // DOMAIN_INFO_H
