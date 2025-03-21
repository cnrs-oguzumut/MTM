#include "../include/geometry/DomainInfo.h"
#include <stdexcept>

DomainInfo compute_domain_size(const std::vector<Point2D>& points) {
    if (points.empty()) {
        throw std::runtime_error("Empty point set");
    }
    
    // Initialize with first point
    double min_x = points[0].coord.x(), max_x = points[0].coord.x();
    double min_y = points[0].coord.y(), max_y = points[0].coord.y();
    
    // Find min and max for each coordinate
    for (const auto& point : points) {
        min_x = std::min(min_x, point.coord.x());
        max_x = std::max(max_x, point.coord.x());
        
        min_y = std::min(min_y, point.coord.y());
        max_y = std::max(max_y, point.coord.y());
    }
    
    // Compute domain sizes
    double domain_x = max_x - min_x;
    double domain_y = max_y - min_y;
    
    // Create and populate domain info
    DomainInfo info;
    info.min_x = min_x;
    info.max_x = max_x;
    info.min_y = min_y;
    info.max_y = max_y;
    
    // Translation vectors are the domain sizes
    info.translation_x = Eigen::Vector2d(domain_x, 0);
    info.translation_y = Eigen::Vector2d(0, domain_y);
    
    return info;
}
