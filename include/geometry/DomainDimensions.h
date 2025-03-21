#pragma once
#ifndef DOMAIN_DIMENSIONS_H
#define DOMAIN_DIMENSIONS_H

struct DomainDimensions {
    double size_x;
    double size_y;
    
    // Default constructor
    DomainDimensions() : size_x(0.0), size_y(0.0) {}
    
    // Constructor with dimensions
    DomainDimensions(double sx, double sy) : size_x(sx), size_y(sy) {}
};

#endif // DOMAIN_DIMENSIONS_H
