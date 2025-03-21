#pragma once
#ifndef TRIANGLE_H
#define TRIANGLE_H

#include <array>

// Structure to represent a triangle
struct Triangle {
    std::array<int, 3> vertex_indices;
    
    Triangle(int v0, int v1, int v2) {
        vertex_indices[0] = v0;
        vertex_indices[1] = v1;
        vertex_indices[2] = v2;
    }
};

#endif // TRIANGLE_H
