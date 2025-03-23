#ifndef TRIANGULATION_VTK_H
#define TRIANGULATION_VTK_H

#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "../include/geometry/Point2D.h"
#include "../include/geometry/DomainDimensions.h"
#include "../include/geometry/LatticeGenerator.h"
#include "../include/geometry/DomainInfo.h"
#include "../include/mesh/MeshGenerator.h"
#include "../include/optimization/LatticeOptimizer.h"

#include "../include/interatomic/inter_atomic.h"

#include "../include/mesh/Triangle.h"
#include "../include/mesh/ElementTriangle2D.h"
#include "../include/mesh/MeshGenerator.h"

void write_triangulation_to_vtk(
    const std::vector<Triangle>& triangles,
    const std::vector<Point2D>& points,
    int iteration_number,
    const std::string& prefix = "triangulation"
);

#endif // TRIANGULATION_VTK_WRITER_H