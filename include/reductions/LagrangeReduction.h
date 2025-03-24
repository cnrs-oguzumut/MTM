#ifndef LAGRANGE_REDUCTION_H
#define LAGRANGE_REDUCTION_H

#include <Eigen/Dense>

namespace lagrange {
    // Define result structure for Lagrange reduction
    struct Result {
        Eigen::Matrix2d C_reduced;         // Reduced metric tensor
        Eigen::Matrix2d m_matrix;          // Transformation matrix
        bool third_condition_satisfied;    // Flag indicating if third condition was met
    };

    // Main Lagrange reduction function - returns both reduced tensor and transformation matrix
    Result reduce(const Eigen::Matrix2d& C);

    // Simple version that only returns the reduced metric tensor
    Eigen::Matrix2d reduceTensor(const Eigen::Matrix2d& C);
} // namespace lagrange

#endif // LAGRANGE_REDUCTION_H