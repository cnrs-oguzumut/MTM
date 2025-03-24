#ifndef LAGRANGE_REDUCTION_H
#define LAGRANGE_REDUCTION_H

#include <Eigen/Dense>

namespace lagrange {
    // Structure to store the reduction result
    struct Result {
        Eigen::Matrix2d C_reduced;       // Reduced metric tensor
        Eigen::Matrix2d m_matrix;        // Transformation matrix
        bool third_condition_satisfied;  // Flag indicating if third condition was met
    };

    // Enhanced structure with count tracking
    struct ResultWithCount {
        Eigen::Matrix2d C_reduced;       // Reduced metric tensor
        Eigen::Matrix2d m_matrix;        // Transformation matrix
        bool third_condition_satisfied;  // Flag indicating if third condition was met at least once
        int third_condition_count;       // Count of how many times the third condition was met
    };

    // Function declarations
    Result reduce(const Eigen::Matrix2d& C);
    ResultWithCount reduceAndCount(const Eigen::Matrix2d& C);
    Eigen::Matrix2d reduceTensor(const Eigen::Matrix2d& C);

} // namespace lagrange

#endif // LAGRANGE_REDUCTION_H
