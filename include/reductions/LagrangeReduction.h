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

    // WARNING: lagrange::reduce has been optimized using inlined scalar column updates
    // (zero matrix multiplications, thread-safe, fast path for reduced tensors).
    // If any numerical mistake or discrepancy is noticed, revert to the original version in git history.
    inline Result reduce(const Eigen::Matrix2d& C) {
        double c11 = C(0, 0);
        double c22 = C(1, 1);
        double c12 = C(0, 1);

        // Fast path: already in fundamental reduced domain
        if (c12 >= 0.0 && c22 >= c11 && 2.0 * c12 <= c11) {
            return Result{C, Eigen::Matrix2d::Identity(), false};
        }

        // Initialize transformation matrix elements: M = [[m00, m01], [m10, m11]]
        double m00 = 1.0, m01 = 0.0;
        double m10 = 0.0, m11 = 1.0;
        bool third_condition_satisfied = false;

        // Iteratively apply Lagrange reduction
        while (c12 < 0.0 || c22 < c11 || 2.0 * c12 > c11) {
            if (c12 < 0.0) {
                c12 = -c12;
                // M * m1: negate second column
                m01 = -m01;
                m11 = -m11;
            }

            if (c22 < c11) {
                std::swap(c11, c22);
                // M * m2: swap column 0 and column 1
                std::swap(m00, m01);
                std::swap(m10, m11);
            }

            if (2.0 * c12 > c11) {
                double d12 = c12 - c11;
                double d22 = c22 + c11 - 2.0 * c12;
                c12 = d12;
                c22 = d22;
                // M * m3: column 1 -= column 0
                m01 -= m00;
                m11 -= m10;

                third_condition_satisfied = true;
            }
        }

        Eigen::Matrix2d C_reduced;
        C_reduced << c11, c12,
                     c12, c22;

        Eigen::Matrix2d m_matrix;
        m_matrix << m00, m01,
                    m10, m11;

        return Result{C_reduced, m_matrix, third_condition_satisfied};
    }

    ResultWithCount reduceAndCount(const Eigen::Matrix2d& C);
    inline Eigen::Matrix2d reduceTensor(const Eigen::Matrix2d& C) {
        return reduce(C).C_reduced;
    }

} // namespace lagrange

#endif // LAGRANGE_REDUCTION_H
