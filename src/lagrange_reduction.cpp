#include "../include/reductions/LagrangeReduction.h"

namespace lagrange {

    ResultWithCount reduceAndCount(const Eigen::Matrix2d& C) {
        double c11 = C(0, 0);
        double c22 = C(1, 1);
        double c12 = C(0, 1);
        
        // Fast path: already in fundamental reduced domain
        if (c12 >= 0.0 && c22 >= c11 && 2.0 * c12 <= c11) {
            return ResultWithCount{C, Eigen::Matrix2d::Identity(), false, 0};
        }

        // Initialize transformation matrix elements directly
        double m00 = 1.0, m01 = 0.0;
        double m10 = 0.0, m11 = 1.0;
        int third_condition_count = 0;
        
        // Iteratively apply Lagrange reduction
        while (c12 < 0.0 || c22 < c11 || 2.0 * c12 > c11) {
            if (c12 < 0.0) {
                c12 = -c12;
                m01 = -m01;
                m11 = -m11;
            }
            
            if (c22 < c11) {
                std::swap(c11, c22);
                std::swap(m00, m01);
                std::swap(m10, m11);
            }
            
            if (2.0 * c12 > c11) {
                double d12 = c12 - c11;
                double d22 = c22 + c11 - 2.0 * c12;
                c12 = d12;
                c22 = d22;
                m01 -= m00;
                m11 -= m10;
                third_condition_count++;
            }
        }
        
        // Create reduced metric tensor
        Eigen::Matrix2d C_reduced;
        C_reduced << c11, c12,
                     c12, c22;
        
        Eigen::Matrix2d m_matrix;
        m_matrix << m00, m01,
                    m10, m11;

        ResultWithCount result;
        result.C_reduced = C_reduced;
        result.m_matrix = m_matrix;
        result.third_condition_count = third_condition_count;
        result.third_condition_satisfied = (third_condition_count > 0);
        return result;
    }

} // namespace lagrange


