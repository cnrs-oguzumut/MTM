#include "../include/reductions/LagrangeReduction.h"

namespace lagrange {
    // Main Lagrange reduction function - returns both reduced tensor, transformation matrix, and flag
    Result reduce(const Eigen::Matrix2d& C) {
        double c11 = C(0, 0);
        double c22 = C(1, 1);
        double c12 = C(0, 1);
        
        // Define the basic transformation matrices
        Eigen::Matrix2d m_identity = Eigen::Matrix2d::Identity();
        static Eigen::Matrix2d m1, m2, m3;
        
        // m1 = [[1, 0], [0, -1]]
        m1 << 1.0, 0.0,
              0.0, -1.0;
        
        // m2 = [[0, 1], [1, 0]]
        m2 << 0.0, 1.0,
              1.0, 0.0;
        
        // m3 = [[1, -1], [0, 1]]
        m3 << 1.0, -1.0,
              0.0, 1.0;
        
        // Initialize transformation matrix as identity
        Eigen::Matrix2d m_matrix = m_identity;
        
        // Flag to track if third condition was satisfied
        bool third_condition_satisfied = false;
        
        // Iteratively apply Lagrange reduction
        while (c12 < 0 || c22 < c11 || 2 * c12 > c11) {
            if (c12 < 0) {
                c12 = -c12;
                // Update transformation matrix
                m_matrix = m_matrix * m1;
            }
            
            if (c22 < c11) {
                double temp = c11;
                c11 = c22;
                c22 = temp;
                // Update transformation matrix
                m_matrix = m_matrix * m2;
            }
            
            if (2 * c12 > c11) {
                double d11 = c11;
                double d12 = c12 - c11;
                double d22 = c22 + c11 - 2 * c12;
                c11 = d11;
                c12 = d12;
                c22 = d22;
                // Update transformation matrix
                m_matrix = m_matrix * m3;
                
                // Set the flag when third condition is satisfied
                third_condition_satisfied = true;
            }
        }
        
        // Create reduced metric tensor
        Eigen::Matrix2d C_reduced;
        C_reduced << c11, c12,
                     c12, c22; // Symmetric matrix
        
        // Prepare the result
        Result result;
        result.C_reduced = C_reduced;
        result.m_matrix = m_matrix;
        result.third_condition_satisfied = third_condition_satisfied;
        
        return result;
    }
    
    // Simple version that only returns the reduced metric tensor
    Eigen::Matrix2d reduceTensor(const Eigen::Matrix2d& C) {
        return reduce(C).C_reduced;
    }

    ResultWithCount reduceAndCount(const Eigen::Matrix2d& C) {
        double c11 = C(0, 0);
        double c22 = C(1, 1);
        double c12 = C(0, 1);
        
        // Define the basic transformation matrices
        Eigen::Matrix2d m_identity = Eigen::Matrix2d::Identity();
        static Eigen::Matrix2d m1, m2, m3;
        
        // m1 = [[1, 0], [0, -1]]
        m1 << 1.0, 0.0,
              0.0, -1.0;
        
        // m2 = [[0, 1], [1, 0]]
        m2 << 0.0, 1.0,
              1.0, 0.0;
        
        // m3 = [[1, -1], [0, 1]]
        m3 << 1.0, -1.0,
              0.0, 1.0;
        
        // Initialize transformation matrix as identity
        Eigen::Matrix2d m_matrix = m_identity;
        
        // Counter for third condition occurrences
        int third_condition_count = 0;
        
        // Iteratively apply Lagrange reduction
        while (c12 < 0 || c22 < c11 || 2 * c12 > c11) {
            if (c12 < 0) {
                c12 = -c12;
                // Update transformation matrix
                m_matrix = m_matrix * m1;
            }
            
            if (c22 < c11) {
                double temp = c11;
                c11 = c22;
                c22 = temp;
                // Update transformation matrix
                m_matrix = m_matrix * m2;
            }
            
            if (2 * c12 > c11) {
                double d11 = c11;
                double d12 = c12 - c11;
                double d22 = c22 + c11 - 2 * c12;
                c11 = d11;
                c12 = d12;
                c22 = d22;
                // Update transformation matrix
                m_matrix = m_matrix * m3;
                
                // Increment counter when third condition is satisfied
                third_condition_count++;
            }
        }
        
        // Create reduced metric tensor
        Eigen::Matrix2d C_reduced;
        C_reduced << c11, c12,
                     c12, c22; // Symmetric matrix
        
        // Prepare the result
        ResultWithCount result;
        result.C_reduced = C_reduced;
        result.m_matrix = m_matrix;
        result.third_condition_count = third_condition_count;
        result.third_condition_satisfied = (third_condition_count > 0);
        return result;
    }

} // namespace lagrange


