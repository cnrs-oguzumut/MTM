#ifndef CG_OPTIMIZER_H
#define CG_OPTIMIZER_H

#include "src/optimization.h" // This path is relative to ALGLIB_DIR
#include <iostream>

class CGOptimizer {
public:
    /**
     * @brief Construct a new Conjugate Gradient Optimizer object
     * 
     * @param epsg Gradient stopping condition (default: 1e-6)
     * @param epsf Function stopping condition (default: 0)
     * @param epsx Parameter stopping condition (default: 0)
     * @param maxits Maximum number of iterations (default: 0 = unlimited)
     */
    CGOptimizer(
        double epsg = 1e-6,
        double epsf = 0.0,
        double epsx = 0.0,
        alglib::ae_int_t maxits = 0
    );

    /**
     * @brief Optimize the given function using Conjugate Gradient method
     * 
     * @param x Initial starting point (will contain result after optimization)
     * @param energy_func Function to optimize (should compute both energy and gradient)
     * @param userData Optional user data to pass to energy function
     */
    void optimize(
        alglib::real_1d_array& x,
        void (*energy_func)(const alglib::real_1d_array& x, double& func, 
                           alglib::real_1d_array& grad, void* ptr),
        void* userData = nullptr
    );

private:
    double epsg_;       // Gradient stopping condition
    double epsf_;       // Function stopping condition
    double epsx_;       // Parameter stopping condition
    alglib::ae_int_t maxits_;  // Maximum number of iterations
};

#endif // CG_OPTIMIZER_H