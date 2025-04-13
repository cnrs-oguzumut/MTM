#ifndef BLEIC_OPTIMIZER_H
#define BLEIC_OPTIMIZER_H

#include "src/optimization.h" // Adjust this path relative to your ALGLIB setup
#include <iostream>

class BLEICOptimizer {
public:
    /**
     * @brief Construct a new BLEIC Optimizer object
     * 
     * @param epsg Gradient norm stopping condition (default: 1e-6)
     * @param epsf Function value change stopping condition (default: 0)
     * @param epsx Parameter change stopping condition (default: 0)
     * @param maxits Maximum number of iterations (default: 0 = unlimited)
     */
    BLEICOptimizer(
        double epsg = 1e-6,
        double epsf = 0.0,
        double epsx = 0.0,
        alglib::ae_int_t maxits = 0
    );

    /**
     * @brief Optimize the given function using the BLEIC method.
     * 
     * @param x Initial guess (input) and result (output)
     * @param energy_func User-defined energy function (must compute gradient)
     * @param userData Optional pointer to user data passed to energy_func
     * @param bndl Lower bounds on variables
     * @param bndu Upper bounds on variables
     */
    void optimize(
        alglib::real_1d_array& x,
        void (*energy_func)(const alglib::real_1d_array& x, double& func, 
                            alglib::real_1d_array& grad, void* ptr),
        void* userData = nullptr,
        const alglib::real_1d_array& bndl = alglib::real_1d_array(),
        const alglib::real_1d_array& bndu = alglib::real_1d_array()
    );

private:
    double epsg_;
    double epsf_;
    double epsx_;
    alglib::ae_int_t maxits_;
};

#endif // BLEIC_OPTIMIZER_H
