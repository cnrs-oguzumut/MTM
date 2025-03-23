#ifndef LBFGS_OPTIMIZER_H
#define LBFGS_OPTIMIZER_H

#include "src/optimization.h" // This path is relative to ALGLIB_DIR

class LBFGSOptimizer {
private:
    int corrections_;
    double epsg_;
    double epsf_;
    double epsx_;
    alglib::ae_int_t maxits_;

public:
    // Constructor with default parameters
    LBFGSOptimizer(
        int corrections = 10,
        double epsg = 0,
        double epsf = 0,
        double epsx = 0,
        alglib::ae_int_t maxits = 0
    );

    // Main optimization method
    void optimize(
        alglib::real_1d_array& x,
        void (*energy_func)(const alglib::real_1d_array& x, double& func, alglib::real_1d_array& grad, void* ptr),
        void* userData
    );
};

#endif // LBFGS_OPTIMIZER_H