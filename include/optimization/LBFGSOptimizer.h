// LBFGSOptimizer.h
#ifndef LBFGS_OPTIMIZER_H
#define LBFGS_OPTIMIZER_H

#include "src/optimization.h" // This path is relative to ALGLIB_DIR
#include "../include/optimization/LatticeOptimizer.h"

class LBFGSOptimizer {
private:
    int corrections_;
    double epsg_;
    double epsf_;
    double epsx_;
    alglib::ae_int_t maxits_;
    alglib::minlbfgsstate state_;
    bool state_created_;

public:
    LBFGSOptimizer(
        int corrections = 12,
        double epsg = 0.0,
        double epsf = 1e-13,
        double epsx = 0.0,
        alglib::ae_int_t maxits = 0
    );

    // Original optimize method (without preconditioner)
    void optimize(
        alglib::real_1d_array& x,
        void (*energy_func)(const alglib::real_1d_array&, double&, alglib::real_1d_array&, void*),
        void* userData
    );

    // New methods for preconditioner support
    void createState(alglib::real_1d_array& x);
    
    void setPreconditioner(const alglib::real_1d_array& diag_precond);
    
    void runOptimization(
        alglib::real_1d_array& x,
        void (*energy_func)(const alglib::real_1d_array&, double&, alglib::real_1d_array&, void*),
        void* userData
    );
    
    // Combined method: create, set preconditioner, and run
    void optimizeWithPreconditioner(
        alglib::real_1d_array& x,
        const alglib::real_1d_array& diag_precond,
        void (*energy_func)(const alglib::real_1d_array&, double&, alglib::real_1d_array&, void*),
        void* userData
    );
};

#endif // LBFGS_OPTIMIZER_H