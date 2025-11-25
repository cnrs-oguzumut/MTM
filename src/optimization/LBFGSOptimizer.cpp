// LBFGSOptimizer.cpp
#include "../include/optimization/LBFGSOptimizer.h"
#include <iostream>

LBFGSOptimizer::LBFGSOptimizer(
    int corrections,
    double epsg,
    double epsf,
    double epsx,
    alglib::ae_int_t maxits
) : corrections_(corrections),
    epsg_(epsg),
    epsf_(epsf),
    epsx_(epsx),
    maxits_(maxits),
    state_created_(false) {}

void LBFGSOptimizer::optimize(
    alglib::real_1d_array& x,
    void (*energy_func)(const alglib::real_1d_array& x, double& func, alglib::real_1d_array& grad, void* ptr),
    void* userData
) {
    // Use local state for non-preconditioned optimization
    alglib::minlbfgsstate state;
    alglib::minlbfgsreport rep;
    
    // Create and configure optimizer
    alglib::minlbfgscreate(corrections_, x, state);
    alglib::minlbfgssetcond(state, epsg_, epsf_, epsx_, maxits_);
    
    // Run optimization
    alglib::minlbfgsoptimize(state, energy_func, nullptr, userData);
    
    // Get results
    alglib::minlbfgsresults(state, x, rep);

    std::cout << "Optimization completed: iterations=" << rep.iterationscount 
              << ", termination type=" << rep.terminationtype << std::endl;
}

void LBFGSOptimizer::createState(alglib::real_1d_array& x) {
    std::cout << "Creating L-BFGS state..." << std::endl;
    
    // Create and configure optimizer state (member variable)
    alglib::minlbfgscreate(corrections_, x, state_);
    alglib::minlbfgssetcond(state_, epsg_, epsf_, epsx_, maxits_);
    
    state_created_ = true;
    
    std::cout << "State created with " << corrections_ << " corrections" << std::endl;
}

void LBFGSOptimizer::setPreconditioner(const alglib::real_1d_array& diag_precond) {
    if (!state_created_) {
        std::cerr << "Error: Must call createState() before setPreconditioner()!" << std::endl;
        throw std::runtime_error("LBFGSOptimizer: state not created");
    }
    
    std::cout << "Setting diagonal preconditioner (length=" << diag_precond.length() << ")..." << std::endl;
    alglib::minlbfgssetprecdiag(state_, diag_precond);
    std::cout << "Preconditioner applied to L-BFGS optimizer" << std::endl;
}

void LBFGSOptimizer::runOptimization(
    alglib::real_1d_array& x,
    void (*energy_func)(const alglib::real_1d_array& x, double& func, alglib::real_1d_array& grad, void* ptr),
    void* userData
) {
    if (!state_created_) {
        std::cerr << "Error: Must call createState() before runOptimization()!" << std::endl;
        throw std::runtime_error("LBFGSOptimizer: state not created");
    }
    
    std::cout << "Running L-BFGS optimization with preconditioner..." << std::endl;
    
    // Run optimization using member state_
    alglib::minlbfgsoptimize(state_, energy_func, nullptr, userData);
    
    // Get results
    alglib::minlbfgsreport rep;
    alglib::minlbfgsresults(state_, x, rep);
    
    std::cout << "Optimization completed: iterations=" << rep.iterationscount 
              << ", termination type=" << rep.terminationtype << std::endl;
    
    state_created_ = false;  // Reset for next use
}

void LBFGSOptimizer::optimizeWithPreconditioner(
    alglib::real_1d_array& x,
    const alglib::real_1d_array& diag_precond,
    void (*energy_func)(const alglib::real_1d_array& x, double& func, alglib::real_1d_array& grad, void* ptr),
    void* userData
) {
    std::cout << "\n=== L-BFGS with Preconditioner ===" << std::endl;
    
    // All-in-one method: create state, set preconditioner, run optimization
    createState(x);
    setPreconditioner(diag_precond);
    runOptimization(x, energy_func, userData);
    
    std::cout << "=== Optimization Complete ===\n" << std::endl;
}