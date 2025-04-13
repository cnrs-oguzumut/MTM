#include "../include/optimization/CGOptimizer.h"

CGOptimizer::CGOptimizer(
    double epsg,
    double epsf,
    double epsx,
    alglib::ae_int_t maxits
) : epsg_(epsg),
    epsf_(epsf),
    epsx_(epsx),
    maxits_(maxits) {}

void CGOptimizer::optimize(
    alglib::real_1d_array& x,
    void (*energy_func)(const alglib::real_1d_array& x, double& func, alglib::real_1d_array& grad, void* ptr),
    void* userData
) {
    // Set up minimization parameters
    alglib::mincgstate state;
    alglib::mincgreport rep;
    
    // Create and configure optimizer
    alglib::mincgcreate(x, state);
    alglib::mincgsetcond(state, epsg_, epsf_, epsx_, maxits_);
    
    // Run optimization
    alglib::mincgoptimize(state, energy_func, nullptr, userData);
    
    // Get results
    alglib::mincgresults(state, x, rep);

    std::cout << "Optimization completed: iterations=" << rep.iterationscount 
              << ", termination type=" << rep.terminationtype << std::endl;
}