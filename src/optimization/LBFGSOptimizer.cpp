#include "../include/optimization/LBFGSOptimizer.h"

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
    maxits_(maxits) {}

void LBFGSOptimizer::optimize(
    alglib::real_1d_array& x,
    void (*energy_func)(const alglib::real_1d_array& x, double& func, alglib::real_1d_array& grad, void* ptr),
    void* userData
) {
    // Set up minimization parameters
    alglib::minlbfgsstate state;
    alglib::minlbfgsreport rep;
    UserData* data = static_cast<UserData*>(userData);
    data->optimizer_state = &state;

    
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