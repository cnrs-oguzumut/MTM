#include "../include/optimization/BLEICOptimizer.h"

using namespace alglib;

BLEICOptimizer::BLEICOptimizer(double epsg, double epsf, double epsx, ae_int_t maxits)
    : epsg_(epsg), epsf_(epsf), epsx_(epsx), maxits_(maxits) {}

void BLEICOptimizer::optimize(
    real_1d_array& x,
    void (*energy_func)(const real_1d_array& x, double& func, real_1d_array& grad, void* ptr),
    void* userData,
    const real_1d_array& bndl,
    const real_1d_array& bndu
) {
    minbleicstate state;
    minbleicreport rep;

    minbleiccreate(x, state);
    minbleicsetcond(state, epsg_, epsf_, epsx_, maxits_);

    if (bndl.length() > 0 && bndu.length() > 0) {
        minbleicsetbc(state, bndl, bndu);
    }

    // Pass function and userData directly (no lambda needed)
    minbleicoptimize(state, energy_func, nullptr, userData);

    minbleicresults(state, x, rep);

    // Optional diagnostics
    std::cout << "BLEIC optimization terminated with code: " << rep.terminationtype << "\n";
    std::cout << "Iterations: " << rep.iterationscount << "\n";
    std::cout << "Function evaluations: " << rep.nfev << "\n";
}
