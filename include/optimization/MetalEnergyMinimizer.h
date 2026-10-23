#pragma once
#ifndef METAL_ENERGY_MINIMIZER_H
#define METAL_ENERGY_MINIMIZER_H

#include <vector>
#include "src/ap.h"

struct UserData;

namespace metal_gpu {

    /// Checks if a compatible Apple Silicon Metal GPU device is available.
    bool is_available();

    /// Invalidates cached GPU element buffers (e.g., after adaptive remeshing).
    void invalidate_mesh_cache();

    /// Computes the total energy and nodal gradient on the GPU using Metal.
    /// Falls back to CPU if Metal is not supported or if the calculator is not
    /// Strain_Energy_LatticeCalculator.
    void compute_energy_and_gradient(
        const alglib::real_1d_array &x,
        double &func,
        alglib::real_1d_array &grad,
        UserData *userData);

} // namespace metal_gpu

#endif // METAL_ENERGY_MINIMIZER_H
