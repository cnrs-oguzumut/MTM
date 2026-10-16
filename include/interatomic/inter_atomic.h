#pragma once

#include <functional>
#include <string>

double lennard_jones_energy_v2(double r);
double lennard_jones_energy_der_v2(double r);
double lennard_jones_energy_sder_v2(double r);

double square_energy(double r);
double square_energy_der(double r);
double square_energy_sder(double r);

double lennard_jones_energy_v3(double r);
double lennard_jones_energy_der_v3(double r);
double lennard_jones_energy_sder_v3(double r);

/**
 * @brief Result of equilibrium lattice parameter calculation for a pair potential V(r)
 */
struct EquilibriumLatticeResult {
    std::string lattice_type = "triangular";
    double lattice_parameter = 1.0; // Equilibrium lattice spacing a_0 (or b) where dE/da = 0
    double cohesive_energy   = 0.0; // Cohesive energy per atom E_coh = E(a_0)
    double energy_density    = 0.0; // Energy per unit area E(a_0) / Omega_0
    double unit_cell_area    = 1.0; // Omega_0 = (sqrt(3)/2)*a_0^2 (triangular) or a_0^2 (square)
    double dE_da             = 0.0; // Residual first derivative dE/da at a_0 (zero-pressure check)
    double d2E_da2           = 0.0; // Second derivative d^2E/da^2 at a_0 (> 0 for stability)
    int iterations           = 0;
    bool converged           = false;

    void printSummary() const;
};

/**
 * @brief Calculate the zero-pressure equilibrium lattice parameter a_0 (where dE/da = 0)
 *        for any interatomic pair potential V(r) on a "triangular" or "square" lattice.
 *
 * Uses a bracketed Newton-Raphson / bisection root finder on dE/da = 0.
 * If dpot or d2pot are omitted (nullptr), accurate central finite differences of pot(r) are used.
 */
EquilibriumLatticeResult compute_equilibrium_lattice_parameter(
    const std::function<double(double)>& pot,
    const std::function<double(double)>& dpot = nullptr,
    const std::function<double(double)>& d2pot = nullptr,
    const std::string& lattice_type = "triangular",
    double r_cutoff = 2.5,
    double a_min = 0.4,
    double a_max = 1.8
);

/**
 * @brief Convenience function returning the equilibrium lattice parameter a_0
 */
double find_optimal_lattice_parameter(
    const std::function<double(double)>& pot,
    const std::string& lattice_type = "triangular",
    double r_cutoff = 2.5,
    const std::function<double(double)>& dpot = nullptr,
    const std::function<double(double)>& d2pot = nullptr
);