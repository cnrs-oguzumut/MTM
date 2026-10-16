#pragma once
#ifndef INTER_ATOMIC_H
#define INTER_ATOMIC_H
//#include "interatomic/inter_atomic.h"
#include "../include/interatomic/inter_atomic.h"

#include "src/optimization.h" // This path is relative to ALGLIB_DIR

#include <cmath>


double lennard_jones_energy_v2(double r){
	 
	
	const double Pi = alglib :: pi ();; // Define Pi
	
	double epsilon = 0.5;
	double sigma = 2 * sin(Pi / 10.0);
	double rin ;// = 2 * sigma;
	double rcut;// = 2.5 * sigma;
	
	
	 rin = 2. ;
	 rcut =2.5 ;
	
	
	
	
	double C1 = 24. * epsilon * pow(sigma, 6.) * (pow(rin, 6.) - 2. * pow(sigma, 6)) / pow(rin, 13.);
	double C2 = 12. * epsilon * pow(sigma, 6.) * (26. * pow(sigma, 6) - 7 * pow(rin, 6)) / pow(rin, 14.);
	double C3 = -(3. * C1 + 4 * C2 * (rcut - rin)) / (3.*pow((rcut - rin), 2));
	double C4 = (C1 + C2 * (rcut - rin)) / (2.*pow((rcut - rin), 3));
		
	double C0 = -(rcut - rin) * (3. * C1 + C2 * (rcut - rin)) / 6.;

	double A = C0 - 4. * epsilon * (pow(sigma / rin, 12.) - pow(sigma / rin, 6.));
	
	double energy=0;
	
	if(r<rin)
		energy =4 * epsilon * (pow(sigma / r, 12) - pow(sigma / r, 6)) + A;
	
	if(r>=rin && r<=rcut )
		energy  =  C0 * pow((r - rin), 0.) + C1 * pow((r - rin), 1.) + C2 * pow((r - rin), 2.) + 
           C3 * pow((r - rin), 3.) + C4 * pow((r - rin), 4.);
    
    if(r>rcut)
    	energy=0;


	return energy;

}

double lennard_jones_energy_der_v2(double r) {
    const double Pi = alglib::pi(); // Define Pi

    double epsilon = 0.5;
    double sigma = 2 * sin(Pi / 10.0);

    double rin = 2.0;
    double rcut = 2.5;

    double C1 = 24.0 * epsilon * pow(sigma, 6.0) * (pow(rin, 6.0) - 2.0 * pow(sigma, 6.0)) / pow(rin, 13.0);
    double C2 = 12.0 * epsilon * pow(sigma, 6.0) * (26.0 * pow(sigma, 6.0) - 7.0 * pow(rin, 6.0)) / pow(rin, 14.0);
    double C3 = -(3.0 * C1 + 4.0 * C2 * (rcut - rin)) / (3.0 * pow((rcut - rin), 2.0));
    double C4 = (C1 + C2 * (rcut - rin)) / (2.0 * pow((rcut - rin), 3.0));

//     double C0 = -(rcut - rin) * (3.0 * C1 + C2 * (rcut - rin)) / 6.0;
//     double A = C0 - 4.0 * epsilon * (pow(sigma / rin, 12.0) - pow(sigma / rin, 6.0));

    double derivative = 0.0;

    if (r < rin) {
        // Derivative of 4 * epsilon * ( (sigma / r)^12 - (sigma / r)^6 ) + A
        double term1 = -48.0 * epsilon * pow(sigma, 12.0) / pow(r, 13.0); // Derivative of (sigma / r)^12
        double term2 = 24.0 * epsilon * pow(sigma, 6.0) / pow(r, 7.0);    // Derivative of (sigma / r)^6
        derivative = term1 + term2;
    } else if (r >= rin && r <= rcut) {
        // Derivative of C0 + C1 * (r - rin) + C2 * (r - rin)^2 + C3 * (r - rin)^3 + C4 * (r - rin)^4
        derivative = C1 + 2.0 * C2 * (r - rin) + 3.0 * C3 * pow((r - rin), 2.0) + 4.0 * C4 * pow((r - rin), 3.0);
    } else if (r > rcut) {
        derivative = 0.0; // The energy is zero for r > rcut, so its derivative is also zero.
    }

    return derivative;
}
// More efficient version with better numerical stability


double lennard_jones_energy_sder_v2(double r){
    const double Pi = alglib::pi(); // Define Pi
    double epsilon = 0.5;
    double sigma = 2 * sin(Pi / 10.0);
    double rin = 2.;
    double rcut = 2.5;
    
    // Calculate constants (same as original function)
    double C1 = 24. * epsilon * pow(sigma, 6.) * (pow(rin, 6.) - 2. * pow(sigma, 6)) / pow(rin, 13.);
    double C2 = 12. * epsilon * pow(sigma, 6.) * (26. * pow(sigma, 6) - 7 * pow(rin, 6)) / pow(rin, 14.);
    double C3 = -(3. * C1 + 4 * C2 * (rcut - rin)) / (3.*pow((rcut - rin), 2));
    double C4 = (C1 + C2 * (rcut - rin)) / (2.*pow((rcut - rin), 3));
    
    double second_derivative = 0;
    
    if(r < rin){
        // Second derivative for r < rin
        second_derivative = 624. * epsilon * pow(sigma, 12.) / pow(r, 14.) 
                          - 168. * epsilon * pow(sigma, 6.) / pow(r, 8.);
    }
    else if(r >= rin && r <= rcut){
        // Second derivative for rin <= r <= rcut (polynomial)
        double dr = r - rin;
        second_derivative = 2. * C2 + 6. * C3 * dr + 12. * C4 * pow(dr, 2.);
    }
    else{ // r > rcut
        second_derivative = 0;
    }
    
    return second_derivative;
}

// double lennard_jones_energy_sder_v2(double r) {
//     // Constants
//     const double Pi = alglib::pi(); // Define Pi
//     const double epsilon = 0.5;
//     const double sigma = 2.0 * sin(Pi / 10.0);
//     const double rin = 2.0;
//     const double rcut = 2.5;
    
//     // Pre-compute powers of sigma
//     const double sigma2 = sigma * sigma;
//     const double sigma6 = sigma2 * sigma2 * sigma2;
//     const double sigma12 = sigma6 * sigma6;
    
//     // Pre-compute rin powers
//     const double rin6 = pow(rin, 6.0);
//     const double rin13 = pow(rin, 13.0);
//     const double rin14 = rin13 * rin;
    
//     // Compute coefficients
//     const double C1 = 24.0 * epsilon * sigma6 * (rin6 - 2.0 * sigma6) / rin13;
//     const double C2 = 12.0 * epsilon * sigma6 * (26.0 * sigma6 - 7.0 * rin6) / rin14;
//     const double delta_r = rcut - rin;
//     const double delta_r2 = delta_r * delta_r;
//     const double delta_r3 = delta_r2 * delta_r;
//     const double C3 = -(3.0 * C1 + 4.0 * C2 * delta_r) / (3.0 * delta_r2);
//     const double C4 = (C1 + C2 * delta_r) / (2.0 * delta_r3);
    
//     double second_derivative = 0.0;
    
//     if (r < rin) {
//         // More stable computation for region 1
//         double r2 = r * r;
//         double r4 = r2 * r2;
//         double r8 = r4 * r4;
//         double r14 = r8 * r4 * r2;
        
//         second_derivative = 624.0 * epsilon * sigma12 / r14 - 168.0 * epsilon * sigma6 / r8;
//     }
//     else if (r <= rcut) {
//         // Region 2
//         double dr = r - rin;
//         second_derivative = 2.0 * C2 + 6.0 * C3 * dr + 12.0 * C4 * dr * dr;
//     }
//     // else: r > rcut, second_derivative remains 0.0
    
//     return second_derivative;
// }

double square_energy(double r) {
    // Constants
    const double e = 2.71828182846;
    const double a = 1.0;
    const double c1 = 2.0 * a;
    const double c2 = 2.0 * a;
    const double b1 = 8.0;
    const double b2 = 8.0;  // Note: b2 isn't actually used in the formula
    const double r1 = 1.0;
    const double r2 = 1.425;
    
    // Simplified energy calculation
    double repulsive_term = a / pow(r, 12);
    double attractive_term1 = c1 * exp(-b1 * pow(r - r1, 2));
    double attractive_term2 = c2 * exp(-b1 * pow(r - r2, 2));
    
    return repulsive_term - attractive_term1 - attractive_term2;
}

double square_energy_der(double r) {
    // Constants
    const double e = 2.71828182846;
    const double a = 1.0;
    const double c1 = 2.0 * a;
    const double c2 = 2.0 * a;
    const double b1 = 8.0;
    const double r1 = 1.0;
    const double r2 = 1.425;
    
    // Derivative terms
    double repulsive_der = -12.0 * a / pow(r, 13);
    double attractive_der1 = 2.0 * b1 * c1 * (r - r1) * exp(-b1 * pow(r - r1, 2));
    double attractive_der2 = 2.0 * b1 * c2 * (r - r2) * exp(-b1 * pow(r - r2, 2));
    
    return repulsive_der + attractive_der1 + attractive_der2;
}

#include <cmath>

// Potential energy function U(r)
double lennard_jones_energy_v3(double r) {
    // Constants
    const double sigma = 1.0;
    const double rc = 1.86602540378444;  // sigma*((sqrt(3)+2)/2)
    const double A = -4.47124703217530;
    const double B = -90.2632082644540;
    const double C2p0 = -5530.88442798764;
    const double C2p1 = 20688.7278150076;
    const double C2p2 = -34486.6459222351;
    const double C2p3 = 32704.6608372787;
    const double C2p4 = -18879.7974113434;
    const double C2p5 = 6595.54667916880;
    const double C2p6 = -1286.11734180526;
    const double C2p7 = 107.717810684010;
    
    double energy = 0.0;
    
    if (r < sigma) {
        // U1 = -2/r^4 + r^(-8)
        double r4 = r * r * r * r;
        double r8 = r4 * r4;
        energy = -2.0 / r4 + 1.0 / r8;
    }
    else if (r >= sigma && r <= rc) {
        // U2 = A/r^8 - B/r^4 + polynomial
        double r2 = r * r;
        double r3 = r2 * r;
        double r4 = r3 * r;
        double r5 = r4 * r;
        double r6 = r5 * r;
        double r7 = r6 * r;
        double r8 = r7 * r;
        
        energy = A / r8 - B / r4 + 
                C2p0 + C2p1 * r + C2p2 * r2 + C2p3 * r3 + 
                C2p4 * r4 + C2p5 * r5 + C2p6 * r6 + C2p7 * r7;
    }
    else if (r > rc) {
        // U3 = 0
        energy = 0.0;
    }
    
    return energy;
}

// Derivative of potential energy function U'(r)
double lennard_jones_energy_der_v3(double r) {
    // Constants
    const double sigma = 1.0;
    const double rc = 1.86602540378444;  // sigma*((sqrt(3)+2)/2)
    const double A = -4.47124703217530;
    const double B = -90.2632082644540;
    const double C2p0 = -5530.88442798764;
    const double C2p1 = 20688.7278150076;
    const double C2p2 = -34486.6459222351;
    const double C2p3 = 32704.6608372787;
    const double C2p4 = -18879.7974113434;
    const double C2p5 = 6595.54667916880;
    const double C2p6 = -1286.11734180526;
    const double C2p7 = 107.717810684010;
    
    double derivative = 0.0;
    
    if (r < sigma) {
        // U1' = 8/r^5 - 8/r^9
        double r5 = r * r * r * r * r;
        double r9 = r5 * r * r * r * r;
        derivative = 8.0 / r5 - 8.0 / r9;
    }
    else if (r >= sigma && r <= rc) {
        // U2' = -8A/r^9 + 4B/r^5 + derivative of polynomial
        double r2 = r * r;
        double r3 = r2 * r;
        double r4 = r3 * r;
        double r5 = r4 * r;
        double r6 = r5 * r;
        double r9 = r5 * r4;
        
        derivative = -8.0 * A / r9 + 4.0 * B / r5 +
                    C2p1 + 2.0 * C2p2 * r + 3.0 * C2p3 * r2 + 
                    4.0 * C2p4 * r3 + 5.0 * C2p5 * r4 + 
                    6.0 * C2p6 * r5 + 7.0 * C2p7 * r6;
    }
    else if (r > rc) {
        // U3' = 0
        derivative = 0.0;
    }
    
    return derivative;
}



double lennard_jones_energy_sder_v3(double r) {
    // Constants
    const double sigma = 1.0;
    const double rc = 1.86602540378444; // sigma*((sqrt(3)+2)/2)
    const double A = -4.47124703217530;
    const double B = -90.2632082644540;
    const double C2p0 = -5530.88442798764;  // Not used in second derivative
    const double C2p1 = 20688.7278150076;   // Not used in second derivative
    const double C2p2 = -34486.6459222351;
    const double C2p3 = 32704.6608372787;
    const double C2p4 = -18879.7974113434;
    const double C2p5 = 6595.54667916880;
    const double C2p6 = -1286.11734180526;
    const double C2p7 = 107.717810684010;
    
    double second_derivative = 0.0;
    
    if (r < sigma) {
        // U1 = -2/r^4 + 1/r^8
        // dU1/dr = 8/r^5 - 8/r^9
        // d²U1/dr² = -40/r^6 + 72/r^10
        double r6 = r * r * r * r * r * r;
        double r10 = r6 * r * r * r * r;
        second_derivative = -40.0 / r6 + 72.0 / r10;
    }
    else if (r >= sigma && r <= rc) {
        // U2 = A/r^8 - B/r^4 + polynomial
        // dU2/dr = -8A/r^9 + 4B/r^5 + C2p1 + 2*C2p2*r + 3*C2p3*r^2 + 4*C2p4*r^3 + 5*C2p5*r^4 + 6*C2p6*r^5 + 7*C2p7*r^6
        // d²U2/dr² = 72A/r^10 - 20B/r^6 + 2*C2p2 + 6*C2p3*r + 12*C2p4*r^2 + 20*C2p5*r^3 + 30*C2p6*r^4 + 42*C2p7*r^5
        
        double r2 = r * r;
        double r3 = r2 * r;
        double r4 = r3 * r;
        double r5 = r4 * r;
        double r6 = r5 * r;
        double r10 = r6 * r4;
        
        second_derivative = 72.0 * A / r10 - 20.0 * B / r6 +
                           2.0 * C2p2 + 6.0 * C2p3 * r + 12.0 * C2p4 * r2 +
                           20.0 * C2p5 * r3 + 30.0 * C2p6 * r4 + 42.0 * C2p7 * r5;
    }
    else if (r > rc) {
        // U3 = 0
        // d²U3/dr² = 0
        second_derivative = 0.0;
    }
    
    return second_derivative;
}

// Alternative implementation with better numerical stability for small r
double lennard_jones_energy_v3_second_derivative_stable(double r) {
    // Constants
    const double sigma = 1.0;
    const double rc = 1.86602540378444;
    const double A = -4.47124703217530;
    const double B = -90.2632082644540;
    const double C2p2 = -34486.6459222351;
    const double C2p3 = 32704.6608372787;
    const double C2p4 = -18879.7974113434;
    const double C2p5 = 6595.54667916880;
    const double C2p6 = -1286.11734180526;
    const double C2p7 = 107.717810684010;
    
    double second_derivative = 0.0;
    
    if (r < sigma) {
        // U1 = -2/r^4 + 1/r^8
        // d²U1/dr² = -40/r^6 + 72/r^10 = (72 - 40*r^4)/r^10
        double r2 = r * r;
        double r4 = r2 * r2;
        double r6 = r4 * r2;
        double r10 = r6 * r4;
        second_derivative = (-40.0 * r4 + 72.0) / r10;
    }
    else if (r >= sigma && r < rc) {
        // Same as above but using powers built incrementally
        double r2 = r * r;
        double r3 = r2 * r;
        double r4 = r3 * r;
        double r5 = r4 * r;
        double r6 = r5 * r;
        double r10 = r4 * r6;
        
        second_derivative = 72.0 * A / r10 - 20.0 * B / r6 +
                           2.0 * C2p2 + 6.0 * C2p3 * r + 12.0 * C2p4 * r2 +
                           20.0 * C2p5 * r3 + 30.0 * C2p6 * r4 + 42.0 * C2p7 * r5;
    }
    else {
        second_derivative = 0.0;
    }
    
    return second_derivative;
}

double square_energy_sder(double r) {
    const double a = 1.0;
    const double c1 = 2.0 * a;
    const double c2 = 2.0 * a;
    const double b1 = 8.0;
    const double r1 = 1.0;
    const double r2 = 1.425;

    double repulsive_sder = 156.0 * a / std::pow(r, 14);
    double dr1 = r - r1;
    double dr2 = r - r2;
    double attractive_sder1 = 2.0 * b1 * c1 * (1.0 - 2.0 * b1 * dr1 * dr1) * std::exp(-b1 * dr1 * dr1);
    double attractive_sder2 = 2.0 * b1 * c2 * (1.0 - 2.0 * b1 * dr2 * dr2) * std::exp(-b1 * dr2 * dr2);

    return repulsive_sder + attractive_sder1 + attractive_sder2;
}

#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>

void EquilibriumLatticeResult::printSummary() const {
    std::cout << "  --- Equilibrium Lattice Parameter (" << lattice_type << " lattice) ---\n";
    std::cout << std::fixed << std::setprecision(15);
    std::cout << "  Equilibrium spacing a_0 = " << lattice_parameter << "\n";
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "  Cohesive energy E_coh   = " << cohesive_energy
              << " | Energy density = " << energy_density << "\n";
    std::cout << "  Unit cell area Omega_0  = " << unit_cell_area
              << " | Curvature d^2E/da^2 = " << d2E_da2 << "\n";
    std::cout << std::scientific << std::setprecision(3);
    std::cout << "  Residual dE/da          = " << dE_da
              << " (" << iterations << " iterations, "
              << (converged ? "CONVERGED" : "NOT CONVERGED") << ")\n"
              << std::defaultfloat;
}

EquilibriumLatticeResult compute_equilibrium_lattice_parameter(
    const std::function<double(double)>& pot,
    const std::function<double(double)>& dpot,
    const std::function<double(double)>& d2pot,
    const std::string& lattice_type,
    double r_cutoff,
    double a_min,
    double a_max
) {
    EquilibriumLatticeResult result;
    result.lattice_type = lattice_type;

    // Fallback to central finite differences if analytical derivatives are not supplied
    std::function<double(double)> dpot_fn = dpot;
    if (!dpot_fn) {
        dpot_fn = [&pot](double r) {
            const double h = 1e-7;
            return (pot(r + h) - pot(r - h)) / (2.0 * h);
        };
    }

    std::function<double(double)> d2pot_fn = d2pot;
    if (!d2pot_fn) {
        d2pot_fn = [&pot, &dpot, &dpot_fn](double r) {
            const double h = 1e-6;
            if (dpot) {
                return (dpot(r + h) - dpot(r - h)) / (2.0 * h);
            }
            return (pot(r + h) - 2.0 * pot(r) + pot(r - h)) / (h * h);
        };
    }

    // Precompute dimensionless lattice neighbor distances d_{mn}
    const int nb_atoms = 10;
    std::vector<double> neighbor_dists;
    neighbor_dists.reserve((2 * nb_atoms + 1) * (2 * nb_atoms + 1));

    const bool is_triangular = (lattice_type == "triangular");
    const double sqrt3_over_2 = std::sqrt(3.0) / 2.0;

    for (int m = -nb_atoms; m <= nb_atoms; ++m) {
        for (int n = -nb_atoms; n <= nb_atoms; ++n) {
            if (m == 0 && n == 0) continue;
            double px = is_triangular ? (static_cast<double>(m) + 0.5 * static_cast<double>(n))
                                      : static_cast<double>(m);
            double py = is_triangular ? (sqrt3_over_2 * static_cast<double>(n))
                                      : static_cast<double>(n);
            double d = std::sqrt(px * px + py * py);
            neighbor_dists.push_back(d);
        }
    }

    // Helper to evaluate E(a), dE/da(a), and d^2E/da^2(a)
    auto eval_lattice_sum = [&](double a, double& E, double& dE, double& d2E) {
        E = 0.0;
        dE = 0.0;
        d2E = 0.0;
        for (double d : neighbor_dists) {
            double r = a * d;
            if (r <= r_cutoff) {
                E   += 0.5 * pot(r);
                dE  += 0.5 * d * dpot_fn(r);
                d2E += 0.5 * d * d * d2pot_fn(r);
            }
        }
    };

    // Step 1: Scan [a_min, a_max] to bracket the equilibrium root where dE/da crosses from < 0 to > 0
    const int n_scan = 250;
    const double da = (a_max - a_min) / static_cast<double>(n_scan);

    double best_aL = a_min;
    double best_aR = a_max;
    double best_E  = std::numeric_limits<double>::max();
    bool found_bracket = false;

    double prev_a = a_min;
    double prev_E = 0.0, prev_dE = 0.0, prev_d2E = 0.0;
    eval_lattice_sum(prev_a, prev_E, prev_dE, prev_d2E);

    for (int i = 1; i <= n_scan; ++i) {
        double cur_a = a_min + i * da;
        double cur_E = 0.0, cur_dE = 0.0, cur_d2E = 0.0;
        eval_lattice_sum(cur_a, cur_E, cur_dE, cur_d2E);

        if (prev_dE <= 0.0 && cur_dE >= 0.0) {
            double mid_E = 0.5 * (prev_E + cur_E);
            if (!found_bracket || mid_E < best_E) {
                best_aL = prev_a;
                best_aR = cur_a;
                best_E  = mid_E;
                found_bracket = true;
            }
        }
        prev_a  = cur_a;
        prev_E  = cur_E;
        prev_dE = cur_dE;
    }

    if (!found_bracket) {
        // Fallback: pick the lowest-energy point in the scan
        double min_E = std::numeric_limits<double>::max();
        double min_a = 1.0;
        for (int i = 0; i <= n_scan; ++i) {
            double cur_a = a_min + i * da;
            double cur_E = 0.0, cur_dE = 0.0, cur_d2E = 0.0;
            eval_lattice_sum(cur_a, cur_E, cur_dE, cur_d2E);
            if (cur_E < min_E) {
                min_E = cur_E;
                min_a = cur_a;
            }
        }
        best_aL = std::max(a_min, min_a - da);
        best_aR = std::min(a_max, min_a + da);
    }

    // Step 2: Safeguarded Newton-Raphson + bisection refinement to machine precision
    double aL = best_aL;
    double aR = best_aR;
    double EL = 0.0, dEL = 0.0, d2EL = 0.0;
    double ER = 0.0, dER = 0.0, d2ER = 0.0;
    eval_lattice_sum(aL, EL, dEL, d2EL);
    eval_lattice_sum(aR, ER, dER, d2ER);

    double a = (std::abs(dER - dEL) > 1e-16)
                   ? (aL - dEL * (aR - aL) / (dER - dEL))
                   : 0.5 * (aL + aR);
    if (a <= aL || a >= aR) {
        a = 0.5 * (aL + aR);
    }

    double E = 0.0, dE = 0.0, d2E = 0.0;
    const int max_iter = 100;
    int iter = 0;

    for (; iter < max_iter; ++iter) {
        eval_lattice_sum(a, E, dE, d2E);

        if (std::abs(dE) < 1e-14 || std::abs(aR - aL) < 1e-15) {
            result.converged = (d2E > 0.0);
            break;
        }

        if (dE < 0.0) {
            aL = a;
            dEL = dE;
        } else {
            aR = a;
            dER = dE;
        }

        // Try Newton step; fall back to bisection if outside bracket or non-convex
        double a_next = 0.5 * (aL + aR);
        if (d2E > 1e-12) {
            double a_newton = a - dE / d2E;
            if (a_newton > aL && a_newton < aR) {
                a_next = a_newton;
            }
        }

        if (std::abs(a_next - a) < 1e-16) {
            a = a_next;
            eval_lattice_sum(a, E, dE, d2E);
            result.converged = (d2E > 0.0);
            break;
        }
        a = a_next;
    }

    result.lattice_parameter = a;
    result.cohesive_energy   = E;
    result.unit_cell_area    = is_triangular ? (sqrt3_over_2 * a * a) : (a * a);
    result.energy_density    = E / result.unit_cell_area;
    result.dE_da             = dE;
    result.d2E_da2           = d2E;
    result.iterations        = iter + 1;
    return result;
}

double find_optimal_lattice_parameter(
    const std::function<double(double)>& pot,
    const std::string& lattice_type,
    double r_cutoff,
    const std::function<double(double)>& dpot,
    const std::function<double(double)>& d2pot
) {
    return compute_equilibrium_lattice_parameter(
        pot, dpot, d2pot, lattice_type, r_cutoff
    ).lattice_parameter;
}

#endif // INTER_ATOMIC_H
