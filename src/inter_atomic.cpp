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


#endif // INTER_ATOMIC_H
