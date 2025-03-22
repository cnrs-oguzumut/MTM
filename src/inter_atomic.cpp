#pragma once
#ifndef INTER_ATOMIC_H
#define INTER_ATOMIC_H
#include "interatomic/inter_atomic.h"
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
#endif // INTER_ATOMIC_H
