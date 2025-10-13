#pragma once

// Remove the circular include if this IS inter_atomic.h
// If you need another header, include it with a different name
// #include "other_header.h"  // Only if needed

double lennard_jones_energy_v2(double r);
double lennard_jones_energy_der_v2(double r);
double lennard_jones_energy_sder_v2(double r);

double square_energy(double r);
double square_energy_der(double r);
// double square_energy_sder(double r);  // Uncomment when implemented

double lennard_jones_energy_v3(double r);
double lennard_jones_energy_der_v3(double r);
double lennard_jones_energy_sder_v3(double r);