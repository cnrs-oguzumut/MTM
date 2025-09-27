#pragma once
#ifndef INTER_ATOMIC_H
#define INTER_ATOMIC_H
#include "../interatomic/inter_atomic.h"


double lennard_jones_energy_v2(double r);
double lennard_jones_energy_der_v2(double r);
double lennard_jones_energy_sder_v2(double r);

double square_energy(double r);
double square_energy_der(double r);
double square_energy_sder(double r);

double lennard_jones_energy_v3(double r);
double lennard_jones_energy_der_v3(double r);
double lennard_jones_energy_sder_v3(double r);

#endif // INTER_ATOMIC_H
