#include "../include/EnergySplineInterpolator.h"
#include "../include/reductions/LagrangeReduction.h"
#include "../include/lattice_energy/Strain_Energy_LatticeCalculator.h"
#include "../include/lattice_energy/EnergyFunctions.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

EnergySplineInterpolator::EnergySplineInterpolator() 
    : is_initialized(false), t_min(0.0), t_max(1.0), 
      p_min(-M_PI), p_max(M_PI) {
}

bool EnergySplineInterpolator::satisfiesLagrangeCondition(
    double c11, double c12, double c22) const {
    double min_val = std::min(c11, c22);
    return (2.0 * fabs(c12) < min_val);
}

void EnergySplineInterpolator::buildSpline(
    std::function<double(double)> potential_func,
    std::function<double(double)> potential_func_der,
    std::function<double(double)> potential_func_sder,
    double scale,
    int nt,
    int np,
    const DomainConstraints& constraints
) {
    std::cout << "\n==================================================" << std::endl;
    std::cout << "Building 2D energy spline (Fixed Large Domain)" << std::endl;
    std::cout << "==================================================" << std::endl;
    
    constraints_ = constraints;
    
    t_min = 0.0;
    t_max = 1.0;
    p_min = -M_PI;
    p_max = M_PI;
    
    std::cout << "Domain Configuration:" << std::endl;
    std::cout << "  t ∈ [" << t_min << ", " << t_max << "]" << std::endl;
    std::cout << "  p ∈ [" << p_min << ", " << p_max << "]" << std::endl;
    std::cout << "  Grid: " << nt << " × " << np << " = " << nt*np << " points" << std::endl;
    
    Strain_Energy_LatticeCalculator strain_calculator(scale);
    double normalisation = strain_calculator.getUnitCellArea();
    
    alglib::real_1d_array t_grid, p_grid;
    t_grid.setlength(nt);
    p_grid.setlength(np);
    
    double dt = (t_max - t_min) / std::max(1, nt - 1);
    double dp = (p_max - p_min) / std::max(1, np - 1);
    
    for (int i = 0; i < nt; i++) {
        t_grid[i] = t_min + i * dt;
    }
    for (int j = 0; j < np; j++) {
        p_grid[j] = p_min + j * dp;
    }
    
    alglib::real_1d_array energy_grid_1d;
    alglib::real_1d_array dE_dc11_grid;
    alglib::real_1d_array dE_dc22_grid;
    alglib::real_1d_array dE_dc12_grid;
    
    energy_grid_1d.setlength(nt * np);
    dE_dc11_grid.setlength(nt * np);
    dE_dc22_grid.setlength(nt * np);
    dE_dc12_grid.setlength(nt * np);
    
    // ADD THESE TWO LINES
    int count = 0;
    int invalid_count = 0;
    
    std::cout << "\nComputing energy grid..." << std::endl;
    
    for (int i = 0; i < nt; i++) {
        for (int j = 0; j < np; j++) {
            double t = t_grid[i];
            double p = p_grid[j];
            
            int idx = j * nt + i;
            
            try {
                double c11 = cosh(t) + sinh(t) * sin(p);
                double c22 = cosh(t) - sinh(t) * sin(p);
                double c12 = sinh(t) * cos(p);
                
                Eigen::Matrix2d C;
                C << c11, c12, c12, c22;
                
                lagrange::Result reduction_result = lagrange::reduce(C);
                Eigen::Matrix2d C_reduced = reduction_result.C_reduced;
                
                double energy = strain_calculator.calculate_energy(C_reduced, potential_func, normalisation);
                energy_grid_1d[idx] = energy;
                
                Eigen::Matrix2d dE = energy_functions::dphi_func(C_reduced);
                dE_dc11_grid[idx] = dE(0, 0);
                dE_dc22_grid[idx] = dE(1, 1);
                dE_dc12_grid[idx] = dE(0, 1);
                
                count++;
                
            } catch (const std::exception& e) {
                energy_grid_1d[idx] = 1e10;
                dE_dc11_grid[idx] = 0.0;
                dE_dc22_grid[idx] = 0.0;
                dE_dc12_grid[idx] = 0.0;
                invalid_count++;
            }
        }
        
        if ((i + 1) % 20 == 0) {
            std::cout << "Progress: " << (i + 1) << "/" << nt 
                      << " rows (" << count << " valid so far)" << std::endl;
        }
    }
    
    std::cout << "\nBuilding bicubic splines (energy + 3 derivatives)..." << std::endl;
    alglib::spline2dbuildbicubicv(t_grid, nt, p_grid, np, energy_grid_1d, 1, energy_spline);
    alglib::spline2dbuildbicubicv(t_grid, nt, p_grid, np, dE_dc11_grid, 1, dE_dc11_spline);
    alglib::spline2dbuildbicubicv(t_grid, nt, p_grid, np, dE_dc22_grid, 1, dE_dc22_spline);
    alglib::spline2dbuildbicubicv(t_grid, nt, p_grid, np, dE_dc12_grid, 1, dE_dc12_spline);
    
    is_initialized = true;
    
    std::cout << "\n==================================================" << std::endl;
    std::cout << "Spline construction complete!" << std::endl;
    std::cout << "==================================================" << std::endl;
    std::cout << "Computed points: " << count << "/" << (nt*np) << std::endl;
    std::cout << "Invalid points:  " << invalid_count << " (set to high energy)" << std::endl;
    std::cout << "==================================================" << std::endl;
}

EnergySplineInterpolator::StrainParams 
EnergySplineInterpolator::C_to_tp(double c11, double c12, double c22) const {
    StrainParams params;
    params.valid = true;
    
    double cosh_t = (c11 + c22) / 2.0;
    
    if (cosh_t < 1.0 - 1e-10) {
        params.valid = false;
        return params;
    }
    
    if (cosh_t < 1.0) cosh_t = 1.0;
    
    if (cosh_t < 1.0001) {
        params.t = sqrt(2.0 * (cosh_t - 1.0));
    } else {
        params.t = acosh(cosh_t);
    }
    
    params.p = atan2((c11 - c22) / 2.0, c12);
    
    if (params.t > t_max) { 
        params.valid = false;
    }
    
    return params;
}

EnergySplineInterpolator::EnergyResult 
EnergySplineInterpolator::evaluate(double c11, double c12, double c22) const {
    EnergyResult result;
    
    if (!is_initialized) {
        std::cerr << "Error: Spline not initialized!" << std::endl;
        return result;
    }
    
    StrainParams tp = C_to_tp(c11, c12, c22);
    
    if (!tp.valid) {
        result.valid = false;
        return result;
    }
    
    double t = tp.t;
    double p = tp.p;
    
    result.energy = alglib::spline2dcalc(energy_spline, t, p);
    result.dE_dc11 = alglib::spline2dcalc(dE_dc11_spline, t, p);
    result.dE_dc22 = alglib::spline2dcalc(dE_dc22_spline, t, p);
    result.dE_dc12 = alglib::spline2dcalc(dE_dc12_spline, t, p);
    
    result.valid = true;
    return result;
}

EnergySplineInterpolator::EnergyResult 
EnergySplineInterpolator::evaluateFromMatrix(const Eigen::Matrix2d& C) const {
    return evaluate(C(0,0), C(0,1), C(1,1));
}

void EnergySplineInterpolator::saveSpline(const std::string& filename) const {
    if (!is_initialized) return;
    alglib::spline2dinterpolant& non_const_spline = 
        const_cast<alglib::spline2dinterpolant&>(energy_spline);
    std::string s;
    alglib::spline2dserialize(non_const_spline, s);
    std::ofstream file(filename, std::ios::binary);
    file.write(s.c_str(), s.size());
    file.close();
    std::cout << "Spline saved to: " << filename << std::endl;
}

void EnergySplineInterpolator::loadSpline(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    std::string s((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    file.close();
    alglib::spline2dunserialize(s, energy_spline);
    is_initialized = true;
    std::cout << "Spline loaded from: " << filename << std::endl;
}