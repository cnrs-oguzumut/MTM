#ifndef ENERGY_SPLINE_INTERPOLATOR_H
#define ENERGY_SPLINE_INTERPOLATOR_H

#include <Eigen/Dense>
#include "src/interpolation.h"  // Alglib
#include <string>
#include <functional>

class EnergySplineInterpolator {
public:
    // Constructor
    EnergySplineInterpolator();
    
    // Result structure for energy and derivatives
    struct EnergyResult {
        double energy;
        double dE_dc11;
        double dE_dc22;
        double dE_dc12;
        bool valid;
        
        EnergyResult() : energy(0.0), dE_dc11(0.0), dE_dc22(0.0), 
                         dE_dc12(0.0), valid(false) {}
    };
    
    // Internal parameter structure (public for debugging)
    struct StrainParams {
        double t;
        double p;
        bool valid;
    };
    
    // Domain constraints configuration
    struct DomainConstraints {
        bool enforce_lagrange_condition;  // Enforce 2*|c12| < min(c11, c22)
        bool c12_positive_only;           // Only c12 > 0
        bool auto_compute_domain;         // Auto-compute safe (t,p) domain
        
        DomainConstraints() 
            : enforce_lagrange_condition(true), 
              c12_positive_only(true),
              auto_compute_domain(true) {}
    };
    
    // Build the spline from parametric data
    void buildSpline(
        std::function<double(double)> potential_func,
        std::function<double(double)> potential_func_der,
        std::function<double(double)> potential_func_sder,
        double scale = 1.0,
        int nt = 100, 
        int np = 256,
        const DomainConstraints& constraints = DomainConstraints()
    );
    
    // Evaluate energy and derivatives from C matrix components
    EnergyResult evaluate(double c11, double c12, double c22) const;
    
    // Evaluate from Eigen matrix
    EnergyResult evaluateFromMatrix(const Eigen::Matrix2d& C) const;
    
    // Save/Load spline
    void saveSpline(const std::string& filename) const;
    void loadSpline(const std::string& filename);
    
    // Check if initialized
    bool isInitialized() const { return is_initialized; }
    
    // Get domain bounds (useful for testing)
    void getDomainBounds(double& t_min_out, double& t_max_out,
                         double& p_min_out, double& p_max_out) const {
        t_min_out = t_min;
        t_max_out = t_max;
        p_min_out = p_min;
        p_max_out = p_max;
    }
    
    // Inverse mapping: C -> (t, p) - public for debugging
    StrainParams C_to_tp(double c11, double c12, double c22) const;
    
private:
    alglib::spline2dinterpolant energy_spline;
    alglib::spline2dinterpolant dE_dc11_spline;
    alglib::spline2dinterpolant dE_dc22_spline;
    alglib::spline2dinterpolant dE_dc12_spline;
    double t_min, t_max, p_min, p_max;
    bool is_initialized;
    DomainConstraints constraints_;
    
    // Lagrange condition checking
    bool satisfiesLagrangeCondition(double c11, double c12, double c22) const;
    
    // Compute maximum valid t for given p (Lagrange constraint)
    static double compute_max_t_for_p(double p);
    
    // Find safe (t,p) domain satisfying Lagrange condition
    static std::pair<double, double> find_safe_tp_domain(
        double p_min, double p_max, int n_samples = 1000);
};

#endif // ENERGY_SPLINE_INTERPOLATOR_H