#pragma once
// Macroscopic nominal-stress control for periodic (PBC) cells (stress_and_force_loading_plan.md, §4).
//
// The cell deformation F_bar (upper triangular, F21 = 0) has stress-driven components that are solved
// for together with the atoms, on the enthalpy
//   H(w, F_bar) = E_el(F_bar w) - A0 * sum_{free ij} Pt_ij (F_ij - delta_ij),
// with internal coordinates w (x = F_bar w). The extra variables live only inside relax_configuration;
// outside it the pipeline keeps physical positions and the loop's F_ext (= *F_bar).
//
// Cauchy (true-stress) presets "cauchy_*" prescribe sigma instead of P. A fixed Cauchy stress is a
// follower load with no potential in general, so each relaxation is a fixed-point loop on the nominal
// target conjugate to the free components, Pt = J sigma_t F_bar^{-T}, re-evaluated at the relaxed box.
// Matching P11, P12, P22 gives the full sigma = sigma_t: rotational invariance makes P F^T symmetric.

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <string>

struct PBCStressControl {
    // ---- options (command line) ----
    // "" (none), "shear" (P12 = alpha), "uniaxial_y" (P22 = alpha), "uniaxial_x" (P11 = alpha),
    // "relax_normal" (F12 = alpha strain-driven, P11 = P22 = 0), and the Cauchy versions
    // "cauchy_shear" (sigma_xy = alpha), "cauchy_uniaxial_y" (sigma_yy = alpha), "cauchy_uniaxial_x"
    // (sigma_xx = alpha); the other in-plane stress components are zero in every stress preset
    std::string preset;
    bool relax_stage1 = false; // zero-stress relaxation of all three components at alpha = 0
    // Shear presets: the relaxation step (step 0) is held at F12 = relax_strain under strain control (affine
    // + noise, remesh cascade), so dislocations nucleate on the unstable branch; the stress schedule then
    // rises from the relaxed stress
    double relax_strain = 0.0;

    // ---- per-step state, set by the loading loop ----
    bool active = false;                  // F_bar components are unknowns in this step
    bool comp_free[3] = {false, false, false}; // (11, 12, 22)
    Eigen::Matrix2d P_target = Eigen::Matrix2d::Zero(); // nominal target on the free components
    bool cauchy = false;                  // P_target follows sigma_target (fixed-point loop)
    Eigen::Matrix2d sigma_target = Eigen::Matrix2d::Zero(); // symmetric Cauchy target
    Eigen::Matrix2d* F_bar = nullptr;     // the loop's F_ext (updated by the relaxation)
    double A0 = 0.0;                      // reference cell area (sum of element reference areas)

    // Limit-load guard: largest change of a stress-driven component within one step
    Eigen::Matrix2d F_step_start = Eigen::Matrix2d::Identity();
    double F_jump_limit = 0.0;            // <= 0: off

    bool enabled() const { return !preset.empty() || relax_stage1; }

    static bool isCauchyPreset(const std::string& p) { return p.rfind("cauchy_", 0) == 0; }
    // Loading pattern without the stress measure: "cauchy_shear" -> "shear"
    static std::string basePreset(const std::string& p) { return isCauchyPreset(p) ? p.substr(7) : p; }
    static bool isForcePreset(const std::string& p) {
        const std::string b = basePreset(p);
        return b == "shear" || b == "uniaxial_y" || b == "uniaxial_x";
    }
    static bool isKnownPreset(const std::string& p) {
        return isForcePreset(p) || p == "relax_normal";
    }
    // Controlled (row, col) of a stress preset: shear (0, 1), uniaxial_y (1, 1), uniaxial_x (0, 0)
    static int targetRow(const std::string& p) { return basePreset(p) == "uniaxial_y" ? 1 : 0; }
    static int targetCol(const std::string& p) { return basePreset(p) == "uniaxial_x" ? 0 : 1; }

    // Pt = J sigma_t F^{-T} at box F
    void setNominalFromCauchy(const Eigen::Matrix2d& F) {
        P_target = F.determinant() * sigma_target * F.inverse().transpose();
    }

    static int compRow(int k) { return k == 2 ? 1 : 0; }
    static int compCol(int k) { return k == 0 ? 0 : 1; }

    // -A0 * sum_free Pt_ij (F_ij - delta_ij)
    double workTerm(const Eigen::Matrix2d& F) const {
        if (!active) return 0.0;
        double w = 0.0;
        for (int k = 0; k < 3; ++k) {
            if (!comp_free[k]) continue;
            const int r = compRow(k), c = compCol(k);
            w -= A0 * P_target(r, c) * (F(r, c) - (r == c ? 1.0 : 0.0));
        }
        return w;
    }

    double maxStressDrivenJump(const Eigen::Matrix2d& F) const {
        double jump = 0.0;
        for (int k = 0; k < 3; ++k) {
            if (comp_free[k]) jump = std::max(jump, std::abs(F(compRow(k), compCol(k)) -
                                                             F_step_start(compRow(k), compCol(k))));
        }
        return jump;
    }
    bool exceedsLimit(const Eigen::Matrix2d& F) const {
        return active && F_jump_limit > 0.0 && maxStressDrivenJump(F) > F_jump_limit;
    }
};

// Options from the command line, read by the PBC loading loop
void set_pbc_stress_control_options(const std::string& preset, bool relax_stage1, double relax_strain = 0.0);
const PBCStressControl& pbc_stress_control_options();
