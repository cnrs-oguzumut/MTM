#pragma once
// Macroscopic nominal-stress control for periodic (PBC) cells (stress_and_force_loading_plan.md, §4).
//
// The cell deformation F_bar (upper triangular, F21 = 0) has stress-driven components that are solved
// for together with the atoms, on the enthalpy
//   H(w, F_bar) = E_el(F_bar w) - A0 * sum_{free ij} Pt_ij (F_ij - delta_ij),
// with internal coordinates w (x = F_bar w). The extra variables live only inside relax_configuration;
// outside it the pipeline keeps physical positions and the loop's F_ext (= *F_bar).

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <string>

struct PBCStressControl {
    // ---- options (command line) ----
    // "" (none), "shear" (P12 = alpha), "uniaxial_y" (P22 = alpha), "uniaxial_x" (P11 = alpha),
    // "relax_normal" (F12 = alpha strain-driven, P11 = P22 = 0)
    std::string preset;
    bool relax_stage1 = false; // zero-stress relaxation of all three components at alpha = 0

    // ---- per-step state, set by the loading loop ----
    bool active = false;                  // F_bar components are unknowns in this step
    bool comp_free[3] = {false, false, false}; // (11, 12, 22)
    Eigen::Matrix2d P_target = Eigen::Matrix2d::Zero(); // nominal target on the free components
    Eigen::Matrix2d* F_bar = nullptr;     // the loop's F_ext (updated by the relaxation)
    double A0 = 0.0;                      // reference cell area (sum of element reference areas)

    // Limit-load guard: largest change of a stress-driven component within one step
    Eigen::Matrix2d F_step_start = Eigen::Matrix2d::Identity();
    double F_jump_limit = 0.0;            // <= 0: off

    bool enabled() const { return !preset.empty() || relax_stage1; }

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
void set_pbc_stress_control_options(const std::string& preset, bool relax_stage1);
const PBCStressControl& pbc_stress_control_options();
