#ifndef REFERENCE_FRAME_H
#define REFERENCE_FRAME_H

#include <Eigen/Dense>
#include <functional>
#include <memory>
#include <string>

#include "BaseLatticeCalculator.h"
#include "../reductions/LagrangeReduction.h"

enum class ReferenceFrame {
    SquareZ2,
    Triangular
};

void set_reference_frame(ReferenceFrame frame);
ReferenceFrame reference_frame();

// Lattice symmetry of the run (triangular or square), recorded by make_reference_setup so output-side
// lattice measures (orientation, rotation) need no separate setup call
bool reference_lattice_is_triangular();

inline bool use_lagrange_reduction() {
    return reference_frame() == ReferenceFrame::SquareZ2;
}

inline lagrange::Result reduce_metric(const Eigen::Matrix2d& C) {
    if (!use_lagrange_reduction()) {
        return lagrange::Result{C, Eigen::Matrix2d::Identity(), false};
    }
    return lagrange::reduce(C);
}

inline lagrange::ResultWithCount reduce_metric_and_count(const Eigen::Matrix2d& C) {
    if (!use_lagrange_reduction()) {
        return lagrange::ResultWithCount{C, Eigen::Matrix2d::Identity(), false, 0};
    }
    return lagrange::reduceAndCount(C);
}

struct ReferenceSetup {
    double calc_scale = 1.0;
    Eigen::Matrix<double, 3, 2> dndx = Eigen::Matrix<double, 3, 2>::Zero();
    double element_area = 0.5;
    Eigen::Matrix2d C_ref = Eigen::Matrix2d::Identity();
    Eigen::Matrix2d F0 = Eigen::Matrix2d::Identity();
    std::unique_ptr<BaseLatticeCalculator> calc;
    double zero = 0.0;
};

ReferenceSetup make_reference_setup(const std::string& lattice_type,
                                    int potential_mode,
                                    double a0,
                                    double r_cutoff,
                                    const std::function<double(double)>& pot);

#endif // REFERENCE_FRAME_H
