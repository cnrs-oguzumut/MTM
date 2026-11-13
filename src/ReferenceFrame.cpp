#include "../include/lattice_energy/ReferenceFrame.h"

#include <cmath>
#include <stdexcept>

#include "../include/experiments/common_simulation_helpers.h"
#include "../include/lattice_energy/SquareLatticeCalculator.h"
#include "../include/lattice_energy/Strain_Energy_LatticeCalculator.h"
#include "../include/lattice_energy/TriangularLatticeCalculator.h"
#include "../include/mesh/ElementTriangle2D.h"

namespace {
ReferenceFrame g_reference_frame = ReferenceFrame::SquareZ2;
bool g_reference_lattice_triangular = false;
} // namespace

bool reference_lattice_is_triangular() {
    return g_reference_lattice_triangular;
}

void set_reference_frame(ReferenceFrame frame) {
    g_reference_frame = frame;
}

ReferenceFrame reference_frame() {
    return g_reference_frame;
}

ReferenceSetup make_reference_setup(const std::string& lattice_type,
                                    int potential_mode,
                                    double a0,
                                    double r_cutoff,
                                    const std::function<double(double)>& pot) {
    g_reference_lattice_triangular = (lattice_type == "triangular");
    if (g_reference_frame == ReferenceFrame::Triangular) {
        if (lattice_type != "triangular") {
            throw std::invalid_argument(
                "Triangular reference frame requires lattice_type == \"triangular\".");
        }
        if (potential_mode == 3) {
            throw std::invalid_argument(
                "Triangular reference frame is not supported with Zanzotto potential (potential_mode == 3).");
        }

        ReferenceSetup setup;
        setup.calc_scale = a0;
        const Eigen::Vector2d p1(0.0, 0.0);
        const Eigen::Vector2d p2(a0, 0.0);
        const Eigen::Vector2d p3(0.5 * a0, a0 * std::sqrt(3.0) / 2.0);
        setup.dndx = calculateShapeDerivatives(p1, p2, p3);
        setup.element_area = ElementTriangle2D::referenceElementArea(setup.dndx);
        setup.C_ref = Eigen::Matrix2d::Identity();
        setup.F0 = Eigen::Matrix2d::Identity();
        setup.calc = std::make_unique<TriangularLatticeCalculator>(a0, r_cutoff);
        setup.zero = setup.calc->calculate_energy(setup.C_ref, pot, 0.0);
        return setup;
    }

    ReferenceSetup setup;
    setup.calc_scale = (potential_mode == 3)
        ? 1.0
        : ((lattice_type == "triangular") ? (a0 / std::pow(4.0 / 3.0, 0.25)) : a0);
    const double h = (potential_mode == 3) ? 1.0 : setup.calc_scale;
    const Eigen::Vector2d p1(0.0, 0.0);
    const Eigen::Vector2d p2(h, 0.0);
    const Eigen::Vector2d p3(0.0, h);
    setup.dndx = calculateShapeDerivatives(p1, p2, p3);
    setup.element_area = ElementTriangle2D::referenceElementArea(setup.dndx);

    if (lattice_type == "triangular") {
        const double a2 = std::sqrt(4.0 / 3.0);
        setup.C_ref << a2, 0.5 * a2,
                       0.5 * a2, a2;
        const double s0 = std::pow(4.0 / 3.0, 0.25);
        setup.F0 << s0, 0.5 * s0,
                    0.0, s0 * std::sqrt(3.0) / 2.0;
    } else {
        setup.C_ref = Eigen::Matrix2d::Identity();
        setup.F0 = Eigen::Matrix2d::Identity();
    }

    if (potential_mode == 3) {
        setup.calc = std::make_unique<Strain_Energy_LatticeCalculator>(1.0);
    } else {
        setup.calc = std::make_unique<SquareLatticeCalculator>(setup.calc_scale, r_cutoff);
    }
    setup.zero = setup.calc->calculate_energy(setup.C_ref, pot, 0.0);
    return setup;
}
