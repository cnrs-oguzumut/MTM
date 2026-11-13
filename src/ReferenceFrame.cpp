#include "../include/lattice_energy/ReferenceFrame.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
#include <stdexcept>

#include "../include/experiments/common_simulation_helpers.h"
#include "../include/lattice_energy/SquareLatticeCalculator.h"
#include "../include/lattice_energy/Strain_Energy_LatticeCalculator.h"
#include "../include/lattice_energy/TriangularLatticeCalculator.h"
#include "../include/mesh/ElementTriangle2D.h"

namespace {
ReferenceFrame g_reference_frame = ReferenceFrame::SquareZ2;
bool g_reference_lattice_triangular = false;

// NOTE: compiled with -ffast-math; no infinity()/NaN sentinels here.
struct GroundState {
    double energy = 0.0;
    Eigen::Matrix2d C_shape = Eigen::Matrix2d::Identity();  // Z^2 metric shape (det 1)
    bool at_search_edge = false;
};

// Lowest energy of the calculator over all 2D Bravais lattices. Shapes: reduced Z^2 metrics with det 1
// (0 <= 2 c12 <= c11 <= c22), plus the exact square and triangular shapes; for each shape the cell size is
// minimised (coarse scan + golden section). The triangular-frame calculator works in triangular lattice
// coordinates, so Z^2 metrics are converted with C = B^-T C_z B^-1.
GroundState find_ground_state(BaseLatticeCalculator& calc, const std::function<double(double)>& pot,
                              bool triangular_coordinates) {
    Eigen::Matrix2d B;
    B << 1.0, 0.5,
         0.0, std::sqrt(3.0) / 2.0;
    const Eigen::Matrix2d B_inv = B.inverse();
    auto energy_at = [&](const Eigen::Matrix2d& Cz, double s) {
        const Eigen::Matrix2d C = triangular_coordinates ? Eigen::Matrix2d(B_inv.transpose() * Cz * B_inv) : Cz;
        return calc.calculate_energy(s * s * C, pot, 0.0);
    };
    const double s_lo = 0.5, s_hi = 2.0;
    auto best_over_size = [&](const Eigen::Matrix2d& Cz, bool& edge) {
        const int n_scan = 30;
        int k_best = 0;
        double e_best = std::numeric_limits<double>::max();
        for (int k = 0; k <= n_scan; ++k) {
            const double e = energy_at(Cz, s_lo + (s_hi - s_lo) * k / n_scan);
            if (e < e_best) { e_best = e; k_best = k; }
        }
        edge = (k_best == 0 || k_best == n_scan);
        const double h = (s_hi - s_lo) / n_scan;
        double a = s_lo + h * std::max(0, k_best - 1), b = s_lo + h * std::min(n_scan, k_best + 1);
        const double g = 0.5 * (std::sqrt(5.0) - 1.0);
        double x1 = b - g * (b - a), x2 = a + g * (b - a);
        double f1 = energy_at(Cz, x1), f2 = energy_at(Cz, x2);
        for (int it = 0; it < 60; ++it) {
            if (f1 < f2) { b = x2; x2 = x1; f2 = f1; x1 = b - g * (b - a); f1 = energy_at(Cz, x1); }
            else         { a = x1; x1 = x2; f1 = f2; x2 = a + g * (b - a); f2 = energy_at(Cz, x2); }
        }
        return std::min(e_best, std::min(f1, f2));
    };

    std::vector<Eigen::Matrix2d> shapes;
    const double a_tri = std::pow(4.0 / 3.0, 0.25);
    Eigen::Matrix2d C_tri;
    C_tri << a_tri * a_tri, 0.5 * a_tri * a_tri,
             0.5 * a_tri * a_tri, a_tri * a_tri;
    shapes.push_back(Eigen::Matrix2d::Identity());
    shapes.push_back(C_tri);
    const int n_a = 24, n_b = 24;
    const double a_min = 0.45;
    for (int i = 0; i <= n_a; ++i) {
        const double c11 = a_min + (a_tri - a_min) * i / n_a;
        for (int j = 0; j <= n_b; ++j) {
            const double c12 = 0.5 * c11 * j / n_b;
            const double c22 = (1.0 + c12 * c12) / c11;
            if (c22 < c11 - 1e-12) continue;
            Eigen::Matrix2d Cz;
            Cz << c11, c12,
                  c12, c22;
            shapes.push_back(Cz);
        }
    }

    GroundState gs;
    gs.energy = std::numeric_limits<double>::max();
    for (const Eigen::Matrix2d& Cz : shapes) {
        bool edge = false;
        const double e = best_over_size(Cz, edge);
        if (e < gs.energy) {
            gs.energy = e;
            gs.C_shape = Cz;
            gs.at_search_edge = edge || (Cz(0, 0) < a_min + 1e-9);
        }
    }
    return gs;
}

// Sets setup.zero to the potential's ground state and setup.reference_excess to the reference lattice's
// energy above it. When the reference lattice is itself the ground state the zero is its exact energy
// (perfect reference lattice = exactly 0, unchanged from before).
void set_energy_zero(ReferenceSetup& setup, const std::function<double(double)>& pot,
                     bool triangular_coordinates, const std::string& lattice_type) {
    const double e_ref = setup.calc->calculate_energy(setup.C_ref, pot, 0.0);
    const GroundState gs = find_ground_state(*setup.calc, pot, triangular_coordinates);
    const double tol = 1e-9 * std::max(1.0, std::abs(e_ref));
    if (e_ref - gs.energy <= tol) {
        setup.zero = e_ref;
        setup.reference_excess = 0.0;
        return;
    }
    setup.zero = gs.energy;
    setup.reference_excess = e_ref - gs.energy;
    const Eigen::Matrix2d& C = gs.C_shape;
    const char* kind = (std::abs(C(0, 1)) < 1e-9 && std::abs(C(0, 0) - C(1, 1)) < 1e-9) ? "square"
                     : (std::abs(C(0, 1) - 0.5 * C(0, 0)) < 1e-9 && std::abs(C(0, 0) - C(1, 1)) < 1e-9) ? "triangular"
                     : "oblique";
    std::cout << "[Reference] Energy zero set at the potential's ground state (" << kind << " lattice). The "
              << lattice_type << " reference lattice lies " << setup.reference_excess
              << " above it (energy per atom, before area normalisation), so a perfect " << lattice_type
              << " crystal no longer reads 0.\n";
    if (gs.at_search_edge) {
        std::cout << "[Reference] WARNING: the lowest energy was found at the edge of the lattice search "
                     "(very elongated cell or extreme size): the energy may be unbounded below for these "
                     "potential parameters.\n";
    }
}
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
        set_energy_zero(setup, pot, /*triangular_coordinates=*/true, lattice_type);
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
    set_energy_zero(setup, pot, /*triangular_coordinates=*/false, lattice_type);
    return setup;
}
