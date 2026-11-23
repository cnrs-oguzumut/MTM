#include <iostream>
#include <fstream>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <array>
#include <vector>
#include "../include/optimization/LatticeOptimizer.h"
#include "../include/geometry/LatticeGenerator.h"
#include "../include/geometry/DomainDimensions.h"
#include "../include/geometry/MacroscopicGeometry.h"
#include "../include/interatomic/inter_atomic.h"
#include "../include/lattice_energy/EnergyFunctions.h"
#include "../include/lattice_energy/SquareLatticeCalculator.h"
#include "../include/mesh/MeshGenerator.h"
#include "../include/mesh/Remesher.h"
#include "../include/output/configuration_saver.h"

static void test_polygon_ray_casting_and_cutouts() {
    MacroscopicGeometryConfig geo;
    geo.enabled = true;
    geo.outer_shape = OuterDomainShapeType::Polygon;
    // Diamond outer polygon with both upward-sloping and downward-sloping edges
    geo.outer_polygon = {
        {0.50, 0.05},
        {0.95, 0.50},
        {0.50, 0.95},
        {0.05, 0.50}
    };

    std::vector<Point2D> pts = LatticeGenerator::generate_2d_lattice(21, 21, 1.0, "square");
    geo.setPhysicalBoundsFromPoints(pts, 1.0);

    // Center (10, 10) must be inside the diamond
    assert(geo.isPointInsideSolid(10.0, 10.0));
    // Points near all four quadrants inside the diamond (testing both upward & downward edges)
    assert(geo.isPointInsideSolid(12.0, 10.0));
    assert(geo.isPointInsideSolid(8.0, 10.0));
    assert(geo.isPointInsideSolid(10.0, 12.0));
    assert(geo.isPointInsideSolid(10.0, 8.0));
    // Corners of the bounding box must be outside the diamond
    assert(!geo.isPointInsideSolid(1.0, 1.0));
    assert(!geo.isPointInsideSolid(19.0, 1.0));
    assert(!geo.isPointInsideSolid(19.0, 19.0));
    assert(!geo.isPointInsideSolid(1.0, 19.0));

    // Add a triangular polygon cutout in the center
    GeometryCutout tri_hole;
    tri_hole.enabled = true;
    tri_hole.type = CutoutShapeType::Polygon;
    tri_hole.vertices = {
        {0.40, 0.40},
        {0.60, 0.40},
        {0.50, 0.60}
    };
    geo.cutouts.push_back(tri_hole);

    // Point (10, 9.5) -> normalized (0.50, 0.475) is inside the triangular hole -> not solid
    assert(!geo.isPointInsideSolid(10.0, 9.5));
    // Point (10, 7.8) -> normalized (0.50, 0.39) is just outside the bottom edge y=0.40 of the hole
    assert(geo.isPointInsideSolid(10.0, 7.8));
    assert(!geo.isPointNearCutout(10.0, 7.8, 0.0));
    // With margin_phys = 0.5 (normalized 0.025), (10, 7.8) IS near the polygon cutout!
    assert(geo.isPointNearCutout(10.0, 7.8, 0.5));

    // Triangle rejection across polygon cutout edge
    assert(geo.shouldRejectTriangle(
        Eigen::Vector2d(7.5, 9.0),
        Eigen::Vector2d(12.5, 9.0),
        Eigen::Vector2d(10.0, 13.0)
    ));

    // Test carving and OuterPerimeter + HoleSurfaces BCs on polygon
    BoundaryZoneRule perim_rule;
    perim_rule.enabled = true;
    perim_rule.region = BoundaryRegionType::OuterPerimeter;
    perim_rule.kind = BoundaryConditionKind::FixedClamped;
    perim_rule.thickness_a0 = 1.5;
    geo.boundary_rules = { perim_rule };

    std::vector<Point2D> carved = geo.carveLatticePoints(pts, 1.0);
    assert(!carved.empty());
    assert(carved.size() < pts.size());

    std::vector<std::pair<int, int>> interior_mapping, full_mapping;
    std::vector<int> node_roles;
    std::vector<Eigen::Vector2d> node_vel;
    geo.buildBoundaryDOFs(carved, interior_mapping, full_mapping, node_roles, node_vel);
    int fixed_nodes = 0;
    for (int r : node_roles) if (r == 1) fixed_nodes++;
    assert(fixed_nodes > 0 && fixed_nodes < static_cast<int>(carved.size()));
    std::cout << "[PASS] Polygon ray-casting, polygon cutouts, and polygon BCs verified." << std::endl;
}

static void test_nanoindenter_shapes_and_penalty() {
    std::vector<Point2D> pts = LatticeGenerator::generate_2d_lattice(15, 15, 1.0, "square");

    for (IndenterShapeType shape : {IndenterShapeType::Cylinder, IndenterShapeType::FlatPunch, IndenterShapeType::Wedge}) {
        MacroscopicGeometryConfig geo;
        geo.enabled = true;
        geo.outer_shape = OuterDomainShapeType::Rectangle;
        geo.indenter.enabled = true;
        geo.indenter.shape = shape;
        geo.indenter.contact_mode = IndenterContactMode::RepulsivePenalty;
        geo.indenter.cx = 0.5;
        geo.indenter.radius_a0 = 4.0;
        geo.indenter.max_depth_a0 = 2.0;
        geo.indenter.penalty_stiffness = 120.0;
        geo.indenter.wedge_half_angle_deg = 65.0;

        BoundaryZoneRule bot_rule;
        bot_rule.enabled = true;
        bot_rule.region = BoundaryRegionType::BottomEdge;
        bot_rule.kind = BoundaryConditionKind::FixedClamped;
        bot_rule.thickness_a0 = 1.5;
        // Also add a LeftEdge FixedClamped rule extending into upper 60% to verify indenter does NOT overwrite it
        BoundaryZoneRule left_rule;
        left_rule.enabled = true;
        left_rule.region = BoundaryRegionType::LeftEdge;
        left_rule.kind = BoundaryConditionKind::FixedClamped;
        left_rule.thickness_a0 = 1.5;
        geo.boundary_rules = { bot_rule, left_rule };

        std::vector<Point2D> carved = geo.carveLatticePoints(pts, 1.0);
        geo.setupIndenterBounds(1.0, &carved);
        assert(std::abs(geo.indenter.physical_y_top - 14.0) < 1e-9);

        std::vector<std::pair<int, int>> interior_mapping, full_mapping;
        std::vector<int> baseline_roles;
        std::vector<Eigen::Vector2d> baseline_vel;
        geo.buildBoundaryDOFs(carved, interior_mapping, full_mapping, baseline_roles, baseline_vel);

        // Verify evaluateIndenterContact derivatives via finite differences at a point under the tip
        Eigen::Vector2d test_pt(7.2, 13.5);
        double alpha = 0.5;
        double kc = geo.indenter.penalty_stiffness;
        auto c0 = geo.evaluateIndenterContact(test_pt, alpha, 0.0);
        assert(c0.in_contact && c0.penetration > 0.0);

        double eps = 1e-6;
        auto cx_plus  = geo.evaluateIndenterContact(test_pt + Eigen::Vector2d(eps, 0.0), alpha, 0.0);
        auto cx_minus = geo.evaluateIndenterContact(test_pt - Eigen::Vector2d(eps, 0.0), alpha, 0.0);
        auto cy_plus  = geo.evaluateIndenterContact(test_pt + Eigen::Vector2d(0.0, eps), alpha, 0.0);
        auto cy_minus = geo.evaluateIndenterContact(test_pt - Eigen::Vector2d(0.0, eps), alpha, 0.0);

        auto pen_energy = [&](const MacroscopicGeometryConfig::IndenterPointContact& c) {
            return c.in_contact ? 0.5 * kc * c.penetration * c.penetration : 0.0;
        };
        auto pen_grad = [&](const MacroscopicGeometryConfig::IndenterPointContact& c) -> Eigen::Vector2d {
            if (!c.in_contact) return Eigen::Vector2d(0.0, 0.0);
            return Eigen::Vector2d(-kc * c.penetration * c.normal);
        };

        Eigen::Vector2d g0 = pen_grad(c0);
        double fd_gx = (pen_energy(cx_plus) - pen_energy(cx_minus)) / (2.0 * eps);
        double fd_gy = (pen_energy(cy_plus) - pen_energy(cy_minus)) / (2.0 * eps);
        assert(std::abs(g0.x() - fd_gx) < 1e-4);
        assert(std::abs(g0.y() - fd_gy) < 1e-4);

        Eigen::Vector2d fd_hx = (pen_grad(cx_plus) - pen_grad(cx_minus)) / (2.0 * eps);
        Eigen::Vector2d fd_hy = (pen_grad(cy_plus) - pen_grad(cy_minus)) / (2.0 * eps);
        assert((c0.hessian.col(0) - fd_hx).norm() < 1e-3);
        assert((c0.hessian.col(1) - fd_hy).norm() < 1e-3);

        // Verify accumulateIndenterPenalty reads trial vector x directly
        int n_vars = static_cast<int>(interior_mapping.size());
        alglib::real_1d_array x_trial, grad;
        x_trial.setlength(2 * n_vars);
        grad.setlength(2 * n_vars);
        for (int i = 0; i < n_vars; ++i) {
            int orig_idx = interior_mapping[i].first;
            x_trial[i] = carved[orig_idx].coord.x();
            x_trial[n_vars + i] = carved[orig_idx].coord.y();
            grad[i] = 0.0;
            grad[n_vars + i] = 0.0;
        }
        double e_pen = geo.accumulateIndenterPenalty(x_trial, interior_mapping, grad, alpha);
        assert(e_pen > 0.0);

        // Verify KinematicProjection preserves baseline fixed BCs in upper 60%
        geo.indenter.contact_mode = IndenterContactMode::KinematicProjection;
        std::vector<Point2D> kin_pts = carved;
        std::vector<int> node_roles = baseline_roles;
        std::vector<int> contact_atoms;
        geo.updateIndenterContactKinematics(
            kin_pts, carved, baseline_roles, node_roles, contact_atoms, alpha
        );
        assert(!contact_atoms.empty());
        for (size_t i = 0; i < baseline_roles.size(); ++i) {
            if (baseline_roles[i] == 1) {
                assert(node_roles[i] == 1);
            }
        }
    }
    std::cout << "[PASS] Nanoindenter (Cylinder, FlatPunch, Wedge) kinematics, penalty & Hessian verified." << std::endl;
}

static void test_periodic_dislocations_setup() {
    MacroscopicGeometryConfig config;
    config.enabled = true;
    config.microstructure.mode = MicrostructureMode::PeriodicDislocations;
    config.microstructure.dislocation_pair_count = 2;
    config.microstructure.dislocation_layout = "custom";
    config.microstructure.enforce_periodic_bc = true;

    // 4 glide planes: rows 10, 20, 30, 40 with signs +1, -1, +1, -1
    config.microstructure.dislocations = {
        {0, 10, 0.35, +1, true},
        {1, 20, 0.65, -1, true},
        {2, 30, 0.40, +1, true},
        {3, 40, 0.70, -1, true}
    };

    int nx = 40;
    int ny = 50;
    double a0 = 1.0;
    
    // Test square lattice
    auto [sq_points, sq_domain] = config.generatePeriodicDislocationLattice("square", nx, ny, a0, 3);
    assert(!sq_points.empty());
    
    // Total atoms:
    // Rows 0..10: 11 rows with 40 atoms = 440
    // Rows 11..20: 10 rows with 41 atoms = 410
    // Rows 21..30: 10 rows with 40 atoms = 400
    // Rows 31..40: 10 rows with 41 atoms = 410
    // Rows 41..49: 9 rows with 40 atoms = 360
    // Total = 440 + 410 + 400 + 410 + 360 = 2020 atoms
    assert(sq_points.size() == 2020);
    assert(std::abs(sq_domain.get_width() - (nx - 1) * a0) < 1e-6);
    assert(std::abs(sq_domain.get_height() - (ny - 1) * a0) < 1e-6);

    // Test triangular lattice
    auto [tri_points, tri_domain] = config.generatePeriodicDislocationLattice("triangular", nx, ny, a0, 3);
    assert(tri_points.size() == 2020);
    double expected_hy = (std::sqrt(3.0) / 2.0) * a0;
    assert(std::abs(tri_domain.get_height() - (ny - 1) * expected_hy) < 1e-6);

    std::cout << "[PASS] Multi-Plane Single-Dislocation Periodic Setup (Square & Triangular, net b = 0) verified." << std::endl;
}

static void test_bicrystal_seam_deadzone_and_triangular_csl() {
    auto check_case = [](const std::string& lattice_type,
                         int nx, int ny,
                         double d_min_a0,
                         double max_void_limit,
                         const std::vector<double>& angles_deg) {
        for (double dth : angles_deg) {
            MacroscopicGeometryConfig cfg;
            cfg.enabled = true;
            cfg.microstructure.mode = MicrostructureMode::Bicrystal;
            cfg.microstructure.gb_cx = 0.5;
            cfg.microstructure.gb_cy = 0.5;
            cfg.microstructure.gb_angle_deg = 0.0;
            cfg.microstructure.theta1_deg = -0.5 * dth;
            cfg.microstructure.theta2_deg = +0.5 * dth;
            cfg.microstructure.min_seam_distance_a0 = d_min_a0;
            cfg.microstructure.seam_resolution = SeamOverlapResolution::MergeMidpoint;
            cfg.microstructure.gb_normal_gap_a0 = 0.0;
            cfg.microstructure.gb_tangential_shift_a0 = 0.0;

            auto [pts, domain] = cfg.generateMicrostructureLattice(lattice_type, nx, ny, 1.0, 0.0, 3);
            assert(!pts.empty());

            // Verify minimum pairwise separation near and away from seam
            const size_t n = pts.size();
            const double seam_x = 0.5 * domain.size_x;
            double min_dist = 1e9;
            std::vector<Eigen::Vector2d> seam_band;
            for (size_t i = 0; i < n; ++i) {
                if (std::abs(pts[i].coord.x() - seam_x) < 3.5) {
                    seam_band.push_back(pts[i].coord);
                }
            }
            for (size_t i = 0; i < seam_band.size(); ++i) {
                for (size_t j = i + 1; j < seam_band.size(); ++j) {
                    if (std::abs(seam_band[i].y() - seam_band[j].y()) > 1.5) continue;
                    const double d = (seam_band[i] - seam_band[j]).norm();
                    if (d < min_dist) min_dist = d;
                }
            }
            assert(min_dist >= d_min_a0 - 1e-3);

            // Verify zero missing-atom monovacancy voids along the interior vertical GB line
            double max_void_r = 0.0;
            const double y_lo = 4.0;
            const double y_hi = domain.size_y - 4.0;
            for (double y = y_lo; y <= y_hi; y += 0.1) {
                double best_r = 1e9;
                const Eigen::Vector2d q(seam_x, y);
                for (const auto& sp : seam_band) {
                    if (std::abs(sp.y() - y) > 2.5) continue;
                    const double r = (sp - q).norm();
                    if (r < best_r) best_r = r;
                }
                if (best_r > max_void_r) max_void_r = best_r;
            }
            if (max_void_r >= max_void_limit) {
                std::cerr << "[FAIL] " << lattice_type << " dTheta=" << dth
                          << " max_void_r=" << max_void_r << " >= " << max_void_limit << std::endl;
            }
            assert(max_void_r < max_void_limit);
        }
    };

    check_case("triangular", 60, 80, 0.52, 0.82, {2.0, 4.0, 6.0, 8.0, 10.0, 21.7868, 32.2042});
    check_case("square", 40, 40, 0.48, 0.88, {6.0, 36.87});
    std::cout << "[PASS] Bicrystal seam dead-zone rescue & triangular CSL snap tests passed." << std::endl;
}

static void test_microstructure_a0_scale_invariance() {
    auto make_ref_dndx = [](const std::string& lattice_type, double a0) -> Eigen::Matrix<double, 3, 2> {
        Eigen::Matrix<double, 3, 2> dNdxi;
        dNdxi << -1.0, -1.0,
                  1.0,  0.0,
                  0.0,  1.0;
        Eigen::Matrix2d J;
        if (lattice_type == "triangular") {
            J << a0, 0.5 * a0,
                 0.0, 0.5 * std::sqrt(3.0) * a0;
        } else {
            J << a0, 0.0,
                 0.0, a0;
        }
        return dNdxi * J.inverse();
    };

    struct CaseSpec {
        std::string name;
        std::string lattice_type;
        MicrostructureMode mode;
        double dtheta_deg;
        int nx;
        int ny;
        int potential_mode;
    };

    const std::vector<CaseSpec> cases = {
        {"tri_bicrystal_21.79", "triangular", MicrostructureMode::Bicrystal, 21.786789, 24, 24, 1},
        {"tri_bicrystal_32.20", "triangular", MicrostructureMode::Bicrystal, 32.204228, 24, 24, 1},
        {"sq_bicrystal_36.87", "square", MicrostructureMode::Bicrystal, 36.869898, 24, 24, 1},
        {"tri_voronoi_4g", "triangular", MicrostructureMode::VoronoiPolycrystal, 0.0, 24, 24, 1},
    };
    const std::array<double, 4> test_a0s = {0.6872044442, 0.99640715, 1.07456993, 1.20933111};

    for (const auto& cs : cases) {
        auto run_at_a0 = [&](double a0, MacroscopicGeometryConfig& cfg, std::vector<Point2D>& pts) -> size_t {
            cfg = MacroscopicGeometryConfig();
            cfg.enabled = true;
            cfg.microstructure.mode = cs.mode;
            if (cs.mode == MicrostructureMode::Bicrystal) {
                cfg.microstructure.gb_cx = 0.5;
                cfg.microstructure.gb_cy = 0.5;
                cfg.microstructure.gb_angle_deg = 0.0;
                cfg.microstructure.theta1_deg = -0.5 * cs.dtheta_deg;
                cfg.microstructure.theta2_deg = +0.5 * cs.dtheta_deg;
                cfg.microstructure.min_seam_distance_a0 = (cs.lattice_type == "triangular") ? 0.52 : 0.50;
                cfg.microstructure.seam_resolution = SeamOverlapResolution::MergeMidpoint;
                cfg.microstructure.gb_normal_gap_a0 = 0.0;
                cfg.microstructure.max_gb_edge_a0 = (cs.lattice_type == "triangular") ? 1.85 : 1.95;
            } else {
                cfg.microstructure.num_grains = 4;
                cfg.microstructure.random_seed = 42;
                cfg.microstructure.min_seam_distance_a0 = 0.52;
                cfg.microstructure.seam_resolution = SeamOverlapResolution::MergeMidpoint;
                cfg.microstructure.max_gb_edge_a0 = 1.85;
            }

            auto [gen_pts, domain] = cfg.generateMicrostructureLattice(
                cs.lattice_type, cs.nx, cs.ny, a0, 0.0, cs.potential_mode);
            pts = std::move(gen_pts);
            assert(!pts.empty());
            assert(std::abs(cfg.a0 - a0) < 1e-12);
            cfg.setPhysicalBoundsFromPoints(pts, a0);

            const std::array<double, 2> offsets = (cs.lattice_type == "square")
                ? std::array<double, 2>{a0, a0}
                : std::array<double, 2>{0.5 * a0, 0.5 * std::sqrt(3.0) * a0};
            auto [orig_map, trans_map] = MeshGenerator::create_domain_maps(pts.size(), domain, offsets);
            auto [interior_map, full_map] = create_dof_mapping_original(pts, 0.25 * a0, false);

            Point2D dims_pt(domain.size_x, domain.size_y);
            AdaptiveMesher mesher(dims_pt, offsets, orig_map, trans_map, full_map, 1e-6 * a0 * a0, false);
            mesher.setUsePeriodicCopies(false);
            AdaptiveMesher::setGlobalGeometryFilter(&cfg, &pts);

            alglib::real_1d_array x;
            const int n_vars = static_cast<int>(interior_map.size());
            x.setlength(2 * n_vars);
            map_points_to_solver_array(x, pts, interior_map, n_vars);

            const auto dndx = make_ref_dndx(cs.lattice_type, a0);
            auto [elems, active] = mesher.createMesh(pts, x, Eigen::Matrix2d::Identity(), &dndx);
            AdaptiveMesher::setGlobalGeometryFilter(nullptr, nullptr);
            assert(!elems.empty());
            assert(active.size() == elems.size());
            return elems.size();
        };

        MacroscopicGeometryConfig ref_cfg;
        std::vector<Point2D> ref_pts;
        const size_t ref_elem_count = run_at_a0(1.0, ref_cfg, ref_pts);

        for (double a0 : test_a0s) {
            MacroscopicGeometryConfig cfg;
            std::vector<Point2D> pts;
            const size_t elem_count = run_at_a0(a0, cfg, pts);
            assert(pts.size() == ref_pts.size());
            for (size_t k = 0; k < pts.size(); ++k) {
                assert(pts[k].grain_id == ref_pts[k].grain_id);
                assert(((pts[k].coord / a0) - ref_pts[k].coord).norm() < 1e-12);
            }
            assert(elem_count == ref_elem_count);
        }
    }
    std::cout << "[PASS] Microstructure a0 scale-invariance (coords + initial AdaptiveMesher) tests passed." << std::endl;
}

static void test_square_calculator_dynamic_cutoff() {
    const double a0 = 0.6872044442;
    const double rcut = 2.5;
    SquareLatticeCalculator calc(a0, rcut);

    for (double gamma : {3.0, 4.0}) {
        Eigen::Matrix2d F;
        F << 1.0, gamma,
             0.0, 1.0;
        const Eigen::Matrix2d C = F.transpose() * F;

        const double e_calc = calc.calculate_energy(C, lennard_jones_energy_v2, 0.0);
        const Eigen::Matrix2d dE_calc = calc.calculate_derivative(C, lennard_jones_energy_der_v2);

        double e_ref = 0.0;
        Eigen::Matrix2d dE_ref = Eigen::Matrix2d::Zero();
        double e_trunc10 = 0.0;
        for (int m = -40; m <= 40; ++m) {
            for (int n = -40; n <= 40; ++n) {
                if (m == 0 && n == 0) continue;
                const Eigen::Vector2d p(static_cast<double>(m), static_cast<double>(n));
                const double r = a0 * std::sqrt(p.transpose() * C * p);
                if (r <= rcut) {
                    const double phi_term = 0.5 * lennard_jones_energy_v2(r);
                    e_ref += phi_term;
                    if (std::abs(m) <= 10 && std::abs(n) <= 10) {
                        e_trunc10 += phi_term;
                    }
                    const double dphi = 0.5 * lennard_jones_energy_der_v2(r) / r;
                    dE_ref(0, 0) += 0.5 * dphi * p.x() * p.x() * a0 * a0;
                    dE_ref(1, 1) += 0.5 * dphi * p.y() * p.y() * a0 * a0;
                    dE_ref(0, 1) += 0.5 * dphi * p.x() * p.y() * a0 * a0;
                    dE_ref(1, 0) += 0.5 * dphi * p.x() * p.y() * a0 * a0;
                }
            }
        }
        assert(std::abs(e_calc - e_ref) < 1e-12);
        assert((dE_calc - dE_ref).cwiseAbs().maxCoeff() < 1e-12);
        if (gamma >= 3.5) {
            assert(std::abs(e_calc - e_trunc10) > 1e-4);
        }
    }
    std::cout << "[PASS] SquareLatticeCalculator dynamic cutoff tests passed." << std::endl;
}

static void test_gb_tangential_shift_candidates() {
    MacroscopicGeometryConfig cfg;
    cfg.enabled = true;
    cfg.microstructure.mode = MicrostructureMode::Bicrystal;
    cfg.microstructure.auto_scan_gb_shift = true;
    cfg.microstructure.gb_shift_scan_steps = 9;
    cfg.microstructure.gb_angle_deg = 0.0;
    cfg.microstructure.gb_tangential_shift_a0 = 0.0;

    // 1) Perfect crystal -> trivial {0.0}
    cfg.microstructure.theta1_deg = 0.0;
    cfg.microstructure.theta2_deg = 0.0;
    auto c0 = cfg.computeTangentialShiftCandidates("triangular");
    assert(c0.size() == 1 && std::abs(c0[0]) < 1e-12);

    cfg.microstructure.theta1_deg = -30.0;
    cfg.microstructure.theta2_deg = 30.0;
    auto c60 = cfg.computeTangentialShiftCandidates("triangular");
    assert(c60.size() == 1 && std::abs(c60[0]) < 1e-12);

    // 2) Symmetric Σ13 (32.2042 deg) -> half period = 0.5 * sqrt(13) ≈ 1.80278
    cfg.microstructure.theta1_deg = -16.102114;
    cfg.microstructure.theta2_deg = +16.102114;
    auto c13 = cfg.computeTangentialShiftCandidates("triangular");
    assert(c13.size() == 9);
    assert(std::abs(c13.front() - 0.0) < 1e-9);
    assert(std::abs(c13.back() - 0.5 * std::sqrt(13.0)) < 1e-6);

    // 3) Symmetric Σ5 square (36.8699 deg) -> half period = 0.5 * sqrt(10) ≈ 1.58114
    cfg.microstructure.theta1_deg = -18.434949;
    cfg.microstructure.theta2_deg = +18.434949;
    auto c5 = cfg.computeTangentialShiftCandidates("square");
    assert(c5.size() == 9);
    assert(std::abs(c5.back() - 0.5 * std::sqrt(10.0)) < 1e-6);

    // 4) Asymmetric tilt -> exact `steps` count, spans negative and positive around 0.0
    cfg.microstructure.theta1_deg = -8.0;
    cfg.microstructure.theta2_deg = +22.0;
    auto casym = cfg.computeTangentialShiftCandidates("triangular");
    assert(casym.size() == 9);
    assert(std::abs(casym[0]) < 1e-12);
    bool has_neg = false, has_pos = false;
    for (double s : casym) {
        if (s < -1e-4) has_neg = true;
        if (s > +1e-4) has_pos = true;
    }
    assert(has_neg && has_pos);

    cfg.microstructure.gb_shift_scan_steps = 8;
    auto casym8 = cfg.computeTangentialShiftCandidates("triangular");
    assert(casym8.size() == 8);
    assert(std::abs(casym8[0]) < 1e-12);
    cfg.microstructure.gb_shift_scan_steps = 9;

    // 5) Horizontal GB (gb_angle_deg = 90 deg, e.g. bicrystalIndentation preset): symmetry-folded
    cfg.microstructure.gb_angle_deg = 90.0;
    cfg.microstructure.theta1_deg = -12.0;
    cfg.microstructure.theta2_deg = +12.0;
    auto ch_sq = cfg.computeTangentialShiftCandidates("square");
    assert(ch_sq.size() == 9);
    assert(std::abs(ch_sq.front() - 0.0) < 1e-9);
    assert(ch_sq.back() > 0.5);
    for (double s : ch_sq) assert(s >= -1e-12); // symmetric -> non-negative [0, max_shift]
    auto ch_tri = cfg.computeTangentialShiftCandidates("triangular");
    assert(ch_tri.size() == 9);
    for (double s : ch_tri) assert(s >= -1e-12);

    // 6) Symmetry-folded angle outside [-30, 30] deg: triangular Σ13 at ±43.8979 deg
    cfg.microstructure.gb_angle_deg = 0.0;
    cfg.microstructure.theta1_deg = -43.897886;
    cfg.microstructure.theta2_deg = +43.897886;
    auto c13_fold = cfg.computeTangentialShiftCandidates("triangular");
    assert(c13_fold.size() == 9);
    assert(std::abs(c13_fold.back() - 0.5 * std::sqrt(13.0)) < 1e-6);

    // 7) Smooth clamp for L_CSL > 5.0 (no cliff drop back to 1.0):
    //    triangular Σ13 (3,1) at ±13.8979 deg -> L_CSL = sqrt(39) ≈ 6.245 -> half = 3.122 -> clamped to 2.5
    cfg.microstructure.theta1_deg = -13.897886;
    cfg.microstructure.theta2_deg = +13.897886;
    auto c13_31 = cfg.computeTangentialShiftCandidates("triangular");
    assert(c13_31.size() == 9);
    assert(std::abs(c13_31.back() - 2.5) < 1e-6);

    //    square Σ13 (5,1) at ±11.3099 deg -> L_CSL = sqrt(26) ≈ 5.099 -> half = 2.550 -> clamped to 2.5
    cfg.microstructure.theta1_deg = -11.309932;
    cfg.microstructure.theta2_deg = +11.309932;
    auto c13_sq = cfg.computeTangentialShiftCandidates("square");
    assert(c13_sq.size() == 9);
    assert(std::abs(c13_sq.back() - 2.5) < 1e-6);

    // 8) Square (1,0) asymmetric bicrystal (0 deg + 26.565 deg -> D1=1, D2=sqrt(5), Δ = sqrt(5)-2 ≈ 0.236 -> clamped 0.5)
    cfg.microstructure.theta1_deg = 0.0;
    cfg.microstructure.theta2_deg = 26.565051;
    auto csq_asym = cfg.computeTangentialShiftCandidates("square");
    assert(csq_asym.size() == 9);
    double max_abs_sq = 0.0;
    for (double s : csq_asym) max_abs_sq = std::max(max_abs_sq, std::abs(s));
    assert(max_abs_sq <= 0.26); // span clamped to 0.5 -> |s| <= 0.25

    // 9) Disabled auto_scan_gb_shift and non-Bicrystal mode return single shift
    cfg.microstructure.gb_tangential_shift_a0 = 0.33;
    cfg.microstructure.auto_scan_gb_shift = false;
    auto c_off = cfg.computeTangentialShiftCandidates("triangular");
    assert(c_off.size() == 1 && std::abs(c_off[0] - 0.33) < 1e-12);
    cfg.microstructure.auto_scan_gb_shift = true;
    cfg.microstructure.mode = MicrostructureMode::VoronoiPolycrystal;
    auto c_vor = cfg.computeTangentialShiftCandidates("triangular");
    assert(c_vor.size() == 1 && std::abs(c_vor[0] - 0.33) < 1e-12);

    // 10) JSON round-trip and clamping for gb_tangential_shift_a0 / auto_scan_gb_shift / gb_shift_scan_steps
    {
        const std::string tmp_json = "/tmp/test_gb_shift_scan_cfg.json";
        {
            std::ofstream out(tmp_json);
            out << R"({
              "enabled": true,
              "microstructure": {
                "mode": "bicrystal",
                "gb_tangential_shift_a0": 1.8028,
                "auto_scan_gb_shift": true,
                "gb_shift_scan_steps": 99
              }
            })";
        }
        MacroscopicGeometryConfig loaded;
        assert(MacroscopicGeometryConfig::loadFromJSON(tmp_json, loaded));
        assert(loaded.microstructure.auto_scan_gb_shift == true);
        assert(loaded.microstructure.gb_shift_scan_steps == 25); // clamped to 25
        assert(std::abs(loaded.microstructure.gb_tangential_shift_a0 - 1.8028) < 1e-6);
        std::remove(tmp_json.c_str());
    }

    // 11) Checkpoint round-trip of gb_tangential_shift_a0
    {
        const std::string tmp_chk = "/tmp/test_gb_shift_roundtrip.chk";
        ConfigurationSaver::CheckpointData chk_out;
        chk_out.nx = 4; chk_out.ny = 4;
        chk_out.gb_tangential_shift_a0 = 0.90138782;
        chk_out.has_gb_tangential_shift = true;
        Point2D p(0.0, 0.0); p.grain_id = 1; p.grain_theta_deg = 16.1021;
        chk_out.points.push_back(p);
        ConfigurationSaver::saveCheckpoint(tmp_chk, chk_out);

        ConfigurationSaver::CheckpointData chk_in;
        Eigen::Matrix<double, 3, 2> dndx = Eigen::Matrix<double, 3, 2>::Zero();
        assert(ConfigurationSaver::loadCheckpoint(tmp_chk, chk_in, dndx));
        assert(chk_in.has_gb_tangential_shift == true);
        assert(std::abs(chk_in.gb_tangential_shift_a0 - 0.90138782) < 1e-12);
        assert(chk_in.points.size() == 1 && chk_in.points[0].grain_id == 1);
        std::remove(tmp_chk.c_str());
    }

    std::cout << "[PASS] GB tangential shift candidate & checkpoint tests passed." << std::endl;
}

void run_triangular_reference_tests();

int main() {
    // Test square lattice generation
    int nx = 2;
    int ny = 2;
    double lattice_constant = 1.0;
    std::string lattice_type = "square";
    
    std::vector<Point2D> points = LatticeGenerator::generate_2d_lattice(
        nx, ny, lattice_constant, lattice_type);
    
    // Check number of points
    assert(points.size() == nx * ny);
    
    // Check a few point coordinates
    assert(points[0].coord.x() == 0.0);
    assert(points[0].coord.y() == 0.0);
    
    assert(points[1].coord.x() == 1.0);
    assert(points[1].coord.y() == 0.0);
    
    assert(points[2].coord.x() == 0.0);
    assert(points[2].coord.y() == 1.0);
    
    assert(points[3].coord.x() == 1.0);
    assert(points[3].coord.y() == 1.0);
    
    test_polygon_ray_casting_and_cutouts();
    test_nanoindenter_shapes_and_penalty();
    test_periodic_dislocations_setup();
    test_bicrystal_seam_deadzone_and_triangular_csl();
    test_microstructure_a0_scale_invariance();
    test_square_calculator_dynamic_cutoff();
    test_gb_tangential_shift_candidates();
    run_triangular_reference_tests();

    std::cout << "All tests passed!" << std::endl;
    return 0;
}
