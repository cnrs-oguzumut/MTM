#include <iostream>
#include <cassert>
#include <cmath>
#include "../include/optimization/LatticeOptimizer.h"
#include "../include/geometry/LatticeGenerator.h"
#include "../include/geometry/DomainDimensions.h"
#include "../include/geometry/MacroscopicGeometry.h"

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

    std::cout << "All tests passed!" << std::endl;
    return 0;
}
