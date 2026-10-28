#pragma once
#ifndef MACROSCOPIC_GEOMETRY_H
#define MACROSCOPIC_GEOMETRY_H

#include <string>
#include <vector>
#include <utility>
#include <Eigen/Dense>
#include "Point2D.h"
#include "DomainDimensions.h"

struct UserData;

enum class OuterDomainShapeType {
    Rectangle = 0,
    Dogbone = 1,
    Disk = 2,
    Polygon = 3
};

enum class CutoutShapeType {
    Circle = 0,
    Ellipse = 1,
    Rectangle = 2,
    EdgeNotch = 3,
    PreCrack = 4,
    Polygon = 5
};

struct GeometryCutout {
    std::string id = "cutout_1";
    CutoutShapeType type = CutoutShapeType::Circle;
    bool enabled = true;

    // Normalized center (cx, cy) in [0, 1] x [0, 1]
    double cx = 0.5;
    double cy = 0.5;

    // Normalized radii / half-dimensions
    // For Circle: radius r is scaled by min(Lx, Ly)
    // For Ellipse/Rectangle: rx is scaled by Lx, ry is scaled by Ly (width = 2*rx, height = 2*ry)
    double rx = 0.15;
    double ry = 0.15;
    double rotation_deg = 0.0;

    // EdgeNotch specific: edge = "left", "right", "top", "bottom"
    // cx or cy acts as position along edge, rx = depth, ry = half-width
    std::string notch_edge = "left";
    std::string notch_profile = "v_notch"; // "v_notch" or "u_notch"

    // PreCrack segment endpoints in normalized [0, 1]^2 coordinates & half-thickness
    double x1 = 0.0;
    double y1 = 0.5;
    double x2 = 0.35;
    double y2 = 0.5;
    double thickness = 0.015; // Normalized thickness

    // Polygon cutout vertices in normalized [0, 1]^2 coordinates
    std::vector<std::pair<double, double>> vertices;
};

enum class BoundaryRegionType {
    TopEdge = 0,
    BottomEdge = 1,
    LeftEdge = 2,
    RightEdge = 3,
    OuterPerimeter = 4,
    HoleSurfaces = 5,
    CustomBox = 6
};

enum class IndenterShapeType {
    Cylinder = 0,   // Spherical / Cylindrical 2D indenter
    FlatPunch = 1,  // Flat punch indenter
    Wedge = 2       // 2D Wedge / Berkovich cross-section
};

enum class IndenterContactMode {
    KinematicProjection = 0, // Option A: Hard projection onto indenter profile + Dirichlet constraint
    RepulsivePenalty = 1     // Option B: Smooth harmonic penalty potential
};

struct IndenterConfig {
    bool enabled = false;
    IndenterShapeType shape = IndenterShapeType::Cylinder;
    IndenterContactMode contact_mode = IndenterContactMode::KinematicProjection;

    // Normalized horizontal tip center in [0, 1]
    double cx = 0.50;

    // Radius in units of lattice spacing a_0 (e.g. 25.0)
    double radius_a0 = 25.0;

    // Maximum indentation stroke delta_max in units of a_0 (e.g. 5.0)
    double max_depth_a0 = 5.0;

    // Sub-surface cushion flag to prevent inverted elements (Option A)
    bool cushion_subsurface = true;
    double cushion_offset_a0 = 0.20;

    // Penalty stiffness k_c for Option B (energy units / a0^2)
    double penalty_stiffness = 100.0;

    // Wedge half-angle in degrees (for IndenterShapeType::Wedge)
    double wedge_half_angle_deg = 70.3;

    // Cached physical parameters computed during bounds setup
    double physical_radius = 25.0;
    double physical_max_depth = 5.0;
    double physical_cx = 50.0;
    double physical_y_top = 100.0;
    double physical_cushion_offset = 0.20;
};

enum class BoundaryConditionKind {
    FreeSurface = 0,            // Unconstrained DOF (natural traction-free boundary)
    FixedClamped = 1,           // Clamped u = 0
    PrescribedDisplacement = 2, // Driven u(alpha) = alpha * (ux * L_ref, uy * L_ref)
    AffineDeformation = 3       // Driven x(alpha) = F(alpha) * X_0
};

struct BoundaryZoneRule {
    std::string id = "bc_1";
    BoundaryRegionType region = BoundaryRegionType::TopEdge;
    BoundaryConditionKind kind = BoundaryConditionKind::FreeSurface;
    bool enabled = true;

    // Grip thickness in units of equilibrium lattice spacing a_0
    double thickness_a0 = 2.0;

    // CustomBox bounds in normalized [0, 1]^2 coordinates
    double box_x_min = 0.0;
    double box_y_min = 0.0;
    double box_x_max = 1.0;
    double box_y_max = 0.1;

    // Prescribed displacement multipliers per unit alpha (scaled by Ly by default)
    double ux = 0.0;
    double uy = 1.0;

    // Affine mode: "simple_shear_x", "simple_shear_y", "uniaxial_y", "uniaxial_x", "pure_shear", "biaxial"
    std::string affine_mode = "simple_shear_x";
};

struct MacroscopicGeometryConfig {
    bool enabled = false;
    std::string preset_name = "custom";

    // Outer domain shape parameters
    OuterDomainShapeType outer_shape = OuterDomainShapeType::Rectangle;
    double dogbone_gauge_width = 0.50;   // Ratio of central gauge width to Lx
    double dogbone_gauge_length = 0.44;  // Ratio of central straight gauge length to Ly
    double dogbone_fillet_radius = 0.12; // Normalized fillet transition radius
    double disk_radius = 0.48;           // Normalized radius in units of min(Lx, Ly)
    std::vector<std::pair<double, double>> outer_polygon;

    // Cutouts (holes, slots, notches, cracks)
    std::vector<GeometryCutout> cutouts;

    // Boundary condition rules (evaluated in order; later rules override earlier rules,
    // or Fixed/Driven rules take precedence over FreeSurface)
    std::vector<BoundaryZoneRule> boundary_rules;

    // Nanoindenter Contact Tool configuration
    IndenterConfig indenter;

    // Primary macroscopic stress component to log in Stress-Strain curves:
    // "auto", "sigma_xy", "sigma_yy", "sigma_xx", "von_mises", "hydrostatic", "indentation_load"
    std::string primary_stress_component = "auto";

    // Cached physical bounding box of the uncarved reference crystal
    double x_min = 0.0;
    double y_min = 0.0;
    double Lx = 100.0;
    double Ly = 100.0;
    double a0 = 1.0;

    // Cached physical bounding box of the carved solid crystal
    double carved_x_min = 0.0;
    double carved_x_max = 100.0;
    double carved_y_min = 0.0;
    double carved_y_max = 100.0;

    // Load from JSON file
    static bool loadFromJSON(const std::string& filepath, MacroscopicGeometryConfig& out_config);

    // Initialize bounding box from uncarved reference crystal
    void setPhysicalBoundsFromPoints(const std::vector<Point2D>& raw_points, double lattice_a0);

    // Test if a physical point (x, y) is inside the macroscopic solid (inside outer shape & outside cutouts)
    bool isPointInsideSolid(double x, double y) const;

    // Test if a physical point (x, y) is inside any active cutout (expanded by optional margin)
    bool isPointNearCutout(double x, double y, double margin_phys) const;

    // Test if a triangle (p0, p1, p2) in reference/initial space should be rejected
    // (centroid or edge midpoint outside solid, or edge crosses a pre-crack slit)
    bool shouldRejectTriangle(const Eigen::Vector2d& p0,
                              const Eigen::Vector2d& p1,
                              const Eigen::Vector2d& p2) const;

    // Carve lattice points: retains only points inside solid
    std::vector<Point2D> carveLatticePoints(const std::vector<Point2D>& raw_points, double lattice_a0);

    // Classify surviving nodes and build DOF mappings + per-node BC roles and velocity vectors
    // bc_role: 0 = Free (interior or traction-free surface), 1 = Fixed (u = 0), 2 = Driven (u(alpha))
    void buildBoundaryDOFs(
        const std::vector<Point2D>& carved_points,
        std::vector<std::pair<int, int>>& interior_mapping,
        std::vector<std::pair<int, int>>& full_mapping,
        std::vector<int>& node_bc_role,
        std::vector<Eigen::Vector2d>& node_bc_velocity) const;

    // Apply load step alpha to points (updates Fixed and Driven nodes exactly,
    // and applies a smooth predictor increment to Free nodes)
    void applyLoadStep(
        std::vector<Point2D>& current_points,
        const std::vector<Point2D>& points_ref,
        const std::vector<int>& node_bc_role,
        const std::vector<Eigen::Vector2d>& node_bc_velocity,
        double alpha,
        double dalpha,
        bool is_first_step,
        unsigned int seed = 42,
        double noise_level = 0.0) const;

    // Compute approximate macroscopic F_ext(alpha) for DefectAnalysis & reference frame
    Eigen::Matrix2d computeMacroscopicF(double alpha) const;

    // Extract the conjugate scalar stress component for logging and avalanche detection
    double extractPrimaryStress(const Eigen::Matrix2d& stress_tensor,
                                const std::vector<Point2D>* current_points = nullptr,
                                double alpha = 0.0,
                                const UserData* userData = nullptr,
                                const std::vector<int>* contact_atoms = nullptr) const;

    struct IndenterPointContact {
        bool in_contact = false;
        double penetration = 0.0;
        Eigen::Vector2d surface_point = Eigen::Vector2d::Zero();
        Eigen::Vector2d normal = Eigen::Vector2d(0.0, -1.0); // Outward from indenter into substrate
        Eigen::Matrix2d hessian = Eigen::Matrix2d::Zero();   // 2x2 nodal Hessian of penalty energy
    };

    // Nanoindenter contact helpers
    bool isIndenterActive() const { return enabled && indenter.enabled; }
    void setupIndenterBounds(double lattice_a0, const std::vector<Point2D>* carved_points = nullptr);
    Eigen::Vector2d getIndenterPosition(double alpha) const;
    IndenterPointContact evaluateIndenterContact(const Eigen::Vector2d& p,
                                                 double alpha,
                                                 double extra_margin = 0.0) const;
    double computeIndenterPenaltyEnergy(const std::vector<Point2D>& points, double alpha) const;
    double accumulateIndenterPenalty(const alglib::real_1d_array& x,
                                     const std::vector<std::pair<int, int>>& interior_mapping,
                                     alglib::real_1d_array& grad,
                                     double alpha) const;
    double computeIndenterVerticalForce(const std::vector<Point2D>& points, double alpha) const;
    void updateIndenterContactKinematics(std::vector<Point2D>& current_points,
                                        const std::vector<Point2D>& points_ref,
                                        const std::vector<int>& baseline_bc_role,
                                        std::vector<int>& node_bc_role,
                                        std::vector<int>& out_contact_atoms,
                                        double alpha) const;
    void rebuildBoundaryDOFs(const std::vector<int>& node_bc_role,
                             std::vector<std::pair<int, int>>& interior_mapping,
                             std::vector<std::pair<int, int>>& full_mapping) const;
    double computeKinematicContactForce(const UserData* userData,
                                        const std::vector<int>& contact_atoms) const;
};

#endif // MACROSCOPIC_GEOMETRY_H
