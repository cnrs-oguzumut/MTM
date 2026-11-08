#include "../include/geometry/MacroscopicGeometry.h"
#include "../include/lattice_energy/ReferenceFrame.h"
#include "../include/geometry/DomainInfo.h"
#include "../include/optimization/LatticeOptimizer.h"
#include "../include/mesh/ElementTriangle2D.h"
#include "../include/utils/dislocation_utils.h"
#include <unordered_set>
#include <unordered_map>
#include <limits>

#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>
#include <cctype>

namespace {

// Centre-symmetric affine modes (biaxial, radial, pure shear, torsion) pivot about the specimen
// centre only when applied on a closed contour (outer rim or hole rims). Edge and custom-box
// rules keep the bottom-left corner origin, so their strain matches a one-sided grip.
bool affineRuleUsesCenterOrigin(const BoundaryZoneRule& rule) {
    if (rule.kind != BoundaryConditionKind::AffineDeformation) return false;
    if (rule.region != BoundaryRegionType::OuterPerimeter &&
        rule.region != BoundaryRegionType::HoleSurfaces) return false;
    const std::string& m = rule.affine_mode;
    return m == "biaxial" || m == "radial_outward" || m == "radial_inward" ||
           m == "biaxial_compression" || m == "pure_shear" ||
           m == "circumferential_shear" || m == "torsion";
}

// Minimal, dependency-free JSON token/object parser helpers for geometry_config.json
std::string trim(const std::string& s) {
    size_t a = 0;
    while (a < s.size() && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    size_t b = s.size();
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
    return s.substr(a, b - a);
}

std::string extractRawValue(const std::string& json, const std::string& key) {
    std::string pattern = "\"" + key + "\"";
    size_t pos = json.find(pattern);
    if (pos == std::string::npos) return "";
    size_t colon = json.find(':', pos + pattern.size());
    if (colon == std::string::npos) return "";
    size_t i = colon + 1;
    while (i < json.size() && std::isspace(static_cast<unsigned char>(json[i]))) ++i;
    if (i >= json.size()) return "";

    if (json[i] == '"') {
        size_t end = i + 1;
        while (end < json.size()) {
            if (json[end] == '"' && json[end - 1] != '\\') break;
            ++end;
        }
        return json.substr(i + 1, end - (i + 1));
    } else if (json[i] == '[' || json[i] == '{') {
        char open_c = json[i];
        char close_c = (open_c == '[') ? ']' : '}';
        int depth = 0;
        size_t end = i;
        bool in_str = false;
        while (end < json.size()) {
            char c = json[end];
            if (c == '"' && (end == 0 || json[end - 1] != '\\')) in_str = !in_str;
            if (!in_str) {
                if (c == open_c) ++depth;
                else if (c == close_c) {
                    --depth;
                    if (depth == 0) {
                        ++end;
                        break;
                    }
                }
            }
            ++end;
        }
        return json.substr(i, end - i);
    } else {
        size_t end = i;
        while (end < json.size() && json[end] != ',' && json[end] != '}' && json[end] != ']' &&
               !std::isspace(static_cast<unsigned char>(json[end]))) {
            ++end;
        }
        return json.substr(i, end - i);
    }
}

bool extractBool(const std::string& json, const std::string& key, bool default_val) {
    std::string raw = extractRawValue(json, key);
    if (raw.empty()) return default_val;
    if (raw == "true" || raw == "1") return true;
    if (raw == "false" || raw == "0") return false;
    return default_val;
}

double extractDouble(const std::string& json, const std::string& key, double default_val) {
    std::string raw = extractRawValue(json, key);
    if (raw.empty()) return default_val;
    try {
        return std::stod(raw);
    } catch (...) {
        return default_val;
    }
}

std::string extractString(const std::string& json, const std::string& key, const std::string& default_val) {
    std::string raw = extractRawValue(json, key);
    if (raw.empty()) return default_val;
    return raw;
}

// Split a JSON array of objects "[ {...}, {...} ]" into individual "{...}" strings
std::vector<std::string> splitJSONObjectArray(const std::string& array_str) {
    std::vector<std::string> items;
    int depth = 0;
    size_t start = std::string::npos;
    bool in_str = false;
    for (size_t i = 0; i < array_str.size(); ++i) {
        char c = array_str[i];
        if (c == '"' && (i == 0 || array_str[i - 1] != '\\')) in_str = !in_str;
        if (in_str) continue;
        if (c == '{') {
            if (depth == 0) start = i;
            ++depth;
        } else if (c == '}') {
            --depth;
            if (depth == 0 && start != std::string::npos) {
                items.push_back(array_str.substr(start, i - start + 1));
                start = std::string::npos;
            }
        }
    }
    return items;
}

// Parse polygon vertices array: either [[x, y], ...] or [{"x": ..., "y": ...}, ...]
std::vector<std::pair<double, double>> parseVerticesArray(const std::string& arr_str) {
    std::vector<std::pair<double, double>> pts;
    if (arr_str.empty()) return pts;

    auto obj_items = splitJSONObjectArray(arr_str);
    if (!obj_items.empty()) {
        for (const auto& item : obj_items) {
            double vx = extractDouble(item, "x", 0.0);
            double vy = extractDouble(item, "y", 0.0);
            pts.emplace_back(vx, vy);
        }
        return pts;
    }

    // Fallback: parse pairs of numbers inside nested brackets [[x,y], [x,y]]
    std::vector<double> nums;
    std::string cur;
    for (char c : arr_str) {
        if (std::isdigit(static_cast<unsigned char>(c)) || c == '.' || c == '-' || c == '+' || c == 'e' || c == 'E') {
            cur.push_back(c);
        } else if (!cur.empty()) {
            try { nums.push_back(std::stod(cur)); } catch (...) {}
            cur.clear();
        }
    }
    if (!cur.empty()) {
        try { nums.push_back(std::stod(cur)); } catch (...) {}
    }
    for (size_t i = 0; i + 1 < nums.size(); i += 2) {
        pts.emplace_back(nums[i], nums[i + 1]);
    }
    return pts;
}

// Point-in-polygon test (ray casting)
bool isPointInPolygon(double x, double y, const std::vector<std::pair<double, double>>& poly) {
    if (poly.size() < 3) return false;
    bool inside = false;
    for (size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++) {
        double xi = poly[i].first, yi = poly[i].second;
        double xj = poly[j].first, yj = poly[j].second;
        if ((yi > y) != (yj > y)) {
            double x_intersect = (xj - xi) * (y - yi) / (yj - yi) + xi;
            if (x < x_intersect) inside = !inside;
        }
    }
    return inside;
}

// Distance from point P to segment AB
double pointToSegmentDistance(const Eigen::Vector2d& p,
                              const Eigen::Vector2d& a,
                              const Eigen::Vector2d& b) {
    Eigen::Vector2d ab = b - a;
    double len2 = ab.squaredNorm();
    if (len2 < 1e-14) return (p - a).norm();
    double t = std::clamp((p - a).dot(ab) / len2, 0.0, 1.0);
    Eigen::Vector2d proj = a + t * ab;
    return (p - proj).norm();
}

// Proper 2D segment-segment intersection test (p1-p2 vs q1-q2)
bool segmentsIntersect(const Eigen::Vector2d& p1, const Eigen::Vector2d& p2,
                       const Eigen::Vector2d& q1, const Eigen::Vector2d& q2) {
    auto ccw = [](const Eigen::Vector2d& A, const Eigen::Vector2d& B, const Eigen::Vector2d& C) {
        return (C.y() - A.y()) * (B.x() - A.x()) > (B.y() - A.y()) * (C.x() - A.x());
    };
    return (ccw(p1, q1, q2) != ccw(p2, q1, q2)) && (ccw(p1, p2, q1) != ccw(p1, p2, q2));
}

// Minimum physical distance from physical point (x, y) to any edge of a normalized [0, 1]^2 polygon
double distanceToPolygonEdges(double x, double y,
                              double x_min, double y_min, double Lx, double Ly,
                              const std::vector<std::pair<double, double>>& poly) {
    if (poly.size() < 2) return 1e30;
    Eigen::Vector2d p(x, y);
    double min_dist = 1e30;
    for (size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++) {
        Eigen::Vector2d a(x_min + poly[j].first * Lx, y_min + poly[j].second * Ly);
        Eigen::Vector2d b(x_min + poly[i].first * Lx, y_min + poly[i].second * Ly);
        min_dist = std::min(min_dist, pointToSegmentDistance(p, a, b));
    }
    return min_dist;
}

// Check if physical segment p0-p1 properly intersects any edge of a normalized [0, 1]^2 polygon
bool segmentIntersectsPolygonEdges(const Eigen::Vector2d& p0, const Eigen::Vector2d& p1,
                                   double x_min, double y_min, double Lx, double Ly,
                                   const std::vector<std::pair<double, double>>& poly) {
    if (poly.size() < 3) return false;
    // Shrink segment slightly so vertices lying right on the polygon boundary are not flagged
    Eigen::Vector2d dir = p1 - p0;
    Eigen::Vector2d s0 = p0 + 1e-4 * dir;
    Eigen::Vector2d s1 = p1 - 1e-4 * dir;
    for (size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++) {
        Eigen::Vector2d a(x_min + poly[j].first * Lx, y_min + poly[j].second * Ly);
        Eigen::Vector2d b(x_min + poly[i].first * Lx, y_min + poly[i].second * Ly);
        if (segmentsIntersect(s0, s1, a, b)) return true;
    }
    return false;
}

} // namespace

bool MacroscopicGeometryConfig::loadFromJSON(const std::string& filepath,
                                             MacroscopicGeometryConfig& out_config) {
    std::ifstream in(filepath);
    if (!in.is_open()) {
        std::cerr << "[MacroscopicGeometry] Error: Could not open JSON file: " << filepath << "\n";
        return false;
    }
    std::stringstream buffer;
    buffer << in.rdbuf();
    std::string json = buffer.str();

    out_config = MacroscopicGeometryConfig();
    out_config.enabled = extractBool(json, "enabled", true);
    out_config.preset_name = extractString(json, "preset_name", "custom");
    out_config.primary_stress_component = extractString(json, "primary_stress_component", "auto");

    std::string outer_str = extractString(json, "outer_shape", "rectangle");
    if (outer_str == "dogbone" || outer_str == "Dogbone") {
        out_config.outer_shape = OuterDomainShapeType::Dogbone;
    } else if (outer_str == "disk" || outer_str == "Disk" || outer_str == "circle") {
        out_config.outer_shape = OuterDomainShapeType::Disk;
    } else if (outer_str == "polygon" || outer_str == "Polygon") {
        out_config.outer_shape = OuterDomainShapeType::Polygon;
    } else {
        out_config.outer_shape = OuterDomainShapeType::Rectangle;
    }

    out_config.dogbone_gauge_width = extractDouble(json, "dogbone_gauge_width", 0.50);
    out_config.dogbone_gauge_length = extractDouble(json, "dogbone_gauge_length", 0.44);
    out_config.dogbone_fillet_radius = extractDouble(json, "dogbone_fillet_radius", 0.12);
    out_config.disk_radius = extractDouble(json, "disk_radius", 0.48);

    std::string outer_poly_raw = extractRawValue(json, "outer_polygon");
    if (!outer_poly_raw.empty()) {
        out_config.outer_polygon = parseVerticesArray(outer_poly_raw);
    }

    // Parse cutouts array
    std::string cutouts_raw = extractRawValue(json, "cutouts");
    for (const auto& item_json : splitJSONObjectArray(cutouts_raw)) {
        GeometryCutout c;
        c.id = extractString(item_json, "id", "cutout");
        c.enabled = extractBool(item_json, "enabled", true);
        std::string type_str = extractString(item_json, "type", "circle");
        if (type_str == "ellipse" || type_str == "Ellipse") c.type = CutoutShapeType::Ellipse;
        else if (type_str == "rectangle" || type_str == "slot" || type_str == "Rectangle") c.type = CutoutShapeType::Rectangle;
        else if (type_str == "edge_notch" || type_str == "edgeNotch" || type_str == "notch") c.type = CutoutShapeType::EdgeNotch;
        else if (type_str == "pre_crack" || type_str == "preCrack" || type_str == "crack") c.type = CutoutShapeType::PreCrack;
        else if (type_str == "polygon" || type_str == "Polygon") c.type = CutoutShapeType::Polygon;
        else c.type = CutoutShapeType::Circle;

        c.cx = extractDouble(item_json, "cx", 0.5);
        c.cy = extractDouble(item_json, "cy", 0.5);
        c.rx = extractDouble(item_json, "rx", 0.15);
        c.ry = extractDouble(item_json, "ry", c.rx);
        c.rotation_deg = extractDouble(item_json, "rotation_deg", 0.0);
        c.notch_edge = extractString(item_json, "notch_edge", "left");
        c.notch_profile = extractString(item_json, "notch_profile", "v_notch");
        c.x1 = extractDouble(item_json, "x1", 0.0);
        c.y1 = extractDouble(item_json, "y1", 0.5);
        c.x2 = extractDouble(item_json, "x2", 0.35);
        c.y2 = extractDouble(item_json, "y2", 0.5);
        c.thickness = extractDouble(item_json, "thickness", 0.015);
        std::string verts_raw = extractRawValue(item_json, "vertices");
        if (!verts_raw.empty()) {
            c.vertices = parseVerticesArray(verts_raw);
        }
        out_config.cutouts.push_back(c);
    }

    auto parseRulesList = [](const std::string& raw_arr) -> std::vector<BoundaryZoneRule> {
        std::vector<BoundaryZoneRule> parsed;
        for (const auto& item_json : splitJSONObjectArray(raw_arr)) {
            BoundaryZoneRule r;
            r.id = extractString(item_json, "id", "bc");
            r.enabled = extractBool(item_json, "enabled", true);
            std::string reg_str = extractString(item_json, "region", "top");
            if (reg_str == "bottom" || reg_str == "bottomEdge" || reg_str == "bottom_edge") r.region = BoundaryRegionType::BottomEdge;
            else if (reg_str == "left" || reg_str == "leftEdge" || reg_str == "left_edge") r.region = BoundaryRegionType::LeftEdge;
            else if (reg_str == "right" || reg_str == "rightEdge" || reg_str == "right_edge") r.region = BoundaryRegionType::RightEdge;
            else if (reg_str == "outer" || reg_str == "outer_perimeter" || reg_str == "outerPerimeter") r.region = BoundaryRegionType::OuterPerimeter;
            else if (reg_str == "hole" || reg_str == "hole_surfaces" || reg_str == "holeSurfaces") r.region = BoundaryRegionType::HoleSurfaces;
            else if (reg_str == "custom_box" || reg_str == "customBox" || reg_str == "box") r.region = BoundaryRegionType::CustomBox;
            else r.region = BoundaryRegionType::TopEdge;

            std::string kind_str = extractString(item_json, "kind", "free");
            if (kind_str == "fixed" || kind_str == "fixedClamped" || kind_str == "fixed_clamped") r.kind = BoundaryConditionKind::FixedClamped;
            else if (kind_str == "prescribed" || kind_str == "prescribedDisplacement" || kind_str == "prescribed_displacement" || kind_str == "driven") r.kind = BoundaryConditionKind::PrescribedDisplacement;
            else if (kind_str == "affine" || kind_str == "affineDeformation" || kind_str == "affine_deformation" || kind_str == "affine_detachment") r.kind = BoundaryConditionKind::AffineDeformation;
            else if (kind_str == "roller_x" || kind_str == "rollerX") r.kind = BoundaryConditionKind::RollerX;
            else if (kind_str == "roller_y" || kind_str == "rollerY") r.kind = BoundaryConditionKind::RollerY;
            else if (kind_str == "traction" || kind_str == "prescribedTraction" || kind_str == "prescribed_traction") r.kind = BoundaryConditionKind::PrescribedTraction;
            else if (kind_str == "force_grip" || kind_str == "forceGrip" || kind_str == "rigid_force_grip" || kind_str == "rigidForceGrip") r.kind = BoundaryConditionKind::RigidForceGrip;
            else r.kind = BoundaryConditionKind::FreeSurface;

            r.thickness_a0 = extractDouble(item_json, "thickness_a0", 2.0);
            r.box_x_min = extractDouble(item_json, "box_x_min", 0.0);
            r.box_y_min = extractDouble(item_json, "box_y_min", 0.0);
            r.box_x_max = extractDouble(item_json, "box_x_max", 1.0);
            r.box_y_max = extractDouble(item_json, "box_y_max", 0.1);
            r.ux = extractDouble(item_json, "ux", 0.0);
            r.uy = extractDouble(item_json, "uy", 0.0);
            r.lock_x = extractBool(item_json, "lock_x", true);
            r.lock_y = extractBool(item_json, "lock_y", true);
            if (r.kind == BoundaryConditionKind::RollerX) {
                r.lock_x = false;
                r.lock_y = true;
            } else if (r.kind == BoundaryConditionKind::RollerY) {
                r.lock_x = true;
                r.lock_y = false;
            }
            r.affine_mode = extractString(item_json, "affine_mode", "simple_shear_x");
            r.traction_mode = extractString(item_json, "traction_mode", "cartesian");
            r.tx = extractDouble(item_json, "tx", 0.0);
            r.ty = extractDouble(item_json, "ty", 0.0);
            parsed.push_back(r);
        }
        return parsed;
    };

    // Parse boundary_rules (Stage 2) and stage1_boundary_rules (Stage 1) arrays
    out_config.boundary_rules = parseRulesList(extractRawValue(json, "boundary_rules"));
    out_config.stage1_boundary_rules = parseRulesList(extractRawValue(json, "stage1_boundary_rules"));

    // Parse indenter object if present
    std::string indenter_raw = extractRawValue(json, "indenter");
    if (!indenter_raw.empty()) {
        out_config.indenter.enabled = extractBool(indenter_raw, "enabled", false);
        std::string shape_str = extractString(indenter_raw, "shape", "cylinder");
        if (shape_str == "flat_punch" || shape_str == "flatPunch") {
            out_config.indenter.shape = IndenterShapeType::FlatPunch;
        } else if (shape_str == "wedge" || shape_str == "Wedge") {
            out_config.indenter.shape = IndenterShapeType::Wedge;
        } else {
            out_config.indenter.shape = IndenterShapeType::Cylinder;
        }

        std::string mode_str = extractString(indenter_raw, "contact_mode", "kinematic");
        if (mode_str == "penalty" || mode_str == "repulsive_penalty" || mode_str == "RepulsivePenalty") {
            out_config.indenter.contact_mode = IndenterContactMode::RepulsivePenalty;
        } else {
            out_config.indenter.contact_mode = IndenterContactMode::KinematicProjection;
        }

        out_config.indenter.cx = extractDouble(indenter_raw, "cx", 0.50);
        out_config.indenter.radius_a0 = extractDouble(indenter_raw, "radius_a0", 25.0);
        out_config.indenter.max_depth_a0 = extractDouble(indenter_raw, "max_depth_a0", 5.0);
        out_config.indenter.cushion_subsurface = extractBool(indenter_raw, "cushion_subsurface", true);
        out_config.indenter.cushion_offset_a0 = extractDouble(indenter_raw, "cushion_offset_a0", 0.20);
        out_config.indenter.penalty_stiffness = extractDouble(indenter_raw, "penalty_stiffness", 100.0);
        out_config.indenter.wedge_half_angle_deg = extractDouble(indenter_raw, "wedge_half_angle_deg", 70.3);
    }

    // Parse microstructure object if present
    std::string micro_raw = extractRawValue(json, "microstructure");
    if (!micro_raw.empty()) {
        std::string mode_str = extractString(micro_raw, "mode", "single_crystal");
        if (mode_str == "bicrystal" || mode_str == "Bicrystal") {
            out_config.microstructure.mode = MicrostructureMode::Bicrystal;
        } else if (mode_str == "voronoi" || mode_str == "voronoi_polycrystal" || mode_str == "VoronoiPolycrystal") {
            out_config.microstructure.mode = MicrostructureMode::VoronoiPolycrystal;
        } else if (mode_str == "periodic_dislocations" || mode_str == "PeriodicDislocations" || mode_str == "dislocations") {
            out_config.microstructure.mode = MicrostructureMode::PeriodicDislocations;
        } else {
            out_config.microstructure.mode = MicrostructureMode::SingleCrystal;
        }

        std::string seam_str = extractString(micro_raw, "seam_resolution", "merge_midpoint");
        if (seam_str == "delete_second" || seam_str == "deleteSecond") {
            out_config.microstructure.seam_resolution = SeamOverlapResolution::DeleteSecond;
        } else {
            out_config.microstructure.seam_resolution = SeamOverlapResolution::MergeMidpoint;
        }

        out_config.microstructure.min_seam_distance_a0 = extractDouble(micro_raw, "min_seam_distance_a0", 0.55);
        out_config.microstructure.max_gb_edge_a0 = std::clamp(extractDouble(micro_raw, "max_gb_edge_a0", 2.60), 1.40, 5.00);
        out_config.microstructure.pre_relax_gb = extractBool(micro_raw, "pre_relax_gb", true);
        out_config.microstructure.relax_only = extractBool(micro_raw, "relax_only", false);
        out_config.microstructure.reanchor_ref_after_relax = extractBool(micro_raw, "reanchor_ref_after_relax", true);
        out_config.microstructure.max_relax_remesh_passes = std::clamp(static_cast<int>(extractDouble(micro_raw, "max_relax_remesh_passes", 5.0)), 1, 25);
        out_config.microstructure.loading_bc_preset = extractString(micro_raw, "loading_bc_preset", "keep_studio");
        out_config.microstructure.stage2_load_type = extractString(micro_raw, "stage2_load_type", "shear_parallel_gb");
        out_config.microstructure.gb_cx = extractDouble(micro_raw, "gb_cx", 0.50);
        out_config.microstructure.gb_cy = extractDouble(micro_raw, "gb_cy", 0.50);
        out_config.microstructure.gb_angle_deg = extractDouble(micro_raw, "gb_angle_deg", 0.0);
        out_config.microstructure.theta1_deg = extractDouble(micro_raw, "theta1_deg", -15.0);
        out_config.microstructure.theta2_deg = extractDouble(micro_raw, "theta2_deg", +15.0);
        out_config.microstructure.gb_normal_gap_a0 = std::clamp(extractDouble(micro_raw, "gb_normal_gap_a0", 0.0), -0.50, 1.50);
        out_config.microstructure.gb_tangential_shift_a0 = std::clamp(extractDouble(micro_raw, "gb_tangential_shift_a0", 0.0), -1.00, 1.00);

        out_config.microstructure.num_grains = static_cast<int>(extractDouble(micro_raw, "num_grains", 6.0));
        out_config.microstructure.random_seed = static_cast<int>(extractDouble(micro_raw, "random_seed", 42.0));
        out_config.microstructure.theta_min_deg = extractDouble(micro_raw, "theta_min_deg", -30.0);
        out_config.microstructure.theta_max_deg = extractDouble(micro_raw, "theta_max_deg", +30.0);

        // Multi-Plane Single-Dislocation Periodic Setup parameters
        out_config.microstructure.dislocation_pair_count = std::clamp(static_cast<int>(extractDouble(micro_raw, "dislocation_pair_count", 2.0)), 1, 16);
        out_config.microstructure.dislocation_layout = extractString(micro_raw, "dislocation_layout", "alternating_random_x");
        out_config.microstructure.min_plane_spacing_rows = std::clamp(static_cast<int>(extractDouble(micro_raw, "min_plane_spacing_rows", 6.0)), 2, 40);
        out_config.microstructure.dislocation_core_radius_a0 = std::clamp(extractDouble(micro_raw, "dislocation_core_radius_a0", 0.80), 0.1, 5.0);
        out_config.microstructure.bound_strip_excess = extractBool(micro_raw, "bound_strip_excess", true);
        out_config.microstructure.enforce_periodic_bc = extractBool(micro_raw, "enforce_periodic_bc", true);

        out_config.microstructure.dislocations.clear();
        std::string dislocs_raw = extractRawValue(micro_raw, "dislocations");
        int disloc_idx = 0;
        for (const auto& d_json : splitJSONObjectArray(dislocs_raw)) {
            SeededDislocation d;
            d.id = static_cast<int>(extractDouble(d_json, "id", disloc_idx));
            d.glide_plane_row = static_cast<int>(extractDouble(d_json, "glide_plane_row", 20.0));
            d.x_frac = std::clamp(extractDouble(d_json, "x_frac", 0.50), 0.01, 0.99);
            d.burgers_sign = (extractDouble(d_json, "burgers_sign", 1.0) < 0) ? -1 : 1;
            d.enabled = extractBool(d_json, "enabled", true);
            out_config.microstructure.dislocations.push_back(d);
            ++disloc_idx;
        }

        out_config.microstructure.seeds.clear();
        std::string seeds_raw = extractRawValue(micro_raw, "seeds");
        int seed_idx = 0;
        for (const auto& s_json : splitJSONObjectArray(seeds_raw)) {
            VoronoiGrainSeed s;
            s.id = static_cast<int>(extractDouble(s_json, "id", seed_idx));
            s.cx = extractDouble(s_json, "cx", 0.5);
            s.cy = extractDouble(s_json, "cy", 0.5);
            s.theta_deg = extractDouble(s_json, "theta_deg", 0.0);
            out_config.microstructure.seeds.push_back(s);
            ++seed_idx;
        }
    }

    std::string micro_label = "SINGLE";
    if (out_config.microstructure.mode == MicrostructureMode::Bicrystal) micro_label = "BICRYSTAL";
    else if (out_config.microstructure.mode == MicrostructureMode::VoronoiPolycrystal) micro_label = "VORONOI";
    else if (out_config.microstructure.mode == MicrostructureMode::PeriodicDislocations) micro_label = "PERIODIC_DISLOCATIONS";

    std::cout << "[MacroscopicGeometry] Loaded '" << out_config.preset_name << "' from " << filepath
              << " (" << out_config.cutouts.size() << " cutouts, "
              << out_config.boundary_rules.size() << " boundary rules, microstructure: " << micro_label
              << ", indenter: "
              << (out_config.indenter.enabled ? (out_config.indenter.contact_mode == IndenterContactMode::RepulsivePenalty ? "PENALTY" : "KINEMATIC") : "OFF")
              << ")\n";
    return true;
}

std::pair<std::vector<Point2D>, DomainDimensions>
MacroscopicGeometryConfig::generateMicrostructureLattice(
    const std::string& lattice_type,
    int nx, int ny,
    double lattice_a0,
    double fallback_theta_deg,
    int potential_mode) {

    if (microstructure.mode == MicrostructureMode::PeriodicDislocations) {
        return generatePeriodicDislocationLattice(lattice_type, nx, ny, lattice_a0, potential_mode);
    }

    const double a = (lattice_a0 > 0.0) ? lattice_a0 : 1.0;
    const double sqrt3 = std::sqrt(3.0);
    const double size_x = nx * a;
    const double size_y = (lattice_type == "triangular") ? (ny * a * sqrt3 * 0.5) : (ny * a);
    const double d_lat = (lattice_type == "triangular") ? (a * sqrt3 * 0.5) : a;
    const double r_box = std::hypot(size_x, size_y);
    const int n_span = static_cast<int>(std::ceil(r_box / d_lat)) + 4;
    const double tol = 1e-5;

    std::vector<Point2D> raw_points;

    auto appendGrainPoints = [&](int gid, double theta_deg, double origin_x, double origin_y,
                                 const std::function<bool(double, double)>& region_pred) {
        const double theta_rad = theta_deg * M_PI / 180.0;
        const double cos_th = std::cos(theta_rad);
        const double sin_th = std::sin(theta_rad);

        const double a1_x = a * cos_th;
        const double a1_y = a * sin_th;
        const double a2_x = (lattice_type == "triangular") ? a * (0.5 * cos_th - 0.5 * sqrt3 * sin_th) : (-a * sin_th);
        const double a2_y = (lattice_type == "triangular") ? a * (0.5 * sin_th + 0.5 * sqrt3 * cos_th) : (a * cos_th);

        for (int i = -n_span; i <= n_span; ++i) {
            for (int j = -n_span; j <= n_span; ++j) {
                double rx = origin_x + i * a1_x + j * a2_x;
                double ry = origin_y + i * a1_y + j * a2_y;
                if (rx >= -tol && rx <= size_x + tol && ry >= -tol && ry <= size_y + tol) {
                    if (region_pred(rx, ry)) {
                        raw_points.emplace_back(rx, ry, gid, theta_deg);
                    }
                }
            }
        }
    };

    if (microstructure.mode == MicrostructureMode::Bicrystal) {
        const double p_gb_x = std::clamp(microstructure.gb_cx, 0.0, 1.0) * size_x;
        const double p_gb_y = std::clamp(microstructure.gb_cy, 0.0, 1.0) * size_y;
        const double phi_rad = microstructure.gb_angle_deg * M_PI / 180.0;
        // Normal vector to the GB plane: phi = 0 deg -> vertical boundary (normal = +x),
        // phi = 90 deg -> horizontal boundary (normal = +y)
        const double nx_gb = std::cos(phi_rad);
        const double ny_gb = std::sin(phi_rad);

        // General (m, n) CSL Rationalizer (gb_code / aimsgb):
        // Snap near-CSL angles (within 0.02 deg of exact atan(n/m) for coprime 1 <= n <= m <= 18)
        // relative to the GB plane so that UI-rounded values produce exact CSL periodicity
        // along the entire boundary without drift between the two halves.
        int matched_csl_m = 0;
        int matched_csl_n = 0;
        auto snapNearCSLAngle = [&](double theta_in_deg) -> double {
            if (lattice_type != "square") return theta_in_deg;
            const double rel_deg = theta_in_deg - microstructure.gb_angle_deg;
            double best_diff = 0.025;
            double best_target_rel = rel_deg;
            for (int m = 1; m <= 18; ++m) {
                for (int n = 1; n <= m; ++n) {
                    if (std::gcd(m, n) != 1) continue;
                    const double exact_abs = std::atan2(static_cast<double>(n), static_cast<double>(m)) * 180.0 / M_PI;
                    for (double sgn : {-1.0, 1.0}) {
                        const double target_rel = sgn * exact_abs;
                        const double diff = std::abs(rel_deg - target_rel);
                        if (diff < best_diff) {
                            best_diff = diff;
                            best_target_rel = target_rel;
                            matched_csl_m = m;
                            matched_csl_n = n;
                        }
                    }
                }
            }
            if (best_diff < 0.025) {
                return microstructure.gb_angle_deg + best_target_rel;
            }
            return theta_in_deg;
        };

        const double th1_exact = snapNearCSLAngle(microstructure.theta1_deg);
        const double th2_exact = snapNearCSLAngle(microstructure.theta2_deg);

        // Include a small coincidence capture band (±0.02 * a) on the GB plane for both grains
        // so coincidence lattice sites on the boundary are never dropped due to floating-point drift,
        // and are then merged symmetrically onto the GB plane by the seam filter.
        const double gb_coincidence_band = 0.02 * a;
        const double tx_gb = -ny_gb;
        const double ty_gb = nx_gb;
        const double delta_par = microstructure.gb_tangential_shift_a0 * a;

        // Grain 0: left/lower half-space (dist <= +gb_coincidence_band)
        appendGrainPoints(0, th1_exact,
                          p_gb_x - 0.5 * delta_par * tx_gb,
                          p_gb_y - 0.5 * delta_par * ty_gb,
                          [&](double x, double y) {
                              double dist = (x - p_gb_x) * nx_gb + (y - p_gb_y) * ny_gb;
                              return dist <= gb_coincidence_band;
                          });

        // Grain 1: right/upper half-space (dist >= -gb_coincidence_band)
        appendGrainPoints(1, th2_exact,
                          p_gb_x + 0.5 * delta_par * tx_gb,
                          p_gb_y + 0.5 * delta_par * ty_gb,
                          [&](double x, double y) {
                              double dist = (x - p_gb_x) * nx_gb + (y - p_gb_y) * ny_gb;
                              return dist >= -gb_coincidence_band;
                          });

        std::cout << "[Microstructure] Generated Bicrystal: GB center=(" << microstructure.gb_cx << ", "
                  << microstructure.gb_cy << "), phi_GB=" << microstructure.gb_angle_deg
                  << "°, theta_1=" << th1_exact << "°, theta_2=" << th2_exact
                  << "°, gap_perp=" << microstructure.gb_normal_gap_a0
                  << " a0, shift_par=" << microstructure.gb_tangential_shift_a0
                  << " a0 (raw points before seam filter: " << raw_points.size() << ")\n";
        if (matched_csl_m > 0 && matched_csl_n > 0) {
            const int sum_sq = matched_csl_m * matched_csl_m + matched_csl_n * matched_csl_n;
            const int sigma = sum_sq / std::gcd(sum_sq, 2);
            const double l_csl = std::sqrt(static_cast<double>(sum_sq));
            const double d_mn = 1.0 / l_csl;
            std::cout << "[Microstructure] CSL Rationalizer: Σ" << sigma << " (" << matched_csl_m
                      << " " << matched_csl_n << " 0), L_CSL = " << l_csl
                      << " a0, d_(mn) = " << d_mn << " a0\n";
        }

    } else if (microstructure.mode == MicrostructureMode::VoronoiPolycrystal) {
        // Ensure physical bounds are initialized before checking solid containment of Voronoi seeds
        x_min = 0.0;
        y_min = 0.0;
        Lx = size_x;
        Ly = size_y;
        a0 = a;

        std::vector<VoronoiGrainSeed> seeds = microstructure.seeds;
        if (seeds.empty()) {
            int ng = std::clamp(microstructure.num_grains, 2, 64);
            std::mt19937 rng(static_cast<unsigned int>(microstructure.random_seed));
            std::uniform_real_distribution<double> u01(0.10, 0.90);
            std::uniform_real_distribution<double> utheta(
                std::min(microstructure.theta_min_deg, microstructure.theta_max_deg),
                std::max(microstructure.theta_min_deg, microstructure.theta_max_deg));

            // Stratified / Poisson-disc-like placement constrained to the solid specimen
            int cols = std::max(1, static_cast<int>(std::ceil(std::sqrt(static_cast<double>(ng)))));
            int rows = std::max(1, static_cast<int>(std::ceil(static_cast<double>(ng) / cols)));
            for (int k = 0; k < ng; ++k) {
                int r = k / cols;
                int c = k % cols;
                std::uniform_real_distribution<double> jx(0.20, 0.80);
                std::uniform_real_distribution<double> jy(0.20, 0.80);
                VoronoiGrainSeed s;
                s.id = k;
                s.cx = std::clamp((c + jx(rng)) / cols, 0.08, 0.92);
                s.cy = std::clamp((r + jy(rng)) / rows, 0.08, 0.92);
                s.theta_deg = utheta(rng);
                seeds.push_back(s);
            }
        }

        // Guarantee every Voronoi seed lies strictly inside the carved solid domain
        // (e.g., never in the empty carved waist of a Dogbone, outside a Disk, or inside a Hole)
        int num_relocated_seeds = 0;
        for (auto& s : seeds) {
            if (!isPointInsideSolid(s.cx * size_x, s.cy * size_y)) {
                ++num_relocated_seeds;
                bool relocated = false;
                // First try pulling horizontally toward x_hat = 0.5 at the seed's y_hat (ideal for Dogbone waist)
                for (int step = 1; step <= 20; ++step) {
                    double t = static_cast<double>(step) / 20.0;
                    double cand_x = s.cx + t * (0.5 - s.cx);
                    double cand_y = s.cy;
                    if (isPointInsideSolid(cand_x * size_x, cand_y * size_y)) {
                        // Step slightly further inward from the boundary rim
                        double inward_t = std::min(1.0, t + 0.12);
                        double inner_x = s.cx + inward_t * (0.5 - s.cx);
                        if (isPointInsideSolid(inner_x * size_x, cand_y * size_y)) {
                            cand_x = inner_x;
                        }
                        s.cx = cand_x;
                        s.cy = cand_y;
                        relocated = true;
                        break;
                    }
                }
                // Next try pulling radially toward (0.5, 0.5) (ideal for Disk or Polygon)
                if (!relocated) {
                    for (int step = 1; step <= 25; ++step) {
                        double t = static_cast<double>(step) / 25.0;
                        double cand_x = s.cx + t * (0.5 - s.cx);
                        double cand_y = s.cy + t * (0.5 - s.cy);
                        if (isPointInsideSolid(cand_x * size_x, cand_y * size_y)) {
                            double inward_t = std::min(1.0, t + 0.10);
                            double inner_x = s.cx + inward_t * (0.5 - s.cx);
                            double inner_y = s.cy + inward_t * (0.5 - s.cy);
                            if (isPointInsideSolid(inner_x * size_x, inner_y * size_y)) {
                                cand_x = inner_x;
                                cand_y = inner_y;
                            }
                            s.cx = cand_x;
                            s.cy = cand_y;
                            relocated = true;
                            break;
                        }
                    }
                }
                // Last resort (e.g. seed inside a central hole): nearest solid point on growing rings
                for (int ring = 1; !relocated && ring <= 40; ++ring) {
                    const double rad = 0.0125 * ring;
                    for (int k = 0; k < 32; ++k) {
                        const double ang = 2.0 * M_PI * k / 32.0;
                        const double tx = s.cx + rad * std::cos(ang);
                        const double ty = s.cy + rad * std::sin(ang);
                        if (isPointInsideSolid(tx * size_x, ty * size_y) &&
                            !isPointNearCutout(tx * size_x, ty * size_y, 1.5 * a)) {
                            s.cx = tx;
                            s.cy = ty;
                            relocated = true;
                            break;
                        }
                    }
                }
                if (!relocated) {
                    --num_relocated_seeds;
                    std::cerr << "[Microstructure] WARNING: Voronoi seed " << s.id
                              << " lies outside the carved solid and could not be relocated\n";
                }
            }
        }
        microstructure.seeds = seeds;
        if (num_relocated_seeds > 0) {
            std::cout << "[Microstructure] Relocated " << num_relocated_seeds
                      << " Voronoi seed(s) from outside the carved solid into the specimen interior\n";
        }

        const int ng = static_cast<int>(seeds.size());
        std::vector<Eigen::Vector2d> seed_phys(ng);
        for (int k = 0; k < ng; ++k) {
            seed_phys[k] = Eigen::Vector2d(seeds[k].cx * size_x, seeds[k].cy * size_y);
        }

        for (int k = 0; k < ng; ++k) {
            appendGrainPoints(k, seeds[k].theta_deg, seed_phys[k].x(), seed_phys[k].y(),
                              [&](double x, double y) {
                                  Eigen::Vector2d p(x, y);
                                  double d_k = (p - seed_phys[k]).squaredNorm();
                                  for (int j = 0; j < ng; ++j) {
                                      if (j == k) continue;
                                      double d_j = (p - seed_phys[j]).squaredNorm();
                                      if (d_j < d_k - 1e-10) return false;
                                      if (std::abs(d_j - d_k) <= 1e-10 && j < k) return false;
                                  }
                                  return true;
                              });
        }

        std::cout << "[Microstructure] Generated 2D Voronoi Polycrystal with " << ng
                  << " grains (raw points before seam filter: " << raw_points.size() << ")\n";
    } else {
        // Fallback SingleCrystal
        appendGrainPoints(0, fallback_theta_deg, 0.5 * size_x, 0.5 * size_y,
                          [](double, double) { return true; });
        return {raw_points, DomainDimensions(size_x, size_y)};
    }

    // ==================== SEAM OVERLAP FILTER (r_ij < d_min * a0) ====================
    // Multi-pass, distance-sorted overlap resolution for Bicrystals and 2D Voronoi Polycrystals
    // (guarantees zero overlapping nodes at grain boundaries and triple junctions).
    const double d_min = std::clamp(microstructure.min_seam_distance_a0, 0.15, 0.98) * a;
    const double d_min_sq = d_min * d_min;
    const double d_same_safe = std::min(d_min, 0.68 * a);
    const double d_same_safe_sq = d_same_safe * d_same_safe;
    const size_t N = raw_points.size();
    std::vector<bool> removed(N, false);
    std::vector<bool> is_merged_seam(N, false);
    int merged_or_removed_count = 0;

    // Spatial hash grid with cell size d_min
    const double cell_h = std::max(1e-4, d_min);
    const double inv_cell = 1.0 / cell_h;
    const int grid_nx = std::max(1, static_cast<int>(std::ceil(size_x * inv_cell)) + 4);
    const int grid_ny = std::max(1, static_cast<int>(std::ceil(size_y * inv_cell)) + 4);
    std::vector<std::vector<size_t>> buckets(static_cast<size_t>(grid_nx * grid_ny));

    auto cellIndex = [&](double x, double y) -> std::pair<int, int> {
        int ix = std::clamp(static_cast<int>(std::floor(x * inv_cell)) + 1, 0, grid_nx - 1);
        int iy = std::clamp(static_cast<int>(std::floor(y * inv_cell)) + 1, 0, grid_ny - 1);
        return {ix, iy};
    };

    auto rebuildBuckets = [&]() {
        for (auto& b : buckets) b.clear();
        for (size_t i = 0; i < N; ++i) {
            if (removed[i]) continue;
            auto [ix, iy] = cellIndex(raw_points[i].coord.x(), raw_points[i].coord.y());
            buckets[static_cast<size_t>(iy * grid_nx + ix)].push_back(i);
        }
    };

    rebuildBuckets();

    // Collect all cross-grain candidate pairs within d_min and sort by distance (closest first)
    struct OverlapPair {
        size_t i;
        size_t j;
        double d2;
    };
    std::vector<OverlapPair> candidate_pairs;
    for (size_t i = 0; i < N; ++i) {
        auto [ix, iy] = cellIndex(raw_points[i].coord.x(), raw_points[i].coord.y());
        for (int dy = -1; dy <= 1; ++dy) {
            int ny_c = iy + dy;
            if (ny_c < 0 || ny_c >= grid_ny) continue;
            for (int dx = -1; dx <= 1; ++dx) {
                int nx_c = ix + dx;
                if (nx_c < 0 || nx_c >= grid_nx) continue;
                const auto& bucket = buckets[static_cast<size_t>(ny_c * grid_nx + nx_c)];
                for (size_t j : bucket) {
                    if (j <= i) continue;
                    if (raw_points[i].grain_id == raw_points[j].grain_id) continue;
                    double d2 = (raw_points[i].coord - raw_points[j].coord).squaredNorm();
                    if (d2 < d_min_sq) {
                        candidate_pairs.push_back({i, j, d2});
                    }
                }
            }
        }
    }
    std::sort(candidate_pairs.begin(), candidate_pairs.end(),
              [](const OverlapPair& a, const OverlapPair& b) { return a.d2 < b.d2; });

    auto isPositionCrowded = [&](const Eigen::Vector2d& pos, size_t ignore_a, size_t ignore_b, int keep_gid) -> bool {
        auto [ix, iy] = cellIndex(pos.x(), pos.y());
        for (int dy = -1; dy <= 1; ++dy) {
            int ny_c = iy + dy;
            if (ny_c < 0 || ny_c >= grid_ny) continue;
            for (int dx = -1; dx <= 1; ++dx) {
                int nx_c = ix + dx;
                if (nx_c < 0 || nx_c >= grid_nx) continue;
                for (size_t k : buckets[static_cast<size_t>(ny_c * grid_nx + nx_c)]) {
                    if (removed[k] || k == ignore_a || k == ignore_b) continue;
                    double d2 = (raw_points[k].coord - pos).squaredNorm();
                    double tol2 = (raw_points[k].grain_id == keep_gid || is_merged_seam[k]) ? d_same_safe_sq : d_min_sq;
                    if (d2 < tol2) return true;
                }
            }
        }
        return false;
    };

    // Pass 1: Resolve sorted cross-grain overlaps (closest pairs first)
    for (const auto& pair : candidate_pairs) {
        size_t i = pair.i;
        size_t j = pair.j;
        if (removed[i] || removed[j]) continue;
        double cur_d2 = (raw_points[i].coord - raw_points[j].coord).squaredNorm();
        if (cur_d2 >= d_min_sq) continue;

        if (microstructure.seam_resolution == SeamOverlapResolution::MergeMidpoint) {
            Eigen::Vector2d mid = 0.5 * (raw_points[i].coord + raw_points[j].coord);
            // Only shift i to midpoint if it doesn't crowd another existing node;
            // otherwise keep i at its lattice position and delete j.
            if (!isPositionCrowded(mid, i, j, raw_points[i].grain_id)) {
                raw_points[i].coord = mid;
                is_merged_seam[i] = true;
                auto [nix, niy] = cellIndex(mid.x(), mid.y());
                buckets[static_cast<size_t>(niy * grid_nx + nix)].push_back(i);
            }
        }
        removed[j] = true;
        ++merged_or_removed_count;
    }

    // Pass 2: Exhaustive cleanup sweep (guarantees zero remaining cross-grain or triple-junction overlaps)
    for (int pass = 0; pass < 3; ++pass) {
        rebuildBuckets();
        int pass_removed = 0;
        for (size_t i = 0; i < N; ++i) {
            if (removed[i]) continue;
            auto [ix, iy] = cellIndex(raw_points[i].coord.x(), raw_points[i].coord.y());
            for (int dy = -1; dy <= 1; ++dy) {
                int ny_c = iy + dy;
                if (ny_c < 0 || ny_c >= grid_ny) continue;
                for (int dx = -1; dx <= 1; ++dx) {
                    int nx_c = ix + dx;
                    if (nx_c < 0 || nx_c >= grid_nx) continue;
                    for (size_t j : buckets[static_cast<size_t>(ny_c * grid_nx + nx_c)]) {
                        if (j <= i || removed[j]) continue;
                        double d2 = (raw_points[i].coord - raw_points[j].coord).squaredNorm();
                        double cutoff2 = (raw_points[i].grain_id != raw_points[j].grain_id && !is_merged_seam[i] && !is_merged_seam[j])
                                         ? d_min_sq
                                         : d_same_safe_sq;
                        if (d2 < cutoff2) {
                            removed[j] = true;
                            ++merged_or_removed_count;
                            ++pass_removed;
                        }
                    }
                }
            }
        }
        if (pass_removed == 0) break;
    }

    const double phi_gb_rad = microstructure.gb_angle_deg * M_PI / 180.0;
    const double nx_gb_sep = std::cos(phi_gb_rad);
    const double ny_gb_sep = std::sin(phi_gb_rad);
    const double delta_perp = (microstructure.mode == MicrostructureMode::Bicrystal)
        ? (microstructure.gb_normal_gap_a0 * a)
        : 0.0;

    std::vector<Point2D> filtered_points;
    filtered_points.reserve(N - static_cast<size_t>(merged_or_removed_count));
    for (size_t i = 0; i < N; ++i) {
        if (!removed[i]) {
            Point2D pt = raw_points[i];
            if (std::abs(delta_perp) > 1e-12 && !is_merged_seam[i]) {
                double sgn = (pt.grain_id == 0) ? -0.5 : +0.5;
                pt.coord.x() += sgn * delta_perp * nx_gb_sep;
                pt.coord.y() += sgn * delta_perp * ny_gb_sep;
            }
            filtered_points.push_back(pt);
        }
    }

    std::cout << "[Microstructure] Seam filter (d_min = " << microstructure.min_seam_distance_a0
              << " a0, mode = "
              << (microstructure.seam_resolution == SeamOverlapResolution::MergeMidpoint ? "MERGE_MIDPOINT" : "DELETE_SECOND")
              << ", gap_perp = " << microstructure.gb_normal_gap_a0
              << " a0): resolved " << merged_or_removed_count << " overlapping seam nodes -> "
              << filtered_points.size() << " atoms retained.\n";

    return {filtered_points, DomainDimensions(size_x + std::abs(delta_perp * nx_gb_sep),
                                              size_y + std::abs(delta_perp * ny_gb_sep))};
}

std::pair<std::vector<Point2D>, DomainDimensions>
MacroscopicGeometryConfig::generatePeriodicDislocationLattice(
    const std::string& lattice_type,
    int nx, int ny,
    double lattice_a0,
    int potential_mode) {

    const double a = (lattice_a0 > 0.0) ? lattice_a0 : 1.0;
    const bool is_triangular = (lattice_type == "triangular");
    const double sqrt3 = std::sqrt(3.0);
    const double h_y = is_triangular ? (a * sqrt3 * 0.5) : a;

    int effective_ny = ny;
    if (is_triangular && (effective_ny % 2 != 0)) {
        effective_ny++;
    }

    const double Lx = nx * a;
    const double Ly = effective_ny * h_y;

    int M = std::max(1, microstructure.dislocation_pair_count);
    std::vector<SeededDislocation> active_dislocs;

    if (!microstructure.dislocations.empty()) {
        for (const auto& d : microstructure.dislocations) {
            if (d.enabled && d.glide_plane_row >= 1 && d.glide_plane_row < effective_ny - 1) {
                active_dislocs.push_back(d);
            }
        }
    }

    int sum_b = 0;
    for (const auto& d : active_dislocs) sum_b += d.burgers_sign;

    if (active_dislocs.size() < 2 || active_dislocs.size() % 2 != 0 || sum_b != 0) {
        active_dislocs.clear();
        int total_planes = 2 * M;

        int min_row = 2;
        int max_row = std::max(min_row + total_planes, effective_ny - 3);

        std::vector<int> plane_rows;
        bool use_equispaced_rows = (
            microstructure.dislocation_layout == "alternating_staggered" ||
            microstructure.dislocation_layout == "alternating_random_x" ||
            microstructure.dislocation_layout == "equispaced_random_signs" ||
            microstructure.dislocation_layout == "grouped_dipole_band"
        );
        if (use_equispaced_rows) {
            for (int k = 0; k < total_planes; ++k) {
                double ideal_row = (static_cast<double>(k) + 0.5) * (static_cast<double>(effective_ny) / static_cast<double>(total_planes)) - 0.5;
                int r = std::clamp(static_cast<int>(std::round(ideal_row)), min_row, max_row);
                if (!plane_rows.empty() && r <= plane_rows.back()) {
                    r = plane_rows.back() + 1;
                }
                plane_rows.push_back(r);
            }
        } else {
            std::mt19937 rng(microstructure.random_seed);
            int min_space = std::max(2, microstructure.min_plane_spacing_rows);
            int attempts = 0;
            bool success = false;
            while (attempts++ < 1000 && !success) {
                std::vector<int> candidates;
                for (int r = min_row; r <= max_row; ++r) candidates.push_back(r);
                std::shuffle(candidates.begin(), candidates.end(), rng);
                plane_rows.clear();
                for (int r : candidates) {
                    bool ok = true;
                    for (int existing : plane_rows) {
                        int diff = std::abs(r - existing);
                        int p_diff = std::min(diff, effective_ny - diff);
                        if (p_diff < min_space) { ok = false; break; }
                    }
                    if (ok) {
                        plane_rows.push_back(r);
                        if (plane_rows.size() == static_cast<size_t>(total_planes)) break;
                    }
                }
                if (plane_rows.size() == static_cast<size_t>(total_planes)) {
                    std::sort(plane_rows.begin(), plane_rows.end());
                    success = true;
                }
            }
            if (!success) {
                plane_rows.clear();
                for (int k = 0; k < total_planes; ++k) {
                    double ideal_row = (static_cast<double>(k) + 0.5) * (static_cast<double>(effective_ny) / static_cast<double>(total_planes)) - 0.5;
                    int r = std::clamp(static_cast<int>(std::round(ideal_row)), min_row, max_row);
                    if (!plane_rows.empty() && r <= plane_rows.back()) r = plane_rows.back() + 1;
                    plane_rows.push_back(r);
                }
            }
        }

        std::mt19937 x_rng(microstructure.random_seed + 101);
        std::uniform_real_distribution<double> x_dist(0.18, 0.82);

        std::vector<int> signs(total_planes);
        if (microstructure.dislocation_layout == "grouped_dipole_band") {
            for (int k = 0; k < total_planes; ++k) signs[k] = (k < M) ? 1 : -1;
        } else if (microstructure.dislocation_layout == "random_planes_random_signs" ||
                   microstructure.dislocation_layout == "equispaced_random_signs") {
            for (int k = 0; k < total_planes; ++k) signs[k] = (k < M) ? 1 : -1;
            std::shuffle(signs.begin(), signs.end(), x_rng);
        } else {
            for (int k = 0; k < total_planes; ++k) signs[k] = (k % 2 == 0) ? 1 : -1;
        }

        for (int k = 0; k < total_planes; ++k) {
            SeededDislocation d;
            d.id = k;
            d.glide_plane_row = plane_rows[k];
            d.burgers_sign = signs[k];
            if (microstructure.dislocation_layout == "alternating_staggered") {
                d.x_frac = (k % 2 == 0) ? 0.35 : 0.65;
            } else {
                d.x_frac = x_dist(x_rng);
            }
            d.enabled = true;
            active_dislocs.push_back(d);
        }
    }

    std::sort(active_dislocs.begin(), active_dislocs.end(), [](const SeededDislocation& a, const SeededDislocation& b) {
        return a.glide_plane_row < b.glide_plane_row;
    });

    std::vector<int> n_strip(effective_ny, 0);
    for (int j = 0; j < effective_ny; ++j) {
        int sum_s = 0;
        for (const auto& d : active_dislocs) {
            if (d.glide_plane_row < j) {
                sum_s += d.burgers_sign;
            }
        }
        n_strip[j] = sum_s;
    }

    DislocationSolutionParams disloc_params =
        (potential_mode == 3)
            ? DislocationSolutionParams::fromCurrentZanzottoPotential(lattice_type)
            : DislocationSolutionParams::fromPotentialMode(potential_mode, lattice_type);

    const double core_r0 = std::clamp(microstructure.dislocation_core_radius_a0, 0.1, 5.0) * a;
    constexpr int Mx = 4;
    constexpr int My = 2;

    // Precompute reference vertical offsets u_ref(mx * Lx, 0; +a) once
    std::array<double, Mx + 1> u_ref_y_pos{};
    for (int mx = 1; mx <= Mx; ++mx) {
        u_ref_y_pos[mx] = NonSingularDisplacement::calculateEdgeDisplacement(
            Eigen::Vector2d(mx * Lx, 0.0), Eigen::Vector2d(a, 0.0), core_r0, disloc_params).y();
    }

    auto evalXPeriodicSinglePlane = [&](double xi, double eta, int b_sign) -> Eigen::Vector2d {
        const Eigen::Vector2d b_vec(b_sign * a, 0.0);
        const double b = b_vec.x();
        const double sgn_eta = (eta >= 0.0) ? 1.0 : -1.0;
        const double shift_x = sgn_eta * 0.25 * b;

        // Primary image (mx = 0)
        Eigen::Vector2d u0 = NonSingularDisplacement::calculateEdgeDisplacement(
            Eigen::Vector2d(xi, eta), b_vec, core_r0, disloc_params);
        double sum_ux = u0.x() - shift_x;
        double sum_uy = u0.y();

        // Symmetric periodic image pairs (+mx, -mx):
        // ux_centered is strictly odd in x; uy is strictly even in x.
        for (int mx = 1; mx <= Mx; ++mx) {
            Eigen::Vector2d u_pos = NonSingularDisplacement::calculateEdgeDisplacement(
                Eigen::Vector2d(xi + mx * Lx, eta), b_vec, core_r0, disloc_params);
            Eigen::Vector2d u_neg = NonSingularDisplacement::calculateEdgeDisplacement(
                Eigen::Vector2d(xi - mx * Lx, eta), b_vec, core_r0, disloc_params);

            sum_ux += (u_pos.x() - shift_x) + (u_neg.x() - shift_x);
            sum_uy += (u_pos.y() + u_neg.y() - 2.0 * b_sign * u_ref_y_pos[mx]);
        }

        // Exact telescoping boundary jump at xi = +Lx/2 vs xi = -Lx/2:
        // Since ux_centered is odd in x, delta_ux = 2 * ux_centered((Mx + 0.5) * Lx, eta), and delta_uy = 0.
        Eigen::Vector2d u_end = NonSingularDisplacement::calculateEdgeDisplacement(
            Eigen::Vector2d((Mx + 0.5) * Lx, eta), b_vec, core_r0, disloc_params);
        double delta_ux = 2.0 * (u_end.x() - shift_x);

        return Eigen::Vector2d(sum_ux - delta_ux * (xi / Lx), sum_uy);
    };

    std::vector<std::vector<Point2D>> row_points(effective_ny);

    #pragma omp parallel for schedule(dynamic)
    for (int j = 0; j < effective_ny; ++j) {
        const int Nj = std::max(4, nx + n_strip[j]);
        const double row_y = j * h_y;
        const double stagger_x = is_triangular ? (0.5 * (j % 2)) : 0.0;
        auto& r_pts = row_points[j];
        r_pts.reserve(Nj);

        for (int i = 0; i < Nj; ++i) {
            // Half-integer offset (i + 0.5) keeps edge atoms away from the x = 0 / Lx branch cut
            double bx = (static_cast<double>(i) + 0.5 + stagger_x) * (Lx / static_cast<double>(Nj));
            while (bx >= Lx) bx -= Lx;
            double by = row_y;

            Eigen::Vector2d total_u = Eigen::Vector2d::Zero();
            Eigen::Vector2d delta_v_y = Eigen::Vector2d::Zero();

            for (const auto& d : active_dislocs) {
                double disloc_x = d.x_frac * Lx;
                double disloc_y = (d.glide_plane_row + 0.5) * h_y;

                double dx_raw = bx - disloc_x;
                double xi = dx_raw - Lx * std::round(dx_raw / Lx);

                for (int my = -My; my <= My; ++my) {
                    double eta = by - (disloc_y + my * Ly);
                    total_u += evalXPeriodicSinglePlane(xi, eta, d.burgers_sign);
                }

                // Exact vertical telescoping residual between y = Ly and y = 0:
                Eigen::Vector2d v_top = evalXPeriodicSinglePlane(xi, (My + 1.0) * Ly - disloc_y, d.burgers_sign);
                Eigen::Vector2d v_bot = evalXPeriodicSinglePlane(xi, -My * Ly - disloc_y, d.burgers_sign);
                delta_v_y += (v_top - v_bot);
            }

            total_u -= delta_v_y * ((by - 0.5 * Ly) / Ly);

            double final_x = bx + total_u.x();
            double final_y = by + total_u.y();

            while (final_y < 0.0) final_y += Ly;
            while (final_y >= Ly) final_y -= Ly;

            r_pts.emplace_back(final_x, final_y, 0, 0.0);
        }
    }

    std::vector<Point2D> points;
    points.reserve(nx * effective_ny + 2 * M);
    for (int j = 0; j < effective_ny; ++j) {
        points.insert(points.end(), row_points[j].begin(), row_points[j].end());
    }

    // Exact periodic domain dimensions so that domain_dims + offsets == (Lx, Ly) to machine precision
    const double offset_x = is_triangular ? (0.5 * a) : a;
    const double offset_y = h_y;
    DomainDimensions dims(Lx - offset_x, Ly - offset_y);

    std::cout << "[PeriodicDislocations] Generated " << points.size() << " atoms with "
              << active_dislocs.size() << " glide-plane dislocations (M = " << (active_dislocs.size() / 2)
              << " pairs, ∑b = 0, lattice = " << lattice_type
              << ", periodic box = " << Lx << " x " << Ly << ")\n";

    return {points, dims};
}

void MacroscopicGeometryConfig::setupIndenterBounds(double lattice_a0,
                                                    const std::vector<Point2D>* carved_points) {
    a0 = (lattice_a0 > 1e-8) ? lattice_a0 : 1.0;
    indenter.physical_radius = std::max(1e-4, indenter.radius_a0 * a0);
    indenter.physical_max_depth = std::max(0.0, indenter.max_depth_a0 * a0);
    indenter.physical_cx = x_min + indenter.cx * Lx;
    indenter.physical_cushion_offset = std::max(0.0, indenter.cushion_offset_a0 * a0);

    // Determine initial top surface height at indenter.physical_cx from actual carved crystal sites
    double y_top = carved_y_max;
    if (carved_points != nullptr && !carved_points->empty()) {
        double search_half_w = std::max(2.5 * a0, 0.5 * indenter.physical_radius);
        double max_y_local = -1e30;
        for (const auto& pt : *carved_points) {
            if (std::abs(pt.coord.x() - indenter.physical_cx) <= search_half_w) {
                max_y_local = std::max(max_y_local, pt.coord.y());
            }
        }
        if (max_y_local > -1e20) {
            y_top = max_y_local;
        }
    }
    indenter.physical_y_top = y_top;
}

Eigen::Vector2d MacroscopicGeometryConfig::getIndenterPosition(double alpha) const {
    double depth = std::max(0.0, alpha) * indenter.physical_max_depth;
    double y_center = indenter.physical_y_top + indenter.physical_radius - depth;
    return Eigen::Vector2d(indenter.physical_cx, y_center);
}

MacroscopicGeometryConfig::IndenterPointContact
MacroscopicGeometryConfig::evaluateIndenterContact(const Eigen::Vector2d& p,
                                                   double alpha,
                                                   double extra_margin) const {
    IndenterPointContact res;
    if (!isIndenterActive()) return res;

    const double depth = std::max(0.0, alpha) * indenter.physical_max_depth;
    const double y_tip = indenter.physical_y_top - depth - extra_margin;
    const double R_eff = indenter.physical_radius + extra_margin;
    const double kc = indenter.penalty_stiffness;
    const double cx = indenter.physical_cx;

    switch (indenter.shape) {
    case IndenterShapeType::Cylinder: {
        Eigen::Vector2d c(cx, y_tip + R_eff);
        if (p.y() < c.y()) {
            Eigen::Vector2d diff = p - c;
            double d = diff.norm();
            if (d < R_eff) {
                res.in_contact = true;
                res.penetration = R_eff - d;
                res.normal = (d > 1e-12) ? (diff / d) : Eigen::Vector2d(0.0, -1.0);
                res.surface_point = c + R_eff * res.normal;
                res.hessian = kc * (res.normal * res.normal.transpose());
            }
        }
        break;
    }
    case IndenterShapeType::FlatPunch: {
        const double r_c = std::min(0.5 * a0, 0.15 * R_eff);
        const double w_flat = std::max(0.0, R_eff - r_c);
        const double dx = p.x() - cx;
        const double adx = std::abs(dx);
        if (p.y() <= indenter.physical_y_top + 2.5 * a0) {
            if (adx <= w_flat) {
                if (p.y() > y_tip) {
                    res.in_contact = true;
                    res.penetration = p.y() - y_tip;
                    res.normal = Eigen::Vector2d(0.0, -1.0);
                    res.surface_point = Eigen::Vector2d(p.x(), y_tip);
                    res.hessian = kc * (res.normal * res.normal.transpose());
                }
            } else if (adx < R_eff) {
                Eigen::Vector2d c_corner(cx + (dx >= 0.0 ? w_flat : -w_flat), y_tip + r_c);
                if (p.y() < c_corner.y()) {
                    Eigen::Vector2d diff = p - c_corner;
                    double d = diff.norm();
                    if (d < r_c) {
                        res.in_contact = true;
                        res.penetration = r_c - d;
                        res.normal = (d > 1e-12) ? (diff / d) : Eigen::Vector2d(0.0, -1.0);
                        res.surface_point = c_corner + r_c * res.normal;
                        res.hessian = kc * (res.normal * res.normal.transpose());
                    }
                } else {
                    double u_c = adx - w_flat;
                    double y_arc = c_corner.y() - std::sqrt(std::max(0.0, r_c * r_c - u_c * u_c));
                    res.in_contact = true;
                    res.penetration = p.y() - y_arc;
                    res.normal = Eigen::Vector2d(0.0, -1.0);
                    res.surface_point = Eigen::Vector2d(p.x(), y_arc);
                    res.hessian = kc * (res.normal * res.normal.transpose());
                }
            }
        }
        break;
    }
    case IndenterShapeType::Wedge: {
        const double theta_rad = std::clamp(indenter.wedge_half_angle_deg, 15.0, 85.0) * M_PI / 180.0;
        const double sin_t = std::sin(theta_rad);
        const double cos_t = std::cos(theta_rad);
        const double r_tip = 0.25 * a0;
        const Eigen::Vector2d c_tip(cx, y_tip + r_tip);
        const double dx = p.x() - cx;
        const double adx = std::abs(dx);
        const double s = (dx >= 0.0) ? 1.0 : -1.0;

        if (adx <= R_eff && p.y() <= indenter.physical_y_top + 2.5 * a0) {
            const double qy = p.y() - c_tip.y();
            const double tau = adx * sin_t + qy * cos_t;
            if (tau <= 0.0) {
                if (qy < 0.0) {
                    Eigen::Vector2d diff(dx, qy);
                    double d = diff.norm();
                    if (d < r_tip) {
                        res.in_contact = true;
                        res.penetration = r_tip - d;
                        res.normal = (d > 1e-12) ? (diff / d) : Eigen::Vector2d(0.0, -1.0);
                        res.surface_point = c_tip + r_tip * res.normal;
                        res.hessian = kc * (res.normal * res.normal.transpose());
                    }
                } else {
                    res.in_contact = true;
                    res.penetration = p.y() - y_tip;
                    res.normal = Eigen::Vector2d(0.0, -1.0);
                    res.surface_point = Eigen::Vector2d(p.x(), y_tip);
                    res.hessian = kc * (res.normal * res.normal.transpose());
                }
            } else {
                const double q_dot_u = adx * cos_t - qy * sin_t;
                const double pen = r_tip - q_dot_u;
                if (pen > 0.0) {
                    res.in_contact = true;
                    res.penetration = pen;
                    res.normal = Eigen::Vector2d(s * cos_t, -sin_t);
                    res.surface_point = p + pen * res.normal;
                    res.hessian = kc * (res.normal * res.normal.transpose());
                }
            }
        }
        break;
    }
    }

    return res;
}

double MacroscopicGeometryConfig::computeIndenterPenaltyEnergy(
    const std::vector<Point2D>& points, double alpha) const {
    if (!isIndenterActive()) return 0.0;
    const double kc = indenter.penalty_stiffness;
    double total_energy = 0.0;

    for (const auto& pt : points) {
        auto contact = evaluateIndenterContact(pt.coord, alpha, 0.0);
        if (contact.in_contact && contact.penetration > 0.0) {
            total_energy += 0.5 * kc * contact.penetration * contact.penetration;
        }
    }
    return total_energy;
}

double MacroscopicGeometryConfig::accumulateIndenterPenalty(
    const alglib::real_1d_array& x,
    const std::vector<std::pair<int, int>>& interior_mapping,
    alglib::real_1d_array& grad,
    double alpha) const {
    if (!isIndenterActive()) return 0.0;
    const double kc = indenter.penalty_stiffness;
    const int n_vars = static_cast<int>(interior_mapping.size());
    double total_energy = 0.0;

    for (int i = 0; i < n_vars; ++i) {
        int dof_idx = interior_mapping[i].second;
        if (dof_idx < 0 || dof_idx >= n_vars) continue;
        Eigen::Vector2d p(x[dof_idx], x[n_vars + dof_idx]);
        auto contact = evaluateIndenterContact(p, alpha, 0.0);
        if (contact.in_contact && contact.penetration > 0.0) {
            total_energy += 0.5 * kc * contact.penetration * contact.penetration;
            // Conservative force: F = -grad(E) -> dE/dx = -kc * penetration * normal
            Eigen::Vector2d dE_dx = -kc * contact.penetration * contact.normal;
            grad[dof_idx] += dE_dx.x();
            grad[n_vars + dof_idx] += dE_dx.y();
        }
    }
    return total_energy;
}

double MacroscopicGeometryConfig::computeIndenterVerticalForce(
    const std::vector<Point2D>& points, double alpha) const {
    if (!isIndenterActive()) return 0.0;
    const double kc = indenter.penalty_stiffness;
    double total_Fy = 0.0;

    for (const auto& pt : points) {
        auto contact = evaluateIndenterContact(pt.coord, alpha, 0.0);
        if (contact.in_contact && contact.penetration > 0.0) {
            // Upward vertical reaction force exerted by substrate on indenter:
            // contact.normal.y() is negative (pointing into substrate), so -normal.y() > 0
            total_Fy += kc * contact.penetration * (-contact.normal.y());
        }
    }
    return total_Fy;
}

void MacroscopicGeometryConfig::updateIndenterContactKinematics(
    std::vector<Point2D>& current_points,
    const std::vector<Point2D>& points_ref,
    const std::vector<int>& baseline_bc_role,
    std::vector<int>& node_bc_role,
    std::vector<int>& out_contact_atoms,
    double alpha) const {
    out_contact_atoms.clear();
    if (!isIndenterActive()) return;

    // Restore baseline BC roles first so user-defined Fixed/Driven BCs are preserved
    // while previously pinned indenter atoms revert to Free before re-detecting contact
    if (baseline_bc_role.size() == current_points.size()) {
        node_bc_role = baseline_bc_role;
    }

    const double y_cutoff = carved_y_min + 0.4 * std::max(1e-8, carved_y_max - carved_y_min);
    std::vector<int> direct_contact;
    for (size_t i = 0; i < current_points.size(); ++i) {
        if (node_bc_role[i] != 0) continue; // Only unconstrained nodes can enter indenter contact
        if (points_ref[i].coord.y() >= y_cutoff) {
            auto contact = evaluateIndenterContact(current_points[i].coord, alpha, 1e-5);
            if (contact.in_contact) {
                // Place atom 0.5e-5 inside the 1e-5 detection band (on the exact surface)
                // so re-evaluating at the same alpha never drops contact due to roundoff.
                current_points[i].coord = contact.surface_point - 0.5e-5 * contact.normal;
                node_bc_role[i] = 1; // Pinned on indenter surface
                direct_contact.push_back(static_cast<int>(i));
            }
        }
    }

    if (indenter.cushion_subsurface && !direct_contact.empty()) {
        double horiz_tol = 1.25 * a0;
        double vert_tol = 1.25 * a0;
        for (int c_idx : direct_contact) {
            const auto& c_pt = current_points[c_idx].coord;
            for (size_t i = 0; i < current_points.size(); ++i) {
                if (static_cast<int>(i) == c_idx || node_bc_role[i] != 0) continue;
                double dx = std::abs(current_points[i].coord.x() - c_pt.x());
                double dy = c_pt.y() - current_points[i].coord.y();
                if (dx <= horiz_tol && dy > 0.0 && dy <= vert_tol) {
                    auto cush = evaluateIndenterContact(
                        current_points[i].coord, alpha, indenter.physical_cushion_offset);
                    if (cush.in_contact) {
                        // Smoothly nudge subsurface atom downward to avoid triangle inversion,
                        // leaving node_bc_role = 0 so it relaxes freely during L-BFGS
                        current_points[i].coord = cush.surface_point;
                    }
                }
            }
        }
    }

    out_contact_atoms = direct_contact;
    std::sort(out_contact_atoms.begin(), out_contact_atoms.end());
    out_contact_atoms.erase(std::unique(out_contact_atoms.begin(), out_contact_atoms.end()), out_contact_atoms.end());
}

void MacroscopicGeometryConfig::rebuildBoundaryDOFs(
    const std::vector<int>& node_bc_role,
    std::vector<std::pair<int, int>>& interior_mapping,
    std::vector<std::pair<int, int>>& full_mapping) const {
    const size_t n_pts = node_bc_role.size();
    interior_mapping.clear();
    full_mapping.clear();
    cached_node_bc_role = node_bc_role;
    has_partial_dof_constraints = false;
    int next_solver_idx = 0;
    for (size_t i = 0; i < n_pts; ++i) {
        int role = node_bc_role[i];
        if (role == 0 || role == 3 || role == 4) {
            interior_mapping.emplace_back(static_cast<int>(i), next_solver_idx);
            full_mapping.emplace_back(static_cast<int>(i), next_solver_idx);
            ++next_solver_idx;
            if (role == 3 || role == 4) {
                has_partial_dof_constraints = true;
            }
        } else {
            full_mapping.emplace_back(static_cast<int>(i), -1);
        }
    }
}

void MacroscopicGeometryConfig::zeroConstrainedComponentGradients(
    const std::vector<std::pair<int, int>>& interior_mapping,
    alglib::real_1d_array& grad) const {
    if (!has_partial_dof_constraints || cached_node_bc_role.empty()) return;
    const int n_vars = static_cast<int>(interior_mapping.size());
    const int n_roles = static_cast<int>(cached_node_bc_role.size());
    for (int i = 0; i < n_vars; ++i) {
        const auto& [point_idx, dof_idx] = interior_mapping[i];
        if (point_idx < 0 || point_idx >= n_roles || dof_idx < 0 || dof_idx >= n_vars) continue;
        const int role = cached_node_bc_role[point_idx];
        if (role == 3) {
            grad[dof_idx] = 0.0; // X-DOF constrained, Y-DOF free
        } else if (role == 4) {
            grad[n_vars + dof_idx] = 0.0; // Y-DOF constrained, X-DOF free
        }
    }
}

double MacroscopicGeometryConfig::computeKinematicContactForce(
    const UserData* userData,
    const std::vector<int>& contact_atoms) const {
    if (!userData || contact_atoms.empty()) return 0.0;
    const double normalisation = userData->calculator.getUnitCellArea();
    std::unordered_set<int> contact_set(contact_atoms.begin(), contact_atoms.end());
    double total_contact_Fy = 0.0;

    for (size_t elem_idx : userData->active_elements) {
        const auto& elem = userData->elements[elem_idx];
        int i0 = elem.getNodeIndex(0);
        int i1 = elem.getNodeIndex(1);
        int i2 = elem.getNodeIndex(2);
        bool touches_contact = (contact_set.count(i0) || contact_set.count(i1) || contact_set.count(i2));
        if (!touches_contact) continue;

        const Eigen::Matrix2d& F = elem.getDeformationGradient();
        Eigen::Matrix2d C = elem.getMetricTensor();
        auto result = reduce_metric(C);
        Eigen::Matrix2d dE_dC = userData->calculator.calculate_derivative(
            result.C_reduced, userData->derivative_function) / normalisation;
        Eigen::Matrix2d P_stress = 2.0 * F * result.m_matrix * dE_dC * result.m_matrix.transpose();

        auto nodal_forces = elem.calculate_nodal_forces(P_stress);
        // Reaction force exerted by crystal substrate UPWARD onto indenter tip:
        // Downward indentation induces negative dE/dy on contact nodes,
        // so upward contact reaction force is -nodal_forces[...].y()
        if (contact_set.count(i0)) total_contact_Fy -= nodal_forces[0].y();
        if (contact_set.count(i1)) total_contact_Fy -= nodal_forces[1].y();
        if (contact_set.count(i2)) total_contact_Fy -= nodal_forces[2].y();
    }
    return std::max(0.0, total_contact_Fy);
}

void MacroscopicGeometryConfig::setPhysicalBoundsFromPoints(
    const std::vector<Point2D>& raw_points, double lattice_a0) {
    a0 = (lattice_a0 > 1e-8) ? lattice_a0 : 1.0;
    if (raw_points.empty()) {
        x_min = 0.0;
        y_min = 0.0;
        Lx = 100.0;
        Ly = 100.0;
        carved_x_min = 0.0;
        carved_x_max = 100.0;
        carved_y_min = 0.0;
        carved_y_max = 100.0;
        return;
    }
    double min_x = raw_points[0].coord.x();
    double max_x = raw_points[0].coord.x();
    double min_y = raw_points[0].coord.y();
    double max_y = raw_points[0].coord.y();
    for (const auto& p : raw_points) {
        min_x = std::min(min_x, p.coord.x());
        max_x = std::max(max_x, p.coord.x());
        min_y = std::min(min_y, p.coord.y());
        max_y = std::max(max_y, p.coord.y());
    }
    x_min = min_x;
    y_min = min_y;
    Lx = std::max(1e-8, max_x - min_x);
    Ly = std::max(1e-8, max_y - min_y);
    carved_x_min = min_x;
    carved_x_max = max_x;
    carved_y_min = min_y;
    carved_y_max = max_y;
    setupIndenterBounds(lattice_a0, &raw_points);
}

int MacroscopicGeometryConfig::findCutoutNear(double x, double y, double margin_phys) const {
    const double min_dim = std::min(Lx, Ly);
    const double x_hat = (x - x_min) / Lx;
    const double y_hat = (y - y_min) / Ly;

    for (size_t ci = 0; ci < cutouts.size(); ++ci) {
        const auto& c = cutouts[ci];
        if (!c.enabled) continue;
        const double cx_phys = x_min + c.cx * Lx;
        const double cy_phys = y_min + c.cy * Ly;
        const double dx = x - cx_phys;
        const double dy = y - cy_phys;

        switch (c.type) {
        case CutoutShapeType::Circle: {
            double r_phys = c.rx * min_dim + margin_phys;
            if (dx * dx + dy * dy < r_phys * r_phys) return static_cast<int>(ci);
            break;
        }
        case CutoutShapeType::Ellipse: {
            double phi = c.rotation_deg * M_PI / 180.0;
            double cos_p = std::cos(-phi);
            double sin_p = std::sin(-phi);
            double lx = cos_p * dx - sin_p * dy;
            double ly = sin_p * dx + cos_p * dy;
            double rx_phys = std::max(1e-8, c.rx * Lx + margin_phys);
            double ry_phys = std::max(1e-8, c.ry * Ly + margin_phys);
            if ((lx * lx) / (rx_phys * rx_phys) + (ly * ly) / (ry_phys * ry_phys) < 1.0) return static_cast<int>(ci);
            break;
        }
        case CutoutShapeType::Rectangle: {
            double phi = c.rotation_deg * M_PI / 180.0;
            double cos_p = std::cos(-phi);
            double sin_p = std::sin(-phi);
            double lx = cos_p * dx - sin_p * dy;
            double ly = sin_p * dx + cos_p * dy;
            double hw_phys = c.rx * Lx + margin_phys;
            double hh_phys = c.ry * Ly + margin_phys;
            if (std::abs(lx) < hw_phys && std::abs(ly) < hh_phys) return static_cast<int>(ci);
            break;
        }
        case CutoutShapeType::EdgeNotch: {
            // Notch from one of the 4 outer edges
            double depth_norm = c.rx;
            double half_w_norm = c.ry;
            double margin_x = margin_phys / Lx;
            double margin_y = margin_phys / Ly;
            if (c.notch_edge == "left") {
                double u = x_hat; // distance from left edge
                double v = std::abs(y_hat - c.cy);
                double d_eff = depth_norm + margin_x;
                double w_eff = std::max(1e-6, half_w_norm + margin_y);
                if (u <= d_eff) {
                    if (c.notch_profile == "u_notch") {
                        if (v < w_eff) return static_cast<int>(ci);
                    } else {
                        // V-notch tapering linearly from half-width w_eff at u=0 to 0 at u=d_eff
                        if (v < w_eff * (1.0 - u / std::max(1e-6, d_eff))) return static_cast<int>(ci);
                    }
                }
            } else if (c.notch_edge == "right") {
                double u = 1.0 - x_hat;
                double v = std::abs(y_hat - c.cy);
                double d_eff = depth_norm + margin_x;
                double w_eff = std::max(1e-6, half_w_norm + margin_y);
                if (u <= d_eff) {
                    if (c.notch_profile == "u_notch") {
                        if (v < w_eff) return static_cast<int>(ci);
                    } else {
                        if (v < w_eff * (1.0 - u / std::max(1e-6, d_eff))) return static_cast<int>(ci);
                    }
                }
            } else if (c.notch_edge == "bottom") {
                double u = y_hat;
                double v = std::abs(x_hat - c.cx);
                double d_eff = depth_norm + margin_y;
                double w_eff = std::max(1e-6, half_w_norm + margin_x);
                if (u <= d_eff) {
                    if (c.notch_profile == "u_notch") {
                        if (v < w_eff) return static_cast<int>(ci);
                    } else {
                        if (v < w_eff * (1.0 - u / std::max(1e-6, d_eff))) return static_cast<int>(ci);
                    }
                }
            } else if (c.notch_edge == "top") {
                double u = 1.0 - y_hat;
                double v = std::abs(x_hat - c.cx);
                double d_eff = depth_norm + margin_y;
                double w_eff = std::max(1e-6, half_w_norm + margin_x);
                if (u <= d_eff) {
                    if (c.notch_profile == "u_notch") {
                        if (v < w_eff) return static_cast<int>(ci);
                    } else {
                        if (v < w_eff * (1.0 - u / std::max(1e-6, d_eff))) return static_cast<int>(ci);
                    }
                }
            }
            break;
        }
        case CutoutShapeType::PreCrack: {
            Eigen::Vector2d a(x_min + c.x1 * Lx, y_min + c.y1 * Ly);
            Eigen::Vector2d b(x_min + c.x2 * Lx, y_min + c.y2 * Ly);
            double half_t = 0.5 * c.thickness * min_dim + margin_phys;
            if (pointToSegmentDistance(Eigen::Vector2d(x, y), a, b) < half_t) return static_cast<int>(ci);
            break;
        }
        case CutoutShapeType::Polygon: {
            if (isPointInPolygon(x_hat, y_hat, c.vertices)) return static_cast<int>(ci);
            if (margin_phys > 0.0 && c.vertices.size() >= 3) {
                double dist = distanceToPolygonEdges(x, y, x_min, y_min, Lx, Ly, c.vertices);
                if (dist <= margin_phys) return static_cast<int>(ci);
            }
            break;
        }
        }
    }
    return -1;
}

bool MacroscopicGeometryConfig::isPointNearCutout(double x, double y, double margin_phys) const {
    return findCutoutNear(x, y, margin_phys) >= 0;
}

Eigen::Vector2d MacroscopicGeometryConfig::cutoutPivot(int idx) const {
    const Eigen::Vector2d specimen_center(x_min + 0.5 * Lx, y_min + 0.5 * Ly);
    if (idx < 0 || idx >= static_cast<int>(cutouts.size())) return specimen_center;
    const auto& c = cutouts[idx];
    switch (c.type) {
    case CutoutShapeType::Circle:
    case CutoutShapeType::Ellipse:
    case CutoutShapeType::Rectangle:
        return Eigen::Vector2d(x_min + c.cx * Lx, y_min + c.cy * Ly);
    case CutoutShapeType::Polygon: {
        if (c.vertices.empty()) return Eigen::Vector2d(x_min + c.cx * Lx, y_min + c.cy * Ly);
        double sx = 0.0, sy = 0.0;
        for (const auto& v : c.vertices) { sx += v.first; sy += v.second; }
        const double n = static_cast<double>(c.vertices.size());
        return Eigen::Vector2d(x_min + (sx / n) * Lx, y_min + (sy / n) * Ly);
    }
    case CutoutShapeType::PreCrack:
        return Eigen::Vector2d(x_min + 0.5 * (c.x1 + c.x2) * Lx, y_min + 0.5 * (c.y1 + c.y2) * Ly);
    case CutoutShapeType::EdgeNotch:
        break;
    }
    return specimen_center;
}

bool MacroscopicGeometryConfig::isPointInsideSolid(double x, double y) const {
    if (!enabled) return true;

    const double x_hat = (x - x_min) / Lx;
    const double y_hat = (y - y_min) / Ly;
    const double min_dim = std::min(Lx, Ly);

    // 1. Check outer domain shape
    switch (outer_shape) {
    case OuterDomainShapeType::Rectangle: {
        if (x_hat < -1e-6 || x_hat > 1.0 + 1e-6 || y_hat < -1e-6 || y_hat > 1.0 + 1e-6) {
            return false;
        }
        break;
    }
    case OuterDomainShapeType::Dogbone: {
        // Vertical tensile dogbone centered at x_hat = 0.5, y_hat = 0.5
        if (x_hat < -1e-6 || x_hat > 1.0 + 1e-6 || y_hat < -1e-6 || y_hat > 1.0 + 1e-6) {
            return false;
        }
        double half_gauge_w = 0.5 * std::clamp(dogbone_gauge_width, 0.15, 0.98);
        double half_gauge_l = 0.5 * std::clamp(dogbone_gauge_length, 0.10, 0.90);
        double r_f = std::max(0.02, dogbone_fillet_radius);
        double dy = std::abs(y_hat - 0.5);
        double half_w_allowed = 0.5;
        if (dy <= half_gauge_l) {
            half_w_allowed = half_gauge_w;
        } else if (dy <= half_gauge_l + r_f) {
            // Circular fillet transition from gauge width to full grip shoulder
            double u = dy - half_gauge_l;
            double fillet_dx = r_f - std::sqrt(std::max(0.0, r_f * r_f - u * u));
            half_w_allowed = std::min(0.5, half_gauge_w + fillet_dx * ((0.5 - half_gauge_w) / r_f));
        } else {
            half_w_allowed = 0.5;
        }
        if (std::abs(x_hat - 0.5) > half_w_allowed + 1e-6) {
            return false;
        }
        break;
    }
    case OuterDomainShapeType::Disk: {
        double cx = x_min + 0.5 * Lx;
        double cy = y_min + 0.5 * Ly;
        double r_phys = disk_radius * min_dim;
        double dx = x - cx;
        double dy = y - cy;
        if (dx * dx + dy * dy > r_phys * r_phys) return false;
        break;
    }
    case OuterDomainShapeType::Polygon: {
        if (outer_polygon.size() >= 3 && !isPointInPolygon(x_hat, y_hat, outer_polygon)) {
            return false;
        }
        break;
    }
    }

    // 2. Check cutouts (holes, notches, cracks)
    if (isPointNearCutout(x, y, 0.0)) {
        return false;
    }

    return true;
}

bool MacroscopicGeometryConfig::shouldRejectTriangle(
    const Eigen::Vector2d& p0,
    const Eigen::Vector2d& p1,
    const Eigen::Vector2d& p2) const {
    if (!enabled) return false;

    // 1. Check triangle centroid
    Eigen::Vector2d centroid = (p0 + p1 + p2) / 3.0;
    if (!isPointInsideSolid(centroid.x(), centroid.y())) {
        return true;
    }

    // 2. Check edge midpoints (prevents bridging across narrow notches, slots, or dogbone fillets)
    // For small crystalline edges (length <= 1.25 * a0), chord sagitta along curved boundaries (e.g. circle)
    // is tiny (~ 0.01 a0). Only check midpoint solid containment for edges > 1.25 * a0 to prevent
    // falsely rejecting genuine crystalline surface facets.
    Eigen::Vector2d m01 = 0.5 * (p0 + p1);
    Eigen::Vector2d m12 = 0.5 * (p1 + p2);
    Eigen::Vector2d m20 = 0.5 * (p2 + p0);
    double l01 = (p1 - p0).norm();
    double l12 = (p2 - p1).norm();
    double l20 = (p0 - p2).norm();

    if (l01 > 1.25 * a0 && !isPointInsideSolid(m01.x(), m01.y())) return true;
    if (l12 > 1.25 * a0 && !isPointInsideSolid(m12.x(), m12.y())) return true;
    if (l20 > 1.25 * a0 && !isPointInsideSolid(m20.x(), m20.y())) return true;

    // 3. Check if any edge crosses a PreCrack slit segment or Polygon boundary
    if (outer_shape == OuterDomainShapeType::Polygon && outer_polygon.size() >= 3) {
        if ((l01 > 1.25 * a0 && segmentIntersectsPolygonEdges(p0, p1, x_min, y_min, Lx, Ly, outer_polygon)) ||
            (l12 > 1.25 * a0 && segmentIntersectsPolygonEdges(p1, p2, x_min, y_min, Lx, Ly, outer_polygon)) ||
            (l20 > 1.25 * a0 && segmentIntersectsPolygonEdges(p2, p0, x_min, y_min, Lx, Ly, outer_polygon))) {
            return true;
        }
    }

    for (const auto& c : cutouts) {
        if (!c.enabled) continue;
        if (c.type == CutoutShapeType::PreCrack) {
            Eigen::Vector2d a(x_min + c.x1 * Lx, y_min + c.y1 * Ly);
            Eigen::Vector2d b(x_min + c.x2 * Lx, y_min + c.y2 * Ly);
            if (segmentsIntersect(p0, p1, a, b) ||
                segmentsIntersect(p1, p2, a, b) ||
                segmentsIntersect(p2, p0, a, b)) {
                return true;
            }
        } else if (c.type == CutoutShapeType::Polygon && c.vertices.size() >= 3) {
            if (segmentIntersectsPolygonEdges(p0, p1, x_min, y_min, Lx, Ly, c.vertices) ||
                segmentIntersectsPolygonEdges(p1, p2, x_min, y_min, Lx, Ly, c.vertices) ||
                segmentIntersectsPolygonEdges(p2, p0, x_min, y_min, Lx, Ly, c.vertices)) {
                return true;
            }
        }
    }

    return false;
}

std::vector<Point2D> MacroscopicGeometryConfig::carveLatticePoints(
    const std::vector<Point2D>& raw_points, double lattice_a0) {
    setPhysicalBoundsFromPoints(raw_points, lattice_a0);
    if (!enabled) return raw_points;

    std::vector<Point2D> carved;
    carved.reserve(raw_points.size());
    for (const auto& pt : raw_points) {
        if (isPointInsideSolid(pt.coord.x(), pt.coord.y())) {
            carved.push_back(pt);
        }
    }

    if (!carved.empty()) {
        double c_min_x = carved[0].coord.x();
        double c_max_x = carved[0].coord.x();
        double c_min_y = carved[0].coord.y();
        double c_max_y = carved[0].coord.y();
        for (const auto& p : carved) {
            c_min_x = std::min(c_min_x, p.coord.x());
            c_max_x = std::max(c_max_x, p.coord.x());
            c_min_y = std::min(c_min_y, p.coord.y());
            c_max_y = std::max(c_max_y, p.coord.y());
        }
        carved_x_min = c_min_x;
        carved_x_max = c_max_x;
        carved_y_min = c_min_y;
        carved_y_max = c_max_y;
    }

    setupIndenterBounds(lattice_a0, &carved);

    std::cout << "[MacroscopicGeometry] Carved discrete crystal lattice (fixed a0 = " << a0 << "): kept "
              << carved.size() << " / " << raw_points.size() << " lattice sites ("
              << (raw_points.size() - carved.size()) << " removed by geometry/holes)\n";
    return carved;
}

void MacroscopicGeometryConfig::buildBoundaryDOFs(
    const std::vector<Point2D>& carved_points,
    std::vector<std::pair<int, int>>& interior_mapping,
    std::vector<std::pair<int, int>>& full_mapping,
    std::vector<int>& node_bc_role,
    std::vector<Eigen::Vector2d>& node_bc_velocity) const {

    const size_t n_pts = carved_points.size();
    interior_mapping.clear();
    full_mapping.clear();
    node_bc_role.assign(n_pts, 0);
    node_bc_velocity.assign(n_pts, Eigen::Vector2d::Zero());
    cached_node_rotation_driven.assign(n_pts, 0);
    cached_node_rotation_pivot.assign(n_pts, Eigen::Vector2d::Zero());
    cached_traction_rule_index.clear();
    cached_traction_membership.clear();
    cached_node_ext_force_unit.assign(n_pts, Eigen::Vector2d::Zero());
    cached_node_ext_force_ref.clear();
    has_external_tractions = false;
    traction_anchor_active = false;
    cached_force_grips.clear();
    cached_node_grip.assign(n_pts, -1);
    for (size_t k = 0; k < boundary_rules.size(); ++k) {
        const auto& rk = boundary_rules[k];
        if (rk.enabled && rk.kind == BoundaryConditionKind::RigidForceGrip) {
            ForceGrip g;
            g.rule_index = static_cast<int>(k);
            // A loaded component (t != 0) is always solved for; an unloaded one is held by its lock
            // (no lateral slip) or slides freely at zero force when unlocked
            g.free_x = !rk.lock_x || std::abs(rk.tx) > 0.0;
            g.free_y = !rk.lock_y || std::abs(rk.ty) > 0.0;
            cached_force_grips.push_back(g);
        }
    }
    for (size_t k = 0; k < boundary_rules.size(); ++k) {
        if (boundary_rules[k].enabled && boundary_rules[k].kind == BoundaryConditionKind::PrescribedTraction) {
            cached_traction_rule_index.push_back(static_cast<int>(k));
            cached_traction_membership.emplace_back(n_pts, 0);
        }
    }

    if (n_pts == 0) return;

    // Use carved solid bounds for Top/Bottom/Left/Right edge grips so inset Custom Polygons
    // and Disks properly match their extreme edges
    double eff_x_min = carved_points[0].coord.x();
    double eff_x_max = carved_points[0].coord.x();
    double eff_y_min = carved_points[0].coord.y();
    double eff_y_max = carved_points[0].coord.y();
    for (const auto& pt : carved_points) {
        eff_x_min = std::min(eff_x_min, pt.coord.x());
        eff_x_max = std::max(eff_x_max, pt.coord.x());
        eff_y_min = std::min(eff_y_min, pt.coord.y());
        eff_y_max = std::max(eff_y_max, pt.coord.y());
    }
    const Eigen::Vector2d center(x_min + 0.5 * Lx, y_min + 0.5 * Ly);

    // Evaluate boundary rules for each node, accumulating independent X/Y locks so
    // intersecting roller edges (e.g. RollerY on Left + RollerX on Bottom) compose at corners.
    for (size_t i = 0; i < n_pts; ++i) {
        const double x = carved_points[i].coord.x();
        const double y = carved_points[i].coord.y();
        const double x_hat = (x - x_min) / Lx;
        const double y_hat = (y - y_min) / Ly;

        bool locked_x = false;
        bool locked_y = false;
        bool driven_x = false;
        bool driven_y = false;
        bool rotation_driven = false;
        Eigen::Vector2d rotation_pivot = Eigen::Vector2d::Zero();
        Eigen::Vector2d vel = Eigen::Vector2d::Zero();
        int grip_of_node = -1;

        for (const auto& rule : boundary_rules) {
            if (!rule.enabled) continue;
            const double strip = std::max(0.35, rule.thickness_a0) * a0;
            bool matches = false;

            const double disk_cap_strip = (outer_shape == OuterDomainShapeType::Disk)
                ? std::max(strip, 0.075 * disk_radius * std::min(Lx, Ly))
                : strip;

            switch (rule.region) {
            case BoundaryRegionType::TopEdge:
                // When the Nanoindenter is active, keep the top surface free for contact
                if (isIndenterActive() && rule.kind != BoundaryConditionKind::FreeSurface) {
                    matches = false;
                } else {
                    matches = (y >= eff_y_max - disk_cap_strip);
                }
                break;
            case BoundaryRegionType::BottomEdge:
                matches = (y <= eff_y_min + disk_cap_strip);
                break;
            case BoundaryRegionType::LeftEdge:
                matches = (x <= eff_x_min + disk_cap_strip);
                break;
            case BoundaryRegionType::RightEdge:
                matches = (x >= eff_x_max - disk_cap_strip);
                break;
            case BoundaryRegionType::OuterPerimeter:
                if (outer_shape == OuterDomainShapeType::Disk) {
                    double r_phys = disk_radius * std::min(Lx, Ly);
                    double dist = (carved_points[i].coord - center).norm();
                    matches = (dist >= r_phys - strip);
                } else if (outer_shape == OuterDomainShapeType::Polygon && outer_polygon.size() >= 3) {
                    double dist = distanceToPolygonEdges(x, y, x_min, y_min, Lx, Ly, outer_polygon);
                    matches = (dist <= strip);
                } else if (outer_shape == OuterDomainShapeType::Dogbone) {
                    double half_gauge_w = 0.5 * std::clamp(dogbone_gauge_width, 0.15, 0.98);
                    double half_gauge_l = 0.5 * std::clamp(dogbone_gauge_length, 0.10, 0.90);
                    double r_f = std::max(0.02, dogbone_fillet_radius);
                    double dy = std::abs(y_hat - 0.5);
                    double half_w_allowed = 0.5;
                    if (dy <= half_gauge_l) {
                        half_w_allowed = half_gauge_w;
                    } else if (dy <= half_gauge_l + r_f) {
                        double u = dy - half_gauge_l;
                        double fillet_dx = r_f - std::sqrt(std::max(0.0, r_f * r_f - u * u));
                        half_w_allowed = std::min(0.5, half_gauge_w + fillet_dx * ((0.5 - half_gauge_w) / r_f));
                    }
                    double x_edge_dist = (half_w_allowed - std::abs(x_hat - 0.5)) * Lx;
                    matches = (x_edge_dist <= strip || y <= eff_y_min + strip || y >= eff_y_max - strip);
                } else {
                    matches = (x <= eff_x_min + strip || x >= eff_x_max - strip ||
                               y <= eff_y_min + strip || y >= eff_y_max - strip);
                }
                break;
            case BoundaryRegionType::HoleSurfaces:
                matches = isPointNearCutout(x, y, strip);
                break;
            case BoundaryRegionType::CustomBox:
                matches = (x_hat >= rule.box_x_min && x_hat <= rule.box_x_max &&
                           y_hat >= rule.box_y_min && y_hat <= rule.box_y_max);
                break;
            }

            if (!matches) continue;

            // Any later constraining rule replaces an earlier torsion trajectory on this node.
            // Tractions are loads, not kinematic constraints, so they leave it untouched.
            if (rule.kind != BoundaryConditionKind::FreeSurface &&
                rule.kind != BoundaryConditionKind::PrescribedTraction) rotation_driven = false;
            // Likewise a later kinematic rule takes the node out of a force grip
            if (rule.kind != BoundaryConditionKind::FreeSurface &&
                rule.kind != BoundaryConditionKind::PrescribedTraction &&
                rule.kind != BoundaryConditionKind::RigidForceGrip) grip_of_node = -1;

            if (rule.kind == BoundaryConditionKind::RigidForceGrip) {
                // Rigid grip: both DOFs follow the shared grip displacement U (not solver DOFs)
                const int rule_k = static_cast<int>(&rule - boundary_rules.data());
                for (size_t g = 0; g < cached_force_grips.size(); ++g) {
                    if (cached_force_grips[g].rule_index == rule_k) grip_of_node = static_cast<int>(g);
                }
                locked_x = true; locked_y = true;
                driven_x = true; driven_y = true;
                vel = Eigen::Vector2d::Zero();
            } else if (rule.kind == BoundaryConditionKind::PrescribedTraction) {
                // Dead load: the node stays a free DOF; finalizeBoundaryTractions assigns its force
                const int rule_k = static_cast<int>(&rule - boundary_rules.data());
                for (size_t m = 0; m < cached_traction_rule_index.size(); ++m) {
                    if (cached_traction_rule_index[m] == rule_k) cached_traction_membership[m][i] = 1;
                }
            } else if (rule.kind == BoundaryConditionKind::FreeSurface) {
                // FreeSurface (Neumann natural BC) marks a traction-free boundary.
                // It MUST NEVER overwrite an already-constrained Dirichlet condition!
            } else if (rule.kind == BoundaryConditionKind::RollerX) {
                // Roller X: Free ux, Fix uy = 0
                locked_y = true;
                driven_y = false;
                vel.y() = 0.0;
            } else if (rule.kind == BoundaryConditionKind::RollerY) {
                // Roller Y: Fix ux = 0, Free uy
                locked_x = true;
                driven_x = false;
                vel.x() = 0.0;
            } else if (rule.kind == BoundaryConditionKind::FixedClamped) {
                if (rule.lock_x) {
                    locked_x = true;
                    driven_x = false;
                    vel.x() = 0.0;
                }
                if (rule.lock_y) {
                    locked_y = true;
                    driven_y = false;
                    vel.y() = 0.0;
                }
            } else if (rule.kind == BoundaryConditionKind::PrescribedDisplacement) {
                // ux, uy are normalized by Ly so alpha = 0.01 means 1% engineering strain
                if (rule.lock_x) {
                    locked_x = true;
                    driven_x = true;
                    vel.x() = rule.ux * Ly;
                }
                if (rule.lock_y) {
                    locked_y = true;
                    driven_y = true;
                    vel.y() = rule.uy * Ly;
                }
            } else if (rule.kind == BoundaryConditionKind::AffineDeformation) {
                locked_x = true;
                locked_y = true;
                driven_x = true;
                driven_y = true;
                Eigen::Vector2d rel_corner(x - x_min, y - y_min);
                // Symmetric modes are centred only on closed contours around the specimen
                // (outer rim, hole rims); edge / custom-box rules keep the bottom-left corner origin.
                // Hole rims pivot about their own hole's centre; the outer rim about the specimen centre.
                Eigen::Vector2d pivot(x_min, y_min);
                if (affineRuleUsesCenterOrigin(rule)) {
                    pivot = (rule.region == BoundaryRegionType::HoleSurfaces)
                        ? cutoutPivot(findCutoutNear(x, y, strip))
                        : center;
                }
                Eigen::Vector2d rel_pivot(x - pivot.x(), y - pivot.y());
                if (rule.affine_mode == "simple_shear_y") {
                    vel = Eigen::Vector2d(0.0, rel_corner.x());
                } else if (rule.affine_mode == "uniaxial_y") {
                    vel = Eigen::Vector2d(0.0, rel_corner.y());
                } else if (rule.affine_mode == "uniaxial_x") {
                    vel = Eigen::Vector2d(rel_corner.x(), 0.0);
                } else if (rule.affine_mode == "pure_shear") {
                    vel = Eigen::Vector2d(rel_pivot.x(), -rel_pivot.y());
                } else if (rule.affine_mode == "biaxial" || rule.affine_mode == "radial_outward") {
                    vel = Eigen::Vector2d(rel_pivot.x(), rel_pivot.y());
                } else if (rule.affine_mode == "radial_inward" || rule.affine_mode == "biaxial_compression") {
                    vel = Eigen::Vector2d(-rel_pivot.x(), -rel_pivot.y());
                } else if (rule.affine_mode == "circumferential_shear" || rule.affine_mode == "torsion") {
                    vel = Eigen::Vector2d(-rel_pivot.y(), rel_pivot.x());
                    rotation_driven = true;
                    rotation_pivot = pivot;
                } else {
                    // default: simple_shear_x -> u_x = alpha * (y - y_min)
                    vel = Eigen::Vector2d(rel_corner.y(), 0.0);
                }
            }
        }

        node_bc_velocity[i] = vel;
        if (grip_of_node >= 0) {
            cached_node_grip[i] = grip_of_node;
            cached_force_grips[grip_of_node].nodes.push_back(static_cast<int>(i));
        }
        if (rotation_driven) {
            cached_node_rotation_driven[i] = 1;
            cached_node_rotation_pivot[i] = rotation_pivot;
        }
        if (locked_x && locked_y) {
            node_bc_role[i] = (driven_x || driven_y) ? 2 : 1;
        } else if (locked_x && !locked_y) {
            node_bc_role[i] = 3; // X-Constrained, Y-Free (Roller Y or partial X lock)
        } else if (!locked_x && locked_y) {
            node_bc_role[i] = 4; // Y-Constrained, X-Free (Roller X or partial Y lock)
        } else {
            node_bc_role[i] = 0; // Free in both X and Y
        }
    }

    // Dead loads with no kinematic constraint at all: pin a minimal isostatic anchor (3 DOFs) so the
    // rigid-body modes are removed. Decided here, before the mesh exists, so the DOF mappings are final.
    if (!cached_traction_rule_index.empty() && cached_force_grips.empty() &&
        std::all_of(node_bc_role.begin(), node_bc_role.end(), [](int r) { return r == 0; })) {
        Eigen::Vector2d centroid = Eigen::Vector2d::Zero();
        for (const auto& pt : carved_points) centroid += pt.coord;
        centroid /= static_cast<double>(n_pts);
        size_t anchor_a = 0;
        double best = std::numeric_limits<double>::max();
        for (size_t i = 0; i < n_pts; ++i) {
            const double d = (carved_points[i].coord - centroid).squaredNorm();
            if (d < best) { best = d; anchor_a = i; }
        }
        // Second atom on the same row, as far as possible horizontally: locking its u_y removes rotation
        size_t anchor_b = anchor_a;
        double best_dx = 0.0;
        const Eigen::Vector2d pa = carved_points[anchor_a].coord;
        for (size_t i = 0; i < n_pts; ++i) {
            const Eigen::Vector2d d = carved_points[i].coord - pa;
            if (std::abs(d.y()) < 0.5 * a0 && std::abs(d.x()) > best_dx) { best_dx = std::abs(d.x()); anchor_b = i; }
        }
        node_bc_role[anchor_a] = 1;
        if (anchor_b != anchor_a) node_bc_role[anchor_b] = 4;
        traction_anchor_active = true;
        std::cout << "[MacroscopicGeometry] Traction-only loading: isostatic anchor at node " << anchor_a
                  << " (ux = uy = 0) and node " << anchor_b << " (uy = 0, dx = " << best_dx << ")\n";
    }

    // Build interior_mapping and full_mapping (nodes with at least 1 free DOF are included)
    cached_node_bc_role = node_bc_role;
    has_partial_dof_constraints = false;
    interior_mapping.reserve(n_pts);
    full_mapping.reserve(n_pts);
    int solver_idx = 0;
    int count_free = 0, count_fixed = 0, count_driven = 0, count_roller_x = 0, count_roller_y = 0;

    for (size_t i = 0; i < n_pts; ++i) {
        int role = node_bc_role[i];
        if (role == 0 || role == 3 || role == 4) {
            interior_mapping.emplace_back(static_cast<int>(i), solver_idx);
            full_mapping.emplace_back(static_cast<int>(i), solver_idx);
            ++solver_idx;
            if (role == 0) {
                ++count_free;
            } else if (role == 3) {
                ++count_roller_y;
                has_partial_dof_constraints = true;
            } else if (role == 4) {
                ++count_roller_x;
                has_partial_dof_constraints = true;
            }
        } else {
            full_mapping.emplace_back(static_cast<int>(i), -1);
            if (role == 1) ++count_fixed;
            else ++count_driven;
        }
    }

    std::cout << "[MacroscopicGeometry] Boundary DOFs classified: "
              << count_free << " Free, "
              << count_fixed << " Fixed (ux=0, uy=0), "
              << count_driven << " Driven (ux, uy), "
              << count_roller_x << " Roller-X (Free ux, Constrained uy), "
              << count_roller_y << " Roller-Y (Constrained ux, Free uy)\n";

    if (hasTorsionDrivenNodes() && count_fixed + count_roller_x + count_roller_y == 0) {
        std::cerr << "[MacroscopicGeometry] WARNING: torsion (circumferential_shear) is active but no node is "
                     "Fixed or on a Roller. A free specimen only rotates rigidly, so the torque stays ~0. "
                     "Clamp an inner boundary (e.g. hole_surfaces fixed).\n";
    }
}

bool MacroscopicGeometryConfig::hasTractionRules() const {
    if (!enabled) return false;
    for (const auto& rule : boundary_rules) {
        if (rule.enabled && rule.kind == BoundaryConditionKind::PrescribedTraction) return true;
    }
    return false;
}

bool MacroscopicGeometryConfig::finalizeBoundaryTractions(
    const std::vector<Point2D>& reference_points,
    const std::vector<ElementTriangle2D>& elements) const {
    has_external_tractions = false;
    const size_t n_pts = reference_points.size();
    cached_node_ext_force_unit.assign(n_pts, Eigen::Vector2d::Zero());
    cached_node_ext_force_ref.resize(n_pts);
    for (size_t i = 0; i < n_pts; ++i) cached_node_ext_force_ref[i] = reference_points[i].coord;
    if (!enabled || (cached_traction_rule_index.empty() && cached_force_grips.empty())) return true;
    for (const auto& mem : cached_traction_membership) {
        if (mem.size() != n_pts) {
            std::cerr << "[Traction] ERROR: boundary DOFs were built for a different point set\n";
            return false;
        }
    }

    // Boundary edges of the reference mesh: edges used by exactly one triangle (keep the third vertex
    // to orient the outward normal). The solid area normalizes the conjugate strain.
    struct EdgeInfo { int count = 0; int a = -1; int b = -1; int third = -1; };
    std::unordered_map<uint64_t, EdgeInfo> edges;
    edges.reserve(elements.size() * 2);
    double solid_area = 0.0;
    for (const auto& el : elements) {
        int nn[3] = {el.getNodeIndex(0), el.getNodeIndex(1), el.getNodeIndex(2)};
        bool valid = true;
        for (int k = 0; k < 3; ++k) valid = valid && nn[k] >= 0 && static_cast<size_t>(nn[k]) < n_pts;
        if (!valid) continue;
        const Eigen::Vector2d e1 = reference_points[nn[1]].coord - reference_points[nn[0]].coord;
        const Eigen::Vector2d e2 = reference_points[nn[2]].coord - reference_points[nn[0]].coord;
        solid_area += 0.5 * std::abs(e1.x() * e2.y() - e1.y() * e2.x());
        for (int k = 0; k < 3; ++k) {
            const int a = nn[k], b = nn[(k + 1) % 3], c = nn[(k + 2) % 3];
            const uint64_t key = static_cast<uint64_t>(std::min(a, b)) * static_cast<uint64_t>(n_pts) +
                                 static_cast<uint64_t>(std::max(a, b));
            EdgeInfo& info = edges[key];
            info.count++;
            info.a = a; info.b = b; info.third = c;
        }
    }
    cached_reference_solid_area = solid_area;

    const Eigen::Vector2d center(x_min + 0.5 * Lx, y_min + 0.5 * Ly);
    int loaded_edges = 0;
    for (const auto& [key, info] : edges) {
        if (info.count != 1) continue;
        const Eigen::Vector2d pa = reference_points[info.a].coord;
        const Eigen::Vector2d pb = reference_points[info.b].coord;
        const Eigen::Vector2d mid = 0.5 * (pa + pb);
        const double len = (pb - pa).norm();
        if (len < 1e-12) continue;
        Eigen::Vector2d n_hat((pb - pa).y() / len, -(pb - pa).x() / len);
        if (n_hat.dot(reference_points[info.third].coord - mid) > 0.0) n_hat = -n_hat; // outward of the solid

        for (size_t m = 0; m < cached_traction_rule_index.size(); ++m) {
            const auto& mem = cached_traction_membership[m];
            if (!mem[info.a] || !mem[info.b]) continue;
            const BoundaryZoneRule& rule = boundary_rules[cached_traction_rule_index[m]];
            // Edge rules load only the faces that look that way: a thick strip also captures the
            // short side edges next to the corners, which must stay traction-free.
            Eigen::Vector2d face_dir = Eigen::Vector2d::Zero();
            switch (rule.region) {
            case BoundaryRegionType::TopEdge: face_dir = Eigen::Vector2d(0.0, 1.0); break;
            case BoundaryRegionType::BottomEdge: face_dir = Eigen::Vector2d(0.0, -1.0); break;
            case BoundaryRegionType::LeftEdge: face_dir = Eigen::Vector2d(-1.0, 0.0); break;
            case BoundaryRegionType::RightEdge: face_dir = Eigen::Vector2d(1.0, 0.0); break;
            default: break;
            }
            if (face_dir.squaredNorm() > 0.0 && n_hat.dot(face_dir) < 0.5) continue;
            const bool on_hole = (rule.region == BoundaryRegionType::HoleSurfaces);
            Eigen::Vector2d t0 = Eigen::Vector2d::Zero();
            if (rule.traction_mode == "radial_outward" || rule.traction_mode == "radial_inward") {
                // Radially away from the contour's pivot: +n on outer contours, -n on hole faces
                const double s = on_hole ? -1.0 : 1.0;
                t0 = (rule.traction_mode == "radial_outward" ? s : -s) * n_hat;
            } else if (rule.traction_mode == "circumferential_shear" || rule.traction_mode == "torsion") {
                const double strip = std::max(0.35, rule.thickness_a0) * a0;
                const Eigen::Vector2d pivot = on_hole ? cutoutPivot(findCutoutNear(mid.x(), mid.y(), strip)) : center;
                const Eigen::Vector2d r = mid - pivot;
                if (r.norm() > 1e-12) t0 = Eigen::Vector2d(-r.y(), r.x()) / r.norm();
            } else {
                t0 = Eigen::Vector2d(rule.tx, rule.ty);
            }
            const Eigen::Vector2d share = 0.5 * len * t0;
            cached_node_ext_force_unit[info.a] += share;
            cached_node_ext_force_unit[info.b] += share;
            ++loaded_edges;
        }
    }

    // Rigid force grips: F0 = (tx, ty) * L_face, with L_face the reference length of the grip's
    // boundary edges that face the region's direction (all boundary edges for non-edge regions)
    for (auto& grip : cached_force_grips) {
        const BoundaryZoneRule& rule = boundary_rules[grip.rule_index];
        Eigen::Vector2d face_dir = Eigen::Vector2d::Zero();
        switch (rule.region) {
        case BoundaryRegionType::TopEdge: face_dir = Eigen::Vector2d(0.0, 1.0); break;
        case BoundaryRegionType::BottomEdge: face_dir = Eigen::Vector2d(0.0, -1.0); break;
        case BoundaryRegionType::LeftEdge: face_dir = Eigen::Vector2d(-1.0, 0.0); break;
        case BoundaryRegionType::RightEdge: face_dir = Eigen::Vector2d(1.0, 0.0); break;
        default: break;
        }
        double face_length = 0.0;
        for (const auto& [key, info] : edges) {
            if (info.count != 1 || !isGripNode(info.a) || !isGripNode(info.b)) continue;
            if (cached_node_grip[info.a] != cached_node_grip[info.b]) continue;
            if (&cached_force_grips[cached_node_grip[info.a]] != &grip) continue;
            const Eigen::Vector2d pa = reference_points[info.a].coord, pb = reference_points[info.b].coord;
            const double len = (pb - pa).norm();
            if (len < 1e-12) continue;
            Eigen::Vector2d n_hat((pb - pa).y() / len, -(pb - pa).x() / len);
            if (n_hat.dot(reference_points[info.third].coord - 0.5 * (pa + pb)) > 0.0) n_hat = -n_hat;
            if (face_dir.squaredNorm() > 0.0 && n_hat.dot(face_dir) < 0.5) continue;
            face_length += len;
        }
        grip.face_length = face_length;
        grip.F0 = face_length * Eigen::Vector2d(rule.tx, rule.ty);
        if (!grip.free_x) grip.F0.x() = 0.0;
        if (!grip.free_y) grip.F0.y() = 0.0;
        std::cout << "[Force Grip] rule '" << rule.id << "': " << grip.nodes.size() << " gripped atoms, face length "
                  << face_length << ", F0 = (" << grip.F0.x() << ", " << grip.F0.y() << ") per unit alpha, U free in "
                  << (grip.free_x ? "x" : "") << (grip.free_y ? "y" : "")
                  << ((grip.free_x || grip.free_y) ? "" : "none") << "\n";
        if (grip.nodes.empty()) {
            std::cerr << "[Force Grip] WARNING: grip rule '" << rule.id << "' matched no atom.\n";
        }
    }
    if (!cached_force_grips.empty()) {
        // A force grip only transmits the load: some other kinematic boundary must hold the specimen
        bool other_support = false;
        for (size_t i = 0; i < n_pts && i < cached_node_bc_role.size(); ++i) {
            if (cached_node_bc_role[i] != 0 && !isGripNode(i)) { other_support = true; break; }
        }
        if (!other_support) {
            std::cerr << "[Force Grip] ERROR: a rigid force grip needs another Fixed, Prescribed, Affine or "
                         "Roller boundary to react the load (e.g. clamp the opposite edge).\n";
            return false;
        }
    }

    double force_scale = 0.0;
    Eigen::Vector2d resultant = Eigen::Vector2d::Zero();
    Eigen::Vector2d centroid = Eigen::Vector2d::Zero();
    for (const auto& pt : reference_points) centroid += pt.coord;
    if (n_pts > 0) centroid /= static_cast<double>(n_pts);
    double moment = 0.0;
    int loaded_nodes = 0;
    for (size_t i = 0; i < n_pts; ++i) {
        const Eigen::Vector2d& f = cached_node_ext_force_unit[i];
        if (f.squaredNorm() == 0.0) continue;
        ++loaded_nodes;
        force_scale += f.norm();
        resultant += f;
        const Eigen::Vector2d r = reference_points[i].coord - centroid;
        moment += r.x() * f.y() - r.y() * f.x();
    }
    has_external_tractions = (force_scale > 0.0);
    if (cached_traction_rule_index.empty()) return true; // grips only

    std::cout << "[Traction] " << cached_traction_rule_index.size() << " traction rule(s): " << loaded_edges
              << " boundary edge share(s) on " << loaded_nodes << " contour atoms | sum|f0| = " << force_scale
              << " | R = (" << resultant.x() << ", " << resultant.y() << "), M = " << moment
              << " (per unit alpha) | A_solid = " << solid_area << "\n";
    if (!has_external_tractions) {
        std::cerr << "[Traction] WARNING: traction rule(s) matched no boundary edge; no load is applied.\n";
        return true;
    }
    if (traction_anchor_active) {
        const double l_char = std::max(Lx, Ly);
        const bool unbalanced = resultant.norm() > 1e-10 * force_scale ||
                                std::abs(moment) > 1e-10 * force_scale * l_char;
        if (unbalanced) {
            std::cerr << "[Traction] ERROR: the dead load is not self-equilibrated (R = (" << resultant.x() << ", "
                      << resultant.y() << "), M = " << moment << " per unit alpha) and no Fixed, Prescribed, "
                      << "Affine or Roller boundary holds the specimen. Add a support or balance the load.\n";
            return false;
        }
    } else if (cached_force_grips.empty()) { // grip loads act through U, not through the nodal forces
        // The supports can leave rigid-body modes free (rollers, a single pinned atom): a translation along x (y) when no
        // node locks u_x (u_y), and a rotation about (x0, y0) when every x-locked node lies on y = y0 and every
        // y-locked node on x = x0. The load must have no component along a free mode.
        const double tol_pos = 1e-3 * a0;
        bool lock_x = false, lock_y = false, rot_free = true, have_y0 = false, have_x0 = false;
        double x0 = centroid.x(), y0 = centroid.y();
        for (size_t i = 0; i < n_pts && i < cached_node_bc_role.size(); ++i) {
            const int role = cached_node_bc_role[i];
            if (role == 0) continue;
            bool hold_x = (role == 1 || role == 2 || role == 3), hold_y = (role == 1 || role == 2 || role == 4);
            const Eigen::Vector2d& p = reference_points[i].coord;
            if (hold_x) { // u_x held
                lock_x = true;
                if (!have_y0) { y0 = p.y(); have_y0 = true; }
                else if (std::abs(p.y() - y0) > tol_pos) rot_free = false;
            }
            if (hold_y) { // u_y held
                lock_y = true;
                if (!have_x0) { x0 = p.x(); have_x0 = true; }
                else if (std::abs(p.x() - x0) > tol_pos) rot_free = false;
            }
        }
        double moment_pivot = 0.0;
        for (size_t i = 0; i < n_pts; ++i) {
            const Eigen::Vector2d& f = cached_node_ext_force_unit[i];
            const Eigen::Vector2d r = reference_points[i].coord - Eigen::Vector2d(x0, y0);
            moment_pivot += r.x() * f.y() - r.y() * f.x();
        }
        const double l_char = std::max(Lx, Ly);
        const bool free_x = !lock_x && std::abs(resultant.x()) > 1e-10 * force_scale;
        const bool free_y = !lock_y && std::abs(resultant.y()) > 1e-10 * force_scale;
        const bool free_rot = rot_free && std::abs(moment_pivot) > 1e-10 * force_scale * l_char;
        if (free_x || free_y || free_rot) {
            std::cerr << "[Traction] ERROR: the supports leave a rigid-body "
                      << (free_rot ? "rotation" : "translation") << " free and the dead load acts along it (R = ("
                      << resultant.x() << ", " << resultant.y() << "), M = " << moment_pivot
                      << " per unit alpha). Add a support that holds it or balance the load.\n";
            return false;
        }
    }
    return true;
}

Eigen::Vector2d MacroscopicGeometryConfig::gripDisplacement(int g, const std::vector<Point2D>& points) const {
    if (g < 0 || g >= static_cast<int>(cached_force_grips.size()) || cached_force_grips[g].nodes.empty()) {
        return Eigen::Vector2d::Zero();
    }
    const size_t n0 = static_cast<size_t>(cached_force_grips[g].nodes.front());
    if (n0 >= points.size() || n0 >= cached_node_ext_force_ref.size()) return Eigen::Vector2d::Zero();
    return points[n0].coord - cached_node_ext_force_ref[n0];
}

void MacroscopicGeometryConfig::setGripPositions(int g, const Eigen::Vector2d& U, std::vector<Point2D>& points) const {
    if (g < 0 || g >= static_cast<int>(cached_force_grips.size())) return;
    for (int n : cached_force_grips[g].nodes) {
        if (n >= 0 && static_cast<size_t>(n) < points.size() && static_cast<size_t>(n) < cached_node_ext_force_ref.size()) {
            points[n].coord = cached_node_ext_force_ref[n] + U;
        }
    }
}

double MacroscopicGeometryConfig::computeExternalForcePotential(
    const std::vector<Point2D>& points, double alpha) const {
    if (!hasForceControl() || alpha == 0.0) return 0.0;
    double energy = 0.0;
    for (int g = 0; g < static_cast<int>(cached_force_grips.size()); ++g) {
        energy -= alpha * cached_force_grips[g].F0.dot(gripDisplacement(g, points));
    }
    if (!hasExternalTractionForces()) return energy;
    const size_t n = std::min({points.size(), cached_node_ext_force_unit.size(),
                               cached_node_ext_force_ref.size(), cached_node_bc_role.size()});
    for (size_t i = 0; i < n; ++i) {
        const int role = cached_node_bc_role[i];
        if (role != 0 && role != 3 && role != 4) continue; // same node set as the optimizer's DOFs
        const Eigen::Vector2d& f0 = cached_node_ext_force_unit[i];
        if (f0.squaredNorm() == 0.0) continue;
        energy -= alpha * f0.dot(points[i].coord - cached_node_ext_force_ref[i]);
    }
    return energy;
}

double MacroscopicGeometryConfig::accumulateExternalForcePotential(
    const alglib::real_1d_array& x,
    const std::vector<std::pair<int, int>>& interior_mapping,
    alglib::real_1d_array& grad,
    double alpha) const {
    if (!hasExternalTractionForces() || alpha == 0.0) return 0.0;
    const int n_vars = static_cast<int>(interior_mapping.size());
    const int n_forces = static_cast<int>(std::min(cached_node_ext_force_unit.size(), cached_node_ext_force_ref.size()));
    double energy = 0.0;
    for (int i = 0; i < n_vars; ++i) {
        const auto& [point_idx, dof_idx] = interior_mapping[i];
        if (point_idx < 0 || point_idx >= n_forces || dof_idx < 0 || dof_idx >= n_vars) continue;
        const Eigen::Vector2d& f0 = cached_node_ext_force_unit[point_idx];
        if (f0.squaredNorm() == 0.0) continue;
        const Eigen::Vector2d f = alpha * f0;
        const Eigen::Vector2d& X = cached_node_ext_force_ref[point_idx];
        energy -= f.x() * (x[dof_idx] - X.x()) + f.y() * (x[n_vars + dof_idx] - X.y());
        grad[dof_idx] -= f.x();
        grad[n_vars + dof_idx] -= f.y();
    }
    return energy;
}

double MacroscopicGeometryConfig::computeTractionConjugateStrain(const std::vector<Point2D>& points) const {
    if (!hasForceControl() || cached_reference_solid_area <= 0.0) return 0.0;
    const size_t n = std::min({points.size(), cached_node_ext_force_unit.size(),
                               cached_node_ext_force_ref.size(), cached_node_bc_role.size()});
    double work = 0.0;
    for (int g = 0; g < static_cast<int>(cached_force_grips.size()); ++g) {
        work += cached_force_grips[g].F0.dot(gripDisplacement(g, points));
    }
    for (size_t i = 0; i < n; ++i) {
        const int role = cached_node_bc_role[i];
        if (role != 0 && role != 3 && role != 4) continue;
        work += cached_node_ext_force_unit[i].dot(points[i].coord - cached_node_ext_force_ref[i]);
    }
    return work / cached_reference_solid_area;
}

void MacroscopicGeometryConfig::setRunawayReference(const std::vector<Point2D>& points, double limit) const {
    runaway_limit = limit;
    runaway_reference_points.clear();
    if (limit <= 0.0) return;
    runaway_reference_points.reserve(points.size());
    for (const auto& pt : points) runaway_reference_points.push_back(pt.coord);
}

double MacroscopicGeometryConfig::maxDisplacementSinceStepStart(const std::vector<Point2D>& points) const {
    double max_jump = 0.0;
    const size_t n = std::min(points.size(), runaway_reference_points.size());
    for (size_t i = 0; i < n; ++i) {
        max_jump = std::max(max_jump, (points[i].coord - runaway_reference_points[i]).norm());
    }
    return max_jump;
}

bool MacroscopicGeometryConfig::hasTorsionDrivenNodes() const {
    for (char flag : cached_node_rotation_driven) {
        if (flag) return true;
    }
    return false;
}

double MacroscopicGeometryConfig::computeRimTorque(
    const UserData* userData,
    const std::vector<Point2D>& current_points) const {
    if (!userData || !hasTorsionDrivenNodes()) return 0.0;
    const double normalisation = userData->calculator.getUnitCellArea();
    const size_t n_flags = cached_node_rotation_driven.size();
    double torque = 0.0;

    for (size_t elem_idx : userData->active_elements) {
        const auto& elem = userData->elements[elem_idx];
        int idx[3] = {elem.getNodeIndex(0), elem.getNodeIndex(1), elem.getNodeIndex(2)};
        bool touches_rim = false;
        for (int k = 0; k < 3; ++k) {
            if (idx[k] >= 0 && static_cast<size_t>(idx[k]) < n_flags && cached_node_rotation_driven[idx[k]]) {
                touches_rim = true;
            }
        }
        if (!touches_rim) continue;

        const Eigen::Matrix2d& F = elem.getDeformationGradient();
        Eigen::Matrix2d C = elem.getMetricTensor();
        auto result = reduce_metric(C);
        Eigen::Matrix2d dE_dC = userData->calculator.calculate_derivative(
            result.C_reduced, userData->derivative_function) / normalisation;
        Eigen::Matrix2d P_stress = 2.0 * F * result.m_matrix * dE_dC * result.m_matrix.transpose();
        auto nodal_forces = elem.calculate_nodal_forces(P_stress);

        for (int k = 0; k < 3; ++k) {
            const int n = idx[k];
            if (n < 0 || static_cast<size_t>(n) >= n_flags || !cached_node_rotation_driven[n] ||
                static_cast<size_t>(n) >= current_points.size()) continue;
            // nodal_forces = dE/dx: the load the grips must apply to hold the rim at R(alpha)
            const Eigen::Vector2d r = current_points[n].coord - cached_node_rotation_pivot[n];
            torque += r.x() * nodal_forces[k].y() - r.y() * nodal_forces[k].x();
        }
    }
    return torque;
}

std::vector<BoundaryZoneRule> MacroscopicGeometryConfig::buildStage1RelaxationRules() const {
    if (!stage1_boundary_rules.empty()) {
        return stage1_boundary_rules;
    }
    if (microstructure.mode == MicrostructureMode::Bicrystal) {
        std::vector<BoundaryZoneRule> rules;
        const double phi_rad = microstructure.gb_angle_deg * M_PI / 180.0;
        const bool is_vertical_gb = std::abs(std::cos(phi_rad)) >= std::abs(std::sin(phi_rad));

        if (is_vertical_gb) {
            // Vertical-dominant GB:
            // Top & Bottom on RollerX (lock uy = 0, free ux) -> keeps GB terminations flat, allows normal breathing
            BoundaryZoneRule r_top;
            r_top.id = "stage1_gb_top_roller";
            r_top.region = BoundaryRegionType::TopEdge;
            r_top.kind = BoundaryConditionKind::RollerX;
            r_top.thickness_a0 = 1.5;
            r_top.lock_x = false;
            r_top.lock_y = true;
            rules.push_back(r_top);

            BoundaryZoneRule r_bot;
            r_bot.id = "stage1_gb_bot_roller";
            r_bot.region = BoundaryRegionType::BottomEdge;
            r_bot.kind = BoundaryConditionKind::RollerX;
            r_bot.thickness_a0 = 1.5;
            r_bot.lock_x = false;
            r_bot.lock_y = true;
            rules.push_back(r_bot);

            // Left Edge on RollerY (lock ux = 0, free uy) -> anchors Grain 1 horizontally while Right Edge is FreeSurface
            // so Grain 2 can freely adjust both normal expansion Δu_⊥ and tangential shift Δu_∥
            BoundaryZoneRule r_left;
            r_left.id = "stage1_gb_left_pin";
            r_left.region = BoundaryRegionType::LeftEdge;
            r_left.kind = BoundaryConditionKind::RollerY;
            r_left.thickness_a0 = 1.5;
            r_left.lock_x = true;
            r_left.lock_y = false;
            rules.push_back(r_left);
        } else {
            // Horizontal-dominant GB:
            // Left & Right on RollerY (lock ux = 0, free uy) -> keeps GB terminations flat
            BoundaryZoneRule r_left;
            r_left.id = "stage1_gb_left_roller";
            r_left.region = BoundaryRegionType::LeftEdge;
            r_left.kind = BoundaryConditionKind::RollerY;
            r_left.thickness_a0 = 1.5;
            r_left.lock_x = true;
            r_left.lock_y = false;
            rules.push_back(r_left);

            BoundaryZoneRule r_right;
            r_right.id = "stage1_gb_right_roller";
            r_right.region = BoundaryRegionType::RightEdge;
            r_right.kind = BoundaryConditionKind::RollerY;
            r_right.thickness_a0 = 1.5;
            r_right.lock_x = true;
            r_right.lock_y = false;
            rules.push_back(r_right);

            // Bottom Edge on RollerX (lock uy = 0, free ux) -> anchors bottom grain vertically while Top Edge breathes
            BoundaryZoneRule r_bot;
            r_bot.id = "stage1_gb_bot_pin";
            r_bot.region = BoundaryRegionType::BottomEdge;
            r_bot.kind = BoundaryConditionKind::RollerX;
            r_bot.thickness_a0 = 1.5;
            r_bot.lock_x = false;
            r_bot.lock_y = true;
            rules.push_back(r_bot);
        }
        return rules;
    }

    // For Voronoi Polycrystal or general cutouts: keep non-driven rules or convert driven rules to rollers at alpha = 0
    return boundary_rules;
}

std::vector<BoundaryZoneRule> MacroscopicGeometryConfig::buildStage2LoadingRules(
    const std::string& fallback_affine_mode) const {

    const std::string& preset = microstructure.loading_bc_preset;
    const std::string& ltype = microstructure.stage2_load_type;
    if (preset.empty() || preset == "keep_studio") {
        return boundary_rules;
    }

    // Map stage2_load_type to affine mode and grip displacement vector
    std::string aff_mode = fallback_affine_mode;
    if (ltype == "shear_parallel_gb") aff_mode = "simple_shear_y";
    else if (ltype == "shear_horizontal") aff_mode = "simple_shear_x";
    else if (ltype == "tension_normal_x" || ltype == "compression_normal_x") aff_mode = "uniaxial_x";
    else if (ltype == "tension_vertical_y" || ltype == "compression_vertical_y") aff_mode = "uniaxial_y";
    else {
        for (const auto& r : boundary_rules) {
            if (r.enabled && r.kind == BoundaryConditionKind::AffineDeformation && !r.affine_mode.empty()) {
                aff_mode = r.affine_mode;
                break;
            }
        }
    }

    if (preset == "affine_from_relaxed") {
        std::vector<BoundaryZoneRule> rules;
        BoundaryZoneRule r_aff;
        r_aff.id = "stage2_affine_outer";
        r_aff.region = BoundaryRegionType::OuterPerimeter;
        r_aff.kind = BoundaryConditionKind::AffineDeformation;
        r_aff.thickness_a0 = 2.0;
        r_aff.lock_x = true;
        r_aff.lock_y = true;
        r_aff.affine_mode = aff_mode;
        rules.push_back(r_aff);
        return rules;
    }

    if (preset == "x_grips_y_roller") {
        // Left & Right edges driven (shear or normal tension/compression across vertical GB), Top & Bottom on RollerX
        std::vector<BoundaryZoneRule> rules;
        double drive_ux = 0.0;
        double drive_uy = 1.0;
        if (ltype == "tension_normal_x") {
            drive_ux = 1.0; drive_uy = 0.0;
        } else if (ltype == "compression_normal_x") {
            drive_ux = -1.0; drive_uy = 0.0;
        } else if (ltype == "shear_parallel_gb") {
            drive_ux = 0.0; drive_uy = 1.0;
        } else if (primary_stress_component == "sigma_xx") {
            drive_ux = 1.0; drive_uy = 0.0;
        } else {
            drive_ux = 0.0; drive_uy = 1.0;
        }

        BoundaryZoneRule r_left;
        r_left.id = "stage2_left_grip";
        r_left.region = BoundaryRegionType::LeftEdge;
        r_left.kind = BoundaryConditionKind::FixedClamped;
        r_left.thickness_a0 = 2.0;
        r_left.lock_x = true;
        r_left.lock_y = true;
        rules.push_back(r_left);

        BoundaryZoneRule r_right;
        r_right.id = "stage2_right_grip";
        r_right.region = BoundaryRegionType::RightEdge;
        r_right.kind = BoundaryConditionKind::PrescribedDisplacement;
        r_right.thickness_a0 = 2.0;
        r_right.ux = drive_ux;
        r_right.uy = drive_uy;
        r_right.lock_x = true;
        r_right.lock_y = true;
        rules.push_back(r_right);

        BoundaryZoneRule r_top;
        r_top.id = "stage2_top_roller";
        r_top.region = BoundaryRegionType::TopEdge;
        r_top.kind = BoundaryConditionKind::RollerX;
        r_top.thickness_a0 = 1.5;
        r_top.lock_x = false;
        r_top.lock_y = true;
        rules.push_back(r_top);

        BoundaryZoneRule r_bot;
        r_bot.id = "stage2_bot_roller";
        r_bot.region = BoundaryRegionType::BottomEdge;
        r_bot.kind = BoundaryConditionKind::RollerX;
        r_bot.thickness_a0 = 1.5;
        r_bot.lock_x = false;
        r_bot.lock_y = true;
        rules.push_back(r_bot);
        return rules;
    }

    if (preset == "y_grips_x_roller") {
        // Top & Bottom edges driven, Left & Right on RollerY
        std::vector<BoundaryZoneRule> rules;
        double drive_ux = 1.0;
        double drive_uy = 0.0;
        if (ltype == "tension_vertical_y") {
            drive_ux = 0.0; drive_uy = 1.0;
        } else if (ltype == "compression_vertical_y") {
            drive_ux = 0.0; drive_uy = -1.0;
        } else if (ltype == "shear_horizontal") {
            drive_ux = 1.0; drive_uy = 0.0;
        } else if (primary_stress_component == "sigma_yy") {
            drive_ux = 0.0; drive_uy = 1.0;
        } else {
            drive_ux = 1.0; drive_uy = 0.0;
        }

        BoundaryZoneRule r_bot;
        r_bot.id = "stage2_bot_grip";
        r_bot.region = BoundaryRegionType::BottomEdge;
        r_bot.kind = BoundaryConditionKind::FixedClamped;
        r_bot.thickness_a0 = 2.0;
        r_bot.lock_x = true;
        r_bot.lock_y = true;
        rules.push_back(r_bot);

        BoundaryZoneRule r_top;
        r_top.id = "stage2_top_grip";
        r_top.region = BoundaryRegionType::TopEdge;
        r_top.kind = BoundaryConditionKind::PrescribedDisplacement;
        r_top.thickness_a0 = 2.0;
        r_top.ux = drive_ux;
        r_top.uy = drive_uy;
        r_top.lock_x = true;
        r_top.lock_y = true;
        rules.push_back(r_top);

        BoundaryZoneRule r_left;
        r_left.id = "stage2_left_roller";
        r_left.region = BoundaryRegionType::LeftEdge;
        r_left.kind = BoundaryConditionKind::RollerY;
        r_left.thickness_a0 = 1.5;
        r_left.lock_x = true;
        r_left.lock_y = false;
        rules.push_back(r_left);

        BoundaryZoneRule r_right;
        r_right.id = "stage2_right_roller";
        r_right.region = BoundaryRegionType::RightEdge;
        r_right.kind = BoundaryConditionKind::RollerY;
        r_right.thickness_a0 = 1.5;
        r_right.lock_x = true;
        r_right.lock_y = false;
        rules.push_back(r_right);
        return rules;
    }

    return boundary_rules;
}

Eigen::Matrix2d MacroscopicGeometryConfig::computeMacroscopicF(double alpha) const {
    Eigen::Matrix2d F = Eigen::Matrix2d::Identity();
    if (!enabled) {
        F(0, 1) = alpha;
        return F;
    }

    // Infer dominant loading mode from active driven boundary rules
    for (const auto& rule : boundary_rules) {
        if (!rule.enabled) continue;
        if (rule.kind == BoundaryConditionKind::AffineDeformation) {
            if (rule.affine_mode == "uniaxial_y") {
                F(1, 1) = 1.0 + alpha;
            } else if (rule.affine_mode == "uniaxial_x") {
                F(0, 0) = 1.0 + alpha;
            } else if (rule.affine_mode == "pure_shear") {
                F(0, 0) = 1.0 + alpha;
                F(1, 1) = 1.0 - alpha;
            } else if (rule.affine_mode == "biaxial" || rule.affine_mode == "radial_outward") {
                F(0, 0) = 1.0 + alpha;
                F(1, 1) = 1.0 + alpha;
            } else if (rule.affine_mode == "radial_inward" || rule.affine_mode == "biaxial_compression") {
                F(0, 0) = 1.0 - alpha;
                F(1, 1) = 1.0 - alpha;
            } else if (rule.affine_mode == "circumferential_shear" || rule.affine_mode == "torsion") {
                // Exact rotation R(alpha), alpha in radians (no spurious sqrt(1 + alpha^2) stretch)
                F(0, 0) = std::cos(alpha);
                F(0, 1) = -std::sin(alpha);
                F(1, 0) = std::sin(alpha);
                F(1, 1) = std::cos(alpha);
            } else if (rule.affine_mode == "simple_shear_y") {
                F(1, 0) = alpha;
            } else {
                F(0, 1) = alpha;
            }
            return F;
        } else if (rule.kind == BoundaryConditionKind::PrescribedDisplacement) {
            if (rule.region == BoundaryRegionType::TopEdge || rule.region == BoundaryRegionType::BottomEdge) {
                double sign = (rule.region == BoundaryRegionType::TopEdge) ? 1.0 : -1.0;
                F(0, 1) += sign * alpha * rule.ux;
                F(1, 1) += sign * alpha * rule.uy;
            } else if (rule.region == BoundaryRegionType::RightEdge || rule.region == BoundaryRegionType::LeftEdge) {
                double sign = (rule.region == BoundaryRegionType::RightEdge) ? 1.0 : -1.0;
                double scale = (Lx > 1e-8) ? (Ly / Lx) : 1.0;
                F(0, 0) += sign * alpha * rule.ux * scale;
                F(1, 0) += sign * alpha * rule.uy * scale;
            }
        }
    }
    return F;
}

void MacroscopicGeometryConfig::applyLoadStep(
    std::vector<Point2D>& current_points,
    const std::vector<Point2D>& points_ref,
    const std::vector<int>& node_bc_role,
    const std::vector<Eigen::Vector2d>& node_bc_velocity,
    double alpha,
    double dalpha,
    bool is_first_step,
    unsigned int seed,
    double noise_level) const {

    const size_t n_pts = current_points.size();
    if (n_pts == 0 || points_ref.size() != n_pts) return;

    Eigen::Matrix2d F_macro_alpha = computeMacroscopicF(alpha);
    // Linear modes: dF = F(dalpha) - I. Torsion's R(alpha) is non-linear, so use the exact
    // increment R(alpha) - R(alpha - dalpha) there.
    Eigen::Matrix2d dF_macro = hasTorsionDrivenNodes()
        ? Eigen::Matrix2d(F_macro_alpha - computeMacroscopicF(alpha - dalpha))
        : Eigen::Matrix2d(computeMacroscopicF(dalpha) - Eigen::Matrix2d::Identity());
    bool use_center_origin = false;
    for (const auto& rule : boundary_rules) {
        if (!rule.enabled) continue;
        if (rule.kind == BoundaryConditionKind::AffineDeformation && affineRuleUsesCenterOrigin(rule)) {
            use_center_origin = true;
            break;
        }
    }
    Eigen::Vector2d origin = use_center_origin
        ? Eigen::Vector2d(x_min + 0.5 * Lx, y_min + 0.5 * Ly)
        : Eigen::Vector2d(x_min, y_min);

    if (is_first_step) {
        if (std::abs(alpha) > 1e-12) {
            std::cout << std::defaultfloat << std::setprecision(6)
                      << "  📐 [Affine Pre-Strain] Deforming full crystal (interior + boundaries) to α_start = "
                      << alpha << " | F_macro = [[" << F_macro_alpha(0,0) << ", " << F_macro_alpha(0,1)
                      << "], [" << F_macro_alpha(1,0) << ", " << F_macro_alpha(1,1) << "]]" << std::endl;
        }
        if (noise_level > 0.0) {
            std::cout << std::scientific << std::setprecision(2)
                      << "  🎲 [Loading Stage Noise] Seeding positional thermal noise (σ = "
                      << noise_level << ", seed = " << seed << ") at loading onset"
                      << std::defaultfloat << std::endl;
        }
    }

    std::mt19937 gen(seed);
    std::normal_distribution<double> noise_dist(0.0, std::max(0.0, noise_level));

    for (size_t j = 0; j < n_pts; ++j) {
        Eigen::Vector2d noise_vec = (is_first_step && noise_level > 0.0)
            ? Eigen::Vector2d(noise_dist(gen), noise_dist(gen))
            : Eigen::Vector2d::Zero();

        // At the first loading step (`is_first_step`), measure `rel` from the Stage 1 relaxed
        // coordinates `current_points[j].coord` so that any Stage 1 grain-boundary or defect
        // relaxation is preserved and affinely mapped by F_macro(alpha), even if
        // `reanchor_ref_after_relax` is false.
        const Eigen::Vector2d base_coord = is_first_step ? current_points[j].coord : points_ref[j].coord;

        if (node_bc_role[j] == 1) {
            // Fixed / Clamped node: stays strictly at reference position X_0.
            // When an indenter is active, only true clamped specimen boundary nodes
            // (e.g. bottom clamped base) must stay at X_0.
            // Upper solid nodes that were in dynamic indenter contact must NEVER be reset to X_0!
            if (isIndenterActive() && points_ref[j].coord.y() >= y_min + 0.4 * Ly) {
                // Keep current deformed position from prior steps
            } else {
                current_points[j].coord = points_ref[j].coord;
            }
        } else if (node_bc_role[j] == 2 && isGripNode(j)) {
            // Rigid force grip: position X + U is a solver unknown; keep the current one as the guess
        } else if (node_bc_role[j] == 2) {
            // Driven / Prescribed node: exact boundary trajectory X_0 + alpha * v_j
            // (torsion rims: exact rotation X_c + R(alpha) (X_0 - X_c) about their pivot)
            if (j < cached_node_rotation_driven.size() && cached_node_rotation_driven[j]) {
                const Eigen::Vector2d& pc = cached_node_rotation_pivot[j];
                const double c = std::cos(alpha), s = std::sin(alpha);
                const Eigen::Vector2d d = points_ref[j].coord - pc;
                current_points[j].coord = pc + Eigen::Vector2d(c * d.x() - s * d.y(), s * d.x() + c * d.y());
            } else {
                current_points[j].coord = points_ref[j].coord + alpha * node_bc_velocity[j];
            }
        } else if (node_bc_role[j] == 3) {
            // Roller-Y / X-Constrained, Y-Free: lock x(alpha) = X_0.x + alpha * vx, free predictor on y
            Eigen::Vector2d rel = base_coord - origin;
            current_points[j].coord.x() = points_ref[j].coord.x() + alpha * node_bc_velocity[j].x();
            if (is_first_step) {
                current_points[j].coord.y() = (origin + F_macro_alpha * rel + noise_vec).y();
            } else {
                current_points[j].coord.y() += (dF_macro * rel).y();
            }
        } else if (node_bc_role[j] == 4) {
            // Roller-X / Y-Constrained, X-Free: lock y(alpha) = X_0.y + alpha * vy, free predictor on x
            Eigen::Vector2d rel = base_coord - origin;
            current_points[j].coord.y() = points_ref[j].coord.y() + alpha * node_bc_velocity[j].y();
            if (is_first_step) {
                current_points[j].coord.x() = (origin + F_macro_alpha * rel + noise_vec).x();
            } else {
                current_points[j].coord.x() += (dF_macro * rel).x();
            }
        } else {
            // Free node (interior or traction-free hole/outer surface):
            // Warm-start predictor using smooth macroscopic affine map at is_first_step, or increment thereafter
            Eigen::Vector2d rel = base_coord - origin;
            if (is_first_step) {
                current_points[j].coord = origin + F_macro_alpha * rel + noise_vec;
            } else {
                current_points[j].coord += dF_macro * rel;
            }
        }
    }
}

double MacroscopicGeometryConfig::extractPrimaryStress(
    const Eigen::Matrix2d& stress_tensor,
    const std::vector<Point2D>* current_points,
    double alpha,
    const UserData* userData,
    const std::vector<int>* contact_atoms) const {
    if (!enabled) return stress_tensor(0, 1);

    if (primary_stress_component == "indentation_load" || isIndenterActive()) {
        if (indenter.contact_mode == IndenterContactMode::KinematicProjection) {
            if (userData != nullptr && contact_atoms != nullptr) {
                return computeKinematicContactForce(userData, *contact_atoms);
            }
        } else if (current_points != nullptr) {
            return computeIndenterVerticalForce(*current_points, alpha);
        }
    }

    const bool torque_available = (userData != nullptr && current_points != nullptr && hasTorsionDrivenNodes());
    if (primary_stress_component == "torque") {
        return torque_available ? computeRimTorque(userData, *current_points) : stress_tensor(0, 1);
    }
    if (primary_stress_component == "sigma_yy") return stress_tensor(1, 1);
    if (primary_stress_component == "sigma_xx") return stress_tensor(0, 0);
    if (primary_stress_component == "sigma_xy") return stress_tensor(0, 1);
    if (primary_stress_component == "hydrostatic") return 0.5 * (stress_tensor(0, 0) + stress_tensor(1, 1));
    if (primary_stress_component == "von_mises") {
        double sxx = stress_tensor(0, 0), syy = stress_tensor(1, 1), sxy = stress_tensor(0, 1);
        return std::sqrt(std::max(0.0, sxx * sxx - sxx * syy + syy * syy + 3.0 * sxy * sxy));
    }

    // "auto": inspect active driven rules to pick the conjugate stress component
    for (const auto& rule : boundary_rules) {
        if (!rule.enabled) continue;
        if (rule.kind == BoundaryConditionKind::RigidForceGrip) {
            const bool horizontal_face = (rule.region == BoundaryRegionType::TopEdge ||
                                          rule.region == BoundaryRegionType::BottomEdge);
            const bool vertical_face = (rule.region == BoundaryRegionType::LeftEdge ||
                                        rule.region == BoundaryRegionType::RightEdge);
            if (horizontal_face && std::abs(rule.ty) > std::abs(rule.tx)) return stress_tensor(1, 1);
            if (vertical_face && std::abs(rule.tx) > std::abs(rule.ty)) return stress_tensor(0, 0);
            return stress_tensor(0, 1);
        }
        if (rule.kind == BoundaryConditionKind::PrescribedTraction) {
            if (rule.traction_mode == "radial_outward" || rule.traction_mode == "radial_inward") {
                return 0.5 * (stress_tensor(0, 0) + stress_tensor(1, 1));
            }
            if (rule.traction_mode == "cartesian") {
                const bool horizontal_face = (rule.region == BoundaryRegionType::TopEdge ||
                                              rule.region == BoundaryRegionType::BottomEdge);
                const bool vertical_face = (rule.region == BoundaryRegionType::LeftEdge ||
                                            rule.region == BoundaryRegionType::RightEdge);
                if (horizontal_face && std::abs(rule.ty) > std::abs(rule.tx)) return stress_tensor(1, 1);
                if (vertical_face && std::abs(rule.tx) > std::abs(rule.ty)) return stress_tensor(0, 0);
            }
            return stress_tensor(0, 1);
        }
        if (rule.kind == BoundaryConditionKind::PrescribedDisplacement) {
            if ((rule.region == BoundaryRegionType::TopEdge || rule.region == BoundaryRegionType::BottomEdge) &&
                std::abs(rule.uy) > std::abs(rule.ux)) {
                return stress_tensor(1, 1);
            }
            if ((rule.region == BoundaryRegionType::RightEdge || rule.region == BoundaryRegionType::LeftEdge) &&
                std::abs(rule.ux) > std::abs(rule.uy)) {
                return stress_tensor(0, 0);
            }
            return stress_tensor(0, 1);
        } else if (rule.kind == BoundaryConditionKind::AffineDeformation) {
            if (rule.affine_mode == "uniaxial_y") return stress_tensor(1, 1);
            if (rule.affine_mode == "uniaxial_x") return stress_tensor(0, 0);
            if (rule.affine_mode == "biaxial" || rule.affine_mode == "radial_outward" ||
                rule.affine_mode == "radial_inward" || rule.affine_mode == "biaxial_compression") {
                return 0.5 * (stress_tensor(0, 0) + stress_tensor(1, 1));
            }
            if (rule.affine_mode == "pure_shear") {
                return 0.5 * (stress_tensor(0, 0) - stress_tensor(1, 1));
            }
            if ((rule.affine_mode == "circumferential_shear" || rule.affine_mode == "torsion") && torque_available) {
                return computeRimTorque(userData, *current_points);
            }
            return stress_tensor(0, 1);
        }
    }
    return stress_tensor(0, 1);
}
