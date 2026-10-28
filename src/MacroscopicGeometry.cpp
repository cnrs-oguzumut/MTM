#include "../include/geometry/MacroscopicGeometry.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <random>
#include <cctype>

namespace {

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
        bool intersect = ((yi > y) != (yj > y)) &&
                         (x < (xj - xi) * (y - yi) / std::max(1e-15, (yj - yi)) + xi);
        if (intersect) inside = !inside;
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

    // Parse boundary_rules array
    std::string rules_raw = extractRawValue(json, "boundary_rules");
    for (const auto& item_json : splitJSONObjectArray(rules_raw)) {
        BoundaryZoneRule r;
        r.id = extractString(item_json, "id", "bc");
        r.enabled = extractBool(item_json, "enabled", true);
        std::string reg_str = extractString(item_json, "region", "top");
        if (reg_str == "bottom" || reg_str == "bottomEdge") r.region = BoundaryRegionType::BottomEdge;
        else if (reg_str == "left" || reg_str == "leftEdge") r.region = BoundaryRegionType::LeftEdge;
        else if (reg_str == "right" || reg_str == "rightEdge") r.region = BoundaryRegionType::RightEdge;
        else if (reg_str == "outer" || reg_str == "outer_perimeter" || reg_str == "outerPerimeter") r.region = BoundaryRegionType::OuterPerimeter;
        else if (reg_str == "hole" || reg_str == "hole_surfaces" || reg_str == "holeSurfaces") r.region = BoundaryRegionType::HoleSurfaces;
        else if (reg_str == "custom_box" || reg_str == "customBox" || reg_str == "box") r.region = BoundaryRegionType::CustomBox;
        else r.region = BoundaryRegionType::TopEdge;

        std::string kind_str = extractString(item_json, "kind", "free");
        if (kind_str == "fixed" || kind_str == "fixedClamped") r.kind = BoundaryConditionKind::FixedClamped;
        else if (kind_str == "prescribed" || kind_str == "prescribedDisplacement" || kind_str == "driven") r.kind = BoundaryConditionKind::PrescribedDisplacement;
        else if (kind_str == "affine" || kind_str == "affineDeformation") r.kind = BoundaryConditionKind::AffineDeformation;
        else r.kind = BoundaryConditionKind::FreeSurface;

        r.thickness_a0 = extractDouble(item_json, "thickness_a0", 2.0);
        r.box_x_min = extractDouble(item_json, "box_x_min", 0.0);
        r.box_y_min = extractDouble(item_json, "box_y_min", 0.0);
        r.box_x_max = extractDouble(item_json, "box_x_max", 1.0);
        r.box_y_max = extractDouble(item_json, "box_y_max", 0.1);
        r.ux = extractDouble(item_json, "ux", 0.0);
        r.uy = extractDouble(item_json, "uy", 0.0);
        r.affine_mode = extractString(item_json, "affine_mode", "simple_shear_x");
        out_config.boundary_rules.push_back(r);
    }

    std::cout << "[MacroscopicGeometry] Loaded '" << out_config.preset_name << "' from " << filepath
              << " (" << out_config.cutouts.size() << " cutouts, "
              << out_config.boundary_rules.size() << " boundary rules)\n";
    return true;
}

void MacroscopicGeometryConfig::setPhysicalBoundsFromPoints(
    const std::vector<Point2D>& raw_points, double lattice_a0) {
    a0 = (lattice_a0 > 1e-8) ? lattice_a0 : 1.0;
    if (raw_points.empty()) {
        x_min = 0.0;
        y_min = 0.0;
        Lx = 100.0;
        Ly = 100.0;
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
}

bool MacroscopicGeometryConfig::isPointNearCutout(double x, double y, double margin_phys) const {
    const double min_dim = std::min(Lx, Ly);
    const double x_hat = (x - x_min) / Lx;
    const double y_hat = (y - y_min) / Ly;

    for (const auto& c : cutouts) {
        if (!c.enabled) continue;
        const double cx_phys = x_min + c.cx * Lx;
        const double cy_phys = y_min + c.cy * Ly;
        const double dx = x - cx_phys;
        const double dy = y - cy_phys;

        switch (c.type) {
        case CutoutShapeType::Circle: {
            double r_phys = c.rx * min_dim + margin_phys;
            if (dx * dx + dy * dy < r_phys * r_phys) return true;
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
            if ((lx * lx) / (rx_phys * rx_phys) + (ly * ly) / (ry_phys * ry_phys) < 1.0) return true;
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
            if (std::abs(lx) < hw_phys && std::abs(ly) < hh_phys) return true;
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
                        if (v < w_eff) return true;
                    } else {
                        // V-notch tapering linearly from half-width w_eff at u=0 to 0 at u=d_eff
                        if (v < w_eff * (1.0 - u / std::max(1e-6, d_eff))) return true;
                    }
                }
            } else if (c.notch_edge == "right") {
                double u = 1.0 - x_hat;
                double v = std::abs(y_hat - c.cy);
                double d_eff = depth_norm + margin_x;
                double w_eff = std::max(1e-6, half_w_norm + margin_y);
                if (u <= d_eff) {
                    if (c.notch_profile == "u_notch") {
                        if (v < w_eff) return true;
                    } else {
                        if (v < w_eff * (1.0 - u / std::max(1e-6, d_eff))) return true;
                    }
                }
            } else if (c.notch_edge == "bottom") {
                double u = y_hat;
                double v = std::abs(x_hat - c.cx);
                double d_eff = depth_norm + margin_y;
                double w_eff = std::max(1e-6, half_w_norm + margin_x);
                if (u <= d_eff) {
                    if (c.notch_profile == "u_notch") {
                        if (v < w_eff) return true;
                    } else {
                        if (v < w_eff * (1.0 - u / std::max(1e-6, d_eff))) return true;
                    }
                }
            } else if (c.notch_edge == "top") {
                double u = 1.0 - y_hat;
                double v = std::abs(x_hat - c.cx);
                double d_eff = depth_norm + margin_y;
                double w_eff = std::max(1e-6, half_w_norm + margin_x);
                if (u <= d_eff) {
                    if (c.notch_profile == "u_notch") {
                        if (v < w_eff) return true;
                    } else {
                        if (v < w_eff * (1.0 - u / std::max(1e-6, d_eff))) return true;
                    }
                }
            }
            break;
        }
        case CutoutShapeType::PreCrack: {
            Eigen::Vector2d a(x_min + c.x1 * Lx, y_min + c.y1 * Ly);
            Eigen::Vector2d b(x_min + c.x2 * Lx, y_min + c.y2 * Ly);
            double half_t = 0.5 * c.thickness * min_dim + margin_phys;
            if (pointToSegmentDistance(Eigen::Vector2d(x, y), a, b) < half_t) return true;
            break;
        }
        case CutoutShapeType::Polygon: {
            if (isPointInPolygon(x_hat, y_hat, c.vertices)) return true;
            break;
        }
        }
    }
    return false;
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

    // 3. Check if any edge crosses a PreCrack slit segment
    for (const auto& c : cutouts) {
        if (!c.enabled || c.type != CutoutShapeType::PreCrack) continue;
        Eigen::Vector2d a(x_min + c.x1 * Lx, y_min + c.y1 * Ly);
        Eigen::Vector2d b(x_min + c.x2 * Lx, y_min + c.y2 * Ly);
        if (segmentsIntersect(p0, p1, a, b) ||
            segmentsIntersect(p1, p2, a, b) ||
            segmentsIntersect(p2, p0, a, b)) {
            return true;
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

    if (n_pts == 0) return;

    const double x_max = x_min + Lx;
    const double y_max = y_min + Ly;
    const Eigen::Vector2d center(x_min + 0.5 * Lx, y_min + 0.5 * Ly);

    // Evaluate boundary rules for each node
    for (size_t i = 0; i < n_pts; ++i) {
        const double x = carved_points[i].coord.x();
        const double y = carved_points[i].coord.y();
        const double x_hat = (x - x_min) / Lx;
        const double y_hat = (y - y_min) / Ly;

        for (const auto& rule : boundary_rules) {
            if (!rule.enabled) continue;
            const double strip = std::max(0.35, rule.thickness_a0) * a0;
            bool matches = false;

            switch (rule.region) {
            case BoundaryRegionType::TopEdge:
                matches = (y >= y_max - strip);
                break;
            case BoundaryRegionType::BottomEdge:
                matches = (y <= y_min + strip);
                break;
            case BoundaryRegionType::LeftEdge:
                matches = (x <= x_min + strip);
                break;
            case BoundaryRegionType::RightEdge:
                matches = (x >= x_max - strip);
                break;
            case BoundaryRegionType::OuterPerimeter:
                if (outer_shape == OuterDomainShapeType::Disk) {
                    double r_phys = disk_radius * std::min(Lx, Ly);
                    double dist = (carved_points[i].coord - center).norm();
                    matches = (dist >= r_phys - strip);
                } else {
                    matches = (x <= x_min + strip || x >= x_max - strip ||
                               y <= y_min + strip || y >= y_max - strip);
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

            if (rule.kind == BoundaryConditionKind::FreeSurface) {
                // FreeSurface (Neumann natural BC) marks a traction-free boundary.
                // It MUST NEVER overwrite an already-constrained Dirichlet condition (Fixed, Prescribed, or Affine)!
                if (node_bc_role[i] == 0) {
                    node_bc_role[i] = 0;
                    node_bc_velocity[i].setZero();
                }
            } else if (rule.kind == BoundaryConditionKind::FixedClamped) {
                node_bc_role[i] = 1;
                node_bc_velocity[i].setZero();
            } else if (rule.kind == BoundaryConditionKind::PrescribedDisplacement) {
                node_bc_role[i] = 2;
                // ux, uy are normalized by Ly so alpha = 0.01 means 1% engineering strain
                node_bc_velocity[i] = Eigen::Vector2d(rule.ux * Ly, rule.uy * Ly);
            } else if (rule.kind == BoundaryConditionKind::AffineDeformation) {
                node_bc_role[i] = 2;
                Eigen::Vector2d rel(x - x_min, y - y_min);
                if (rule.affine_mode == "simple_shear_y") {
                    node_bc_velocity[i] = Eigen::Vector2d(0.0, rel.x());
                } else if (rule.affine_mode == "uniaxial_y") {
                    node_bc_velocity[i] = Eigen::Vector2d(0.0, rel.y());
                } else if (rule.affine_mode == "uniaxial_x") {
                    node_bc_velocity[i] = Eigen::Vector2d(rel.x(), 0.0);
                } else if (rule.affine_mode == "pure_shear") {
                    node_bc_velocity[i] = Eigen::Vector2d(rel.x(), -rel.y());
                } else if (rule.affine_mode == "biaxial") {
                    node_bc_velocity[i] = Eigen::Vector2d(rel.x(), rel.y());
                } else {
                    // default: simple_shear_x -> u_x = alpha * (y - y_min)
                    node_bc_velocity[i] = Eigen::Vector2d(rel.y(), 0.0);
                }
            }
        }
    }

    // Build interior_mapping and full_mapping
    interior_mapping.reserve(n_pts);
    full_mapping.reserve(n_pts);
    int solver_idx = 0;
    int count_free = 0, count_fixed = 0, count_driven = 0;

    for (size_t i = 0; i < n_pts; ++i) {
        if (node_bc_role[i] == 0) {
            interior_mapping.emplace_back(static_cast<int>(i), solver_idx);
            full_mapping.emplace_back(static_cast<int>(i), solver_idx);
            ++solver_idx;
            ++count_free;
        } else {
            full_mapping.emplace_back(static_cast<int>(i), -1);
            if (node_bc_role[i] == 1) ++count_fixed;
            else ++count_driven;
        }
    }

    std::cout << "[MacroscopicGeometry] Boundary DOFs classified: "
              << count_free << " Free (interior + traction-free surfaces), "
              << count_fixed << " Fixed (clamped), "
              << count_driven << " Driven (prescribed BC)\n";
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
            } else if (rule.affine_mode == "biaxial") {
                F(0, 0) = 1.0 + alpha;
                F(1, 1) = 1.0 + alpha;
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
                F(0, 0) += sign * alpha * rule.ux;
                F(1, 0) += sign * alpha * rule.uy;
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
    Eigen::Matrix2d dF_macro = computeMacroscopicF(dalpha) - Eigen::Matrix2d::Identity();
    Eigen::Vector2d origin(x_min, y_min);

    std::mt19937 gen(seed);
    std::normal_distribution<double> noise_dist(0.0, std::max(0.0, noise_level));

    for (size_t j = 0; j < n_pts; ++j) {
        Eigen::Vector2d noise_vec = (is_first_step && noise_level > 0.0)
            ? Eigen::Vector2d(noise_dist(gen), noise_dist(gen))
            : Eigen::Vector2d::Zero();

        if (node_bc_role[j] == 1) {
            // Fixed / Clamped node: stays strictly at reference position X_0
            current_points[j].coord = points_ref[j].coord;
        } else if (node_bc_role[j] == 2) {
            // Driven / Prescribed node: exact boundary trajectory X_0 + alpha * v_j
            current_points[j].coord = points_ref[j].coord + alpha * node_bc_velocity[j];
        } else {
            // Free node (interior or traction-free hole/outer surface):
            // Warm-start predictor using smooth macroscopic increment
            Eigen::Vector2d rel = points_ref[j].coord - origin;
            if (is_first_step) {
                current_points[j].coord = origin + F_macro_alpha * rel + noise_vec;
            } else {
                current_points[j].coord += dF_macro * rel;
            }
        }
    }
}

double MacroscopicGeometryConfig::extractPrimaryStress(const Eigen::Matrix2d& stress_tensor) const {
    if (!enabled) return stress_tensor(0, 1);

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
            if (rule.affine_mode == "biaxial") return 0.5 * (stress_tensor(0, 0) + stress_tensor(1, 1));
            return stress_tensor(0, 1);
        }
    }
    return stress_tensor(0, 1);
}
