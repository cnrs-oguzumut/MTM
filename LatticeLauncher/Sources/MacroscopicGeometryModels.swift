import Foundation
import CoreGraphics

public enum OuterDomainShape: String, CaseIterable, Codable, Identifiable {
    case rectangle = "rectangle"
    case dogbone = "dogbone"
    case disk = "disk"
    case polygon = "polygon"
    
    public var id: String { rawValue }
    public var displayName: String {
        switch self {
        case .rectangle: return "Rectangular Plate"
        case .dogbone: return "Dogbone Tensile Specimen"
        case .disk: return "Circular Disk"
        case .polygon: return "Custom Polygon"
        }
    }
    public var icon: String {
        switch self {
        case .rectangle: return "rectangle"
        case .dogbone: return "figure.strengthtraining.traditional"
        case .disk: return "circle"
        case .polygon: return "skew"
        }
    }
}

public enum CutoutPrimitiveType: String, CaseIterable, Codable, Identifiable {
    case circle = "circle"
    case ellipse = "ellipse"
    case rectangle = "rectangle"
    case edgeNotch = "edge_notch"
    case preCrack = "pre_crack"
    case polygon = "polygon"
    
    public var id: String { rawValue }
    public var displayName: String {
        switch self {
        case .circle: return "Circular Hole"
        case .ellipse: return "Elliptical Hole"
        case .rectangle: return "Rectangular Slot"
        case .edgeNotch: return "Edge Notch (V / U)"
        case .preCrack: return "Pre-Crack Slit"
        case .polygon: return "Polygon Cutout"
        }
    }
    public var icon: String {
        switch self {
        case .circle: return "circle.dashed"
        case .ellipse: return "oval"
        case .rectangle: return "rectangle.dashed"
        case .edgeNotch: return "triangle"
        case .preCrack: return "bolt.horizontal"
        case .polygon: return "pentagon"
        }
    }
}

public struct NormalizedVertex: Codable, Equatable, Hashable, Identifiable {
    public var id: UUID = UUID()
    public var x: Double
    public var y: Double
    
    public init(x: Double, y: Double) {
        self.x = x
        self.y = y
    }
    
    enum CodingKeys: String, CodingKey {
        case x, y
    }
}

public struct GeometryCutoutItem: Codable, Equatable, Identifiable {
    public var id: String = UUID().uuidString
    public var type: CutoutPrimitiveType = .circle
    public var enabled: Bool = true
    public var cx: Double = 0.50
    public var cy: Double = 0.50
    public var rx: Double = 0.16
    public var ry: Double = 0.16
    public var rotation_deg: Double = 0.0
    public var notch_edge: String = "left"      // "left", "right", "top", "bottom"
    public var notch_profile: String = "v_notch" // "v_notch" or "u_notch"
    public var x1: Double = 0.0
    public var y1: Double = 0.50
    public var x2: Double = 0.35
    public var y2: Double = 0.50
    public var thickness: Double = 0.015
    public var vertices: [NormalizedVertex] = [
        NormalizedVertex(x: 0.40, y: 0.40),
        NormalizedVertex(x: 0.60, y: 0.40),
        NormalizedVertex(x: 0.55, y: 0.60),
        NormalizedVertex(x: 0.42, y: 0.58)
    ]
    
    public init(
        id: String = UUID().uuidString,
        type: CutoutPrimitiveType = .circle,
        enabled: Bool = true,
        cx: Double = 0.50,
        cy: Double = 0.50,
        rx: Double = 0.16,
        ry: Double = 0.16,
        rotation_deg: Double = 0.0,
        notch_edge: String = "left",
        notch_profile: String = "v_notch",
        x1: Double = 0.0,
        y1: Double = 0.50,
        x2: Double = 0.35,
        y2: Double = 0.50,
        thickness: Double = 0.015,
        vertices: [NormalizedVertex] = []
    ) {
        self.id = id
        self.type = type
        self.enabled = enabled
        self.cx = cx
        self.cy = cy
        self.rx = rx
        self.ry = ry
        self.rotation_deg = rotation_deg
        self.notch_edge = notch_edge
        self.notch_profile = notch_profile
        self.x1 = x1
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2
        self.thickness = thickness
        if !vertices.isEmpty {
            self.vertices = vertices
        }
    }
    
    public var summaryLabel: String {
        switch type {
        case .circle:
            return String(format: "Circle (%.2f, %.2f) R=%.2f", cx, cy, rx)
        case .ellipse:
            return String(format: "Ellipse (%.2f, %.2f) %.2f×%.2f ∠%.0f°", cx, cy, rx, ry, rotation_deg)
        case .rectangle:
            return String(format: "Slot (%.2f, %.2f) %.2f×%.2f ∠%.0f°", cx, cy, 2 * rx, 2 * ry, rotation_deg)
        case .edgeNotch:
            return String(format: "%@ %@ d=%.2f w=%.2f", notch_edge.capitalized, notch_profile == "u_notch" ? "U-Notch" : "V-Notch", rx, 2 * ry)
        case .preCrack:
            return String(format: "Crack (%.2f,%.2f)→(%.2f,%.2f)", x1, y1, x2, y2)
        case .polygon:
            return "Polygon (\(vertices.count) pts)"
        }
    }
}

public enum BoundaryRegionSelector: String, CaseIterable, Codable, Identifiable {
    case top = "top"
    case bottom = "bottom"
    case left = "left"
    case right = "right"
    case outerPerimeter = "outer_perimeter"
    case holeSurfaces = "hole_surfaces"
    case customBox = "custom_box"
    
    public var id: String { rawValue }
    public var displayName: String {
        switch self {
        case .top: return "Top Edge Grip"
        case .bottom: return "Bottom Edge Grip"
        case .left: return "Left Edge Grip"
        case .right: return "Right Edge Grip"
        case .outerPerimeter: return "All Outer Perimeter"
        case .holeSurfaces: return "Internal Hole Surfaces"
        case .customBox: return "Custom Box Zone"
        }
    }
}

public enum BoundaryConditionKindModel: String, CaseIterable, Codable, Identifiable {
    case free = "free"
    case fixed = "fixed"
    case prescribed = "prescribed"
    case affine = "affine"
    
    public var id: String { rawValue }
    public var displayName: String {
        switch self {
        case .free: return "Free Surface (Traction-Free)"
        case .fixed: return "Fixed / Clamped (u = 0)"
        case .prescribed: return "Prescribed Grip Δu(α)"
        case .affine: return "Affine Macro-Gradient F(α)"
        }
    }
    public var shortBadge: String {
        switch self {
        case .free: return "FREE"
        case .fixed: return "FIXED"
        case .prescribed: return "DRIVEN Δu"
        case .affine: return "AFFINE F(α)"
        }
    }
}

public struct BoundaryZoneRuleModel: Codable, Equatable, Identifiable {
    public var id: String = UUID().uuidString
    public var region: BoundaryRegionSelector = .top
    public var kind: BoundaryConditionKindModel = .prescribed
    public var enabled: Bool = true
    public var thickness_a0: Double = 2.0
    public var box_x_min: Double = 0.0
    public var box_y_min: Double = 0.90
    public var box_x_max: Double = 1.0
    public var box_y_max: Double = 1.0
    public var ux: Double = 0.0
    public var uy: Double = 1.0
    public var affine_mode: String = "simple_shear_x"
    
    public init(
        id: String = UUID().uuidString,
        region: BoundaryRegionSelector = .top,
        kind: BoundaryConditionKindModel = .prescribed,
        enabled: Bool = true,
        thickness_a0: Double = 2.0,
        box_x_min: Double = 0.0,
        box_y_min: Double = 0.90,
        box_x_max: Double = 1.0,
        box_y_max: Double = 1.0,
        ux: Double = 0.0,
        uy: Double = 1.0,
        affine_mode: String = "simple_shear_x"
    ) {
        self.id = id
        self.region = region
        self.kind = kind
        self.enabled = enabled
        self.thickness_a0 = thickness_a0
        self.box_x_min = box_x_min
        self.box_y_min = box_y_min
        self.box_x_max = box_x_max
        self.box_y_max = box_y_max
        self.ux = ux
        self.uy = uy
        self.affine_mode = affine_mode
    }
}

public enum IndenterShapeModel: String, CaseIterable, Codable, Identifiable {
    case cylinder = "cylinder"
    case flatPunch = "flat_punch"
    case wedge = "wedge"
    
    public var id: String { rawValue }
    public var displayName: String {
        switch self {
        case .cylinder: return "Cylindrical Tip (2D Circular)"
        case .flatPunch: return "Flat Punch Tip"
        case .wedge: return "Wedge / Berkovich Tip"
        }
    }
}

public enum IndenterContactModeModel: String, CaseIterable, Codable, Identifiable {
    case kinematic = "kinematic"
    case penalty = "penalty"
    
    public var id: String { rawValue }
    public var displayName: String {
        switch self {
        case .kinematic: return "Kinematic Projection (Option A - Hard)"
        case .penalty: return "Repulsive Penalty (Option B - Soft)"
        }
    }
    public var shortBadge: String {
        switch self {
        case .kinematic: return "HARD"
        case .penalty: return "PENALTY"
        }
    }
}

public struct IndenterConfigModel: Codable, Equatable {
    public var enabled: Bool = false
    public var shape: IndenterShapeModel = .cylinder
    public var contact_mode: IndenterContactModeModel = .kinematic
    public var cx: Double = 0.50
    public var radius_a0: Double = 25.0
    public var max_depth_a0: Double = 5.0
    public var cushion_subsurface: Bool = true
    public var cushion_offset_a0: Double = 0.20
    public var penalty_stiffness: Double = 100.0
    public var wedge_half_angle_deg: Double = 70.3
    
    public init(
        enabled: Bool = false,
        shape: IndenterShapeModel = .cylinder,
        contact_mode: IndenterContactModeModel = .kinematic,
        cx: Double = 0.50,
        radius_a0: Double = 25.0,
        max_depth_a0: Double = 5.0,
        cushion_subsurface: Bool = true,
        cushion_offset_a0: Double = 0.20,
        penalty_stiffness: Double = 100.0,
        wedge_half_angle_deg: Double = 70.3
    ) {
        self.enabled = enabled
        self.shape = shape
        self.contact_mode = contact_mode
        self.cx = cx
        self.radius_a0 = radius_a0
        self.max_depth_a0 = max_depth_a0
        self.cushion_subsurface = cushion_subsurface
        self.cushion_offset_a0 = cushion_offset_a0
        self.penalty_stiffness = penalty_stiffness
        self.wedge_half_angle_deg = wedge_half_angle_deg
    }
}

public enum GeometryEngineeringPreset: String, CaseIterable, Identifiable {
    case plateWithHoleTension = "plate_hole_tension"
    case plateWithHoleShear = "plate_hole_shear"
    case dogboneTension = "dogbone_tension"
    case doubleNotchTension = "dent_tension"
    case preCrackedSlab = "pre_cracked_mode1"
    case tiltedEllipseTension = "tilted_ellipse_tension"
    case porousVoidArray = "porous_void_array"
    case nanoIndentation = "nano_indentation"
    
    public var id: String { rawValue }
    public var title: String {
        switch self {
        case .plateWithHoleTension: return "Plate with Hole (Uniaxial Tension)"
        case .plateWithHoleShear: return "Plate with Hole (Simple Shear)"
        case .dogboneTension: return "Dogbone Specimen (Uniaxial Tension)"
        case .doubleNotchTension: return "Double-Edge Notched Tension (DENT)"
        case .preCrackedSlab: return "Pre-Cracked Mode-I Fracture Slab"
        case .tiltedEllipseTension: return "Tilted 45° Elliptical Slot (Tension)"
        case .porousVoidArray: return "Porous 2×2 Void Metamaterial (Shear)"
        case .nanoIndentation: return "Nanoindentation on Crystal Substrate"
        }
    }
    public var subtitle: String {
        switch self {
        case .plateWithHoleTension:
            return "Central circular void (R=0.16 L), bottom clamped, top pulled (+α Ly), left/right/hole traction-free."
        case .plateWithHoleShear:
            return "Central circular void (R=0.16 L), outer perimeter under affine simple shear, hole surface traction-free."
        case .dogboneTension:
            return "Filleted tensile gauge bar (w=0.50 Lx), top & bottom shoulders gripped under uniaxial tension."
        case .doubleNotchTension:
            return "Symmetric left & right V-notches focusing plastic slip bands at the specimen waist."
        case .preCrackedSlab:
            return "Edge slit crack (a=0.35 Lx) under vertical Mode-I opening tension on top/bottom grips."
        case .tiltedEllipseTension:
            return "Inclined elliptical cavity (45°) inducing mixed-mode shear & tensile dislocation emission."
        case .porousVoidArray:
            return "2×2 periodic-like pattern of circular cavities subject to macroscopic shear."
        case .nanoIndentation:
            return "Rigid cylindrical tip (R=25 a₀) penetrating top surface, clamped base, traction-free sides, P-δ curve."
        }
    }
}

public struct MacroscopicGeometryConfigModel: Codable, Equatable {
    public var enabled: Bool = false
    public var preset_name: String = GeometryEngineeringPreset.plateWithHoleTension.rawValue
    public var outer_shape: OuterDomainShape = .rectangle
    public var dogbone_gauge_width: Double = 0.50
    public var dogbone_gauge_length: Double = 0.44
    public var dogbone_fillet_radius: Double = 0.12
    public var disk_radius: Double = 0.48
    public var outer_polygon: [NormalizedVertex] = [
        NormalizedVertex(x: 0.05, y: 0.05),
        NormalizedVertex(x: 0.95, y: 0.05),
        NormalizedVertex(x: 0.95, y: 0.95),
        NormalizedVertex(x: 0.05, y: 0.95)
    ]
    public var cutouts: [GeometryCutoutItem] = [
        GeometryCutoutItem(type: .circle, cx: 0.50, cy: 0.50, rx: 0.16, ry: 0.16)
    ]
    public var boundary_rules: [BoundaryZoneRuleModel] = [
        BoundaryZoneRuleModel(region: .bottom, kind: .fixed, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
        BoundaryZoneRuleModel(region: .top, kind: .prescribed, thickness_a0: 2.0, ux: 0.0, uy: 1.0),
        BoundaryZoneRuleModel(region: .left, kind: .free, thickness_a0: 1.5, ux: 0.0, uy: 0.0),
        BoundaryZoneRuleModel(region: .right, kind: .free, thickness_a0: 1.5, ux: 0.0, uy: 0.0),
        BoundaryZoneRuleModel(region: .holeSurfaces, kind: .free, thickness_a0: 1.5, ux: 0.0, uy: 0.0)
    ]
    public var indenter: IndenterConfigModel = IndenterConfigModel()
    public var primary_stress_component: String = "auto"
    
    public init() {}
    
    public mutating func applyPreset(_ preset: GeometryEngineeringPreset) {
        self.enabled = true
        self.preset_name = preset.rawValue
        self.indenter.enabled = false
        switch preset {
        case .plateWithHoleTension:
            outer_shape = .rectangle
            cutouts = [
                GeometryCutoutItem(type: .circle, cx: 0.50, cy: 0.50, rx: 0.16, ry: 0.16)
            ]
            boundary_rules = [
                BoundaryZoneRuleModel(region: .bottom, kind: .fixed, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                BoundaryZoneRuleModel(region: .top, kind: .prescribed, thickness_a0: 2.0, ux: 0.0, uy: 1.0),
                BoundaryZoneRuleModel(region: .left, kind: .free, thickness_a0: 1.5),
                BoundaryZoneRuleModel(region: .right, kind: .free, thickness_a0: 1.5),
                BoundaryZoneRuleModel(region: .holeSurfaces, kind: .free, thickness_a0: 1.5)
            ]
            primary_stress_component = "sigma_yy"
            
        case .plateWithHoleShear:
            outer_shape = .rectangle
            cutouts = [
                GeometryCutoutItem(type: .circle, cx: 0.50, cy: 0.50, rx: 0.16, ry: 0.16)
            ]
            boundary_rules = [
                BoundaryZoneRuleModel(region: .outerPerimeter, kind: .affine, thickness_a0: 2.0, affine_mode: "simple_shear_x"),
                BoundaryZoneRuleModel(region: .holeSurfaces, kind: .free, thickness_a0: 1.5)
            ]
            primary_stress_component = "sigma_xy"
            
        case .dogboneTension:
            outer_shape = .dogbone
            dogbone_gauge_width = 0.50
            dogbone_gauge_length = 0.44
            dogbone_fillet_radius = 0.12
            cutouts = []
            boundary_rules = [
                BoundaryZoneRuleModel(region: .bottom, kind: .fixed, thickness_a0: 2.5, ux: 0.0, uy: 0.0),
                BoundaryZoneRuleModel(region: .top, kind: .prescribed, thickness_a0: 2.5, ux: 0.0, uy: 1.0)
            ]
            primary_stress_component = "sigma_yy"
            
        case .doubleNotchTension:
            outer_shape = .rectangle
            cutouts = [
                GeometryCutoutItem(type: .edgeNotch, cx: 0.0, cy: 0.50, rx: 0.18, ry: 0.10, notch_edge: "left", notch_profile: "v_notch"),
                GeometryCutoutItem(type: .edgeNotch, cx: 1.0, cy: 0.50, rx: 0.18, ry: 0.10, notch_edge: "right", notch_profile: "v_notch")
            ]
            boundary_rules = [
                BoundaryZoneRuleModel(region: .bottom, kind: .fixed, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                BoundaryZoneRuleModel(region: .top, kind: .prescribed, thickness_a0: 2.0, ux: 0.0, uy: 1.0)
            ]
            primary_stress_component = "sigma_yy"
            
        case .preCrackedSlab:
            outer_shape = .rectangle
            cutouts = [
                GeometryCutoutItem(type: .preCrack, x1: 0.0, y1: 0.50, x2: 0.36, y2: 0.50, thickness: 0.022)
            ]
            boundary_rules = [
                BoundaryZoneRuleModel(region: .bottom, kind: .fixed, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                BoundaryZoneRuleModel(region: .top, kind: .prescribed, thickness_a0: 2.0, ux: 0.0, uy: 1.0)
            ]
            primary_stress_component = "sigma_yy"
            
        case .tiltedEllipseTension:
            outer_shape = .rectangle
            cutouts = [
                GeometryCutoutItem(type: .ellipse, cx: 0.50, cy: 0.50, rx: 0.20, ry: 0.08, rotation_deg: 45.0)
            ]
            boundary_rules = [
                BoundaryZoneRuleModel(region: .bottom, kind: .fixed, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                BoundaryZoneRuleModel(region: .top, kind: .prescribed, thickness_a0: 2.0, ux: 0.0, uy: 1.0)
            ]
            primary_stress_component = "sigma_yy"
            
        case .porousVoidArray:
            outer_shape = .rectangle
            cutouts = [
                GeometryCutoutItem(type: .circle, cx: 0.32, cy: 0.32, rx: 0.11, ry: 0.11),
                GeometryCutoutItem(type: .circle, cx: 0.68, cy: 0.32, rx: 0.11, ry: 0.11),
                GeometryCutoutItem(type: .circle, cx: 0.32, cy: 0.68, rx: 0.11, ry: 0.11),
                GeometryCutoutItem(type: .circle, cx: 0.68, cy: 0.68, rx: 0.11, ry: 0.11)
            ]
            boundary_rules = [
                BoundaryZoneRuleModel(region: .bottom, kind: .fixed, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                BoundaryZoneRuleModel(region: .top, kind: .prescribed, thickness_a0: 2.0, ux: 1.0, uy: 0.0)
            ]
            primary_stress_component = "sigma_xy"
            
        case .nanoIndentation:
            outer_shape = .rectangle
            cutouts = []
            boundary_rules = [
                BoundaryZoneRuleModel(region: .bottom, kind: .fixed, thickness_a0: 2.5, ux: 0.0, uy: 0.0),
                BoundaryZoneRuleModel(region: .left, kind: .free, thickness_a0: 1.5),
                BoundaryZoneRuleModel(region: .right, kind: .free, thickness_a0: 1.5),
                BoundaryZoneRuleModel(region: .top, kind: .free, thickness_a0: 1.5)
            ]
            indenter.enabled = true
            indenter.shape = .cylinder
            indenter.contact_mode = .kinematic
            indenter.cx = 0.50
            indenter.radius_a0 = 25.0
            indenter.max_depth_a0 = 5.0
            indenter.cushion_subsurface = true
            indenter.cushion_offset_a0 = 0.20
            indenter.penalty_stiffness = 100.0
            primary_stress_component = "indentation_load"
        }
    }
    
    // MARK: - Live Client-Side Discrete Lattice Masking (Matches C++ MacroscopicGeometry)
    public func isPointInsideCutout(xHat: Double, yHat: Double, margin: Double = 0.0) -> Bool {
        for c in cutouts where c.enabled {
            let dx = xHat - c.cx
            let dy = yHat - c.cy
            switch c.type {
            case .circle:
                let r = c.rx + margin
                if dx * dx + dy * dy < r * r { return true }
            case .ellipse:
                let phi = c.rotation_deg * .pi / 180.0
                let cosP = cos(-phi), sinP = sin(-phi)
                let lx = cosP * dx - sinP * dy
                let ly = sinP * dx + cosP * dy
                let rx = max(1e-5, c.rx + margin)
                let ry = max(1e-5, c.ry + margin)
                if (lx * lx) / (rx * rx) + (ly * ly) / (ry * ry) < 1.0 { return true }
            case .rectangle:
                let phi = c.rotation_deg * .pi / 180.0
                let cosP = cos(-phi), sinP = sin(-phi)
                let lx = cosP * dx - sinP * dy
                let ly = sinP * dx + cosP * dy
                if abs(lx) < c.rx + margin && abs(ly) < c.ry + margin { return true }
            case .edgeNotch:
                let dEff = max(1e-5, c.rx + margin)
                let wEff = max(1e-5, c.ry + margin)
                var u = 1.0, v = 1.0
                if c.notch_edge == "left" { u = xHat; v = abs(yHat - c.cy) }
                else if c.notch_edge == "right" { u = 1.0 - xHat; v = abs(yHat - c.cy) }
                else if c.notch_edge == "bottom" { u = yHat; v = abs(xHat - c.cx) }
                else { u = 1.0 - yHat; v = abs(xHat - c.cx) }
                if u <= dEff {
                    if c.notch_profile == "u_notch" {
                        if v < wEff { return true }
                    } else {
                        if v < wEff * (1.0 - u / dEff) { return true }
                    }
                }
            case .preCrack:
                let ax = c.x1, ay = c.y1, bx = c.x2, by = c.y2
                let abx = bx - ax, aby = by - ay
                let len2 = abx * abx + aby * aby
                let t = len2 > 1e-12 ? min(1.0, max(0.0, ((xHat - ax) * abx + (yHat - ay) * aby) / len2)) : 0.0
                let px = ax + t * abx, py = ay + t * aby
                let dist = hypot(xHat - px, yHat - py)
                if dist < 0.5 * c.thickness + margin { return true }
            case .polygon:
                if Self.pointInPolygon(x: xHat, y: yHat, poly: c.vertices) { return true }
            }
        }
        return false
    }
    
    public func isPointInsideSolid(xHat: Double, yHat: Double) -> Bool {
        switch outer_shape {
        case .rectangle:
            if xHat < -1e-6 || xHat > 1.0 + 1e-6 || yHat < -1e-6 || yHat > 1.0 + 1e-6 { return false }
        case .dogbone:
            if xHat < -1e-6 || xHat > 1.0 + 1e-6 || yHat < -1e-6 || yHat > 1.0 + 1e-6 { return false }
            let halfGW = 0.5 * min(0.98, max(0.15, dogbone_gauge_width))
            let halfGL = 0.5 * min(0.90, max(0.10, dogbone_gauge_length))
            let rf = max(0.02, dogbone_fillet_radius)
            let dy = abs(yHat - 0.5)
            let halfWAllowed: Double
            if dy <= halfGL {
                halfWAllowed = halfGW
            } else if dy <= halfGL + rf {
                let u = dy - halfGL
                let filletDx = rf - sqrt(max(0.0, rf * rf - u * u))
                halfWAllowed = min(0.5, halfGW + filletDx * ((0.5 - halfGW) / rf))
            } else {
                halfWAllowed = 0.5
            }
            if abs(xHat - 0.5) > halfWAllowed + 1e-6 { return false }
        case .disk:
            let dx = xHat - 0.5, dy = yHat - 0.5
            if dx * dx + dy * dy > disk_radius * disk_radius { return false }
        case .polygon:
            if outer_polygon.count >= 3 && !Self.pointInPolygon(x: xHat, y: yHat, poly: outer_polygon) {
                return false
            }
        }
        if isPointNearCutoutNormalized(xHat: xHat, yHat: yHat) { return false }
        return true
    }
    
    private func isPointNearCutoutNormalized(xHat: Double, yHat: Double) -> Bool {
        return isPointInsideCutout(xHat: xHat, yHat: yHat, margin: 0.0)
    }
    
    /// Returns 0 = Free, 1 = Fixed, 2 = Driven
    public func classifyNodeRole(xHat: Double, yHat: Double, nx: Int, ny: Int) -> Int {
        var role = 0
        let a0x = 1.0 / Double(max(1, nx - 1))
        let a0y = 1.0 / Double(max(1, ny - 1))
        for rule in boundary_rules where rule.enabled {
            let stripX = max(0.35, rule.thickness_a0) * a0x
            let stripY = max(0.35, rule.thickness_a0) * a0y
            var matches = false
            switch rule.region {
            case .top: matches = (yHat >= 1.0 - stripY)
            case .bottom: matches = (yHat <= stripY)
            case .left: matches = (xHat <= stripX)
            case .right: matches = (xHat >= 1.0 - stripX)
            case .outerPerimeter:
                if outer_shape == .disk {
                    let dist = hypot(xHat - 0.5, yHat - 0.5)
                    matches = (dist >= disk_radius - max(stripX, stripY))
                } else {
                    matches = (xHat <= stripX || xHat >= 1.0 - stripX || yHat <= stripY || yHat >= 1.0 - stripY)
                }
            case .holeSurfaces:
                matches = isPointInsideCutout(xHat: xHat, yHat: yHat, margin: max(stripX, stripY))
            case .customBox:
                matches = (xHat >= rule.box_x_min && xHat <= rule.box_x_max &&
                           yHat >= rule.box_y_min && yHat <= rule.box_y_max)
            }
            if matches {
                switch rule.kind {
                case .free: role = 0
                case .fixed: role = 1
                case .prescribed, .affine: role = 2
                }
            }
        }
        if indenter.enabled {
            let rNorm = indenter.radius_a0 * a0x
            let cx = indenter.cx
            let cy = 1.0 + rNorm
            let dx = xHat - cx
            let dy = yHat - cy
            let dist = sqrt(dx * dx + dy * dy)
            if dist <= rNorm + 1.25 * a0y && yHat >= 1.0 - 2.5 * a0y {
                role = 2 // Contact zone highlight
            }
        }
        return role
    }
    
    public static func pointInPolygon(x: Double, y: Double, poly: [NormalizedVertex]) -> Bool {
        guard poly.count >= 3 else { return false }
        var inside = false
        var j = poly.count - 1
        for i in 0..<poly.count {
            let xi = poly[i].x, yi = poly[i].y
            let xj = poly[j].x, yj = poly[j].y
            let intersect = ((yi > y) != (yj > y)) &&
                (x < (xj - xi) * (y - yi) / max(1e-15, (yj - yi)) + xi)
            if intersect { inside.toggle() }
            j = i
        }
        return inside
    }
    
    public func writeJSON(to path: String) throws {
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.prettyPrinted, .sortedKeys]
        let data = try encoder.encode(self)
        let url = URL(fileURLWithPath: path)
        try FileManager.default.createDirectory(at: url.deletingLastPathComponent(), withIntermediateDirectories: true)
        try data.write(to: url)
    }
}
