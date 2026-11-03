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
    case rollerX = "roller_x"
    case rollerY = "roller_y"
    case prescribed = "prescribed"
    case affine = "affine"
    
    public var id: String { rawValue }
    public var displayName: String {
        switch self {
        case .free: return "Free Surface (Traction-Free)"
        case .fixed: return "Fixed / Clamped (u = 0)"
        case .rollerX: return "Roller X (Free ux, Fix uy = 0)"
        case .rollerY: return "Roller Y (Fix ux = 0, Free uy)"
        case .prescribed: return "Prescribed Grip Δu(α)"
        case .affine: return "Affine Macro-Gradient F(α)"
        }
    }
    public var shortBadge: String {
        switch self {
        case .free: return "FREE"
        case .fixed: return "FIXED"
        case .rollerX: return "ROLLER-X (uy=0)"
        case .rollerY: return "ROLLER-Y (ux=0)"
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
    public var lock_x: Bool = true
    public var lock_y: Bool = true
    public var affine_mode: String = "simple_shear_x"
    
    enum CodingKeys: String, CodingKey {
        case id, region, kind, enabled, thickness_a0
        case box_x_min, box_y_min, box_x_max, box_y_max
        case ux, uy, lock_x, lock_y, affine_mode
    }

    public static func snapUnitDirection(_ val: Double) -> Double {
        if val <= -0.25 { return -1.0 }
        if val >= 0.25 { return 1.0 }
        return 0.0
    }

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
        lock_x: Bool = true,
        lock_y: Bool = true,
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
        self.ux = Self.snapUnitDirection(ux)
        self.uy = Self.snapUnitDirection(uy)
        if kind == .rollerX {
            self.lock_x = false
            self.lock_y = true
        } else if kind == .rollerY {
            self.lock_x = true
            self.lock_y = false
        } else {
            self.lock_x = lock_x
            self.lock_y = lock_y
        }
        self.affine_mode = affine_mode
    }

    public init(from decoder: Decoder) throws {
        let c = try decoder.container(keyedBy: CodingKeys.self)
        self.id = try c.decodeIfPresent(String.self, forKey: .id) ?? UUID().uuidString
        self.region = try c.decodeIfPresent(BoundaryRegionSelector.self, forKey: .region) ?? .top
        self.kind = try c.decodeIfPresent(BoundaryConditionKindModel.self, forKey: .kind) ?? .prescribed
        self.enabled = try c.decodeIfPresent(Bool.self, forKey: .enabled) ?? true
        self.thickness_a0 = try c.decodeIfPresent(Double.self, forKey: .thickness_a0) ?? 2.0
        self.box_x_min = try c.decodeIfPresent(Double.self, forKey: .box_x_min) ?? 0.0
        self.box_y_min = try c.decodeIfPresent(Double.self, forKey: .box_y_min) ?? 0.90
        self.box_x_max = try c.decodeIfPresent(Double.self, forKey: .box_x_max) ?? 1.0
        self.box_y_max = try c.decodeIfPresent(Double.self, forKey: .box_y_max) ?? 1.0
        self.ux = try c.decodeIfPresent(Double.self, forKey: .ux) ?? 0.0
        self.uy = try c.decodeIfPresent(Double.self, forKey: .uy) ?? 1.0
        self.lock_x = try c.decodeIfPresent(Bool.self, forKey: .lock_x) ?? (self.kind != .rollerX)
        self.lock_y = try c.decodeIfPresent(Bool.self, forKey: .lock_y) ?? (self.kind != .rollerY)
        self.affine_mode = try c.decodeIfPresent(String.self, forKey: .affine_mode) ?? "simple_shear_x"
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

public enum MicrostructureModeModel: String, CaseIterable, Codable, Identifiable {
    case singleCrystal = "single_crystal"
    case bicrystal = "bicrystal"
    case voronoi = "voronoi_polycrystal"
    case periodicDislocations = "periodic_dislocations"
    
    public var id: String { rawValue }
    public var displayName: String {
        switch self {
        case .singleCrystal: return "Single Crystal"
        case .bicrystal: return "Bicrystal (Planar GB)"
        case .voronoi: return "2D Voronoi Polycrystal"
        case .periodicDislocations: return "Periodic Dislocations (2M Glide Planes)"
        }
    }
    public var shortBadge: String {
        switch self {
        case .singleCrystal: return "SINGLE"
        case .bicrystal: return "BICRYSTAL"
        case .voronoi: return "VORONOI"
        case .periodicDislocations: return "DISLOC (2M)"
        }
    }
}

public enum SeamOverlapResolutionModel: String, CaseIterable, Codable, Identifiable {
    case mergeMidpoint = "merge_midpoint"
    case deleteSecond = "delete_second"
    
    public var id: String { rawValue }
    public var displayName: String {
        switch self {
        case .mergeMidpoint: return "Merge Pair to Midpoint (Recommended)"
        case .deleteSecond: return "Delete Grain 2 Atom"
        }
    }
}

public struct VoronoiGrainSeedModel: Codable, Equatable, Identifiable {
    public var id: Int
    public var cx: Double
    public var cy: Double
    public var theta_deg: Double
    
    public init(id: Int = 0, cx: Double = 0.5, cy: Double = 0.5, theta_deg: Double = 0.0) {
        self.id = id
        self.cx = cx
        self.cy = cy
        self.theta_deg = theta_deg
    }
}

public struct SeededDislocationModel: Codable, Equatable, Identifiable {
    public var id: Int
    public var glide_plane_row: Int
    public var x_frac: Double
    public var burgers_sign: Int
    public var enabled: Bool
    
    public init(id: Int = 0, glide_plane_row: Int = 25, x_frac: Double = 0.50, burgers_sign: Int = 1, enabled: Bool = true) {
        self.id = id
        self.glide_plane_row = glide_plane_row
        self.x_frac = x_frac
        self.burgers_sign = burgers_sign
        self.enabled = enabled
    }
}

public enum DislocationLayoutPreset: String, CaseIterable, Codable, Identifiable {
    case alternatingRandomX = "alternating_random_x"
    case alternatingStaggered = "alternating_staggered"
    case randomBalanced = "random_planes_random_signs"
    case custom = "custom"
    
    public var id: String { rawValue }
    public var displayName: String {
        switch self {
        case .alternatingRandomX: return "Alternating ±b (Random X)"
        case .alternatingStaggered: return "Alternating ±b (Staggered X)"
        case .randomBalanced: return "Random Planes & ±b (∑b = 0)"
        case .custom: return "Custom Dislocation Table"
        }
    }
}

public enum SimulationProtocolMode: String, CaseIterable, Codable, Identifiable {
    case relaxOnly = "relax_only"
    case relaxThenLoad = "relax_then_load"
    case sweepMisorientation = "sweep_misorientation"
    
    public var id: String { rawValue }
    public var shortTitle: String {
        switch self {
        case .relaxOnly: return "1. Relax Only"
        case .relaxThenLoad: return "2. Relax → Load"
        case .sweepMisorientation: return "3. Sweep Δθ"
        }
    }
}

public enum Stage2LoadingBCPreset: String, CaseIterable, Codable, Identifiable {
    case xGripsYRoller = "x_grips_y_roller"
    case affineFromRelaxed = "affine_from_relaxed"
    case yGripsXRoller = "y_grips_x_roller"
    case keepStudio = "keep_studio"
    
    public var id: String { rawValue }
    public var displayName: String {
        switch self {
        case .xGripsYRoller: return "X-Grips Driven + Y-Roller (Left Fixed, Right Driven)"
        case .affineFromRelaxed: return "Affine Driven from Relaxed X₀ (All 4 Walls)"
        case .yGripsXRoller: return "Y-Grips Driven + X-Roller (Bottom Fixed, Top Driven)"
        case .keepStudio: return "Keep Custom Studio Boundary Rules"
        }
    }
}

public enum Stage2LoadTypePreset: String, CaseIterable, Codable, Identifiable {
    case shearParallelGB = "shear_parallel_gb"
    case tensionNormalX = "tension_normal_x"
    case compressionNormalX = "compression_normal_x"
    case shearHorizontal = "shear_horizontal"
    case tensionVerticalY = "tension_vertical_y"
    case compressionVerticalY = "compression_vertical_y"
    
    public var id: String { rawValue }
    public var displayName: String {
        switch self {
        case .shearParallelGB: return "Shear Parallel to Vertical GB (+uy on Right Grip)"
        case .tensionNormalX: return "Tension Normal to Vertical GB (+ux on Right Grip)"
        case .compressionNormalX: return "Compression Normal to Vertical GB (-ux on Right Grip)"
        case .shearHorizontal: return "Horizontal Simple Shear (+ux on Top Grip)"
        case .tensionVerticalY: return "Vertical Tension (+uy on Top Grip)"
        case .compressionVerticalY: return "Vertical Compression (-uy on Top Grip)"
        }
    }
    
    public var recommendedGripPreset: Stage2LoadingBCPreset {
        switch self {
        case .shearParallelGB, .tensionNormalX, .compressionNormalX:
            return .xGripsYRoller
        case .shearHorizontal, .tensionVerticalY, .compressionVerticalY:
            return .yGripsXRoller
        }
    }
    
    public var primaryStressComponent: String {
        switch self {
        case .shearParallelGB, .shearHorizontal:
            return "sigma_xy"
        case .tensionNormalX, .compressionNormalX:
            return "sigma_xx"
        case .tensionVerticalY, .compressionVerticalY:
            return "sigma_yy"
        }
    }
}

public struct MicrostructureConfigModel: Codable, Equatable {
    public var mode: MicrostructureModeModel = .singleCrystal
    public var seam_resolution: SeamOverlapResolutionModel = .mergeMidpoint
    public var min_seam_distance_a0: Double = 0.55
    public var max_gb_edge_a0: Double = 2.60
    public var pre_relax_gb: Bool = true
    public var relax_only: Bool = false
    public var reanchor_ref_after_relax: Bool = true
    public var max_relax_remesh_passes: Int = 5
    public var loading_bc_preset: String = Stage2LoadingBCPreset.xGripsYRoller.rawValue
    public var stage2_load_type: String = Stage2LoadTypePreset.shearParallelGB.rawValue
    
    // Bicrystal parameters
    public var gb_cx: Double = 0.50
    public var gb_cy: Double = 0.50
    public var gb_angle_deg: Double = 0.0 // 0° = Vertical GB (x = cx), 90° = Horizontal GB (y = cy)
    public var theta1_deg: Double = -15.0
    public var theta2_deg: Double = +15.0
    public var gb_normal_gap_a0: Double = 0.0
    public var gb_tangential_shift_a0: Double = 0.0
    
    // Voronoi Polycrystal parameters
    public var num_grains: Int = 6
    public var random_seed: Int = 42
    public var theta_min_deg: Double = -30.0
    public var theta_max_deg: Double = +30.0
    public var seeds: [VoronoiGrainSeedModel] = []
    
    // Multi-Plane Single-Dislocation Periodic Setup (2M active glide planes, 1 dislocation per plane, sum(b) = 0)
    public var dislocation_pair_count: Int = 2
    public var dislocation_layout: String = DislocationLayoutPreset.alternatingRandomX.rawValue
    public var min_plane_spacing_rows: Int = 6
    public var dislocation_core_radius_a0: Double = 0.8
    public var bound_strip_excess: Bool = true
    public var enforce_periodic_bc: Bool = true
    public var dislocations: [SeededDislocationModel] = []
    
    enum CodingKeys: String, CodingKey {
        case mode
        case seam_resolution
        case min_seam_distance_a0
        case max_gb_edge_a0
        case pre_relax_gb
        case relax_only
        case reanchor_ref_after_relax
        case max_relax_remesh_passes
        case loading_bc_preset
        case stage2_load_type
        case gb_cx
        case gb_cy
        case gb_angle_deg
        case theta1_deg
        case theta2_deg
        case gb_normal_gap_a0
        case gb_tangential_shift_a0
        case num_grains
        case random_seed
        case theta_min_deg
        case theta_max_deg
        case seeds
        case dislocation_pair_count
        case dislocation_layout
        case min_plane_spacing_rows
        case dislocation_core_radius_a0
        case bound_strip_excess
        case enforce_periodic_bc
        case dislocations
    }
    
    public init() {
        self.regenerateSeeds()
        self.regeneratePeriodicDislocations()
    }
    
    public init(from decoder: Decoder) throws {
        let c = try decoder.container(keyedBy: CodingKeys.self)
        self.mode = try c.decodeIfPresent(MicrostructureModeModel.self, forKey: .mode) ?? .singleCrystal
        self.seam_resolution = try c.decodeIfPresent(SeamOverlapResolutionModel.self, forKey: .seam_resolution) ?? .mergeMidpoint
        self.min_seam_distance_a0 = try c.decodeIfPresent(Double.self, forKey: .min_seam_distance_a0) ?? 0.55
        self.max_gb_edge_a0 = try c.decodeIfPresent(Double.self, forKey: .max_gb_edge_a0) ?? 2.60
        self.pre_relax_gb = try c.decodeIfPresent(Bool.self, forKey: .pre_relax_gb) ?? true
        self.relax_only = try c.decodeIfPresent(Bool.self, forKey: .relax_only) ?? false
        self.reanchor_ref_after_relax = try c.decodeIfPresent(Bool.self, forKey: .reanchor_ref_after_relax) ?? true
        self.max_relax_remesh_passes = try c.decodeIfPresent(Int.self, forKey: .max_relax_remesh_passes) ?? 5
        self.loading_bc_preset = try c.decodeIfPresent(String.self, forKey: .loading_bc_preset) ?? Stage2LoadingBCPreset.xGripsYRoller.rawValue
        self.stage2_load_type = try c.decodeIfPresent(String.self, forKey: .stage2_load_type) ?? Stage2LoadTypePreset.shearParallelGB.rawValue
        self.gb_cx = try c.decodeIfPresent(Double.self, forKey: .gb_cx) ?? 0.50
        self.gb_cy = try c.decodeIfPresent(Double.self, forKey: .gb_cy) ?? 0.50
        self.gb_angle_deg = try c.decodeIfPresent(Double.self, forKey: .gb_angle_deg) ?? 0.0
        self.theta1_deg = try c.decodeIfPresent(Double.self, forKey: .theta1_deg) ?? -15.0
        self.theta2_deg = try c.decodeIfPresent(Double.self, forKey: .theta2_deg) ?? +15.0
        self.gb_normal_gap_a0 = try c.decodeIfPresent(Double.self, forKey: .gb_normal_gap_a0) ?? 0.0
        self.gb_tangential_shift_a0 = try c.decodeIfPresent(Double.self, forKey: .gb_tangential_shift_a0) ?? 0.0
        self.num_grains = try c.decodeIfPresent(Int.self, forKey: .num_grains) ?? 6
        self.random_seed = try c.decodeIfPresent(Int.self, forKey: .random_seed) ?? 42
        self.theta_min_deg = try c.decodeIfPresent(Double.self, forKey: .theta_min_deg) ?? -30.0
        self.theta_max_deg = try c.decodeIfPresent(Double.self, forKey: .theta_max_deg) ?? +30.0
        if let decodedSeeds = try c.decodeIfPresent([VoronoiGrainSeedModel].self, forKey: .seeds), !decodedSeeds.isEmpty {
            self.seeds = decodedSeeds
        } else {
            self.regenerateSeeds()
        }
        self.dislocation_pair_count = try c.decodeIfPresent(Int.self, forKey: .dislocation_pair_count) ?? 2
        self.dislocation_layout = try c.decodeIfPresent(String.self, forKey: .dislocation_layout) ?? DislocationLayoutPreset.alternatingRandomX.rawValue
        self.min_plane_spacing_rows = try c.decodeIfPresent(Int.self, forKey: .min_plane_spacing_rows) ?? 6
        self.dislocation_core_radius_a0 = try c.decodeIfPresent(Double.self, forKey: .dislocation_core_radius_a0) ?? 0.8
        self.bound_strip_excess = try c.decodeIfPresent(Bool.self, forKey: .bound_strip_excess) ?? true
        self.enforce_periodic_bc = try c.decodeIfPresent(Bool.self, forKey: .enforce_periodic_bc) ?? true
        if let decodedDislocs = try c.decodeIfPresent([SeededDislocationModel].self, forKey: .dislocations), !decodedDislocs.isEmpty {
            self.dislocations = decodedDislocs
        } else {
            self.regeneratePeriodicDislocations()
        }
    }
    
    public mutating func regenerateSeeds() {
        let ng = max(2, min(32, num_grains))
        let cols = max(1, Int(ceil(sqrt(Double(ng)))))
        let rows = max(1, Int(ceil(Double(ng) / Double(cols))))
        let tMin = min(theta_min_deg, theta_max_deg)
        let tMax = max(theta_min_deg, theta_max_deg)
        var state = UInt64(bitPattern: Int64(random_seed == 0 ? 42 : random_seed))
        func nextRand() -> Double {
            state = state &* 6364136223846793005 &+ 1442695040888963407
            return Double((state >> 11) & 0x1FFFFF) / Double(0x1FFFFF)
        }
        var out: [VoronoiGrainSeedModel] = []
        for k in 0..<ng {
            let r = k / cols
            let c = k % cols
            let jx = 0.20 + 0.60 * nextRand()
            let jy = 0.20 + 0.60 * nextRand()
            let cx = min(0.94, max(0.06, (Double(c) + jx) / Double(cols)))
            let cy = min(0.94, max(0.06, (Double(r) + jy) / Double(rows)))
            let th = tMin + (tMax - tMin) * nextRand()
            out.append(VoronoiGrainSeedModel(
                id: k,
                cx: (cx * 1000).rounded() / 1000,
                cy: (cy * 1000).rounded() / 1000,
                theta_deg: (th * 10).rounded() / 10
            ))
        }
        self.seeds = out
    }
    
    public mutating func regeneratePeriodicDislocations(nx: Int = 60, ny: Int = 60) {
        let M = max(1, dislocation_pair_count)
        let totalPlanes = 2 * M
        let minRow = 1
        let maxRow = max(minRow + totalPlanes, ny - 2)
        
        var planeRows: [Int] = []
        if dislocation_layout == DislocationLayoutPreset.alternatingStaggered.rawValue || dislocation_layout == DislocationLayoutPreset.alternatingRandomX.rawValue {
            for k in 0..<totalPlanes {
                let idealRow = (Double(k) + 0.5) * (Double(ny) / Double(totalPlanes)) - 0.5
                var r = min(maxRow, max(minRow, Int(idealRow.rounded())))
                if let last = planeRows.last, r <= last {
                    r = last + 1
                }
                planeRows.append(r)
            }
        } else {
            let candidates = Array(minRow...maxRow).shuffled()
            let minSpace = max(2, min_plane_spacing_rows)
            var chosen: [Int] = []
            for r in candidates {
                let tooClose = chosen.contains { existing in
                    let diff = abs(existing - r)
                    let pDiff = min(diff, max(1, ny - diff))
                    return pDiff < minSpace
                }
                if !tooClose {
                    chosen.append(r)
                    if chosen.count == totalPlanes { break }
                }
            }
            if chosen.count == totalPlanes {
                planeRows = chosen.sorted()
            } else {
                for k in 0..<totalPlanes {
                    let idealRow = (Double(k) + 0.5) * (Double(ny) / Double(totalPlanes)) - 0.5
                    var r = min(maxRow, max(minRow, Int(idealRow.rounded())))
                    if let last = planeRows.last, r <= last { r = last + 1 }
                    planeRows.append(r)
                }
            }
        }
        
        var signs = [Int](repeating: 1, count: totalPlanes)
        if dislocation_layout == DislocationLayoutPreset.randomBalanced.rawValue && !bound_strip_excess {
            for k in 0..<totalPlanes { signs[k] = (k < M) ? 1 : -1 }
            signs.shuffle()
        } else {
            for k in 0..<totalPlanes { signs[k] = (k % 2 == 0) ? 1 : -1 }
        }
        
        var newDislocs: [SeededDislocationModel] = []
        for k in 0..<totalPlanes {
            let xFrac: Double
            if dislocation_layout == DislocationLayoutPreset.alternatingStaggered.rawValue {
                xFrac = (k % 2 == 0) ? 0.35 : 0.65
            } else {
                xFrac = Double.random(in: 0.20...0.80)
            }
            newDislocs.append(SeededDislocationModel(
                id: k,
                glide_plane_row: planeRows[k],
                x_frac: (xFrac * 100).rounded() / 100.0,
                burgers_sign: signs[k],
                enabled: true
            ))
        }
        self.dislocations = newDislocs
    }
    
    /// Returns (grainID, thetaDeg) at normalized coordinate (xHat, yHat)
    public func grainAt(xHat: Double, yHat: Double, fallbackTheta: Double = 0.0) -> (Int, Double) {
        switch mode {
        case .singleCrystal, .periodicDislocations:
            return (0, fallbackTheta)
        case .bicrystal:
            let phi = gb_angle_deg * .pi / 180.0
            let nx = cos(phi)
            let ny = sin(phi)
            let dist = (xHat - gb_cx) * nx + (yHat - gb_cy) * ny
            return dist < 0.0 ? (0, theta1_deg) : (1, theta2_deg)
        case .voronoi:
            guard !seeds.isEmpty else { return (0, fallbackTheta) }
            var bestIdx = 0
            var bestD2 = Double.infinity
            for (idx, s) in seeds.enumerated() {
                let d2 = (xHat - s.cx) * (xHat - s.cx) + (yHat - s.cy) * (yHat - s.cy)
                if d2 < bestD2 {
                    bestD2 = d2
                    bestIdx = idx
                }
            }
            return (bestIdx, seeds[bestIdx].theta_deg)
        }
    }
}

public enum GeometryEngineeringPreset: String, CaseIterable, Identifiable {
    case perfectCrystalSlab = "perfect_crystal_slab"
    case plateWithHoleTension = "plate_hole_tension"
    case plateWithHoleShear = "plate_hole_shear"
    case dogboneTension = "dogbone_tension"
    case doubleNotchTension = "dent_tension"
    case preCrackedSlab = "pre_cracked_mode1"
    case tiltedEllipseTension = "tilted_ellipse_tension"
    case porousVoidArray = "porous_void_array"
    case nanoIndentation = "nano_indentation"
    case bicrystalTiltTension = "bicrystal_tilt_tension"
    case bicrystalIndentation = "bicrystal_indentation"
    case voronoiPolycrystalDogbone = "voronoi_polycrystal_dogbone"
    
    public var id: String { rawValue }
    public var title: String {
        switch self {
        case .perfectCrystalSlab: return "Perfect Crystal Slab (No Cutouts)"
        case .plateWithHoleTension: return "Plate with Hole (Uniaxial Tension)"
        case .plateWithHoleShear: return "Plate with Hole (Simple Shear)"
        case .dogboneTension: return "Dogbone Specimen (Uniaxial Tension)"
        case .doubleNotchTension: return "Double-Edge Notched Tension (DENT)"
        case .preCrackedSlab: return "Pre-Cracked Mode-I Fracture Slab"
        case .tiltedEllipseTension: return "Tilted 45° Elliptical Slot (Tension)"
        case .porousVoidArray: return "Porous 2×2 Void Metamaterial (Shear)"
        case .nanoIndentation: return "Nanoindentation on Crystal Substrate"
        case .bicrystalTiltTension: return "Symmetric Tilt Bicrystal (Uniaxial Tension)"
        case .bicrystalIndentation: return "Nanoindentation Near Bicrystal Grain Boundary"
        case .voronoiPolycrystalDogbone: return "2D Voronoi Polycrystal Dogbone (Tension)"
        }
    }
    public var subtitle: String {
        switch self {
        case .perfectCrystalSlab:
            return "Pristine rectangular single crystal without holes or notches; bottom clamped, top driven, left/right traction-free."
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
        case .bicrystalTiltTension:
            return "Vertical symmetric tilt grain boundary (θ₁=-15°, θ₂=+15°) under vertical tension with traction-free sides."
        case .bicrystalIndentation:
            return "Cylindrical indenter pressing near a horizontal/vertical tilt grain boundary to probe slip transmission."
        case .voronoiPolycrystalDogbone:
            return "Multi-grain 2D Voronoi polycrystal carved into a tensile dogbone specimen."
        }
    }
}

public struct MacroscopicGeometryConfigModel: Codable, Equatable {
    public var enabled: Bool = false
    public var preset_name: String = GeometryEngineeringPreset.perfectCrystalSlab.rawValue
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
    public var cutouts: [GeometryCutoutItem] = []
    public var boundary_rules: [BoundaryZoneRuleModel] = [
        BoundaryZoneRuleModel(region: .bottom, kind: .fixed, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
        BoundaryZoneRuleModel(region: .top, kind: .prescribed, thickness_a0: 2.0, ux: 0.0, uy: 1.0),
        BoundaryZoneRuleModel(region: .left, kind: .free, thickness_a0: 1.5, ux: 0.0, uy: 0.0),
        BoundaryZoneRuleModel(region: .right, kind: .free, thickness_a0: 1.5, ux: 0.0, uy: 0.0)
    ]
    public var stage1_boundary_rules: [BoundaryZoneRuleModel] = [
        BoundaryZoneRuleModel(region: .top, kind: .rollerX, thickness_a0: 1.5, ux: 0.0, uy: 0.0, lock_x: false, lock_y: true),
        BoundaryZoneRuleModel(region: .bottom, kind: .rollerX, thickness_a0: 1.5, ux: 0.0, uy: 0.0, lock_x: false, lock_y: true),
        BoundaryZoneRuleModel(region: .left, kind: .rollerY, thickness_a0: 1.5, ux: 0.0, uy: 0.0, lock_x: true, lock_y: false),
        BoundaryZoneRuleModel(region: .right, kind: .free, thickness_a0: 1.5, ux: 0.0, uy: 0.0)
    ]
    public var indenter: IndenterConfigModel = IndenterConfigModel()
    public var microstructure: MicrostructureConfigModel = MicrostructureConfigModel()
    public var primary_stress_component: String = "auto"
    
    public init() {}
    
    public mutating func applyPreset(_ preset: GeometryEngineeringPreset) {
        self.enabled = true
        self.preset_name = preset.rawValue
        self.indenter.enabled = false
        self.microstructure.mode = .singleCrystal
        switch preset {
        case .perfectCrystalSlab:
            outer_shape = .rectangle
            cutouts = []
            boundary_rules = [
                BoundaryZoneRuleModel(region: .bottom, kind: .fixed, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                BoundaryZoneRuleModel(region: .top, kind: .prescribed, thickness_a0: 2.0, ux: 0.0, uy: 1.0),
                BoundaryZoneRuleModel(region: .left, kind: .free, thickness_a0: 1.5),
                BoundaryZoneRuleModel(region: .right, kind: .free, thickness_a0: 1.5)
            ]
            primary_stress_component = "sigma_yy"
            
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
            
        case .bicrystalTiltTension:
            outer_shape = .rectangle
            cutouts = []
            boundary_rules = [
                BoundaryZoneRuleModel(region: .bottom, kind: .fixed, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                BoundaryZoneRuleModel(region: .top, kind: .prescribed, thickness_a0: 2.0, ux: 0.0, uy: 1.0),
                BoundaryZoneRuleModel(region: .left, kind: .free, thickness_a0: 1.5),
                BoundaryZoneRuleModel(region: .right, kind: .free, thickness_a0: 1.5)
            ]
            microstructure.mode = .bicrystal
            microstructure.gb_cx = 0.50
            microstructure.gb_cy = 0.50
            microstructure.gb_angle_deg = 0.0
            microstructure.theta1_deg = -15.0
            microstructure.theta2_deg = +15.0
            primary_stress_component = "sigma_yy"
            
        case .bicrystalIndentation:
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
            microstructure.mode = .bicrystal
            microstructure.gb_cx = 0.50
            microstructure.gb_cy = 0.68
            microstructure.gb_angle_deg = 90.0
            microstructure.theta1_deg = -12.0
            microstructure.theta2_deg = +12.0
            primary_stress_component = "indentation_load"
            
        case .voronoiPolycrystalDogbone:
            outer_shape = .dogbone
            dogbone_gauge_width = 0.52
            dogbone_gauge_length = 0.44
            dogbone_fillet_radius = 0.12
            cutouts = []
            boundary_rules = [
                BoundaryZoneRuleModel(region: .bottom, kind: .fixed, thickness_a0: 2.5, ux: 0.0, uy: 0.0),
                BoundaryZoneRuleModel(region: .top, kind: .prescribed, thickness_a0: 2.5, ux: 0.0, uy: 1.0)
            ]
            microstructure.mode = .voronoi
            microstructure.num_grains = 6
            microstructure.theta_min_deg = -25.0
            microstructure.theta_max_deg = +25.0
            microstructure.regenerateSeeds()
            primary_stress_component = "sigma_yy"
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
                if margin > 0 && Self.distanceToPolygonEdges(x: xHat, y: yHat, poly: c.vertices) <= margin {
                    return true
                }
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
    
    /// Highest solid yHat at the given xHat (used for positioning nanoindenter on carved top surface)
    public func topSurfaceYHat(at xHat: Double) -> Double {
        let steps = 200
        let searchWidths = [0.0, 0.015, 0.03]
        for halfW in searchWidths {
            var bestY: Double? = nil
            let xSamples = halfW == 0.0 ? [xHat] : [max(0.0, xHat - halfW), xHat, min(1.0, xHat + halfW)]
            for xs in xSamples {
                for k in 0...steps {
                    let yh = 1.0 - Double(k) / Double(steps)
                    if isPointInsideSolid(xHat: xs, yHat: yh) {
                        if bestY == nil || yh > bestY! {
                            bestY = yh
                        }
                        break
                    }
                }
            }
            if let found = bestY {
                return found
            }
        }
        return 1.0
    }
    
    /// Returns the Stage 1 (α = 0) relaxation boundary rules (user-editable stage1_boundary_rules)
    public func stage1RelaxationRules() -> [BoundaryZoneRuleModel] {
        if !stage1_boundary_rules.isEmpty {
            return stage1_boundary_rules
        }
        return Self.defaultStage1RollerRules(forAngleDeg: microstructure.gb_angle_deg)
    }

    public static func defaultStage1RollerRules(forAngleDeg gbAngleDeg: Double) -> [BoundaryZoneRuleModel] {
        let phiRad = gbAngleDeg * .pi / 180.0
        let isVerticalGB = abs(cos(phiRad)) >= abs(sin(phiRad))
        if isVerticalGB {
            return [
                BoundaryZoneRuleModel(region: .top, kind: .rollerX, thickness_a0: 1.5, ux: 0.0, uy: 0.0, lock_x: false, lock_y: true),
                BoundaryZoneRuleModel(region: .bottom, kind: .rollerX, thickness_a0: 1.5, ux: 0.0, uy: 0.0, lock_x: false, lock_y: true),
                BoundaryZoneRuleModel(region: .left, kind: .rollerY, thickness_a0: 1.5, ux: 0.0, uy: 0.0, lock_x: true, lock_y: false),
                BoundaryZoneRuleModel(region: .right, kind: .free, thickness_a0: 1.5)
            ]
        } else {
            return [
                BoundaryZoneRuleModel(region: .left, kind: .rollerY, thickness_a0: 1.5, ux: 0.0, uy: 0.0, lock_x: true, lock_y: false),
                BoundaryZoneRuleModel(region: .right, kind: .rollerY, thickness_a0: 1.5, ux: 0.0, uy: 0.0, lock_x: true, lock_y: false),
                BoundaryZoneRuleModel(region: .bottom, kind: .rollerX, thickness_a0: 1.5, ux: 0.0, uy: 0.0, lock_x: false, lock_y: true),
                BoundaryZoneRuleModel(region: .top, kind: .free, thickness_a0: 1.5)
            ]
        }
    }

    public mutating func resetStage1ToGBFrameRollerBCs(forAngleDeg gbAngleDeg: Double) {
        stage1_boundary_rules = Self.defaultStage1RollerRules(forAngleDeg: gbAngleDeg)
    }

    /// Returns 0 = Free, 1 = Fixed, 2 = Driven, 3 = RollerY (X-locked), 4 = RollerX (Y-locked)
    public func classifyNodeRole(
        xHat: Double,
        yHat: Double,
        nx: Int,
        ny: Int,
        overrideRules: [BoundaryZoneRuleModel]? = nil,
        overrideIndenterEnabled: Bool? = nil
    ) -> Int {
        var lockedX = false
        var lockedY = false
        var drivenX = false
        var drivenY = false
        let a0x = 1.0 / Double(max(1, nx - 1))
        let a0y = 1.0 / Double(max(1, ny - 1))
        let activeRules = overrideRules ?? boundary_rules
        let indenterActive = overrideIndenterEnabled ?? indenter.enabled
        
        var domMinX = 0.0, domMaxX = 1.0, domMinY = 0.0, domMaxY = 1.0
        if outer_shape == .polygon && outer_polygon.count >= 3 {
            domMinX = outer_polygon.map(\.x).min() ?? 0.0
            domMaxX = outer_polygon.map(\.x).max() ?? 1.0
            domMinY = outer_polygon.map(\.y).min() ?? 0.0
            domMaxY = outer_polygon.map(\.y).max() ?? 1.0
        } else if outer_shape == .disk {
            domMinX = max(0.0, 0.5 - disk_radius)
            domMaxX = min(1.0, 0.5 + disk_radius)
            domMinY = max(0.0, 0.5 - disk_radius)
            domMaxY = min(1.0, 0.5 + disk_radius)
        }
        
        for rule in activeRules where rule.enabled {
            let stripX = max(0.35, rule.thickness_a0) * a0x
            let stripY = max(0.35, rule.thickness_a0) * a0y
            var matches = false
            switch rule.region {
            case .top:
                if indenterActive && rule.kind != .free {
                    matches = false
                } else {
                    matches = (yHat >= domMaxY - stripY)
                }
            case .bottom: matches = (yHat <= domMinY + stripY)
            case .left: matches = (xHat <= domMinX + stripX)
            case .right: matches = (xHat >= domMaxX - stripX)
            case .outerPerimeter:
                if outer_shape == .disk {
                    let dist = hypot(xHat - 0.5, yHat - 0.5)
                    matches = (dist >= disk_radius - max(stripX, stripY))
                } else if outer_shape == .polygon && outer_polygon.count >= 3 {
                    let dist = Self.distanceToPolygonEdges(x: xHat, y: yHat, poly: outer_polygon)
                    matches = (dist <= max(stripX, stripY))
                } else if outer_shape == .dogbone {
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
                    matches = (yHat <= stripY || yHat >= 1.0 - stripY || (halfWAllowed - abs(xHat - 0.5)) <= stripX)
                } else {
                    matches = (xHat <= domMinX + stripX || xHat >= domMaxX - stripX ||
                               yHat <= domMinY + stripY || yHat >= domMaxY - stripY)
                }
            case .holeSurfaces:
                matches = isPointInsideCutout(xHat: xHat, yHat: yHat, margin: max(stripX, stripY))
            case .customBox:
                matches = (xHat >= rule.box_x_min && xHat <= rule.box_x_max &&
                           yHat >= rule.box_y_min && yHat <= rule.box_y_max)
            }
            if matches {
                switch rule.kind {
                case .free:
                    break
                case .rollerX:
                    lockedY = true
                    drivenY = false
                case .rollerY:
                    lockedX = true
                    drivenX = false
                case .fixed:
                    if rule.lock_x { lockedX = true; drivenX = false }
                    if rule.lock_y { lockedY = true; drivenY = false }
                case .prescribed:
                    if rule.lock_x { lockedX = true; drivenX = true }
                    if rule.lock_y { lockedY = true; drivenY = true }
                case .affine:
                    lockedX = true; lockedY = true
                    drivenX = true; drivenY = true
                }
            }
        }
        var role = 0
        if lockedX && lockedY {
            role = (drivenX || drivenY) ? 2 : 1
        } else if lockedX && !lockedY {
            role = 3
        } else if !lockedX && lockedY {
            role = 4
        } else {
            role = 0
        }
        if indenterActive && role == 0 {
            let rNorm = indenter.radius_a0 * a0x
            let cx = indenter.cx
            let yTop = topSurfaceYHat(at: cx)
            let dx = abs(xHat - cx)
            switch indenter.shape {
            case .cylinder:
                let cy = yTop + rNorm
                let dy = yHat - cy
                let dist = hypot(xHat - cx, dy)
                if dist <= rNorm + 1.25 * a0y && yHat >= yTop - 2.5 * a0y {
                    role = 2
                }
            case .flatPunch:
                if dx <= rNorm + 0.75 * a0x && yHat >= yTop - 1.75 * a0y {
                    role = 2
                }
            case .wedge:
                let theta = min(85.0, max(15.0, indenter.wedge_half_angle_deg)) * .pi / 180.0
                let halfSpan = max(2.5 * a0x, min(rNorm, 2.5 * a0y * tan(theta)))
                if dx <= halfSpan && yHat >= yTop - 2.0 * a0y {
                    role = 2
                }
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
                (x < (xj - xi) * (y - yi) / (yj - yi) + xi)
            if intersect { inside.toggle() }
            j = i
        }
        return inside
    }
    
    public static func distanceToPolygonEdges(x: Double, y: Double, poly: [NormalizedVertex]) -> Double {
        guard poly.count >= 2 else { return .infinity }
        var minDist = Double.infinity
        var j = poly.count - 1
        for i in 0..<poly.count {
            let ax = poly[j].x, ay = poly[j].y
            let bx = poly[i].x, by = poly[i].y
            let abx = bx - ax, aby = by - ay
            let len2 = abx * abx + aby * aby
            let t = len2 > 1e-15 ? min(1.0, max(0.0, ((x - ax) * abx + (y - ay) * aby) / len2)) : 0.0
            let px = ax + t * abx, py = ay + t * aby
            let d = hypot(x - px, y - py)
            if d < minDist { minDist = d }
            j = i
        }
        return minDist
    }
    
    public mutating func applyGBFrameRollerBCs(forAngleDeg gbAngleDeg: Double) {
        let phiRad = gbAngleDeg * .pi / 180.0
        let isMostlyVertical = abs(cos(phiRad)) >= abs(sin(phiRad))
        let rollerKind: BoundaryConditionKindModel = isMostlyVertical ? .rollerX : .rollerY
        let lockX = !isMostlyVertical
        let lockY = isMostlyVertical
        boundary_rules = [
            BoundaryZoneRuleModel(region: .bottom, kind: rollerKind, thickness_a0: 2.0, ux: 0.0, uy: 0.0, lock_x: lockX, lock_y: lockY),
            BoundaryZoneRuleModel(region: .top, kind: rollerKind, thickness_a0: 2.0, ux: 0.0, uy: 0.0, lock_x: lockX, lock_y: lockY),
            BoundaryZoneRuleModel(region: .left, kind: rollerKind, thickness_a0: 2.0, ux: 0.0, uy: 0.0, lock_x: lockX, lock_y: lockY),
            BoundaryZoneRuleModel(region: .right, kind: rollerKind, thickness_a0: 2.0, ux: 0.0, uy: 0.0, lock_x: lockX, lock_y: lockY)
        ]
        primary_stress_component = isMostlyVertical ? "sigma_xx" : "sigma_yy"
    }

    public mutating func rotateBicrystalOrientation(toAngleDeg newPhiDeg: Double) {
        let oldPhi = microstructure.gb_angle_deg
        let rel1 = microstructure.theta1_deg - oldPhi
        let rel2 = microstructure.theta2_deg - oldPhi
        microstructure.gb_angle_deg = newPhiDeg
        microstructure.theta1_deg = newPhiDeg + rel1
        microstructure.theta2_deg = newPhiDeg + rel2
        resetStage1ToGBFrameRollerBCs(forAngleDeg: newPhiDeg)
        let allRollers = !boundary_rules.isEmpty && boundary_rules.allSatisfy { $0.kind == .rollerX || $0.kind == .rollerY }
        if allRollers {
            applyGBFrameRollerBCs(forAngleDeg: newPhiDeg)
        }
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

// MARK: - General (m, n) Coincidence Site Lattice (CSL) Rationalizer (gb_code / aimsgb)
public struct CSLRationalMatch: Equatable, Identifiable {
    public let m: Int
    public let n: Int
    public let sigma: Int
    public let halfAngleDeg: Double
    public let fullMisoriDeg: Double
    public let periodA0: Double
    public let layerSpacingA0: Double
    public let recommendedDMinA0: Double
    public let recommendedGapA0: Double

    public var id: String { String(format: "%.3f_%d_%d", fullMisoriDeg, m, n) }

    public var shortLabel: String {
        if sigma > 0 {
            return String(format: "Σ%d (%d %d 0) • ±%.2f° (Δθ=%.2f°)", sigma, m, n, halfAngleDeg, fullMisoriDeg)
        }
        return String(format: "Δθ = %.2f° (θ = ±%.2f°)", fullMisoriDeg, halfAngleDeg)
    }

    public func commensurateSize(periods k: Int) -> Int {
        max(16, Int(round(Double(max(1, k)) * periodA0)))
    }

    public func suggestedPeriodCounts(around targetN: Int) -> [Int] {
        let p = max(1.0, periodA0)
        let kRaw = max(2, Int(round(Double(targetN) / p)))
        // Quantize the 3-button window into fixed blocks of 6 periods
        // (e.g. [43, 45, 47] for kRaw in 43...48) so clicking any of the 3
        // commensurate buttons never shifts the button list under the cursor.
        let block = max(0, (kRaw - 1) / 6)
        let kCenter = max(3, block * 6 + 3)
        var counts: [Int] = []
        for dk in [-2, 0, 2] {
            let k = max(1, kCenter + dk)
            if !counts.contains(k) {
                counts.append(k)
            }
        }
        return counts
    }
}

public enum CSLRationalizer {
    public static func gcd(_ a: Int, _ b: Int) -> Int {
        var x = abs(a)
        var y = abs(b)
        while y != 0 {
            let r = x % y
            x = y
            y = r
        }
        return max(1, x)
    }

    public static func makeMatch(m: Int, n: Int) -> CSLRationalMatch {
        if n <= 0 {
            return CSLRationalMatch(
                m: 1, n: 0, sigma: 1,
                halfAngleDeg: 0.0, fullMisoriDeg: 0.0,
                periodA0: 1.0, layerSpacingA0: 1.0,
                recommendedDMinA0: 0.48, recommendedGapA0: 0.0
            )
        }
        if m == 1 && n == 1 {
            let s2 = sqrt(2.0)
            return CSLRationalMatch(
                m: 1, n: 1, sigma: 1,
                halfAngleDeg: 45.0, fullMisoriDeg: 90.0,
                periodA0: s2, layerSpacingA0: 1.0 / s2,
                recommendedDMinA0: 0.48, recommendedGapA0: 0.0
            )
        }
        let sumSq = m * m + n * n
        let g2 = gcd(sumSq, 2)
        let sigma = sumSq / g2
        let halfDeg = atan2(Double(n), Double(m)) * 180.0 / .pi
        let period = sqrt(Double(sumSq))
        let dLayer = 1.0 / period
        let recDMin: Double
        if m == 3 && n == 1 {
            recDMin = 0.48
        } else if m == 2 && n == 1 {
            recDMin = 0.50
        } else {
            recDMin = min(0.60, max(0.48, (round((1.12 * dLayer) * 100.0) / 100.0)))
        }
        let recGap: Double
        if m == 3 && n == 1 {
            recGap = 0.26
        } else if m == 2 && n == 1 {
            recGap = 0.0
        } else if m - n == 1 && m >= 4 {
            // High-angle (m, m-1) branch approaching (110) single crystal at Δθ = 90°
            recGap = min(0.16, max(0.0, (round((0.68 / Double(m)) * 50.0) / 50.0)))
        } else {
            let kFirst = max(1, Int(ceil((recDMin - 1e-6) / (2.0 * dLayer))))
            let rFirst = 2.0 * Double(kFirst) * dLayer
            let rawGap = max(0.0, 0.94 - rFirst)
            recGap = min(0.38, max(0.0, (round(rawGap * 50.0) / 50.0)))
        }
        return CSLRationalMatch(
            m: m,
            n: n,
            sigma: sigma,
            halfAngleDeg: halfDeg,
            fullMisoriDeg: 2.0 * halfDeg,
            periodA0: period,
            layerSpacingA0: dLayer,
            recommendedDMinA0: recDMin,
            recommendedGapA0: recGap
        )
    }

    public static let allCandidates: [CSLRationalMatch] = {
        var list: [CSLRationalMatch] = [makeMatch(m: 1, n: 0)]
        for m in 1...18 {
            for n in 1...m {
                if gcd(m, n) == 1 {
                    list.append(makeMatch(m: m, n: n))
                }
            }
        }
        list.sort { $0.halfAngleDeg < $1.halfAngleDeg }
        return list
    }()

    /// Curated CSL sequence spanning the full [0°, 90°] misorientation arch on a 2D Square Lattice
    public static let curatedCatalog: [CSLRationalMatch] = {
        let pairs: [(Int, Int)] = [
            (1, 0),   // ±0.00°  (Δθ = 0.00°,  Σ1 (100) Single Crystal)
            (12, 1),  // ±4.76°  (Δθ = 9.53°,  Σ145)
            (10, 1),  // ±5.71°  (Δθ = 11.42°, Σ101)
            (9, 1),   // ±6.34°  (Δθ = 12.68°, Σ41)
            (8, 1),   // ±7.13°  (Δθ = 14.25°, Σ65)
            (7, 1),   // ±8.13°  (Δθ = 16.26°, Σ25)
            (6, 1),   // ±9.46°  (Δθ = 18.92°, Σ37)
            (5, 1),   // ±11.31° (Δθ = 22.62°, Σ13)
            (4, 1),   // ±14.04° (Δθ = 28.07°, Σ17)
            (3, 1),   // ±18.43° (Δθ = 36.87°, Σ5 (310) Cusp)
            (5, 2),   // ±21.80° (Δθ = 43.60°, Σ29)
            (2, 1),   // ±26.57° (Δθ = 53.13°, Σ5 (210) Deep Cusp)
            (5, 3),   // ±30.96° (Δθ = 61.93°, Σ17)
            (3, 2),   // ±33.69° (Δθ = 67.38°, Σ13)
            (4, 3),   // ±36.87° (Δθ = 73.74°, Σ25)
            (5, 4),   // ±38.66° (Δθ = 77.32°, Σ41)
            (6, 5),   // ±39.81° (Δθ = 79.61°, Σ61)
            (8, 7),   // ±41.19° (Δθ = 82.37°, Σ113)
            (10, 9),  // ±41.99° (Δθ = 83.97°, Σ181)
            (12, 11), // ±42.51° (Δθ = 85.01°, Σ265)
            (1, 1)    // ±45.00° (Δθ = 90.00°, Σ1 (110) Single Crystal)
        ]
        return pairs.map { makeMatch(m: $0.0, n: $0.1) }
    }()

    public static func closestMatch(forHalfAngleDeg targetDeg: Double, maxIndex: Int = 18) -> CSLRationalMatch {
        let absTarget = max(0.0, min(45.0, abs(targetDeg)))
        if absTarget < 0.25 {
            return makeMatch(m: 1, n: 0)
        }
        var best = makeMatch(m: 3, n: 1)
        var bestErr = Double.infinity
        for cand in allCandidates where cand.m <= maxIndex && cand.n > 0 {
            // Slightly favor lower-index (higher-coincidence) boundaries when within 0.15 deg
            let angleErr = abs(cand.halfAngleDeg - absTarget)
            let complexityPenalty = 0.008 * Double(cand.m)
            let score = angleErr + complexityPenalty
            if score < bestErr {
                bestErr = score
                best = cand
            }
        }
        return best
    }

    public static func makeContinuousCandidate(
        fullMisoriDeg: Double,
        isTriangular: Bool,
        sigma: Int = 0,
        m: Int = 0,
        n: Int = 0,
        periodA0: Double = 1.0
    ) -> CSLRationalMatch {
        let symDeg = isTriangular ? 60.0 : 90.0
        let clampedDTheta = min(symDeg, max(0.0, fullMisoriDeg))
        let halfDeg = 0.5 * clampedDTheta
        let isEndpoint = (clampedDTheta <= 0.05 || clampedDTheta >= symDeg - 0.05)
        let dMin = isTriangular ? 0.52 : 0.50
        let gap: Double = isEndpoint
            ? 0.0
            : (round((0.08 * sin(.pi * clampedDTheta / symDeg)) * 100.0) / 100.0)
        let effPeriod: Double = (periodA0 > 1.01 || isEndpoint)
            ? max(1.0, periodA0)
            : closestMatch(forHalfAngleDeg: halfDeg).periodA0
        return CSLRationalMatch(
            m: m,
            n: n,
            sigma: sigma,
            halfAngleDeg: halfDeg,
            fullMisoriDeg: clampedDTheta,
            periodA0: effPeriod,
            layerSpacingA0: 1.0 / max(1.0, effPeriod),
            recommendedDMinA0: dMin,
            recommendedGapA0: gap
        )
    }

    public static func buildContinuousSweep(
        thetaMinDeg: Double,
        thetaMaxDeg: Double,
        stepDeltaThetaDeg: Double,
        isTriangular: Bool,
        includeExactCusps: Bool
    ) -> [CSLRationalMatch] {
        let symDeg = isTriangular ? 60.0 : 90.0
        let dThetaMin = min(symDeg, max(0.0, 2.0 * min(thetaMinDeg, thetaMaxDeg)))
        let dThetaMax = min(symDeg, max(0.0, 2.0 * max(thetaMinDeg, thetaMaxDeg)))
        let step = max(0.5, min(30.0, stepDeltaThetaDeg))

        // Major cusps for Square (0-90°) vs Triangular (0-60°)
        let exactCusps: [(dtheta: Double, sigma: Int, m: Int, n: Int, per: Double)] = isTriangular
            ? [
                (0.0, 1, 1, 0, 1.0),
                (21.786789, 7, 2, 1, sqrt(7.0)),
                (27.795772, 13, 3, 1, sqrt(13.0)),
                (32.204228, 13, 3, 2, sqrt(13.0)),
                (38.213211, 7, 3, 2, sqrt(7.0)),
                (60.0, 1, 1, 1, 1.0)
            ]
            : [
                (0.0, 1, 1, 0, 1.0),
                (22.619865, 13, 5, 1, sqrt(26.0)),
                (28.072487, 17, 4, 1, sqrt(17.0)),
                (36.869898, 5, 3, 1, sqrt(10.0)),
                (53.130102, 5, 2, 1, sqrt(5.0)),
                (67.380135, 13, 3, 2, sqrt(13.0)),
                (90.0, 1, 1, 1, sqrt(2.0))
            ]

        var angles: [Double] = []
        var cur = dThetaMin
        while cur <= dThetaMax + 1e-6 {
            angles.append(cur)
            cur += step
        }
        if let last = angles.last, abs(last - dThetaMax) > 0.15 {
            angles.append(dThetaMax)
        }

        if includeExactCusps {
            for cusp in exactCusps where cusp.dtheta >= dThetaMin - 1e-4 && cusp.dtheta <= dThetaMax + 1e-4 {
                // Replace any grid point within 0.45 * step of an exact cusp so the cusp isn't duplicated
                angles.removeAll { abs($0 - cusp.dtheta) < min(0.85, 0.42 * step) }
                angles.append(cusp.dtheta)
            }
        }
        angles.sort()

        return angles.map { ang in
            if let cusp = exactCusps.first(where: { abs($0.dtheta - ang) < 0.05 }) {
                return makeContinuousCandidate(
                    fullMisoriDeg: cusp.dtheta,
                    isTriangular: isTriangular,
                    sigma: cusp.sigma,
                    m: cusp.m,
                    n: cusp.n,
                    periodA0: cusp.per
                )
            }
            return makeContinuousCandidate(fullMisoriDeg: ang, isTriangular: isTriangular)
        }
    }
}
