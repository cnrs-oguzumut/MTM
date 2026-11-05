import SwiftUI
import AppKit
import simd

// MARK: - Colormaps

public enum MeshColormap: String, CaseIterable, Identifiable {
    case coolwarm = "Cool–Warm"
    case viridis = "Viridis"
    case plasma = "Plasma"
    case inferno = "Inferno"
    case turbo = "Turbo"
    case grayscale = "Grayscale"
    
    public var id: String { rawValue }
    
    private var stops: [(Float, UInt32)] {
        switch self {
        case .coolwarm:
            return [(0, 0x3B4CC0), (0.25, 0x7B9FF9), (0.5, 0xDDDDDD), (0.75, 0xF49A7B), (1, 0xB40426)]
        case .viridis:
            return [(0, 0x440154), (0.1, 0x482475), (0.2, 0x414487), (0.3, 0x355F8D), (0.4, 0x2A788E),
                    (0.5, 0x21918C), (0.6, 0x22A884), (0.7, 0x44BF70), (0.8, 0x7AD151), (0.9, 0xBDDF26), (1, 0xFDE725)]
        case .plasma:
            return [(0, 0x0D0887), (0.1, 0x41049D), (0.2, 0x6A00A8), (0.3, 0x8F0DA4), (0.4, 0xB12A90),
                    (0.5, 0xCC4778), (0.6, 0xE16462), (0.7, 0xF2844B), (0.8, 0xFCA636), (0.9, 0xFCCE25), (1, 0xF0F921)]
        case .inferno:
            return [(0, 0x000004), (0.1, 0x160B39), (0.2, 0x420A68), (0.3, 0x6A176E), (0.4, 0x932667),
                    (0.5, 0xBC3754), (0.6, 0xDD513A), (0.7, 0xF37819), (0.8, 0xFCA50A), (0.9, 0xF6D746), (1, 0xFCFFA4)]
        case .turbo:
            return [(0, 0x30123B), (0.13, 0x4662D7), (0.25, 0x36AAF9), (0.38, 0x1AE4B6), (0.5, 0x72FE5E),
                    (0.63, 0xC8EF34), (0.75, 0xFABA39), (0.88, 0xF66B19), (1, 0x7A0403)]
        case .grayscale:
            return [(0, 0x000000), (1, 0xFFFFFF)]
        }
    }
    
    private static func rgb(_ hex: UInt32) -> SIMD3<Float> {
        SIMD3(Float((hex >> 16) & 0xFF), Float((hex >> 8) & 0xFF), Float(hex & 0xFF)) / 255
    }
    
    public func sample(_ t: Float) -> SIMD3<Float> {
        let u = min(1, max(0, t.isFinite ? t : 0))
        let s = stops
        for k in 1..<s.count where u <= s[k].0 {
            let (t0, c0) = s[k - 1], (t1, c1) = s[k]
            let f = (u - t0) / max(1e-6, t1 - t0)
            return MeshColormap.rgb(c0) * (1 - f) + MeshColormap.rgb(c1) * f
        }
        return MeshColormap.rgb(s.last!.1)
    }
    
    /// 256-entry RGBA8 lookup table for the GPU.
    func lut() -> [UInt8] {
        var out = [UInt8](repeating: 255, count: 256 * 4)
        for i in 0..<256 {
            let c = sample(Float(i) / 255)
            out[4 * i] = UInt8(min(255, max(0, c.x * 255)))
            out[4 * i + 1] = UInt8(min(255, max(0, c.y * 255)))
            out[4 * i + 2] = UInt8(min(255, max(0, c.z * 255)))
        }
        return out
    }
    
    public func gradientColors(_ n: Int = 12) -> [Color] {
        (0...n).map { i in
            let c = sample(Float(i) / Float(n))
            return Color(red: Double(c.x), green: Double(c.y), blue: Double(c.z))
        }
    }
}

// MARK: - Geometry passed to the renderer

/// Triangle "soup" (3 vertices per triangle) plus optional node and outline layers, in world coordinates.
public struct MeshGeometry {
    public let id = UUID()
    public var triPositions: [SIMD2<Float>]
    public var nodePositions: [SIMD2<Float>] = []
    /// Optional per-node colours (same count as `nodePositions`).
    public var nodeColors: [SIMD4<Float>]? = nil
    /// Line segments as point pairs (domain borders etc.).
    public var outlineSegments: [SIMD2<Float>] = []
    
    public var triangleCount: Int { triPositions.count / 3 }
    public var nodeCount: Int { nodePositions.count }
    
    public var bounds: (min: SIMD2<Float>, max: SIMD2<Float>) {
        let all = triPositions.isEmpty ? (nodePositions.isEmpty ? outlineSegments : nodePositions) : triPositions
        guard var lo = all.first else { return (.zero, SIMD2(1, 1)) }
        var hi = lo
        for p in all { lo = simd_min(lo, p); hi = simd_max(hi, p) }
        for p in outlineSegments { lo = simd_min(lo, p); hi = simd_max(hi, p) }
        return (lo, hi)
    }
    
    public init(triPositions: [SIMD2<Float>], nodePositions: [SIMD2<Float>] = [], nodeColors: [SIMD4<Float>]? = nil, outlineSegments: [SIMD2<Float>] = []) {
        self.triPositions = triPositions
        self.nodePositions = nodePositions
        self.nodeColors = nodeColors
        self.outlineSegments = outlineSegments
    }
    
    public init(dataset: VTKDataset) {
        var tri = [SIMD2<Float>]()
        tri.reserveCapacity(dataset.triangles.count)
        for i in dataset.triangles { tri.append(dataset.points[Int(i)]) }
        self.init(triPositions: tri, nodePositions: dataset.points)
    }
}

/// One value per triangle vertex (`triPositions.count` entries). `version` must change whenever the values do.
public struct MeshValues {
    public let version: Int
    public let perVertex: [Float]
}

public struct MeshAppearance: Equatable {
    public var colormap: MeshColormap = .coolwarm
    public var rangeMin: Float = 0
    public var rangeMax: Float = 1
    public var showFill = true
    public var showEdges = false
    public var showNodes = false
    public var showOutline = true
    public var solidFill = SIMD4<Float>(0.31, 0.55, 0.85, 0.30)
    public var edgeColor = SIMD4<Float>(0, 0, 0, 0.45)
    public var nodeColor = SIMD4<Float>(0.35, 0.55, 1.0, 1.0)
    public var outlineColor = SIMD4<Float>(0.40, 0.75, 1.0, 0.9)
    public var nodeSizePx: Float = 4
    public var background = SIMD4<Float>(0.043, 0.059, 0.094, 1)
    
    public init() {}
}

// MARK: - Spatial lookup for hover

/// Uniform-grid bucket index over the triangles, for O(1) "which triangle is under the cursor".
public final class TriangleLocator {
    private let tri: [SIMD2<Float>]
    private let lo: SIMD2<Float>
    private let cell: SIMD2<Float>
    private let nx: Int, ny: Int
    private var buckets: [[Int32]]
    
    public init(triPositions: [SIMD2<Float>]) {
        tri = triPositions
        let nT = triPositions.count / 3
        var mn = triPositions.first ?? .zero, mx = mn
        for p in triPositions { mn = simd_min(mn, p); mx = simd_max(mx, p) }
        lo = mn
        let side = max(1, Int(Double(max(1, nT)).squareRoot() / 1.5))
        let ext = simd_max(mx - mn, SIMD2(1e-6, 1e-6))
        let aspect = ext.x / ext.y
        nx = max(1, Int(Float(side) * aspect.squareRoot()))
        ny = max(1, Int(Float(side) / aspect.squareRoot()))
        cell = ext / SIMD2(Float(nx), Float(ny))
        buckets = Array(repeating: [], count: nx * ny)
        for t in 0..<nT {
            let a = tri[3 * t], b = tri[3 * t + 1], c = tri[3 * t + 2]
            let bmin = simd_min(a, simd_min(b, c)), bmax = simd_max(a, simd_max(b, c))
            let i0 = cellX(bmin.x), i1 = cellX(bmax.x), j0 = cellY(bmin.y), j1 = cellY(bmax.y)
            if (i1 - i0 + 1) * (j1 - j0 + 1) > 64 { continue } // degenerate / periodic-wrap triangle
            for j in j0...j1 { for i in i0...i1 { buckets[j * nx + i].append(Int32(t)) } }
        }
    }
    
    @inline(__always) private func cellX(_ x: Float) -> Int { min(nx - 1, max(0, Int((x - lo.x) / cell.x))) }
    @inline(__always) private func cellY(_ y: Float) -> Int { min(ny - 1, max(0, Int((y - lo.y) / cell.y))) }
    
    public func triangle(at p: SIMD2<Float>) -> Int? {
        for t in buckets[cellY(p.y) * nx + cellX(p.x)] {
            let a = tri[3 * Int(t)], b = tri[3 * Int(t) + 1], c = tri[3 * Int(t) + 2]
            let d1 = cross2(p - a, b - a), d2 = cross2(p - b, c - b), d3 = cross2(p - c, a - c)
            let hasNeg = d1 < 0 || d2 < 0 || d3 < 0, hasPos = d1 > 0 || d2 > 0 || d3 > 0
            if !(hasNeg && hasPos) { return Int(t) }
        }
        return nil
    }
    
    @inline(__always) private func cross2(_ u: SIMD2<Float>, _ v: SIMD2<Float>) -> Float { u.x * v.y - u.y * v.x }
}
