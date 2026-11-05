import SwiftUI
import simd

/// Builds and caches GPU geometry for the FEM Mesh preview, so it is rebuilt only when
/// the mesh, the 1×1 / 3×3 mode, or the colour field changes – never on pan / zoom.
final class FEMPreviewGPUCache {
    private var geometryKey = ""
    private(set) var geometry: MeshGeometry?
    private(set) var locator: TriangleLocator?
    private var valuesKey = ""
    private(set) var values: MeshValues?
    private var valuesVersion = 0
    
    static func usesTiling(_ mesh: InitialMeshModel, periodic: Bool) -> Bool {
        periodic && !mesh.elements_3x3.isEmpty
    }
    
    static func triangleCount(_ mesh: InitialMeshModel, periodic: Bool) -> Int {
        usesTiling(mesh, periodic: periodic) ? mesh.elements_3x3.count : mesh.elements.count
    }
    
    static func nodeCount(_ mesh: InitialMeshModel, periodic: Bool) -> Int {
        periodic ? mesh.translated_domains.reduce(0) { $0 + $1.nodes.count } : mesh.primary_nodes.count
    }
    
    func geometry(for mesh: InitialMeshModel, version: Int, periodic: Bool) -> MeshGeometry {
        let key = "\(version)|\(periodic)"
        if key == geometryKey, let g = geometry { return g }
        
        var tri: [SIMD2<Float>] = []
        func add(_ c: [[Double]]) {
            guard c.count >= 3, c[0].count >= 2, c[1].count >= 2, c[2].count >= 2 else {
                tri.append(contentsOf: [.zero, .zero, .zero]); return   // keep 3 vertices per element
            }
            for k in 0..<3 { tri.append(SIMD2(Float(c[k][0]), Float(c[k][1]))) }
        }
        if FEMPreviewGPUCache.usesTiling(mesh, periodic: periodic) {
            tri.reserveCapacity(mesh.elements_3x3.count * 3)
            for e in mesh.elements_3x3 { add(e.node_coords) }
        } else {
            tri.reserveCapacity(mesh.elements.count * 3)
            for e in mesh.elements { add(e.node_coords) }
        }
        
        var nodes: [SIMD2<Float>] = []
        var colors: [SIMD4<Float>] = []
        let free = SIMD4<Float>(0.30, 0.50, 1.00, 1), fixed = SIMD4<Float>(1.00, 0.42, 0.42, 1)
        let driven = SIMD4<Float>(0.31, 0.80, 0.77, 1), ghost = SIMD4<Float>(0.55, 0.55, 0.60, 0.65)
        if periodic {
            for d in mesh.translated_domains {
                let primary = (d.tx == 0 && d.ty == 0)
                for n in d.nodes { nodes.append(SIMD2(Float(n.x), Float(n.y))); colors.append(primary ? free : ghost) }
            }
        } else {
            for n in mesh.primary_nodes {
                nodes.append(SIMD2(Float(n.x), Float(n.y)))
                let role = n.bc_role ?? ((n.is_boundary == true) ? 1 : 0)
                colors.append(role == 2 ? driven : (role == 1 ? fixed : free))
            }
        }
        
        var outline: [SIMD2<Float>] = []
        let lx = Float(mesh.Lx), ly = Float(mesh.Ly)
        let tiles: [(Float, Float)] = periodic
            ? [(-1, -1), (0, -1), (1, -1), (-1, 0), (0, 0), (1, 0), (-1, 1), (0, 1), (1, 1)]
            : [(0, 0)]
        for (tx, ty) in tiles {
            let a = SIMD2(tx * lx, ty * ly), b = SIMD2((tx + 1) * lx, ty * ly)
            let c = SIMD2((tx + 1) * lx, (ty + 1) * ly), d = SIMD2(tx * lx, (ty + 1) * ly)
            outline.append(contentsOf: [a, b, b, c, c, d, d, a])
        }
        
        let g = MeshGeometry(triPositions: tri, nodePositions: nodes, nodeColors: colors, outlineSegments: outline)
        geometry = g
        locator = TriangleLocator(triPositions: tri)
        geometryKey = key
        valuesKey = ""
        return g
    }
    
    func values(for mesh: InitialMeshModel, version: Int, periodic: Bool, field: InitialMeshColorField) -> MeshValues? {
        guard field != .uniform else { return nil }
        let key = "\(version)|\(periodic)|\(field.rawValue)"
        if key == valuesKey { return values }
        func pick(_ e: InitialMeshElement) -> Float { Float(FEMPreviewGPUCache.value(field, e.strain_energy_density, e.von_mises, e.sigma_xy, e.pressure, e.sigma_xx, e.sigma_yy)) }
        func pick3(_ e: InitialMesh3x3Element) -> Float { Float(FEMPreviewGPUCache.value(field, e.strain_energy_density, e.von_mises, e.sigma_xy, e.pressure, e.sigma_xx, e.sigma_yy)) }
        var per: [Float] = []
        if FEMPreviewGPUCache.usesTiling(mesh, periodic: periodic) {
            per.reserveCapacity(mesh.elements_3x3.count * 3)
            for e in mesh.elements_3x3 { let v = pick3(e); per.append(contentsOf: [v, v, v]) }
        } else {
            per.reserveCapacity(mesh.elements.count * 3)
            for e in mesh.elements { let v = pick(e); per.append(contentsOf: [v, v, v]) }
        }
        valuesVersion += 1
        values = MeshValues(version: valuesVersion, perVertex: per)
        valuesKey = key
        return values
    }
    
    static func value(_ f: InitialMeshColorField, _ w: Double, _ vm: Double, _ sxy: Double, _ p: Double, _ sxx: Double, _ syy: Double) -> Double {
        switch f {
        case .uniform: return 0
        case .energyDensity: return w
        case .vonMises: return vm
        case .shearStress: return sxy
        case .pressure: return p
        case .sigmaXX: return sxx
        case .sigmaYY: return syy
        }
    }
}

extension StiffnessColorPalette {
    var meshColormap: MeshColormap {
        switch self {
        case .coolwarm: return .coolwarm
        case .viridis: return .viridis
        case .plasma: return .plasma
        case .grayscale: return .grayscale
        }
    }
}
