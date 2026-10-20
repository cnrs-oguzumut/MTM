import Foundation

public struct InitialMeshNode: Codable, Identifiable {
    public var id: Int
    public var x: Double
    public var y: Double
    public var is_boundary: Bool?
    
    public init(id: Int, x: Double, y: Double, is_boundary: Bool? = false) {
        self.id = id
        self.x = x
        self.y = y
        self.is_boundary = is_boundary
    }
}

public struct InitialMeshDomain: Codable {
    public var tx: Int
    public var ty: Int
    public var offset_x: Double
    public var offset_y: Double
    public var nodes: [InitialMeshNode]
}

public struct InitialMeshElement: Codable, Identifiable {
    public var id: Int
    public var nodes: [Int]
    public var node_coords: [[Double]] // 3 pairs of [x, y]
    public var translations: [[Double]]?
    public var crosses_boundary: Bool?
    public var area: Double
    public var strain_energy_density: Double
    public var sigma_xx: Double
    public var sigma_yy: Double
    public var sigma_xy: Double
    public var von_mises: Double
    public var pressure: Double
}

public struct InitialMesh3x3Element: Codable, Identifiable {
    public var id: String { "\(domain_tx)_\(domain_ty)_\(orig_id)" }
    public var domain_tx: Int
    public var domain_ty: Int
    public var orig_id: Int
    public var node_coords: [[Double]]
    public var strain_energy_density: Double
    public var sigma_xx: Double
    public var sigma_yy: Double
    public var sigma_xy: Double
    public var von_mises: Double
    public var pressure: Double
}

public struct InitialMeshStatistics: Codable {
    public var min_energy: Double = 0.0
    public var max_energy: Double = 0.0
    public var min_von_mises: Double = 0.0
    public var max_von_mises: Double = 0.0
    public var min_shear: Double = 0.0
    public var max_shear: Double = 0.0
    public var min_pressure: Double = 0.0
    public var max_pressure: Double = 0.0
    
    public init() {}
}

public struct InitialMeshModel: Codable {
    public var lattice_type: String = "square"
    public var bc: String = "periodic"
    public var rotation_angle: Double = 0.0
    public var is_csl: Bool = false
    public var m: Int = 1
    public var n: Int = 0
    public var sigma: Int = 1
    public var kx: Int = 1
    public var ky: Int = 1
    public var Lx: Double = 1.0
    public var Ly: Double = 1.0
    public var num_primary_nodes: Int = 0
    public var num_elements: Int = 0
    public var primary_nodes: [InitialMeshNode] = []
    public var translated_domains: [InitialMeshDomain] = []
    public var elements: [InitialMeshElement] = []
    public var elements_3x3: [InitialMesh3x3Element] = []
    public var statistics: InitialMeshStatistics = InitialMeshStatistics()
    
    public init() {}
}

public struct CSLAngleOption: Identifiable, Hashable {
    public var id: String { "\(m)_\(n)" }
    public var m: Int
    public var n: Int
    public var sigma: Int
    public var theta_deg: Double
    public var Lx_ratio: Double
    public var Ly_ratio: Double
    public var atoms_cell: Int
    public var is_reduced: Bool
    public var lattice_type: String
    
    public var displayName: String {
        let degStr = String(format: "%.2f°", theta_deg)
        let redStr = is_reduced ? " (reduced)" : ""
        return "θ=\(degStr) [m=\(m), n=\(n), Σ=\(sigma), N=\(atoms_cell)]\(redStr)"
    }
    
    public static func allowedAngles(for lattice: LatticeSymmetry) -> [CSLAngleOption] {
        if lattice == .square {
            return [
                CSLAngleOption(m: 1, n: 0, sigma: 1, theta_deg: 0.0, Lx_ratio: 1.0, Ly_ratio: 1.0, atoms_cell: 1, is_reduced: false, lattice_type: "square"),
                CSLAngleOption(m: 3, n: 1, sigma: 10, theta_deg: 18.435, Lx_ratio: 3.162, Ly_ratio: 3.162, atoms_cell: 10, is_reduced: false, lattice_type: "square"),
                CSLAngleOption(m: 7, n: 1, sigma: 50, theta_deg: 22.620, Lx_ratio: 7.071, Ly_ratio: 7.071, atoms_cell: 50, is_reduced: false, lattice_type: "square"),
                CSLAngleOption(m: 2, n: 1, sigma: 5, theta_deg: 26.565, Lx_ratio: 2.236, Ly_ratio: 2.236, atoms_cell: 5, is_reduced: false, lattice_type: "square"),
                CSLAngleOption(m: 5, n: 2, sigma: 29, theta_deg: 28.072, Lx_ratio: 5.385, Ly_ratio: 5.385, atoms_cell: 29, is_reduced: false, lattice_type: "square"),
                CSLAngleOption(m: 8, n: 5, sigma: 89, theta_deg: 31.891, Lx_ratio: 9.434, Ly_ratio: 9.434, atoms_cell: 89, is_reduced: false, lattice_type: "square"),
                CSLAngleOption(m: 3, n: 2, sigma: 13, theta_deg: 33.690, Lx_ratio: 3.606, Ly_ratio: 3.606, atoms_cell: 13, is_reduced: false, lattice_type: "square"),
                CSLAngleOption(m: 4, n: 3, sigma: 25, theta_deg: 36.870, Lx_ratio: 5.0, Ly_ratio: 5.0, atoms_cell: 25, is_reduced: false, lattice_type: "square"),
                CSLAngleOption(m: 6, n: 5, sigma: 61, theta_deg: 39.806, Lx_ratio: 7.810, Ly_ratio: 7.810, atoms_cell: 61, is_reduced: false, lattice_type: "square"),
                CSLAngleOption(m: 1, n: 1, sigma: 2, theta_deg: 45.000, Lx_ratio: 1.414, Ly_ratio: 1.414, atoms_cell: 2, is_reduced: false, lattice_type: "square")
            ]
        } else {
            return [
                CSLAngleOption(m: 1, n: 0, sigma: 1, theta_deg: 0.0, Lx_ratio: 1.0, Ly_ratio: 1.732, atoms_cell: 2, is_reduced: false, lattice_type: "triangular"),
                CSLAngleOption(m: 4, n: 1, sigma: 21, theta_deg: 10.893, Lx_ratio: 4.583, Ly_ratio: 2.646, atoms_cell: 14, is_reduced: true, lattice_type: "triangular"),
                CSLAngleOption(m: 3, n: 1, sigma: 13, theta_deg: 13.898, Lx_ratio: 3.606, Ly_ratio: 6.245, atoms_cell: 26, is_reduced: false, lattice_type: "triangular"),
                CSLAngleOption(m: 5, n: 2, sigma: 39, theta_deg: 16.102, Lx_ratio: 6.245, Ly_ratio: 3.606, atoms_cell: 26, is_reduced: true, lattice_type: "triangular"),
                CSLAngleOption(m: 2, n: 1, sigma: 7, theta_deg: 19.107, Lx_ratio: 2.646, Ly_ratio: 4.583, atoms_cell: 14, is_reduced: false, lattice_type: "triangular"),
                CSLAngleOption(m: 5, n: 3, sigma: 49, theta_deg: 21.787, Lx_ratio: 7.000, Ly_ratio: 12.124, atoms_cell: 98, is_reduced: false, lattice_type: "triangular"),
                CSLAngleOption(m: 3, n: 2, sigma: 19, theta_deg: 23.413, Lx_ratio: 4.359, Ly_ratio: 7.550, atoms_cell: 38, is_reduced: false, lattice_type: "triangular"),
                CSLAngleOption(m: 4, n: 3, sigma: 37, theta_deg: 25.285, Lx_ratio: 6.083, Ly_ratio: 10.536, atoms_cell: 74, is_reduced: false, lattice_type: "triangular"),
                CSLAngleOption(m: 1, n: 1, sigma: 3, theta_deg: 30.000, Lx_ratio: 1.732, Ly_ratio: 1.0, atoms_cell: 2, is_reduced: true, lattice_type: "triangular"),
                CSLAngleOption(m: 3, n: 4, sigma: 37, theta_deg: 34.715, Lx_ratio: 6.083, Ly_ratio: 10.536, atoms_cell: 74, is_reduced: false, lattice_type: "triangular"),
                CSLAngleOption(m: 2, n: 3, sigma: 19, theta_deg: 36.587, Lx_ratio: 4.359, Ly_ratio: 7.550, atoms_cell: 38, is_reduced: false, lattice_type: "triangular"),
                CSLAngleOption(m: 1, n: 2, sigma: 7, theta_deg: 40.893, Lx_ratio: 2.646, Ly_ratio: 4.583, atoms_cell: 14, is_reduced: false, lattice_type: "triangular"),
                CSLAngleOption(m: 2, n: 5, sigma: 39, theta_deg: 43.898, Lx_ratio: 6.245, Ly_ratio: 3.606, atoms_cell: 26, is_reduced: true, lattice_type: "triangular"),
                CSLAngleOption(m: 1, n: 3, sigma: 13, theta_deg: 46.102, Lx_ratio: 3.606, Ly_ratio: 6.245, atoms_cell: 26, is_reduced: false, lattice_type: "triangular"),
                CSLAngleOption(m: 1, n: 4, sigma: 21, theta_deg: 49.107, Lx_ratio: 4.583, Ly_ratio: 2.646, atoms_cell: 14, is_reduced: true, lattice_type: "triangular"),
                CSLAngleOption(m: 0, n: 1, sigma: 1, theta_deg: 60.000, Lx_ratio: 1.0, Ly_ratio: 1.732, atoms_cell: 2, is_reduced: false, lattice_type: "triangular"),
                CSLAngleOption(m: -1, n: 5, sigma: 21, theta_deg: 70.893, Lx_ratio: 4.583, Ly_ratio: 2.646, atoms_cell: 14, is_reduced: true, lattice_type: "triangular"),
                CSLAngleOption(m: -1, n: 4, sigma: 13, theta_deg: 73.898, Lx_ratio: 3.606, Ly_ratio: 6.245, atoms_cell: 26, is_reduced: false, lattice_type: "triangular"),
                CSLAngleOption(m: -2, n: 7, sigma: 39, theta_deg: 76.102, Lx_ratio: 6.245, Ly_ratio: 3.606, atoms_cell: 26, is_reduced: true, lattice_type: "triangular"),
                CSLAngleOption(m: -1, n: 3, sigma: 7, theta_deg: 79.107, Lx_ratio: 2.646, Ly_ratio: 4.583, atoms_cell: 14, is_reduced: false, lattice_type: "triangular"),
                CSLAngleOption(m: -2, n: 5, sigma: 19, theta_deg: 83.413, Lx_ratio: 4.359, Ly_ratio: 7.550, atoms_cell: 38, is_reduced: false, lattice_type: "triangular"),
                CSLAngleOption(m: -3, n: 7, sigma: 37, theta_deg: 85.285, Lx_ratio: 6.083, Ly_ratio: 10.536, atoms_cell: 74, is_reduced: false, lattice_type: "triangular"),
                CSLAngleOption(m: -1, n: 2, sigma: 3, theta_deg: 90.000, Lx_ratio: 1.732, Ly_ratio: 1.0, atoms_cell: 2, is_reduced: true, lattice_type: "triangular"),
                CSLAngleOption(m: -4, n: 7, sigma: 37, theta_deg: 94.715, Lx_ratio: 6.083, Ly_ratio: 10.536, atoms_cell: 74, is_reduced: false, lattice_type: "triangular"),
                CSLAngleOption(m: -3, n: 5, sigma: 19, theta_deg: 96.587, Lx_ratio: 4.359, Ly_ratio: 7.550, atoms_cell: 38, is_reduced: false, lattice_type: "triangular"),
                CSLAngleOption(m: -2, n: 3, sigma: 7, theta_deg: 100.893, Lx_ratio: 2.646, Ly_ratio: 4.583, atoms_cell: 14, is_reduced: false, lattice_type: "triangular"),
                CSLAngleOption(m: -5, n: 7, sigma: 39, theta_deg: 103.898, Lx_ratio: 6.245, Ly_ratio: 3.606, atoms_cell: 26, is_reduced: true, lattice_type: "triangular"),
                CSLAngleOption(m: -3, n: 4, sigma: 13, theta_deg: 106.102, Lx_ratio: 3.606, Ly_ratio: 6.245, atoms_cell: 26, is_reduced: false, lattice_type: "triangular"),
                CSLAngleOption(m: -4, n: 5, sigma: 21, theta_deg: 109.107, Lx_ratio: 4.583, Ly_ratio: 2.646, atoms_cell: 14, is_reduced: true, lattice_type: "triangular"),
                CSLAngleOption(m: -1, n: 1, sigma: 1, theta_deg: 120.000, Lx_ratio: 1.0, Ly_ratio: 1.732, atoms_cell: 2, is_reduced: false, lattice_type: "triangular")
            ]
        }
    }
}
