import Foundation

public struct DislocationProfileJSONRoot: Codable {
    public let status: String
    public let latticeType: String
    public let potentialMode: Int
    public let beta: Double
    public let bulkK: Double
    public let rFree: Double
    public let latticeConstant: Double
    public let burgersVector: Double
    public let coreX: Double
    public let coreY: Double
    public let elasticConstants: DislocationElasticConstantsModel
    public let unrelaxed: DislocationStateDataModel
    public let relaxed: DislocationStateDataModel?

    enum CodingKeys: String, CodingKey {
        case status
        case latticeType = "lattice_type"
        case potentialMode = "potential_mode"
        case beta
        case bulkK = "bulk_K"
        case rFree = "R_free"
        case latticeConstant = "lattice_constant"
        case burgersVector = "burgers_vector"
        case coreX = "core_x"
        case coreY = "core_y"
        case elasticConstants = "elastic_constants"
        case unrelaxed
        case relaxed
    }

    public init(from decoder: Decoder) throws {
        let c = try decoder.container(keyedBy: CodingKeys.self)
        status = (try? c.decode(String.self, forKey: .status)) ?? "completed"
        latticeType = (try? c.decode(String.self, forKey: .latticeType)) ?? "triangular"
        potentialMode = (try? c.decode(Int.self, forKey: .potentialMode)) ?? 2
        beta = (try? c.decode(Double.self, forKey: .beta)) ?? -0.25
        bulkK = (try? c.decode(Double.self, forKey: .bulkK)) ?? 5.0
        rFree = (try? c.decode(Double.self, forKey: .rFree)) ?? 20.0
        latticeConstant = (try? c.decode(Double.self, forKey: .latticeConstant)) ?? 1.0
        burgersVector = (try? c.decode(Double.self, forKey: .burgersVector)) ?? latticeConstant
        coreX = (try? c.decode(Double.self, forKey: .coreX)) ?? 0.0
        coreY = (try? c.decode(Double.self, forKey: .coreY)) ?? 0.0
        elasticConstants = (try? c.decode(DislocationElasticConstantsModel.self, forKey: .elasticConstants))
            ?? DislocationElasticConstantsModel(c11: 21.7, c12: 18.2, c66: 6.0, poissonRatio: 0.25)
        unrelaxed = try c.decode(DislocationStateDataModel.self, forKey: .unrelaxed)
        relaxed = try? c.decode(DislocationStateDataModel.self, forKey: .relaxed)
    }
}

public struct DislocationElasticConstantsModel: Codable {
    public let c11: Double
    public let c12: Double
    public let c66: Double
    public let poissonRatio: Double

    public init(c11: Double = 0, c12: Double = 0, c66: Double = 0, poissonRatio: Double = 0.25) {
        self.c11 = c11
        self.c12 = c12
        self.c66 = c66
        self.poissonRatio = poissonRatio
    }

    enum CodingKeys: String, CodingKey {
        case c11 = "C11"
        case c12 = "C12"
        case c66 = "C66"
        case poissonRatio = "poisson_ratio"
    }
}

public struct DislocationCoreFitModel: Codable {
    public let coreEnergy: Double
    public let elasticK: Double
    public let r2: Double

    enum CodingKeys: String, CodingKey {
        case coreEnergy = "core_energy"
        case elasticK = "elastic_K"
        case r2 = "R2"
    }
}

public struct DislocationStateDataModel: Codable {
    public let totalEnergy: Double
    public let totalStress: Double
    public let coreFit: DislocationCoreFitModel
    public let radial: [DislocationRadialPointModel]
    public let slipPlane: [DislocationSlipPointModel]
    public let slipAvg: [DislocationSlipAvgPointModel]

    enum CodingKeys: String, CodingKey {
        case totalEnergy = "total_energy"
        case totalStress = "total_stress"
        case coreFit = "core_fit"
        case radial
        case slipPlane = "slip_plane"
        case slipAvg = "slip_avg"
    }
}

public struct DislocationRadialPointModel: Codable, Identifiable {
    public var id: Double { r }
    public let r: Double
    public let logR: Double
    public let cumulativeEnergy: Double
    public let energyDensity: Double
    public let avgAtomEnergy: Double
    public let sXx: Double
    public let sYy: Double
    public let sXy: Double
    public let sHydro: Double
    public let sVm: Double
    public let atomCount: Int

    enum CodingKeys: String, CodingKey {
        case r
        case logR = "log_r"
        case cumulativeEnergy = "cumulative_energy"
        case energyDensity = "energy_density"
        case avgAtomEnergy = "avg_atom_energy"
        case sXx = "s_xx"
        case sYy = "s_yy"
        case sXy = "s_xy"
        case sHydro = "s_hydro"
        case sVm = "s_vm"
        case atomCount = "atom_count"
    }
}

public struct DislocationSlipPointModel: Codable, Identifiable {
    public var id: String { "\(dx)_\(dy)_\(row)" }
    public let dx: Double
    public let dy: Double
    public let x: Double
    public let y: Double
    public let r: Double
    public let energy: Double
    public let sXx: Double
    public let sYy: Double
    public let sXy: Double
    public let sHydro: Double
    public let sVm: Double
    public let row: String

    enum CodingKeys: String, CodingKey {
        case dx, dy, x, y, r, energy
        case sXx = "s_xx"
        case sYy = "s_yy"
        case sXy = "s_xy"
        case sHydro = "s_hydro"
        case sVm = "s_vm"
        case row
    }
}

public struct DislocationSlipAvgPointModel: Codable, Identifiable {
    public var id: Double { dx }
    public let dx: Double
    public let energy: Double
    public let sXx: Double
    public let sYy: Double
    public let sXy: Double
    public let sHydro: Double
    public let sVm: Double

    enum CodingKeys: String, CodingKey {
        case dx, energy
        case sXx = "s_xx"
        case sYy = "s_yy"
        case sXy = "s_xy"
        case sHydro = "s_hydro"
        case sVm = "s_vm"
    }
}
