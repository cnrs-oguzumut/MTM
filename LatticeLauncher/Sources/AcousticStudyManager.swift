import Foundation
import Combine
import SwiftUI
import AppKit

// Note: AcousticPotentialMode is defined and aliased to PotentialMode in SimulationConfig.swift

public enum ColorFieldMode: String, CaseIterable, Identifiable {
    case acousticStability = "Acoustic Stability (min det Q)"
    case strainEnergy = "Strain Energy Density Φ(C)"
    
    public var id: String { rawValue }
    public var shortLabel: String {
        switch self {
        case .acousticStability: return "Stability"
        case .strainEnergy: return "Energy"
        }
    }
}

public enum RenderStyleMode: String, CaseIterable, Identifiable {
    case smoothSurface = "Smooth Surface"
    case scatterPoints = "Point Scatter"
    
    public var id: String { rawValue }
}

public enum EnergyColormapPreset: String, CaseIterable, Identifiable {
    case inferno = "Inferno"
    case viridis = "Viridis"
    case plasma = "Plasma"
    case turbo = "Turbo"
    case spectral = "Spectral"
    case coolwarm = "CoolWarm"
    
    public var id: String { rawValue }
}

public struct ContourSegment: Identifiable {
    public let id = UUID()
    public let startX: CGFloat
    public let startY: CGFloat
    public let endX: CGFloat
    public let endY: CGFloat
    public let startXRed: CGFloat
    public let startYRed: CGFloat
    public let endXRed: CGFloat
    public let endYRed: CGFloat
}

public struct SurfaceCell: Identifiable {
    public let id = UUID()
    public let ptsOrig: [CGPoint] // 3 or 4 points in Poincaré space [-1, 1]
    public let ptsRed: [CGPoint]
    public let avgDetAc: Double
    public let avgEnergy: Double
    public let maxEnergy: Double
}

public struct AcousticPointModel: Identifiable, Codable {
    public var id: UUID = UUID()
    public var t: Double = 0.0
    public var p: Double = 0.0
    public var c11: Double = 1.0
    public var c22: Double = 1.0
    public var c12: Double = 0.0
    public var c11_red: Double = 1.0
    public var c22_red: Double = 1.0
    public var c12_red: Double = 0.0
    public var min_detAc: Double = 0.0
    public var energy: Double = 0.0
    public var angle_deg: Double = 0.0
    public var third_condition: Bool = false
    public var x: Double = 0.0
    public var y: Double = 0.0
    public var x_red: Double = 0.0
    public var y_red: Double = 0.0
    public var is_stable: Bool = true
    
    enum CodingKeys: String, CodingKey {
        case t, p, c11, c22, c12, c11_red, c22_red, c12_red
        case min_detAc, energy, angle_deg, third_condition, x, y, x_red, y_red, is_stable
    }
    
    public init(t: Double = 0, p: Double = 0, c11: Double = 1, c22: Double = 1, c12: Double = 0,
                c11_red: Double = 1, c22_red: Double = 1, c12_red: Double = 0, min_detAc: Double = 0,
                energy: Double = 0.0, angle_deg: Double = 0, third_condition: Bool = false,
                x: Double = 0, y: Double = 0, x_red: Double = 0, y_red: Double = 0, is_stable: Bool = true) {
        self.id = UUID()
        self.t = t
        self.p = p
        self.c11 = c11
        self.c22 = c22
        self.c12 = c12
        self.c11_red = c11_red
        self.c22_red = c22_red
        self.c12_red = c12_red
        self.min_detAc = min_detAc
        self.energy = energy
        self.angle_deg = angle_deg
        self.third_condition = third_condition
        self.x = x
        self.y = y
        self.x_red = x_red
        self.y_red = y_red
        self.is_stable = is_stable
    }
    
    public init(from decoder: Decoder) throws {
        let container = try decoder.container(keyedBy: CodingKeys.self)
        self.id = UUID()
        self.t = try container.decodeIfPresent(Double.self, forKey: .t) ?? 0.0
        self.p = try container.decodeIfPresent(Double.self, forKey: .p) ?? 0.0
        self.c11 = try container.decodeIfPresent(Double.self, forKey: .c11) ?? 1.0
        self.c22 = try container.decodeIfPresent(Double.self, forKey: .c22) ?? 1.0
        self.c12 = try container.decodeIfPresent(Double.self, forKey: .c12) ?? 0.0
        self.c11_red = try container.decodeIfPresent(Double.self, forKey: .c11_red) ?? 1.0
        self.c22_red = try container.decodeIfPresent(Double.self, forKey: .c22_red) ?? 1.0
        self.c12_red = try container.decodeIfPresent(Double.self, forKey: .c12_red) ?? 0.0
        self.min_detAc = try container.decodeIfPresent(Double.self, forKey: .min_detAc) ?? 0.0
        self.energy = try container.decodeIfPresent(Double.self, forKey: .energy) ?? 0.0
        self.angle_deg = try container.decodeIfPresent(Double.self, forKey: .angle_deg) ?? 0.0
        self.third_condition = try container.decodeIfPresent(Bool.self, forKey: .third_condition) ?? false
        self.x = try container.decodeIfPresent(Double.self, forKey: .x) ?? 0.0
        self.y = try container.decodeIfPresent(Double.self, forKey: .y) ?? 0.0
        self.x_red = try container.decodeIfPresent(Double.self, forKey: .x_red) ?? 0.0
        self.y_red = try container.decodeIfPresent(Double.self, forKey: .y_red) ?? 0.0
        self.is_stable = try container.decodeIfPresent(Bool.self, forKey: .is_stable) ?? true
    }
}

public struct AcousticElasticConstantsModel: Codable {
    public var equilibrium_lattice_parameter: Double = 1.0
    public var C11: Double = 0.0
    public var C22: Double = 0.0
    public var C12: Double = 0.0
    public var C66: Double = 0.0
    public var Zener_A: Double = 1.0
    public var bulk_modulus: Double = 0.0
    public var shear_modulus: Double = 0.0
    public var poisson_ratio: Double = 0.0
    public var stroh_kappa: Double = 1.0
    public var stroh_xi: Double = 0.0
    public var stroh_zeta: Double = 0.0
    public var stroh_g1: Double = 0.0
    public var stroh_g2: Double = 0.0
    public var stroh_g3: Double = 0.0
    public var foreman_K: Double = 0.0
    public var born_stable: Bool = true
    public var is_isotropic: Bool = true
    public var anisotropy_comment: String = ""
    
    enum CodingKeys: String, CodingKey {
        case equilibrium_lattice_parameter
        case C11, C22, C12, C66, Zener_A
        case bulk_modulus, shear_modulus
        case poisson_ratio, stroh_kappa, stroh_xi, stroh_zeta
        case stroh_g1, stroh_g2, stroh_g3, foreman_K
        case born_stable, is_isotropic, anisotropy_comment
    }
    
    public init() {}
    
    public init(from decoder: Decoder) throws {
        let container = try decoder.container(keyedBy: CodingKeys.self)
        self.equilibrium_lattice_parameter = try container.decodeIfPresent(Double.self, forKey: .equilibrium_lattice_parameter) ?? 1.0
        self.C11 = try container.decodeIfPresent(Double.self, forKey: .C11) ?? 0.0
        self.C22 = try container.decodeIfPresent(Double.self, forKey: .C22) ?? 0.0
        self.C12 = try container.decodeIfPresent(Double.self, forKey: .C12) ?? 0.0
        self.C66 = try container.decodeIfPresent(Double.self, forKey: .C66) ?? 0.0
        self.Zener_A = try container.decodeIfPresent(Double.self, forKey: .Zener_A) ?? 1.0
        self.bulk_modulus = try container.decodeIfPresent(Double.self, forKey: .bulk_modulus) ?? 0.0
        self.shear_modulus = try container.decodeIfPresent(Double.self, forKey: .shear_modulus) ?? 0.0
        self.born_stable = try container.decodeIfPresent(Bool.self, forKey: .born_stable) ?? true
        self.is_isotropic = try container.decodeIfPresent(Bool.self, forKey: .is_isotropic) ?? true
        self.anisotropy_comment = try container.decodeIfPresent(String.self, forKey: .anisotropy_comment) ?? ""
        
        if let nu = try container.decodeIfPresent(Double.self, forKey: .poisson_ratio) {
            self.poisson_ratio = nu
            self.stroh_kappa = try container.decodeIfPresent(Double.self, forKey: .stroh_kappa) ?? 1.0
            self.stroh_xi = try container.decodeIfPresent(Double.self, forKey: .stroh_xi) ?? 0.0
            self.stroh_zeta = try container.decodeIfPresent(Double.self, forKey: .stroh_zeta) ?? 0.0
            self.stroh_g1 = try container.decodeIfPresent(Double.self, forKey: .stroh_g1) ?? 0.0
            self.stroh_g2 = try container.decodeIfPresent(Double.self, forKey: .stroh_g2) ?? 0.0
            self.stroh_g3 = try container.decodeIfPresent(Double.self, forKey: .stroh_g3) ?? 0.0
            self.foreman_K = try container.decodeIfPresent(Double.self, forKey: .foreman_K) ?? 0.0
        } else {
            self.populateDislocationParamsFallback()
        }
    }
    
    public mutating func populateDislocationParamsFallback() {
        let c11Safe = abs(C11) > 1e-12 ? C11 : 1.0
        let c66Safe = abs(C66) > 1e-12 ? abs(C66) : 1.0
        self.poisson_ratio = C12 / c11Safe
        let num = C11 * C11 - C12 * C12 - 2.0 * C12 * C66
        let den = 2.0 * abs(C11) * c66Safe
        self.stroh_kappa = den > 1e-12 ? (num / den) : 1.0
        if stroh_kappa > 1.0 + 1e-7 {
            let s = sqrt(max(0.0, stroh_kappa * stroh_kappa - 1.0))
            self.stroh_xi = sqrt(max(0.0, stroh_kappa + s))
            self.stroh_zeta = sqrt(max(0.0, stroh_kappa - s))
        } else if stroh_kappa < 1.0 - 1e-7 {
            self.stroh_xi = sqrt(max(0.0, 0.5 * (1.0 + stroh_kappa)))
            self.stroh_zeta = sqrt(max(0.0, 0.5 * (1.0 - stroh_kappa)))
        } else {
            self.stroh_xi = 1.0
            self.stroh_zeta = 0.0
        }
        let d = C12 + C66
        if abs(d) > 1e-12 {
            self.stroh_g1 = (C11 - C66 * stroh_kappa) / d
            self.stroh_g2 = C11 / d
            self.stroh_g3 = (C11 + C12) * C66 / (C11 * d)
        }
        let twoPlus2k = max(0.0, 2.0 + 2.0 * stroh_kappa)
        if abs(C11) > 1e-12 && twoPlus2k > 1e-12 {
            self.foreman_K = (C11 * C11 - C12 * C12) * sqrt(C66 / (C11 * twoPlus2k))
        } else {
            self.foreman_K = C66 * (1.0 + poisson_ratio)
        }
    }
}

public struct AcousticStudySummary: Codable {
    public var total_points: Int = 0
    public var stable_points: Int = 0
    public var unstable_points: Int = 0
    public var min_detAc_min: Double = 0.0
    public var min_detAc_max: Double = 0.0
    public var energy_min: Double = 0.0
    public var energy_max: Double = 0.0
}

public struct AcousticStudyJSONRoot: Codable {
    public var mode: Int?
    public var mode_name: String?
    public var lattice_type: String?
    public var beta: Double?
    public var K: Double?
    public var t_max: Double?
    public var dt: Double?
    public var dp: Double?
    public var elastic_constants: AcousticElasticConstantsModel?
    public var summary: AcousticStudySummary?
    public var points: [AcousticPointModel]?
    public var shear_path: ShearPathModel?
}

public struct ShearPathPointModel: Codable, Identifiable {
    public var id: Double { alpha }
    public var alpha: Double = 0.0
    public var x: Double = 0.0
    public var y: Double = 0.0
    public var x_red: Double = 0.0
    public var y_red: Double = 0.0
    public var c11: Double = 1.0
    public var c22: Double = 1.0
    public var c12: Double = 0.0
    public var min_detAc: Double? = nil
    public var energy: Double = 0.0
    public var is_stable: Bool? = nil
    public var angle_deg: Double? = nil
    public var angle_euler_deg: Double? = nil
    public var lambda_min_Q: Double? = nil
    public var lambda_max_Q: Double? = nil
    public var mode_mx: Double? = nil
    public var mode_my: Double? = nil
    public var f11: Double? = nil
    public var f12: Double? = nil
    public var f21: Double? = nil
    public var f22: Double? = nil
    public var q11_xx: Double? = nil
    public var q11_yy: Double? = nil
    public var q11_xy: Double? = nil
    public var q22_xx: Double? = nil
    public var q22_yy: Double? = nil
    public var q22_xy: Double? = nil
    public var q12_xx: Double? = nil
    public var q12_yy: Double? = nil
    public var q12_xy: Double? = nil
}

public struct ShearPathModel: Codable {
    public var path_type: String? = "simple_shear"
    public var theta_deg: Double? = 0.0
    public var period_alpha: Double? = 1.0
    public var alpha_min: Double = -1.5
    public var alpha_max: Double = 1.5
    public var dalpha: Double = 0.005
    public var critical_alphas: [Double] = []
    public var critical_points: [ShearPathPointModel] = []
    public var points: [ShearPathPointModel] = []
}

public class AcousticStudyManager: ObservableObject {
    // Parameters
    @Published public var potentialMode: AcousticPotentialMode = .zanzotto
    @Published public var latticeSymmetry: LatticeSymmetry = .square
    @Published public var beta: Double = -0.25
    @Published public var bulkK: Double = 5.0
    
    // Sweep Grid Settings
    @Published public var tMax: Double = 1.4
    @Published public var dt: Double = 0.02
    @Published public var dp: Double = 0.05
    
    // Deformation path settings
    @Published public var showShearPath: Bool = true
    @Published public var pathType: String = "simple_shear" // "simple_shear", "pure_shear", "diagonal_shear"
    @Published public var shearTheta: Double = 0.0 // Angle in degrees
    @Published public var shearAlphaMin: Double = -1.5
    @Published public var shearAlphaMax: Double = 1.5
    @Published public var shearDAlpha: Double = 0.005
    @Published public var shearPath: ShearPathModel? = nil
    
    // Bidirectional Cursor Synchronization
    @Published public var hoveredAlpha: Double? = nil
    @Published public var hoveredPathPoint: ShearPathPointModel? = nil
    @Published public var selectedCriticalAlpha: Double? = nil
    
    public var estimatedPoints: Int {
        let tSteps = (tMax / max(0.001, dt)) + 1.0
        let pSteps = (2.0 * Double.pi / max(0.001, dp)) + 1.0
        return Int(tSteps * pSteps)
    }
    
    // View Options
    @Published public var useReducedCoordinates: Bool = false
    @Published public var colorField: ColorFieldMode = .acousticStability
    @Published public var renderStyle: RenderStyleMode = .smoothSurface
    @Published public var highlightZeroBoundary: Bool = true
    @Published public var pointRadius: Double = 3.5
    @Published public var colormapPreset: String = "RdYlBu"
    
    // Poincaré Energy Color Scale & Percentile Filtering
    @Published public var energyColormap: EnergyColormapPreset = .inferno
    @Published public var invertEnergyColormap: Bool = false
    @Published public var energyUpperPercentile: Double = 95.0
    @Published public var energyLowerPercentile: Double = 0.0
    @Published public var discardAbovePercentile: Bool = true
    @Published public var energyGamma: Double = 0.65
    @Published public var useCustomEnergyRange: Bool = false
    @Published public var customEnergyMin: Double = 0.0
    @Published public var customEnergyMax: Double = 10.0
    @Published public private(set) var sortedEnergies: [Double] = []
    
    public var effectiveEnergyBounds: (min: Double, max: Double, keptCount: Int, totalCount: Int) {
        let total = sortedEnergies.count
        guard total > 0 else {
            return (min: summary.energy_min, max: max(summary.energy_max, summary.energy_min + 1e-6), keptCount: 0, totalCount: 0)
        }
        if useCustomEnergyRange {
            let lo = customEnergyMin
            let hi = max(customEnergyMax, lo + 1e-9)
            // Binary search or filter count
            let kept = sortedEnergies.filter { $0 <= hi }.count
            return (min: lo, max: hi, keptCount: kept, totalCount: total)
        }
        let loPct = max(0.0, min(99.0, energyLowerPercentile)) / 100.0
        let hiPct = max(loPct + 0.001, min(100.0, energyUpperPercentile)) / 100.0
        let loIdx = min(total - 1, max(0, Int(Double(total - 1) * loPct)))
        let hiIdx = min(total - 1, max(loIdx, Int(Double(total - 1) * hiPct)))
        let loVal = sortedEnergies[loIdx]
        let hiVal = max(sortedEnergies[hiIdx], loVal + 1e-9)
        return (min: loVal, max: hiVal, keptCount: hiIdx + 1, totalCount: total)
    }
    
    // Results
    @Published public var elasticConstants: AcousticElasticConstantsModel = AcousticElasticConstantsModel()
    @Published public var points: [AcousticPointModel] = []
    @Published public var summary: AcousticStudySummary = AcousticStudySummary()
    
    // Precalculated Visual Geometry (Marching Squares Contours + Filled Surface Mesh)
    @Published public var zeroContourSegments: [ContourSegment] = []
    @Published public var surfaceCells: [SurfaceCell] = []
    
    // Process State
    @Published public var isRunning: Bool = false
    @Published public var progress: Double = 0.0
    @Published public var statusMessage: String = "Ready"
    @Published public var logLines: [String] = []
    @Published public var lastLoadedFile: String = ""
    
    private var process: Process?
    private var cancellables = Set<AnyCancellable>()
    
    public init() {
        // Automatically recompute linear elastic constants when inputs change
        Publishers.CombineLatest4($potentialMode, $latticeSymmetry, $beta, $bulkK)
            .debounce(for: .milliseconds(150), scheduler: RunLoop.main)
            .sink { [weak self] _ in
                self?.refreshElasticConstants()
            }
            .store(in: &cancellables)
            
        // Initial computation
        refreshElasticConstants()
        
        // Auto-load previous or default results if present
        let root = SimulationConfig.detectProjectRoot()
        let guiJSON = "\(root)/scratch/acoustic_study_gui.json"
        let defaultJSON = "\(root)/scratch/acoustic_study_results.json"
        if FileManager.default.fileExists(atPath: guiJSON) {
            loadResultsFromJSON(path: guiJSON)
        } else if FileManager.default.fileExists(atPath: defaultJSON) {
            loadResultsFromJSON(path: defaultJSON)
        }
    }
    
    public func refreshElasticConstants() {
        // Compute natively or via binary
        computeElasticConstantsFast()
    }
    
    public func computeElasticConstantsFast() {
        // Quick native evaluation for Zanzotto or background CLI query
        let root = SimulationConfig.detectProjectRoot()
        let binary = SimulationConfig.detectBinaryPath(projectRoot: root)
        
        guard FileManager.default.fileExists(atPath: binary) else {
            // Fallback native calculation
            computeNativeElasticConstants()
            return
        }
        
        let p = Process()
        p.executableURL = URL(fileURLWithPath: binary)
        p.arguments = [
            "--acoustic-study",
            "--print-elastic-constants",
            "--potential-mode=\(potentialMode.rawValue)",
            "--lattice=\(latticeSymmetry.rawValue)",
            "--beta=\(beta)",
            "--K=\(bulkK)"
        ]
        
        let pipe = Pipe()
        p.standardOutput = pipe
        p.standardError = Pipe()
        
        DispatchQueue.global(qos: .userInitiated).async {
            do {
                try p.run()
                p.waitUntilExit()
                let data = pipe.fileHandleForReading.readDataToEndOfFile()
                if let str = String(data: data, encoding: .utf8) {
                    // Extract JSON block between { and }
                    if let startIdx = str.range(of: "{\n")?.lowerBound,
                       let endIdx = str.range(of: "\n}", options: .backwards)?.upperBound {
                        let jsonSub = String(str[startIdx..<endIdx])
                        if let jsonData = jsonSub.data(using: .utf8),
                           let parsed = try? JSONDecoder().decode(AcousticElasticConstantsModel.self, from: jsonData) {
                            DispatchQueue.main.async {
                                self.elasticConstants = parsed
                            }
                            return
                        }
                    }
                }
            } catch {
                // Ignore and fall back
            }
            DispatchQueue.main.async {
                self.computeNativeElasticConstants()
            }
        }
    }
    
    private func computeNativeElasticConstants() {
        var ec = AcousticElasticConstantsModel()
        
        if potentialMode == .zanzotto {
            ec.equilibrium_lattice_parameter = (latticeSymmetry == .triangular) ? pow(4.0 / 3.0, 0.25) : 1.0
            // Zanzotto square ground state
            if latticeSymmetry == .square {
                let C11 = 4.0 * (1.0 + 3.0 * beta + 0.25 * bulkK)
                let C22 = C11
                let C12 = 4.0 * (0.25 * bulkK - 1.0 - beta)
                let C66 = 1.0 + beta
                let diff = C11 - C12
                ec.C11 = C11
                ec.C22 = C22
                ec.C12 = C12
                ec.C66 = C66
                ec.Zener_A = abs(diff) > 1e-10 ? (2.0 * C66 / diff) : 1.0
            } else {
                // Zanzotto triangular ground state
                ec.C11 = 20.0 - 0.05
                ec.C22 = ec.C11
                ec.C12 = 20.0 + 0.05
                ec.C66 = -0.05
                ec.Zener_A = 1.0
            }
        } else if potentialMode == .ljV2 {
            if latticeSymmetry == .triangular {
                ec.equilibrium_lattice_parameter = 0.687204426704008
                ec.C11 = 110.936
                ec.C22 = 110.936
                ec.C12 = 36.9788
                ec.C66 = 36.9788
                ec.Zener_A = 1.0
            } else {
                ec.equilibrium_lattice_parameter = 0.692978556224715
                ec.C11 = 83.1502
                ec.C22 = 83.1502
                ec.C12 = -5.8797
                ec.C66 = -5.8797
                ec.Zener_A = -0.132
            }
        } else if potentialMode == .square {
            // Mode 4 Square Gaussian
            ec.equilibrium_lattice_parameter = (latticeSymmetry == .triangular) ? 1.0 : 1.065846222482338
            ec.C11 = 79.875
            ec.C22 = 79.875
            ec.C12 = 10.874
            ec.C66 = 10.874
            ec.Zener_A = 0.315
        } else {
            // Mode 2 LJ v3
            ec.equilibrium_lattice_parameter = (latticeSymmetry == .triangular) ? 0.996407177045697 : 0.985415026120263
            ec.C11 = 42.0
            ec.C22 = 42.0
            ec.C12 = 14.0
            ec.C66 = 14.0
            ec.Zener_A = 1.0
        }
        
        ec.bulk_modulus = 0.25 * (ec.C11 + ec.C22 + 2.0 * ec.C12)
        ec.shear_modulus = ec.C66
        ec.born_stable = (ec.C11 > 0) && (ec.C22 > 0) && (ec.C11 * ec.C22 - ec.C12 * ec.C12 > 0) && (ec.C66 > 0)
        ec.is_isotropic = (abs(ec.Zener_A - 1.0) < 0.01)
        ec.populateDislocationParamsFallback()
        
        if ec.is_isotropic {
            ec.anisotropy_comment = "Elastically isotropic in linear regime (Zener ratio A ≈ 1.000). Wave speeds are invariant to planar propagation angle."
        } else if ec.Zener_A > 1.0 {
            ec.anisotropy_comment = String(format: "Elastically anisotropic (Zener ratio A = %.3f > 1). Crystal is softer in shear along principal <10> axes.", ec.Zener_A)
        } else {
            ec.anisotropy_comment = String(format: "Elastically anisotropic (Zener ratio A = %.3f < 1). Crystal is softer in shear along diagonal <11> directions.", ec.Zener_A)
        }
        
        self.elasticConstants = ec
    }
    
    public func runParametricStudy() {
        guard !isRunning else { return }
        
        let root = SimulationConfig.detectProjectRoot()
        let binary = SimulationConfig.detectBinaryPath(projectRoot: root)
        guard FileManager.default.fileExists(atPath: binary) else {
            statusMessage = "Binary not found at \(binary). Build project first."
            return
        }
        
        try? FileManager.default.createDirectory(atPath: "\(root)/scratch", withIntermediateDirectories: true)
        let datOut = "\(root)/scratch/acoustic_study_gui.dat"
        let jsonOut = "\(root)/scratch/acoustic_study_gui.json"
        let csvOut = "\(root)/scratch/energy_path_gui.csv"
        
        let p = Process()
        p.executableURL = URL(fileURLWithPath: binary)
        p.currentDirectoryURL = URL(fileURLWithPath: root)
        p.arguments = [
            "--acoustic-study",
            "--potential-mode=\(potentialMode.rawValue)",
            "--lattice=\(latticeSymmetry.rawValue)",
            "--beta=\(beta)",
            "--K=\(bulkK)",
            "--t-max=\(tMax)",
            "--dt=\(dt)",
            "--dp=\(dp)",
            "--shear-path-type=\(pathType)",
            "--shear-theta=\(shearTheta)",
            "--shear-alpha-min=\(shearAlphaMin)",
            "--shear-alpha-max=\(shearAlphaMax)",
            "--shear-dalpha=\(shearDAlpha)",
            "--dat-output=\(datOut)",
            "--json-output=\(jsonOut)",
            "--csv-output=\(csvOut)"
        ]
        
        let pipe = Pipe()
        p.standardOutput = pipe
        p.standardError = pipe
        
        self.isRunning = true
        self.progress = 0.0
        self.statusMessage = "Running parametric acoustic sweep..."
        self.logLines.removeAll()
        self.appendLog("🚀 Starting Parametric Acoustic Tensor Study...")
        self.appendLog("Parameters: Mode=\(potentialMode.shortName), Lattice=\(latticeSymmetry.rawValue), β=\(beta), K=\(bulkK)")
        self.appendLog("Grid: t_max=\(tMax), dt=\(dt), dp=\(dp)")
        
        self.process = p
        
        pipe.fileHandleForReading.readabilityHandler = { [weak self] handle in
            let data = handle.availableData
            guard !data.isEmpty, let text = String(data: data, encoding: .utf8) else { return }
            
            DispatchQueue.main.async {
                self?.handleProcessOutput(text)
            }
        }
        
        p.terminationHandler = { [weak self] proc in
            DispatchQueue.main.async {
                self?.isRunning = false
                pipe.fileHandleForReading.readabilityHandler = nil
                self?.process = nil
                
                if proc.terminationStatus == 0 {
                    self?.statusMessage = "Sweep completed successfully!"
                    self?.appendLog("✓ Sweep completed with exit code 0.")
                    self?.loadResultsFromJSON(path: jsonOut)
                } else if proc.terminationReason == .uncaughtSignal {
                    self?.statusMessage = "Stopped by user"
                    self?.appendLog("⏹ Process terminated.")
                } else {
                    self?.statusMessage = "Process exited with code \(proc.terminationStatus)"
                    self?.appendLog("❌ Failed with code \(proc.terminationStatus)")
                }
            }
        }
        
        do {
            try p.run()
        } catch {
            self.isRunning = false
            self.statusMessage = "Failed to launch: \(error.localizedDescription)"
            self.appendLog("❌ Error: \(error.localizedDescription)")
        }
    }
    
    public func stopStudy() {
        process?.terminate()
        isRunning = false
        statusMessage = "Stopped by user"
        appendLog("⏹ Process stopped by user.")
    }
    
    private func handleProcessOutput(_ text: String) {
        let lines = text.split(separator: "\n").map(String.init)
        for line in lines {
            appendLog(line)
            
            // Check for progress line: "Computed 250 points..."
            if line.contains("Computed") && line.contains("points") {
                let parts = line.split(separator: " ")
                if parts.count >= 2, let current = Double(parts[1]) {
                    let estTotal = (tMax / dt + 1.0) * (2.0 * .pi / dp + 1.0)
                    self.progress = min(0.95, current / max(1.0, estTotal))
                    self.statusMessage = "Progress: \(Int(current)) points calculated"
                }
            }
        }
    }
    
    private func appendLog(_ text: String) {
        logLines.append(text)
        if logLines.count > 1000 {
            logLines.removeFirst(100)
        }
    }
    
    public func loadResultsFromJSON(path: String) {
        guard let data = try? Data(contentsOf: URL(fileURLWithPath: path)) else {
            appendLog("⚠️ Could not read JSON results from \(path)")
            return
        }
        
        do {
            let root = try JSONDecoder().decode(AcousticStudyJSONRoot.self, from: data)
            if let ec = root.elastic_constants {
                self.elasticConstants = ec
            }
            if let sum = root.summary {
                self.summary = sum
            }
            self.shearPath = root.shear_path
            if let sp = root.shear_path {
                let list = sp.critical_alphas.map { String(format: "%.5f", $0) }.joined(separator: ", ")
                self.appendLog("✓ Simple shear path: \(sp.points.count) points, critical α = [\(list)]")
            }
            if let pts = root.points {
                self.points = pts
                self.progress = 1.0
                self.lastLoadedFile = (path as NSString).lastPathComponent
                self.generateMeshAndContours()
                self.appendLog("✓ Successfully loaded \(pts.count) Poincaré stability points from JSON (Mesh: \(surfaceCells.count) cells, Contours: \(zeroContourSegments.count) segments).")
            }
        } catch {
            self.appendLog("❌ JSON decode error: \(error.localizedDescription)")
        }
    }
    
    public func loadDatFile(url: URL) {
        do {
            let text = try String(contentsOf: url, encoding: .utf8)
            let lines = text.split(separator: "\n")
            
            var loadedPoints: [AcousticPointModel] = []
            var stableCount = 0
            var unstableCount = 0
            var minDet = 1e20
            var maxDet = -1e20
            var minEnergy = 1e20
            var maxEnergy = -1e20
            
            for line in lines {
                let trimmed = line.trimmingCharacters(in: .whitespaces)
                if trimmed.isEmpty || trimmed.hasPrefix("#") { continue }
                
                let tokens = trimmed.split(whereSeparator: { $0.isWhitespace }).map(String.init)
                guard tokens.count >= 9 else { continue }
                
                guard let t = Double(tokens[0]),
                      let p = Double(tokens[1]),
                      let c11 = Double(tokens[2]),
                      let c22 = Double(tokens[3]),
                      let c12 = Double(tokens[4]),
                      let c11_red = Double(tokens[5]),
                      let c22_red = Double(tokens[6]),
                      let c12_red = Double(tokens[7]),
                      let min_detAc = Double(tokens[8]) else { continue }
                      
                let angle_deg = tokens.count > 9 ? (Double(tokens[9]) ?? 0.0) : 0.0
                let third_cond = tokens.count > 10 ? (tokens[10] == "1") : false
                let strainEnergy = tokens.count > 15 ? (Double(tokens[15]) ?? 0.0) : 0.0
                
                // Stereographic projection formula requested by user:
                // t_proj = 2. / (2. + c11 + c22)
                // x = t_proj * (c11 - c22) / 2.
                // y = t_proj * c12
                let t_proj = 2.0 / (2.0 + c11 + c22)
                let x = t_proj * (c11 - c22) * 0.5
                let y = t_proj * c12
                
                let t_proj_red = 2.0 / (2.0 + c11_red + c22_red)
                let x_red = t_proj_red * (c11_red - c22_red) * 0.5
                let y_red = t_proj_red * c12_red
                
                let isStable = min_detAc > 0.0
                if isStable { stableCount += 1 } else { unstableCount += 1 }
                minDet = min(minDet, min_detAc)
                maxDet = max(maxDet, min_detAc)
                minEnergy = min(minEnergy, strainEnergy)
                maxEnergy = max(maxEnergy, strainEnergy)
                
                let pt = AcousticPointModel(
                    t: t, p: p,
                    c11: c11, c22: c22, c12: c12,
                    c11_red: c11_red, c22_red: c22_red, c12_red: c12_red,
                    min_detAc: min_detAc, energy: strainEnergy, angle_deg: angle_deg, third_condition: third_cond,
                    x: x, y: y, x_red: x_red, y_red: y_red, is_stable: isStable
                )
                loadedPoints.append(pt)
            }
            
            self.points = loadedPoints
            var sum = AcousticStudySummary()
            sum.total_points = loadedPoints.count
            sum.stable_points = stableCount
            sum.unstable_points = unstableCount
            sum.min_detAc_min = minDet
            sum.min_detAc_max = maxDet
            sum.energy_min = minEnergy < 1e19 ? minEnergy : 0.0
            sum.energy_max = maxEnergy > -1e19 ? maxEnergy : 0.0
            self.summary = sum
            self.lastLoadedFile = url.lastPathComponent
            self.generateMeshAndContours()
            self.statusMessage = "Loaded \(loadedPoints.count) points from \(url.lastPathComponent)"
            self.appendLog("✓ Loaded \(loadedPoints.count) points from DAT file: \(url.path)")
        } catch {
            self.statusMessage = "Error opening file: \(error.localizedDescription)"
            self.appendLog("❌ Failed to open DAT file: \(error.localizedDescription)")
        }
    }
    
    // MARK: - Smooth Surface Mesh & Marching Squares Zero-Boundary Extraction
    public func generateMeshAndContours() {
        self.sortedEnergies = points.map { $0.energy }.filter { $0.isFinite }.sorted()
        if !sortedEnergies.isEmpty {
            self.customEnergyMin = sortedEnergies.first ?? 0.0
            let idx95 = min(sortedEnergies.count - 1, Int(Double(sortedEnergies.count - 1) * 0.95))
            self.customEnergyMax = sortedEnergies[idx95]
        }
        
        guard !points.isEmpty else {
            self.zeroContourSegments = []
            self.surfaceCells = []
            return
        }
        
        // Group points by radius t (using rounding to 4 decimals to avoid float jitter)
        var tGroups = [Int: [AcousticPointModel]]()
        for pt in points {
            let tKey = Int((pt.t * 10000.0).rounded())
            tGroups[tKey, default: []].append(pt)
        }
        
        let sortedTKeys = tGroups.keys.sorted()
        guard sortedTKeys.count >= 2 else {
            self.zeroContourSegments = []
            self.surfaceCells = []
            return
        }
        
        // Sort each radial ring by angle p
        var rings: [[AcousticPointModel]] = []
        for k in sortedTKeys {
            let ring = tGroups[k]!.sorted { $0.p < $1.p }
            rings.append(ring)
        }
        
        var newSegments: [ContourSegment] = []
        var newCells: [SurfaceCell] = []
        
        // Helper: linear interpolation between two sample points for exact zero-crossing
        func interpEdge(_ a: AcousticPointModel, _ b: AcousticPointModel) -> (orig: CGPoint, red: CGPoint) {
            let denom = b.min_detAc - a.min_detAc
            let w = abs(denom) > 1e-12 ? max(0.0, min(1.0, (0.0 - a.min_detAc) / denom)) : 0.5
            let ox = CGFloat((1.0 - w) * a.x + w * b.x)
            let oy = CGFloat((1.0 - w) * a.y + w * b.y)
            let rx = CGFloat((1.0 - w) * a.x_red + w * b.x_red)
            let ry = CGFloat((1.0 - w) * a.y_red + w * b.y_red)
            return (CGPoint(x: ox, y: oy), CGPoint(x: rx, y: ry))
        }
        
        // Connect cells between ring i and ring i + 1
        for r in 0 ..< (rings.count - 1) {
            let r0 = rings[r]
            let r1 = rings[r + 1]
            guard !r0.isEmpty, !r1.isEmpty else { continue }
            
            let isCenter = (sortedTKeys[r] == 0)
            
            if isCenter {
                // Center triangle fan: from origin to first ring
                let centerPt = r0[0]
                for j in 0 ..< r1.count {
                    let jNext = (j + 1) % r1.count
                    let p1 = r1[j]
                    let p2 = r1[jNext]
                    
                    let avgDet = (centerPt.min_detAc + p1.min_detAc + p2.min_detAc) / 3.0
                    let avgE = (centerPt.energy + p1.energy + p2.energy) / 3.0
                    let maxE = max(centerPt.energy, max(p1.energy, p2.energy))
                    
                    newCells.append(SurfaceCell(
                        ptsOrig: [CGPoint(x: centerPt.x, y: centerPt.y),
                                  CGPoint(x: p1.x, y: p1.y),
                                  CGPoint(x: p2.x, y: p2.y)],
                        ptsRed: [CGPoint(x: centerPt.x_red, y: centerPt.y_red),
                                 CGPoint(x: p1.x_red, y: p1.y_red),
                                 CGPoint(x: p2.x_red, y: p2.y_red)],
                        avgDetAc: avgDet,
                        avgEnergy: avgE,
                        maxEnergy: maxE
                    ))
                    
                    let z0 = centerPt.min_detAc
                    let z1 = p1.min_detAc
                    let z2 = p2.min_detAc
                    
                    var cuts: [(orig: CGPoint, red: CGPoint)] = []
                    if (z0 > 0) != (z1 > 0) { cuts.append(interpEdge(centerPt, p1)) }
                    if (z1 > 0) != (z2 > 0) { cuts.append(interpEdge(p1, p2)) }
                    if (z2 > 0) != (z0 > 0) { cuts.append(interpEdge(p2, centerPt)) }
                    
                    if cuts.count == 2 {
                        newSegments.append(ContourSegment(
                            startX: cuts[0].orig.x, startY: cuts[0].orig.y,
                            endX: cuts[1].orig.x, endY: cuts[1].orig.y,
                            startXRed: cuts[0].red.x, startYRed: cuts[0].red.y,
                            endXRed: cuts[1].red.x, endYRed: cuts[1].red.y
                        ))
                    }
                }
            } else {
                // Quad cells between ring r0 and ring r1
                let count0 = r0.count
                let count1 = r1.count
                let nSlices = min(count0, count1)
                
                for j in 0 ..< nSlices {
                    let jNext0 = (j + 1) % count0
                    let jNext1 = (j + 1) % count1
                    
                    let p0 = r0[j]
                    let p1 = r0[jNext0]
                    let p2 = r1[jNext1]
                    let p3 = r1[j]
                    
                    let avgDet = (p0.min_detAc + p1.min_detAc + p2.min_detAc + p3.min_detAc) * 0.25
                    let avgE = (p0.energy + p1.energy + p2.energy + p3.energy) * 0.25
                    let maxE = max(max(p0.energy, p1.energy), max(p2.energy, p3.energy))
                    
                    newCells.append(SurfaceCell(
                        ptsOrig: [CGPoint(x: p0.x, y: p0.y),
                                  CGPoint(x: p1.x, y: p1.y),
                                  CGPoint(x: p2.x, y: p2.y),
                                  CGPoint(x: p3.x, y: p3.y)],
                        ptsRed: [CGPoint(x: p0.x_red, y: p0.y_red),
                                 CGPoint(x: p1.x_red, y: p1.y_red),
                                 CGPoint(x: p2.x_red, y: p2.y_red),
                                 CGPoint(x: p3.x_red, y: p3.y_red)],
                        avgDetAc: avgDet,
                        avgEnergy: avgE,
                        maxEnergy: maxE
                    ))
                    
                    // Marching Squares on the 4 corners:
                    let z0 = p0.min_detAc
                    let z1 = p1.min_detAc
                    let z2 = p2.min_detAc
                    let z3 = p3.min_detAc
                    
                    let mask = (z0 > 0 ? 1 : 0) | (z1 > 0 ? 2 : 0) | (z2 > 0 ? 4 : 0) | (z3 > 0 ? 8 : 0)
                    guard mask > 0 && mask < 15 else { continue }
                    
                    var e0: (orig: CGPoint, red: CGPoint)? = nil
                    var e1: (orig: CGPoint, red: CGPoint)? = nil
                    var e2: (orig: CGPoint, red: CGPoint)? = nil
                    var e3: (orig: CGPoint, red: CGPoint)? = nil
                    
                    if (z0 > 0) != (z1 > 0) { e0 = interpEdge(p0, p1) }
                    if (z1 > 0) != (z2 > 0) { e1 = interpEdge(p1, p2) }
                    if (z2 > 0) != (z3 > 0) { e2 = interpEdge(p2, p3) }
                    if (z3 > 0) != (z0 > 0) { e3 = interpEdge(p3, p0) }
                    
                    func addSeg(_ a: (orig: CGPoint, red: CGPoint)?, _ b: (orig: CGPoint, red: CGPoint)?) {
                        guard let a = a, let b = b else { return }
                        newSegments.append(ContourSegment(
                            startX: a.orig.x, startY: a.orig.y,
                            endX: b.orig.x, endY: b.orig.y,
                            startXRed: a.red.x, startYRed: a.red.y,
                            endXRed: b.red.x, endYRed: b.red.y
                        ))
                    }
                    
                    switch mask {
                    case 1, 14: addSeg(e3, e0)
                    case 2, 13: addSeg(e0, e1)
                    case 3, 12: addSeg(e3, e1)
                    case 4, 11: addSeg(e1, e2)
                    case 5:
                        addSeg(e3, e0)
                        addSeg(e1, e2)
                    case 6, 9:  addSeg(e0, e2)
                    case 7, 8:  addSeg(e2, e3)
                    case 10:
                        addSeg(e0, e1)
                        addSeg(e2, e3)
                    default: break
                    }
                }
            }
        }
        
        self.surfaceCells = newCells
        self.zeroContourSegments = newSegments
    }
    
    // MARK: - Data Export Generators
    
    /// Exports the exact zero-contour stability boundary line segments and near-boundary points as CSV.
    public func exportStabilityBoundariesCSV() -> String {
        var csv = "# Acoustic Stability Boundary (min det Ac = 0 Loss of Strong Ellipticity)\n"
        csv += "# Symmetry: \(latticeSymmetry.displayName), Potential Mode: \(potentialMode.displayName)\n"
        csv += "# Elastic Moduli: C11=\(elasticConstants.C11), C22=\(elasticConstants.C22), C12=\(elasticConstants.C12), C66=\(elasticConstants.C66)\n"
        csv += "type,segment_id,point_id,poincare_x,poincare_y,poincare_x_red,poincare_y_red,t,p,c11,c22,c12,min_detAc,critical_angle_deg,is_stable\n"
        
        // 1. Boundary Contour Line Vertices
        for (idx, seg) in zeroContourSegments.enumerated() {
            csv += "contour_vertex,\(idx),start,\(seg.startX),\(seg.startY),\(seg.startXRed),\(seg.startYRed),0,0,0,0,0,0.0,0,false\n"
            csv += "contour_vertex,\(idx),end,\(seg.endX),\(seg.endY),\(seg.endXRed),\(seg.endYRed),0,0,0,0,0,0.0,0,false\n"
        }
        
        // 2. Near Boundary Sampling Points (|min_detAc| <= 0.08)
        let nearPoints = points.filter { abs($0.min_detAc) <= 0.08 }
        for (idx, pt) in nearPoints.enumerated() {
            csv += "sampled_boundary,\(idx),\(pt.id),\(pt.x),\(pt.y),\(pt.x_red),\(pt.y_red),\(pt.t),\(pt.p),\(pt.c11),\(pt.c22),\(pt.c12),\(pt.min_detAc),\(pt.angle_deg),\(pt.is_stable)\n"
        }
        
        return csv
    }
    
    /// Exports the full parametric acoustic stability scan dataset as CSV.
    public func exportFullParametricScanCSV() -> String {
        var csv = "point_id,poincare_x,poincare_y,poincare_x_red,poincare_y_red,t,p,c11,c22,c12,c11_red,c22_red,c12_red,min_detAc,critical_angle_deg,energy,is_stable\n"
        for (idx, pt) in points.enumerated() {
            csv += "\(idx),\(pt.x),\(pt.y),\(pt.x_red),\(pt.y_red),\(pt.t),\(pt.p),\(pt.c11),\(pt.c22),\(pt.c12),\(pt.c11_red),\(pt.c22_red),\(pt.c12_red),\(pt.min_detAc),\(pt.angle_deg),\(pt.energy),\(pt.is_stable)\n"
        }
        return csv
    }
    
    /// Exports the simple shear path and critical transition points as CSV.
    public func exportShearPathCSV() -> String {
        guard let sp = shearPath else { return "No shear path computed.\n" }
        var csv = "# Simple Shear Path Stability\n"
        csv += "# Critical points count: \(sp.critical_points.count)\n"
        csv += "type,alpha,poincare_x,poincare_y,poincare_x_red,poincare_y_red,c11,c22,c12,min_detAc,energy,is_stable\n"
        for pt in sp.critical_points {
            let detStr = pt.min_detAc.map { String($0) } ?? "nan"
            let stabStr = pt.is_stable.map { String($0) } ?? "unknown"
            csv += "critical_transition,\(pt.alpha),\(pt.x),\(pt.y),\(pt.x_red),\(pt.y_red),\(pt.c11),\(pt.c22),\(pt.c12),\(detStr),\(pt.energy),\(stabStr)\n"
        }
        for pt in sp.points {
            let detStr = pt.min_detAc.map { String($0) } ?? "nan"
            let stabStr = pt.is_stable.map { String($0) } ?? "unknown"
            csv += "shear_step,\(pt.alpha),\(pt.x),\(pt.y),\(pt.x_red),\(pt.y_red),\(pt.c11),\(pt.c22),\(pt.c12),\(detStr),\(pt.energy),\(stabStr)\n"
        }
        return csv
    }
}

