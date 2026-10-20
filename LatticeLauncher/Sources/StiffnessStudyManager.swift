import Foundation
import Combine
import SwiftUI
import AppKit

public enum BoundaryConditionMode: String, CaseIterable, Identifiable {
    case periodic = "periodic"
    case fixed = "fixed"
    
    public var id: String { rawValue }
    public var displayName: String {
        switch self {
        case .periodic: return "Periodic (PBC)"
        case .fixed: return "Fixed Boundary"
        }
    }
    public var icon: String {
        switch self {
        case .periodic: return "arrow.triangle.2.circlepath"
        case .fixed: return "lock.fill"
        }
    }
}

public struct StiffnessMeshModel: Codable {
    public var a_0: Double = 1.0
    public var nodes_x: [Double] = []
    public var nodes_y: [Double] = []
    public var is_free: [Int] = []
    public var triangles: [[Int]] = []
    
    public init() {}
}

public struct StiffnessPointModel: Identifiable, Codable {
    public var id: UUID = UUID()
    public var alpha: Double = 0.0
    public var lambda_min: Double = 0.0
    public var lambda_2: Double = 0.0
    public var lambda_3: Double = 0.0
    public var is_stable: Bool = true
    public var energy: Double = 0.0
    public var energy_density: Double = 0.0
    public var continuum_energy: Double = 0.0
    public var mode1_ux: [Double]? = nil
    public var mode1_uy: [Double]? = nil
    public var mode2_ux: [Double]? = nil
    public var mode2_uy: [Double]? = nil
    public var mode3_ux: [Double]? = nil
    public var mode3_uy: [Double]? = nil
    
    enum CodingKeys: String, CodingKey {
        case alpha, lambda_min, lambda_2, lambda_3, is_stable, energy, energy_density, continuum_energy
        case mode1_ux, mode1_uy, mode2_ux, mode2_uy, mode3_ux, mode3_uy
    }
    
    public init() {}
    
    public init(from decoder: Decoder) throws {
        let container = try decoder.container(keyedBy: CodingKeys.self)
        self.id = UUID()
        self.alpha = try container.decode(Double.self, forKey: .alpha)
        self.lambda_min = try container.decode(Double.self, forKey: .lambda_min)
        self.lambda_2 = try container.decode(Double.self, forKey: .lambda_2)
        self.lambda_3 = try container.decode(Double.self, forKey: .lambda_3)
        self.is_stable = try container.decode(Bool.self, forKey: .is_stable)
        self.energy = try container.decode(Double.self, forKey: .energy)
        self.energy_density = try container.decodeIfPresent(Double.self, forKey: .energy_density) ?? energy
        self.continuum_energy = try container.decodeIfPresent(Double.self, forKey: .continuum_energy) ?? energy_density
        self.mode1_ux = try container.decodeIfPresent([Double].self, forKey: .mode1_ux)
        self.mode1_uy = try container.decodeIfPresent([Double].self, forKey: .mode1_uy)
        self.mode2_ux = try container.decodeIfPresent([Double].self, forKey: .mode2_ux)
        self.mode2_uy = try container.decodeIfPresent([Double].self, forKey: .mode2_uy)
        self.mode3_ux = try container.decodeIfPresent([Double].self, forKey: .mode3_ux)
        self.mode3_uy = try container.decodeIfPresent([Double].self, forKey: .mode3_uy)
    }
}

public struct StiffnessStudySummaryModel: Codable {
    public var total_nodes: Int = 0
    public var free_nodes: Int = 0
    public var total_dofs: Int = 0
    public var total_elements: Int = 0
    public var alpha_c_pos: Double = 0.0
    public var alpha_c_neg: Double = 0.0
    public var min_lambda: Double = 0.0
    public var max_lambda: Double = 0.0
    public var execution_seconds: Double = 0.0
}

public struct StiffnessStudyConfigModel: Codable {
    public var nx: Int?
    public var ny: Int?
    public var lattice_type: String?
    public var potential_mode: Int?
    public var bc: String?
    public var alpha_min: Double?
    public var alpha_max: Double?
    public var dalpha: Double?
    public var beta: Double?
    public var K: Double?
}

public struct StiffnessFullScaleSnapshot: Codable {
    public var alpha: Double = 0.0
    public var lambda_min: Double = 0.0
    public var lambda_2: Double = 0.0
    public var lambda_3: Double = 0.0
    public var is_stable: Bool = true
    public var mode1_ux: [Double]? = nil
    public var mode1_uy: [Double]? = nil
    public var mode2_ux: [Double]? = nil
    public var mode2_uy: [Double]? = nil
    public var mode3_ux: [Double]? = nil
    public var mode3_uy: [Double]? = nil
}

public struct StiffnessFullScaleWrapper: Codable {
    public var nx: Int = 10
    public var ny: Int = 10
    public var a_0: Double = 1.0
    public var lattice_type: String = "square"
    public var bc: String = "periodic"
    public var mesh: StiffnessMeshModel? = nil
    public var snapshots: [StiffnessFullScaleSnapshot] = []
}

public struct StiffnessStudyJSONWrapper: Codable {
    public var config: StiffnessStudyConfigModel?
    public var summary: StiffnessStudySummaryModel
    public var mesh: StiffnessMeshModel?
    public var points: [StiffnessPointModel]
}

public class StiffnessStudyManager: ObservableObject {
    // Inputs
    @Published public var nx: Int = 10
    @Published public var ny: Int = 10
    @Published public var latticeSymmetry: LatticeSymmetry = .square
    @Published public var potentialMode: AcousticPotentialMode = .zanzotto
    @Published public var boundaryCondition: BoundaryConditionMode = .periodic
    
    // Shear Range
    @Published public var alphaMin: Double = -0.3
    @Published public var alphaMax: Double = 0.3
    @Published public var dalpha: Double = 0.005
    
    // Zanzotto Mode 3 parameters
    @Published public var beta: Double = -0.25
    @Published public var bulkK: Double = 5.0
    
    // Process State
    @Published public var isRunning: Bool = false
    @Published public var isComputingFullScale: Bool = false
    @Published public var progress: Double = 0.0
    @Published public var completedSteps: Int = 0
    @Published public var totalSteps: Int = 0
    @Published public var elapsedSeconds: Double = 0.0
    @Published public var etaSeconds: Double? = nil
    @Published public var statusMessage: String = "Ready"
    @Published public var logLines: [String] = []
    
    // Results
    @Published public var points: [StiffnessPointModel] = []
    @Published public var summary: StiffnessStudySummaryModel = StiffnessStudySummaryModel()
    @Published public var mesh: StiffnessMeshModel? = nil
    @Published public var fullScaleData: StiffnessFullScaleWrapper? = nil
    @Published public var selectedAlpha: Double? = nil
    
    private var process: Process?
    
    public init() {
        let root = SimulationConfig.detectProjectRoot()
        let jsonOut = "\(root)/scratch/stiffness_shear_study.json"
        if FileManager.default.fileExists(atPath: jsonOut) {
            loadResultsJSON(jsonPath: jsonOut, restoreConfig: true)
        }
        let fsOut = "\(root)/scratch/stiffness_fullscale.json"
        if FileManager.default.fileExists(atPath: fsOut) {
            loadFullScaleJSON(jsonPath: fsOut, mergeExisting: false)
        }
    }
    
    public func runStudy() {
        guard !isRunning else { return }
        
        let root = SimulationConfig.detectProjectRoot()
        let binary = "\(root)/build/lattice_triangulation"
        guard FileManager.default.fileExists(atPath: binary) else {
            statusMessage = "Binary not found at \(binary). Build project first."
            return
        }
        
        guard dalpha > 1e-9 else {
            statusMessage = "Error: Step Δα must be greater than 0"
            appendLog("[ERROR] Step Δα must be greater than 0")
            return
        }
        guard alphaMax >= alphaMin else {
            statusMessage = "Error: α Max must be greater than or equal to α Min"
            appendLog("[ERROR] α Max must be greater than or equal to α Min")
            return
        }
        
        let scratchDir = "\(root)/scratch"
        try? FileManager.default.createDirectory(atPath: scratchDir, withIntermediateDirectories: true)
        let csvOut = "\(scratchDir)/stiffness_shear_study.csv"
        let jsonOut = "\(scratchDir)/stiffness_shear_study.json"
        let fsOut = "\(scratchDir)/stiffness_fullscale.json"
        
        let p = Process()
        p.executableURL = URL(fileURLWithPath: binary)
        p.currentDirectoryURL = URL(fileURLWithPath: root)
        var args = [
            "--stiffness-study",
            "--nx=\(nx)",
            "--ny=\(ny)",
            "--lattice=\(latticeSymmetry.rawValue)",
            "--potential-mode=\(potentialMode.rawValue)",
            "--bc=\(boundaryCondition.rawValue)",
            "--alpha-min=\(alphaMin)",
            "--alpha-max=\(alphaMax)",
            "--dalpha=\(dalpha)",
            "--stiffness-csv=\(csvOut)",
            "--stiffness-json=\(jsonOut)",
            "--stiffness-fullscale-json=\(fsOut)"
        ]
        if potentialMode == .zanzotto {
            args.append("--beta=\(beta)")
            args.append("--K=\(bulkK)")
        }
        p.arguments = args
        
        let pipe = Pipe()
        p.standardOutput = pipe
        p.standardError = pipe
        
        let estSteps = max(1, Int(round((alphaMax - alphaMin) / dalpha))) + 1
        self.isRunning = true
        self.progress = 0.0
        self.completedSteps = 0
        self.totalSteps = estSteps
        self.elapsedSeconds = 0.0
        self.etaSeconds = nil
        self.statusMessage = "Step 0 / \(estSteps) (assembling stiffness matrix)..."
        self.logLines.removeAll()
        self.points.removeAll()
        self.summary = StiffnessStudySummaryModel()
        
        self.appendLog("🚀 Starting Stiffness Matrix K(α) Shear Study...")
        self.appendLog("Grid: \(nx)x\(ny), Symmetry: \(latticeSymmetry.displayName), BC: \(boundaryCondition.displayName)")
        self.appendLog("Potential: \(potentialMode.displayName)")
        self.appendLog("Shear sweep: [\(alphaMin), \(alphaMax)], step \(dalpha) (\(estSteps) steps)")
        
        self.process = p
        
        pipe.fileHandleForReading.readabilityHandler = { [weak self] handle in
            let data = handle.availableData
            guard !data.isEmpty, let str = String(data: data, encoding: .utf8) else { return }
            DispatchQueue.main.async {
                let lines = str.components(separatedBy: .newlines)
                for line in lines where !line.trimmingCharacters(in: .whitespaces).isEmpty {
                    self?.handleOutputLine(line)
                }
            }
        }
        
        p.terminationHandler = { [weak self] proc in
            DispatchQueue.main.async {
                self?.isRunning = false
                pipe.fileHandleForReading.readabilityHandler = nil
                
                if proc.terminationStatus == 0 {
                    self?.progress = 1.0
                    self?.etaSeconds = 0.0
                    self?.statusMessage = "Complete! Data parsed successfully."
                    self?.loadResultsJSON(jsonPath: jsonOut)
                    self?.loadFullScaleJSON(jsonPath: fsOut, mergeExisting: false)
                } else {
                    self?.statusMessage = "Process exited with error code \(proc.terminationStatus)"
                }
            }
        }
        
        do {
            try p.run()
        } catch {
            self.isRunning = false
            self.statusMessage = "Failed to launch process: \(error.localizedDescription)"
        }
    }
    
    public func stopStudy() {
        process?.terminate()
        process = nil
        isRunning = false
        statusMessage = "Cancelled by user"
    }
    
    private func handleOutputLine(_ line: String) {
        if line.hasPrefix("[PROGRESS]") {
            // Format: "[PROGRESS] 12/121 alpha=-0.2400 elapsed=1.42"
            let parts = line.split(separator: " ")
            if parts.count >= 2 {
                let stepParts = parts[1].split(separator: "/")
                if stepParts.count == 2,
                   let cur = Int(stepParts[0]),
                   let tot = Int(stepParts[1]), tot > 0 {
                    self.completedSteps = cur
                    self.totalSteps = tot
                    self.progress = min(0.99, Double(cur) / Double(tot))
                }
            }
            for p in parts {
                if p.hasPrefix("elapsed="), let el = Double(p.dropFirst("elapsed=".count)) {
                    self.elapsedSeconds = el
                    if self.completedSteps > 0 && self.totalSteps > self.completedSteps {
                        let rem = el * Double(self.totalSteps - self.completedSteps) / Double(self.completedSteps)
                        self.etaSeconds = rem
                    } else {
                        self.etaSeconds = 0.0
                    }
                } else if p.hasPrefix("pt=") {
                    let vals = p.dropFirst("pt=".count).split(separator: ",")
                    if vals.count >= 8,
                       let a = Double(vals[0]),
                       let l1 = Double(vals[1]),
                       let l2 = Double(vals[2]),
                       let l3 = Double(vals[3]),
                       let st = Int(vals[4]),
                       let en = Double(vals[5]),
                       let ed = Double(vals[6]),
                       let ce = Double(vals[7]) {
                        var pt = StiffnessPointModel()
                        pt.alpha = a
                        pt.lambda_min = l1
                        pt.lambda_2 = l2
                        pt.lambda_3 = l3
                        pt.is_stable = st != 0
                        pt.energy = en
                        pt.energy_density = ed
                        pt.continuum_energy = ce
                        self.points.append(pt)
                    }
                }
            }
            let pct = Int(round(self.progress * 100.0))
            if let eta = self.etaSeconds, eta > 0.2 {
                self.statusMessage = String(format: "Step %d/%d (%d%%) • %.1fs elapsed • ~%.1fs left", self.completedSteps, self.totalSteps, pct, self.elapsedSeconds, eta)
            } else {
                self.statusMessage = String(format: "Step %d/%d (%d%%) • %.1fs elapsed", self.completedSteps, self.totalSteps, pct, self.elapsedSeconds)
            }
            return
        }
        appendLog(line)
    }
    
    private func appendLog(_ line: String) {
        logLines.append(line)
        if logLines.count > 1000 {
            logLines.removeFirst(logLines.count - 1000)
        }
    }
    
    private func loadResultsJSON(jsonPath: String, restoreConfig: Bool = false) {
        guard FileManager.default.fileExists(atPath: jsonPath),
              let data = try? Data(contentsOf: URL(fileURLWithPath: jsonPath)) else {
            statusMessage = "Could not read output JSON at \(jsonPath)"
            return
        }
        
        do {
            let decoder = JSONDecoder()
            let res = try decoder.decode(StiffnessStudyJSONWrapper.self, from: data)
            self.points = res.points
            self.summary = res.summary
            self.mesh = res.mesh
            if self.selectedAlpha == nil {
                if res.summary.alpha_c_pos > 0 {
                    self.selectedAlpha = res.summary.alpha_c_pos
                } else {
                    self.selectedAlpha = res.points.first(where: { !$0.is_stable })?.alpha ?? res.points.first?.alpha
                }
            }
            
            // If high-precision CSV is alongside JSON, merge its 10-digit precision scalars while keeping eigenvectors
            let csvPath = jsonPath.replacingOccurrences(of: ".json", with: ".csv")
            if FileManager.default.fileExists(atPath: csvPath) {
                loadResultsCSV(csvPath: csvPath)
            }
            
            if restoreConfig, let cfg = res.config {
                if let nxVal = cfg.nx { self.nx = nxVal }
                if let nyVal = cfg.ny { self.ny = nyVal }
                if let lat = cfg.lattice_type, let sym = LatticeSymmetry(rawValue: lat) {
                    self.latticeSymmetry = sym
                }
                if let pm = cfg.potential_mode, let mode = AcousticPotentialMode(rawValue: pm) {
                    self.potentialMode = mode
                }
                if let bcStr = cfg.bc, let bcMode = BoundaryConditionMode(rawValue: bcStr) {
                    self.boundaryCondition = bcMode
                }
                if let aMin = cfg.alpha_min { self.alphaMin = aMin }
                if let aMax = cfg.alpha_max { self.alphaMax = aMax }
                if let da = cfg.dalpha { self.dalpha = da }
                if let b = cfg.beta { self.beta = b }
                if let k = cfg.K { self.bulkK = k }
                self.statusMessage = "Loaded previous study (\(self.nx)×\(self.ny), \(self.points.count) steps)"
            }
            self.appendLog("✅ Loaded \(self.points.count) shear steps. Critical α_c: [\(String(format: "%.4f", res.summary.alpha_c_neg)), \(String(format: "%.4f", res.summary.alpha_c_pos))]")
        } catch {
            self.appendLog("⚠️ JSON decode error: \(error.localizedDescription). Parsing CSV fallback...")
            loadResultsCSV(csvPath: jsonPath.replacingOccurrences(of: ".json", with: ".csv"))
        }
    }
    
    private func loadResultsCSV(csvPath: String) {
        guard let content = try? String(contentsOfFile: csvPath, encoding: .utf8) else { return }
        var parsed: [StiffnessPointModel] = []
        let lines = content.components(separatedBy: .newlines)
        for line in lines.dropFirst() {
            let parts = line.components(separatedBy: ",").map { $0.trimmingCharacters(in: .whitespaces) }
            guard parts.count >= 6,
                  let a = Double(parts[0]),
                  let lmin = Double(parts[1]),
                  let l2 = Double(parts[2]),
                  let l3 = Double(parts[3]),
                  let st = Int(parts[4]),
                  let en = Double(parts[5]) else { continue }
            var pt = StiffnessPointModel()
            pt.alpha = a
            pt.lambda_min = lmin
            pt.lambda_2 = l2
            pt.lambda_3 = l3
            pt.is_stable = (st == 1)
            pt.energy = en
            pt.energy_density = (parts.count >= 7) ? (Double(parts[6]) ?? en) : en
            pt.continuum_energy = (parts.count >= 8) ? (Double(parts[7]) ?? pt.energy_density) : pt.energy_density
            parsed.append(pt)
        }
        if self.points.count == parsed.count {
            for i in 0..<parsed.count {
                parsed[i].mode1_ux = self.points[i].mode1_ux
                parsed[i].mode1_uy = self.points[i].mode1_uy
                parsed[i].mode2_ux = self.points[i].mode2_ux
                parsed[i].mode2_uy = self.points[i].mode2_uy
                parsed[i].mode3_ux = self.points[i].mode3_ux
                parsed[i].mode3_uy = self.points[i].mode3_uy
            }
        }
        self.points = parsed
    }
    
    public func loadFullScaleJSON(jsonPath: String, mergeExisting: Bool) {
        guard FileManager.default.fileExists(atPath: jsonPath),
              let data = try? Data(contentsOf: URL(fileURLWithPath: jsonPath)),
              let decoded = try? JSONDecoder().decode(StiffnessFullScaleWrapper.self, from: data) else {
            return
        }
        if mergeExisting, var existing = self.fullScaleData,
           existing.nx == decoded.nx, existing.ny == decoded.ny,
           existing.lattice_type == decoded.lattice_type {
            for snap in decoded.snapshots {
                existing.snapshots.removeAll(where: { abs($0.alpha - snap.alpha) < 1e-6 })
                existing.snapshots.append(snap)
            }
            existing.snapshots.sort(by: { $0.alpha < $1.alpha })
            self.fullScaleData = existing
        } else {
            self.fullScaleData = decoded
        }
    }
    
    public func computeFullScaleAtAlpha(_ alpha: Double) {
        guard !isRunning && !isComputingFullScale else { return }
        let root = SimulationConfig.detectProjectRoot()
        let binary = "\(root)/build/lattice_triangulation"
        guard FileManager.default.fileExists(atPath: binary) else { return }
        
        let scratchDir = "\(root)/scratch"
        try? FileManager.default.createDirectory(atPath: scratchDir, withIntermediateDirectories: true)
        let fsOut = "\(scratchDir)/stiffness_fullscale_single.json"
        
        let p = Process()
        p.executableURL = URL(fileURLWithPath: binary)
        p.currentDirectoryURL = URL(fileURLWithPath: root)
        var args = [
            "--stiffness-study",
            "--nx=\(nx)",
            "--ny=\(ny)",
            "--lattice=\(latticeSymmetry.rawValue)",
            "--potential-mode=\(potentialMode.rawValue)",
            "--bc=\(boundaryCondition.rawValue)",
            "--single-alpha=\(String(format: "%.6f", alpha))",
            "--stiffness-fullscale-json=\(fsOut)"
        ]
        if potentialMode == .zanzotto {
            args.append("--beta=\(beta)")
            args.append("--K=\(bulkK)")
        }
        p.arguments = args
        
        self.isComputingFullScale = true
        self.statusMessage = String(format: "Computing full %d×%d mode at α = %+.4f...", nx, ny, alpha)
        
        p.terminationHandler = { [weak self] proc in
            DispatchQueue.main.async {
                self?.isComputingFullScale = false
                if proc.terminationStatus == 0 {
                    self?.loadFullScaleJSON(jsonPath: fsOut, mergeExisting: true)
                    self?.statusMessage = String(format: "Full %d×%d mode ready at α = %+.4f", self?.nx ?? 0, self?.ny ?? 0, alpha)
                } else {
                    self?.statusMessage = "Full-scale mode computation failed (code \(proc.terminationStatus))"
                }
            }
        }
        try? p.run()
    }
}
