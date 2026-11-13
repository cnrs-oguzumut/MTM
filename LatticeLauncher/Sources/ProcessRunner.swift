import Foundation
import Combine
import AppKit

public enum ProcessStatus: Equatable {
    case idle
    case running(pid: Int32)
    case completed(exitCode: Int32, duration: TimeInterval)
    case failed(error: String)
    case terminated
    
    public var isRunning: Bool {
        if case .running = self { return true }
        return false
    }
    
    public var label: String {
        switch self {
        case .idle: return "Ready"
        case .running(let pid): return "Running (PID: \(pid))"
        case .completed(let code, _): return code == 0 ? "Completed (Success)" : "Finished with code \(code)"
        case .failed(let err): return "Failed: \(err)"
        case .terminated: return "Stopped by User"
        }
    }
}

public struct LogMessage: Identifiable {
    public let id = UUID()
    public let timestamp: Date
    public let text: String
    public let isError: Bool
    
    public init(text: String, isError: Bool = false) {
        self.timestamp = Date()
        self.text = text
        self.isError = isError
    }
}

public class ProcessRunner: ObservableObject {
    @Published public var status: ProcessStatus = .idle
    @Published public var logLines: [LogMessage] = []
    @Published public var currentStep: Int = 0
    @Published public var currentAlpha: Double = 0.0
    @Published public var currentStress: Double = 0.0
    @Published public var currentEnergy: Double = 0.0
    @Published public var avalancheCount: Int = 0
    /// Force / stress control: conjugate strain ε of the last point (nil for strain-driven runs)
    @Published public var currentConjugateStrain: Double? = nil
    @Published public var elapsedTime: TimeInterval = 0
    
    // Live Stress-Strain & Eigenvalue Curve points
    @Published public var curvePoints: [StressStrainPoint] = []
    @Published public var eigenPoints: [EigenPoint] = []
    @Published public var latticeRotationPoints: [LatticeRotationPoint] = []
    @Published public var currentLambdaMin: Double? = nil
    @Published public var peakStress: Double = 0.0
    @Published public var yieldStrain: Double = 0.0
    @Published public var activeConfig: SimulationConfig?
    
    // Grain Boundary Energy (Bicrystal alpha = 0 summary & cumulative curve vs misorientation)
    @Published public var currentGBEnergy: GBEnergyPoint? = nil
    @Published public var gbEnergyPoints: [GBEnergyPoint] = []
    @Published public var currentStageLabel: String = ""
    @Published public var isSweepRunning: Bool = false
    @Published public var sweepProgressText: String = ""
    @Published public var requestedChartMode: ChartMetricMode? = nil
    @Published public var sweepOutputDirectory: String = ""
    @Published public var sweepAngleVsEnergyFilePath: String = ""
    @Published public var sweepRelaxedVTKDirectory: String = ""
    
    // Dislocation Profiles (for Volterra mode)
    @Published public var dislocationData: DislocationProfileJSONRoot? = nil
    
    private var process: Process?
    private var timer: Timer?
    private var csvPollingTimer: Timer?
    private var startTime: Date?
    private let maxLogLines = 10_000
    @Published public var outputDirectory: String = ""
    private var pendingSweepCandidates: [CSLRationalMatch] = []
    private var currentSweepCandidate: CSLRationalMatch? = nil
    private var totalSweepCount: Int = 0
    private var currentSweepIndex: Int = 0
    private var sweepBaseTargetDir: String = ""
    private var sweepSavedBaseNx: Int = 100
    private var sweepSavedBaseNy: Int = 100
    private var pendingOutputChunks: [String] = []
    private var isOutputFlushScheduled: Bool = false
    private var lastDislocationProfileKey: String = ""
    
    public init() {
        self.pollGBEnergyCSV()
    }
    
    public func clearLogs() {
        logLines.removeAll()
        curvePoints.removeAll()
        eigenPoints.removeAll()
        latticeRotationPoints.removeAll()
        currentLambdaMin = nil
        currentGBEnergy = nil
        peakStress = 0.0
        yieldStrain = 0.0
        currentStep = 0
        currentAlpha = 0.0
        currentStress = 0.0
        currentEnergy = 0.0
        avalancheCount = 0
        currentConjugateStrain = nil
        dislocationData = nil
        lastDislocationProfileKey = ""
    }
    
    public func loadExistingIfAvailable(dir: String) {
        self.outputDirectory = dir
        let csvPath = "\(dir)/energy_stress_log.csv"
        if FileManager.default.fileExists(atPath: csvPath) {
            self.curvePoints.removeAll()
            self.eigenPoints.removeAll()
            self.latticeRotationPoints.removeAll()
            self.pollCSV()
        } else {
            self.curvePoints.removeAll()
            self.eigenPoints.removeAll()
            self.latticeRotationPoints.removeAll()
            self.currentStep = 0
            self.currentAlpha = 0.0
            self.currentStress = 0.0
            self.currentEnergy = 0.0
            self.avalancheCount = 0
            self.currentConjugateStrain = nil
        }
        self.pollGBEnergyCSV()
        self.pollDislocationProfiles(dir: dir)
    }
    
    public func openOutputFolder(fallbackDirectory: String? = nil) {
        let fm = FileManager.default
        var candidate = ""
        if !sweepOutputDirectory.isEmpty && fm.fileExists(atPath: sweepOutputDirectory) {
            candidate = sweepOutputDirectory
        } else if !outputDirectory.isEmpty && fm.fileExists(atPath: outputDirectory) {
            candidate = outputDirectory
        } else if let fb = fallbackDirectory, !fb.isEmpty {
            candidate = fb
        } else if !outputDirectory.isEmpty {
            candidate = outputDirectory
        } else {
            candidate = SimulationConfig.detectProjectRoot()
        }

        if !fm.fileExists(atPath: candidate) {
            try? fm.createDirectory(atPath: candidate, withIntermediateDirectories: true)
        }
        if !fm.fileExists(atPath: candidate) {
            candidate = SimulationConfig.detectProjectRoot()
        }

        let folderURL = URL(fileURLWithPath: candidate, isDirectory: true)
        if !NSWorkspace.shared.open(folderURL) {
            NSWorkspace.shared.selectFile(nil, inFileViewerRootedAtPath: candidate)
        }
    }

    public func openInParaView(fallbackDirectory: String? = nil, specificVTKURL: URL? = nil) -> Bool {
        let fm = FileManager.default
        var appURL: URL? = NSWorkspace.shared.urlForApplication(withBundleIdentifier: "org.paraview.ParaView")
        if appURL == nil {
            for dir in ["/Applications", NSHomeDirectory() + "/Applications"] {
                if let names = try? fm.contentsOfDirectory(atPath: dir),
                   let n = names.filter({ $0.lowercased().hasPrefix("paraview") && $0.hasSuffix(".app") }).sorted().last {
                    appURL = URL(fileURLWithPath: dir).appendingPathComponent(n)
                    break
                }
            }
        }
        guard let app = appURL else { return false }

        // Resolve target VTK file if none explicitly provided
        var targetVTK: URL? = specificVTKURL
        if targetVTK == nil {
            let baseDir: String = {
                if !sweepRelaxedVTKDirectory.isEmpty && fm.fileExists(atPath: sweepRelaxedVTKDirectory) {
                    return sweepRelaxedVTKDirectory
                }
                if !outputDirectory.isEmpty && fm.fileExists(atPath: outputDirectory) {
                    return outputDirectory
                }
                if let fb = fallbackDirectory, !fb.isEmpty && fm.fileExists(atPath: fb) {
                    return fb
                }
                return SimulationConfig.detectProjectRoot()
            }()
            let vtkDir = fm.fileExists(atPath: "\(baseDir)/vtk_output") ? "\(baseDir)/vtk_output" : baseDir
            if let items = try? fm.contentsOfDirectory(atPath: vtkDir) {
                let vtkFiles = items.filter { $0.lowercased().hasSuffix(".vtk") }.sorted()
                if let latest = vtkFiles.last {
                    targetVTK = URL(fileURLWithPath: vtkDir).appendingPathComponent(latest)
                }
            }
        }

        let execURL = app.appendingPathComponent("Contents/MacOS/paraview")
        if fm.isExecutableFile(atPath: execURL.path) {
            let proc = Process()
            proc.executableURL = execURL
            if let vtk = targetVTK {
                proc.arguments = [vtk.path]
                proc.currentDirectoryURL = vtk.deletingLastPathComponent()
            }
            do {
                try proc.run()
                return true
            } catch {
                // Fall back to NSWorkspace below
            }
        }

        let config = NSWorkspace.OpenConfiguration()
        config.activates = true
        if let vtk = targetVTK {
            config.arguments = [vtk.path]
            NSWorkspace.shared.open([vtk], withApplicationAt: app, configuration: config)
        } else {
            NSWorkspace.shared.openApplication(at: app, configuration: config)
        }
        return true
    }
    
    public func openSweepTextFile() {
        let fm = FileManager.default
        if !sweepAngleVsEnergyFilePath.isEmpty && fm.fileExists(atPath: sweepAngleVsEnergyFilePath) {
            NSWorkspace.shared.open(URL(fileURLWithPath: sweepAngleVsEnergyFilePath))
            return
        }
        let fallback1 = "\(outputDirectory)/angle_vs_energy.txt"
        if fm.fileExists(atPath: fallback1) {
            NSWorkspace.shared.open(URL(fileURLWithPath: fallback1))
            return
        }
        openOutputFolder()
    }
    
    public func openSweepRelaxedVTKFolder() {
        let fm = FileManager.default
        if !sweepRelaxedVTKDirectory.isEmpty && fm.fileExists(atPath: sweepRelaxedVTKDirectory) {
            NSWorkspace.shared.selectFile(nil, inFileViewerRootedAtPath: sweepRelaxedVTKDirectory)
            return
        }
        openOutputFolder()
    }
    
    public func start(
        config: SimulationConfig,
        overrideRunDirectory: String? = nil,
        overrideNx: Int? = nil,
        overrideNy: Int? = nil,
        overrideGeometryConfig: MacroscopicGeometryConfigModel? = nil,
        overrideProtocolMode: SimulationProtocolMode? = nil
    ) {
        guard !status.isRunning else { return }
        
        // Auto-select appropriate live chart tab
        let effProto = overrideProtocolMode ?? config.protocolMode
        if isSweepRunning || effProto == .sweepMisorientation || effProto == .relaxOnly {
            self.requestedChartMode = .gbEnergy
        }
        
        // 1. Verify binary exists
        let fm = FileManager.default
        let binary = config.binaryPath
        guard fm.fileExists(atPath: binary) else {
            self.status = .failed(error: "Binary not found at: \(binary)\nPlease compile first with cmake or make.")
            self.appendLog(text: "❌ Error: Executable does not exist at '\(binary)'", isError: true)
            return
        }
        
        // 2. Prepare target directory
        let runDir = overrideRunDirectory ?? config.targetDirectory
        self.outputDirectory = runDir
        do {
            try fm.createDirectory(atPath: runDir, withIntermediateDirectories: true, attributes: nil)
            // Remove stale CSV logs before fresh start so live chart starts clean
            try? fm.removeItem(atPath: "\(runDir)/energy_stress_log.csv")
            try? fm.removeItem(atPath: "\(runDir)/eigen_log.csv")
            try? fm.removeItem(atPath: "\(runDir)/gb_energy_summary.csv")
        } catch {
            self.status = .failed(error: "Cannot create output directory: \(error.localizedDescription)")
            return
        }
        
        // 3. Clear stats and set running state
        self.activeConfig = config
        self.curvePoints.removeAll()
        self.eigenPoints.removeAll()
        self.latticeRotationPoints.removeAll()
        self.currentLambdaMin = nil
        if !isSweepRunning {
            self.currentGBEnergy = nil
        }
        self.peakStress = 0.0
        self.yieldStrain = 0.0
        self.currentStep = 0
        self.currentAlpha = config.alphaStart
        self.currentStress = 0.0
        self.currentEnergy = 0.0
        self.avalancheCount = 0
        self.currentConjugateStrain = nil
        if !isSweepRunning {
            self.elapsedTime = 0
            self.startTime = Date()
        } else if self.startTime == nil {
            self.startTime = Date()
        }
        
        let effNx = overrideNx ?? config.nx
        let effNy = overrideNy ?? config.ny
        self.appendLog(text: "==========================================================")
        self.appendLog(text: "🚀 Launching Simulation: \(effNx)x\(effNy), mode=\(config.mode.rawValue), seed=\(config.seed)")
        self.appendLog(text: "📁 Working Directory: \(runDir)")
        self.appendLog(text: "⚙️ Command: \(config.fullCommandLine())")
        self.appendLog(text: "==========================================================")
        
        // 4. Configure Process
        let proc = Process()
        proc.executableURL = URL(fileURLWithPath: binary)
        proc.arguments = config.buildArguments(
            overrideRunDirectory: runDir,
            overrideNx: effNx,
            overrideNy: effNy,
            overrideGeometryConfig: overrideGeometryConfig,
            overrideProtocolMode: overrideProtocolMode
        )
        proc.currentDirectoryURL = URL(fileURLWithPath: runDir)
        
        proc.environment = SimulationConfig.subprocessEnvironment(
            ompThreads: config.ompThreads,
            blasThreads: config.blasThreads
        )
        
        // 5. Pipes
        let pipe = Pipe()
        proc.standardOutput = pipe
        proc.standardError = pipe
        
        pipe.fileHandleForReading.readabilityHandler = { [weak self] handle in
            let data = handle.availableData
            guard !data.isEmpty, let text = String(data: data, encoding: .utf8) else { return }
            
            DispatchQueue.main.async {
                self?.enqueueIncomingOutput(text)
            }
        }
        
        proc.terminationHandler = { [weak self] terminatedProc in
            DispatchQueue.main.async {
                self?.flushPendingOutput()
                self?.timer?.invalidate()
                self?.timer = nil
                self?.csvPollingTimer?.invalidate()
                self?.csvPollingTimer = nil
                self?.pollCSV() // Final flush
                if self?.isSweepRunning == true {
                    self?.harvestSweepCandidateResult()
                } else {
                    self?.pollGBEnergyCSV()
                    if let pt = self?.currentGBEnergy {
                        self?.upsertGBEnergyPointAndWriteFiles(pt)
                    }
                }
                self?.pollDislocationProfiles()
                
                let duration = self?.startTime.map { Date().timeIntervalSince($0) } ?? 0
                let exitCode = terminatedProc.terminationStatus
                
                pipe.fileHandleForReading.readabilityHandler = nil
                self?.process = nil
                
                if exitCode == 0 {
                    self?.status = .completed(exitCode: exitCode, duration: duration)
                    self?.appendLog(text: "🎉 Simulation completed successfully in \(String(format: "%.1f", duration))s")
                    if let strongSelf = self, strongSelf.isSweepRunning {
                        if !strongSelf.pendingSweepCandidates.isEmpty {
                            strongSelf.runNextSweepCandidate()
                            return
                        } else {
                            // Sweep finished! Restore base nx, ny and summarize outputs
                            if let cfg = strongSelf.activeConfig, cfg.sweepAutoCommensurateBox {
                                cfg.nx = strongSelf.sweepSavedBaseNx
                                cfg.ny = strongSelf.sweepSavedBaseNy
                            }
                            strongSelf.outputDirectory = strongSelf.sweepOutputDirectory
                            strongSelf.appendLog(text: "==========================================================")
                            strongSelf.appendLog(text: "✅ CSL Misorientation Sweep Complete (\(strongSelf.gbEnergyPoints.count) angles)!")
                            strongSelf.appendLog(text: "📄 Angle vs Energy Table: \(strongSelf.sweepAngleVsEnergyFilePath)")
                            strongSelf.appendLog(text: "📦 Relaxed VTK Configurations: \(strongSelf.sweepRelaxedVTKDirectory)")
                            strongSelf.appendLog(text: "🎬 ParaView Series File: \(strongSelf.sweepRelaxedVTKDirectory)/csl_sweep_angles.vtk.series")
                            strongSelf.appendLog(text: "==========================================================")
                        }
                    }
                } else if terminatedProc.terminationReason == .uncaughtSignal {
                    self?.status = .terminated
                    self?.appendLog(text: "🛑 Process was terminated by signal", isError: true)
                    if let strongSelf = self, strongSelf.isSweepRunning, let cfg = strongSelf.activeConfig, cfg.sweepAutoCommensurateBox {
                        cfg.nx = strongSelf.sweepSavedBaseNx
                        cfg.ny = strongSelf.sweepSavedBaseNy
                    }
                } else {
                    self?.status = .completed(exitCode: exitCode, duration: duration)
                    self?.appendLog(text: "⚠️ Simulation ended with exit code \(exitCode)", isError: true)
                    if let strongSelf = self, strongSelf.isSweepRunning, let cfg = strongSelf.activeConfig, cfg.sweepAutoCommensurateBox {
                        cfg.nx = strongSelf.sweepSavedBaseNx
                        cfg.ny = strongSelf.sweepSavedBaseNy
                    }
                }
                self?.isSweepRunning = false
                self?.currentSweepCandidate = nil
                self?.sweepProgressText = ""
                self?.currentStageLabel = ""
            }
        }
        
        do {
            try proc.run()
            self.process = proc
            self.status = .running(pid: proc.processIdentifier)
            
            // Start elapsed time timer
            self.timer = Timer.scheduledTimer(withTimeInterval: 1.0, repeats: true) { [weak self] _ in
                guard let self = self, let start = self.startTime else { return }
                self.elapsedTime = Date().timeIntervalSince(start)
            }
            
            // Start CSV polling timer for live Stress-Strain & GB Energy curves
            self.csvPollingTimer = Timer.scheduledTimer(withTimeInterval: 0.35, repeats: true) { [weak self] _ in
                self?.pollCSV()
                if self?.isSweepRunning != true {
                    self?.pollGBEnergyCSV()
                }
            }
        } catch {
            self.status = .failed(error: error.localizedDescription)
            self.appendLog(text: "❌ Failed to start process: \(error.localizedDescription)", isError: true)
        }
    }
    
    public func startCSLSweep(config: SimulationConfig, candidates: [CSLRationalMatch]) {
        guard !status.isRunning, !candidates.isEmpty else { return }
        self.activeConfig = config
        self.sweepSavedBaseNx = config.nx
        self.sweepSavedBaseNy = config.ny
        self.sweepBaseTargetDir = config.targetDirectory
        self.sweepOutputDirectory = "\(sweepBaseTargetDir)/csl_sweep"
        self.sweepRelaxedVTKDirectory = "\(sweepOutputDirectory)/relaxed_configurations"
        self.sweepAngleVsEnergyFilePath = "\(sweepOutputDirectory)/angle_vs_energy.txt"
        
        let fm = FileManager.default
        try? fm.createDirectory(atPath: sweepOutputDirectory, withIntermediateDirectories: true)
        try? fm.createDirectory(atPath: sweepRelaxedVTKDirectory, withIntermediateDirectories: true)
        try? fm.removeItem(atPath: "\(sweepOutputDirectory)/gb_energy_vs_misorientation.csv")
        try? fm.removeItem(atPath: sweepAngleVsEnergyFilePath)
        
        // Reset live chart points so user watches the new sweep plot live point-by-point
        self.gbEnergyPoints.removeAll()
        self.currentGBEnergy = nil
        self.requestedChartMode = .gbEnergy
        self.startTime = Date()
        self.elapsedTime = 0
        
        self.pendingSweepCandidates = candidates
        self.totalSweepCount = candidates.count
        self.currentSweepIndex = 0
        self.isSweepRunning = true
        self.appendLog(text: "==========================================================")
        self.appendLog(text: "🔄 Starting CSL Misorientation Sweep (\(candidates.count) angles)")
        self.appendLog(text: "📂 Sweep Output Folder: \(sweepOutputDirectory)")
        self.appendLog(text: "📦 Relaxed VTK Folder:  \(sweepRelaxedVTKDirectory)")
        self.appendLog(text: "📄 Angle vs Energy Txt: \(sweepAngleVsEnergyFilePath)")
        self.appendLog(text: "==========================================================")
        runNextSweepCandidate()
    }
    
    private func runNextSweepCandidate() {
        guard isSweepRunning, !pendingSweepCandidates.isEmpty, let cfg = activeConfig else {
            isSweepRunning = false
            currentSweepCandidate = nil
            sweepProgressText = ""
            return
        }
        let cand = pendingSweepCandidates.removeFirst()
        currentSweepCandidate = cand
        currentSweepIndex += 1
        let sigSuffix = cand.sigma > 0 ? " (Σ\(cand.sigma))" : ""
        sweepProgressText = "Sweep \(currentSweepIndex)/\(totalSweepCount): Δθ=\(String(format: "%.2f°", cand.fullMisoriDeg))\(sigSuffix)"
        currentStageLabel = sweepProgressText
        
        // Commensurate box size K * L_CSL closest to user's System Geometry NX/NY
        let runNx: Int
        let runNy: Int
        if cfg.sweepAutoCommensurateBox && cand.periodA0 > 1.01 {
            let kPeriodsX = max(2, Int(round(Double(sweepSavedBaseNx) / max(2.0, cand.periodA0))))
            let kPeriodsY = max(2, Int(round(Double(sweepSavedBaseNy) / max(2.0, cand.periodA0))))
            runNx = cand.commensurateSize(periods: kPeriodsX)
            runNy = cand.commensurateSize(periods: kPeriodsY)
        } else {
            runNx = sweepSavedBaseNx
            runNy = sweepSavedBaseNy
        }
        
        // Configure Bicrystal for this CSL candidate on a value-type copy
        // so the live UI SimulationConfig does not jitter or trigger recursive .onChange handlers!
        var sweepGeom = cfg.geometryConfig
        sweepGeom.enabled = true
        sweepGeom.microstructure.mode = .bicrystal
        sweepGeom.microstructure.relax_only = true
        sweepGeom.microstructure.theta1_deg = -cand.halfAngleDeg
        sweepGeom.microstructure.theta2_deg = +cand.halfAngleDeg
        sweepGeom.microstructure.min_seam_distance_a0 = cand.recommendedDMinA0
        sweepGeom.microstructure.gb_normal_gap_a0 = cand.recommendedGapA0
        sweepGeom.microstructure.gb_tangential_shift_a0 = 0.0
        
        let angleSubfolder = String(
            format: "angle_%02d_dtheta_%05.2fdeg_Sigma%d_%d%d0",
            currentSweepIndex, cand.fullMisoriDeg, cand.sigma, cand.m, cand.n
        )
        let angleRunDir = "\(sweepOutputDirectory)/\(angleSubfolder)"
        start(
            config: cfg,
            overrideRunDirectory: angleRunDir,
            overrideNx: runNx,
            overrideNy: runNy,
            overrideGeometryConfig: sweepGeom,
            overrideProtocolMode: .sweepMisorientation
        )
    }
    
    private func harvestSweepCandidateResult() {
        let fm = FileManager.default
        let sumPath = "\(outputDirectory)/gb_energy_summary.csv"
        guard fm.fileExists(atPath: sumPath),
              let content = try? String(contentsOfFile: sumPath, encoding: .utf8) else {
            return
        }
        var parsedPt: GBEnergyPoint? = nil
        for line in content.components(separatedBy: .newlines) {
            if let pt = parseGBEnergyLine(line) {
                parsedPt = pt
                break
            }
        }
        guard var pt = parsedPt else { return }
        
        if let cand = currentSweepCandidate {
            pt.sigma = cand.sigma
            pt.m = cand.m
            pt.n = cand.n
        }
        
        // Copy the relaxed VTK configuration into relaxed_configurations/ and baseTargetDir
        let srcPostVTK = "\(outputDirectory)/vtk_files/configuration_0_post.vtk"
        let srcFallbackVTK = "\(outputDirectory)/vtk_files/configuration_0.vtk"
        let srcVTK = fm.fileExists(atPath: srcPostVTK) ? srcPostVTK : (fm.fileExists(atPath: srcFallbackVTK) ? srcFallbackVTK : "")
        
        if !srcVTK.isEmpty {
            let vtkFileName = String(
                format: "relaxed_angle_%02d_dtheta_%05.2fdeg_Sigma%d_%d%d0.vtk",
                currentSweepIndex, pt.misorientationDeg, pt.sigma, pt.m, pt.n
            )
            let dstSweepVTK = "\(sweepRelaxedVTKDirectory)/\(vtkFileName)"
            try? fm.removeItem(atPath: dstSweepVTK)
            try? fm.copyItem(atPath: srcVTK, toPath: dstSweepVTK)
            pt.relaxedVTKPath = dstSweepVTK
            
            if !sweepBaseTargetDir.isEmpty {
                let dstBaseVTK = "\(sweepBaseTargetDir)/\(vtkFileName)"
                try? fm.removeItem(atPath: dstBaseVTK)
                try? fm.copyItem(atPath: srcVTK, toPath: dstBaseVTK)
            }
        }
        
        upsertGBEnergyPointAndWriteFiles(pt)
        let vtkShort = (pt.relaxedVTKPath as NSString).lastPathComponent
        appendLog(text: String(
            format: "✅ [Sweep %d/%d] Δθ = %.2f° (Σ%d, %d %d 0) | γ_GB(int) = %.6f | γ_GB(tot) = %.6f | VTK: %@",
            currentSweepIndex, totalSweepCount, pt.misorientationDeg, pt.sigma, pt.m, pt.n,
            pt.gbEnergyInterior, pt.gbEnergyTotal, vtkShort
        ))
    }
    
    private func upsertGBEnergyPointAndWriteFiles(_ pt: GBEnergyPoint) {
        self.currentGBEnergy = pt
        var pts = self.gbEnergyPoints
        if let existingIdx = pts.firstIndex(where: { abs($0.misorientationDeg - pt.misorientationDeg) < 0.05 }) {
            pts[existingIdx] = pt
        } else {
            pts.append(pt)
        }
        pts.sort { $0.misorientationDeg < $1.misorientationDeg }
        self.gbEnergyPoints = pts
        writeAngleVsEnergyFiles(points: pts)
    }
    
    private func writeAngleVsEnergyFiles(points: [GBEnergyPoint]) {
        guard !points.isEmpty else { return }
        let fm = FileManager.default
        
        // 1. Human & Gnuplot/Python-readable ASCII table: angle_vs_energy.txt
        var txt = "# =====================================================================================================================\n"
        txt += "# Symmetric Tilt Grain Boundary Energy vs. Misorientation Angle (Relaxed Equilibrium at alpha = 0)\n"
        txt += "# Columns:\n"
        txt += "#  1:DeltaTheta_deg  2:HalfAngle_deg  3:Sigma  4:m  5:n  6:Gamma_GB_interior  7:Gamma_GB_total  8:E_total_relaxed  9:E_interior_relaxed  10:L_GB_a0  11:NX  12:NY  13:Gap_a0  14:d_min_a0  15:Relaxed_VTK_File\n"
        txt += "# =====================================================================================================================\n"
        for p in points {
            let vtkName = p.relaxedVTKPath.isEmpty ? "configuration_0_post.vtk" : (p.relaxedVTKPath as NSString).lastPathComponent
            txt += String(
                format: "%10.4f %10.4f %6d %4d %4d %16.8f %16.8f %16.8f %16.8f %10.2f %5d %5d %9.4f %9.4f  %@\n",
                p.misorientationDeg, p.halfAngleDeg, p.sigma, p.m, p.n,
                p.gbEnergyInterior, p.gbEnergyTotal, p.totalEnergy, p.interiorEnergy,
                p.lengthGBA0, p.nx, p.ny, p.normalGapA0, p.seamDMinA0, vtkName
            )
        }
        
        // 2. Standard CSV: gb_energy_vs_misorientation.csv
        var csv = "MisorientationDeg,HalfAngleDeg,Theta1Deg,Theta2Deg,GBAngleDeg,Nx,Ny,L_GB_a0,L_Int_a0,TotalEnergy,InteriorEnergy,GBEnergyTotal,GBEnergyInterior,NormalGapA0,SeamDMinA0,Sigma,m,n,RelaxedVTKPath\n"
        for p in points {
            csv += String(
                format: "%.6f,%.6f,%.6f,%.6f,%.6f,%d,%d,%.4f,%.4f,%.8f,%.8f,%.8f,%.8f,%.4f,%.4f,%d,%d,%d,%@\n",
                p.misorientationDeg, p.halfAngleDeg, p.theta1Deg, p.theta2Deg, p.gbAngleDeg,
                p.nx, p.ny, p.lengthGBA0, p.lengthIntA0, p.totalEnergy, p.interiorEnergy,
                p.gbEnergyTotal, p.gbEnergyInterior, p.normalGapA0, p.seamDMinA0,
                p.sigma, p.m, p.n, p.relaxedVTKPath
            )
        }
        
        var targetDirs: [String] = []
        if !sweepOutputDirectory.isEmpty { targetDirs.append(sweepOutputDirectory) }
        if !sweepBaseTargetDir.isEmpty { targetDirs.append(sweepBaseTargetDir) }
        if !outputDirectory.isEmpty { targetDirs.append(outputDirectory) }
        let root = SimulationConfig.detectProjectRoot()
        if !root.isEmpty { targetDirs.append(root) }
        
        for dir in Set(targetDirs) {
            try? txt.write(toFile: "\(dir)/angle_vs_energy.txt", atomically: true, encoding: .utf8)
            try? csv.write(toFile: "\(dir)/gb_energy_vs_misorientation.csv", atomically: true, encoding: .utf8)
        }
        if sweepAngleVsEnergyFilePath.isEmpty, let firstDir = targetDirs.first {
            sweepAngleVsEnergyFilePath = "\(firstDir)/angle_vs_energy.txt"
        }
        
        // 3. ParaView .vtk.series JSON in sweepRelaxedVTKDirectory so ParaView can play through all angles Δθ
        if !sweepRelaxedVTKDirectory.isEmpty && fm.fileExists(atPath: sweepRelaxedVTKDirectory) {
            var fileEntries: [String] = []
            for p in points where !p.relaxedVTKPath.isEmpty {
                let baseName = (p.relaxedVTKPath as NSString).lastPathComponent
                fileEntries.append(String(format: "    { \"name\": \"%@\", \"time\": %.4f }", baseName, p.misorientationDeg))
            }
            if !fileEntries.isEmpty {
                let seriesJSON = "{\n  \"file-series-version\": \"1.0\",\n  \"files\": [\n" + fileEntries.joined(separator: ",\n") + "\n  ]\n}\n"
                try? seriesJSON.write(
                    toFile: "\(sweepRelaxedVTKDirectory)/csl_sweep_angles.vtk.series",
                    atomically: true,
                    encoding: .utf8
                )
            }
        }
    }

    public func stop() {
        if isSweepRunning, let cfg = activeConfig, cfg.sweepAutoCommensurateBox {
            cfg.nx = sweepSavedBaseNx
            cfg.ny = sweepSavedBaseNy
        }
        isSweepRunning = false
        currentSweepCandidate = nil
        pendingSweepCandidates.removeAll()
        sweepProgressText = ""
        currentStageLabel = ""
        guard let proc = process, proc.isRunning else { return }
        appendLog(text: "🛑 Sending SIGTERM to process (PID: \(proc.processIdentifier))...", isError: true)
        proc.terminate()
        csvPollingTimer?.invalidate()
        csvPollingTimer = nil
        pollCSV()
    }
    
    private var isPolling = false
    public func pollCSV() {
        guard !isPolling else { return }
        let csvPath = "\(outputDirectory)/energy_stress_log.csv"
        guard FileManager.default.fileExists(atPath: csvPath) else { return }
        
        isPolling = true
        DispatchQueue.global(qos: .utility).async { [weak self] in
            guard let self = self,
                  let content = try? String(contentsOfFile: csvPath, encoding: .utf8) else {
                DispatchQueue.main.async { self?.isPolling = false }
                return
            }
            
            // Force / stress control logs the enthalpy, which falls at every load step: count the
            // solver's own avalanche flags (fractional_drops.csv) instead of energy drops
            var solverAvalanches: Int? = nil
            if content.contains("StrainBurstFraction"),
               let drops = try? String(contentsOfFile: "\(self.outputDirectory)/fractional_drops.csv", encoding: .utf8) {
                var count = 0
                for row in drops.components(separatedBy: .newlines).dropFirst() {
                    let last = row.trimmingCharacters(in: .whitespaces).components(separatedBy: ",").last
                    if last == "1" { count += 1 }
                }
                solverAvalanches = count
            }

            let lines = content.components(separatedBy: .newlines)
            var allPoints: [StressStrainPoint] = []
            allPoints.reserveCapacity(lines.count)

            // Columns 0–22 are fixed; the optional trailing columns differ per run type (traction
            // runs: EffectiveStrain, ElasticEnergy; Volterra cylinder loads: AppliedSigmaXY, Fbar…,
            // Core…, CoreShift, CoreJump), so look those up by header name
            var effectiveStrainCol: Int? = nil
            var elasticEnergyCol: Int? = nil
            var coreShiftCol: Int? = nil
            var coreShiftRelaxedCol: Int? = nil
            var appliedSigmaCol: Int? = nil
            var burgersCol: Int? = nil
            var coreStartCol: Int? = nil
            if let header = lines.first(where: { $0.hasPrefix("Iteration") }) {
                let names = header.components(separatedBy: ",").map { $0.trimmingCharacters(in: .whitespacesAndNewlines) }
                effectiveStrainCol = names.firstIndex(of: "EffectiveStrain")
                elasticEnergyCol = names.firstIndex(of: "ElasticEnergy")
                coreShiftCol = names.firstIndex(of: "CoreShift")
                coreShiftRelaxedCol = names.firstIndex(of: "CoreShiftRelaxed")
                appliedSigmaCol = names.firstIndex(of: "AppliedSigmaXY")
                burgersCol = names.firstIndex(of: "Burgers")
                coreStartCol = names.firstIndex(of: "CoreStart")
            }

            for line in lines {
                let trimmed = line.trimmingCharacters(in: .whitespacesAndNewlines)
                if trimmed.isEmpty || trimmed.hasPrefix("Iteration") { continue }
                
                let parts = trimmed.components(separatedBy: ",")
                guard parts.count >= 11 else { continue }
                
                guard let iter = Int(parts[0].trimmingCharacters(in: .whitespaces)),
                      let alpha = Double(parts[1].trimmingCharacters(in: .whitespaces)),
                      let preEnergy = Double(parts[2].trimmingCharacters(in: .whitespaces)),
                      let preStress = Double(parts[3].trimmingCharacters(in: .whitespaces)),
                      let postEnergy = Double(parts[4].trimmingCharacters(in: .whitespaces)),
                      let postStress = Double(parts[5].trimmingCharacters(in: .whitespaces)) else {
                    continue
                }
                
                let remeshStr = parts[10].trimmingCharacters(in: .whitespacesAndNewlines)
                let shouldRemesh = (remeshStr == "1" || remeshStr.lowercased() == "true")

                func parseCol(_ idx: Int, fallback: Double = 0.0) -> Double {
                    guard idx < parts.count,
                          let val = Double(parts[idx].trimmingCharacters(in: .whitespacesAndNewlines)),
                          val.isFinite else {
                        return fallback
                    }
                    return val
                }
                func parseOptional(_ idx: Int?) -> Double? {
                    guard let idx = idx, idx < parts.count,
                          let val = Double(parts[idx].trimmingCharacters(in: .whitespacesAndNewlines)),
                          val.isFinite else {
                        return nil
                    }
                    return val
                }

                allPoints.append(StressStrainPoint(
                    iteration: iter,
                    alpha: alpha,
                    preStress: preStress,
                    postStress: postStress,
                    preEnergy: preEnergy,
                    postEnergy: postEnergy,
                    shouldRemesh: shouldRemesh,
                    cauchyXX: parseCol(11),
                    cauchyXY: parseCol(12, fallback: postStress),
                    cauchyYY: parseCol(13),
                    hydrostaticPressure: parseCol(14),
                    piola11: parseCol(15),
                    piola12: parseCol(16),
                    piola22: parseCol(17),
                    vonMisesAvg: parseCol(18),
                    vonMisesMax: parseCol(19),
                    indenterLoad: parseCol(20),
                    meanContactPressure: parseCol(21),
                    contactAtoms: parseCol(22),
                    effectiveStrain: parseOptional(effectiveStrainCol) ?? 0.0, // traction runs only
                    elasticEnergy: parseOptional(elasticEnergyCol), // same runs
                    coreShift: parseOptional(coreShiftCol), // Volterra cylinder load runs
                    coreShiftRelaxed: parseOptional(coreShiftRelaxedCol),
                    appliedSigmaXY: parseOptional(appliedSigmaCol),
                    burgers: parseOptional(burgersCol),
                    coreStart: parseOptional(coreStartCol)
                ))
            }
            
            DispatchQueue.main.async {
                self.isPolling = false
                
                if allPoints.count != self.curvePoints.count || allPoints.last?.iteration != self.curvePoints.last?.iteration {
                    self.curvePoints = allPoints
                    if let last = allPoints.last {
                        self.currentStep = last.iteration
                        self.currentAlpha = last.alpha
                        self.currentStress = last.postStress
                        self.currentEnergy = last.postEnergy
                        self.currentConjugateStrain = allPoints.contains { $0.effectiveStrain != 0 }
                            ? last.effectiveStrain : nil
                        if abs(last.postStress) > self.peakStress {
                            self.peakStress = abs(last.postStress)
                            self.yieldStrain = last.alpha
                        }
                    }
                    self.updateAvalancheStats(solverCount: solverAvalanches)
                }
                self.pollEigenCSV()
                self.pollLatticeRotationCSV()
                if !self.isSweepRunning {
                    self.pollGBEnergyCSV()
                }
            }
        }
    }
    
    private func parseGBEnergyLine(_ line: String) -> GBEnergyPoint? {
        let trimmed = line.trimmingCharacters(in: .whitespacesAndNewlines)
        if trimmed.isEmpty || trimmed.hasPrefix("MisorientationDeg") { return nil }
        let parts = trimmed.components(separatedBy: ",")
        guard parts.count >= 13 else { return nil }
        func d(_ idx: Int) -> Double? {
            guard idx < parts.count else { return nil }
            return Double(parts[idx].trimmingCharacters(in: .whitespacesAndNewlines))
        }
        guard let misori = d(0),
              let halfA = d(1),
              let th1 = d(2),
              let th2 = d(3),
              let gbAng = d(4),
              let nx = Int(parts[5].trimmingCharacters(in: .whitespacesAndNewlines)),
              let ny = Int(parts[6].trimmingCharacters(in: .whitespacesAndNewlines)),
              let lGB = d(7),
              let lInt = d(8),
              let totE = d(9),
              let intE = d(10),
              let gbTot = d(11),
              let gbInt = d(12) else {
            return nil
        }
        let gap = d(13) ?? 0.0
        let dmin = d(14) ?? 0.0
        let sig = parts.count > 15 ? Int(parts[15].trimmingCharacters(in: .whitespacesAndNewlines)) : nil
        let mVal = parts.count > 16 ? Int(parts[16].trimmingCharacters(in: .whitespacesAndNewlines)) : nil
        let nVal = parts.count > 17 ? Int(parts[17].trimmingCharacters(in: .whitespacesAndNewlines)) : nil
        let vtkPath = parts.count > 18 ? parts[18].trimmingCharacters(in: .whitespacesAndNewlines) : ""
        return GBEnergyPoint(
            misorientationDeg: misori,
            halfAngleDeg: halfA,
            theta1Deg: th1,
            theta2Deg: th2,
            gbAngleDeg: gbAng,
            nx: nx,
            ny: ny,
            lengthGBA0: lGB,
            lengthIntA0: lInt,
            totalEnergy: totE,
            interiorEnergy: intE,
            gbEnergyTotal: gbTot,
            gbEnergyInterior: gbInt,
            normalGapA0: gap,
            seamDMinA0: dmin,
            sigma: sig,
            m: mVal,
            n: nVal,
            relaxedVTKPath: vtkPath
        )
    }
    
    public func pollGBEnergyCSV() {
        let fm = FileManager.default
        if !outputDirectory.isEmpty {
            let sumPath = "\(outputDirectory)/gb_energy_summary.csv"
            if fm.fileExists(atPath: sumPath),
               let content = try? String(contentsOfFile: sumPath, encoding: .utf8) {
                for line in content.components(separatedBy: .newlines) {
                    if let pt = parseGBEnergyLine(line) {
                        if self.currentGBEnergy != pt {
                            self.currentGBEnergy = pt
                        }
                        break
                    }
                }
            }
        }
        
        var cumPaths: [String] = []
        if !outputDirectory.isEmpty {
            cumPaths.append("\(outputDirectory)/csl_sweep/gb_energy_vs_misorientation.csv")
            cumPaths.append("\(outputDirectory)/gb_energy_vs_misorientation.csv")
        } else {
            let root = SimulationConfig.detectProjectRoot()
            cumPaths.append("\(root)/gb_energy_vs_misorientation.csv")
        }
        
        var foundPoints = false
        for cPath in cumPaths where !cPath.isEmpty {
            if fm.fileExists(atPath: cPath),
               let content = try? String(contentsOfFile: cPath, encoding: .utf8) {
                var pts: [GBEnergyPoint] = []
                for line in content.components(separatedBy: .newlines) {
                    if let pt = parseGBEnergyLine(line) {
                        pts.append(pt)
                    }
                }
                pts.sort { $0.misorientationDeg < $1.misorientationDeg }
                if pts != self.gbEnergyPoints {
                    self.gbEnergyPoints = pts
                }
                foundPoints = true
                break
            }
        }
        if !foundPoints && !isSweepRunning && !self.gbEnergyPoints.isEmpty {
            self.gbEnergyPoints = []
        }
    }
    
    public func clearGBEnergyHistory() {
        let fm = FileManager.default
        if !outputDirectory.isEmpty {
            let parentDir = (outputDirectory as NSString).deletingLastPathComponent
            try? fm.removeItem(atPath: "\(parentDir)/gb_energy_vs_misorientation.csv")
        }
        let root = SimulationConfig.detectProjectRoot()
        try? fm.removeItem(atPath: "\(root)/gb_energy_vs_misorientation.csv")
        self.gbEnergyPoints.removeAll()
    }
    
    /// Lattice rotation log written by the solver every load step (lattice_rotation_log.csv)
    public func pollLatticeRotationCSV() {
        let path = "\(outputDirectory)/lattice_rotation_log.csv"
        guard FileManager.default.fileExists(atPath: path),
              let content = try? String(contentsOfFile: path, encoding: .utf8) else { return }
        let parsed = LatticeRotationLog.parse(content)
        if parsed.count != latticeRotationPoints.count || parsed.last != latticeRotationPoints.last {
            latticeRotationPoints = parsed
        }
    }

    public func pollEigenCSV() {
        let eigPath = "\(outputDirectory)/eigen_log.csv"
        guard FileManager.default.fileExists(atPath: eigPath) else { return }
        guard let content = try? String(contentsOfFile: eigPath, encoding: .utf8) else { return }
        
        let lines = content.components(separatedBy: .newlines)
        
        var parsed: [EigenPoint] = []
        // Optional stress-control columns, found by name (they follow the lambda_k columns)
        var strainCol = -1, boxCol = -1, tangentCol = -1
        for line in lines {
            let trimmed = line.trimmingCharacters(in: .whitespacesAndNewlines)
            if trimmed.hasPrefix("Iteration") {
                let names = trimmed.components(separatedBy: ",").map { $0.trimmingCharacters(in: .whitespaces) }
                strainCol = names.firstIndex(of: "Strain") ?? -1
                boxCol = names.firstIndex(of: "lambda_box") ?? -1
                tangentCol = names.firstIndex(of: "box_tangent") ?? -1
                continue
            }
            if trimmed.isEmpty { continue }
            
            let parts = trimmed.components(separatedBy: ",")
            guard parts.count >= 6 else { continue }
            
            guard let iter = Int(parts[0].trimmingCharacters(in: .whitespaces)),
                  let alpha = Double(parts[1].trimmingCharacters(in: .whitespaces)),
                  let failedShifts = Int(parts[3].trimmingCharacters(in: .whitespaces)),
                  let lambda1 = Double(parts[5].trimmingCharacters(in: .whitespaces)),
                  lambda1.isFinite else {
                continue
            }
            
            let trigger = parts[2].trimmingCharacters(in: .whitespaces)
            var lambda2: Double? = nil
            if parts.count >= 7, let l2 = Double(parts[6].trimmingCharacters(in: .whitespaces)), l2.isFinite {
                lambda2 = l2
            }
            
            func column(_ k: Int) -> Double? {
                guard k >= 0, k < parts.count,
                      let v = Double(parts[k].trimmingCharacters(in: .whitespaces)), v.isFinite else { return nil }
                return v
            }
            
            parsed.append(EigenPoint(
                iteration: iter,
                alpha: alpha,
                trigger: trigger,
                lambdaMin: lambda1,
                lambda2: lambda2,
                failedShifts: failedShifts,
                strain: column(strainCol),
                lambdaBox: column(boxCol),
                boxTangent: column(tangentCol)
            ))
        }
        
        if parsed.count != self.eigenPoints.count || parsed.last?.iteration != self.eigenPoints.last?.iteration {
            self.eigenPoints = parsed
            self.currentLambdaMin = parsed.last?.lambdaMin
        }
        
        self.pollDislocationProfiles()
    }
    
    public func pollDislocationProfiles(dir: String? = nil) {
        let outDir = outputDirectory
        let targetDir = activeConfig?.targetDirectory ?? ""
        DispatchQueue.global(qos: .utility).async { [weak self] in
            guard let self = self else { return }
            let fm = FileManager.default
            var pathsToTry: [String] = []
            if let d = dir, !d.isEmpty {
                pathsToTry.append("\(d)/dislocation_profiles.json")
            }
            if !outDir.isEmpty {
                pathsToTry.append("\(outDir)/dislocation_profiles.json")
            }
            if !targetDir.isEmpty {
                pathsToTry.append("\(targetDir)/dislocation_profiles.json")
            }
            
            let root = SimulationConfig.detectProjectRoot()
            pathsToTry.append("\(root)/dislocation_profiles.json")
            
            // Also scan for recent run_volterra_* directories in project root
            if let contents = try? fm.contentsOfDirectory(atPath: root) {
                let sorted = contents.filter { $0.hasPrefix("run_volterra") }.compactMap { sub -> (path: String, date: Date)? in
                    let fullPath = "\(root)/\(sub)/dislocation_profiles.json"
                    if let attr = try? fm.attributesOfItem(atPath: fullPath),
                       let date = attr[.modificationDate] as? Date {
                        return (fullPath, date)
                    }
                    return nil
                }.sorted { $0.date > $1.date }
                
                for item in sorted {
                    pathsToTry.append(item.path)
                }
            }
            
            pathsToTry.append("\(fm.currentDirectoryPath)/dislocation_profiles.json")
            
            for jsonPath in pathsToTry where !jsonPath.isEmpty {
                guard fm.fileExists(atPath: jsonPath),
                      let attr = try? fm.attributesOfItem(atPath: jsonPath) else { continue }
                let modTime = (attr[.modificationDate] as? Date)?.timeIntervalSince1970 ?? 0
                let size = (attr[.size] as? NSNumber)?.int64Value ?? 0
                let cacheKey = "\(jsonPath)|\(modTime)|\(size)"
                if cacheKey == self.lastDislocationProfileKey {
                    return
                }
                guard let data = try? Data(contentsOf: URL(fileURLWithPath: jsonPath)) else { continue }
                do {
                    let decoded = try JSONDecoder().decode(DislocationProfileJSONRoot.self, from: data)
                    self.lastDislocationProfileKey = cacheKey
                    DispatchQueue.main.async {
                        self.dislocationData = decoded
                    }
                    return
                } catch {
                    print("Error decoding dislocation_profiles.json at \(jsonPath): \(error)")
                }
            }
        }
    }
    
    /// Computes true physical energy drops (E_{i} < E_{i-1}) across the entire relaxed curve
    private func updateAvalancheStats(solverCount: Int? = nil) {
        if let n = solverCount {
            self.avalancheCount = n
            return
        }
        guard curvePoints.count >= 2 else {
            self.avalancheCount = 0
            return
        }
        
        var drops = 0
        for i in 1..<curvePoints.count {
            let prevE = curvePoints[i - 1].postEnergy
            let currE = curvePoints[i].postEnergy
            if (prevE - currE) > 1e-8 {
                drops += 1
            }
        }
        self.avalancheCount = drops
    }
    
    private func enqueueIncomingOutput(_ text: String) {
        pendingOutputChunks.append(text)
        guard !isOutputFlushScheduled else { return }
        isOutputFlushScheduled = true
        DispatchQueue.main.asyncAfter(deadline: .now() + 0.10) { [weak self] in
            self?.flushPendingOutput()
        }
    }
    
    private func flushPendingOutput() {
        isOutputFlushScheduled = false
        guard !pendingOutputChunks.isEmpty else { return }
        let combined = pendingOutputChunks.joined()
        pendingOutputChunks.removeAll(keepingCapacity: true)
        processIncomingOutput(combined)
    }
    
    private func processIncomingOutput(_ text: String) {
        let lines = text.components(separatedBy: .newlines)
        var newLogs: [LogMessage] = []
        var latestStageLabel = self.currentStageLabel
        var latestStep = self.currentStep
        var latestAlpha = self.currentAlpha
        var latestStress = self.currentStress
        var latestLambdaMin = self.currentLambdaMin
        
        for line in lines where !line.isEmpty {
            newLogs.append(LogMessage(text: line, isError: false))
            parseMetricsBatch(
                from: line,
                stageLabel: &latestStageLabel,
                step: &latestStep,
                alpha: &latestAlpha,
                stress: &latestStress,
                lambdaMin: &latestLambdaMin
            )
        }
        
        if !newLogs.isEmpty {
            var updated = self.logLines
            updated.append(contentsOf: newLogs)
            if updated.count > maxLogLines + 200 {
                updated.removeFirst(updated.count - maxLogLines)
            }
            self.logLines = updated
        }
        if latestStageLabel != self.currentStageLabel {
            self.currentStageLabel = latestStageLabel
        }
        if latestStep != self.currentStep {
            self.currentStep = latestStep
        }
        if latestAlpha != self.currentAlpha {
            self.currentAlpha = latestAlpha
        }
        if latestStress != self.currentStress {
            self.currentStress = latestStress
        }
        if latestLambdaMin != self.currentLambdaMin {
            self.currentLambdaMin = latestLambdaMin
        }
    }
    
    private func parseMetricsBatch(
        from line: String,
        stageLabel: inout String,
        step: inout Int,
        alpha: inout Double,
        stress: inout Double,
        lambdaMin: inout Double?
    ) {
        if line.contains("[Stage 1/2]") {
            stageLabel = "Stage 1/2: Equilibrating GB (α = 0)"
        } else if line.contains("[Stage 2/2]") {
            stageLabel = "Stage 2/2: Mechanical Loading (α > 0)"
        } else if line.contains("[Stage 1/1]") || line.contains("Relax-Only mode active") {
            if !isSweepRunning {
                stageLabel = "Stage 1/1: Relax Only (α = 0)"
            }
        }
        // Parse step progress: "Step   42/1250 | alpha = 0.033600 | E = 1.47195e-01 | stress = 1.86084e+00 | mesh: unchanged (39204 tris) | 0.04s"
        if line.contains("Step") || line.contains("step") {
            let parts = line
                .components(separatedBy: CharacterSet(charactersIn: ",| \t=/"))
                .filter { !$0.isEmpty }
            for (idx, p) in parts.enumerated() {
                if p.lowercased() == "step" && idx + 1 < parts.count, let s = Int(parts[idx + 1]) {
                    step = s
                }
                if p.lowercased() == "alpha" && idx + 1 < parts.count, let a = Double(parts[idx + 1]) {
                    alpha = a
                }
                if (p.lowercased() == "sigma" || p.lowercased() == "stress") && idx + 1 < parts.count,
                   let val = Double(parts[idx + 1]) {
                    stress = val
                }
            }
        }
        if line.contains("lambda_min") {
            if let range = line.range(of: "lambda_min = ") {
                let substr = line[range.upperBound...]
                let token = substr.prefix(while: { !$0.isWhitespace && $0 != "," && $0 != "|" })
                if let val = Double(token), val.isFinite {
                    lambdaMin = val
                }
            }
        }
    }
    
    private func appendLog(text: String, isError: Bool = false) {
        logLines.append(LogMessage(text: text, isError: isError))
        if logLines.count > maxLogLines + 200 {
            logLines.removeFirst(200)
        }
    }
}
