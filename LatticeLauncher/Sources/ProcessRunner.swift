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
    @Published public var elapsedTime: TimeInterval = 0
    
    // Live Stress-Strain & Eigenvalue Curve points
    @Published public var curvePoints: [StressStrainPoint] = []
    @Published public var eigenPoints: [EigenPoint] = []
    @Published public var currentLambdaMin: Double? = nil
    @Published public var peakStress: Double = 0.0
    @Published public var yieldStrain: Double = 0.0
    @Published public var activeConfig: SimulationConfig?
    
    private var process: Process?
    private var timer: Timer?
    private var csvPollingTimer: Timer?
    private var startTime: Date?
    private let maxLogLines = 10_000
    private var outputDirectory: String = ""
    private var lastReadCSVLineCount: Int = 0
    
    public init() {}
    
    public func clearLogs() {
        logLines.removeAll()
        curvePoints.removeAll()
        eigenPoints.removeAll()
        currentLambdaMin = nil
        peakStress = 0.0
        yieldStrain = 0.0
    }
    
    public func loadExistingIfAvailable(dir: String) {
        let csvPath = "\(dir)/energy_stress_log.csv"
        guard FileManager.default.fileExists(atPath: csvPath) else { return }
        self.outputDirectory = dir
        self.pollCSV()
    }
    
    public func openOutputFolder() {
        let dir = outputDirectory.isEmpty ? FileManager.default.currentDirectoryPath : outputDirectory
        NSWorkspace.shared.selectFile(nil, inFileViewerRootedAtPath: dir)
    }
    
    public func start(config: SimulationConfig) {
        guard !status.isRunning else { return }
        
        // 1. Verify binary exists
        let fm = FileManager.default
        let binary = config.binaryPath
        guard fm.fileExists(atPath: binary) else {
            self.status = .failed(error: "Binary not found at: \(binary)\nPlease compile first with cmake or make.")
            self.appendLog(text: "❌ Error: Executable does not exist at '\(binary)'", isError: true)
            return
        }
        
        // 2. Prepare target directory
        let runDir = config.targetDirectory
        self.outputDirectory = runDir
        do {
            try fm.createDirectory(atPath: runDir, withIntermediateDirectories: true, attributes: nil)
            // Remove stale CSV logs before fresh start so live chart starts clean
            try? fm.removeItem(atPath: "\(runDir)/energy_stress_log.csv")
            try? fm.removeItem(atPath: "\(runDir)/eigen_log.csv")
        } catch {
            self.status = .failed(error: "Cannot create output directory: \(error.localizedDescription)")
            return
        }
        
        // 3. Clear stats and set running state
        self.activeConfig = config
        self.curvePoints.removeAll()
        self.eigenPoints.removeAll()
        self.currentLambdaMin = nil
        self.peakStress = 0.0
        self.yieldStrain = 0.0
        self.currentStep = 0
        self.currentAlpha = config.alphaStart
        self.currentStress = 0.0
        self.currentEnergy = 0.0
        self.avalancheCount = 0
        self.elapsedTime = 0
        self.startTime = Date()
        
        self.appendLog(text: "==========================================================")
        self.appendLog(text: "🚀 Launching Simulation: \(config.nx)x\(config.ny), mode=\(config.mode.rawValue), seed=\(config.seed)")
        self.appendLog(text: "📁 Working Directory: \(runDir)")
        self.appendLog(text: "⚙️ Command: \(config.fullCommandLine())")
        self.appendLog(text: "==========================================================")
        
        // 4. Configure Process
        let proc = Process()
        proc.executableURL = URL(fileURLWithPath: binary)
        proc.arguments = config.buildArguments()
        proc.currentDirectoryURL = URL(fileURLWithPath: runDir)
        
        var env = ProcessInfo.processInfo.environment
        env["OMP_NUM_THREADS"] = String(config.ompThreads)
        env["MKL_NUM_THREADS"] = String(config.blasThreads)
        env["OPENBLAS_NUM_THREADS"] = String(config.blasThreads)
        env["VECLIB_MAXIMUM_THREADS"] = String(config.blasThreads)
        
        // Ensure Homebrew and Apple libraries are found
        let extraPaths = "/opt/homebrew/lib:/opt/homebrew/opt/libomp/lib:/usr/local/lib"
        if let existingDyld = env["DYLD_LIBRARY_PATH"] {
            env["DYLD_LIBRARY_PATH"] = "\(extraPaths):\(existingDyld)"
        } else {
            env["DYLD_LIBRARY_PATH"] = extraPaths
        }
        proc.environment = env
        
        // 5. Pipes
        let pipe = Pipe()
        proc.standardOutput = pipe
        proc.standardError = pipe
        
        pipe.fileHandleForReading.readabilityHandler = { [weak self] handle in
            let data = handle.availableData
            guard !data.isEmpty, let text = String(data: data, encoding: .utf8) else { return }
            
            DispatchQueue.main.async {
                self?.processIncomingOutput(text)
            }
        }
        
        proc.terminationHandler = { [weak self] terminatedProc in
            DispatchQueue.main.async {
                self?.timer?.invalidate()
                self?.timer = nil
                self?.csvPollingTimer?.invalidate()
                self?.csvPollingTimer = nil
                self?.pollCSV() // Final flush
                
                let duration = self?.startTime.map { Date().timeIntervalSince($0) } ?? 0
                let exitCode = terminatedProc.terminationStatus
                
                if exitCode == 0 {
                    self?.status = .completed(exitCode: exitCode, duration: duration)
                    self?.appendLog(text: "🎉 Simulation completed successfully in \(String(format: "%.1f", duration))s")
                } else if terminatedProc.terminationReason == .uncaughtSignal {
                    self?.status = .terminated
                    self?.appendLog(text: "🛑 Process was terminated by signal", isError: true)
                } else {
                    self?.status = .completed(exitCode: exitCode, duration: duration)
                    self?.appendLog(text: "⚠️ Simulation ended with exit code \(exitCode)", isError: true)
                }
                pipe.fileHandleForReading.readabilityHandler = nil
                self?.process = nil
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
            
            // Start CSV polling timer for live Stress-Strain curve
            self.csvPollingTimer = Timer.scheduledTimer(withTimeInterval: 0.4, repeats: true) { [weak self] _ in
                self?.pollCSV()
            }
        } catch {
            self.status = .failed(error: error.localizedDescription)
            self.appendLog(text: "❌ Failed to start process: \(error.localizedDescription)", isError: true)
        }
    }
    
    public func stop() {
        guard let proc = process, proc.isRunning else { return }
        appendLog(text: "🛑 Sending SIGTERM to process (PID: \(proc.processIdentifier))...", isError: true)
        proc.terminate()
        csvPollingTimer?.invalidate()
        csvPollingTimer = nil
        pollCSV()
    }
    
    public func pollCSV() {
        self.pollEigenCSV()
        
        let csvPath = "\(outputDirectory)/energy_stress_log.csv"
        guard FileManager.default.fileExists(atPath: csvPath) else { return }
        guard let content = try? String(contentsOfFile: csvPath, encoding: .utf8) else { return }
        
        var lines = content.components(separatedBy: "\n")
        if !content.hasSuffix("\n") && !lines.isEmpty {
            // Drop incomplete last line currently being written
            lines.removeLast()
        }
        
        var newPoints: [StressStrainPoint] = []
        let currentLastIter = self.curvePoints.last?.iteration ?? -1
        
        for line in lines {
            let trimmed = line.trimmingCharacters(in: .whitespacesAndNewlines)
            if trimmed.isEmpty || trimmed.hasPrefix("Iteration") { continue }
            
            let parts = trimmed.components(separatedBy: ",")
            guard parts.count >= 11 else { continue }
            
            guard let iter = Int(parts[0].trimmingCharacters(in: .whitespaces)),
                  iter > currentLastIter,
                  let alpha = Double(parts[1].trimmingCharacters(in: .whitespaces)),
                  let preEnergy = Double(parts[2].trimmingCharacters(in: .whitespaces)),
                  let preStress = Double(parts[3].trimmingCharacters(in: .whitespaces)),
                  let postEnergy = Double(parts[4].trimmingCharacters(in: .whitespaces)),
                  let postStress = Double(parts[5].trimmingCharacters(in: .whitespaces)) else {
                continue
            }
            
            let remeshStr = parts[10].trimmingCharacters(in: .whitespacesAndNewlines)
            let shouldRemesh = (remeshStr == "1" || remeshStr.lowercased() == "true")
            
            let pt = StressStrainPoint(
                iteration: iter,
                alpha: alpha,
                preStress: preStress,
                postStress: postStress,
                preEnergy: preEnergy,
                postEnergy: postEnergy,
                shouldRemesh: shouldRemesh
            )
            newPoints.append(pt)
        }
        
        if !newPoints.isEmpty {
            self.curvePoints.append(contentsOf: newPoints)
            if let last = self.curvePoints.last {
                self.currentStep = last.iteration
                self.currentAlpha = last.alpha
                self.currentStress = last.postStress
                self.currentEnergy = last.postEnergy
                if abs(last.postStress) > self.peakStress {
                    self.peakStress = abs(last.postStress)
                    self.yieldStrain = last.alpha
                }
            }
            self.updateAvalancheStats()
        }
        
        self.pollEigenCSV()
    }
    
    public func pollEigenCSV() {
        let eigPath = "\(outputDirectory)/eigen_log.csv"
        guard FileManager.default.fileExists(atPath: eigPath) else { return }
        guard let content = try? String(contentsOfFile: eigPath, encoding: .utf8) else { return }
        
        var lines = content.components(separatedBy: "\n")
        if !content.hasSuffix("\n") && !lines.isEmpty {
            lines.removeLast()
        }
        
        var parsed: [EigenPoint] = []
        for line in lines {
            let trimmed = line.trimmingCharacters(in: .whitespacesAndNewlines)
            if trimmed.isEmpty || trimmed.hasPrefix("Iteration") { continue }
            
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
            
            parsed.append(EigenPoint(
                iteration: iter,
                alpha: alpha,
                trigger: trigger,
                lambdaMin: lambda1,
                lambda2: lambda2,
                failedShifts: failedShifts
            ))
        }
        
        if parsed.count != self.eigenPoints.count {
            self.eigenPoints = parsed
            self.currentLambdaMin = parsed.last?.lambdaMin
        }
    }
    
    /// Computes true physical energy drops (E_{i} < E_{i-1}) across the entire relaxed curve
    private func updateAvalancheStats() {
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
    
    private func processIncomingOutput(_ text: String) {
        let lines = text.components(separatedBy: .newlines)
        for line in lines where !line.isEmpty {
            appendLog(text: line)
            parseMetrics(from: line)
        }
    }
    
    private func parseMetrics(from line: String) {
        // Parse step progress: "Step 00120 | alpha = 0.1472 | E = -1234.56 | sigma = 0.281"
        if line.contains("Step") || line.contains("step") {
            let parts = line.components(separatedBy: CharacterSet(charactersIn: ",| \t"))
            for (idx, p) in parts.enumerated() {
                if p.lowercased() == "step" && idx + 1 < parts.count, let s = Int(parts[idx + 1]) {
                    self.currentStep = s
                }
                if p.lowercased().contains("alpha") && idx + 1 < parts.count, let a = Double(parts[idx + 1].replacingOccurrences(of: "=", with: "")) {
                    self.currentAlpha = a
                }
                if (p.lowercased().contains("sigma") || p.lowercased().contains("stress")) && idx + 1 < parts.count,
                   let val = Double(parts[idx + 1].replacingOccurrences(of: "=", with: "")) {
                    self.currentStress = val
                }
            }
        }
        if line.contains("lambda_min") {
            if let range = line.range(of: "lambda_min = ") {
                let substr = line[range.upperBound...]
                let token = substr.prefix(while: { !$0.isWhitespace && $0 != "," && $0 != "|" })
                if let val = Double(token), val.isFinite {
                    self.currentLambdaMin = val
                }
            }
        }
    }
    
    private func appendLog(text: String, isError: Bool = false) {
        logLines.append(LogMessage(text: text, isError: isError))
        if logLines.count > maxLogLines {
            logLines.removeFirst(logLines.count - maxLogLines)
        }
    }
}
