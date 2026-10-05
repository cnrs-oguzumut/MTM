import Foundation
import Combine

public enum PreconditionerType: String, CaseIterable, Identifiable {
    case stiffness = "stiffness"
    case none = "none"
    case laplacian = "laplacian"
    case diag = "diag"
    
    public var id: String { rawValue }
    public var displayName: String {
        switch self {
        case .stiffness: return "Stiffness (CHOLMOD / Supernodal)"
        case .none: return "None (Standard L-BFGS)"
        case .laplacian: return "Laplacian Metric"
        case .diag: return "Diagonal Metric"
        }
    }
}

public enum LoadingMode: String, CaseIterable, Identifiable {
    case positive = "positive"
    case negative = "negative"
    case shift = "shift"
    case volterra = "volterra"
    
    public var id: String { rawValue }
    public var displayName: String {
        switch self {
        case .positive: return "Positive Shear (+α)"
        case .negative: return "Negative Shear (-α)"
        case .shift: return "Shift / Staircase"
        case .volterra: return "Volterra Dislocation"
        }
    }
}

public enum EigenSolverType: String, CaseIterable, Identifiable {
    case fast = "fast"
    case legacy = "legacy"
    
    public var id: String { rawValue }
    public var displayName: String {
        switch self {
        case .fast: return "Fast (Shift-Invert + Cholesky)"
        case .legacy: return "Legacy (Spectra / SparseLU)"
        }
    }
}

public class SimulationConfig: ObservableObject {
    // 1. System Size
    @Published public var nx: Int = 100
    @Published public var ny: Int = 100
    
    // 2. Loading & Mode
    @Published public var mode: LoadingMode = .positive
    @Published public var seed: Int = 42
    @Published public var enableRemeshing: Bool = true
    
    // 3. Loading Schedule
    @Published public var alphaStart: Double = 0.14
    @Published public var alphaEnd: Double = 1.0
    @Published public var stepSize: Double = 0.00006
    
    // 4. Preconditioning
    @Published public var preconditioner: PreconditionerType = .stiffness
    @Published public var precondFromStep: Int = 1
    @Published public var precondRefresh: Int = 0
    @Published public var precondTol: Double = 1e-6
    
    // 5. Stability & Eigenvalues
    @Published public var enableEigenvalues: Bool = false
    @Published public var eigEvery: Int = 5
    @Published public var eigAtAvalanche: Bool = true
    @Published public var eigModes: Int = 5
    @Published public var eigVectors: Int = 2
    @Published public var eigSolver: EigenSolverType = .fast
    
    // 6. Data Saving & Output Controls
    @Published public var traceInnerIterations: Bool = true
    @Published public var saveSurgeryVTK: Bool = true
    @Published public var saveTriangleData: Bool = false
    @Published public var checkpointInterval: Int = 500
    @Published public var stressDropThreshold: Double = 0.10
    @Published public var maxAvalanches: Int = 0
    
    // 7. Hardware & Threads
    @Published public var ompThreads: Int = max(2, ProcessInfo.processInfo.activeProcessorCount - 2)
    @Published public var blasThreads: Int = max(2, ProcessInfo.processInfo.activeProcessorCount - 2)
    
    // 8. Paths & Directory
    @Published public var projectRoot: String = SimulationConfig.detectProjectRoot()
    @Published public var customOutputDirectory: String = ""
    @Published public var useCustomDirectory: Bool = false
    
    public static func detectProjectRoot() -> String {
        let hardcodedDefault = "/Users/usalman/programming/FEM_2D/factorized/lattice_triangulation"
        let bundleURL = Bundle.main.bundleURL
        
        // 1. If running inside LatticeLauncher/Lattice Launcher.app
        let grandParent = bundleURL.deletingLastPathComponent().deletingLastPathComponent().path
        if FileManager.default.fileExists(atPath: "\(grandParent)/src/main.cpp") {
            return grandParent
        }
        
        // 2. If running inside project root
        let parent = bundleURL.deletingLastPathComponent().path
        if FileManager.default.fileExists(atPath: "\(parent)/src/main.cpp") {
            return parent
        }
        
        // 3. Current directory (if not root /)
        let cwd = FileManager.default.currentDirectoryPath
        if cwd != "/" && FileManager.default.fileExists(atPath: "\(cwd)/src/main.cpp") {
            return cwd
        }
        
        // 4. Default project directory
        if FileManager.default.fileExists(atPath: hardcodedDefault) {
            return hardcodedDefault
        }
        
        // 5. Fallback to user Documents
        return "\(FileManager.default.homeDirectoryForCurrentUser.path)/Documents/LatticeRuns"
    }
    
    public var binaryPath: String {
        let path1 = "\(projectRoot)/build/lattice_triangulation"
        let path2 = "\(projectRoot)/lattice_triangulation"
        if FileManager.default.fileExists(atPath: path1) { return path1 }
        if FileManager.default.fileExists(atPath: path2) { return path2 }
        return path1
    }
    
    public init() {
        self.projectRoot = SimulationConfig.detectProjectRoot()
    }
    
    public var defaultRunDirName: String {
        let remeshTag = enableRemeshing ? "remesh" : "noremesh"
        return "run_\(mode.rawValue)_\(nx)x\(ny)_seed\(seed)_\(remeshTag)"
    }
    
    public var targetDirectory: String {
        if useCustomDirectory && !customOutputDirectory.isEmpty {
            return customOutputDirectory
        }
        return "\(projectRoot)/\(defaultRunDirName)"
    }
    
    public func randomizeSeed() {
        seed = Int.random(in: 1...999999)
    }
    
    public func setSize(nx: Int, ny: Int) {
        self.nx = nx
        self.ny = ny
    }
    
    /// Constructs the array of command-line arguments for the executable
    public func buildArguments() -> [String] {
        var args: [String] = []
        
        // Positional arguments: nx ny mode seed remesh
        args.append(String(nx))
        args.append(String(ny))
        args.append(mode.rawValue)
        args.append(String(seed))
        args.append(enableRemeshing ? "1" : "0")
        
        // Preconditioner
        args.append("--precond=\(preconditioner.rawValue)")
        if preconditioner != .none {
            if precondFromStep > 0 {
                args.append("--precond-from-step=\(precondFromStep)")
            }
            if precondRefresh > 0 {
                args.append("--precond-refresh=\(precondRefresh)")
            }
            if precondTol != 1e-6 {
                args.append("--precond-tol=\(precondTol)")
            }
        }
        
        // Loading Schedule
        args.append(String(format: "--alpha-start=%.6f", alphaStart))
        args.append(String(format: "--alpha-end=%.6f", alphaEnd))
        args.append(String(format: "--step-size=%.8f", stepSize))
        
        // Data & Checkpoints
        if checkpointInterval > 0 {
            args.append("--checkpoint-interval=\(checkpointInterval)")
        } else {
            args.append("--checkpoint-interval=0")
        }
        args.append(String(format: "--stress-drop-threshold=%.4f", stressDropThreshold))
        
        // Tracing and VTK
        if traceInnerIterations {
            args.append("--trace=1")
        } else {
            args.append("--trace=0")
        }
        
        if saveSurgeryVTK {
            args.append("--surgery-vtk=1")
        } else {
            args.append("--surgery-vtk=0")
        }
        
        if saveTriangleData {
            args.append("--save-triangle-data=1")
        } else {
            args.append("--no-triangle-data")
        }
        
        if maxAvalanches > 0 {
            args.append("--max-avalanches=\(maxAvalanches)")
        }
        
        // Stability Monitor
        if enableEigenvalues {
            args.append("--eig-every=\(eigEvery)")
            args.append("--eig-at-avalanche=\(eigAtAvalanche ? 1 : 0)")
            args.append("--eig-modes=\(eigModes)")
            args.append("--eig-vectors=\(eigVectors)")
            args.append("--eig-solver=\(eigSolver.rawValue)")
        } else {
            args.append("--eig-every=0")
        }
        
        return args
    }
    
    /// Returns the exact shell command string for copying or display
    public func fullCommandLine() -> String {
        let envs = "OMP_NUM_THREADS=\(ompThreads) MKL_NUM_THREADS=\(blasThreads) OPENBLAS_NUM_THREADS=\(blasThreads)"
        let argsStr = buildArguments().joined(separator: " ")
        return "\(envs) \(binaryPath) \(argsStr)"
    }
}
