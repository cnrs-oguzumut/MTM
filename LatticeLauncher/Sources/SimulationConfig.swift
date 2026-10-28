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
        case .volterra: return "Edge Dislocation (Volterra / Stroh)"
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

public enum PotentialMode: Int, CaseIterable, Identifiable {
    case zanzotto = 3
    case ljV2 = 1
    case ljV3 = 2
    case square = 4
    
    public var id: Int { rawValue }
    public var displayName: String {
        switch self {
        case .zanzotto: return "Mode 3: Zanzotto Continuum Φ(C)"
        case .ljV2: return "Mode 1: Lennard-Jones v2 V(r)"
        case .ljV3: return "Mode 2: Lennard-Jones v3 V(r)"
        case .square: return "Mode 4: Square Gaussian Potential V(r)"
        }
    }
    public var shortName: String {
        switch self {
        case .zanzotto: return "Zanzotto Φ(C)"
        case .ljV2: return "LJ v2"
        case .ljV3: return "LJ v3"
        case .square: return "Square V(r)"
        }
    }
}

public typealias AcousticPotentialMode = PotentialMode

public enum LatticeSymmetry: String, CaseIterable, Identifiable {
    case square = "square"
    case triangular = "triangular"
    
    public var id: String { rawValue }
    public var displayName: String {
        switch self {
        case .square: return "Square Lattice"
        case .triangular: return "Triangular Lattice"
        }
    }
}

public class SimulationConfig: ObservableObject {
    // 1. System Size
    @Published public var nx: Int = 100
    @Published public var ny: Int = 100
    
    // 2. Lattice & Potential Model
    @Published public var potentialMode: PotentialMode = .zanzotto
    @Published public var latticeSymmetry: LatticeSymmetry = .square
    @Published public var beta: Double = -0.25
    @Published public var bulkK: Double = 5.0
    
    // 3. Loading & Mode
    @Published public var mode: LoadingMode = .positive
    @Published public var seed: Int = 42
    @Published public var noiseLevel: Double = 0.04
    @Published public var enableRemeshing: Bool = true
    @Published public var rFree: Double = 35.0
    @Published public var exportFullMesh: Bool = true
    
    // Dislocation Stress Relaxation (Full Cauchy Tensor)
    @Published public var relaxMacroStress: Bool = false
    @Published public var stressTol: Double = 1e-5
    @Published public var maxMacroIters: Int = 30
    
    // 3. Loading Schedule (Shear modes)
    @Published public var boundaryCondition: BoundaryConditionMode = .periodic
    @Published public var alphaStart: Double = 0.14
    @Published public var alphaEnd: Double = 1.0
    @Published public var stepSize: Double = 0.00006
    
    // 3b. Staircase Shift Schedule (Shift mode)
    @Published public var shiftAmplitude: Double = 2.0
    @Published public var shiftHorizontalSteps: Int = 2
    @Published public var shiftVerticalSteps: Int = 2
    @Published public var shiftPerturbMesh: Bool = false
    
    // Rotated Lattice & CSL inputs
    @Published public var useRotation: Bool = false
    @Published public var rotationAngle: Double = 0.0
    @Published public var cslM: Int = 2
    @Published public var cslN: Int = 1
    @Published public var cslKx: Int = 45
    @Published public var cslKy: Int = 45
    @Published public var autoScaleSupercell: Bool = true
    
    public func updateCSLSupercellAuto() {
        guard autoScaleSupercell else { return }
        let allowed = CSLAngleOption.allowedAngles(for: latticeSymmetry)
        let opt = allowed.first(where: { $0.m == cslM && $0.n == cslN }) ?? allowed[0]
        cslKx = max(1, Int(round(Double(nx) / opt.Lx_ratio)))
        cslKy = max(1, Int(round(Double(ny) / opt.Ly_ratio)))
    }
    
    public var cslEstimatedAtoms: Int {
        let allowed = CSLAngleOption.allowedAngles(for: latticeSymmetry)
        let opt = allowed.first(where: { $0.m == cslM && $0.n == cslN }) ?? allowed[0]
        return cslKx * cslKy * opt.atoms_cell
    }
    
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
    @Published public var traceInnerIterations: Bool = false
    @Published public var saveSurgeryVTK: Bool = false
    @Published public var saveTriangleData: Bool = false
    @Published public var saveVTKInterval: Int = 50
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
    
    // 9. Macroscopic Geometry & Boundary Condition Studio
    @Published public var geometryConfig: MacroscopicGeometryConfigModel = MacroscopicGeometryConfigModel()
    
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
        
        // 5. Fallback to user Documents (auto-create for standalone .app / .dmg users)
        let docsFallback = "\(FileManager.default.homeDirectoryForCurrentUser.path)/Documents/LatticeRuns"
        try? FileManager.default.createDirectory(atPath: "\(docsFallback)/scratch", withIntermediateDirectories: true)
        return docsFallback
    }
    
    public static func detectBinaryPath(projectRoot: String? = nil) -> String {
        let root = projectRoot ?? detectProjectRoot()
        let path1 = "\(root)/build/lattice_triangulation"
        let path2 = "\(root)/lattice_triangulation"
        let bundledPath1 = Bundle.main.bundleURL.appendingPathComponent("Contents/MacOS/lattice_triangulation").path
        let bundledPath2 = Bundle.main.executableURL?.deletingLastPathComponent().appendingPathComponent("lattice_triangulation").path ?? ""
        
        if FileManager.default.fileExists(atPath: path1) { return path1 }
        if FileManager.default.fileExists(atPath: path2) { return path2 }
        if FileManager.default.fileExists(atPath: bundledPath1) { return bundledPath1 }
        if !bundledPath2.isEmpty && FileManager.default.fileExists(atPath: bundledPath2) { return bundledPath2 }
        return bundledPath1
    }
    
    public var binaryPath: String {
        return SimulationConfig.detectBinaryPath(projectRoot: projectRoot)
    }
    
    public init() {
        self.projectRoot = SimulationConfig.detectProjectRoot()
        updateCSLSupercellAuto()
    }
    
    public var defaultRunDirName: String {
        let remeshTag = enableRemeshing ? "remesh" : "noremesh"
        let latTag = latticeSymmetry == .triangular ? "_tri" : ""
        let bcTag = (mode != .shift && mode != .volterra && boundaryCondition == .fixed) ? "_fixed" : ""
        let potTag = potentialMode == .zanzotto ? "" : "_pot\(potentialMode.rawValue)"
        let rotTag = useRotation ? ((boundaryCondition == .fixed || geometryConfig.enabled || mode == .volterra || mode == .shift) ? String(format: "_rot%.1f", rotationAngle) : "_csl\(cslM)_\(cslN)_K\(cslKx)x\(cslKy)") : ""
        let relaxTag = (mode == .volterra && relaxMacroStress) ? "_relax" : ""
        let geomTag = (mode != .volterra && mode != .shift && geometryConfig.enabled) ? "_\(geometryConfig.preset_name)" : ""
        if mode == .shift {
            let pertTag = shiftPerturbMesh ? "_pert" : ""
            return String(format: "run_shift_%dx%d_amp%.1f_%@%@%@%@%@", nx, ny, shiftAmplitude, remeshTag, pertTag, latTag, potTag, rotTag)
        }
        return "run_\(mode.rawValue)_\(nx)x\(ny)\(bcTag)\(geomTag)_seed\(seed)_\(remeshTag)\(latTag)\(potTag)\(rotTag)\(relaxTag)"
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
        updateCSLSupercellAuto()
    }
    
    /// Constructs the array of command-line arguments for the executable
    public func buildArguments() -> [String] {
        var args: [String] = []
        
        // Positional arguments: nx ny mode seed remesh
        args.append(String(nx))
        let effectiveNy = (latticeSymmetry == .triangular && ny % 2 != 0) ? (ny + 1) : ny
        args.append(String(effectiveNy))
        args.append(mode.rawValue)
        args.append(String(seed))
        args.append(enableRemeshing ? "1" : "0")
        
        // Lattice Symmetry & Potential Parameters
        args.append("--lattice=\(latticeSymmetry.rawValue)")
        args.append("--potential-mode=\(potentialMode.rawValue)")
        if potentialMode == .zanzotto {
            args.append(String(format: "--beta=%.4f", beta))
            args.append(String(format: "--K=%.4f", bulkK))
        }
        args.append(String(format: "--noise-level=%.6f", noiseLevel))
        
        
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
        
        // Loading Schedule vs Dislocation Cylinder vs Staircase Shift
        if mode == .volterra {
            args.append(String(format: "--r-free=%.2f", rFree))
            args.append(exportFullMesh ? "--full-mesh" : "--circle")
            if relaxMacroStress {
                args.append("--relax-stress=1")
                args.append(String(format: "--stress-tol=%.1e", stressTol))
                args.append("--max-macro-iters=\(maxMacroIters)")
            }
        } else if mode == .shift {
            args.append(String(format: "--shift-amplitude=%.4f", shiftAmplitude))
            args.append("--shift-h-steps=\(shiftHorizontalSteps)")
            args.append("--shift-v-steps=\(shiftVerticalSteps)")
            args.append("--shift-perturb=\(shiftPerturbMesh ? 1 : 0)")
        } else {
            args.append(String(format: "--alpha-start=%.6f", alphaStart))
            args.append(String(format: "--alpha-end=%.6f", alphaEnd))
            args.append(String(format: "--step-size=%.8f", stepSize))
        }
        
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
        
        if saveVTKInterval > 0 {
            args.append("--save-vtk-every=\(saveVTKInterval)")
        } else {
            args.append("--save-vtk-every=0")
        }
        
        if maxAvalanches > 0 {
            args.append("--max-avalanches=\(maxAvalanches)")
        }
        
        // Stability Monitor
        if enableEigenvalues && mode != .shift {
            let everyVal = (mode == .volterra) ? max(1, eigEvery) : eigEvery
            args.append("--eig-every=\(everyVal)")
            args.append("--eig-at-avalanche=\(eigAtAvalanche ? 1 : 0)")
            args.append("--eig-modes=\(eigModes)")
            args.append("--eig-vectors=\(eigVectors)")
            args.append("--eig-solver=\(eigSolver.rawValue)")
        } else {
            args.append("--eig-every=0")
        }
        
        // Boundary condition & Macroscopic Geometry Studio
        if mode != .shift && mode != .volterra {
            if geometryConfig.enabled {
                args.append("--bc=fixed")
                let scratchDir = "\(projectRoot)/scratch"
                try? FileManager.default.createDirectory(atPath: scratchDir, withIntermediateDirectories: true)
                let scratchGeomPath = "\(scratchDir)/geometry_config.json"
                _ = try? geometryConfig.writeJSON(to: scratchGeomPath)
                
                let targetDir = targetDirectory
                try? FileManager.default.createDirectory(atPath: targetDir, withIntermediateDirectories: true)
                let runGeomPath = "\(targetDir)/geometry_config.json"
                if (try? geometryConfig.writeJSON(to: runGeomPath)) != nil {
                    args.append("--geometry-json=\(runGeomPath)")
                } else {
                    args.append("--geometry-json=\(scratchGeomPath)")
                }
            } else {
                args.append("--bc=\(boundaryCondition.rawValue)")
            }
        }
        
        // Rotated Crystal Orientation / CSL
        if mode != .volterra && useRotation {
            args.append("--use-rotation")
            if boundaryCondition == .fixed || geometryConfig.enabled || mode == .shift {
                args.append(String(format: "--rotation-angle=%.4f", rotationAngle))
            } else {
                args.append("--csl-m=\(cslM)")
                args.append("--csl-n=\(cslN)")
                args.append("--csl-kx=\(cslKx)")
                args.append("--csl-ky=\(cslKy)")
                if rotationAngle > 0.001 {
                    args.append(String(format: "--rotation-angle=%.4f", rotationAngle))
                }
            }
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
