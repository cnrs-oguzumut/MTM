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
    
    // 3. Loading Schedule & Two-Stage Protocol (Shear & Microstructure modes)
    @Published public var protocolMode: SimulationProtocolMode = .relaxThenLoad
    @Published public var sweepContinuousMode: Bool = true
    @Published public var sweepDeltaThetaStepDeg: Double = 2.0
    @Published public var sweepIncludeExactCusps: Bool = true
    @Published public var sweepThetaMinDeg: Double = 0.0
    @Published public var sweepThetaMaxDeg: Double = 45.0
    @Published public var sweepUseCuratedCSLOnly: Bool = true
    @Published public var sweepAutoCommensurateBox: Bool = true
    @Published public var boundaryCondition: BoundaryConditionMode = .periodic
    @Published public var alphaStart: Double = 0.0
    @Published public var alphaEnd: Double = 0.25
    @Published public var stepSize: Double = 0.0005

    public var sweepCSLCandidates: [CSLRationalMatch] {
        if sweepContinuousMode {
            return CSLRationalizer.buildContinuousSweep(
                thetaMinDeg: sweepThetaMinDeg,
                thetaMaxDeg: sweepThetaMaxDeg,
                stepDeltaThetaDeg: sweepDeltaThetaStepDeg,
                isTriangular: latticeSymmetry == .triangular,
                includeExactCusps: sweepIncludeExactCusps
            )
        }
        let lo = min(sweepThetaMinDeg, sweepThetaMaxDeg)
        let hi = max(sweepThetaMinDeg, sweepThetaMaxDeg)
        let source = sweepUseCuratedCSLOnly ? CSLRationalizer.curatedCatalog : CSLRationalizer.allCandidates
        return source.filter { $0.halfAngleDeg >= lo - 1e-4 && $0.halfAngleDeg <= hi + 1e-4 }
    }
    
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
        let isPeriodicDislocs = (geometryConfig.enabled && geometryConfig.microstructure.mode == .periodicDislocations)
        let bcTag = (mode != .shift && mode != .volterra && boundaryCondition == .fixed && !isPeriodicDislocs) ? "_fixed" : ""
        let potTag = potentialMode == .zanzotto ? "" : "_pot\(potentialMode.rawValue)"
        let rotTag = useRotation ? ((boundaryCondition == .fixed || geometryConfig.enabled || mode == .volterra || mode == .shift) ? String(format: "_rot%.1f", rotationAngle) : "_csl\(cslM)_\(cslN)_K\(cslKx)x\(cslKy)") : ""
        let relaxTag = (mode == .volterra && relaxMacroStress) ? "_relax" : ""
        let geomTag = (mode != .volterra && mode != .shift && geometryConfig.enabled) ? (isPeriodicDislocs ? "_pbc_dislocs_\(2 * geometryConfig.microstructure.dislocation_pair_count)D" : "_\(geometryConfig.preset_name)") : ""
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
        if geometryConfig.microstructure.mode == .periodicDislocations &&
            geometryConfig.microstructure.dislocation_layout != DislocationLayoutPreset.custom.rawValue {
            geometryConfig.microstructure.regeneratePeriodicDislocations(nx: nx, ny: ny)
        }
    }
    
    /// Constructs the array of command-line arguments for the executable
    public func buildArguments(
        overrideRunDirectory: String? = nil,
        overrideNx: Int? = nil,
        overrideNy: Int? = nil,
        overrideGeometryConfig: MacroscopicGeometryConfigModel? = nil,
        overrideProtocolMode: SimulationProtocolMode? = nil
    ) -> [String] {
        var args: [String] = []
        let activeProtocol = overrideProtocolMode ?? protocolMode
        let activeGeometry = overrideGeometryConfig ?? geometryConfig
        
        // Positional arguments: nx ny mode seed remesh
        let useNx = overrideNx ?? nx
        let baseNy = overrideNy ?? ny
        args.append(String(useNx))
        let effectiveNy = (latticeSymmetry == .triangular && baseNy % 2 != 0) ? (baseNy + 1) : baseNy
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
        } else if activeProtocol == .relaxOnly || activeProtocol == .sweepMisorientation {
            args.append("--alpha-start=0.000000")
            args.append("--alpha-end=0.000000")
            args.append(String(format: "--step-size=%.8f", max(1e-6, stepSize)))
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
            let isPeriodicDislocs = (activeGeometry.enabled && activeGeometry.microstructure.mode == .periodicDislocations)
            if activeGeometry.enabled && (isPeriodicDislocs || boundaryCondition != .periodic) {
                args.append(isPeriodicDislocs ? "--bc=periodic" : "--bc=fixed")
                var effectiveGeom = activeGeometry
                effectiveGeom.microstructure.relax_only = (activeProtocol != .relaxThenLoad)
                let scratchDir = "\(projectRoot)/scratch"
                try? FileManager.default.createDirectory(atPath: scratchDir, withIntermediateDirectories: true)
                let scratchGeomPath = "\(scratchDir)/geometry_config.json"
                _ = try? effectiveGeom.writeJSON(to: scratchGeomPath)
                
                let targetDir = overrideRunDirectory ?? targetDirectory
                try? FileManager.default.createDirectory(atPath: targetDir, withIntermediateDirectories: true)
                let runGeomPath = "\(targetDir)/geometry_config.json"
                if (try? effectiveGeom.writeJSON(to: runGeomPath)) != nil {
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
            if boundaryCondition == .fixed || activeGeometry.enabled || mode == .shift {
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

// MARK: - Unified Boundary Regime & Mechanical Loading Truth Model

public enum ActiveBoundaryRegime: String, CaseIterable, Identifiable {
    case periodicPBC = "periodic_pbc"
    case finiteSurfaceBCs = "finite_surface_bcs"
    
    public var id: String { rawValue }
    
    public var title: String {
        switch self {
        case .periodicPBC: return "Full 2D Periodic (PBC)"
        case .finiteSurfaceBCs: return "Surface Grips / Custom BCs"
        }
    }
    
    public var shortBadge: String {
        switch self {
        case .periodicPBC: return "2D Periodic (PBC)"
        case .finiteSurfaceBCs: return "Surface BCs"
        }
    }
    
    public var icon: String {
        switch self {
        case .periodicPBC: return "repeat.circle.fill"
        case .finiteSurfaceBCs: return "rectangle.inset.filled"
        }
    }
}

public enum ActiveMechanicalLoading: String, CaseIterable, Identifiable {
    case simpleShearPos = "simple_shear_pos"
    case simpleShearNeg = "simple_shear_neg"
    case uniaxialTensionY = "uniaxial_tension_y"
    case uniaxialCompY = "uniaxial_comp_y"
    case uniaxialTensionX = "uniaxial_tension_x"
    case shearParallelGB = "shear_parallel_gb"
    case affineShearBox = "affine_shear_box"
    case radialExpansionOutward = "radial_expansion_outward"
    case radialCompressionInward = "radial_compression_inward"
    case nanoindentation = "nanoindentation"
    case staircaseShift = "staircase_shift"
    case volterraDislocation = "volterra"
    case customRegional = "custom_regional"
    
    public var id: String { rawValue }
    
    /// Primary mechanical loading modes shown in the visual grid.
    public static var primaryGridModes: [ActiveMechanicalLoading] {
        [
            .simpleShearPos,
            .simpleShearNeg,
            .uniaxialTensionY,
            .uniaxialCompY,
            .uniaxialTensionX,
            .shearParallelGB,
            .radialExpansionOutward,
            .radialCompressionInward,
            .affineShearBox,
            .nanoindentation
        ]
    }
    
    public var shortTitle: String {
        switch self {
        case .simpleShearPos: return "Simple Shear (+γ_xy)"
        case .simpleShearNeg: return "Reverse Shear (-γ_xy)"
        case .uniaxialTensionY: return "Uniaxial Tension Y (+ε_yy)"
        case .uniaxialCompY: return "Compression Y (-ε_yy)"
        case .uniaxialTensionX: return "Uniaxial Tension X (+ε_xx)"
        case .shearParallelGB: return "Shear ∥ Vertical GB (+Y)"
        case .affineShearBox: return "Affine Shear Box (4 Walls)"
        case .radialExpansionOutward: return "Radial Outward Pull (+ε_r)"
        case .radialCompressionInward: return "Radial Inward Push (-ε_r)"
        case .nanoindentation: return "Nanoindentation (Top Tip)"
        case .staircaseShift: return "Staircase Shift"
        case .volterraDislocation: return "Volterra Dislocation Core"
        case .customRegional: return "Custom Regional BCs"
        }
    }
    
    public var displayName: String { shortTitle }
    
    public var detailSubtitle: String {
        switch self {
        case .simpleShearPos: return "Top driven +X (or PBC F_xy = +α) • Measures σ_xy"
        case .simpleShearNeg: return "Top driven -X (or PBC F_xy = -α) • Measures σ_xy"
        case .uniaxialTensionY: return "Top pulled +Y, Bottom fixed, Sides free • Measures σ_yy"
        case .uniaxialCompY: return "Top pushed -Y, Bottom fixed, Sides free • Measures σ_yy"
        case .uniaxialTensionX: return "Right pulled +X, Left roller, Top/Bot roller • Measures σ_xx"
        case .shearParallelGB: return "Right driven +Y, Left fixed, Top/Bot free • Measures σ_xy"
        case .affineShearBox: return "All outer walls driven by F_xy = +α • Measures σ_xy"
        case .radialExpansionOutward: return "Outer rim pulled radially outward u = +α(r - r_c) • Measures ½(σ_xx+σ_yy)"
        case .radialCompressionInward: return "Outer rim pushed radially inward u = -α(r - r_c) • Measures ½(σ_xx+σ_yy)"
        case .nanoindentation: return "Rigid tip pressing into top free surface • Measures Load P(δ)"
        case .staircaseShift: return "Alternating horizontal & vertical block shifts"
        case .volterraDislocation: return "Single edge dislocation inside cylindrical R_free"
        case .customRegional: return "User-defined regional boundary rules in Step 2"
        }
    }
    
    public var icon: String {
        switch self {
        case .simpleShearPos: return "arrow.right.square.fill"
        case .simpleShearNeg: return "arrow.left.square.fill"
        case .uniaxialTensionY: return "arrow.up.and.down.square.fill"
        case .uniaxialCompY: return "arrow.down.to.line.compact"
        case .uniaxialTensionX: return "arrow.left.and.right.square.fill"
        case .shearParallelGB: return "arrow.up.right.and.arrow.down.left.rectangle.fill"
        case .affineShearBox: return "skew"
        case .radialExpansionOutward: return "arrow.up.left.and.down.right.and.arrow.up.right.and.down.left"
        case .radialCompressionInward: return "arrow.down.right.and.arrow.up.left"
        case .nanoindentation: return "Scope"
        case .staircaseShift: return "stairs"
        case .volterraDislocation: return " hurricane"
        case .customRegional: return "slider.horizontal.3"
        }
    }
    
    public var sfSymbol: String {
        switch self {
        case .simpleShearPos: return "arrow.right.circle.fill"
        case .simpleShearNeg: return "arrow.left.circle.fill"
        case .uniaxialTensionY: return "arrow.up.circle.fill"
        case .uniaxialCompY: return "arrow.down.circle.fill"
        case .uniaxialTensionX: return "arrow.left.and.right.circle.fill"
        case .shearParallelGB: return "arrow.up.and.down.circle.fill"
        case .affineShearBox: return "skew"
        case .radialExpansionOutward: return "arrow.up.left.and.down.right.and.arrow.up.right.and.down.left"
        case .radialCompressionInward: return "arrow.down.right.and.arrow.up.left"
        case .nanoindentation: return "arrow.down.to.line.circle.fill"
        case .staircaseShift: return "stairs"
        case .volterraDislocation: return "scope"
        case .customRegional: return "slider.horizontal.3"
        }
    }
    
    /// Whether this loading mode is compatible with Full 2D Periodic (PBC) torus boundary conditions.
    public var supportsPeriodicPBC: Bool {
        switch self {
        case .simpleShearPos, .simpleShearNeg:
            return true
        default:
            return false
        }
    }
}

public extension SimulationConfig {
    /// True iff the C++ solver will execute with Full 2D Periodic Boundary Conditions (T² torus).
    var isEffectivePeriodicPBC: Bool {
        guard mode != .shift && mode != .volterra else { return false }
        if geometryConfig.enabled {
            return geometryConfig.microstructure.mode == .periodicDislocations
        }
        return boundaryCondition == .periodic
    }
    
    /// The active boundary condition regime that the C++ solver will actually run.
    var activeBoundaryRegime: ActiveBoundaryRegime {
        isEffectivePeriodicPBC ? .periodicPBC : .finiteSurfaceBCs
    }
    
    /// The exact mechanical loading mode that the C++ solver will actually execute.
    var activeMechanicalLoading: ActiveMechanicalLoading {
        if mode == .volterra { return .volterraDislocation }
        if mode == .shift { return .staircaseShift }
        
        // When in 2D PBC or when custom geometry rules are disabled, C++ runs homogeneous F_xy = ±α simple shear.
        if isEffectivePeriodicPBC || !geometryConfig.enabled {
            return (mode == .negative) ? .simpleShearNeg : .simpleShearPos
        }
        
        // Custom geometry rules are active (finite specimen)
        if geometryConfig.indenter.enabled {
            return .nanoindentation
        }
        
        let rules = geometryConfig.boundary_rules.filter(\.enabled)
        func ruleOf(_ reg: BoundaryRegionSelector) -> BoundaryZoneRuleModel? {
            rules.first(where: { $0.region == reg })
        }
        let sign: Double = (mode == .negative) ? -1.0 : 1.0
        
        if let perim = ruleOf(.outerPerimeter), perim.kind == .affine {
            if perim.affine_mode == "uniaxial_y" {
                return sign >= 0 ? .uniaxialTensionY : .uniaxialCompY
            } else if perim.affine_mode == "uniaxial_x" {
                return .uniaxialTensionX
            } else if perim.affine_mode == "biaxial" || perim.affine_mode == "radial_outward" {
                return sign >= 0 ? .radialExpansionOutward : .radialCompressionInward
            } else if perim.affine_mode == "radial_inward" || perim.affine_mode == "biaxial_compression" {
                return sign >= 0 ? .radialCompressionInward : .radialExpansionOutward
            }
            return .affineShearBox
        }
        
        if let bot = ruleOf(.bottom), bot.kind == .fixed,
           let top = ruleOf(.top), top.kind == .prescribed {
            let ux = sign * top.ux
            let uy = sign * top.uy
            if abs(ux - 1.0) < 0.05 && abs(uy) < 0.05 {
                return .simpleShearPos
            }
            if abs(ux + 1.0) < 0.05 && abs(uy) < 0.05 {
                return .simpleShearNeg
            }
            if abs(ux) < 0.05 && abs(uy - 1.0) < 0.05 {
                return .uniaxialTensionY
            }
            if abs(ux) < 0.05 && abs(uy + 1.0) < 0.05 {
                return .uniaxialCompY
            }
        }
        
        if let left = ruleOf(.left), (left.kind == .fixed || left.kind == .rollerY),
           let right = ruleOf(.right), right.kind == .prescribed {
            let ux = sign * right.ux
            let uy = sign * right.uy
            if left.kind == .fixed && abs(ux) < 0.05 && abs(uy - 1.0) < 0.05 {
                return .shearParallelGB
            }
            if abs(ux - 1.0) < 0.05 && abs(uy) < 0.05 {
                return .uniaxialTensionX
            }
        }
        
        return .customRegional
    }
    
    /// Switches the global Boundary Regime while keeping the loading direction consistent.
    func setBoundaryRegime(_ regime: ActiveBoundaryRegime) {
        if mode == .shift || mode == .volterra {
            mode = .positive
        }
        switch regime {
        case .periodicPBC:
            let wasNegShear = (activeMechanicalLoading == .simpleShearNeg)
            boundaryCondition = .periodic
            geometryConfig.indenter.enabled = false
            if geometryConfig.microstructure.mode == .periodicDislocations {
                geometryConfig.enabled = true
                geometryConfig.microstructure.enforce_periodic_bc = true
                geometryConfig.outer_shape = .rectangle
                geometryConfig.cutouts = []
            } else {
                geometryConfig.enabled = false
                if geometryConfig.microstructure.mode != .singleCrystal {
                    geometryConfig.microstructure.mode = .singleCrystal
                }
            }
            mode = wasNegShear ? .negative : .positive
            
        case .finiteSurfaceBCs:
            let prevLoading = activeMechanicalLoading
            boundaryCondition = .fixed
            if geometryConfig.microstructure.mode == .periodicDislocations {
                geometryConfig.microstructure.mode = .singleCrystal
            }
            geometryConfig.enabled = true
            // Preserve shear if the user was previously in Simple Shear
            if prevLoading == .simpleShearPos {
                applyMechanicalLoading(.simpleShearPos)
            } else if prevLoading == .simpleShearNeg {
                applyMechanicalLoading(.simpleShearNeg)
            }
        }
    }
    
    /// Atomically applies the chosen mechanical loading mode and synchronizes all boundary rules and CLI flags.
    func applyMechanicalLoading(_ load: ActiveMechanicalLoading) {
        switch load {
        case .simpleShearPos:
            mode = .positive
            geometryConfig.indenter.enabled = false
            if isEffectivePeriodicPBC {
                // Already in 2D PBC: mode = .positive drives F_xy = +α
                return
            }
            if geometryConfig.enabled {
                geometryConfig.boundary_rules = [
                    BoundaryZoneRuleModel(region: .bottom, kind: .fixed, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                    BoundaryZoneRuleModel(region: .top, kind: .prescribed, thickness_a0: 2.0, ux: 1.0, uy: 0.0),
                    BoundaryZoneRuleModel(region: .left, kind: .free, thickness_a0: 1.5),
                    BoundaryZoneRuleModel(region: .right, kind: .free, thickness_a0: 1.5),
                    BoundaryZoneRuleModel(region: .holeSurfaces, kind: .free, thickness_a0: 1.5)
                ]
                geometryConfig.primary_stress_component = "sigma_xy"
            }
            
        case .simpleShearNeg:
            geometryConfig.indenter.enabled = false
            if isEffectivePeriodicPBC || !geometryConfig.enabled {
                mode = .negative
                return
            }
            mode = .positive
            geometryConfig.boundary_rules = [
                BoundaryZoneRuleModel(region: .bottom, kind: .fixed, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                BoundaryZoneRuleModel(region: .top, kind: .prescribed, thickness_a0: 2.0, ux: -1.0, uy: 0.0),
                BoundaryZoneRuleModel(region: .left, kind: .free, thickness_a0: 1.5),
                BoundaryZoneRuleModel(region: .right, kind: .free, thickness_a0: 1.5),
                BoundaryZoneRuleModel(region: .holeSurfaces, kind: .free, thickness_a0: 1.5)
            ]
            geometryConfig.primary_stress_component = "sigma_xy"
            
        case .uniaxialTensionY:
            mode = .positive
            boundaryCondition = .fixed
            if geometryConfig.microstructure.mode == .periodicDislocations {
                geometryConfig.microstructure.mode = .singleCrystal
            }
            geometryConfig.enabled = true
            geometryConfig.indenter.enabled = false
            geometryConfig.boundary_rules = [
                BoundaryZoneRuleModel(region: .bottom, kind: .fixed, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                BoundaryZoneRuleModel(region: .top, kind: .prescribed, thickness_a0: 2.0, ux: 0.0, uy: 1.0),
                BoundaryZoneRuleModel(region: .left, kind: .free, thickness_a0: 1.5),
                BoundaryZoneRuleModel(region: .right, kind: .free, thickness_a0: 1.5),
                BoundaryZoneRuleModel(region: .holeSurfaces, kind: .free, thickness_a0: 1.5)
            ]
            geometryConfig.primary_stress_component = "sigma_yy"
            
        case .uniaxialCompY:
            mode = .positive
            boundaryCondition = .fixed
            if geometryConfig.microstructure.mode == .periodicDislocations {
                geometryConfig.microstructure.mode = .singleCrystal
            }
            geometryConfig.enabled = true
            geometryConfig.indenter.enabled = false
            geometryConfig.boundary_rules = [
                BoundaryZoneRuleModel(region: .bottom, kind: .fixed, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                BoundaryZoneRuleModel(region: .top, kind: .prescribed, thickness_a0: 2.0, ux: 0.0, uy: -1.0),
                BoundaryZoneRuleModel(region: .left, kind: .free, thickness_a0: 1.5),
                BoundaryZoneRuleModel(region: .right, kind: .free, thickness_a0: 1.5),
                BoundaryZoneRuleModel(region: .holeSurfaces, kind: .free, thickness_a0: 1.5)
            ]
            geometryConfig.primary_stress_component = "sigma_yy"
            
        case .uniaxialTensionX:
            mode = .positive
            boundaryCondition = .fixed
            if geometryConfig.microstructure.mode == .periodicDislocations {
                geometryConfig.microstructure.mode = .singleCrystal
            }
            geometryConfig.enabled = true
            geometryConfig.indenter.enabled = false
            geometryConfig.boundary_rules = [
                BoundaryZoneRuleModel(region: .left, kind: .rollerY, thickness_a0: 2.0),
                BoundaryZoneRuleModel(region: .right, kind: .prescribed, thickness_a0: 2.0, ux: 1.0, uy: 0.0, lock_x: true, lock_y: false),
                BoundaryZoneRuleModel(region: .bottom, kind: .rollerX, thickness_a0: 2.0),
                BoundaryZoneRuleModel(region: .top, kind: .rollerX, thickness_a0: 2.0),
                BoundaryZoneRuleModel(region: .holeSurfaces, kind: .free, thickness_a0: 1.5)
            ]
            geometryConfig.primary_stress_component = "sigma_xx"
            
        case .shearParallelGB:
            mode = .positive
            boundaryCondition = .fixed
            if geometryConfig.microstructure.mode == .periodicDislocations {
                geometryConfig.microstructure.mode = .bicrystal
            }
            geometryConfig.enabled = true
            geometryConfig.indenter.enabled = false
            geometryConfig.boundary_rules = [
                BoundaryZoneRuleModel(region: .left, kind: .fixed, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                BoundaryZoneRuleModel(region: .right, kind: .prescribed, thickness_a0: 2.0, ux: 0.0, uy: 1.0, lock_x: true, lock_y: true),
                BoundaryZoneRuleModel(region: .bottom, kind: .free, thickness_a0: 1.5),
                BoundaryZoneRuleModel(region: .top, kind: .free, thickness_a0: 1.5),
                BoundaryZoneRuleModel(region: .holeSurfaces, kind: .free, thickness_a0: 1.5)
            ]
            geometryConfig.primary_stress_component = "sigma_xy"
            
        case .affineShearBox:
            mode = .positive
            boundaryCondition = .fixed
            if geometryConfig.microstructure.mode == .periodicDislocations {
                geometryConfig.microstructure.mode = .singleCrystal
            }
            geometryConfig.enabled = true
            geometryConfig.indenter.enabled = false
            geometryConfig.boundary_rules = [
                BoundaryZoneRuleModel(region: .outerPerimeter, kind: .affine, thickness_a0: 2.0, affine_mode: "simple_shear_x"),
                BoundaryZoneRuleModel(region: .holeSurfaces, kind: .free, thickness_a0: 1.5)
            ]
            geometryConfig.primary_stress_component = "sigma_xy"
            
        case .radialExpansionOutward:
            mode = .positive
            boundaryCondition = .fixed
            if geometryConfig.microstructure.mode == .periodicDislocations {
                geometryConfig.microstructure.mode = .singleCrystal
            }
            geometryConfig.enabled = true
            geometryConfig.indenter.enabled = false
            geometryConfig.boundary_rules = [
                BoundaryZoneRuleModel(region: .outerPerimeter, kind: .affine, thickness_a0: 2.0, affine_mode: "radial_outward"),
                BoundaryZoneRuleModel(region: .holeSurfaces, kind: .free, thickness_a0: 1.5)
            ]
            geometryConfig.primary_stress_component = "hydrostatic"
            
        case .radialCompressionInward:
            mode = .positive
            boundaryCondition = .fixed
            if geometryConfig.microstructure.mode == .periodicDislocations {
                geometryConfig.microstructure.mode = .singleCrystal
            }
            geometryConfig.enabled = true
            geometryConfig.indenter.enabled = false
            geometryConfig.boundary_rules = [
                BoundaryZoneRuleModel(region: .outerPerimeter, kind: .affine, thickness_a0: 2.0, affine_mode: "radial_inward"),
                BoundaryZoneRuleModel(region: .holeSurfaces, kind: .free, thickness_a0: 1.5)
            ]
            geometryConfig.primary_stress_component = "hydrostatic"
            
        case .nanoindentation:
            mode = .positive
            boundaryCondition = .fixed
            if geometryConfig.microstructure.mode == .periodicDislocations {
                geometryConfig.microstructure.mode = .singleCrystal
            }
            geometryConfig.enabled = true
            geometryConfig.indenter.enabled = true
            geometryConfig.boundary_rules = [
                BoundaryZoneRuleModel(region: .bottom, kind: .fixed, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                BoundaryZoneRuleModel(region: .left, kind: .rollerY, thickness_a0: 2.0),
                BoundaryZoneRuleModel(region: .right, kind: .rollerY, thickness_a0: 2.0),
                BoundaryZoneRuleModel(region: .top, kind: .free, thickness_a0: 1.5),
                BoundaryZoneRuleModel(region: .holeSurfaces, kind: .free, thickness_a0: 1.5)
            ]
            geometryConfig.primary_stress_component = "indentation_load"
            if protocolMode == .relaxOnly {
                protocolMode = .relaxThenLoad
                geometryConfig.microstructure.relax_only = false
            }
            
        case .staircaseShift:
            mode = .shift
            
        case .volterraDislocation:
            mode = .volterra
            
        case .customRegional:
            mode = .positive
            boundaryCondition = .fixed
            geometryConfig.enabled = true
        }
        // Stage 1 follows these rules automatically while it is still linked to them
        // (see MacroscopicGeometryConfigModel.boundary_rules didSet); custom Stage 1 rules are kept.
    }
    
    /// Affine mode of the enabled outer-perimeter rule, if any.
    private var outerPerimeterAffineMode: String? {
        geometryConfig.boundary_rules.first { $0.enabled && $0.region == .outerPerimeter && $0.kind == .affine }?.affine_mode
    }

    var activeStressComponentLabel: String {
        if activeMechanicalLoading == .affineShearBox, let m = outerPerimeterAffineMode {
            if m == "circumferential_shear" || m == "torsion" { return "Rim Torque T = Σ r × f" }
            if m == "pure_shear" { return "Pure-Shear Stress ½(σ_xx−σ_yy)" }
        }
        switch activeMechanicalLoading {
        case .simpleShearPos, .simpleShearNeg, .shearParallelGB, .affineShearBox:
            return "Shear Stress σ_xy"
        case .uniaxialTensionY:
            return "Tensile Stress σ_yy"
        case .uniaxialCompY:
            return "Compressive Stress σ_yy"
        case .uniaxialTensionX:
            return "Axial Stress σ_xx"
        case .radialExpansionOutward, .radialCompressionInward:
            return "Mean Radial Stress ½(σ_xx+σ_yy)"
        case .nanoindentation:
            return "Indentation Load P(δ)"
        case .staircaseShift:
            return "Shear Stress σ_xy"
        case .volterraDislocation:
            return "Core Energy E(r)"
        case .customRegional:
            return "Stress (\(geometryConfig.primary_stress_component))"
        }
    }
    
    var activeStrainParameterNoun: String {
        if activeMechanicalLoading == .affineShearBox, let m = outerPerimeterAffineMode {
            if m == "circumferential_shear" || m == "torsion" { return "rim rotation θ = α (rad)" }
            if m == "pure_shear" { return "pure-shear strain ε" }
        }
        switch activeMechanicalLoading {
        case .simpleShearPos, .simpleShearNeg, .shearParallelGB, .affineShearBox:
            return "shear strain γ_xy"
        case .uniaxialTensionY:
            return "tensile strain ε_yy"
        case .uniaxialCompY:
            return "compressive strain |ε_yy|"
        case .uniaxialTensionX:
            return "axial strain ε_xx"
        case .radialExpansionOutward:
            return "radial outward strain +ε_r"
        case .radialCompressionInward:
            return "radial inward strain |ε_r|"
        case .nanoindentation:
            return "depth ratio α = δ/δ_max"
        default:
            return "load parameter α"
        }
    }
    
    /// Concise, unmistakable one-line summary of the active BC + Loading physics.
    var activePhysicsSummaryBadge: String {
        if mode == .volterra {
            return "Fixed Cylinder (R=\(Int(rFree))a₀) • Single Edge Dislocation"
        }
        if mode == .shift {
            return "Staircase BCs • Block Shift (A=\(String(format: "%.1f", shiftAmplitude))a₀)"
        }
        if protocolMode == .sweepMisorientation {
            return "\(activeBoundaryRegime.shortBadge) • Δθ Misorientation Sweep (\(sweepCSLCandidates.count) runs)"
        }
        if protocolMode == .relaxOnly {
            return "\(activeBoundaryRegime.shortBadge) • Equilibrate Only (α = 0, No External Load)"
        }
        return "\(activeBoundaryRegime.shortBadge) • \(activeMechanicalLoading.shortTitle)"
    }
    
    /// Explicit title for the primary Launch button so the user always sees what will run.
    var launchButtonSummaryTitle: String {
        if mode == .volterra {
            return "Launch Volterra Dislocation"
        }
        if mode == .shift {
            return "Launch Staircase Shift"
        }
        switch protocolMode {
        case .relaxOnly:
            return "Launch Relax Only (α = 0)"
        case .sweepMisorientation:
            return "Launch Δθ Sweep (\(sweepCSLCandidates.count) runs)"
        case .relaxThenLoad:
            let bcShort = isEffectivePeriodicPBC ? "PBC" : "Surface"
            return "Launch: \(bcShort) \(activeMechanicalLoading.shortTitle)"
        }
    }
}
