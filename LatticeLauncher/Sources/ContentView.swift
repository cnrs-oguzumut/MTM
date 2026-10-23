import SwiftUI
import AppKit

enum AppNavigationTab: String, CaseIterable, Identifiable {
    case simulation = "Lattice Simulation"
    case acousticStudy = "Parametric Acoustic Study"
    case stiffnessStudy = "Stiffness Spectrum K(α)"
    
    var id: String { rawValue }
    var icon: String {
        switch self {
        case .simulation: return "square.grid.3x3.fill"
        case .acousticStudy: return "waveform.path.ecg"
        case .stiffnessStudy: return "point.3.filled.connected.trianglepath.dotted"
        }
    }
}

enum MainDisplayTab: String, CaseIterable, Identifiable {
    case curve = "Curve"
    case initialMesh = "Initial Mesh"
    case terminal = "Terminal Log"
    
    var id: String { rawValue }
    
    func title(mode: LoadingMode) -> String {
        switch self {
        case .curve:
            return mode == .volterra ? "Dislocation Profiles" : "Stress–Strain Curve"
        case .initialMesh:
            return "Initial Mesh"
        case .terminal:
            return "Terminal Log"
        }
    }
    
    var icon: String {
        switch self {
        case .curve: return "chart.xyaxis.line"
        case .initialMesh: return "square.stack.3d.down.right.fill"
        case .terminal: return "terminal"
        }
    }
}

struct ContentView: View {
    @ObservedObject var config: SimulationConfig
    @EnvironmentObject var runner: ProcessRunner
    @StateObject private var previewManager = StiffnessStudyManager()
    @State private var appTab: AppNavigationTab = .simulation
    @State private var displayTab: MainDisplayTab = .curve
    @State private var filterText: String = ""
    @State private var autoScroll: Bool = true
    
    var body: some View {
        VStack(spacing: 0) {
            // Global App Top Bar
            HStack(spacing: 12) {
                ScrollView(.horizontal, showsIndicators: false) {
                    Picker("Navigation Mode", selection: $appTab) {
                        ForEach(AppNavigationTab.allCases) { tab in
                            Label(tab.rawValue, systemImage: tab.icon).tag(tab)
                        }
                    }
                    .labelsHidden()
                    .pickerStyle(.segmented)
                    .fixedSize()
                }
                .frame(maxWidth: .infinity, alignment: .leading)
                
                // Active Compute Device Badge
                HStack(spacing: 6) {
                    Image(systemName: config.computeDevice == .gpu ? "bolt.circle.fill" : "cpu")
                        .foregroundColor(config.computeDevice == .gpu ? .purple : .blue)
                        .imageScale(.small)
                    Text(config.computeDevice == .gpu ? "Metal GPU" : "CPU (OpenMP)")
                        .font(.caption.bold())
                        .foregroundColor(.secondary)
                }
                .padding(.horizontal, 10)
                .padding(.vertical, 4)
                .background(Color(NSColor.controlBackgroundColor), in: Capsule())
                .fixedSize()

                // Active Lattice Badge
                HStack(spacing: 6) {
                    Circle()
                        .fill(config.latticeSymmetry == .square ? Color.blue : Color.orange)
                        .frame(width: 8, height: 8)
                    Text(config.latticeSymmetry == .square ? "Square" : "Triangular")
                        .font(.caption.bold())
                        .foregroundColor(.secondary)
                }
                .padding(.horizontal, 10)
                .padding(.vertical, 4)
                .background(Color(NSColor.controlBackgroundColor), in: Capsule())
                .fixedSize()
            }
            .padding(.horizontal, 16)
            .padding(.vertical, 8)
            .background(Color(NSColor.windowBackgroundColor))
            
            Divider()
            
            // Tab Contents
            if appTab == .simulation {
                HSplitView {
                    // LEFT: Configuration Sidebar (Width ~420px)
                    SidebarControlsView(config: config, previewManager: previewManager, displayTab: $displayTab)
                        .frame(minWidth: 380, idealWidth: 420, maxWidth: 480)
                    
                    // RIGHT: Monitor & Live Terminal / Curve
                    VStack(spacing: 0) {
                        // Unified Top Header: Status, Live Metrics & View Switcher
                        LiveStatusHeaderView(config: config, displayTab: $displayTab)
                            .padding(.horizontal, 14)
                            .padding(.vertical, 10)
                            .background(Color(NSColor.controlBackgroundColor))
                        
                        Divider()
                        
                        // Main Content View: Live Curve, Initial Mesh, or Terminal
                        if displayTab == .curve {
                            if config.mode == .volterra {
                                DislocationProfileChartView()
                            } else {
                                StressStrainChartView()
                            }
                        } else if displayTab == .initialMesh {
                            InitialLatticeVisualizerView(manager: previewManager)
                        } else {
                            TerminalView(
                                logLines: runner.logLines,
                                filterText: $filterText,
                                autoScroll: $autoScroll
                            )
                        }
                        
                        Divider()
                        
                        // Bottom Action Bar
                        BottomActionBar(config: config)
                            .padding(.horizontal, 16)
                            .padding(.vertical, 12)
                            .background(Color(NSColor.windowBackgroundColor))
                    }
                    .frame(minWidth: 580)
                }
                .onAppear {
                    runner.loadExistingIfAvailable(dir: config.targetDirectory)
                }
                .onChange(of: config.targetDirectory) { newDir in
                    if !runner.status.isRunning {
                        runner.loadExistingIfAvailable(dir: newDir)
                    }
                }
            } else if appTab == .acousticStudy {
                AcousticStudyView()
            } else if appTab == .stiffnessStudy {
                StiffnessStudyView()
            }
        }
    }
}

// MARK: - Sidebar Controls
struct SidebarControlsView: View {
    @ObservedObject var config: SimulationConfig
    @ObservedObject var previewManager: StiffnessStudyManager
    @Binding var displayTab: MainDisplayTab
    
    var body: some View {
        ScrollView {
            VStack(alignment: .leading, spacing: 18) {
                // 1. System Size
                GroupBox(label: Label("System Geometry", systemImage: "grid")) {
                    VStack(alignment: .leading, spacing: 10) {
                        HStack {
                            VStack(alignment: .leading, spacing: 4) {
                                Text("NX (Cols)")
                                    .font(.caption)
                                    .foregroundColor(.secondary)
                                IntField(title: "100", value: $config.nx, width: 80)
                            }
                            VStack(alignment: .leading, spacing: 4) {
                                Text("NY (Rows)")
                                    .font(.caption)
                                    .foregroundColor(.secondary)
                                IntField(title: "100", value: $config.ny, width: 80)
                            }
                            Spacer()
                        }
                        
                        // Preset buttons
                        ScrollView(.horizontal, showsIndicators: false) {
                            HStack(spacing: 6) {
                                Text("Presets:")
                                    .font(.caption2)
                                    .foregroundColor(.secondary)
                                ForEach([(50, 50), (100, 100), (150, 150), (200, 200), (300, 300)], id: \.0) { size in
                                    Button("\(size.0)×\(size.1)") {
                                        config.setSize(nx: size.0, ny: size.1)
                                    }
                                    .buttonStyle(.bordered)
                                    .controlSize(.mini)
                                }
                            }
                        }
                    }
                    .frame(maxWidth: .infinity, alignment: .leading)
                    .padding(8)
                }
                
                // 2. Loading Mode & Remeshing
                GroupBox(label: Label("Mode & Triangulation", systemImage: "arrow.triangle.swap")) {
                    VStack(alignment: .leading, spacing: 10) {
                        Picker("Loading Mode", selection: $config.mode) {
                            ForEach(LoadingMode.allCases) { m in
                                Text(m.displayName).tag(m)
                            }
                        }
                        
                        if config.mode == .volterra {
                            Divider()
                            HStack {
                                Text("Free Core Radius R_free:")
                                    .font(.subheadline)
                                Spacer()
                                DoubleField(title: "35.0", value: $config.rFree, width: 75)
                            }
                            Toggle("Export Full Mesh (Unclipped VTK)", isOn: $config.exportFullMesh)
                                .font(.caption)
                            
                            Text(config.latticeSymmetry == .triangular ?
                                 "Uses Isotropic 2D Volterra solution (Eq. 15) with ν = C₁₂/C₁₁ and equilibrium spacing a₀ computed from active potential." :
                                 "Uses Anisotropic Stroh solution (Eqs. A5–A16) with (κ, ξ, ζ, g₁, g₂, g₃) and equilibrium spacing a₀ computed from active potential.")
                                .font(.caption2)
                                .foregroundColor(.accentColor)
                        }
                        
                        HStack {
                            Text("Random Seed:")
                                .font(.subheadline)
                            IntField(title: "Seed", value: $config.seed, width: 90)
                            Button(action: { config.randomizeSeed() }) {
                                Label("Random", systemImage: "dice")
                            }
                            .buttonStyle(.bordered)
                            .controlSize(.small)
                        }
                        
                        VStack(alignment: .leading, spacing: 4) {
                            HStack {
                                Text("Initial Noise (σ):")
                                    .font(.subheadline)
                                Spacer()
                                DoubleField(title: "0.04", value: $config.noiseLevel, width: 75, minValue: 0.0, maxValue: 1.0)
                            }
                            Slider(value: Binding(
                                get: { config.noiseLevel },
                                set: { config.noiseLevel = (round($0 * 1000) / 1000) }
                            ), in: 0.0...0.20)
                        }
                        
                        Toggle("Adaptive Remeshing (Delaunay Surgery)", isOn: $config.enableRemeshing)
                            .toggleStyle(.switch)
                        
                        Divider().padding(.vertical, 2)
                        
                        // Boundary Conditions Selector
                        VStack(alignment: .leading, spacing: 4) {
                            Text("Boundary Condition:")
                                .font(.subheadline)
                            Picker("Boundary Condition", selection: $config.boundaryCondition) {
                                ForEach(BoundaryConditionMode.allCases) { bc in
                                    Text(bc.displayName).tag(bc)
                                }
                            }
                            .pickerStyle(.segmented)
                        }
                    }
                    .frame(maxWidth: .infinity, alignment: .leading)
                    .padding(8)
                }
                
                // 3. Lattice & Potential Model
                GroupBox(label: Label("Lattice & Potential Model", systemImage: "atom")) {
                    VStack(alignment: .leading, spacing: 10) {
                        Picker("Potential", selection: $config.potentialMode) {
                            ForEach(PotentialMode.allCases) { pot in
                                Text(pot.displayName).tag(pot)
                            }
                        }
                        
                        Picker("Symmetry", selection: $config.latticeSymmetry) {
                            ForEach(LatticeSymmetry.allCases) { sym in
                                Text(sym.displayName).tag(sym)
                            }
                        }
                        
                        if config.potentialMode == .zanzotto {
                            HStack {
                                Text("Potential β:")
                                    .font(.subheadline)
                                Spacer()
                                DoubleField(title: "β", value: $config.beta, width: 80)
                            }
                            
                            HStack {
                                Text("Bulk Modulus K:")
                                    .font(.subheadline)
                                Spacer()
                                DoubleField(title: "K", value: $config.bulkK, width: 80)
                            }
                            
                            // Presets
                            HStack(spacing: 6) {
                                Text("Presets:")
                                    .font(.caption2)
                                    .foregroundColor(.secondary)
                                Button("Square (β=-0.25)") {
                                    config.latticeSymmetry = .square
                                    config.beta = -0.25
                                    config.bulkK = 5.0
                                }
                                .buttonStyle(.bordered)
                                .controlSize(.mini)
                                
                                Button("Triangular (β=4.0)") {
                                    config.latticeSymmetry = .triangular
                                    config.beta = 4.0
                                    config.bulkK = 5.0
                                }
                                .buttonStyle(.bordered)
                                .controlSize(.mini)
                            }
                            
                            Text("Reference element remains strictly square; deformation gradient F is measured relative to square reference.")
                                .font(.caption2)
                                .foregroundColor(.secondary)
                        } else {
                            VStack(alignment: .leading, spacing: 4) {
                                if config.potentialMode == .ljV2 {
                                    Text("Lennard-Jones v2: V(r) = 4[(1/r)¹² - (1/r)⁶], r_cut = 2.50")
                                        .font(.caption)
                                        .foregroundColor(.primary)
                                    Text("Ground state: b = 0.6872, mesh size h = b.")
                                        .font(.caption2)
                                        .foregroundColor(.secondary)
                                } else {
                                    Text("Lennard-Jones v3: V(r) = (1/r)¹² - 2(1/r)⁶, r_cut = 1.866")
                                        .font(.caption)
                                        .foregroundColor(.primary)
                                    Text("Ground state: b = 0.9964, mesh size h = b.")
                                        .font(.caption2)
                                        .foregroundColor(.secondary)
                                }
                                Text("Coarse-grained via Cauchy-Born with Lagrange reduction.")
                                    .font(.caption2)
                                    .foregroundColor(.secondary)
                            }
                            .padding(.vertical, 2)
                        }
                    }
                    .frame(maxWidth: .infinity, alignment: .leading)
                    .padding(8)
                }
                
                // 3b. Rotated Crystal Orientation / CSL
                GroupBox(label: Label("Rotated Crystal Orientation", systemImage: "rotate.right")) {
                    VStack(alignment: .leading, spacing: 10) {
                        Toggle(isOn: $config.useRotation) {
                            Text("Enable Rotated Lattice")
                                .font(.caption.bold())
                        }
                        .toggleStyle(.checkbox)
                        
                        if config.useRotation {
                            if config.boundaryCondition == .fixed || config.mode == .volterra {
                                VStack(alignment: .leading, spacing: 6) {
                                    Text("Arbitrary Continuous Rotation (Fixed BC):")
                                        .font(.caption2)
                                        .foregroundColor(.secondary)
                                    HStack {
                                        Text("Angle θ:")
                                            .font(.caption2.bold())
                                        Spacer()
                                        DoubleField(title: "0.0", value: $config.rotationAngle, width: 65, minValue: 0.0, maxValue: 360.0)
                                        Text("°").font(.caption2)
                                    }
                                    Slider(value: $config.rotationAngle, in: 0.0...360.0, step: 0.5)
                                    Text("Under fixed boundary conditions, the lattice can be rotated to any continuous angle θ ∈ [0, 360°].")
                                        .font(.caption2)
                                        .foregroundColor(.secondary)
                                }
                            } else {
                                // Periodic BC: CSL Discrete Angles
                                VStack(alignment: .leading, spacing: 8) {
                                    Text("Coincident Site Lattice (CSL) for PBC:")
                                        .font(.caption2)
                                        .foregroundColor(.secondary)
                                    
                                    // CSL Presets Picker
                                    let allowedCSL = CSLAngleOption.allowedAngles(for: config.latticeSymmetry)
                                    Picker("CSL Angle", selection: Binding(
                                        get: {
                                            allowedCSL.first(where: { $0.m == config.cslM && $0.n == config.cslN }) ?? allowedCSL[0]
                                        },
                                        set: { opt in
                                            config.cslM = opt.m
                                            config.cslN = opt.n
                                            config.rotationAngle = opt.theta_deg
                                            config.updateCSLSupercellAuto()
                                        }
                                    )) {
                                        ForEach(allowedCSL) { opt in
                                            Text(opt.displayName).tag(opt)
                                        }
                                    }
                                    .pickerStyle(.menu)
                                    .controlSize(.small)
                                    
                                    // Coprime (m, n) manual overrides
                                    HStack(spacing: 12) {
                                        HStack {
                                            Text("m:")
                                                .font(.caption2.bold())
                                            IntField(title: "m", value: $config.cslM, width: 40)
                                            Stepper("", value: Binding(
                                                get: { config.cslM },
                                                set: { config.cslM = $0; config.updateCSLSupercellAuto() }
                                            ), in: -20...20)
                                                .labelsHidden()
                                        }
                                        HStack {
                                            Text("n:")
                                                .font(.caption2.bold())
                                            IntField(title: "n", value: $config.cslN, width: 40)
                                            Stepper("", value: Binding(
                                                get: { config.cslN },
                                                set: { config.cslN = $0; config.updateCSLSupercellAuto() }
                                            ), in: 0...20)
                                                .labelsHidden()
                                        }
                                    }
                                    
                                    // Supercell Auto-scaling to Nx x Ny
                                    Toggle("Auto-scale Supercell to match \(config.nx)×\(config.ny)", isOn: $config.autoScaleSupercell)
                                        .font(.caption2)
                                        .onChange(of: config.autoScaleSupercell) { auto in
                                            if auto { config.updateCSLSupercellAuto() }
                                        }
                                    
                                    HStack {
                                        Text("Supercell: \(config.cslKx) × \(config.cslKy)")
                                            .font(.caption2.bold())
                                        Spacer()
                                        Text("~\(config.cslEstimatedAtoms) atoms")
                                            .font(.caption2)
                                            .foregroundColor(.accentColor)
                                    }
                                    
                                    if !config.autoScaleSupercell {
                                        // Manual Supercell Multipliers (Kx, Ky)
                                        HStack(spacing: 12) {
                                            HStack {
                                                Text("Kx:")
                                                    .font(.caption2.bold())
                                                IntField(title: "Kx", value: $config.cslKx, width: 40)
                                                Stepper("", value: $config.cslKx, in: 1...100)
                                                    .labelsHidden()
                                            }
                                            HStack {
                                                Text("Ky:")
                                                    .font(.caption2.bold())
                                                IntField(title: "Ky", value: $config.cslKy, width: 40)
                                                Stepper("", value: $config.cslKy, in: 1...100)
                                                    .labelsHidden()
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        
                        Button(action: {
                            previewManager.nx = config.nx
                            previewManager.ny = config.ny
                            previewManager.latticeSymmetry = config.latticeSymmetry
                            previewManager.potentialMode = config.potentialMode
                            previewManager.boundaryCondition = (config.mode == .volterra) ? .fixed : config.boundaryCondition
                            previewManager.beta = config.beta
                            previewManager.bulkK = config.bulkK
                            previewManager.useRotation = config.useRotation
                            previewManager.rotationAngle = config.rotationAngle
                            previewManager.cslM = config.cslM
                            previewManager.cslN = config.cslN
                            previewManager.cslKx = config.cslKx
                            previewManager.cslKy = config.cslKy
                            previewManager.generateInitialMeshPreview()
                            displayTab = .initialMesh
                        }) {
                            HStack {
                                if previewManager.isGeneratingPreview {
                                    ProgressView().scaleEffect(0.6).frame(width: 14, height: 14)
                                } else {
                                    Image(systemName: "eye.fill")
                                }
                                Text("Preview Initial Mesh")
                            }
                            .frame(maxWidth: .infinity)
                        }
                        .buttonStyle(.borderedProminent)
                        .controlSize(.small)
                        .padding(.top, 2)
                    }
                    .frame(maxWidth: .infinity, alignment: .leading)
                    .padding(8)
                }
                
                // 4. Preconditioning & Energy Relaxation
                GroupBox(label: Label("Compute Device & Preconditioner", systemImage: "bolt.fill")) {
                    VStack(alignment: .leading, spacing: 10) {
                        VStack(alignment: .leading, spacing: 4) {
                            Text("Compute Device:")
                                .font(.subheadline)
                            Picker("Device", selection: $config.computeDevice) {
                                ForEach(ComputeDevice.allCases) { dev in
                                    Text(dev.displayName).tag(dev)
                                }
                            }
                            .pickerStyle(.segmented)
                        }
                        
                        if config.computeDevice == .gpu {
                            HStack {
                                Text("GPU Gradient Tol:")
                                    .font(.subheadline)
                                Spacer()
                                DoubleField(title: "1e-4", value: $config.gpuTol, width: 80)
                            }
                            .help("Stopping tolerance for GPU L-BFGS iterations.")
                            .padding(.top, 2)
                        }
                        
                        Divider()
                        
                        Picker("Metric", selection: $config.preconditioner) {
                            ForEach(PreconditionerType.allCases) { p in
                                Text(p.displayName).tag(p)
                            }
                        }
                        
                        if config.preconditioner != .none {
                            HStack {
                                Text("Active from step:")
                                    .font(.subheadline)
                                Spacer()
                                IntField(title: "1", value: $config.precondFromStep, width: 70)
                            }
                            HStack {
                                Text("Refresh every N steps (0=adaptive):")
                                    .font(.subheadline)
                                Spacer()
                                IntField(title: "0", value: $config.precondRefresh, width: 70)
                            }
                        }
                    }
                    .padding(8)
                }
                
                // 5. Loading Schedule vs Dislocation Cylinder Settings
                if config.mode == .volterra {
                    GroupBox(label: Label("Dislocation Cylinder Settings", systemImage: "circle.dashed")) {
                        VStack(alignment: .leading, spacing: 10) {
                            HStack {
                                Text("Relaxation radius R_free (in h):")
                                    .font(.subheadline)
                                Spacer()
                                DoubleField(title: "10.0", value: $config.rFree, width: 75)
                            }
                            
                            Toggle("Export Full Mesh (all atoms)", isOn: $config.exportFullMesh)
                                .toggleStyle(.switch)
                            
                            Text("Note: Loading schedule is inactive in Volterra mode. The Volterra field is installed directly, and atoms within R_free relax by energy minimization.")
                                .font(.caption2)
                                .foregroundColor(.secondary)
                        }
                        .padding(8)
                    }
                } else {
                    GroupBox(label: Label("Loading Schedule (|α|)", systemImage: "chart.line.uptrend.xyaxis")) {
                        VStack(alignment: .leading, spacing: 8) {
                            HStack {
                                Text("Start strain |α_start|:")
                                    .font(.subheadline)
                                Spacer()
                                DoubleField(title: "0.14", value: $config.alphaStart, width: 80)
                            }
                            HStack {
                                Text("End strain |α_end|:")
                                    .font(.subheadline)
                                Spacer()
                                DoubleField(title: "1.0", value: $config.alphaEnd, width: 80)
                            }
                            HStack {
                                Text("Step size |dα|:")
                                    .font(.subheadline)
                                Spacer()
                                DoubleField(title: "0.00006", value: $config.stepSize, width: 90)
                            }
                        }
                        .padding(8)
                    }
                }
                
                // 6. Stability Monitor & Eigenvalues
                GroupBox(label: Label("Stability & Eigenvalues", systemImage: "waveform.path.ecg")) {
                    VStack(alignment: .leading, spacing: 10) {
                        Toggle("Compute Lowest Stiffness Eigenvalues", isOn: $config.enableEigenvalues)
                            .toggleStyle(.switch)
                        
                        if config.enableEigenvalues {
                            HStack {
                                Text("Every N steps:")
                                Spacer()
                                IntField(title: "5", value: $config.eigEvery, width: 60)
                            }
                            Toggle("Compute at avalanches / stress jumps", isOn: $config.eigAtAvalanche)
                            HStack {
                                Text("Modes to compute:")
                                Spacer()
                                Stepper("\(config.eigModes)", value: $config.eigModes, in: 1...20)
                            }
                            HStack {
                                Text("Soft modes exported:")
                                Spacer()
                                Stepper("\(config.eigVectors)", value: $config.eigVectors, in: 0...10)
                            }
                            Picker("Solver", selection: $config.eigSolver) {
                                ForEach(EigenSolverType.allCases) { s in
                                    Text(s.displayName).tag(s)
                                }
                            }
                        }
                    }
                    .padding(8)
                }
                
                // 7. Data Saving & Output Controls
                GroupBox(label: Label("Output & Checkpoint Controls", systemImage: "externaldrive")) {
                    VStack(alignment: .leading, spacing: 8) {
                        Toggle("Inner iterations trace (--trace=1)", isOn: $config.traceInnerIterations)
                        Toggle("Save Surgery VTK Snapshots (--surgery-vtk)", isOn: $config.saveSurgeryVTK)
                        Toggle("Save Legacy triangle_data", isOn: $config.saveTriangleData)
                        
                        HStack {
                            Text("Checkpoint every N steps:")
                            Spacer()
                            IntField(title: "500", value: $config.checkpointInterval, width: 70)
                        }
                        HStack {
                            Text("Stress drop threshold Δσ/σ:")
                            Spacer()
                            DoubleField(title: "0.10", value: $config.stressDropThreshold, width: 70)
                        }
                    }
                    .padding(8)
                }
                
                // 7. Hardware & Threads
                GroupBox(label: Label("Hardware & CPU Threads", systemImage: "cpu")) {
                    VStack(alignment: .leading, spacing: 8) {
                        HStack {
                            Text("OpenMP threads (FEM loops):")
                            Spacer()
                            Stepper("\(config.ompThreads)", value: $config.ompThreads, in: 1...64)
                        }
                        HStack {
                            Text("BLAS threads (CHOLMOD):")
                            Spacer()
                            Stepper("\(config.blasThreads)", value: $config.blasThreads, in: 1...64)
                        }
                    }
                    .padding(8)
                }
                
                // 8. Working & Output Directory
                GroupBox(label: Label("Output Directory", systemImage: "folder.fill")) {
                    VStack(alignment: .leading, spacing: 6) {
                        Text(config.targetDirectory)
                            .font(.system(size: 10, design: .monospaced))
                            .lineLimit(2)
                            .truncationMode(.middle)
                            .foregroundColor(.secondary)
                        
                        HStack {
                            Button("Change Base Folder...") {
                                let panel = NSOpenPanel()
                                panel.canChooseFiles = false
                                panel.canChooseDirectories = true
                                panel.allowsMultipleSelection = false
                                if panel.runModal() == .OK, let url = panel.url {
                                    config.projectRoot = url.path
                                }
                            }
                            .buttonStyle(.bordered)
                            .controlSize(.small)
                            
                            Spacer()
                            
                            Button("Reset") {
                                config.projectRoot = SimulationConfig.detectProjectRoot()
                            }
                            .buttonStyle(.borderless)
                            .controlSize(.small)
                        }
                    }
                    .padding(8)
                }
            }
            .frame(maxWidth: .infinity, alignment: .leading)
            .padding(14)
        }
        .background(Color(NSColor.controlBackgroundColor).opacity(0.5))
    }
}

// MARK: - Unified Live Status & View Switcher Header
struct LiveStatusHeaderView: View {
    @ObservedObject var config: SimulationConfig
    @Binding var displayTab: MainDisplayTab
    @EnvironmentObject var runner: ProcessRunner
    
    var body: some View {
        HStack(spacing: 10) {
            // 1. Pinned View Switcher (Curve vs Terminal) — never moves or clips
            Picker("View", selection: $displayTab) {
                ForEach(MainDisplayTab.allCases) { tab in
                    Label(tab.rawValue, systemImage: tab.icon).tag(tab)
                }
            }
            .labelsHidden()
            .pickerStyle(.segmented)
            .controlSize(.small)
            .fixedSize()
            
            // 2. Pinned Concise Status Pill
            HStack(spacing: 5) {
                Circle()
                    .fill(statusColor)
                    .frame(width: 8, height: 8)
                Text(shortStatusLabel)
                    .font(.system(size: 11, weight: .semibold))
                    .lineLimit(1)
            }
            .padding(.horizontal, 8)
            .padding(.vertical, 5)
            .background(Color(NSColor.quaternaryLabelColor).opacity(0.6))
            .cornerRadius(6)
            .fixedSize()
            .help(runner.status.label)
            
            Divider()
                .frame(height: 18)
            
            // 3. Horizontally scrollable Live Metrics Pills so they never push or clip the header
            ScrollView(.horizontal, showsIndicators: false) {
                HStack(spacing: 6) {
                    MetricPill(label: "Step", value: "\(runner.currentStep)")
                    MetricPill(label: "Strain |α|", value: String(format: "%.4f", runner.currentAlpha))
                    MetricPill(label: "Stress σ", value: String(format: "%.4f", runner.currentStress))
                    if let lmin = runner.currentLambdaMin {
                        MetricPill(
                            label: "λ_min",
                            value: String(format: "%.2e", lmin),
                            valueColor: lmin <= 0 ? .red : .teal
                        )
                    }
                    MetricPill(label: "Avalanches", value: "\(runner.avalancheCount)")
                    MetricPill(label: "Elapsed", value: formatTime(runner.elapsedTime))
                }
                .padding(.vertical, 1)
            }
            .frame(maxWidth: .infinity, alignment: .leading)
        }
        .frame(maxWidth: .infinity, alignment: .leading)
    }
    
    private var shortStatusLabel: String {
        switch runner.status {
        case .idle: return "Ready"
        case .running: return "Running"
        case .completed(let code, _): return code == 0 ? "Completed" : "Exit \(code)"
        case .failed: return "Failed"
        case .terminated: return "Stopped"
        }
    }
    
    private var statusColor: Color {
        switch runner.status {
        case .idle: return .gray
        case .running: return .green
        case .completed(let code, _): return code == 0 ? .blue : .orange
        case .failed, .terminated: return .red
        }
    }
    
    private func formatTime(_ sec: TimeInterval) -> String {
        let m = Int(sec) / 60
        let s = Int(sec) % 60
        return String(format: "%02d:%02d", m, s)
    }
}

struct MetricPill: View {
    let label: String
    let value: String
    var valueColor: Color = .primary
    
    var body: some View {
        VStack(spacing: 1) {
            Text(label)
                .font(.system(size: 9, weight: .semibold))
                .foregroundColor(.secondary)
            Text(value)
                .font(.system(size: 12, weight: .bold, design: .monospaced))
                .foregroundColor(valueColor)
        }
        .padding(.horizontal, 8)
        .padding(.vertical, 4)
        .background(Color(NSColor.quaternaryLabelColor))
        .cornerRadius(6)
        .lineLimit(1)
        .fixedSize()
    }
}

// MARK: - Terminal View
struct TerminalView: View {
    let logLines: [LogMessage]
    @Binding var filterText: String
    @Binding var autoScroll: Bool
    
    var filteredLines: [LogMessage] {
        if filterText.isEmpty { return logLines }
        return logLines.filter { $0.text.localizedCaseInsensitiveContains(filterText) }
    }
    
    var body: some View {
        VStack(spacing: 0) {
            // Search & Tools Bar
            HStack {
                Image(systemName: "magnifyingglass")
                    .foregroundColor(.secondary)
                TextField("Filter output...", text: $filterText)
                    .textFieldStyle(.plain)
                    .frame(maxWidth: 220)
                
                Spacer()
                
                Toggle("Auto-scroll", isOn: $autoScroll)
                    .toggleStyle(.checkbox)
                    .controlSize(.small)
                
                Text("\(logLines.count) lines")
                    .font(.caption2)
                    .foregroundColor(.secondary)
            }
            .padding(.horizontal, 12)
            .padding(.vertical, 6)
            .background(Color(NSColor.controlBackgroundColor))
            
            Divider()
            
            // Console Log ScrollView
            ScrollViewReader { proxy in
                ScrollView {
                    LazyVStack(alignment: .leading, spacing: 2) {
                        ForEach(filteredLines) { msg in
                            Text(msg.text)
                                .font(.system(size: 11, design: .monospaced))
                                .foregroundColor(msg.isError ? Color.red : Color(NSColor.textColor))
                                .textSelection(.enabled)
                                .id(msg.id)
                        }
                    }
                    .padding(8)
                }
                .background(Color(NSColor.textBackgroundColor))
                .onChange(of: logLines.count) { _ in
                    if autoScroll, let last = filteredLines.last {
                        proxy.scrollTo(last.id, anchor: .bottom)
                    }
                }
            }
        }
    }
}

// MARK: - Bottom Action Bar
struct BottomActionBar: View {
    @ObservedObject var config: SimulationConfig
    @EnvironmentObject var runner: ProcessRunner
    @State private var copiedCLI: Bool = false
    
    var body: some View {
        HStack(spacing: 10) {
            // Clear Console
            Button(action: { runner.clearLogs() }) {
                Label("Clear Log", systemImage: "trash")
            }
            .controlSize(.regular)
            
            // Reveal in Finder
            Button(action: { runner.openOutputFolder() }) {
                Label("Open Output Folder", systemImage: "folder")
            }
            .controlSize(.regular)
            
            // Copy CLI Command
            Button(action: {
                NSPasteboard.general.clearContents()
                NSPasteboard.general.setString(config.fullCommandLine(), forType: .string)
                copiedCLI = true
                DispatchQueue.main.asyncAfter(deadline: .now() + 1.5) {
                    copiedCLI = false
                }
            }) {
                Label(copiedCLI ? "Copied CLI!" : "Copy CLI", systemImage: copiedCLI ? "checkmark" : "terminal")
            }
            .controlSize(.regular)
            .help(config.fullCommandLine())
            
            Spacer()
            
            if runner.status.isRunning {
                Button(role: .destructive, action: { runner.stop() }) {
                    Label("Stop Simulation", systemImage: "stop.fill")
                        .font(.headline)
                        .padding(.horizontal, 8)
                }
                .buttonStyle(.borderedProminent)
                .tint(.red)
                .controlSize(.large)
            } else {
                Button(action: { runner.start(config: config) }) {
                    Label("Launch Simulation", systemImage: "play.fill")
                        .font(.headline)
                        .padding(.horizontal, 12)
                }
                .buttonStyle(.borderedProminent)
                .tint(.accentColor)
                .controlSize(.large)
            }
        }
    }
}
