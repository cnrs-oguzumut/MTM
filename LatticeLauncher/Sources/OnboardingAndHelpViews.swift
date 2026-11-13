import SwiftUI
import AppKit

// MARK: - UI Complexity Mode (Progressive Disclosure)
public enum UIComplexityMode: String, CaseIterable, Identifiable {
    case guided = "guided"
    case expert = "expert"
    
    public var id: String { rawValue }
    
    public var title: String {
        switch self {
        case .guided: return "Guided"
        case .expert: return "Expert"
        }
    }
    
    public var badgeLabel: String {
        switch self {
        case .guided: return "🌱 Guided"
        case .expert: return "⚙️ Expert"
        }
    }
    
    public var subtitle: String {
        switch self {
        case .guided:
            return "Step-by-step workflow with solver internals tucked in an expandable drawer."
        case .expert:
            return "All 35+ research parameters, solver metrics, and thread controls expanded."
        }
    }
}

// MARK: - Notification Names for Menu Bar Commands
public extension Notification.Name {
    static let showWelcomeTour = Notification.Name("LatticeLauncher.showWelcomeTour")
    static let showDocumentationCenter = Notification.Name("LatticeLauncher.showDocumentationCenter")
    static let applySampleExperiment = Notification.Name("LatticeLauncher.applySampleExperiment")
    static let switchMainDisplayTab = Notification.Name("LatticeLauncher.switchMainDisplayTab")
    static let toggleConsoleDrawer = Notification.Name("LatticeLauncher.toggleConsoleDrawer")
    static let openVTKFile = Notification.Name("LatticeLauncher.openVTKFile")
    static let openVTKFilePanel = Notification.Name("LatticeLauncher.openVTKFilePanel")
    static let toggleComplexityMode = Notification.Name("LatticeLauncher.toggleComplexityMode")
    static let runSimulationCommand = Notification.Name("LatticeLauncher.runSimulationCommand")
    static let stopSimulationCommand = Notification.Name("LatticeLauncher.stopSimulationCommand")
}

// MARK: - Reusable Inline '?' Contextual Help Popover Button
public struct HelpPopoverButton: View {
    public let title: String
    public let summary: String
    public var formulaOrTip: String? = nil
    public var docTopic: DocumentationTopic? = nil
    public var onOpenDocs: ((DocumentationTopic) -> Void)? = nil
    
    @State private var isPresented: Bool = false
    
    public init(
        title: String,
        summary: String,
        formulaOrTip: String? = nil,
        docTopic: DocumentationTopic? = nil,
        onOpenDocs: ((DocumentationTopic) -> Void)? = nil
    ) {
        self.title = title
        self.summary = summary
        self.formulaOrTip = formulaOrTip
        self.docTopic = docTopic
        self.onOpenDocs = onOpenDocs
    }
    
    public var body: some View {
        Button(action: { isPresented.toggle() }) {
            Image(systemName: "questionmark.circle")
                .font(.system(size: 12, weight: .semibold))
                .foregroundColor(BlueHour.ice.opacity(0.85))
        }
        .buttonStyle(.plain)
        .help("Click for explanation: \(title)")
        .popover(isPresented: $isPresented, arrowEdge: .trailing) {
            VStack(alignment: .leading, spacing: 8) {
                HStack(spacing: 6) {
                    Image(systemName: "lightbulb.fill")
                        .foregroundColor(Color(red: 1.0, green: 0.82, blue: 0.35))
                    Text(title)
                        .font(.system(size: 12.5, weight: .bold))
                        .foregroundColor(BlueHour.ice)
                }
                
                Text(summary)
                    .font(.system(size: 11.5))
                    .foregroundColor(.primary)
                    .fixedSize(horizontal: false, vertical: true)
                
                if let tip = formulaOrTip, !tip.isEmpty {
                    Text(tip)
                        .font(.system(size: 10.5, weight: .medium, design: .monospaced))
                        .foregroundColor(Color(red: 0.45, green: 0.92, blue: 0.62))
                        .padding(7)
                        .frame(maxWidth: .infinity, alignment: .leading)
                        .background(BlueHour.midnight.opacity(0.75), in: RoundedRectangle(cornerRadius: 6))
                }
                
                if let topic = docTopic, let openAction = onOpenDocs {
                    Divider()
                    Button(action: {
                        isPresented = false
                        openAction(topic)
                    }) {
                        Label("Read Full Guide: \(topic.title)", systemImage: "book.fill")
                            .font(.system(size: 11, weight: .semibold))
                    }
                    .buttonStyle(.link)
                }
            }
            .padding(12)
            .frame(width: 320)
            .background(BlueHour.windowBackground)
        }
    }
}

// MARK: - Curated 1-Click Sample Experiments Catalog
public enum SampleExperimentPreset: String, CaseIterable, Identifiable {
    case homogeneousShearAvalanches = "homogeneous_shear_avalanches"
    case periodicDislocations6D = "periodic_dislocations_6d"
    case periodicDislocations60D = "periodic_dislocations_60d"
    case sigma5SymmetricTiltGB = "sigma5_symmetric_tilt_gb"
    case gbEnergyReadShockleySweep = "gb_energy_sweep"
    case nanoindentationCrystal = "nanoindentation_crystal"
    case plateWithHoleTension = "plate_with_hole_tension"
    case singleVolterraDislocation = "single_volterra_dislocation"
    
    public var id: String { rawValue }
    
    public var subtitle: String { difficultyBadge }
    
    public var title: String {
        switch self {
        case .homogeneousShearAvalanches:
            return "1. Crystal Shear & Plastic Avalanches"
        case .periodicDislocations6D:
            return "2. Periodic Dislocation Dipoles (6 Planes)"
        case .periodicDislocations60D:
            return "3. Dense DDD Dislocation Ensemble (60 Dislocations)"
        case .sigma5SymmetricTiltGB:
            return "4. Special Σ5 (310) Symmetric Tilt Grain Boundary"
        case .gbEnergyReadShockleySweep:
            return "5. Grain Boundary Energy Curve γ_GB(Δθ)"
        case .nanoindentationCrystal:
            return "6. Cylindrical Nanoindentation on Free Surface"
        case .plateWithHoleTension:
            return "7. Plate with Circular Hole under Tension"
        case .singleVolterraDislocation:
            return "8. Single Edge Dislocation Core (Volterra / Stroh)"
        }
    }
    
    public var difficultyBadge: String {
        switch self {
        case .homogeneousShearAvalanches, .singleVolterraDislocation:
            return "Beginner"
        case .periodicDislocations6D, .sigma5SymmetricTiltGB, .plateWithHoleTension:
            return "Intermediate"
        case .periodicDislocations60D, .gbEnergyReadShockleySweep, .nanoindentationCrystal:
            return "Advanced"
        }
    }
    
    public var difficultyColor: Color {
        switch difficultyBadge {
        case "Beginner": return Color(red: 0.45, green: 0.90, blue: 0.62)
        case "Intermediate": return Color(red: 0.31, green: 0.85, blue: 0.78)
        default: return Color(red: 1.0, green: 0.78, blue: 0.35)
        }
    }
    
    public var icon: String {
        switch self {
        case .homogeneousShearAvalanches: return "waveform.path.ecg"
        case .periodicDislocations6D: return "arrow.left.and.right.square"
        case .periodicDislocations60D: return "square.grid.3x3.middle.filled"
        case .sigma5SymmetricTiltGB: return "square.split.diagonal.2x2"
        case .gbEnergyReadShockleySweep: return "chart.xyaxis.line"
        case .nanoindentationCrystal: return "arrow.down.to.line.compact"
        case .plateWithHoleTension: return "circle.dashed.inset.filled"
        case .singleVolterraDislocation: return "scope"
        }
    }
    
    public var physicsSummary: String {
        switch self {
        case .homogeneousShearAvalanches:
            return "Shears a 100×100 periodic Square crystal with weak thermal disorder (σ=0.04) and adaptive Delaunay surgery to observe elastic loading, homogeneous dislocation nucleation, and serrated stress-drop avalanches."
        case .periodicDislocations6D:
            return "Seeds M=3 balanced ±a₀ edge dislocation pairs (6 active glide planes, 1 dislocation per plane) on a 100×100 periodic torus. Stage 1 relaxes cores at α=0; Stage 2 shears the periodic box."
        case .periodicDislocations60D:
            return "Seeds M=30 balanced ±a₀ pairs (60 dislocations on 60 distinct glide planes) inside a 240×240 periodic crystal with randomly shuffled ±b signs."
        case .sigma5SymmetricTiltGB:
            return "Constructs an exact Σ5 (310) (±18.43°, Δθ=36.87°) vertical bicrystal in a commensurate box, merges overlapping seam atoms, relaxes the GB core at α=0, and shears parallel to the boundary."
        case .gbEnergyReadShockleySweep:
            return "Automatically sweeps misorientation angles Δθ across low-Σ cusps (Σ5, Σ13, Σ17) at α=0 to plot the Read–Shockley grain boundary energy curve γ_GB(Δθ)."
        case .nanoindentationCrystal:
            return "Presses a rigid cylindrical nanoindenter tip (R=12a₀) into the traction-free top surface of a 100×100 crystal after Stage 1 equilibrium relaxation."
        case .plateWithHoleTension:
            return "Carves a central circular hole (R=0.18L) inside a 100×100 plate and applies uniaxial horizontal tension to nucleate dislocations from the hole rim."
        case .singleVolterraDislocation:
            return "Installs a single anisotropic Stroh edge dislocation at (0,0), relaxes the core inside R_free=35a₀, and plots radial logarithmic energy E(r) ~ A ln(r/r₀)."
        }
    }
    
    public var whatToObserve: String {
        switch self {
        case .homogeneousShearAvalanches:
            return "Watch the live Stress–Strain curve drop sharply during plastic avalanches and open VTK snapshots in ParaView."
        case .periodicDislocations6D, .periodicDislocations60D:
            return "Inspect the ⊥/⊤ glide planes in Geometry & BCs Studio, then watch dipolar passing and annihilation under shear."
        case .sigma5SymmetricTiltGB:
            return "See the periodic kite-shaped structural units along the vertical GB seam in Geometry & BCs Studio."
        case .gbEnergyReadShockleySweep:
            return "Observe the live γ_GB(Δθ) curve populate with energy cusps at special Σ5 (36.87° & 53.13°) orientations."
        case .nanoindentationCrystal:
            return "Inspect the yellow indenter contact arc in Geometry Studio and watch load–depth P(δ) bursts."
        case .plateWithHoleTension:
            return "Observe stress concentration σ_xx at the top/bottom poles of the circular hole."
        case .singleVolterraDislocation:
            return "Inspect the 4-panel Dislocation Profiles view comparing unrelaxed vs. relaxed core energy."
        }
    }
    
    public func apply(to config: SimulationConfig) -> Bool {
        switch self {
        case .homogeneousShearAvalanches:
            config.setSize(nx: 100, ny: 100)
            config.latticeSymmetry = .square
            config.potentialMode = .zanzotto
            config.beta = -0.25
            config.bulkK = 5.0
            config.mode = .positive
            config.noiseLevel = 0.04
            config.enableRemeshing = true
            config.useRotation = false
            config.boundaryCondition = .periodic
            config.geometryConfig.enabled = false
            config.geometryConfig.microstructure.mode = .singleCrystal
            config.protocolMode = .relaxThenLoad
            config.alphaStart = 0.0
            config.alphaEnd = 0.20
            config.stepSize = 0.0005
            return false // show curve tab
            
        case .periodicDislocations6D:
            config.setSize(nx: 100, ny: 100)
            config.latticeSymmetry = .square
            config.potentialMode = .zanzotto
            config.beta = -0.25
            config.bulkK = 5.0
            config.mode = .positive
            config.noiseLevel = 0.0
            config.enableRemeshing = false
            config.useRotation = false
            config.boundaryCondition = .periodic
            config.geometryConfig.enabled = true
            config.geometryConfig.outer_shape = .rectangle
            config.geometryConfig.cutouts = []
            config.geometryConfig.indenter.enabled = false
            config.geometryConfig.microstructure.mode = .periodicDislocations
            config.geometryConfig.microstructure.dislocation_pair_count = 3
            config.geometryConfig.microstructure.dislocation_layout = DislocationLayoutPreset.equispacedRandomSigns.rawValue
            config.geometryConfig.microstructure.regeneratePeriodicDislocations(nx: 100, ny: 100)
            config.protocolMode = .relaxThenLoad
            config.alphaStart = 0.0
            config.alphaEnd = 0.10
            config.stepSize = 0.0005
            return true // show geometry studio
            
        case .periodicDislocations60D:
            config.setSize(nx: 240, ny: 240)
            config.latticeSymmetry = .square
            config.potentialMode = .zanzotto
            config.beta = -0.25
            config.bulkK = 5.0
            config.mode = .positive
            config.noiseLevel = 0.0
            config.enableRemeshing = false
            config.useRotation = false
            config.boundaryCondition = .periodic
            config.geometryConfig.enabled = true
            config.geometryConfig.outer_shape = .rectangle
            config.geometryConfig.cutouts = []
            config.geometryConfig.indenter.enabled = false
            config.geometryConfig.microstructure.mode = .periodicDislocations
            config.geometryConfig.microstructure.dislocation_pair_count = 30
            config.geometryConfig.microstructure.dislocation_layout = DislocationLayoutPreset.equispacedRandomSigns.rawValue
            config.geometryConfig.microstructure.regeneratePeriodicDislocations(nx: 240, ny: 240)
            config.protocolMode = .relaxThenLoad
            config.alphaStart = 0.0
            config.alphaEnd = 0.10
            config.stepSize = 0.0005
            return true
            
        case .sigma5SymmetricTiltGB:
            let match = CSLRationalizer.makeMatch(m: 3, n: 1)
            let boxN = match.commensurateSize(periods: 30)
            config.setSize(nx: boxN, ny: boxN)
            config.latticeSymmetry = .square
            config.potentialMode = .zanzotto
            config.beta = -0.25
            config.bulkK = 5.0
            config.mode = .positive
            config.noiseLevel = 0.0
            config.enableRemeshing = true
            config.geometryConfig.enabled = true
            config.geometryConfig.outer_shape = .rectangle
            config.geometryConfig.cutouts = []
            config.geometryConfig.indenter.enabled = false
            config.geometryConfig.microstructure.mode = .bicrystal
            config.geometryConfig.microstructure.gb_angle_deg = 0.0
            config.geometryConfig.microstructure.theta1_deg = -match.halfAngleDeg
            config.geometryConfig.microstructure.theta2_deg = +match.halfAngleDeg
            config.geometryConfig.microstructure.min_seam_distance_a0 = match.recommendedDMinA0
            config.geometryConfig.microstructure.gb_normal_gap_a0 = match.recommendedGapA0
            config.geometryConfig.boundary_rules = [
                BoundaryZoneRuleModel(region: .left, kind: .fixed, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                BoundaryZoneRuleModel(region: .right, kind: .prescribed, thickness_a0: 2.0, ux: 0.0, uy: 1.0, lock_x: true, lock_y: true),
                BoundaryZoneRuleModel(region: .bottom, kind: .free, thickness_a0: 1.5),
                BoundaryZoneRuleModel(region: .top, kind: .free, thickness_a0: 1.5)
            ]
            config.geometryConfig.primary_stress_component = "sigma_xy"
            config.protocolMode = .relaxThenLoad
            config.alphaStart = 0.0
            config.alphaEnd = 0.12
            config.stepSize = 0.0005
            return true
            
        case .gbEnergyReadShockleySweep:
            config.setSize(nx: 95, ny: 95)
            config.latticeSymmetry = .square
            config.potentialMode = .zanzotto
            config.beta = -0.25
            config.bulkK = 5.0
            config.mode = .positive
            config.noiseLevel = 0.0
            config.enableRemeshing = true
            config.geometryConfig.enabled = true
            config.geometryConfig.outer_shape = .rectangle
            config.geometryConfig.cutouts = []
            config.geometryConfig.indenter.enabled = false
            config.geometryConfig.microstructure.mode = .bicrystal
            config.geometryConfig.microstructure.gb_angle_deg = 0.0
            config.protocolMode = .sweepMisorientation
            config.sweepContinuousMode = true
            config.sweepDeltaThetaStepDeg = 2.5
            config.sweepIncludeExactCusps = true
            config.sweepAutoCommensurateBox = true
            return false
            
        case .nanoindentationCrystal:
            config.setSize(nx: 100, ny: 100)
            config.latticeSymmetry = .square
            config.potentialMode = .zanzotto
            config.beta = -0.25
            config.bulkK = 5.0
            config.mode = .positive
            config.noiseLevel = 0.0
            config.enableRemeshing = true
            config.geometryConfig.enabled = true
            config.geometryConfig.microstructure.mode = .singleCrystal
            config.geometryConfig.applyPreset(.nanoIndentation)
            config.protocolMode = .relaxThenLoad
            config.alphaStart = 0.0
            config.alphaEnd = 1.0
            config.stepSize = 0.01
            return true
            
        case .plateWithHoleTension:
            config.setSize(nx: 100, ny: 100)
            config.latticeSymmetry = .square
            config.potentialMode = .zanzotto
            config.beta = -0.25
            config.bulkK = 5.0
            config.mode = .positive
            config.noiseLevel = 0.0
            config.enableRemeshing = true
            config.geometryConfig.enabled = true
            config.geometryConfig.microstructure.mode = .singleCrystal
            config.geometryConfig.applyPreset(.plateWithHoleTension)
            config.protocolMode = .relaxThenLoad
            config.alphaStart = 0.0
            config.alphaEnd = 0.12
            config.stepSize = 0.0005
            return true
            
        case .singleVolterraDislocation:
            config.setSize(nx: 100, ny: 100)
            config.latticeSymmetry = .square
            config.potentialMode = .zanzotto
            config.beta = -0.25
            config.bulkK = 5.0
            config.mode = .volterra
            config.rFree = 35.0
            config.exportFullMesh = true
            config.relaxMacroStress = true
            config.enableEigenvalues = true
            return false
        }
    }
}

// MARK: - Documentation Center Topics
public enum DocumentationTopic: String, CaseIterable, Identifiable {
    case quickStart = "quick_start"
    case installingDislocations = "installing_dislocations"
    case grainBoundariesCSL = "grain_boundaries_csl"
    case boundaryConditionsIndenter = "bcs_and_indenter"
    case acousticAndStiffness = "acoustic_and_stiffness"
    case parameterGlossary = "parameter_glossary"
    
    public static var dislocationsGuide: DocumentationTopic { .installingDislocations }
    public static var grainBoundaries: DocumentationTopic { .grainBoundariesCSL }
    public static var geometryAndBCs: DocumentationTopic { .boundaryConditionsIndenter }
    public static var potentialsAndRemeshing: DocumentationTopic { .parameterGlossary }
    public static var solverAndEigenvalues: DocumentationTopic { .parameterGlossary }
    
    public var id: String { rawValue }
    
    public var title: String {
        switch self {
        case .quickStart: return "1. Quick Start & 3-Step Workflow"
        case .installingDislocations: return "2. Installing Edge Dislocations (2M Periodic & Volterra)"
        case .grainBoundariesCSL: return "3. Bicrystals, Voronoi & CSL Rationalizer"
        case .boundaryConditionsIndenter: return "4. Two-Stage Boundary Conditions & Nanoindenter"
        case .acousticAndStiffness: return "5. Acoustic Tensor & Stiffness Spectrum K(α)"
        case .parameterGlossary: return "6. Glossary of Physical Symbols & Parameters"
        }
    }
    
    public var icon: String {
        switch self {
        case .quickStart: return "sparkles"
        case .installingDislocations: return "arrow.left.and.right.square"
        case .grainBoundariesCSL: return "square.split.diagonal.2x2"
        case .boundaryConditionsIndenter: return "slider.horizontal.below.square.and.square.filled"
        case .acousticAndStiffness: return "waveform.path.ecg"
        case .parameterGlossary: return "character.book.closed.fill"
        }
    }
}

// MARK: - Welcome Onboarding & Sample Experiments Sheet
public struct WelcomeOnboardingSheet: View {
    @ObservedObject var config: SimulationConfig
    @Binding var isPresented: Bool
    @AppStorage("showWelcomeOnStartup_v2") private var showOnStartup: Bool = true
    var onSelectExperiment: (SampleExperimentPreset, Bool) -> Void
    var onOpenDocs: (DocumentationTopic) -> Void
    
    @State private var selectedTab: Int = 0 // 0 = Interactive Tour + Quick Examples, 1 = Full 8-Experiment Gallery
    
    public var body: some View {
        VStack(spacing: 0) {
            // Top Hero Banner
            HStack(spacing: 14) {
                ZStack {
                    RoundedRectangle(cornerRadius: 12)
                        .fill(
                            LinearGradient(
                                colors: [Color(red: 0.12, green: 0.64, blue: 0.54), BlueHour.midnight],
                                startPoint: .topLeading,
                                endPoint: .bottomTrailing
                            )
                        )
                        .frame(width: 52, height: 52)
                    Image(systemName: "atom")
                        .font(.system(size: 28, weight: .bold))
                        .foregroundColor(.white)
                }
                
                VStack(alignment: .leading, spacing: 3) {
                    HStack(spacing: 8) {
                        Text("Welcome to Lattice Launcher")
                            .font(.system(size: 20, weight: .bold))
                            .foregroundColor(BlueHour.ice)
                        Text("v2.0 Studio")
                            .font(.system(size: 10, weight: .bold, design: .monospaced))
                            .padding(.horizontal, 7)
                            .padding(.vertical, 2)
                            .background(Color(red: 0.12, green: 0.64, blue: 0.54), in: Capsule())
                            .foregroundColor(.white)
                    }
                    Text("Mesoscopic Tensorial Model (MTM) for Crystal Plasticity, Dislocations, Grain Boundaries & Phase Transformations")
                        .font(.system(size: 12))
                        .foregroundColor(BlueHour.mist)
                }
                
                Spacer()
                
                BlueHourSegmentedPicker(
                    selection: $selectedTab,
                    items: [
                        BlueHourTabItem(0, title: "Quick-Start Tour", icon: "sparkles"),
                        BlueHourTabItem(1, title: "Sample Experiments (\(SampleExperimentPreset.allCases.count))", icon: "flask")
                    ],
                    scale: .regular
                )
            }
            .padding(.horizontal, 20)
            .padding(.vertical, 14)
            .background(BlueHour.controlSurface)
            
            Divider()
            
            ScrollView {
                if selectedTab == 0 {
                    quickTourContent
                } else {
                    experimentsGalleryGrid
                }
            }
            .background(BlueHour.windowBackground)
            
            Divider()
            
            // Footer Bar
            HStack(spacing: 12) {
                Toggle("Show this Welcome Guide on startup", isOn: $showOnStartup)
                    .toggleStyle(.checkbox)
                    .font(.caption)
                
                Button(action: {
                    isPresented = false
                    onOpenDocs(.quickStart)
                }) {
                    Label("Open Full Physics & UI Documentation", systemImage: "book.fill")
                }
                .buttonStyle(.bordered)
                .controlSize(.small)
                
                Spacer()
                
                Button("Start Exploring Studio") {
                    isPresented = false
                }
                .buttonStyle(.borderedProminent)
                .controlSize(.regular)
            }
            .padding(.horizontal, 20)
            .padding(.vertical, 12)
            .background(BlueHour.controlSurface)
        }
        .frame(width: 860, height: 610)
    }
    
    private var quickTourContent: some View {
        VStack(alignment: .leading, spacing: 18) {
            // 3-Step Visual Workflow Cards
            Text("How Lattice Launcher Works in 3 Simple Steps")
                .font(.system(size: 14, weight: .bold))
                .foregroundColor(BlueHour.ice)
            
            HStack(alignment: .top, spacing: 12) {
                tourStepCard(
                    stepNumber: "STEP 1",
                    title: "Choose Crystal & Microstructure",
                    icon: "cube.transparent.fill",
                    accent: Color(red: 0.31, green: 0.85, blue: 0.78),
                    points: [
                        "Pick Square or Triangular lattice symmetry and system size (e.g. 100×100).",
                        "Choose a defect mode: Perfect Crystal, Bicrystal Grain Boundary (Σ5, Σ13…), 2D Voronoi Polycrystal, or 2M Periodic Dislocations.",
                        "Use 🌱 Guided Mode for a clean 3-step sidebar, or switch to ⚙️ Expert Mode anytime for all 35+ research controls."
                    ]
                )
                
                tourStepCard(
                    stepNumber: "STEP 2",
                    title: "Customize Geometry & BCs Visually",
                    icon: "hand.draw.fill",
                    accent: Color(red: 1.0, green: 0.78, blue: 0.35),
                    points: [
                        "Open the interactive 'Geometry & BCs' tab to see every discrete atom and boundary arrow.",
                        "Drag dislocation cores (⊥/⊤), grain boundary planes, Voronoi seeds, or circular/crack cutouts directly on the canvas.",
                        "Inspect both Stage 1 (Relaxation at α=0) and Stage 2 (Mechanical Loading / Nanoindentation) boundary conditions."
                    ]
                )
                
                tourStepCard(
                    stepNumber: "STEP 3",
                    title: "Run & Analyze Live Mechanics",
                    icon: "chart.xyaxis.line",
                    accent: Color(red: 0.45, green: 0.90, blue: 0.62),
                    points: [
                        "Click 'Launch Simulation' (⌘R) to run the C++ CHOLMOD-preconditioned L-BFGS solver.",
                        "Watch real-time Stress–Strain σ(α), Grain Boundary Energy γ_GB(Δθ), and plastic avalanche drops.",
                        "Click 'Open Output Folder' to inspect configuration_*.vtk snapshots in ParaView."
                    ]
                )
            }
            
            Divider().padding(.vertical, 2)
            
            HStack {
                Text("Try a Featured Experiment in One Click")
                    .font(.system(size: 14, weight: .bold))
                    .foregroundColor(BlueHour.ice)
                Spacer()
                Button("View All \(SampleExperimentPreset.allCases.count) Sample Experiments →") {
                    selectedTab = 1
                }
                .buttonStyle(.link)
                .font(.caption.bold())
            }
            
            LazyVGrid(columns: [GridItem(.flexible(), spacing: 12), GridItem(.flexible(), spacing: 12)], spacing: 12) {
                ForEach([
                    SampleExperimentPreset.homogeneousShearAvalanches,
                    SampleExperimentPreset.periodicDislocations6D,
                    SampleExperimentPreset.sigma5SymmetricTiltGB,
                    SampleExperimentPreset.nanoindentationCrystal
                ]) { preset in
                    experimentCard(preset)
                }
            }
        }
        .padding(20)
    }
    
    private var experimentsGalleryGrid: some View {
        VStack(alignment: .leading, spacing: 14) {
            Text("Curated Sample Experiments — Click 'Load & Preview' to configure all parameters automatically:")
                .font(.system(size: 12.5))
                .foregroundColor(BlueHour.mist)
            
            LazyVGrid(columns: [GridItem(.flexible(), spacing: 12), GridItem(.flexible(), spacing: 12)], spacing: 12) {
                ForEach(SampleExperimentPreset.allCases) { preset in
                    experimentCard(preset)
                }
            }
        }
        .padding(20)
    }
    
    private func tourStepCard(stepNumber: String, title: String, icon: String, accent: Color, points: [String]) -> some View {
        VStack(alignment: .leading, spacing: 8) {
            HStack {
                Text(stepNumber)
                    .font(.system(size: 9.5, weight: .heavy, design: .monospaced))
                    .padding(.horizontal, 6)
                    .padding(.vertical, 2)
                    .background(accent.opacity(0.22), in: Capsule())
                    .foregroundColor(accent)
                Spacer()
                Image(systemName: icon)
                    .font(.system(size: 15, weight: .bold))
                    .foregroundColor(accent)
            }
            
            Text(title)
                .font(.system(size: 13, weight: .bold))
                .foregroundColor(BlueHour.ice)
            
            ForEach(points, id: \.self) { pt in
                HStack(alignment: .top, spacing: 6) {
                    Circle()
                        .fill(accent)
                        .frame(width: 5, height: 5)
                        .padding(.top, 5)
                    Text(pt)
                        .font(.system(size: 11))
                        .foregroundColor(.secondary)
                        .fixedSize(horizontal: false, vertical: true)
                }
            }
            Spacer(minLength: 0)
        }
        .padding(12)
        .frame(maxWidth: .infinity, minHeight: 175, alignment: .topLeading)
        .background(BlueHour.midnight.opacity(0.65), in: RoundedRectangle(cornerRadius: 8))
        .overlay(RoundedRectangle(cornerRadius: 8).strokeBorder(accent.opacity(0.45), lineWidth: 1))
    }
    
    private func experimentCard(_ preset: SampleExperimentPreset) -> some View {
        VStack(alignment: .leading, spacing: 8) {
            HStack(spacing: 8) {
                Image(systemName: preset.icon)
                    .font(.system(size: 15, weight: .bold))
                    .foregroundColor(preset.difficultyColor)
                    .frame(width: 24, height: 24)
                    .background(preset.difficultyColor.opacity(0.15), in: RoundedRectangle(cornerRadius: 6))
                
                Text(preset.title)
                    .font(.system(size: 12.5, weight: .bold))
                    .foregroundColor(BlueHour.ice)
                    .lineLimit(1)
                
                Spacer()
                
                Text(preset.difficultyBadge)
                    .font(.system(size: 9.5, weight: .bold, design: .monospaced))
                    .padding(.horizontal, 6)
                    .padding(.vertical, 2)
                    .background(preset.difficultyColor.opacity(0.18), in: Capsule())
                    .foregroundColor(preset.difficultyColor)
            }
            
            Text(preset.physicsSummary)
                .font(.system(size: 11))
                .foregroundColor(.primary.opacity(0.88))
                .lineLimit(3)
                .fixedSize(horizontal: false, vertical: true)
            
            HStack(alignment: .top, spacing: 5) {
                Image(systemName: "eye.fill")
                    .font(.system(size: 9.5))
                    .foregroundColor(BlueHour.mist)
                    .padding(.top, 2)
                Text(preset.whatToObserve)
                    .font(.system(size: 10.5).italic())
                    .foregroundColor(BlueHour.mist)
                    .lineLimit(2)
            }
            
            Spacer(minLength: 2)
            
            HStack(spacing: 8) {
                Button(action: {
                    onSelectExperiment(preset, false)
                    isPresented = false
                }) {
                    Label("Load & Inspect Setup", systemImage: "slider.horizontal.3")
                        .font(.system(size: 11, weight: .semibold))
                        .frame(maxWidth: .infinity)
                }
                .buttonStyle(.bordered)
                .controlSize(.small)
                
                Button(action: {
                    onSelectExperiment(preset, true)
                    isPresented = false
                }) {
                    Label("Load & Run", systemImage: "play.fill")
                        .font(.system(size: 11, weight: .bold))
                        .frame(maxWidth: .infinity)
                }
                .buttonStyle(.borderedProminent)
                .tint( Color(red: 0.12, green: 0.64, blue: 0.54) )
                .controlSize(.small)
            }
        }
        .padding(12)
        .frame(maxWidth: .infinity, minHeight: 158, alignment: .topLeading)
        .background(BlueHour.midnight.opacity(0.55), in: RoundedRectangle(cornerRadius: 8))
        .overlay(RoundedRectangle(cornerRadius: 8).strokeBorder(BlueHour.border, lineWidth: 1))
    }
}

// MARK: - Built-in Documentation & Physics Reference Center
public struct DocumentationCenterSheet: View {
    @Binding var isPresented: Bool
    @Binding var selectedTopic: DocumentationTopic
    @State private var searchText: String = ""
    
    public var body: some View {
        VStack(spacing: 0) {
            // Header
            HStack(spacing: 12) {
                Image(systemName: "book.closed.fill")
                    .font(.system(size: 18, weight: .bold))
                    .foregroundColor(Color(red: 0.31, green: 0.85, blue: 0.78))
                VStack(alignment: .leading, spacing: 1) {
                    Text("Lattice Launcher Documentation & Physics Reference")
                        .font(.system(size: 16, weight: .bold))
                        .foregroundColor(BlueHour.ice)
                    Text("Interactive guide to MTM crystal mechanics, dislocations, grain boundaries & boundary conditions")
                        .font(.system(size: 11))
                        .foregroundColor(BlueHour.mist)
                }
                Spacer()
                
                Button("Done") {
                    isPresented = false
                }
                .buttonStyle(.borderedProminent)
                .controlSize(.small)
            }
            .padding(.horizontal, 18)
            .padding(.vertical, 12)
            .background(BlueHour.controlSurface)
            
            Divider()
            
            HSplitView {
                // Left Chapter List
                VStack(alignment: .leading, spacing: 4) {
                    Text("CHAPTERS")
                        .font(.system(size: 10, weight: .bold, design: .monospaced))
                        .foregroundColor(BlueHour.mist)
                        .padding(.horizontal, 10)
                        .padding(.top, 10)
                    
                    ForEach(DocumentationTopic.allCases) { topic in
                        let isSel = (selectedTopic == topic)
                        Button(action: { selectedTopic = topic }) {
                            HStack(spacing: 8) {
                                Image(systemName: topic.icon)
                                    .font(.system(size: 12, weight: .semibold))
                                    .frame(width: 18)
                                Text(topic.title)
                                    .font(.system(size: 11.5, weight: isSel ? .bold : .medium))
                                    .multilineTextAlignment(.leading)
                                Spacer()
                            }
                            .foregroundColor(isSel ? .white : BlueHour.ice.opacity(0.88))
                            .padding(.horizontal, 10)
                            .padding(.vertical, 7)
                            .background(
                                isSel ? Color(red: 0.12, green: 0.64, blue: 0.54) : Color.clear,
                                in: RoundedRectangle(cornerRadius: 6)
                            )
                        }
                        .buttonStyle(.plain)
                    }
                    
                    Divider().padding(.vertical, 6)
                    
                    Text("TECHNICAL PDF / MD GUIDES")
                        .font(.system(size: 10, weight: .bold, design: .monospaced))
                        .foregroundColor(BlueHour.mist)
                        .padding(.horizontal, 10)
                    
                    docFileButton(title: "How to Install Dislocations (.md)", fileName: "how_to_install_dislocations.md")
                    docFileButton(title: "Special GBs & BCs Math (.pdf)", fileName: "special_grain_boundaries_and_bcs_math.pdf")
                    docFileButton(title: "Discrete Burgers Circuit (.pdf)", fileName: "discrete_burgers_circuit_method.pdf")
                    docFileButton(title: "Plasticity & Remeshing (.pdf)", fileName: "plasticity_memory_remeshing.pdf")
                    docFileButton(title: "Aligned Deformation Gradient (.pdf)", fileName: "aligned_deformation_gradient.pdf")
                    docFileButton(title: "Stress & Traction Control (.pdf)", fileName: "stress_and_traction_control.pdf")
                    docFileButton(title: "Periodic Rotated Lattices (.pdf)", fileName: "periodic_rotated_lattices_recipe.pdf")
                    docFileButton(title: "Preconditioned L-BFGS (.pdf)", fileName: "preconditioned_lbfgs.pdf")
                    
                    Spacer()
                }
                .padding(8)
                .frame(minWidth: 235, maxWidth: 255)
                .background(BlueHour.sidebarBackground)
                
                // Right Chapter Content
                ScrollView {
                    VStack(alignment: .leading, spacing: 14) {
                        topicDetailView(selectedTopic)
                    }
                    .padding(22)
                    .frame(maxWidth: .infinity, alignment: .leading)
                }
                .background(BlueHour.windowBackground)
            }
        }
        .frame(width: 880, height: 600)
    }
    
    private func docFileButton(title: String, fileName: String) -> some View {
        Button(action: { openDocumentationFile(fileName) }) {
            HStack(spacing: 6) {
                Image(systemName: fileName.hasSuffix(".pdf") ? "doc.richtext.fill" : "doc.text.fill")
                    .font(.system(size: 11))
                    .foregroundColor(Color(red: 1.0, green: 0.82, blue: 0.35))
                Text(title)
                    .font(.system(size: 11))
                    .foregroundColor(BlueHour.ice)
                    .lineLimit(1)
                Spacer()
                Image(systemName: "arrow.up.right.square")
                    .font(.system(size: 10))
                    .foregroundColor(BlueHour.mist)
            }
            .padding(.horizontal, 10)
            .padding(.vertical, 5)
        }
        .buttonStyle(.plain)
    }
    
    private func openDocumentationFile(_ fileName: String) {
        let root = SimulationConfig.detectProjectRoot()
        let candidates = [
            Bundle.main.resourceURL?.appendingPathComponent("Docs/\(fileName)").path,
            Bundle.main.resourceURL?.appendingPathComponent(fileName).path,
            "\(root)/docs/\(fileName)"
        ].compactMap { $0 }
        
        for path in candidates {
            if FileManager.default.fileExists(atPath: path) {
                NSWorkspace.shared.open(URL(fileURLWithPath: path))
                return
            }
        }
    }
    
    @ViewBuilder
    private func topicDetailView(_ topic: DocumentationTopic) -> some View {
        switch topic {
        case .quickStart:
            docSectionHeader("Quick Start & 3-Step Workflow", subtitle: "How to set up, preview, and run your first MTM crystal simulation")
            docCard(
                title: "1. 🌱 Guided Mode vs. ⚙️ Expert Mode",
                lines: [
                    "• 🌱 Guided Mode (default) organizes the left sidebar into 3 clear steps: (1) Crystal & Defect Setup, (2) Macroscopic Geometry & Boundary Studio, and (3) Loading Protocol & Schedule, while placing low-level solver settings inside an expandable '4. Advanced Solver & Hardware' drawer.",
                    "• ⚙️ Expert Mode expands every solver card immediately so all 35+ research parameters are visible at a glance.",
                    "• You can switch between Guided and Expert Mode anytime from the top window bar or View > Toggle Guided / Expert Mode (⌥⌘E)."
                ]
            )
            docCard(
                title: "2. Main Workspace Tabs (⌘1 – ⌘4)",
                lines: [
                    "• 📈 Curve (⌘1): Plots live Stress–Strain σ(α), Grain Boundary Energy γ_GB(Δθ), Lattice Rotation ω(α), or radial Volterra Dislocation Profiles.",
                    "• ω(α) is the rotation of the crystal lattice measured from the current atom positions (per grain and mean ± spread; optionally net of the applied rigid rotation). Per-element fields ElementRotationAngleDeg, ElementLatticeOrientationDeg and ElementLatticeDistortion are in every VTK snapshot; cells distorted by more than 0.35 shear (dislocation cores, transformed regions) are excluded from the statistics.",
                    "• 📐 Geometry & BCs (⌘2): Interactive 2D studio showing discrete lattice atoms, grain boundaries, dislocation glide planes, cutouts, nanoindenter, and boundary condition vectors.",
                    "• 🕸 Initial Mesh (⌘3): Generates and renders the actual Delaunay triangulation and per-element elastic energy prior to running.",
                    "• 💻 Terminal Log (⌘4): Streams real-time solver iterations, CHOLMOD factorization metrics, and avalanche timestamps."
                ]
            )
            
        case .installingDislocations:
            docSectionHeader("Installing Edge Dislocations", subtitle: "Mode A: Multi-Plane Periodic Torus (2M Planes) vs. Mode B: Single Volterra Core")
            docCard(
                title: "Mode A: Multi-Plane Single-Dislocation Periodic Torus (Periodic Dislocs)",
                lines: [
                    "• On a 2D periodic torus [0, Lx) × [0, Ly), Stokes' theorem requires the net Burgers vector sum to vanish: ∑ b_k = 0.",
                    "• To install isolated single dislocations on glide planes under full 2D PBCs, choose M balanced (±a₀) pairs across 2M distinct horizontal glide planes — strictly 1 dislocation per active glide plane.",
                    "• Across glide plane y_k = (j_k + 0.5)h_y, the integer number of atoms per row jumps by N_{j_k+1} - N_{j_k} = s_k ∈ {+1, -1}.",
                    "• Vertical ±b sign order does NOT need to be alternating (+ - + -): use '🎲 Shuffle ±b' for random sign order, '+ − Alternating' for minimal strip pre-strain, or '+ + − − Grouped' for a dipolar tilt band.",
                    "• For large ensembles (e.g. M=30 pairs = 60 dislocations), use a box with NY ≥ 200 (such as 240×240 or 300×300) so adjacent glide planes are separated by at least 4–5 atomic rows."
                ]
            )
            docCard(
                title: "Mode B: Isolated Single Edge Dislocation (Volterra / Stroh)",
                lines: [
                    "• Select Loading Mode = 'Edge Dislocation (Volterra / Stroh)' to place 1 dislocation at the center (0,0) of a finite circular cylinder of radius R_free.",
                    "• Atoms outside r > R_free are locked to the analytical anisotropic Stroh (Square) or isotropic Volterra (Triangular) displacement field, while interior atoms r ≤ R_free relax freely."
                ]
            )
            
        case .grainBoundariesCSL:
            docSectionHeader("Bicrystals, 2D Voronoi & CSL Rationalizer", subtitle: "Constructing commensurate symmetric tilt grain boundaries and polycrystals")
            docCard(
                title: "Coincidence Site Lattice (CSL) Rationalizer (gb_code / aimsgb)",
                lines: [
                    "• In 'Bicrystal (GB)' mode, a planar interface divides Grain 0 (θ₁) and Grain 1 (θ₂) at inclination angle ϕ_GB.",
                    "• For a Symmetric Tilt Grain Boundary (STGB) generated by coprime integers (m, n), the exact half-tilt angle is θ = arctan(n/m) with Σ = (m² + n²) / gcd(m² + n², 2).",
                    "• Clicking 'Snap to Exact (m,n) CSL + Box & BCs' automatically sets θ₁, θ₂, the commensurate box size K·L_CSL, the seam overlap cutoff d_min*, and the normal expansion gap δ_⊥*."
                ]
            )
            docCard(
                title: "Atomic Seam Overlap Resolution",
                lines: [
                    "• Where two misoriented lattices meet at the GB plane, pairs of cross-grain atoms closer than d_min·a₀ are resolved via 'Merge Pair to Midpoint' (recommended) or 'Delete Grain 2 Atom'.",
                    "• Merged seam-core atoms are highlighted with a gold ring on the Geometry Studio canvas."
                ]
            )
            docCard(
                title: "2D Voronoi Polycrystals (Finite or Periodic)",
                lines: [
                    "• Each grain is a lattice rotated by θ_k and anchored at its seed; an atom belongs to the grain with the nearest seed (physical distance). Seeds come from a stratified layout driven by Grains N_g and Seed, or are dragged on the canvas.",
                    "• Min. Neighbour Misorientation redraws θ_k when it lies within that angle of an adjacent grain (modulo 60° triangular / 90° square), avoiding accidental near-identical neighbours.",
                    "• Finite specimens relax in Stage 1 under a minimal anchor by default (Bottom RollerX + Left RollerY, top and right free), so the grain network relaxes without the loading grips clamping it.",
                    "• 'Periodic (T² torus)' tiles a fully periodic box: no free surfaces or grips, every atom free in Stage 1, the periodic box sheared in Stage 2. Grains that cross the box edge stay continuous, so the only seams are real grain boundaries.",
                    "• On the torus a grain must never border its own periodic copy (its two rotated copies would not line up). That happens with 2–3 grains or aligned seeds; the canvas draws such edges in red and the solver refuses the run."
                ]
            )
            
        case .boundaryConditionsIndenter:
            docSectionHeader("Two-Stage Boundary Conditions & Nanoindenter", subtitle: "Decoupling Stage 1 equilibrium relaxation (α = 0) from Stage 2 mechanical loading (α > 0)")
            docCard(
                title: "Why Two Stages? (Stage 1 Relax → Stage 2 Load)",
                lines: [
                    "• Stage 1 — Relaxation Step (α = 0): Equilibrates grain boundary cores, cutouts, or periodic dislocations at zero external load. For bicrystals, roller grips allow Poisson lateral breathing without artificial clamping stress.",
                    "• Stage 2 — Loading Step (α > 0): Re-anchors the reference frame (X₀ ← x_eq) and applies shear grips, uniaxial tension/compression, or nanoindentation.",
                    "• In Geometry & BCs Studio, use the top-left canvas badge button ('← Relaxation Step BCs' / 'Loading Step BCs →') to inspect and edit both stages."
                ]
            )
            
        case .acousticAndStiffness:
            docSectionHeader("Acoustic Tensor & Stiffness Spectrum K(α)", subtitle: "Continuum rank-1 convexity vs. discrete Hessian eigenvalue stability")
            docCard(
                title: "Parametric Acoustic Study & Stiffness Spectrum Tabs",
                lines: [
                    "• Parametric Acoustic Study (Tab 2): Evaluates the continuum acoustic tensor Q(n) = min_n det A(F, n) along homogeneous deformation paths to map rank-1 ellipticity loss.",
                    "• Stiffness Spectrum K(α) (Tab 3): Builds the discrete MTM tangent stiffness matrix K(α) on an Nx×Ny mesh and computes the lowest eigenmodes (λ_min → 0) via shift-invert Lanczos."
                ]
            )
            
        case .parameterGlossary:
            docSectionHeader("Glossary of Physical Symbols & Parameters", subtitle: "Quick reference for all symbols used across the interface")
            docCard(
                title: "Crystal & Potential Parameters",
                lines: [
                    "• a₀ : Equilibrium nearest-neighbor lattice spacing of the active potential.",
                    "• β (Potential Beta) : Shear anisotropy / barrier height parameter in the Zanzotto continuum potential Φ(C) (standard Square: β = -0.25; Triangular: β = 4.0).",
                    "• K (Bulk Modulus) : Volumetric compressibility modulus in Φ(C) (standard: K = 5.0).",
                    "• σ (Initial Noise) : Standard deviation of Gaussian random displacement added to initial positions (use σ = 0.0 for clean GBs and periodic dislocations; σ = 0.04 for homogeneous shear)."
                ]
            )
            docCard(
                title: "Microstructure, Dislocation & Solver Symbols",
                lines: [
                    "• M (Dislocation Pairs) : Number of balanced (+a₀, -a₀) edge dislocation dipoles across 2M distinct glide planes.",
                    "• r₀ (Core Radius) : Cai non-singular core regularization radius (in units of a₀, default 0.8 a₀).",
                    "• d_min : Cross-grain seam merge cutoff distance (in units of a₀, typically 0.55–0.70 a₀).",
                    "• δ_⊥* : Normal rigid-body expansion across a grain boundary (in units of a₀).",
                    "• α : Loading parameter (engineering shear strain γ, axial strain ε, or normalized indenter depth δ/δ_max).",
                    "• λ_min : Smallest eigenvalue of the free-DOF Hessian stiffness matrix K; λ_min ≤ 0 signals plastic instability (avalanche)."
                ]
            )
        }
    }
    
    private func docSectionHeader(_ title: String, subtitle: String) -> some View {
        VStack(alignment: .leading, spacing: 3) {
            Text(title)
                .font(.system(size: 18, weight: .bold))
                .foregroundColor(BlueHour.ice)
            Text(subtitle)
                .font(.system(size: 12))
                .foregroundColor(BlueHour.mist)
        }
    }
    
    private func docCard(title: String, lines: [String]) -> some View {
        VStack(alignment: .leading, spacing: 7) {
            Text(title)
                .font(.system(size: 13, weight: .bold))
                .foregroundColor(Color(red: 0.31, green: 0.85, blue: 0.78))
            ForEach(lines, id: \.self) { line in
                Text(line)
                    .font(.system(size: 12))
                    .foregroundColor(.primary.opacity(0.92))
                    .fixedSize(horizontal: false, vertical: true)
            }
        }
        .padding(14)
        .frame(maxWidth: .infinity, alignment: .leading)
        .background(BlueHour.midnight.opacity(0.60), in: RoundedRectangle(cornerRadius: 8))
        .overlay(RoundedRectangle(cornerRadius: 8).strokeBorder(BlueHour.border, lineWidth: 1))
    }
}
