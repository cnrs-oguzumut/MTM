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

/// Navigation intents used by the sidebar ("show me the blueprint / the mesh / the results").
/// The actual on-screen layout is driven by `WorkspaceMode` + `SpecimenSubView`.
enum MainDisplayTab: String, CaseIterable, Identifiable {
    case geometryStudio = "blueprint"
    case initialMesh = "femMesh"
    case curve = "results"
    
    var id: String { rawValue }
}

/// The three top-level workspaces of the main viewport.
enum WorkspaceMode: String, CaseIterable, Identifiable {
    case specimen = "Specimen & Mesh"
    case results = "Results & Plots"
    case split = "Split View"
    
    var id: String { rawValue }
    var icon: String {
        switch self {
        case .specimen: return "square.on.square.dashed"
        case .results: return "chart.xyaxis.line"
        case .split: return "rectangle.split.2x1"
        }
    }
}

/// Two views of the same initial specimen.
enum SpecimenSubView: String, CaseIterable, Identifiable {
    case blueprint = "Blueprint"
    case femMesh = "FEM Mesh"
    
    var id: String { rawValue }
    var icon: String {
        switch self {
        case .blueprint: return "hand.draw"
        case .femMesh: return "triangle"
        }
    }
}

/// Results workspace: live curves or the built-in VTK snapshot viewer.
enum ResultsSubView: String { case curves, snapshots }

/// The four steps of the unified left setup inspector.
enum SetupInspectorStep: String, CaseIterable, Identifiable {
    case crystal
    case shapeAndBCs
    case protocolAndSchedule
    case solver
    
    var id: String { rawValue }
    var shortTitle: String {
        switch self {
        case .crystal: return "1 Crystal"
        case .shapeAndBCs: return "2 Shape & BCs"
        case .protocolAndSchedule: return "3 Protocol"
        case .solver: return "4 Solver"
        }
    }
    var icon: String {
        switch self {
        case .crystal: return "atom"
        case .shapeAndBCs: return "skew"
        case .protocolAndSchedule: return "slider.horizontal.below.square.filled.and.square"
        case .solver: return "gearshape.2"
        }
    }
    var helpText: String {
        switch self {
        case .crystal: return "Size, lattice, potential, grains & dislocations"
        case .shapeAndBCs: return "Specimen shape, holes, boundary conditions & indenter"
        case .protocolAndSchedule: return "Loading mode, noise, remeshing & strain schedule"
        case .solver: return "Preconditioner, eigenvalues, output & threads"
        }
    }
}

struct ContentView: View {
    @ObservedObject var config: SimulationConfig
    @EnvironmentObject var runner: ProcessRunner
    @StateObject private var previewManager = StiffnessStudyManager()
    @State private var appTab: AppNavigationTab = .simulation
    @State private var filterText: String = ""
    @State private var autoScroll: Bool = true
    
    // Layout state (one sidebar, one viewport, one console drawer)
    @State private var workspace: WorkspaceMode = .specimen
    @State private var specimenSub: SpecimenSubView = .blueprint
    @State private var resultsSub: ResultsSubView = .curves
    @State private var inspectorStep: SetupInspectorStep = .crystal
    @AppStorage("showConsoleDrawer") private var showConsole: Bool = false
    
    // Onboarding & Help
    @AppStorage("uiComplexityMode") private var uiComplexityModeRaw: String = UIComplexityMode.guided.rawValue
    @AppStorage("showWelcomeOnStartup_v2") private var showWelcomeOnStartup: Bool = true
    @State private var didCheckStartupWelcome: Bool = false
    @State private var showWelcomeSheet: Bool = false
    @State private var showDocsSheet: Bool = false
    @State private var selectedDocTopic: DocumentationTopic = .quickStart
    
    private var uiComplexityMode: UIComplexityMode {
        UIComplexityMode(rawValue: uiComplexityModeRaw) ?? .guided
    }
    
    // MARK: Navigation
    
    /// Sidebar-facing binding: maps "show blueprint / mesh / results" intents onto the workspace layout.
    private var displayTabBinding: Binding<MainDisplayTab> {
        Binding(
            get: {
                if workspace == .results { return .curve }
                return specimenSub == .femMesh ? .initialMesh : .geometryStudio
            },
            set: { navigate(to: $0) }
        )
    }
    
    private func navigate(to tab: MainDisplayTab) {
        appTab = .simulation
        switch tab {
        case .geometryStudio:
            specimenSub = .blueprint
            if workspace == .results { workspace = .specimen }
        case .initialMesh:
            specimenSub = .femMesh
            if workspace == .results { workspace = .specimen }
        case .curve:
            if workspace == .specimen { workspace = .results }
        }
    }
    
    private func applyPreset(_ preset: SampleExperimentPreset, runNow: Bool = false) {
        appTab = .simulation
        let openStudio = preset.apply(to: config)
        inspectorStep = .crystal
        navigate(to: openStudio ? .geometryStudio : .curve)
        if runNow && !runner.status.isRunning {
            navigate(to: .curve)
            launchRun()
        }
    }
    
    private func launchRun() {
        guard !runner.status.isRunning else { return }
        if config.mode != .volterra && config.mode != .shift && config.protocolMode == .sweepMisorientation {
            runner.startCSLSweep(config: config, candidates: config.sweepCSLCandidates)
        } else {
            runner.start(config: config)
        }
    }
    
    /// Syncs the preview manager with the current configuration and regenerates the FEM mesh.
    private func refreshMeshPreview() {
        previewManager.nx = config.nx
        previewManager.ny = config.ny
        previewManager.latticeSymmetry = config.latticeSymmetry
        previewManager.potentialMode = config.potentialMode
        previewManager.boundaryCondition = (config.mode == .shift || config.geometryConfig.enabled) ? .fixed : config.boundaryCondition
        previewManager.beta = config.beta
        previewManager.bulkK = config.bulkK
        previewManager.useRotation = config.useRotation
        previewManager.rotationAngle = config.rotationAngle
        previewManager.cslM = config.cslM
        previewManager.cslN = config.cslN
        previewManager.cslKx = config.cslKx
        previewManager.cslKy = config.cslKy
        if config.geometryConfig.enabled {
            let root = SimulationConfig.detectProjectRoot()
            let scratchDir = "\(root)/scratch"
            try? FileManager.default.createDirectory(atPath: scratchDir, withIntermediateDirectories: true)
            let geomPath = "\(scratchDir)/geometry_config.json"
            try? config.geometryConfig.writeJSON(to: geomPath)
            previewManager.geometryJSONPath = geomPath
        } else {
            previewManager.geometryJSONPath = ""
        }
        previewManager.generateInitialMeshPreview()
    }
    
    // MARK: Top Bar
    
    private var globalTopHeaderBar: some View {
        HStack(spacing: 12) {
            // Left: module selector (compact menu instead of 3 big tabs)
            Menu {
                ForEach(AppNavigationTab.allCases) { tab in
                    Button(action: { appTab = tab }) {
                        Label(tab.rawValue, systemImage: tab.icon)
                    }
                }
            } label: {
                Label(appTab.rawValue, systemImage: appTab.icon)
                    .font(.system(size: 12, weight: .semibold))
            }
            .menuStyle(.borderlessButton)
            .fixedSize()
            .help("Switch tool: Lattice Simulation, Parametric Acoustic Study or Stiffness Spectrum K(α)")
            
            Spacer(minLength: 8)
            
            // Center: the ONE viewport switcher
            if appTab == .simulation {
                Picker("Workspace", selection: $workspace) {
                    ForEach(WorkspaceMode.allCases) { mode in
                        Label(mode.rawValue, systemImage: mode.icon).tag(mode)
                    }
                }
                .labelsHidden()
                .pickerStyle(.segmented)
                .fixedSize()
            }
            
            Spacer(minLength: 8)
            
            // Right: examples & help
            Menu {
                Section("1-Click Sample Experiments") {
                    ForEach(SampleExperimentPreset.allCases) { preset in
                        Button(action: { applyPreset(preset) }) {
                            Label(preset.title, systemImage: preset.icon)
                        }
                    }
                }
                Divider()
                Button(action: { showWelcomeSheet = true }) {
                    Label("Browse Visual Experiment Gallery…", systemImage: "sparkles.rectangle.stack")
                }
            } label: {
                Label("Examples", systemImage: "flask")
                    .font(.system(size: 11, weight: .semibold))
            }
            .menuStyle(.borderedButton)
            .controlSize(.small)
            .fixedSize()
            
            Button(action: { showWelcomeSheet = true }) {
                Label("Tour", systemImage: "sparkles")
            }
            .buttonStyle(.bordered)
            .controlSize(.small)
            .help("Welcome tour & sample gallery")
            
            Button(action: { showDocsSheet = true }) {
                Label("Help", systemImage: "questionmark.circle")
            }
            .buttonStyle(.bordered)
            .controlSize(.small)
            .help("Physics & user guide (⌘?)")
        }
        .padding(.horizontal, 14)
        .frame(height: 44)
        .background(BlueHour.windowBackground)
    }
    
    // MARK: Viewport Panes
    
    private var specimenPane: some View {
        VStack(spacing: 0) {
            HStack(spacing: 10) {
                Text("Specimen")
                    .font(.system(size: 11.5, weight: .bold))
                    .foregroundColor(BlueHour.ice)
                Picker("View", selection: $specimenSub) {
                    ForEach(SpecimenSubView.allCases) { sub in
                        Label(sub.rawValue, systemImage: sub.icon).tag(sub)
                    }
                }
                .labelsHidden()
                .pickerStyle(.segmented)
                .controlSize(.small)
                .fixedSize()
                
                if specimenSub == .femMesh {
                    Button(action: refreshMeshPreview) {
                        if previewManager.isGeneratingPreview {
                            ProgressView().scaleEffect(0.5).frame(width: 14, height: 14)
                        } else {
                            Label("Regenerate Mesh", systemImage: "arrow.clockwise")
                        }
                    }
                    .controlSize(.small)
                    .help("Re-mesh using the current crystal, shape and BC settings")
                }
                Spacer()
                Text("\(config.nx)×\(config.ny) \(config.latticeSymmetry == .square ? "square" : "triangular") lattice")
                    .font(.system(size: 10, design: .monospaced))
                    .foregroundColor(BlueHour.mist)
            }
            .padding(.horizontal, 12)
            .frame(height: 34)
            .background(BlueHour.controlSurface)
            
            Divider()
            
            if specimenSub == .blueprint {
                GeometryBuilderView(
                    config: config,
                    previewManager: previewManager,
                    displayMode: .canvasOnly,
                    onSwitchToCarvedMesh: { specimenSub = .femMesh },
                    onOpenDocs: { topic in
                        selectedDocTopic = topic
                        showDocsSheet = true
                    }
                )
            } else {
                InitialLatticeVisualizerView(manager: previewManager)
            }
        }
        .onChange(of: specimenSub) { sub in
            if sub == .femMesh { refreshMeshPreview() }
        }
    }
    
    private var resultsPane: some View {
        VStack(spacing: 0) {
            HStack(spacing: 10) {
                Text("Results")
                    .font(.system(size: 11.5, weight: .bold))
                    .foregroundColor(BlueHour.ice)
                Picker("Results view", selection: $resultsSub) {
                    Label("Curves", systemImage: "chart.xyaxis.line").tag(ResultsSubView.curves)
                    Label("Snapshots", systemImage: "cube.transparent").tag(ResultsSubView.snapshots)
                }
                .labelsHidden()
                .pickerStyle(.segmented)
                .controlSize(.small)
                .fixedSize()
                Divider().frame(height: 18)
                LiveStatusHeaderView(config: config)
            }
            .padding(.horizontal, 12)
            .frame(height: 44)
            .background(BlueHour.controlSurface)
            Divider()
            if resultsSub == .snapshots {
                VTKSnapshotViewer(runDirectory: config.targetDirectory, isRunning: runner.status.isRunning)
            } else if config.mode == .volterra {
                DislocationProfileChartView()
            } else {
                StressStrainChartView()
            }
        }
    }
    
    @ViewBuilder
    private var viewport: some View {
        switch workspace {
        case .specimen:
            specimenPane
        case .results:
            resultsPane
        case .split:
            HSplitView {
                specimenPane.frame(minWidth: 360)
                resultsPane.frame(minWidth: 360)
            }
        }
    }
    
    private var consoleDrawer: some View {
        VStack(spacing: 0) {
            HStack(spacing: 8) {
                Image(systemName: "terminal")
                    .foregroundColor(BlueHour.ice)
                Text("Console")
                    .font(.system(size: 11, weight: .bold))
                    .foregroundColor(BlueHour.ice)
                Spacer()
                Button(action: { runner.clearLogs() }) {
                    Label("Clear", systemImage: "trash")
                }
                .controlSize(.small)
                Button(action: { withAnimation(.easeInOut(duration: 0.15)) { showConsole = false } }) {
                    Image(systemName: "chevron.down")
                }
                .buttonStyle(.borderless)
                .help("Hide console (⌘J)")
            }
            .padding(.horizontal, 12)
            .frame(height: 28)
            .background(BlueHour.controlSurface)
            Divider()
            TerminalView(filterText: $filterText, autoScroll: $autoScroll)
        }
        .frame(minHeight: 140, idealHeight: 220, maxHeight: 420)
    }
    
    // MARK: Body
    
    var body: some View {
        VStack(spacing: 0) {
            globalTopHeaderBar
            Divider()
            
            if appTab == .simulation {
                HSplitView {
                    SidebarControlsView(
                        config: config,
                        previewManager: previewManager,
                        displayTab: displayTabBinding,
                        inspectorStep: $inspectorStep,
                        uiComplexityMode: uiComplexityMode,
                        onOpenWelcomeGallery: { showWelcomeSheet = true },
                        onOpenDocs: { topic in
                            selectedDocTopic = topic
                            showDocsSheet = true
                        }
                    )
                    .frame(minWidth: 360, idealWidth: 410, maxWidth: 480)
                    
                    VStack(spacing: 0) {
                        if showConsole {
                            VSplitView {
                                viewport.frame(minHeight: 260)
                                consoleDrawer
                            }
                        } else {
                            viewport
                        }
                        Divider()
                        BottomActionBar(config: config, showConsole: $showConsole, onLaunch: launchRun)
                            .padding(.horizontal, 14)
                            .padding(.vertical, 8)
                            .background(BlueHour.windowBackground)
                    }
                    .frame(minWidth: 580)
                }
            } else if appTab == .acousticStudy {
                AcousticStudyView()
            } else if appTab == .stiffnessStudy {
                StiffnessStudyView()
            }
        }
        .onAppear {
            if !didCheckStartupWelcome {
                didCheckStartupWelcome = true
                if showWelcomeOnStartup {
                    DispatchQueue.main.asyncAfter(deadline: .now() + 0.35) {
                        showWelcomeSheet = true
                    }
                }
            }
        }
        // When a run starts, bring the results into view.
        .onChange(of: runner.status.isRunning) { running in
            if running && workspace == .specimen { workspace = .results }
        }
        // Editing the crystal or the shape? Show the specimen.
        .onChange(of: inspectorStep) { step in
            if (step == .crystal || step == .shapeAndBCs) && workspace == .results && !runner.status.isRunning {
                workspace = .specimen
            }
        }
        .sheet(isPresented: $showWelcomeSheet) {
            WelcomeOnboardingSheet(
                config: config,
                isPresented: $showWelcomeSheet,
                onSelectExperiment: { preset, runImmediately in
                    applyPreset(preset, runNow: runImmediately)
                },
                onOpenDocs: { topic in
                    showWelcomeSheet = false
                    selectedDocTopic = topic
                    DispatchQueue.main.asyncAfter(deadline: .now() + 0.2) {
                        showDocsSheet = true
                    }
                }
            )
        }
        .sheet(isPresented: $showDocsSheet) {
            DocumentationCenterSheet(isPresented: $showDocsSheet, selectedTopic: $selectedDocTopic)
        }
        .onReceive(NotificationCenter.default.publisher(for: .showWelcomeTour)) { _ in
            showWelcomeSheet = true
        }
        .onReceive(NotificationCenter.default.publisher(for: .showDocumentationCenter)) { note in
            if let raw = note.object as? String, let topic = DocumentationTopic(rawValue: raw) {
                selectedDocTopic = topic
            }
            showDocsSheet = true
        }
        .onReceive(NotificationCenter.default.publisher(for: .applySampleExperiment)) { note in
            if let raw = note.object as? String, let preset = SampleExperimentPreset(rawValue: raw) {
                applyPreset(preset)
            }
        }
        .onReceive(NotificationCenter.default.publisher(for: .switchMainDisplayTab)) { note in
            guard let raw = note.object as? String else { return }
            appTab = .simulation
            switch raw {
            case "blueprint": workspace = .specimen; specimenSub = .blueprint
            case "femMesh": workspace = .specimen; specimenSub = .femMesh
            case "results": workspace = .results; resultsSub = .curves
            case "split": workspace = .split
            case "snapshots":
                if workspace == .specimen { workspace = .results }
                resultsSub = .snapshots
            default: break
            }
        }
        .onReceive(NotificationCenter.default.publisher(for: .openVTKFile)) { _ in
            appTab = .simulation
            if workspace == .specimen { workspace = .results }
            resultsSub = .snapshots
            // Give the viewer a moment to appear, then ask it to show the open panel
            DispatchQueue.main.asyncAfter(deadline: .now() + 0.25) {
                NotificationCenter.default.post(name: .openVTKFilePanel, object: nil)
            }
        }
        .onReceive(NotificationCenter.default.publisher(for: .toggleConsoleDrawer)) { _ in
            withAnimation(.easeInOut(duration: 0.15)) { showConsole.toggle() }
        }
        .onReceive(NotificationCenter.default.publisher(for: .toggleComplexityMode)) { _ in
            uiComplexityModeRaw = (uiComplexityMode == .guided)
                ? UIComplexityMode.expert.rawValue
                : UIComplexityMode.guided.rawValue
        }
        .onReceive(NotificationCenter.default.publisher(for: .runSimulationCommand)) { _ in
            appTab = .simulation
            launchRun()
        }
        .onReceive(NotificationCenter.default.publisher(for: .stopSimulationCommand)) { _ in
            if runner.status.isRunning {
                runner.stop()
            }
        }
    }
}

// MARK: - Sidebar Controls
struct SidebarControlsView: View {
    @ObservedObject var config: SimulationConfig
    @ObservedObject var previewManager: StiffnessStudyManager
    @Binding var displayTab: MainDisplayTab
    @Binding var inspectorStep: SetupInspectorStep
    var uiComplexityMode: UIComplexityMode = .guided
    var onOpenWelcomeGallery: () -> Void = {}
    var onOpenDocs: (DocumentationTopic) -> Void = { _ in }
    @EnvironmentObject var runner: ProcessRunner
    
    @ViewBuilder
    private var quickStartExperimentBanner: some View {
        VStack(alignment: .leading, spacing: 8) {
            HStack(spacing: 8) {
                Image(systemName: "flask.fill")
                    .font(.system(size: 13, weight: .bold))
                    .foregroundColor(Color(red: 0.31, green: 0.85, blue: 0.78))
                VStack(alignment: .leading, spacing: 1) {
                    Text("Quick-Start Sample Experiments")
                        .font(.system(size: 12, weight: .bold))
                        .foregroundColor(BlueHour.ice)
                    Text("New to Lattice Launcher? Load a ready-to-run crystal or defect setup.")
                        .font(.system(size: 10))
                        .foregroundColor(BlueHour.mist)
                }
                Spacer()
                Button(action: onOpenWelcomeGallery) {
                    Text("Gallery…")
                        .font(.system(size: 10, weight: .bold))
                }
                .buttonStyle(.borderedProminent)
                .tint(Color(red: 0.14, green: 0.62, blue: 0.56))
                .controlSize(.mini)
            }
            
            Menu {
                ForEach(SampleExperimentPreset.allCases) { preset in
                    Button(action: {
                        let openStudio = preset.apply(to: config)
                        displayTab = openStudio ? .geometryStudio : .curve
                    }) {
                        Label(preset.title, systemImage: preset.icon)
                    }
                }
            } label: {
                HStack {
                    Image(systemName: "sparkles")
                    Text("Choose a 1-Click Sample Experiment…")
                        .font(.system(size: 11, weight: .semibold))
                    Spacer()
                }
                .frame(maxWidth: .infinity)
            }
            .menuStyle(.borderedButton)
            .controlSize(.small)
        }
        .padding(10)
        .background(
            LinearGradient(
                colors: [Color(red: 0.10, green: 0.25, blue: 0.34).opacity(0.65), BlueHour.pillSurface.opacity(0.85)],
                startPoint: .topLeading,
                endPoint: .bottomTrailing
            ),
            in: RoundedRectangle(cornerRadius: 9)
        )
        .overlay(
            RoundedRectangle(cornerRadius: 9)
                .strokeBorder(Color(red: 0.31, green: 0.85, blue: 0.78).opacity(0.45), lineWidth: 1)
        )
    }
    
    var body: some View {
        VStack(spacing: 0) {
            setupStepTabStrip
            Divider()
            ScrollView {
                VStack(alignment: .leading, spacing: 16) {
                    switch inspectorStep {
                    case .crystal:
                        quickStartExperimentBanner
                        systemGeometryCard
                        boundaryAndLoadingSelectorCard
                        latticePotentialCard
                        if config.mode != .volterra && config.mode != .shift {
                            GeometryBuilderView(
                                config: config,
                                previewManager: previewManager,
                                displayMode: .inspectorMicrostructure,
                                onSwitchToCarvedMesh: { displayTab = .initialMesh },
                                onOpenDocs: onOpenDocs
                            )
                        }
                        crystalOrientationCard
                    case .shapeAndBCs:
                        if config.mode == .volterra || config.mode == .shift {
                            stepNotApplicableNote(
                                "Shape & boundary conditions are fixed by the selected loading mode (\(config.mode.displayName)). Switch the loading mode in step 3 to carve a custom specimen."
                            )
                        } else {
                            boundaryAndLoadingSelectorCard
                            specimenShapeCard
                            if config.geometryConfig.enabled && !config.isEffectivePeriodicPBC {
                                GeometryBuilderView(
                                    config: config,
                                    previewManager: previewManager,
                                    displayMode: .inspectorShapeAndBCs,
                                    onSwitchToCarvedMesh: { displayTab = .initialMesh },
                                    onOpenDocs: onOpenDocs
                                )
                            }
                        }
                    case .protocolAndSchedule:
                        boundaryAndLoadingSelectorCard
                        loadingModeCard
                        protocolScheduleCard
                    case .solver:
                        if uiComplexityMode == .guided {
                            stepNotApplicableNote(
                                "These solver settings already have sensible defaults. You only need to change them for very large crystals or for stability (eigenvalue) studies."
                            )
                        }
                        preconditionerCard
                        solverAndOutputCards
                    }
                }
                .frame(maxWidth: .infinity, alignment: .leading)
                .padding(14)
            }
            
            Divider()
            stepFooterNavigation
        }
        .background(BlueHour.sidebarBackground)
    }
    
    // MARK: Step Tab Strip
    private var setupStepTabStrip: some View {
        HStack(spacing: 4) {
            ForEach(SetupInspectorStep.allCases) { step in
                let isSel = (inspectorStep == step)
                Button(action: { inspectorStep = step }) {
                    VStack(spacing: 2) {
                        Image(systemName: step.icon)
                            .font(.system(size: 13, weight: .semibold))
                        Text(step.shortTitle)
                            .font(.system(size: 10, weight: isSel ? .bold : .medium))
                            .lineLimit(1)
                    }
                    .foregroundColor(isSel ? .white : BlueHour.ice.opacity(0.8))
                    .frame(maxWidth: .infinity, minHeight: 40)
                    .background(
                        isSel ? Color(red: 0.12, green: 0.64, blue: 0.54) : BlueHour.midnight.opacity(0.55),
                        in: RoundedRectangle(cornerRadius: 7)
                    )
                    .overlay(
                        RoundedRectangle(cornerRadius: 7)
                            .strokeBorder(isSel ? Color(red: 0.42, green: 0.96, blue: 0.84) : BlueHour.steel.opacity(0.4), lineWidth: isSel ? 1.2 : 0.8)
                    )
                    .contentShape(Rectangle())
                }
                .buttonStyle(.plain)
                .help(step.helpText)
            }
        }
        .padding(.horizontal, 10)
        .padding(.vertical, 8)
        .background(BlueHour.controlSurface)
    }
    
    // MARK: Back / Next footer
    private var stepFooterNavigation: some View {
        let all = SetupInspectorStep.allCases
        let idx = all.firstIndex(of: inspectorStep) ?? 0
        return HStack {
            if idx > 0 {
                Button(action: { inspectorStep = all[idx - 1] }) {
                    Label(all[idx - 1].shortTitle, systemImage: "chevron.left")
                }
                .controlSize(.small)
            }
            Spacer()
            Text("Step \(idx + 1) of \(all.count)")
                .font(.system(size: 10))
                .foregroundColor(BlueHour.mist)
            Spacer()
            if idx < all.count - 1 {
                Button(action: { inspectorStep = all[idx + 1] }) {
                    HStack(spacing: 4) {
                        Text("Next: \(all[idx + 1].shortTitle)")
                        Image(systemName: "chevron.right")
                    }
                }
                .buttonStyle(.borderedProminent)
                .controlSize(.small)
            }
        }
        .padding(.horizontal, 12)
        .padding(.vertical, 8)
        .background(BlueHour.controlSurface)
    }
    
    private func stepNotApplicableNote(_ text: String) -> some View {
        HStack(alignment: .top, spacing: 8) {
            Image(systemName: "info.circle.fill")
                .foregroundColor(BlueHour.ice)
            Text(text)
                .font(.system(size: 11))
                .foregroundColor(BlueHour.mist)
                .fixedSize(horizontal: false, vertical: true)
        }
        .padding(10)
        .frame(maxWidth: .infinity, alignment: .leading)
        .background(BlueHour.pillSurface, in: RoundedRectangle(cornerRadius: 8))
    }
    
    @ViewBuilder
    private var systemGeometryCard: some View {
                // 1. System Size
                GroupBox(label: HStack {
                    Label("System Size", systemImage: "grid")
                    Spacer()
                    HelpPopoverButton(
                        title: "System Size (NX × NY)",
                        summary: "Sets the number of discrete crystal unit cells along X (columns) and Y (rows). Each node represents a discrete atom or lattice site.",
                        formulaOrTip: "• 100×100 (10,000 atoms): Fast interactive runs (< 1 min).\n• 200×200 to 400×400: High-resolution dislocation ensembles & multi-plane DDD.",
                        docTopic: .quickStart,
                        onOpenDocs: onOpenDocs
                    )
                }) {
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
                                ForEach([(50, 50), (100, 100), (150, 150), (200, 200), (300, 300), (400, 400)], id: \.0) { size in
                                    let isSelected = (config.nx == size.0 && config.ny == size.1)
                                    Button(action: {
                                        config.setSize(nx: size.0, ny: size.1)
                                    }) {
                                        HStack(spacing: 3) {
                                            if isSelected {
                                                Image(systemName: "checkmark.circle.fill")
                                                    .font(.system(size: 8.5, weight: .bold))
                                            }
                                            Text("\(size.0)×\(size.1)")
                                                .font(.system(size: 10, weight: isSelected ? .bold : .medium, design: .monospaced))
                                        }
                                        .foregroundColor(isSelected ? .white : BlueHour.ice.opacity(0.88))
                                        .padding(.horizontal, 7)
                                        .padding(.vertical, 3.5)
                                        .background(
                                            isSelected
                                                ? Color(red: 0.12, green: 0.64, blue: 0.54)
                                                : BlueHour.midnight.opacity(0.85),
                                            in: RoundedRectangle(cornerRadius: 5)
                                        )
                                        .overlay(
                                            RoundedRectangle(cornerRadius: 5)
                                                .strokeBorder(
                                                    isSelected
                                                        ? Color(red: 0.42, green: 0.96, blue: 0.84)
                                                        : BlueHour.steel.opacity(0.55),
                                                    lineWidth: isSelected ? 1.3 : 1.0
                                                )
                                        )
                                    }
                                    .buttonStyle(.plain)
                                }
                            }
                        }
                    }
                    .frame(maxWidth: .infinity, alignment: .leading)
                    .padding(8)
                }
    }
    
    // MARK: - Unified Boundary Regime & Loading Selector Card
    @ViewBuilder
    private var boundaryAndLoadingSelectorCard: some View {
        GroupBox(label: HStack {
            Label("Boundary Regime & Mechanical Loading", systemImage: "arrow.triangle.swap")
            Spacer()
            HelpPopoverButton(
                title: "Boundary Conditions & Loading Mode",
                summary: "Controls whether the crystal is wrapped as a continuous 2D periodic torus (PBC) or treated as a finite specimen with surface grips/rollers.\n\n• Full 2D Periodic (PBC): Enforces periodic boundary conditions across X and Y. Deforms via macroscopic simple shear (F_xy = ±α).\n• Surface Grips / Custom BCs: Fixes or pulls specific specimen surfaces (e.g. top edge, bottom edge, indenter).",
                formulaOrTip: "• For pure shear on an infinite crystal: Select 'Full 2D Periodic (PBC)'.\n• For tension/compression or cutouts: Select 'Surface Grips / Custom BCs'.",
                docTopic: .geometryAndBCs,
                onOpenDocs: onOpenDocs
            )
        }) {
            VStack(alignment: .leading, spacing: 12) {
                // 1. Boundary Regime Selector
                VStack(alignment: .leading, spacing: 5) {
                    Text("Boundary Condition Regime:")
                        .font(.caption.bold())
                        .foregroundColor(BlueHour.ice)
                    
                    HStack(spacing: 6) {
                        ForEach(ActiveBoundaryRegime.allCases) { regime in
                            let isSel = (config.activeBoundaryRegime == regime)
                            Button(action: {
                                config.setBoundaryRegime(regime)
                            }) {
                                HStack(spacing: 5) {
                                    Image(systemName: isSel ? "checkmark.circle.fill" : regime.icon)
                                        .font(.system(size: 11, weight: .bold))
                                    Text(regime.title)
                                        .font(.system(size: 11, weight: isSel ? .bold : .medium))
                                }
                                .foregroundColor(isSel ? .white : BlueHour.ice.opacity(0.85))
                                .frame(maxWidth: .infinity)
                                .padding(.vertical, 6)
                                .background(
                                    isSel
                                        ? (regime == .periodicPBC ? Color(red: 0.12, green: 0.64, blue: 0.54) : Color(red: 0.20, green: 0.50, blue: 0.85))
                                        : BlueHour.midnight.opacity(0.75),
                                    in: RoundedRectangle(cornerRadius: 6)
                                )
                                .overlay(
                                    RoundedRectangle(cornerRadius: 6)
                                        .strokeBorder(
                                            isSel
                                                ? (regime == .periodicPBC ? Color(red: 0.45, green: 0.95, blue: 0.80) : Color(red: 0.50, green: 0.80, blue: 1.0))
                                                : BlueHour.steel.opacity(0.55),
                                            lineWidth: isSel ? 1.4 : 0.9
                                        )
                                )
                            }
                            .buttonStyle(.plain)
                        }
                    }
                    
                    Text(config.isEffectivePeriodicPBC
                         ? "✓ Full 2D Periodic (T² torus): X ∈ [0, Lx) and Y ∈ [0, Ly). No artificial surfaces. Driven by macroscopic simple shear F_xy = ±α."
                         : "Specimen has physical outer surfaces (clamped bottom, driven top, or traction-free sides/cutouts).")
                        .font(.caption2)
                        .foregroundColor(config.isEffectivePeriodicPBC ? Color(red: 0.45, green: 0.90, blue: 0.62) : .secondary)
                }
                
                Divider()
                
                // 2. Mechanical Loading Type Selector
                VStack(alignment: .leading, spacing: 6) {
                    HStack {
                        Text("Mechanical Deformation / Loading:")
                            .font(.caption.bold())
                            .foregroundColor(BlueHour.ice)
                        Spacer()
                        Text("Measures: \(config.activeStressComponentLabel)")
                            .font(.system(size: 10, weight: .semibold, design: .monospaced))
                            .foregroundColor(Color(red: 1.0, green: 0.82, blue: 0.35))
                    }
                    
                    let active = config.activeMechanicalLoading
                    LazyVGrid(columns: [GridItem(.flexible()), GridItem(.flexible())], spacing: 6) {
                        ForEach(ActiveMechanicalLoading.primaryGridModes) { load in
                            let isSel = (active == load)
                            Button(action: {
                                config.applyMechanicalLoading(load)
                            }) {
                                HStack(spacing: 5) {
                                    Image(systemName: isSel ? "checkmark.circle.fill" : load.sfSymbol)
                                        .font(.system(size: 11, weight: .bold))
                                        .foregroundColor(isSel ? .white : Color(red: 0.31, green: 0.85, blue: 0.78))
                                    VStack(alignment: .leading, spacing: 1) {
                                        Text(load.shortTitle)
                                            .font(.system(size: 10.5, weight: isSel ? .bold : .medium))
                                            .lineLimit(1)
                                        Text(load.supportsPeriodicPBC ? "PBC or Surface" : "Surface BCs")
                                            .font(.system(size: 8.5))
                                            .foregroundColor(isSel ? .white.opacity(0.85) : BlueHour.mist)
                                    }
                                    Spacer(minLength: 0)
                                }
                                .foregroundColor(isSel ? .white : BlueHour.ice.opacity(0.90))
                                .padding(.horizontal, 7)
                                .padding(.vertical, 5)
                                .background(
                                    isSel ? Color(red: 0.14, green: 0.58, blue: 0.52) : BlueHour.midnight.opacity(0.80),
                                    in: RoundedRectangle(cornerRadius: 6)
                                )
                                .overlay(
                                    RoundedRectangle(cornerRadius: 6)
                                        .strokeBorder(
                                            isSel ? Color(red: 0.45, green: 0.95, blue: 0.80) : BlueHour.steel.opacity(0.50),
                                            lineWidth: isSel ? 1.3 : 0.8
                                        )
                                )
                            }
                            .buttonStyle(.plain)
                            .help(load.detailSubtitle)
                        }
                    }
                }
                
                // 3. Prominent Live Confirmation Banner
                HStack(alignment: .center, spacing: 8) {
                    Circle()
                        .fill(config.isEffectivePeriodicPBC ? Color(red: 0.45, green: 0.90, blue: 0.62) : Color(red: 0.31, green: 0.85, blue: 0.78))
                        .frame(width: 8, height: 8)
                    
                    VStack(alignment: .leading, spacing: 2) {
                        HStack(spacing: 6) {
                            Text("Active Physics:")
                                .font(.system(size: 9.5, weight: .bold))
                                .foregroundColor(BlueHour.mist)
                            Text(config.activePhysicsSummaryBadge)
                                .font(.system(size: 10.5, weight: .bold, design: .monospaced))
                                .foregroundColor(BlueHour.ice)
                        }
                        Text(config.activeMechanicalLoading.detailSubtitle)
                            .font(.system(size: 9))
                            .foregroundColor(BlueHour.mist)
                            .lineLimit(1)
                    }
                    Spacer()
                }
                .padding(8)
                .background(
                    (config.isEffectivePeriodicPBC
                        ? Color(red: 0.12, green: 0.64, blue: 0.54).opacity(0.25)
                        : Color(red: 0.20, green: 0.50, blue: 0.85).opacity(0.25)),
                    in: RoundedRectangle(cornerRadius: 6)
                )
                .overlay(
                    RoundedRectangle(cornerRadius: 6)
                        .strokeBorder(
                            (config.isEffectivePeriodicPBC
                                ? Color(red: 0.45, green: 0.90, blue: 0.62).opacity(0.55)
                                : Color(red: 0.40, green: 0.75, blue: 1.0).opacity(0.55)),
                            lineWidth: 1
                        )
                )
            }
            .frame(maxWidth: .infinity, alignment: .leading)
            .padding(8)
        }
    }
    
    @ViewBuilder
    private var loadingModeCard: some View {
                // 2. Loading Mode & Remeshing
                GroupBox(label: HStack {
                    Label("Loading Mode & Remeshing", systemImage: "arrow.triangle.swap")
                    Spacer()
                    HelpPopoverButton(
                        title: "Loading Mode & Delaunay Remeshing",
                        summary: "Selects how the crystal is deformed and whether topological Delaunay flips (surgery) are active.",
                        formulaOrTip: "• Adaptive Remeshing (Delaunay Surgery) MUST be ON for plastic slip and dislocation nucleation.\n• Turn it OFF only to study purely hyperelastic (fixed-connectivity) response.",
                        docTopic: .potentialsAndRemeshing,
                        onOpenDocs: onOpenDocs
                    )
                }) {
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
                            
                            Toggle("Relax Macroscopic Stress to Zero (Full Cauchy)", isOn: $config.relaxMacroStress)
                                .font(.caption)
                            
                            if config.relaxMacroStress {
                                HStack {
                                    Text("Stress Tolerance (max |σ|):")
                                        .font(.caption)
                                    Spacer()
                                    DoubleField(title: "1e-5", value: $config.stressTol, width: 75)
                                }
                                HStack {
                                    Text("Max Macro Iterations:")
                                        .font(.caption)
                                    Spacer()
                                    IntField(title: "30", value: $config.maxMacroIters, width: 75)
                                }
                            }
                            
                            Text(config.latticeSymmetry == .triangular ?
                                 "Uses Isotropic 2D Volterra solution (Eq. 15) with ν = C₁₂/C₁₁ and equilibrium spacing a₀ computed from active potential." :
                                 "Uses Anisotropic Stroh solution (Eqs. A5–A16) with (κ, ξ, ζ, g₁, g₂, g₃) and equilibrium spacing a₀ computed from active potential.")
                                .font(.caption2)
                                .foregroundColor(.accentColor)
                        }
                        
                        if config.mode != .volterra && config.mode != .shift {
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
                        }
                        
                        Toggle("Adaptive Remeshing (Delaunay Surgery)", isOn: $config.enableRemeshing)
                            .toggleStyle(.switch)
                        
                        if config.mode == .shift {
                            Toggle("Perturb Initial Mesh Diagonal (-1e-7)", isOn: $config.shiftPerturbMesh)
                                .font(.caption)
                            Text("Special staircase boundary conditions: Left & bottom edges are supported; upper-left window is shifted horizontally and lower-right window is pushed vertically in alternating steps.")
                                .font(.caption2)
                                .foregroundColor(.accentColor)
                        } else if config.mode != .volterra {
                            Divider().padding(.vertical, 2)
                            
                            HStack {
                                Text("Active Boundary Regime:")
                                    .font(.caption)
                                    .foregroundColor(.secondary)
                                Spacer()
                                HStack(spacing: 4) {
                                    Circle()
                                        .fill(config.isEffectivePeriodicPBC ? Color(red: 0.45, green: 0.90, blue: 0.62) : Color(red: 0.31, green: 0.85, blue: 0.78))
                                        .frame(width: 6, height: 6)
                                    Text(config.isEffectivePeriodicPBC ? "2D Periodic (PBC)" : "Surface Grips / Custom BCs")
                                        .font(.system(size: 10, weight: .bold))
                                        .foregroundColor(config.isEffectivePeriodicPBC ? Color(red: 0.45, green: 0.90, blue: 0.62) : BlueHour.ice)
                                }
                                .padding(.horizontal, 7)
                                .padding(.vertical, 3)
                                .background(BlueHour.pillSurface, in: RoundedRectangle(cornerRadius: 5))
                            }
                        }
                    }
                    .frame(maxWidth: .infinity, alignment: .leading)
                    .padding(8)
                }
    }
    
    @ViewBuilder
    private var specimenShapeCard: some View {
                if config.mode != .volterra && config.mode != .shift {
                    // 2b. Macroscopic Geometry & Regional Boundary Conditions Studio
                    GroupBox(label: HStack {
                        Label("Specimen Shape & Custom BCs", systemImage: "skew")
                        Spacer()
                        HelpPopoverButton(
                            title: "Macroscopic Geometry, Dislocations & BCs",
                            summary: "Enables custom specimen shapes (notches, holes, pillars), polycrystals (Bicrystal, Tri-Crystal, Voronoi), Periodic Dislocations, and regional boundary conditions.",
                            formulaOrTip: "• Click 'Open Geometry Studio' to visually edit grain boundaries, seed up to 100 edge dislocations, or configure Nanoindentation.",
                            docTopic: .geometryAndBCs,
                            onOpenDocs: onOpenDocs
                        )
                    }) {
                        VStack(alignment: .leading, spacing: 10) {
                            Toggle(isOn: Binding(
                                get: { config.geometryConfig.enabled },
                                set: { enabled in
                                    config.geometryConfig.enabled = enabled
                                    if enabled {
                                        config.boundaryCondition = .fixed
                                    }
                                }
                            )) {
                                Text("Carve Macroscopic Shape & Custom BCs")
                                    .font(.caption.bold())
                            }
                            .toggleStyle(.switch)
                            
                            if config.geometryConfig.enabled {
                                // Preset quick-selector
                                VStack(alignment: .leading, spacing: 4) {
                                    Text("Engineering Specimen Presets:")
                                        .font(.caption2)
                                        .foregroundColor(.secondary)
                                    ScrollView(.horizontal, showsIndicators: false) {
                                        HStack(spacing: 5) {
                                            ForEach(GeometryEngineeringPreset.allCases) { preset in
                                                let isSel = (config.geometryConfig.preset_name == preset.rawValue)
                                                Button(action: {
                                                    config.geometryConfig.applyPreset(preset)
                                                }) {
                                                    Text(preset.title)
                                                        .font(.system(size: 10, weight: isSel ? .bold : .medium))
                                                }
                                                .buttonStyle(.bordered)
                                                .tint(isSel ? BlueHour.accent : nil)
                                                .controlSize(.mini)
                                            }
                                        }
                                    }
                                }
                                
                                HStack(spacing: 8) {
                                    Label(config.geometryConfig.outer_shape.displayName, systemImage: config.geometryConfig.outer_shape.icon)
                                        .font(.system(size: 10, weight: .semibold))
                                        .foregroundColor(BlueHour.ice)
                                    Text("• \(config.geometryConfig.cutouts.count) cutout(s)")
                                        .font(.system(size: 10))
                                        .foregroundColor(BlueHour.mist)
                                    Text("• \(config.geometryConfig.boundary_rules.count) BC rule(s)")
                                        .font(.system(size: 10))
                                        .foregroundColor(BlueHour.mist)
                                }
                                
                                Text("Discrete \(config.latticeSymmetry.displayName) unit cell size a₀ is strictly preserved.")
                                    .font(.caption2)
                                    .foregroundColor(.secondary)
                            }
                            
                            HStack(spacing: 8) {
                                Button(action: {
                                    if !config.geometryConfig.enabled {
                                        config.geometryConfig.enabled = true
                                        config.boundaryCondition = .fixed
                                    }
                                    displayTab = .geometryStudio
                                }) {
                                    Label("Show Specimen on Canvas", systemImage: "slider.horizontal.below.square.and.square.filled")
                                        .frame(maxWidth: .infinity)
                                }
                                .buttonStyle(.borderedProminent)
                                .controlSize(.small)
                                
                                if config.geometryConfig.enabled {
                                    Button(action: {
                                        let root = SimulationConfig.detectProjectRoot()
                                        let scratchDir = "\(root)/scratch"
                                        try? FileManager.default.createDirectory(atPath: scratchDir, withIntermediateDirectories: true)
                                        let geomPath = "\(scratchDir)/geometry_config.json"
                                        try? config.geometryConfig.writeJSON(to: geomPath)
                                        
                                        previewManager.nx = config.nx
                                        previewManager.ny = config.ny
                                        previewManager.latticeSymmetry = config.latticeSymmetry
                                        previewManager.potentialMode = config.potentialMode
                                        previewManager.boundaryCondition = .fixed
                                        previewManager.beta = config.beta
                                        previewManager.bulkK = config.bulkK
                                        previewManager.useRotation = config.useRotation
                                        previewManager.rotationAngle = config.rotationAngle
                                        previewManager.geometryJSONPath = geomPath
                                        previewManager.generateInitialMeshPreview()
                                        displayTab = .initialMesh
                                    }) {
                                        Label("Carved Mesh", systemImage: "eye.fill")
                                    }
                                    .buttonStyle(.bordered)
                                    .controlSize(.small)
                                }
                            }
                        }
                        .frame(maxWidth: .infinity, alignment: .leading)
                        .padding(8)
                    }
                }
    }
    
    @ViewBuilder
    private var latticePotentialCard: some View {
                // 3. Lattice & Potential Model
                GroupBox(label: HStack {
                    Label("Lattice & Potential", systemImage: "atom")
                    Spacer()
                    HelpPopoverButton(
                        title: "Lattice Symmetry & Interatomic Potential",
                        summary: "Selects the Bravais lattice symmetry (Square or Hexagonal/Triangular) and the strain-energy density W(C).",
                        formulaOrTip: "• Zanzotto GL Strain Energy: Exact invariant polynomial with GL3 reduction.\n  – Square preset: β = -0.25, K = 5.0\n  – Triangular preset: β = 4.0, K = 5.0\n• Lennard-Jones v2/v3: Pairwise potentials coarse-grained via Cauchy-Born.",
                        docTopic: .potentialsAndRemeshing,
                        onOpenDocs: onOpenDocs
                    )
                }) {
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
    }
    
    @ViewBuilder
    private var crystalOrientationCard: some View {
                if config.mode != .volterra {
                    // 3b. Rotated Crystal Orientation / CSL
                    GroupBox(label: HStack {
                        Label("Crystal Orientation & Mesh Preview", systemImage: "rotate.right")
                        Spacer()
                        HelpPopoverButton(
                            title: "Rotated Crystal & Initial Mesh Preview",
                            summary: "Optionally rotates the reference crystal lattice by an angle θ (or a discrete Coincident Site Lattice Σ angle under PBC) and previews the discrete finite-element triangulation.",
                            formulaOrTip: "• Click 'Preview Initial Mesh' at any time to inspect the discrete lattice before launching the C++ solver.",
                            docTopic: .grainBoundaries,
                            onOpenDocs: onOpenDocs
                        )
                    }) {
                        VStack(alignment: .leading, spacing: 10) {
                            Toggle(isOn: $config.useRotation) {
                                Text("Enable Rotated Lattice")
                                    .font(.caption.bold())
                            }
                            .toggleStyle(.checkbox)
                            
                            if config.useRotation {
                                if config.boundaryCondition == .fixed || config.geometryConfig.enabled || config.mode == .shift {
                                    VStack(alignment: .leading, spacing: 6) {
                                        Text(config.geometryConfig.enabled
                                             ? "Arbitrary Continuous Rotation (Finite Specimen in Geometry Studio):"
                                             : (config.mode == .shift ? "Arbitrary Continuous Rotation (Staircase BC):" : "Arbitrary Continuous Rotation (Fixed BC):"))
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
                                        Text("The finite crystal block can be rotated to any continuous angle θ ∈ [0, 360°].")
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
                                previewManager.boundaryCondition = (config.mode == .shift || config.geometryConfig.enabled) ? .fixed : config.boundaryCondition
                                previewManager.beta = config.beta
                                previewManager.bulkK = config.bulkK
                                previewManager.useRotation = config.useRotation
                                previewManager.rotationAngle = config.rotationAngle
                                previewManager.cslM = config.cslM
                                previewManager.cslN = config.cslN
                                previewManager.cslKx = config.cslKx
                                previewManager.cslKy = config.cslKy
                                if config.geometryConfig.enabled {
                                    let root = SimulationConfig.detectProjectRoot()
                                    let scratchDir = "\(root)/scratch"
                                    try? FileManager.default.createDirectory(atPath: scratchDir, withIntermediateDirectories: true)
                                    let geomPath = "\(scratchDir)/geometry_config.json"
                                    try? config.geometryConfig.writeJSON(to: geomPath)
                                    previewManager.geometryJSONPath = geomPath
                                } else {
                                    previewManager.geometryJSONPath = ""
                                }
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
                }
    }
    
    @ViewBuilder
    private var preconditionerCard: some View {
                // 4. Preconditioning & Energy Relaxation
                VStack(alignment: .leading, spacing: 16) {
                    GroupBox(label: HStack {
                        Label("L-BFGS Preconditioner", systemImage: "bolt.fill")
                        Spacer()
                        HelpPopoverButton(
                            title: "L-BFGS Metric Preconditioner",
                            summary: "Accelerates energy minimization by preconditioning the L-BFGS search direction with a sparse Cholesky factorization (CHOLMOD).",
                            formulaOrTip: "• 'Laplacian (Graph)' or 'Elastic Tangent' dramatically reduces iterations for large crystals (N ≥ 150).\n• Refresh = 0 re-factorizes adaptively when convergence slows.",
                            docTopic: .solverAndEigenvalues,
                            onOpenDocs: onOpenDocs
                        )
                    }) {
                        VStack(alignment: .leading, spacing: 10) {
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
                }
    }
    
    @ViewBuilder
    private var protocolScheduleCard: some View {
                // 5. Loading Schedule
                if config.mode == .shift {
                    GroupBox(label: Label("Staircase Shift Schedule", systemImage: "stairs")) {
                        VStack(alignment: .leading, spacing: 8) {
                            HStack {
                                Text("Shift Amplitude (in a₀):")
                                    .font(.subheadline)
                                Spacer()
                                DoubleField(title: "2.0", value: $config.shiftAmplitude, width: 80)
                            }
                            HStack {
                                Text("Horizontal Steps (per a₀):")
                                    .font(.subheadline)
                                Spacer()
                                IntField(title: "2", value: $config.shiftHorizontalSteps, width: 80)
                            }
                            HStack {
                                Text("Vertical Steps (per a₀):")
                                    .font(.subheadline)
                                Spacer()
                                IntField(title: "2", value: $config.shiftVerticalSteps, width: 80)
                            }
                        }
                        .padding(8)
                    }
                } else if config.mode != .volterra {
                    GroupBox(label: HStack {
                        Label("Simulation Protocol & Schedule", systemImage: "slider.horizontal.below.square.filled.and.square")
                        Spacer()
                        HelpPopoverButton(
                            title: "Simulation Protocol & Loading Schedule",
                            summary: "Choose between 3 protocols:\n1. Relax Only (α = 0): Equilibrates defects or grain boundaries with zero external load.\n2. Relax → Load: Two-stage handover (equilibrates at α = 0, then shears/loads from α_start to α_end).\n3. Sweep Δθ: Automated multi-angle grain boundary energy curve builder γ_GB(Δθ).",
                            formulaOrTip: "• Step size |dα| = 0.0005–0.001 captures sharp plastic avalanches.\n• Re-anchor reference frame (X₀ ← x_eq) prevents boundary snap at Stage 1 → Stage 2 handover.",
                            docTopic: .quickStart,
                            onOpenDocs: onOpenDocs
                        )
                    }) {
                        VStack(alignment: .leading, spacing: 10) {
                            // 3-Mode Protocol Switcher: 1. Relax Only | 2. Relax -> Load | 3. Sweep Δθ
                            Picker("Protocol", selection: $config.protocolMode) {
                                ForEach(SimulationProtocolMode.allCases) { mode in
                                    Text(mode.shortTitle).tag(mode)
                                }
                            }
                            .pickerStyle(.segmented)
                            .labelsHidden()
                            .onChange(of: config.protocolMode) { newMode in
                                config.geometryConfig.microstructure.relax_only = (newMode != .relaxThenLoad)
                                if newMode == .relaxThenLoad {
                                    config.geometryConfig.microstructure.pre_relax_gb = true
                                } else {
                                    runner.requestedChartMode = .gbEnergy
                                    displayTab = .curve
                                }
                            }
                            
                            if config.protocolMode == .relaxOnly {
                                // MODE 1: RELAX ONLY
                                VStack(alignment: .leading, spacing: 8) {
                                    HStack(spacing: 6) {
                                        Image(systemName: "checkmark.seal.fill")
                                            .foregroundColor(.green)
                                            .font(.caption)
                                        Text("Relaxation Step")
                                            .font(.caption.bold())
                                    }
                                    Text("Relaxes the grain boundary / crystal under the Relaxation Step BCs, performs Delaunay core reconstruction until ΔE converges, saves the equilibrium VTK, and records γ_GB(Δθ).")
                                        .font(.caption2)
                                        .foregroundColor(.secondary)
                                }
                                .padding(8)
                                .background(Color.green.opacity(0.08), in: RoundedRectangle(cornerRadius: 6))
                                
                            } else if config.protocolMode == .relaxThenLoad {
                                // MODE 2: RELAX -> LOAD (Two-Stage Handover)
                                VStack(alignment: .leading, spacing: 7) {
                                    HStack(spacing: 5) {
                                        Circle().fill(Color.green).frame(width: 7, height: 7)
                                        Text("1. Relaxation Step (α = 0, Noise OFF)")
                                            .font(.caption.bold())
                                    }
                                    Toggle("Equilibrate crystal / GB before loading", isOn: $config.geometryConfig.microstructure.pre_relax_gb)
                                        .font(.caption)
                                    if config.geometryConfig.microstructure.pre_relax_gb {
                                        Text("Relaxation Step BCs: \(config.geometryConfig.stage1_boundary_rules.filter(\.enabled).count) active rules • Indenter OFF • Thermal noise OFF")
                                            .font(.caption2)
                                            .foregroundColor(Color(red: 0.45, green: 0.92, blue: 0.62))
                                        Toggle("Re-anchor reference frame (X₀ ← x_eq)", isOn: $config.geometryConfig.microstructure.reanchor_ref_after_relax)
                                            .font(.caption)
                                            .help("Prevents boundary snap by measuring Loading Step boundary displacements and indenter surface contact from the relaxed positions.")
                                    }
                                    
                                    Divider().padding(.vertical, 2)
                                    
                                    HStack(spacing: 5) {
                                        Circle().fill(Color.blue).frame(width: 7, height: 7)
                                        Text("2. Loading Step (Affine to α_start, Noise ON)")
                                            .font(.caption.bold())
                                        Spacer()
                                        Button("BCs & Indenter →") {
                                            inspectorStep = .shapeAndBCs
                                            displayTab = .geometryStudio
                                        }
                                        .buttonStyle(.bordered)
                                        .controlSize(.mini)
                                    }
                                    
                                    let indenterActive = config.geometryConfig.indenter.enabled
                                    let hasDrivenBC = config.geometryConfig.boundary_rules.contains {
                                        $0.enabled && ($0.kind == .prescribed || $0.kind == .affine)
                                    }
                                    if indenterActive {
                                        HStack(spacing: 5) {
                                            Image(systemName: "arrow.down.to.line.compact")
                                                .foregroundColor(.yellow)
                                                .font(.caption2)
                                            Text("Nanoindenter (\(config.geometryConfig.indenter.shape.displayName), δ_max = \(String(format: "%.1f", config.geometryConfig.indenter.max_depth_a0))a₀) + \(config.geometryConfig.boundary_rules.filter(\.enabled).count) BCs")
                                                .font(.caption2)
                                                .foregroundColor(BlueHour.ice)
                                        }
                                    } else if hasDrivenBC {
                                        HStack(spacing: 5) {
                                            Image(systemName: "arrow.up.and.down.and.sparkles")
                                                .foregroundColor(BlueHour.ice)
                                                .font(.caption2)
                                            Text("Using BCs & Indenter rules (\(config.geometryConfig.boundary_rules.filter(\.enabled).count) active, \(config.geometryConfig.primary_stress_component))")
                                                .font(.caption2)
                                                .foregroundColor(BlueHour.ice)
                                        }
                                    } else {
                                        HStack(spacing: 5) {
                                            Image(systemName: "exclamationmark.triangle.fill")
                                                .foregroundColor(.orange)
                                                .font(.caption2)
                                            Text("No driven load or Indenter enabled in BCs & Indenter")
                                                .font(.caption2)
                                                .foregroundColor(.orange)
                                        }
                                    }
                                    
                                    VStack(alignment: .leading, spacing: 6) {
                                        HStack {
                                            Text(indenterActive ? "Start depth ratio α_start:" : "Start strain |α_start|:")
                                                .font(.subheadline)
                                            Spacer()
                                            DoubleField(title: "0.0", value: $config.alphaStart, width: 80)
                                        }
                                        HStack {
                                            Text(indenterActive ? "End depth ratio α_end:" : "End strain |α_end|:")
                                                .font(.subheadline)
                                            Spacer()
                                            DoubleField(title: "0.25", value: $config.alphaEnd, width: 80)
                                        }
                                        HStack {
                                            Text("Step size |dα|:")
                                                .font(.subheadline)
                                            Spacer()
                                            DoubleField(title: "0.0005", value: $config.stepSize, width: 90)
                                        }
                                        Divider()
                                        HStack {
                                            Text("Save VTK every N steps:")
                                                .font(.subheadline)
                                            Spacer()
                                            IntField(title: "50", value: $config.saveVTKInterval, width: 70)
                                        }
                                    }
                                }
                                .padding(8)
                                .background(BlueHour.pillSurface.opacity(0.55), in: RoundedRectangle(cornerRadius: 6))
                                
                            } else {
                                // MODE 3: SWEEP Δθ (Continuous or CSL Grain Boundary Energy Curve Builder)
                                VStack(alignment: .leading, spacing: 8) {
                                    HStack(spacing: 6) {
                                        Image(systemName: "chart.xyaxis.line")
                                            .foregroundColor(.orange)
                                            .font(.caption)
                                        Text("Misorientation Sweep γ_GB(Δθ) (α = 0)")
                                            .font(.caption.bold())
                                    }
                                    
                                    Picker("Sweep Sampling", selection: $config.sweepContinuousMode) {
                                        Text("Continuous Δθ Grid").tag(true)
                                        Text("Discrete Σ Catalog").tag(false)
                                    }
                                    .pickerStyle(.segmented)
                                    .labelsHidden()
                                    .controlSize(.small)
                                    
                                    HStack {
                                        Text("Half-Angle θ_min (°):")
                                            .font(.caption)
                                        Spacer()
                                        DoubleField(title: "0.0", value: $config.sweepThetaMinDeg, width: 65, minValue: 0.0, maxValue: 44.5)
                                    }
                                    HStack {
                                        Text("Half-Angle θ_max (°):")
                                            .font(.caption)
                                        Spacer()
                                        DoubleField(title: config.latticeSymmetry == .triangular ? "30.0" : "45.0",
                                                    value: $config.sweepThetaMaxDeg,
                                                    width: 65,
                                                    minValue: 0.5,
                                                    maxValue: config.latticeSymmetry == .triangular ? 30.0 : 45.0)
                                    }
                                    
                                    if config.sweepContinuousMode {
                                        HStack(spacing: 6) {
                                            Text("Step d(Δθ) (°):")
                                                .font(.caption)
                                            Spacer()
                                            ForEach([1.0, 2.0, 3.0], id: \.self) { st in
                                                Button("\(Int(st))°") {
                                                    config.sweepDeltaThetaStepDeg = st
                                                }
                                                .buttonStyle(.bordered)
                                                .controlSize(.mini)
                                            }
                                            DoubleField(title: "2.0", value: $config.sweepDeltaThetaStepDeg, width: 52, minValue: 0.5, maxValue: 15.0)
                                        }
                                        Toggle("Include exact low-Σ cusp angles (Σ5, Σ7, Σ13)", isOn: $config.sweepIncludeExactCusps)
                                            .font(.caption2)
                                    } else {
                                        Toggle("Curated CSL catalog (0°–90° arch, 21 angles)", isOn: $config.sweepUseCuratedCSLOnly)
                                            .font(.caption2)
                                    }
                                    
                                    Toggle("Auto-adjust box to commensurate K·L_CSL per angle", isOn: $config.sweepAutoCommensurateBox)
                                        .font(.caption2)
                                        .help("Automatically sets each angle's box size to K·round(L_CSL) closest to your System Geometry NX×NY.")
                                    
                                    let cands = config.sweepCSLCandidates
                                    let previewStr = cands.prefix(8).map {
                                        $0.sigma > 0 ? String(format: "%.1f°(Σ%d)", $0.fullMisoriDeg, $0.sigma) : String(format: "%.1f°", $0.fullMisoriDeg)
                                    }.joined(separator: ", ")
                                    Text(config.sweepAutoCommensurateBox
                                         ? "Target box ≈ \(config.nx)×\(config.ny) a₀ (each angle snaps to K·L_CSL near \(config.nx)a₀)"
                                         : "Uniform box size: \(config.nx)×\(config.ny) a₀ across all \(cands.count) angles")
                                        .font(.system(size: 10, weight: .semibold, design: .monospaced))
                                        .foregroundColor(BlueHour.ice)
                                    Text("Scheduled \(cands.count) angles (Δθ = 2θ): \(previewStr)\(cands.count > 8 ? "..." : "")")
                                        .font(.system(size: 10, design: .monospaced))
                                        .foregroundColor(.secondary)
                                        .fixedSize(horizontal: false, vertical: true)
                                    
                                    Button(action: {
                                        displayTab = .curve
                                        runner.requestedChartMode = .gbEnergy
                                        runner.startCSLSweep(config: config, candidates: cands)
                                    }) {
                                        Label("Launch Δθ Angle Sweep (\(cands.count) runs)", systemImage: "play.circle.fill")
                                            .frame(maxWidth: .infinity)
                                    }
                                    .buttonStyle(.borderedProminent)
                                    .tint(.orange)
                                    .controlSize(.small)
                                    .disabled(runner.status.isRunning || cands.isEmpty)
                                    
                                    if !runner.gbEnergyPoints.isEmpty {
                                        HStack(spacing: 6) {
                                            Button(action: { runner.openSweepTextFile() }) {
                                                Label("Open angle_vs_energy.txt", systemImage: "doc.text")
                                                    .font(.system(size: 10, weight: .semibold))
                                                    .frame(maxWidth: .infinity)
                                            }
                                            .buttonStyle(.bordered)
                                            .controlSize(.mini)
                                            
                                            Button(action: { runner.openSweepRelaxedVTKFolder() }) {
                                                Label("Open Relaxed VTKs", systemImage: "folder")
                                                    .font(.system(size: 10, weight: .semibold))
                                                    .frame(maxWidth: .infinity)
                                            }
                                            .buttonStyle(.bordered)
                                            .controlSize(.mini)
                                        }
                                    }
                                }
                                .padding(8)
                                .background(Color.orange.opacity(0.08), in: RoundedRectangle(cornerRadius: 6))
                            }
                        }
                        .padding(8)
                    }
                }
                
    }
    
    @ViewBuilder
    private var solverAndOutputCards: some View {
                VStack(alignment: .leading, spacing: 16) {
                    // 6. Stability Monitor & Eigenvalues
                    if config.mode != .shift {
                        GroupBox(label: HStack {
                            Label("Stability & Eigenvalues", systemImage: "waveform.path.ecg")
                            Spacer()
                            HelpPopoverButton(
                                title: "Hessian Stability & Soft Modes",
                                summary: "Computes the lowest eigenvalues λ_min of the tangent stiffness matrix H = ∂²E/∂u² via shift-invert Lanczos to detect incipient plastic instabilities and soft dislocation core modes.",
                                formulaOrTip: "• λ_min → 0 signals a saddle-node bifurcation (plastic avalanche).",
                                docTopic: .solverAndEigenvalues,
                                onOpenDocs: onOpenDocs
                            )
                        }) {
                            VStack(alignment: .leading, spacing: 10) {
                                Toggle(config.mode == .volterra ? "Analyze Dislocation Core Stability" : "Compute Lowest Stiffness Eigenvalues", isOn: $config.enableEigenvalues)
                                    .toggleStyle(.switch)
                                
                                if config.enableEigenvalues {
                                    if config.mode == .volterra {
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
                                        Text("Computes the lowest eigenvalues and soft deformation modes of the relaxed dislocation core via Cholesky shift-invert Lanczos.")
                                            .font(.caption2)
                                            .foregroundColor(.secondary)
                                    } else {
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
                            }
                            .padding(8)
                        }
                    }
                    
                    // 7. Data Saving & Output Controls
                    GroupBox(label: HStack {
                        Label("Output & Checkpoint Controls", systemImage: "externaldrive")
                        Spacer()
                        HelpPopoverButton(
                            title: "VTK Snapshots & Checkpoints",
                            summary: "Controls how frequently Legacy VTK meshes (.vtk) for ParaView and binary simulation checkpoints are written to disk.",
                            formulaOrTip: "• 'Save Surgery VTK Snapshots' writes an extra VTK frame whenever Delaunay flips occur.",
                            docTopic: .solverAndEigenvalues,
                            onOpenDocs: onOpenDocs
                        )
                    }) {
                        VStack(alignment: .leading, spacing: 8) {
                            Toggle("Inner iterations trace (--trace=1)", isOn: $config.traceInnerIterations)
                            Toggle("Save Surgery VTK Snapshots (--surgery-vtk)", isOn: $config.saveSurgeryVTK)
                            Toggle("Save Legacy triangle_data", isOn: $config.saveTriangleData)
                            
                            if config.mode != .shift && config.mode != .volterra {
                                if config.protocolMode != .relaxThenLoad {
                                    HStack {
                                        Text("Save VTK every N steps:")
                                        Spacer()
                                        IntField(title: "50", value: $config.saveVTKInterval, width: 70)
                                    }
                                }
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
                        }
                        .padding(8)
                    }
                    
                    // 8. Hardware & Threads
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
                    
                    // 9. Working & Output Directory
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
    }
}

// MARK: - Unified Live Status & View Switcher Header
struct LiveStatusHeaderView: View {
    @ObservedObject var config: SimulationConfig
    @EnvironmentObject var runner: ProcessRunner
    
    var body: some View {
        HStack(spacing: 10) {
            // 1. Concise Status Pill
            HStack(spacing: 5) {
                Circle()
                    .fill(statusColor)
                    .frame(width: 8, height: 8)
                Text(shortStatusLabel)
                    .font(.system(size: 11, weight: .semibold))
                    .foregroundColor(BlueHour.ice)
                    .lineLimit(1)
            }
            .padding(.horizontal, 8)
            .padding(.vertical, 5)
            .background(BlueHour.pillSurface)
            .overlay(RoundedRectangle(cornerRadius: 6).strokeBorder(BlueHour.border, lineWidth: 1))
            .cornerRadius(6)
            .fixedSize()
            .help(runner.status.label)
            
            Divider()
                .frame(height: 18)
            
            // 3. Horizontally scrollable Live Metrics Pills so they never push or clip the header
            ScrollView(.horizontal, showsIndicators: false) {
                HStack(spacing: 6) {
                    if config.mode == .volterra, let root = runner.dislocationData {
                        let activeState = root.relaxed ?? root.unrelaxed
                        let maxMacroStress: Double = {
                            if let ms = root.finalMacroStress {
                                return max(abs(ms.xx), max(abs(ms.yy), abs(ms.xy)))
                            }
                            return abs(activeState.totalStress)
                        }()
                        MetricPill(
                            label: "State",
                            value: root.relaxed != nil ? "Relaxed" : "Unrelaxed",
                            valueColor: root.relaxed != nil ? .green : .orange
                        )
                        MetricPill(
                            label: "Total E",
                            value: String(format: "%.4f", activeState.totalEnergy)
                        )
                        MetricPill(
                            label: "E_core",
                            value: String(format: "%.4f", activeState.coreFit.coreEnergy)
                        )
                        MetricPill(
                            label: "max |σ|",
                            value: String(format: "%.2e", maxMacroStress),
                            valueColor: root.macroStressRelaxed ? .green : .primary
                        )
                    } else {
                        if !runner.currentStageLabel.isEmpty {
                            MetricPill(
                                label: "Protocol Stage",
                                value: runner.currentStageLabel,
                                valueColor: runner.currentStageLabel.contains("Stage 1") || runner.currentStageLabel.contains("Sweep") ? .green : .cyan
                            )
                        }
                        MetricPill(label: "Step", value: "\(runner.currentStep)")
                        MetricPill(label: config.mode == .shift ? "Shift (u/a₀)" : "Strain |α|", value: String(format: "%.4f", runner.currentAlpha))
                        MetricPill(label: "Stress σ", value: String(format: "%.4f", runner.currentStress))
                    }
                    if let lmin = runner.currentLambdaMin ?? runner.dislocationData?.lambdaMin {
                        MetricPill(
                            label: "λ_min",
                            value: String(format: "%.2e", lmin),
                            valueColor: lmin <= 0 ? .red : .teal
                        )
                    }
                    if config.mode != .volterra {
                        MetricPill(label: "Avalanches", value: "\(runner.avalancheCount)")
                    }
                    MetricPill(label: "Elapsed", value: formatTime(runner.elapsedTime))
                }
                .padding(.vertical, 1)
            }
            .frame(height: 34)
            .frame(maxWidth: .infinity, alignment: .leading)
        }
        .frame(height: 34)
        .frame(maxWidth: .infinity, alignment: .leading)
        .onAppear {
            runner.loadExistingIfAvailable(dir: config.targetDirectory)
        }
        .onChange(of: config.targetDirectory) { newDir in
            if !runner.status.isRunning {
                runner.loadExistingIfAvailable(dir: newDir)
            }
        }
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
    var valueColor: Color = BlueHour.ice
    
    var body: some View {
        VStack(spacing: 1) {
            Text(label)
                .font(.system(size: 9, weight: .semibold))
                .foregroundColor(BlueHour.mist)
            Text(value)
                .font(.system(size: 11.5, weight: .bold, design: .monospaced))
                .monospacedDigit()
                .foregroundColor(valueColor)
        }
        .frame(minWidth: 56, minHeight: 26, maxHeight: 26)
        .padding(.horizontal, 8)
        .padding(.vertical, 3)
        .background(BlueHour.pillSurface)
        .overlay(RoundedRectangle(cornerRadius: 6).strokeBorder(BlueHour.border, lineWidth: 1))
        .cornerRadius(6)
        .lineLimit(1)
        .fixedSize(horizontal: true, vertical: false)
    }
}

// MARK: - Terminal View
struct TerminalView: View {
    @EnvironmentObject var runner: ProcessRunner
    @Binding var filterText: String
    @Binding var autoScroll: Bool
    
    var filteredLines: [LogMessage] {
        if filterText.isEmpty { return runner.logLines }
        return runner.logLines.filter { $0.text.localizedCaseInsensitiveContains(filterText) }
    }
    
    private func colorForLine(_ msg: LogMessage) -> Color {
        if msg.isError { return Color.red }
        let t = msg.text
        if t.contains("[AVALANCHE") || t.contains("⚡") {
            return Color.orange
        }
        if t.contains("[Remesh") || t.contains("↻") || t.contains("mesh: remeshed") {
            return Color.cyan
        }
        if t.contains("[VTK Snapshot]") || t.contains("💾") {
            return Color.green
        }
        if t.contains("[Triangulation]") || t.contains("[MacroscopicGeometry]") || t.contains("[POTENTIAL]") || t.contains("[Setup]") {
            return Color.yellow.opacity(0.9)
        }
        return BlueHour.ice
    }

    var body: some View {
        VStack(spacing: 0) {
            // Search & Tools Bar
            HStack {
                Image(systemName: "magnifyingglass")
                    .foregroundColor(BlueHour.mist)
                TextField("Filter output...", text: $filterText)
                    .textFieldStyle(.plain)
                    .frame(maxWidth: 220)
                
                Spacer()
                
                Toggle("Auto-scroll", isOn: $autoScroll)
                    .toggleStyle(.checkbox)
                    .controlSize(.small)
                
                Text("\(runner.logLines.count) lines")
                    .font(.caption2)
                    .foregroundColor(BlueHour.mist)
            }
            .padding(.horizontal, 12)
            .padding(.vertical, 6)
            .background(BlueHour.controlSurface)
            
            Divider()
            
            // Console Log ScrollView
            ScrollViewReader { proxy in
                ScrollView {
                    LazyVStack(alignment: .leading, spacing: 2) {
                        ForEach(filteredLines) { msg in
                            Text(msg.text)
                                .font(.system(size: 11, weight: (msg.text.contains("[AVALANCHE") || msg.text.contains("[Triangulation]")) ? .semibold : .regular, design: .monospaced))
                                .foregroundColor(colorForLine(msg))
                                .textSelection(.enabled)
                                .id(msg.id)
                        }
                    }
                    .padding(8)
                }
                .background(BlueHour.canvasBackground)
                .onChange(of: runner.logLines.count) { _ in
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
    @Binding var showConsole: Bool
    var onLaunch: () -> Void
    @EnvironmentObject var runner: ProcessRunner
    @State private var copiedCLI: Bool = false
    
    private var tickerText: String {
        if runner.status.isRunning, let last = runner.logLines.last?.text.trimmingCharacters(in: .whitespacesAndNewlines), !last.isEmpty {
            return last
        }
        return "▶ \(config.activePhysicsSummaryBadge)"
    }
    
    private var tickerColor: Color {
        switch runner.status {
        case .running: return .green
        case .failed, .terminated: return .red
        case .completed(let code, _): return code == 0 ? .blue : .orange
        case .idle: return Color(red: 0.45, green: 0.90, blue: 0.62)
        }
    }
    
    var body: some View {
        HStack(spacing: 10) {
            // Console drawer toggle
            Button(action: { withAnimation(.easeInOut(duration: 0.15)) { showConsole.toggle() } }) {
                HStack(spacing: 5) {
                    Image(systemName: "terminal")
                    Text("Console")
                    if !runner.logLines.isEmpty {
                        Text("\(runner.logLines.count)")
                            .font(.system(size: 9.5, weight: .bold, design: .monospaced))
                            .padding(.horizontal, 5)
                            .padding(.vertical, 1)
                            .background(BlueHour.pillSurface, in: Capsule())
                    }
                    Image(systemName: showConsole ? "chevron.down" : "chevron.up")
                        .font(.system(size: 9, weight: .bold))
                }
            }
            .controlSize(.regular)
            .help("Show / hide the solver console (⌘J)")
            
            // Live one-line status ticker
            HStack(spacing: 6) {
                Circle().fill(tickerColor).frame(width: 7, height: 7)
                Text(tickerText)
                    .font(.system(size: 10.5, design: .monospaced))
                    .foregroundColor(BlueHour.mist)
                    .lineLimit(1)
                    .truncationMode(.tail)
            }
            .frame(maxWidth: .infinity, alignment: .leading)
            .help(tickerText)
            
            // Reveal in Finder
            Button(action: { runner.openOutputFolder() }) {
                Label("Output", systemImage: "folder")
            }
            .controlSize(.regular)
            .help("Open output folder")
            
            // Copy CLI Command
            Button(action: {
                NSPasteboard.general.clearContents()
                NSPasteboard.general.setString(config.fullCommandLine(), forType: .string)
                copiedCLI = true
                DispatchQueue.main.asyncAfter(deadline: .now() + 1.5) {
                    copiedCLI = false
                }
            }) {
                Label(copiedCLI ? "Copied!" : "CLI", systemImage: copiedCLI ? "checkmark" : "doc.on.doc")
            }
            .controlSize(.regular)
            .help("Copy the equivalent command line:\n" + config.fullCommandLine())
            
            if runner.status.isRunning {
                Button(role: .destructive, action: { runner.stop() }) {
                    Label(runner.isSweepRunning ? "Stop CSL Sweep" : "Stop Simulation", systemImage: "stop.fill")
                        .font(.headline)
                        .padding(.horizontal, 8)
                }
                .buttonStyle(.borderedProminent)
                .tint(.red)
                .controlSize(.large)
            } else {
                Button(action: onLaunch) {
                    let btnTitle: String = config.launchButtonSummaryTitle
                    Label(btnTitle, systemImage: "play.fill")
                        .font(.headline)
                        .padding(.horizontal, 12)
                }
                .buttonStyle(.borderedProminent)
                .tint(config.protocolMode == .sweepMisorientation ? .orange : .accentColor)
                .controlSize(.large)
            }
        }
    }
}
