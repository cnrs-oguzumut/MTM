import SwiftUI
import AppKit

enum MainDisplayTab: String, CaseIterable, Identifiable {
    case curve = "Stress–Strain Curve"
    case terminal = "Terminal Log"
    
    var id: String { rawValue }
    var icon: String {
        switch self {
        case .curve: return "chart.xyaxis.line"
        case .terminal: return "terminal"
        }
    }
}

struct ContentView: View {
    @ObservedObject var config: SimulationConfig
    @EnvironmentObject var runner: ProcessRunner
    @State private var displayTab: MainDisplayTab = .curve
    @State private var filterText: String = ""
    @State private var autoScroll: Bool = true
    
    var body: some View {
        HSplitView {
            // LEFT: Configuration Sidebar (Width ~420px)
            SidebarControlsView(config: config)
                .frame(minWidth: 380, idealWidth: 420, maxWidth: 480)
            
            // RIGHT: Monitor & Live Terminal / Curve
            VStack(spacing: 0) {
                // Top Header: Live Status & Metrics
                LiveStatusHeaderView()
                    .padding(.horizontal, 16)
                    .padding(.top, 14)
                    .padding(.bottom, 10)
                    .background(Color(NSColor.controlBackgroundColor))
                
                Divider()
                
                // View Selector Bar (Curve vs Terminal)
                HStack {
                    Picker("View", selection: $displayTab) {
                        ForEach(MainDisplayTab.allCases) { tab in
                            Label(tab.rawValue, systemImage: tab.icon).tag(tab)
                        }
                    }
                    .pickerStyle(.segmented)
                    .frame(width: 320)
                    
                    Spacer()
                    
                    // Quick stats tag
                    if !runner.curvePoints.isEmpty {
                        Text("\(runner.curvePoints.count) data points plotted")
                            .font(.caption)
                            .foregroundColor(.secondary)
                    }
                }
                .padding(.horizontal, 16)
                .padding(.vertical, 8)
                .background(Color(NSColor.controlBackgroundColor))
                
                Divider()
                
                // Command Preview Banner
                CommandPreviewBar(config: config)
                
                Divider()
                
                // Main Content View: Live Curve or Terminal
                if displayTab == .curve {
                    StressStrainChartView()
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
            .frame(minWidth: 550)
        }
        .onAppear {
            runner.loadExistingIfAvailable(dir: config.targetDirectory)
        }
    }
}

// MARK: - Sidebar Controls
struct SidebarControlsView: View {
    @ObservedObject var config: SimulationConfig
    
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
                                TextField("100", value: $config.nx, formatter: NumberFormatter())
                                    .textFieldStyle(.roundedBorder)
                            }
                            VStack(alignment: .leading, spacing: 4) {
                                Text("NY (Rows)")
                                    .font(.caption)
                                    .foregroundColor(.secondary)
                                TextField("100", value: $config.ny, formatter: NumberFormatter())
                                    .textFieldStyle(.roundedBorder)
                            }
                        }
                        
                        // Preset buttons
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
                        
                        HStack {
                            Text("Random Seed:")
                                .font(.subheadline)
                            TextField("Seed", value: $config.seed, formatter: NumberFormatter())
                                .textFieldStyle(.roundedBorder)
                                .frame(width: 90)
                            Button(action: { config.randomizeSeed() }) {
                                Label("Random", systemImage: "dice")
                            }
                            .buttonStyle(.bordered)
                            .controlSize(.small)
                        }
                        
                        Toggle("Adaptive Remeshing (Delaunay Surgery)", isOn: $config.enableRemeshing)
                            .toggleStyle(.switch)
                    }
                    .padding(8)
                }
                
                // 3. Preconditioning & Energy Relaxation
                GroupBox(label: Label("L-BFGS Preconditioner", systemImage: "bolt.fill")) {
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
                                TextField("1", value: $config.precondFromStep, formatter: NumberFormatter())
                                    .textFieldStyle(.roundedBorder)
                                    .frame(width: 70)
                            }
                            HStack {
                                Text("Refresh every N steps (0=adaptive):")
                                    .font(.subheadline)
                                Spacer()
                                TextField("0", value: $config.precondRefresh, formatter: NumberFormatter())
                                    .textFieldStyle(.roundedBorder)
                                    .frame(width: 70)
                            }
                        }
                    }
                    .padding(8)
                }
                
                // 4. Loading Schedule
                GroupBox(label: Label("Loading Schedule (|α|)", systemImage: "chart.line.uptrend.xyaxis")) {
                    VStack(alignment: .leading, spacing: 8) {
                        HStack {
                            Text("Start strain |α_start|:")
                                .font(.subheadline)
                            Spacer()
                            TextField("0.14", value: $config.alphaStart, format: .number.precision(.fractionLength(4)))
                                .textFieldStyle(.roundedBorder)
                                .frame(width: 80)
                        }
                        HStack {
                            Text("End strain |α_end|:")
                                .font(.subheadline)
                            Spacer()
                            TextField("1.0", value: $config.alphaEnd, format: .number.precision(.fractionLength(2)))
                                .textFieldStyle(.roundedBorder)
                                .frame(width: 80)
                        }
                        HStack {
                            Text("Step size |dα|:")
                                .font(.subheadline)
                            Spacer()
                            TextField("0.00006", value: $config.stepSize, format: .number.precision(.fractionLength(6)))
                                .textFieldStyle(.roundedBorder)
                                .frame(width: 90)
                        }
                    }
                    .padding(8)
                }
                
                // 5. Stability Monitor & Eigenvalues
                GroupBox(label: Label("Stability & Eigenvalues", systemImage: "waveform.path.ecg")) {
                    VStack(alignment: .leading, spacing: 10) {
                        Toggle("Compute Lowest Stiffness Eigenvalues", isOn: $config.enableEigenvalues)
                            .toggleStyle(.switch)
                        
                        if config.enableEigenvalues {
                            HStack {
                                Text("Every N steps:")
                                Spacer()
                                TextField("5", value: $config.eigEvery, formatter: NumberFormatter())
                                    .textFieldStyle(.roundedBorder)
                                    .frame(width: 60)
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
                
                // 6. Data Saving & Output Controls
                GroupBox(label: Label("Output & Checkpoint Controls", systemImage: "externaldrive")) {
                    VStack(alignment: .leading, spacing: 8) {
                        Toggle("Inner iterations trace (--trace=1)", isOn: $config.traceInnerIterations)
                        Toggle("Save Surgery VTK Snapshots (--surgery-vtk)", isOn: $config.saveSurgeryVTK)
                        Toggle("Save Legacy triangle_data", isOn: $config.saveTriangleData)
                        
                        HStack {
                            Text("Checkpoint every N steps:")
                            Spacer()
                            TextField("500", value: $config.checkpointInterval, formatter: NumberFormatter())
                                .textFieldStyle(.roundedBorder)
                                .frame(width: 70)
                        }
                        HStack {
                            Text("Stress drop threshold Δσ/σ:")
                            Spacer()
                            TextField("0.10", value: $config.stressDropThreshold, format: .number.precision(.fractionLength(2)))
                                .textFieldStyle(.roundedBorder)
                                .frame(width: 70)
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
            .padding(14)
        }
        .background(Color(NSColor.controlBackgroundColor).opacity(0.5))
    }
}

// MARK: - Live Status Header
struct LiveStatusHeaderView: View {
    @EnvironmentObject var runner: ProcessRunner
    
    var body: some View {
        HStack(spacing: 16) {
            // Status Badge
            HStack(spacing: 6) {
                Circle()
                    .fill(statusColor)
                    .frame(width: 10, height: 10)
                Text(runner.status.label)
                    .font(.headline)
            }
            
            Spacer()
            
            // Metrics Cards
            MetricPill(label: "Step", value: "\(runner.currentStep)")
            MetricPill(label: "Strain |α|", value: String(format: "%.4f", runner.currentAlpha))
            MetricPill(label: "Stress σ", value: String(format: "%.4f", runner.currentStress))
            MetricPill(label: "Avalanches", value: "\(runner.avalancheCount)")
            MetricPill(label: "Elapsed", value: formatTime(runner.elapsedTime))
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
    
    var body: some View {
        VStack(spacing: 2) {
            Text(label)
                .font(.system(size: 10, weight: .semibold))
                .foregroundColor(.secondary)
            Text(value)
                .font(.system(size: 13, weight: .bold, design: .monospaced))
        }
        .padding(.horizontal, 10)
        .padding(.vertical, 5)
        .background(Color(NSColor.quaternaryLabelColor))
        .cornerRadius(6)
    }
}

// MARK: - Command Preview Bar
struct CommandPreviewBar: View {
    @ObservedObject var config: SimulationConfig
    @State private var copied: Bool = false
    
    var body: some View {
        HStack {
            Image(systemName: "terminal")
                .foregroundColor(.secondary)
            Text(config.fullCommandLine())
                .font(.system(size: 11, design: .monospaced))
                .lineLimit(1)
                .truncationMode(.middle)
                .foregroundColor(.primary)
            Spacer()
            Button(action: {
                NSPasteboard.general.clearContents()
                NSPasteboard.general.setString(config.fullCommandLine(), forType: .string)
                copied = true
                DispatchQueue.main.asyncAfter(deadline: .now() + 1.5) {
                    copied = false
                }
            }) {
                Label(copied ? "Copied!" : "Copy", systemImage: copied ? "checkmark" : "doc.on.doc")
            }
            .buttonStyle(.borderless)
            .controlSize(.small)
        }
        .padding(.horizontal, 16)
        .padding(.vertical, 6)
        .background(Color(NSColor.textBackgroundColor).opacity(0.4))
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
    
    var body: some View {
        HStack(spacing: 12) {
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
