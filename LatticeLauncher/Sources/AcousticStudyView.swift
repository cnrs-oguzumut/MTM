import SwiftUI
import AppKit

public enum AcousticViewMode: String, CaseIterable, Identifiable {
    case split = "Split View"
    case poincare = "Poincaré Disk"
    case energyPath = "Deformation Energy E(α)"
    case instabilityMode = "Mode Q(n)"
    
    public var id: String { rawValue }
    public var iconName: String {
        switch self {
        case .split: return "rectangle.split.2x1"
        case .poincare: return "circle"
        case .energyPath: return "chart.xyaxis.line"
        case .instabilityMode: return "scope"
        }
    }
}

public struct AcousticStudyView: View {
    @StateObject private var manager = AcousticStudyManager()
    @State private var viewMode: AcousticViewMode = .split
    @State private var showLogs: Bool = false
    
    public init() {}
    
    public var body: some View {
        HSplitView {
            // LEFT SIDEBAR: Scrollable Controls + Pinned Action Footer
            VStack(spacing: 0) {
                ScrollView {
                    VStack(alignment: .leading, spacing: 14) {
                        // Header
                        HStack {
                            Image(systemName: "waveform.path.ecg")
                                .font(.title2)
                                .foregroundColor(.accentColor)
                            VStack(alignment: .leading, spacing: 2) {
                                Text("Parametric Acoustic Study")
                                    .font(.headline)
                                Text("Acoustic tensor stability & Poincaré projection")
                                    .font(.caption2)
                                    .foregroundColor(.secondary)
                            }
                        }
                        .padding(.bottom, 2)
                        
                        // 1. Potential Selection
                        GroupBox(label: Label("Constitutive Potential", systemImage: "atom")) {
                            VStack(alignment: .leading, spacing: 10) {
                                Picker("", selection: $manager.potentialMode) {
                                    ForEach(AcousticPotentialMode.allCases) { mode in
                                        Text(mode.displayName).tag(mode)
                                    }
                                }
                                .labelsHidden()
                                .pickerStyle(.radioGroup)
                                
                                if manager.potentialMode == .zanzotto {
                                    Divider()
                                    
                                    // Beta input
                                    VStack(alignment: .leading, spacing: 4) {
                                        HStack {
                                            Text("β (shear coefficient):")
                                                .font(.caption)
                                            Spacer()
                                            DoubleField(title: "β", value: $manager.beta, width: 80)
                                        }
                                        Slider(value: $manager.beta, in: -1.0...1.0, step: 0.01)
                                    }
                                    
                                    // K input
                                    VStack(alignment: .leading, spacing: 4) {
                                        HStack {
                                            Text("K (bulk modulus):")
                                                .font(.caption)
                                            Spacer()
                                            DoubleField(title: "K", value: $manager.bulkK, width: 80)
                                        }
                                        Slider(value: $manager.bulkK, in: 0.5...20.0, step: 0.25)
                                    }
                                } else if manager.potentialMode == .square {
                                    Text("Mode 4: 2-well Gaussian potential V(r) designed to stabilize square crystal lattice against close-packed collapse.")
                                        .font(.caption2)
                                        .foregroundColor(.secondary)
                                        .padding(.top, 4)
                                } else {
                                    Text(manager.potentialMode == .ljV2 ? 
                                         "Interatomic pair potential V(r) with cutoff r_cut = 2.5, σ = 2·sin(π/10)." : 
                                         "Interatomic pair potential V(r) with cutoff r_cut = 1.866.")
                                        .font(.caption2)
                                        .foregroundColor(.secondary)
                                        .padding(.top, 4)
                                }
                            }
                            .padding(6)
                        }
                        
                        // 2. Lattice Symmetry
                        GroupBox(label: Label("Lattice Symmetry", systemImage: "square.grid.3x3")) {
                            VStack(alignment: .leading, spacing: 8) {
                                Picker("", selection: $manager.latticeSymmetry) {
                                    ForEach(LatticeSymmetry.allCases) { sym in
                                        Text(sym.displayName).tag(sym)
                                    }
                                }
                                .labelsHidden()
                                .pickerStyle(.radioGroup)
                            }
                            .padding(6)
                        }
                        
                        // Deformation Path & Stability Limits
                        GroupBox(label: Label("Deformation Path & Stability Limits", systemImage: "point.topleft.down.curvedto.point.bottomright.up")) {
                            VStack(alignment: .leading, spacing: 8) {
                                Toggle("Show deformation path on Poincaré disk", isOn: $manager.showShearPath)
                                    .font(.caption)
                                
                                Picker("Path Family:", selection: $manager.pathType) {
                                    Text("Simple Shear F(α, θ)").tag("simple_shear")
                                    Text("Pure Shear F□(α)").tag("pure_shear")
                                    Text("Diagonal Shear F◇(α)").tag("diagonal_shear")
                                }
                                .pickerStyle(.menu)
                                .font(.caption)
                                
                                if manager.pathType == "simple_shear" {
                                    VStack(alignment: .leading, spacing: 3) {
                                        HStack {
                                            Text("Slip Angle θ:")
                                                .font(.caption)
                                            Spacer()
                                            Text(String(format: "%.1f°", manager.shearTheta))
                                                .font(.caption.monospaced().bold())
                                        }
                                        Slider(value: $manager.shearTheta, in: 0.0...180.0, step: 1.0)
                                    }
                                }
                                
                                HStack {
                                    Text("α min:").font(.caption)
                                    Spacer()
                                    DoubleField(title: "-1.5", value: $manager.shearAlphaMin, width: 80)
                                }
                                HStack {
                                    Text("α max:").font(.caption)
                                    Spacer()
                                    DoubleField(title: "1.5", value: $manager.shearAlphaMax, width: 80)
                                }
                                HStack {
                                    Text("Step dα:").font(.caption)
                                    Spacer()
                                    DoubleField(title: "0.005", value: $manager.shearDAlpha, width: 80)
                                }
                                
                                if let sp = manager.shearPath, !sp.critical_alphas.isEmpty {
                                    Divider()
                                    VStack(alignment: .leading, spacing: 5) {
                                        HStack {
                                            Text("Stability Limits (α_c):")
                                                .font(.caption.bold())
                                            Spacer()
                                            Text("\(sp.critical_alphas.count) critical")
                                                .font(.system(size: 10, design: .monospaced))
                                                .foregroundColor(.secondary)
                                        }
                                        ForEach(sp.critical_alphas, id: \.self) { ac in
                                            HStack {
                                                Circle().fill(Color.yellow).frame(width: 5, height: 5)
                                                Text(String(format: "α_c = %+.4f", ac))
                                                    .font(.system(size: 10, weight: .semibold, design: .monospaced))
                                                Spacer()
                                                Button("Inspect") {
                                                    manager.hoveredAlpha = ac
                                                    manager.selectedCriticalAlpha = ac
                                                }
                                                .buttonStyle(.borderless)
                                                .font(.system(size: 10))
                                                .foregroundColor(.accentColor)
                                                
                                                Button("Mode") {
                                                    manager.hoveredAlpha = ac
                                                    manager.selectedCriticalAlpha = ac
                                                    viewMode = .instabilityMode
                                                }
                                                .buttonStyle(.borderless)
                                                .font(.system(size: 10, weight: .semibold))
                                                .foregroundColor(.orange)
                                            }
                                        }
                                    }
                                }
                            }
                            .padding(6)
                        }

                        // 3. Elastic Constants, Equilibrium Lattice & Dislocation Solution Card
                        ElasticConstantsCard(constants: manager.elasticConstants, symmetry: manager.latticeSymmetry)
                        
                        // 3B. Poincaré Energy Scale & Percentile Cutoff Card
                        GroupBox(label: Label("Poincaré Energy Scale & Cutoff", systemImage: "paintpalette")) {
                            VStack(alignment: .leading, spacing: 8) {
                                HStack {
                                    Text("Active Field:")
                                        .font(.caption)
                                    Spacer()
                                    Picker("", selection: $manager.colorField) {
                                        ForEach(ColorFieldMode.allCases) { mode in
                                            Text(mode.shortLabel).tag(mode)
                                        }
                                    }
                                    .labelsHidden()
                                    .pickerStyle(.segmented)
                                    .frame(width: 200)
                                }
                                
                                Divider()
                                
                                HStack {
                                    Text("Colormap:")
                                        .font(.caption)
                                    Picker("", selection: $manager.energyColormap) {
                                        ForEach(EnergyColormapPreset.allCases) { cmap in
                                            Text(cmap.rawValue).tag(cmap)
                                        }
                                    }
                                    .labelsHidden()
                                    .frame(width: 110)
                                    Spacer()
                                    Toggle("Invert", isOn: $manager.invertEnergyColormap)
                                        .toggleStyle(.checkbox)
                                        .font(.caption)
                                }
                                
                                VStack(alignment: .leading, spacing: 4) {
                                    HStack {
                                        Text("Upper Percentile Cutoff (P_max):")
                                            .font(.caption)
                                        Spacer()
                                        Text(String(format: "%.1f%%", manager.energyUpperPercentile))
                                            .font(.caption.monospaced().bold())
                                            .foregroundColor(.orange)
                                    }
                                    Slider(value: $manager.energyUpperPercentile, in: 10.0...100.0, step: 0.5)
                                    
                                    ScrollView(.horizontal, showsIndicators: false) {
                                        HStack(spacing: 4) {
                                            ForEach([100.0, 99.0, 95.0, 90.0, 85.0, 75.0, 50.0], id: \.self) { pct in
                                                Button(String(format: "%.0f%%", pct)) {
                                                    manager.useCustomEnergyRange = false
                                                    manager.energyUpperPercentile = pct
                                                }
                                                .buttonStyle(.bordered)
                                                .controlSize(.mini)
                                            }
                                        }
                                    }
                                }
                                
                                VStack(alignment: .leading, spacing: 4) {
                                    HStack {
                                        Text("Lower Percentile (P_min):")
                                            .font(.caption)
                                        Spacer()
                                        Text(String(format: "%.1f%%", manager.energyLowerPercentile))
                                            .font(.caption.monospaced())
                                    }
                                    Slider(value: $manager.energyLowerPercentile, in: 0.0...40.0, step: 0.5)
                                }
                                
                                Toggle("Discard / Hide Values Above Cutoff (Clean Plot)", isOn: $manager.discardAbovePercentile)
                                    .font(.caption.bold())
                                
                                VStack(alignment: .leading, spacing: 4) {
                                    HStack {
                                        Text("Contrast Power (γ):")
                                            .font(.caption)
                                        Spacer()
                                        Text(String(format: "%.2f", manager.energyGamma))
                                            .font(.caption.monospaced().bold())
                                    }
                                    Slider(value: $manager.energyGamma, in: 0.15...2.50, step: 0.05)
                                }
                                
                                let b = manager.effectiveEnergyBounds
                                if b.totalCount > 0 {
                                    Text(String(format: "Effective ΔΦ ∈ [%.4f, %.4f] (%d/%d pts kept)", b.min, b.max, b.keptCount, b.totalCount))
                                        .font(.caption2.monospaced())
                                        .foregroundColor(.secondary)
                                }
                            }
                            .padding(6)
                        }
                        
                        // 4. Sweep Parameters & Refinement Settings
                        GroupBox(label: Label("Parametric Sweep & Refinement", systemImage: "slider.horizontal.3")) {
                            VStack(alignment: .leading, spacing: 10) {
                                HStack {
                                    Text("Max Radius t_max:")
                                        .font(.caption)
                                    Spacer()
                                    Text(String(format: "%.2f", manager.tMax))
                                        .font(.caption.monospaced().bold())
                                }
                                Slider(value: $manager.tMax, in: 0.2...1.5, step: 0.05)
                                
                                VStack(alignment: .leading, spacing: 6) {
                                    Text("Resolution Presets:")
                                        .font(.caption)
                                        .foregroundColor(.secondary)
                                    
                                    ScrollView(.horizontal, showsIndicators: false) {
                                        HStack(spacing: 6) {
                                            Button("Fast") {
                                                manager.dt = 0.04
                                                manager.dp = 0.08
                                            }
                                            .buttonStyle(.bordered)
                                            .controlSize(.small)
                                            
                                            Button("Medium") {
                                                manager.dt = 0.02
                                                manager.dp = 0.05
                                            }
                                            .buttonStyle(.bordered)
                                            .controlSize(.small)
                                            
                                            Button("Fine") {
                                                manager.dt = 0.015
                                                manager.dp = 0.03
                                            }
                                            .buttonStyle(.bordered)
                                            .controlSize(.small)
                                            
                                            Button("Ultra-Fine") {
                                                manager.dt = 0.01
                                                manager.dp = 0.02
                                            }
                                            .buttonStyle(.borderedProminent)
                                            .tint(.purple)
                                            .controlSize(.small)
                                        }
                                    }
                                }
                                
                                Divider()
                                
                                // Fine Refinement Sliders
                                VStack(alignment: .leading, spacing: 6) {
                                    HStack {
                                        Text("Radial Step Δt:")
                                            .font(.caption)
                                        Spacer()
                                        Text(String(format: "%.3f", manager.dt))
                                            .font(.caption.monospaced())
                                    }
                                    Slider(value: $manager.dt, in: 0.005...0.05, step: 0.005)
                                    
                                    HStack {
                                        Text("Angular Step Δp:")
                                            .font(.caption)
                                        Spacer()
                                        Text(String(format: "%.3f", manager.dp))
                                            .font(.caption.monospaced())
                                    }
                                    Slider(value: $manager.dp, in: 0.01...0.10, step: 0.005)
                                    
                                    HStack(spacing: 4) {
                                        Image(systemName: "number")
                                            .font(.caption2)
                                            .foregroundColor(.secondary)
                                        Text("Est. Points: ~\(manager.estimatedPoints) sample nodes")
                                            .font(.caption2.monospaced())
                                            .foregroundColor(.secondary)
                                    }
                                    .padding(.top, 2)
                                }
                            }
                            .padding(6)
                        }
                    }
                    .frame(maxWidth: .infinity, alignment: .leading)
                    .padding(12)
                }
                
                Divider()
                
                // 5. Pinned Sticky Action Footer (Always visible at bottom of sidebar)
                VStack(spacing: 8) {
                    if manager.isRunning {
                        VStack(spacing: 6) {
                            ProgressView(value: manager.progress)
                                .progressViewStyle(.linear)
                            HStack {
                                Text(manager.statusMessage)
                                    .font(.caption)
                                    .foregroundColor(.secondary)
                                    .lineLimit(1)
                                Spacer()
                                Button("Cancel", role: .destructive, action: manager.stopStudy)
                                    .buttonStyle(.borderedProminent)
                                    .tint(.red)
                                    .controlSize(.small)
                            }
                        }
                        .padding(8)
                        .background(Color.accentColor.opacity(0.1), in: RoundedRectangle(cornerRadius: 8))
                    } else {
                        Button(action: manager.runParametricStudy) {
                            HStack {
                                Image(systemName: "play.fill")
                                Text("Run Parametric Acoustic Study")
                                    .fontWeight(.semibold)
                            }
                            .frame(maxWidth: .infinity)
                        }
                        .buttonStyle(.borderedProminent)
                        .controlSize(.large)
                    }
                    
                    HStack(spacing: 8) {
                        Button(action: openDatFileDialog) {
                            Label("Load .dat/.json...", systemImage: "folder")
                                .frame(maxWidth: .infinity)
                        }
                        .buttonStyle(.bordered)
                        .controlSize(.small)
                        
                        Button(action: { showLogs.toggle() }) {
                            Label(showLogs ? "Hide Logs" : "Logs", systemImage: "terminal")
                        }
                        .buttonStyle(.bordered)
                        .controlSize(.small)
                    }
                }
                .padding(.horizontal, 12)
                .padding(.vertical, 10)
                .background(Color(NSColor.controlBackgroundColor))
            }
            .frame(minWidth: 360, idealWidth: 390, maxWidth: 440)
            
            // RIGHT MAIN: Poincaré Stability Map & Console
            VStack(spacing: 0) {
                // Top Header Bar: Pinned View Mode Switcher + Horizontally Scrollable Status & Stats
                HStack(spacing: 10) {
                    // 1. Pinned View Mode Switcher on the left — never moves or clips
                    Picker("", selection: $viewMode) {
                        Label("Poincaré", systemImage: "circle").tag(AcousticViewMode.poincare)
                        Label("Split", systemImage: "rectangle.split.2x1").tag(AcousticViewMode.split)
                        Label("Energy E(α)", systemImage: "chart.xyaxis.line").tag(AcousticViewMode.energyPath)
                        Label("Mode Q(n)", systemImage: "scope").tag(AcousticViewMode.instabilityMode)
                    }
                    .labelsHidden()
                    .pickerStyle(.segmented)
                    .controlSize(.small)
                    .fixedSize()
                    
                    Divider()
                        .frame(height: 18)
                    
                    // 2. Horizontally scrollable stats & subtitle
                    ScrollView(.horizontal, showsIndicators: false) {
                        HStack(spacing: 10) {
                            if !manager.points.isEmpty {
                                HStack(spacing: 5) {
                                    Circle().fill(Color.blue).frame(width: 7, height: 7)
                                    Text("\(manager.summary.stable_points) stable")
                                        .font(.caption.monospaced().bold())
                                    Text("/")
                                        .foregroundColor(.secondary)
                                        .font(.caption2)
                                    Circle().fill(Color.red).frame(width: 7, height: 7)
                                    Text("\(manager.summary.unstable_points) unstable")
                                        .font(.caption.monospaced().bold())
                                    
                                    let pct = manager.summary.total_points > 0 ? (Double(manager.summary.stable_points) / Double(manager.summary.total_points) * 100.0) : 0
                                    Text(String(format: "(%.0f%%)", pct))
                                        .font(.caption2.monospaced())
                                        .foregroundColor(.secondary)
                                }
                                .padding(.horizontal, 8)
                                .padding(.vertical, 4)
                                .background(Color.secondary.opacity(0.12), in: Capsule())
                                .fixedSize()
                            }
                            
                            if let sp = manager.shearPath, !sp.critical_alphas.isEmpty {
                                HStack(spacing: 5) {
                                    Circle().fill(Color.yellow).frame(width: 6, height: 6)
                                    Text("α_c: " + sp.critical_alphas.map { String(format: "%+.4f", $0) }.joined(separator: ", "))
                                        .font(.system(size: 11, weight: .semibold, design: .monospaced))
                                }
                                .padding(.horizontal, 8)
                                .padding(.vertical, 4)
                                .background(Color.yellow.opacity(0.12), in: Capsule())
                                .fixedSize()
                            }
                            
                            Text(manager.lastLoadedFile.isEmpty ?
                                 "Stereographic coordinates x = t·(c₁₁−c₂₂)/2, y = t·c₁₂" :
                                 "File: \(manager.lastLoadedFile)")
                                .font(.caption2)
                                .foregroundColor(.secondary)
                                .lineLimit(1)
                                .fixedSize()
                        }
                        .padding(.vertical, 1)
                    }
                    .frame(maxWidth: .infinity, alignment: .leading)
                    
                    // Export Data Menu
                    if !manager.points.isEmpty {
                        Menu {
                            Button(action: {
                                let csv = manager.exportStabilityBoundariesCSV()
                                ExportDataHelper.saveStringToFile(
                                    content: csv,
                                    defaultFileName: "stability_boundary_points.csv"
                                )
                            }) {
                                Label("Save Stability Boundary Points (CSV)...", systemImage: "circle.dashed")
                            }
                            
                            Button(action: {
                                let csv = manager.exportFullParametricScanCSV()
                                ExportDataHelper.saveStringToFile(
                                    content: csv,
                                    defaultFileName: "acoustic_parametric_scan.csv"
                                )
                            }) {
                                Label("Save Full Stability Scan (CSV)...", systemImage: "tablecells")
                            }
                            
                            if manager.shearPath != nil {
                                Button(action: {
                                    let csv = manager.exportShearPathCSV()
                                    ExportDataHelper.saveStringToFile(
                                        content: csv,
                                        defaultFileName: "shear_path_stability.csv"
                                    )
                                }) {
                                    Label("Save Simple Shear Path (CSV)...", systemImage: "point.topleft.down.curvedto.point.bottomright.up")
                                }
                            }
                            
                            Divider()
                            
                            Button(action: {
                                let csv = manager.exportStabilityBoundariesCSV()
                                ExportDataHelper.copyToClipboard(text: csv)
                            }) {
                                Label("Copy Boundary Points (CSV)", systemImage: "doc.on.doc")
                            }
                            
                            Button(action: {
                                let csv = manager.exportFullParametricScanCSV()
                                ExportDataHelper.copyToClipboard(text: csv)
                            }) {
                                Label("Copy Full Scan (CSV)", systemImage: "doc.on.clipboard")
                            }
                        } label: {
                            Label("Export Data", systemImage: "square.and.arrow.up")
                        }
                        .menuStyle(.borderedButton)
                        .controlSize(.regular)
                    }
                }
                .frame(maxWidth: .infinity, alignment: .leading)
                .padding(.horizontal, 14)
                .padding(.vertical, 8)
                .background(Color(NSColor.controlBackgroundColor))
                
                Divider()
                
                // Main Content Canvas based on viewMode
                Group {
                    switch viewMode {
                    case .poincare:
                        PoincareStabilityMapView(manager: manager)
                    case .energyPath:
                        PathEnergyChartView(manager: manager)
                    case .instabilityMode:
                        AcousticInstabilityModeView(manager: manager)
                    case .split:
                        HSplitView {
                            PoincareStabilityMapView(manager: manager)
                                .frame(minWidth: 340)
                            PathEnergyChartView(manager: manager)
                                .frame(minWidth: 340)
                        }
                    }
                }
                .frame(maxWidth: .infinity, maxHeight: .infinity)
                
                // Collapsible Log Drawer
                if showLogs {
                    Divider()
                    VStack(alignment: .leading, spacing: 4) {
                        HStack {
                            Text("Execution Log")
                                .font(.caption.bold())
                            Spacer()
                            Button("Clear") { manager.logLines.removeAll() }
                                .font(.caption2)
                                .buttonStyle(.borderless)
                        }
                        .padding(.horizontal, 10)
                        .padding(.top, 6)
                        
                        ScrollViewReader { proxy in
                            ScrollView {
                                LazyVStack(alignment: .leading, spacing: 2) {
                                    ForEach(Array(manager.logLines.enumerated()), id: \.offset) { item in
                                        Text(item.element)
                                            .font(.system(size: 11, design: .monospaced))
                                            .foregroundColor(.secondary)
                                            .id(item.offset)
                                    }
                                }
                                .padding(8)
                            }
                            .onChange(of: manager.logLines.count) { count in
                                if count > 0 {
                                    proxy.scrollTo(count - 1, anchor: .bottom)
                                }
                            }
                        }
                    }
                    .frame(height: 160)
                    .background(Color(NSColor.textBackgroundColor))
                }
            }
        }
    }
    
    private func openDatFileDialog() {
        let panel = NSOpenPanel()
        panel.canChooseFiles = true
        panel.canChooseDirectories = false
        panel.allowsMultipleSelection = false
        panel.allowedContentTypes = []
        
        if panel.runModal() == .OK, let url = panel.url {
            if url.pathExtension.lowercased() == "json" {
                manager.loadResultsFromJSON(path: url.path)
            } else {
                manager.loadDatFile(url: url)
            }
        }
    }
}

// MARK: - Elastic Constants Card
struct ElasticConstantsCard: View {
    let constants: AcousticElasticConstantsModel
    let symmetry: LatticeSymmetry
    
    var body: some View {
        GroupBox(label: Label("Elastic Constants & Dislocation Parameters", systemImage: "sparkles")) {
            VStack(alignment: .leading, spacing: 10) {
                // Equilibrium Lattice Parameter a0 (b)
                HStack {
                    Label("Equilibrium Lattice Spacing a₀ (b):", systemImage: "ruler")
                        .font(.caption2.bold())
                        .foregroundColor(.secondary)
                    Spacer()
                    Text(String(format: "%.9f", constants.equilibrium_lattice_parameter))
                        .font(.system(size: 11, weight: .bold, design: .monospaced))
                        .foregroundColor(.accentColor)
                }
                .padding(.horizontal, 8)
                .padding(.vertical, 5)
                .background(Color.accentColor.opacity(0.10), in: RoundedRectangle(cornerRadius: 6))
                
                // Moduli Grid
                Grid(alignment: .leading, horizontalSpacing: 12, verticalSpacing: 6) {
                    GridRow {
                        VStack(alignment: .leading, spacing: 1) {
                            Text("C₁₁").font(.caption2).foregroundColor(.secondary)
                            Text(String(format: "%.3f", constants.C11))
                                .font(.caption.monospaced().bold())
                        }
                        VStack(alignment: .leading, spacing: 1) {
                            Text("C₂₂").font(.caption2).foregroundColor(.secondary)
                            Text(String(format: "%.3f", constants.C22))
                                .font(.caption.monospaced().bold())
                        }
                        VStack(alignment: .leading, spacing: 1) {
                            Text("C₁₂").font(.caption2).foregroundColor(.secondary)
                            Text(String(format: "%.3f", constants.C12))
                                .font(.caption.monospaced().bold())
                        }
                        VStack(alignment: .leading, spacing: 1) {
                            Text("C₆₆ (Shear)").font(.caption2).foregroundColor(.secondary)
                            Text(String(format: "%.3f", constants.C66))
                                .font(.caption.monospaced().bold())
                        }
                    }
                    GridRow {
                        VStack(alignment: .leading, spacing: 1) {
                            Text("Bulk K_2D").font(.caption2).foregroundColor(.secondary)
                            Text(String(format: "%.3f", constants.bulk_modulus))
                                .font(.caption.monospaced())
                        }
                        VStack(alignment: .leading, spacing: 1) {
                            Text("Shear μ").font(.caption2).foregroundColor(.secondary)
                            Text(String(format: "%.3f", constants.shear_modulus))
                                .font(.caption.monospaced())
                        }
                        VStack(alignment: .leading, spacing: 1) {
                            Text("Born Stable").font(.caption2).foregroundColor(.secondary)
                            Text(constants.born_stable ? "STABLE" : "UNSTABLE")
                                .font(.caption2.bold())
                                .foregroundColor(constants.born_stable ? .green : .red)
                        }
                        VStack(alignment: .leading, spacing: 1) {
                            Text("Zener Ratio A").font(.caption2).foregroundColor(.secondary)
                            Text(String(format: "%.3f", constants.Zener_A))
                                .font(.caption.monospaced().bold())
                                .foregroundColor(constants.is_isotropic ? .blue : (constants.Zener_A > 1 ? .purple : .orange))
                        }
                    }
                }
                .padding(6)
                .background(Color.secondary.opacity(0.08), in: RoundedRectangle(cornerRadius: 6))
                
                // Dynamic Dislocation Solution Parameters (Volterra vs Stroh)
                VStack(alignment: .leading, spacing: 4) {
                    if symmetry == .triangular {
                        HStack {
                            Text("Triangular Dislocation: Isotropic 2D Volterra (Eq. 15)")
                                .font(.caption2.bold())
                                .foregroundColor(.orange)
                            Spacer()
                        }
                        HStack(spacing: 14) {
                            VStack(alignment: .leading, spacing: 1) {
                                Text("Poisson ν = C₁₂/C₁₁").font(.caption2).foregroundColor(.secondary)
                                Text(String(format: "%.4f", constants.poisson_ratio))
                                    .font(.caption.monospaced().bold())
                            }
                            VStack(alignment: .leading, spacing: 1) {
                                Text("(1−ν)/(2π)").font(.caption2).foregroundColor(.secondary)
                                Text(String(format: "%.4f", (1.0 - constants.poisson_ratio) / (2.0 * .pi)))
                                    .font(.caption.monospaced())
                            }
                            VStack(alignment: .leading, spacing: 1) {
                                Text("Factor K = μ(1+ν)").font(.caption2).foregroundColor(.secondary)
                                Text(String(format: "%.3f", constants.C66 * (1.0 + constants.poisson_ratio)))
                                    .font(.caption.monospaced().bold())
                            }
                        }
                    } else {
                        HStack {
                            Text("Square Dislocation: Anisotropic Stroh (Eqs. A5–A16)")
                                .font(.caption2.bold())
                                .foregroundColor(.blue)
                            Spacer()
                            Text(constants.stroh_kappa > 1.0001 ? "κ > 1" : (constants.stroh_kappa < 0.9999 ? "κ < 1" : "κ = 1"))
                                .font(.caption2.monospaced().bold())
                                .foregroundColor(.secondary)
                        }
                        Grid(alignment: .leading, horizontalSpacing: 12, verticalSpacing: 4) {
                            GridRow {
                                VStack(alignment: .leading, spacing: 1) {
                                    Text("κ").font(.caption2).foregroundColor(.secondary)
                                    Text(String(format: "%.4f", constants.stroh_kappa)).font(.caption.monospaced().bold())
                                }
                                VStack(alignment: .leading, spacing: 1) {
                                    Text("ξ").font(.caption2).foregroundColor(.secondary)
                                    Text(String(format: "%.4f", constants.stroh_xi)).font(.caption.monospaced())
                                }
                                VStack(alignment: .leading, spacing: 1) {
                                    Text("ζ").font(.caption2).foregroundColor(.secondary)
                                    Text(String(format: "%.4f", constants.stroh_zeta)).font(.caption.monospaced())
                                }
                                VStack(alignment: .leading, spacing: 1) {
                                    Text("Foreman K").font(.caption2).foregroundColor(.secondary)
                                    Text(String(format: "%.3f", constants.foreman_K)).font(.caption.monospaced().bold())
                                }
                            }
                            GridRow {
                                VStack(alignment: .leading, spacing: 1) {
                                    Text("g₁").font(.caption2).foregroundColor(.secondary)
                                    Text(String(format: "%.4f", constants.stroh_g1)).font(.caption.monospaced())
                                }
                                VStack(alignment: .leading, spacing: 1) {
                                    Text("g₂").font(.caption2).foregroundColor(.secondary)
                                    Text(String(format: "%.4f", constants.stroh_g2)).font(.caption.monospaced())
                                }
                                VStack(alignment: .leading, spacing: 1) {
                                    Text("g₃").font(.caption2).foregroundColor(.secondary)
                                    Text(String(format: "%.4f", constants.stroh_g3)).font(.caption.monospaced())
                                }
                                VStack(alignment: .leading, spacing: 1) {
                                    Text("ν = C₁₂/C₁₁").font(.caption2).foregroundColor(.secondary)
                                    Text(String(format: "%.4f", constants.poisson_ratio)).font(.caption.monospaced())
                                }
                            }
                        }
                    }
                }
                .padding(6)
                .background(Color.secondary.opacity(0.06), in: RoundedRectangle(cornerRadius: 6))
                
                // Anisotropy Badge
                HStack(spacing: 8) {
                    if constants.is_isotropic {
                        Label("Isotropic (A ≈ 1.0)", systemImage: "circle")
                            .font(.caption2.bold())
                            .foregroundColor(.blue)
                    } else if constants.Zener_A > 1.0 {
                        Label("Anisotropic: <10> Shear Soft", systemImage: "arrow.up.and.down.and.arrow.left.and.right")
                            .font(.caption2.bold())
                            .foregroundColor(.purple)
                    } else {
                        Label("Anisotropic: <11> Diagonal Soft", systemImage: "arrow.up.left.and.down.right.and.arrow.up.right.and.down.left")
                            .font(.caption2.bold())
                            .foregroundColor(.orange)
                    }
                    
                    Spacer()
                    
                    Image(systemName: constants.born_stable ? "checkmark.seal.fill" : "exclamationmark.triangle.fill")
                        .foregroundColor(constants.born_stable ? .green : .red)
                }
                
                // Commentary Box
                Text(constants.anisotropy_comment)
                    .font(.system(size: 11))
                    .foregroundColor(.secondary)
                    .lineLimit(4)
                    .fixedSize(horizontal: false, vertical: true)
                    .padding(8)
                    .background(Color(NSColor.textBackgroundColor).opacity(0.6), in: RoundedRectangle(cornerRadius: 6))
            }
            .padding(6)
        }
    }
}

struct BadgeItem: View {
    let title: String
    let value: String
    let color: Color
    
    var body: some View {
        HStack(spacing: 4) {
            Text(title + ":")
                .font(.caption2)
                .foregroundColor(.secondary)
            Text(value)
                .font(.caption2.bold())
                .foregroundColor(color)
        }
        .padding(.horizontal, 8)
        .padding(.vertical, 4)
        .background(color.opacity(0.12), in: Capsule())
    }
}
