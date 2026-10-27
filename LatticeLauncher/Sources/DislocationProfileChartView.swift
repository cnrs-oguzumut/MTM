import SwiftUI
import Charts
import AppKit

public enum DislocationPlotMode: String, CaseIterable, Identifiable {
    case radialEnergy = "E(r) Radial"
    case slipEnergy = "E(x–x_c) Slip Line"
    case slipStress = "σ(x–x_c) Stresses"
    case radialStress = "σ(r) Stresses"
    
    public var id: String { rawValue }
    
    public var icon: String {
        switch self {
        case .radialEnergy: return "circle.dashed"
        case .slipEnergy: return "waveform.path.ecg"
        case .slipStress: return "arrow.left.and.right"
        case .radialStress: return "chart.bar.xaxis"
        }
    }
    
    public var filePrefix: String {
        switch self {
        case .radialEnergy: return "radial_energy"
        case .slipEnergy: return "slip_line_energy"
        case .slipStress: return "slip_line_stresses"
        case .radialStress: return "radial_stresses"
        }
    }
}

public enum DislocationStateSelection: String, CaseIterable, Identifiable {
    case relaxed = "Relaxed"
    case unrelaxed = "Unrelaxed"
    case compare = "Compare"
    
    public var id: String { rawValue }
}

public enum RadialXAxisScale: String, CaseIterable, Identifiable {
    case linear = "Linear (r/b)"
    case log = "Log ln(r/b)"
    
    public var id: String { rawValue }
}

public enum SlipRowDisplay: String, CaseIterable, Identifiable {
    case average = "Average Plane"
    case upper = "Upper Row"
    case lower = "Lower Row"
    case bothRows = "Both Rows"
    
    public var id: String { rawValue }
    
    public var shortTitle: String {
        switch self {
        case .average: return "Avg Plane"
        case .upper: return "Upper"
        case .lower: return "Lower"
        case .bothRows: return "Both"
        }
    }
}

public enum SlipXRange: String, CaseIterable, Identifiable {
    case full = "Full Crystal"
    case rFree = "Within R_free"
    case nearCore = "Near Core"
    
    public var id: String { rawValue }
    
    public var shortTitle: String {
        switch self {
        case .full: return "Full Crystal"
        case .rFree: return "R_free"
        case .nearCore: return "Near Core"
        }
    }
}

public struct DislocationProfileChartView: View {
    @EnvironmentObject var runner: ProcessRunner
    
    @State private var plotMode: DislocationPlotMode = .radialEnergy
    @State private var stateSelection: DislocationStateSelection = .relaxed
    @State private var radialScale: RadialXAxisScale = .linear
    @State private var slipRow: SlipRowDisplay = .average
    @State private var slipRange: SlipXRange = .full
    @State private var nearCoreSpan: Double = 5.0
    
    // Stress component toggles
    @State private var showSxx: Bool = true
    @State private var showSyy: Bool = true
    @State private var showSxy: Bool = true
    @State private var showSVM: Bool = false
    
    // Crosshair hover state
    @State private var hoveredX: Double?
    
    private var data: DislocationProfileJSONRoot? {
        runner.dislocationData
    }
    
    public var body: some View {
        VStack(spacing: 0) {
            if let root = data {
                // 1. Controls Header
                controlsHeader
                
                Divider()
                
                // 2. KPI Summary Cards
                kpiBar(root: root)
                    .padding(.horizontal, 16)
                    .padding(.vertical, 8)
                    .background(BlueHour.controlSurface.opacity(0.5))
                
                Divider()
                
                // 3. Chart Area
                chartView(root: root)
                    .padding(16)
                    .frame(maxWidth: .infinity, maxHeight: .infinity)
            } else {
                emptyStateView
            }
        }
        .background(BlueHour.canvasBackground)
        .onAppear {
            runner.pollDislocationProfiles()
        }
    }
    
    // MARK: - Controls Header
    private var controlsHeader: some View {
        VStack(spacing: 8) {
            tier1Header
            tier2Header
        }
        .padding(.horizontal, 16)
        .padding(.vertical, 8)
        .background(BlueHour.windowBackground)
    }
    
    // MARK: - Tier 1 Header
    private var tier1Header: some View {
        HStack(spacing: 10) {
            ScrollView(.horizontal, showsIndicators: false) {
                Picker("", selection: $plotMode) {
                    ForEach(DislocationPlotMode.allCases) { mode in
                        Label(mode.rawValue, systemImage: mode.icon).tag(mode)
                    }
                }
                .pickerStyle(.segmented)
                .labelsHidden()
                .fixedSize()
            }
            .frame(maxWidth: .infinity, alignment: .leading)
            
            Picker("", selection: $stateSelection) {
                ForEach(DislocationStateSelection.allCases) { state in
                    Text(state.rawValue).tag(state)
                }
            }
            .pickerStyle(.segmented)
            .labelsHidden()
            .fixedSize()
            
            exportMenu
        }
    }
    
    // MARK: - Export Menu
    private var exportMenu: some View {
        Menu {
            Button(action: exportActivePlotCSV) {
                Label("Save Active Plot Data (CSV)...", systemImage: "chart.xyaxis.line")
            }
            Button(action: exportAllProfilesCSV) {
                Label("Save All Profiles (Full Crystal CSV)...", systemImage: "tablecells")
            }
            Button(action: exportAllProfilesJSON) {
                Label("Save Full Dataset (JSON)...", systemImage: "curlybraces")
            }
            Divider()
            Button(action: copyActivePlotCSVToClipboard) {
                Label("Copy Active Plot (CSV)", systemImage: "doc.on.doc")
            }
            Button(action: copyAllProfilesCSVToClipboard) {
                Label("Copy All Profiles (CSV)", systemImage: "doc.on.clipboard")
            }
        } label: {
            HStack(spacing: 4) {
                Image(systemName: "square.and.arrow.up")
                Text("Export")
            }
            .font(.caption.weight(.medium))
            .padding(.horizontal, 8)
            .padding(.vertical, 4)
            .background(BlueHour.controlSurface)
            .cornerRadius(5)
            .overlay(
                RoundedRectangle(cornerRadius: 5)
                    .stroke(Color.secondary.opacity(0.35), lineWidth: 1)
            )
        }
        .menuStyle(.borderlessButton)
        .fixedSize()
        .help("Export Data (CSV / JSON)")
    }
    
    // MARK: - Tier 2 Header
    @ViewBuilder
    private var tier2Header: some View {
        switch plotMode {
        case .radialEnergy:
            radialEnergyOptions
        case .slipEnergy:
            slipEnergyOptions
        case .slipStress:
            slipStressOptions
        case .radialStress:
            radialStressOptions
        }
    }
    
    private var radialEnergyOptions: some View {
        ScrollView(.horizontal, showsIndicators: false) {
            HStack(spacing: 8) {
                Text("Scale:")
                    .font(.caption.bold())
                    .foregroundColor(.secondary)
                    .lineLimit(1)
                    .fixedSize()
                Picker("", selection: $radialScale) {
                    ForEach(RadialXAxisScale.allCases) { s in
                        Text(s.rawValue).tag(s)
                    }
                }
                .pickerStyle(.segmented)
                .labelsHidden()
                .frame(width: 190)
            }
            .font(.caption)
        }
        .frame(maxWidth: .infinity, alignment: .leading)
    }
    
    private var slipEnergyOptions: some View {
        ScrollView(.horizontal, showsIndicators: false) {
            HStack(spacing: 12) {
                HStack(spacing: 5) {
                    Text("Atoms:")
                        .font(.caption.bold())
                        .foregroundColor(.secondary)
                        .lineLimit(1)
                        .fixedSize()
                    Picker("", selection: $slipRow) {
                        ForEach(SlipRowDisplay.allCases) { r in
                            Text(r.shortTitle).tag(r)
                        }
                    }
                    .pickerStyle(.segmented)
                    .labelsHidden()
                    .frame(width: 220)
                }
                
                Divider().frame(height: 14)
                
                HStack(spacing: 5) {
                    Text("Domain:")
                        .font(.caption.bold())
                        .foregroundColor(.secondary)
                        .lineLimit(1)
                        .fixedSize()
                    Picker("", selection: $slipRange) {
                        ForEach(SlipXRange.allCases) { r in
                            Text(r.shortTitle).tag(r)
                        }
                    }
                    .pickerStyle(.segmented)
                    .labelsHidden()
                    .frame(width: 200)
                }
                
                if slipRange == .nearCore {
                    Divider().frame(height: 14)
                    
                    HStack(spacing: 5) {
                        Text("Span:")
                            .font(.caption.bold())
                            .foregroundColor(.secondary)
                            .lineLimit(1)
                            .fixedSize()
                        
                        nearCoreSpanSelector
                    }
                }
            }
            .font(.caption)
        }
        .frame(maxWidth: .infinity, alignment: .leading)
    }
    
    private var slipStressOptions: some View {
        ScrollView(.horizontal, showsIndicators: false) {
            HStack(spacing: 12) {
                HStack(spacing: 5) {
                    Text("Stress:")
                        .font(.caption.bold())
                        .foregroundColor(.secondary)
                        .lineLimit(1)
                        .fixedSize()
                    Toggle("σ_xx", isOn: $showSxx).toggleStyle(.checkbox).fixedSize()
                    Toggle("σ_yy", isOn: $showSyy).toggleStyle(.checkbox).fixedSize()
                    Toggle("σ_xy", isOn: $showSxy).toggleStyle(.checkbox).fixedSize()
                    Toggle("VM", isOn: $showSVM).toggleStyle(.checkbox).fixedSize()
                }
                
                Divider().frame(height: 14)
                
                HStack(spacing: 5) {
                    Text("Plane:")
                        .font(.caption.bold())
                        .foregroundColor(.secondary)
                        .lineLimit(1)
                        .fixedSize()
                    Picker("", selection: $slipRow) {
                        Text("Avg Plane").tag(SlipRowDisplay.average)
                        Text("Upper").tag(SlipRowDisplay.upper)
                        Text("Lower").tag(SlipRowDisplay.lower)
                    }
                    .pickerStyle(.segmented)
                    .labelsHidden()
                    .frame(width: 175)
                }
                
                Divider().frame(height: 14)
                
                HStack(spacing: 5) {
                    Text("Domain:")
                        .font(.caption.bold())
                        .foregroundColor(.secondary)
                        .lineLimit(1)
                        .fixedSize()
                    Picker("", selection: $slipRange) {
                        ForEach(SlipXRange.allCases) { r in
                            Text(r.shortTitle).tag(r)
                        }
                    }
                    .pickerStyle(.segmented)
                    .labelsHidden()
                    .frame(width: 200)
                }
                
                if slipRange == .nearCore {
                    Divider().frame(height: 14)
                    
                    HStack(spacing: 5) {
                        Text("Span:")
                            .font(.caption.bold())
                            .foregroundColor(.secondary)
                            .lineLimit(1)
                            .fixedSize()
                        
                        nearCoreSpanSelector
                    }
                }
            }
            .font(.caption)
        }
        .frame(maxWidth: .infinity, alignment: .leading)
    }
    
    private var radialStressOptions: some View {
        ScrollView(.horizontal, showsIndicators: false) {
            HStack(spacing: 12) {
                HStack(spacing: 6) {
                    Text("Stress:")
                        .font(.caption.bold())
                        .foregroundColor(.secondary)
                        .lineLimit(1)
                        .fixedSize()
                    Toggle("σ_xx", isOn: $showSxx).toggleStyle(.checkbox).fixedSize()
                    Toggle("σ_yy", isOn: $showSyy).toggleStyle(.checkbox).fixedSize()
                    Toggle("σ_xy", isOn: $showSxy).toggleStyle(.checkbox).fixedSize()
                    Toggle("VM", isOn: $showSVM).toggleStyle(.checkbox).fixedSize()
                }
                
                Divider().frame(height: 14)
                
                HStack(spacing: 6) {
                    Text("Scale:")
                        .font(.caption.bold())
                        .foregroundColor(.secondary)
                        .lineLimit(1)
                        .fixedSize()
                    Picker("", selection: $radialScale) {
                        ForEach(RadialXAxisScale.allCases) { s in
                            Text(s.rawValue).tag(s)
                        }
                    }
                    .pickerStyle(.segmented)
                    .labelsHidden()
                    .frame(width: 190)
                }
            }
            .font(.caption)
        }
        .frame(maxWidth: .infinity, alignment: .leading)
    }
    
    private var nearCoreSpanSelector: some View {
        HStack(spacing: 4) {
            Menu {
                ForEach([2, 3, 5, 8, 10, 15, 20, 30, 50], id: \.self) { val in
                    Button("±\(val)b") { nearCoreSpan = Double(val) }
                }
            } label: {
                HStack(spacing: 3) {
                    Text("±\(Int(nearCoreSpan))b")
                        .font(.caption.monospacedDigit().bold())
                    Image(systemName: "chevron.up.chevron.down")
                        .font(.system(size: 8))
                }
                .padding(.horizontal, 6)
                .padding(.vertical, 3)
                .background(BlueHour.controlSurface)
                .cornerRadius(4)
                .overlay(
                    RoundedRectangle(cornerRadius: 4)
                        .stroke(Color.secondary.opacity(0.35), lineWidth: 1)
                )
            }
            .menuStyle(.borderlessButton)
            .fixedSize()
            
            Stepper("", value: $nearCoreSpan, in: 1.0...100.0, step: 1.0)
                .labelsHidden()
                .controlSize(.mini)
                .fixedSize()
        }
    }
    
    // MARK: - KPI Bar
    private func kpiBar(root: DislocationProfileJSONRoot) -> some View {
        let activeState = (stateSelection == .unrelaxed) ? root.unrelaxed : (root.relaxed ?? root.unrelaxed)
        let coreFit = activeState.coreFit
        let eCore = coreFit.coreEnergy
        let kFactor = coreFit.elasticK
        let r2 = coreFit.r2
        
        let unrelCore = root.unrelaxed.coreFit.coreEnergy
        let relCore = root.relaxed?.coreFit.coreEnergy ?? unrelCore
        let drop = unrelCore - relCore
        // Macro Stress KPI
        let maxMacroStress: Double = {
            if let ms = root.finalMacroStress {
                return max(abs(ms.xx), max(abs(ms.yy), abs(ms.xy)))
            }
            return abs(activeState.totalStress)
        }()
        let isZeroStress = root.macroStressRelaxed

        return ScrollView(.horizontal, showsIndicators: false) {
            HStack(spacing: 10) {
                DislocKpiPill(
                    title: "Core Position (x_c, y_c)",
                    value: String(format: "(%.2f, %.2f)", root.coreX, root.coreY),
                    icon: "scope",
                    color: .blue
                )
                
                DislocKpiPill(
                    title: "Fitted Core Energy E_core",
                    value: String(format: "%.4f", eCore),
                    subtitle: String(format: "R² = %.3f", r2),
                    icon: "atom",
                    color: .purple
                )
                
                DislocKpiPill(
                    title: "Elastic Factor K",
                    value: String(format: "%.4f", kFactor),
                    subtitle: "μb²/[4π(1-ν)]",
                    icon: "function",
                    color: .orange
                )
                
                DislocKpiPill(
                    title: "Core Relaxation Drop ΔE",
                    value: String(format: "%.4f", drop),
                    subtitle: String(format: "%.1f%% drop", unrelCore != 0 ? (drop / unrelCore) * 100 : 0),
                    icon: "arrow.down.forward.circle.fill",
                    color: .green
                )
                
                DislocKpiPill(
                    title: "Macro Cauchy Stress",
                    value: String(format: "%.2e", maxMacroStress),
                    subtitle: isZeroStress ? "Relaxed to 0 (Full Cauchy)" : "Unrelaxed Boundary",
                    icon: isZeroStress ? "checkmark.circle.fill" : "exclamationmark.triangle.fill",
                    color: isZeroStress ? .green : .red
                )
                
                if let f = root.finalF, f.count >= 2, f[0].count >= 2, f[1].count >= 2 {
                    DislocKpiPill(
                        title: "Macro Deformation F̄",
                        value: String(format: "F̄₁₁=%.5f, F̄₂₂=%.5f", f[0][0], f[1][1]),
                        subtitle: String(format: "F̄₁₂=%.2e", f[0][1]),
                        icon: "arrow.up.left.and.arrow.down.right",
                        color: .teal
                    )
                }
                
                if let lmin = root.lambdaMin {
                    let isStable = root.isStable ?? (lmin > 0)
                    DislocKpiPill(
                        title: "Core Stability λ_min",
                        value: String(format: "%.4e", lmin),
                        subtitle: isStable ? "Stable Core (λ > 0)" : "Unstable (λ ≤ 0)",
                        icon: isStable ? "shield.lefthalf.filled" : "exclamationmark.shield.fill",
                        color: isStable ? .mint : .red
                    )
                }
            }
            .padding(.vertical, 2)
        }
        .frame(maxWidth: .infinity, alignment: .leading)
    }
    
    // MARK: - Chart View Dispatch
    @ViewBuilder
    private func chartView(root: DislocationProfileJSONRoot) -> some View {
        switch plotMode {
        case .radialEnergy:
            radialEnergyChart(root: root)
        case .slipEnergy:
            slipEnergyChart(root: root)
        case .slipStress:
            slipStressChart(root: root)
        case .radialStress:
            radialStressChart(root: root)
        }
    }
    
    // MARK: - Color Scales
    private static let radialEnergyColorScale: KeyValuePairs<String, Color> = [
        "Relaxed (State 1)": Color.purple,
        "Unrelaxed Volterra (State 0)": Color.orange,
        "Linear Fit (Relaxed)": Color.cyan
    ]

    private static let slipEnergyColorScale: KeyValuePairs<String, Color> = [
        "Slip Plane Center": Color.purple,
        "Upper Row (y > y_c)": Color.blue,
        "Lower Row (y < y_c)": Color.red,
        "Relaxed (State 1)": Color.purple,
        "Unrelaxed Volterra (State 0)": Color.orange,
        "Upper Row (Relaxed)": Color.blue,
        "Lower Row (Relaxed)": Color.red,
        "Upper Row (Unrelaxed)": Color.cyan,
        "Lower Row (Unrelaxed)": Color.orange,
        "Slip Plane Center (Unrelaxed)": Color.orange
    ]

    private static let slipStressColorScale: KeyValuePairs<String, Color> = [
        "σ_xx (Normal)": Color.blue,
        "σ_yy (Normal)": Color.green,
        "σ_xy (Shear)": Color.red,
        "von Mises σ_vm": Color.purple,
        "σ_xx (Relaxed)": Color.blue,
        "σ_yy (Relaxed)": Color.green,
        "σ_xy (Relaxed)": Color.red,
        "σ_vm (Relaxed)": Color.purple,
        "σ_xx (Unrelaxed)": Color.cyan,
        "σ_yy (Unrelaxed)": Color.mint,
        "σ_xy (Unrelaxed)": Color.orange,
        "σ_vm (Unrelaxed)": Color.pink
    ]

    private static let radialStressColorScale: KeyValuePairs<String, Color> = [
        "σ_xx": Color.blue,
        "σ_yy": Color.green,
        "σ_xy": Color.red,
        "von Mises": Color.purple,
        "σ_xx (Relaxed)": Color.blue,
        "σ_yy (Relaxed)": Color.green,
        "σ_xy (Relaxed)": Color.red,
        "von Mises (Relaxed)": Color.purple,
        "σ_xx (Unrelaxed)": Color.cyan,
        "σ_yy (Unrelaxed)": Color.mint,
        "σ_xy (Unrelaxed)": Color.orange,
        "von Mises (Unrelaxed)": Color.pink
    ]

    // MARK: - Internal Chart Plot Data Models
    private struct SlipEnergyPoint: Identifiable {
        let id = UUID()
        let seriesName: String
        let color: Color
        let isDashed: Bool
        let dx: Double
        let energy: Double
    }
    
    private struct SlipStressPoint: Identifiable {
        let id = UUID()
        let seriesName: String
        let color: Color
        let isDashed: Bool
        let dx: Double
        let stress: Double
    }
    
    private struct RadialStressPoint: Identifiable {
        let id = UUID()
        let seriesName: String
        let color: Color
        let isDashed: Bool
        let xVal: Double
        let stress: Double
    }

    // MARK: - Hover Overlay Helper
    private func hoverOverlay(proxy: ChartProxy) -> some View {
        GeometryReader { geo in
            Rectangle().fill(Color.clear).contentShape(Rectangle())
                .onContinuousHover { phase in
                    switch phase {
                    case .active(let loc):
                        let origin = geo[proxy.plotAreaFrame].origin
                        let plotX = loc.x - origin.x
                        if let x: Double = proxy.value(atX: plotX) {
                            hoveredX = x
                        }
                    case .ended:
                        hoveredX = nil
                    }
                }
        }
    }

    // MARK: - 1. Radial Energy Chart
    private func radialEnergyChart(root: DislocationProfileJSONRoot) -> some View {
        let b = root.burgersVector
        let showRelaxed = (stateSelection == .relaxed || stateSelection == .compare)
        let showUnrelaxed = (stateSelection == .unrelaxed || stateSelection == .compare)
        let isCompare = (stateSelection == .compare)
        
        let rel = root.relaxed
        let unrel = root.unrelaxed
        
        var asymptoteStart: (x: Double, y: Double)? = nil
        var asymptoteEnd: (x: Double, y: Double)? = nil
        if radialScale == .log, let rel = rel {
            let k = rel.coreFit.elasticK
            let e0 = rel.coreFit.coreEnergy
            let posLogPoints = rel.radial.filter { $0.logR >= 0.0 }
            if let first = posLogPoints.first, let last = posLogPoints.last {
                asymptoteStart = (first.logR, k * first.logR + e0)
                asymptoteEnd = (last.logR, k * last.logR + e0)
            }
        }
        
        return Chart {
            if showRelaxed, let rel = rel {
                ForEach(rel.radial) { pt in
                    let xVal = radialScale == .linear ? (pt.r / b) : pt.logR
                    LineMark(
                        x: .value("x", xVal),
                        y: .value("Cumulative Energy E(r)", pt.cumulativeEnergy),
                        series: .value("State", "Relaxed (State 1)")
                    )
                    .foregroundStyle(Color.purple)
                    .lineStyle(StrokeStyle(lineWidth: 2.5))
                    
                    PointMark(
                        x: .value("x", xVal),
                        y: .value("Cumulative Energy E(r)", pt.cumulativeEnergy)
                    )
                    .foregroundStyle(Color.purple)
                    .symbolSize(14)
                }
                
                if let start = asymptoteStart, let end = asymptoteEnd {
                    LineMark(
                        x: .value("x", start.x),
                        y: .value("Cumulative Energy E(r)", start.y),
                        series: .value("State", "Linear Fit (Relaxed)")
                    )
                    .foregroundStyle(Color.cyan)
                    .lineStyle(StrokeStyle(lineWidth: 1.5, dash: [4, 4]))
                    
                    LineMark(
                        x: .value("x", end.x),
                        y: .value("Cumulative Energy E(r)", end.y),
                        series: .value("State", "Linear Fit (Relaxed)")
                    )
                    .foregroundStyle(Color.cyan)
                    .lineStyle(StrokeStyle(lineWidth: 1.5, dash: [4, 4]))
                }
            }
            
            if showUnrelaxed {
                ForEach(unrel.radial) { pt in
                    let xVal = radialScale == .linear ? (pt.r / b) : pt.logR
                    LineMark(
                        x: .value("x", xVal),
                        y: .value("Cumulative Energy E(r)", pt.cumulativeEnergy),
                        series: .value("State", "Unrelaxed Volterra (State 0)")
                    )
                    .foregroundStyle(Color.orange)
                    .lineStyle(StrokeStyle(lineWidth: 2.0, dash: isCompare ? [5, 3] : []))
                    
                    PointMark(
                        x: .value("x", xVal),
                        y: .value("Cumulative Energy E(r)", pt.cumulativeEnergy)
                    )
                    .foregroundStyle(Color.orange)
                    .symbolSize(12)
                }
            }
            
            if let hx = hoveredX {
                RuleMark(x: .value("Hover", hx))
                    .lineStyle(StrokeStyle(lineWidth: 1, dash: [3, 3]))
                    .foregroundStyle(Color.secondary.opacity(0.7))
            }
        }
        .chartForegroundStyleScale(Self.radialEnergyColorScale)
        .chartXAxisLabel(radialScale == .linear ? "Radial Distance from Core r / b" : "Logarithmic Distance ln(r / b)")
        .chartYAxisLabel("Cumulative Strain Energy E(r)")
        .chartOverlay { proxy in hoverOverlay(proxy: proxy) }
    }
    
    // MARK: - Filtering Helper for Slip Plane Domain
    private func filterSlipPoints<T>(_ list: [T], rMax: Double, b: Double) -> [T] {
        switch slipRange {
        case .full:
            return list
        case .rFree:
            if let pts = list as? [DislocationSlipPointModel] {
                return pts.filter { abs($0.dx) <= rMax } as! [T]
            }
            if let pts = list as? [DislocationSlipAvgPointModel] {
                return pts.filter { abs($0.dx) <= rMax } as! [T]
            }
            return list
        case .nearCore:
            let limit = nearCoreSpan * b
            if let pts = list as? [DislocationSlipPointModel] {
                return pts.filter { abs($0.dx) <= limit } as! [T]
            }
            if let pts = list as? [DislocationSlipAvgPointModel] {
                return pts.filter { abs($0.dx) <= limit } as! [T]
            }
            return list
        }
    }
    
    // MARK: - 2. Slip Plane Energy Chart E(x - xc)
    private func slipEnergyChart(root: DislocationProfileJSONRoot) -> some View {
        let rMax = root.rFree * root.latticeConstant
        let b = root.burgersVector
        let showRelaxed = (stateSelection == .relaxed || stateSelection == .compare)
        let showUnrelaxed = (stateSelection == .unrelaxed || stateSelection == .compare)
        let isCompare = (stateSelection == .compare)
        
        var plotPoints: [SlipEnergyPoint] = []
        
        if showRelaxed, let rel = root.relaxed ?? (stateSelection == .relaxed ? root.unrelaxed : nil) {
            if slipRow == .average || slipRow == .bothRows {
                let avgPts = filterSlipPoints(rel.slipAvg.sorted { $0.dx < $1.dx }, rMax: rMax, b: b)
                for pt in avgPts {
                    plotPoints.append(SlipEnergyPoint(
                        seriesName: isCompare ? "Relaxed (State 1)" : "Slip Plane Center",
                        color: .purple,
                        isDashed: false,
                        dx: pt.dx,
                        energy: pt.energy
                    ))
                }
            }
            if slipRow == .upper || slipRow == .bothRows {
                let upperPts = filterSlipPoints(rel.slipPlane.filter { $0.row == "upper" }.sorted { $0.dx < $1.dx }, rMax: rMax, b: b)
                for pt in upperPts {
                    plotPoints.append(SlipEnergyPoint(
                        seriesName: isCompare ? "Upper Row (Relaxed)" : "Upper Row (y > y_c)",
                        color: .blue,
                        isDashed: false,
                        dx: pt.dx,
                        energy: pt.energy
                    ))
                }
            }
            if slipRow == .lower || slipRow == .bothRows {
                let lowerPts = filterSlipPoints(rel.slipPlane.filter { $0.row == "lower" }.sorted { $0.dx < $1.dx }, rMax: rMax, b: b)
                for pt in lowerPts {
                    plotPoints.append(SlipEnergyPoint(
                        seriesName: isCompare ? "Lower Row (Relaxed)" : "Lower Row (y < y_c)",
                        color: .red,
                        isDashed: false,
                        dx: pt.dx,
                        energy: pt.energy
                    ))
                }
            }
        }
        
        if showUnrelaxed {
            let unrel = root.unrelaxed
            if slipRow == .average || slipRow == .bothRows {
                let avgPts = filterSlipPoints(unrel.slipAvg.sorted { $0.dx < $1.dx }, rMax: rMax, b: b)
                for pt in avgPts {
                    plotPoints.append(SlipEnergyPoint(
                        seriesName: isCompare ? "Unrelaxed Volterra (State 0)" : "Slip Plane Center (Unrelaxed)",
                        color: .orange,
                        isDashed: true,
                        dx: pt.dx,
                        energy: pt.energy
                    ))
                }
            }
            if slipRow == .upper || slipRow == .bothRows {
                let upperPts = filterSlipPoints(unrel.slipPlane.filter { $0.row == "upper" }.sorted { $0.dx < $1.dx }, rMax: rMax, b: b)
                for pt in upperPts {
                    plotPoints.append(SlipEnergyPoint(
                        seriesName: isCompare ? "Upper Row (Unrelaxed)" : "Upper Row (Unrelaxed)",
                        color: .cyan,
                        isDashed: true,
                        dx: pt.dx,
                        energy: pt.energy
                    ))
                }
            }
            if slipRow == .lower || slipRow == .bothRows {
                let lowerPts = filterSlipPoints(unrel.slipPlane.filter { $0.row == "lower" }.sorted { $0.dx < $1.dx }, rMax: rMax, b: b)
                for pt in lowerPts {
                    plotPoints.append(SlipEnergyPoint(
                        seriesName: isCompare ? "Lower Row (Unrelaxed)" : "Lower Row (Unrelaxed)",
                        color: .orange,
                        isDashed: true,
                        dx: pt.dx,
                        energy: pt.energy
                    ))
                }
            }
        }
        
        return Chart {
            ForEach(plotPoints) { pt in
                LineMark(
                    x: .value("x - x_c", pt.dx),
                    y: .value("Energy E", pt.energy),
                    series: .value("Row", pt.seriesName)
                )
                .foregroundStyle(pt.color)
                .lineStyle(StrokeStyle(lineWidth: pt.isDashed ? 1.8 : 2.2, dash: pt.isDashed ? [5, 3] : []))
                
                PointMark(x: .value("x - x_c", pt.dx), y: .value("Energy E", pt.energy))
                    .foregroundStyle(pt.color)
                    .symbolSize(pt.isDashed ? 10 : 12)
            }
            
            // Core vertical reference line at x - xc = 0
            RuleMark(x: .value("Core x_c", 0.0))
                .lineStyle(StrokeStyle(lineWidth: 1.5, dash: [4, 4]))
                .foregroundStyle(Color.red.opacity(0.8))
                .annotation(position: .top, alignment: .leading) {
                    Text("Core x_c")
                        .font(.caption2.bold())
                        .foregroundColor(.red)
                        .padding(2)
                        .background(BlueHour.controlSurface.opacity(0.8))
                }
            
            if let hx = hoveredX {
                RuleMark(x: .value("Hover", hx))
                    .lineStyle(StrokeStyle(lineWidth: 1, dash: [3, 3]))
                    .foregroundStyle(Color.secondary.opacity(0.7))
            }
        }
        .chartForegroundStyleScale(Self.slipEnergyColorScale)
        .chartXAxisLabel("Distance along Slip Plane (x - x_c)")
        .chartYAxisLabel("Atomic Energy E(x)")
        .chartOverlay { proxy in hoverOverlay(proxy: proxy) }
    }
    
    // MARK: - 3. Slip Plane Stresses σ(x - xc)
    private func slipStressChart(root: DislocationProfileJSONRoot) -> some View {
        let rMax = root.rFree * root.latticeConstant
        let b = root.burgersVector
        let showRelaxed = (stateSelection == .relaxed || stateSelection == .compare)
        let showUnrelaxed = (stateSelection == .unrelaxed || stateSelection == .compare)
        let isCompare = (stateSelection == .compare)
        
        func extractPoints(from state: DislocationStateDataModel) -> [DislocationSlipAvgPointModel] {
            let raw: [DislocationSlipAvgPointModel]
            if slipRow == .upper {
                raw = state.slipPlane.filter { $0.row == "upper" }.map {
                    DislocationSlipAvgPointModel(dx: $0.dx, energy: $0.energy, sXx: $0.sXx, sYy: $0.sYy, sXy: $0.sXy, sHydro: $0.sHydro, sVm: $0.sVm)
                }
            } else if slipRow == .lower {
                raw = state.slipPlane.filter { $0.row == "lower" }.map {
                    DislocationSlipAvgPointModel(dx: $0.dx, energy: $0.energy, sXx: $0.sXx, sYy: $0.sYy, sXy: $0.sXy, sHydro: $0.sHydro, sVm: $0.sVm)
                }
            } else {
                raw = state.slipAvg
            }
            return filterSlipPoints(raw.sorted { $0.dx < $1.dx }, rMax: rMax, b: b)
        }
        
        var plotPoints: [SlipStressPoint] = []
        if showRelaxed, let rel = root.relaxed ?? (stateSelection == .relaxed ? root.unrelaxed : nil) {
            let pts = extractPoints(from: rel)
            for p in pts {
                if showSxx {
                    plotPoints.append(SlipStressPoint(seriesName: isCompare ? "σ_xx (Relaxed)" : "σ_xx (Normal)", color: .blue, isDashed: false, dx: p.dx, stress: p.sXx))
                }
                if showSyy {
                    plotPoints.append(SlipStressPoint(seriesName: isCompare ? "σ_yy (Relaxed)" : "σ_yy (Normal)", color: .green, isDashed: false, dx: p.dx, stress: p.sYy))
                }
                if showSxy {
                    plotPoints.append(SlipStressPoint(seriesName: isCompare ? "σ_xy (Relaxed)" : "σ_xy (Shear)", color: .red, isDashed: false, dx: p.dx, stress: p.sXy))
                }
                if showSVM {
                    plotPoints.append(SlipStressPoint(seriesName: isCompare ? "σ_vm (Relaxed)" : "von Mises σ_vm", color: .purple, isDashed: false, dx: p.dx, stress: p.sVm))
                }
            }
        }
        if showUnrelaxed {
            let pts = extractPoints(from: root.unrelaxed)
            for p in pts {
                if showSxx {
                    plotPoints.append(SlipStressPoint(seriesName: isCompare ? "σ_xx (Unrelaxed)" : "σ_xx (Normal)", color: .cyan, isDashed: true, dx: p.dx, stress: p.sXx))
                }
                if showSyy {
                    plotPoints.append(SlipStressPoint(seriesName: isCompare ? "σ_yy (Unrelaxed)" : "σ_yy (Normal)", color: .mint, isDashed: true, dx: p.dx, stress: p.sYy))
                }
                if showSxy {
                    plotPoints.append(SlipStressPoint(seriesName: isCompare ? "σ_xy (Unrelaxed)" : "σ_xy (Shear)", color: .orange, isDashed: true, dx: p.dx, stress: p.sXy))
                }
                if showSVM {
                    plotPoints.append(SlipStressPoint(seriesName: isCompare ? "σ_vm (Unrelaxed)" : "von Mises σ_vm", color: .pink, isDashed: true, dx: p.dx, stress: p.sVm))
                }
            }
        }
        
        return Chart {
            ForEach(plotPoints) { pt in
                LineMark(
                    x: .value("x - x_c", pt.dx),
                    y: .value("Stress", pt.stress),
                    series: .value("Component", pt.seriesName)
                )
                .foregroundStyle(pt.color)
                .lineStyle(StrokeStyle(lineWidth: pt.isDashed ? 1.8 : 2.0, dash: pt.isDashed ? [5, 3] : []))
                
                PointMark(x: .value("x - x_c", pt.dx), y: .value("Stress", pt.stress))
                    .foregroundStyle(pt.color)
                    .symbolSize(pt.isDashed ? 10 : 12)
            }
            
            // Core line marker
            RuleMark(x: .value("Core x_c", 0.0))
                .lineStyle(StrokeStyle(lineWidth: 1.5, dash: [4, 4]))
                .foregroundStyle(Color.secondary.opacity(0.7))
                .annotation(position: .top) {
                    Text("x_c = 0").font(.caption2.bold()).foregroundColor(.secondary)
                }
            
            RuleMark(y: .value("Zero", 0.0))
                .lineStyle(StrokeStyle(lineWidth: 0.8))
                .foregroundStyle(Color.secondary.opacity(0.4))
            
            if let hx = hoveredX {
                RuleMark(x: .value("Hover", hx))
                    .lineStyle(StrokeStyle(lineWidth: 1, dash: [3, 3]))
                    .foregroundStyle(Color.secondary.opacity(0.7))
            }
        }
        .chartForegroundStyleScale(Self.slipStressColorScale)
        .chartXAxisLabel("Distance along Slip Plane (x - x_c)")
        .chartYAxisLabel("Cauchy Stress Components")
        .chartOverlay { proxy in hoverOverlay(proxy: proxy) }
    }
    
    // MARK: - 4. Radial Stresses Chart σ(r)
    private func radialStressChart(root: DislocationProfileJSONRoot) -> some View {
        let b = root.burgersVector
        let showRelaxed = (stateSelection == .relaxed || stateSelection == .compare)
        let showUnrelaxed = (stateSelection == .unrelaxed || stateSelection == .compare)
        let isCompare = (stateSelection == .compare)
        
        var plotPoints: [RadialStressPoint] = []
        if showRelaxed, let rel = root.relaxed ?? (stateSelection == .relaxed ? root.unrelaxed : nil) {
            for pt in rel.radial {
                let x = (radialScale == .linear) ? (pt.r / b) : pt.logR
                if showSxx { plotPoints.append(RadialStressPoint(seriesName: isCompare ? "σ_xx (Relaxed)" : "σ_xx", color: .blue, isDashed: false, xVal: x, stress: pt.sXx)) }
                if showSyy { plotPoints.append(RadialStressPoint(seriesName: isCompare ? "σ_yy (Relaxed)" : "σ_yy", color: .green, isDashed: false, xVal: x, stress: pt.sYy)) }
                if showSxy { plotPoints.append(RadialStressPoint(seriesName: isCompare ? "σ_xy (Relaxed)" : "σ_xy", color: .red, isDashed: false, xVal: x, stress: pt.sXy)) }
                if showSVM { plotPoints.append(RadialStressPoint(seriesName: isCompare ? "von Mises (Relaxed)" : "von Mises", color: .purple, isDashed: false, xVal: x, stress: pt.sVm)) }
            }
        }
        if showUnrelaxed {
            for pt in root.unrelaxed.radial {
                let x = (radialScale == .linear) ? (pt.r / b) : pt.logR
                if showSxx { plotPoints.append(RadialStressPoint(seriesName: isCompare ? "σ_xx (Unrelaxed)" : "σ_xx", color: .cyan, isDashed: true, xVal: x, stress: pt.sXx)) }
                if showSyy { plotPoints.append(RadialStressPoint(seriesName: isCompare ? "σ_yy (Unrelaxed)" : "σ_yy", color: .mint, isDashed: true, xVal: x, stress: pt.sYy)) }
                if showSxy { plotPoints.append(RadialStressPoint(seriesName: isCompare ? "σ_xy (Unrelaxed)" : "σ_xy", color: .orange, isDashed: true, xVal: x, stress: pt.sXy)) }
                if showSVM { plotPoints.append(RadialStressPoint(seriesName: isCompare ? "von Mises (Unrelaxed)" : "von Mises", color: .pink, isDashed: true, xVal: x, stress: pt.sVm)) }
            }
        }
        
        return Chart {
            ForEach(plotPoints) { pt in
                LineMark(
                    x: .value(radialScale == .linear ? "Radius r/b" : "ln(r/b)", pt.xVal),
                    y: .value("Stress", pt.stress),
                    series: .value("Component", pt.seriesName)
                )
                .foregroundStyle(pt.color)
                .lineStyle(StrokeStyle(lineWidth: pt.isDashed ? 1.8 : 2.0, dash: pt.isDashed ? [5, 3] : []))
            }
            
            RuleMark(y: .value("Zero", 0.0))
                .lineStyle(StrokeStyle(lineWidth: 0.8))
                .foregroundStyle(Color.secondary.opacity(0.4))
        }
        .chartForegroundStyleScale(Self.radialStressColorScale)
        .chartXAxisLabel(radialScale == .linear ? "Radial Distance r / b" : "Logarithmic Distance ln(r / b)")
        .chartYAxisLabel("Average Shell Stress")
        .chartOverlay { proxy in hoverOverlay(proxy: proxy) }
    }
    
    // MARK: - Empty State View
    private var emptyStateView: some View {
        VStack(spacing: 16) {
            Image(systemName: "circle.dashed")
                .font(.system(size: 48))
                .foregroundColor(.purple)
            
            Text("Volterra Dislocation Analysis Ready")
                .font(.title2.bold())
            
            Text("Click 'Launch Simulation' to relax the dislocation core and inspect:\n• Radial Energy E(r) and Dislocation Core Energy E_core\n• Core Line Energy profile E(x - x_c)\n• Cauchy stress components (σ_xx, σ_yy, σ_xy) along the slip plane")
                .font(.subheadline)
                .multilineTextAlignment(.center)
                .foregroundColor(.secondary)
                .frame(maxWidth: 500)
            
            if runner.status.isRunning {
                ProgressView("Simulating dislocation core relaxation...")
                    .progressViewStyle(.circular)
                    .padding(.top, 8)
            } else {
                Button(action: { runner.pollDislocationProfiles() }) {
                    Label("Reload Output Files", systemImage: "arrow.clockwise")
                }
                .buttonStyle(.bordered)
                .controlSize(.regular)
                .padding(.top, 4)
            }
        }
        .frame(maxWidth: .infinity, maxHeight: .infinity)
        .padding(32)
    }
    
    // MARK: - Data Export Methods
    
    private func generateActivePlotCSV() -> String {
        guard let root = data else { return "" }
        var csv = ""
        
        let states: [(name: String, profile: DislocationStateDataModel)]
        switch stateSelection {
        case .relaxed:
            if let rel = root.relaxed { states = [("relaxed", rel)] }
            else { states = [("unrelaxed", root.unrelaxed)] }
        case .unrelaxed:
            states = [("unrelaxed", root.unrelaxed)]
        case .compare:
            var list: [(String, DislocationStateDataModel)] = [("unrelaxed", root.unrelaxed)]
            if let rel = root.relaxed { list.append(("relaxed", rel)) }
            states = list
        }
        
        switch plotMode {
        case .radialEnergy:
            csv += "state,r,log_r_over_b,cumulative_energy,energy_density,avg_atom_energy\n"
            for s in states {
                for pt in s.profile.radial {
                    csv += "\(s.name),\(pt.r),\(pt.logR),\(pt.cumulativeEnergy),\(pt.energyDensity),\(pt.avgAtomEnergy)\n"
                }
            }
        case .slipEnergy:
            csv += "state,row,x,dx,energy\n"
            for s in states {
                switch slipRow {
                case .average:
                    for pt in s.profile.slipAvg {
                        csv += "\(s.name),average,\(pt.dx + root.coreX),\(pt.dx),\(pt.energy)\n"
                    }
                case .upper:
                    for pt in s.profile.slipPlane.filter({ $0.row == "upper" }) {
                        csv += "\(s.name),upper,\(pt.x),\(pt.dx),\(pt.energy)\n"
                    }
                case .lower:
                    for pt in s.profile.slipPlane.filter({ $0.row == "lower" }) {
                        csv += "\(s.name),lower,\(pt.x),\(pt.dx),\(pt.energy)\n"
                    }
                case .bothRows:
                    for pt in s.profile.slipPlane {
                        csv += "\(s.name),\(pt.row),\(pt.x),\(pt.dx),\(pt.energy)\n"
                    }
                }
            }
        case .slipStress:
            csv += "state,row,x,dx,s_xx,s_yy,s_xy,s_hydro,s_vm\n"
            for s in states {
                switch slipRow {
                case .average:
                    for pt in s.profile.slipAvg {
                        csv += "\(s.name),average,\(pt.dx + root.coreX),\(pt.dx),\(pt.sXx),\(pt.sYy),\(pt.sXy),\(pt.sHydro),\(pt.sVm)\n"
                    }
                case .upper:
                    for pt in s.profile.slipPlane.filter({ $0.row == "upper" }) {
                        csv += "\(s.name),upper,\(pt.x),\(pt.dx),\(pt.sXx),\(pt.sYy),\(pt.sXy),\(pt.sHydro),\(pt.sVm)\n"
                    }
                case .lower:
                    for pt in s.profile.slipPlane.filter({ $0.row == "lower" }) {
                        csv += "\(s.name),lower,\(pt.x),\(pt.dx),\(pt.sXx),\(pt.sYy),\(pt.sXy),\(pt.sHydro),\(pt.sVm)\n"
                    }
                case .bothRows:
                    for pt in s.profile.slipPlane {
                        csv += "\(s.name),\(pt.row),\(pt.x),\(pt.dx),\(pt.sXx),\(pt.sYy),\(pt.sXy),\(pt.sHydro),\(pt.sVm)\n"
                    }
                }
            }
        case .radialStress:
            csv += "state,r,log_r_over_b,s_xx,s_yy,s_xy,s_hydro,s_vm\n"
            for s in states {
                for pt in s.profile.radial {
                    csv += "\(s.name),\(pt.r),\(pt.logR),\(pt.sXx),\(pt.sYy),\(pt.sXy),\(pt.sHydro),\(pt.sVm)\n"
                }
            }
        }
        return csv
    }
    
    private func exportActivePlotCSV() {
        let csv = generateActivePlotCSV()
        ExportDataHelper.saveStringToFile(content: csv, defaultFileName: "dislocation_\(plotMode.filePrefix).csv")
    }
    
    private func exportAllProfilesCSV() {
        let dir = runner.activeConfig?.targetDirectory ?? FileManager.default.currentDirectoryPath
        let csvPath = "\(dir)/dislocation_profiles.csv"
        if let content = try? String(contentsOfFile: csvPath, encoding: .utf8) {
            ExportDataHelper.saveStringToFile(content: content, defaultFileName: "dislocation_profiles_full.csv")
        } else {
            let csv = generateActivePlotCSV()
            ExportDataHelper.saveStringToFile(content: csv, defaultFileName: "dislocation_profiles.csv")
        }
    }
    
    private func exportAllProfilesJSON() {
        let dir = runner.activeConfig?.targetDirectory ?? FileManager.default.currentDirectoryPath
        let jsonPath = "\(dir)/dislocation_profiles.json"
        if let content = try? String(contentsOfFile: jsonPath, encoding: .utf8) {
            ExportDataHelper.saveStringToFile(content: content, defaultFileName: "dislocation_profiles.json", allowedExtensions: ["json", "txt"])
        }
    }
    
    private func copyActivePlotCSVToClipboard() {
        ExportDataHelper.copyToClipboard(text: generateActivePlotCSV())
    }
    
    private func copyAllProfilesCSVToClipboard() {
        let dir = runner.activeConfig?.targetDirectory ?? FileManager.default.currentDirectoryPath
        let csvPath = "\(dir)/dislocation_profiles.csv"
        if let content = try? String(contentsOfFile: csvPath, encoding: .utf8) {
            ExportDataHelper.copyToClipboard(text: content)
        } else {
            ExportDataHelper.copyToClipboard(text: generateActivePlotCSV())
        }
    }
}

// MARK: - KPI Pill Helper
private struct DislocKpiPill: View {
    let title: String
    let value: String
    var subtitle: String? = nil
    let icon: String
    let color: Color
    
    var body: some View {
        HStack(spacing: 8) {
            Circle()
                .fill(color.opacity(0.15))
                .frame(width: 30, height: 30)
                .overlay(Image(systemName: icon).foregroundColor(color).font(.system(size: 13)))
            
            VStack(alignment: .leading, spacing: 2) {
                Text(title)
                    .font(.system(size: 10, weight: .semibold))
                    .foregroundColor(.secondary)
                    .lineLimit(1)
                HStack(spacing: 4) {
                    Text(value)
                        .font(.system(size: 12.5, weight: .bold, design: .monospaced))
                        .lineLimit(1)
                    if let sub = subtitle {
                        Text("(\(sub))")
                            .font(.system(size: 10))
                            .foregroundColor(.secondary)
                            .lineLimit(1)
                    }
                }
            }
        }
        .padding(.horizontal, 10)
        .padding(.vertical, 6)
        .background(BlueHour.pillSurface.opacity(0.6))
        .cornerRadius(8)
        .fixedSize(horizontal: true, vertical: false)
    }
}
