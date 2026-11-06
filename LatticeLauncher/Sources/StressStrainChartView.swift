import SwiftUI
import Charts

public enum ChartMetricMode: String, CaseIterable, Identifiable {
    case stress = "Stress σ(α)"
    case energy = "Energy E(α)"
    case eigenvalue = "λ_min(α)"
    case gbEnergy = "γ_GB(Δθ)"
    
    public var id: String { rawValue }
}

public enum ChartScaleMode: String, CaseIterable, Identifiable {
    case autoFit = "Auto-Fit"
    case fullSchedule = "Full α"
    case fromZero = "From 0"
    
    public var id: String { rawValue }
    
    public var tooltip: String {
        switch self {
        case .autoFit:
            return "Clever Auto-Fit: dynamically zooms to active data range"
        case .fullSchedule:
            return "Full Schedule: starts at α_start up to α_end"
        case .fromZero:
            return "Origin at 0: shows strain from 0.0"
        }
    }
}

struct StressStrainChartView: View {
    @EnvironmentObject var runner: ProcessRunner
    @State private var selectedAlpha: Double?
    @State private var chartMode: ChartMetricMode = .stress
    @State private var stressComponent: StressComponentType = .primary
    @State private var scaleMode: ChartScaleMode = .autoFit
    @State private var showRemeshEvents: Bool = true
    @State private var showYieldLine: Bool = true
    @State private var useInteriorGBGauge: Bool = true
    @State private var showGBSigmaLabels: Bool = true
    
    /// Torsion runs report the rim torque T = Σ r × f as their primary quantity.
    private var primaryIsTorque: Bool {
        guard let g = runner.activeConfig?.geometryConfig else { return false }
        if g.primary_stress_component == "torque" { return true }
        return g.primary_stress_component == "auto" && g.boundary_rules.contains {
            $0.enabled && $0.kind == .affine && ($0.affine_mode == "circumferential_shear" || $0.affine_mode == "torsion")
        }
    }

    private func stressValue(for pt: StressStrainPoint) -> Double {
        pt.value(for: stressComponent)
    }
    
    private func gbEnergyValue(for pt: GBEnergyPoint) -> Double {
        useInteriorGBGauge ? pt.gbEnergyInterior : pt.gbEnergyTotal
    }
    
    private var displayGBPoints: [GBEnergyPoint] {
        if !runner.gbEnergyPoints.isEmpty {
            return runner.gbEnergyPoints
        }
        if let cur = runner.currentGBEnergy {
            return [cur]
        }
        return []
    }
    
    private var currentPeakStress: Double {
        runner.curvePoints.map { abs(stressValue(for: $0)) }.max() ?? 0.0
    }
    
    // Remesh / surgery points (excluding initial setup mesh at iteration 0)
    private var remeshEvents: [StressStrainPoint] {
        runner.curvePoints.filter { $0.iteration > 0 && $0.shouldRemesh }
    }
    
    // MARK: - Dynamic Clever Domain Calculations
    
    private var activeAlphas: [Double] {
        if chartMode == .gbEnergy {
            return displayGBPoints.map(\.misorientationDeg)
        }
        if chartMode == .eigenvalue && !runner.eigenPoints.isEmpty {
            return runner.eigenPoints.map(\.alpha)
        }
        return runner.curvePoints.map(\.alpha)
    }
    
    private var xDomain: ClosedRange<Double> {
        if chartMode == .gbEnergy {
            let angles = displayGBPoints.map(\.misorientationDeg)
            guard !angles.isEmpty else { return 0.0...60.0 }
            let minA = angles.min() ?? 0.0
            let maxA = angles.max() ?? 60.0
            if scaleMode == .fromZero {
                return 0.0...max(10.0, maxA * 1.08)
            }
            let span = max(4.0, maxA - minA)
            let lo = max(0.0, minA - span * 0.08)
            let hi = maxA + span * 0.08
            return lo < hi ? lo...hi : 0.0...60.0
        }
        
        let alphas = activeAlphas
        let firstAlpha = alphas.first ?? (runner.activeConfig?.alphaStart ?? 0.14)
        
        guard !alphas.isEmpty else {
            let endAlpha = runner.activeConfig?.alphaEnd ?? 1.0
            switch scaleMode {
            case .autoFit:
                return (firstAlpha - 0.001)...(firstAlpha + 0.005)
            case .fullSchedule:
                return firstAlpha...endAlpha
            case .fromZero:
                return 0.0...endAlpha
            }
        }
        
        let minA = alphas.min() ?? firstAlpha
        let maxA = alphas.max() ?? (minA + 0.001)
        let span = maxA - minA
        
        switch scaleMode {
        case .autoFit:
            let effectiveSpan = max(0.0008, span)
            let padLeft = effectiveSpan * 0.025
            let padRight = effectiveSpan * 0.065
            let xMin = minA - padLeft
            let xMax = maxA + padRight
            return xMin <= xMax ? xMin...xMax : (minA - 0.001)...(minA + 0.005)
            
        case .fullSchedule:
            let schedStart = runner.activeConfig?.alphaStart ?? minA
            let schedEnd = max(maxA + 0.005, runner.activeConfig?.alphaEnd ?? 1.0)
            let fullSpan = max(0.01, schedEnd - schedStart)
            let xMin = schedStart - fullSpan * 0.01
            let xMax = schedEnd + fullSpan * 0.02
            return xMin <= xMax ? xMin...xMax : schedStart...schedEnd
            
        case .fromZero:
            let schedEnd = max(0.05, max(maxA * 1.05, runner.activeConfig?.alphaEnd ?? 1.0))
            return 0.0...schedEnd
        }
    }
    
    private var yDomain: ClosedRange<Double> {
        switch chartMode {
        case .gbEnergy:
            let pts = displayGBPoints
            guard !pts.isEmpty else { return 0.0...1.0 }
            let vals = pts.map { gbEnergyValue(for: $0) }
            let minY = vals.min() ?? 0.0
            let maxY = vals.max() ?? 0.5
            let span = max(0.05, maxY - minY)
            let lo = scaleMode == .fromZero ? 0.0 : max(0.0, minY - span * 0.12)
            let hi = maxY + span * 0.14
            return lo < hi ? lo...hi : 0.0...max(0.5, maxY * 1.2)
            
        case .stress:
            guard !runner.curvePoints.isEmpty else { return 0.0...0.05 }
            let vals = runner.curvePoints.map { stressValue(for: $0) }
            let minY = vals.min() ?? 0.0
            let maxY = vals.max() ?? 0.01
            let span = maxY - minY
            
            if span < 1e-6 {
                let pad = max(abs(minY) * 0.1, 0.002)
                let yMin = min(0.0, minY - pad)
                return yMin...(maxY + pad)
            }
            
            let pad = span * 0.08
            let yMin: Double
            if scaleMode == .fromZero || (minY >= 0 && minY < span * 0.3) {
                yMin = min(0.0, minY - pad * 0.5)
            } else {
                yMin = minY - pad
            }
            
            var yMax = maxY + pad
            let peakVal = currentPeakStress
            if showYieldLine && peakVal > 0 && maxY > 0 {
                yMax = max(yMax, peakVal + pad * 0.5)
            }
            return yMin <= yMax ? yMin...yMax : 0.0...0.05
            
        case .energy:
            guard !runner.curvePoints.isEmpty else { return 0.0...100.0 }
            let vals = runner.curvePoints.map(\.postEnergy)
            let minY = vals.min() ?? 0.0
            let maxY = vals.max() ?? 0.01
            let span = maxY - minY
            
            if scaleMode == .fromZero {
                let lo = min(0.0, minY * 1.1)
                let hi = max(0.01, maxY * 1.1)
                return lo...hi
            }
            
            if span < 1e-6 {
                let pad = max(abs(minY) * 0.02, 0.01)
                return (minY - pad)...(maxY + pad)
            }
            
            let pad = span * 0.08
            let yMin = minY - pad
            let yMax = maxY + pad
            return yMin <= yMax ? yMin...yMax : (minY - 0.1)...(maxY + 0.1)
            
        case .eigenvalue:
            guard !runner.eigenPoints.isEmpty else { return -0.05...0.5 }
            let vals = runner.eigenPoints.map(\.lambdaMin)
            let minY = vals.min() ?? 0.0
            let maxY = vals.max() ?? 0.1
            let span = max(maxY - minY, 1e-5)
            let pad = span * 0.12
            let yMin = min(0.0, minY) - pad
            let yMax = max(0.0, maxY) + pad
            return yMin < yMax ? yMin...yMax : -0.05...0.5
        }
    }
    
    // MARK: - Axis Formatter Helpers
    
    private func formatXValue(_ val: Double) -> String {
        if chartMode == .gbEnergy {
            return String(format: "%.1f°", val)
        }
        let span = xDomain.upperBound - xDomain.lowerBound
        if span <= 0.002 {
            return String(format: "%.5f", val)
        } else if span <= 0.02 {
            return String(format: "%.4f", val)
        } else if span <= 0.2 {
            return String(format: "%.3f", val)
        } else {
            return String(format: "%.2f", val)
        }
    }
    
    private func formatYValue(_ val: Double) -> String {
        let span = yDomain.upperBound - yDomain.lowerBound
        switch chartMode {
        case .stress, .eigenvalue, .gbEnergy:
            if span <= 0.002 {
                return String(format: "%.5f", val)
            } else if span <= 0.02 {
                return String(format: "%.4f", val)
            } else if span <= 0.2 {
                return String(format: "%.3f", val)
            } else {
                return String(format: "%.2f", val)
            }
        case .energy:
            if abs(val) >= 10_000 || (abs(val) < 0.001 && abs(val) > 0) {
                return String(format: "%.3e", val)
            } else if span <= 0.01 {
                return String(format: "%.4f", val)
            } else if span <= 0.1 {
                return String(format: "%.3f", val)
            } else {
                return String(format: "%.2f", val)
            }
        }
    }
    
    // MARK: - View Body
    
    var body: some View {
        VStack(spacing: 0) {
            // Chart Control Toolbar: Pinned left pickers + horizontally scrollable options/badges
            HStack(spacing: 10) {
                // 1. Pinned Metric Selector: Stress vs Energy vs Eigenvalue vs GB Energy
                Picker("", selection: $chartMode) {
                    ForEach(ChartMetricMode.allCases) { mode in
                        Text(mode.rawValue).tag(mode)
                    }
                }
                .pickerStyle(.segmented)
                .labelsHidden()
                .controlSize(.small)
                .fixedSize()
                .help("Select plotted quantity: Stress / Load σ(α), Internal Energy E(α), Smallest Stiffness Eigenvalue λ_min(α), or Grain Boundary Energy γ_GB(Δθ)")
                
                // 1b. Stress / Nanoindentation Component Dropdown Selector
                if chartMode == .stress {
                    Picker("", selection: $stressComponent) {
                        Section(header: Text("Primary")) {
                            Text(StressComponentType.primary.rawValue).tag(StressComponentType.primary)
                        }
                        Section(header: Text("Cauchy Stress & Pressure")) {
                            Text(StressComponentType.cauchyXY.rawValue).tag(StressComponentType.cauchyXY)
                            Text(StressComponentType.cauchyXX.rawValue).tag(StressComponentType.cauchyXX)
                            Text(StressComponentType.cauchyYY.rawValue).tag(StressComponentType.cauchyYY)
                            Text(StressComponentType.hydrostaticPressure.rawValue).tag(StressComponentType.hydrostaticPressure)
                        }
                        Section(header: Text("Aligned Piola / Projected Stress")) {
                            Text(StressComponentType.piola12.rawValue).tag(StressComponentType.piola12)
                            Text(StressComponentType.piola11.rawValue).tag(StressComponentType.piola11)
                            Text(StressComponentType.piola22.rawValue).tag(StressComponentType.piola22)
                        }
                        Section(header: Text("Von Mises Equivalent Stress (2D)")) {
                            Text(StressComponentType.vonMisesAvg.rawValue).tag(StressComponentType.vonMisesAvg)
                            Text(StressComponentType.vonMisesMax.rawValue).tag(StressComponentType.vonMisesMax)
                        }
                        Section(header: Text("Nanoindentation Metrics")) {
                            Text(StressComponentType.indenterLoad.rawValue).tag(StressComponentType.indenterLoad)
                            Text(StressComponentType.meanContactPressure.rawValue).tag(StressComponentType.meanContactPressure)
                            Text(StressComponentType.contactAtoms.rawValue).tag(StressComponentType.contactAtoms)
                        }
                    }
                    .pickerStyle(.menu)
                    .labelsHidden()
                    .controlSize(.small)
                    .frame(width: 195)
                    .help("Select which stress component (Cauchy, Aligned Piola / Projected, Von Mises) or nanoindentation metric to monitor live")
                }
                
                // 2. Pinned Clever Auto-Scaling Selector
                Picker("", selection: $scaleMode) {
                    ForEach(ChartScaleMode.allCases) { mode in
                        Text(mode.rawValue).tag(mode)
                    }
                }
                .pickerStyle(.segmented)
                .labelsHidden()
                .controlSize(.small)
                .fixedSize()
                .help(scaleMode.tooltip)
                
                Divider()
                    .frame(height: 14)
                
                // 3. Horizontally scrollable toggles & badges so left pickers never shift or clip
                ScrollView(.horizontal, showsIndicators: false) {
                    HStack(spacing: 10) {
                        if chartMode == .gbEnergy {
                            Toggle("Interior Gauge (exclude 4a₀ collar)", isOn: $useInteriorGBGauge)
                                .toggleStyle(.checkbox)
                                .controlSize(.small)
                                .lineLimit(1)
                                .fixedSize()
                                .help("When checked, uses E_int / L_int (excluding 4a₀ near the outer boundaries). Uncheck for total E / L_GB.")
                            Toggle("Σ Labels", isOn: $showGBSigmaLabels)
                                .toggleStyle(.checkbox)
                                .controlSize(.small)
                                .lineLimit(1)
                                .fixedSize()
                                .help("Show Σ cusp labels on the γ_GB(Δθ) curve")
                        } else if chartMode != .eigenvalue {
                            Toggle("Remesh (\(remeshEvents.count))", isOn: $showRemeshEvents)
                                .toggleStyle(.checkbox)
                                .controlSize(.small)
                                .lineLimit(1)
                                .fixedSize()
                                .help("Highlight Delaunay remeshing & plastic surgery avalanche events")
                        }
                        
                        if chartMode == .stress {
                            Toggle("Peak Yield", isOn: $showYieldLine)
                                .toggleStyle(.checkbox)
                                .controlSize(.small)
                                .lineLimit(1)
                                .fixedSize()
                                .help("Display peak stress horizontal guide line")
                        }
                        
                        // Live Grain Boundary Energy badge (always visible when a Bicrystal run produced γ_GB)
                        if let curGB = runner.currentGBEnergy {
                            let gVal = useInteriorGBGauge ? curGB.gbEnergyInterior : curGB.gbEnergyTotal
                            Button(action: { chartMode = .gbEnergy }) {
                                HStack(spacing: 4) {
                                    Image(systemName: "square.split.2x1.fill")
                                        .font(.system(size: 9))
                                        .foregroundColor(.green)
                                    Text("γ_GB: \(String(format: "%.5f", gVal)) (Δθ=\(String(format: "%.2f°", curGB.misorientationDeg)))")
                                        .font(.system(size: 10, weight: .bold, design: .monospaced))
                                        .foregroundColor(.green)
                                }
                                .padding(.horizontal, 6)
                                .padding(.vertical, 3)
                                .background(BlueHour.controlSurface)
                                .cornerRadius(4)
                            }
                            .buttonStyle(.plain)
                            .help("Grain Boundary Energy at α = 0. Click to view the cumulative γ_GB(Δθ) curve across misorientations.")
                            .lineLimit(1)
                            .fixedSize()
                        }
                        
                        let peakVal = currentPeakStress
                        if chartMode == .stress && peakVal > 0 {
                            HStack(spacing: 4) {
                                Image(systemName: "crown.fill")
                                    .font(.system(size: 9))
                                    .foregroundColor(.orange)
                                Text("|\(stressComponent.shortLabel)|_max: \(String(format: "%.4f", peakVal))")
                                    .font(.system(size: 10, weight: .bold, design: .monospaced))
                                    .foregroundColor(.orange)
                            }
                            .padding(.horizontal, 6)
                            .padding(.vertical, 3)
                            .background(BlueHour.controlSurface)
                            .cornerRadius(4)
                            .lineLimit(1)
                            .fixedSize()
                        } else if chartMode == .eigenvalue, let lastEig = runner.eigenPoints.last {
                            HStack(spacing: 4) {
                                Circle()
                                    .fill(lastEig.isUnstable ? Color.red : Color.green)
                                    .frame(width: 6, height: 6)
                                Text("λ_min: \(String(format: "%.4e", lastEig.lambdaMin))")
                                    .font(.system(size: 10, weight: .bold, design: .monospaced))
                                    .foregroundColor(lastEig.isUnstable ? .red : .green)
                            }
                            .padding(.horizontal, 6)
                            .padding(.vertical, 3)
                            .background(BlueHour.controlSurface)
                            .cornerRadius(4)
                            .lineLimit(1)
                            .fixedSize()
                        }
                        
                        let ptCount: Int = {
                            switch chartMode {
                            case .eigenvalue: return runner.eigenPoints.count
                            case .gbEnergy: return displayGBPoints.count
                            default: return runner.curvePoints.count
                            }
                        }()
                        Text("\(ptCount) pts")
                            .font(.system(size: 10, weight: .medium, design: .monospaced))
                            .foregroundColor(.secondary)
                            .padding(.horizontal, 6)
                            .padding(.vertical, 3)
                            .background(BlueHour.controlSurface)
                            .cornerRadius(4)
                            .lineLimit(1)
                            .fixedSize()
                    }
                }
                .frame(maxWidth: .infinity, alignment: .leading)
                
                // Active Strain Range Badge
                if chartMode != .gbEnergy, let firstA = runner.curvePoints.first?.alpha, let lastA = runner.curvePoints.last?.alpha {
                    Text("α: [\(String(format: "%.3f", firstA)) → \(String(format: "%.3f", lastA))]")
                        .font(.system(size: 10, weight: .medium, design: .monospaced))
                        .foregroundColor(.secondary)
                        .padding(.horizontal, 6)
                        .padding(.vertical, 3)
                        .background(BlueHour.controlSurface)
                        .cornerRadius(4)
                        .lineLimit(1)
                        .fixedSize()
                }
                
                // Export Data Menu
                if !runner.curvePoints.isEmpty || !displayGBPoints.isEmpty {
                    Menu {
                        if !runner.curvePoints.isEmpty {
                            Button(action: exportCurveCSV) {
                                Label("Save Stress–Strain Curve (CSV)...", systemImage: "chart.xyaxis.line")
                            }
                            if !remeshEvents.isEmpty {
                                Button(action: exportAvalanchesCSV) {
                                    Label("Save Avalanches / Remesh Events (CSV)...", systemImage: "burst")
                                }
                            }
                            Button(action: copyCurveCSVToClipboard) {
                                Label("Copy Curve Data (CSV)", systemImage: "doc.on.doc")
                            }
                        }
                        if !displayGBPoints.isEmpty {
                            if !runner.curvePoints.isEmpty {
                                Divider()
                            }
                            Button(action: exportGBEnergyText) {
                                Label("Save Angle vs. Energy Table (.txt)...", systemImage: "doc.text")
                            }
                            Button(action: exportGBEnergyCSV) {
                                Label("Save γ_GB vs Misorientation Δθ (CSV)...", systemImage: "square.split.2x1")
                            }
                            Button(action: copyGBEnergyCSVToClipboard) {
                                Label("Copy γ_GB(Δθ) Table (CSV)", systemImage: "doc.on.doc")
                            }
                            Button(action: { runner.openSweepTextFile() }) {
                                Label("Open angle_vs_energy.txt", systemImage: "arrow.up.forward.app")
                            }
                            Button(action: { runner.openSweepRelaxedVTKFolder() }) {
                                Label("Open Relaxed VTK Configurations Folder", systemImage: "folder")
                            }
                            Button(role: .destructive, action: { runner.clearGBEnergyHistory() }) {
                                Label("Clear Cumulative γ_GB(Δθ) History", systemImage: "trash")
                            }
                        }
                    } label: {
                        Label("Export Data", systemImage: "square.and.arrow.up")
                    }
                    .menuStyle(.borderedButton)
                    .controlSize(.small)
                }
            }
            .frame(maxWidth: .infinity, alignment: .leading)
            .padding(.horizontal, 12)
            .padding(.vertical, 6)
            .background(BlueHour.controlSurface.opacity(0.85))
            
            Divider()
            
            // Main Chart Canvas
            if chartMode == .gbEnergy && displayGBPoints.isEmpty {
                VStack(spacing: 12) {
                    Image(systemName: "square.split.2x1")
                        .font(.system(size: 44))
                        .foregroundColor(.secondary.opacity(0.5))
                    Text("No Grain Boundary Energy Data Yet")
                        .font(.headline)
                        .foregroundColor(.secondary)
                    Text("Run a Bicrystal simulation (at α = 0) to compute the specific grain boundary energy γ_GB = E_GB / L_GB. Each angle Δθ you run is automatically appended to gb_energy_vs_misorientation.csv and plotted here.")
                        .font(.subheadline)
                        .foregroundColor(.secondary.opacity(0.8))
                        .multilineTextAlignment(.center)
                        .frame(maxWidth: 460)
                }
                .frame(maxWidth: .infinity, maxHeight: .infinity)
                .background(BlueHour.canvasBackground)
            } else if chartMode == .eigenvalue && runner.eigenPoints.isEmpty {
                VStack(spacing: 12) {
                    Image(systemName: "waveform.path.ecg")
                        .font(.system(size: 44))
                        .foregroundColor(.secondary.opacity(0.5))
                    Text("No Eigenvalue Data Recorded Yet")
                        .font(.headline)
                        .foregroundColor(.secondary)
                    Text(runner.activeConfig?.enableEigenvalues == true || runner.status == .idle
                         ? "Enable 'Compute Lowest Stiffness Eigenvalues' in the Stability & Eigenvalues sidebar section and launch the simulation to plot λ_min(α) live."
                         : "Eigenvalue computation was disabled for this run. Enable 'Compute Lowest Stiffness Eigenvalues' in the sidebar for the next run.")
                        .font(.subheadline)
                        .foregroundColor(.secondary.opacity(0.8))
                        .multilineTextAlignment(.center)
                        .frame(maxWidth: 440)
                }
                .frame(maxWidth: .infinity, maxHeight: .infinity)
                .background(BlueHour.canvasBackground)
            } else if chartMode != .gbEnergy && runner.curvePoints.isEmpty {
                VStack(spacing: 14) {
                    Image(systemName: "chart.xyaxis.line")
                        .font(.system(size: 48))
                        .foregroundColor(.secondary.opacity(0.5))
                    Text("No Simulation Data Yet")
                        .font(.headline)
                        .foregroundColor(.secondary)
                    Text("Click 'Launch Simulation' to stream the live stress–strain, energy, and lowest eigenvalue curves in real-time.")
                        .font(.subheadline)
                        .foregroundColor(.secondary.opacity(0.8))
                        .multilineTextAlignment(.center)
                        .frame(maxWidth: 420)
                }
                .frame(maxWidth: .infinity, maxHeight: .infinity)
                .background(BlueHour.canvasBackground)
            } else {
                Chart {
                    if chartMode == .gbEnergy {
                        ForEach(displayGBPoints) { pt in
                            let yVal = gbEnergyValue(for: pt)
                            LineMark(
                                x: .value("Misorientation Δθ (°)", pt.misorientationDeg),
                                y: .value("Grain Boundary Energy γ_GB", yVal)
                            )
                            .foregroundStyle(Color.green)
                            .interpolationMethod(.linear)
                            
                            PointMark(
                                x: .value("Misorientation Δθ (°)", pt.misorientationDeg),
                                y: .value("Grain Boundary Energy γ_GB", yVal)
                            )
                            .symbol(Circle())
                            .symbolSize(displayGBPoints.count > 25 ? 36 : 52)
                            .foregroundStyle(Color.green)
                            .annotation(position: .top) {
                                let isMajorCusp = (pt.sigma > 1 && pt.sigma <= 17)
                                if showGBSigmaLabels && (displayGBPoints.count <= 22 || isMajorCusp) {
                                    Text(pt.sigma > 0
                                         ? String(format: "Σ%d (%.1f°)", pt.sigma, pt.misorientationDeg)
                                         : String(format: "%.1f°", pt.misorientationDeg))
                                        .font(.system(size: 9, weight: .medium, design: .monospaced))
                                        .foregroundColor(.secondary)
                                }
                            }
                        }
                        
                        // Highlight current run's point if present
                        if let cur = runner.currentGBEnergy {
                            let curY = gbEnergyValue(for: cur)
                            PointMark(
                                x: .value("Current Δθ (°)", cur.misorientationDeg),
                                y: .value("Current γ_GB", curY)
                            )
                            .symbol(Circle().strokeBorder(lineWidth: 2.0))
                            .symbolSize(105)
                            .foregroundStyle(Color.orange)
                        }
                        
                        if let selAngle = selectedAlpha, let gbPt = findClosestGBPoint(to: selAngle) {
                            RuleMark(x: .value("Selected Δθ", gbPt.misorientationDeg))
                                .foregroundStyle(Color.secondary.opacity(0.4))
                                .lineStyle(StrokeStyle(lineWidth: 1, dash: [2, 2]))
                            
                            PointMark(
                                x: .value("Selected Δθ", gbPt.misorientationDeg),
                                y: .value("γ_GB", gbEnergyValue(for: gbPt))
                            )
                            .symbolSize(90)
                            .foregroundStyle(Color.cyan)
                        }
                    } else if chartMode == .eigenvalue {
                        // Zero-eigenvalue stability threshold line
                        RuleMark(y: .value("Stability Threshold", 0.0))
                            .foregroundStyle(Color.red.opacity(0.65))
                            .lineStyle(StrokeStyle(lineWidth: 1, dash: [4, 4]))
                            .annotation(position: .top, alignment: .leading) {
                                Text("λ = 0 (Instability Threshold)")
                                    .font(.system(size: 9, weight: .semibold, design: .monospaced))
                                    .foregroundColor(.red.opacity(0.85))
                                    .padding(2)
                            }
                        
                        ForEach(displayEigenPoints) { pt in
                            LineMark(
                                x: .value("Shear Strain (α)", pt.alpha),
                                y: .value("Smallest Eigenvalue (λ_min)", pt.lambdaMin)
                            )
                            .foregroundStyle(Color.teal)
                            .interpolationMethod(.linear)
                            
                            PointMark(
                                x: .value("Shear Strain (α)", pt.alpha),
                                y: .value("Smallest Eigenvalue (λ_min)", pt.lambdaMin)
                            )
                            .symbol(Circle())
                            .symbolSize(pt.isUnstable ? 55 : 32)
                            .foregroundStyle(pt.isUnstable ? Color.red : Color.teal)
                        }
                        
                        if let selAlpha = selectedAlpha, let eigPt = findClosestEigenPoint(to: selAlpha) {
                            RuleMark(x: .value("Selected α", eigPt.alpha))
                                .foregroundStyle(Color.secondary.opacity(0.4))
                                .lineStyle(StrokeStyle(lineWidth: 1, dash: [2, 2]))
                            
                            PointMark(
                                x: .value("Selected α", eigPt.alpha),
                                y: .value("Smallest Eigenvalue (λ_min)", eigPt.lambdaMin)
                            )
                            .symbolSize(85)
                            .foregroundStyle(Color.cyan)
                        }
                    } else {
                        // 1. Primary Stress / Energy Line Mark (+ PointMark when sparse)
                        ForEach(displayCurvePoints) { pt in
                            let yVal = chartMode == .stress ? stressValue(for: pt) : pt.postEnergy
                            let yLabel = chartMode == .stress ? ((stressComponent == .primary && primaryIsTorque) ? "Rim Torque T" : stressComponent.rawValue) : "Internal Energy (E)"
                            LineMark(
                                x: .value("Shear Strain (α)", pt.alpha),
                                y: .value(yLabel, yVal)
                            )
                            .foregroundStyle(chartMode == .stress ? Color.blue : Color.purple)
                            .interpolationMethod(.linear)
                            
                            if runner.curvePoints.count <= 50 {
                                PointMark(
                                    x: .value("Shear Strain (α)", pt.alpha),
                                    y: .value(yLabel, yVal)
                                )
                                .symbol(Circle())
                                .symbolSize(30)
                                .foregroundStyle(chartMode == .stress ? Color.blue : Color.purple)
                            }
                        }
                        
                        // 2. Plasticity / Remesh Events Overlay
                        if showRemeshEvents {
                            ForEach(displayRemeshEvents) { pt in
                                let yVal = chartMode == .stress ? stressValue(for: pt) : pt.postEnergy
                                let yLabel = chartMode == .stress ? ((stressComponent == .primary && primaryIsTorque) ? "Rim Torque T" : stressComponent.rawValue) : "Internal Energy (E)"
                                PointMark(
                                    x: .value("Shear Strain (α)", pt.alpha),
                                    y: .value(yLabel, yVal)
                                )
                                .symbol(Circle())
                                .symbolSize(48)
                                .foregroundStyle(Color.red)
                                
                                PointMark(
                                    x: .value("Shear Strain (α)", pt.alpha),
                                    y: .value(yLabel, yVal)
                                )
                                .symbol(Circle().strokeBorder(lineWidth: 1.5))
                                .symbolSize(80)
                                .foregroundStyle(Color.red.opacity(0.85))
                            }
                        }
                        
                        // 3. Peak / Yield Stress Guide Line
                        let peakVal = currentPeakStress
                        if showYieldLine && chartMode == .stress && peakVal > 0 {
                            let lastStress = runner.curvePoints.last.map { stressValue(for: $0) } ?? 0.0
                            let peakY = lastStress < 0 ? -peakVal : peakVal
                            RuleMark(
                                y: .value("Peak Stress", peakY)
                            )
                            .foregroundStyle(Color.orange.opacity(0.7))
                            .lineStyle(StrokeStyle(lineWidth: 1, dash: [4, 4]))
                            .annotation(position: .top, alignment: .trailing) {
                                Text("|\(stressComponent.shortLabel)|_max: \(String(format: "%.4f", peakVal))")
                                    .font(.system(size: 9, weight: .bold, design: .monospaced))
                                    .foregroundColor(.orange)
                                    .padding(2)
                                    .background(BlueHour.windowBackground.opacity(0.8))
                            }
                        }
                        
                        // 4. Hover Cursor / Selection Indicator
                        if let selAlpha = selectedAlpha, let pt = findClosestPoint(to: selAlpha) {
                            RuleMark(x: .value("Selected α", pt.alpha))
                                .foregroundStyle(Color.secondary.opacity(0.4))
                                .lineStyle(StrokeStyle(lineWidth: 1, dash: [2, 2]))
                            
                            PointMark(
                                x: .value("Selected α", pt.alpha),
                                y: .value("Value", chartMode == .stress ? stressValue(for: pt) : pt.postEnergy)
                            )
                            .symbolSize(80)
                            .foregroundStyle(Color.cyan)
                        }
                    }
                }
                .chartXScale(domain: xDomain)
                .chartYScale(domain: yDomain)
                .chartXAxis {
                    AxisMarks(values: .automatic(desiredCount: 7)) { value in
                        AxisGridLine(stroke: StrokeStyle(lineWidth: 0.5, dash: [2, 2]))
                            .foregroundStyle(Color.secondary.opacity(0.25))
                        AxisTick()
                        if let d = value.as(Double.self) {
                            AxisValueLabel {
                                Text(formatXValue(d))
                                    .font(.system(size: 10, design: .monospaced))
                                    .foregroundColor(.secondary)
                            }
                        }
                    }
                }
                .chartYAxis {
                    AxisMarks(values: .automatic(desiredCount: 6)) { value in
                        AxisGridLine(stroke: StrokeStyle(lineWidth: 0.5, dash: [2, 2]))
                            .foregroundStyle(Color.secondary.opacity(0.25))
                        AxisTick()
                        if let d = value.as(Double.self) {
                            AxisValueLabel {
                                Text(formatYValue(d))
                                    .font(.system(size: 10, design: .monospaced))
                                    .foregroundColor(.secondary)
                            }
                        }
                    }
                }
                .chartOverlay { proxy in
                    GeometryReader { geo in
                        Rectangle()
                            .fill(Color.clear)
                            .contentShape(Rectangle())
                            .onContinuousHover { phase in
                                switch phase {
                                case .active(let location):
                                    if let x: Double = proxy.value(atX: location.x) {
                                        selectedAlpha = x
                                    }
                                case .ended:
                                    selectedAlpha = nil
                                }
                            }
                    }
                }
                .padding(.horizontal, 16)
                .padding(.top, 16)
                .padding(.bottom, 8)
                .background(BlueHour.canvasBackground)
                
                // Hovered point details bar
                if chartMode == .gbEnergy, let gbPt = (selectedAlpha.flatMap { findClosestGBPoint(to: $0) } ?? runner.currentGBEnergy ?? displayGBPoints.last) {
                    HStack(spacing: 12) {
                        Text(gbPt.sigma > 0
                             ? String(format: "Δθ: %.2f° (Σ%d, %d %d 0)", gbPt.misorientationDeg, gbPt.sigma, gbPt.m, gbPt.n)
                             : String(format: "Δθ: %.3f° (θ = ±%.3f°)", gbPt.misorientationDeg, gbPt.halfAngleDeg))
                            .font(.system(size: 11, weight: .semibold, design: .monospaced))
                        Text("γ_GB (Int): \(String(format: "%.6f", gbPt.gbEnergyInterior))")
                            .font(.system(size: 11, weight: .bold, design: .monospaced))
                            .foregroundColor(.green)
                        Text("γ_GB (Tot): \(String(format: "%.6f", gbPt.gbEnergyTotal))")
                            .font(.system(size: 11, design: .monospaced))
                            .foregroundColor(.teal)
                        Text("L_GB: \(String(format: "%.0f", gbPt.lengthGBA0))a₀ (\(gbPt.nx)×\(gbPt.ny))")
                            .font(.system(size: 10, design: .monospaced))
                            .foregroundColor(.secondary)
                        
                        Spacer()
                        
                        if !gbPt.relaxedVTKPath.isEmpty && FileManager.default.fileExists(atPath: gbPt.relaxedVTKPath) {
                            Button(action: {
                                NSWorkspace.shared.open(URL(fileURLWithPath: gbPt.relaxedVTKPath))
                            }) {
                                Label("Open Angle VTK", systemImage: "cube.transparent")
                                    .font(.system(size: 10, weight: .semibold))
                            }
                            .buttonStyle(.bordered)
                            .controlSize(.mini)
                            .help("Open relaxed configuration VTK for Δθ = \(String(format: "%.2f°", gbPt.misorientationDeg)) in ParaView")
                        }
                        
                        Button(action: { runner.openSweepTextFile() }) {
                            Label("angle_vs_energy.txt", systemImage: "doc.text")
                                .font(.system(size: 10, weight: .semibold))
                        }
                        .buttonStyle(.bordered)
                        .controlSize(.mini)
                        
                        Button(action: { runner.openSweepRelaxedVTKFolder() }) {
                            Label("VTK Folder", systemImage: "folder")
                                .font(.system(size: 10, weight: .semibold))
                        }
                        .buttonStyle(.bordered)
                        .controlSize(.mini)
                    }
                    .padding(.horizontal, 14)
                    .padding(.vertical, 6)
                    .background(BlueHour.controlSurface)
                } else if chartMode == .eigenvalue, let selAlpha = selectedAlpha, let eigPt = findClosestEigenPoint(to: selAlpha) {
                    HStack(spacing: 18) {
                        Text("Step: \(eigPt.iteration)")
                            .font(.system(size: 11, weight: .semibold, design: .monospaced))
                        Text("Strain α: \(String(format: "%.6f", eigPt.alpha))")
                            .font(.system(size: 11, design: .monospaced))
                        Text("λ₁ (min): \(String(format: "%.6e", eigPt.lambdaMin))")
                            .font(.system(size: 11, weight: .semibold, design: .monospaced))
                            .foregroundColor(eigPt.isUnstable ? .red : .teal)
                        if let l2 = eigPt.lambda2 {
                            Text("λ₂: \(String(format: "%.6e", l2))")
                                .font(.system(size: 11, design: .monospaced))
                                .foregroundColor(.secondary)
                        }
                        Text("Trigger: \(eigPt.trigger)")
                            .font(.system(size: 10, design: .monospaced))
                            .foregroundColor(.secondary)
                    }
                    .padding(.horizontal, 14)
                    .padding(.vertical, 6)
                    .background(BlueHour.controlSurface)
                } else if let selAlpha = selectedAlpha, let pt = findClosestPoint(to: selAlpha) {
                    let isIndentation = (runner.activeConfig?.geometryConfig.indenter.enabled == true ||
                                         runner.activeConfig?.geometryConfig.primary_stress_component == "indentation_load")
                    let curStress = stressValue(for: pt)
                    let stressTitle: String = {
                        if stressComponent == .primary {
                            return isIndentation ? "Load P" : (primaryIsTorque ? "Torque T" : "Stress σ")
                        }
                        return stressComponent.shortLabel
                    }()
                    HStack(spacing: 18) {
                        Text("Step: \(pt.iteration)")
                            .font(.system(size: 11, weight: .semibold, design: .monospaced))
                        Text(isIndentation ? "Depth δ: \(String(format: "%.4f", pt.alpha))" : "Strain α: \(String(format: "%.6f", pt.alpha))")
                            .font(.system(size: 11, design: .monospaced))
                        Text("\(stressTitle): \(String(format: "%.6f", curStress))")
                            .font(.system(size: 11, weight: .semibold, design: .monospaced))
                            .foregroundColor(.blue)
                        if stressComponent != .vonMisesAvg {
                            Text("⟨σ_vM⟩: \(String(format: "%.4f", pt.vonMisesAvg))")
                                .font(.system(size: 11, design: .monospaced))
                                .foregroundColor(.orange)
                        }
                        Text("Energy E: \(String(format: "%.6e", pt.postEnergy))")
                            .font(.system(size: 11, design: .monospaced))
                            .foregroundColor(.purple)
                        if pt.shouldRemesh && pt.iteration > 0 {
                            Text("⚡ Delaunay Surgery Event")
                                .font(.system(size: 10, weight: .bold))
                                .foregroundColor(.red)
                        }
                    }
                    .padding(.horizontal, 14)
                    .padding(.vertical, 6)
                    .background(BlueHour.controlSurface)
                }
            }
        }
        .onAppear {
            if let req = runner.requestedChartMode {
                chartMode = req
            } else if runner.isSweepRunning || runner.activeConfig?.protocolMode == .sweepMisorientation || runner.activeConfig?.protocolMode == .relaxOnly {
                chartMode = .gbEnergy
            }
        }
        .onReceive(runner.$requestedChartMode) { newMode in
            if let m = newMode {
                chartMode = m
            }
        }
    }
    
    private var displayCurvePoints: [StressStrainPoint] {
        let pts = runner.curvePoints
        guard pts.count > 1000 else { return pts }
        let stride = max(1, pts.count / 800)
        var sampled: [StressStrainPoint] = []
        sampled.reserveCapacity(800 + 2)
        for i in Swift.stride(from: 0, to: pts.count, by: stride) {
            sampled.append(pts[i])
        }
        if let last = pts.last, sampled.last?.iteration != last.iteration {
            sampled.append(last)
        }
        return sampled
    }
    
    private var displayRemeshEvents: [StressStrainPoint] {
        let events = remeshEvents
        guard events.count > 400 else { return events }
        let stride = max(1, events.count / 300)
        var sampled: [StressStrainPoint] = []
        for i in Swift.stride(from: 0, to: events.count, by: stride) {
            sampled.append(events[i])
        }
        return sampled
    }
    
    private var displayEigenPoints: [EigenPoint] {
        let pts = runner.eigenPoints
        guard pts.count > 1000 else { return pts }
        let stride = max(1, pts.count / 800)
        var sampled: [EigenPoint] = []
        for i in Swift.stride(from: 0, to: pts.count, by: stride) {
            sampled.append(pts[i])
        }
        if let last = pts.last, sampled.last?.iteration != last.iteration {
            sampled.append(last)
        }
        return sampled
    }
    
    private func findClosestPoint(to alpha: Double) -> StressStrainPoint? {
        guard !runner.curvePoints.isEmpty else { return nil }
        return runner.curvePoints.min(by: { abs($0.alpha - alpha) < abs($1.alpha - alpha) })
    }
    
    private func findClosestEigenPoint(to alpha: Double) -> EigenPoint? {
        guard !runner.eigenPoints.isEmpty else { return nil }
        return runner.eigenPoints.min(by: { abs($0.alpha - alpha) < abs($1.alpha - alpha) })
    }
    
    private func findClosestGBPoint(to angle: Double) -> GBEnergyPoint? {
        let pts = displayGBPoints
        guard !pts.isEmpty else { return nil }
        return pts.min(by: { abs($0.misorientationDeg - angle) < abs($1.misorientationDeg - angle) })
    }
    
    // MARK: - Export Methods
    
    private func generateCurveCSV() -> String {
        var csv = "step,alpha_strain,stress_sigma,energy_E,pre_stress,pre_energy,is_remesh,cauchy_xx,cauchy_xy,cauchy_yy,hydrostatic_pressure,piola_11,piola_12,piola_22,von_mises_avg,von_mises_max,indenter_load,mean_contact_pressure,contact_atoms\n"
        for pt in runner.curvePoints {
            csv += "\(pt.iteration),\(pt.alpha),\(pt.postStress),\(pt.postEnergy),\(pt.preStress),\(pt.preEnergy),\(pt.shouldRemesh),\(pt.cauchyXX),\(pt.cauchyXY),\(pt.cauchyYY),\(pt.hydrostaticPressure),\(pt.piola11),\(pt.piola12),\(pt.piola22),\(pt.vonMisesAvg),\(pt.vonMisesMax),\(pt.indenterLoad),\(pt.meanContactPressure),\(pt.contactAtoms)\n"
        }
        return csv
    }
    
    private func generateGBEnergyCSV() -> String {
        var csv = "MisorientationDeg,HalfAngleDeg,Theta1Deg,Theta2Deg,GBAngleDeg,Nx,Ny,L_GB_a0,L_Int_a0,TotalEnergy,InteriorEnergy,GBEnergyTotal,GBEnergyInterior,NormalGapA0,SeamDMinA0,Sigma,m,n,RelaxedVTKPath\n"
        for pt in displayGBPoints {
            csv += "\(pt.misorientationDeg),\(pt.halfAngleDeg),\(pt.theta1Deg),\(pt.theta2Deg),\(pt.gbAngleDeg),\(pt.nx),\(pt.ny),\(pt.lengthGBA0),\(pt.lengthIntA0),\(pt.totalEnergy),\(pt.interiorEnergy),\(pt.gbEnergyTotal),\(pt.gbEnergyInterior),\(pt.normalGapA0),\(pt.seamDMinA0),\(pt.sigma),\(pt.m),\(pt.n),\(pt.relaxedVTKPath)\n"
        }
        return csv
    }
    
    private func generateGBEnergyText() -> String {
        var txt = "# Symmetric Tilt Grain Boundary Energy vs. Misorientation Angle (alpha = 0)\n"
        txt += "#  1:DeltaTheta_deg  2:HalfAngle_deg  3:Sigma  4:m  5:n  6:Gamma_GB_interior  7:Gamma_GB_total  8:E_total_relaxed  9:E_interior_relaxed  10:L_GB_a0  11:NX  12:NY  13:Gap_a0  14:d_min_a0  15:Relaxed_VTK_File\n"
        for p in displayGBPoints {
            let vtkName = p.relaxedVTKPath.isEmpty ? "configuration_0_post.vtk" : (p.relaxedVTKPath as NSString).lastPathComponent
            txt += String(
                format: "%10.4f %10.4f %6d %4d %4d %16.8f %16.8f %16.8f %16.8f %10.2f %5d %5d %9.4f %9.4f  %@\n",
                p.misorientationDeg, p.halfAngleDeg, p.sigma, p.m, p.n,
                p.gbEnergyInterior, p.gbEnergyTotal, p.totalEnergy, p.interiorEnergy,
                p.lengthGBA0, p.nx, p.ny, p.normalGapA0, p.seamDMinA0, vtkName
            )
        }
        return txt
    }
    
    private func exportCurveCSV() {
        let csv = generateCurveCSV()
        ExportDataHelper.saveStringToFile(content: csv, defaultFileName: "stress_strain_curve.csv")
    }
    
    private func exportGBEnergyText() {
        ExportDataHelper.saveStringToFile(content: generateGBEnergyText(), defaultFileName: "angle_vs_energy.txt")
    }
    
    private func exportGBEnergyCSV() {
        let csv = generateGBEnergyCSV()
        ExportDataHelper.saveStringToFile(content: csv, defaultFileName: "gb_energy_vs_misorientation.csv")
    }
    
    private func copyGBEnergyCSVToClipboard() {
        ExportDataHelper.copyToClipboard(text: generateGBEnergyCSV())
    }
    
    private func exportAvalanchesCSV() {
        var csv = "step,alpha_strain,stress_before,stress_after,stress_drop,energy_before,energy_after,energy_drop\n"
        for pt in remeshEvents {
            let sDrop = pt.preStress - pt.postStress
            let eDrop = pt.preEnergy - pt.postEnergy
            csv += "\(pt.iteration),\(pt.alpha),\(pt.preStress),\(pt.postStress),\(sDrop),\(pt.preEnergy),\(pt.postEnergy),\(eDrop)\n"
        }
        ExportDataHelper.saveStringToFile(content: csv, defaultFileName: "plastic_avalanches_remesh_events.csv")
    }
    
    private func copyCurveCSVToClipboard() {
        ExportDataHelper.copyToClipboard(text: generateCurveCSV())
    }
}
