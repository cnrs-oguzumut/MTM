import SwiftUI
import Charts

public enum ChartMetricMode: String, CaseIterable, Identifiable {
    case stress = "Stress σ(α)"
    case energy = "Energy E(α)"
    case eigenvalue = "λ_min(α)"
    
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
    @State private var scaleMode: ChartScaleMode = .autoFit
    @State private var showRemeshEvents: Bool = true
    @State private var showYieldLine: Bool = true
    
    // Remesh / surgery points (excluding initial setup mesh at iteration 0)
    private var remeshEvents: [StressStrainPoint] {
        runner.curvePoints.filter { $0.iteration > 0 && $0.shouldRemesh }
    }
    
    // MARK: - Dynamic Clever Domain Calculations
    
    private var activeAlphas: [Double] {
        if chartMode == .eigenvalue && !runner.eigenPoints.isEmpty {
            return runner.eigenPoints.map(\.alpha)
        }
        return runner.curvePoints.map(\.alpha)
    }
    
    private var xDomain: ClosedRange<Double> {
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
        case .stress:
            guard !runner.curvePoints.isEmpty else { return 0.0...0.05 }
            let vals = runner.curvePoints.map(\.postStress)
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
            if showYieldLine && runner.peakStress > 0 && maxY > 0 {
                yMax = max(yMax, runner.peakStress + pad * 0.5)
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
        case .stress, .eigenvalue:
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
                // 1. Pinned Metric Selector: Stress vs Energy vs Eigenvalue
                Picker("", selection: $chartMode) {
                    ForEach(ChartMetricMode.allCases) { mode in
                        Text(mode.rawValue).tag(mode)
                    }
                }
                .pickerStyle(.segmented)
                .labelsHidden()
                .controlSize(.small)
                .fixedSize()
                .help("Select plotted quantity: Shear Stress σ(α), Internal Energy E(α), or Smallest Stiffness Eigenvalue λ_min(α)")
                
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
                        if chartMode != .eigenvalue {
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
                        
                        if chartMode == .stress && runner.peakStress > 0 {
                            HStack(spacing: 4) {
                                Image(systemName: "crown.fill")
                                    .font(.system(size: 9))
                                    .foregroundColor(.orange)
                                Text("|σ|_max: \(String(format: "%.4f", runner.peakStress))")
                                    .font(.system(size: 10, weight: .bold, design: .monospaced))
                                    .foregroundColor(.orange)
                            }
                            .padding(.horizontal, 6)
                            .padding(.vertical, 3)
                            .background(Color(NSColor.controlBackgroundColor))
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
                            .background(Color(NSColor.controlBackgroundColor))
                            .cornerRadius(4)
                            .lineLimit(1)
                            .fixedSize()
                        }
                        
                        let ptCount = chartMode == .eigenvalue ? runner.eigenPoints.count : runner.curvePoints.count
                        Text("\(ptCount) pts")
                            .font(.system(size: 10, weight: .medium, design: .monospaced))
                            .foregroundColor(.secondary)
                            .padding(.horizontal, 6)
                            .padding(.vertical, 3)
                            .background(Color(NSColor.controlBackgroundColor))
                            .cornerRadius(4)
                            .lineLimit(1)
                            .fixedSize()
                    }
                    .padding(.vertical, 1)
                }
                .frame(maxWidth: .infinity, alignment: .leading)
            }
            .frame(maxWidth: .infinity, alignment: .leading)
            .padding(.horizontal, 12)
            .padding(.vertical, 6)
            .background(Color(NSColor.controlBackgroundColor).opacity(0.85))
            
            Divider()
            
            // Main Chart Canvas
            if chartMode == .eigenvalue && runner.eigenPoints.isEmpty {
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
                .background(Color(NSColor.textBackgroundColor))
            } else if runner.curvePoints.isEmpty {
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
                .background(Color(NSColor.textBackgroundColor))
            } else {
                Chart {
                    if chartMode == .eigenvalue {
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
                        
                        ForEach(runner.eigenPoints) { pt in
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
                        ForEach(runner.curvePoints) { pt in
                            LineMark(
                                x: .value("Shear Strain (α)", pt.alpha),
                                y: .value(chartMode == .stress ? "Shear Stress (σ)" : "Internal Energy (E)",
                                         chartMode == .stress ? pt.postStress : pt.postEnergy)
                            )
                            .foregroundStyle(chartMode == .stress ? Color.blue : Color.purple)
                            .interpolationMethod(.linear)
                            
                            if runner.curvePoints.count <= 50 {
                                PointMark(
                                    x: .value("Shear Strain (α)", pt.alpha),
                                    y: .value(chartMode == .stress ? "Shear Stress (σ)" : "Internal Energy (E)",
                                             chartMode == .stress ? pt.postStress : pt.postEnergy)
                                )
                                .symbol(Circle())
                                .symbolSize(30)
                                .foregroundStyle(chartMode == .stress ? Color.blue : Color.purple)
                            }
                        }
                        
                        // 2. Plasticity / Remesh Events Overlay
                        if showRemeshEvents {
                            ForEach(remeshEvents) { pt in
                                PointMark(
                                    x: .value("Shear Strain (α)", pt.alpha),
                                    y: .value(chartMode == .stress ? "Shear Stress (σ)" : "Internal Energy (E)",
                                             chartMode == .stress ? pt.postStress : pt.postEnergy)
                                )
                                .symbol(Circle())
                                .symbolSize(48)
                                .foregroundStyle(Color.red)
                                
                                PointMark(
                                    x: .value("Shear Strain (α)", pt.alpha),
                                    y: .value(chartMode == .stress ? "Shear Stress (σ)" : "Internal Energy (E)",
                                             chartMode == .stress ? pt.postStress : pt.postEnergy)
                                )
                                .symbol(Circle().strokeBorder(lineWidth: 1.5))
                                .symbolSize(80)
                                .foregroundStyle(Color.red.opacity(0.85))
                            }
                        }
                        
                        // 3. Peak / Yield Stress Guide Line
                        if showYieldLine && chartMode == .stress && runner.peakStress > 0 {
                            let peakY = (runner.curvePoints.last?.postStress ?? 0) < 0 ? -runner.peakStress : runner.peakStress
                            RuleMark(
                                y: .value("Peak Stress", peakY)
                            )
                            .foregroundStyle(Color.orange.opacity(0.7))
                            .lineStyle(StrokeStyle(lineWidth: 1, dash: [4, 4]))
                            .annotation(position: .top, alignment: .trailing) {
                                Text("|σ|_max: \(String(format: "%.4f", runner.peakStress))")
                                    .font(.system(size: 9, weight: .bold, design: .monospaced))
                                    .foregroundColor(.orange)
                                    .padding(2)
                                    .background(Color(NSColor.windowBackgroundColor).opacity(0.8))
                            }
                        }
                        
                        // 4. Hover Cursor / Selection Indicator
                        if let selAlpha = selectedAlpha, let pt = findClosestPoint(to: selAlpha) {
                            RuleMark(x: .value("Selected α", pt.alpha))
                                .foregroundStyle(Color.secondary.opacity(0.4))
                                .lineStyle(StrokeStyle(lineWidth: 1, dash: [2, 2]))
                            
                            PointMark(
                                x: .value("Selected α", pt.alpha),
                                y: .value("Value", chartMode == .stress ? pt.postStress : pt.postEnergy)
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
                .background(Color(NSColor.textBackgroundColor))
                
                // Hovered point details bar
                if chartMode == .eigenvalue, let selAlpha = selectedAlpha, let eigPt = findClosestEigenPoint(to: selAlpha) {
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
                    .background(Color(NSColor.controlBackgroundColor))
                } else if let selAlpha = selectedAlpha, let pt = findClosestPoint(to: selAlpha) {
                    HStack(spacing: 20) {
                        Text("Step: \(pt.iteration)")
                            .font(.system(size: 11, weight: .semibold, design: .monospaced))
                        Text("Strain α: \(String(format: "%.6f", pt.alpha))")
                            .font(.system(size: 11, design: .monospaced))
                        Text("Stress σ: \(String(format: "%.6f", pt.postStress))")
                            .font(.system(size: 11, design: .monospaced))
                            .foregroundColor(.blue)
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
                    .background(Color(NSColor.controlBackgroundColor))
                }
            }
        }
    }
    
    private func findClosestPoint(to alpha: Double) -> StressStrainPoint? {
        guard !runner.curvePoints.isEmpty else { return nil }
        return runner.curvePoints.min(by: { abs($0.alpha - alpha) < abs($1.alpha - alpha) })
    }
    
    private func findClosestEigenPoint(to alpha: Double) -> EigenPoint? {
        guard !runner.eigenPoints.isEmpty else { return nil }
        return runner.eigenPoints.min(by: { abs($0.alpha - alpha) < abs($1.alpha - alpha) })
    }
}
