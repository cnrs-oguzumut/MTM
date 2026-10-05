import SwiftUI
import Charts

public enum ChartMetricMode: String, CaseIterable, Identifiable {
    case stress = "Stress σ(α)"
    case energy = "Energy E(α)"
    
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
            return "Full Schedule: starts at α_start (0.14) up to α_end (1.0)"
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
    
    private var xDomain: ClosedRange<Double> {
        let firstAlpha = runner.curvePoints.first?.alpha ?? (runner.activeConfig?.alphaStart ?? 0.14)
        
        guard !runner.curvePoints.isEmpty else {
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
        
        let alphas = runner.curvePoints.map(\.alpha)
        let minA = alphas.min() ?? firstAlpha
        let maxA = alphas.max() ?? (minA + 0.001)
        let span = maxA - minA
        
        switch scaleMode {
        case .autoFit:
            // Clever tight auto-fit around the active data
            let effectiveSpan = max(0.0008, span)
            let padLeft = effectiveSpan * 0.025
            let padRight = effectiveSpan * 0.065 // extra room on right for incoming points
            let xMin = minA - padLeft
            let xMax = maxA + padRight
            return xMin <= xMax ? xMin...xMax : (minA - 0.001)...(minA + 0.005)
            
        case .fullSchedule:
            // Starts at simulation start strain (e.g. 0.14), NEVER at 0
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
        guard !runner.curvePoints.isEmpty else {
            return chartMode == .stress ? 0.0...0.05 : 0.0...100.0
        }
        
        let vals = runner.curvePoints.map { chartMode == .stress ? $0.postStress : $0.postEnergy }
        let minY = vals.min() ?? 0.0
        let maxY = vals.max() ?? 0.01
        let span = maxY - minY
        
        switch chartMode {
        case .stress:
            if span < 1e-6 {
                let pad = max(abs(minY) * 0.1, 0.002)
                let yMin = max(0.0, minY - pad)
                return yMin...(maxY + pad)
            }
            
            let pad = span * 0.08
            let yMin: Double
            if scaleMode == .fromZero || (minY >= 0 && minY < span * 0.3) {
                // If baseline is naturally close to 0, anchor cleanly
                yMin = min(0.0, minY - pad * 0.5)
            } else {
                // Auto-scale tight to stress range
                yMin = minY - pad
            }
            
            var yMax = maxY + pad
            if showYieldLine && runner.peakStress > 0 {
                yMax = max(yMax, runner.peakStress + pad * 0.5)
            }
            return yMin <= yMax ? yMin...yMax : 0.0...0.05
            
        case .energy:
            // Energy levels (e.g. ~51.8) MUST be tightly zoomed; otherwise it looks completely flat
            if scaleMode == .fromZero {
                return 0.0...(maxY * 1.1)
            }
            
            if span < 1e-6 {
                let pad = max(abs(minY) * 0.02, 0.01)
                return (minY - pad)...(maxY + pad)
            }
            
            let pad = span * 0.08
            let yMin = minY - pad
            let yMax = maxY + pad
            return yMin <= yMax ? yMin...yMax : (minY - 0.1)...(maxY + 0.1)
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
        case .stress:
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
            // Chart Control & Metrics Sub-header with clean single-line alignment
            HStack(spacing: 10) {
                // 1. Metric Selector: Stress vs Energy
                Picker("", selection: $chartMode) {
                    ForEach(ChartMetricMode.allCases) { mode in
                        Text(mode.rawValue).tag(mode)
                    }
                }
                .pickerStyle(.segmented)
                .labelsHidden()
                .frame(width: 175)
                .help("Select plotted quantity: Shear Stress σ(α) or Internal Energy E(α)")
                
                // 2. Clever Auto-Scaling Selector
                Picker("", selection: $scaleMode) {
                    ForEach(ChartScaleMode.allCases) { mode in
                        Text(mode.rawValue).tag(mode)
                    }
                }
                .pickerStyle(.segmented)
                .labelsHidden()
                .frame(width: 195)
                .help(scaleMode.tooltip)
                
                Divider()
                    .frame(height: 16)
                
                // 3. Toggles with live remesh count and single-line protection
                Toggle("Remesh (\(remeshEvents.count))", isOn: $showRemeshEvents)
                    .toggleStyle(.checkbox)
                    .controlSize(.small)
                    .lineLimit(1)
                    .fixedSize()
                    .help("Highlight Delaunay remeshing & plastic surgery avalanche events")
                
                Toggle("Peak Yield", isOn: $showYieldLine)
                    .toggleStyle(.checkbox)
                    .controlSize(.small)
                    .lineLimit(1)
                    .fixedSize()
                    .help("Display peak stress horizontal guide line")
                
                Spacer(minLength: 8)
                
                // 4. Yield Stats Badge (Single un-wrapped line)
                if runner.peakStress > 0 {
                    HStack(spacing: 6) {
                        Image(systemName: "crown.fill")
                            .font(.system(size: 9))
                            .foregroundColor(.orange)
                        Text("σ_max: \(String(format: "%.4f", runner.peakStress))")
                            .font(.system(size: 10, weight: .bold, design: .monospaced))
                            .foregroundColor(.orange)
                        Text("α: \(String(format: "%.4f", runner.yieldStrain))")
                            .font(.system(size: 10, weight: .medium, design: .monospaced))
                            .foregroundColor(.secondary)
                    }
                    .padding(.horizontal, 7)
                    .padding(.vertical, 3)
                    .background(Color(NSColor.controlBackgroundColor))
                    .cornerRadius(5)
                    .lineLimit(1)
                    .fixedSize()
                }
                
                // 5. Active Strain Range Badge
                if let firstA = runner.curvePoints.first?.alpha, let lastA = runner.curvePoints.last?.alpha {
                    Text("α: [\(String(format: "%.3f", firstA)) → \(String(format: "%.3f", lastA))]")
                        .font(.system(size: 10, weight: .medium, design: .monospaced))
                        .foregroundColor(.secondary)
                        .padding(.horizontal, 6)
                        .padding(.vertical, 3)
                        .background(Color(NSColor.controlBackgroundColor))
                        .cornerRadius(4)
                        .lineLimit(1)
                        .fixedSize()
                }
                
                // 6. Point Counter
                Text("\(runner.curvePoints.count) pts")
                    .font(.system(size: 10, weight: .medium, design: .monospaced))
                    .foregroundColor(.secondary)
                    .padding(.horizontal, 6)
                    .padding(.vertical, 3)
                    .background(Color(NSColor.controlBackgroundColor))
                    .cornerRadius(4)
                    .lineLimit(1)
                    .fixedSize()
            }
            .padding(.horizontal, 12)
            .padding(.vertical, 6)
            .background(Color(NSColor.controlBackgroundColor).opacity(0.85))
            
            Divider()
            
            // Main Chart Canvas
            if runner.curvePoints.isEmpty {
                VStack(spacing: 14) {
                    Image(systemName: "chart.xyaxis.line")
                        .font(.system(size: 48))
                        .foregroundColor(.secondary.opacity(0.5))
                    Text("No Simulation Data Yet")
                        .font(.headline)
                        .foregroundColor(.secondary)
                    Text("Click 'Launch Simulation' to stream the live stress–strain curve in real-time. Clever auto-scaling dynamically fits your data.")
                        .font(.subheadline)
                        .foregroundColor(.secondary.opacity(0.8))
                        .multilineTextAlignment(.center)
                        .frame(maxWidth: 420)
                }
                .frame(maxWidth: .infinity, maxHeight: .infinity)
                .background(Color(NSColor.textBackgroundColor))
            } else {
                Chart {
                    // 1. Primary Line Mark
                    ForEach(runner.curvePoints) { pt in
                        LineMark(
                            x: .value("Shear Strain (α)", pt.alpha),
                            y: .value(chartMode == .stress ? "Shear Stress (σ)" : "Internal Energy (E)",
                                     chartMode == .stress ? pt.postStress : pt.postEnergy)
                        )
                        .foregroundStyle(chartMode == .stress ? Color.blue : Color.purple)
                        .interpolationMethod(.linear)
                    }
                    
                    // 2. Plasticity / Remesh Events Overlay (Concentric Target Ring + Solid Core)
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
                        RuleMark(
                            y: .value("Peak Stress", runner.peakStress)
                        )
                        .foregroundStyle(Color.orange.opacity(0.7))
                        .lineStyle(StrokeStyle(lineWidth: 1, dash: [4, 4]))
                        .annotation(position: .top, alignment: .trailing) {
                            Text("σ_max: \(String(format: "%.4f", runner.peakStress))")
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
                if let selAlpha = selectedAlpha, let pt = findClosestPoint(to: selAlpha) {
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
}
