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
    case fullSchedule = "Full"
    case fromZero = "From 0"
    
    public var id: String { rawValue }
    
    public var tooltip: String {
        switch self {
        case .autoFit:
            return "Clever Auto-Fit: dynamically zooms to active data range"
        case .fullSchedule:
            return "Full Schedule: spans the full configured loading schedule"
        case .fromZero:
            return "Origin at 0: anchors horizontal axis at 0.0"
        }
    }
}

public enum ChartPointSymbol: String, CaseIterable, Identifiable {
    case circle = "● Circle"
    case square = "■ Square"
    case diamond = "◆ Diamond"
    case triangle = "▲ Triangle"
    case cross = "✕ Cross"
    case plus = "+ Plus"
    
    public var id: String { rawValue }
    
    public var chartSymbol: BasicChartSymbolShape {
        switch self {
        case .circle: return .circle
        case .square: return .square
        case .diamond: return .diamond
        case .triangle: return .triangle
        case .cross: return .cross
        case .plus: return .plus
        }
    }
}

struct StressStrainChartView: View {
    @EnvironmentObject var runner: ProcessRunner
    @State private var selectedAlpha: Double?
    @State private var chartMode: ChartMetricMode = .stress
    @State private var stressComponent: StressComponentType = .primary
    @State private var scaleMode: ChartScaleMode = .autoFit
    /// x-axis zoom (drag on the chart): lower bound, and upper bound or nil = follow the newest data
    @State private var xZoomLower: Double? = nil
    @State private var xZoomUpper: Double? = nil
    @State private var zoomDragStart: CGFloat? = nil
    @State private var zoomDragCurrent: CGFloat? = nil
    @State private var showRemeshEvents: Bool = true
    @State private var showYieldLine: Bool = true
    @State private var useInteriorGBGauge: Bool = true
    @State private var showGBSigmaLabels: Bool = true
    @State private var showStylePopover: Bool = false
    // Line & Point appearance controls
    @State private var showLine: Bool = true
    @State private var showPoints: Bool = false
    @State private var lineWidth: Double = 2.25
    @State private var pointSize: Double = 5.0
    @State private var pointSymbol: ChartPointSymbol = .circle
    @State private var stressLineColor: Color = Color(red: 0.25, green: 0.62, blue: 1.00)
    @State private var stressPointColor: Color = Color(red: 0.35, green: 0.70, blue: 1.00)
    @State private var energyLineColor: Color = Color(red: 0.72, green: 0.48, blue: 1.00)
    @State private var energyPointColor: Color = Color(red: 0.78, green: 0.58, blue: 1.00)
    @State private var eigenLineColor: Color = Color(red: 0.20, green: 0.86, blue: 0.78)
    @State private var eigenPointColor: Color = Color(red: 0.28, green: 0.92, blue: 0.84)
    @State private var gbLineColor: Color = Color(red: 0.28, green: 0.85, blue: 0.52)
    @State private var gbPointColor: Color = Color(red: 0.36, green: 0.90, blue: 0.58)
    /// Force / stress control: plot against the conjugate strain (bursts show as jumps along x)
    @State private var plotAgainstConjugateStrain: Bool = true
    /// Volterra cylinder load: plot against the applied far-field stress σ_xy instead of the load α
    /// Volterra cylinder load: x axis = applied far-field stress σ_xy, load α (γ or τ), or dislocation position x_d / b
    enum VolterraXAxis: Hashable { case appliedStress, load, corePosition }
    @State private var volterraXAxis: VolterraXAxis = .corePosition
    /// Force / stress control: the logged energy is the enthalpy H = E_el - W_ext; plot E_el unless asked
    @State private var plotEnthalpy: Bool = false
    /// PBC stress control: λ_min chart shows the free-cell eigenvalue λ_box (true) or the fixed-cell λ_min(K)
    @State private var eigenShowsCell: Bool = true

    private var isShiftMode: Bool {
        runner.activeConfig?.mode == .shift
    }

    private var isIndentationMode: Bool {
        runner.activeConfig?.geometryConfig.indenter.enabled == true ||
        runner.activeConfig?.geometryConfig.primary_stress_component == "indentation_load"
    }

    private var hasElasticEnergy: Bool {
        runner.curvePoints.contains { $0.elasticEnergy != nil }
    }
    private func energyValue(_ pt: StressStrainPoint) -> Double {
        (hasElasticEnergy && !plotEnthalpy) ? (pt.elasticEnergy ?? pt.postEnergy) : pt.postEnergy
    }
    private var energyLabel: String {
        hasElasticEnergy ? (plotEnthalpy ? "Enthalpy H = E_el − W_ext" : "Elastic Energy E_el") : "Internal Energy (E)"
    }

    /// The run logged a conjugate strain (dead load, force grip or PBC stress control)
    private var hasConjugateStrain: Bool {
        runner.curvePoints.contains { $0.effectiveStrain != 0 }
    }
    /// Volterra cylinder load run: the log carries the core glide x_c − x₀
    private var hasCoreShift: Bool {
        runner.curvePoints.contains { $0.coreShift != nil }
    }
    private var conjugateStrainAxis: Bool {
        // Core glide is always plotted against the load α (γ or τ)
        hasConjugateStrain && plotAgainstConjugateStrain && chartMode != .eigenvalue && chartMode != .gbEnergy
            && !(chartMode == .stress && stressComponent == .coreShift)
    }
    /// Volterra cylinder load: the log carries the ring's applied stress and the core before remeshing
    private var hasAppliedStress: Bool {
        runner.curvePoints.contains { $0.appliedSigmaXY != nil }
    }
    private var hasCoreShiftRelaxed: Bool {
        runner.curvePoints.contains { $0.coreShiftRelaxed != nil }
    }
    private var volterraAxisActive: Bool {
        hasAppliedStress && chartMode != .eigenvalue && chartMode != .gbEnergy
    }
    private var appliedStressAxis: Bool { volterraAxisActive && volterraXAxis == .appliedStress }
    /// Dislocation position x_d / b (core glide in Burgers vectors), as in Fig. 8 of the MTM cylinder paper
    private var corePositionAxis: Bool { volterraAxisActive && volterraXAxis == .corePosition }
    /// x coordinate of a curve point: the load parameter α, the conjugate strain ε under force control, or the
    /// applied far-field stress in a Volterra cylinder load
    private func xValue(_ pt: StressStrainPoint) -> Double {
        if corePositionAxis { return (pt.coreShift ?? 0.0) / max(pt.burgers ?? 1.0, 1e-12) }
        if appliedStressAxis { return pt.appliedSigmaXY ?? pt.alpha }
        return conjugateStrainAxis ? pt.effectiveStrain : pt.alpha
    }
    private var xAxisLabel: String {
        if corePositionAxis { return "Dislocation Position (x_d / b)" }
        if appliedStressAxis { return "Applied Stress (σ_xy)" }
        if conjugateStrainAxis { return "Conjugate Strain (ε)" }
        if isShiftMode { return "Block Shift (u / a₀)" }
        if isIndentationMode { return "Indentation Depth δ (a₀)" }
        return "Shear Strain (α)"
    }
    private var xSymbolShort: String {
        if chartMode == .gbEnergy { return "Δθ" }
        if chartMode == .eigenvalue {
            if eigenStrainAxis { return "ε" }
            if isShiftMode { return "u/a₀" }
            if isIndentationMode { return "δ" }
            return "α"
        }
        if corePositionAxis { return "x_d/b" }
        if appliedStressAxis { return "σ_xy" }
        if conjugateStrainAxis { return "ε" }
        if isShiftMode { return "u/a₀" }
        if isIndentationMode { return "δ" }
        return "α"
    }

    /// PBC stress control: the eigen log carries the free-cell stability λ_box and the strain
    private var hasCellEigen: Bool {
        runner.eigenPoints.contains { $0.lambdaBox != nil }
    }
    private var showCellEigen: Bool { hasCellEigen && eigenShowsCell }
    private var eigenStrainAxis: Bool {
        runner.eigenPoints.contains { $0.strain != nil } && plotAgainstConjugateStrain
    }
    private func eigenX(_ pt: EigenPoint) -> Double {
        eigenStrainAxis ? (pt.strain ?? pt.alpha) : pt.alpha
    }
    private func eigenY(_ pt: EigenPoint) -> Double? {
        showCellEigen ? pt.lambdaBox : pt.lambdaMin
    }
    private func eigenUnstable(_ pt: EigenPoint) -> Bool {
        showCellEigen ? pt.isCellUnstable : pt.isUnstable
    }
    private var eigenXLabel: String {
        if eigenStrainAxis { return "Conjugate Strain (ε)" }
        if isShiftMode { return "Block Shift (u / a₀)" }
        if isIndentationMode { return "Indentation Depth δ (a₀)" }
        return hasConjugateStrain ? "Load (α)" : "Shear Strain (α)"
    }
    private var eigenYLabel: String {
        showCellEigen ? "Cell Eigenvalue (λ_box)" : "Smallest Eigenvalue (λ_min)"
    }

    private var activeXAxisTitle: String {
        switch chartMode {
        case .gbEnergy:   return "Misorientation Angle Δθ (°)"
        case .eigenvalue: return eigenXLabel
        case .stress, .energy: return xAxisLabel
        }
    }

    private var activeYAxisTitle: String {
        switch chartMode {
        case .gbEnergy:
            return useInteriorGBGauge ? "Interior GB Energy γ_GB (E_int / L_int)" : "Total GB Energy γ_GB (E / L_GB)"
        case .eigenvalue:
            return eigenYLabel
        case .energy:
            return energyLabel
        case .stress:
            if stressComponent == .primary {
                if isIndentationMode { return "Indentation Load P" }
                if primaryIsTorque { return "Rim Torque T" }
                if isShiftMode { return "Shear Stress σ" }
                return "Primary Stress σ"
            }
            return stressComponent.rawValue
        }
    }

    private var chartModeTabItems: [BlueHourTabItem<ChartMetricMode>] {
        let paramSym = isShiftMode ? "u" : (isIndentationMode ? "δ" : "α")
        let stressTitle = isIndentationMode ? "Load P(δ)" : "Stress σ(\(paramSym))"
        let energyTitle = "Energy E(\(paramSym))"
        let eigenTitle  = "λ_min(\(paramSym))"
        let hasCurve    = !runner.curvePoints.isEmpty
        let hasEigen    = !runner.eigenPoints.isEmpty
        let hasGB       = !displayGBPoints.isEmpty

        return [
            BlueHourTabItem(
                .stress,
                title: stressTitle,
                hasData: hasCurve || runner.status == .idle,
                liveDotColor: hasCurve ? stressLineColor : nil,
                activeFill: BlueHour.stressTabFill,
                activeBorder: BlueHour.stressTabBorder,
                helpText: "Plot mechanical response: \(stressTitle)"
            ),
            BlueHourTabItem(
                .energy,
                title: energyTitle,
                hasData: hasCurve || runner.status == .idle,
                liveDotColor: hasCurve ? energyLineColor : nil,
                activeFill: BlueHour.energyTabFill,
                activeBorder: BlueHour.energyTabBorder,
                helpText: "Plot stored elastic / internal energy: \(energyTitle)"
            ),
            BlueHourTabItem(
                .eigenvalue,
                title: eigenTitle,
                hasData: hasEigen,
                liveDotColor: hasEigen ? eigenLineColor : nil,
                activeFill: BlueHour.eigenTabFill,
                activeBorder: BlueHour.eigenTabBorder,
                helpText: hasEigen
                    ? "Plot smallest Hessian eigenvalue \(eigenTitle)"
                    : "Smallest Hessian eigenvalue \(eigenTitle) — no eigenvalues recorded in this run (click for instructions)"
            ),
            BlueHourTabItem(
                .gbEnergy,
                title: "γ_GB(Δθ)",
                badge: hasGB ? "\(displayGBPoints.count)" : nil,
                hasData: hasGB,
                liveDotColor: hasGB ? gbLineColor : nil,
                activeFill: BlueHour.gbTabFill,
                activeBorder: BlueHour.gbTabBorder,
                helpText: hasGB
                    ? "Grain Boundary Energy γ_GB vs. Misorientation Δθ (\(displayGBPoints.count) recorded angles)"
                    : "Grain Boundary Energy γ_GB(Δθ) — run a Bicrystal relaxation or Δθ sweep to populate"
            )
        ]
    }
    
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
    
    // MARK: - x-axis zoom

    private var isXZoomed: Bool { xZoomLower != nil && chartMode != .gbEnergy }
    /// Largest x of all data in the current mode (not only the visible part)
    private var dataXMax: Double? {
        chartMode == .eigenvalue
            ? runner.eigenPoints.filter { eigenY($0) != nil }.map { eigenX($0) }.max()
            : runner.curvePoints.map { xValue($0) }.max()
    }
    private func inXZoom(_ x: Double) -> Bool {
        guard isXZoomed, let lo = xZoomLower else { return true }
        return x >= lo && (xZoomUpper.map { x <= $0 } ?? true)
    }
    private func resetXZoom() {
        xZoomLower = nil
        xZoomUpper = nil
        zoomDragStart = nil
        zoomDragCurrent = nil
    }

    // MARK: - Dynamic Clever Domain Calculations
    
    private var activeAlphas: [Double] {
        if chartMode == .gbEnergy {
            return displayGBPoints.map(\.misorientationDeg)
        }
        if chartMode == .eigenvalue && !runner.eigenPoints.isEmpty {
            return displayEigenPoints.map { eigenX($0) }
        }
        return runner.curvePoints.map { xValue($0) }
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
        
        if isXZoomed, let lo = xZoomLower {
            let hi = xZoomUpper ?? max(dataXMax ?? lo, lo)
            let span = max(hi - lo, 1e-9)
            return lo...(hi + (xZoomUpper == nil ? span * 0.05 : 0.0))
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

        // Conjugate-strain axis: the schedule (α_start → α_end) is a load range, not a strain range
        if conjugateStrainAxis {
            let effectiveSpan = max(0.0008, span)
            let lo = scaleMode == .autoFit ? minA - effectiveSpan * 0.025 : min(0.0, minA)
            return lo...(maxA + effectiveSpan * 0.065)
        }
        
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
            let zoomed = runner.curvePoints.filter { inXZoom(xValue($0)) }
            let vals = (zoomed.isEmpty ? runner.curvePoints : zoomed).map { stressValue(for: $0) }
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
            let zoomed = runner.curvePoints.filter { inXZoom(xValue($0)) }
            let vals = (zoomed.isEmpty ? runner.curvePoints : zoomed).map { energyValue($0) }
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
            let zoomedEig = runner.eigenPoints.filter { inXZoom(eigenX($0)) }
            let vals = (zoomedEig.isEmpty ? runner.eigenPoints : zoomedEig).compactMap { eigenY($0) }
            guard !vals.isEmpty else { return -0.05...0.5 }
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
    
    // MARK: - Line & Point Style Helpers
    
    private var currentLineColor: Color {
        switch chartMode {
        case .stress: return stressLineColor
        case .energy: return energyLineColor
        case .eigenvalue: return eigenLineColor
        case .gbEnergy: return gbLineColor
        }
    }
    
    private var currentPointColor: Color {
        switch chartMode {
        case .stress: return stressPointColor
        case .energy: return energyPointColor
        case .eigenvalue: return eigenPointColor
        case .gbEnergy: return gbPointColor
        }
    }
    
    private var activeLineColorBinding: Binding<Color> {
        Binding(
            get: { currentLineColor },
            set: { newColor in
                switch chartMode {
                case .stress: stressLineColor = newColor
                case .energy: energyLineColor = newColor
                case .eigenvalue: eigenLineColor = newColor
                case .gbEnergy: gbLineColor = newColor
                }
                if !showLine { showLine = true }
            }
        )
    }
    
    private var activePointColorBinding: Binding<Color> {
        Binding(
            get: { currentPointColor },
            set: { newColor in
                switch chartMode {
                case .stress: stressPointColor = newColor
                case .energy: energyPointColor = newColor
                case .eigenvalue: eigenPointColor = newColor
                case .gbEnergy: gbPointColor = newColor
                }
                if !showPoints { showPoints = true }
            }
        )
    }
    
    private var lineToggleBinding: Binding<Bool> {
        Binding(
            get: { showLine },
            set: { newValue in
                showLine = newValue
                if !newValue && !showPoints {
                    showPoints = true
                }
            }
        )
    }
    
    private var pointsToggleBinding: Binding<Bool> {
        Binding(
            get: { showPoints },
            set: { newValue in
                showPoints = newValue
                if !newValue && !showLine {
                    showLine = true
                }
            }
        )
    }
    
    private var pointSymbolArea: CGFloat {
        CGFloat(pointSize * pointSize)
    }
    
    private func resetPlotStyle() {
        showLine = true
        showPoints = false
        lineWidth = 2.25
        pointSize = 5.0
        pointSymbol = .circle
        stressLineColor = Color(red: 0.25, green: 0.62, blue: 1.00)
        stressPointColor = Color(red: 0.35, green: 0.70, blue: 1.00)
        energyLineColor = Color(red: 0.72, green: 0.48, blue: 1.00)
        energyPointColor = Color(red: 0.78, green: 0.58, blue: 1.00)
        eigenLineColor = Color(red: 0.20, green: 0.86, blue: 0.78)
        eigenPointColor = Color(red: 0.28, green: 0.92, blue: 0.84)
        gbLineColor = Color(red: 0.28, green: 0.85, blue: 0.52)
        gbPointColor = Color(red: 0.36, green: 0.90, blue: 0.58)
    }
    
    @ViewBuilder
    private var plotStylePopoverContent: some View {
        VStack(alignment: .leading, spacing: 12) {
            HStack {
                Label("Curve Appearance", systemImage: "paintbrush.pointed.fill")
                    .font(.system(size: 12, weight: .bold))
                    .foregroundColor(BlueHour.ice)
                Spacer()
                Button(action: resetPlotStyle) {
                    Label("Reset", systemImage: "arrow.counterclockwise")
                        .font(.system(size: 10.5, weight: .semibold))
                }
                .buttonStyle(.bordered)
                .controlSize(.mini)
                .help("Reset line & point colors, widths, and markers to defaults")
            }

            Divider()

            // Line controls
            VStack(alignment: .leading, spacing: 6) {
                HStack {
                    Toggle("Connecting Line", isOn: lineToggleBinding)
                        .toggleStyle(.checkbox)
                        .font(.system(size: 11.5, weight: .semibold))
                    Spacer()
                    ColorPicker("", selection: activeLineColorBinding, supportsOpacity: true)
                        .labelsHidden()
                        .controlSize(.small)
                        .frame(width: 28)
                }
                HStack(spacing: 8) {
                    Text("Width")
                        .font(.system(size: 10.5))
                        .foregroundColor(BlueHour.mist)
                        .frame(width: 40, alignment: .leading)
                    Slider(
                        value: Binding(
                            get: { lineWidth },
                            set: { v in
                                lineWidth = v
                                if !showLine { showLine = true }
                            }
                        ),
                        in: 0.5...6.0,
                        step: 0.25
                    )
                    .controlSize(.small)
                    Text(String(format: "%.2g pt", lineWidth))
                        .font(.system(size: 10.5, design: .monospaced))
                        .foregroundColor(BlueHour.ice)
                        .frame(width: 42, alignment: .trailing)
                }
                .opacity(showLine ? 1.0 : 0.45)
            }

            Divider()

            // Point marker controls
            VStack(alignment: .leading, spacing: 6) {
                HStack {
                    Toggle("Data Point Markers", isOn: pointsToggleBinding)
                        .toggleStyle(.checkbox)
                        .font(.system(size: 11.5, weight: .semibold))
                    Spacer()
                    ColorPicker("", selection: activePointColorBinding, supportsOpacity: true)
                        .labelsHidden()
                        .controlSize(.small)
                        .frame(width: 28)
                }
                HStack(spacing: 8) {
                    Text("Size")
                        .font(.system(size: 10.5))
                        .foregroundColor(BlueHour.mist)
                        .frame(width: 40, alignment: .leading)
                    Slider(
                        value: Binding(
                            get: { pointSize },
                            set: { v in
                                pointSize = v
                                if !showPoints { showPoints = true }
                            }
                        ),
                        in: 2.0...14.0,
                        step: 0.5
                    )
                    .controlSize(.small)
                    Text(String(format: "%.2g pt", pointSize))
                        .font(.system(size: 10.5, design: .monospaced))
                        .foregroundColor(BlueHour.ice)
                        .frame(width: 42, alignment: .trailing)
                }
                .opacity(showPoints ? 1.0 : 0.45)

                HStack(spacing: 8) {
                    Text("Shape")
                        .font(.system(size: 10.5))
                        .foregroundColor(BlueHour.mist)
                        .frame(width: 40, alignment: .leading)
                    Picker("", selection: $pointSymbol) {
                        ForEach(ChartPointSymbol.allCases) { sym in
                            Text(sym.rawValue).tag(sym)
                        }
                    }
                    .labelsHidden()
                    .pickerStyle(.menu)
                    .controlSize(.small)
                }
                .opacity(showPoints ? 1.0 : 0.45)
            }
        }
        .padding(14)
        .frame(width: 278)
        .background(BlueHour.cardSurface)
    }

    @ViewBuilder
    private var chartTitleBanner: some View {
        let ptCount: Int = {
            switch chartMode {
            case .eigenvalue: return runner.eigenPoints.count
            case .gbEnergy:   return displayGBPoints.count
            default:          return runner.curvePoints.count
            }
        }()

        HStack(spacing: 10) {
            // Active curve visual swatch
            ZStack {
                if showLine {
                    Capsule()
                        .fill(currentLineColor)
                        .frame(width: 20, height: max(2.0, min(4.5, lineWidth)))
                }
                if showPoints {
                    Circle()
                        .fill(currentPointColor)
                        .frame(width: 7, height: 7)
                }
            }
            .frame(width: 22, height: 12)

            // Self-describing plot title: Y vs X
            HStack(spacing: 5) {
                Text(activeYAxisTitle)
                    .font(.system(size: 12, weight: .bold))
                    .foregroundColor(BlueHour.ice)
                Text("vs.")
                    .font(.system(size: 11, weight: .regular))
                    .foregroundColor(BlueHour.mist.opacity(0.8))
                Text(activeXAxisTitle)
                    .font(.system(size: 12, weight: .semibold))
                    .foregroundColor(BlueHour.mist)
            }
            .lineLimit(1)

            Spacer(minLength: 8)

            // Live GB Energy shortcut pill (when viewing another tab during a Bicrystal run)
            if chartMode != .gbEnergy, let curGB = runner.currentGBEnergy {
                let gVal = useInteriorGBGauge ? curGB.gbEnergyInterior : curGB.gbEnergyTotal
                Button(action: { chartMode = .gbEnergy }) {
                    HStack(spacing: 4) {
                        Image(systemName: "square.split.2x1.fill")
                            .font(.system(size: 9))
                        Text("γ_GB = \(String(format: "%.5f", gVal)) ↗")
                            .font(.system(size: 10, weight: .bold, design: .monospaced))
                    }
                    .foregroundColor(Color(red: 0.43, green: 0.91, blue: 0.61))
                    .padding(.horizontal, 7)
                    .padding(.vertical, 2.5)
                    .background(BlueHour.gbTabFill.opacity(0.32))
                    .overlay(RoundedRectangle(cornerRadius: 5).strokeBorder(BlueHour.gbTabBorder.opacity(0.55), lineWidth: 1))
                    .cornerRadius(5)
                }
                .buttonStyle(.plain)
                .help("Click to switch to Grain Boundary Energy γ_GB(Δθ) plot")
            }

            // Peak Stress or Stability summary pill
            let peakVal = currentPeakStress
            if chartMode == .stress && stressComponent != .coreShift && peakVal > 0 {
                HStack(spacing: 4) {
                    Image(systemName: "crown.fill")
                        .font(.system(size: 9))
                    Text("|\(stressComponent.shortLabel)|_max = \(String(format: "%.4f", peakVal))")
                        .font(.system(size: 10, weight: .bold, design: .monospaced))
                }
                .foregroundColor(.orange)
                .padding(.horizontal, 7)
                .padding(.vertical, 2.5)
                .background(Color.orange.opacity(0.12))
                .overlay(RoundedRectangle(cornerRadius: 5).strokeBorder(Color.orange.opacity(0.42), lineWidth: 1))
                .cornerRadius(5)
                .fixedSize()
            } else if chartMode == .eigenvalue,
                      let lastEig = runner.eigenPoints.last(where: { eigenY($0) != nil }),
                      let lastY = eigenY(lastEig) {
                let unstable = eigenUnstable(lastEig)
                HStack(spacing: 4) {
                    Circle()
                        .fill(unstable ? Color.red : Color.green)
                        .frame(width: 6, height: 6)
                    Text("\(showCellEigen ? "λ_box" : "λ_min") = \(String(format: "%.3e", lastY))")
                        .font(.system(size: 10, weight: .bold, design: .monospaced))
                }
                .foregroundColor(unstable ? .red : Color(red: 0.43, green: 0.91, blue: 0.61))
                .padding(.horizontal, 7)
                .padding(.vertical, 2.5)
                .background((unstable ? Color.red : Color.green).opacity(0.12))
                .overlay(RoundedRectangle(cornerRadius: 5).strokeBorder((unstable ? Color.red : Color.green).opacity(0.42), lineWidth: 1))
                .cornerRadius(5)
                .fixedSize()
            }

            // Active X range badge
            if let firstX = activeAlphas.first, let lastX = activeAlphas.last {
                Text("\(xSymbolShort) ∈ [\(formatXValue(firstX)) → \(formatXValue(lastX))]")
                    .font(.system(size: 10, weight: .medium, design: .monospaced))
                    .foregroundColor(BlueHour.mist)
                    .padding(.horizontal, 7)
                    .padding(.vertical, 2.5)
                    .background(BlueHour.controlSurface.opacity(0.85))
                    .overlay(RoundedRectangle(cornerRadius: 5).strokeBorder(BlueHour.border, lineWidth: 1))
                    .cornerRadius(5)
                    .fixedSize()
            }

            // Point count pill
            Text("\(ptCount) pts")
                .font(.system(size: 10, weight: .semibold, design: .monospaced))
                .foregroundColor(BlueHour.ice)
                .padding(.horizontal, 7)
                .padding(.vertical, 2.5)
                .background(BlueHour.controlSurface.opacity(0.85))
                .overlay(RoundedRectangle(cornerRadius: 5).strokeBorder(BlueHour.border, lineWidth: 1))
                .cornerRadius(5)
                .fixedSize()
        }
        .padding(.horizontal, 16)
        .padding(.top, 10)
        .padding(.bottom, 2)
        .background(BlueHour.canvasBackground)
    }

    // MARK: - View Body
    
    var body: some View {
        VStack(spacing: 0) {
            // Single clean, high-contrast Chart Control Toolbar (scrollable if narrowed in Split View so nothing ever clips)
            ScrollView(.horizontal, showsIndicators: false) {
                HStack(spacing: 8) {
                    // 1. Color-coded High-Contrast Metric Tabs
                    BlueHourSegmentedPicker(
                        selection: $chartMode,
                        items: chartModeTabItems,
                        scale: .compact
                    )

                    // 1b. Contextual Metric Sub-Selector
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
                            Section(header: Text("Force / Stress Control")) {
                                Text(StressComponentType.effectiveStrain.rawValue).tag(StressComponentType.effectiveStrain)
                            }
                            if hasCoreShift {
                                Section(header: Text("Dislocation Core (Volterra Load)")) {
                                    Text(StressComponentType.appliedSigmaXY.rawValue).tag(StressComponentType.appliedSigmaXY)
                                    Text(StressComponentType.coreShift.rawValue).tag(StressComponentType.coreShift)
                                }
                            }
                        }
                        .onChange(of: hasCoreShift, initial: true) { _, present in
                            if present && stressComponent == .primary { stressComponent = .cauchyXY }
                            if !present && (stressComponent == .coreShift || stressComponent == .appliedSigmaXY) {
                                stressComponent = .primary
                            }
                        }
                        .pickerStyle(.menu)
                        .labelsHidden()
                        .controlSize(.small)
                        .frame(width: 162)
                        .help("Select which stress component (Cauchy, Aligned Piola / Projected, Von Mises) or nanoindentation metric to plot")
                    } else if chartMode == .energy && hasElasticEnergy {
                        BlueHourSegmentedPicker(
                            selection: $plotEnthalpy,
                            items: [
                                BlueHourTabItem(false, title: "E_el", activeFill: BlueHour.energyTabFill, activeBorder: BlueHour.energyTabBorder, helpText: "Stored elastic energy E_el"),
                                BlueHourTabItem(true, title: "H", activeFill: BlueHour.energyTabFill, activeBorder: BlueHour.energyTabBorder, helpText: "Enthalpy H = E_el − W_ext")
                            ],
                            scale: .compact
                        )
                    } else if chartMode == .eigenvalue && hasCellEigen {
                        BlueHourSegmentedPicker(
                            selection: $eigenShowsCell,
                            items: [
                                BlueHourTabItem(true, title: "Cell", activeFill: BlueHour.eigenTabFill, activeBorder: BlueHour.eigenTabBorder, helpText: "Relaxed free-cell eigenvalue λ_box"),
                                BlueHourTabItem(false, title: "Atoms", activeFill: BlueHour.eigenTabFill, activeBorder: BlueHour.eigenTabBorder, helpText: "Fixed-cell atomic stiffness λ_min(K)")
                            ],
                            scale: .compact
                        )
                    }

                    // 1c. Contextual X-Axis Variable Switcher (Volterra or Force/Stress Control)
                    if hasAppliedStress && chartMode != .gbEnergy && chartMode != .eigenvalue {
                        BlueHourSegmentedPicker(
                            selection: $volterraXAxis,
                            items: [
                                BlueHourTabItem(.corePosition, title: "x: x_d/b", helpText: "Dislocation position x_d / b"),
                                BlueHourTabItem(.appliedStress, title: "x: σ_xy", helpText: "Far-field applied shear stress σ_xy"),
                                BlueHourTabItem(.load, title: "x: α", helpText: "Applied load parameter α")
                            ],
                            scale: .compact
                        )
                    } else if (hasConjugateStrain && chartMode != .gbEnergy && chartMode != .eigenvalue) ||
                              (chartMode == .eigenvalue && runner.eigenPoints.contains { $0.strain != nil }) {
                        BlueHourSegmentedPicker(
                            selection: $plotAgainstConjugateStrain,
                            items: [
                                BlueHourTabItem(true, title: "x: ε", helpText: "Plot against conjugate strain ε (shows bursts as jumps along x)"),
                                BlueHourTabItem(false, title: "x: α", helpText: "Plot against applied load parameter α")
                            ],
                            scale: .compact
                        )
                    }

                    Divider()
                        .frame(height: 16)

                    // 2. X-Axis Domain Mode
                    BlueHourSegmentedPicker(
                        selection: $scaleMode,
                        items: ChartScaleMode.allCases.map { m in
                            BlueHourTabItem(m, title: m.rawValue, helpText: m.tooltip)
                        },
                        scale: .compact
                    )

                    if isXZoomed {
                        Button(action: resetXZoom) {
                            Label("Reset Zoom", systemImage: "arrow.uturn.backward")
                                .font(.system(size: 10.5, weight: .semibold))
                        }
                        .buttonStyle(.bordered)
                        .controlSize(.small)
                        .tint(.cyan)
                        .help("Reset custom X-axis zoom (or double-click on the plot)")
                    }

                    Divider()
                        .frame(height: 16)

                    // 3. Explicit Overlay Toggle Chips (never clipped!)
                    if chartMode == .gbEnergy {
                        BlueHourToggleChip(
                            "Interior Gauge",
                            isOn: $useInteriorGBGauge,
                            activeFill: BlueHour.gbTabFill,
                            activeBorder: BlueHour.gbTabBorder,
                            scale: .compact,
                            helpText: "Use interior gauge E_int / L_int (excluding 4a₀ outer boundary collar)"
                        )
                        BlueHourToggleChip(
                            "Σ Labels",
                            isOn: $showGBSigmaLabels,
                            activeFill: BlueHour.gbTabFill,
                            activeBorder: BlueHour.gbTabBorder,
                            scale: .compact,
                            helpText: "Display Σ CSL cusp labels above points"
                        )
                    } else if chartMode != .eigenvalue {
                        BlueHourToggleChip(
                            "Remesh",
                            badge: "\(remeshEvents.count)",
                            isOn: $showRemeshEvents,
                            activeFill: Color(red: 0.62, green: 0.16, blue: 0.22),
                            activeBorder: Color(red: 1.0, green: 0.48, blue: 0.52),
                            scale: .compact,
                            helpText: "Highlight Delaunay remeshing & plastic slip surgery events on the curve"
                        )
                    }

                    if chartMode == .stress && stressComponent != .coreShift {
                        BlueHourToggleChip(
                            "Peak",
                            isOn: $showYieldLine,
                            activeFill: Color(red: 0.58, green: 0.34, blue: 0.08),
                            activeBorder: Color.orange,
                            scale: .compact,
                            helpText: "Show dashed horizontal peak stress guide line"
                        )
                    }

                    Divider()
                        .frame(height: 16)

                    // 4. Curve Style Popover Button + Export Menu
                    Button {
                        showStylePopover.toggle()
                    } label: {
                        HStack(spacing: 5) {
                            Circle()
                                .fill(currentLineColor)
                                .frame(width: 7, height: 7)
                            Text("Style")
                                .font(.system(size: 10.5, weight: .semibold))
                            Image(systemName: "chevron.down")
                                .font(.system(size: 8.5, weight: .bold))
                                .foregroundColor(BlueHour.mist)
                        }
                        .foregroundColor(BlueHour.ice)
                        .padding(.horizontal, 8)
                        .padding(.vertical, 4)
                        .background(
                            RoundedRectangle(cornerRadius: 7, style: .continuous)
                                .fill(showStylePopover ? BlueHour.elevatedSurface : BlueHour.tabTrack)
                        )
                        .overlay(
                            RoundedRectangle(cornerRadius: 7, style: .continuous)
                                .strokeBorder(BlueHour.steel.opacity(0.5), lineWidth: 1)
                        )
                    }
                    .buttonStyle(.plain)
                    .popover(isPresented: $showStylePopover, arrowEdge: .bottom) {
                        plotStylePopoverContent
                    }
                    .help("Customize curve line width, point markers, colors, and symbols")

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
                            Label("Export", systemImage: "square.and.arrow.up")
                        }
                        .menuStyle(.borderedButton)
                        .controlSize(.small)
                        .fixedSize()
                    }
                }
                .padding(.horizontal, 12)
                .padding(.vertical, 6)
            }
            .background(BlueHour.controlSurface.opacity(0.9))
            
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
                         ? "Enable 'Compute Lowest Stiffness Eigenvalues' in the 4 Solver sidebar tab and launch the simulation to plot λ_min live."
                         : "Eigenvalue computation was disabled for this run. Enable 'Compute Lowest Stiffness Eigenvalues' in 4 Solver for the next run.")
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
                    Text("Click 'Launch' in the bottom-right corner to stream live stress, energy, and stability curves in real time.")
                        .font(.subheadline)
                        .foregroundColor(.secondary.opacity(0.8))
                        .multilineTextAlignment(.center)
                        .frame(maxWidth: 420)
                }
                .frame(maxWidth: .infinity, maxHeight: .infinity)
                .background(BlueHour.canvasBackground)
            } else {
                chartTitleBanner

                Chart {
                    if chartMode == .gbEnergy {
                        ForEach(displayGBPoints) { pt in
                            let yVal = gbEnergyValue(for: pt)
                            if showLine {
                                LineMark(
                                    x: .value("Misorientation Δθ (°)", pt.misorientationDeg),
                                    y: .value("Grain Boundary Energy γ_GB", yVal)
                                )
                                .foregroundStyle(currentLineColor)
                                .lineStyle(StrokeStyle(lineWidth: lineWidth))
                                .interpolationMethod(.linear)
                            }
                            
                            if showPoints {
                                PointMark(
                                    x: .value("Misorientation Δθ (°)", pt.misorientationDeg),
                                    y: .value("Grain Boundary Energy γ_GB", yVal)
                                )
                                .symbol(pointSymbol.chartSymbol)
                                .symbolSize(pointSymbolArea)
                                .foregroundStyle(currentPointColor)
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
                            if let yVal = eigenY(pt) {
                                if showLine {
                                    LineMark(
                                        x: .value(eigenXLabel, eigenX(pt)),
                                        y: .value(eigenYLabel, yVal)
                                    )
                                    .foregroundStyle(currentLineColor)
                                    .lineStyle(StrokeStyle(lineWidth: lineWidth))
                                    .interpolationMethod(.linear)
                                }
                                
                                if showPoints {
                                    PointMark(
                                        x: .value(eigenXLabel, eigenX(pt)),
                                        y: .value(eigenYLabel, yVal)
                                    )
                                    .symbol(pointSymbol.chartSymbol)
                                    .symbolSize(eigenUnstable(pt) ? max(pointSymbolArea * 1.4, 55) : pointSymbolArea)
                                    .foregroundStyle(eigenUnstable(pt) ? Color.red : currentPointColor)
                                }
                            }
                        }
                    } else {
                        // 1. Primary Stress / Energy Line Mark & PointMark
                        ForEach(displayCurvePoints) { pt in
                            let yVal = chartMode == .stress ? stressValue(for: pt) : energyValue(pt)
                            let yLabel = chartMode == .stress ? ((stressComponent == .primary && primaryIsTorque) ? "Rim Torque T" : stressComponent.rawValue) : energyLabel
                            if showLine {
                                LineMark(
                                    x: .value(xAxisLabel, xValue(pt)),
                                    y: .value(yLabel, yVal),
                                    series: .value("Curve", "main")
                                )
                                .foregroundStyle(currentLineColor)
                                .lineStyle(StrokeStyle(lineWidth: lineWidth))
                                .interpolationMethod(.linear)
                            }
                            
                            if showPoints {
                                PointMark(
                                    x: .value(xAxisLabel, xValue(pt)),
                                    y: .value(yLabel, yVal)
                                )
                                .symbol(pointSymbol.chartSymbol)
                                .symbolSize(pointSymbolArea)
                                .foregroundStyle(currentPointColor)
                            }
                        }
                        
                        // 1b. Volterra cylinder load: core glide before remeshing (dashed)
                        if chartMode == .stress && stressComponent == .coreShift && hasCoreShiftRelaxed {
                            ForEach(displayCurvePoints) { pt in
                                if let yRelaxed = pt.coreShiftRelaxed {
                                    LineMark(
                                        x: .value(xAxisLabel, xValue(pt)),
                                        y: .value(stressComponent.rawValue, yRelaxed),
                                        series: .value("Curve", "before remesh")
                                    )
                                    .foregroundStyle(Color.orange)
                                    .lineStyle(StrokeStyle(lineWidth: 1.2, dash: [4, 3]))
                                    .interpolationMethod(.linear)
                                }
                            }
                        }

                        // 2. Plasticity / Remesh Events Overlay
                        if showRemeshEvents {
                            ForEach(displayRemeshEvents) { pt in
                                let yVal = chartMode == .stress ? stressValue(for: pt) : energyValue(pt)
                                let yLabel = chartMode == .stress ? ((stressComponent == .primary && primaryIsTorque) ? "Rim Torque T" : stressComponent.rawValue) : energyLabel
                                PointMark(
                                    x: .value(xAxisLabel, xValue(pt)),
                                    y: .value(yLabel, yVal)
                                )
                                .symbol(Circle())
                                .symbolSize(48)
                                .foregroundStyle(Color.red)
                                
                                PointMark(
                                    x: .value(xAxisLabel, xValue(pt)),
                                    y: .value(yLabel, yVal)
                                )
                                .symbol(Circle().strokeBorder(lineWidth: 1.5))
                                .symbolSize(80)
                                .foregroundStyle(Color.red.opacity(0.85))
                            }
                        }
                        
                        // 3. Peak / Yield Stress Guide Line
                        let peakVal = currentPeakStress
                        // (not for the core glide, where a "peak" line has no meaning)
                        if showYieldLine && chartMode == .stress && stressComponent != .coreShift && peakVal > 0 {
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
                    }
                }
                .chartXScale(domain: xDomain)
                .chartYScale(domain: yDomain)
                // Volterra load: legend for the before-remesh curve (on the plot, so the toolbar stays one row)
                .overlay(alignment: .topLeading) {
                    if chartMode == .stress && stressComponent == .coreShift && hasCoreShiftRelaxed {
                        Text("solid: final · dashed: before remesh")
                            .font(.system(size: 10, design: .monospaced))
                            .foregroundColor(.secondary)
                            .lineLimit(1)
                            .fixedSize()
                            .padding(.horizontal, 6)
                            .padding(.vertical, 3)
                            .background(BlueHour.controlSurface.opacity(0.85))
                            .cornerRadius(4)
                            .padding(6)
                            .allowsHitTesting(false)
                            .help("Core glide after each step (solid) and after the fixed-mesh relaxation, before that step's remeshing (dashed). A jump only on the solid curve came from the reconnection; a jump on the dashed curve is an instability of the lattice energy.")
                    }
                }
                .chartXAxis {
                    AxisMarks(values: .automatic(desiredCount: 7)) { value in
                        AxisGridLine(stroke: StrokeStyle(lineWidth: 0.5, dash: [2, 2]))
                            .foregroundStyle(Color.secondary.opacity(0.25))
                        AxisTick()
                        if let d = value.as(Double.self) {
                            AxisValueLabel {
                                Text(formatXValue(d))
                                    .font(.system(size: 10, design: .monospaced))
                                    .foregroundColor(BlueHour.mist)
                            }
                        }
                    }
                }
                .chartYAxis {
                    AxisMarks(position: .leading, values: .automatic(desiredCount: 6)) { value in
                        AxisGridLine(stroke: StrokeStyle(lineWidth: 0.5, dash: [2, 2]))
                            .foregroundStyle(Color.secondary.opacity(0.25))
                        AxisTick()
                        if let d = value.as(Double.self) {
                            AxisValueLabel {
                                Text(formatYValue(d))
                                    .font(.system(size: 10, design: .monospaced))
                                    .foregroundColor(BlueHour.mist)
                            }
                        }
                    }
                }
                .chartXAxisLabel(position: .bottom, alignment: .center) {
                    Text(activeXAxisTitle)
                        .font(.system(size: 11, weight: .semibold))
                        .foregroundColor(BlueHour.mist)
                        .padding(.top, 2)
                }
                .chartYAxisLabel(position: .leading, alignment: .center) {
                    Text(activeYAxisTitle)
                        .font(.system(size: 11, weight: .semibold))
                        .foregroundColor(BlueHour.mist)
                        .padding(.bottom, 2)
                }
                .chartOverlay { proxy in
                    GeometryReader { geo in
                        let plotFrame = proxy.plotFrame.map { geo[$0] } ?? CGRect(origin: .zero, size: geo.size)
                        Rectangle()
                            .fill(Color.clear)
                            .contentShape(Rectangle())
                            .onContinuousHover { phase in
                                switch phase {
                                case .active(let location):
                                    let relX = location.x - plotFrame.origin.x
                                    if let rawX: Double = proxy.value(atX: relX) {
                                        let snapped: Double? = {
                                            switch chartMode {
                                            case .gbEnergy:
                                                return findClosestGBPoint(to: rawX)?.misorientationDeg
                                            case .eigenvalue:
                                                return findClosestEigenPoint(to: rawX).map { eigenX($0) }
                                            default:
                                                return findClosestPoint(to: rawX).map { xValue($0) }
                                            }
                                        }()
                                        if selectedAlpha != snapped {
                                            selectedAlpha = snapped
                                        }
                                    }
                                case .ended:
                                    if selectedAlpha != nil {
                                        selectedAlpha = nil
                                    }
                                }
                            }
                            // Drag across the chart to zoom the x axis; double-click to show all points
                            .gesture(
                                DragGesture(minimumDistance: 6)
                                    .onChanged { value in
                                        guard chartMode != .gbEnergy else { return }
                                        zoomDragStart = value.startLocation.x
                                        zoomDragCurrent = value.location.x
                                    }
                                    .onEnded { value in
                                        defer { zoomDragStart = nil; zoomDragCurrent = nil }
                                        guard chartMode != .gbEnergy,
                                              let a: Double = proxy.value(atX: value.startLocation.x - plotFrame.origin.x),
                                              let b: Double = proxy.value(atX: value.location.x - plotFrame.origin.x),
                                              abs(value.location.x - value.startLocation.x) > 6 else { return }
                                        let lo = min(a, b), hi = max(a, b)
                                        // A selection that reaches the newest point keeps following the live data
                                        let follow = (dataXMax.map { hi >= $0 } ?? true)
                                        xZoomLower = lo
                                        xZoomUpper = follow ? nil : hi
                                    }
                            )
                            .onTapGesture(count: 2) { resetXZoom() }
                        
                        // Lightweight hover cursor rendered in the overlay (avoids rebuilding Chart marks on hover)
                        if let selX = selectedAlpha,
                           let hovered = hoveredCoordinates(for: selX),
                           let px = proxy.position(forX: hovered.x),
                           let py = proxy.position(forY: hovered.y) {
                            let screenX = plotFrame.origin.x + px
                            let screenY = plotFrame.origin.y + py
                            if screenX >= plotFrame.minX && screenX <= plotFrame.maxX {
                                Path { path in
                                    path.move(to: CGPoint(x: screenX, y: plotFrame.minY))
                                    path.addLine(to: CGPoint(x: screenX, y: plotFrame.maxY))
                                }
                                .stroke(Color.cyan.opacity(0.55), style: StrokeStyle(lineWidth: 1, dash: [3, 2]))
                                .allowsHitTesting(false)
                                
                                if screenY >= plotFrame.minY && screenY <= plotFrame.maxY {
                                    Circle()
                                        .fill(Color.cyan)
                                        .overlay(Circle().strokeBorder(Color.white, lineWidth: 1.5))
                                        .frame(width: 9, height: 9)
                                        .position(x: screenX, y: screenY)
                                        .allowsHitTesting(false)
                                }
                            }
                        }
                        
                        if let s0 = zoomDragStart, let s1 = zoomDragCurrent {
                            Rectangle()
                                .fill(BlueHour.activeAccent.opacity(0.16))
                                .overlay(Rectangle().stroke(BlueHour.activeAccent.opacity(0.7), lineWidth: 1))
                                .frame(width: abs(s1 - s0), height: geo.size.height)
                                .position(x: (s0 + s1) / 2, y: geo.size.height / 2)
                                .allowsHitTesting(false)
                        }
                    }
                }
                .help("Drag across the chart to zoom the x axis; double-click to show all points.")
                .onChange(of: chartMode) { _, _ in resetXZoom() }
                .onChange(of: plotAgainstConjugateStrain) { _, _ in resetXZoom() }
                .onChange(of: volterraXAxis) { _, _ in resetXZoom() }
                .onChange(of: scaleMode) { _, _ in resetXZoom() }
                .onChange(of: runner.outputDirectory) { _, _ in resetXZoom() }
                .padding(.horizontal, 16)
                .padding(.top, 8)
                .padding(.bottom, 8)
                .background(BlueHour.canvasBackground)
                
                Divider()

                // Fixed-height (30 pt) non-jittering bottom readout strip (Cursor when hovering, Latest otherwise)
                bottomReadoutBar
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

    @ViewBuilder
    private var bottomReadoutBar: some View {
        let isHovering = (selectedAlpha != nil)

        HStack(spacing: 12) {
            // Mode badge: CURSOR vs LATEST
            Text(isHovering ? "CURSOR" : "LATEST")
                .font(.system(size: 9, weight: .bold, design: .monospaced))
                .foregroundColor(isHovering ? .black : BlueHour.mist)
                .padding(.horizontal, 5.5)
                .padding(.vertical, 2)
                .background(
                    Capsule()
                        .fill(isHovering ? Color.cyan : BlueHour.slate.opacity(0.65))
                )

            if chartMode == .gbEnergy,
               let gbPt = (selectedAlpha.flatMap { findClosestGBPoint(to: $0) } ?? runner.currentGBEnergy ?? displayGBPoints.last) {
                Text(gbPt.sigma > 0
                     ? String(format: "Δθ = %.2f° (Σ%d, %d %d 0)", gbPt.misorientationDeg, gbPt.sigma, gbPt.m, gbPt.n)
                     : String(format: "Δθ = %.3f° (θ = ±%.3f°)", gbPt.misorientationDeg, gbPt.halfAngleDeg))
                    .font(.system(size: 11, weight: .semibold, design: .monospaced))
                    .foregroundColor(BlueHour.ice)
                Text("γ_GB(Int) = \(String(format: "%.6f", gbPt.gbEnergyInterior))")
                    .font(.system(size: 11, weight: .bold, design: .monospaced))
                    .foregroundColor(gbLineColor)
                Text("γ_GB(Tot) = \(String(format: "%.6f", gbPt.gbEnergyTotal))")
                    .font(.system(size: 11, design: .monospaced))
                    .foregroundColor(.teal)
                Text("L_GB = \(String(format: "%.0f", gbPt.lengthGBA0))a₀ (\(gbPt.nx)×\(gbPt.ny))")
                    .font(.system(size: 10, design: .monospaced))
                    .foregroundColor(BlueHour.mist)

                Spacer(minLength: 6)

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
            } else if chartMode == .eigenvalue,
                      let eigPt = (selectedAlpha.flatMap { findClosestEigenPoint(to: $0) } ?? runner.eigenPoints.last) {
                let paramLabel = isShiftMode ? "Shift u/a₀" : (isIndentationMode ? "Depth δ" : (hasConjugateStrain ? "Load α" : "Strain α"))
                Text("Step \(eigPt.iteration)")
                    .font(.system(size: 11, weight: .semibold, design: .monospaced))
                    .foregroundColor(BlueHour.ice)
                Text("\(paramLabel) = \(String(format: "%.6f", eigPt.alpha))")
                    .font(.system(size: 11, design: .monospaced))
                    .foregroundColor(BlueHour.mist)
                if let eps = eigPt.strain {
                    Text("ε = \(String(format: "%.6f", eps))")
                        .font(.system(size: 11, design: .monospaced))
                        .foregroundColor(BlueHour.mist)
                }
                if let lb = eigPt.lambdaBox {
                    Text("λ_box = \(String(format: "%.4e", lb))")
                        .font(.system(size: 11, weight: .semibold, design: .monospaced))
                        .foregroundColor(eigPt.isCellUnstable ? .red : eigenLineColor)
                }
                Text("λ₁(min) = \(String(format: "%.5e", eigPt.lambdaMin))")
                    .font(.system(size: 11, weight: .semibold, design: .monospaced))
                    .foregroundColor(eigPt.isUnstable ? .red : eigenLineColor)
                if let l2 = eigPt.lambda2 {
                    Text("λ₂ = \(String(format: "%.5e", l2))")
                        .font(.system(size: 10.5, design: .monospaced))
                        .foregroundColor(BlueHour.mist)
                }
                Spacer(minLength: 6)
                if !isHovering {
                    Text("Hover curve to inspect • Drag to zoom X")
                        .font(.system(size: 10))
                        .foregroundColor(BlueHour.mist.opacity(0.75))
                        .lineLimit(1)
                }
            } else if let pt = (selectedAlpha.flatMap { findClosestPoint(to: $0) } ?? runner.curvePoints.last) {
                let curStress = stressValue(for: pt)
                let stressTitle: String = {
                    if stressComponent == .primary {
                        return isIndentationMode ? "Load P" : (primaryIsTorque ? "Torque T" : "Stress σ")
                    }
                    return stressComponent.shortLabel
                }()
                let xReadout: String = {
                    if isIndentationMode {
                        return "Depth δ = \(String(format: "%.4f", pt.alpha))"
                    } else if isShiftMode {
                        return "Shift u/a₀ = \(String(format: "%.5f", pt.alpha))"
                    } else if hasConjugateStrain {
                        return "Load α = \(String(format: "%.5f", pt.alpha))  ε = \(String(format: "%.5f", pt.effectiveStrain))"
                    } else {
                        return "Strain α = \(String(format: "%.5f", pt.alpha))"
                    }
                }()

                Text("Step \(pt.iteration)")
                    .font(.system(size: 11, weight: .semibold, design: .monospaced))
                    .foregroundColor(BlueHour.ice)
                Text(xReadout)
                    .font(.system(size: 11, design: .monospaced))
                    .foregroundColor(BlueHour.mist)
                Text("\(stressTitle) = \(String(format: "%.6f", curStress))")
                    .font(.system(size: 11, weight: .semibold, design: .monospaced))
                    .foregroundColor(stressLineColor)
                if stressComponent != .vonMisesAvg {
                    Text("⟨σ_vM⟩ = \(String(format: "%.4f", pt.vonMisesAvg))")
                        .font(.system(size: 11, design: .monospaced))
                        .foregroundColor(.orange)
                }
                Text(pt.elasticEnergy.map { "E_el = \(String(format: "%.4e", $0))" }
                     ?? "E = \(String(format: "%.4e", pt.postEnergy))")
                    .font(.system(size: 11, design: .monospaced))
                    .foregroundColor(energyLineColor)
                if pt.shouldRemesh && pt.iteration > 0 {
                    Text("⚡ Surgery")
                        .font(.system(size: 10, weight: .bold))
                        .foregroundColor(.red)
                }

                Spacer(minLength: 6)

                if !isHovering {
                    Text("Hover curve to inspect • Drag to zoom X")
                        .font(.system(size: 10))
                        .foregroundColor(BlueHour.mist.opacity(0.75))
                        .lineLimit(1)
                }
            } else {
                Spacer()
            }
        }
        .lineLimit(1)
        .padding(.horizontal, 14)
        .frame(height: 30)
        .background(BlueHour.controlSurface)
    }
    
    private func hoveredCoordinates(for selX: Double) -> (x: Double, y: Double)? {
        switch chartMode {
        case .gbEnergy:
            guard let pt = findClosestGBPoint(to: selX) else { return nil }
            return (pt.misorientationDeg, gbEnergyValue(for: pt))
        case .eigenvalue:
            guard let pt = findClosestEigenPoint(to: selX), let y = eigenY(pt) else { return nil }
            return (eigenX(pt), y)
        default:
            guard let pt = findClosestPoint(to: selX) else { return nil }
            let y = chartMode == .stress ? stressValue(for: pt) : energyValue(pt)
            return (xValue(pt), y)
        }
    }
    
    private var displayCurvePoints: [StressStrainPoint] {
        let pts = isXZoomed ? runner.curvePoints.filter { inXZoom(xValue($0)) } : runner.curvePoints
        let maxPoints = showPoints ? 450 : 800
        guard pts.count > maxPoints else { return pts }
        let stride = max(1, pts.count / maxPoints)
        var sampled: [StressStrainPoint] = []
        sampled.reserveCapacity(maxPoints + 2)
        for i in Swift.stride(from: 0, to: pts.count, by: stride) {
            sampled.append(pts[i])
        }
        if let last = pts.last, sampled.last?.iteration != last.iteration {
            sampled.append(last)
        }
        return sampled
    }
    
    private var displayRemeshEvents: [StressStrainPoint] {
        let events = isXZoomed ? remeshEvents.filter { inXZoom(xValue($0)) } : remeshEvents
        guard events.count > 400 else { return events }
        let stride = max(1, events.count / 300)
        var sampled: [StressStrainPoint] = []
        for i in Swift.stride(from: 0, to: events.count, by: stride) {
            sampled.append(events[i])
        }
        return sampled
    }
    
    private var displayEigenPoints: [EigenPoint] {
        let pts = isXZoomed ? runner.eigenPoints.filter { inXZoom(eigenX($0)) } : runner.eigenPoints
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
        return runner.curvePoints.min(by: { abs(xValue($0) - alpha) < abs(xValue($1) - alpha) })
    }
    
    private func findClosestEigenPoint(to alpha: Double) -> EigenPoint? {
        guard !runner.eigenPoints.isEmpty else { return nil }
        return runner.eigenPoints.filter { eigenY($0) != nil }
            .min(by: { abs(eigenX($0) - alpha) < abs(eigenX($1) - alpha) })
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
