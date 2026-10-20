import SwiftUI
import AppKit

public struct PathEnergyChartView: View {
    @ObservedObject var manager: AcousticStudyManager
    
    @State private var hoveredScreenX: CGFloat? = nil
    @State private var showAcousticTensorCurve: Bool = true
    
    private var pathPoints: [ShearPathPointModel] {
        manager.shearPath?.points ?? []
    }
    
    private var criticalAlphas: [Double] {
        manager.shearPath?.critical_alphas ?? []
    }
    
    private var criticalPoints: [ShearPathPointModel] {
        manager.shearPath?.critical_points ?? []
    }
    
    private var periodAlpha: Double {
        manager.shearPath?.period_alpha ?? ((manager.latticeSymmetry == .triangular) ? sqrt(4.0 / 3.0) : 1.0)
    }
    
    public init(manager: AcousticStudyManager) {
        self.manager = manager
    }
    
    public var body: some View {
        VStack(spacing: 0) {
            // =========================================================================
            // 1. TOP RESPONSIVE TOOLBAR (Horizontally scrollable so it never clips in Split View)
            // =========================================================================
            ScrollView(.horizontal, showsIndicators: false) {
                HStack(spacing: 8) {
                    VStack(alignment: .leading, spacing: 1) {
                        HStack(spacing: 5) {
                            Image(systemName: "chart.xyaxis.line")
                                .foregroundColor(.accentColor)
                                .font(.system(size: 13, weight: .bold))
                            Text("Path Energy E(α)")
                                .font(.system(size: 13, weight: .bold))
                        }
                        Text(pathSubtitle)
                            .font(.system(size: 10, design: .monospaced))
                            .foregroundColor(.secondary)
                            .lineLimit(1)
                    }
                    .fixedSize()
                    
                    Divider().frame(height: 16)
                    
                    // Dropdown menu for selecting alpha / critical alpha / wells
                    if !criticalAlphas.isEmpty || periodAlpha > 0 {
                        Menu {
                            if !criticalAlphas.isEmpty {
                                Section("Critical Stability Limits (α_c)") {
                                    ForEach(criticalAlphas, id: \.self) { aC in
                                        Button(action: {
                                            manager.hoveredAlpha = aC
                                            manager.selectedCriticalAlpha = aC
                                        }) {
                                            Text(String(format: "α_c = %+.4f (Loss of Ellipticity)", aC))
                                        }
                                    }
                                }
                            }
                            
                            Section("Symmetry Wells") {
                                Button(action: { manager.hoveredAlpha = 0.0 }) {
                                    Text("α = 0.0000 (Reference State)")
                                }
                                let p = periodAlpha
                                let minK = Int(floor(alphaBounds.min / p))
                                let maxK = Int(ceil(alphaBounds.max / p))
                                ForEach(minK...maxK, id: \.self) { k in
                                    if k != 0 {
                                        let wellA = Double(k) * p
                                        if wellA >= alphaBounds.min && wellA <= alphaBounds.max {
                                            Button(action: { manager.hoveredAlpha = wellA }) {
                                                Text(String(format: "α = %+.4f (%dΔα Well)", wellA, k))
                                            }
                                        }
                                    }
                                }
                            }
                        } label: {
                            HStack(spacing: 4) {
                                Image(systemName: "list.bullet")
                                if let hA = manager.hoveredAlpha {
                                    Text(String(format: "α = %+.3f", hA))
                                        .font(.system(size: 11, weight: .semibold, design: .monospaced))
                                } else {
                                    Text("Jump to α")
                                        .font(.system(size: 11, weight: .medium))
                                }
                            }
                        }
                        .menuStyle(.borderedButton)
                        .controlSize(.small)
                        .fixedSize()
                        .help("Jump directly to any critical stability limit α_c or periodic symmetry well")
                    }
                    
                    Divider().frame(height: 14)
                    
                    // Compact Legend
                    HStack(spacing: 8) {
                        HStack(spacing: 3) {
                            Circle().fill(Color(red: 0.2, green: 0.95, blue: 0.4)).frame(width: 6, height: 6)
                            Text("Stable").font(.system(size: 10))
                        }
                        HStack(spacing: 3) {
                            Circle().fill(Color(red: 1.0, green: 0.2, blue: 0.85)).frame(width: 6, height: 6)
                            Text("Unstable").font(.system(size: 10))
                        }
                        HStack(spacing: 3) {
                            Rectangle().fill(Color.yellow).frame(width: 8, height: 2)
                            Text("α_c").font(.system(size: 10))
                        }
                    }
                    .fixedSize()
                    
                    Divider().frame(height: 14)
                    
                    // Toggle Acoustic Tensor Curve
                    Button(action: { showAcousticTensorCurve.toggle() }) {
                        Image(systemName: showAcousticTensorCurve ? "waveform.path.ecg" : "waveform")
                            .foregroundColor(showAcousticTensorCurve ? .accentColor : .secondary)
                    }
                    .buttonStyle(.bordered)
                    .controlSize(.small)
                    .fixedSize()
                    .help(showAcousticTensorCurve ? "Hide min det Ac curve" : "Show min det Ac curve")
                    
                    // Export CSV Button
                    Button(action: exportCSV) {
                        Image(systemName: "square.and.arrow.up")
                    }
                    .buttonStyle(.bordered)
                    .controlSize(.small)
                    .fixedSize()
                    .help("Export energy path data as CSV")
                }
                .padding(.vertical, 1)
            }
            .frame(maxWidth: .infinity, alignment: .leading)
            .padding(.horizontal, 10)
            .padding(.vertical, 6)
            .background(Color(NSColor.controlBackgroundColor))
            
            Divider()
            
            // =========================================================================
            // 3. REAL-TIME TELEMETRY BAR (Horizontally scrollable, zero plot occlusion)
            // =========================================================================
            ScrollView(.horizontal, showsIndicators: false) {
                HStack(spacing: 8) {
                    if let pt = currentInspectedPoint {
                        let isStable = pt.is_stable ?? ((pt.min_detAc ?? 0.0) > 0)
                        HStack(spacing: 4) {
                            Text(String(format: "α = %+.4f", pt.alpha))
                                .font(.system(size: 11, weight: .bold, design: .monospaced))
                                .foregroundColor(.cyan)
                            
                            Text(isStable ? "STABLE" : "UNSTABLE")
                                .font(.system(size: 9, weight: .black))
                                .padding(.horizontal, 4)
                                .padding(.vertical, 1)
                                .background(isStable ? Color.green.opacity(0.25) : Color.red.opacity(0.25), in: Capsule())
                                .foregroundColor(isStable ? .green : .red)
                        }
                        .fixedSize()
                        
                        Divider().frame(height: 12)
                        
                        HStack(spacing: 2) {
                            Text("E:").foregroundColor(.secondary)
                            Text(String(format: "%.4f", pt.energy)).bold()
                        }
                        .font(.system(size: 10, design: .monospaced))
                        .fixedSize()
                        
                        if let det = pt.min_detAc {
                            Divider().frame(height: 12)
                            HStack(spacing: 2) {
                                Text("det Ac:").foregroundColor(.secondary)
                                Text(String(format: "%+.2f", det)).bold()
                                    .foregroundColor(det > 0 ? .green : .red)
                            }
                            .font(.system(size: 10, design: .monospaced))
                            .fixedSize()
                        }
                        
                        Divider().frame(height: 12)
                        
                        HStack(spacing: 2) {
                            Text("Poincaré:").foregroundColor(.secondary)
                            Text(String(format: "(%.2f, %.2f)", pt.x, pt.y))
                        }
                        .font(.system(size: 10, design: .monospaced))
                        .fixedSize()
                    } else {
                        Image(systemName: "hand.point.up.left")
                            .foregroundColor(.secondary)
                            .font(.caption2)
                        Text("Hover or drag across the chart to inspect local stability & Poincaré mapping.")
                            .font(.caption2)
                            .foregroundColor(.secondary)
                            .lineLimit(1)
                            .fixedSize()
                    }
                }
                .padding(.vertical, 1)
            }
            .frame(maxWidth: .infinity, alignment: .leading)
            .padding(.horizontal, 10)
            .padding(.vertical, 4)
            .background(Color(NSColor.controlBackgroundColor).opacity(0.85))
            
            Divider()
            
            // =========================================================================
            // 4. CHART CANVASES
            // =========================================================================
            GeometryReader { geo in
                let width = geo.size.width
                let height = geo.size.height
                
                if pathPoints.isEmpty {
                    VStack(spacing: 10) {
                        Image(systemName: "chart.xyaxis.line")
                            .font(.system(size: 36))
                            .foregroundColor(.secondary)
                        Text("No deformation path calculated.")
                            .font(.headline)
                            .foregroundColor(.secondary)
                        Text("Run study with simple shear to view E(α) and stability limits.")
                            .font(.caption)
                            .foregroundColor(.secondary)
                    }
                    .frame(maxWidth: .infinity, maxHeight: .infinity)
                } else {
                    let eHeight = showAcousticTensorCurve ? max(120, height * 0.64) : height
                    let acHeight = showAcousticTensorCurve ? max(80, height - eHeight - 1) : 0
                    
                    VStack(spacing: 0) {
                        // Energy Chart Canvas
                        Canvas { context, size in
                            drawEnergyChart(context: context, size: size, showBottomLabels: !showAcousticTensorCurve)
                        }
                        .frame(width: width, height: eHeight)
                        
                        if showAcousticTensorCurve {
                            Divider()
                            // Acoustic Tensor min det Ac Canvas
                            Canvas { context, size in
                                drawAcousticTensorChart(context: context, size: size)
                            }
                            .frame(width: width, height: acHeight)
                        }
                    }
                    .contentShape(Rectangle())
                    .gesture(
                        DragGesture(minimumDistance: 0)
                            .onChanged { val in
                                hoveredScreenX = val.location.x
                                updateHoveredAlpha(screenX: val.location.x, width: width)
                            }
                    )
                }
            }
        }
        .background(Color(NSColor.textBackgroundColor))
    }
    
    // MARK: - Subtitle
    private var pathSubtitle: String {
        let ptype = manager.pathType
        let pName = ptype == "pure_shear" ? "Pure Shear F□" : (ptype == "diagonal_shear" ? "Diagonal Shear F◇" : String(format: "Simple Shear (θ=%.0f°)", manager.shearTheta))
        let latName = manager.latticeSymmetry == .triangular ? "triangular" : "square"
        return String(format: "%@ · Δα=%.3f (%@)", pName, periodAlpha, latName)
    }
    
    private var alphaBounds: (min: Double, max: Double) {
        guard !pathPoints.isEmpty else { return (-1.0, 1.0) }
        let aMin = pathPoints.map { $0.alpha }.min() ?? -1.0
        let aMax = pathPoints.map { $0.alpha }.max() ?? 1.0
        return (min: aMin, max: aMax)
    }
    
    private var energyBounds: (min: Double, max: Double) {
        guard !pathPoints.isEmpty else { return (0.0, 1.0) }
        let energies = pathPoints.map { $0.energy }.filter { $0.isFinite }
        let eMin = energies.min() ?? 0.0
        let eMax = energies.max() ?? 1.0
        let pad = max(0.01, (eMax - eMin) * 0.08)
        return (min: max(0.0, eMin - pad * 0.5), max: eMax + pad)
    }
    
    private var detAcBounds: (min: Double, max: Double) {
        guard !pathPoints.isEmpty else { return (-10.0, 10.0) }
        let dets = pathPoints.compactMap { $0.min_detAc }.filter { $0.isFinite }
        let dMin = dets.min() ?? -1.0
        let dMax = dets.max() ?? 1.0
        let pad = max(0.1, (dMax - dMin) * 0.08)
        return (min: dMin - pad, max: dMax + pad)
    }
    
    private var currentInspectedPoint: ShearPathPointModel? {
        guard let hAlpha = manager.hoveredAlpha, !pathPoints.isEmpty else {
            return pathPoints.min(by: { abs($0.alpha) < abs($1.alpha) }) // Default to alpha = 0
        }
        return pathPoints.min(by: { abs($0.alpha - hAlpha) < abs($1.alpha - hAlpha) })
    }
    
    private func updateHoveredAlpha(screenX: CGFloat, width: CGFloat) {
        let marginL: CGFloat = 46
        let marginR: CGFloat = 16
        let plotW = max(10, width - marginL - marginR)
        let clampedX = max(marginL, min(width - marginR, screenX))
        let normX = (clampedX - marginL) / plotW
        let bounds = alphaBounds
        let alpha = bounds.min + Double(normX) * (bounds.max - bounds.min)
        manager.hoveredAlpha = alpha
        if let pt = pathPoints.min(by: { abs($0.alpha - alpha) < abs($1.alpha - alpha) }) {
            manager.hoveredPathPoint = pt
        }
    }
    
    // MARK: - Canvas Rendering: Energy Chart E(α)
    private func drawEnergyChart(context: GraphicsContext, size: CGSize, showBottomLabels: Bool) {
        let marginL: CGFloat = 46
        let marginR: CGFloat = 16
        let marginT: CGFloat = 24
        let marginB: CGFloat = showBottomLabels ? 22 : 8
        
        let plotW = size.width - marginL - marginR
        let plotH = size.height - marginT - marginB
        guard plotW > 20, plotH > 20 else { return }
        
        let aBounds = alphaBounds
        let eBounds = energyBounds
        let aSpan = max(1e-9, aBounds.max - aBounds.min)
        let eSpan = max(1e-9, eBounds.max - eBounds.min)
        
        func xPos(_ a: Double) -> CGFloat {
            marginL + CGFloat((a - aBounds.min) / aSpan) * plotW
        }
        
        func yPos(_ e: Double) -> CGFloat {
            marginT + plotH - CGFloat((e - eBounds.min) / eSpan) * plotH
        }
        
        // Plot Area Frame
        var bgPath = Path()
        bgPath.addRect(CGRect(x: marginL, y: marginT, width: plotW, height: plotH))
        context.fill(bgPath, with: .color(Color(NSColor.textBackgroundColor)))
        context.stroke(bgPath, with: .color(Color.secondary.opacity(0.2)), lineWidth: 1)
        
        // Horizontal grid lines (Energy) - Safe count
        let targetETicks = max(2, min(5, Int(plotH / 38.0)))
        let eStep = niceStep(eSpan, targetSteps: targetETicks)
        var eVal = (ceil((eBounds.min - 1e-9) / eStep)) * eStep
        while eVal <= eBounds.max + 1e-9 {
            let y = yPos(eVal)
            if y >= marginT && y <= marginT + plotH {
                var line = Path()
                line.move(to: CGPoint(x: marginL, y: y))
                line.addLine(to: CGPoint(x: marginL + plotW, y: y))
                let isZero = abs(eVal) < 1e-6
                context.stroke(line, with: .color(Color.secondary.opacity(isZero ? 0.45 : 0.12)),
                               style: StrokeStyle(lineWidth: isZero ? 1.2 : 0.6, dash: isZero ? [] : [3, 3]))
                
                let labelStr = abs(eVal) < 1e-6 ? "0.00" : String(format: "%.2f", eVal)
                let label = Text(labelStr)
                    .font(.system(size: 9, design: .monospaced))
                    .foregroundColor(.secondary)
                context.draw(label, at: CGPoint(x: marginL - 18, y: y))
            }
            eVal += eStep
        }
        
        // Vertical grid lines (Alpha) - Spaced by at least 65px so labels never collide!
        let targetTicks = max(2, min(6, Int(plotW / 65.0)))
        let aStep = niceStep(aSpan, targetSteps: targetTicks)
        var aVal = (ceil((aBounds.min - 1e-9) / aStep)) * aStep
        while aVal <= aBounds.max + 1e-9 {
            let x = xPos(aVal)
            if x >= marginL && x <= marginL + plotW {
                var line = Path()
                line.move(to: CGPoint(x: x, y: marginT))
                line.addLine(to: CGPoint(x: x, y: marginT + plotH))
                let isZero = abs(aVal) < 1e-6
                context.stroke(line, with: .color(Color.secondary.opacity(isZero ? 0.45 : 0.12)),
                               style: StrokeStyle(lineWidth: isZero ? 1.2 : 0.6, dash: isZero ? [] : [3, 3]))
                
                if showBottomLabels {
                    let labelStr = abs(aVal) < 1e-6 ? "0.0" : String(format: "%.2f", aVal)
                    let label = Text(labelStr)
                        .font(.system(size: 9, design: .monospaced))
                        .foregroundColor(.secondary)
                    context.draw(label, at: CGPoint(x: x, y: marginT + plotH + 10))
                }
            }
            aVal += aStep
        }
        
        // Symmetry Well Markers (dotted blue lines at k * periodAlpha)
        let p = periodAlpha
        var k = floor(aBounds.min / p)
        while k * p <= aBounds.max {
            let wellAlpha = k * p
            if wellAlpha >= aBounds.min && wellAlpha <= aBounds.max {
                let x = xPos(wellAlpha)
                var line = Path()
                line.move(to: CGPoint(x: x, y: marginT))
                line.addLine(to: CGPoint(x: x, y: marginT + plotH))
                context.stroke(line, with: .color(Color.blue.opacity(0.4)),
                               style: StrokeStyle(lineWidth: 1.0, dash: [2, 3]))
                let tag = Text(k == 0 ? "Well (0)" : String(format: "%dΔα", Int(k)))
                    .font(.system(size: 8, weight: .bold, design: .monospaced))
                    .foregroundColor(Color.blue.opacity(0.75))
                context.draw(tag, at: CGPoint(x: x, y: marginT - 8))
            }
            k += 1.0
        }
        
        // Critical Alpha lines α_c (Loss of Ellipticity)
        for (idx, alphaC) in criticalAlphas.enumerated() {
            if alphaC >= aBounds.min && alphaC <= aBounds.max {
                let x = xPos(alphaC)
                var line = Path()
                line.move(to: CGPoint(x: x, y: marginT))
                line.addLine(to: CGPoint(x: x, y: marginT + plotH))
                context.stroke(line, with: .color(Color.yellow),
                               style: StrokeStyle(lineWidth: 1.5, dash: [4, 3]))
                
                // Stagger badge height alternating to prevent collision if two alpha_c are close!
                let badgeY = (idx % 2 == 0) ? (marginT - 6) : (marginT - 15)
                let badge = Text(String(format: "αc=%.3f", alphaC))
                    .font(.system(size: 8, weight: .heavy, design: .monospaced))
                    .foregroundColor(.yellow)
                context.draw(badge, at: CGPoint(x: x, y: badgeY))
            }
        }
        
        // Energy Path Curves colored by acoustic stability
        guard pathPoints.count > 1 else { return }
        for i in 1..<pathPoints.count {
            let p1 = pathPoints[i - 1]
            let p2 = pathPoints[i]
            
            let isStable = (p1.is_stable ?? true) && (p2.is_stable ?? true)
            let color = isStable ? Color(red: 0.2, green: 0.95, blue: 0.4) : Color(red: 1.0, green: 0.2, blue: 0.85)
            
            var seg = Path()
            seg.move(to: CGPoint(x: xPos(p1.alpha), y: yPos(p1.energy)))
            seg.addLine(to: CGPoint(x: xPos(p2.alpha), y: yPos(p2.energy)))
            
            // Crisp segment line
            context.stroke(seg, with: .color(color),
                           style: StrokeStyle(lineWidth: 2.5, lineCap: .round))
        }
        
        // Scrubber Line
        if let hAlpha = manager.hoveredAlpha, hAlpha >= aBounds.min && hAlpha <= aBounds.max {
            let hx = xPos(hAlpha)
            var scrubLine = Path()
            scrubLine.move(to: CGPoint(x: hx, y: marginT))
            scrubLine.addLine(to: CGPoint(x: hx, y: marginT + plotH))
            context.stroke(scrubLine, with: .color(Color.cyan),
                           style: StrokeStyle(lineWidth: 1.5, dash: [3, 2]))
            
            if let hPt = currentInspectedPoint {
                let hy = yPos(hPt.energy)
                let r: CGFloat = 4.5
                context.fill(Path(ellipseIn: CGRect(x: hx - r, y: hy - r, width: 2 * r, height: 2 * r)), with: .color(Color.cyan))
                context.stroke(Path(ellipseIn: CGRect(x: hx - r, y: hy - r, width: 2 * r, height: 2 * r)), with: .color(.black), lineWidth: 1.5)
            }
        }
        
        // Small Y-axis Tag
        let yTag = Text("E(α)")
            .font(.system(size: 9, weight: .bold))
            .foregroundColor(.secondary)
        context.draw(yTag, at: CGPoint(x: marginL - 18, y: marginT - 6))
    }
    
    // MARK: - Canvas Rendering: Acoustic Tensor min det Ac Chart
    private func drawAcousticTensorChart(context: GraphicsContext, size: CGSize) {
        let marginL: CGFloat = 46
        let marginR: CGFloat = 16
        let marginT: CGFloat = 10
        let marginB: CGFloat = 22
        
        let plotW = size.width - marginL - marginR
        let plotH = size.height - marginT - marginB
        guard plotW > 20, plotH > 20 else { return }
        
        let aBounds = alphaBounds
        let dBounds = detAcBounds
        let aSpan = max(1e-9, aBounds.max - aBounds.min)
        let dSpan = max(1e-9, dBounds.max - dBounds.min)
        
        func xPos(_ a: Double) -> CGFloat {
            marginL + CGFloat((a - aBounds.min) / aSpan) * plotW
        }
        
        func yPos(_ d: Double) -> CGFloat {
            marginT + plotH - CGFloat((d - dBounds.min) / dSpan) * plotH
        }
        
        // Plot Area Frame
        var bgPath = Path()
        bgPath.addRect(CGRect(x: marginL, y: marginT, width: plotW, height: plotH))
        context.fill(bgPath, with: .color(Color(NSColor.textBackgroundColor)))
        context.stroke(bgPath, with: .color(Color.secondary.opacity(0.2)), lineWidth: 1)
        
        // Zero Line for Loss of Ellipticity
        if 0.0 >= dBounds.min && 0.0 <= dBounds.max {
            let y0 = yPos(0.0)
            var zeroLine = Path()
            zeroLine.move(to: CGPoint(x: marginL, y: y0))
            zeroLine.addLine(to: CGPoint(x: marginL + plotW, y: y0))
            context.stroke(zeroLine, with: .color(Color.red.opacity(0.8)), lineWidth: 1.2)
            
            let tag = Text("det Ac = 0")
                .font(.system(size: 8, weight: .bold))
                .foregroundColor(.red)
            context.draw(tag, at: CGPoint(x: marginL + plotW - 32, y: y0 - 6))
        }
        
        // Vertical grid lines (Alpha) with safe non-overlapping spacing
        let targetTicks = max(2, min(6, Int(plotW / 65.0)))
        let aStep = niceStep(aSpan, targetSteps: targetTicks)
        var aVal = (ceil((aBounds.min - 1e-9) / aStep)) * aStep
        while aVal <= aBounds.max + 1e-9 {
            let x = xPos(aVal)
            if x >= marginL && x <= marginL + plotW {
                var line = Path()
                line.move(to: CGPoint(x: x, y: marginT))
                line.addLine(to: CGPoint(x: x, y: marginT + plotH))
                let isZero = abs(aVal) < 1e-6
                context.stroke(line, with: .color(Color.secondary.opacity(isZero ? 0.45 : 0.12)),
                               style: StrokeStyle(lineWidth: isZero ? 1.2 : 0.6, dash: isZero ? [] : [3, 3]))
                
                let labelStr = abs(aVal) < 1e-6 ? "0.0" : String(format: "%.2f", aVal)
                let label = Text(labelStr)
                    .font(.system(size: 9, design: .monospaced))
                    .foregroundColor(.secondary)
                context.draw(label, at: CGPoint(x: x, y: marginT + plotH + 10))
            }
            aVal += aStep
        }
        
        // Critical Alpha vertical lines
        for alphaC in criticalAlphas {
            if alphaC >= aBounds.min && alphaC <= aBounds.max {
                let x = xPos(alphaC)
                var line = Path()
                line.move(to: CGPoint(x: x, y: marginT))
                line.addLine(to: CGPoint(x: x, y: marginT + plotH))
                context.stroke(line, with: .color(Color.yellow),
                               style: StrokeStyle(lineWidth: 1.2, dash: [3, 3]))
            }
        }
        
        // min det Ac Curve
        guard pathPoints.count > 1 else { return }
        for i in 1..<pathPoints.count {
            let p1 = pathPoints[i - 1]
            let p2 = pathPoints[i]
            guard let d1 = p1.min_detAc, let d2 = p2.min_detAc else { continue }
            
            let isStable = d1 > 0 && d2 > 0
            let color = isStable ? Color(red: 0.2, green: 0.85, blue: 0.4) : Color(red: 1.0, green: 0.2, blue: 0.4)
            
            var seg = Path()
            seg.move(to: CGPoint(x: xPos(p1.alpha), y: yPos(d1)))
            seg.addLine(to: CGPoint(x: xPos(p2.alpha), y: yPos(d2)))
            context.stroke(seg, with: .color(color), style: StrokeStyle(lineWidth: 2.0, lineCap: .round))
        }
        
        // Scrubber line
        if let hAlpha = manager.hoveredAlpha, hAlpha >= aBounds.min && hAlpha <= aBounds.max {
            let hx = xPos(hAlpha)
            var scrubLine = Path()
            scrubLine.move(to: CGPoint(x: hx, y: marginT))
            scrubLine.addLine(to: CGPoint(x: hx, y: marginT + plotH))
            context.stroke(scrubLine, with: .color(Color.cyan),
                           style: StrokeStyle(lineWidth: 1.2, dash: [3, 2]))
        }
        
        // Y-axis Label
        let yLabel = Text("det Ac")
            .font(.system(size: 9, weight: .bold))
            .foregroundColor(.secondary)
        context.draw(yLabel, at: CGPoint(x: marginL - 18, y: marginT + 6))
    }
    
    // MARK: - CSV Export Action
    private func exportCSV() {
        let panel = NSSavePanel()
        panel.title = "Save Deformation Energy Path CSV"
        panel.nameFieldStringValue = "energy_path.csv"
        panel.allowedContentTypes = [.commaSeparatedText, .plainText]
        if panel.runModal() == .OK, let url = panel.url {
            var text = "# alpha,energy,min_detAc,is_stable,x,y,x_red,y_red,c11,c22,c12\n"
            for pt in pathPoints {
                text += String(format: "%.8e,%.8e,%.8e,%d,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e\n",
                               pt.alpha, pt.energy, pt.min_detAc ?? 0.0, (pt.is_stable ?? true) ? 1 : 0,
                               pt.x, pt.y, pt.x_red, pt.y_red, pt.c11, pt.c22, pt.c12)
            }
            try? text.write(to: url, atomically: true, encoding: .utf8)
        }
    }
    
    private func niceStep(_ span: Double, targetSteps: Int) -> Double {
        guard targetSteps > 0, span > 1e-12 else { return 1.0 }
        let rawStep = span / Double(targetSteps)
        let mag = pow(10.0, floor(log10(rawStep)))
        let norm = rawStep / mag
        let stepMultiplier: Double
        if norm < 1.5 { stepMultiplier = 1.0 }
        else if norm < 3.5 { stepMultiplier = 2.0 }
        else if norm < 7.5 { stepMultiplier = 5.0 }
        else { stepMultiplier = 10.0 }
        return stepMultiplier * mag
    }
}
