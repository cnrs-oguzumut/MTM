import SwiftUI
import AppKit

public struct PoincareStabilityMapView: View {
    @ObservedObject var manager: AcousticStudyManager
    
    // Zoom and pan state
    @State private var scale: CGFloat = 1.0
    @State private var lastScale: CGFloat = 1.0
    @State private var offset: CGSize = .zero
    @State private var lastOffset: CGSize = .zero
    
    // Hover inspection state
    @State private var hoverLocation: CGPoint? = nil
    @State private var hoveredPoint: AcousticPointModel? = nil
    
    // Color scaling
    private var minDetVal: Double { manager.summary.min_detAc_min }
    private var maxDetVal: Double { manager.summary.min_detAc_max }
    private var energyBounds: (min: Double, max: Double, keptCount: Int, totalCount: Int) {
        manager.effectiveEnergyBounds
    }
    
    public var body: some View {
        VStack(spacing: 0) {
            // =========================================================================
            // 1. TOP DEDICATED TOOLBAR (Fixed height, completely ABOVE canvas)
            // =========================================================================
            HStack(spacing: 6) {
                // Coordinate Mode Switcher
                Picker("", selection: $manager.useReducedCoordinates) {
                    Text("Original").tag(false)
                    Text("Reduced").tag(true)
                }
                .labelsHidden()
                .pickerStyle(.segmented)
                .controlSize(.small)
                .frame(width: 130)
                
                // Color Field Mode (Stability vs Energy)
                Picker("", selection: $manager.colorField) {
                    ForEach(ColorFieldMode.allCases) { mode in
                        Text(mode.shortLabel).tag(mode)
                    }
                }
                .labelsHidden()
                .pickerStyle(.segmented)
                .controlSize(.small)
                .frame(width: 125)
                
                // Render Style Mode (Smooth Surface vs Scatter)
                Picker("", selection: $manager.renderStyle) {
                    Text("Surface").tag(RenderStyleMode.smoothSurface)
                    Text("Scatter").tag(RenderStyleMode.scatterPoints)
                }
                .labelsHidden()
                .pickerStyle(.segmented)
                .controlSize(.small)
                .frame(width: 115)
                
                Spacer(minLength: 4)
                
                // Display Options Menu
                Menu {
                    Toggle("Zero Boundary Contour", isOn: $manager.highlightZeroBoundary)
                    Toggle("Deformation Path", isOn: $manager.showShearPath)
                    if manager.renderStyle == .scatterPoints {
                        Divider()
                        Text("Point Radius: \(manager.pointRadius, specifier: "%.1f")")
                    }
                } label: {
                    Image(systemName: "slider.horizontal.2.square")
                }
                .menuStyle(.borderlessButton)
                .frame(width: 26)
                .help("Display options (Zero contour, Deformation path)")
                
                // Reset View Button
                Button(action: resetView) {
                    Image(systemName: "arrow.counterclockwise")
                }
                .buttonStyle(.bordered)
                .controlSize(.small)
                .help("Reset Zoom & Pan")
            }
            .padding(.horizontal, 10)
            .padding(.vertical, 6)
            .background(Color(NSColor.controlBackgroundColor))
            
            // =========================================================================
            // 1B. ENERGY COLOR SCALE & PERCENTILE FILTER TOOLBAR (When Energy Φ(C) is active)
            // =========================================================================
            if manager.colorField == .strainEnergy {
                Divider()
                VStack(spacing: 5) {
                    HStack(spacing: 8) {
                        // Colormap Picker
                        HStack(spacing: 4) {
                            Image(systemName: "paintpalette.fill")
                                .foregroundColor(.orange)
                                .font(.caption)
                            Picker("", selection: $manager.energyColormap) {
                                ForEach(EnergyColormapPreset.allCases) { cmap in
                                    Text(cmap.rawValue).tag(cmap)
                                }
                            }
                            .labelsHidden()
                            .controlSize(.small)
                            .frame(width: 95)
                            
                            Toggle("Invert", isOn: $manager.invertEnergyColormap)
                                .toggleStyle(.checkbox)
                                .controlSize(.small)
                                .font(.caption2)
                        }
                        
                        Divider().frame(height: 16)
                        
                        // Percentile Cutoff Slider + Presets
                        if !manager.useCustomEnergyRange {
                            HStack(spacing: 5) {
                                Text("P_max:")
                                    .font(.caption.bold())
                                Slider(value: $manager.energyUpperPercentile, in: 10.0...100.0, step: 0.5)
                                    .controlSize(.small)
                                    .frame(width: 100)
                                Text(String(format: "%.1f%%", manager.energyUpperPercentile))
                                    .font(.system(size: 11, weight: .bold, design: .monospaced))
                                    .foregroundColor(.orange)
                                    .frame(width: 44, alignment: .trailing)
                                
                                ForEach([100.0, 99.0, 95.0, 90.0, 80.0, 65.0], id: \.self) { pct in
                                    Button(String(format: "%.0f%%", pct)) {
                                        manager.energyUpperPercentile = pct
                                    }
                                    .buttonStyle(.bordered)
                                    .controlSize(.mini)
                                }
                            }
                        } else {
                            HStack(spacing: 6) {
                                Text("Min:")
                                    .font(.caption)
                                DoubleField(title: "0.0", value: $manager.customEnergyMin, width: 60)
                                Text("Max:")
                                    .font(.caption)
                                DoubleField(title: "10.0", value: $manager.customEnergyMax, width: 60)
                            }
                        }
                        
                        Toggle("Custom", isOn: $manager.useCustomEnergyRange)
                            .toggleStyle(.checkbox)
                            .controlSize(.small)
                            .font(.caption2)
                        
                        Spacer(minLength: 4)
                    }
                    
                    HStack(spacing: 8) {
                        // Outlier Filtering Mode: Hide/Discard vs Clamp
                        Picker("", selection: $manager.discardAbovePercentile) {
                            Text("Hide > Cutoff").tag(true)
                            Text("Clamp Color").tag(false)
                        }
                        .labelsHidden()
                        .pickerStyle(.segmented)
                        .controlSize(.small)
                        .frame(width: 195)
                        
                        Divider().frame(height: 16)
                        
                        // Gamma / Contrast Slider
                        HStack(spacing: 5) {
                            Text("γ:")
                                .font(.caption.bold())
                            Slider(value: $manager.energyGamma, in: 0.15...2.50, step: 0.05)
                                .controlSize(.small)
                                .frame(width: 85)
                            Text(String(format: "%.2f", manager.energyGamma))
                                .font(.system(size: 11, weight: .semibold, design: .monospaced))
                                .frame(width: 34)
                            Button("0.5") { manager.energyGamma = 0.50 }
                                .buttonStyle(.bordered).controlSize(.mini)
                            Button("0.7") { manager.energyGamma = 0.70 }
                                .buttonStyle(.bordered).controlSize(.mini)
                            Button("1.0") { manager.energyGamma = 1.00 }
                                .buttonStyle(.bordered).controlSize(.mini)
                        }
                        
                        Spacer(minLength: 4)
                        
                        // Live Cutoff Readout Badge
                        let b = energyBounds
                        Text(String(format: "ΔΦ ∈ [%.2f, %.2f] (%d/%d)", b.min, b.max, b.keptCount, b.totalCount))
                            .font(.system(size: 10, weight: .medium, design: .monospaced))
                            .foregroundColor(.secondary)
                            .lineLimit(1)
                    }
                }
                .padding(.horizontal, 12)
                .padding(.vertical, 6)
                .background(Color(NSColor.controlBackgroundColor).opacity(0.95))
            }
            
            Divider()
            
            // =========================================================================
            // 2. CENTER: POINCARÉ DISK CANVAS (Clean, dedicated viewport)
            // =========================================================================
            GeometryReader { geo in
                let center = CGPoint(x: geo.size.width / 2 + offset.width, y: geo.size.height / 2 + offset.height)
                // Use 0.40 of the smaller dimension so the disk never touches the borders
                let baseRadius = min(geo.size.width, geo.size.height) * 0.40 * scale
                
                ZStack {
                    // Dark Canvas Backdrop
                    Color(NSColor.underPageBackgroundColor)
                        .edgesIgnoringSafeArea(.all)
                    
                    Canvas { context, size in
                        drawBackgroundGrid(context: context, center: center, diskRadius: baseRadius)
                        drawSurfaceOrPoints(context: context, center: center, diskRadius: baseRadius)
                        drawZeroBoundaryContours(context: context, center: center, diskRadius: baseRadius)
                        drawShearPath(context: context, center: center, diskRadius: baseRadius)
                        drawHoverIndicator(context: context, center: center, diskRadius: baseRadius)
                    }
                    .gesture(
                        DragGesture()
                            .onChanged { val in
                                offset = CGSize(
                                    width: lastOffset.width + val.translation.width,
                                    height: lastOffset.height + val.translation.height
                                )
                            }
                            .onEnded { _ in
                                lastOffset = offset
                            }
                    )
                    .simultaneousGesture(
                        MagnificationGesture()
                            .onChanged { val in
                                scale = max(0.5, min(10.0, lastScale * val))
                            }
                            .onEnded { _ in
                                lastScale = scale
                            }
                    )
                    .onContinuousHover { phase in
                        switch phase {
                        case .active(let loc):
                            hoverLocation = loc
                            updateHoveredPoint(at: loc, center: center, diskRadius: baseRadius)
                        case .ended:
                            hoverLocation = nil
                            hoveredPoint = nil
                        }
                    }
                    
                }
            }
            
            Divider()
            
            // =========================================================================
            // 3. BOTTOM DEDICATED FOOTER (Fixed height, completely BELOW canvas)
            // =========================================================================
            HStack(spacing: 16) {
                let b = energyBounds
                // Colorbar Legend
                ColorbarLegendView(
                    colorField: manager.colorField,
                    minDet: minDetVal,
                    maxDet: maxDetVal,
                    minEnergy: b.min,
                    maxEnergy: b.max,
                    colormap: manager.energyColormap,
                    invertColormap: manager.invertEnergyColormap,
                    gamma: manager.energyGamma,
                    upperPercentile: manager.useCustomEnergyRange ? nil : manager.energyUpperPercentile
                )
                .frame(width: 310)
                
                Divider()
                    .frame(height: 38)
                
                // Hover Inspector or Summary Info
                if let pt = hoveredPoint {
                    HoverInfoInlineView(point: pt, useReduced: manager.useReducedCoordinates)
                } else {
                    HStack(spacing: 12) {
                        Image(systemName: "hand.point.up.left")
                            .foregroundColor(.secondary)
                        VStack(alignment: .leading, spacing: 2) {
                            Text("\(manager.points.count) data points plotted (\(manager.renderStyle == .smoothSurface ? "Smooth Mesh" : "Scatter Points"))")
                                .font(.caption.bold())
                            Text("Hover cursor over Poincaré disk to inspect local metric, energy & stability")
                                .font(.caption2)
                                .foregroundColor(.secondary)
                        }
                        Spacer()
                    }
                }
            }
            .padding(.horizontal, 16)
            .padding(.vertical, 8)
            .background(Color(NSColor.controlBackgroundColor))
        }
    }
    
    private func resetView() {
        withAnimation(.easeOut(duration: 0.25)) {
            scale = 1.0
            lastScale = 1.0
            offset = .zero
            lastOffset = .zero
        }
    }
    
    // MARK: - Drawing Routines
    
    private func drawBackgroundGrid(context: GraphicsContext, center: CGPoint, diskRadius: CGFloat) {
        // Disk Interior
        let diskRect = CGRect(x: center.x - diskRadius, y: center.y - diskRadius, width: diskRadius * 2, height: diskRadius * 2)
        let diskPath = Path(ellipseIn: diskRect)
        
        context.fill(diskPath, with: .color(Color(NSColor.textBackgroundColor).opacity(0.85)))
        
        // Concentric Hyperbolic Circles: r = tanh(t/2) for t = 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4
        let tValues: [(t: Double, showLabel: Bool)] = [
            (0.2, false), (0.4, true), (0.6, false), (0.8, true), (1.0, false), (1.2, true), (1.4, false)
        ]
        let diag = CGFloat(0.70710678)
        for item in tValues {
            let r_hyp = CGFloat(tanh(item.t / 2.0)) * diskRadius
            if r_hyp < diskRadius {
                let rRect = CGRect(x: center.x - r_hyp, y: center.y - r_hyp, width: r_hyp * 2, height: r_hyp * 2)
                context.stroke(Path(ellipseIn: rRect), with: .color(Color.gray.opacity(0.22)), style: StrokeStyle(lineWidth: 1, dash: [4, 4]))
                
                // Label only alternate circles along the 45° diagonal so they never overlap
                if item.showLabel && diskRadius > 110 {
                    let text = Text(String(format: "t=%.1f", item.t))
                        .font(.system(size: 9, weight: .regular, design: .monospaced))
                        .foregroundColor(.gray.opacity(0.75))
                    let lx = center.x - r_hyp * diag
                    let ly = center.y - r_hyp * diag
                    context.draw(text, at: CGPoint(x: lx, y: ly))
                }
            }
        }
        
        // Coordinate Axes: x (horizontal) and y (vertical)
        var xAxes = Path()
        xAxes.move(to: CGPoint(x: center.x - diskRadius, y: center.y))
        xAxes.addLine(to: CGPoint(x: center.x + diskRadius, y: center.y))
        context.stroke(xAxes, with: .color(Color.gray.opacity(0.35)), lineWidth: 1)
        
        var yAxes = Path()
        yAxes.move(to: CGPoint(x: center.x, y: center.y - diskRadius))
        yAxes.addLine(to: CGPoint(x: center.x, y: center.y + diskRadius))
        context.stroke(yAxes, with: .color(Color.gray.opacity(0.35)), lineWidth: 1)
        
        // Axis Labels (properly positioned with clear margins)
        let xLabel = Text("x = t·(c₁₁−c₂₂)/2").font(.system(size: 11, weight: .semibold)).foregroundColor(.secondary)
        context.draw(xLabel, at: CGPoint(x: center.x + diskRadius - 55, y: center.y + 14))
        
        let yLabel = Text("y = t·c₁₂").font(.system(size: 11, weight: .semibold)).foregroundColor(.secondary)
        context.draw(yLabel, at: CGPoint(x: center.x + 30, y: center.y - diskRadius + 14))
        
        // Unit Disk Boundary (r = 1.0)
        context.stroke(diskPath, with: .color(Color.blue.opacity(0.9)), lineWidth: 2)
    }
    
    private func drawSurfaceOrPoints(context: GraphicsContext, center: CGPoint, diskRadius: CGFloat) {
        guard !manager.points.isEmpty else { return }
        
        let minZ = min(0.0, minDetVal)
        let maxZ = max(1.0, maxDetVal)
        let spanZ = max(1e-6, max(abs(minZ), abs(maxZ)))
        
        let eBounds = energyBounds
        let minE = eBounds.min
        let maxE = max(minE + 1e-6, eBounds.max)
        let shouldFilterHighEnergy = (manager.colorField == .strainEnergy) && manager.discardAbovePercentile && (manager.useCustomEnergyRange || manager.energyUpperPercentile < 99.99)
        
        if manager.renderStyle == .smoothSurface && !manager.surfaceCells.isEmpty {
            // Draw filled quad surface mesh
            for cell in manager.surfaceCells {
                if shouldFilterHighEnergy && cell.avgEnergy > maxE {
                    continue
                }
                let pts = manager.useReducedCoordinates ? cell.ptsRed : cell.ptsOrig
                guard pts.count >= 3 else { continue }
                
                var path = Path()
                let p0 = CGPoint(x: center.x + pts[0].x * diskRadius, y: center.y - pts[0].y * diskRadius)
                path.move(to: p0)
                for k in 1..<pts.count {
                    let pk = CGPoint(x: center.x + pts[k].x * diskRadius, y: center.y - pts[k].y * diskRadius)
                    path.addLine(to: pk)
                }
                path.closeSubpath()
                
                let color: Color
                if manager.colorField == .acousticStability {
                    color = colorForDetAc(cell.avgDetAc, minVal: minZ, maxVal: maxZ, span: spanZ)
                } else {
                    color = colorForEnergy(cell.avgEnergy, minVal: minE, maxVal: maxE)
                }
                
                context.fill(path, with: .color(color))
            }
        } else {
            // Draw discrete point scatter
            let ptRadius = CGFloat(manager.pointRadius)
            
            for pt in manager.points {
                if shouldFilterHighEnergy && pt.energy > maxE {
                    continue
                }
                let px = manager.useReducedCoordinates ? pt.x_red : pt.x
                let py = manager.useReducedCoordinates ? pt.y_red : pt.y
                
                let sx = center.x + CGFloat(px) * diskRadius
                let sy = center.y - CGFloat(py) * diskRadius
                
                let color: Color
                if manager.colorField == .acousticStability {
                    color = colorForDetAc(pt.min_detAc, minVal: minZ, maxVal: maxZ, span: spanZ)
                } else {
                    color = colorForEnergy(pt.energy, minVal: minE, maxVal: maxE)
                }
                
                let rect = CGRect(x: sx - ptRadius, y: sy - ptRadius, width: ptRadius * 2, height: ptRadius * 2)
                context.fill(Path(ellipseIn: rect), with: .color(color))
            }
        }
    }
    
    private func drawZeroBoundaryContours(context: GraphicsContext, center: CGPoint, diskRadius: CGFloat) {
        guard manager.highlightZeroBoundary, !manager.zeroContourSegments.isEmpty else { return }
        
        var contourPath = Path()
        for seg in manager.zeroContourSegments {
            let sx = manager.useReducedCoordinates ? seg.startXRed : seg.startX
            let sy = manager.useReducedCoordinates ? seg.startYRed : seg.startY
            let ex = manager.useReducedCoordinates ? seg.endXRed : seg.endX
            let ey = manager.useReducedCoordinates ? seg.endYRed : seg.endY
            
            let pA = CGPoint(x: center.x + sx * diskRadius, y: center.y - sy * diskRadius)
            let pB = CGPoint(x: center.x + ex * diskRadius, y: center.y - ey * diskRadius)
            
            contourPath.move(to: pA)
            contourPath.addLine(to: pB)
        }
        
        // High contrast understroke (dark outline)
        context.stroke(contourPath, with: .color(Color.black.opacity(0.85)), lineWidth: 3.5)
        // Crisp bright white zero contour stroke
        context.stroke(contourPath, with: .color(Color.white), lineWidth: 2.0)
    }
    
    private func drawShearPath(context: GraphicsContext, center: CGPoint, diskRadius: CGFloat) {
        guard manager.showShearPath, let sp = manager.shearPath, sp.points.count > 1 else { return }
        let useRed = manager.useReducedCoordinates
        
        func screen(_ p: ShearPathPointModel) -> CGPoint {
            let px = useRed ? p.x_red : p.x
            let py = useRed ? p.y_red : p.y
            return CGPoint(x: center.x + CGFloat(px) * diskRadius, y: center.y - CGFloat(py) * diskRadius)
        }
        
        // Segment-by-segment so color follows local stability; in reduced coordinates the
        // path folds back at fundamental-domain boundaries, so skip large jumps.
        for i in 1..<sp.points.count {
            let a = sp.points[i - 1]
            let b = sp.points[i]
            let dx = (useRed ? b.x_red - a.x_red : b.x - a.x)
            let dy = (useRed ? b.y_red - a.y_red : b.y - a.y)
            if dx * dx + dy * dy > 0.01 { continue }
            var seg = Path()
            seg.move(to: screen(a))
            seg.addLine(to: screen(b))
            let stable = (a.is_stable ?? true) && (b.is_stable ?? true)
            context.stroke(seg, with: .color(.black.opacity(0.8)), style: StrokeStyle(lineWidth: 5, lineCap: .round))
            context.stroke(seg, with: .color(stable ? Color(red: 0.2, green: 0.95, blue: 0.4) : Color(red: 1.0, green: 0.2, blue: 0.85)),
                           style: StrokeStyle(lineWidth: 2.8, lineCap: .round))
        }
        
        // Reference state (α = 0)
        if let ref = sp.points.min(by: { abs($0.alpha) < abs($1.alpha) }) {
            let s = screen(ref)
            let r: CGFloat = 5
            context.fill(Path(ellipseIn: CGRect(x: s.x - r, y: s.y - r, width: 2 * r, height: 2 * r)), with: .color(.white))
            context.stroke(Path(ellipseIn: CGRect(x: s.x - r, y: s.y - r, width: 2 * r, height: 2 * r)), with: .color(.black), lineWidth: 1.5)
            context.draw(Text("α=0").font(.system(size: 10, weight: .bold)).foregroundColor(.white), at: CGPoint(x: s.x + 18, y: s.y - 10))
        }
        
        // Critical points α_c (loss of ellipticity)
        for cp in sp.critical_points {
            let s = screen(cp)
            let r: CGFloat = 7
            let rect = CGRect(x: s.x - r, y: s.y - r, width: 2 * r, height: 2 * r)
            context.fill(Path(ellipseIn: rect), with: .color(.yellow))
            context.stroke(Path(ellipseIn: rect), with: .color(.black), lineWidth: 2)
            let label = Text(String(format: "αc=%.4f", cp.alpha))
                .font(.system(size: 11, weight: .heavy, design: .monospaced))
                .foregroundColor(.yellow)
            context.draw(label, at: CGPoint(x: s.x + 46, y: s.y + (cp.alpha < 0 ? 14 : -14)))
        }
        
        // Synchronized hover cursor from 1D energy chart
        if let hAlpha = manager.hoveredAlpha, let pt = sp.points.min(by: { abs($0.alpha - hAlpha) < abs($1.alpha - hAlpha) }) {
            let s = screen(pt)
            let r: CGFloat = 8
            let haloRect = CGRect(x: s.x - r - 4, y: s.y - r - 4, width: 2 * (r + 4), height: 2 * (r + 4))
            context.fill(Path(ellipseIn: haloRect), with: .color(Color.cyan.opacity(0.35)))
            context.fill(Path(ellipseIn: CGRect(x: s.x - r, y: s.y - r, width: 2 * r, height: 2 * r)), with: .color(.cyan))
            context.stroke(Path(ellipseIn: CGRect(x: s.x - r, y: s.y - r, width: 2 * r, height: 2 * r)), with: .color(.black), lineWidth: 2)
            let tag = Text(String(format: "α=%.3f (E=%.3f)", pt.alpha, pt.energy))
                .font(.system(size: 11, weight: .bold, design: .monospaced))
                .foregroundColor(.cyan)
            context.draw(tag, at: CGPoint(x: s.x + 55, y: s.y - 14))
        }
    }
    
    private func drawHoverIndicator(context: GraphicsContext, center: CGPoint, diskRadius: CGFloat) {
        guard let pt = hoveredPoint else { return }
        let px = manager.useReducedCoordinates ? pt.x_red : pt.x
        let py = manager.useReducedCoordinates ? pt.y_red : pt.y
        let sx = center.x + CGFloat(px) * diskRadius
        let sy = center.y - CGFloat(py) * diskRadius
        
        let ringRect = CGRect(x: sx - 10, y: sy - 10, width: 20, height: 20)
        context.stroke(Path(ellipseIn: ringRect), with: .color(Color.accentColor), lineWidth: 3.0)
        
        let dotRect = CGRect(x: sx - 3, y: sy - 3, width: 6, height: 6)
        context.fill(Path(ellipseIn: dotRect), with: .color(.white))
    }
    
    private func updateHoveredPoint(at screenPos: CGPoint, center: CGPoint, diskRadius: CGFloat) {
        guard !manager.points.isEmpty else { return }
        
        var closestPt: AcousticPointModel? = nil
        var minDist: CGFloat = 28.0
        
        for pt in manager.points {
            let px = manager.useReducedCoordinates ? pt.x_red : pt.x
            let py = manager.useReducedCoordinates ? pt.y_red : pt.y
            let sx = center.x + CGFloat(px) * diskRadius
            let sy = center.y - CGFloat(py) * diskRadius
            
            let dx = sx - screenPos.x
            let dy = sy - screenPos.y
            let dist = sqrt(dx * dx + dy * dy)
            if dist < minDist {
                minDist = dist
                closestPt = pt
            }
        }
        
        self.hoveredPoint = closestPt
    }
    
    // MARK: - Colormap Helpers
    
    private func colorForDetAc(_ val: Double, minVal: Double, maxVal: Double, span: Double) -> Color {
        // Diverging colormap:
        // val < 0: Red / Crimson (loss of ellipticity)
        // val = 0: Yellow / Gold boundary
        // val > 0: Cyan / Blue (stable)
        if val < 0.0 {
            let norm = min(1.0, abs(val) / span)
            return Color(red: 0.95 + 0.05 * norm,
                         green: max(0.08, 0.65 * (1.0 - norm)),
                         blue: max(0.08, 0.2 * (1.0 - norm)),
                         opacity: 0.90)
        } else {
            let norm = min(1.0, val / span)
            return Color(red: max(0.05, 0.3 * (1.0 - norm)),
                         green: 0.5 + 0.3 * (1.0 - norm),
                         blue: 0.85 + 0.15 * norm,
                         opacity: 0.90)
        }
    }
    
    private func colorForEnergy(_ val: Double, minVal: Double, maxVal: Double) -> Color {
        return Self.evaluateEnergyColormap(
            val: val,
            minVal: minVal,
            maxVal: maxVal,
            preset: manager.energyColormap,
            invert: manager.invertEnergyColormap,
            gamma: manager.energyGamma
        )
    }
    
    public static func evaluateEnergyColormap(
        val: Double,
        minVal: Double,
        maxVal: Double,
        preset: EnergyColormapPreset,
        invert: Bool,
        gamma: Double
    ) -> Color {
        let span = max(1e-9, maxVal - minVal)
        let rawU = max(0.0, min(1.0, (val - minVal) / span))
        let gammaSafe = max(0.1, min(5.0, gamma))
        var u = pow(rawU, gammaSafe)
        if invert {
            u = 1.0 - u
        }
        
        let stops: [(Double, Double, Double, Double)]
        switch preset {
        case .inferno:
            stops = [
                (0.00, 0.00, 0.00, 0.05),
                (0.25, 0.32, 0.07, 0.48),
                (0.50, 0.73, 0.22, 0.33),
                (0.75, 0.97, 0.58, 0.11),
                (1.00, 0.99, 0.99, 0.65)
            ]
        case .viridis:
            stops = [
                (0.00, 0.267, 0.005, 0.329),
                (0.25, 0.230, 0.322, 0.546),
                (0.50, 0.128, 0.567, 0.551),
                (0.75, 0.369, 0.788, 0.383),
                (1.00, 0.993, 0.906, 0.144)
            ]
        case .plasma:
            stops = [
                (0.00, 0.050, 0.030, 0.528),
                (0.25, 0.494, 0.012, 0.658),
                (0.50, 0.798, 0.280, 0.470),
                (0.75, 0.973, 0.580, 0.254),
                (1.00, 0.940, 0.975, 0.131)
            ]
        case .turbo:
            stops = [
                (0.00, 0.190, 0.072, 0.232),
                (0.20, 0.157, 0.565, 0.933),
                (0.40, 0.180, 0.925, 0.600),
                (0.60, 0.733, 0.949, 0.204),
                (0.80, 0.984, 0.502, 0.133),
                (1.00, 0.478, 0.016, 0.012)
            ]
        case .spectral:
            stops = [
                (0.00, 0.196, 0.533, 0.741),
                (0.25, 0.400, 0.761, 0.647),
                (0.50, 0.902, 0.961, 0.596),
                (0.75, 0.992, 0.682, 0.380),
                (1.00, 0.835, 0.243, 0.310)
            ]
        case .coolwarm:
            stops = [
                (0.00, 0.230, 0.299, 0.754),
                (0.25, 0.553, 0.690, 0.996),
                (0.50, 0.865, 0.865, 0.865),
                (0.75, 0.957, 0.604, 0.482),
                (1.00, 0.706, 0.016, 0.150)
            ]
        }
        
        for i in 0 ..< (stops.count - 1) {
            let s0 = stops[i]
            let s1 = stops[i + 1]
            if u >= s0.0 && u <= s1.0 {
                let f = (u - s0.0) / max(1e-9, s1.0 - s0.0)
                let r = s0.1 + f * (s1.1 - s0.1)
                let g = s0.2 + f * (s1.2 - s0.2)
                let b = s0.3 + f * (s1.3 - s0.3)
                return Color(red: r, green: g, blue: b, opacity: 0.92)
            }
        }
        let last = stops.last!
        return Color(red: last.1, green: last.2, blue: last.3, opacity: 0.92)
    }
}

// MARK: - Colorbar Legend
struct ColorbarLegendView: View {
    let colorField: ColorFieldMode
    let minDet: Double
    let maxDet: Double
    let minEnergy: Double
    let maxEnergy: Double
    var colormap: EnergyColormapPreset = .inferno
    var invertColormap: Bool = false
    var gamma: Double = 1.0
    var upperPercentile: Double? = nil
    
    var body: some View {
        VStack(alignment: .leading, spacing: 4) {
            HStack {
                Text(colorField == .acousticStability ? "min det Q(n)" : "Strain Energy Density ΔΦ(C)")
                    .font(.system(size: 10, weight: .bold))
                Spacer()
                if colorField == .acousticStability {
                    Text("Acoustic Stability")
                        .font(.system(size: 9))
                        .foregroundColor(.secondary)
                } else {
                    let pctText = upperPercentile != nil ? String(format: "P%.0f", upperPercentile!) : "Custom"
                    Text("\(colormap.rawValue) (\(pctText), γ=\(String(format: "%.2f", gamma)))")
                        .font(.system(size: 9))
                        .foregroundColor(.secondary)
                }
            }
            
            // Gradient Strip
            if colorField == .acousticStability {
                LinearGradient(
                    colors: [
                        Color(red: 0.95, green: 0.1, blue: 0.1),  // Red (unstable)
                        Color(red: 0.95, green: 0.7, blue: 0.2),  // Orange/Gold (near 0)
                        Color(red: 0.2, green: 0.75, blue: 0.85), // Cyan
                        Color(red: 0.05, green: 0.3, blue: 0.95)  // Deep Blue (stable)
                    ],
                    startPoint: .leading,
                    endPoint: .trailing
                )
                .frame(height: 10)
                .cornerRadius(3)
                
                HStack {
                    Text(String(format: "%.1f (Loss)", minDet))
                        .font(.system(size: 8))
                        .foregroundColor(.red)
                    Spacer()
                    Text("0 (Boundary)")
                        .font(.system(size: 8, weight: .semibold))
                        .foregroundColor(.primary)
                    Spacer()
                    Text(String(format: "%.1f (Stable)", maxDet))
                        .font(.system(size: 8))
                        .foregroundColor(.blue)
                }
            } else {
                let sampledColors = (0...16).map { i -> Color in
                    let frac = Double(i) / 16.0
                    let v = minEnergy + frac * (maxEnergy - minEnergy)
                    return PoincareStabilityMapView.evaluateEnergyColormap(
                        val: v,
                        minVal: minEnergy,
                        maxVal: maxEnergy,
                        preset: colormap,
                        invert: invertColormap,
                        gamma: gamma
                    )
                }
                LinearGradient(
                    colors: sampledColors,
                    startPoint: .leading,
                    endPoint: .trailing
                )
                .frame(height: 10)
                .cornerRadius(3)
                
                HStack {
                    Text(String(format: "%.3f (Min)", minEnergy))
                        .font(.system(size: 8))
                        .foregroundColor(.secondary)
                    Spacer()
                    Text(String(format: "%.3f", (minEnergy + maxEnergy) * 0.5))
                        .font(.system(size: 8))
                        .foregroundColor(.secondary)
                    Spacer()
                    Text(String(format: "%.3f (Cutoff)", maxEnergy))
                        .font(.system(size: 8, weight: .semibold))
                        .foregroundColor(.orange)
                }
            }
        }
    }
}

// MARK: - Hover Inspection Inline View (inside footer bar)
struct HoverInfoInlineView: View {
    let point: AcousticPointModel
    let useReduced: Bool
    
    var body: some View {
        HStack(spacing: 16) {
            // Status Tag
            Text(point.is_stable ? "STABLE" : "UNSTABLE")
                .font(.system(size: 10, weight: .bold))
                .padding(.horizontal, 8)
                .padding(.vertical, 3)
                .background(point.is_stable ? Color.blue.opacity(0.2) : Color.red.opacity(0.2))
                .foregroundColor(point.is_stable ? .blue : .red)
                .cornerRadius(4)
            
            VStack(alignment: .leading, spacing: 2) {
                HStack(spacing: 8) {
                    Text("min det Q:").font(.caption2).foregroundColor(.secondary)
                    Text(String(format: "%.3f", point.min_detAc))
                        .font(.caption.monospaced().bold())
                        .foregroundColor(point.is_stable ? .primary : .red)
                    
                    Text("• Energy ΔΦ:").font(.caption2).foregroundColor(.secondary)
                    Text(String(format: "%.4f", point.energy))
                        .font(.caption.monospaced().bold())
                        .foregroundColor(.orange)
                    
                    Text("• Angle:").font(.caption2).foregroundColor(.secondary)
                    Text(String(format: "%.1f°", point.angle_deg)).font(.caption.monospaced())
                }
                
                HStack(spacing: 8) {
                    Text("Poincaré:").font(.caption2).foregroundColor(.secondary)
                    Text(String(format: "(%.3f, %.3f)", useReduced ? point.x_red : point.x, useReduced ? point.y_red : point.y))
                        .font(.caption.monospaced())
                    
                    Text("• C:").font(.caption2).foregroundColor(.secondary)
                    Text(String(format: "[%.2f, %.2f; %.2f]",
                                useReduced ? point.c11_red : point.c11,
                                useReduced ? point.c22_red : point.c22,
                                useReduced ? point.c12_red : point.c12))
                        .font(.caption.monospaced())
                }
            }
            
            Spacer()
        }
    }
}

// MARK: - Simple Shear Critical Strain Overlay

struct ShearCriticalOverlay: View {
    let path: ShearPathModel
    let lattice: String
    
    var body: some View {
        VStack(alignment: .leading, spacing: 4) {
            HStack(spacing: 6) {
                Image(systemName: "arrow.left.and.right.righttriangle.left.righttriangle.right")
                    .foregroundColor(.green)
                Text("Simple shear  F = [1 α; 0 1]·F₀")
                    .font(.caption.bold())
            }
            Text("Reference: \(lattice) ground state  ·  α ∈ [\(String(format: "%.2f", path.alpha_min)), \(String(format: "%.2f", path.alpha_max))]")
                .font(.caption2)
                .foregroundColor(.secondary)
            Divider()
            if path.critical_alphas.isEmpty {
                Text("No sign change of min det Q(n) in range")
                    .font(.caption2)
                    .foregroundColor(.orange)
            } else {
                ForEach(Array(path.critical_alphas.enumerated()), id: \.offset) { idx, a in
                    HStack(spacing: 6) {
                        Circle().fill(Color.yellow).frame(width: 8, height: 8)
                        Text(String(format: "α_c%d = %+.5f", idx + 1, a))
                            .font(.system(size: 12, weight: .semibold, design: .monospaced))
                    }
                }
            }
            HStack(spacing: 10) {
                HStack(spacing: 4) {
                    RoundedRectangle(cornerRadius: 1).fill(Color(red: 0.2, green: 0.95, blue: 0.4)).frame(width: 14, height: 3)
                    Text("stable").font(.caption2)
                }
                HStack(spacing: 4) {
                    RoundedRectangle(cornerRadius: 1).fill(Color(red: 1.0, green: 0.2, blue: 0.85)).frame(width: 14, height: 3)
                    Text("unstable").font(.caption2)
                }
            }
            .foregroundColor(.secondary)
        }
        .padding(10)
        .background(.ultraThinMaterial, in: RoundedRectangle(cornerRadius: 10))
        .overlay(RoundedRectangle(cornerRadius: 10).stroke(Color.white.opacity(0.15)))
    }
}
