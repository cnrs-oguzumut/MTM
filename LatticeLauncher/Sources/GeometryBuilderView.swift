import SwiftUI
import AppKit

public struct GeometryBuilderView: View {
    @ObservedObject var config: SimulationConfig
    @ObservedObject var previewManager: StiffnessStudyManager
    var onSwitchToCarvedMesh: () -> Void
    
    @State private var selectedCutoutID: String? = nil
    @State private var selectedRuleID: String? = nil
    @State private var showDiscreteLatticePreview: Bool = true
    @State private var showBCArrows: Bool = true
    @State private var activeDragMode: DragHandleMode = .none
    @State private var dragStartCutout: GeometryCutoutItem? = nil
    @State private var inspectorSection: InspectorSection = .microstructure
    
    private enum InspectorSection: String, CaseIterable, Identifiable {
        case microstructure = "Microstructure & GB"
        case shapeAndCutouts = "Shape & Cutouts"
        case bcsAndIndenter = "BCs & Indenter"
        case all = "All"
        var id: String { rawValue }
    }
    
    private enum DragHandleMode {
        case none
        case moveCenter(String)
        case resizeRadius(String)
        case crackStart(String)
        case crackEnd(String)
        case polygonVertex(String, Int)
        case movePolygon(String, Double, Double)
        case outerPolygonVertex(Int)
        case moveBicrystalGB
        case moveVoronoiSeed(Int)
    }
    
    public init(
        config: SimulationConfig,
        previewManager: StiffnessStudyManager,
        onSwitchToCarvedMesh: @escaping () -> Void
    ) {
        self.config = config
        self.previewManager = previewManager
        self.onSwitchToCarvedMesh = onSwitchToCarvedMesh
    }
    
    public var body: some View {
        VStack(spacing: 0) {
            studioTopBar
            Divider()
            HStack(spacing: 0) {
                // LEFT: Interactive Macroscopic Geometry & Discrete Lattice Canvas
                VStack(spacing: 0) {
                    canvasHeaderHintBar
                    Divider()
                    GeometryReader { geo in
                        interactiveCanvas(size: geo.size)
                    }
                    Divider()
                    canvasFooterLegend
                }
                .frame(minWidth: 380, maxWidth: .infinity, maxHeight: .infinity)
                .background(BlueHour.canvasBackground)
                
                Divider()
                
                // RIGHT: Organized Inspector Panel with Section Switcher (Fixed Width)
                VStack(spacing: 0) {
                    VStack(alignment: .leading, spacing: 6) {
                        HStack {
                            Label("Geometry & GB Inspector", systemImage: "slider.horizontal.3")
                                .font(.caption.bold())
                                .foregroundColor(BlueHour.ice)
                            Spacer()
                            Text("\(config.nx)×\(config.ny)")
                                .font(.system(size: 10, weight: .bold, design: .monospaced))
                                .foregroundColor(Color(red: 0.31, green: 0.85, blue: 0.78))
                                .padding(.horizontal, 6)
                                .padding(.vertical, 2)
                                .background(BlueHour.midnight.opacity(0.65), in: Capsule())
                        }
                        Picker("", selection: $inspectorSection) {
                            Text("Grains & GB").tag(InspectorSection.microstructure)
                            Text("Shape & Holes").tag(InspectorSection.shapeAndCutouts)
                            Text("BCs & Indenter").tag(InspectorSection.bcsAndIndenter)
                            Text("All").tag(InspectorSection.all)
                        }
                        .labelsHidden()
                        .pickerStyle(.segmented)
                        .controlSize(.small)
                    }
                    .padding(.horizontal, 12)
                    .padding(.vertical, 8)
                    .background(BlueHour.controlSurface)
                    
                    Divider()
                    
                    ScrollView {
                        VStack(alignment: .leading, spacing: 12) {
                            if inspectorSection == .microstructure || inspectorSection == .all {
                                microstructureCard
                                if config.geometryConfig.microstructure.mode == .singleCrystal {
                                    crystalLatticeCard
                                }
                            }
                            if inspectorSection == .shapeAndCutouts || inspectorSection == .all {
                                presetsCard
                                outerDomainCard
                                cutoutsCard
                            }
                            if inspectorSection == .bcsAndIndenter || inspectorSection == .all {
                                boundaryConditionsCard
                                nanoindenterCard
                            }
                        }
                        .frame(maxWidth: .infinity, alignment: .leading)
                        .padding(12)
                    }
                }
                .frame(width: 400)
                .background(BlueHour.sidebarBackground)
            }
        }
        .onAppear {
            if selectedCutoutID == nil {
                selectedCutoutID = config.geometryConfig.cutouts.first?.id
            }
        }
    }
    
    private var canvasHeaderHintBar: some View {
        HStack(spacing: 8) {
            Image(systemName: "hand.draw.fill")
                .font(.system(size: 10))
                .foregroundColor(BlueHour.ice)
            Text("Interactive Macroscopic Geometry Canvas:")
                .font(.system(size: 10.5, weight: .bold))
                .foregroundColor(BlueHour.ice)
            Text("Drag GB/seed/cutout handles directly on canvas • Click 'Generate & View Carved Triangulation' to mesh")
                .font(.system(size: 10))
                .foregroundColor(BlueHour.mist)
                .lineLimit(1)
            Spacer()
        }
        .padding(.horizontal, 12)
        .padding(.vertical, 5)
        .background(BlueHour.midnight.opacity(0.65))
    }
    
    // MARK: - Top Toolbar
    private var studioTopBar: some View {
        ScrollView(.horizontal, showsIndicators: false) {
            HStack(spacing: 12) {
                Toggle(isOn: $config.geometryConfig.enabled) {
                    Label("Enable Custom Macroscopic Geometry & BCs", systemImage: "DirectSelection")
                        .font(.caption.bold())
                }
                .toggleStyle(.switch)
                .controlSize(.small)
                
                Divider().frame(height: 18)
                
                // Discrete Crystal Symmetry Picker (Synced with SimulationConfig)
                HStack(spacing: 6) {
                    Text("Symmetry:")
                        .font(.caption2.bold())
                        .foregroundColor(.secondary)
                    Picker("", selection: $config.latticeSymmetry) {
                        ForEach(LatticeSymmetry.allCases) { sym in
                            Text(sym.displayName).tag(sym)
                        }
                    }
                    .pickerStyle(.segmented)
                    .controlSize(.small)
                    .fixedSize()
                }
                
                // Rotated Crystal Controls in Studio Top Bar
                HStack(spacing: 6) {
                    Toggle(isOn: $config.useRotation) {
                        Text("Rotate θ:")
                            .font(.caption2.bold())
                    }
                    .toggleStyle(.checkbox)
                    
                    if config.useRotation {
                        DoubleField(title: "0.0", value: $config.rotationAngle, width: 44, minValue: 0.0, maxValue: 360.0)
                        Text("°").font(.caption2)
                        Slider(value: $config.rotationAngle, in: 0...360, step: 1.0)
                            .frame(width: 80)
                    }
                }
                
                Divider().frame(height: 18)
                
                Toggle("Show Discrete Lattice Sites", isOn: $showDiscreteLatticePreview)
                    .font(.caption)
                    .toggleStyle(.checkbox)
                
                Toggle("Show BC Vectors", isOn: $showBCArrows)
                    .font(.caption)
                    .toggleStyle(.checkbox)
                
                Divider().frame(height: 18)
                
                Button(action: {
                    config.geometryConfig.enabled = true
                    triggerCarvedMeshPreview()
                    onSwitchToCarvedMesh()
                }) {
                    if previewManager.isGeneratingPreview {
                        ProgressView().scaleEffect(0.6).frame(width: 14, height: 14)
                    } else {
                        Label("Generate & View Carved Triangulation", systemImage: "triangle.circle.fill")
                            .font(.caption.bold())
                    }
                }
                .buttonStyle(.borderedProminent)
                .controlSize(.small)
                
                // Discrete Lattice Fixed-Scale Badge
                HStack(spacing: 6) {
                    Image(systemName: "atom")
                        .font(.system(size: 10))
                    Text("\(config.nx)×\(config.ny) \(config.latticeSymmetry == .square ? "Square" : "Triangular") | Fixed Unit Cell a₀")
                        .font(.system(size: 10, weight: .semibold, design: .monospaced))
                }
                .padding(.horizontal, 8)
                .padding(.vertical, 4)
                .background(BlueHour.pillSurface, in: Capsule())
                .overlay(Capsule().strokeBorder(BlueHour.border, lineWidth: 1))
            }
            .padding(.horizontal, 12)
            .padding(.vertical, 7)
        }
        .background(BlueHour.controlSurface)
    }
    
    // MARK: - Interactive Canvas
    private func interactiveCanvas(size: CGSize) -> some View {
        let isIndenterActive = config.geometryConfig.indenter.enabled
        let padSide: CGFloat = 44.0
        let padTop: CGFloat = isIndenterActive ? 75.0 : 44.0
        let padBottom: CGFloat = 44.0
        let availW = max(100.0, size.width - 2 * padSide)
        let availH = max(100.0, size.height - (padTop + padBottom))
        let aspect = (config.latticeSymmetry == .triangular)
            ? (Double(config.nx) / (Double(config.ny) * sqrt(3.0) / 2.0))
            : (Double(config.nx) / Double(max(1, config.ny)))
        let boxW: CGFloat
        let boxH: CGFloat
        if availW / availH > CGFloat(aspect) {
            boxH = availH
            boxW = availH * CGFloat(aspect)
        } else {
            boxW = availW
            boxH = availW / CGFloat(max(0.2, aspect))
        }
        let originX = (size.width - boxW) * 0.5
        let originY = padTop + (availH - boxH) * 0.5
        
        func toScreen(_ xHat: Double, _ yHat: Double) -> CGPoint {
            CGPoint(
                x: originX + CGFloat(xHat) * boxW,
                y: originY + CGFloat(1.0 - yHat) * boxH
            )
        }
        
        func toNormalized(_ pt: CGPoint) -> (Double, Double) {
            let xh = Double((pt.x - originX) / max(1.0, boxW))
            let yh = 1.0 - Double((pt.y - originY) / max(1.0, boxH))
            return (min(1.0, max(0.0, xh)), min(1.0, max(0.0, yh)))
        }
        
        return ZStack {
            Canvas { ctx, canvasSize in
                let boxRect = CGRect(x: originX, y: originY, width: boxW, height: boxH)
                
                // 1. Reference bounding box outline (dashed)
                var boxPath = Path()
                boxPath.addRect(boxRect)
                ctx.stroke(
                    boxPath,
                    with: .color(BlueHour.steel.opacity(0.45)),
                    style: StrokeStyle(lineWidth: 1.0, dash: [5, 4])
                )
                
                // 2. Outer Macroscopic Domain Path
                let outerPath = buildOuterDomainPath(originX: originX, originY: originY, boxW: boxW, boxH: boxH)
                ctx.fill(outerPath, with: .color(BlueHour.midnight.opacity(0.85)))
                ctx.stroke(outerPath, with: .color(BlueHour.mist.opacity(0.9)), lineWidth: 2.0)
                
                // 3. Cutout / Hole Shapes (Filled with Abyss background & highlighted rim)
                for cutout in config.geometryConfig.cutouts where cutout.enabled {
                    let isSelected = (cutout.id == selectedCutoutID)
                    let cPath = buildCutoutPath(cutout: cutout, originX: originX, originY: originY, boxW: boxW, boxH: boxH)
                    ctx.fill(cPath, with: .color(BlueHour.abyss))
                    ctx.stroke(
                        cPath,
                        with: .color(isSelected ? Color(red: 0.31, green: 0.85, blue: 0.78) : BlueHour.ice.opacity(0.75)),
                        style: StrokeStyle(lineWidth: isSelected ? 2.4 : 1.5, dash: cutout.type == .preCrack ? [4, 2] : [])
                    )
                }
                
                // 4. Discrete Crystal Lattice Sites (Square or Triangular, Color-Coded by Grain & BC Role)
                if showDiscreteLatticePreview {
                    drawDiscreteLatticeSitesPreview(ctx: &ctx, toScreen: toScreen, boxW: boxW, boxH: boxH)
                }
                
                // 4b. Bicrystal GB Plane / Voronoi Grain Cell Boundaries & Seed Handles
                if config.geometryConfig.microstructure.mode != .singleCrystal {
                    drawMicrostructureOverlay(ctx: &ctx, toScreen: toScreen, originX: originX, originY: originY, boxW: boxW, boxH: boxH)
                }
                
                // 5. Boundary Condition Velocity Arrows on Driven Edges
                if showBCArrows {
                    drawBoundaryArrows(ctx: &ctx, originX: originX, originY: originY, boxW: boxW, boxH: boxH)
                }
                
                // 6. Interactive Drag Handles for Selected Cutout
                if let selID = selectedCutoutID,
                   let cutout = config.geometryConfig.cutouts.first(where: { $0.id == selID && $0.enabled }) {
                    drawCutoutHandles(ctx: &ctx, cutout: cutout, toScreen: toScreen, boxW: boxW, boxH: boxH)
                }
                
                // 7. Outer Polygon Vertex Handles
                if config.geometryConfig.outer_shape == .polygon {
                    for vert in config.geometryConfig.outer_polygon {
                        let p = toScreen(vert.x, vert.y)
                        let hRect = CGRect(x: p.x - 5, y: p.y - 5, width: 10, height: 10)
                        ctx.fill(Path(ellipseIn: hRect), with: .color(BlueHour.ice))
                        ctx.stroke(Path(ellipseIn: hRect), with: .color(BlueHour.abyss), lineWidth: 1.5)
                    }
                }
                
                // 8. Nanoindenter Tip Visual Preview
                if config.geometryConfig.indenter.enabled {
                    drawNanoindenterTip(ctx: &ctx, toScreen: toScreen, boxW: boxW, boxH: boxH)
                }
            }
            .gesture(
                DragGesture(minimumDistance: 1)
                    .onChanged { value in
                        handleCanvasDragChanged(
                            startLocation: value.startLocation,
                            currentLocation: value.location,
                            toScreen: toScreen,
                            toNormalized: toNormalized,
                            boxW: boxW,
                            boxH: boxH
                        )
                    }
                    .onEnded { _ in
                        activeDragMode = .none
                        dragStartCutout = nil
                    }
            )
        }
    }
    
    // MARK: - Canvas Path & Handle Helpers
    private func buildOuterDomainPath(originX: CGFloat, originY: CGFloat, boxW: CGFloat, boxH: CGFloat) -> Path {
        var path = Path()
        switch config.geometryConfig.outer_shape {
        case .rectangle:
            path.addRect(CGRect(x: originX, y: originY, width: boxW, height: boxH))
        case .disk:
            let r = CGFloat(config.geometryConfig.disk_radius) * min(boxW, boxH)
            let cx = originX + 0.5 * boxW
            let cy = originY + 0.5 * boxH
            path.addEllipse(in: CGRect(x: cx - r, y: cy - r, width: 2 * r, height: 2 * r))
        case .dogbone:
            let steps = 120
            var pts: [CGPoint] = []
            let halfGW = 0.5 * min(0.98, max(0.15, config.geometryConfig.dogbone_gauge_width))
            let halfGL = 0.5 * min(0.90, max(0.10, config.geometryConfig.dogbone_gauge_length))
            let rf = max(0.02, config.geometryConfig.dogbone_fillet_radius)
            func halfWidthAt(_ yHat: Double) -> Double {
                let dy = abs(yHat - 0.5)
                if dy <= halfGL { return halfGW }
                if dy <= halfGL + rf {
                    let u = dy - halfGL
                    let filletDx = rf - sqrt(max(0.0, rf * rf - u * u))
                    return min(0.5, halfGW + filletDx * ((0.5 - halfGW) / rf))
                }
                return 0.5
            }
            for k in 0...steps {
                let yh = Double(k) / Double(steps)
                let xh = 0.5 + halfWidthAt(yh)
                pts.append(CGPoint(x: originX + CGFloat(xh) * boxW, y: originY + CGFloat(1.0 - yh) * boxH))
            }
            for k in (0...steps).reversed() {
                let yh = Double(k) / Double(steps)
                let xh = 0.5 - halfWidthAt(yh)
                pts.append(CGPoint(x: originX + CGFloat(xh) * boxW, y: originY + CGFloat(1.0 - yh) * boxH))
            }
            if let first = pts.first {
                path.move(to: first)
                for p in pts.dropFirst() { path.addLine(to: p) }
                path.closeSubpath()
            }
        case .polygon:
            let verts = config.geometryConfig.outer_polygon
            if let first = verts.first {
                path.move(to: CGPoint(x: originX + CGFloat(first.x) * boxW, y: originY + CGFloat(1.0 - first.y) * boxH))
                for v in verts.dropFirst() {
                    path.addLine(to: CGPoint(x: originX + CGFloat(v.x) * boxW, y: originY + CGFloat(1.0 - v.y) * boxH))
                }
                path.closeSubpath()
            }
        }
        return path
    }
    
    private func buildCutoutPath(cutout: GeometryCutoutItem, originX: CGFloat, originY: CGFloat, boxW: CGFloat, boxH: CGFloat) -> Path {
        var path = Path()
        let cx = originX + CGFloat(cutout.cx) * boxW
        let cy = originY + CGFloat(1.0 - cutout.cy) * boxH
        let minDim = min(boxW, boxH)
        
        switch cutout.type {
        case .circle:
            let r = CGFloat(cutout.rx) * minDim
            path.addEllipse(in: CGRect(x: cx - r, y: cy - r, width: 2 * r, height: 2 * r))
        case .ellipse:
            let rx = CGFloat(cutout.rx) * boxW
            let ry = CGFloat(cutout.ry) * boxH
            var local = Path(ellipseIn: CGRect(x: -rx, y: -ry, width: 2 * rx, height: 2 * ry))
            let tr = CGAffineTransform(translationX: cx, y: cy)
                .rotated(by: CGFloat(-cutout.rotation_deg * .pi / 180.0))
            path.addPath(local, transform: tr)
        case .rectangle:
            let rx = CGFloat(cutout.rx) * boxW
            let ry = CGFloat(cutout.ry) * boxH
            let local = Path(CGRect(x: -rx, y: -ry, width: 2 * rx, height: 2 * ry))
            let tr = CGAffineTransform(translationX: cx, y: cy)
                .rotated(by: CGFloat(-cutout.rotation_deg * .pi / 180.0))
            path.addPath(local, transform: tr)
        case .edgeNotch:
            let d = cutout.rx
            let w = cutout.ry
            if cutout.notch_edge == "left" {
                let p1 = CGPoint(x: originX, y: originY + CGFloat(1.0 - (cutout.cy - w)) * boxH)
                let pTip = CGPoint(x: originX + CGFloat(d) * boxW, y: originY + CGFloat(1.0 - cutout.cy) * boxH)
                let p2 = CGPoint(x: originX, y: originY + CGFloat(1.0 - (cutout.cy + w)) * boxH)
                path.move(to: p1); path.addLine(to: pTip); path.addLine(to: p2); path.closeSubpath()
            } else if cutout.notch_edge == "right" {
                let p1 = CGPoint(x: originX + boxW, y: originY + CGFloat(1.0 - (cutout.cy - w)) * boxH)
                let pTip = CGPoint(x: originX + CGFloat(1.0 - d) * boxW, y: originY + CGFloat(1.0 - cutout.cy) * boxH)
                let p2 = CGPoint(x: originX + boxW, y: originY + CGFloat(1.0 - (cutout.cy + w)) * boxH)
                path.move(to: p1); path.addLine(to: pTip); path.addLine(to: p2); path.closeSubpath()
            } else if cutout.notch_edge == "bottom" {
                let p1 = CGPoint(x: originX + CGFloat(cutout.cx - w) * boxW, y: originY + boxH)
                let pTip = CGPoint(x: originX + CGFloat(cutout.cx) * boxW, y: originY + CGFloat(1.0 - d) * boxH)
                let p2 = CGPoint(x: originX + CGFloat(cutout.cx + w) * boxW, y: originY + boxH)
                path.move(to: p1); path.addLine(to: pTip); path.addLine(to: p2); path.closeSubpath()
            } else {
                let p1 = CGPoint(x: originX + CGFloat(cutout.cx - w) * boxW, y: originY)
                let pTip = CGPoint(x: originX + CGFloat(cutout.cx) * boxW, y: originY + CGFloat(d) * boxH)
                let p2 = CGPoint(x: originX + CGFloat(cutout.cx + w) * boxW, y: originY)
                path.move(to: p1); path.addLine(to: pTip); path.addLine(to: p2); path.closeSubpath()
            }
        case .preCrack:
            let p1 = CGPoint(x: originX + CGFloat(cutout.x1) * boxW, y: originY + CGFloat(1.0 - cutout.y1) * boxH)
            let p2 = CGPoint(x: originX + CGFloat(cutout.x2) * boxW, y: originY + CGFloat(1.0 - cutout.y2) * boxH)
            path.move(to: p1)
            path.addLine(to: p2)
        case .polygon:
            if let first = cutout.vertices.first {
                path.move(to: CGPoint(x: originX + CGFloat(first.x) * boxW, y: originY + CGFloat(1.0 - first.y) * boxH))
                for v in cutout.vertices.dropFirst() {
                    path.addLine(to: CGPoint(x: originX + CGFloat(v.x) * boxW, y: originY + CGFloat(1.0 - v.y) * boxH))
                }
                path.closeSubpath()
            }
        }
        return path
    }
    
    private func drawCutoutHandles(
        ctx: inout GraphicsContext,
        cutout: GeometryCutoutItem,
        toScreen: (Double, Double) -> CGPoint,
        boxW: CGFloat,
        boxH: CGFloat
    ) {
        let cyan = Color(red: 0.31, green: 0.85, blue: 0.78)
        let gold = Color(red: 1.0, green: 0.82, blue: 0.40)
        
        switch cutout.type {
        case .circle, .ellipse, .rectangle, .edgeNotch:
            let centerPt: CGPoint
            if cutout.type == .edgeNotch {
                if cutout.notch_edge == "left" { centerPt = toScreen(0.0, cutout.cy) }
                else if cutout.notch_edge == "right" { centerPt = toScreen(1.0, cutout.cy) }
                else if cutout.notch_edge == "bottom" { centerPt = toScreen(cutout.cx, 0.0) }
                else { centerPt = toScreen(cutout.cx, 1.0) }
            } else {
                centerPt = toScreen(cutout.cx, cutout.cy)
            }
            let cRect = CGRect(x: centerPt.x - 6, y: centerPt.y - 6, width: 12, height: 12)
            ctx.fill(Path(ellipseIn: cRect), with: .color(gold))
            ctx.stroke(Path(ellipseIn: cRect), with: .color(BlueHour.abyss), lineWidth: 1.5)
            
            // Resize handle
            let rimPt: CGPoint
            if cutout.type == .circle {
                let rPx = CGFloat(cutout.rx) * min(boxW, boxH)
                rimPt = CGPoint(x: centerPt.x + rPx, y: centerPt.y)
            } else if cutout.type == .edgeNotch {
                if cutout.notch_edge == "left" { rimPt = toScreen(cutout.rx, cutout.cy) }
                else if cutout.notch_edge == "right" { rimPt = toScreen(1.0 - cutout.rx, cutout.cy) }
                else if cutout.notch_edge == "bottom" { rimPt = toScreen(cutout.cx, cutout.rx) }
                else { rimPt = toScreen(cutout.cx, 1.0 - cutout.rx) }
            } else {
                rimPt = toScreen(min(1.0, cutout.cx + cutout.rx), min(1.0, cutout.cy + cutout.ry))
            }
            let rRect = CGRect(x: rimPt.x - 5.5, y: rimPt.y - 5.5, width: 11, height: 11)
            ctx.fill(Path(ellipseIn: rRect), with: .color(cyan))
            ctx.stroke(Path(ellipseIn: rRect), with: .color(BlueHour.abyss), lineWidth: 1.5)
            
        case .preCrack:
            let p1 = toScreen(cutout.x1, cutout.y1)
            let p2 = toScreen(cutout.x2, cutout.y2)
            for p in [p1, p2] {
                let r = CGRect(x: p.x - 6, y: p.y - 6, width: 12, height: 12)
                ctx.fill(Path(ellipseIn: r), with: .color(cyan))
                ctx.stroke(Path(ellipseIn: r), with: .color(BlueHour.abyss), lineWidth: 1.5)
            }
        case .polygon:
            if !cutout.vertices.isEmpty {
                let cx = cutout.vertices.map(\.x).reduce(0.0, +) / Double(cutout.vertices.count)
                let cy = cutout.vertices.map(\.y).reduce(0.0, +) / Double(cutout.vertices.count)
                let centerPt = toScreen(cx, cy)
                let cRect = CGRect(x: centerPt.x - 6, y: centerPt.y - 6, width: 12, height: 12)
                ctx.fill(Path(ellipseIn: cRect), with: .color(gold))
                ctx.stroke(Path(ellipseIn: cRect), with: .color(BlueHour.abyss), lineWidth: 1.5)
            }
            for v in cutout.vertices {
                let p = toScreen(v.x, v.y)
                let r = CGRect(x: p.x - 5.5, y: p.y - 5.5, width: 11, height: 11)
                ctx.fill(Path(ellipseIn: r), with: .color(cyan))
                ctx.stroke(Path(ellipseIn: r), with: .color(BlueHour.abyss), lineWidth: 1.5)
            }
        }
    }
    
    private func drawBoundaryArrows(
        ctx: inout GraphicsContext,
        originX: CGFloat,
        originY: CGFloat,
        boxW: CGFloat,
        boxH: CGFloat
    ) {
        let drivenColor = Color(red: 0.31, green: 0.85, blue: 0.78)
        let fixedColor = Color(red: 1.0, green: 0.42, blue: 0.42)
        let rollerColor = Color(red: 0.98, green: 0.72, blue: 0.24)
        
        for rule in config.geometryConfig.boundary_rules where rule.enabled {
            if rule.region == .customBox {
                let bx = originX + CGFloat(rule.box_x_min) * boxW
                let by = originY + CGFloat(1.0 - rule.box_y_max) * boxH
                let bw = CGFloat(max(0.01, rule.box_x_max - rule.box_x_min)) * boxW
                let bh = CGFloat(max(0.01, rule.box_y_max - rule.box_y_min)) * boxH
                let rect = CGRect(x: bx, y: by, width: bw, height: bh)
                let col: Color
                switch rule.kind {
                case .fixed: col = fixedColor
                case .rollerX, .rollerY: col = rollerColor
                case .free: col = Color.gray
                case .prescribed, .affine: col = drivenColor
                }
                ctx.fill(Path(rect), with: .color(col.opacity(0.12)))
                ctx.stroke(Path(rect), with: .color(col.opacity(0.85)), style: StrokeStyle(lineWidth: 1.5, dash: [4, 3]))
            }

            if rule.kind == .fixed || rule.kind == .rollerX || rule.kind == .rollerY {
                let isRoller = (rule.kind == .rollerX || rule.kind == .rollerY || !rule.lock_x || !rule.lock_y)
                let markerColor = isRoller ? rollerColor : fixedColor
                let count = 7
                for k in 1..<count {
                    let t = CGFloat(k) / CGFloat(count)
                    var pt = CGPoint.zero
                    if rule.region == .bottom { pt = CGPoint(x: originX + t * boxW, y: originY + boxH + 8) }
                    else if rule.region == .top { pt = CGPoint(x: originX + t * boxW, y: originY - 8) }
                    else if rule.region == .left { pt = CGPoint(x: originX - 8, y: originY + t * boxH) }
                    else if rule.region == .right { pt = CGPoint(x: originX + boxW + 8, y: originY + t * boxH) }
                    else { continue }
                    let tri = Path(ellipseIn: CGRect(x: pt.x - 3.5, y: pt.y - 3.5, width: 7, height: 7))
                    ctx.stroke(tri, with: .color(markerColor), lineWidth: 1.6)
                    if isRoller {
                        var slideLine = Path()
                        if rule.kind == .rollerX || (!rule.lock_x && rule.lock_y) {
                            slideLine.move(to: CGPoint(x: pt.x - 5.5, y: pt.y))
                            slideLine.addLine(to: CGPoint(x: pt.x + 5.5, y: pt.y))
                        } else {
                            slideLine.move(to: CGPoint(x: pt.x, y: pt.y - 5.5))
                            slideLine.addLine(to: CGPoint(x: pt.x, y: pt.y + 5.5))
                        }
                        ctx.stroke(slideLine, with: .color(markerColor), lineWidth: 1.4)
                    }
                }
            } else if rule.kind == .prescribed || rule.kind == .affine {
                var vx = rule.ux
                var vy = rule.uy
                if rule.kind == .affine {
                    if rule.affine_mode == "uniaxial_y" { vx = 0.0; vy = 1.0 }
                    else if rule.affine_mode == "uniaxial_x" { vx = 1.0; vy = 0.0 }
                    else { vx = 1.0; vy = 0.0 }
                }
                let mag = hypot(vx, vy)
                guard mag > 1e-4 else { continue }
                let dx = CGFloat(vx / mag) * 18.0
                let dy = CGFloat(-vy / mag) * 18.0
                
                let count = 5
                for k in 1...count {
                    let t = CGFloat(k) / CGFloat(count + 1)
                    var start = CGPoint.zero
                    if rule.region == .top || rule.region == .outerPerimeter {
                        start = CGPoint(x: originX + t * boxW, y: originY - 6)
                    } else if rule.region == .bottom {
                        start = CGPoint(x: originX + t * boxW, y: originY + boxH + 6)
                    } else if rule.region == .right {
                        start = CGPoint(x: originX + boxW + 6, y: originY + t * boxH)
                    } else if rule.region == .left {
                        start = CGPoint(x: originX - 6, y: originY + t * boxH)
                    } else {
                        continue
                    }
                    let end = CGPoint(x: start.x + dx, y: start.y + dy)
                    var arrow = Path()
                    arrow.move(to: start)
                    arrow.addLine(to: end)
                    ctx.stroke(arrow, with: .color(drivenColor), lineWidth: 2.0)
                    let head = CGRect(x: end.x - 3, y: end.y - 3, width: 6, height: 6)
                    ctx.fill(Path(ellipseIn: head), with: .color(drivenColor))
                }
            }
        }
    }
    
    // MARK: - Drag Gesture Logic
    private func handleCanvasDragChanged(
        startLocation: CGPoint,
        currentLocation: CGPoint,
        toScreen: (Double, Double) -> CGPoint,
        toNormalized: (CGPoint) -> (Double, Double),
        boxW: CGFloat,
        boxH: CGFloat
    ) {
        if case .none = activeDragMode {
            // Hit-test Bicrystal GB center handle or Voronoi grain seed handles first
            let micro = config.geometryConfig.microstructure
            if micro.mode == .bicrystal {
                let gbPt = toScreen(micro.gb_cx, micro.gb_cy)
                if hypot(startLocation.x - gbPt.x, startLocation.y - gbPt.y) < 16 {
                    activeDragMode = .moveBicrystalGB
                }
            } else if micro.mode == .voronoi {
                for (sIdx, s) in micro.seeds.enumerated() {
                    let sPt = toScreen(s.cx, s.cy)
                    if hypot(startLocation.x - sPt.x, startLocation.y - sPt.y) < 16 {
                        activeDragMode = .moveVoronoiSeed(sIdx)
                        break
                    }
                }
            }

            // Hit-test outer polygon vertices
            if case .none = activeDragMode, config.geometryConfig.outer_shape == .polygon {
                for (idx, v) in config.geometryConfig.outer_polygon.enumerated() {
                    let p = toScreen(v.x, v.y)
                    if hypot(startLocation.x - p.x, startLocation.y - p.y) < 14 {
                        activeDragMode = .outerPolygonVertex(idx)
                        break
                    }
                }
            }
            
            // Hit-test cutout handles
            if case .none = activeDragMode {
                for cutout in config.geometryConfig.cutouts where cutout.enabled {
                    switch cutout.type {
                    case .circle, .ellipse, .rectangle, .edgeNotch:
                        let centerPt: CGPoint
                        if cutout.type == .edgeNotch {
                            if cutout.notch_edge == "left" { centerPt = toScreen(0.0, cutout.cy) }
                            else if cutout.notch_edge == "right" { centerPt = toScreen(1.0, cutout.cy) }
                            else if cutout.notch_edge == "bottom" { centerPt = toScreen(cutout.cx, 0.0) }
                            else { centerPt = toScreen(cutout.cx, 1.0) }
                        } else {
                            centerPt = toScreen(cutout.cx, cutout.cy)
                        }
                        let rimPt: CGPoint
                        if cutout.type == .circle {
                            let rPx = CGFloat(cutout.rx) * min(boxW, boxH)
                            rimPt = CGPoint(x: centerPt.x + rPx, y: centerPt.y)
                        } else if cutout.type == .edgeNotch {
                            if cutout.notch_edge == "left" { rimPt = toScreen(cutout.rx, cutout.cy) }
                            else if cutout.notch_edge == "right" { rimPt = toScreen(1.0 - cutout.rx, cutout.cy) }
                            else if cutout.notch_edge == "bottom" { rimPt = toScreen(cutout.cx, cutout.rx) }
                            else { rimPt = toScreen(cutout.cx, 1.0 - cutout.rx) }
                        } else {
                            rimPt = toScreen(min(1.0, cutout.cx + cutout.rx), min(1.0, cutout.cy + cutout.ry))
                        }
                        if hypot(startLocation.x - rimPt.x, startLocation.y - rimPt.y) < 15 {
                            selectedCutoutID = cutout.id
                            dragStartCutout = cutout
                            activeDragMode = .resizeRadius(cutout.id)
                            break
                        }
                        if hypot(startLocation.x - centerPt.x, startLocation.y - centerPt.y) < 22 {
                            selectedCutoutID = cutout.id
                            dragStartCutout = cutout
                            activeDragMode = .moveCenter(cutout.id)
                            break
                        }
                    case .preCrack:
                        let p1 = toScreen(cutout.x1, cutout.y1)
                        let p2 = toScreen(cutout.x2, cutout.y2)
                        if hypot(startLocation.x - p1.x, startLocation.y - p1.y) < 16 {
                            selectedCutoutID = cutout.id
                            activeDragMode = .crackStart(cutout.id)
                            break
                        }
                        if hypot(startLocation.x - p2.x, startLocation.y - p2.y) < 16 {
                            selectedCutoutID = cutout.id
                            activeDragMode = .crackEnd(cutout.id)
                            break
                        }
                    case .polygon:
                        var matchedVertex = false
                        for (vIdx, v) in cutout.vertices.enumerated() {
                            let p = toScreen(v.x, v.y)
                            if hypot(startLocation.x - p.x, startLocation.y - p.y) < 14 {
                                selectedCutoutID = cutout.id
                                activeDragMode = .polygonVertex(cutout.id, vIdx)
                                matchedVertex = true
                                break
                            }
                        }
                        if matchedVertex { break }
                        if !cutout.vertices.isEmpty {
                            let cx = cutout.vertices.map(\.x).reduce(0.0, +) / Double(cutout.vertices.count)
                            let cy = cutout.vertices.map(\.y).reduce(0.0, +) / Double(cutout.vertices.count)
                            let centerPt = toScreen(cx, cy)
                            if hypot(startLocation.x - centerPt.x, startLocation.y - centerPt.y) < 20 {
                                selectedCutoutID = cutout.id
                                dragStartCutout = cutout
                                let (sxh, syh) = toNormalized(startLocation)
                                activeDragMode = .movePolygon(cutout.id, sxh, syh)
                                break
                            }
                        }
                    }
                }
            }
        }
        
        let (xh, yh) = toNormalized(currentLocation)
        switch activeDragMode {
        case .none:
            break
        case .moveBicrystalGB:
            config.geometryConfig.microstructure.gb_cx = (round(xh * 200) / 200)
            config.geometryConfig.microstructure.gb_cy = (round(yh * 200) / 200)
        case .moveVoronoiSeed(let sIdx):
            if sIdx >= 0 && sIdx < config.geometryConfig.microstructure.seeds.count {
                config.geometryConfig.microstructure.seeds[sIdx].cx = (round(xh * 200) / 200)
                config.geometryConfig.microstructure.seeds[sIdx].cy = (round(yh * 200) / 200)
            }
        case .outerPolygonVertex(let idx):
            if idx < config.geometryConfig.outer_polygon.count {
                config.geometryConfig.outer_polygon[idx].x = xh
                config.geometryConfig.outer_polygon[idx].y = yh
            }
        case .moveCenter(let id):
            if let idx = config.geometryConfig.cutouts.firstIndex(where: { $0.id == id }) {
                config.geometryConfig.cutouts[idx].cx = (round(xh * 200) / 200)
                config.geometryConfig.cutouts[idx].cy = (round(yh * 200) / 200)
            }
        case .resizeRadius(let id):
            if let idx = config.geometryConfig.cutouts.firstIndex(where: { $0.id == id }) {
                let c = config.geometryConfig.cutouts[idx]
                if c.type == .circle {
                    let r = max(0.03, min(0.46, hypot(xh - c.cx, yh - c.cy)))
                    config.geometryConfig.cutouts[idx].rx = (round(r * 200) / 200)
                    config.geometryConfig.cutouts[idx].ry = config.geometryConfig.cutouts[idx].rx
                } else if c.type == .edgeNotch {
                    if c.notch_edge == "left" { config.geometryConfig.cutouts[idx].rx = max(0.03, min(0.45, xh)) }
                    else if c.notch_edge == "right" { config.geometryConfig.cutouts[idx].rx = max(0.03, min(0.45, 1.0 - xh)) }
                    else if c.notch_edge == "bottom" { config.geometryConfig.cutouts[idx].rx = max(0.03, min(0.45, yh)) }
                    else { config.geometryConfig.cutouts[idx].rx = max(0.03, min(0.45, 1.0 - yh)) }
                } else {
                    config.geometryConfig.cutouts[idx].rx = max(0.03, min(0.46, abs(xh - c.cx)))
                    config.geometryConfig.cutouts[idx].ry = max(0.03, min(0.46, abs(yh - c.cy)))
                }
            }
        case .crackStart(let id):
            if let idx = config.geometryConfig.cutouts.firstIndex(where: { $0.id == id }) {
                config.geometryConfig.cutouts[idx].x1 = xh
                config.geometryConfig.cutouts[idx].y1 = yh
            }
        case .crackEnd(let id):
            if let idx = config.geometryConfig.cutouts.firstIndex(where: { $0.id == id }) {
                config.geometryConfig.cutouts[idx].x2 = xh
                config.geometryConfig.cutouts[idx].y2 = yh
            }
        case .polygonVertex(let id, let vIdx):
            if let idx = config.geometryConfig.cutouts.firstIndex(where: { $0.id == id }),
               vIdx < config.geometryConfig.cutouts[idx].vertices.count {
                config.geometryConfig.cutouts[idx].vertices[vIdx].x = xh
                config.geometryConfig.cutouts[idx].vertices[vIdx].y = yh
            }
        case .movePolygon(let id, let startXH, let startYH):
            if let idx = config.geometryConfig.cutouts.firstIndex(where: { $0.id == id }),
               let orig = dragStartCutout, !orig.vertices.isEmpty {
                var dx = xh - startXH
                var dy = yh - startYH
                let minX = orig.vertices.map(\.x).min() ?? 0.0
                let maxX = orig.vertices.map(\.x).max() ?? 1.0
                let minY = orig.vertices.map(\.y).min() ?? 0.0
                let maxY = orig.vertices.map(\.y).max() ?? 1.0
                dx = max(-minX, min(1.0 - maxX, dx))
                dy = max(-minY, min(1.0 - maxY, dy))
                for k in 0..<orig.vertices.count {
                    config.geometryConfig.cutouts[idx].vertices[k].x = (round((orig.vertices[k].x + dx) * 200) / 200)
                    config.geometryConfig.cutouts[idx].vertices[k].y = (round((orig.vertices[k].y + dy) * 200) / 200)
                }
                let verts = config.geometryConfig.cutouts[idx].vertices
                config.geometryConfig.cutouts[idx].cx = verts.map(\.x).reduce(0.0, +) / Double(verts.count)
                config.geometryConfig.cutouts[idx].cy = verts.map(\.y).reduce(0.0, +) / Double(verts.count)
            }
        }
    }
    
    // MARK: - Footer Legend
    private var canvasFooterLegend: some View {
        HStack(spacing: 16) {
            HStack(spacing: 5) {
                Circle().fill(BlueHour.ice).frame(width: 8, height: 8)
                Text("Free / Traction-Free Surface Atoms")
                    .font(.caption2)
            }
            HStack(spacing: 5) {
                Circle().fill(Color(red: 1.0, green: 0.42, blue: 0.42)).frame(width: 8, height: 8)
                Text("Fixed (u = 0)")
                    .font(.caption2)
            }
            HStack(spacing: 5) {
                Circle().fill(Color(red: 0.98, green: 0.72, blue: 0.24)).frame(width: 8, height: 8)
                Text("Roller (1-DOF Locked)")
                    .font(.caption2)
            }
            HStack(spacing: 5) {
                Circle().fill(Color(red: 0.31, green: 0.85, blue: 0.78)).frame(width: 8, height: 8)
                Text("Driven Δu(α)")
                    .font(.caption2)
            }
            if config.geometryConfig.indenter.enabled {
                HStack(spacing: 5) {
                    Circle().strokeBorder(Color(red: 1.0, green: 0.82, blue: 0.22), lineWidth: 1.5).frame(width: 8, height: 8)
                    Text("Rigid Indenter Tip (R=\(Int(config.geometryConfig.indenter.radius_a0))a₀)")
                        .font(.caption2)
                        .foregroundColor(Color(red: 1.0, green: 0.82, blue: 0.22))
                }
            }
            Spacer()
            Text("Discrete MTM Lattice: unit-cell spacing a₀ is invariant")
                .font(.caption2.italic())
                .foregroundColor(.secondary)
        }
        .padding(.horizontal, 14)
        .padding(.vertical, 7)
        .background(BlueHour.controlSurface)
    }
    
    // MARK: - Right Inspector Cards
    private var crystalLatticeCard: some View {
        GroupBox(label: Label("1. Discrete Crystal & Rotation (Fills Box)", systemImage: "atom")) {
            VStack(alignment: .leading, spacing: 8) {
                Text("Step 1: The discrete crystal lattice completely fills the rectangular bounding box [0, Lx] × [0, Ly] at rotation angle θ. Step 2: Macroscopic geometry then carves the outer shape, cutouts, and regional boundary conditions.")
                    .font(.caption2)
                    .foregroundColor(.secondary)
                
                HStack {
                    Text("Symmetry:")
                        .font(.caption.bold())
                    Spacer()
                    Picker("", selection: $config.latticeSymmetry) {
                        ForEach(LatticeSymmetry.allCases) { sym in
                            Text(sym.displayName).tag(sym)
                        }
                    }
                    .pickerStyle(.segmented)
                    .controlSize(.small)
                }
                
                Toggle(isOn: $config.useRotation) {
                    Text("Rotate Crystal Lattice (Angle θ)")
                        .font(.caption.bold())
                }
                .toggleStyle(.checkbox)
                
                if config.useRotation {
                    VStack(alignment: .leading, spacing: 6) {
                        HStack {
                            Text("Rotation Angle θ:")
                                .font(.caption2.bold())
                            Spacer()
                            DoubleField(title: "0.0", value: $config.rotationAngle, width: 60, minValue: 0.0, maxValue: 360.0)
                            Text("°").font(.caption2)
                        }
                        Slider(value: $config.rotationAngle, in: 0.0...360.0, step: 0.5)
                        
                        HStack(spacing: 5) {
                            Text("Angles:")
                                .font(.caption2)
                                .foregroundColor(.secondary)
                            ForEach([0.0, 15.0, 30.0, 45.0, 60.0, 90.0], id: \.self) { deg in
                                Button("\(Int(deg))°") {
                                    config.rotationAngle = deg
                                }
                                .buttonStyle(.bordered)
                                .controlSize(.mini)
                            }
                        }
                    }
                    .padding(6)
                    .background(BlueHour.midnight.opacity(0.45), in: RoundedRectangle(cornerRadius: 6))
                }
                
                HStack {
                    Image(systemName: "checkmark.shield")
                        .foregroundColor(Color(red: 0.31, green: 0.85, blue: 0.78))
                    Text("Discrete unit-cell spacing a₀ is strictly preserved.")
                        .font(.system(size: 9.5))
                        .foregroundColor(BlueHour.mist)
                }
            }
            .padding(6)
        }
    }
    
    // MARK: - Microstructure Card (Bicrystal & 2D Voronoi Polycrystal)
    private var microstructureCard: some View {
        GroupBox(label: Label("1b. Grain Boundaries & Polycrystal", systemImage: "square.split.diagonal.2x2")) {
            VStack(alignment: .leading, spacing: 8) {
                Picker("", selection: $config.geometryConfig.microstructure.mode) {
                    Text("Single Crystal").tag(MicrostructureModeModel.singleCrystal)
                    Text("Bicrystal (GB)").tag(MicrostructureModeModel.bicrystal)
                    Text("2D Voronoi").tag(MicrostructureModeModel.voronoi)
                }
                .labelsHidden()
                .pickerStyle(.segmented)
                .controlSize(.small)
                .frame(maxWidth: .infinity)
                .onChange(of: config.geometryConfig.microstructure.mode) { newMode in
                    if newMode != .singleCrystal {
                        config.geometryConfig.enabled = true
                        if newMode == .voronoi && config.geometryConfig.microstructure.seeds.isEmpty {
                            config.geometryConfig.microstructure.regenerateSeeds()
                        }
                    }
                }
                
                if config.geometryConfig.microstructure.mode == .bicrystal {
                    VStack(alignment: .leading, spacing: 6) {
                        HStack(spacing: 4) {
                            Button("Vertical GB") {
                                config.geometryConfig.microstructure.gb_cx = 0.50
                                config.geometryConfig.microstructure.gb_cy = 0.50
                                config.geometryConfig.microstructure.gb_angle_deg = 0.0
                            }
                            .buttonStyle(.bordered)
                            .controlSize(.mini)
                            
                            Button("Horizontal GB") {
                                config.geometryConfig.microstructure.gb_cx = 0.50
                                config.geometryConfig.microstructure.gb_cy = 0.50
                                config.geometryConfig.microstructure.gb_angle_deg = 90.0
                            }
                            .buttonStyle(.bordered)
                            .controlSize(.mini)
                            
                            Button("Sym ±15°") {
                                config.geometryConfig.microstructure.theta1_deg = -15.0
                                config.geometryConfig.microstructure.theta2_deg = +15.0
                            }
                            .buttonStyle(.bordered)
                            .controlSize(.mini)
                        }
                        
                        HStack(spacing: 4) {
                            Text("CSL Presets:")
                                .font(.system(size: 9.5, weight: .semibold))
                                .foregroundColor(.secondary)
                            Button("Σ5 (310) ±18.43°") {
                                let exact = atan2(1.0, 3.0) * 180.0 / .pi
                                config.geometryConfig.microstructure.gb_cx = 0.50
                                config.geometryConfig.microstructure.gb_cy = 0.50
                                config.geometryConfig.microstructure.gb_angle_deg = 0.0
                                config.geometryConfig.microstructure.theta1_deg = -exact
                                config.geometryConfig.microstructure.theta2_deg = +exact
                                config.geometryConfig.microstructure.seam_resolution = .mergeMidpoint
                                config.geometryConfig.microstructure.min_seam_distance_a0 = 0.68
                                config.geometryConfig.microstructure.gb_normal_gap_a0 = 0.0
                                config.geometryConfig.microstructure.gb_tangential_shift_a0 = 0.0
                                config.geometryConfig.microstructure.max_gb_edge_a0 = 2.60
                                config.geometryConfig.boundary_rules = [
                                    BoundaryZoneRuleModel(region: .bottom, kind: .rollerX, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                                    BoundaryZoneRuleModel(region: .top, kind: .rollerX, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                                    BoundaryZoneRuleModel(region: .left, kind: .rollerX, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                                    BoundaryZoneRuleModel(region: .right, kind: .rollerX, thickness_a0: 2.0, ux: 0.0, uy: 0.0)
                                ]
                                config.geometryConfig.primary_stress_component = "sigma_xx"
                                config.noiseLevel = 0.0
                                config.nx = 76
                                config.ny = 76
                            }
                            .buttonStyle(.borderedProminent)
                            .controlSize(.mini)
                            
                            Button("Σ5 (210) ±26.57°") {
                                let exact = atan2(1.0, 2.0) * 180.0 / .pi
                                config.geometryConfig.microstructure.gb_cx = 0.50
                                config.geometryConfig.microstructure.gb_cy = 0.50
                                config.geometryConfig.microstructure.gb_angle_deg = 0.0
                                config.geometryConfig.microstructure.theta1_deg = -exact
                                config.geometryConfig.microstructure.theta2_deg = +exact
                                config.geometryConfig.microstructure.seam_resolution = .mergeMidpoint
                                config.geometryConfig.microstructure.min_seam_distance_a0 = 0.50
                                config.geometryConfig.microstructure.gb_normal_gap_a0 = 0.0
                                config.geometryConfig.microstructure.gb_tangential_shift_a0 = 0.0
                                config.geometryConfig.microstructure.max_gb_edge_a0 = 2.60
                                config.geometryConfig.boundary_rules = [
                                    BoundaryZoneRuleModel(region: .bottom, kind: .rollerX, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                                    BoundaryZoneRuleModel(region: .top, kind: .rollerX, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                                    BoundaryZoneRuleModel(region: .left, kind: .rollerX, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                                    BoundaryZoneRuleModel(region: .right, kind: .rollerX, thickness_a0: 2.0, ux: 0.0, uy: 0.0)
                                ]
                                config.geometryConfig.primary_stress_component = "sigma_xx"
                                config.noiseLevel = 0.0
                                config.nx = 76
                                config.ny = 76
                            }
                            .buttonStyle(.bordered)
                            .controlSize(.mini)
                        }
                        
                        Button(action: {
                            config.geometryConfig.boundary_rules = [
                                BoundaryZoneRuleModel(region: .bottom, kind: .rollerX, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                                BoundaryZoneRuleModel(region: .top, kind: .rollerX, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                                BoundaryZoneRuleModel(region: .left, kind: .rollerX, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                                BoundaryZoneRuleModel(region: .right, kind: .rollerX, thickness_a0: 2.0, ux: 0.0, uy: 0.0)
                            ]
                            config.geometryConfig.primary_stress_component = "sigma_xx"
                        }) {
                            Label("Apply GB Roller BCs (Fix uy = 0, Free ux on All 4 Edges)", systemImage: "arrow.left.and.right.square")
                                .font(.system(size: 10, weight: .semibold))
                                .frame(maxWidth: .infinity)
                        }
                        .buttonStyle(.bordered)
                        .tint(Color(red: 0.98, green: 0.72, blue: 0.24))
                        .controlSize(.mini)
                        
                        HStack(spacing: 4) {
                            Text("CSL-Commensurate N×N:")
                                .font(.system(size: 9.5, weight: .semibold))
                                .foregroundColor(.secondary)
                            Button("64×64 (20L)") {
                                config.nx = 64
                                config.ny = 64
                            }
                            .buttonStyle(.bordered)
                            .controlSize(.mini)
                            Button("76×76 (24L)") {
                                config.nx = 76
                                config.ny = 76
                            }
                            .buttonStyle(.bordered)
                            .controlSize(.mini)
                            Button("95×95 (30L)") {
                                config.nx = 95
                                config.ny = 95
                            }
                            .buttonStyle(.bordered)
                            .controlSize(.mini)
                            if config.noiseLevel > 0.0 {
                                Button("σ=0") {
                                    config.noiseLevel = 0.0
                                }
                                .buttonStyle(.borderedProminent)
                                .tint(.orange)
                                .controlSize(.mini)
                            }
                        }
                        
                        HStack {
                            Text("GB Center (cx, cy):")
                                .font(.caption2.bold())
                            Spacer()
                            DoubleField(title: "0.50", value: $config.geometryConfig.microstructure.gb_cx, width: 46, minValue: 0.05, maxValue: 0.95)
                            DoubleField(title: "0.50", value: $config.geometryConfig.microstructure.gb_cy, width: 46, minValue: 0.05, maxValue: 0.95)
                        }
                        
                        HStack {
                            Text("GB Plane Angle ϕ (°):")
                                .font(.caption2.bold())
                            Spacer()
                            Slider(value: $config.geometryConfig.microstructure.gb_angle_deg, in: -90.0...90.0, step: 1.0)
                                .frame(width: 85)
                            DoubleField(title: "0.0", value: $config.geometryConfig.microstructure.gb_angle_deg, width: 48, minValue: -180.0, maxValue: 180.0)
                        }
                        
                        HStack {
                            Text("Grain 1 Angle θ₁ (°):")
                                .font(.caption2.bold())
                            Spacer()
                            Slider(value: $config.geometryConfig.microstructure.theta1_deg, in: -90.0...90.0, step: 0.5)
                                .frame(width: 85)
                            DoubleField(title: "-15.0", value: $config.geometryConfig.microstructure.theta1_deg, width: 48, minValue: -180.0, maxValue: 180.0)
                        }
                        
                        HStack {
                            Text("Grain 2 Angle θ₂ (°):")
                                .font(.caption2.bold())
                            Spacer()
                            Slider(value: $config.geometryConfig.microstructure.theta2_deg, in: -90.0...90.0, step: 0.5)
                                .frame(width: 85)
                            DoubleField(title: "+15.0", value: $config.geometryConfig.microstructure.theta2_deg, width: 48, minValue: -180.0, maxValue: 180.0)
                        }
                        
                        Divider().padding(.vertical, 1)
                        
                        HStack(spacing: 4) {
                            Text("Grain Separation δ_⊥ (in a₀):")
                                .font(.caption2.bold())
                                .foregroundColor(Color(red: 0.31, green: 0.85, blue: 0.78))
                            Spacer()
                            Slider(value: $config.geometryConfig.microstructure.gb_normal_gap_a0, in: -0.40...1.00, step: 0.02)
                                .frame(width: 78)
                            DoubleField(title: "0.0", value: $config.geometryConfig.microstructure.gb_normal_gap_a0, width: 44, minValue: -0.60, maxValue: 1.50)
                            Button("0") {
                                config.geometryConfig.microstructure.gb_normal_gap_a0 = 0.0
                            }
                            .buttonStyle(.bordered)
                            .controlSize(.mini)
                        }
                        
                        HStack(spacing: 4) {
                            Text("GB Parallel Shift δ_∥ (in a₀):")
                                .font(.caption2.bold())
                            Spacer()
                            Slider(value: $config.geometryConfig.microstructure.gb_tangential_shift_a0, in: -0.50...0.50, step: 0.02)
                                .frame(width: 78)
                            DoubleField(title: "0.0", value: $config.geometryConfig.microstructure.gb_tangential_shift_a0, width: 44, minValue: -1.00, maxValue: 1.00)
                            Button("0") {
                                config.geometryConfig.microstructure.gb_tangential_shift_a0 = 0.0
                            }
                            .buttonStyle(.bordered)
                            .controlSize(.mini)
                        }
                        
                        let misori = abs(config.geometryConfig.microstructure.theta2_deg - config.geometryConfig.microstructure.theta1_deg)
                        HStack {
                            Image(systemName: "angle")
                                .foregroundColor(Color(red: 1.0, green: 0.82, blue: 0.35))
                            Text(String(format: "Misorientation Δθ = %.2f° • Gap δ_⊥ = %.2f a₀", misori, config.geometryConfig.microstructure.gb_normal_gap_a0))
                                .font(.system(size: 9.5, weight: .medium).monospacedDigit())
                                .foregroundColor(Color(red: 1.0, green: 0.82, blue: 0.35))
                        }
                    }
                    .frame(maxWidth: .infinity, alignment: .leading)
                    .padding(6)
                    .background(BlueHour.midnight.opacity(0.45), in: RoundedRectangle(cornerRadius: 6))
                } else if config.geometryConfig.microstructure.mode == .voronoi {
                    VStack(alignment: .leading, spacing: 6) {
                        HStack {
                            Text("Grains N_g:")
                                .font(.caption2.bold())
                            Stepper(
                                "\(config.geometryConfig.microstructure.num_grains)",
                                value: $config.geometryConfig.microstructure.num_grains,
                                in: 2...32,
                                onEditingChanged: { _ in
                                    config.geometryConfig.microstructure.regenerateSeeds()
                                }
                            )
                            .font(.caption2.monospacedDigit())
                            Spacer()
                            Text("Seed:")
                                .font(.caption2)
                            IntField(title: "42", value: $config.geometryConfig.microstructure.random_seed, width: 44)
                            Button("Regenerate") {
                                config.geometryConfig.microstructure.regenerateSeeds()
                            }
                            .buttonStyle(.borderedProminent)
                            .controlSize(.mini)
                        }
                        
                        HStack {
                            Text("Orientation Range [θ_min, θ_max] (°):")
                                .font(.caption2)
                            Spacer()
                            DoubleField(title: "-30", value: $config.geometryConfig.microstructure.theta_min_deg, width: 46, minValue: -180.0, maxValue: 180.0)
                            Text("to").font(.caption2)
                            DoubleField(title: "+30", value: $config.geometryConfig.microstructure.theta_max_deg, width: 46, minValue: -180.0, maxValue: 180.0)
                        }
                        
                        Text("Drag grain seeds (G0, G1, ...) directly on canvas or edit below:")
                            .font(.system(size: 9.5))
                            .foregroundColor(.secondary)
                        
                        ScrollView {
                            VStack(spacing: 4) {
                                ForEach(config.geometryConfig.microstructure.seeds.indices, id: \.self) { idx in
                                    HStack(spacing: 6) {
                                        Text("G\(idx)")
                                            .font(.system(size: 10, weight: .bold, design: .monospaced))
                                            .foregroundColor(grainPaletteColor(idx))
                                            .frame(width: 24, alignment: .leading)
                                        Text("cx:")
                                            .font(.system(size: 9.5))
                                        DoubleField(title: "0.5", value: $config.geometryConfig.microstructure.seeds[idx].cx, width: 42, minValue: 0.02, maxValue: 0.98)
                                        Text("cy:")
                                            .font(.system(size: 9.5))
                                        DoubleField(title: "0.5", value: $config.geometryConfig.microstructure.seeds[idx].cy, width: 42, minValue: 0.02, maxValue: 0.98)
                                        Text("θ:")
                                            .font(.system(size: 9.5))
                                        DoubleField(title: "0.0", value: $config.geometryConfig.microstructure.seeds[idx].theta_deg, width: 46, minValue: -180.0, maxValue: 180.0)
                                        Text("°").font(.system(size: 9.5))
                                    }
                                }
                            }
                        }
                        .frame(maxHeight: 110)
                    }
                    .frame(maxWidth: .infinity, alignment: .leading)
                    .padding(6)
                    .background(BlueHour.midnight.opacity(0.45), in: RoundedRectangle(cornerRadius: 6))
                }
                
                if config.geometryConfig.microstructure.mode != .singleCrystal {
                    VStack(alignment: .leading, spacing: 5) {
                        HStack {
                            Text("Seam Overlap (r_ij < d_min):")
                                .font(.caption2.bold())
                            Spacer()
                            Picker("", selection: $config.geometryConfig.microstructure.seam_resolution) {
                                ForEach(SeamOverlapResolutionModel.allCases) { s in
                                    Text(s.displayName).tag(s)
                                }
                            }
                            .labelsHidden()
                            .controlSize(.mini)
                            .frame(width: 170)
                        }
                        HStack {
                            Text("Min Seam Distance d_min (in a₀):")
                                .font(.caption2)
                            Spacer()
                            Slider(value: $config.geometryConfig.microstructure.min_seam_distance_a0, in: 0.25...0.90, step: 0.05)
                                .frame(width: 85)
                            DoubleField(title: "0.55", value: $config.geometryConfig.microstructure.min_seam_distance_a0, width: 46, minValue: 0.20, maxValue: 0.95)
                        }
                        HStack {
                            Text("Max GB Element Edge L_max (in a₀):")
                                .font(.caption2)
                            Spacer()
                            Slider(value: $config.geometryConfig.microstructure.max_gb_edge_a0, in: 1.50...3.50, step: 0.05)
                                .frame(width: 85)
                            DoubleField(title: "2.60", value: $config.geometryConfig.microstructure.max_gb_edge_a0, width: 46, minValue: 1.40, maxValue: 5.00)
                        }
                        let dMin = config.geometryConfig.microstructure.min_seam_distance_a0
                        let lMax = config.geometryConfig.microstructure.max_gb_edge_a0
                        let isMerge = (config.geometryConfig.microstructure.seam_resolution == .mergeMidpoint)
                        if isMerge && dMin <= 0.62 && lMax >= 2.25 {
                            Text("✓ Defect-Free Kite Seam: merges coincidence apices (d=0) while preserving kite-wing pairs; use Grain Separation δ_⊥ to relieve core compression.")
                                .font(.system(size: 9))
                                .foregroundColor(Color(red: 0.45, green: 0.90, blue: 0.62))
                        } else if isMerge && dMin > 0.62 && lMax >= 2.25 {
                            Text("⚠ Large d_min (> 0.63 a₀) merges kite-wing atoms onto the GB line, which can buckle into opposite zig-zag phases (disconnection). Use d_min ≈ 0.45–0.50 a₀.")
                                .font(.system(size: 9))
                                .foregroundColor(.orange)
                        } else {
                            Text("⚠ Small L_max (< 2.1 a₀) or Delete Mode with large d_min (> 0.70 a₀) removes stretched GB elements to form core voids.")
                                .font(.system(size: 9))
                                .foregroundColor(.orange)
                        }
                    }
                    .frame(maxWidth: .infinity, alignment: .leading)
                    .padding(6)
                    .background(BlueHour.midnight.opacity(0.35), in: RoundedRectangle(cornerRadius: 6))
                }
            }
            .frame(maxWidth: .infinity, alignment: .leading)
            .padding(6)
        }
    }
    
    private func grainPaletteColor(_ gid: Int) -> Color {
        let palette: [Color] = [
            Color(red: 0.56, green: 0.78, blue: 1.00),
            Color(red: 1.00, green: 0.76, blue: 0.42),
            Color(red: 0.62, green: 0.92, blue: 0.68),
            Color(red: 0.92, green: 0.62, blue: 0.95),
            Color(red: 1.00, green: 0.60, blue: 0.62),
            Color(red: 0.52, green: 0.92, blue: 0.90),
            Color(red: 0.95, green: 0.90, blue: 0.52),
            Color(red: 0.75, green: 0.68, blue: 1.00)
        ]
        return palette[abs(gid) % palette.count]
    }
    
    private struct CanvasLatticeSite {
        var xPhys: Double
        var yPhys: Double
        var gid: Int
        var thetaDeg: Double
        var isMergedSeam: Bool = false
    }

    private func snapNearCSLAngle(_ deg: Double) -> Double {
        let r2d = 180.0 / Double.pi
        let cslHalfAngles: [Double] = [
            atan2(1.0, 3.0) * r2d,
            atan2(1.0, 2.0) * r2d,
            atan2(2.0, 3.0) * r2d,
            atan2(1.0, 4.0) * r2d,
            atan2(1.0, 5.0) * r2d,
            atan2(2.0, 5.0) * r2d,
            atan2(3.0, 5.0) * r2d,
            0.5 * atan2(3.0, 4.0) * r2d,
            0.5 * atan2(4.0, 3.0) * r2d
        ]
        for exact in cslHalfAngles {
            if abs(abs(deg) - exact) < 0.02 {
                return (deg < 0.0) ? -exact : exact
            }
        }
        return deg
    }

    private func drawDiscreteLatticeSitesPreview(
        ctx: inout GraphicsContext,
        toScreen: (Double, Double) -> CGPoint,
        boxW: CGFloat,
        boxH: CGFloat
    ) {
        let stepStride = max(1, Int(ceil(Double(max(config.nx, config.ny)) / 125.0)))
        let isTri = (config.latticeSymmetry == .triangular)
        let dotR: CGFloat = max(1.6, min(4.2, min(boxW, boxH) / CGFloat(max(config.nx, config.ny) / stepStride) * 0.30))
        let sqrt3 = sqrt(3.0)
        let a0 = 1.0
        let Lx = max(1e-4, Double(config.nx) * a0)
        let Ly = max(1e-4, isTri ? (Double(config.ny) * a0 * sqrt3 * 0.5) : (Double(config.ny) * a0))
        let scaleStep = Double(stepStride)
        let aCanvas = scaleStep * a0
        let dMinGrid = isTri ? (aCanvas * sqrt3 * 0.5) : aCanvas
        let rBox = hypot(Lx, Ly)
        let nSpan = Int(ceil(rBox / dMinGrid)) + 4
        let micro = config.geometryConfig.microstructure

        var rawSites: [CanvasLatticeSite] = []
        rawSites.reserveCapacity(6400)

        func appendGrainLattice(gid: Int, thetaDeg: Double, originXHat: Double, originYHat: Double, pred: (Double, Double) -> Bool) {
            let snappedDeg = snapNearCSLAngle(thetaDeg)
            let thetaRad = snappedDeg * .pi / 180.0
            let cosT = cos(thetaRad), sinT = sin(thetaRad)
            let a1X = aCanvas * cosT
            let a1Y = aCanvas * sinT
            let a2X = isTri ? (aCanvas * (0.5 * cosT - 0.5 * sqrt3 * sinT)) : (-aCanvas * sinT)
            let a2Y = isTri ? (aCanvas * (0.5 * sinT + 0.5 * sqrt3 * cosT)) : (aCanvas * cosT)
            let ox = originXHat * Lx
            let oy = originYHat * Ly

            for i in -nSpan...nSpan {
                for j in -nSpan...nSpan {
                    let rx = ox + Double(i) * a1X + Double(j) * a2X
                    let ry = oy + Double(i) * a1Y + Double(j) * a2Y
                    guard rx >= -1e-4 && rx <= Lx + 1e-4 && ry >= -1e-4 && ry <= Ly + 1e-4 else { continue }
                    let xHat = min(1.0, max(0.0, rx / Lx))
                    let yHat = min(1.0, max(0.0, ry / Ly))
                    guard pred(xHat, yHat) else { continue }
                    guard config.geometryConfig.isPointInsideSolid(xHat: xHat, yHat: yHat) else { continue }
                    rawSites.append(CanvasLatticeSite(xPhys: xHat * Lx, yPhys: yHat * Ly, gid: gid, thetaDeg: snappedDeg))
                }
            }
        }

        switch micro.mode {
        case .singleCrystal:
            let theta = config.useRotation ? config.rotationAngle : 0.0
            appendGrainLattice(gid: 0, thetaDeg: theta, originXHat: 0.5, originYHat: 0.5) { _, _ in true }
        case .bicrystal:
            let phi = micro.gb_angle_deg * .pi / 180.0
            let nx = cos(phi), ny = sin(phi)
            let tx = -ny, ty = nx
            let deltaPar = micro.gb_tangential_shift_a0 * aCanvas
            let x0 = micro.gb_cx * Lx
            let y0 = micro.gb_cy * Ly
            let cslBand = 0.02 * aCanvas
            let ox0 = (x0 - 0.5 * deltaPar * tx) / Lx
            let oy0 = (y0 - 0.5 * deltaPar * ty) / Ly
            let ox1 = (x0 + 0.5 * deltaPar * tx) / Lx
            let oy1 = (y0 + 0.5 * deltaPar * ty) / Ly
            appendGrainLattice(gid: 0, thetaDeg: micro.theta1_deg, originXHat: ox0, originYHat: oy0) { x, y in
                ((x * Lx - x0) * nx + (y * Ly - y0) * ny) <= cslBand
            }
            appendGrainLattice(gid: 1, thetaDeg: micro.theta2_deg, originXHat: ox1, originYHat: oy1) { x, y in
                ((x * Lx - x0) * nx + (y * Ly - y0) * ny) >= -cslBand
            }
        case .voronoi:
            for (k, s) in micro.seeds.enumerated() {
                appendGrainLattice(gid: k, thetaDeg: s.theta_deg, originXHat: s.cx, originYHat: s.cy) { x, y in
                    micro.grainAt(xHat: x, yHat: y).0 == k
                }
            }
        }

        // Apply live Seam Overlap Filter using fast 2D spatial grid (O(N) on UI thread)
        var filteredSites = rawSites
        if micro.mode != .singleCrystal && filteredSites.count > 1 {
            let dMinPhys = max(0.15, min(0.98, micro.min_seam_distance_a0)) * aCanvas
            let dMinSq = dMinPhys * dMinPhys
            let dSameSafeSq = pow(min(dMinPhys, 0.68 * aCanvas), 2)
            let n = filteredSites.count
            var removed = [Bool](repeating: false, count: n)

            let cellH = max(1e-3, 2.5 * aCanvas)
            let invCell = 1.0 / cellH
            let gNx = max(1, Int(ceil(Lx * invCell)) + 3)
            let gNy = max(1, Int(ceil(Ly * invCell)) + 3)
            var buckets = [[Int]](repeating: [], count: gNx * gNy)

            @inline(__always)
            func cellCoord(_ x: Double, _ y: Double) -> (Int, Int) {
                let ix = min(gNx - 1, max(0, Int(floor(x * invCell)) + 1))
                let iy = min(gNy - 1, max(0, Int(floor(y * invCell)) + 1))
                return (ix, iy)
            }

            for i in 0..<n {
                let (ix, iy) = cellCoord(filteredSites[i].xPhys, filteredSites[i].yPhys)
                buckets[iy * gNx + ix].append(i)
            }

            var crossPairs: [(Int, Int, Double)] = []
            for i in 0..<n {
                let xi = filteredSites[i].xPhys, yi = filteredSites[i].yPhys
                let gi = filteredSites[i].gid
                let (ix, iy) = cellCoord(xi, yi)
                for dy in -1...1 {
                    let ny = iy + dy
                    if ny < 0 || ny >= gNy { continue }
                    for dx in -1...1 {
                        let nx = ix + dx
                        if nx < 0 || nx >= gNx { continue }
                        for j in buckets[ny * gNx + nx] where j > i {
                            if filteredSites[j].gid == gi { continue }
                            let ddx = filteredSites[j].xPhys - xi
                            let ddy = filteredSites[j].yPhys - yi
                            let d2 = ddx * ddx + ddy * ddy
                            if d2 < dMinSq {
                                crossPairs.append((i, j, d2))
                            }
                        }
                    }
                }
            }
            crossPairs.sort { $0.2 < $1.2 }

            for (i, j, _) in crossPairs {
                if removed[i] || removed[j] { continue }
                let dx = filteredSites[j].xPhys - filteredSites[i].xPhys
                let dy = filteredSites[j].yPhys - filteredSites[i].yPhys
                if dx * dx + dy * dy >= dMinSq { continue }

                if micro.seam_resolution == .mergeMidpoint {
                    let mx = 0.5 * (filteredSites[i].xPhys + filteredSites[j].xPhys)
                    let my = 0.5 * (filteredSites[i].yPhys + filteredSites[j].yPhys)
                    let (mix, miy) = cellCoord(mx, my)
                    var crowded = false
                    outerCheck: for dy in -1...1 {
                        let ny = miy + dy
                        if ny < 0 || ny >= gNy { continue }
                        for dx in -1...1 {
                            let nx = mix + dx
                            if nx < 0 || nx >= gNx { continue }
                            for k in buckets[ny * gNx + nx] where !removed[k] && k != i && k != j {
                                let kx = filteredSites[k].xPhys - mx
                                let ky = filteredSites[k].yPhys - my
                                let sameOrSeam = (filteredSites[k].gid == filteredSites[i].gid) || filteredSites[k].isMergedSeam
                                let tol2 = sameOrSeam ? dSameSafeSq : dMinSq
                                if kx * kx + ky * ky < tol2 {
                                    crowded = true
                                    break outerCheck
                                }
                            }
                        }
                    }
                    if !crowded {
                        filteredSites[i].xPhys = mx
                        filteredSites[i].yPhys = my
                        buckets[miy * gNx + mix].append(i)
                    }
                    filteredSites[i].isMergedSeam = true
                }
                removed[j] = true
            }

            // Secondary cleanup pass on spatial buckets to guarantee zero overlapping nodes at Voronoi triple junctions
            for idx in 0..<buckets.count { buckets[idx].removeAll(keepingCapacity: true) }
            for i in 0..<n where !removed[i] {
                let (ix, iy) = cellCoord(filteredSites[i].xPhys, filteredSites[i].yPhys)
                buckets[iy * gNx + ix].append(i)
            }
            for i in 0..<n where !removed[i] {
                let xi = filteredSites[i].xPhys, yi = filteredSites[i].yPhys
                let gi = filteredSites[i].gid
                let (ix, iy) = cellCoord(xi, yi)
                for dy in -1...1 {
                    let ny = iy + dy
                    if ny < 0 || ny >= gNy { continue }
                    for dx in -1...1 {
                        let nx = ix + dx
                        if nx < 0 || nx >= gNx { continue }
                        for j in buckets[ny * gNx + nx] where j > i && !removed[j] {
                            let ddx = filteredSites[j].xPhys - xi
                            let ddy = filteredSites[j].yPhys - yi
                            let bothMergedSeam = filteredSites[i].isMergedSeam && filteredSites[j].isMergedSeam
                            let cutoff2 = (filteredSites[j].gid != gi && !bothMergedSeam) ? dMinSq : dSameSafeSq
                            if ddx * ddx + ddy * ddy < cutoff2 {
                                removed[j] = true
                            }
                        }
                    }
                }
            }

            let phiSep = micro.gb_angle_deg * .pi / 180.0
            let nxSep = cos(phiSep), nySep = sin(phiSep)
            let deltaPerp = (micro.mode == .bicrystal) ? (micro.gb_normal_gap_a0 * aCanvas) : 0.0

            var survivors: [CanvasLatticeSite] = []
            survivors.reserveCapacity(n)
            for i in 0..<n where !removed[i] {
                var site = filteredSites[i]
                if abs(deltaPerp) > 1e-12 && !site.isMergedSeam {
                    let sgn = (site.gid == 0) ? -0.5 : 0.5
                    site.xPhys += sgn * deltaPerp * nxSep
                    site.yPhys += sgn * deltaPerp * nySep
                }
                if site.xPhys >= -1e-4 && site.xPhys <= Lx + 1e-4 && site.yPhys >= -1e-4 && site.yPhys <= Ly + 1e-4 {
                    survivors.append(site)
                }
            }
            filteredSites = survivors

            // Rebuild buckets on survivors for O(N) live cross-seam GB bond preview
            for idx in 0..<buckets.count { buckets[idx].removeAll(keepingCapacity: true) }
            let m = filteredSites.count
            for i in 0..<m {
                let (ix, iy) = cellCoord(filteredSites[i].xPhys, filteredSites[i].yPhys)
                buckets[iy * gNx + ix].append(i)
            }

            let lMaxPhys = micro.max_gb_edge_a0 * aCanvas
            let maxSearchPhys = 2.40 * aCanvas
            let maxSearchSq = maxSearchPhys * maxSearchPhys
            var connectedPath = Path()
            var brokenPath = Path()

            for i in 0..<m {
                let si = filteredSites[i]
                let (ix, iy) = cellCoord(si.xPhys, si.yPhys)
                var nearestCrossIdx = -1
                var nearestD2 = maxSearchSq
                for dy in -1...1 {
                    let ny = iy + dy
                    if ny < 0 || ny >= gNy { continue }
                    for dx in -1...1 {
                        let nx = ix + dx
                        if nx < 0 || nx >= gNx { continue }
                        for j in buckets[ny * gNx + nx] where filteredSites[j].gid != si.gid {
                            let ddx = filteredSites[j].xPhys - si.xPhys
                            let ddy = filteredSites[j].yPhys - si.yPhys
                            let d2 = ddx * ddx + ddy * ddy
                            if d2 < nearestD2 {
                                nearestD2 = d2
                                nearestCrossIdx = j
                            }
                        }
                    }
                }
                if nearestCrossIdx >= 0 {
                    let sj = filteredSites[nearestCrossIdx]
                    let dist = sqrt(nearestD2)
                    let p1 = toScreen(si.xPhys / Lx, si.yPhys / Ly)
                    let p2 = toScreen(sj.xPhys / Lx, sj.yPhys / Ly)
                    if dist <= lMaxPhys {
                        connectedPath.move(to: p1)
                        connectedPath.addLine(to: p2)
                    } else {
                        brokenPath.move(to: p1)
                        brokenPath.addLine(to: p2)
                    }
                }
            }
            ctx.stroke(connectedPath, with: .color(Color(red: 0.45, green: 0.90, blue: 0.62).opacity(0.42)), lineWidth: 1.1)
            ctx.stroke(brokenPath, with: .color(Color(red: 1.0, green: 0.38, blue: 0.36).opacity(0.90)), style: StrokeStyle(lineWidth: 1.6, dash: [3, 2]))
        }

        // Render all surviving lattice sites (highlighting merged seam-core atoms)
        let goldRing = Color(red: 1.0, green: 0.85, blue: 0.35)
        for site in filteredSites {
            let xHat = min(1.0, max(0.0, site.xPhys / Lx))
            let yHat = min(1.0, max(0.0, site.yPhys / Ly))
            let role = config.geometryConfig.classifyNodeRole(
                xHat: xHat, yHat: yHat, nx: config.nx, ny: config.ny
            )
            let pt = toScreen(xHat, yHat)
            let r = site.isMergedSeam ? (dotR * 1.18) : dotR
            let rect = CGRect(x: pt.x - r, y: pt.y - r, width: r * 2, height: r * 2)
            let baseColor = (micro.mode == .singleCrystal) ? BlueHour.ice.opacity(0.62) : grainPaletteColor(site.gid).opacity(0.82)
            let dotColor: Color
            if role == 1 {
                dotColor = Color(red: 1.0, green: 0.42, blue: 0.42)
            } else if role == 2 {
                dotColor = Color(red: 0.31, green: 0.85, blue: 0.78)
            } else if role == 3 || role == 4 {
                dotColor = Color(red: 0.98, green: 0.72, blue: 0.24)
            } else {
                dotColor = baseColor
            }
            ctx.fill(Path(ellipseIn: rect), with: .color(dotColor))
            if site.isMergedSeam {
                ctx.stroke(Path(ellipseIn: rect.insetBy(dx: -0.8, dy: -0.8)), with: .color(goldRing.opacity(0.95)), lineWidth: 1.1)
            }
        }
    }
    
    private func drawMicrostructureOverlay(
        ctx: inout GraphicsContext,
        toScreen: (Double, Double) -> CGPoint,
        originX: CGFloat,
        originY: CGFloat,
        boxW: CGFloat,
        boxH: CGFloat
    ) {
        let micro = config.geometryConfig.microstructure
        let gold = Color(red: 1.0, green: 0.82, blue: 0.35)
        
        if micro.mode == .bicrystal {
            let phi = micro.gb_angle_deg * .pi / 180.0
            // Tangent vector along the GB plane: (-sin(phi), cos(phi))
            let tx = -sin(phi)
            let ty = cos(phi)
            let p1 = toScreen(micro.gb_cx - 1.5 * tx, micro.gb_cy - 1.5 * ty)
            let p2 = toScreen(micro.gb_cx + 1.5 * tx, micro.gb_cy + 1.5 * ty)
            
            var clipCtx = ctx
            clipCtx.clip(to: Path(CGRect(x: originX, y: originY, width: boxW, height: boxH)))
            var gbLine = Path()
            gbLine.move(to: p1)
            gbLine.addLine(to: p2)
            clipCtx.stroke(gbLine, with: .color(gold), style: StrokeStyle(lineWidth: 2.2, dash: [6, 4]))
            
            // Draggable GB center handle
            let cPt = toScreen(micro.gb_cx, micro.gb_cy)
            let hRect = CGRect(x: cPt.x - 6, y: cPt.y - 6, width: 12, height: 12)
            ctx.fill(Path(ellipseIn: hRect), with: .color(gold))
            ctx.stroke(Path(ellipseIn: hRect), with: .color(BlueHour.abyss), lineWidth: 1.5)
            
            let nx = cos(phi), ny = sin(phi)
            let g0Pt = toScreen(min(0.88, max(0.12, micro.gb_cx - 0.18 * nx)), min(0.88, max(0.12, micro.gb_cy - 0.18 * ny)))
            let g1Pt = toScreen(min(0.88, max(0.12, micro.gb_cx + 0.18 * nx)), min(0.88, max(0.12, micro.gb_cy + 0.18 * ny)))
            let pill0 = CGRect(x: g0Pt.x - 42, y: g0Pt.y - 10, width: 84, height: 20)
            let pill1 = CGRect(x: g1Pt.x - 42, y: g1Pt.y - 10, width: 84, height: 20)
            ctx.fill(Path(roundedRect: pill0, cornerRadius: 5), with: .color(BlueHour.abyss.opacity(0.86)))
            ctx.stroke(Path(roundedRect: pill0, cornerRadius: 5), with: .color(grainPaletteColor(0).opacity(0.65)), lineWidth: 1.0)
            ctx.fill(Path(roundedRect: pill1, cornerRadius: 5), with: .color(BlueHour.abyss.opacity(0.86)))
            ctx.stroke(Path(roundedRect: pill1, cornerRadius: 5), with: .color(grainPaletteColor(1).opacity(0.65)), lineWidth: 1.0)
            ctx.draw(
                Text(String(format: "G0 (%.1f°)", micro.theta1_deg))
                    .font(.system(size: 10, weight: .bold, design: .monospaced))
                    .foregroundColor(grainPaletteColor(0)),
                at: g0Pt
            )
            ctx.draw(
                Text(String(format: "G1 (%.1f°)", micro.theta2_deg))
                    .font(.system(size: 10, weight: .bold, design: .monospaced))
                    .foregroundColor(grainPaletteColor(1)),
                at: g1Pt
            )
        } else if micro.mode == .voronoi {
            // Draw analytical Voronoi grain-boundary network (dashed gold bisector edges)
            let seeds = micro.seeds
            let ng = seeds.count
            if ng >= 2 {
                var clipCtx = ctx
                clipCtx.clip(to: Path(CGRect(x: originX, y: originY, width: boxW, height: boxH)))
                var voronoiPath = Path()
                for i in 0..<ng {
                    for j in (i + 1)..<ng {
                        let sx_i = seeds[i].cx, sy_i = seeds[i].cy
                        let sx_j = seeds[j].cx, sy_j = seeds[j].cy
                        let dx = sx_j - sx_i, dy = sy_j - sy_i
                        let dist = hypot(dx, dy)
                        if dist < 1e-6 { continue }
                        let mx = 0.5 * (sx_i + sx_j)
                        let my = 0.5 * (sy_i + sy_j)
                        let tx = -dy / dist
                        let ty = dx / dist
                        var tMin = -2.0
                        var tMax = 2.0

                        // Clip to [0, 1]^2
                        if abs(tx) > 1e-9 {
                            let t0 = (0.0 - mx) / tx
                            let t1 = (1.0 - mx) / tx
                            tMin = max(tMin, min(t0, t1))
                            tMax = min(tMax, max(t0, t1))
                        }
                        if abs(ty) > 1e-9 {
                            let t0 = (0.0 - my) / ty
                            let t1 = (1.0 - my) / ty
                            tMin = max(tMin, min(t0, t1))
                            tMax = min(tMax, max(t0, t1))
                        }

                        // Clip against all other Voronoi seeds k != i, j
                        let r_i2 = (mx - sx_i) * (mx - sx_i) + (my - sy_i) * (my - sy_i)
                        for k in 0..<ng where k != i && k != j {
                            let sx_k = seeds[k].cx, sy_k = seeds[k].cy
                            let a = 2.0 * (tx * (sx_k - sx_i) + ty * (sy_k - sy_i))
                            let b = ((mx - sx_k) * (mx - sx_k) + (my - sy_k) * (my - sy_k)) - r_i2
                            if abs(a) < 1e-9 {
                                if b < 0.0 { tMin = 1.0; tMax = -1.0; break }
                            } else if a > 0.0 {
                                tMax = min(tMax, b / a)
                            } else {
                                tMin = max(tMin, b / a)
                            }
                        }

                        if tMin < tMax {
                            let p1 = toScreen(mx + tMin * tx, my + tMin * ty)
                            let p2 = toScreen(mx + tMax * tx, my + tMax * ty)
                            voronoiPath.move(to: p1)
                            voronoiPath.addLine(to: p2)
                        }
                    }
                }
                clipCtx.stroke(voronoiPath, with: .color(gold.opacity(0.82)), style: StrokeStyle(lineWidth: 1.8, dash: [5, 3]))
            }

            for (idx, s) in micro.seeds.enumerated() {
                let p = toScreen(s.cx, s.cy)
                let col = grainPaletteColor(idx)
                let hRect = CGRect(x: p.x - 5.5, y: p.y - 5.5, width: 11, height: 11)
                ctx.fill(Path(ellipseIn: hRect), with: .color(col))
                ctx.stroke(Path(ellipseIn: hRect), with: .color(BlueHour.abyss), lineWidth: 1.5)
                let labelPt = CGPoint(x: p.x, y: p.y - 15)
                let pillRect = CGRect(x: labelPt.x - 34, y: labelPt.y - 8, width: 68, height: 16)
                ctx.fill(Path(roundedRect: pillRect, cornerRadius: 4), with: .color(BlueHour.abyss.opacity(0.86)))
                ctx.draw(
                    Text(String(format: "G%d (%.0f°)", idx, s.theta_deg))
                        .font(.system(size: 9.5, weight: .bold, design: .monospaced))
                        .foregroundColor(col),
                    at: labelPt
                )
            }
        }
    }
    
    private var presetsCard: some View {
        GroupBox(label: Label("Engineering Benchmark Presets", systemImage: "sparkles.rectangle.stack")) {
            VStack(alignment: .leading, spacing: 8) {
                ForEach(GeometryEngineeringPreset.allCases) { preset in
                    Button(action: {
                        config.geometryConfig.enabled = true
                        config.geometryConfig.applyPreset(preset)
                        config.alphaStart = 0.0
                        config.noiseLevel = 0.0
                        selectedCutoutID = config.geometryConfig.cutouts.first?.id
                    }) {
                        HStack {
                            VStack(alignment: .leading, spacing: 2) {
                                Text(preset.title)
                                    .font(.caption.bold())
                                    .foregroundColor(BlueHour.ice)
                                Text(preset.subtitle)
                                    .font(.system(size: 9.5))
                                    .foregroundColor(BlueHour.mist)
                                    .lineLimit(2)
                            }
                            Spacer()
                            if config.geometryConfig.preset_name == preset.rawValue {
                                Image(systemName: "checkmark.circle.fill")
                                    .foregroundColor(Color(red: 0.31, green: 0.85, blue: 0.78))
                            }
                        }
                        .padding(6)
                        .background(
                            config.geometryConfig.preset_name == preset.rawValue
                                ? BlueHour.slate.opacity(0.65)
                                : BlueHour.midnight.opacity(0.45),
                            in: RoundedRectangle(cornerRadius: 6)
                        )
                    }
                    .buttonStyle(.plain)
                }
            }
            .padding(6)
        }
    }
    
    private var outerDomainCard: some View {
        GroupBox(label: Label("Outer Macroscopic Shape", systemImage: "square.on.circle")) {
            VStack(alignment: .leading, spacing: 8) {
                Picker("Shape", selection: $config.geometryConfig.outer_shape) {
                    ForEach(OuterDomainShape.allCases) { shape in
                        Label(shape.displayName, systemImage: shape.icon).tag(shape)
                    }
                }
                .controlSize(.small)
                
                if config.geometryConfig.outer_shape == .dogbone {
                    HStack {
                        Text("Gauge Width Ratio (w/Lx):")
                            .font(.caption2)
                        Spacer()
                        Slider(value: $config.geometryConfig.dogbone_gauge_width, in: 0.25...0.85)
                            .frame(width: 110)
                        Text(String(format: "%.2f", config.geometryConfig.dogbone_gauge_width))
                            .font(.caption2.monospacedDigit())
                    }
                    HStack {
                        Text("Gauge Length Ratio (l/Ly):")
                            .font(.caption2)
                        Spacer()
                        Slider(value: $config.geometryConfig.dogbone_gauge_length, in: 0.20...0.75)
                            .frame(width: 110)
                        Text(String(format: "%.2f", config.geometryConfig.dogbone_gauge_length))
                            .font(.caption2.monospacedDigit())
                    }
                } else if config.geometryConfig.outer_shape == .disk {
                    HStack {
                        Text("Disk Radius (R/L):")
                            .font(.caption2)
                        Spacer()
                        Slider(value: $config.geometryConfig.disk_radius, in: 0.25...0.50)
                            .frame(width: 110)
                        Text(String(format: "%.2f", config.geometryConfig.disk_radius))
                            .font(.caption2.monospacedDigit())
                    }
                } else if config.geometryConfig.outer_shape == .polygon {
                    VStack(alignment: .leading, spacing: 6) {
                        ScrollView(.horizontal, showsIndicators: false) {
                            HStack(spacing: 5) {
                                Button("Hexagon") {
                                    config.geometryConfig.outer_polygon = (0..<6).map { k in
                                        let a = Double(k) * .pi / 3.0
                                        return NormalizedVertex(x: round((0.5 + 0.46 * cos(a)) * 100) / 100,
                                                                y: round((0.5 + 0.46 * sin(a)) * 100) / 100)
                                    }
                                }
                                Button("Octagon") {
                                    config.geometryConfig.outer_polygon = (0..<8).map { k in
                                        let a = Double(k) * .pi / 4.0 + .pi / 8.0
                                        return NormalizedVertex(x: round((0.5 + 0.48 * cos(a)) * 100) / 100,
                                                                y: round((0.5 + 0.48 * sin(a)) * 100) / 100)
                                    }
                                }
                                Button("Trapezoid") {
                                    config.geometryConfig.outer_polygon = [
                                        NormalizedVertex(x: 0.05, y: 0.02),
                                        NormalizedVertex(x: 0.95, y: 0.02),
                                        NormalizedVertex(x: 0.72, y: 0.98),
                                        NormalizedVertex(x: 0.28, y: 0.98)
                                    ]
                                }
                                Button("Hourglass") {
                                    config.geometryConfig.outer_polygon = [
                                        NormalizedVertex(x: 0.08, y: 0.02),
                                        NormalizedVertex(x: 0.92, y: 0.02),
                                        NormalizedVertex(x: 0.68, y: 0.50),
                                        NormalizedVertex(x: 0.92, y: 0.98),
                                        NormalizedVertex(x: 0.08, y: 0.98),
                                        NormalizedVertex(x: 0.32, y: 0.50)
                                    ]
                                }
                                Button("Reset Box") {
                                    config.geometryConfig.outer_polygon = [
                                        NormalizedVertex(x: 0.02, y: 0.02),
                                        NormalizedVertex(x: 0.98, y: 0.02),
                                        NormalizedVertex(x: 0.98, y: 0.98),
                                        NormalizedVertex(x: 0.02, y: 0.98)
                                    ]
                                }
                            }
                            .buttonStyle(.bordered)
                            .controlSize(.mini)
                        }
                        
                        HStack {
                            Text("Vertices (\(config.geometryConfig.outer_polygon.count)):")
                                .font(.caption2.bold())
                            Spacer()
                            Button("+ Add Vertex") {
                                insertVertexOnLongestEdge(vertices: &config.geometryConfig.outer_polygon)
                            }
                            .buttonStyle(.bordered)
                            .controlSize(.mini)
                        }
                        
                        ForEach(config.geometryConfig.outer_polygon.indices, id: \.self) { idx in
                            HStack(spacing: 6) {
                                Text("V\(idx + 1)")
                                    .font(.caption2.monospacedDigit())
                                    .foregroundColor(.secondary)
                                    .frame(width: 22, alignment: .leading)
                                Text("x:").font(.caption2)
                                DoubleField(title: "0.5", value: $config.geometryConfig.outer_polygon[idx].x, width: 48, minValue: 0.0, maxValue: 1.0)
                                Text("y:").font(.caption2)
                                DoubleField(title: "0.5", value: $config.geometryConfig.outer_polygon[idx].y, width: 48, minValue: 0.0, maxValue: 1.0)
                                Spacer()
                                if config.geometryConfig.outer_polygon.count > 3 {
                                    Button(role: .destructive) {
                                        config.geometryConfig.outer_polygon.remove(at: idx)
                                    } label: {
                                        Image(systemName: "minus.circle")
                                            .font(.caption2)
                                    }
                                    .buttonStyle(.borderless)
                                }
                            }
                        }
                    }
                }
            }
            .padding(6)
        }
    }
    
    private var cutoutsCard: some View {
        GroupBox(label: Label("Holes, Notches & Cracks (\(config.geometryConfig.cutouts.count))", systemImage: "circle.dashed.inset.filled")) {
            VStack(alignment: .leading, spacing: 8) {
                // Add Cutout Buttons
                ScrollView(.horizontal, showsIndicators: false) {
                    HStack(spacing: 6) {
                        Button("+ Circle") { addCutout(.circle) }
                        Button("+ Ellipse") { addCutout(.ellipse) }
                        Button("+ Slot") { addCutout(.rectangle) }
                        Button("+ Notch") { addCutout(.edgeNotch) }
                        Button("+ Crack") { addCutout(.preCrack) }
                        Button("+ Polygon") { addCutout(.polygon) }
                    }
                    .buttonStyle(.bordered)
                    .controlSize(.mini)
                }
                
                if config.geometryConfig.cutouts.isEmpty {
                    Text("No holes or cutouts active (solid outer domain).")
                        .font(.caption2)
                        .foregroundColor(.secondary)
                } else {
                    ForEach($config.geometryConfig.cutouts) { $cutout in
                        VStack(alignment: .leading, spacing: 6) {
                            HStack {
                                Toggle("", isOn: $cutout.enabled)
                                    .labelsHidden()
                                    .controlSize(.mini)
                                Text(cutout.summaryLabel)
                                    .font(.caption.bold())
                                    .foregroundColor(selectedCutoutID == cutout.id ? Color(red: 0.31, green: 0.85, blue: 0.78) : BlueHour.ice)
                                    .onTapGesture { selectedCutoutID = cutout.id }
                                Spacer()
                                Button(role: .destructive, action: {
                                    let idToRemove = cutout.id
                                    config.geometryConfig.cutouts.removeAll { $0.id == idToRemove }
                                }) {
                                    Image(systemName: "trash")
                                        .font(.caption2)
                                }
                                .buttonStyle(.borderless)
                            }
                            
                            if selectedCutoutID == cutout.id {
                                Divider()
                                cutoutDetailEditor(cutout: $cutout)
                            }
                        }
                        .padding(6)
                        .background(
                            selectedCutoutID == cutout.id ? BlueHour.slate.opacity(0.45) : BlueHour.midnight.opacity(0.4),
                            in: RoundedRectangle(cornerRadius: 6)
                        )
                    }
                }
            }
            .padding(6)
        }
    }
    
    @ViewBuilder
    private func cutoutDetailEditor(cutout: Binding<GeometryCutoutItem>) -> some View {
        VStack(alignment: .leading, spacing: 6) {
            Picker("Type", selection: cutout.type) {
                ForEach(CutoutPrimitiveType.allCases) { t in
                    Text(t.displayName).tag(t)
                }
            }
            .controlSize(.mini)
            
            if cutout.wrappedValue.type == .preCrack {
                HStack {
                    Text("Start (x₁, y₁):").font(.caption2)
                    Spacer()
                    Slider(value: cutout.x1, in: 0.0...1.0).frame(width: 65)
                    Slider(value: cutout.y1, in: 0.0...1.0).frame(width: 65)
                }
                HStack {
                    Text("Tip (x₂, y₂):").font(.caption2)
                    Spacer()
                    Slider(value: cutout.x2, in: 0.0...1.0).frame(width: 65)
                    Slider(value: cutout.y2, in: 0.0...1.0).frame(width: 65)
                }
                HStack {
                    Text("Slit Opening:").font(.caption2)
                    Spacer()
                    Slider(value: cutout.thickness, in: 0.002...0.08).frame(width: 100)
                    Text(String(format: "%.3f", cutout.wrappedValue.thickness)).font(.caption2.monospacedDigit())
                }
            } else if cutout.wrappedValue.type == .edgeNotch {
                Picker("Edge", selection: cutout.notch_edge) {
                    Text("Left").tag("left")
                    Text("Right").tag("right")
                    Text("Top").tag("top")
                    Text("Bottom").tag("bottom")
                }
                .pickerStyle(.segmented)
                .controlSize(.mini)
                
                HStack {
                    Text("Notch Depth:").font(.caption2)
                    Spacer()
                    Slider(value: cutout.rx, in: 0.03...0.45).frame(width: 110)
                    Text(String(format: "%.2f", cutout.wrappedValue.rx)).font(.caption2.monospacedDigit())
                }
                HStack {
                    Text("Half-Width:").font(.caption2)
                    Spacer()
                    Slider(value: cutout.ry, in: 0.02...0.30).frame(width: 110)
                    Text(String(format: "%.2f", cutout.wrappedValue.ry)).font(.caption2.monospacedDigit())
                }
            } else if cutout.wrappedValue.type != .polygon {
                HStack {
                    Text("Center X:").font(.caption2)
                    Spacer()
                    Slider(value: cutout.cx, in: 0.05...0.95).frame(width: 110)
                    Text(String(format: "%.2f", cutout.wrappedValue.cx)).font(.caption2.monospacedDigit())
                }
                HStack {
                    Text("Center Y:").font(.caption2)
                    Spacer()
                    Slider(value: cutout.cy, in: 0.05...0.95).frame(width: 110)
                    Text(String(format: "%.2f", cutout.wrappedValue.cy)).font(.caption2.monospacedDigit())
                }
                HStack {
                    Text(cutout.wrappedValue.type == .circle ? "Radius R:" : "Semi-Axis Rx:").font(.caption2)
                    Spacer()
                    Slider(value: Binding(
                        get: { cutout.wrappedValue.rx },
                        set: { val in
                            cutout.wrappedValue.rx = val
                            if cutout.wrappedValue.type == .circle {
                                cutout.wrappedValue.ry = val
                            }
                        }
                    ), in: 0.03...0.45).frame(width: 110)
                    Text(String(format: "%.2f", cutout.wrappedValue.rx)).font(.caption2.monospacedDigit())
                }
                if cutout.wrappedValue.type != .circle {
                    HStack {
                        Text("Semi-Axis Ry:").font(.caption2)
                        Spacer()
                        Slider(value: cutout.ry, in: 0.02...0.45).frame(width: 110)
                        Text(String(format: "%.2f", cutout.wrappedValue.ry)).font(.caption2.monospacedDigit())
                    }
                    HStack {
                        Text("Rotation φ:").font(.caption2)
                        Spacer()
                        Slider(value: cutout.rotation_deg, in: -90.0...90.0).frame(width: 110)
                        Text(String(format: "%.0f°", cutout.wrappedValue.rotation_deg)).font(.caption2.monospacedDigit())
                    }
                }
            } else {
                // Polygon cutout inspector controls
                ScrollView(.horizontal, showsIndicators: false) {
                    HStack(spacing: 5) {
                        Button("Triangle") {
                            let cx = cutout.wrappedValue.cx
                            let cy = cutout.wrappedValue.cy
                            let r = max(0.08, min(0.28, cutout.wrappedValue.rx))
                            cutout.wrappedValue.vertices = (0..<3).map { k in
                                let a = Double(k) * 2.0 * .pi / 3.0 + .pi / 2.0
                                return NormalizedVertex(
                                    x: min(0.98, max(0.02, round((cx + r * cos(a)) * 100) / 100)),
                                    y: min(0.98, max(0.02, round((cy + r * sin(a)) * 100) / 100))
                                )
                            }
                        }
                        Button("Diamond") {
                            let cx = cutout.wrappedValue.cx
                            let cy = cutout.wrappedValue.cy
                            let r = max(0.08, min(0.28, cutout.wrappedValue.rx))
                            cutout.wrappedValue.vertices = [
                                NormalizedVertex(x: cx, y: min(0.98, cy + r)),
                                NormalizedVertex(x: min(0.98, cx + r), y: cy),
                                NormalizedVertex(x: cx, y: max(0.02, cy - r)),
                                NormalizedVertex(x: max(0.02, cx - r), y: cy)
                            ]
                        }
                        Button("Hexagon") {
                            let cx = cutout.wrappedValue.cx
                            let cy = cutout.wrappedValue.cy
                            let r = max(0.08, min(0.28, cutout.wrappedValue.rx))
                            cutout.wrappedValue.vertices = (0..<6).map { k in
                                let a = Double(k) * .pi / 3.0
                                return NormalizedVertex(
                                    x: min(0.98, max(0.02, round((cx + r * cos(a)) * 100) / 100)),
                                    y: min(0.98, max(0.02, round((cy + r * sin(a)) * 100) / 100))
                                )
                            }
                        }
                    }
                    .buttonStyle(.bordered)
                    .controlSize(.mini)
                }
                
                HStack {
                    Text("Vertices (\(cutout.wrappedValue.vertices.count)):")
                        .font(.caption2.bold())
                    Spacer()
                    Button("+ Add Vertex") {
                        insertVertexOnLongestEdge(vertices: &cutout.wrappedValue.vertices)
                    }
                    .buttonStyle(.bordered)
                    .controlSize(.mini)
                }
                
                ForEach(cutout.wrappedValue.vertices.indices, id: \.self) { vIdx in
                    HStack(spacing: 6) {
                        Text("V\(vIdx + 1)")
                            .font(.caption2.monospacedDigit())
                            .foregroundColor(.secondary)
                            .frame(width: 22, alignment: .leading)
                        Text("x:").font(.caption2)
                        DoubleField(title: "0.5", value: cutout.vertices[vIdx].x, width: 48, minValue: 0.0, maxValue: 1.0)
                        Text("y:").font(.caption2)
                        DoubleField(title: "0.5", value: cutout.vertices[vIdx].y, width: 48, minValue: 0.0, maxValue: 1.0)
                        Spacer()
                        if cutout.wrappedValue.vertices.count > 3 {
                            Button(role: .destructive) {
                                cutout.wrappedValue.vertices.remove(at: vIdx)
                            } label: {
                                Image(systemName: "minus.circle")
                                    .font(.caption2)
                            }
                            .buttonStyle(.borderless)
                        }
                    }
                }
            }
        }
    }
    
    private var boundaryConditionsCard: some View {
        GroupBox(label: Label("Boundary Conditions on Geometry", systemImage: "arrow.up.and.down.and.sparkles")) {
            VStack(alignment: .leading, spacing: 8) {
                // Quick BC Templates
                ScrollView(.horizontal, showsIndicators: false) {
                    HStack(spacing: 6) {
                        Button("GB Roller (Fix uy=0, Free ux)") {
                            config.geometryConfig.boundary_rules = [
                                BoundaryZoneRuleModel(region: .bottom, kind: .rollerX, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                                BoundaryZoneRuleModel(region: .top, kind: .rollerX, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                                BoundaryZoneRuleModel(region: .left, kind: .rollerX, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                                BoundaryZoneRuleModel(region: .right, kind: .rollerX, thickness_a0: 2.0, ux: 0.0, uy: 0.0)
                            ]
                            config.geometryConfig.primary_stress_component = "sigma_xx"
                        }
                        Button("Tension (+Y)") {
                            config.geometryConfig.boundary_rules = [
                                BoundaryZoneRuleModel(region: .bottom, kind: .fixed, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                                BoundaryZoneRuleModel(region: .top, kind: .prescribed, thickness_a0: 2.0, ux: 0.0, uy: 1.0),
                                BoundaryZoneRuleModel(region: .left, kind: .free, thickness_a0: 1.5),
                                BoundaryZoneRuleModel(region: .right, kind: .free, thickness_a0: 1.5),
                                BoundaryZoneRuleModel(region: .holeSurfaces, kind: .free, thickness_a0: 1.5)
                            ]
                            config.geometryConfig.primary_stress_component = "sigma_yy"
                        }
                        Button("Compression (-Y)") {
                            config.geometryConfig.boundary_rules = [
                                BoundaryZoneRuleModel(region: .bottom, kind: .fixed, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                                BoundaryZoneRuleModel(region: .top, kind: .prescribed, thickness_a0: 2.0, ux: 0.0, uy: -1.0),
                                BoundaryZoneRuleModel(region: .left, kind: .free, thickness_a0: 1.5),
                                BoundaryZoneRuleModel(region: .right, kind: .free, thickness_a0: 1.5),
                                BoundaryZoneRuleModel(region: .holeSurfaces, kind: .free, thickness_a0: 1.5)
                            ]
                            config.geometryConfig.primary_stress_component = "sigma_yy"
                        }
                        Button("Shear Grips (Top +X)") {
                            config.geometryConfig.boundary_rules = [
                                BoundaryZoneRuleModel(region: .bottom, kind: .fixed, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                                BoundaryZoneRuleModel(region: .top, kind: .prescribed, thickness_a0: 2.0, ux: 1.0, uy: 0.0),
                                BoundaryZoneRuleModel(region: .holeSurfaces, kind: .free, thickness_a0: 1.5)
                            ]
                            config.geometryConfig.primary_stress_component = "sigma_xy"
                        }
                        Button("Affine Shear Box") {
                            config.geometryConfig.boundary_rules = [
                                BoundaryZoneRuleModel(region: .outerPerimeter, kind: .affine, thickness_a0: 2.0, affine_mode: "simple_shear_x"),
                                BoundaryZoneRuleModel(region: .holeSurfaces, kind: .free, thickness_a0: 1.5)
                            ]
                            config.geometryConfig.primary_stress_component = "sigma_xy"
                        }
                    }
                    .buttonStyle(.bordered)
                    .controlSize(.mini)
                }
                
                ForEach($config.geometryConfig.boundary_rules) { $rule in
                    VStack(alignment: .leading, spacing: 5) {
                        HStack {
                            Toggle("", isOn: $rule.enabled)
                                .labelsHidden()
                                .controlSize(.mini)
                            Picker("", selection: $rule.region) {
                                ForEach(BoundaryRegionSelector.allCases) { reg in
                                    Text(reg.displayName).tag(reg)
                                }
                            }
                            .labelsHidden()
                            .controlSize(.mini)
                            
                            Picker("", selection: $rule.kind) {
                                ForEach(BoundaryConditionKindModel.allCases) { k in
                                    Text(k.shortBadge).tag(k)
                                }
                            }
                            .labelsHidden()
                            .controlSize(.mini)
                            .frame(width: 125)
                            .onChange(of: rule.kind) { newKind in
                                if newKind == .rollerX {
                                    rule.lock_x = false
                                    rule.lock_y = true
                                } else if newKind == .rollerY {
                                    rule.lock_x = true
                                    rule.lock_y = false
                                } else if newKind == .fixed || newKind == .prescribed {
                                    if !rule.lock_x && !rule.lock_y {
                                        rule.lock_x = true
                                        rule.lock_y = true
                                    }
                                }
                            }
                            
                            Button(role: .destructive, action: {
                                let idToRemove = rule.id
                                config.geometryConfig.boundary_rules.removeAll { $0.id == idToRemove }
                            }) {
                                Image(systemName: "xmark.circle")
                                    .font(.caption2)
                            }
                            .buttonStyle(.borderless)
                        }

                        if rule.region == .customBox {
                            VStack(alignment: .leading, spacing: 3) {
                                Text("Custom Box Bounds in [0, 1]:").font(.caption2).foregroundColor(BlueHour.ice)
                                HStack(spacing: 6) {
                                    Text("X:").font(.caption2.bold())
                                    DoubleField(title: "0.0", value: $rule.box_x_min, width: 44)
                                    Text("→").font(.caption2).foregroundColor(.secondary)
                                    DoubleField(title: "1.0", value: $rule.box_x_max, width: 44)
                                    
                                    Text("Y:").font(.caption2.bold())
                                    DoubleField(title: "0.0", value: $rule.box_y_min, width: 44)
                                    Text("→").font(.caption2).foregroundColor(.secondary)
                                    DoubleField(title: "1.0", value: $rule.box_y_max, width: 44)
                                }
                            }
                        }
                        
                        if rule.kind == .prescribed {
                            VStack(alignment: .leading, spacing: 4) {
                                HStack(spacing: 6) {
                                    Text("Grip Δu/α:").font(.caption2).foregroundColor(.secondary)
                                    Text("ux:").font(.caption2.bold())
                                    Picker("", selection: Binding(
                                        get: { BoundaryZoneRuleModel.snapUnitDirection(rule.ux) },
                                        set: { rule.ux = $0 }
                                    )) {
                                        Text("-1").tag(-1.0)
                                        Text("0").tag(0.0)
                                        Text("+1").tag(1.0)
                                    }
                                    .pickerStyle(.segmented)
                                    .labelsHidden()
                                    .controlSize(.mini)
                                    .frame(width: 72)
                                    .disabled(!rule.lock_x)

                                    Text("uy:").font(.caption2.bold())
                                    Picker("", selection: Binding(
                                        get: { BoundaryZoneRuleModel.snapUnitDirection(rule.uy) },
                                        set: { rule.uy = $0 }
                                    )) {
                                        Text("-1").tag(-1.0)
                                        Text("0").tag(0.0)
                                        Text("+1").tag(1.0)
                                    }
                                    .pickerStyle(.segmented)
                                    .labelsHidden()
                                    .controlSize(.mini)
                                    .frame(width: 72)
                                    .disabled(!rule.lock_y)

                                    Text("Strip (a₀):").font(.caption2)
                                    DoubleField(title: "2.0", value: $rule.thickness_a0, width: 42)
                                }
                                HStack(spacing: 12) {
                                    Text("Active DOFs:").font(.caption2).foregroundColor(.secondary)
                                    Toggle("Lock ux", isOn: $rule.lock_x)
                                        .font(.caption2)
                                        .controlSize(.mini)
                                    Toggle("Lock uy", isOn: $rule.lock_y)
                                        .font(.caption2)
                                        .controlSize(.mini)
                                    Spacer()
                                    Text(rule.lock_x && rule.lock_y ? "Both Prescribed" : (rule.lock_x ? "Free uy (Roller)" : (rule.lock_y ? "Free ux (Roller)" : "All Free")))
                                        .font(.system(size: 9.5, weight: .medium))
                                        .foregroundColor(Color(red: 0.98, green: 0.72, blue: 0.24))
                                }
                            }
                        } else if rule.kind == .fixed {
                            VStack(alignment: .leading, spacing: 4) {
                                HStack(spacing: 8) {
                                    Text("Clamped Strip Thickness (in a₀):").font(.caption2).foregroundColor(.secondary)
                                    Spacer()
                                    DoubleField(title: "2.0", value: $rule.thickness_a0, width: 52)
                                }
                                HStack(spacing: 12) {
                                    Text("Clamped DOFs:").font(.caption2).foregroundColor(.secondary)
                                    Toggle("Lock ux = 0", isOn: $rule.lock_x)
                                        .font(.caption2)
                                        .controlSize(.mini)
                                    Toggle("Lock uy = 0", isOn: $rule.lock_y)
                                        .font(.caption2)
                                        .controlSize(.mini)
                                    Spacer()
                                }
                            }
                        } else if rule.kind == .rollerX || rule.kind == .rollerY {
                            HStack(spacing: 8) {
                                Text(rule.kind == .rollerX ? "Free ux, Fix uy = 0 (Strip in a₀):" : "Fix ux = 0, Free uy (Strip in a₀):")
                                    .font(.caption2)
                                    .foregroundColor(Color(red: 0.98, green: 0.72, blue: 0.24))
                                Spacer()
                                DoubleField(title: "2.0", value: $rule.thickness_a0, width: 52)
                            }
                        } else if rule.kind == .affine {
                            HStack(spacing: 8) {
                                Picker("Affine Mode", selection: $rule.affine_mode) {
                                    Text("Simple Shear X").tag("simple_shear_x")
                                    Text("Uniaxial Y").tag("uniaxial_y")
                                    Text("Uniaxial X").tag("uniaxial_x")
                                    Text("Pure Shear").tag("pure_shear")
                                    Text("Biaxial").tag("biaxial")
                                }
                                .controlSize(.mini)
                                Text("Strip:").font(.caption2)
                                DoubleField(title: "2.0", value: $rule.thickness_a0, width: 42)
                            }
                        }
                    }
                    .padding(6)
                    .background(BlueHour.midnight.opacity(0.45), in: RoundedRectangle(cornerRadius: 6))
                }
                
                Button("+ Add Boundary Zone Rule") {
                    config.geometryConfig.boundary_rules.append(
                        BoundaryZoneRuleModel(region: .top, kind: .prescribed, thickness_a0: 2.0, ux: 0.0, uy: 1.0)
                    )
                }
                .buttonStyle(.bordered)
                .controlSize(.mini)
            }
            .padding(6)
        }
    }
    
    // MARK: - Helpers
    private func addCutout(_ type: CutoutPrimitiveType) {
        config.geometryConfig.enabled = true
        let item = GeometryCutoutItem(type: type)
        config.geometryConfig.cutouts.append(item)
        selectedCutoutID = item.id
    }
    
    private func insertVertexOnLongestEdge(vertices: inout [NormalizedVertex]) {
        guard vertices.count >= 2 else { return }
        var bestIdx = 0
        var bestLen2 = -1.0
        for i in 0..<vertices.count {
            let j = (i + 1) % vertices.count
            let dx = vertices[j].x - vertices[i].x
            let dy = vertices[j].y - vertices[i].y
            let len2 = dx * dx + dy * dy
            if len2 > bestLen2 {
                bestLen2 = len2
                bestIdx = i
            }
        }
        let nextIdx = (bestIdx + 1) % vertices.count
        let midX = round(0.5 * (vertices[bestIdx].x + vertices[nextIdx].x) * 200) / 200
        let midY = round(0.5 * (vertices[bestIdx].y + vertices[nextIdx].y) * 200) / 200
        vertices.insert(NormalizedVertex(x: midX, y: midY), at: bestIdx + 1)
    }
    
    public func triggerCarvedMeshPreview() {
        let root = config.projectRoot
        let scratchGeomPath = "\(root)/scratch/geometry_config.json"
        if config.geometryConfig.enabled {
            try? config.geometryConfig.writeJSON(to: scratchGeomPath)
            previewManager.geometryJSONPath = scratchGeomPath
        } else {
            previewManager.geometryJSONPath = ""
        }
        previewManager.nx = config.nx
        previewManager.ny = config.ny
        previewManager.latticeSymmetry = config.latticeSymmetry
        previewManager.potentialMode = config.potentialMode
        previewManager.boundaryCondition = config.geometryConfig.enabled ? .fixed : config.boundaryCondition
        previewManager.beta = config.beta
        previewManager.bulkK = config.bulkK
        previewManager.useRotation = config.useRotation
        previewManager.rotationAngle = config.rotationAngle
        previewManager.cslM = config.cslM
        previewManager.cslN = config.cslN
        previewManager.cslKx = config.cslKx
        previewManager.cslKy = config.cslKy
        previewManager.generateInitialMeshPreview()
    }
    
    private var nanoindenterCard: some View {
        NanoindenterCardView(config: config)
    }
    
    private func drawNanoindenterTip(ctx: inout GraphicsContext, toScreen: (Double, Double) -> CGPoint, boxW: CGFloat, boxH: CGFloat) {
        let indenter = config.geometryConfig.indenter
        let cx = indenter.cx
        let isTri = (config.latticeSymmetry == .triangular)
        let sqrt3 = sqrt(3.0)
        let Ly = isTri ? (Double(config.ny) * sqrt3 * 0.5) : Double(config.ny)
        let Lx = Double(config.nx)
        
        let rxNorm = indenter.radius_a0 / max(1.0, Lx)
        let rx = CGFloat(rxNorm) * boxW
        let ry = rx
        let yTop = config.geometryConfig.topSurfaceYHat(at: cx)
        let contactPt = toScreen(cx, yTop)
        
        let depthNormY = indenter.max_depth_a0 / max(1.0, Ly)
        let depthPx = CGFloat(depthNormY) * boxH
        
        let goldFill = Color(red: 1.0, green: 0.82, blue: 0.22).opacity(0.22)
        let goldStroke = Color(red: 1.0, green: 0.82, blue: 0.22)
        let cyanDash = Color(red: 0.31, green: 0.85, blue: 0.78).opacity(0.85)
        
        switch indenter.shape {
        case .cylinder:
            let circleRect = CGRect(x: contactPt.x - rx, y: contactPt.y - 2 * ry, width: 2 * rx, height: 2 * ry)
            var tipPath = Path()
            tipPath.addEllipse(in: circleRect)
            ctx.fill(tipPath, with: .color(goldFill))
            ctx.stroke(tipPath, with: .color(goldStroke), lineWidth: 2.2)
            
            if depthPx > 1 {
                let strokeRect = CGRect(x: contactPt.x - rx, y: contactPt.y - 2 * ry + depthPx, width: 2 * rx, height: 2 * ry)
                var strokePath = Path()
                strokePath.addEllipse(in: strokeRect)
                ctx.stroke(strokePath, with: .color(cyanDash), style: StrokeStyle(lineWidth: 1.5, dash: [4, 3]))
            }
            
            let arrowTop = contactPt.y - min(ry * 1.3, 50.0)
            let arrowBottom = contactPt.y - min(ry * 0.15, 8.0)
            var arrowPath = Path()
            arrowPath.move(to: CGPoint(x: contactPt.x, y: arrowTop))
            arrowPath.addLine(to: CGPoint(x: contactPt.x, y: arrowBottom))
            arrowPath.move(to: CGPoint(x: contactPt.x - 5, y: arrowBottom - 7))
            arrowPath.addLine(to: CGPoint(x: contactPt.x, y: arrowBottom))
            arrowPath.addLine(to: CGPoint(x: contactPt.x + 5, y: arrowBottom - 7))
            ctx.stroke(arrowPath, with: .color(goldStroke), style: StrokeStyle(lineWidth: 2.4, lineCap: .round, lineJoin: .round))
            
        case .flatPunch:
            let h = min(22.0, max(12.0, ry * 0.5))
            let barRect = CGRect(x: contactPt.x - rx, y: contactPt.y - h, width: 2 * rx, height: h)
            let tipPath = Path(roundedRect: barRect, cornerRadius: 3)
            ctx.fill(tipPath, with: .color(goldFill))
            ctx.stroke(tipPath, with: .color(goldStroke), lineWidth: 2.2)
            
            if depthPx > 1 {
                let strokeRect = CGRect(x: contactPt.x - rx, y: contactPt.y - h + depthPx, width: 2 * rx, height: h)
                let strokePath = Path(roundedRect: strokeRect, cornerRadius: 3)
                ctx.stroke(strokePath, with: .color(cyanDash), style: StrokeStyle(lineWidth: 1.5, dash: [4, 3]))
            }
            
            var arrowPath = Path()
            arrowPath.move(to: CGPoint(x: contactPt.x, y: contactPt.y - h - 16))
            arrowPath.addLine(to: CGPoint(x: contactPt.x, y: contactPt.y - h + 2))
            arrowPath.move(to: CGPoint(x: contactPt.x - 5, y: contactPt.y - h - 5))
            arrowPath.addLine(to: CGPoint(x: contactPt.x, y: contactPt.y - h + 2))
            arrowPath.addLine(to: CGPoint(x: contactPt.x + 5, y: contactPt.y - h - 5))
            ctx.stroke(arrowPath, with: .color(goldStroke), style: StrokeStyle(lineWidth: 2.4, lineCap: .round, lineJoin: .round))
            
        case .wedge:
            let theta = min(85.0, max(15.0, indenter.wedge_half_angle_deg)) * .pi / 180.0
            let h = min(50.0, max(20.0, ry / max(0.35, CGFloat(tan(theta)))))
            let halfSpan = min(rx, max(12.0, h * CGFloat(tan(theta))))
            var tipPath = Path()
            tipPath.move(to: contactPt)
            tipPath.addLine(to: CGPoint(x: contactPt.x - halfSpan, y: contactPt.y - h))
            tipPath.addLine(to: CGPoint(x: contactPt.x + halfSpan, y: contactPt.y - h))
            tipPath.closeSubpath()
            ctx.fill(tipPath, with: .color(goldFill))
            ctx.stroke(tipPath, with: .color(goldStroke), lineWidth: 2.2)
            
            if depthPx > 1 {
                var strokePath = Path()
                strokePath.move(to: CGPoint(x: contactPt.x, y: contactPt.y + depthPx))
                strokePath.addLine(to: CGPoint(x: contactPt.x - halfSpan, y: contactPt.y - h + depthPx))
                strokePath.addLine(to: CGPoint(x: contactPt.x + halfSpan, y: contactPt.y - h + depthPx))
                strokePath.closeSubpath()
                ctx.stroke(strokePath, with: .color(cyanDash), style: StrokeStyle(lineWidth: 1.5, dash: [4, 3]))
            }
            
            var arrowPath = Path()
            arrowPath.move(to: CGPoint(x: contactPt.x, y: contactPt.y - h - 14))
            arrowPath.addLine(to: CGPoint(x: contactPt.x, y: contactPt.y - h * 0.35))
            arrowPath.move(to: CGPoint(x: contactPt.x - 5, y: contactPt.y - h * 0.35 - 7))
            arrowPath.addLine(to: CGPoint(x: contactPt.x, y: contactPt.y - h * 0.35))
            arrowPath.addLine(to: CGPoint(x: contactPt.x + 5, y: contactPt.y - h * 0.35 - 7))
            ctx.stroke(arrowPath, with: .color(goldStroke), style: StrokeStyle(lineWidth: 2.4, lineCap: .round, lineJoin: .round))
        }
    }
}

// MARK: - Nanoindenter Card View
struct NanoindenterCardView: View {
    @ObservedObject var config: SimulationConfig
    
    var body: some View {
        GroupBox(label: Label("Nanoindenter Contact Tool", systemImage: "arrow.down.to.line.compact")) {
            VStack(alignment: .leading, spacing: 8) {
                Toggle("Enable Nanoindenter Tool", isOn: $config.geometryConfig.indenter.enabled)
                    .font(.caption.bold())
                    .onChange(of: config.geometryConfig.indenter.enabled) { enabled in
                        if enabled {
                            config.geometryConfig.enabled = true
                            config.geometryConfig.primary_stress_component = "indentation_load"
                        }
                    }
                
                if config.geometryConfig.indenter.enabled {
                    VStack(alignment: .leading, spacing: 8) {
                        Picker("Tip Shape", selection: $config.geometryConfig.indenter.shape) {
                            ForEach(IndenterShapeModel.allCases) { s in
                                Text(s.displayName).tag(s)
                            }
                        }
                        .pickerStyle(.menu)
                        
                        Picker("Contact Enforcement", selection: $config.geometryConfig.indenter.contact_mode) {
                            ForEach(IndenterContactModeModel.allCases) { m in
                                Text(m.displayName).tag(m)
                            }
                        }
                        .pickerStyle(.menu)
                        
                        if config.geometryConfig.indenter.contact_mode == .penalty {
                            HStack {
                                Text("Penalty Stiffness kc:")
                                    .font(.caption2.bold())
                                Spacer()
                                DoubleField(title: "100.0", value: $config.geometryConfig.indenter.penalty_stiffness, width: 60)
                            }
                        }
                        
                        HStack {
                            Text(config.geometryConfig.indenter.shape == .flatPunch ? "Half-Width R (in a₀):" : "Tip Radius / Cap R (in a₀):")
                                .font(.caption2.bold())
                            Spacer()
                            DoubleField(title: "25.0", value: $config.geometryConfig.indenter.radius_a0, width: 60)
                        }
                        
                        if config.geometryConfig.indenter.shape == .wedge {
                            HStack {
                                Text("Wedge Half-Angle θ (°):")
                                    .font(.caption2.bold())
                                Spacer()
                                Slider(value: $config.geometryConfig.indenter.wedge_half_angle_deg, in: 15.0...85.0)
                                    .frame(width: 75)
                                DoubleField(title: "70.3", value: $config.geometryConfig.indenter.wedge_half_angle_deg, width: 48, minValue: 15.0, maxValue: 85.0)
                            }
                        }
                        
                        HStack {
                            Text("Indenter Center X:")
                                .font(.caption2.bold())
                            Spacer()
                            Slider(value: $config.geometryConfig.indenter.cx, in: 0.1...0.9)
                                .frame(width: 80)
                            DoubleField(title: "0.50", value: $config.geometryConfig.indenter.cx, width: 50)
                        }
                        
                        HStack {
                            Text("Max Depth δ_max (in a₀):")
                                .font(.caption2.bold())
                            Spacer()
                            DoubleField(title: "5.0", value: $config.geometryConfig.indenter.max_depth_a0, width: 60)
                        }
                        
                        Toggle("Subsurface Cushion Buffer", isOn: $config.geometryConfig.indenter.cushion_subsurface)
                            .font(.caption2)
                        
                        if config.geometryConfig.indenter.cushion_subsurface {
                            HStack {
                                Text("Cushion Offset (in a₀):")
                                    .font(.caption2)
                                Spacer()
                                DoubleField(title: "0.20", value: $config.geometryConfig.indenter.cushion_offset_a0, width: 60)
                            }
                        }
                        
                        Divider()
                        
                        HStack {
                            VStack(alignment: .leading, spacing: 2) {
                                Text("Save VTK every N steps:")
                                    .font(.caption2.bold())
                                Text("Snapshots for defect nucleation")
                                    .font(.system(size: 8.5))
                                    .foregroundColor(.secondary)
                            }
                            Spacer()
                            IntField(title: "50", value: $config.saveVTKInterval, width: 60)
                        }
                    }
                    .padding(6)
                    .background(BlueHour.midnight.opacity(0.45), in: RoundedRectangle(cornerRadius: 6))
                }
            }
            .padding(6)
        }
    }
}
