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
    
    private enum DragHandleMode {
        case none
        case moveCenter(String)
        case resizeRadius(String)
        case crackStart(String)
        case crackEnd(String)
        case polygonVertex(String, Int)
        case outerPolygonVertex(Int)
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
            HSplitView {
                // LEFT: Interactive Macroscopic Geometry & Discrete Lattice Canvas
                VStack(spacing: 0) {
                    GeometryReader { geo in
                        interactiveCanvas(size: geo.size)
                    }
                    Divider()
                    canvasFooterLegend
                }
                .frame(minWidth: 380, maxWidth: .infinity, maxHeight: .infinity)
                .background(BlueHour.canvasBackground)
                
                // RIGHT: Cutouts & Boundary Condition Inspector Panel
                ScrollView {
                    VStack(alignment: .leading, spacing: 14) {
                        crystalLatticeCard
                        presetsCard
                        outerDomainCard
                        cutoutsCard
                        boundaryConditionsCard
                    }
                    .padding(12)
                }
                .frame(minWidth: 320, idealWidth: 360, maxWidth: 420)
                .background(BlueHour.sidebarBackground)
            }
        }
        .onAppear {
            if selectedCutoutID == nil {
                selectedCutoutID = config.geometryConfig.cutouts.first?.id
            }
        }
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
        let pad: CGFloat = 44.0
        let availW = max(100.0, size.width - 2 * pad)
        let availH = max(100.0, size.height - 2 * pad)
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
        let originY = (size.height - boxH) * 0.5
        
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
                
                // 4. Discrete Crystal Lattice Sites (Square or Triangular, Color-Coded by BC Role)
                // Step 1: Crystal lattice (at rotation theta) completely fills the rectangular box [0, 1] x [0, 1].
                // Step 2: Solid mask carves macroscopic outer boundary & cutouts from this filled crystal.
                if showDiscreteLatticePreview {
                    let stepStride = max(1, Int(ceil(Double(max(config.nx, config.ny)) / 75.0)))
                    let isTri = (config.latticeSymmetry == .triangular)
                    let dotR: CGFloat = max(1.6, min(4.2, min(boxW, boxH) / CGFloat(max(config.nx, config.ny) / stepStride) * 0.30))
                    
                    let thetaRad = (config.useRotation ? config.rotationAngle : 0.0) * .pi / 180.0
                    let cosT = cos(thetaRad), sinT = sin(thetaRad)
                    let sqrt3 = sqrt(3.0)
                    
                    let a0 = 1.0
                    let Lx = Double(config.nx) * a0
                    let Ly = isTri ? (Double(config.ny) * a0 * sqrt3 * 0.5) : (Double(config.ny) * a0)
                    let dMin = isTri ? (a0 * sqrt3 * 0.5) : a0
                    
                    // Stride-scaled spacing in normalized [0, 1] x [0, 1] coordinates
                    let scaleStep = Double(stepStride)
                    let a1HatX = (scaleStep * a0 * cosT) / max(1e-4, Lx)
                    let a1HatY = (scaleStep * a0 * sinT) / max(1e-4, Ly)
                    let a2HatX = (isTri ? (scaleStep * a0 * (0.5 * cosT - 0.5 * sqrt3 * sinT)) : (-scaleStep * a0 * sinT)) / max(1e-4, Lx)
                    let a2HatY = (isTri ? (scaleStep * a0 * (0.5 * sinT + 0.5 * sqrt3 * cosT)) : (scaleStep * a0 * cosT)) / max(1e-4, Ly)
                    
                    let rBox = hypot(Lx * 0.5, Ly * 0.5)
                    let nSpan = Int(ceil(rBox / (dMin * scaleStep))) + 4
                    
                    for i in -nSpan...nSpan {
                        for j in -nSpan...nSpan {
                            let xHat = 0.5 + Double(i) * a1HatX + Double(j) * a2HatX
                            let yHat = 0.5 + Double(i) * a1HatY + Double(j) * a2HatY
                            
                            guard xHat >= -1e-5 && xHat <= 1.0 + 1e-5 && yHat >= -1e-5 && yHat <= 1.0 + 1e-5 else { continue }
                            let clampedX = min(1.0, max(0.0, xHat))
                            let clampedY = min(1.0, max(0.0, yHat))
                            guard config.geometryConfig.isPointInsideSolid(xHat: clampedX, yHat: clampedY) else { continue }
                            
                            let role = config.geometryConfig.classifyNodeRole(
                                xHat: clampedX, yHat: clampedY, nx: config.nx, ny: config.ny
                            )
                            let pt = toScreen(clampedX, clampedY)
                            let rect = CGRect(x: pt.x - dotR, y: pt.y - dotR, width: dotR * 2, height: dotR * 2)
                            let dotColor: Color
                            if role == 1 {
                                dotColor = Color(red: 1.0, green: 0.42, blue: 0.42) // Fixed Clamped
                            } else if role == 2 {
                                dotColor = Color(red: 0.31, green: 0.85, blue: 0.78) // Driven Prescribed
                            } else {
                                dotColor = BlueHour.ice.opacity(0.62) // Free Interior / Free Surface
                            }
                            ctx.fill(Path(ellipseIn: rect), with: .color(dotColor))
                        }
                    }
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
            
            // Top-Left Overlay HUD Badge
            VStack(alignment: .leading, spacing: 4) {
                HStack(spacing: 6) {
                    Image(systemName: "hand.draw.fill")
                        .foregroundColor(BlueHour.ice)
                    Text("Interactive Macroscopic Geometry Canvas")
                        .font(.system(size: 11, weight: .bold))
                        .foregroundColor(BlueHour.ice)
                }
                Text("Drag center handle to move hole • Drag cyan rim handle to resize • Click 'Generate & View Carved Triangulation' to mesh")
                    .font(.system(size: 9.5))
                    .foregroundColor(BlueHour.mist)
            }
            .padding(.horizontal, 10)
            .padding(.vertical, 6)
            .background(BlueHour.midnight.opacity(0.88), in: RoundedRectangle(cornerRadius: 8))
            .overlay(RoundedRectangle(cornerRadius: 8).strokeBorder(BlueHour.border, lineWidth: 1))
            .frame(maxWidth: .infinity, maxHeight: .infinity, alignment: .topLeading)
            .padding(10)
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
        
        for rule in config.geometryConfig.boundary_rules where rule.enabled {
            if rule.region == .customBox {
                let bx = originX + CGFloat(rule.box_x_min) * boxW
                let by = originY + CGFloat(1.0 - rule.box_y_max) * boxH
                let bw = CGFloat(max(0.01, rule.box_x_max - rule.box_x_min)) * boxW
                let bh = CGFloat(max(0.01, rule.box_y_max - rule.box_y_min)) * boxH
                let rect = CGRect(x: bx, y: by, width: bw, height: bh)
                let col = (rule.kind == .fixed) ? fixedColor : (rule.kind == .free ? Color.gray : drivenColor)
                ctx.fill(Path(rect), with: .color(col.opacity(0.12)))
                ctx.stroke(Path(rect), with: .color(col.opacity(0.85)), style: StrokeStyle(lineWidth: 1.5, dash: [4, 3]))
            }

            if rule.kind == .fixed {
                // Draw subtle clamped hash marks along fixed edge
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
                    ctx.stroke(tri, with: .color(fixedColor), lineWidth: 1.6)
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
            // Hit-test outer polygon vertices first
            if config.geometryConfig.outer_shape == .polygon {
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
                        for (vIdx, v) in cutout.vertices.enumerated() {
                            let p = toScreen(v.x, v.y)
                            if hypot(startLocation.x - p.x, startLocation.y - p.y) < 14 {
                                selectedCutoutID = cutout.id
                                activeDragMode = .polygonVertex(cutout.id, vIdx)
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
                Text("Fixed Clamped Grip Atoms (u = 0)")
                    .font(.caption2)
            }
            HStack(spacing: 5) {
                Circle().fill(Color(red: 0.31, green: 0.85, blue: 0.78)).frame(width: 8, height: 8)
                Text("Driven Prescribed Grip Atoms Δu(α)")
                    .font(.caption2)
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
            }
        }
    }
    
    private var boundaryConditionsCard: some View {
        GroupBox(label: Label("Boundary Conditions on Geometry", systemImage: "arrow.up.and.down.and.sparkles")) {
            VStack(alignment: .leading, spacing: 8) {
                // Quick BC Templates
                ScrollView(.horizontal, showsIndicators: false) {
                    HStack(spacing: 6) {
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
                            .frame(width: 95)
                            
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
                            HStack(spacing: 8) {
                                Text("Grip Δu/α:").font(.caption2).foregroundColor(.secondary)
                                Text("ux:").font(.caption2.bold())
                                DoubleField(title: "0.0", value: $rule.ux, width: 48)
                                Text("uy:").font(.caption2.bold())
                                DoubleField(title: "1.0", value: $rule.uy, width: 48)
                                Text("Strip (a₀):").font(.caption2)
                                DoubleField(title: "2.0", value: $rule.thickness_a0, width: 42)
                            }
                        } else if rule.kind == .fixed {
                            HStack(spacing: 8) {
                                Text("Clamped Strip Thickness (in a₀):").font(.caption2).foregroundColor(.secondary)
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
}
