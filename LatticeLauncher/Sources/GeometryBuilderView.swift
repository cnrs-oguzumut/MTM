import SwiftUI
import AppKit

public enum GeometryBuilderDisplayMode {
    case canvasOnly
    case inspectorMicrostructure
    case inspectorShapeAndBCs
    case splitCanvasAndInspector
}

public struct GeometryBuilderView: View {
    @ObservedObject var config: SimulationConfig
    @ObservedObject var previewManager: StiffnessStudyManager
    var displayMode: GeometryBuilderDisplayMode = .canvasOnly
    var onSwitchToCarvedMesh: () -> Void
    var onOpenDocs: (DocumentationTopic) -> Void = { _ in }
    
    @State private var selectedCutoutID: String? = nil
    @State private var selectedRuleID: String? = nil
    @State private var showDiscreteLatticePreview: Bool = true
    @State private var showBCArrows: Bool = true
    @State private var activeDragMode: DragHandleMode = .none
    @State private var dragStartCutout: GeometryCutoutItem? = nil
    @State private var inspectorSection: InspectorSection = .microstructure
    @State private var autoUpdateGBSeamParams: Bool = true
    @State private var approvedLatticeSiteCount: Int? = nil
    
    /// Rough number of lattice sites the blueprint would draw.
    private var estimatedLatticeSiteCount: Int { max(0, config.nx * config.ny) }
    
    /// Large crystal and the user has not asked to see every site → draw no individual sites.
    private var latticeSitesAutoHidden: Bool {
        LargeMeshPolicy.isLarge(estimatedLatticeSiteCount) && approvedLatticeSiteCount != estimatedLatticeSiteCount
    }
    @AppStorage("geometryBCPreviewStage") private var previewStageRaw: String = BCStagePreviewMode.stage2Load.rawValue
    
    private var previewStage: BCStagePreviewMode {
        get { BCStagePreviewMode(rawValue: previewStageRaw) ?? .stage2Load }
        nonmutating set { previewStageRaw = newValue.rawValue }
    }
    
    private enum BCStagePreviewMode: String, CaseIterable, Identifiable {
        case stage1Relax = "Relaxation Step BCs"
        case stage2Load = "Loading Step BCs & Indenter"
        var id: String { rawValue }
    }
    
    private var isViewingStage1RelaxBCs: Bool {
        if config.protocolMode != .relaxThenLoad {
            return true
        }
        return previewStage == .stage1Relax
    }
    
    private var activeEditingRulesBinding: Binding<[BoundaryZoneRuleModel]> {
        isViewingStage1RelaxBCs
            ? $config.geometryConfig.stage1_boundary_rules
            : $config.geometryConfig.boundary_rules
    }
    
    private var activeCanvasBoundaryRules: [BoundaryZoneRuleModel] {
        if isViewingStage1RelaxBCs {
            return config.geometryConfig.stage1RelaxationRules()
        }
        if config.geometryConfig.indenter.enabled {
            return config.geometryConfig.boundary_rules.map { rule in
                if rule.region == .top && rule.kind != .free {
                    var copy = rule
                    copy.kind = .free
                    return copy
                }
                return rule
            }
        }
        return config.geometryConfig.boundary_rules
    }
    
    private var activeCanvasIndenterEnabled: Bool {
        if isViewingStage1RelaxBCs {
            return false
        }
        return config.geometryConfig.indenter.enabled
    }
    
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
        case moveDislocationCore(Int)
    }
    
    public init(
        config: SimulationConfig,
        previewManager: StiffnessStudyManager,
        displayMode: GeometryBuilderDisplayMode = .canvasOnly,
        onSwitchToCarvedMesh: @escaping () -> Void,
        onOpenDocs: @escaping (DocumentationTopic) -> Void = { _ in }
    ) {
        self.config = config
        self.previewManager = previewManager
        self.displayMode = displayMode
        self.onSwitchToCarvedMesh = onSwitchToCarvedMesh
        self.onOpenDocs = onOpenDocs
    }
    
    public var body: some View {
        Group {
            switch displayMode {
            case .canvasOnly:
                VStack(spacing: 0) {
                    studioTopBar
                    Divider()
                    VStack(spacing: 0) {
                        canvasHeaderHintBar
                        Divider()
                        GeometryReader { geo in
                            interactiveCanvas(size: geo.size)
                        }
                        .overlay(alignment: .bottom) {
                            if showDiscreteLatticePreview && latticeSitesAutoHidden {
                                NodesHiddenNote(nodeCount: estimatedLatticeSiteCount, noun: "lattice sites") {
                                    approvedLatticeSiteCount = estimatedLatticeSiteCount
                                }
                                .padding(.bottom, 10)
                            }
                        }
                        Divider()
                        canvasFooterLegend
                    }
                    .frame(minWidth: 380, maxWidth: .infinity, maxHeight: .infinity)
                    .background(BlueHour.canvasBackground)
                }
            case .inspectorMicrostructure:
                VStack(alignment: .leading, spacing: 12) {
                    microstructureCard
                }
                .frame(maxWidth: .infinity, alignment: .leading)
            case .inspectorShapeAndBCs:
                VStack(alignment: .leading, spacing: 12) {
                    presetsCard
                    outerDomainCard
                    cutoutsCard
                    stageProtocolBannerCard
                    boundaryConditionsCard
                    nanoindenterCard
                        .onChange(of: config.geometryConfig.indenter.enabled) { enabled in
                            if enabled {
                                if config.protocolMode == .relaxOnly {
                                    config.protocolMode = .relaxThenLoad
                                    config.geometryConfig.microstructure.relax_only = false
                                }
                                previewStage = .stage2Load
                            }
                        }
                }
                .frame(maxWidth: .infinity, alignment: .leading)
            case .splitCanvasAndInspector:
                VStack(spacing: 0) {
                    studioTopBar
                    Divider()
                    HStack(spacing: 0) {
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
                    }
                }
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
    
    private var canvasSurfaceBadgeTitle: String {
        if isViewingStage1RelaxBCs {
            return "STAGE 1 RELAXATION (Free/Roller BCs, α = 0)"
        }
        switch config.activeMechanicalLoading {
        case .simpleShearPos:
            return "STAGE 2: SIMPLE SHEAR (+γ_xy) • Driven Top ux → • Fixed Base"
        case .simpleShearNeg:
            return "STAGE 2: REVERSE SHEAR (-γ_xy) • Driven Top ux ← • Fixed Base"
        case .uniaxialTensionY:
            return "STAGE 2: UNIAXIAL TENSION Y (+ε_yy) • Pulling Top uy ↑ • Fixed Base"
        case .uniaxialCompY:
            return "STAGE 2: COMPRESSION Y (-ε_yy) • Pushing Top uy ↓ • Fixed Base"
        case .uniaxialTensionX:
            return "STAGE 2: UNIAXIAL TENSION X (+ε_xx) • Pulling Right ux → • Rollers"
        case .shearParallelGB:
            return "STAGE 2: SHEAR ∥ VERTICAL GB (+uy) • Driven Right Wall"
        case .affineShearBox:
            return "STAGE 2: AFFINE SHEAR BOX (4 WALLS) • Linear Shear Gradient F_xy"
        case .nanoindentation:
            return "STAGE 2: NANOINDENTATION • Tip δ_max = \(String(format: "%.1f", config.geometryConfig.indenter.max_depth_a0)) a₀"
        default:
            return "SURFACE BCs · \(config.activeMechanicalLoading.displayName.uppercased()) (\(config.geometryConfig.primary_stress_component))"
        }
    }

    private var canvasSurfaceBadgeColor: Color {
        if isViewingStage1RelaxBCs {
            return Color(red: 0.45, green: 0.92, blue: 0.62)
        }
        if config.activeMechanicalLoading == .uniaxialTensionY || config.activeMechanicalLoading == .uniaxialCompY {
            return Color(red: 1.0, green: 0.80, blue: 0.40)
        }
        if config.activeMechanicalLoading == .nanoindentation {
            return Color(red: 1.0, green: 0.55, blue: 0.35)
        }
        return Color(red: 0.31, green: 0.85, blue: 0.78)
    }

    // MARK: - Top Toolbar
    private var studioTopBar: some View {
        ScrollView(.horizontal, showsIndicators: false) {
            HStack(spacing: 12) {
                Toggle(isOn: Binding(
                    get: { !config.isEffectivePeriodicPBC },
                    set: { isSurface in
                        config.setBoundaryRegime(isSurface ? .finiteSurfaceBCs : .periodicPBC)
                    }
                )) {
                    Label("Surface Grips / Custom BCs", systemImage: "slider.horizontal.2.square.on.square")
                        .font(.caption.bold())
                }
                .toggleStyle(.switch)
                .controlSize(.small)
                
                Divider().frame(height: 18)
                
                Toggle("Show Discrete Lattice Sites", isOn: $showDiscreteLatticePreview)
                    .font(.caption)
                    .toggleStyle(.checkbox)
                
                Toggle("Show BC Vectors", isOn: $showBCArrows)
                    .font(.caption)
                    .toggleStyle(.checkbox)
                
                Divider().frame(height: 18)
                
                Button(action: {
                    if !config.isEffectivePeriodicPBC {
                        config.geometryConfig.enabled = true
                    }
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
                
                // Active Physics / BC Loading Interactive Menu
                Menu {
                    Section("Boundary Regime") {
                        Button(action: { config.setBoundaryRegime(.periodicPBC) }) {
                            if config.isEffectivePeriodicPBC {
                                Label("Full 2D Periodic (PBC)", systemImage: "checkmark")
                            } else {
                                Label("Full 2D Periodic (PBC)", systemImage: "arrow.triangle.2.circlepath")
                            }
                        }
                        Button(action: { config.setBoundaryRegime(.finiteSurfaceBCs) }) {
                            if !config.isEffectivePeriodicPBC {
                                Label("Surface Grips / Custom BCs", systemImage: "checkmark")
                            } else {
                                Label("Surface Grips / Custom BCs", systemImage: "slider.horizontal.2.square.on.square")
                            }
                        }
                    }
                    Section("Mechanical Deformation / Loading") {
                        ForEach(ActiveMechanicalLoading.primaryGridModes) { load in
                            Button(action: { config.applyMechanicalLoading(load) }) {
                                if config.activeMechanicalLoading == load {
                                    Label(load.shortTitle, systemImage: "checkmark")
                                } else {
                                    Label(load.shortTitle, systemImage: load.sfSymbol)
                                }
                            }
                        }
                    }
                } label: {
                    HStack(spacing: 5) {
                        Circle()
                            .fill(config.isEffectivePeriodicPBC ? Color(red: 0.45, green: 0.92, blue: 0.62) : Color(red: 0.31, green: 0.85, blue: 0.78))
                            .frame(width: 7, height: 7)
                        Text(config.activePhysicsSummaryBadge)
                            .font(.system(size: 10, weight: .bold, design: .monospaced))
                            .foregroundColor(config.isEffectivePeriodicPBC ? Color(red: 0.45, green: 0.92, blue: 0.62) : BlueHour.ice)
                        Image(systemName: "chevron.up.chevron.down")
                            .font(.system(size: 8))
                            .foregroundColor(BlueHour.mist)
                    }
                    .padding(.horizontal, 8)
                    .padding(.vertical, 4)
                    .background(BlueHour.pillSurface, in: Capsule())
                    .overlay(Capsule().strokeBorder(config.isEffectivePeriodicPBC ? Color.green.opacity(0.4) : BlueHour.border, lineWidth: 1))
                }
                .menuStyle(.borderlessButton)
                
                HelpPopoverButton(
                    title: "Macroscopic Geometry & BC Studio",
                    summary: "Design your specimen in 3 steps using the right-hand inspector:\n1. Grains & Dislocs: Single Crystal, Bicrystal, Tri-Crystal, Voronoi, or Periodic Dislocations.\n2. Shape & Holes: Carve notches, dogbones, circular/elliptical holes, or cracks.\n3. BCs & Indenter: Configure Stage 1 (Relaxation) and Stage 2 (Loading) boundary conditions or Nanoindentation.",
                    formulaOrTip: "• Drag handles directly on the canvas to move grain boundaries, Voronoi seeds, cutouts, or dislocation cores.\n• Click 'Generate & View Carved Triangulation' to inspect the exact finite-element mesh.",
                    docTopic: .geometryAndBCs,
                    onOpenDocs: onOpenDocs
                )
            }
            .padding(.horizontal, 12)
            .padding(.vertical, 6)
        }
        .frame(height: 36)
        .background(BlueHour.controlSurface)
    }
    
    // MARK: - Reusable High-Contrast Choice Chip
    private func selectableChipButton(
        title: String,
        isSelected: Bool,
        activeFill: Color = Color(red: 0.12, green: 0.64, blue: 0.54),
        activeBorder: Color = Color(red: 0.42, green: 0.96, blue: 0.84),
        action: @escaping () -> Void
    ) -> some View {
        Button(action: action) {
            HStack(spacing: 3.5) {
                if isSelected {
                    Image(systemName: "checkmark.circle.fill")
                        .font(.system(size: 8.5, weight: .bold))
                        .foregroundColor(.white)
                }
                Text(title)
                    .font(.system(size: 9.5, weight: isSelected ? .bold : .medium, design: .monospaced))
                    .foregroundColor(isSelected ? .white : BlueHour.ice.opacity(0.90))
                    .lineLimit(1)
            }
            .padding(.horizontal, 6.5)
            .padding(.vertical, 3.5)
            .background(
                isSelected ? activeFill : BlueHour.midnight.opacity(0.88),
                in: RoundedRectangle(cornerRadius: 5)
            )
            .overlay(
                RoundedRectangle(cornerRadius: 5)
                    .strokeBorder(
                        isSelected ? activeBorder : BlueHour.steel.opacity(0.58),
                        lineWidth: isSelected ? 1.35 : 0.95
                    )
            )
        }
        .buttonStyle(.plain)
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
                //    When large crystals hide interior sites for speed, we still draw the perimeter boundary sites
                //    so boundary conditions (fixed, driven, roller, free) remain immediately visible on every wall!
                if showDiscreteLatticePreview {
                    if !latticeSitesAutoHidden {
                        drawDiscreteLatticeSitesPreview(ctx: &ctx, toScreen: toScreen, boxW: boxW, boxH: boxH)
                    } else {
                        drawPerimeterLatticeSitesPreview(ctx: &ctx, toScreen: toScreen, boxW: boxW, boxH: boxH)
                    }
                }
                
                // 4b. Bicrystal GB Plane / Voronoi Grain Cell Boundaries & Seed Handles
                if config.geometryConfig.microstructure.mode != .singleCrystal {
                    drawMicrostructureOverlay(ctx: &ctx, toScreen: toScreen, originX: originX, originY: originY, boxW: boxW, boxH: boxH)
                }
                
                // 5. Boundary Condition Velocity Arrows on Driven Edges (suppressed in Full 2D Periodic Dislocations mode)
                if showBCArrows && (config.isEffectivePeriodicPBC || config.geometryConfig.microstructure.mode != .periodicDislocations) {
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
                
                // 8. Nanoindenter Tip Visual Preview (always visible when enabled so user can position it in either tab)
                if config.geometryConfig.indenter.enabled && config.geometryConfig.microstructure.mode != .periodicDislocations {
                    drawNanoindenterTip(ctx: &ctx, toScreen: toScreen, boxW: boxW, boxH: boxH)
                }
            }
            .overlay(alignment: .topLeading) {
                let isPeriodic = config.isEffectivePeriodicPBC
                let isPeriodicDislocs = (config.geometryConfig.microstructure.mode == .periodicDislocations)
                HStack(spacing: 6) {
                    Circle()
                        .fill((isPeriodic || isViewingStage1RelaxBCs) ? Color.green : canvasSurfaceBadgeColor)
                        .frame(width: 7, height: 7)
                    if isPeriodic {
                        if isPeriodicDislocs {
                            Text(config.protocolMode == .relaxOnly
                                 ? "FULL 2D PERIODIC TORUS (PBC) · Stage 1 Core Relaxation Only (α = 0)"
                                 : "FULL 2D PERIODIC TORUS (PBC) · Stage 1 Relax (α = 0) → Stage 2 Shear (\(config.mode == .negative ? "-α" : "+α"))")
                                .font(.system(size: 10, weight: .bold, design: .monospaced))
                                .foregroundColor(Color(red: 0.45, green: 0.92, blue: 0.62))
                        } else {
                            Text("FULL 2D PERIODIC TORUS (PBC) · SIMPLE SHEAR γ_xy (\(config.mode == .negative ? "-α" : "+α"))")
                                .font(.system(size: 10, weight: .bold, design: .monospaced))
                                .foregroundColor(Color(red: 0.45, green: 0.92, blue: 0.62))
                        }
                    } else {
                        Text(canvasSurfaceBadgeTitle)
                            .font(.system(size: 10, weight: .bold, design: .monospaced))
                            .foregroundColor(canvasSurfaceBadgeColor)
                        if config.protocolMode == .relaxThenLoad {
                            Button(isViewingStage1RelaxBCs ? "Loading Step BCs →" : "← Relaxation Step BCs") {
                                previewStage = isViewingStage1RelaxBCs ? .stage2Load : .stage1Relax
                            }
                            .buttonStyle(.bordered)
                            .controlSize(.mini)
                        }
                    }
                }
                .padding(.horizontal, 10)
                .padding(.vertical, 5)
                .background(BlueHour.abyss.opacity(0.88), in: Capsule())
                .overlay(Capsule().strokeBorder((isPeriodic || isViewingStage1RelaxBCs) ? Color.green.opacity(0.55) : BlueHour.ice.opacity(0.45), lineWidth: 1))
                .padding(10)
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
    
    // MARK: - Crisp Vector Arrow with Filled Head
    private func drawVectorArrow(
        ctx: inout GraphicsContext,
        start: CGPoint,
        end: CGPoint,
        color: Color,
        lineWidth: CGFloat = 2.2,
        headLength: CGFloat = 6.5
    ) {
        let dx = end.x - start.x
        let dy = end.y - start.y
        let len = hypot(dx, dy)
        guard len > 1e-4 else { return }
        let ux = dx / len
        let uy = dy / len
        let nx = -uy
        let ny = ux
        
        var shaft = Path()
        shaft.move(to: start)
        shaft.addLine(to: end)
        ctx.stroke(shaft, with: .color(color), lineWidth: lineWidth)
        
        var head = Path()
        let p1 = end
        let p2 = CGPoint(x: end.x - headLength * ux + (headLength * 0.45) * nx,
                         y: end.y - headLength * uy + (headLength * 0.45) * ny)
        let p3 = CGPoint(x: end.x - headLength * ux - (headLength * 0.45) * nx,
                         y: end.y - headLength * uy - (headLength * 0.45) * ny)
        head.move(to: p1)
        head.addLine(to: p2)
        head.addLine(to: p3)
        head.closeSubpath()
        ctx.fill(head, with: .color(color))
    }

    // MARK: - Structural Mechanics Fixed-Base Clamp Marker
    private func drawFixedClampSupport(
        ctx: inout GraphicsContext,
        at pt: CGPoint,
        tangent: CGPoint,
        normalOut: CGPoint,
        color: Color = Color(red: 1.0, green: 0.42, blue: 0.42)
    ) {
        let r: CGFloat = 3.5
        let nodeRect = CGRect(x: pt.x - r, y: pt.y - r, width: r * 2, height: r * 2)
        ctx.fill(Path(ellipseIn: nodeRect), with: .color(color))
        ctx.stroke(Path(ellipseIn: nodeRect), with: .color(BlueHour.abyss), lineWidth: 1.0)
        
        let barHalf: CGFloat = 6.5
        let b1 = CGPoint(x: pt.x - barHalf * tangent.x + 3.0 * normalOut.x,
                         y: pt.y - barHalf * tangent.y + 3.0 * normalOut.y)
        let b2 = CGPoint(x: pt.x + barHalf * tangent.x + 3.0 * normalOut.x,
                         y: pt.y + barHalf * tangent.y + 3.0 * normalOut.y)
        var bar = Path()
        bar.move(to: b1)
        bar.addLine(to: b2)
        ctx.stroke(bar, with: .color(color), lineWidth: 1.8)
        
        let hatchLen: CGFloat = 4.5
        let hatchSpacing: CGFloat = 4.5
        for s in [-hatchSpacing, 0.0, hatchSpacing] {
            let hp = CGPoint(x: pt.x + s * tangent.x + 3.0 * normalOut.x,
                             y: pt.y + s * tangent.y + 3.0 * normalOut.y)
            let he = CGPoint(x: hp.x + (hatchLen * 0.7) * normalOut.x - (hatchLen * 0.7) * tangent.x,
                             y: hp.y + (hatchLen * 0.7) * normalOut.y - (hatchLen * 0.7) * tangent.y)
            var hLine = Path()
            hLine.move(to: hp)
            hLine.addLine(to: he)
            ctx.stroke(hLine, with: .color(color.opacity(0.85)), lineWidth: 1.2)
        }
    }

    // MARK: - Structural Mechanics Roller Support Marker
    private func drawRollerSupport(
        ctx: inout GraphicsContext,
        at pt: CGPoint,
        freeDir: CGPoint,
        normalOut: CGPoint,
        color: Color = Color(red: 0.98, green: 0.72, blue: 0.24)
    ) {
        let r: CGFloat = 3.5
        let rollerCenter = CGPoint(x: pt.x + 4.5 * normalOut.x, y: pt.y + 4.5 * normalOut.y)
        let nodeRect = CGRect(x: rollerCenter.x - r, y: rollerCenter.y - r, width: r * 2, height: r * 2)
        ctx.stroke(Path(ellipseIn: nodeRect), with: .color(color), lineWidth: 1.6)
        
        let trackHalf: CGFloat = 7.0
        let t1 = CGPoint(x: rollerCenter.x - trackHalf * freeDir.x + 3.5 * normalOut.x,
                         y: rollerCenter.y - trackHalf * freeDir.y + 3.5 * normalOut.y)
        let t2 = CGPoint(x: rollerCenter.x + trackHalf * freeDir.x + 3.5 * normalOut.x,
                         y: rollerCenter.y + trackHalf * freeDir.y + 3.5 * normalOut.y)
        var track = Path()
        track.move(to: t1)
        track.addLine(to: t2)
        ctx.stroke(track, with: .color(color.opacity(0.85)), lineWidth: 1.4)
    }

    // MARK: - Perimeter Boundary Atoms (Drawn Even When Large Meshes Hide Interior Sites)
    private func drawPerimeterLatticeSitesPreview(
        ctx: inout GraphicsContext,
        toScreen: (Double, Double) -> CGPoint,
        boxW: CGFloat,
        boxH: CGFloat
    ) {
        let micro = config.geometryConfig.microstructure
        let nx = max(2, config.nx)
        let ny = max(2, config.ny)
        let r: CGFloat = 2.8
        
        let countX = min(nx, 80)
        let countY = min(ny, 60)
        
        var coords: [(Double, Double)] = []
        coords.reserveCapacity(countX * 2 + countY * 2)
        
        // Bottom edge: y = 0
        for i in 0..<countX {
            let x = Double(i) / Double(countX - 1)
            coords.append((x, 0.0))
        }
        // Top edge: y = 1
        for i in 0..<countX {
            let x = Double(i) / Double(countX - 1)
            coords.append((x, 1.0))
        }
        // Left edge: x = 0
        for j in 1..<(countY - 1) {
            let y = Double(j) / Double(countY - 1)
            coords.append((0.0, y))
        }
        // Right edge: x = 1
        for j in 1..<(countY - 1) {
            let y = Double(j) / Double(countY - 1)
            coords.append((1.0, y))
        }
        
        for (xHat, yHat) in coords {
            guard config.geometryConfig.isPointInsideSolid(xHat: xHat, yHat: yHat) else { continue }
            let role = (micro.mode == .periodicDislocations)
                ? 0
                : config.geometryConfig.classifyNodeRole(
                    xHat: xHat,
                    yHat: yHat,
                    nx: nx,
                    ny: ny,
                    overrideRules: activeCanvasBoundaryRules,
                    overrideIndenterEnabled: activeCanvasIndenterEnabled
                )
            let pt = toScreen(xHat, yHat)
            let rect = CGRect(x: pt.x - r, y: pt.y - r, width: r * 2, height: r * 2)
            let dotColor: Color
            switch role {
            case 1: dotColor = Color(red: 1.0, green: 0.42, blue: 0.42) // Fixed (red)
            case 2: dotColor = Color(red: 0.31, green: 0.85, blue: 0.78) // Driven (cyan)
            case 3, 4: dotColor = Color(red: 0.98, green: 0.72, blue: 0.24) // Roller (yellow)
            default: dotColor = BlueHour.ice.opacity(0.85) // Free / Interior (ice)
            }
            ctx.fill(Path(ellipseIn: rect), with: .color(dotColor))
            ctx.stroke(Path(ellipseIn: rect), with: .color(BlueHour.abyss.opacity(0.8)), lineWidth: 0.8)
        }
    }
    
    private func drawBoundaryArrows(
        ctx: inout GraphicsContext,
        originX: CGFloat,
        originY: CGFloat,
        boxW: CGFloat,
        boxH: CGFloat
    ) {
        // -------------------------------------------------------------
        // CASE 1: FULL 2D PERIODIC BOUNDARY CONDITIONS (PBC)
        // -------------------------------------------------------------
        if config.isEffectivePeriodicPBC {
            let shearSign: CGFloat = (config.mode == .negative) ? -1.0 : 1.0
            let pbcColor = Color(red: 0.45, green: 0.90, blue: 0.62)
            
            // Top Edge: Macroscopic shear vectors pointing right (+X) or left (-X)
            let topY = originY - 14
            for k in 1...6 {
                let t = CGFloat(k) / 7.0
                let x = originX + t * boxW
                let start = CGPoint(x: x - shearSign * 12, y: topY)
                let end = CGPoint(x: x + shearSign * 12, y: topY)
                drawVectorArrow(ctx: &ctx, start: start, end: end, color: pbcColor, lineWidth: 2.2, headLength: 6.5)
            }
            
            // Bottom Edge: Complementary macroscopic shear vectors
            let botY = originY + boxH + 14
            for k in 1...6 {
                let t = CGFloat(k) / 7.0
                let x = originX + t * boxW
                let start = CGPoint(x: x + shearSign * 12, y: botY)
                let end = CGPoint(x: x - shearSign * 12, y: botY)
                drawVectorArrow(ctx: &ctx, start: start, end: end, color: pbcColor.opacity(0.85), lineWidth: 1.8, headLength: 6.0)
            }
            
            // Left & Right Edges: Periodic Torus wrap dashed lines
            var leftWrap = Path()
            leftWrap.move(to: CGPoint(x: originX, y: originY))
            leftWrap.addLine(to: CGPoint(x: originX, y: originY + boxH))
            ctx.stroke(leftWrap, with: .color(pbcColor.opacity(0.7)), style: StrokeStyle(lineWidth: 1.8, dash: [6, 4]))
            
            var rightWrap = Path()
            rightWrap.move(to: CGPoint(x: originX + boxW, y: originY))
            rightWrap.addLine(to: CGPoint(x: originX + boxW, y: originY + boxH))
            ctx.stroke(rightWrap, with: .color(pbcColor.opacity(0.7)), style: StrokeStyle(lineWidth: 1.8, dash: [6, 4]))
            
            // Informative Labels on Boundaries
            ctx.draw(
                Text("Top: Macroscopic Shear \(shearSign > 0 ? "+γ_xy" : "-γ_xy")")
                    .font(.system(size: 9.5, weight: .bold, design: .monospaced))
                    .foregroundColor(pbcColor),
                at: CGPoint(x: originX + boxW * 0.5, y: originY - 26)
            )
            ctx.draw(
                Text("Bottom: Wrap to Top (-γ_xy)")
                    .font(.system(size: 9.5, weight: .bold, design: .monospaced))
                    .foregroundColor(pbcColor.opacity(0.8)),
                at: CGPoint(x: originX + boxW * 0.5, y: originY + boxH + 26)
            )
            return
        }

        // -------------------------------------------------------------
        // CASE 2: SURFACE BOUNDARY CONDITIONS (FINITE SPECIMEN)
        // -------------------------------------------------------------
        let drivenColor = Color(red: 0.31, green: 0.85, blue: 0.78)
        let fixedColor = Color(red: 1.0, green: 0.42, blue: 0.42)
        let rollerColor = Color(red: 0.98, green: 0.72, blue: 0.24)
        
        let hasOuterAffine = activeCanvasBoundaryRules.contains { $0.enabled && ($0.region == .outerPerimeter || $0.kind == .affine) }
        
        if hasOuterAffine {
            // =========================================================
            // AFFINE SHEAR BOX (4 WALLS): ALL 4 BOUNDARY WALLS DRIVEN
            // =========================================================
            // Top Wall: Full rightward shear displacement
            let topY = originY - 14
            for k in 1...6 {
                let t = CGFloat(k) / 7.0
                let x = originX + t * boxW
                let start = CGPoint(x: x - 12, y: topY)
                let end = CGPoint(x: x + 12, y: topY)
                drawVectorArrow(ctx: &ctx, start: start, end: end, color: drivenColor, lineWidth: 2.2, headLength: 6.5)
            }
            ctx.draw(
                Text("Top Wall Driven: +ux = +γ·Ly (Max Shear)")
                    .font(.system(size: 9.5, weight: .bold, design: .monospaced))
                    .foregroundColor(drivenColor),
                at: CGPoint(x: originX + boxW * 0.5, y: originY - 26)
            )
            
            // Bottom Wall: Fixed Base u = 0 (clamped support)
            for k in 1...7 {
                let t = CGFloat(k) / 8.0
                let x = originX + t * boxW
                let pt = CGPoint(x: x, y: originY + boxH)
                drawFixedClampSupport(ctx: &ctx, at: pt, tangent: CGPoint(x: 1, y: 0), normalOut: CGPoint(x: 0, y: 1))
            }
            ctx.draw(
                Text("Bottom Wall: Fixed Base u = 0 (Clamped)")
                    .font(.system(size: 9.5, weight: .bold, design: .monospaced))
                    .foregroundColor(fixedColor),
                at: CGPoint(x: originX + boxW * 0.5, y: originY + boxH + 22)
            )
            
            // Left Wall: Affine shear gradient ux(y) = γ·y
            for k in 1...5 {
                let t = CGFloat(k) / 6.0
                let y = originY + boxH * (1.0 - t)
                let arrowLen = 5.0 + 17.0 * t
                let end = CGPoint(x: originX - 4, y: y)
                let start = CGPoint(x: originX - 4 - arrowLen, y: y)
                drawVectorArrow(ctx: &ctx, start: start, end: end, color: drivenColor, lineWidth: 1.8, headLength: 5.5)
            }
            ctx.draw(
                Text("Left: ux(y) = γ·y")
                    .font(.system(size: 8.5, weight: .bold, design: .monospaced))
                    .foregroundColor(drivenColor),
                at: CGPoint(x: originX - 35, y: originY + boxH * 0.5)
            )
            
            // Right Wall: Affine shear gradient ux(y) = γ·y
            for k in 1...5 {
                let t = CGFloat(k) / 6.0
                let y = originY + boxH * (1.0 - t)
                let arrowLen = 5.0 + 17.0 * t
                let start = CGPoint(x: originX + boxW + 4, y: y)
                let end = CGPoint(x: originX + boxW + 4 + arrowLen, y: y)
                drawVectorArrow(ctx: &ctx, start: start, end: end, color: drivenColor, lineWidth: 1.8, headLength: 5.5)
            }
            ctx.draw(
                Text("Right: ux(y) = γ·y")
                    .font(.system(size: 8.5, weight: .bold, design: .monospaced))
                    .foregroundColor(drivenColor),
                at: CGPoint(x: originX + boxW + 35, y: originY + boxH * 0.5)
            )
            
            // Outer Glowing Frame
            let frameRect = CGRect(x: originX, y: originY, width: boxW, height: boxH)
            ctx.stroke(Path(frameRect), with: .color(drivenColor.opacity(0.85)), style: StrokeStyle(lineWidth: 1.8, dash: [6, 4]))
            return
        }

        // =============================================================
        // REGIONAL RULES (TOP, BOTTOM, LEFT, RIGHT, CUSTOM BOX)
        // =============================================================
        for rule in activeCanvasBoundaryRules where rule.enabled {
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

            if rule.kind == .fixed {
                let count = 7
                for k in 1..<count {
                    let t = CGFloat(k) / CGFloat(count)
                    var pt = CGPoint.zero
                    var tang = CGPoint(x: 1, y: 0)
                    var nOut = CGPoint(x: 0, y: 1)
                    if rule.region == .bottom {
                        pt = CGPoint(x: originX + t * boxW, y: originY + boxH)
                        tang = CGPoint(x: 1, y: 0)
                        nOut = CGPoint(x: 0, y: 1)
                    } else if rule.region == .top {
                        pt = CGPoint(x: originX + t * boxW, y: originY)
                        tang = CGPoint(x: 1, y: 0)
                        nOut = CGPoint(x: 0, y: -1)
                    } else if rule.region == .left {
                        pt = CGPoint(x: originX, y: originY + t * boxH)
                        tang = CGPoint(x: 0, y: 1)
                        nOut = CGPoint(x: -1, y: 0)
                    } else if rule.region == .right {
                        pt = CGPoint(x: originX + boxW, y: originY + t * boxH)
                        tang = CGPoint(x: 0, y: 1)
                        nOut = CGPoint(x: 1, y: 0)
                    } else { continue }
                    drawFixedClampSupport(ctx: &ctx, at: pt, tangent: tang, normalOut: nOut, color: fixedColor)
                }
                if rule.region == .bottom {
                    ctx.draw(
                        Text("Bottom Fixed Base: u = 0 (Clamped)")
                            .font(.system(size: 9.5, weight: .bold, design: .monospaced))
                            .foregroundColor(fixedColor),
                        at: CGPoint(x: originX + boxW * 0.5, y: originY + boxH + 22)
                    )
                }
            } else if rule.kind == .rollerX || rule.kind == .rollerY {
                let count = 6
                for k in 1...count {
                    let t = CGFloat(k) / CGFloat(count + 1)
                    var pt = CGPoint.zero
                    var freeDir = CGPoint(x: 1, y: 0)
                    var nOut = CGPoint(x: 0, y: 1)
                    if rule.region == .bottom {
                        pt = CGPoint(x: originX + t * boxW, y: originY + boxH)
                        freeDir = CGPoint(x: 1, y: 0)
                        nOut = CGPoint(x: 0, y: 1)
                    } else if rule.region == .top {
                        pt = CGPoint(x: originX + t * boxW, y: originY)
                        freeDir = CGPoint(x: 1, y: 0)
                        nOut = CGPoint(x: 0, y: -1)
                    } else if rule.region == .left {
                        pt = CGPoint(x: originX, y: originY + t * boxH)
                        freeDir = CGPoint(x: 0, y: 1)
                        nOut = CGPoint(x: -1, y: 0)
                    } else if rule.region == .right {
                        pt = CGPoint(x: originX + boxW, y: originY + t * boxH)
                        freeDir = CGPoint(x: 0, y: 1)
                        nOut = CGPoint(x: 1, y: 0)
                    } else { continue }
                    drawRollerSupport(ctx: &ctx, at: pt, freeDir: freeDir, normalOut: nOut, color: rollerColor)
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
                let dx = CGFloat(vx / mag) * 20.0
                let dy = CGFloat(-vy / mag) * 20.0
                
                let count = 6
                for k in 1...count {
                    let t = CGFloat(k) / CGFloat(count + 1)
                    var edgePt = CGPoint.zero
                    if rule.region == .top {
                        edgePt = CGPoint(x: originX + t * boxW, y: originY - 4)
                    } else if rule.region == .bottom {
                        edgePt = CGPoint(x: originX + t * boxW, y: originY + boxH + 4)
                    } else if rule.region == .right {
                        edgePt = CGPoint(x: originX + boxW + 4, y: originY + t * boxH)
                    } else if rule.region == .left {
                        edgePt = CGPoint(x: originX - 4, y: originY + t * boxH)
                    } else {
                        continue
                    }
                    let start = edgePt
                    let end = CGPoint(x: edgePt.x + dx, y: edgePt.y + dy)
                    drawVectorArrow(ctx: &ctx, start: start, end: end, color: drivenColor, lineWidth: 2.2, headLength: 6.5)
                }
                
                let modeNoun: String
                if abs(vx) > 0.5 && abs(vy) < 0.5 {
                    modeNoun = (vx > 0 ? "Driven Grip: +ux (Shear)" : "Driven Grip: -ux (Reverse Shear)")
                } else if abs(vy) > 0.5 && abs(vx) < 0.5 {
                    modeNoun = (vy > 0 ? "Driven Grip: +uy (Tension)" : "Driven Grip: -uy (Compression)")
                } else {
                    modeNoun = String(format: "Driven Grip: (%.1f, %.1f)", vx, vy)
                }
                if rule.region == .top {
                    ctx.draw(
                        Text("Top \(modeNoun)")
                            .font(.system(size: 9.5, weight: .bold, design: .monospaced))
                            .foregroundColor(drivenColor),
                        at: CGPoint(x: originX + boxW * 0.5, y: originY - 26)
                    )
                } else if rule.region == .right {
                    ctx.draw(
                        Text("Right \(modeNoun)")
                            .font(.system(size: 9.5, weight: .bold, design: .monospaced))
                            .foregroundColor(drivenColor),
                        at: CGPoint(x: originX + boxW + 35, y: originY + boxH * 0.5)
                    )
                }
            } else if rule.kind == .free {
                if rule.region == .left {
                    var line = Path()
                    line.move(to: CGPoint(x: originX, y: originY))
                    line.addLine(to: CGPoint(x: originX, y: originY + boxH))
                    ctx.stroke(line, with: .color(BlueHour.ice.opacity(0.35)), style: StrokeStyle(lineWidth: 1.2, dash: [4, 4]))
                    ctx.draw(
                        Text("Free Surface")
                            .font(.system(size: 8.5, weight: .medium, design: .monospaced))
                            .foregroundColor(BlueHour.mist.opacity(0.85)),
                        at: CGPoint(x: originX - 35, y: originY + boxH * 0.5)
                    )
                } else if rule.region == .right {
                    var line = Path()
                    line.move(to: CGPoint(x: originX + boxW, y: originY))
                    line.addLine(to: CGPoint(x: originX + boxW, y: originY + boxH))
                    ctx.stroke(line, with: .color(BlueHour.ice.opacity(0.35)), style: StrokeStyle(lineWidth: 1.2, dash: [4, 4]))
                    ctx.draw(
                        Text("Free Surface")
                            .font(.system(size: 8.5, weight: .medium, design: .monospaced))
                            .foregroundColor(BlueHour.mist.opacity(0.85)),
                        at: CGPoint(x: originX + boxW + 35, y: originY + boxH * 0.5)
                    )
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
            // Hit-test Bicrystal GB center handle, Voronoi grain seeds, or Periodic Dislocation cores
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
            } else if micro.mode == .periodicDislocations {
                let ny = max(4, config.ny)
                for (dIdx, d) in micro.dislocations.enumerated() where d.enabled {
                    let yFrac = (Double(d.glide_plane_row) + 0.5) / Double(ny - 1)
                    let dPt = toScreen(d.x_frac, yFrac)
                    if hypot(startLocation.x - dPt.x, startLocation.y - dPt.y) < 18 {
                        activeDragMode = .moveDislocationCore(dIdx)
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
        case .moveDislocationCore(let dIdx):
            if dIdx >= 0 && dIdx < config.geometryConfig.microstructure.dislocations.count {
                config.geometryConfig.microstructure.dislocations[dIdx].x_frac = min(0.95, max(0.05, round(xh * 100) / 100.0))
                let ny = max(4, config.ny)
                let targetRow = Int(round(yh * Double(ny - 1) - 0.5))
                config.geometryConfig.microstructure.dislocations[dIdx].glide_plane_row = min(ny - 3, max(2, targetRow))
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
        GroupBox(label: Label("Discrete Crystal & Rotation (Fills Box)", systemImage: "atom")) {
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
                                selectableChipButton(
                                    title: "\(Int(deg))°",
                                    isSelected: abs(config.rotationAngle - deg) < 0.05
                                ) {
                                    config.rotationAngle = deg
                                }
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
        GroupBox(label: HStack {
            Label("Grains & Dislocations", systemImage: "square.split.diagonal.2x2")
            Spacer()
            HelpPopoverButton(
                title: "Microstructure, Grain Boundaries & Periodic Dislocations",
                summary: "Choose between:\n• Single Crystal: Defect-free crystal block.\n• Bicrystal (GB): Symmetric tilt or general grain boundary with Σ CSL presets.\n• 2D Voronoi: Multi-grain polycrystal or Tri-Crystal Y-junction.\n• Periodic Dislocs: Multi-plane edge dislocation dipoles (up to 100 dislocations) on a 2D periodic torus.",
                formulaOrTip: "• For Periodic Dislocs: Net Burgers vector Σ b_k = 0 is required for 2D PBC.\n• Keep vertical plane spacing Δr ≥ 5 rows so cores do not annihilate at α = 0.",
                docTopic: config.geometryConfig.microstructure.mode == .periodicDislocations ? .dislocationsGuide : .grainBoundaries,
                onOpenDocs: onOpenDocs
            )
        }) {
            VStack(alignment: .leading, spacing: 8) {
                Picker("", selection: $config.geometryConfig.microstructure.mode) {
                    Text("Single Crystal").tag(MicrostructureModeModel.singleCrystal)
                    Text("Bicrystal (GB)").tag(MicrostructureModeModel.bicrystal)
                    Text("2D Voronoi").tag(MicrostructureModeModel.voronoi)
                    Text("Periodic Dislocs").tag(MicrostructureModeModel.periodicDislocations)
                }
                .labelsHidden()
                .pickerStyle(.segmented)
                .controlSize(.small)
                .frame(maxWidth: .infinity)
                .onChange(of: config.geometryConfig.microstructure.mode) { newMode in
                    if newMode != .singleCrystal {
                        config.geometryConfig.enabled = true
                        if newMode == .bicrystal {
                            config.boundaryCondition = .fixed
                            // Ensure default loading is Simple Shear, never accidental Tension Y
                            if config.activeMechanicalLoading == .uniaxialTensionY {
                                config.applyMechanicalLoading(.simpleShearPos)
                            }
                        } else if newMode == .voronoi && config.geometryConfig.microstructure.seeds.isEmpty {
                            config.geometryConfig.microstructure.regenerateSeeds()
                        } else if newMode == .periodicDislocations {
                            let hasOutOfBounds = config.geometryConfig.microstructure.dislocations.contains {
                                $0.glide_plane_row < 1 || $0.glide_plane_row >= config.ny - 1
                            }
                            if config.geometryConfig.microstructure.dislocations.isEmpty || hasOutOfBounds {
                                config.geometryConfig.microstructure.regeneratePeriodicDislocations(nx: config.nx, ny: config.ny)
                            }
                            config.boundaryCondition = .periodic
                            config.noiseLevel = 0.0
                        }
                    }
                }
                .onChange(of: config.ny) { newNy in
                    if config.geometryConfig.microstructure.mode == .periodicDislocations &&
                        config.geometryConfig.microstructure.dislocation_layout != DislocationLayoutPreset.custom.rawValue {
                        config.geometryConfig.microstructure.regeneratePeriodicDislocations(nx: config.nx, ny: newNy)
                    }
                }
                
                if config.geometryConfig.microstructure.mode == .bicrystal {
                    let phiGB = config.geometryConfig.microstructure.gb_angle_deg
                    let phiRad = phiGB * .pi / 180.0
                    let isMostlyVertical = abs(cos(phiRad)) >= abs(sin(phiRad))
                    let currentHalfTilt = max(0.5, 0.5 * abs(config.geometryConfig.microstructure.theta2_deg - config.geometryConfig.microstructure.theta1_deg))
                    let cslMatch = CSLRationalizer.closestMatch(forHalfAngleDeg: currentHalfTilt, maxIndex: 18)
                    let isExactCSL = abs(abs(config.geometryConfig.microstructure.theta1_deg - phiGB) - cslMatch.halfAngleDeg) < 0.025 &&
                                     abs(abs(config.geometryConfig.microstructure.theta2_deg - phiGB) - cslMatch.halfAngleDeg) < 0.025

                    let isVerticalGBSelected = abs(phiGB) < 0.5 || abs(abs(phiGB) - 180.0) < 0.5
                    let isHorizontalGBSelected = abs(abs(phiGB) - 90.0) < 0.5

                    VStack(alignment: .leading, spacing: 6) {
                        HStack(spacing: 5) {
                            selectableChipButton(
                                title: "Vertical GB (ϕ=0°)",
                                isSelected: isVerticalGBSelected
                            ) {
                                config.geometryConfig.microstructure.gb_cx = 0.50
                                config.geometryConfig.microstructure.gb_cy = 0.50
                                config.geometryConfig.rotateBicrystalOrientation(toAngleDeg: 0.0)
                            }
                            
                            selectableChipButton(
                                title: "Horizontal GB (ϕ=90°)",
                                isSelected: isHorizontalGBSelected
                            ) {
                                config.geometryConfig.microstructure.gb_cx = 0.50
                                config.geometryConfig.microstructure.gb_cy = 0.50
                                config.geometryConfig.rotateBicrystalOrientation(toAngleDeg: 90.0)
                            }
                            
                            Menu {
                                ForEach(CSLRationalizer.curatedCatalog) { item in
                                    Button(item.shortLabel) {
                                        applyCSLMatch(item, applyAnalyticGap: autoUpdateGBSeamParams)
                                    }
                                }
                            } label: {
                                HStack(spacing: 3) {
                                    Image(systemName: "list.bullet.rectangle")
                                        .font(.system(size: 9))
                                    Text("CSL Catalog")
                                        .font(.system(size: 9.5, weight: .semibold))
                                }
                                .foregroundColor(Color(red: 1.0, green: 0.85, blue: 0.40))
                                .padding(.horizontal, 6)
                                .padding(.vertical, 3)
                            }
                            .menuStyle(.borderlessButton)
                            .background(BlueHour.midnight.opacity(0.90), in: RoundedRectangle(cornerRadius: 5))
                            .overlay(
                                RoundedRectangle(cornerRadius: 5)
                                    .strokeBorder(Color(red: 1.0, green: 0.82, blue: 0.35).opacity(0.75), lineWidth: 1.0)
                            )
                        }
                        
                        HStack(spacing: 4) {
                            Text("Quick CSL:")
                                .font(.system(size: 9.5, weight: .semibold))
                                .foregroundColor(BlueHour.mist)
                            selectableChipButton(
                                title: "Σ5 (310) ±18.43°",
                                isSelected: isExactCSL && cslMatch.m == 3 && cslMatch.n == 1
                            ) {
                                applyCSLMatch(CSLRationalizer.makeMatch(m: 3, n: 1), applyAnalyticGap: autoUpdateGBSeamParams)
                            }
                            
                            selectableChipButton(
                                title: "Σ5 (210) ±26.57°",
                                isSelected: isExactCSL && cslMatch.m == 2 && cslMatch.n == 1
                            ) {
                                applyCSLMatch(CSLRationalizer.makeMatch(m: 2, n: 1), applyAnalyticGap: autoUpdateGBSeamParams)
                            }

                            selectableChipButton(
                                title: "Σ13 (320) ±33.69°",
                                isSelected: isExactCSL && cslMatch.m == 3 && cslMatch.n == 2
                            ) {
                                applyCSLMatch(CSLRationalizer.makeMatch(m: 3, n: 2), applyAnalyticGap: autoUpdateGBSeamParams)
                            }
                        }
                        
                        // General (m, n) CSL Rationalizer Box (gb_code / aimsgb)
                        VStack(alignment: .leading, spacing: 5) {
                            HStack(spacing: 4) {
                                Image(systemName: isExactCSL ? "checkmark.seal.fill" : "wand.and.stars")
                                    .foregroundColor(isExactCSL ? Color(red: 0.45, green: 0.90, blue: 0.62) : Color(red: 0.31, green: 0.85, blue: 0.78))
                                    .font(.system(size: 10))
                                Text("CSL Rationalizer (\(cslMatch.m),\(cslMatch.n)): Σ\(cslMatch.sigma) (\(cslMatch.m) \(cslMatch.n) 0)")
                                    .font(.system(size: 10, weight: .bold))
                                    .foregroundColor(BlueHour.ice)
                                Spacer()
                                Text(String(format: "±%.3f° (Δθ=%.2f°)", cslMatch.halfAngleDeg, cslMatch.fullMisoriDeg))
                                    .font(.system(size: 9.5, weight: .semibold, design: .monospaced))
                                    .foregroundColor(Color(red: 1.0, green: 0.82, blue: 0.35))
                            }
                            
                            Text(String(format: "L_CSL = √%d a₀ = %.3f a₀ • d_(%d%d) = %.3f a₀ • d_min = %.2f a₀ • δ_⊥* = %.2f a₀",
                                        cslMatch.m * cslMatch.m + cslMatch.n * cslMatch.n,
                                        cslMatch.periodA0,
                                        cslMatch.m, cslMatch.n,
                                        cslMatch.layerSpacingA0,
                                        cslMatch.recommendedDMinA0,
                                        cslMatch.recommendedGapA0))
                                .font(.system(size: 9, design: .monospaced))
                                .foregroundColor(.secondary)
                            
                            HStack(spacing: 4) {
                                Button(action: {
                                    applyCSLMatch(cslMatch, applyAnalyticGap: autoUpdateGBSeamParams)
                                }) {
                                    Label(isExactCSL ? "Re-Sync (\(cslMatch.m),\(cslMatch.n)) Box & BCs" : "Snap to Exact (\(cslMatch.m),\(cslMatch.n)) CSL + Box & BCs",
                                          systemImage: "scope")
                                        .font(.system(size: 9.5, weight: .semibold))
                                        .frame(maxWidth: .infinity)
                                }
                                .buttonStyle(.borderedProminent)
                                .tint(Color(red: 0.24, green: 0.68, blue: 0.64))
                                .controlSize(.mini)
                            }
                            
                            HStack(spacing: 4) {
                                Text("Commensurate K·L_CSL:")
                                    .font(.system(size: 9, weight: .semibold))
                                    .foregroundColor(BlueHour.mist)
                                ForEach(cslMatch.suggestedPeriodCounts(around: config.nx), id: \.self) { k in
                                    let nBox = cslMatch.commensurateSize(periods: k)
                                    selectableChipButton(
                                        title: "\(nBox)×\(nBox) (\(k)L)",
                                        isSelected: config.nx == nBox && config.ny == nBox
                                    ) {
                                        config.nx = nBox
                                        config.ny = nBox
                                    }
                                }
                                if config.noiseLevel > 0.0 {
                                    Button("σ=0") {
                                        config.noiseLevel = 0.0
                                    }
                                    .buttonStyle(.borderedProminent)
                                    .tint(.orange)
                                    .controlSize(.mini)
                                }
                            }
                        }
                        .padding(6)
                        .background(BlueHour.slate.opacity(0.35), in: RoundedRectangle(cornerRadius: 6))
                        
                        // Applied Loading Quick-View & Switcher for Bicrystal
                        HStack(spacing: 6) {
                            Circle()
                                .fill(Color(red: 0.31, green: 0.85, blue: 0.78))
                                .frame(width: 7, height: 7)
                            Text("Applied Loading:")
                                .font(.system(size: 9.5, weight: .bold))
                                .foregroundColor(BlueHour.mist)
                            Text(config.activePhysicsSummaryBadge)
                                .font(.system(size: 10, weight: .bold, design: .monospaced))
                                .foregroundColor(BlueHour.ice)
                            Spacer()
                            Menu("Change Loading ▾") {
                                ForEach(ActiveMechanicalLoading.primaryGridModes) { load in
                                    Button(action: { config.applyMechanicalLoading(load) }) {
                                        if config.activeMechanicalLoading == load {
                                            Label(load.shortTitle, systemImage: "checkmark")
                                        } else {
                                            Label(load.shortTitle, systemImage: load.sfSymbol)
                                        }
                                    }
                                }
                            }
                            .font(.system(size: 9, weight: .semibold))
                            .menuStyle(.borderlessButton)
                        }
                        .padding(.horizontal, 8)
                        .padding(.vertical, 4)
                        .background(BlueHour.midnight.opacity(0.85), in: RoundedRectangle(cornerRadius: 6))
                        .overlay(RoundedRectangle(cornerRadius: 6).strokeBorder(BlueHour.border.opacity(0.6), lineWidth: 0.8))
                        
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
                            Slider(value: Binding(
                                get: { config.geometryConfig.microstructure.gb_angle_deg },
                                set: { newPhi in config.geometryConfig.rotateBicrystalOrientation(toAngleDeg: newPhi) }
                            ), in: -90.0...90.0, step: 1.0)
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
                        .onChange(of: config.geometryConfig.microstructure.theta1_deg) { _ in
                            syncAutoSeamParamsFromCurrentAngles()
                        }
                        
                        HStack {
                            Text("Grain 2 Angle θ₂ (°):")
                                .font(.caption2.bold())
                            Spacer()
                            Slider(value: $config.geometryConfig.microstructure.theta2_deg, in: -90.0...90.0, step: 0.5)
                                .frame(width: 85)
                            DoubleField(title: "+15.0", value: $config.geometryConfig.microstructure.theta2_deg, width: 48, minValue: -180.0, maxValue: 180.0)
                        }
                        .onChange(of: config.geometryConfig.microstructure.theta2_deg) { _ in
                            syncAutoSeamParamsFromCurrentAngles()
                        }
                        
                        Divider().padding(.vertical, 1)
                        
                        HStack {
                            Toggle(isOn: $autoUpdateGBSeamParams) {
                                Text("Auto-set δ_⊥, δ_∥ & d_min when angles change (editable below)")
                                    .font(.system(size: 9.5, weight: .semibold))
                                    .foregroundColor(Color(red: 0.45, green: 0.90, blue: 0.62))
                            }
                            .toggleStyle(.checkbox)
                            .controlSize(.mini)
                            Spacer()
                        }
                        
                        HStack(spacing: 4) {
                            Text("Grain Separation δ_⊥ (in a₀):")
                                .font(.caption2.bold())
                                .foregroundColor(Color(red: 0.31, green: 0.85, blue: 0.78))
                            Spacer()
                            Slider(value: $config.geometryConfig.microstructure.gb_normal_gap_a0, in: -0.40...1.00, step: 0.02)
                                .frame(width: 78)
                            DoubleField(title: "0.0", value: $config.geometryConfig.microstructure.gb_normal_gap_a0, width: 44, minValue: -0.60, maxValue: 1.50)
                            Button("Auto") {
                                config.geometryConfig.microstructure.gb_normal_gap_a0 = cslMatch.recommendedGapA0
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
                } else if config.geometryConfig.microstructure.mode == .periodicDislocations {
                    periodicDislocationsSection
                }
                
                if config.geometryConfig.microstructure.mode == .bicrystal || config.geometryConfig.microstructure.mode == .voronoi {
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
    
    // MARK: - Periodic Dislocations Section (2M Planes, 1 Dislocation Per Plane, Net b = 0)
    private var periodicDislocationsSection: some View {
        let micro = config.geometryConfig.microstructure
        let dislocations = micro.dislocations
        let totalDislocs = dislocations.filter(\.enabled).count
        let sumB = dislocations.filter(\.enabled).reduce(0) { $0 + $1.burgers_sign }
        let posCount = dislocations.filter { $0.enabled && $0.burgers_sign > 0 }.count
        let negCount = dislocations.filter { $0.enabled && $0.burgers_sign < 0 }.count
        let isBalanced = (sumB == 0) && (totalDislocs > 0)
        let isTri = (config.latticeSymmetry == .triangular)

        return VStack(alignment: .leading, spacing: 8) {
            // Pair Count & Presets (supports up to M = 150 pairs = 300 dislocations)
            HStack(spacing: 6) {
                Text("Pairs M:")
                    .font(.caption2.bold())
                IntField(
                    title: "30",
                    value: Binding(
                        get: { config.geometryConfig.microstructure.dislocation_pair_count },
                        set: { newM in
                            let clamped = max(1, min(150, newM))
                            config.geometryConfig.microstructure.dislocation_pair_count = clamped
                            config.geometryConfig.microstructure.regeneratePeriodicDislocations(nx: config.nx, ny: config.ny)
                        }
                    ),
                    width: 44
                )
                Stepper(
                    "= \(config.geometryConfig.microstructure.dislocation_pair_count * 2) dislocs",
                    value: $config.geometryConfig.microstructure.dislocation_pair_count,
                    in: 1...150,
                    onEditingChanged: { _ in
                        config.geometryConfig.microstructure.regeneratePeriodicDislocations(nx: config.nx, ny: config.ny)
                    }
                )
                .font(.caption2.monospacedDigit())
                Spacer()
                Button("Randomize") {
                    config.geometryConfig.microstructure.regeneratePeriodicDislocations(nx: config.nx, ny: config.ny)
                }
                .buttonStyle(.borderedProminent)
                .controlSize(.mini)
            }

            // Quick Pair Selection Chips (including M=10, M=15, M=30 (60D), M=50 (100D))
            ScrollView(.horizontal, showsIndicators: false) {
                HStack(spacing: 4) {
                    Text("Quick M:")
                        .font(.system(size: 9.5, weight: .semibold))
                        .foregroundColor(BlueHour.mist)
                    ForEach([1, 2, 3, 5, 10, 15, 30, 50], id: \.self) { mVal in
                        let isSel = (config.geometryConfig.microstructure.dislocation_pair_count == mVal)
                        selectableChipButton(
                            title: "M=\(mVal) (\(2 * mVal)D)",
                            isSelected: isSel
                        ) {
                            config.geometryConfig.microstructure.dislocation_pair_count = mVal
                            if config.ny < 4 * mVal {
                                let recN = max(config.nx, ((6 * mVal + 49) / 50) * 50)
                                config.setSize(nx: recN, ny: recN)
                            } else {
                                config.geometryConfig.microstructure.regeneratePeriodicDislocations(nx: config.nx, ny: config.ny)
                            }
                        }
                    }
                }
            }

            let targetPlanes = 2 * config.geometryConfig.microstructure.dislocation_pair_count
            if config.ny < 3 * targetPlanes {
                let recNy = ((4 * targetPlanes + 49) / 50) * 50
                HStack(spacing: 5) {
                    Image(systemName: "arrow.up.left.and.arrow.down.right")
                        .foregroundColor(.orange)
                    Text("\(targetPlanes) planes on Ny=\(config.ny) gives Δj ≈ \(String(format: "%.1f", Double(config.ny) / Double(targetPlanes))) rows.")
                        .font(.system(size: 9.5))
                        .foregroundColor(.orange)
                    Spacer()
                    Button("Set \(recNy)×\(recNy)") {
                        config.setSize(nx: recNy, ny: recNy)
                    }
                    .buttonStyle(.borderedProminent)
                    .tint(.orange)
                    .controlSize(.mini)
                }
                .padding(4)
                .background(Color.orange.opacity(0.12), in: RoundedRectangle(cornerRadius: 5))
            }

            // Layout Preset & Spacing
            HStack(spacing: 6) {
                Text("Layout:")
                    .font(.caption2.bold())
                Picker("", selection: $config.geometryConfig.microstructure.dislocation_layout) {
                    ForEach(DislocationLayoutPreset.allCases) { preset in
                        Text(preset.displayName).tag(preset.rawValue)
                    }
                }
                .labelsHidden()
                .controlSize(.mini)
                .frame(maxWidth: .infinity)
                .onChange(of: config.geometryConfig.microstructure.dislocation_layout) { _ in
                    config.geometryConfig.microstructure.regeneratePeriodicDislocations(nx: config.nx, ny: config.ny)
                }
            }

            // Vertical ±b Sign Order Quick Controls
            let currentSigns = config.geometryConfig.microstructure.dislocations.map(\.burgers_sign)
            let isAlternatingOrder = !currentSigns.isEmpty && currentSigns.enumerated().allSatisfy { $0.element == (($0.offset % 2 == 0) ? 1 : -1) }
            let isGroupedOrder = !currentSigns.isEmpty && currentSigns.enumerated().allSatisfy { $0.element == (($0.offset < currentSigns.count / 2) ? 1 : -1) }
            HStack(spacing: 4) {
                Text("±b Order in Y:")
                    .font(.system(size: 9.5, weight: .semibold))
                    .foregroundColor(BlueHour.mist)
                selectableChipButton(
                    title: "🎲 Shuffle ±b",
                    isSelected: !isAlternatingOrder && !isGroupedOrder && isBalanced
                ) {
                    shuffleBurgersSigns()
                }
                selectableChipButton(
                    title: "+ − Alternating",
                    isSelected: isAlternatingOrder
                ) {
                    setAlternatingBurgersSigns()
                }
                selectableChipButton(
                    title: "+ + − − Grouped",
                    isSelected: isGroupedOrder && currentSigns.count >= 4
                ) {
                    setGroupedBurgersSigns()
                }
            }

            HStack(spacing: 8) {
                Text("Min Spacing Δj:")
                    .font(.caption2)
                Stepper("\(config.geometryConfig.microstructure.min_plane_spacing_rows) rows",
                        value: $config.geometryConfig.microstructure.min_plane_spacing_rows,
                        in: 2...20)
                    .font(.caption2.monospacedDigit())
                Spacer()
                Text("Core r₀:")
                    .font(.caption2)
                DoubleField(title: "0.8", value: $config.geometryConfig.microstructure.dislocation_core_radius_a0, width: 42, minValue: 0.2, maxValue: 3.0)
                Text("a₀").font(.caption2)
            }

            // Topological Balance Indicator
            HStack(spacing: 6) {
                if isBalanced {
                    Image(systemName: "checkmark.seal.fill")
                        .foregroundColor(Color(red: 0.45, green: 0.90, blue: 0.62))
                    Text("Topologically Balanced: ∑b = 0 (\(posCount)×[+a₀] ⊥, \(negCount)×[-a₀] ⊤) · Full 2D PBC Valid")
                        .font(.system(size: 9.5, weight: .semibold))
                        .foregroundColor(Color(red: 0.45, green: 0.90, blue: 0.62))
                } else {
                    Image(systemName: "exclamationmark.triangle.fill")
                        .foregroundColor(.orange)
                    Text("Net ∑b = \(sumB > 0 ? "+\(sumB)" : "\(sumB)") a₀ ≠ 0! Periodic boundary conditions require ∑b = 0.")
                        .font(.system(size: 9.5, weight: .semibold))
                        .foregroundColor(.orange)
                    Spacer()
                    Button("Auto-Balance") {
                        balanceBurgersVector()
                    }
                    .buttonStyle(.bordered)
                    .controlSize(.mini)
                }
            }
            .padding(5)
            .background(BlueHour.midnight.opacity(0.55), in: RoundedRectangle(cornerRadius: 5))

            if isTri && (config.ny % 2 != 0) {
                HStack(spacing: 4) {
                    Image(systemName: "info.circle")
                        .foregroundColor(.orange)
                    Text("Triangular lattice vertical periodicity requires even Ny (Ny = \(config.ny) is odd). Auto-adjust to \(config.ny + 1)?")
                        .font(.system(size: 9.5))
                        .foregroundColor(.orange)
                    Spacer()
                    Button("Set Ny = \(config.ny + 1)") {
                        config.ny += 1
                        config.geometryConfig.microstructure.regeneratePeriodicDislocations(nx: config.nx, ny: config.ny)
                    }
                    .buttonStyle(.bordered)
                    .controlSize(.mini)
                }
                .padding(4)
                .background(Color.orange.opacity(0.12), in: RoundedRectangle(cornerRadius: 5))
            }

            // Per-Plane Dislocation Table
            Text("Single Dislocation per Glide Plane (Drag ⊥/⊤ on canvas or edit below):")
                .font(.system(size: 9.5))
                .foregroundColor(.secondary)

            ScrollView {
                VStack(spacing: 4) {
                    ForEach(config.geometryConfig.microstructure.dislocations.indices, id: \.self) { idx in
                        let disloc = config.geometryConfig.microstructure.dislocations[idx]
                        HStack(spacing: 5) {
                            Text("D\(idx)")
                                .font(.system(size: 10, weight: .bold, design: .monospaced))
                                .foregroundColor(disloc.burgers_sign > 0 ? Color(red: 0.35, green: 0.90, blue: 0.85) : Color(red: 1.0, green: 0.65, blue: 0.45))
                                .frame(width: 24, alignment: .leading)

                            Text("Row j:")
                                .font(.system(size: 9.5))
                            IntField(title: "20", value: $config.geometryConfig.microstructure.dislocations[idx].glide_plane_row, width: 38)

                            Text("x/L:")
                                .font(.system(size: 9.5))
                            DoubleField(title: "0.50", value: $config.geometryConfig.microstructure.dislocations[idx].x_frac, width: 44, minValue: 0.02, maxValue: 0.98)

                            Button(action: {
                                config.geometryConfig.microstructure.dislocations[idx].burgers_sign *= -1
                            }) {
                                Text(disloc.burgers_sign > 0 ? "⊥ +a₀" : "⊤ -a₀")
                                    .font(.system(size: 10, weight: .bold, design: .monospaced))
                                    .foregroundColor(disloc.burgers_sign > 0 ? Color(red: 0.20, green: 1.0, blue: 0.85) : Color(red: 1.0, green: 0.70, blue: 0.50))
                                    .padding(.horizontal, 6)
                                    .padding(.vertical, 2)
                                    .background(
                                        disloc.burgers_sign > 0
                                            ? Color(red: 0.12, green: 0.60, blue: 0.50).opacity(0.35)
                                            : Color(red: 0.80, green: 0.40, blue: 0.20).opacity(0.35),
                                        in: RoundedRectangle(cornerRadius: 4)
                                    )
                                    .overlay(
                                        RoundedRectangle(cornerRadius: 4)
                                            .strokeBorder(
                                                disloc.burgers_sign > 0
                                                    ? Color(red: 0.20, green: 1.0, blue: 0.85).opacity(0.6)
                                                    : Color(red: 1.0, green: 0.70, blue: 0.50).opacity(0.6),
                                                lineWidth: 1
                                            )
                                    )
                            }
                            .buttonStyle(.plain)

                            Spacer()

                            Toggle("", isOn: $config.geometryConfig.microstructure.dislocations[idx].enabled)
                                .labelsHidden()
                                .controlSize(.mini)

                            if config.geometryConfig.microstructure.dislocations.count > 2 {
                                Button(role: .destructive, action: {
                                    config.geometryConfig.microstructure.dislocations.remove(at: idx)
                                    config.geometryConfig.microstructure.dislocation_pair_count = max(1, config.geometryConfig.microstructure.dislocations.count / 2)
                                }) {
                                    Image(systemName: "trash")
                                        .font(.system(size: 9))
                                        .foregroundColor(.red.opacity(0.8))
                                }
                                .buttonStyle(.plain)
                            }
                        }
                        .padding(.horizontal, 4)
                        .padding(.vertical, 2)
                        .background(BlueHour.midnight.opacity(0.35), in: RoundedRectangle(cornerRadius: 4))
                    }
                }
            }
            .frame(maxHeight: 195)

            HStack {
                Button("+ Add Balanced (±b) Pair") {
                    addDislocationPair()
                }
                .buttonStyle(.bordered)
                .controlSize(.mini)
                
                Spacer()
                
                Text("\(totalDislocs) dislocations across \(totalDislocs) distinct glide planes")
                    .font(.system(size: 9))
                    .foregroundColor(BlueHour.mist)
            }
        }
        .frame(maxWidth: .infinity, alignment: .leading)
        .padding(6)
        .background(BlueHour.midnight.opacity(0.45), in: RoundedRectangle(cornerRadius: 6))
    }

    private func balanceBurgersVector() {
        var dislocs = config.geometryConfig.microstructure.dislocations
        let total = dislocs.count
        guard total >= 2 else { return }
        let half = total / 2
        for i in 0..<total {
            dislocs[i].burgers_sign = (i < half) ? 1 : -1
        }
        config.geometryConfig.microstructure.dislocations = dislocs
    }

    private func shuffleBurgersSigns() {
        var dislocs = config.geometryConfig.microstructure.dislocations
        let total = dislocs.count
        guard total >= 2 else { return }
        let half = total / 2
        let oldSigns = dislocs.map(\.burgers_sign)
        var signs = (0..<total).map { ($0 < half) ? 1 : -1 }
        for _ in 0..<16 {
            signs.shuffle()
            let isAlt = (0..<total).allSatisfy { signs[$0] == (($0 % 2 == 0) ? 1 : -1) }
            if (!isAlt || total < 4) && signs != oldSigns { break }
        }
        for i in 0..<total {
            dislocs[i].burgers_sign = signs[i]
        }
        config.geometryConfig.microstructure.dislocations = dislocs
    }

    private func setAlternatingBurgersSigns() {
        var dislocs = config.geometryConfig.microstructure.dislocations
        for i in 0..<dislocs.count {
            dislocs[i].burgers_sign = (i % 2 == 0) ? 1 : -1
        }
        config.geometryConfig.microstructure.dislocations = dislocs
    }

    private func setGroupedBurgersSigns() {
        var dislocs = config.geometryConfig.microstructure.dislocations
        let half = dislocs.count / 2
        for i in 0..<dislocs.count {
            dislocs[i].burgers_sign = (i < half) ? 1 : -1
        }
        config.geometryConfig.microstructure.dislocations = dislocs
    }

    private func addDislocationPair() {
        let ny = max(6, config.ny)
        var dislocs = config.geometryConfig.microstructure.dislocations
        let usedRows = Set(dislocs.map(\.glide_plane_row))
        var candidate1 = 3
        while candidate1 < ny - 3 && usedRows.contains(candidate1) {
            candidate1 += 1
        }
        var candidate2 = candidate1 + 4
        while candidate2 < ny - 3 && (usedRows.contains(candidate2) || candidate2 == candidate1) {
            candidate2 += 1
        }
        let nextID = (dislocs.map(\.id).max() ?? 0) + 1
        dislocs.append(SeededDislocationModel(id: nextID, glide_plane_row: candidate1, x_frac: 0.35, burgers_sign: 1, enabled: true))
        dislocs.append(SeededDislocationModel(id: nextID + 1, glide_plane_row: candidate2, x_frac: 0.65, burgers_sign: -1, enabled: true))
        config.geometryConfig.microstructure.dislocations = dislocs
        config.geometryConfig.microstructure.dislocation_pair_count = dislocs.count / 2
    }
    
    private struct CanvasLatticeSite {
        var xPhys: Double
        var yPhys: Double
        var gid: Int
        var thetaDeg: Double
        var isMergedSeam: Bool = false
    }

    private func syncAutoSeamParamsFromCurrentAngles() {
        guard autoUpdateGBSeamParams else { return }
        let halfTilt = max(0.5, 0.5 * abs(config.geometryConfig.microstructure.theta2_deg - config.geometryConfig.microstructure.theta1_deg))
        let match = CSLRationalizer.closestMatch(forHalfAngleDeg: halfTilt, maxIndex: 18)
        config.geometryConfig.microstructure.min_seam_distance_a0 = match.recommendedDMinA0
        config.geometryConfig.microstructure.gb_normal_gap_a0 = match.recommendedGapA0
        config.geometryConfig.microstructure.gb_tangential_shift_a0 = 0.0
    }

    private func applyCSLMatch(_ match: CSLRationalMatch, applyAnalyticGap: Bool) {
        let phiGB = config.geometryConfig.microstructure.gb_angle_deg
        config.geometryConfig.microstructure.gb_cx = 0.50
        config.geometryConfig.microstructure.gb_cy = 0.50
        config.geometryConfig.microstructure.theta1_deg = phiGB - match.halfAngleDeg
        config.geometryConfig.microstructure.theta2_deg = phiGB + match.halfAngleDeg
        config.geometryConfig.microstructure.seam_resolution = .mergeMidpoint
        config.geometryConfig.microstructure.min_seam_distance_a0 = match.recommendedDMinA0
        config.geometryConfig.microstructure.gb_normal_gap_a0 = applyAnalyticGap ? match.recommendedGapA0 : 0.0
        config.geometryConfig.microstructure.gb_tangential_shift_a0 = 0.0
        config.geometryConfig.microstructure.max_gb_edge_a0 = 2.60
        if config.protocolMode == .relaxOnly && !config.geometryConfig.indenter.enabled {
            config.geometryConfig.applyGBFrameRollerBCs(forAngleDeg: phiGB)
        } else {
            config.geometryConfig.primary_stress_component = "sigma_xy"
            config.geometryConfig.microstructure.stage2_load_type = "shear_parallel_gb"
            // When snapping to CSL or re-syncing, ensure loading is Simple Shear along boundary unless user chose indenter/affine
            if config.geometryConfig.boundary_rules.contains(where: { $0.uy > 0.5 && $0.region == .top && $0.kind == .prescribed }) ||
               config.activeMechanicalLoading == .uniaxialTensionY {
                config.geometryConfig.boundary_rules = [
                    BoundaryZoneRuleModel(region: .bottom, kind: .fixed, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                    BoundaryZoneRuleModel(region: .top, kind: .prescribed, thickness_a0: 2.0, ux: 1.0, uy: 0.0),
                    BoundaryZoneRuleModel(region: .left, kind: .free, thickness_a0: 1.5),
                    BoundaryZoneRuleModel(region: .right, kind: .free, thickness_a0: 1.5),
                    BoundaryZoneRuleModel(region: .holeSurfaces, kind: .free, thickness_a0: 1.5)
                ]
            }
        }
        config.noiseLevel = 0.0
        let kPeriodsX = max(2, Int(round(Double(max(config.nx, 32)) / match.periodA0)))
        let kPeriodsY = max(2, Int(round(Double(max(config.ny, 32)) / match.periodA0)))
        config.nx = match.commensurateSize(periods: kPeriodsX)
        config.ny = match.commensurateSize(periods: kPeriodsY)
    }

    private func snapNearCSLAngle(_ deg: Double) -> Double {
        let phiGB = config.geometryConfig.microstructure.gb_angle_deg
        let rel = deg - phiGB
        for cand in CSLRationalizer.allCandidates {
            if abs(abs(rel) - cand.halfAngleDeg) < 0.025 {
                return phiGB + ((rel < 0.0) ? -cand.halfAngleDeg : cand.halfAngleDeg)
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
        let stepStride = max(1, Int(ceil(Double(max(config.nx, config.ny)) / 350.0)))
        let isTri = (config.latticeSymmetry == .triangular)
        let dotR: CGFloat = max(0.85, min(4.2, min(boxW, boxH) / CGFloat(max(config.nx, config.ny) / stepStride) * 0.30))
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
        case .periodicDislocations:
            let ny = config.ny
            let nx = config.nx
            let hy = isTri ? (a0 * sqrt3 * 0.5) : a0
            let activeDislocs = micro.dislocations.filter(\.enabled)
            let rc2 = pow(max(0.4, micro.dislocation_core_radius_a0 * a0), 2)
            let nu = 0.30
            let pStride = max(1, Int(ceil(sqrt(Double(nx * ny) / 9000.0))))

            for j in stride(from: 0, to: ny, by: pStride) {
                var nStrip = 0
                for d in activeDislocs where d.glide_plane_row < j {
                    nStrip += d.burgers_sign
                }
                let Nj = max(2, nx + nStrip)
                let rowStagger = isTri ? (0.5 * Double(j % 2)) : 0.0
                let y0 = Double(j) * hy

                for i in stride(from: 0, to: Nj, by: pStride) {
                    var x0 = (Double(i) + 0.5 + rowStagger) * (Lx / Double(Nj))
                    x0 = fmod(x0, Lx)
                    if x0 < 0 { x0 += Lx }

                    var ux = 0.0
                    var uy = 0.0

                    for d in activeDislocs {
                        let yk = (Double(d.glide_plane_row) + 0.5) * hy
                        let xk = d.x_frac * Lx
                        let bk = Double(d.burgers_sign) * a0

                        var xi = x0 - xk
                        while xi < -0.5 * Lx { xi += Lx }
                        while xi >= 0.5 * Lx { xi -= Lx }
                        var eta = y0 - yk
                        while eta < -0.5 * Ly { eta += Ly }
                        while eta >= 0.5 * Ly { eta -= Ly }
                        let r2 = xi * xi + eta * eta + rc2

                        let sgnEta = eta >= 0 ? 1.0 : -1.0
                        let thetaK = atan2(eta, xi)
                        let thetaCentered = thetaK - sgnEta * (.pi / 2.0)
                        let termX = (bk / (2.0 * .pi)) * (thetaCentered + (xi * eta) / (2.0 * (1.0 - nu) * r2))
                        let termY = -(bk / (2.0 * .pi)) * (((1.0 - 2.0 * nu) / (4.0 * (1.0 - nu))) * log(r2) + (xi * xi - eta * eta) / (4.0 * (1.0 - nu) * r2))
                        let closureTilt = -sgnEta * (bk / (2.0 * Lx)) * xi

                        ux += termX + closureTilt
                        uy += termY
                    }

                    let xFinal = x0 + ux
                    var yFinal = fmod(y0 + uy, Ly)
                    if yFinal < 0 { yFinal += Ly }
                    let xHat = min(1.0, max(0.0, xFinal / Lx))
                    let yHat = min(1.0, max(0.0, yFinal / Ly))
                    if config.geometryConfig.isPointInsideSolid(xHat: xHat, yHat: yHat) {
                        rawSites.append(CanvasLatticeSite(xPhys: xHat * Lx, yPhys: yHat * Ly, gid: 0, thetaDeg: 0.0))
                    }
                }
            }
        }

        // Apply live Seam Overlap Filter using fast 2D spatial grid (O(N) on UI thread)
        var filteredSites = rawSites
        if (micro.mode == .bicrystal || micro.mode == .voronoi) && filteredSites.count > 1 {
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
            let role = (micro.mode == .periodicDislocations)
                ? 0
                : config.geometryConfig.classifyNodeRole(
                    xHat: xHat,
                    yHat: yHat,
                    nx: config.nx,
                    ny: config.ny,
                    overrideRules: activeCanvasBoundaryRules,
                    overrideIndenterEnabled: activeCanvasIndenterEnabled
                )
            let pt = toScreen(xHat, yHat)
            let r = site.isMergedSeam ? (dotR * 1.18) : dotR
            let rect = CGRect(x: pt.x - r, y: pt.y - r, width: r * 2, height: r * 2)
            let baseColor = (micro.mode == .singleCrystal || micro.mode == .periodicDislocations)
                ? BlueHour.ice.opacity(0.68)
                : grainPaletteColor(site.gid).opacity(0.82)
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
        } else if micro.mode == .periodicDislocations {
            let ny = max(4, config.ny)
            let h_y_norm = 1.0 / Double(ny - 1)
            let activeCount = micro.dislocations.filter(\.enabled).count
            let isDense = (activeCount > 14)
            
            for (dIdx, d) in micro.dislocations.enumerated() where d.enabled {
                let yFrac = (Double(d.glide_plane_row) + 0.5) * h_y_norm
                let pLeft = toScreen(0.0, yFrac)
                let pRight = toScreen(1.0, yFrac)
                let pCore = toScreen(d.x_frac, yFrac)
                
                let isPos = (d.burgers_sign > 0)
                let planeColor = isPos
                    ? Color(red: 0.20, green: 0.88, blue: 0.85)
                    : Color(red: 1.0, green: 0.62, blue: 0.42)
                
                // Horizontal dashed glide plane line
                var planePath = Path()
                planePath.move(to: pLeft)
                planePath.addLine(to: pRight)
                ctx.stroke(
                    planePath,
                    with: .color(planeColor.opacity(isDense ? 0.38 : 0.70)),
                    style: StrokeStyle(lineWidth: isDense ? 0.8 : 1.4, dash: [4, 3])
                )
                
                // Glide motion arrow under shear tau_xy > 0 (+b moves right, -b moves left)
                let arrowDir: CGFloat = isPos ? 1.0 : -1.0
                let arrowHalfLen: CGFloat = isDense ? 13 : 22
                let arrowStartX = pCore.x - arrowDir * arrowHalfLen
                let arrowEndX = pCore.x + arrowDir * arrowHalfLen
                var arrowPath = Path()
                arrowPath.move(to: CGPoint(x: arrowStartX, y: pCore.y))
                arrowPath.addLine(to: CGPoint(x: arrowEndX, y: pCore.y))
                arrowPath.addLine(to: CGPoint(x: arrowEndX - arrowDir * 4, y: pCore.y - 3.0))
                arrowPath.move(to: CGPoint(x: arrowEndX, y: pCore.y))
                arrowPath.addLine(to: CGPoint(x: arrowEndX - arrowDir * 4, y: pCore.y + 3.0))
                ctx.stroke(arrowPath, with: .color(planeColor.opacity(isDense ? 0.65 : 0.85)), lineWidth: isDense ? 0.9 : 1.2)
                
                // Core handle with ⊥ or ⊤ symbol
                let coreR: CGFloat = isDense ? 6.0 : 8.5
                let hRect = CGRect(x: pCore.x - coreR, y: pCore.y - coreR, width: 2 * coreR, height: 2 * coreR)
                ctx.fill(Path(ellipseIn: hRect), with: .color(BlueHour.midnight))
                ctx.stroke(Path(ellipseIn: hRect), with: .color(planeColor), lineWidth: isDense ? 1.3 : 1.8)
                
                ctx.draw(
                    Text(isPos ? "⊥" : "⊤")
                        .font(.system(size: isDense ? 8.5 : 11, weight: .heavy))
                        .foregroundColor(planeColor),
                    at: CGPoint(x: pCore.x, y: pCore.y)
                )
                
                if !isDense {
                    // Core tag pill (clamp y inside the canvas box so top/bottom planes never clip)
                    let clampedLabelY = max(originY + 12, min(originY + boxH - 12, pCore.y - 15))
                    let clampedLabelX = max(originX + 40, min(originX + boxW - 40, pCore.x))
                    let labelPt = CGPoint(x: clampedLabelX, y: clampedLabelY)
                    let pillW: CGFloat = 76
                    let pillH: CGFloat = 16
                    let pillRect = CGRect(x: labelPt.x - pillW * 0.5, y: labelPt.y - pillH * 0.5, width: pillW, height: pillH)
                    ctx.fill(Path(roundedRect: pillRect, cornerRadius: 4), with: .color(BlueHour.abyss.opacity(0.88)))
                    ctx.stroke(Path(roundedRect: pillRect, cornerRadius: 4), with: .color(planeColor.opacity(0.60)), lineWidth: 0.8)
                    ctx.draw(
                        Text("D\(dIdx) \(isPos ? "+a₀" : "-a₀") j=\(d.glide_plane_row)")
                            .font(.system(size: 8.5, weight: .bold, design: .monospaced))
                            .foregroundColor(planeColor),
                        at: labelPt
                    )
                }
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
    
    private var stageProtocolBannerCard: some View {
        GroupBox(label: Label("Select Step to Edit in BCs & Indenter", systemImage: "arrow.triangle.branch")) {
            VStack(alignment: .leading, spacing: 8) {
                if config.protocolMode == .relaxThenLoad {
                    HStack(spacing: 4) {
                        ForEach(BCStagePreviewMode.allCases) { mode in
                            let isSelected = previewStage == mode
                            let activeColor = mode == .stage1Relax
                                ? Color(red: 0.16, green: 0.72, blue: 0.52)
                                : Color(red: 0.18, green: 0.62, blue: 0.78)
                            let borderColor = mode == .stage1Relax
                                ? Color(red: 0.45, green: 0.96, blue: 0.78)
                                : Color(red: 0.48, green: 0.90, blue: 0.98)
                            Button {
                                previewStage = mode
                            } label: {
                                HStack(spacing: 4) {
                                    Image(systemName: isSelected ? "checkmark.circle.fill" : "circle")
                                        .font(.system(size: 9, weight: .bold))
                                    Text(mode.rawValue)
                                        .font(.system(size: 10.5, weight: isSelected ? .bold : .medium))
                                        .lineLimit(1)
                                }
                                .foregroundColor(isSelected ? .white : BlueHour.ice.opacity(0.88))
                                .frame(maxWidth: .infinity)
                                .padding(.vertical, 5)
                                .padding(.horizontal, 6)
                                .background(
                                    RoundedRectangle(cornerRadius: 6)
                                        .fill(isSelected ? activeColor : BlueHour.midnight.opacity(0.75))
                                )
                                .overlay(
                                    RoundedRectangle(cornerRadius: 6)
                                        .strokeBorder(
                                            isSelected ? borderColor : BlueHour.steel.opacity(0.65),
                                            lineWidth: isSelected ? 1.3 : 0.8
                                        )
                                )
                            }
                            .buttonStyle(.plain)
                        }
                    }
                    .padding(3)
                    .background(
                        RoundedRectangle(cornerRadius: 8)
                            .fill(BlueHour.midnight.opacity(0.92))
                    )
                    .overlay(
                        RoundedRectangle(cornerRadius: 8)
                            .strokeBorder(BlueHour.steel.opacity(0.55), lineWidth: 1)
                    )
                } else {
                    HStack {
                        Label("Protocol Mode: \(config.protocolMode.shortTitle)", systemImage: "checkmark.seal.fill")
                            .font(.caption.bold())
                            .foregroundColor(Color(red: 0.45, green: 0.92, blue: 0.62))
                        Spacer()
                        Button("Switch to 2. Relax → Load") {
                            config.protocolMode = .relaxThenLoad
                            config.geometryConfig.microstructure.relax_only = false
                            previewStage = .stage2Load
                        }
                        .buttonStyle(.borderedProminent)
                        .controlSize(.mini)
                    }
                }
                
                if isViewingStage1RelaxBCs {
                    HStack(alignment: .top, spacing: 6) {
                        Circle().fill(Color.green).frame(width: 7, height: 7).padding(.top, 3)
                        VStack(alignment: .leading, spacing: 2) {
                            Text("Editing Relaxation Step BCs (\(config.geometryConfig.stage1_boundary_rules.filter(\.enabled).count) active rules)")
                                .font(.caption.bold())
                                .foregroundColor(Color(red: 0.45, green: 0.92, blue: 0.62))
                            Text("Boundary conditions used during the initial relaxation step before mechanical loading.")
                                .font(.system(size: 9.5))
                                .foregroundColor(.secondary)
                        }
                    }
                    .padding(7)
                    .frame(maxWidth: .infinity, alignment: .leading)
                    .background(Color.green.opacity(0.12), in: RoundedRectangle(cornerRadius: 6))
                    .overlay(RoundedRectangle(cornerRadius: 6).strokeBorder(Color.green.opacity(0.55), lineWidth: 1))
                } else {
                    HStack(alignment: .top, spacing: 6) {
                        Circle().fill(Color(red: 0.31, green: 0.85, blue: 0.78)).frame(width: 7, height: 7).padding(.top, 3)
                        VStack(alignment: .leading, spacing: 2) {
                            Text("Editing Loading Step BCs & Indenter (\(config.geometryConfig.boundary_rules.filter(\.enabled).count) active rules)")
                                .font(.caption.bold())
                                .foregroundColor(BlueHour.ice)
                            Text("Boundary conditions and Nanoindenter used during the mechanical loading step.")
                                .font(.system(size: 9.5))
                                .foregroundColor(.secondary)
                        }
                    }
                    .padding(7)
                    .frame(maxWidth: .infinity, alignment: .leading)
                    .background(BlueHour.slate.opacity(0.45), in: RoundedRectangle(cornerRadius: 6))
                    .overlay(RoundedRectangle(cornerRadius: 6).strokeBorder(BlueHour.ice.opacity(0.55), lineWidth: 1))
                }
            }
            .padding(6)
        }
    }
    
    private var boundaryConditionsCard: some View {
        GroupBox(label: Label(
            isViewingStage1RelaxBCs
                ? "Relaxation Step BCs"
                : "Loading Step BCs",
            systemImage: isViewingStage1RelaxBCs ? "checkmark.seal.fill" : "arrow.up.and.down.and.sparkles"
        )) {
            VStack(alignment: .leading, spacing: 8) {
                // Quick BC Templates (Stage-aware)
                ScrollView(.horizontal, showsIndicators: false) {
                    HStack(spacing: 6) {
                        if isViewingStage1RelaxBCs {
                            let s1 = config.geometryConfig.stage1_boundary_rules
                            let s2 = config.geometryConfig.boundary_rules
                            let kindOf = { (reg: BoundaryRegionSelector) -> BoundaryConditionKindModel? in
                                s1.first(where: { $0.enabled && $0.region == reg })?.kind
                            }
                            let isSameAsLoading = !s1.isEmpty && s1.count == s2.count && zip(s1, s2).allSatisfy { a, b in
                                a.region == b.region && a.kind == b.kind && a.affine_mode == b.affine_mode
                            }
                            let isGBRollerDefault = !isSameAsLoading && ((kindOf(.top) == .rollerX && kindOf(.bottom) == .rollerX && kindOf(.left) == .rollerY && kindOf(.right) == .free)
                                || (kindOf(.left) == .rollerY && kindOf(.right) == .rollerY && kindOf(.bottom) == .rollerX && kindOf(.top) == .free))
                            let is4EdgeRoller = !isSameAsLoading && ((kindOf(.top) == .rollerX && kindOf(.bottom) == .rollerX && kindOf(.left) == .rollerX && kindOf(.right) == .rollerX)
                                || (kindOf(.top) == .rollerY && kindOf(.bottom) == .rollerY && kindOf(.left) == .rollerY && kindOf(.right) == .rollerY))
                            let isIndentBase = !isSameAsLoading && (kindOf(.bottom) == .fixed && kindOf(.left) == .rollerY && kindOf(.right) == .rollerY && kindOf(.top) == .free)
                            let isAllFree = !isSameAsLoading && (kindOf(.bottom) == .free && kindOf(.top) == .free && kindOf(.left) == .free && kindOf(.right) == .free)
                            
                            selectableChipButton(title: "GB Roller (Default)", isSelected: isGBRollerDefault) {
                                config.geometryConfig.resetStage1ToGBFrameRollerBCs(forAngleDeg: config.geometryConfig.microstructure.gb_angle_deg)
                            }
                            selectableChipButton(title: "Same as Loading (α = 0)", isSelected: isSameAsLoading) {
                                config.geometryConfig.stage1_boundary_rules = config.geometryConfig.boundary_rules
                            }
                            selectableChipButton(title: "4-Edge Roller (Auto X/Y)", isSelected: is4EdgeRoller) {
                                let phiRad = config.geometryConfig.microstructure.gb_angle_deg * .pi / 180.0
                                let isVert = abs(cos(phiRad)) >= abs(sin(phiRad))
                                let rk: BoundaryConditionKindModel = isVert ? .rollerX : .rollerY
                                config.geometryConfig.stage1_boundary_rules = [
                                    BoundaryZoneRuleModel(region: .bottom, kind: rk, thickness_a0: 2.0, lock_x: !isVert, lock_y: isVert),
                                    BoundaryZoneRuleModel(region: .top, kind: rk, thickness_a0: 2.0, lock_x: !isVert, lock_y: isVert),
                                    BoundaryZoneRuleModel(region: .left, kind: rk, thickness_a0: 2.0, lock_x: !isVert, lock_y: isVert),
                                    BoundaryZoneRuleModel(region: .right, kind: rk, thickness_a0: 2.0, lock_x: !isVert, lock_y: isVert)
                                ]
                            }
                            selectableChipButton(title: "Indent Base (Bottom Fixed)", isSelected: isIndentBase) {
                                config.geometryConfig.stage1_boundary_rules = [
                                    BoundaryZoneRuleModel(region: .bottom, kind: .fixed, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                                    BoundaryZoneRuleModel(region: .left, kind: .rollerY, thickness_a0: 2.0),
                                    BoundaryZoneRuleModel(region: .right, kind: .rollerY, thickness_a0: 2.0),
                                    BoundaryZoneRuleModel(region: .top, kind: .free, thickness_a0: 1.5)
                                ]
                            }
                            selectableChipButton(title: "All Free Surface", isSelected: isAllFree) {
                                config.geometryConfig.stage1_boundary_rules = [
                                    BoundaryZoneRuleModel(region: .bottom, kind: .free, thickness_a0: 1.5),
                                    BoundaryZoneRuleModel(region: .top, kind: .free, thickness_a0: 1.5),
                                    BoundaryZoneRuleModel(region: .left, kind: .free, thickness_a0: 1.5),
                                    BoundaryZoneRuleModel(region: .right, kind: .free, thickness_a0: 1.5)
                                ]
                            }
                        } else {
                            let s2 = config.geometryConfig.boundary_rules
                            let ruleOf = { (reg: BoundaryRegionSelector) -> BoundaryZoneRuleModel? in
                                s2.first(where: { $0.enabled && $0.region == reg })
                            }
                            let isIndentBase = ruleOf(.bottom)?.kind == .fixed && ruleOf(.left)?.kind == .rollerY && ruleOf(.right)?.kind == .rollerY && ruleOf(.top)?.kind == .free
                            let isShearGB = ruleOf(.left)?.kind == .fixed && ruleOf(.right)?.kind == .prescribed && abs((ruleOf(.right)?.uy ?? 0) - 1.0) < 0.01
                            let isTensionGB = ruleOf(.left)?.kind == .rollerY && ruleOf(.right)?.kind == .prescribed && abs((ruleOf(.right)?.ux ?? 0) - 1.0) < 0.01
                            let isTensionY = ruleOf(.bottom)?.kind == .fixed && ruleOf(.top)?.kind == .prescribed && abs((ruleOf(.top)?.uy ?? 0) - 1.0) < 0.01
                            let isCompY = ruleOf(.bottom)?.kind == .fixed && ruleOf(.top)?.kind == .prescribed && abs((ruleOf(.top)?.uy ?? 0) + 1.0) < 0.01
                            let isShearTop = ruleOf(.bottom)?.kind == .fixed && ruleOf(.top)?.kind == .prescribed && abs((ruleOf(.top)?.ux ?? 0) - 1.0) < 0.01
                            let isAffine = ruleOf(.outerPerimeter)?.kind == .affine
                            let isGBRoller2 = (ruleOf(.top)?.kind == .rollerX && ruleOf(.bottom)?.kind == .rollerX && ruleOf(.left)?.kind == .rollerY && ruleOf(.right)?.kind == .free)
                                || (ruleOf(.left)?.kind == .rollerY && ruleOf(.right)?.kind == .rollerY && ruleOf(.bottom)?.kind == .rollerX && ruleOf(.top)?.kind == .free)
                            
                            selectableChipButton(title: "Indent Base (Bottom Fixed)", isSelected: isIndentBase) {
                                config.geometryConfig.boundary_rules = [
                                    BoundaryZoneRuleModel(region: .bottom, kind: .fixed, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                                    BoundaryZoneRuleModel(region: .left, kind: .rollerY, thickness_a0: 2.0),
                                    BoundaryZoneRuleModel(region: .right, kind: .rollerY, thickness_a0: 2.0),
                                    BoundaryZoneRuleModel(region: .top, kind: .free, thickness_a0: 1.5),
                                    BoundaryZoneRuleModel(region: .holeSurfaces, kind: .free, thickness_a0: 1.5)
                                ]
                                config.geometryConfig.primary_stress_component = "indentation_load"
                            }
                            selectableChipButton(title: "Shear ∥ GB (Right +Y)", isSelected: isShearGB) {
                                config.geometryConfig.boundary_rules = [
                                    BoundaryZoneRuleModel(region: .left, kind: .fixed, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                                    BoundaryZoneRuleModel(region: .right, kind: .prescribed, thickness_a0: 2.0, ux: 0.0, uy: 1.0, lock_x: true, lock_y: true),
                                    BoundaryZoneRuleModel(region: .bottom, kind: .free, thickness_a0: 1.5),
                                    BoundaryZoneRuleModel(region: .top, kind: .free, thickness_a0: 1.5),
                                    BoundaryZoneRuleModel(region: .holeSurfaces, kind: .free, thickness_a0: 1.5)
                                ]
                                config.geometryConfig.primary_stress_component = "sigma_xy"
                            }
                            selectableChipButton(title: "Tension ⊥ GB (Right +X)", isSelected: isTensionGB) {
                                config.geometryConfig.boundary_rules = [
                                    BoundaryZoneRuleModel(region: .left, kind: .rollerY, thickness_a0: 2.0),
                                    BoundaryZoneRuleModel(region: .right, kind: .prescribed, thickness_a0: 2.0, ux: 1.0, uy: 0.0, lock_x: true, lock_y: false),
                                    BoundaryZoneRuleModel(region: .bottom, kind: .rollerX, thickness_a0: 2.0),
                                    BoundaryZoneRuleModel(region: .top, kind: .rollerX, thickness_a0: 2.0),
                                    BoundaryZoneRuleModel(region: .holeSurfaces, kind: .free, thickness_a0: 1.5)
                                ]
                                config.geometryConfig.primary_stress_component = "sigma_xx"
                            }
                            selectableChipButton(title: "Tension (+Y)", isSelected: isTensionY) {
                                config.geometryConfig.boundary_rules = [
                                    BoundaryZoneRuleModel(region: .bottom, kind: .fixed, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                                    BoundaryZoneRuleModel(region: .top, kind: .prescribed, thickness_a0: 2.0, ux: 0.0, uy: 1.0),
                                    BoundaryZoneRuleModel(region: .left, kind: .free, thickness_a0: 1.5),
                                    BoundaryZoneRuleModel(region: .right, kind: .free, thickness_a0: 1.5),
                                    BoundaryZoneRuleModel(region: .holeSurfaces, kind: .free, thickness_a0: 1.5)
                                ]
                                config.geometryConfig.primary_stress_component = "sigma_yy"
                            }
                            selectableChipButton(title: "Compression (-Y)", isSelected: isCompY) {
                                config.geometryConfig.boundary_rules = [
                                    BoundaryZoneRuleModel(region: .bottom, kind: .fixed, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                                    BoundaryZoneRuleModel(region: .top, kind: .prescribed, thickness_a0: 2.0, ux: 0.0, uy: -1.0),
                                    BoundaryZoneRuleModel(region: .left, kind: .free, thickness_a0: 1.5),
                                    BoundaryZoneRuleModel(region: .right, kind: .free, thickness_a0: 1.5),
                                    BoundaryZoneRuleModel(region: .holeSurfaces, kind: .free, thickness_a0: 1.5)
                                ]
                                config.geometryConfig.primary_stress_component = "sigma_yy"
                            }
                            selectableChipButton(title: "Shear Grips (Top +X)", isSelected: isShearTop) {
                                config.geometryConfig.boundary_rules = [
                                    BoundaryZoneRuleModel(region: .bottom, kind: .fixed, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                                    BoundaryZoneRuleModel(region: .top, kind: .prescribed, thickness_a0: 2.0, ux: 1.0, uy: 0.0),
                                    BoundaryZoneRuleModel(region: .holeSurfaces, kind: .free, thickness_a0: 1.5)
                                ]
                                config.geometryConfig.primary_stress_component = "sigma_xy"
                            }
                            selectableChipButton(title: "Affine Shear Box", isSelected: isAffine) {
                                config.geometryConfig.boundary_rules = [
                                    BoundaryZoneRuleModel(region: .outerPerimeter, kind: .affine, thickness_a0: 2.0, affine_mode: "simple_shear_x"),
                                    BoundaryZoneRuleModel(region: .holeSurfaces, kind: .free, thickness_a0: 1.5)
                                ]
                                config.geometryConfig.primary_stress_component = "sigma_xy"
                            }
                            selectableChipButton(title: "GB Roller (Auto X/Y)", isSelected: isGBRoller2) {
                                config.geometryConfig.applyGBFrameRollerBCs(forAngleDeg: config.geometryConfig.microstructure.gb_angle_deg)
                            }
                        }
                    }
                    .padding(.vertical, 2)
                }
                
                ForEach(activeEditingRulesBinding) { $rule in
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
                                if isViewingStage1RelaxBCs {
                                    config.geometryConfig.stage1_boundary_rules.removeAll { $0.id == idToRemove }
                                } else {
                                    config.geometryConfig.boundary_rules.removeAll { $0.id == idToRemove }
                                }
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
                    if isViewingStage1RelaxBCs {
                        config.geometryConfig.stage1_boundary_rules.append(
                            BoundaryZoneRuleModel(region: .top, kind: .rollerX, thickness_a0: 1.5, ux: 0.0, uy: 0.0, lock_x: false, lock_y: true)
                        )
                    } else {
                        config.geometryConfig.boundary_rules.append(
                            BoundaryZoneRuleModel(region: .top, kind: .prescribed, thickness_a0: 2.0, ux: 0.0, uy: 1.0)
                        )
                    }
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
        GroupBox(label: Label("Loading Step: Nanoindenter Contact Tool", systemImage: "arrow.down.to.line.compact")) {
            VStack(alignment: .leading, spacing: 8) {
                Toggle("Enable Nanoindenter Tool", isOn: $config.geometryConfig.indenter.enabled)
                    .font(.caption.bold())
                    .onChange(of: config.geometryConfig.indenter.enabled) { enabled in
                        if enabled {
                            config.geometryConfig.enabled = true
                            config.boundaryCondition = .fixed
                            config.geometryConfig.primary_stress_component = "indentation_load"
                            config.geometryConfig.boundary_rules = [
                                BoundaryZoneRuleModel(region: .bottom, kind: .fixed, thickness_a0: 2.0, ux: 0.0, uy: 0.0),
                                BoundaryZoneRuleModel(region: .left, kind: .rollerY, thickness_a0: 2.0),
                                BoundaryZoneRuleModel(region: .right, kind: .rollerY, thickness_a0: 2.0),
                                BoundaryZoneRuleModel(region: .top, kind: .free, thickness_a0: 1.5),
                                BoundaryZoneRuleModel(region: .holeSurfaces, kind: .free, thickness_a0: 1.5)
                            ]
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
