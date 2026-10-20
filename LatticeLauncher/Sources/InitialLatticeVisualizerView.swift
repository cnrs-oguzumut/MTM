import SwiftUI
import AppKit

public enum InitialMeshColorField: String, CaseIterable, Identifiable {
    case uniform = "Uniform / Mesh"
    case energyDensity = "Strain Energy (W_e)"
    case vonMises = "von Mises (σ_vM)"
    case shearStress = "Shear Stress (σ_xy)"
    case pressure = "Pressure (p)"
    case sigmaXX = "Normal σ_xx"
    case sigmaYY = "Normal σ_yy"
    
    public var id: String { rawValue }
}

public enum InitialMeshViewMode: String, CaseIterable, Identifiable {
    case primary = "Primary Cell (1×1)"
    case periodic3x3 = "Periodic Tiling (3×3)"
    
    public var id: String { rawValue }
}

public struct InitialLatticeVisualizerView: View {
    @ObservedObject var manager: StiffnessStudyManager
    
    @State private var viewMode: InitialMeshViewMode = .periodic3x3
    @State private var colorField: InitialMeshColorField = .uniform
    @State private var colorPalette: StiffnessColorPalette = .coolwarm
    
    @State private var showNodes: Bool = true
    @State private var showNodeLabels: Bool = false
    @State private var showWireframe: Bool = true
    @State private var showElementLabels: Bool = false
    @State private var highlightCrossing: Bool = true
    @State private var showDomainBorders: Bool = true
    @State private var nodeScale: Double = 1.0
    
    // Pan & Zoom
    @State private var zoomScale: Double = 1.0
    @State private var panOffset: CGSize = .zero
    @State private var dragStartOffset: CGSize = .zero
    
    // Hovered element
    @State private var hoveredElementId: Int? = nil
    @State private var hoveredElementCoords: [[Double]]? = nil
    @State private var hoveredElementEnergy: Double? = nil
    @State private var hoveredElementVM: Double? = nil
    @State private var hoveredElementShear: Double? = nil
    @State private var hoveredElementPress: Double? = nil
    @State private var hoveredElementNodes: [Int]? = nil
    
    public init(manager: StiffnessStudyManager) {
        self.manager = manager
    }
    
    private var mesh: InitialMeshModel? {
        manager.initialMeshData
    }
    
    public var body: some View {
        VStack(spacing: 0) {
            // TOP TOOLBAR
            topToolbar
            
            Divider()
            
            // MAIN CANVAS AREA
            ZStack(alignment: .bottomLeading) {
                if let mesh = mesh, !mesh.primary_nodes.isEmpty {
                    GeometryReader { geo in
                        canvasView(mesh: mesh, size: geo.size)
                            .gesture(
                                DragGesture()
                                    .onChanged { val in
                                        panOffset = CGSize(
                                            width: dragStartOffset.width + val.translation.width,
                                            height: dragStartOffset.height + val.translation.height
                                        )
                                    }
                                    .onEnded { _ in
                                        dragStartOffset = panOffset
                                    }
                            )
                    }
                    
                    // BOTTOM OVERLAYS: HUD INFO & COLORBAR LEGEND
                    HStack(alignment: .bottom, spacing: 12) {
                        // Left: Hover Inspector HUD
                        if let elemId = hoveredElementId {
                            hoverInspectorCard(elemId: elemId)
                                .transition(.opacity)
                        }
                        
                        Spacer()
                        
                        // Right: Colorbar Legend (if field is not uniform)
                        if colorField != .uniform {
                            colorbarLegend(mesh: mesh)
                                .padding(10)
                        }
                    }
                    .padding(14)
                } else {
                    emptyPlaceholder
                }
            }
            .frame(maxWidth: .infinity, maxHeight: .infinity)
            .background(Color(NSColor.underPageBackgroundColor))
        }
        .onAppear {
            if manager.initialMeshData == nil {
                manager.generateInitialMeshPreview()
            }
        }
    }
    
    // MARK: - Top Toolbar
    private var topToolbar: some View {
        ScrollView(.horizontal, showsIndicators: false) {
            HStack(spacing: 12) {
                // View Mode Picker (1x1 vs 3x3)
                Picker("", selection: $viewMode) {
                    ForEach(InitialMeshViewMode.allCases) { mode in
                        Text(mode.rawValue).tag(mode)
                    }
                }
                .pickerStyle(.segmented)
                .controlSize(.small)
                .fixedSize()
                
                Divider().frame(height: 18)
                
                // Color Field Picker
                HStack(spacing: 4) {
                    Text("Color:")
                        .font(.caption2.bold())
                        .foregroundColor(.secondary)
                    Picker("", selection: $colorField) {
                        ForEach(InitialMeshColorField.allCases) { field in
                            Text(field.rawValue).tag(field)
                        }
                    }
                    .pickerStyle(.menu)
                    .controlSize(.small)
                    .frame(width: 170)
                }
                
                // Palette Picker (when colored)
                if colorField != .uniform {
                    Picker("", selection: $colorPalette) {
                        ForEach(StiffnessColorPalette.allCases) { pal in
                            Text(pal.rawValue).tag(pal)
                        }
                    }
                    .pickerStyle(.menu)
                    .controlSize(.small)
                    .frame(width: 110)
                }
                
                Divider().frame(height: 18)
                
                // Display Toggles Menu
                Menu {
                    Toggle("Show Nodes", isOn: $showNodes)
                    Toggle("Node ID Numbers", isOn: $showNodeLabels)
                    Toggle("Element Wireframe", isOn: $showWireframe)
                    Toggle("Element ID Numbers", isOn: $showElementLabels)
                    Toggle("Highlight Boundary Crossing", isOn: $highlightCrossing)
                    Toggle("Domain Box Borders", isOn: $showDomainBorders)
                } label: {
                    Label("Layers", systemImage: "square.3.layers.3d")
                        .font(.caption)
                }
                .menuStyle(.borderlessButton)
                .fixedSize()
                
                Divider().frame(height: 18)
                
                // Zoom & Pan Actions
                HStack(spacing: 4) {
                    Button(action: { zoomScale = max(0.2, zoomScale - 0.25) }) {
                        Image(systemName: "minus.magnifyingglass")
                    }
                    .buttonStyle(.bordered)
                    .controlSize(.small)
                    
                    Text("\(Int(round(zoomScale * 100)))%")
                        .font(.system(size: 11, design: .monospaced))
                        .frame(width: 42)
                    
                    Button(action: { zoomScale = min(15.0, zoomScale + 0.25) }) {
                        Image(systemName: "plus.magnifyingglass")
                    }
                    .buttonStyle(.bordered)
                    .controlSize(.small)
                    
                    Button("Reset") {
                        withAnimation(.easeOut(duration: 0.2)) {
                            zoomScale = 1.0
                            panOffset = .zero
                            dragStartOffset = .zero
                        }
                    }
                    .buttonStyle(.bordered)
                    .controlSize(.small)
                }
                
                Divider().frame(height: 18)
                
                // Refresh / Preview Button
                Button(action: { manager.generateInitialMeshPreview() }) {
                    if manager.isGeneratingPreview {
                        ProgressView().scaleEffect(0.6).frame(width: 16, height: 16)
                    } else {
                        Label("Preview Mesh", systemImage: "arrow.triangle.2.circlepath")
                            .font(.caption)
                    }
                }
                .buttonStyle(.borderedProminent)
                .controlSize(.small)
                .disabled(manager.isGeneratingPreview)
                
                // Lattice Info Badge
                if let mesh = mesh {
                    HStack(spacing: 6) {
                        Text("\(mesh.lattice_type.capitalized) | θ=\(String(format: "%.2f°", mesh.rotation_angle))")
                            .font(.system(size: 10, weight: .bold, design: .monospaced))
                        if mesh.is_csl {
                            Text("Σ=\(mesh.sigma) (m=\(mesh.m), n=\(mesh.n))")
                                .font(.system(size: 10, weight: .medium, design: .monospaced))
                        }
                        Text("Box: \(String(format: "%.2f×%.2f", mesh.Lx, mesh.Ly)) | \(mesh.num_primary_nodes) atoms | \(mesh.num_elements) elements")
                            .font(.system(size: 10, design: .monospaced))
                            .foregroundColor(.secondary)
                    }
                    .padding(.horizontal, 8)
                    .padding(.vertical, 4)
                    .background(Color.accentColor.opacity(0.12), in: RoundedRectangle(cornerRadius: 6))
                    .fixedSize()
                }
            }
            .padding(.horizontal, 12)
            .padding(.vertical, 6)
        }
        .background(Color(NSColor.controlBackgroundColor))
    }
    
    // MARK: - Canvas Rendering View
    private func canvasView(mesh: InitialMeshModel, size: CGSize) -> some View {
        Canvas { ctx, canvasSize in
            let center = CGPoint(x: canvasSize.width / 2.0, y: canvasSize.height / 2.0)
            
            // Domain box size
            let boxLx = mesh.Lx
            let boxLy = mesh.Ly
            let extent = max(boxLx, boxLy) * (viewMode == .periodic3x3 ? 3.4 : 1.4)
            let baseScale = min(canvasSize.width, canvasSize.height) / max(extent, 1e-4)
            let scale = baseScale * zoomScale
            
            // Map global (x, y) to canvas (X, Y)
            let mapPoint: (Double, Double) -> CGPoint = { gx, gy in
                let normX = viewMode == .periodic3x3 ? (gx - boxLx / 2.0) : (gx - boxLx / 2.0)
                let normY = viewMode == .periodic3x3 ? (gy - boxLy / 2.0) : (gy - boxLy / 2.0)
                let px = center.x + panOffset.width + CGFloat(normX) * CGFloat(scale)
                let py = center.y + panOffset.height - CGFloat(normY) * CGFloat(scale)
                return CGPoint(x: px, y: py)
            }
            
            // 1. Draw Domain Bounding Boxes
            if showDomainBorders {
                let translations: [(Int, Int)] = (viewMode == .periodic3x3)
                    ? [(-1, 1), (0, 1), (1, 1), (-1, 0), (0, 0), (1, 0), (-1, -1), (0, -1), (1, -1)]
                    : [(0, 0)]
                
                for (tx, ty) in translations {
                    let minP = mapPoint(Double(tx) * boxLx, Double(ty) * boxLy)
                    let maxP = mapPoint(Double(tx + 1) * boxLx, Double(ty + 1) * boxLy)
                    let rect = CGRect(
                        x: min(minP.x, maxP.x),
                        y: min(minP.y, maxP.y),
                        width: abs(maxP.x - minP.x),
                        height: abs(maxP.y - minP.y)
                    )
                    
                    let isPrimary = (tx == 0 && ty == 0)
                    var path = Path()
                    path.addRect(rect)
                    
                    if isPrimary {
                        ctx.stroke(path, with: .color(Color.accentColor.opacity(0.9)), style: StrokeStyle(lineWidth: 2.0))
                        ctx.fill(path, with: .color(Color.accentColor.opacity(0.04)))
                    } else {
                        ctx.stroke(path, with: .color(Color.secondary.opacity(0.45)), style: StrokeStyle(lineWidth: 1.0, dash: [4, 4]))
                    }
                    
                    // Domain tag badge (e.g. [0, 0], [+1, 0])
                    if viewMode == .periodic3x3 {
                        let tagStr = isPrimary ? "Primary [0, 0]" : "[\(tx > 0 ? "+\(tx)" : "\(tx)"), \(ty > 0 ? "+\(ty)" : "\(ty)")]"
                        let text = Text(tagStr)
                            .font(.system(size: 9, weight: isPrimary ? .bold : .regular, design: .monospaced))
                            .foregroundColor(isPrimary ? Color.accentColor : Color.secondary)
                        ctx.draw(text, at: CGPoint(x: rect.minX + 6, y: rect.minY + 8), anchor: .topLeading)
                    }
                }
            }
            
            // 2. Draw Elements (Triangles)
            let fieldMin = colorFieldMin(mesh: mesh)
            let fieldMax = colorFieldMax(mesh: mesh)
            let fieldRange = max(1e-12, fieldMax - fieldMin)
            
            if viewMode == .periodic3x3 && !mesh.elements_3x3.isEmpty {
                for tri3 in mesh.elements_3x3 {
                    guard tri3.node_coords.count >= 3 else { continue }
                    let p0 = mapPoint(tri3.node_coords[0][0], tri3.node_coords[0][1])
                    let p1 = mapPoint(tri3.node_coords[1][0], tri3.node_coords[1][1])
                    let p2 = mapPoint(tri3.node_coords[2][0], tri3.node_coords[2][1])
                    
                    var path = Path()
                    path.move(to: p0)
                    path.addLine(to: p1)
                    path.addLine(to: p2)
                    path.closeSubpath()
                    
                    // Fill Color
                    if colorField == .uniform {
                        let isCentral = (tri3.domain_tx == 0 && tri3.domain_ty == 0)
                        ctx.fill(path, with: .color(isCentral ? Color.accentColor.opacity(0.18) : Color.gray.opacity(0.08)))
                    } else {
                        let val = valueForField(tri3: tri3)
                        let t = max(0.0, min(1.0, (val - fieldMin) / fieldRange))
                        let color = sampleColormap(t: t, palette: colorPalette)
                        ctx.fill(path, with: .color(color.opacity(0.85)))
                    }
                    
                    // Wireframe Edges
                    if showWireframe {
                        let isCentral = (tri3.domain_tx == 0 && tri3.domain_ty == 0)
                        let edgeColor = isCentral ? Color.black.opacity(0.4) : Color.gray.opacity(0.2)
                        ctx.stroke(path, with: .color(edgeColor), style: StrokeStyle(lineWidth: isCentral ? 0.8 : 0.4))
                    }
                }
            } else {
                for elem in mesh.elements {
                    guard elem.node_coords.count >= 3 else { continue }
                    let p0 = mapPoint(elem.node_coords[0][0], elem.node_coords[0][1])
                    let p1 = mapPoint(elem.node_coords[1][0], elem.node_coords[1][1])
                    let p2 = mapPoint(elem.node_coords[2][0], elem.node_coords[2][1])
                    
                    var path = Path()
                    path.move(to: p0)
                    path.addLine(to: p1)
                    path.addLine(to: p2)
                    path.closeSubpath()
                    
                    // Fill Color
                    if colorField == .uniform {
                        let crosses = (elem.crosses_boundary == true)
                        if highlightCrossing && crosses {
                            ctx.fill(path, with: .color(Color.orange.opacity(0.25)))
                        } else {
                            ctx.fill(path, with: .color(Color.accentColor.opacity(0.18)))
                        }
                    } else {
                        let val = valueForField(elem: elem)
                        let t = max(0.0, min(1.0, (val - fieldMin) / fieldRange))
                        let color = sampleColormap(t: t, palette: colorPalette)
                        ctx.fill(path, with: .color(color.opacity(0.85)))
                    }
                    
                    // Wireframe Edges
                    if showWireframe {
                        let crosses = (elem.crosses_boundary == true)
                        let edgeColor = (highlightCrossing && crosses) ? Color.orange.opacity(0.8) : Color.black.opacity(0.5)
                        let width: CGFloat = (highlightCrossing && crosses) ? 1.2 : 0.8
                        ctx.stroke(path, with: .color(edgeColor), style: StrokeStyle(lineWidth: width))
                    }
                    
                    // Element Centroid Label
                    if showElementLabels {
                        let cx = (p0.x + p1.x + p2.x) / 3.0
                        let cy = (p0.y + p1.y + p2.y) / 3.0
                        let label = Text("\(elem.id)")
                            .font(.system(size: 8, weight: .bold, design: .monospaced))
                            .foregroundColor(.primary)
                        ctx.draw(label, at: CGPoint(x: cx, y: cy), anchor: .center)
                    }
                }
            }
            
            // 3. Draw Nodes (Dots & Labels)
            if showNodes {
                let r: CGFloat = CGFloat(2.5 * nodeScale)
                
                if viewMode == .periodic3x3 {
                    // Draw translated nodes
                    for dom in mesh.translated_domains {
                        let isPrimary = (dom.tx == 0 && dom.ty == 0)
                        for node in dom.nodes {
                            let pt = mapPoint(node.x, node.y)
                            let circleRect = CGRect(x: pt.x - r, y: pt.y - r, width: 2 * r, height: 2 * r)
                            
                            if isPrimary {
                                ctx.fill(Path(ellipseIn: circleRect), with: .color(Color.blue))
                                ctx.stroke(Path(ellipseIn: circleRect), with: .color(Color.white), lineWidth: 0.8)
                            } else {
                                ctx.fill(Path(ellipseIn: circleRect), with: .color(Color.gray.opacity(0.65)))
                            }
                            
                            if showNodeLabels && isPrimary && zoomScale >= 1.5 {
                                let label = Text("\(node.id)")
                                    .font(.system(size: 8, weight: .semibold, design: .monospaced))
                                    .foregroundColor(.primary)
                                ctx.draw(label, at: CGPoint(x: pt.x, y: pt.y - r - 4), anchor: .bottom)
                            }
                        }
                    }
                } else {
                    for node in mesh.primary_nodes {
                        let pt = mapPoint(node.x, node.y)
                        let circleRect = CGRect(x: pt.x - r, y: pt.y - r, width: 2 * r, height: 2 * r)
                        let isBnd = (node.is_boundary == true)
                        
                        let fillCol = isBnd ? Color.red : Color.blue
                        ctx.fill(Path(ellipseIn: circleRect), with: .color(fillCol))
                        ctx.stroke(Path(ellipseIn: circleRect), with: .color(Color.white), lineWidth: 0.8)
                        
                        if showNodeLabels && zoomScale >= 1.5 {
                            let label = Text("\(node.id)")
                                .font(.system(size: 8, weight: .semibold, design: .monospaced))
                                .foregroundColor(.primary)
                            ctx.draw(label, at: CGPoint(x: pt.x, y: pt.y - r - 4), anchor: .bottom)
                        }
                    }
                }
            }
        }
    }
    
    // MARK: - Field Value Helpers
    private func valueForField(elem: InitialMeshElement) -> Double {
        switch colorField {
        case .uniform: return 0.0
        case .energyDensity: return elem.strain_energy_density
        case .vonMises: return elem.von_mises
        case .shearStress: return elem.sigma_xy
        case .pressure: return elem.pressure
        case .sigmaXX: return elem.sigma_xx
        case .sigmaYY: return elem.sigma_yy
        }
    }
    
    private func valueForField(tri3: InitialMesh3x3Element) -> Double {
        switch colorField {
        case .uniform: return 0.0
        case .energyDensity: return tri3.strain_energy_density
        case .vonMises: return tri3.von_mises
        case .shearStress: return tri3.sigma_xy
        case .pressure: return tri3.pressure
        case .sigmaXX: return tri3.sigma_xx
        case .sigmaYY: return tri3.sigma_yy
        }
    }
    
    private func colorFieldMin(mesh: InitialMeshModel) -> Double {
        switch colorField {
        case .uniform: return 0.0
        case .energyDensity: return mesh.statistics.min_energy
        case .vonMises: return mesh.statistics.min_von_mises
        case .shearStress: return mesh.statistics.min_shear
        case .pressure: return mesh.statistics.min_pressure
        case .sigmaXX, .sigmaYY: return min(mesh.statistics.min_shear, 0.0)
        }
    }
    
    private func colorFieldMax(mesh: InitialMeshModel) -> Double {
        switch colorField {
        case .uniform: return 1.0
        case .energyDensity: return max(mesh.statistics.max_energy, 1e-6)
        case .vonMises: return max(mesh.statistics.max_von_mises, 1e-6)
        case .shearStress: return max(mesh.statistics.max_shear, 1e-6)
        case .pressure: return max(mesh.statistics.max_pressure, 1e-6)
        case .sigmaXX, .sigmaYY: return max(mesh.statistics.max_shear, 1e-6)
        }
    }
    
    // MARK: - Colormap Sampling
    private func sampleColormap(t: Double, palette: StiffnessColorPalette) -> Color {
        let u = max(0.0, min(1.0, t))
        switch palette {
        case .coolwarm:
            let r = min(1.0, max(0.0, 0.23 + 0.77 * u))
            let b = min(1.0, max(0.0, 0.95 - 0.75 * u))
            let g = min(1.0, max(0.0, 1.0 - 2.0 * abs(u - 0.5)))
            return Color(red: r, green: g * 0.85, blue: b)
        case .viridis:
            let r = 0.267 + 0.733 * u * u
            let g = 0.004 + 0.896 * u
            let b = 0.329 + 0.5 * sin(u * .pi)
            return Color(red: min(1.0, r), green: min(1.0, g), blue: min(1.0, b))
        case .plasma:
            let r = min(1.0, 0.05 + 0.95 * u)
            let g = min(1.0, max(0.0, sin(u * .pi) * 0.9))
            let b = min(1.0, 0.53 * (1.0 - u) + 0.1)
            return Color(red: r, green: g, blue: b)
        case .grayscale:
            return Color(white: u)
        }
    }
    
    // MARK: - Colorbar Legend
    private func colorbarLegend(mesh: InitialMeshModel) -> some View {
        let fMin = colorFieldMin(mesh: mesh)
        let fMax = colorFieldMax(mesh: mesh)
        
        return VStack(alignment: .leading, spacing: 4) {
            Text(colorField.rawValue)
                .font(.system(size: 10, weight: .bold))
                .foregroundColor(.primary)
            
            HStack(spacing: 6) {
                Text(String(format: "%.3e", fMin))
                    .font(.system(size: 9, design: .monospaced))
                    .foregroundColor(.secondary)
                
                LinearGradient(
                    gradient: Gradient(colors: (0...10).map { sampleColormap(t: Double($0) / 10.0, palette: colorPalette) }),
                    startPoint: .leading,
                    endPoint: .trailing
                )
                .frame(width: 130, height: 10)
                .cornerRadius(3)
                .overlay(RoundedRectangle(cornerRadius: 3).stroke(Color.primary.opacity(0.2), lineWidth: 0.5))
                
                Text(String(format: "%.3e", fMax))
                    .font(.system(size: 9, design: .monospaced))
                    .foregroundColor(.secondary)
            }
        }
        .padding(8)
        .background(.ultraThinMaterial, in: RoundedRectangle(cornerRadius: 8))
        .shadow(color: Color.black.opacity(0.12), radius: 3, x: 0, y: 1)
    }
    
    // MARK: - Hover Inspector HUD
    private func hoverInspectorCard(elemId: Int) -> some View {
        VStack(alignment: .leading, spacing: 4) {
            HStack {
                Text("Element #\(elemId)")
                    .font(.caption.bold())
                Spacer()
                if let nodes = hoveredElementNodes {
                    Text("Nodes: \(nodes.map(String.init).joined(separator: ", "))")
                        .font(.system(size: 9, design: .monospaced))
                        .foregroundColor(.secondary)
                }
            }
            
            Divider().padding(.vertical, 2)
            
            if let e = hoveredElementEnergy {
                HStack {
                    Text("Strain Energy Density (W_e):")
                        .font(.caption2)
                        .foregroundColor(.secondary)
                    Spacer()
                    Text(String(format: "%.5e", e))
                        .font(.system(size: 10, weight: .semibold, design: .monospaced))
                }
            }
            
            if let vm = hoveredElementVM {
                HStack {
                    Text("von Mises Stress (σ_vM):")
                        .font(.caption2)
                        .foregroundColor(.secondary)
                    Spacer()
                    Text(String(format: "%.5e", vm))
                        .font(.system(size: 10, weight: .semibold, design: .monospaced))
                }
            }
            
            if let s = hoveredElementShear {
                HStack {
                    Text("Cauchy Shear Stress (σ_xy):")
                        .font(.caption2)
                        .foregroundColor(.secondary)
                    Spacer()
                    Text(String(format: "%.5e", s))
                        .font(.system(size: 10, weight: .semibold, design: .monospaced))
                }
            }
            
            if let p = hoveredElementPress {
                HStack {
                    Text("Hydrostatic Pressure (p):")
                        .font(.caption2)
                        .foregroundColor(.secondary)
                    Spacer()
                    Text(String(format: "%.5e", p))
                        .font(.system(size: 10, weight: .semibold, design: .monospaced))
                }
            }
        }
        .padding(10)
        .frame(width: 280)
        .background(.ultraThinMaterial, in: RoundedRectangle(cornerRadius: 8))
        .shadow(color: Color.black.opacity(0.15), radius: 4, x: 0, y: 2)
    }
    
    // MARK: - Empty Placeholder
    private var emptyPlaceholder: some View {
        VStack(spacing: 12) {
            Image(systemName: "square.grid.3x3.fill")
                .font(.system(size: 40))
                .foregroundColor(.secondary)
            Text("No Initial Mesh Data Available")
                .font(.headline)
            Text("Click 'Preview Mesh' above to generate and inspect the initial crystal triangulation.")
                .font(.caption)
                .foregroundColor(.secondary)
                .multilineTextAlignment(.center)
                .frame(maxWidth: 320)
            
            Button("Generate Preview") {
                manager.generateInitialMeshPreview()
            }
            .buttonStyle(.borderedProminent)
            .controlSize(.regular)
        }
    }
}
