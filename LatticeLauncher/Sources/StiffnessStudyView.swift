import SwiftUI
import AppKit

public struct StiffnessStudyView: View {
    @StateObject private var manager = StiffnessStudyManager()
    @State private var selectedDisplayTab: Int = 0 // 0: Chart, 1: Table, 2: Terminal
    @State private var hoveredPoint: StiffnessPointModel? = nil
    @State private var showHigherModes: Bool = false
    
    public init() {}
    
    public var body: some View {
        HSplitView {
            // LEFT SIDEBAR: Controls & Parameter Setup
            VStack(spacing: 0) {
                ScrollView {
                    VStack(alignment: .leading, spacing: 14) {
                        // Title Header
                        HStack(spacing: 10) {
                            Image(systemName: "point.3.filled.connected.trianglepath.dotted")
                                .font(.title)
                                .foregroundColor(.accentColor)
                            VStack(alignment: .leading, spacing: 2) {
                                Text("Stiffness Spectrum K(α)")
                                    .font(.headline)
                                Text("Smallest non-rigid eigenvalues of FEM stiffness matrix")
                                    .font(.caption2)
                                    .foregroundColor(.secondary)
                            }
                        }
                        .padding(.bottom, 2)
                        
                        // 1. Grid Size (Node Numbers)
                        GroupBox(label: Label("Lattice Grid & Symmetry", systemImage: "grid")) {
                            VStack(alignment: .leading, spacing: 10) {
                                Picker("Symmetry", selection: $manager.latticeSymmetry) {
                                    Text("Square").tag(LatticeSymmetry.square)
                                    Text("Triangular").tag(LatticeSymmetry.triangular)
                                }
                                .pickerStyle(.segmented)
                                
                                HStack(spacing: 16) {
                                    HStack {
                                        Text("Nx:")
                                            .font(.caption.bold())
                                            .foregroundColor(.secondary)
                                        IntField(title: "Nx", value: $manager.nx, width: 50)
                                        Stepper("", value: $manager.nx, in: 4...200, step: 2)
                                            .labelsHidden()
                                    }
                                    
                                    HStack {
                                        Text("Ny:")
                                            .font(.caption.bold())
                                            .foregroundColor(.secondary)
                                        IntField(title: "Ny", value: $manager.ny, width: 50)
                                        Stepper("", value: $manager.ny, in: 4...200, step: 2)
                                            .labelsHidden()
                                    }
                                }
                                
                                Text("Total Nodes: \(manager.nx * manager.ny)  |  Est. DOFs: \(2 * manager.nx * manager.ny)")
                                    .font(.caption2)
                                    .foregroundColor(.secondary)
                            }
                            .padding(.top, 4)
                        }
                        
                        // 2. Boundary Conditions
                        GroupBox(label: Label("Boundary Conditions", systemImage: "rectangle.dashed")) {
                            VStack(alignment: .leading, spacing: 8) {
                                Picker("Boundary Condition", selection: $manager.boundaryCondition) {
                                    ForEach(BoundaryConditionMode.allCases) { bc in
                                        Label(bc.displayName, systemImage: bc.icon).tag(bc)
                                    }
                                }
                                .pickerStyle(.segmented)
                                
                                if manager.boundaryCondition == .periodic {
                                    Text("Periodic cell across boundaries. Rigid translations (u_x, u_y) are automatically deflated.")
                                        .font(.caption2)
                                        .foregroundColor(.secondary)
                                } else {
                                    Text("Perimeter nodes are held fixed at sheared positions. Internal DOFs are unconstrained.")
                                        .font(.caption2)
                                        .foregroundColor(.secondary)
                                }
                            }
                            .padding(.top, 4)
                        }
                        
                        // 3. Constitutive Potential
                        GroupBox(label: Label("Underlying Energy Potential", systemImage: "atom")) {
                            VStack(alignment: .leading, spacing: 8) {
                                Picker("Potential Mode", selection: $manager.potentialMode) {
                                    ForEach(AcousticPotentialMode.allCases) { mode in
                                        Text(mode.displayName).tag(mode)
                                    }
                                }
                                .labelsHidden()
                                .pickerStyle(.radioGroup)
                                
                                if manager.potentialMode == .zanzotto {
                                    Divider().padding(.vertical, 2)
                                    VStack(alignment: .leading, spacing: 6) {
                                        HStack {
                                            Text("Cubic Well Parameter β:")
                                                .font(.caption2.bold())
                                            Spacer()
                                            DoubleField(title: "-0.25", value: $manager.beta, width: 70, minValue: -0.25, maxValue: 4.0)
                                        }
                                        Slider(value: $manager.beta, in: -0.25...4.0, step: 0.01)
                                        
                                        // Quick Presets
                                        ScrollView(.horizontal, showsIndicators: false) {
                                            HStack(spacing: 4) {
                                                Text("Presets:")
                                                    .font(.caption2)
                                                    .foregroundColor(.secondary)
                                                ForEach([-0.25, -0.1, 0.0, 0.5, 1.0, 2.0, 4.0], id: \.self) { val in
                                                    Button(String(format: "%.2f", val)) {
                                                        manager.beta = val
                                                    }
                                                    .buttonStyle(.bordered)
                                                    .controlSize(.mini)
                                                }
                                            }
                                        }
                                        
                                        Divider().padding(.vertical, 2)
                                        
                                        HStack {
                                            Text("Bulk Modulus K:")
                                                .font(.caption2.bold())
                                            Spacer()
                                            DoubleField(title: "5.0", value: $manager.bulkK, width: 70, minValue: 0.01)
                                        }
                                        Slider(value: $manager.bulkK, in: 1.0...20.0, step: 0.5)
                                    }
                                }
                            }
                            .padding(.top, 4)
                        }
                        
                        // 4. Simple Shear Path Sweep
                        GroupBox(label: Label("Simple Shear Path (α)", systemImage: "arrow.left.and.right")) {
                            VStack(alignment: .leading, spacing: 8) {
                                HStack(spacing: 12) {
                                    VStack(alignment: .leading, spacing: 2) {
                                        Text("α Min:")
                                            .font(.caption2.bold())
                                            .foregroundColor(.secondary)
                                        DoubleField(title: "-0.3", value: $manager.alphaMin, width: 75)
                                    }
                                    VStack(alignment: .leading, spacing: 2) {
                                        Text("α Max:")
                                            .font(.caption2.bold())
                                            .foregroundColor(.secondary)
                                        DoubleField(title: "0.3", value: $manager.alphaMax, width: 75)
                                    }
                                    VStack(alignment: .leading, spacing: 2) {
                                        Text("Step Δα:")
                                            .font(.caption2.bold())
                                            .foregroundColor(.secondary)
                                        DoubleField(title: "0.005", value: $manager.dalpha, width: 75, minValue: 1e-5, maxValue: 1.0)
                                    }
                                }
                                
                                let nSteps: Int = {
                                    guard manager.dalpha > 1e-9, !manager.dalpha.isNaN, !manager.dalpha.isInfinite else { return 0 }
                                    let span = manager.alphaMax - manager.alphaMin
                                    guard span >= 0, !span.isNaN, !span.isInfinite else { return 0 }
                                    let raw = span / manager.dalpha
                                    guard !raw.isNaN, !raw.isInfinite, raw >= 0, raw < 1_000_000 else { return 0 }
                                    return max(1, Int(round(raw))) + 1
                                }()
                                
                                if nSteps > 0 {
                                    Text("Total Shear Steps: \(nSteps)")
                                        .font(.caption2)
                                        .foregroundColor(.secondary)
                                } else {
                                    Text("Total Shear Steps: — (enter Δα > 0)")
                                        .font(.caption2)
                                        .foregroundColor(.orange)
                                }
                            }
                            .padding(.top, 4)
                        }
                        
                        // 5. Summary Metrics (if completed)
                        if !manager.points.isEmpty {
                            GroupBox(label: Label("Stability Summary", systemImage: "checkmark.seal.fill")) {
                                VStack(alignment: .leading, spacing: 6) {
                                    HStack {
                                        Text("Dislocation Threshold α_c:")
                                            .font(.caption2)
                                        Spacer()
                                        Text("[\(manager.summary.alpha_c_neg, specifier: "%.4f"), +\(manager.summary.alpha_c_pos, specifier: "%.4f")]")
                                            .font(.caption2.bold().monospaced())
                                            .foregroundColor(.red)
                                    }
                                    HStack {
                                        Text("Min Eigenvalue min λ_min:")
                                            .font(.caption2)
                                        Spacer()
                                        Text("\(manager.summary.min_lambda, specifier: "%.3e")")
                                            .font(.caption2.bold().monospaced())
                                            .foregroundColor(manager.summary.min_lambda < 0 ? .red : .green)
                                    }
                                    HStack {
                                        Text("Free DOFs / Elements:")
                                            .font(.caption2)
                                        Spacer()
                                        Text("\(manager.summary.total_dofs) DOFs / \(manager.summary.total_elements) elems")
                                            .font(.caption2.monospaced())
                                            .foregroundColor(.secondary)
                                    }
                                    HStack {
                                        Text("Compute Time:")
                                            .font(.caption2)
                                        Spacer()
                                        Text("\(manager.summary.execution_seconds, specifier: "%.3f") s")
                                            .font(.caption2.monospaced())
                                            .foregroundColor(.secondary)
                                    }
                                }
                                .padding(.top, 2)
                            }
                        }
                    }
                    .frame(maxWidth: .infinity, alignment: .leading)
                    .padding(14)
                }
                
                Divider()
                
                // Bottom Action Footer
                VStack(spacing: 8) {
                    if manager.isRunning {
                        VStack(spacing: 6) {
                            ProgressView(value: manager.progress)
                                .progressViewStyle(.linear)
                            HStack {
                                Text(manager.statusMessage)
                                    .font(.caption)
                                    .foregroundColor(.secondary)
                                    .lineLimit(1)
                                Spacer()
                                Button("Cancel", role: .destructive, action: { manager.stopStudy() })
                                    .buttonStyle(.borderedProminent)
                                    .tint(.red)
                                    .controlSize(.small)
                            }
                        }
                        .padding(8)
                        .background(Color.accentColor.opacity(0.1), in: RoundedRectangle(cornerRadius: 8))
                    } else {
                        Button(action: { manager.runStudy() }) {
                            Label("Run Stiffness Study", systemImage: "play.fill")
                                .frame(maxWidth: .infinity)
                                .padding(.vertical, 4)
                        }
                        .buttonStyle(.borderedProminent)
                        .keyboardShortcut(.defaultAction)
                        
                        if !manager.statusMessage.isEmpty && manager.statusMessage != "Ready" {
                            HStack(spacing: 4) {
                                Image(systemName: manager.statusMessage.contains("failed") || manager.statusMessage.contains("Error") ? "exclamationmark.triangle.fill" : "checkmark.circle.fill")
                                    .foregroundColor(manager.statusMessage.contains("failed") || manager.statusMessage.contains("Error") ? .orange : .green)
                                    .font(.caption2)
                                Text(manager.statusMessage)
                                    .font(.caption2)
                                    .foregroundColor(.secondary)
                                    .lineLimit(1)
                                Spacer()
                            }
                        }
                    }
                }
                .padding(12)
                .background(Color(NSColor.controlBackgroundColor))
            }
            .frame(minWidth: 350, idealWidth: 380, maxWidth: 430)
            
            // RIGHT MAIN AREA: Visualizations & Data
            VStack(spacing: 0) {
                // Header Bar: Pinned Tab Picker + Horizontally Scrollable Options & Badges
                HStack(spacing: 10) {
                    Picker("", selection: $selectedDisplayTab) {
                        Label("λ_min(α)", systemImage: "chart.xyaxis.line").tag(0)
                        Label("Energy E(α)", systemImage: "waveform.path.ecg").tag(1)
                        Label("Split", systemImage: "rectangle.split.2x1").tag(2)
                        Label("Table", systemImage: "tablecells").tag(3)
                        Label("Log", systemImage: "terminal").tag(4)
                    }
                    .labelsHidden()
                    .pickerStyle(.segmented)
                    .controlSize(.small)
                    .fixedSize()
                    
                    Divider()
                        .frame(height: 18)
                    
                    ScrollView(.horizontal, showsIndicators: false) {
                        HStack(spacing: 10) {
                            if (selectedDisplayTab == 0 || selectedDisplayTab == 2) && !manager.points.isEmpty {
                                Toggle(isOn: $showHigherModes) {
                                    Text("Show λ₂, λ₃")
                                        .font(.caption)
                                }
                                .toggleStyle(.checkbox)
                                .controlSize(.small)
                                .fixedSize()
                            }
                            
                            if !manager.points.isEmpty {
                                HStack(spacing: 5) {
                                    Text(String(format: "α_c: [%.4f, +%.4f]", manager.summary.alpha_c_neg, manager.summary.alpha_c_pos))
                                        .font(.system(size: 11, weight: .semibold, design: .monospaced))
                                        .foregroundColor(.red)
                                }
                                .padding(.horizontal, 8)
                                .padding(.vertical, 3)
                                .background(Color.red.opacity(0.10), in: Capsule())
                                .fixedSize()
                                
                                Button(action: exportCSV) {
                                    Label("Export CSV", systemImage: "square.and.arrow.up")
                                        .font(.caption)
                                }
                                .buttonStyle(.bordered)
                                .controlSize(.small)
                                .fixedSize()
                            }
                        }
                        .padding(.vertical, 1)
                    }
                    .frame(maxWidth: .infinity, alignment: .leading)
                }
                .frame(maxWidth: .infinity, alignment: .leading)
                .padding(.horizontal, 14)
                .padding(.vertical, 8)
                .background(Color(NSColor.controlBackgroundColor))
                
                Divider()
                
                // Content Views
                Group {
                    if selectedDisplayTab == 0 {
                        StiffnessChartView(
                            points: manager.points,
                            summary: manager.summary,
                            showHigherModes: showHigherModes,
                            hoveredPoint: $hoveredPoint
                        )
                    } else if selectedDisplayTab == 1 {
                        StiffnessEnergyChartView(
                            points: manager.points,
                            summary: manager.summary,
                            hoveredPoint: $hoveredPoint
                        )
                    } else if selectedDisplayTab == 2 {
                        HSplitView {
                            StiffnessChartView(
                                points: manager.points,
                                summary: manager.summary,
                                showHigherModes: showHigherModes,
                                hoveredPoint: $hoveredPoint
                            )
                            .frame(minWidth: 350)
                            
                            StiffnessEnergyChartView(
                                points: manager.points,
                                summary: manager.summary,
                                hoveredPoint: $hoveredPoint
                            )
                            .frame(minWidth: 350)
                        }
                    } else if selectedDisplayTab == 3 {
                        StiffnessTableView(points: manager.points)
                    } else {
                        StiffnessLogView(logLines: manager.logLines, statusMessage: manager.statusMessage)
                    }
                }
                .frame(maxWidth: .infinity, maxHeight: .infinity)
            }
        }
    }
    
    private func exportCSV() {
        let panel = NSSavePanel()
        panel.allowedContentTypes = [.commaSeparatedText]
        panel.nameFieldStringValue = "stiffness_shear_spectrum.csv"
        if panel.runModal() == .OK, let url = panel.url {
            var csv = "alpha,lambda_min,lambda_2,lambda_3,is_stable,energy,energy_density,continuum_energy\n"
            for pt in manager.points {
                csv += "\(pt.alpha),\(pt.lambda_min),\(pt.lambda_2),\(pt.lambda_3),\(pt.is_stable ? 1 : 0),\(pt.energy),\(pt.energy_density),\(pt.continuum_energy)\n"
            }
            try? csv.write(to: url, atomically: true, encoding: .utf8)
        }
    }
}

// MARK: - Interactive Canvas Chart for λ_min(α)
// MARK: - Interactive Canvas Chart for λ_min(α)
struct StiffnessChartView: View {
    let points: [StiffnessPointModel]
    let summary: StiffnessStudySummaryModel
    let showHigherModes: Bool
    @Binding var hoveredPoint: StiffnessPointModel?
    
    private let padLeft: CGFloat = 70
    private let padRight: CGFloat = 30
    private let padTop: CGFloat = 30
    private let padBottom: CGFloat = 50
    
    private var minAlpha: Double {
        points.map(\.alpha).min() ?? -0.3
    }
    private var maxAlpha: Double {
        points.map(\.alpha).max() ?? 0.3
    }
    private var rangeAlpha: Double {
        max(1e-6, maxAlpha - minAlpha)
    }
    
    private var minL: Double {
        let raw = points.map(\.lambda_min).min() ?? -1.0
        return min(raw * 1.1, -0.1)
    }
    private var maxL: Double {
        let raw = points.map(\.lambda_min).max() ?? 1.0
        return max(raw * 1.1, 0.1)
    }
    private var rangeL: Double {
        max(1e-6, maxL - minL)
    }
    
    private func plotWidth(for size: CGSize) -> CGFloat {
        max(50, size.width - padLeft - padRight)
    }
    
    private func plotHeight(for size: CGSize) -> CGFloat {
        max(50, size.height - padTop - padBottom)
    }
    
    private func toScreen(alpha: Double, lambda: Double, size: CGSize) -> CGPoint {
        let pw = plotWidth(for: size)
        let ph = plotHeight(for: size)
        let normX = (alpha - minAlpha) / rangeAlpha
        let normY = (lambda - minL) / rangeL
        let x = padLeft + CGFloat(normX) * pw
        let y = padTop + CGFloat(1.0 - normY) * ph
        return CGPoint(x: x, y: y)
    }
    
    var body: some View {
        GeometryReader { geo in
            if points.isEmpty {
                emptyPlaceholder
            } else {
                chartContainer(size: geo.size)
            }
        }
        .padding(8)
    }
    
    private var emptyPlaceholder: some View {
        VStack(spacing: 12) {
            Image(systemName: "point.3.filled.connected.trianglepath.dotted")
                .font(.system(size: 44))
                .foregroundColor(.secondary.opacity(0.5))
            Text("No Stiffness Spectrum Data")
                .font(.headline)
                .foregroundColor(.secondary)
            Text("Select your parameters on the left and click 'Run Stiffness Study'")
                .font(.caption)
                .foregroundColor(.secondary)
        }
        .frame(maxWidth: .infinity, maxHeight: .infinity)
    }
    
    @ViewBuilder
    private func chartContainer(size: CGSize) -> some View {
        let pw = plotWidth(for: size)
        
        ZStack(alignment: .topLeading) {
            Canvas { ctx, cSize in
                drawChart(ctx: &ctx, size: cSize)
            }
            
            // Top Legend Labels
            VStack {
                ScrollView(.horizontal, showsIndicators: false) {
                    HStack(spacing: 10) {
                        Text("λ_min > 0 : ELASTICALLY STABLE")
                            .font(.caption2.bold())
                            .foregroundColor(.green)
                            .padding(.horizontal, 6)
                            .padding(.vertical, 2)
                            .background(Color.green.opacity(0.12), in: RoundedRectangle(cornerRadius: 4))
                            .fixedSize()
                        
                        Text("λ_min < 0 : UNSTABLE / DISLOCATION NUCLEATION")
                            .font(.caption2.bold())
                            .foregroundColor(.red)
                            .padding(.horizontal, 6)
                            .padding(.vertical, 2)
                            .background(Color.red.opacity(0.12), in: RoundedRectangle(cornerRadius: 4))
                            .fixedSize()
                    }
                }
                .padding(.leading, padLeft + 8)
                .padding(.trailing, padRight + 8)
                .padding(.top, padTop + 6)
                Spacer()
            }
            
            // Critical Alpha Labels
            if summary.alpha_c_pos != 0.0 {
                let ptPos = toScreen(alpha: summary.alpha_c_pos, lambda: 0, size: size)
                Text(String(format: "+α_c = %.3f", summary.alpha_c_pos))
                    .font(.system(size: 10, weight: .bold, design: .monospaced))
                    .foregroundColor(.red)
                    .padding(.horizontal, 4)
                    .padding(.vertical, 2)
                    .background(Color(NSColor.windowBackgroundColor).opacity(0.85), in: RoundedRectangle(cornerRadius: 3))
                    .position(x: ptPos.x + 35, y: ptPos.y - 12)
            }
            if summary.alpha_c_neg != 0.0 {
                let ptNeg = toScreen(alpha: summary.alpha_c_neg, lambda: 0, size: size)
                Text(String(format: "-α_c = %.3f", summary.alpha_c_neg))
                    .font(.system(size: 10, weight: .bold, design: .monospaced))
                    .foregroundColor(.red)
                    .padding(.horizontal, 4)
                    .padding(.vertical, 2)
                    .background(Color(NSColor.windowBackgroundColor).opacity(0.85), in: RoundedRectangle(cornerRadius: 3))
                    .position(x: ptNeg.x - 35, y: ptNeg.y - 12)
            }
            
            // Hover Tooltip Overlay
            if let hovered = hoveredPoint {
                let scr = toScreen(alpha: hovered.alpha, lambda: hovered.lambda_min, size: size)
                tooltipView(for: hovered)
                    .position(
                        x: min(max(scr.x, padLeft + 60), padLeft + pw - 60),
                        y: max(scr.y - 50, padTop + 30)
                    )
            }
        }
        .contentShape(Rectangle())
        .gesture(
            DragGesture(minimumDistance: 0)
                .onChanged { val in
                    let locX = val.location.x
                    guard locX >= padLeft && locX <= padLeft + pw else {
                        hoveredPoint = nil
                        return
                    }
                    let normX = Double((locX - padLeft) / pw)
                    let targetAlpha = minAlpha + normX * rangeAlpha
                    hoveredPoint = points.min(by: { abs($0.alpha - targetAlpha) < abs($1.alpha - targetAlpha) })
                }
                .onEnded { _ in
                    hoveredPoint = nil
                }
        )
    }
    
    @ViewBuilder
    private func tooltipView(for hovered: StiffnessPointModel) -> some View {
        VStack(alignment: .leading, spacing: 3) {
            Text(String(format: "α = %.4f", hovered.alpha))
                .font(.caption.bold())
            Text(String(format: "λ_min = %.4e", hovered.lambda_min))
                .font(.caption2.monospaced())
                .foregroundColor(hovered.is_stable ? .green : .red)
            if showHigherModes {
                Text(String(format: "λ_2 = %.4e", hovered.lambda_2))
                    .font(.caption2.monospaced())
                    .foregroundColor(.orange)
                Text(String(format: "λ_3 = %.4e", hovered.lambda_3))
                    .font(.caption2.monospaced())
                    .foregroundColor(.purple)
            }
            Text(String(format: "Energy = %.4e", hovered.energy))
                .font(.caption2.monospaced())
                .foregroundColor(.secondary)
        }
        .padding(8)
        .background(Color(NSColor.controlBackgroundColor).opacity(0.95), in: RoundedRectangle(cornerRadius: 6))
        .shadow(radius: 3)
    }
    
    private func drawChart(ctx: inout GraphicsContext, size: CGSize) {
        let pw = plotWidth(for: size)
        let ph = plotHeight(for: size)
        let zeroY = toScreen(alpha: 0, lambda: 0, size: size).y
        
        // Unstable Region Fill (below zero line)
        if zeroY < padTop + ph {
            let rectUnstable = CGRect(x: padLeft, y: zeroY, width: pw, height: padTop + ph - zeroY)
            ctx.fill(Path(rectUnstable), with: .color(Color.red.opacity(0.08)))
        }
        // Stable Region Fill (above zero line)
        if zeroY > padTop {
            let rectStable = CGRect(x: padLeft, y: padTop, width: pw, height: zeroY - padTop)
            ctx.fill(Path(rectStable), with: .color(Color.green.opacity(0.06)))
        }
        
        // Zero Axis Y = 0
        var zeroPath = Path()
        zeroPath.move(to: CGPoint(x: padLeft, y: zeroY))
        zeroPath.addLine(to: CGPoint(x: padLeft + pw, y: zeroY))
        ctx.stroke(zeroPath, with: .color(Color.primary.opacity(0.6)), style: StrokeStyle(lineWidth: 1.5, dash: [5, 3]))
        
        // Vertical Axis at α = 0
        let zeroX = toScreen(alpha: 0, lambda: 0, size: size).x
        var zeroXAxis = Path()
        zeroXAxis.move(to: CGPoint(x: zeroX, y: padTop))
        zeroXAxis.addLine(to: CGPoint(x: zeroX, y: padTop + ph))
        ctx.stroke(zeroXAxis, with: .color(Color.primary.opacity(0.3)), style: StrokeStyle(lineWidth: 1.0, dash: [2, 2]))
        
        // Critical Alpha Marker Lines
        if summary.alpha_c_neg != 0.0 {
            let ptNeg = toScreen(alpha: summary.alpha_c_neg, lambda: 0, size: size)
            var pLine = Path()
            pLine.move(to: CGPoint(x: ptNeg.x, y: padTop))
            pLine.addLine(to: CGPoint(x: ptNeg.x, y: padTop + ph))
            ctx.stroke(pLine, with: .color(Color.red.opacity(0.7)), style: StrokeStyle(lineWidth: 1.0, dash: [4, 4]))
        }
        if summary.alpha_c_pos != 0.0 {
            let ptPos = toScreen(alpha: summary.alpha_c_pos, lambda: 0, size: size)
            var pLine = Path()
            pLine.move(to: CGPoint(x: ptPos.x, y: padTop))
            pLine.addLine(to: CGPoint(x: ptPos.x, y: padTop + ph))
            ctx.stroke(pLine, with: .color(Color.red.opacity(0.7)), style: StrokeStyle(lineWidth: 1.0, dash: [4, 4]))
        }
        
        // Higher modes λ_2, λ_3 (if enabled)
        if showHigherModes {
            var p2Path = Path()
            var p3Path = Path()
            for (idx, pt) in points.enumerated() {
                let scr2 = toScreen(alpha: pt.alpha, lambda: pt.lambda_2, size: size)
                let scr3 = toScreen(alpha: pt.alpha, lambda: pt.lambda_3, size: size)
                if idx == 0 {
                    p2Path.move(to: scr2)
                    p3Path.move(to: scr3)
                } else {
                    p2Path.addLine(to: scr2)
                    p3Path.addLine(to: scr3)
                }
            }
            ctx.stroke(p2Path, with: .color(Color.orange.opacity(0.6)), lineWidth: 1.2)
            ctx.stroke(p3Path, with: .color(Color.purple.opacity(0.6)), lineWidth: 1.2)
        }
        
        // Main Curve: λ_min(α)
        var linePath = Path()
        for (idx, pt) in points.enumerated() {
            let scr = toScreen(alpha: pt.alpha, lambda: pt.lambda_min, size: size)
            if idx == 0 {
                linePath.move(to: scr)
            } else {
                linePath.addLine(to: scr)
            }
        }
        ctx.stroke(linePath, with: .color(Color.blue), lineWidth: 2.5)
        
        // Plot Points
        for pt in points {
            let scr = toScreen(alpha: pt.alpha, lambda: pt.lambda_min, size: size)
            let ptColor = pt.is_stable ? Color.green : Color.red
            let circleRect = CGRect(x: scr.x - 3.5, y: scr.y - 3.5, width: 7, height: 7)
            ctx.fill(Path(ellipseIn: circleRect), with: .color(ptColor))
        }
    }
}

// MARK: - Interactive Canvas Chart for Strain Energy E(α) vs Continuum Theory
struct StiffnessEnergyChartView: View {
    let points: [StiffnessPointModel]
    let summary: StiffnessStudySummaryModel
    @Binding var hoveredPoint: StiffnessPointModel?
    @State private var energyMode: Int = 0 // 0: Density (per unit area), 1: Total Integrated
    
    private let padLeft: CGFloat = 70
    private let padRight: CGFloat = 30
    private let padTop: CGFloat = 35
    private let padBottom: CGFloat = 50
    
    private var minAlpha: Double {
        points.map(\.alpha).min() ?? -0.3
    }
    private var maxAlpha: Double {
        points.map(\.alpha).max() ?? 0.3
    }
    private var rangeAlpha: Double {
        max(1e-6, maxAlpha - minAlpha)
    }
    
    private func getFEMEnergy(_ pt: StiffnessPointModel) -> Double {
        energyMode == 0 ? pt.energy_density : pt.energy
    }
    
    private func getContEnergy(_ pt: StiffnessPointModel) -> Double {
        if energyMode == 0 {
            return pt.continuum_energy
        } else {
            let scale = (pt.energy_density > 1e-12) ? (pt.energy / pt.energy_density) : 1.0
            return pt.continuum_energy * scale
        }
    }
    
    private var minE: Double {
        let allValues = points.flatMap { [getFEMEnergy($0), getContEnergy($0)] }
        let raw = allValues.min() ?? 0.0
        return min(raw, 0.0)
    }
    
    private var maxE: Double {
        let allValues = points.flatMap { [getFEMEnergy($0), getContEnergy($0)] }
        let raw = allValues.max() ?? 1.0
        return max(raw * 1.1, 0.001)
    }
    
    private var rangeE: Double {
        max(1e-6, maxE - minE)
    }
    
    private func plotWidth(for size: CGSize) -> CGFloat {
        max(50, size.width - padLeft - padRight)
    }
    
    private func plotHeight(for size: CGSize) -> CGFloat {
        max(50, size.height - padTop - padBottom)
    }
    
    private func toScreen(alpha: Double, energy: Double, size: CGSize) -> CGPoint {
        let pw = plotWidth(for: size)
        let ph = plotHeight(for: size)
        let normX = (alpha - minAlpha) / rangeAlpha
        let normY = (energy - minE) / rangeE
        let x = padLeft + CGFloat(normX) * pw
        let y = padTop + CGFloat(1.0 - normY) * ph
        return CGPoint(x: x, y: y)
    }
    
    var body: some View {
        GeometryReader { geo in
            if points.isEmpty {
                emptyPlaceholder
            } else {
                chartContainer(size: geo.size)
            }
        }
        .padding(8)
    }
    
    private var emptyPlaceholder: some View {
        VStack(spacing: 12) {
            Image(systemName: "waveform.path.ecg")
                .font(.system(size: 44))
                .foregroundColor(.secondary.opacity(0.5))
            Text("No Energy Data Available")
                .font(.headline)
                .foregroundColor(.secondary)
            Text("Run a simulation on the left to compare FEM strain energy with continuum theory.")
                .font(.caption)
                .foregroundColor(.secondary)
        }
        .frame(maxWidth: .infinity, maxHeight: .infinity)
    }
    
    @ViewBuilder
    private func chartContainer(size: CGSize) -> some View {
        let pw = plotWidth(for: size)
        
        ZStack(alignment: .topLeading) {
            Canvas { ctx, cSize in
                drawChart(ctx: &ctx, size: cSize)
            }
            
            // Top Legend & Energy Mode Picker
            VStack {
                ScrollView(.horizontal, showsIndicators: false) {
                    HStack(spacing: 12) {
                        Picker("", selection: $energyMode) {
                            Text("Per Triangle").tag(0)
                            Text("Total Integrated").tag(1)
                        }
                        .labelsHidden()
                        .pickerStyle(.segmented)
                        .controlSize(.small)
                        .fixedSize()
                        
                        HStack(spacing: 5) {
                            Circle().fill(Color.blue).frame(width: 8, height: 8)
                            Text("Lattice FEM")
                                .font(.caption.bold())
                                .foregroundColor(.blue)
                        }
                        .fixedSize()
                        
                        HStack(spacing: 5) {
                            Rectangle().fill(Color.orange).frame(width: 14, height: 2)
                            Text("Continuum W(C)")
                                .font(.caption.bold())
                                .foregroundColor(.orange)
                        }
                        .fixedSize()
                    }
                }
                .padding(.leading, padLeft + 8)
                .padding(.trailing, padRight + 8)
                .padding(.top, padTop - 26)
                
                Spacer()
            }
            
            // Critical Alpha Labels
            if summary.alpha_c_pos != 0.0 {
                let ptPos = toScreen(alpha: summary.alpha_c_pos, energy: minE, size: size)
                Text(String(format: "+α_c = %.3f", summary.alpha_c_pos))
                    .font(.system(size: 10, weight: .bold, design: .monospaced))
                    .foregroundColor(.red)
                    .padding(.horizontal, 4)
                    .padding(.vertical, 2)
                    .background(Color(NSColor.windowBackgroundColor).opacity(0.85), in: RoundedRectangle(cornerRadius: 3))
                    .position(x: ptPos.x + 35, y: padTop + 20)
            }
            if summary.alpha_c_neg != 0.0 {
                let ptNeg = toScreen(alpha: summary.alpha_c_neg, energy: minE, size: size)
                Text(String(format: "-α_c = %.3f", summary.alpha_c_neg))
                    .font(.system(size: 10, weight: .bold, design: .monospaced))
                    .foregroundColor(.red)
                    .padding(.horizontal, 4)
                    .padding(.vertical, 2)
                    .background(Color(NSColor.windowBackgroundColor).opacity(0.85), in: RoundedRectangle(cornerRadius: 3))
                    .position(x: ptNeg.x - 35, y: padTop + 20)
            }
            
            // Hover Tooltip Overlay
            if let hovered = hoveredPoint {
                let scr = toScreen(alpha: hovered.alpha, energy: getFEMEnergy(hovered), size: size)
                tooltipView(for: hovered)
                    .position(
                        x: min(max(scr.x, padLeft + 90), padLeft + pw - 90),
                        y: max(scr.y - 65, padTop + 45)
                    )
            }
        }
        .contentShape(Rectangle())
        .gesture(
            DragGesture(minimumDistance: 0)
                .onChanged { val in
                    let locX = val.location.x
                    guard locX >= padLeft && locX <= padLeft + pw else {
                        hoveredPoint = nil
                        return
                    }
                    let normX = Double((locX - padLeft) / pw)
                    let targetAlpha = minAlpha + normX * rangeAlpha
                    hoveredPoint = points.min(by: { abs($0.alpha - targetAlpha) < abs($1.alpha - targetAlpha) })
                }
                .onEnded { _ in
                    hoveredPoint = nil
                }
        )
    }
    
    @ViewBuilder
    private func tooltipView(for hovered: StiffnessPointModel) -> some View {
        let eFem = getFEMEnergy(hovered)
        let eCont = getContEnergy(hovered)
        let diff = abs(eFem - eCont)
        
        VStack(alignment: .leading, spacing: 3) {
            Text(String(format: "α = %.4f", hovered.alpha))
                .font(.caption.bold())
            Text(String(format: "E_FEM  = %.5e", eFem))
                .font(.caption2.monospaced())
                .foregroundColor(.blue)
            Text(String(format: "W_cont = %.5e", eCont))
                .font(.caption2.monospaced())
                .foregroundColor(.orange)
            Text(String(format: "|ΔE|   = %.3e", diff))
                .font(.caption2.monospaced())
                .foregroundColor(.secondary)
            Text(String(format: "min λ(K) = %.4e", hovered.lambda_min))
                .font(.caption2.monospaced())
                .foregroundColor(hovered.is_stable ? .green : .red)
        }
        .padding(8)
        .background(Color(NSColor.controlBackgroundColor).opacity(0.95), in: RoundedRectangle(cornerRadius: 6))
        .shadow(radius: 3)
    }
    
    private func drawChart(ctx: inout GraphicsContext, size: CGSize) {
        let pw = plotWidth(for: size)
        let ph = plotHeight(for: size)
        let zeroY = toScreen(alpha: 0, energy: 0, size: size).y
        
        // Zero Axis Line
        if zeroY >= padTop && zeroY <= padTop + ph {
            var zeroPath = Path()
            zeroPath.move(to: CGPoint(x: padLeft, y: zeroY))
            zeroPath.addLine(to: CGPoint(x: padLeft + pw, y: zeroY))
            ctx.stroke(zeroPath, with: .color(Color.primary.opacity(0.4)), style: StrokeStyle(lineWidth: 1.0, dash: [4, 4]))
        }
        
        // Vertical Axis at α = 0
        let zeroX = toScreen(alpha: 0, energy: 0, size: size).x
        if zeroX >= padLeft && zeroX <= padLeft + pw {
            var zeroXAxis = Path()
            zeroXAxis.move(to: CGPoint(x: zeroX, y: padTop))
            zeroXAxis.addLine(to: CGPoint(x: zeroX, y: padTop + ph))
            ctx.stroke(zeroXAxis, with: .color(Color.primary.opacity(0.3)), style: StrokeStyle(lineWidth: 1.0, dash: [2, 2]))
        }
        
        // Critical Alpha Marker Lines
        if summary.alpha_c_neg != 0.0 {
            let ptNeg = toScreen(alpha: summary.alpha_c_neg, energy: 0, size: size)
            var pLine = Path()
            pLine.move(to: CGPoint(x: ptNeg.x, y: padTop))
            pLine.addLine(to: CGPoint(x: ptNeg.x, y: padTop + ph))
            ctx.stroke(pLine, with: .color(Color.red.opacity(0.7)), style: StrokeStyle(lineWidth: 1.0, dash: [4, 4]))
        }
        if summary.alpha_c_pos != 0.0 {
            let ptPos = toScreen(alpha: summary.alpha_c_pos, energy: 0, size: size)
            var pLine = Path()
            pLine.move(to: CGPoint(x: ptPos.x, y: padTop))
            pLine.addLine(to: CGPoint(x: ptPos.x, y: padTop + ph))
            ctx.stroke(pLine, with: .color(Color.red.opacity(0.7)), style: StrokeStyle(lineWidth: 1.0, dash: [4, 4]))
        }
        
        // 1. Analytical Continuum Energy Curve (dashed orange)
        var contPath = Path()
        for (idx, pt) in points.enumerated() {
            let scr = toScreen(alpha: pt.alpha, energy: getContEnergy(pt), size: size)
            if idx == 0 {
                contPath.move(to: scr)
            } else {
                contPath.addLine(to: scr)
            }
        }
        ctx.stroke(contPath, with: .color(Color.orange), style: StrokeStyle(lineWidth: 2.2, dash: [5, 3]))
        
        // 2. Simulated Lattice Energy Curve (solid blue)
        var simPath = Path()
        for (idx, pt) in points.enumerated() {
            let scr = toScreen(alpha: pt.alpha, energy: getFEMEnergy(pt), size: size)
            if idx == 0 {
                simPath.move(to: scr)
            } else {
                simPath.addLine(to: scr)
            }
        }
        ctx.stroke(simPath, with: .color(Color.blue), lineWidth: 2.5)
        
        // 3. Points on Simulated Curve
        for pt in points {
            let scr = toScreen(alpha: pt.alpha, energy: getFEMEnergy(pt), size: size)
            let circleRect = CGRect(x: scr.x - 3.0, y: scr.y - 3.0, width: 6, height: 6)
            ctx.fill(Path(ellipseIn: circleRect), with: .color(Color.blue))
        }
    }
}

// MARK: - Data Table View
struct StiffnessTableView: View {
    let points: [StiffnessPointModel]
    
    var body: some View {
        Table(points) {
            TableColumn("Shear α") { pt in
                Text(String(format: "%.4f", pt.alpha))
                    .monospaced()
            }
            .width(min: 75, ideal: 85)
            
            TableColumn("min λ(K)") { pt in
                Text(String(format: "%.4e", pt.lambda_min))
                    .monospaced()
                    .foregroundColor(pt.is_stable ? .primary : .red)
            }
            .width(min: 110, ideal: 120)
            
            TableColumn("λ₂") { pt in
                Text(String(format: "%.4e", pt.lambda_2))
                    .monospaced()
                    .foregroundColor(.secondary)
            }
            .width(min: 100, ideal: 110)
            
            TableColumn("Status") { pt in
                HStack(spacing: 4) {
                    Circle()
                        .fill(pt.is_stable ? Color.green : Color.red)
                        .frame(width: 8, height: 8)
                    Text(pt.is_stable ? "STABLE" : "UNSTABLE")
                        .font(.caption2.bold())
                        .foregroundColor(pt.is_stable ? .green : .red)
                }
            }
            .width(min: 90, ideal: 100)
            
            TableColumn("Energy / Triangle") { pt in
                Text(String(format: "%.4e", pt.energy_density))
                    .monospaced()
            }
            .width(min: 110, ideal: 120)
            
            TableColumn("Continuum W(C)") { pt in
                Text(String(format: "%.4e", pt.continuum_energy))
                    .monospaced()
                    .foregroundColor(.orange)
            }
            .width(min: 115, ideal: 125)
            
            TableColumn("Total Energy") { pt in
                Text(String(format: "%.4e", pt.energy))
                    .monospaced()
                    .foregroundColor(.secondary)
            }
            .width(min: 100, ideal: 110)
        }
    }
}

// MARK: - Console Log View
struct StiffnessLogView: View {
    let logLines: [String]
    let statusMessage: String
    
    var body: some View {
        VStack(alignment: .leading, spacing: 0) {
            HStack {
                Text("Status: \(statusMessage)")
                    .font(.caption.bold())
                Spacer()
                Text("\(logLines.count) log lines")
                    .font(.caption2)
                    .foregroundColor(.secondary)
            }
            .padding(.horizontal, 12)
            .padding(.vertical, 6)
            .background(Color(NSColor.controlBackgroundColor))
            
            Divider()
            
            ScrollView {
                LazyVStack(alignment: .leading, spacing: 2) {
                    ForEach(logLines.indices, id: \.self) { idx in
                        Text(logLines[idx])
                            .font(.system(size: 11, design: .monospaced))
                            .textSelection(.enabled)
                    }
                }
                .padding(10)
            }
            .background(Color(NSColor.textBackgroundColor))
        }
    }
}
