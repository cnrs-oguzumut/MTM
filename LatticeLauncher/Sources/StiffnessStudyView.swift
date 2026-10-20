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
                        
                        // 2b. Rotated Crystal Orientation
                        GroupBox(label: Label("Rotated Crystal Orientation", systemImage: "rotate.right")) {
                            VStack(alignment: .leading, spacing: 10) {
                                Toggle(isOn: $manager.useRotation) {
                                    Text("Enable Rotated Lattice")
                                        .font(.caption.bold())
                                }
                                .toggleStyle(.checkbox)
                                
                                if manager.useRotation {
                                    if manager.boundaryCondition == .fixed {
                                        VStack(alignment: .leading, spacing: 6) {
                                            Text("Arbitrary Rotation (Fixed BC):")
                                                .font(.caption2)
                                                .foregroundColor(.secondary)
                                            HStack {
                                                Text("Angle θ:")
                                                    .font(.caption2.bold())
                                                Spacer()
                                                DoubleField(title: "0.0", value: $manager.rotationAngle, width: 65, minValue: 0.0, maxValue: 360.0)
                                                Text("°").font(.caption2)
                                            }
                                            Slider(value: $manager.rotationAngle, in: 0.0...360.0, step: 0.5)
                                        }
                                    } else {
                                        // Periodic BC: CSL Discrete Angles
                                        VStack(alignment: .leading, spacing: 8) {
                                            Text("Coincident Site Lattice (CSL) for PBC:")
                                                .font(.caption2)
                                                .foregroundColor(.secondary)
                                            
                                            // CSL Presets Picker
                                            let allowedCSL = CSLAngleOption.allowedAngles(for: manager.latticeSymmetry)
                                            Picker("CSL Angle", selection: Binding(
                                                get: {
                                                    allowedCSL.first(where: { $0.m == manager.cslM && $0.n == manager.cslN }) ?? allowedCSL[0]
                                                },
                                                set: { opt in
                                                    manager.cslM = opt.m
                                                    manager.cslN = opt.n
                                                    manager.rotationAngle = opt.theta_deg
                                                }
                                            )) {
                                                ForEach(allowedCSL) { opt in
                                                    Text(opt.displayName).tag(opt)
                                                }
                                            }
                                            .pickerStyle(.menu)
                                            .controlSize(.small)
                                            
                                            // Coprime (m, n) manual overrides
                                            HStack(spacing: 12) {
                                                HStack {
                                                    Text("m:")
                                                        .font(.caption2.bold())
                                                    IntField(title: "m", value: $manager.cslM, width: 40)
                                                    Stepper("", value: $manager.cslM, in: -20...20)
                                                        .labelsHidden()
                                                }
                                                HStack {
                                                    Text("n:")
                                                        .font(.caption2.bold())
                                                    IntField(title: "n", value: $manager.cslN, width: 40)
                                                    Stepper("", value: $manager.cslN, in: 0...20)
                                                        .labelsHidden()
                                                }
                                            }
                                            
                                            // Supercell Multipliers (Kx, Ky)
                                            HStack(spacing: 12) {
                                                HStack {
                                                    Text("Kx:")
                                                        .font(.caption2.bold())
                                                    IntField(title: "Kx", value: $manager.cslKx, width: 40)
                                                    Stepper("", value: $manager.cslKx, in: 1...10)
                                                        .labelsHidden()
                                                }
                                                HStack {
                                                    Text("Ky:")
                                                        .font(.caption2.bold())
                                                    IntField(title: "Ky", value: $manager.cslKy, width: 40)
                                                    Stepper("", value: $manager.cslKy, in: 1...10)
                                                        .labelsHidden()
                                                }
                                            }
                                        }
                                    }
                                    
                                    Button(action: {
                                        manager.generateInitialMeshPreview()
                                        selectedDisplayTab = 6
                                    }) {
                                        HStack {
                                            if manager.isGeneratingPreview {
                                                ProgressView().scaleEffect(0.6).frame(width: 14, height: 14)
                                            } else {
                                                Image(systemName: "eye.fill")
                                            }
                                            Text("Preview Initial Mesh")
                                        }
                                        .frame(maxWidth: .infinity)
                                    }
                                    .buttonStyle(.borderedProminent)
                                    .controlSize(.small)
                                    .padding(.top, 2)
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
                        Label("Mode v(α)", systemImage: "scope").tag(5)
                        Label("Initial Mesh", systemImage: "square.grid.3x3.fill").tag(6)
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
                                    Text("Show Higher Modes (λ₂…λ₅)")
                                        .font(.caption)
                                }
                                .toggleStyle(.checkbox)
                                .controlSize(.small)
                                .fixedSize()
                                
                                Toggle(isOn: $manager.uniqueLevelsOnly) {
                                    Text("Unique Levels (no doublets)")
                                        .font(.caption)
                                }
                                .toggleStyle(.checkbox)
                                .controlSize(.small)
                                .fixedSize()
                                .help("Filter out sin/cos Bloch doublet partners under PBC to compute 5 distinct physical energy levels")
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
                    } else if selectedDisplayTab == 5 {
                        StiffnessModeInspectorView(manager: manager)
                    } else if selectedDisplayTab == 6 {
                        InitialLatticeVisualizerView(manager: manager)
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
            var csv = "alpha,lambda_min,lambda_2,lambda_3,lambda_4,lambda_5,is_stable,energy,energy_density,continuum_energy\n"
            for pt in manager.points {
                csv += "\(pt.alpha),\(pt.lambda_min),\(pt.lambda_2),\(pt.lambda_3),\(pt.lambda_4),\(pt.lambda_5),\(pt.is_stable ? 1 : 0),\(pt.energy),\(pt.energy_density),\(pt.continuum_energy)\n"
            }
            try? csv.write(to: url, atomically: true, encoding: .utf8)
        }
    }
}

// MARK: - Interactive Canvas Chart for λ_min(α)
struct StiffnessChartView: View {
    let points: [StiffnessPointModel]
    let summary: StiffnessStudySummaryModel
    let showHigherModes: Bool
    @Binding var hoveredPoint: StiffnessPointModel?
    
    // 0: Normalized (±1 split scale so positive & negative branches each use full half-height)
    // 1: Zoom Stable (uniform scale normalized by max positive λ, negative clamped near -1.2)
    // 2: Raw (unscaled eigenvalues)
    @State private var scaleMode: Int = 0
    
    private let padLeft: CGFloat = 78
    private let padRight: CGFloat = 30
    private let padTop: CGFloat = 40
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
    
    // Reference positive scale Λ_+ (largest positive eigenvalue among displayed modes)
    private var posScale: Double {
        var vals: [Double] = points.map(\.lambda_min).filter { $0 > 0 }
        if showHigherModes {
            vals.append(contentsOf: points.map(\.lambda_2).filter { $0 > 0 })
            vals.append(contentsOf: points.map(\.lambda_3).filter { $0 > 0 })
            vals.append(contentsOf: points.map(\.lambda_4).filter { $0 > 0 })
            vals.append(contentsOf: points.map(\.lambda_5).filter { $0 > 0 })
        }
        return max(1e-12, vals.max() ?? 1.0)
    }
    
    // Reference negative scale |Λ_-| (largest negative magnitude among displayed modes)
    private var negScale: Double {
        var vals: [Double] = points.map { -$0.lambda_min }.filter { $0 > 0 }
        if showHigherModes {
            vals.append(contentsOf: points.map { -$0.lambda_2 }.filter { $0 > 0 })
            vals.append(contentsOf: points.map { -$0.lambda_3 }.filter { $0 > 0 })
            vals.append(contentsOf: points.map { -$0.lambda_4 }.filter { $0 > 0 })
            vals.append(contentsOf: points.map { -$0.lambda_5 }.filter { $0 > 0 })
        }
        return max(posScale, vals.max() ?? posScale)
    }
    
    private func transformLambda(_ raw: Double) -> Double {
        switch scaleMode {
        case 0: // Normalized split scale: positive -> [0, 1], negative -> [-1, 0]
            if raw >= 0 {
                return raw / posScale
            } else {
                return raw / negScale
            }
        case 1: // Zoom Stable: uniform scale raw / posScale, clamped at -1.25 below zero
            return max(-1.25, raw / posScale)
        default: // Raw
            return raw
        }
    }
    
    private var minL: Double {
        switch scaleMode {
        case 0:
            let hasNeg = points.contains { $0.lambda_min < 0 }
            return hasNeg ? -1.12 : -0.15
        case 1:
            let rawMin = points.map { transformLambda($0.lambda_min) }.min() ?? -0.2
            return min(rawMin * 1.08, -0.15)
        default:
            var allVals = points.map(\.lambda_min)
            if showHigherModes {
                allVals.append(contentsOf: points.map(\.lambda_2))
                allVals.append(contentsOf: points.map(\.lambda_3))
                allVals.append(contentsOf: points.map(\.lambda_4))
                allVals.append(contentsOf: points.map(\.lambda_5))
            }
            let raw = allVals.min() ?? -1.0
            return min(raw * 1.1, -0.05)
        }
    }
    
    private var maxL: Double {
        switch scaleMode {
        case 0, 1:
            return 1.15
        default:
            var allVals = points.map(\.lambda_min)
            if showHigherModes {
                allVals.append(contentsOf: points.map(\.lambda_2))
                allVals.append(contentsOf: points.map(\.lambda_3))
                allVals.append(contentsOf: points.map(\.lambda_4))
                allVals.append(contentsOf: points.map(\.lambda_5))
            }
            let raw = allVals.max() ?? 1.0
            return max(raw * 1.15, 1e-4)
        }
    }
    
    private var rangeL: Double {
        max(1e-9, maxL - minL)
    }
    
    // Effective critical alpha (from summary or live points fallback)
    private var effectiveAlphaCNeg: Double {
        if summary.alpha_c_neg != 0.0 { return summary.alpha_c_neg }
        for i in 1..<points.count {
            let p0 = points[i - 1], p1 = points[i]
            if (p0.lambda_min <= 0 && p1.lambda_min > 0) || (p0.lambda_min > 0 && p1.lambda_min <= 0) {
                let frac = abs(p0.lambda_min) / max(1e-15, abs(p0.lambda_min) + abs(p1.lambda_min))
                let cross = p0.alpha + frac * (p1.alpha - p0.alpha)
                if cross < 0 { return cross }
            }
        }
        return 0.0
    }
    
    private var effectiveAlphaCPos: Double {
        if summary.alpha_c_pos != 0.0 { return summary.alpha_c_pos }
        for i in 1..<points.count {
            let p0 = points[i - 1], p1 = points[i]
            if (p0.lambda_min <= 0 && p1.lambda_min > 0) || (p0.lambda_min > 0 && p1.lambda_min <= 0) {
                let frac = abs(p0.lambda_min) / max(1e-15, abs(p0.lambda_min) + abs(p1.lambda_min))
                let cross = p0.alpha + frac * (p1.alpha - p0.alpha)
                if cross > 0 { return cross }
            }
        }
        return 0.0
    }
    
    // Check if λ_2 is degenerate with λ_min in the stable region (typical for PBC Bloch pairs)
    private var isLambda2DegenerateWithMin: Bool {
        let stablePts = points.filter(\.is_stable)
        guard !stablePts.isEmpty else { return false }
        let maxDiff = stablePts.map { abs($0.lambda_2 - $0.lambda_min) }.max() ?? 0.0
        return maxDiff < 1e-5 * posScale
    }
    
    // Check if λ_4 is degenerate with λ_3 in the stable region (second PBC Bloch doublet)
    private var isLambda4DegenerateWith3: Bool {
        let stablePts = points.filter(\.is_stable)
        guard !stablePts.isEmpty else { return false }
        let maxDiff = stablePts.map { abs($0.lambda_4 - $0.lambda_3) }.max() ?? 0.0
        return maxDiff < 1e-5 * posScale
    }
    
    private func plotWidth(for size: CGSize) -> CGFloat {
        max(50, size.width - padLeft - padRight)
    }
    
    private func plotHeight(for size: CGSize) -> CGFloat {
        max(50, size.height - padTop - padBottom)
    }
    
    private func toScreen(alpha: Double, rawLambda: Double, size: CGSize) -> CGPoint {
        let pw = plotWidth(for: size)
        let ph = plotHeight(for: size)
        let normX = (alpha - minAlpha) / rangeAlpha
        let tLambda = transformLambda(rawLambda)
        let normY = (tLambda - minL) / rangeL
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
        let acPos = effectiveAlphaCPos
        let acNeg = effectiveAlphaCNeg
        
        ZStack(alignment: .topLeading) {
            Canvas { ctx, cSize in
                drawChart(ctx: &ctx, size: cSize)
            }
            
            // Top Controls & Mode Legend Bar
            VStack {
                ScrollView(.horizontal, showsIndicators: false) {
                    HStack(spacing: 10) {
                        Picker("", selection: $scaleMode) {
                            Text("Normalized (±1)").tag(0)
                            Text("Zoom Stable").tag(1)
                            Text("Raw λ").tag(2)
                        }
                        .labelsHidden()
                        .pickerStyle(.segmented)
                        .controlSize(.small)
                        .fixedSize()
                        
                        HStack(spacing: 4) {
                            Circle().fill(Color.green).frame(width: 7, height: 7)
                            Text("λ_min > 0 (Stable)")
                                .font(.caption2.bold())
                                .foregroundColor(.green)
                        }
                        .padding(.horizontal, 6)
                        .padding(.vertical, 2)
                        .background(Color.green.opacity(0.12), in: RoundedRectangle(cornerRadius: 4))
                        .fixedSize()
                        
                        HStack(spacing: 4) {
                            Circle().fill(Color.red).frame(width: 7, height: 7)
                            Text("λ_min < 0 (Unstable)")
                                .font(.caption2.bold())
                                .foregroundColor(.red)
                        }
                        .padding(.horizontal, 6)
                        .padding(.vertical, 2)
                        .background(Color.red.opacity(0.12), in: RoundedRectangle(cornerRadius: 4))
                        .fixedSize()
                        
                        if showHigherModes {
                            HStack(spacing: 5) {
                                Circle()
                                    .stroke(Color.orange, lineWidth: 2)
                                    .frame(width: 8, height: 8)
                                Text(isLambda2DegenerateWithMin ? "λ₂ (= λ_min doublet)" : "λ₂ (2nd mode)")
                                    .font(.caption2.bold())
                                    .foregroundColor(.orange)
                            }
                            .padding(.horizontal, 6)
                            .padding(.vertical, 2)
                            .background(Color.orange.opacity(0.14), in: RoundedRectangle(cornerRadius: 4))
                            .fixedSize()
                            
                            HStack(spacing: 5) {
                                Rectangle()
                                    .fill(Color.cyan)
                                    .frame(width: 7, height: 7)
                                    .rotationEffect(.degrees(45))
                                Text("λ₃ (3rd mode)")
                                    .font(.caption2.bold())
                                    .foregroundColor(.cyan)
                            }
                            .padding(.horizontal, 6)
                            .padding(.vertical, 2)
                            .background(Color.cyan.opacity(0.14), in: RoundedRectangle(cornerRadius: 4))
                            .fixedSize()
                            
                            HStack(spacing: 5) {
                                Circle()
                                    .stroke(Color.purple, lineWidth: 2)
                                    .frame(width: 8, height: 8)
                                Text(isLambda4DegenerateWith3 ? "λ₄ (= λ₃ doublet)" : "λ₄ (4th mode)")
                                    .font(.caption2.bold())
                                    .foregroundColor(.purple)
                            }
                            .padding(.horizontal, 6)
                            .padding(.vertical, 2)
                            .background(Color.purple.opacity(0.14), in: RoundedRectangle(cornerRadius: 4))
                            .fixedSize()
                            
                            HStack(spacing: 5) {
                                Circle()
                                    .fill(Color.mint)
                                    .frame(width: 7, height: 7)
                                Text("λ₅ (5th mode)")
                                    .font(.caption2.bold())
                                    .foregroundColor(.mint)
                            }
                            .padding(.horizontal, 6)
                            .padding(.vertical, 2)
                            .background(Color.mint.opacity(0.14), in: RoundedRectangle(cornerRadius: 4))
                            .fixedSize()
                        }
                    }
                }
                .padding(.leading, padLeft + 4)
                .padding(.trailing, padRight + 4)
                .padding(.top, 4)
                Spacer()
            }
            
            // Critical Alpha Labels
            if acPos != 0.0 {
                let ptPos = toScreen(alpha: acPos, rawLambda: 0, size: size)
                Text(String(format: "+α_c = %.3f", acPos))
                    .font(.system(size: 10, weight: .bold, design: .monospaced))
                    .foregroundColor(.red)
                    .padding(.horizontal, 4)
                    .padding(.vertical, 2)
                    .background(Color(NSColor.windowBackgroundColor).opacity(0.85), in: RoundedRectangle(cornerRadius: 3))
                    .position(x: min(padLeft + pw - 35, ptPos.x + 36), y: max(padTop + 14, ptPos.y - 12))
            }
            if acNeg != 0.0 {
                let ptNeg = toScreen(alpha: acNeg, rawLambda: 0, size: size)
                Text(String(format: "-α_c = %.3f", acNeg))
                    .font(.system(size: 10, weight: .bold, design: .monospaced))
                    .foregroundColor(.red)
                    .padding(.horizontal, 4)
                    .padding(.vertical, 2)
                    .background(Color(NSColor.windowBackgroundColor).opacity(0.85), in: RoundedRectangle(cornerRadius: 3))
                    .position(x: max(padLeft + 35, ptNeg.x - 36), y: max(padTop + 14, ptNeg.y - 12))
            }
            
            // Hover Tooltip Overlay
            if let hovered = hoveredPoint {
                let scr = toScreen(alpha: hovered.alpha, rawLambda: hovered.lambda_min, size: size)
                tooltipView(for: hovered)
                    .position(
                        x: min(max(scr.x, padLeft + 80), padLeft + pw - 80),
                        y: max(min(scr.y - 55, size.height - padBottom - 45), padTop + 45)
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
            HStack(spacing: 6) {
                Text(String(format: "α = %.4f", hovered.alpha))
                    .font(.caption.bold())
                Text(hovered.is_stable ? "STABLE" : "UNSTABLE")
                    .font(.system(size: 9, weight: .bold, design: .monospaced))
                    .padding(.horizontal, 4)
                    .padding(.vertical, 1)
                    .background((hovered.is_stable ? Color.green : Color.red).opacity(0.2), in: RoundedRectangle(cornerRadius: 3))
                    .foregroundColor(hovered.is_stable ? .green : .red)
            }
            Text(String(format: "λ_min = %.4e", hovered.lambda_min))
                .font(.caption2.monospaced())
                .foregroundColor(hovered.is_stable ? .green : .red)
            if showHigherModes {
                Text(String(format: "λ_2   = %.4e", hovered.lambda_2))
                    .font(.caption2.monospaced())
                    .foregroundColor(.orange)
                Text(String(format: "λ_3   = %.4e", hovered.lambda_3))
                    .font(.caption2.monospaced())
                    .foregroundColor(.cyan)
                Text(String(format: "λ_4   = %.4e", hovered.lambda_4))
                    .font(.caption2.monospaced())
                    .foregroundColor(.purple)
                Text(String(format: "λ_5   = %.4e", hovered.lambda_5))
                    .font(.caption2.monospaced())
                    .foregroundColor(.mint)
            }
            Text(String(format: "E_FEM = %.4e", hovered.energy_density))
                .font(.caption2.monospaced())
                .foregroundColor(.secondary)
        }
        .padding(8)
        .background(Color(NSColor.controlBackgroundColor).opacity(0.95), in: RoundedRectangle(cornerRadius: 6))
        .shadow(radius: 3)
    }
    
    private func formatLambdaLabel(_ rawVal: Double) -> String {
        if abs(rawVal) < 1e-14 { return "0" }
        if abs(rawVal) >= 0.1 && abs(rawVal) < 1000 {
            return String(format: "%.2f", rawVal)
        }
        return String(format: "%.2e", rawVal)
    }
    
    private func drawChart(ctx: inout GraphicsContext, size: CGSize) {
        let pw = plotWidth(for: size)
        let ph = plotHeight(for: size)
        let zeroY = toScreen(alpha: 0, rawLambda: 0, size: size).y
        let acNeg = effectiveAlphaCNeg
        let acPos = effectiveAlphaCPos
        
        // 1. Shaded Stability Background Regions
        // Horizontal split at λ = 0 AND vertical interval shading between [-α_c, +α_c]
        if acNeg < 0 && acPos > 0 {
            let xNeg = toScreen(alpha: acNeg, rawLambda: 0, size: size).x
            let xPos = toScreen(alpha: acPos, rawLambda: 0, size: size).x
            
            // Left Unstable Band (α < -α_c)
            if xNeg > padLeft {
                let rectLeft = CGRect(x: padLeft, y: padTop, width: xNeg - padLeft, height: ph)
                ctx.fill(Path(rectLeft), with: .color(Color.red.opacity(0.08)))
            }
            // Center Stable Band (-α_c <= α <= +α_c)
            if xPos > xNeg {
                let rectCenter = CGRect(x: max(padLeft, xNeg), y: padTop, width: min(padLeft + pw, xPos) - max(padLeft, xNeg), height: ph)
                ctx.fill(Path(rectCenter), with: .color(Color.green.opacity(0.07)))
            }
            // Right Unstable Band (α > +α_c)
            if xPos < padLeft + pw {
                let rectRight = CGRect(x: xPos, y: padTop, width: padLeft + pw - xPos, height: ph)
                ctx.fill(Path(rectRight), with: .color(Color.red.opacity(0.08)))
            }
        } else {
            if zeroY < padTop + ph {
                let rectUnstable = CGRect(x: padLeft, y: max(padTop, zeroY), width: pw, height: padTop + ph - max(padTop, zeroY))
                ctx.fill(Path(rectUnstable), with: .color(Color.red.opacity(0.08)))
            }
            if zeroY > padTop {
                let rectStable = CGRect(x: padLeft, y: padTop, width: pw, height: min(padTop + ph, zeroY) - padTop)
                ctx.fill(Path(rectStable), with: .color(Color.green.opacity(0.06)))
            }
        }
        
        // 2. Y-Axis Gridlines & Labels
        if scaleMode == 0 {
            // Normalized split scale ticks: +1.0 (posScale), +0.5, 0, -0.5, -1.0 (-negScale)
            let ticks: [(normYVal: Double, rawEquiv: Double, label: String)] = [
                (1.0, posScale, "+\(formatLambdaLabel(posScale))"),
                (0.5, 0.5 * posScale, "+\(formatLambdaLabel(0.5 * posScale))"),
                (0.0, 0.0, "0 (Threshold)"),
                (-0.5, -0.5 * negScale, "-\(formatLambdaLabel(0.5 * negScale))"),
                (-1.0, -negScale, "-\(formatLambdaLabel(negScale))")
            ]
            for tick in ticks where tick.normYVal >= minL && tick.normYVal <= maxL {
                let y = padTop + CGFloat(1.0 - (tick.normYVal - minL) / rangeL) * ph
                var gPath = Path()
                gPath.move(to: CGPoint(x: padLeft, y: y))
                gPath.addLine(to: CGPoint(x: padLeft + pw, y: y))
                ctx.stroke(gPath, with: .color(Color.primary.opacity(tick.normYVal == 0 ? 0.55 : 0.12)),
                           style: StrokeStyle(lineWidth: tick.normYVal == 0 ? 1.5 : 0.8, dash: tick.normYVal == 0 ? [5, 3] : [2, 3]))
                
                let txt = Text(tick.label)
                    .font(.system(size: 9, weight: tick.normYVal == 0 ? .bold : .regular, design: .monospaced))
                    .foregroundColor(tick.normYVal > 0 ? .green : (tick.normYVal < 0 ? .red : .primary))
                ctx.draw(txt, at: CGPoint(x: padLeft - 6, y: y), anchor: .trailing)
            }
        } else {
            let nTicks = 5
            for i in 0...nTicks {
                let frac = Double(i) / Double(nTicks)
                let tVal = minL + frac * rangeL
                let rawVal = (scaleMode == 1) ? (tVal * posScale) : tVal
                let y = padTop + CGFloat(1.0 - frac) * ph
                var gPath = Path()
                gPath.move(to: CGPoint(x: padLeft, y: y))
                gPath.addLine(to: CGPoint(x: padLeft + pw, y: y))
                ctx.stroke(gPath, with: .color(Color.primary.opacity(0.12)), style: StrokeStyle(lineWidth: 0.8, dash: [2, 3]))
                
                let txt = Text(formatLambdaLabel(rawVal))
                    .font(.system(size: 9, design: .monospaced))
                    .foregroundColor(.secondary)
                ctx.draw(txt, at: CGPoint(x: padLeft - 6, y: y), anchor: .trailing)
            }
            if zeroY >= padTop && zeroY <= padTop + ph {
                var zeroPath = Path()
                zeroPath.move(to: CGPoint(x: padLeft, y: zeroY))
                zeroPath.addLine(to: CGPoint(x: padLeft + pw, y: zeroY))
                ctx.stroke(zeroPath, with: .color(Color.primary.opacity(0.6)), style: StrokeStyle(lineWidth: 1.5, dash: [5, 3]))
            }
        }
        
        // 3. X-Axis Ticks & Vertical Reference Lines
        let xTicks = 6
        for i in 0...xTicks {
            let frac = Double(i) / Double(xTicks)
            let aVal = minAlpha + frac * rangeAlpha
            let x = padLeft + CGFloat(frac) * pw
            let txt = Text(String(format: "%.2f", aVal))
                .font(.system(size: 9, design: .monospaced))
                .foregroundColor(.secondary)
            ctx.draw(txt, at: CGPoint(x: x, y: padTop + ph + 14), anchor: .center)
        }
        let xAxisTitle = Text("Simple Shear Amplitude α")
            .font(.system(size: 10, weight: .medium))
            .foregroundColor(.secondary)
        ctx.draw(xAxisTitle, at: CGPoint(x: padLeft + pw / 2, y: padTop + ph + 32), anchor: .center)
        
        let zeroX = toScreen(alpha: 0, rawLambda: 0, size: size).x
        if zeroX >= padLeft && zeroX <= padLeft + pw {
            var zeroXAxis = Path()
            zeroXAxis.move(to: CGPoint(x: zeroX, y: padTop))
            zeroXAxis.addLine(to: CGPoint(x: zeroX, y: padTop + ph))
            ctx.stroke(zeroXAxis, with: .color(Color.primary.opacity(0.28)), style: StrokeStyle(lineWidth: 1.0, dash: [2, 2]))
        }
        
        // Critical Alpha Marker Lines
        if acNeg != 0.0 {
            let ptNeg = toScreen(alpha: acNeg, rawLambda: 0, size: size)
            var pLine = Path()
            pLine.move(to: CGPoint(x: ptNeg.x, y: padTop))
            pLine.addLine(to: CGPoint(x: ptNeg.x, y: padTop + ph))
            ctx.stroke(pLine, with: .color(Color.red.opacity(0.75)), style: StrokeStyle(lineWidth: 1.2, dash: [4, 4]))
        }
        if acPos != 0.0 {
            let ptPos = toScreen(alpha: acPos, rawLambda: 0, size: size)
            var pLine = Path()
            pLine.move(to: CGPoint(x: ptPos.x, y: padTop))
            pLine.addLine(to: CGPoint(x: ptPos.x, y: padTop + ph))
            ctx.stroke(pLine, with: .color(Color.red.opacity(0.75)), style: StrokeStyle(lineWidth: 1.2, dash: [4, 4]))
        }
        
        // 4. Main Curve: λ_min(α) colored by stability segment
        if points.count > 1 {
            for i in 1..<points.count {
                let p0 = points[i - 1]
                let p1 = points[i]
                let s0 = toScreen(alpha: p0.alpha, rawLambda: p0.lambda_min, size: size)
                let s1 = toScreen(alpha: p1.alpha, rawLambda: p1.lambda_min, size: size)
                var seg = Path()
                seg.move(to: s0)
                seg.addLine(to: s1)
                let segStable = p0.is_stable && p1.is_stable
                ctx.stroke(seg, with: .color(segStable ? Color.green : Color.red), lineWidth: 2.4)
            }
        }
        
        // Plot Points for λ_min(α)
        for pt in points {
            let scr = toScreen(alpha: pt.alpha, rawLambda: pt.lambda_min, size: size)
            let ptColor = pt.is_stable ? Color.green : Color.red
            let circleRect = CGRect(x: scr.x - 3.2, y: scr.y - 3.2, width: 6.4, height: 6.4)
            ctx.fill(Path(ellipseIn: circleRect), with: .color(ptColor))
        }
        
        // 5. Higher modes λ_2 ... λ_5 drawn ON TOP with distinct styles and markers
        if showHigherModes {
            var p5Path = Path()
            var p4Path = Path()
            var p3Path = Path()
            var p2Path = Path()
            for (idx, pt) in points.enumerated() {
                let scr2 = toScreen(alpha: pt.alpha, rawLambda: pt.lambda_2, size: size)
                let scr3 = toScreen(alpha: pt.alpha, rawLambda: pt.lambda_3, size: size)
                let scr4 = toScreen(alpha: pt.alpha, rawLambda: pt.lambda_4, size: size)
                let scr5 = toScreen(alpha: pt.alpha, rawLambda: pt.lambda_5, size: size)
                if idx == 0 {
                    p2Path.move(to: scr2)
                    p3Path.move(to: scr3)
                    p4Path.move(to: scr4)
                    p5Path.move(to: scr5)
                } else {
                    p2Path.addLine(to: scr2)
                    p3Path.addLine(to: scr3)
                    p4Path.addLine(to: scr4)
                    p5Path.addLine(to: scr5)
                }
            }
            // Draw λ_5 in mint with round dot markers
            ctx.stroke(p5Path, with: .color(Color.mint), style: StrokeStyle(lineWidth: 1.8, dash: [4, 3]))
            for (idx, pt) in points.enumerated() where idx % 2 == 0 {
                let scr5 = toScreen(alpha: pt.alpha, rawLambda: pt.lambda_5, size: size)
                let r: CGFloat = 2.8
                ctx.fill(Path(ellipseIn: CGRect(x: scr5.x - r, y: scr5.y - r, width: 2 * r, height: 2 * r)), with: .color(Color.mint))
            }

            // Draw λ_4 in purple with square markers
            ctx.stroke(p4Path, with: .color(Color.purple), style: StrokeStyle(lineWidth: 1.8, dash: [6, 3]))
            for (idx, pt) in points.enumerated() where idx % 2 == 0 {
                let scr4 = toScreen(alpha: pt.alpha, rawLambda: pt.lambda_4, size: size)
                let r: CGFloat = 2.8
                ctx.stroke(Path(roundedRect: CGRect(x: scr4.x - r, y: scr4.y - r, width: 2 * r, height: 2 * r), cornerRadius: 1), with: .color(Color.purple), lineWidth: 1.5)
            }

            // Draw λ_3 in bright cyan with diamond markers
            ctx.stroke(p3Path, with: .color(Color.cyan), style: StrokeStyle(lineWidth: 2.0))
            for (idx, pt) in points.enumerated() where idx % 2 == 0 || pt.is_stable {
                let scr3 = toScreen(alpha: pt.alpha, rawLambda: pt.lambda_3, size: size)
                var diamond = Path()
                let r: CGFloat = 3.6
                diamond.move(to: CGPoint(x: scr3.x, y: scr3.y - r))
                diamond.addLine(to: CGPoint(x: scr3.x + r, y: scr3.y))
                diamond.addLine(to: CGPoint(x: scr3.x, y: scr3.y + r))
                diamond.addLine(to: CGPoint(x: scr3.x - r, y: scr3.y))
                diamond.closeSubpath()
                ctx.fill(diamond, with: .color(Color.cyan))
            }
            
            // Draw λ_2 in dashed orange with hollow ring markers (visible even when degenerate with λ_min!)
            ctx.stroke(p2Path, with: .color(Color.orange), style: StrokeStyle(lineWidth: 2.0, dash: [5, 4]))
            for (idx, pt) in points.enumerated() where idx % 2 == 0 {
                let scr2 = toScreen(alpha: pt.alpha, rawLambda: pt.lambda_2, size: size)
                let ringRect = CGRect(x: scr2.x - 4.2, y: scr2.y - 4.2, width: 8.4, height: 8.4)
                ctx.stroke(Path(ellipseIn: ringRect), with: .color(Color.orange), lineWidth: 1.5)
            }
        }
    }
}

// MARK: - Interactive Canvas Chart for Strain Energy E(α) vs Continuum Theory
struct StiffnessEnergyChartView: View {
    let points: [StiffnessPointModel]
    let summary: StiffnessStudySummaryModel
    @Binding var hoveredPoint: StiffnessPointModel?
    @State private var energyMode: Int = 0 // 0: Density (per unit area), 1: Total Integrated
    
    private let padLeft: CGFloat = 78
    private let padRight: CGFloat = 30
    private let padTop: CGFloat = 40
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
    
    private var effectiveAlphaCNeg: Double {
        if summary.alpha_c_neg != 0.0 { return summary.alpha_c_neg }
        for i in 1..<points.count {
            let p0 = points[i - 1], p1 = points[i]
            if (p0.lambda_min <= 0 && p1.lambda_min > 0) || (p0.lambda_min > 0 && p1.lambda_min <= 0) {
                let frac = abs(p0.lambda_min) / max(1e-15, abs(p0.lambda_min) + abs(p1.lambda_min))
                let cross = p0.alpha + frac * (p1.alpha - p0.alpha)
                if cross < 0 { return cross }
            }
        }
        return 0.0
    }
    
    private var effectiveAlphaCPos: Double {
        if summary.alpha_c_pos != 0.0 { return summary.alpha_c_pos }
        for i in 1..<points.count {
            let p0 = points[i - 1], p1 = points[i]
            if (p0.lambda_min <= 0 && p1.lambda_min > 0) || (p0.lambda_min > 0 && p1.lambda_min <= 0) {
                let frac = abs(p0.lambda_min) / max(1e-15, abs(p0.lambda_min) + abs(p1.lambda_min))
                let cross = p0.alpha + frac * (p1.alpha - p0.alpha)
                if cross > 0 { return cross }
            }
        }
        return 0.0
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
        let acPos = effectiveAlphaCPos
        let acNeg = effectiveAlphaCNeg
        
        ZStack(alignment: .topLeading) {
            Canvas { ctx, cSize in
                drawChart(ctx: &ctx, size: cSize)
            }
            
            // Top Legend & Energy Mode Picker
            VStack {
                ScrollView(.horizontal, showsIndicators: false) {
                    HStack(spacing: 10) {
                        Picker("", selection: $energyMode) {
                            Text("Per Triangle").tag(0)
                            Text("Total Integrated").tag(1)
                        }
                        .labelsHidden()
                        .pickerStyle(.segmented)
                        .controlSize(.small)
                        .fixedSize()
                        
                        HStack(spacing: 4) {
                            Circle().fill(Color.green).frame(width: 8, height: 8)
                            Text("Stable FEM (λ_min > 0)")
                                .font(.caption2.bold())
                                .foregroundColor(.green)
                        }
                        .padding(.horizontal, 6)
                        .padding(.vertical, 2)
                        .background(Color.green.opacity(0.12), in: RoundedRectangle(cornerRadius: 4))
                        .fixedSize()
                        
                        HStack(spacing: 4) {
                            Circle().fill(Color.red).frame(width: 8, height: 8)
                            Text("Unstable FEM (λ_min < 0)")
                                .font(.caption2.bold())
                                .foregroundColor(.red)
                        }
                        .padding(.horizontal, 6)
                        .padding(.vertical, 2)
                        .background(Color.red.opacity(0.12), in: RoundedRectangle(cornerRadius: 4))
                        .fixedSize()
                        
                        HStack(spacing: 5) {
                            Rectangle().fill(Color.orange).frame(width: 14, height: 2)
                            Text("Continuum W(C)")
                                .font(.caption2.bold())
                                .foregroundColor(.orange)
                        }
                        .padding(.horizontal, 6)
                        .padding(.vertical, 2)
                        .background(Color.orange.opacity(0.12), in: RoundedRectangle(cornerRadius: 4))
                        .fixedSize()
                    }
                }
                .padding(.leading, padLeft + 4)
                .padding(.trailing, padRight + 4)
                .padding(.top, 4)
                
                Spacer()
            }
            
            // Critical Alpha Labels
            if acPos != 0.0 {
                let ptPos = toScreen(alpha: acPos, energy: minE, size: size)
                Text(String(format: "+α_c = %.3f", acPos))
                    .font(.system(size: 10, weight: .bold, design: .monospaced))
                    .foregroundColor(.red)
                    .padding(.horizontal, 4)
                    .padding(.vertical, 2)
                    .background(Color(NSColor.windowBackgroundColor).opacity(0.85), in: RoundedRectangle(cornerRadius: 3))
                    .position(x: min(padLeft + pw - 35, ptPos.x + 36), y: padTop + 18)
            }
            if acNeg != 0.0 {
                let ptNeg = toScreen(alpha: acNeg, energy: minE, size: size)
                Text(String(format: "-α_c = %.3f", acNeg))
                    .font(.system(size: 10, weight: .bold, design: .monospaced))
                    .foregroundColor(.red)
                    .padding(.horizontal, 4)
                    .padding(.vertical, 2)
                    .background(Color(NSColor.windowBackgroundColor).opacity(0.85), in: RoundedRectangle(cornerRadius: 3))
                    .position(x: max(padLeft + 35, ptNeg.x - 36), y: padTop + 18)
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
            HStack(spacing: 6) {
                Text(String(format: "α = %.4f", hovered.alpha))
                    .font(.caption.bold())
                Text(hovered.is_stable ? "STABLE" : "UNSTABLE")
                    .font(.system(size: 9, weight: .bold, design: .monospaced))
                    .padding(.horizontal, 4)
                    .padding(.vertical, 1)
                    .background((hovered.is_stable ? Color.green : Color.red).opacity(0.2), in: RoundedRectangle(cornerRadius: 3))
                    .foregroundColor(hovered.is_stable ? .green : .red)
            }
            Text(String(format: "E_FEM  = %.5e", eFem))
                .font(.caption2.monospaced())
                .foregroundColor(hovered.is_stable ? .green : .red)
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
        let acNeg = effectiveAlphaCNeg
        let acPos = effectiveAlphaCPos
        
        // 1. Shaded Stability Intervals along α
        if acNeg < 0 && acPos > 0 {
            let xNeg = toScreen(alpha: acNeg, energy: 0, size: size).x
            let xPos = toScreen(alpha: acPos, energy: 0, size: size).x
            
            // Left Unstable Band (α < -α_c)
            if xNeg > padLeft {
                let rectLeft = CGRect(x: padLeft, y: padTop, width: xNeg - padLeft, height: ph)
                ctx.fill(Path(rectLeft), with: .color(Color.red.opacity(0.09)))
            }
            // Center Stable Band (-α_c <= α <= +α_c)
            if xPos > xNeg {
                let rectCenter = CGRect(x: max(padLeft, xNeg), y: padTop, width: min(padLeft + pw, xPos) - max(padLeft, xNeg), height: ph)
                ctx.fill(Path(rectCenter), with: .color(Color.green.opacity(0.08)))
            }
            // Right Unstable Band (α > +α_c)
            if xPos < padLeft + pw {
                let rectRight = CGRect(x: xPos, y: padTop, width: padLeft + pw - xPos, height: ph)
                ctx.fill(Path(rectRight), with: .color(Color.red.opacity(0.09)))
            }
        }
        
        // 2. Y-Axis Gridlines & Labels
        let nTicks = 5
        for i in 0...nTicks {
            let frac = Double(i) / Double(nTicks)
            let eVal = minE + frac * rangeE
            let y = padTop + CGFloat(1.0 - frac) * ph
            var gPath = Path()
            gPath.move(to: CGPoint(x: padLeft, y: y))
            gPath.addLine(to: CGPoint(x: padLeft + pw, y: y))
            ctx.stroke(gPath, with: .color(Color.primary.opacity(0.12)), style: StrokeStyle(lineWidth: 0.8, dash: [2, 3]))
            
            let labelStr = (maxE >= 100.0) ? String(format: "%.1e", eVal) : String(format: "%.3f", eVal)
            let txt = Text(labelStr)
                .font(.system(size: 9, design: .monospaced))
                .foregroundColor(.secondary)
            ctx.draw(txt, at: CGPoint(x: padLeft - 6, y: y), anchor: .trailing)
        }
        
        // Zero Axis Line
        if zeroY >= padTop && zeroY <= padTop + ph {
            var zeroPath = Path()
            zeroPath.move(to: CGPoint(x: padLeft, y: zeroY))
            zeroPath.addLine(to: CGPoint(x: padLeft + pw, y: zeroY))
            ctx.stroke(zeroPath, with: .color(Color.primary.opacity(0.4)), style: StrokeStyle(lineWidth: 1.0, dash: [4, 4]))
        }
        
        // X-Axis Ticks
        let xTicks = 6
        for i in 0...xTicks {
            let frac = Double(i) / Double(xTicks)
            let aVal = minAlpha + frac * rangeAlpha
            let x = padLeft + CGFloat(frac) * pw
            let txt = Text(String(format: "%.2f", aVal))
                .font(.system(size: 9, design: .monospaced))
                .foregroundColor(.secondary)
            ctx.draw(txt, at: CGPoint(x: x, y: padTop + ph + 14), anchor: .center)
        }
        let xAxisTitle = Text("Simple Shear Amplitude α")
            .font(.system(size: 10, weight: .medium))
            .foregroundColor(.secondary)
        ctx.draw(xAxisTitle, at: CGPoint(x: padLeft + pw / 2, y: padTop + ph + 32), anchor: .center)
        
        // Vertical Axis at α = 0
        let zeroX = toScreen(alpha: 0, energy: 0, size: size).x
        if zeroX >= padLeft && zeroX <= padLeft + pw {
            var zeroXAxis = Path()
            zeroXAxis.move(to: CGPoint(x: zeroX, y: padTop))
            zeroXAxis.addLine(to: CGPoint(x: zeroX, y: padTop + ph))
            ctx.stroke(zeroXAxis, with: .color(Color.primary.opacity(0.3)), style: StrokeStyle(lineWidth: 1.0, dash: [2, 2]))
        }
        
        // Critical Alpha Marker Lines
        if acNeg != 0.0 {
            let ptNeg = toScreen(alpha: acNeg, energy: 0, size: size)
            var pLine = Path()
            pLine.move(to: CGPoint(x: ptNeg.x, y: padTop))
            pLine.addLine(to: CGPoint(x: ptNeg.x, y: padTop + ph))
            ctx.stroke(pLine, with: .color(Color.red.opacity(0.75)), style: StrokeStyle(lineWidth: 1.2, dash: [4, 4]))
        }
        if acPos != 0.0 {
            let ptPos = toScreen(alpha: acPos, energy: 0, size: size)
            var pLine = Path()
            pLine.move(to: CGPoint(x: ptPos.x, y: padTop))
            pLine.addLine(to: CGPoint(x: ptPos.x, y: padTop + ph))
            ctx.stroke(pLine, with: .color(Color.red.opacity(0.75)), style: StrokeStyle(lineWidth: 1.2, dash: [4, 4]))
        }
        
        // 3. Analytical Continuum Energy Curve (dashed orange)
        var contPath = Path()
        for (idx, pt) in points.enumerated() {
            let scr = toScreen(alpha: pt.alpha, energy: getContEnergy(pt), size: size)
            if idx == 0 {
                contPath.move(to: scr)
            } else {
                contPath.addLine(to: scr)
            }
        }
        ctx.stroke(contPath, with: .color(Color.orange.opacity(0.9)), style: StrokeStyle(lineWidth: 2.0, dash: [6, 4]))
        
        // 4. Simulated Lattice Energy Curve colored by stability (Green = Stable, Red = Unstable)
        if points.count > 1 {
            for i in 1..<points.count {
                let p0 = points[i - 1]
                let p1 = points[i]
                let s0 = toScreen(alpha: p0.alpha, energy: getFEMEnergy(p0), size: size)
                let s1 = toScreen(alpha: p1.alpha, energy: getFEMEnergy(p1), size: size)
                var seg = Path()
                seg.move(to: s0)
                seg.addLine(to: s1)
                let segStable = p0.is_stable && p1.is_stable
                ctx.stroke(seg, with: .color(segStable ? Color.green : Color.red), lineWidth: 2.5)
            }
        }
        
        // 5. Points on Simulated Curve colored by stability
        for pt in points {
            let scr = toScreen(alpha: pt.alpha, energy: getFEMEnergy(pt), size: size)
            let circleRect = CGRect(x: scr.x - 3.2, y: scr.y - 3.2, width: 6.4, height: 6.4)
            ctx.fill(Path(ellipseIn: circleRect), with: .color(pt.is_stable ? Color.green : Color.red))
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
