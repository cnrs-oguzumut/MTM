import SwiftUI
import Combine
import AppKit
import UniformTypeIdentifiers

public enum StiffnessEigenmodeChoice: Int, CaseIterable, Identifiable {
    case mode1 = 1
    case mode2 = 2
    case mode3 = 3
    case mode4 = 4
    case mode5 = 5
    
    public var id: Int { rawValue }
    public var label: String {
        switch self {
        case .mode1: return "Mode 1 (λ₁)"
        case .mode2: return "Mode 2 (λ₂)"
        case .mode3: return "Mode 3 (λ₃)"
        case .mode4: return "Mode 4 (λ₄)"
        case .mode5: return "Mode 5 (λ₅)"
        }
    }
    
    public func dynamicLabel(isUniqueLevels: Bool) -> String {
        let prefix = isUniqueLevels ? "Level" : "Mode"
        switch self {
        case .mode1: return "\(prefix) 1 (λ₁)"
        case .mode2: return "\(prefix) 2 (λ₂)"
        case .mode3: return "\(prefix) 3 (λ₃)"
        case .mode4: return "\(prefix) 4 (λ₄)"
        case .mode5: return "\(prefix) 5 (λ₅)"
        }
    }
}

public enum StiffnessColorField: String, CaseIterable, Identifiable {
    case magnitude = "|δu| Mag"
    case elementShear = "Elem Shear"
    case direction = "Phase arg(δu)"
    
    public var id: String { rawValue }
}

public enum StiffnessColorPalette: String, CaseIterable, Identifiable {
    case coolwarm = "Cool-Warm"
    case viridis = "Viridis"
    case plasma = "Plasma"
    case grayscale = "Grayscale"
    
    public var id: String { rawValue }
}

public struct StiffnessModeInspectorView: View {
    @ObservedObject var manager: StiffnessStudyManager
    
    @State private var selectedMode: StiffnessEigenmodeChoice = .mode1
    @State private var frameMode: AcousticFrameMode = .eulerian
    @State private var colorField: StiffnessColorField = .magnitude
    @State private var colorPalette: StiffnessColorPalette = .coolwarm
    @State private var nodeScale: Double = 1.0
    @State private var toastMessage: String? = nil
    @State private var modeAmplitude: Double = 0.35
    @State private var showQuiver: Bool = true
    @State private var showMeshEdges: Bool = true
    @State private var isAnimating: Bool = false
    @State private var animPhase: Double = 0.0
    @State private var zoomScale: Double = 1.0
    @State private var panOffset: CGSize = .zero
    @State private var dragStartOffset: CGSize = .zero
    @State private var hoveredAngleDeg: Double? = nil
    
    private let animTimer = Timer.publish(every: 1.0 / 30.0, on: .main, in: .common).autoconnect()
    
    public init(manager: StiffnessStudyManager) {
        self.manager = manager
    }
    
    private var currentPoint: StiffnessPointModel? {
        let pts = manager.points
        guard !pts.isEmpty else { return nil }
        if let sel = manager.selectedAlpha {
            return pts.min(by: { abs($0.alpha - sel) < abs($1.alpha - sel) })
        }
        return pts.first(where: { !$0.is_stable }) ?? pts[pts.count / 2]
    }
    
    private var alphaRange: ClosedRange<Double> {
        guard let first = manager.points.first?.alpha,
              let last = manager.points.last?.alpha,
              last > first else {
            return -0.3...0.3
        }
        return first...last
    }
    
    private func stepAlpha(by deltaIndex: Int) {
        let pts = manager.points
        guard !pts.isEmpty else { return }
        let curAlpha = currentPoint?.alpha ?? pts[0].alpha
        if let idx = pts.enumerated().min(by: { abs($0.element.alpha - curAlpha) < abs($1.element.alpha - curAlpha) })?.offset {
            let nextIdx = max(0, min(pts.count - 1, idx + deltaIndex))
            manager.selectedAlpha = pts[nextIdx].alpha
        }
    }
    
    public var body: some View {
        VStack(spacing: 0) {
            // =========================================================================
            // 1. TOP RESPONSIVE CONTROL BAR (Horizontally scrollable, never clips)
            // =========================================================================
            ScrollView(.horizontal, showsIndicators: false) {
                HStack(spacing: 10) {
                    // A. Shear α Scrubber + Step Buttons
                    HStack(spacing: 5) {
                        Button(action: { stepAlpha(by: -1) }) {
                            Image(systemName: "chevron.left")
                                .font(.system(size: 10, weight: .bold))
                        }
                        .buttonStyle(.bordered)
                        .controlSize(.mini)
                        .help("Previous shear step α")
                        
                        Text("Scrub α:")
                            .font(.caption2)
                            .foregroundColor(.secondary)
                        
                        Slider(
                            value: Binding(
                                get: { currentPoint?.alpha ?? 0.0 },
                                set: { newVal in manager.selectedAlpha = newVal }
                            ),
                            in: alphaRange
                        )
                        .frame(width: 110)
                        .controlSize(.small)
                        
                        Button(action: { stepAlpha(by: 1) }) {
                            Image(systemName: "chevron.right")
                                .font(.system(size: 10, weight: .bold))
                        }
                        .buttonStyle(.bordered)
                        .controlSize(.mini)
                        .help("Next shear step α")
                    }
                    .fixedSize()
                    
                    // B. Critical α_c Quick Jump Pills
                    if manager.summary.alpha_c_neg < -1e-6 || manager.summary.alpha_c_pos > 1e-6 {
                        Divider().frame(height: 16)
                        HStack(spacing: 4) {
                            if manager.summary.alpha_c_neg < -1e-6 {
                                let acNeg = manager.summary.alpha_c_neg
                                let isSel = abs((currentPoint?.alpha ?? -999) - acNeg) < max(0.003, manager.dalpha * 0.75)
                                Button(action: { manager.selectedAlpha = acNeg }) {
                                    HStack(spacing: 3) {
                                        Circle().fill(Color.yellow).frame(width: 5, height: 5)
                                        Text(String(format: "α_c(-) = %.4f", acNeg))
                                            .font(.system(size: 10, weight: .semibold, design: .monospaced))
                                    }
                                    .padding(.horizontal, 6)
                                    .padding(.vertical, 3)
                                    .background(isSel ? Color.accentColor.opacity(0.25) : Color.secondary.opacity(0.12), in: Capsule())
                                }
                                .buttonStyle(.plain)
                            }
                            if manager.summary.alpha_c_pos > 1e-6 {
                                let acPos = manager.summary.alpha_c_pos
                                let isSel = abs((currentPoint?.alpha ?? -999) - acPos) < max(0.003, manager.dalpha * 0.75)
                                Button(action: { manager.selectedAlpha = acPos }) {
                                    HStack(spacing: 3) {
                                        Circle().fill(Color.yellow).frame(width: 5, height: 5)
                                        Text(String(format: "α_c(+) = +%.4f", acPos))
                                            .font(.system(size: 10, weight: .semibold, design: .monospaced))
                                    }
                                    .padding(.horizontal, 6)
                                    .padding(.vertical, 3)
                                    .background(isSel ? Color.accentColor.opacity(0.25) : Color.secondary.opacity(0.12), in: Capsule())
                                }
                                .buttonStyle(.plain)
                            }
                        }
                        .fixedSize()
                    }
                    
                    Divider().frame(height: 16)
                    
                    // C. Resolution Mode Picker: Subgrid (Live) vs Full N×N (Static)
                    if let pt = currentPoint {
                        if manager.isComputingFullScale {
                            HStack(spacing: 4) {
                                ProgressView().controlSize(.small)
                                Text("Solving \(manager.nx)×\(manager.ny)...")
                                    .font(.caption2.bold())
                                    .foregroundColor(.accentColor)
                            }
                            .fixedSize()
                        } else {
                            Button(action: { manager.computeFullScaleAtAlpha(pt.alpha) }) {
                                Label(String(format: "Compute %d×%d at α=%+.4f", manager.nx, manager.ny, pt.alpha),
                                      systemImage: "cpu")
                                    .font(.system(size: 10, weight: .semibold))
                            }
                            .buttonStyle(.borderedProminent)
                            .controlSize(.mini)
                            .fixedSize()
                        }
                    }
                    
                    Divider().frame(height: 16)
                    
                    // D1. Unique Levels Toggle (No Doublets)
                    Button(action: {
                        manager.uniqueLevelsOnly.toggle()
                        if let curAlpha = currentPoint?.alpha {
                            manager.computeFullScaleAtAlpha(curAlpha)
                        }
                    }) {
                        HStack(spacing: 4) {
                            Image(systemName: manager.uniqueLevelsOnly ? "checkmark.circle.fill" : "circle")
                                .font(.system(size: 11))
                                .foregroundColor(manager.uniqueLevelsOnly ? .accentColor : .secondary)
                            Text(manager.uniqueLevelsOnly ? "Unique Levels" : "All Modes")
                                .font(.caption2.weight(.medium))
                        }
                        .padding(.horizontal, 6)
                        .padding(.vertical, 3)
                        .background(manager.uniqueLevelsOnly ? Color.accentColor.opacity(0.18) : Color.secondary.opacity(0.1), in: RoundedRectangle(cornerRadius: 5))
                    }
                    .buttonStyle(.plain)
                    .help("Toggle between unique physical energy levels (skipping sin/cos doublets under PBC) and all raw eigenmodes")
                    .fixedSize()
                    
                    Divider().frame(height: 16)
                    
                    // D2. Eigenmode / Level Picker (Level 1 ... Level 5)
                    Picker("", selection: $selectedMode) {
                        ForEach(StiffnessEigenmodeChoice.allCases) { m in
                            Text(m.dynamicLabel(isUniqueLevels: manager.uniqueLevelsOnly)).tag(m)
                        }
                    }
                    .labelsHidden()
                    .pickerStyle(.segmented)
                    .controlSize(.small)
                    .fixedSize()
                    
                    Divider().frame(height: 16)
                    
                    // D. Eulerian vs Lagrangian Frame Picker
                    Picker("", selection: $frameMode) {
                        ForEach(AcousticFrameMode.allCases) { mode in
                            Text(mode.label).tag(mode)
                        }
                    }
                    .labelsHidden()
                    .pickerStyle(.segmented)
                    .controlSize(.small)
                    .fixedSize()
                    
                    Divider().frame(height: 16)
                    
                    // E. Color Field Picker
                    Picker("", selection: $colorField) {
                        ForEach(StiffnessColorField.allCases) { cf in
                            Text(cf.rawValue).tag(cf)
                        }
                    }
                    .labelsHidden()
                    .pickerStyle(.segmented)
                    .controlSize(.small)
                    .fixedSize()
                    
                    Divider().frame(height: 16)
                    
                    // F. Publication Color Palette Picker
                    Picker("", selection: $colorPalette) {
                        ForEach(StiffnessColorPalette.allCases) { cp in
                            Text(cp.rawValue).tag(cp)
                        }
                    }
                    .labelsHidden()
                    .pickerStyle(.segmented)
                    .controlSize(.small)
                    .fixedSize()
                    
                    Divider().frame(height: 16)
                    
                    // G. Publication Export Menu
                    Menu {
                        Section("Publication High-Res (300 DPI PNG)") {
                            Button(action: { exportFigure(target: .lattice, bg: .white) }) {
                                Label("Save Lattice Figure (White PNG)...", systemImage: "photo")
                            }
                            Button(action: { exportFigure(target: .rosette, bg: .white) }) {
                                Label("Save Polar Rosette (White PNG)...", systemImage: "scope")
                            }
                            Button(action: { exportFigure(target: .combined, bg: .white) }) {
                                Label("Save Dual Publication Figure (White PNG)...", systemImage: "rectangle.split.2x1")
                            }
                        }
                        Section("Screen Theme (Dark PNG)") {
                            Button(action: { exportFigure(target: .lattice, bg: .dark) }) {
                                Label("Save Lattice (Dark PNG)...", systemImage: "moon.fill")
                            }
                            Button(action: { exportFigure(target: .combined, bg: .dark) }) {
                                Label("Save Dual Figure (Dark PNG)...", systemImage: "moon.stars.fill")
                            }
                        }
                        Section("Clipboard (Direct Paste into Overleaf / Keynote)") {
                            Button(action: { copyFigureToClipboard(target: .lattice) }) {
                                Label("Copy Lattice Figure (White)", systemImage: "doc.on.clipboard")
                            }
                            Button(action: { copyFigureToClipboard(target: .combined) }) {
                                Label("Copy Dual Figure (White)", systemImage: "doc.on.clipboard.fill")
                            }
                        }
                    } label: {
                        Label("Export Figure", systemImage: "square.and.arrow.up")
                            .font(.caption.bold())
                    }
                    .menuStyle(.borderedButton)
                    .controlSize(.small)
                    .fixedSize()
                }
                .padding(.vertical, 1)
            }
            .frame(maxWidth: .infinity, alignment: .leading)
            .padding(.horizontal, 10)
            .padding(.vertical, 6)
            .background(BlueHour.controlSurface)
            
            Divider()
            
            // =========================================================================
            // 2. REAL-TIME TELEMETRY BAR (Horizontally scrollable, zero plot occlusion)
            // =========================================================================
            ScrollView(.horizontal, showsIndicators: false) {
                HStack(spacing: 8) {
                    if let pt = currentPoint {
                        let resolved = resolveLatticeMode(for: pt)
                        
                        HStack(spacing: 4) {
                            Text(String(format: "α = %+.4f", pt.alpha))
                                .font(.system(size: 11, weight: .bold, design: .monospaced))
                                .foregroundColor(.cyan)
                            
                            Text(pt.is_stable ? "STABLE" : "UNSTABLE")
                                .font(.system(size: 9, weight: .black))
                                .padding(.horizontal, 4)
                                .padding(.vertical, 1)
                                .background(pt.is_stable ? Color.green.opacity(0.25) : Color.red.opacity(0.25), in: Capsule())
                                .foregroundColor(pt.is_stable ? .green : .red)
                        }
                        .fixedSize()
                        
                        Divider().frame(height: 12)
                        
                        if manager.uniqueLevelsOnly {
                            HStack(spacing: 3) {
                                Image(systemName: "sparkles")
                                    .font(.system(size: 8))
                                Text("UNIQUE LEVELS (NO DOUBLETS)")
                                    .font(.system(size: 8.5, weight: .bold))
                            }
                            .padding(.horizontal, 5)
                            .padding(.vertical, 1.5)
                            .background(Color.accentColor.opacity(0.18), in: Capsule())
                            .foregroundColor(.accentColor)
                            .fixedSize()
                            
                            Divider().frame(height: 12)
                        }
                        
                        HStack(spacing: 3) {
                            Text("λ₁:").foregroundColor(.secondary)
                            Text(String(format: "%+.3e", pt.lambda_min)).bold()
                                .foregroundColor(pt.lambda_min >= 0 ? .green : .red)
                            Text("λ₂:").foregroundColor(.secondary)
                            Text(String(format: "%+.2e", pt.lambda_2))
                                .foregroundColor(.orange)
                            Text("λ₃:").foregroundColor(.secondary)
                            Text(String(format: "%+.2e", pt.lambda_3))
                                .foregroundColor(.cyan)
                            Text("λ₄:").foregroundColor(.secondary)
                            Text(String(format: "%+.2e", pt.lambda_4))
                                .foregroundColor(.purple)
                            Text("λ₅:").foregroundColor(.secondary)
                            Text(String(format: "%+.2e", pt.lambda_5))
                                .foregroundColor(.mint)
                        }
                        .font(.system(size: 10, design: .monospaced))
                        .fixedSize()
                        
                        Divider().frame(height: 12)
                        
                        HStack(spacing: 3) {
                            Text("Eulerian:")
                                .foregroundColor(.orange)
                                .bold()
                            Text(String(format: "n* = %.1f°", resolved.eulerNormalDeg))
                                .foregroundColor(.orange)
                            Text(String(format: "m* = (%+.2f, %+.2f) [%.1f°]", resolved.eulerMx, resolved.eulerMy, resolved.eulerPolDeg))
                                .foregroundColor(.yellow)
                        }
                        .font(.system(size: 10, design: .monospaced))
                        .padding(.horizontal, 6)
                        .padding(.vertical, 2)
                        .background(frameMode != .lagrangian ? Color.orange.opacity(0.14) : Color.clear, in: RoundedRectangle(cornerRadius: 4))
                        .fixedSize()
                        
                        Divider().frame(height: 12)
                        
                        HStack(spacing: 3) {
                            Text("Lagrangian:")
                                .foregroundColor(.purple)
                                .bold()
                            Text(String(format: "N* = %.1f°", resolved.lagrNormalDeg))
                                .foregroundColor(.purple)
                            Text(String(format: "M* = (%+.2f, %+.2f) [%.1f°]", resolved.lagrMx, resolved.lagrMy, resolved.lagrPolDeg))
                                .foregroundColor(.pink)
                        }
                        .font(.system(size: 10, design: .monospaced))
                        .padding(.horizontal, 6)
                        .padding(.vertical, 2)
                        .background(frameMode != .eulerian ? Color.purple.opacity(0.14) : Color.clear, in: RoundedRectangle(cornerRadius: 4))
                        .fixedSize()
                        
                        Divider().frame(height: 12)
                        
                        HStack(spacing: 4) {
                            Text(String(format: "IPR = %.2f (%@)", resolved.ipr, resolved.ipr > 2.2 ? "Localized Band" : "Extended Wave"))
                                .font(.system(size: 9.5, weight: .bold, design: .monospaced))
                                .padding(.horizontal, 5)
                                .padding(.vertical, 1)
                                .background(Color.cyan.opacity(0.18), in: Capsule())
                                .foregroundColor(.cyan)
                            
                            if !resolved.isExactFEM {
                                Text("Synthetic Preview (Click 'Run Stiffness Study' for exact K(α) eigenvectors)")
                                    .font(.system(size: 9.5, weight: .semibold))
                                    .foregroundColor(.orange)
                            }
                        }
                        .fixedSize()
                    } else {
                        Image(systemName: "info.circle")
                            .foregroundColor(.secondary)
                            .font(.caption2)
                        Text("Run Stiffness Study to inspect discrete N×N eigenvectors v_k(α).")
                            .font(.caption2)
                            .foregroundColor(.secondary)
                    }
                }
                .padding(.vertical, 1)
            }
            .frame(maxWidth: .infinity, alignment: .leading)
            .padding(.horizontal, 10)
            .padding(.vertical, 4)
            .background(BlueHour.controlSurface.opacity(0.85))
            
            Divider()
            
            // =========================================================================
            // 3. DUAL VISUALIZATION WORKSPACE (Left: N×N Lattice Mode | Right: Fourier Spectrum)
            // =========================================================================
            if let pt = currentPoint {
                let resolved = resolveLatticeMode(for: pt)
                HSplitView {
                    // LEFT PANE: Discrete N×N Lattice & Eigenvector Field x_i(α) + A * v_k(α)
                    VStack(spacing: 0) {
                        ScrollView(.horizontal, showsIndicators: false) {
                            HStack(spacing: 10) {
                                Label("Lattice \(resolved.gridNx)×\(resolved.gridNy) (\(resolved.nodesRefX.count) nodes) — \(selectedMode.label)",
                                      systemImage: "square.grid.4x3.fill")
                                    .font(.caption.bold())
                                    .foregroundColor(.secondary)
                                    .lineLimit(1)
                                    .fixedSize()
                                
                                Divider().frame(height: 12)
                                
                                HStack(spacing: 5) {
                                    Text(String(format: "Zoom %.1f×", zoomScale))
                                        .font(.caption2.monospaced())
                                        .foregroundColor(.secondary)
                                        .lineLimit(1)
                                        .fixedSize()
                                    Slider(value: $zoomScale, in: 1.0...10.0)
                                        .frame(width: 65)
                                        .controlSize(.mini)
                                    if zoomScale > 1.01 || abs(panOffset.width) > 1 || abs(panOffset.height) > 1 {
                                        Button("Reset") {
                                            zoomScale = 1.0
                                            panOffset = .zero
                                            dragStartOffset = .zero
                                        }
                                        .buttonStyle(.bordered)
                                        .controlSize(.mini)
                                        .font(.caption2)
                                        .fixedSize()
                                    }
                                }
                                .fixedSize()
                                
                                Toggle("Mesh", isOn: $showMeshEdges)
                                    .toggleStyle(.checkbox)
                                    .controlSize(.mini)
                                    .font(.caption2)
                                    .fixedSize()
                                
                                Toggle("Vectors", isOn: $showQuiver)
                                    .toggleStyle(.checkbox)
                                    .controlSize(.mini)
                                    .font(.caption2)
                                    .fixedSize()
                                
                                Button(action: { isAnimating.toggle() }) {
                                    Label(isAnimating ? "Pause" : "Animate",
                                          systemImage: isAnimating ? "pause.fill" : "play.fill")
                                        .font(.caption2)
                                }
                                .buttonStyle(.bordered)
                                .controlSize(.mini)
                                .fixedSize()
                                
                                HStack(spacing: 4) {
                                    Text("Amp:")
                                        .font(.caption2)
                                        .foregroundColor(.secondary)
                                        .lineLimit(1)
                                        .fixedSize()
                                    Slider(value: $modeAmplitude, in: 0.0...1.2)
                                        .frame(width: 65)
                                        .controlSize(.mini)
                                }
                                .fixedSize()
                                
                                HStack(spacing: 4) {
                                    Text("Nodes:")
                                        .font(.caption2)
                                        .foregroundColor(.secondary)
                                        .lineLimit(1)
                                        .fixedSize()
                                    Slider(value: $nodeScale, in: 0.0...3.0)
                                        .frame(width: 55)
                                        .controlSize(.mini)
                                    Text(nodeScale < 0.05 ? "Off" : String(format: "%.1f×", nodeScale))
                                        .font(.caption2.monospaced())
                                        .foregroundColor(.secondary)
                                        .lineLimit(1)
                                        .fixedSize()
                                }
                                .fixedSize()
                            }
                            .padding(.horizontal, 12)
                            .padding(.vertical, 6)
                        }
                        .background(BlueHour.controlSurface.opacity(0.35))
                        
                        ZStack(alignment: .top) {
                            GeometryReader { _ in
                                Canvas { ctx, size in
                                    drawDiscreteLatticeCanvas(ctx: &ctx, size: size, resolved: resolved)
                                }
                                .contentShape(Rectangle())
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
                            
                            if let msg = toastMessage {
                                HStack(spacing: 6) {
                                    Image(systemName: "checkmark.circle.fill").foregroundColor(.green)
                                    Text(msg).font(.caption.bold())
                                }
                                .padding(.horizontal, 14)
                                .padding(.vertical, 8)
                                .background(.ultraThinMaterial, in: RoundedRectangle(cornerRadius: 8))
                                .shadow(color: .black.opacity(0.3), radius: 6)
                                .padding(.top, 12)
                                .transition(.opacity.combined(with: .move(edge: .top)))
                            }
                        }
                        
                        // Footer Bar for Left Lattice Canvas
                        ScrollView(.horizontal, showsIndicators: false) {
                            HStack(spacing: 12) {
                                Text("Mode Character:")
                                    .font(.caption2.bold())
                                    .foregroundColor(.secondary)
                                    .fixedSize()
                                
                                GeometryReader { barGeo in
                                    let w = barGeo.size.width
                                    let sFrac = CGFloat(max(0.0, min(1.0, resolved.shearPercent / 100.0)))
                                    HStack(spacing: 0) {
                                        Rectangle()
                                            .fill(Color.cyan)
                                            .frame(width: w * sFrac)
                                        Rectangle()
                                            .fill(Color.pink.opacity(0.8))
                                            .frame(width: max(0, w * (1.0 - sFrac)))
                                    }
                                    .clipShape(Capsule())
                                }
                                .frame(width: 100, height: 7)
                                
                                Text(String(format: "Shear |m*×n*| = %.1f%%  •  Dilatation |m*·n*| = %.1f%%",
                                            resolved.shearPercent, 100.0 - resolved.shearPercent))
                                    .font(.caption2.monospaced())
                                    .foregroundColor(.secondary)
                                    .fixedSize()
                            }
                            .padding(.horizontal, 12)
                            .padding(.vertical, 6)
                        }
                        .background(BlueHour.controlSurface.opacity(0.5))
                    }
                    .frame(minWidth: 360, maxWidth: .infinity, maxHeight: .infinity)
                    
                    // RIGHT PANE: Directional Fourier Wave-Normal Spectrum |v̂_k(ξ)|² & Nodal Polarization
                    VStack(spacing: 0) {
                        ScrollView(.horizontal, showsIndicators: false) {
                            HStack(spacing: 10) {
                                Label("Directional Mode Spectrum |v̂(ξ)|² & (n*, m*) / (N*, M*)",
                                      systemImage: "scope")
                                    .font(.caption.bold())
                                    .foregroundColor(.secondary)
                                    .lineLimit(1)
                                    .fixedSize()
                                if let hDeg = hoveredAngleDeg {
                                    Text(String(format: "Cursor ξ = %.1f°", hDeg))
                                        .font(.caption2.monospaced())
                                        .foregroundColor(.accentColor)
                                        .fixedSize()
                                }
                            }
                            .padding(.horizontal, 12)
                            .padding(.vertical, 6)
                        }
                        .background(BlueHour.controlSurface.opacity(0.35))
                        
                        GeometryReader { geo in
                            Canvas { ctx, size in
                                drawDirectionalSpectrumRosette(ctx: &ctx, size: size, resolved: resolved)
                            }
                            .contentShape(Rectangle())
                            .gesture(
                                DragGesture(minimumDistance: 0)
                                    .onChanged { val in
                                        let center = CGPoint(x: geo.size.width / 2, y: geo.size.height / 2)
                                        let dx = Double(val.location.x - center.x)
                                        let dy = Double(center.y - val.location.y)
                                        if hypot(dx, dy) > 8 {
                                            var deg = atan2(dy, dx) * 180.0 / .pi
                                            if deg < 0 { deg += 360.0 }
                                            hoveredAngleDeg = deg
                                        }
                                    }
                                    .onEnded { _ in
                                        hoveredAngleDeg = nil
                                    }
                            )
                        }
                        
                        // Legend Footer for Right Polar Spectrum
                        ScrollView(.horizontal, showsIndicators: false) {
                            HStack(spacing: 12) {
                                if frameMode != .lagrangian {
                                    HStack(spacing: 4) {
                                        Circle().fill(Color.orange).frame(width: 7, height: 7)
                                        Text(String(format: "Eulerian n* = %.1f°", resolved.eulerNormalDeg))
                                            .font(.caption2.monospaced().bold())
                                            .foregroundColor(.orange)
                                    }
                                    HStack(spacing: 4) {
                                        Circle().fill(Color.yellow).frame(width: 7, height: 7)
                                        Text(String(format: "m* = %.1f°", resolved.eulerPolDeg))
                                            .font(.caption2.monospaced().bold())
                                            .foregroundColor(.yellow)
                                    }
                                }
                                if frameMode != .eulerian {
                                    HStack(spacing: 4) {
                                        Circle().fill(Color.purple).frame(width: 7, height: 7)
                                        Text(String(format: "Lagrangian N* = %.1f°", resolved.lagrNormalDeg))
                                            .font(.caption2.monospaced().bold())
                                            .foregroundColor(.purple)
                                    }
                                    HStack(spacing: 4) {
                                        Circle().fill(Color.pink).frame(width: 7, height: 7)
                                        Text(String(format: "M* = %.1f°", resolved.lagrPolDeg))
                                            .font(.caption2.monospaced().bold())
                                            .foregroundColor(.pink)
                                    }
                                }
                                HStack(spacing: 4) {
                                    Rectangle().fill(Color.cyan).frame(width: 10, height: 2)
                                    Text(String(format: "Slip Plane t* = %.1f°",
                                                frameMode == .lagrangian ? resolved.lagrSlipDeg : resolved.eulerSlipDeg))
                                        .font(.caption2.monospaced())
                                        .foregroundColor(.cyan)
                                }
                            }
                            .padding(.horizontal, 12)
                            .padding(.vertical, 6)
                        }
                        .background(BlueHour.controlSurface.opacity(0.5))
                    }
                    .frame(minWidth: 320, maxWidth: .infinity, maxHeight: .infinity)
                }
            } else {
                VStack(spacing: 10) {
                    Image(systemName: "scope")
                        .font(.system(size: 38))
                        .foregroundColor(.secondary.opacity(0.6))
                    Text("No Stiffness Eigenvector Data Loaded")
                        .font(.headline)
                        .foregroundColor(.secondary)
                    Text("Click 'Run Stiffness Study' on the left to compute N×N stiffness eigenvectors v₁(α), v₂(α), v₃(α).")
                        .font(.caption)
                        .foregroundColor(.secondary)
                }
                .frame(maxWidth: .infinity, maxHeight: .infinity)
            }
        }
        .background(BlueHour.canvasBackground)
        .onReceive(animTimer) { _ in
            if isAnimating {
                animPhase += 0.16
                if animPhase > 2.0 * .pi { animPhase -= 2.0 * .pi }
            }
        }
    }
    
    // MARK: - Resolved Discrete Lattice & Mode Data
    fileprivate struct ResolvedLatticeMode {
        var alpha: Double
        var lambdaMin: Double
        var lambda2: Double
        var lambda3: Double
        var lambda4: Double
        var lambda5: Double
        var isExactFEM: Bool
        var gridNx: Int
        var gridNy: Int
        var nodesRefX: [Double]
        var nodesRefY: [Double]
        var nodesEulerX: [Double]
        var nodesEulerY: [Double]
        var ux: [Double]
        var uy: [Double]
        var lagrUx: [Double]
        var lagrUy: [Double]
        var isFree: [Bool]
        var triangles: [(Int, Int, Int)]
        
        var eulerNormalDeg: Double
        var eulerSlipDeg: Double
        var eulerMx: Double
        var eulerMy: Double
        var eulerPolDeg: Double
        
        var lagrNormalDeg: Double
        var lagrSlipDeg: Double
        var lagrMx: Double
        var lagrMy: Double
        var lagrPolDeg: Double
        
        var ipr: Double
        var shearPercent: Double
        
        // 180-angle directional structure factor S(ξ) in [0, 360°)
        var anglesDeg: [Double]
        var eulerSpectrum: [Double]
        var lagrSpectrum: [Double]
    }
    
    private func resolveLatticeMode(for pt: StiffnessPointModel) -> ResolvedLatticeMode {
        var refX: [Double] = []
        var refY: [Double] = []
        var freeFlags: [Bool] = []
        var tris: [(Int, Int, Int)] = []
        var gridNx = max(4, manager.nx)
        var gridNy = max(4, manager.ny)
        var fullScaleSnap: StiffnessFullScaleSnapshot? = nil
        
        if let fs = manager.fullScaleData, fs.nx > 0, fs.ny > 0 {
            fullScaleSnap = fs.snapshots.min(by: { abs($0.alpha - pt.alpha) < abs($1.alpha - pt.alpha) })
            if let snap = fullScaleSnap, abs(snap.alpha - pt.alpha) > max(0.015, manager.dalpha * 2.5) {
                fullScaleSnap = nil
            }
        }
        
        var ux: [Double] = []
        var uy: [Double] = []
        var isExact = false
        
        // 1. USE THE EXACT CALCULATION MESH DIRECTLY (NO REMESHING!)
        let calcMesh = manager.fullScaleData?.mesh ?? manager.mesh
        if let m = calcMesh, !m.nodes_x.isEmpty, m.nodes_x.count == m.nodes_y.count {
            refX = m.nodes_x
            refY = m.nodes_y
            freeFlags = m.is_free.map { $0 != 0 }
            if freeFlags.count != refX.count {
                freeFlags = Array(repeating: true, count: refX.count)
            }
            tris = m.triangles.compactMap { t in
                guard t.count >= 3, t[0] < refX.count, t[1] < refX.count, t[2] < refX.count else { return nil }
                return (t[0], t[1], t[2])
            }
        } else {
            // Pristine generator fallback only if no calculation mesh file was loaded
            let isTri = (manager.latticeSymmetry == .triangular)
            let a0 = manager.fullScaleData?.a_0 ?? (isTri ? pow(4.0 / 3.0, 0.25) : 1.0)
            let isPBC = (manager.boundaryCondition == .periodic)
            refX.reserveCapacity(gridNx * gridNy)
            refY.reserveCapacity(gridNx * gridNy)
            freeFlags.reserveCapacity(gridNx * gridNy)
            for j in 0..<gridNy {
                for i in 0..<gridNx {
                    let x = isTri ? (Double(i) * a0 + (j % 2 == 1 ? 0.5 * a0 : 0.0)) : (Double(i) * a0)
                    let y = isTri ? (Double(j) * a0 * sqrt(3.0) / 2.0) : (Double(j) * a0)
                    refX.append(x)
                    refY.append(y)
                    let isPerim = (i == 0 || i == gridNx - 1 || j == 0 || j == gridNy - 1)
                    freeFlags.append(isPBC ? true : !isPerim)
                }
            }
            for j in 0..<(gridNy - 1) {
                for i in 0..<(gridNx - 1) {
                    let n00 = j * gridNx + i
                    let n10 = j * gridNx + (i + 1)
                    let n01 = (j + 1) * gridNx + i
                    let n11 = (j + 1) * gridNx + (i + 1)
                    tris.append((n00, n10, n11))
                    tris.append((n00, n11, n01))
                }
            }
        }
        
        let nNodes = refX.count
        let meanX = refX.reduce(0.0, +) / Double(max(1, nNodes))
        let meanY = refY.reduce(0.0, +) / Double(max(1, nNodes))
        let cRefX = refX.map { $0 - meanX }
        let cRefY = refY.map { $0 - meanY }
        
        let cEulerX = zip(cRefX, cRefY).map { $0.0 + pt.alpha * $0.1 }
        let cEulerY = cRefY
        
        // 2. Extract eigenvector (ux, uy) for selectedMode
        let chosenUx: [Double]?
        let chosenUy: [Double]?
        if let snap = fullScaleSnap {
            switch selectedMode {
            case .mode1:
                chosenUx = snap.mode1_ux; chosenUy = snap.mode1_uy
            case .mode2:
                chosenUx = snap.mode2_ux ?? snap.mode1_ux; chosenUy = snap.mode2_uy ?? snap.mode1_uy
            case .mode3:
                chosenUx = snap.mode3_ux ?? snap.mode1_ux; chosenUy = snap.mode3_uy ?? snap.mode1_uy
            case .mode4:
                chosenUx = snap.mode4_ux ?? snap.mode1_ux; chosenUy = snap.mode4_uy ?? snap.mode1_uy
            case .mode5:
                chosenUx = snap.mode5_ux ?? snap.mode1_ux; chosenUy = snap.mode5_uy ?? snap.mode1_uy
            }
        } else {
            switch selectedMode {
            case .mode1:
                chosenUx = pt.mode1_ux; chosenUy = pt.mode1_uy
            case .mode2:
                chosenUx = pt.mode2_ux ?? pt.mode1_ux; chosenUy = pt.mode2_uy ?? pt.mode1_uy
            case .mode3:
                chosenUx = pt.mode3_ux ?? pt.mode1_ux; chosenUy = pt.mode3_uy ?? pt.mode1_uy
            case .mode4:
                chosenUx = pt.mode4_ux ?? pt.mode1_ux; chosenUy = pt.mode4_uy ?? pt.mode1_uy
            case .mode5:
                chosenUx = pt.mode5_ux ?? pt.mode1_ux; chosenUy = pt.mode5_uy ?? pt.mode1_uy
            }
        }
        
        if let vx = chosenUx, let vy = chosenUy, vx.count == nNodes, vy.count == nNodes {
            ux = vx
            uy = vy
            isExact = true
        } else {
            // Analytical preview fallback only if no mode was computed at this step
            let spanX = max(1.0, (cRefX.max() ?? 1.0) - (cRefX.min() ?? -1.0))
            let spanY = max(1.0, (cRefY.max() ?? 1.0) - (cRefY.min() ?? -1.0))
            for i in 0..<nNodes {
                guard freeFlags[i] else {
                    ux.append(0.0); uy.append(0.0)
                    continue
                }
                let phase: Double
                switch selectedMode {
                case .mode1:
                    phase = 2.0 * .pi * cRefY[i] / spanY
                    ux.append(sin(phase))
                    uy.append(0.0)
                case .mode2:
                    phase = 2.0 * .pi * cRefX[i] / spanX
                    ux.append(0.0)
                    uy.append(sin(phase))
                case .mode3:
                    phase = 2.0 * .pi * (cRefX[i] / spanX + cRefY[i] / spanY)
                    ux.append(0.707 * sin(phase))
                    uy.append(-0.707 * sin(phase))
                case .mode4:
                    phase = 2.0 * .pi * (cRefX[i] / spanX - cRefY[i] / spanY)
                    ux.append(0.707 * sin(phase))
                    uy.append(0.707 * sin(phase))
                case .mode5:
                    phase = 4.0 * .pi * cRefY[i] / spanY
                    ux.append(sin(phase))
                    uy.append(0.0)
                }
            }
        }
        
        // Pull back nodal displacements to Lagrangian material frame: δU = F^{-1} δu = (δu_x - α δu_y, δu_y)
        var lagrUx: [Double] = Array(repeating: 0.0, count: nNodes)
        var lagrUy: [Double] = Array(repeating: 0.0, count: nNodes)
        var maxLagr = 0.0
        for i in 0..<nNodes {
            let lux = ux[i] - pt.alpha * uy[i]
            let luy = uy[i]
            lagrUx[i] = lux
            lagrUy[i] = luy
            maxLagr = max(maxLagr, hypot(lux, luy))
        }
        if maxLagr > 1e-12 {
            for i in 0..<nNodes {
                lagrUx[i] /= maxLagr
                lagrUy[i] /= maxLagr
            }
        }
        
        // 3. Compute Inverse Participation Ratio (IPR) and Dominant Polarization m*
        var sumSq = 0.0
        var sumQuad = 0.0
        // Structure tensor of nodal displacements M_ij = sum_k du_i(k) du_j(k) to get principal polarization m*
        var mxx = 0.0, myy = 0.0, mxy = 0.0
        for i in 0..<nNodes {
            let u2 = ux[i] * ux[i] + uy[i] * uy[i]
            sumSq += u2
            sumQuad += u2 * u2
            mxx += ux[i] * ux[i]
            myy += uy[i] * uy[i]
            mxy += ux[i] * uy[i]
        }
        let ipr = (sumSq > 1e-12) ? (Double(nNodes) * sumQuad / (sumSq * sumSq)) : 1.0
        
        // Principal eigenvector of 2x2 polarization covariance [mxx, mxy; mxy, myy]
        let trM = 0.5 * (mxx + myy)
        let discM = sqrt(max(0.0, 0.25 * (mxx - myy) * (mxx - myy) + mxy * mxy))
        let lMaxM = trM + discM
        var eulerMx = 1.0, eulerMy = 0.0
        if abs(mxy) > 1e-10 {
            let vx = mxy, vy = lMaxM - mxx
            let vn = max(1e-12, hypot(vx, vy))
            eulerMx = vx / vn; eulerMy = vy / vn
        } else if myy > mxx {
            eulerMx = 0.0; eulerMy = 1.0
        }
        
        // 4. Compute Directional Wave-Normal Spectrum S(ξ) via Element Gradient Tensor (∇δu)
        // Average element displacement gradient G = ∇_X δu across triangles weighted by |∇δu|^2
        var nxxL = 0.0, nyyL = 0.0, nxyL = 0.0
        var nxxE = 0.0, nyyE = 0.0, nxyE = 0.0
        for tri in tris {
            let i0 = tri.0, i1 = tri.1, i2 = tri.2
            guard i0 < nNodes, i1 < nNodes, i2 < nNodes else { continue }
            let dx1 = cRefX[i1] - cRefX[i0], dy1 = cRefY[i1] - cRefY[i0]
            let dx2 = cRefX[i2] - cRefX[i0], dy2 = cRefY[i2] - cRefY[i0]
            let detJ = dx1 * dy2 - dx2 * dy1
            guard abs(detJ) > 1e-9 else { continue }
            let dux1 = ux[i1] - ux[i0], duy1 = uy[i1] - uy[i0]
            let dux2 = ux[i2] - ux[i0], duy2 = uy[i2] - uy[i0]
            // Lagrangian gradient ∂(δu)/∂X
            let gu_xX = ( dy2 * dux1 - dy1 * dux2) / detJ
            let gu_xY = (-dx2 * dux1 + dx1 * dux2) / detJ
            let gu_yX = ( dy2 * duy1 - dy1 * duy2) / detJ
            let gu_yY = (-dx2 * duy1 + dx1 * duy2) / detJ
            
            // G^T G in Lagrangian frame: normal direction N maximizes ||G N||^2 = N^T (G^T G) N
            let gtg_xx = gu_xX * gu_xX + gu_yX * gu_yX
            let gtg_yy = gu_xY * gu_xY + gu_yY * gu_yY
            let gtg_xy = gu_xX * gu_xY + gu_yX * gu_yY
            nxxL += gtg_xx; nyyL += gtg_yy; nxyL += gtg_xy
            
            // Eulerian gradient ∂(δu)/∂x = (∂δu/∂X) F^{-1}, where F^{-1} = [1, -α; 0, 1]
            let gu_xx = gu_xX
            let gu_xy = gu_xY - pt.alpha * gu_xX
            let gu_yx = gu_yX
            let gu_yy = gu_yY - pt.alpha * gu_yX
            nxxE += gu_xx * gu_xx + gu_yx * gu_yx
            nyyE += gu_xy * gu_xy + gu_yy * gu_yy
            nxyE += gu_xx * gu_xy + gu_yx * gu_yy
        }
        
        // Evaluate directional spectrum S_Euler(ξ) = n(ξ)^T N_E n(ξ) and S_Lagr(ξ) = N(ξ)^T N_L N(ξ)
        let numAngles = 180
        var angles: [Double] = []
        var eSpec: [Double] = []
        var lSpec: [Double] = []
        var bestEDeg = 90.0, maxEVal = -1.0
        var bestLDeg = 90.0, maxLVal = -1.0
        
        for k in 0..<numAngles {
            let deg = Double(k) * (360.0 / Double(numAngles))
            let rad = deg * .pi / 180.0
            let c = cos(rad), s = sin(rad)
            angles.append(deg)
            let ev = c * c * nxxE + 2.0 * c * s * nxyE + s * s * nyyE
            let lv = c * c * nxxL + 2.0 * c * s * nxyL + s * s * nyyL
            eSpec.append(ev)
            lSpec.append(lv)
            if k < numAngles / 2 && ev > maxEVal {
                maxEVal = ev; bestEDeg = deg
            }
            if k < numAngles / 2 && lv > maxLVal {
                maxLVal = lv; bestLDeg = deg
            }
        }
        
        // Orient eulerMx, eulerMy with slip tangent t* = (-sin bestEDeg, cos bestEDeg)
        let eRad = bestEDeg * .pi / 180.0
        let tx = -sin(eRad), ty = cos(eRad)
        if eulerMx * tx + eulerMy * ty < 0 {
            eulerMx = -eulerMx; eulerMy = -eulerMy
        }
        var eulerPolDeg = atan2(eulerMy, eulerMx) * 180.0 / .pi
        if eulerPolDeg < 0 { eulerPolDeg += 360.0 }
        
        // Pull back m* to Lagrangian material polarization M* = F^{-1} m* / ||F^{-1} m*||
        let rawLMx = eulerMx - pt.alpha * eulerMy
        let rawLMy = eulerMy
        let normLM = max(1e-12, hypot(rawLMx, rawLMy))
        let lagrMx = rawLMx / normLM
        let lagrMy = rawLMy / normLM
        var lagrPolDeg = atan2(lagrMy, lagrMx) * 180.0 / .pi
        if lagrPolDeg < 0 { lagrPolDeg += 360.0 }
        
        // Shear vs Dilatation decomposition
        let nx = cos(eRad), ny = sin(eRad)
        let shearComp = abs(eulerMx * ny - eulerMy * nx)
        let dilatComp = abs(eulerMx * nx + eulerMy * ny)
        let totSq = max(1e-12, shearComp * shearComp + dilatComp * dilatComp)
        let shearPct = (shearComp * shearComp / totSq) * 100.0
        
        return ResolvedLatticeMode(
            alpha: pt.alpha,
            lambdaMin: pt.lambda_min,
            lambda2: pt.lambda_2,
            lambda3: pt.lambda_3,
            lambda4: pt.lambda_4,
            lambda5: pt.lambda_5,
            isExactFEM: isExact,
            gridNx: gridNx,
            gridNy: gridNy,
            nodesRefX: cRefX,
            nodesRefY: cRefY,
            nodesEulerX: cEulerX,
            nodesEulerY: cEulerY,
            ux: ux,
            uy: uy,
            lagrUx: lagrUx,
            lagrUy: lagrUy,
            isFree: freeFlags,
            triangles: tris,
            eulerNormalDeg: bestEDeg,
            eulerSlipDeg: fmod(bestEDeg + 90.0, 180.0),
            eulerMx: eulerMx,
            eulerMy: eulerMy,
            eulerPolDeg: eulerPolDeg,
            lagrNormalDeg: bestLDeg,
            lagrSlipDeg: fmod(bestLDeg + 90.0, 180.0),
            lagrMx: lagrMx,
            lagrMy: lagrMy,
            lagrPolDeg: lagrPolDeg,
            ipr: ipr,
            shearPercent: shearPct,
            anglesDeg: angles,
            eulerSpectrum: eSpec,
            lagrSpectrum: lSpec
        )
    }
    
    // MARK: - Publication Colormap Interpolation
    private func samplePalette(t: Double, palette: StiffnessColorPalette) -> Color {
        let u = max(0.0, min(1.0, t))
        switch palette {
        case .coolwarm:
            if u < 0.25 {
                let f = u / 0.25
                return Color(red: 0.05 + 0.05 * f, green: 0.25 + 0.60 * f, blue: 0.85 + 0.10 * f)
            } else if u < 0.5 {
                let f = (u - 0.25) / 0.25
                return Color(red: 0.10 + 0.10 * f, green: 0.85 + 0.00 * f, blue: 0.95 - 0.65 * f)
            } else if u < 0.75 {
                let f = (u - 0.5) / 0.25
                return Color(red: 0.20 + 0.75 * f, green: 0.85 + 0.05 * f, blue: 0.30 - 0.20 * f)
            } else {
                let f = (u - 0.75) / 0.25
                return Color(red: 0.95 + 0.00 * f, green: 0.90 - 0.75 * f, blue: 0.10 - 0.00 * f)
            }
        case .viridis:
            if u < 0.333 {
                let f = u / 0.333
                return Color(red: 0.267 - 0.076 * f, green: 0.004 + 0.403 * f, blue: 0.329 + 0.227 * f)
            } else if u < 0.666 {
                let f = (u - 0.333) / 0.333
                return Color(red: 0.191 - 0.063 * f, green: 0.407 + 0.226 * f, blue: 0.556 - 0.028 * f)
            } else {
                let f = (u - 0.666) / 0.334
                return Color(red: 0.128 + 0.865 * f, green: 0.633 + 0.273 * f, blue: 0.528 - 0.384 * f)
            }
        case .plasma:
            if u < 0.333 {
                let f = u / 0.333
                return Color(red: 0.050 + 0.565 * f, green: 0.029 + 0.061 * f, blue: 0.529 + 0.055 * f)
            } else if u < 0.666 {
                let f = (u - 0.333) / 0.333
                return Color(red: 0.615 + 0.275 * f, green: 0.090 + 0.302 * f, blue: 0.584 - 0.306 * f)
            } else {
                let f = (u - 0.666) / 0.334
                return Color(red: 0.890 + 0.051 * f, green: 0.392 + 0.584 * f, blue: 0.278 - 0.149 * f)
            }
        case .grayscale:
            let g = 0.12 + 0.85 * u
            return Color(red: g, green: g, blue: g)
        }
    }
    
    // MARK: - Left Canvas: Discrete N×N Lattice & Eigenvector Field
    fileprivate func drawDiscreteLatticeCanvas(
        ctx: inout GraphicsContext,
        size: CGSize,
        resolved: ResolvedLatticeMode,
        isPublicationWhite: Bool = false
    ) {
        let nNodes = resolved.nodesRefX.count
        guard nNodes > 0 else { return }
        
        let useLagrCoords = (frameMode == .lagrangian)
        let baseX = useLagrCoords ? resolved.nodesRefX : resolved.nodesEulerX
        let baseY = useLagrCoords ? resolved.nodesRefY : resolved.nodesEulerY
        let dispUx = useLagrCoords ? resolved.lagrUx : resolved.ux
        let dispUy = useLagrCoords ? resolved.lagrUy : resolved.uy
        
        let minX = (resolved.nodesEulerX + resolved.nodesRefX).min() ?? -5.0
        let maxX = (resolved.nodesEulerX + resolved.nodesRefX).max() ?? 5.0
        let minY = (resolved.nodesEulerY + resolved.nodesRefY).min() ?? -5.0
        let maxY = (resolved.nodesEulerY + resolved.nodesRefY).max() ?? -5.0
        let span = max(1e-3, max(maxX - minX, maxY - minY))
        
        let effZoom = CGFloat(zoomScale)
        let effPan  = panOffset
        let center  = CGPoint(x: size.width / 2 + effPan.width,
                              y: size.height / 2 + 8 + effPan.height)
        let scale   = (min(size.width, size.height) * 0.72 / CGFloat(span)) * effZoom
        
        func toScreen(_ x: Double, _ y: Double) -> CGPoint {
            CGPoint(x: center.x + CGFloat(x) * scale,
                    y: center.y - CGFloat(y) * scale)
        }
        
        let animFactor = isAnimating ? sin(animPhase) : 1.0
        let amp = modeAmplitude * animFactor
        
        var defX = Array(repeating: 0.0, count: nNodes)
        var defY = Array(repeating: 0.0, count: nNodes)
        for i in 0..<nNodes {
            defX[i] = baseX[i] + amp * dispUx[i]
            defY[i] = baseY[i] + amp * dispUy[i]
        }
        
        // 1. Draw Calculation FEM Triangles
        if showMeshEdges || nNodes <= 900 || effZoom >= 2.0 {
            for tri in resolved.triangles {
                let i0 = tri.0, i1 = tri.1, i2 = tri.2
                guard i0 < nNodes, i1 < nNodes, i2 < nNodes else { continue }
                let p0 = toScreen(defX[i0], defY[i0])
                let p1 = toScreen(defX[i1], defY[i1])
                let p2 = toScreen(defX[i2], defY[i2])
                
                // Viewport culling
                if (p0.x < -20 && p1.x < -20 && p2.x < -20) ||
                   (p0.x > size.width + 20 && p1.x > size.width + 20 && p2.x > size.width + 20) ||
                   (p0.y < -20 && p1.y < -20 && p2.y < -20) ||
                   (p0.y > size.height + 20 && p1.y > size.height + 20 && p2.y > size.height + 20) {
                    continue
                }
                
                var path = Path()
                path.move(to: p0)
                path.addLine(to: p1)
                path.addLine(to: p2)
                path.closeSubpath()
                
                // 1. Element Color Filling (Always active for rich continuous colormaps at any zoom level)
                let mag0 = hypot(dispUx[i0], dispUy[i0])
                let mag1 = hypot(dispUx[i1], dispUy[i1])
                let mag2 = hypot(dispUx[i2], dispUy[i2])
                let avgMag = (mag0 + mag1 + mag2) / 3.0
                
                let fillCol: Color
                switch colorField {
                case .magnitude:
                    fillCol = samplePalette(t: min(1.0, avgMag), palette: colorPalette).opacity(0.85)
                case .elementShear:
                    let du1 = hypot(dispUx[i1] - dispUx[i0], dispUy[i1] - dispUy[i0])
                    let du2 = hypot(dispUx[i2] - dispUx[i0], dispUy[i2] - dispUy[i0])
                    let shearInt = min(1.0, 0.8 * (du1 + du2))
                    fillCol = samplePalette(t: shearInt, palette: colorPalette).opacity(0.85)
                case .direction:
                    let uxAvg = (dispUx[i0] + dispUx[i1] + dispUx[i2]) / 3.0
                    let uyAvg = (dispUy[i0] + dispUy[i1] + dispUy[i2]) / 3.0
                    let angleNorm = fmod((atan2(uyAvg, uxAvg) / (2.0 * .pi)) + 1.0, 1.0)
                    fillCol = samplePalette(t: angleNorm, palette: colorPalette).opacity(0.85)
                }
                ctx.fill(path, with: .color(fillCol))
                
                if showMeshEdges {
                    let edgeCol = isPublicationWhite
                        ? Color.black.opacity(effZoom > 2.0 ? 0.35 : 0.08)
                        : Color.white.opacity(effZoom > 2.0 ? 0.35 : 0.08)
                    ctx.stroke(path, with: .color(edgeCol), lineWidth: effZoom > 3.0 ? 0.85 : 0.45)
                }
            }
        }
        
        // 2. If in .both mode, also draw the Lagrangian undeformed boundary box dashed in purple
        if frameMode == .both {
            let rMinX = resolved.nodesRefX.min() ?? 0, rMaxX = resolved.nodesRefX.max() ?? 0
            let rMinY = resolved.nodesRefY.min() ?? 0, rMaxY = resolved.nodesRefY.max() ?? 0
            var refBox = Path()
            refBox.move(to: toScreen(rMinX, rMinY))
            refBox.addLine(to: toScreen(rMaxX, rMinY))
            refBox.addLine(to: toScreen(rMaxX, rMaxY))
            refBox.addLine(to: toScreen(rMinX, rMaxY))
            refBox.closeSubpath()
            ctx.stroke(refBox, with: .color(Color.purple.opacity(isPublicationWhite ? 0.75 : 0.55)),
                       style: StrokeStyle(lineWidth: 1.4, dash: [5, 4]))
        }
        
        // 3. Draw Displaced Lattice Nodes & Decimated Quiver Arrows
        let cellPix = max(1.3, min(7.0, scale * 0.48))
        let baseNodeR: CGFloat = nNodes > 900 ? cellPix : (nNodes > 400 ? 2.4 : (nNodes > 150 ? 3.2 : 4.2))
        let nodeR: CGFloat = baseNodeR * CGFloat(nodeScale)
        let quiverStepX = max(1, resolved.gridNx / 16)
        let quiverStepY = max(1, resolved.gridNy / 16)
        let arrowScale = max(1.0, Double(quiverStepX) * 0.65)
        
        if nodeScale > 0.01 {
            for i in 0..<nNodes {
                let pDef = toScreen(defX[i], defY[i])
                // Viewport culling for zoomed full-scale mode
                if pDef.x < -12 || pDef.x > size.width + 12 || pDef.y < -12 || pDef.y > size.height + 12 {
                    continue
                }
                let mag = hypot(dispUx[i], dispUy[i])
                
                let dotColor: Color
                if !resolved.isFree[i] {
                    dotColor = .gray.opacity(0.6)
                } else {
                    switch colorField {
                    case .magnitude, .elementShear:
                        let t = min(1.0, max(0.0, mag))
                        dotColor = samplePalette(t: t, palette: colorPalette)
                    case .direction:
                        let angleNorm = fmod((atan2(dispUy[i], dispUx[i]) / (2.0 * .pi)) + 1.0, 1.0)
                        dotColor = samplePalette(t: angleNorm, palette: colorPalette)
                    }
                }
                let rect = CGRect(x: pDef.x - nodeR, y: pDef.y - nodeR, width: 2 * nodeR, height: 2 * nodeR)
                if nNodes > 2500 && effZoom < 2.5 && nodeScale <= 1.0 {
                    ctx.fill(Path(rect), with: .color(dotColor))
                } else {
                    ctx.fill(Path(ellipseIn: rect), with: .color(dotColor))
                }
            }
        }
        
        if showQuiver {
            let arrowColor = isPublicationWhite ? Color.black.opacity(0.85) : Color.white.opacity(0.85)
            for i in 0..<nNodes {
                if nNodes > 900 {
                    let ix = i % max(1, resolved.gridNx)
                    let iy = i / max(1, resolved.gridNx)
                    if (ix % quiverStepX != 0) || (iy % quiverStepY != 0) { continue }
                }
                let mag = hypot(dispUx[i], dispUy[i])
                guard mag > 0.04 else { continue }
                let pBase = toScreen(baseX[i], baseY[i])
                if pBase.x < -20 || pBase.x > size.width + 20 || pBase.y < -20 || pBase.y > size.height + 20 {
                    continue
                }
                let pTip = toScreen(baseX[i] + arrowScale * dispUx[i], baseY[i] + arrowScale * dispUy[i])
                var qLine = Path()
                qLine.move(to: pBase)
                qLine.addLine(to: pTip)
                ctx.stroke(qLine, with: .color(arrowColor), lineWidth: 1.15)
                drawArrowHead(ctx: &ctx, from: pBase, to: pTip, color: arrowColor, size: 4.5)
            }
        }
        
        // 4. Overlay Extracted Macroscopic Wave-Normal & Polarization Vectors at Center
        let vecR = span * 0.34
        if frameMode != .eulerian {
            let lRad = resolved.lagrNormalDeg * .pi / 180.0
            let NLx = cos(lRad), NLy = sin(lRad)
            let NTip = toScreen(vecR * NLx, vecR * NLy)
            var nPath = Path()
            nPath.move(to: center); nPath.addLine(to: NTip)
            ctx.stroke(nPath, with: .color(.purple),
                       style: StrokeStyle(lineWidth: 2.4, dash: frameMode == .both ? [4, 3] : []))
            drawArrowHead(ctx: &ctx, from: center, to: NTip, color: .purple, size: 8)
            
            let MTip = toScreen(vecR * resolved.lagrMx, vecR * resolved.lagrMy)
            var mPath = Path()
            mPath.move(to: center); mPath.addLine(to: MTip)
            ctx.stroke(mPath, with: .color(.pink),
                       style: StrokeStyle(lineWidth: 2.4, dash: frameMode == .both ? [4, 3] : []))
            drawArrowHead(ctx: &ctx, from: center, to: MTip, color: .pink, size: 8)
        }
        
        if frameMode != .lagrangian {
            let eRad = resolved.eulerNormalDeg * .pi / 180.0
            let nx = cos(eRad), ny = sin(eRad)
            let nTip = toScreen(vecR * nx, vecR * ny)
            var nPath = Path()
            nPath.move(to: center); nPath.addLine(to: nTip)
            ctx.stroke(nPath, with: .color(.orange), lineWidth: 2.5)
            drawArrowHead(ctx: &ctx, from: center, to: nTip, color: .orange, size: 8)
            
            let mTip = toScreen(vecR * resolved.eulerMx, vecR * resolved.eulerMy)
            var mPath = Path()
            mPath.move(to: center); mPath.addLine(to: mTip)
            let mColor = isPublicationWhite ? Color(red: 0.85, green: 0.72, blue: 0.0) : Color.yellow
            ctx.stroke(mPath, with: .color(mColor), lineWidth: 2.5)
            drawArrowHead(ctx: &ctx, from: center, to: mTip, color: mColor, size: 8)
        }
    }
    
    // MARK: - Right Canvas: Directional Wave-Normal Spectrum Rosette |v̂(ξ)|²
    fileprivate func drawDirectionalSpectrumRosette(
        ctx: inout GraphicsContext,
        size: CGSize,
        resolved: ResolvedLatticeMode,
        isPublicationWhite: Bool = false
    ) {
        let center = CGPoint(x: size.width / 2, y: size.height / 2)
        let maxR = max(40.0, min(size.width, size.height) * 0.39)
        let maxVal = max(1e-9, (resolved.eulerSpectrum + resolved.lagrSpectrum).max() ?? 1.0)
        
        let ringColor = isPublicationWhite ? Color.black.opacity(0.15) : Color.primary.opacity(0.10)
        let spokeColor = isPublicationWhite ? Color.black.opacity(0.12) : Color.primary.opacity(0.10)
        let labelColor = isPublicationWhite ? Color(white: 0.30) : Color.secondary
        
        // 1. Polar Grid Circles & Spokes
        for frac in [0.25, 0.5, 0.75, 1.0] {
            let r = maxR * CGFloat(frac)
            let rect = CGRect(x: center.x - r, y: center.y - r, width: 2 * r, height: 2 * r)
            ctx.stroke(Path(ellipseIn: rect), with: .color(ringColor), lineWidth: 0.8)
        }
        
        for deg in stride(from: 0, to: 360, by: 30) {
            let rad = Double(deg) * .pi / 180.0
            let pEnd = CGPoint(x: center.x + maxR * CGFloat(cos(rad)),
                               y: center.y - maxR * CGFloat(sin(rad)))
            var spoke = Path()
            spoke.move(to: center); spoke.addLine(to: pEnd)
            ctx.stroke(spoke, with: .color(spokeColor), lineWidth: 0.8)
            
            let lblR = maxR + 15
            let pLbl = CGPoint(x: center.x + lblR * CGFloat(cos(rad)),
                               y: center.y - lblR * CGFloat(sin(rad)))
            ctx.draw(Text("\(deg)°").font(.system(size: 9, design: .monospaced)).foregroundColor(labelColor),
                     at: pLbl, anchor: .center)
        }
        
        // 2. Draw Directional Gradient Power Lobes ||(∇δu) · n(ξ)||²
        func drawLobe(spec: [Double], color: Color, isDashed: Bool) {
            guard spec.count == resolved.anglesDeg.count, !spec.isEmpty else { return }
            var path = Path()
            for i in 0..<spec.count {
                let rad = resolved.anglesDeg[i] * .pi / 180.0
                let r = maxR * CGFloat(max(0.0, min(1.0, spec[i] / maxVal)))
                let pt = CGPoint(x: center.x + r * CGFloat(cos(rad)),
                                 y: center.y - r * CGFloat(sin(rad)))
                if i == 0 { path.move(to: pt) } else { path.addLine(to: pt) }
            }
            path.closeSubpath()
            if !isDashed {
                ctx.fill(path, with: .color(color.opacity(isPublicationWhite ? 0.18 : 0.12)))
            }
            ctx.stroke(path, with: .color(color),
                       style: StrokeStyle(lineWidth: isDashed ? 1.8 : 2.2, dash: isDashed ? [5, 4] : []))
        }
        
        if frameMode != .eulerian {
            drawLobe(spec: resolved.lagrSpectrum, color: .purple, isDashed: frameMode == .both)
        }
        if frameMode != .lagrangian {
            let eulerCol = isPublicationWhite ? Color(red: 0.05, green: 0.45, blue: 0.85) : Color.cyan
            drawLobe(spec: resolved.eulerSpectrum, color: eulerCol, isDashed: false)
        }
        
        // 3. Draw Slip Plane & Dominant (n*, m*) / (N*, M*) Rays
        if frameMode != .eulerian {
            let lRad = resolved.lagrNormalDeg * .pi / 180.0
            let NLx = cos(lRad), NLy = sin(lRad)
            let NTip = CGPoint(x: center.x + maxR * CGFloat(NLx), y: center.y - maxR * CGFloat(NLy))
            var nLine = Path(); nLine.move(to: center); nLine.addLine(to: NTip)
            ctx.stroke(nLine, with: .color(.purple),
                       style: StrokeStyle(lineWidth: 2.2, dash: frameMode == .both ? [4, 3] : []))
            drawArrowHead(ctx: &ctx, from: center, to: NTip, color: .purple, size: 8)
            
            let MTip = CGPoint(x: center.x + maxR * 0.85 * CGFloat(resolved.lagrMx),
                               y: center.y - maxR * 0.85 * CGFloat(resolved.lagrMy))
            var mLine = Path(); mLine.move(to: center); mLine.addLine(to: MTip)
            ctx.stroke(mLine, with: .color(.pink),
                       style: StrokeStyle(lineWidth: 2.2, dash: frameMode == .both ? [4, 3] : []))
            drawArrowHead(ctx: &ctx, from: center, to: MTip, color: .pink, size: 8)
        }
        
        if frameMode != .lagrangian {
            let eRad = resolved.eulerNormalDeg * .pi / 180.0
            let slipRad = eRad + .pi / 2.0
            var slipLine = Path()
            slipLine.move(to: CGPoint(x: center.x - maxR * CGFloat(cos(slipRad)),
                                      y: center.y + maxR * CGFloat(sin(slipRad))))
            slipLine.addLine(to: CGPoint(x: center.x + maxR * CGFloat(cos(slipRad)),
                                         y: center.y - maxR * CGFloat(sin(slipRad))))
            let slipCol = isPublicationWhite ? Color(red: 0.0, green: 0.60, blue: 0.70) : Color.cyan
            ctx.stroke(slipLine, with: .color(slipCol.opacity(0.85)),
                       style: StrokeStyle(lineWidth: 1.5, dash: [6, 4]))
            
            let nx = cos(eRad), ny = sin(eRad)
            let nTip = CGPoint(x: center.x + maxR * CGFloat(nx), y: center.y - maxR * CGFloat(ny))
            var nLine = Path(); nLine.move(to: center); nLine.addLine(to: nTip)
            ctx.stroke(nLine, with: .color(.orange), lineWidth: 2.4)
            drawArrowHead(ctx: &ctx, from: center, to: nTip, color: .orange, size: 8)
            
            let mTip = CGPoint(x: center.x + maxR * 0.85 * CGFloat(resolved.eulerMx),
                               y: center.y - maxR * 0.85 * CGFloat(resolved.eulerMy))
            var mLine = Path(); mLine.move(to: center); mLine.addLine(to: mTip)
            let mColor = isPublicationWhite ? Color(red: 0.85, green: 0.72, blue: 0.0) : Color.yellow
            ctx.stroke(mLine, with: .color(mColor), lineWidth: 2.4)
            drawArrowHead(ctx: &ctx, from: center, to: mTip, color: mColor, size: 8)
        }
    }
    
    private func drawArrowHead(ctx: inout GraphicsContext, from p0: CGPoint, to p1: CGPoint, color: Color, size: CGFloat) {
        let dx = p1.x - p0.x
        let dy = p1.y - p0.y
        let len = hypot(dx, dy)
        guard len > 1e-3 else { return }
        let ux = dx / len, uy = dy / len
        let px = -uy, py = ux
        let base = CGPoint(x: p1.x - ux * size, y: p1.y - uy * size)
        let left = CGPoint(x: base.x + px * (size * 0.55), y: base.y + py * (size * 0.55))
        let right = CGPoint(x: base.x - px * (size * 0.55), y: base.y - py * (size * 0.55))
        var head = Path()
        head.move(to: p1)
        head.addLine(to: left)
        head.addLine(to: right)
        head.closeSubpath()
        ctx.fill(head, with: .color(color))
    }
    
    // MARK: - Publication Export Actions
    fileprivate enum ExportTarget {
        case lattice
        case rosette
        case combined
    }
    
    fileprivate enum ExportBackground {
        case white
        case dark
    }
    
    @MainActor
    private func exportFigure(target: ExportTarget, bg: ExportBackground) {
        guard let pt = currentPoint else { return }
        let resolved = resolveLatticeMode(for: pt)
        
        let panel = NSSavePanel()
        panel.title = "Save Publication Figure (300 DPI PNG)"
        panel.allowedContentTypes = [UTType.png]
        
        let bgTag = (bg == .white) ? "white" : "dark"
        let tgtTag: String
        switch target {
        case .lattice: tgtTag = "lattice"
        case .rosette: tgtTag = "rosette"
        case .combined: tgtTag = "dual_figure"
        }
        panel.nameFieldStringValue = String(format: "%@_mode%d_alpha_%+.4f_%@.png",
                                            tgtTag, selectedMode.rawValue, pt.alpha, bgTag)
        
        if panel.runModal() == .OK, let url = panel.url {
            let width: CGFloat = (target == .combined) ? 1400 : 800
            let height: CGFloat = (target == .combined) ? 700 : 750
            
            let exportView = PublicationFigureView(
                target: target,
                bg: bg,
                resolved: resolved,
                selectedMode: selectedMode,
                frameMode: frameMode,
                colorField: colorField,
                colorPalette: colorPalette,
                showQuiver: showQuiver,
                showMeshEdges: showMeshEdges,
                modeAmplitude: modeAmplitude,
                zoomScale: zoomScale,
                panOffset: panOffset,
                nodeScale: nodeScale,
                inspector: self
            )
            .frame(width: width, height: height)
            
            let renderer = ImageRenderer(content: exportView)
            renderer.scale = 2.5
            if let nsImage = renderer.nsImage,
               let tiff = nsImage.tiffRepresentation,
               let rep = NSBitmapImageRep(data: tiff),
               let png = rep.representation(using: .png, properties: [:]) {
                try? png.write(to: url)
                showToast("Saved 300 DPI figure to \(url.lastPathComponent)")
            }
        }
    }
    
    @MainActor
    private func copyFigureToClipboard(target: ExportTarget) {
        guard let pt = currentPoint else { return }
        let resolved = resolveLatticeMode(for: pt)
        let width: CGFloat = (target == .combined) ? 1400 : 800
        let height: CGFloat = (target == .combined) ? 700 : 750
        
        let exportView = PublicationFigureView(
            target: target,
            bg: .white,
            resolved: resolved,
            selectedMode: selectedMode,
            frameMode: frameMode,
            colorField: colorField,
            colorPalette: colorPalette,
            showQuiver: showQuiver,
            showMeshEdges: showMeshEdges,
            modeAmplitude: modeAmplitude,
            zoomScale: zoomScale,
            panOffset: panOffset,
            nodeScale: nodeScale,
            inspector: self
        )
        .frame(width: width, height: height)
        
        let renderer = ImageRenderer(content: exportView)
        renderer.scale = 2.5
        if let nsImage = renderer.nsImage {
            NSPasteboard.general.clearContents()
            NSPasteboard.general.writeObjects([nsImage])
            showToast("Copied 300 DPI figure to clipboard! (Ready to paste)")
        }
    }
    
    private func showToast(_ msg: String) {
        withAnimation { toastMessage = msg }
        DispatchQueue.main.asyncAfter(deadline: .now() + 3.0) {
            withAnimation {
                if toastMessage == msg { toastMessage = nil }
            }
        }
    }
}

// MARK: - Dedicated Publication Figure View (Supports High-Res 300 DPI White/Dark)
fileprivate struct PublicationFigureView: View {
    let target: StiffnessModeInspectorView.ExportTarget
    let bg: StiffnessModeInspectorView.ExportBackground
    let resolved: StiffnessModeInspectorView.ResolvedLatticeMode
    let selectedMode: StiffnessEigenmodeChoice
    let frameMode: AcousticFrameMode
    let colorField: StiffnessColorField
    let colorPalette: StiffnessColorPalette
    let showQuiver: Bool
    let showMeshEdges: Bool
    let modeAmplitude: Double
    let zoomScale: Double
    let panOffset: CGSize
    let nodeScale: Double
    let inspector: StiffnessModeInspectorView
    
    var isWhite: Bool { bg == .white }
    var bgColor: Color { isWhite ? .white : Color(red: 0.11, green: 0.11, blue: 0.12) }
    var fgColor: Color { isWhite ? .black : .white }
    var subColor: Color { isWhite ? Color(white: 0.35) : Color(white: 0.70) }
    
    var body: some View {
        VStack(spacing: 8) {
            // Header Title
            HStack {
                let curLambda: Double = {
                    switch selectedMode {
                    case .mode1: return resolved.lambdaMin
                    case .mode2: return resolved.lambda2
                    case .mode3: return resolved.lambda3
                    case .mode4: return resolved.lambda4
                    case .mode5: return resolved.lambda5
                    }
                }()
                Text(String(format: "Discrete Stiffness Eigenmode %@   (α = %+.4f,  λ = %+.3e)",
                            selectedMode.label, resolved.alpha, curLambda))
                    .font(.system(size: 15, weight: .bold, design: .serif))
                    .foregroundColor(fgColor)
                Spacer()
                Text(String(format: "Lattice %d×%d (%d nodes)", resolved.gridNx, resolved.gridNy, resolved.nodesRefX.count))
                    .font(.system(size: 12, weight: .medium, design: .monospaced))
                    .foregroundColor(subColor)
            }
            .padding(.horizontal, 20)
            .padding(.top, 14)
            
            Divider()
            
            // Canvas View(s)
            HStack(spacing: 12) {
                if target == .lattice || target == .combined {
                    Canvas { ctx, size in
                        inspector.drawDiscreteLatticeCanvas(ctx: &ctx, size: size, resolved: resolved, isPublicationWhite: isWhite)
                    }
                }
                if target == .rosette || target == .combined {
                    Canvas { ctx, size in
                        inspector.drawDirectionalSpectrumRosette(ctx: &ctx, size: size, resolved: resolved, isPublicationWhite: isWhite)
                    }
                }
            }
            .padding(.horizontal, 14)
            
            Divider()
            
            // Footer Caption
            HStack(spacing: 16) {
                Text(String(format: "Eulerian: n* = %.1f°  m* = (%+.2f, %+.2f) [%.1f°]",
                            resolved.eulerNormalDeg, resolved.eulerMx, resolved.eulerMy, resolved.eulerPolDeg))
                    .font(.system(size: 11, weight: .semibold, design: .monospaced))
                    .foregroundColor(isWhite ? Color(red: 0.8, green: 0.4, blue: 0.0) : .orange)
                
                Text(String(format: "Lagrangian: N* = %.1f°  M* = (%+.2f, %+.2f) [%.1f°]",
                            resolved.lagrNormalDeg, resolved.lagrMx, resolved.lagrMy, resolved.lagrPolDeg))
                    .font(.system(size: 11, weight: .semibold, design: .monospaced))
                    .foregroundColor(isWhite ? Color(red: 0.5, green: 0.1, blue: 0.6) : .purple)
                
                Spacer()
                
                Text(String(format: "Slip t* = %.1f°  •  Shear = %.1f%%", resolved.eulerSlipDeg, resolved.shearPercent))
                    .font(.system(size: 11, weight: .medium, design: .monospaced))
                    .foregroundColor(subColor)
            }
            .padding(.horizontal, 20)
            .padding(.bottom, 12)
        }
        .background(bgColor)
    }
}
