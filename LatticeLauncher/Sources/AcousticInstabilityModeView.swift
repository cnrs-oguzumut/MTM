import SwiftUI
import AppKit

public enum AcousticFrameMode: Int, CaseIterable, Identifiable {
    case eulerian = 0
    case lagrangian = 1
    case both = 2
    
    public var id: Int { rawValue }
    public var label: String {
        switch self {
        case .eulerian: return "Eulerian n (Deformed)"
        case .lagrangian: return "Lagrangian N (Ref)"
        case .both: return "Both (n & N)"
        }
    }
}

public enum AcousticPolarMetric: Int, CaseIterable, Identifiable {
    case detQ = 0
    case lambdaMin = 1
    
    public var id: Int { rawValue }
    public var label: String {
        switch self {
        case .detQ: return "det Q(ξ)"
        case .lambdaMin: return "λ_min(Q(ξ))"
        }
    }
}

public struct AcousticInstabilityModeView: View {
    @ObservedObject var manager: AcousticStudyManager
    
    @State private var frameMode: AcousticFrameMode = .eulerian
    @State private var polarMetric: AcousticPolarMetric = .detQ
    @State private var modeAmplitude: Double = 0.22
    @State private var hoveredAngleDeg: Double? = nil
    
    private var pathPoints: [ShearPathPointModel] {
        manager.shearPath?.points ?? []
    }
    
    private var criticalAlphas: [Double] {
        manager.shearPath?.critical_alphas ?? []
    }
    
    private var periodAlpha: Double {
        manager.shearPath?.period_alpha ?? ((manager.latticeSymmetry == .triangular) ? sqrt(4.0 / 3.0) : 1.0)
    }
    
    private var alphaBounds: (min: Double, max: Double) {
        guard !pathPoints.isEmpty else { return (-1.5, 1.5) }
        let aMin = pathPoints.map(\.alpha).min() ?? -1.5
        let aMax = pathPoints.map(\.alpha).max() ?? 1.5
        return (aMin, aMax)
    }
    
    private var activeAlphaBinding: Binding<Double> {
        Binding<Double>(
            get: {
                if let h = manager.hoveredAlpha { return h }
                if let c = manager.selectedCriticalAlpha { return c }
                if let firstCrit = criticalAlphas.first { return firstCrit }
                return 0.0
            },
            set: { newVal in
                manager.hoveredAlpha = newVal
            }
        )
    }
    
    private var currentPoint: ShearPathPointModel? {
        guard !pathPoints.isEmpty else { return nil }
        let targetAlpha = activeAlphaBinding.wrappedValue
        if let crits = manager.shearPath?.critical_points,
           let cp = crits.first(where: { abs($0.alpha - targetAlpha) < 1e-4 }) {
            return cp
        }
        return pathPoints.min(by: { abs($0.alpha - targetAlpha) < abs($1.alpha - targetAlpha) })
    }
    
    public init(manager: AcousticStudyManager) {
        self.manager = manager
    }
    
    public var body: some View {
        VStack(spacing: 0) {
            // =========================================================================
            // 1. TOP RESPONSIVE TOOLBAR (Pinned grammar matching PathEnergyChartView)
            // =========================================================================
            ScrollView(.horizontal, showsIndicators: false) {
                HStack(spacing: 10) {
                    VStack(alignment: .leading, spacing: 1) {
                        HStack(spacing: 5) {
                            Image(systemName: "scope")
                                .foregroundColor(.accentColor)
                                .font(.system(size: 13, weight: .bold))
                            Text("Acoustic Instability Mode Q(n)")
                                .font(.system(size: 13, weight: .bold))
                        }
                        Text("Directional wave propagation normal & slip polarization eigenvector")
                            .font(.system(size: 10, design: .monospaced))
                            .foregroundColor(.secondary)
                            .lineLimit(1)
                    }
                    .fixedSize()
                    
                    Divider().frame(height: 16)
                    
                    // Jump to α menu + compact scrubber
                    if !pathPoints.isEmpty {
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
                                Text(String(format: "α = %+.4f", activeAlphaBinding.wrappedValue))
                                    .font(.system(size: 11, weight: .semibold, design: .monospaced))
                            }
                        }
                        .menuStyle(.borderedButton)
                        .controlSize(.small)
                        .fixedSize()
                        
                        HStack(spacing: 4) {
                            Text("Scrub α:")
                                .font(.caption2)
                                .foregroundColor(.secondary)
                            Slider(value: activeAlphaBinding, in: alphaBounds.min...alphaBounds.max)
                                .frame(width: 110)
                                .controlSize(.small)
                        }
                        .fixedSize()
                    }
                    
                    Divider().frame(height: 16)
                    
                    // Frame Picker: Eulerian n vs Lagrangian N vs Both
                    Picker("", selection: $frameMode) {
                        ForEach(AcousticFrameMode.allCases) { fm in
                            Text(fm.label).tag(fm)
                        }
                    }
                    .labelsHidden()
                    .pickerStyle(.segmented)
                    .controlSize(.small)
                    .fixedSize()
                    
                    Divider().frame(height: 16)
                    
                    // Polar Metric Picker: det Q(ξ) vs λ_min(Q(ξ))
                    Picker("", selection: $polarMetric) {
                        ForEach(AcousticPolarMetric.allCases) { pm in
                            Text(pm.label).tag(pm)
                        }
                    }
                    .labelsHidden()
                    .pickerStyle(.segmented)
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
                        let resolved = resolveModeAnalysis(for: pt)
                        let isStable = (pt.is_stable ?? (resolved.minDetQ > 0))
                        
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
                            Text("min det Q:").foregroundColor(.secondary)
                            Text(String(format: "%+.3f", resolved.minDetQ)).bold()
                                .foregroundColor(resolved.minDetQ > 0 ? .green : .red)
                        }
                        .font(.system(size: 10, design: .monospaced))
                        .fixedSize()
                        
                        Divider().frame(height: 12)
                        
                        HStack(spacing: 2) {
                            Text("λ_min(Q):").foregroundColor(.secondary)
                            Text(String(format: "%+.3f", resolved.lambdaMinQ)).bold()
                                .foregroundColor(resolved.lambdaMinQ > 0 ? .green : .red)
                        }
                        .font(.system(size: 10, design: .monospaced))
                        .fixedSize()
                        
                        Divider().frame(height: 12)
                        
                        HStack(spacing: 3) {
                            Text("Eulerian:")
                                .foregroundColor(.orange)
                                .bold()
                            Text(String(format: "n* = %.1f°", resolved.eulerAngleDeg))
                                .foregroundColor(.orange)
                            Text(String(format: "m* = (%+.2f, %+.2f) [%.1f°]", resolved.mx, resolved.my, resolved.eulerPolDeg))
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
                            Text(String(format: "N* = %.1f°", resolved.lagrAngleDeg))
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
                            Text(String(format: "%.0f%% Shear Slip", resolved.shearPercent))
                                .font(.system(size: 9, weight: .bold, design: .monospaced))
                                .padding(.horizontal, 5)
                                .padding(.vertical, 1)
                                .background(Color.cyan.opacity(0.18), in: Capsule())
                                .foregroundColor(.cyan)
                        }
                        .fixedSize()
                    } else {
                        Image(systemName: "info.circle")
                            .foregroundColor(.secondary)
                            .font(.caption2)
                        Text("Run Parametric Acoustic Study to inspect directional acoustic tensor modes Q(n).")
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
            // 3. DUAL VISUALIZATION WORKSPACE (Left: Polar Rosette | Right: Slip Band Cell)
            // =========================================================================
            if let pt = currentPoint {
                let resolved = resolveModeAnalysis(for: pt)
                HSplitView {
                    // LEFT PANE: Directional Polar Stability Rosette Q(ξ)
                    VStack(spacing: 0) {
                        ScrollView(.horizontal, showsIndicators: false) {
                            HStack(spacing: 10) {
                                Label(polarMetric == .detQ ? "Directional Stability Profile det Q(ξ)" : "Directional Soft-Mode Profile λ_min(Q(ξ))",
                                      systemImage: "circle.dashed")
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
                            ZStack {
                                Canvas { ctx, size in
                                    drawPolarRosette(ctx: &ctx, size: size, resolved: resolved)
                                }
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
                        
                        // Legend footer for Polar Plot
                        ScrollView(.horizontal, showsIndicators: false) {
                            HStack(spacing: 12) {
                                if frameMode == .eulerian || frameMode == .both {
                                    HStack(spacing: 4) {
                                        Circle().fill(Color.orange).frame(width: 7, height: 7)
                                        Text(String(format: "Eulerian n* = %.1f°", resolved.eulerAngleDeg))
                                            .font(.caption2.monospaced().bold())
                                            .foregroundColor(.orange)
                                    }
                                }
                                if frameMode == .lagrangian || frameMode == .both {
                                    HStack(spacing: 4) {
                                        Circle().fill(Color.purple).frame(width: 7, height: 7)
                                        Text(String(format: "Lagrangian N* = %.1f°", resolved.lagrAngleDeg))
                                            .font(.caption2.monospaced().bold())
                                            .foregroundColor(.purple)
                                    }
                                }
                                HStack(spacing: 4) {
                                    Rectangle().fill(Color.cyan).frame(width: 10, height: 2)
                                    Text(String(format: "Slip Plane t* = %.1f°",
                                                frameMode == .lagrangian ? resolved.slipPlaneLagrDeg : resolved.slipPlaneEulerDeg))
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
                    
                    // RIGHT PANE: 4-Atom Unit Cell (Undeformed Lagrangian vs Deformed Eulerian) + (n*, m*) / (N*, M*)
                    VStack(spacing: 0) {
                        ScrollView(.horizontal, showsIndicators: false) {
                            HStack(spacing: 10) {
                                Label(frameMode == .eulerian
                                      ? "4-Atom Unit Cell (Eulerian x = F·X) & (n*, m*)"
                                      : (frameMode == .lagrangian
                                         ? "4-Atom Unit Cell (Lagrangian X) & (N*, M*)"
                                         : "4-Atom Unit Cell: Eulerian (n*, m*) vs Lagrangian (N*, M*)"),
                                      systemImage: "square.grid.2x2")
                                    .font(.caption.bold())
                                    .foregroundColor(.secondary)
                                    .lineLimit(1)
                                    .fixedSize()
                                
                                Divider().frame(height: 12)
                                
                                HStack(spacing: 4) {
                                    Text("Mode Disp:")
                                        .font(.caption2)
                                        .foregroundColor(.secondary)
                                        .lineLimit(1)
                                        .fixedSize()
                                    Slider(value: $modeAmplitude, in: 0.0...0.45)
                                        .frame(width: 75)
                                        .controlSize(.mini)
                                }
                                .fixedSize()
                            }
                            .padding(.horizontal, 12)
                            .padding(.vertical, 6)
                        }
                        .background(BlueHour.controlSurface.opacity(0.35))
                        
                        GeometryReader { geo in
                            Canvas { ctx, size in
                                drawCrystalModeCell(ctx: &ctx, size: size, resolved: resolved)
                            }
                        }
                        
                        // Mode Decomposition Bar Footer
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
                                .frame(width: 110, height: 7)
                                
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
                    .frame(minWidth: 320, maxWidth: .infinity, maxHeight: .infinity)
                }
            } else {
                VStack(spacing: 10) {
                    Image(systemName: "scope")
                        .font(.system(size: 38))
                        .foregroundColor(.secondary.opacity(0.6))
                    Text("No Acoustic Mode Data Loaded")
                        .font(.headline)
                        .foregroundColor(.secondary)
                    Text("Click 'Run Parametric Acoustic Study' on the left to compute directional acoustic modes Q(n).")
                        .font(.caption)
                        .foregroundColor(.secondary)
                }
                .frame(maxWidth: .infinity, maxHeight: .infinity)
            }
        }
        .background(BlueHour.canvasBackground)
    }
    
    // MARK: - Mode Analysis Resolution (Exact C++ Quadratic Form + Fallback for Legacy JSON)
    private struct ResolvedModeData {
        var alpha: Double
        var minDetQ: Double
        var lambdaMinQ: Double
        var lambdaMaxQ: Double
        var eulerAngleDeg: Double
        var lagrAngleDeg: Double
        var slipPlaneEulerDeg: Double
        var slipPlaneLagrDeg: Double
        var mx: Double
        var my: Double
        var eulerPolDeg: Double
        var lagrMx: Double
        var lagrMy: Double
        var lagrPolDeg: Double
        var shearPercent: Double
        var f11: Double, f12: Double, f21: Double, f22: Double
        // 180-point polar profiles on [0, 360°)
        var eulerAngles: [Double]
        var eulerDetQ: [Double]
        var eulerLambdaMin: [Double]
        var lagrDetQ: [Double]
        var lagrLambdaMin: [Double]
    }
    
    private func resolveModeAnalysis(for pt: ShearPathPointModel) -> ResolvedModeData {
        // Reconstruct macroscopic deformation gradient F_path(alpha) from undeformed B_0 to deformed B_t
        let f11: Double
        let f12: Double
        let f21: Double
        let f22: Double
        if let a = pt.f11, let b = pt.f12, let c = pt.f21, let d = pt.f22 {
            f11 = a; f12 = b; f21 = c; f22 = d
        } else {
            var fp11 = 1.0, fp12 = 0.0, fp21 = 0.0, fp22 = 1.0
            if manager.pathType == "pure_shear" {
                let lam = exp(0.5 * pt.alpha)
                fp11 = 1.0 / lam; fp22 = lam
            } else if manager.pathType == "diagonal_shear" {
                let ch = cosh(pt.alpha), sh = sinh(pt.alpha), sq = sqrt(ch)
                fp11 = ch / sq; fp12 = sh / sq; fp22 = 1.0 / sq
            } else {
                let th = manager.shearTheta * .pi / 180.0
                let c = cos(th), s = sin(th)
                fp11 = 1.0 + pt.alpha * c * s
                fp12 = pt.alpha * c * c
                fp21 = -pt.alpha * s * s
                fp22 = 1.0 - pt.alpha * c * s
            }
            f11 = fp11; f12 = fp12; f21 = fp21; f22 = fp22
        }
        
        // Basis matrices Q11, Q22, Q12 in Lagrangian N coordinates (undeformed crystal B_0)
        let q11_xx: Double, q11_yy: Double, q11_xy: Double
        let q22_xx: Double, q22_yy: Double, q22_xy: Double
        let q12_xx: Double, q12_yy: Double, q12_xy: Double
        
        if let k11x = pt.q11_xx, let k11y = pt.q11_yy, let k11xy = pt.q11_xy,
           let k22x = pt.q22_xx, let k22y = pt.q22_yy, let k22xy = pt.q22_xy,
           let k12x = pt.q12_xx, let k12y = pt.q12_yy, let k12xy = pt.q12_xy,
           (abs(k11x) + abs(k22y) > 1e-12) {
            q11_xx = k11x; q11_yy = k11y; q11_xy = k11xy
            q22_xx = k22x; q22_yy = k22y; q22_xy = k22xy
            q12_xx = k12x; q12_yy = k12y; q12_xy = k12xy
        } else {
            // Fallback for legacy JSON files: construct anisotropic acoustic tensor pushed by F_path
            let ec = manager.elasticConstants
            let c11 = max(1.0, ec.C11)
            let c66 = max(0.2, ec.C66)
            let c12 = ec.C12
            let nearestGrid = manager.points.min(by: {
                hypot($0.x - pt.x, $0.y - pt.y) < hypot($1.x - pt.x, $1.y - pt.y)
            })
            let critDeg = pt.angle_deg ?? nearestGrid?.angle_deg ?? 90.0
            let critRad = critDeg * .pi / 180.0
            let targetDet = pt.min_detAc ?? 0.0
            let baseQ11_xx = c11, baseQ11_yy = c66
            let baseQ22_xx = c66, baseQ22_yy = c11
            let baseQ12_xy = c12 + c66
            let c = cos(critRad), s = sin(critRad)
            let qxx0 = c * c * baseQ11_xx + s * s * baseQ22_xx
            let qyy0 = c * c * baseQ11_yy + s * s * baseQ22_yy
            let qxy0 = c * s * baseQ12_xy
            let det0 = max(1e-6, qxx0 * qyy0 - qxy0 * qxy0)
            let scale = targetDet / det0
            q11_xx = baseQ11_xx; q11_yy = baseQ11_yy * scale; q11_xy = 0.0
            q22_xx = baseQ22_xx * scale; q22_yy = baseQ22_yy; q22_xy = 0.0
            q12_xx = 0.0; q12_yy = 0.0; q12_xy = baseQ12_xy
        }
        
        // Evaluate Q^Lagr(N) and Q^Euler(n) around the full circle [0, 360°)
        let numSamples = 180
        var angles: [Double] = []
        var eDet: [Double] = []
        var eLMin: [Double] = []
        var lDet: [Double] = []
        var lLMin: [Double] = []
        
        var bestEulerIdx = 0
        var bestEulerDet = Double.greatestFiniteMagnitude
        var bestLagrIdx = 0
        var bestLagrDet = Double.greatestFiniteMagnitude
        
        for i in 0..<numSamples {
            let deg = Double(i) * (360.0 / Double(numSamples))
            let rad = deg * .pi / 180.0
            let cn = cos(rad), sn = sin(rad)
            angles.append(deg)
            
            // 1. Lagrangian Q(N) at N = (cos ξ, sin ξ) in undeformed crystal B_0
            let l_xx = cn * cn * q11_xx + cn * sn * q12_xx + sn * sn * q22_xx
            let l_yy = cn * cn * q11_yy + cn * sn * q12_yy + sn * sn * q22_yy
            let l_xy = cn * cn * q11_xy + cn * sn * q12_xy + sn * sn * q22_xy
            let ld = l_xx * l_yy - l_xy * l_xy
            let lTr = 0.5 * (l_xx + l_yy)
            let lDisc = sqrt(max(0.0, 0.25 * (l_xx - l_yy) * (l_xx - l_yy) + l_xy * l_xy))
            let lMin = lTr - lDisc
            lDet.append(ld)
            lLMin.append(lMin)
            if i < numSamples / 2 && ld < bestLagrDet {
                bestLagrDet = ld
                bestLagrIdx = i
            }
            
            // 2. Eulerian Q(n): map spatial n = (cos ξ, sin ξ) back to N = F^T n / ||F^T n||
            let nt_x = f11 * cn + f21 * sn
            let nt_y = f12 * cn + f22 * sn
            let nt_norm = max(1e-12, hypot(nt_x, nt_y))
            let nx = nt_x / nt_norm
            let ny = nt_y / nt_norm
            let e_xx = nx * nx * q11_xx + nx * ny * q12_xx + ny * ny * q22_xx
            let e_yy = nx * ny * q12_yy + nx * nx * q11_yy + ny * ny * q22_yy
            let e_xy = nx * nx * q11_xy + nx * ny * q12_xy + ny * ny * q22_xy
            let ed = e_xx * e_yy - e_xy * e_xy
            let eTr = 0.5 * (e_xx + e_yy)
            let eDisc = sqrt(max(0.0, 0.25 * (e_xx - e_yy) * (e_xx - e_yy) + e_xy * e_xy))
            let eMin = eTr - eDisc
            eDet.append(ed)
            eLMin.append(eMin)
            if i < numSamples / 2 && ed < bestEulerDet {
                bestEulerDet = ed
                bestEulerIdx = i
            }
        }
        
        let lagrDeg = pt.angle_deg ?? angles[bestLagrIdx]
        // Compute exact Eulerian normal n* = F^{-T} N* / ||F^{-T} N*||
        let eulerDeg: Double
        if let eDeg = pt.angle_euler_deg {
            eulerDeg = eDeg
        } else {
            let lRad = lagrDeg * .pi / 180.0
            let NLx = cos(lRad), NLy = sin(lRad)
            let detF = f11 * f22 - f12 * f21
            if abs(detF) > 1e-12 {
                // F^{-T} = (1/detF) * [ f22, -f21 ; -f12, f11 ]
                let nex = ( f22 * NLx - f21 * NLy) / detF
                let ney = (-f12 * NLx + f11 * NLy) / detF
                var deg = atan2(ney, nex) * 180.0 / .pi
                deg = fmod(deg, 180.0)
                if deg < 0 { deg += 180.0 }
                eulerDeg = deg
            } else {
                eulerDeg = angles[bestEulerIdx]
            }
        }
        
        // Eigenvector m* at critical normal
        let eRad = eulerDeg * .pi / 180.0
        let cn = cos(eRad), sn = sin(eRad)
        let nt_x = f11 * cn + f21 * sn
        let nt_y = f12 * cn + f22 * sn
        let nt_norm = max(1e-12, hypot(nt_x, nt_y))
        let nx = nt_x / nt_norm, ny = nt_y / nt_norm
        let qxx = nx * nx * q11_xx + nx * ny * q12_xx + ny * ny * q22_xx
        let qyy = nx * nx * q11_yy + nx * ny * q12_yy + ny * ny * q22_yy
        let qxy = nx * nx * q11_xy + nx * ny * q12_xy + ny * ny * q22_xy
        let tr = 0.5 * (qxx + qyy)
        let disc = sqrt(max(0.0, 0.25 * (qxx - qyy) * (qxx - qyy) + qxy * qxy))
        let lMinVal = pt.lambda_min_Q ?? (tr - disc)
        let lMaxVal = pt.lambda_max_Q ?? (tr + disc)
        
        var mx: Double
        var my: Double
        if let pmx = pt.mode_mx, let pmy = pt.mode_my, hypot(pmx, pmy) > 1e-6 {
            mx = pmx; my = pmy
        } else {
            if abs(qxy) > 1e-10 {
                let vx = qxy
                let vy = (tr - disc) - qxx
                let vn = max(1e-12, hypot(vx, vy))
                mx = vx / vn; my = vy / vn
            } else if qxx <= qyy {
                mx = 1.0; my = 0.0
            } else {
                mx = 0.0; my = 1.0
            }
            let tx = -sn, ty = cn
            if mx * tx + my * ty < 0 {
                mx = -mx; my = -my
            }
        }
        
        // Pull Eulerian polarization m* back to Lagrangian material frame via M* = F^{-1} m* / ||F^{-1} m*||
        let detF = f11 * f22 - f12 * f21
        var lagrMx = mx
        var lagrMy = my
        if abs(detF) > 1e-12 {
            // F^{-1} = (1/detF) * [ f22, -f12 ; -f21, f11 ]
            let rawMx = ( f22 * mx - f12 * my) / detF
            let rawMy = (-f21 * mx + f11 * my) / detF
            let normM = hypot(rawMx, rawMy)
            if normM > 1e-12 {
                lagrMx = rawMx / normM
                lagrMy = rawMy / normM
            }
        }
        var eulerPolDeg = atan2(my, mx) * 180.0 / .pi
        if eulerPolDeg < 0 { eulerPolDeg += 360.0 }
        var lagrPolDeg = atan2(lagrMy, lagrMx) * 180.0 / .pi
        if lagrPolDeg < 0 { lagrPolDeg += 360.0 }
        
        // Shear vs Dilatation decomposition in Eulerian frame:
        // n* = (cos eRad, sin eRad), shear component = |m* x n*|, dilatation = |m* . n*|
        let shearComp = abs(mx * sn - my * cn)
        let dilatComp = abs(mx * cn + my * sn)
        let totSq = max(1e-12, shearComp * shearComp + dilatComp * dilatComp)
        let shearPct = (shearComp * shearComp / totSq) * 100.0
        
        let slipEuler = fmod(eulerDeg + 90.0, 180.0)
        let slipLagr = fmod(lagrDeg + 90.0, 180.0)
        
        return ResolvedModeData(
            alpha: pt.alpha,
            minDetQ: pt.min_detAc ?? bestEulerDet,
            lambdaMinQ: lMinVal,
            lambdaMaxQ: lMaxVal,
            eulerAngleDeg: eulerDeg,
            lagrAngleDeg: lagrDeg,
            slipPlaneEulerDeg: slipEuler,
            slipPlaneLagrDeg: slipLagr,
            mx: mx,
            my: my,
            eulerPolDeg: eulerPolDeg,
            lagrMx: lagrMx,
            lagrMy: lagrMy,
            lagrPolDeg: lagrPolDeg,
            shearPercent: shearPct,
            f11: f11, f12: f12, f21: f21, f22: f22,
            eulerAngles: angles,
            eulerDetQ: eDet,
            eulerLambdaMin: eLMin,
            lagrDetQ: lDet,
            lagrLambdaMin: lLMin
        )
    }
    
    // MARK: - Left Canvas: Polar Stability Rosette
    private func drawPolarRosette(ctx: inout GraphicsContext, size: CGSize, resolved: ResolvedModeData) {
        let center = CGPoint(x: size.width / 2, y: size.height / 2)
        let maxR = max(40.0, min(size.width, size.height) * 0.39)
        let zeroR = maxR * 0.36 // Radius representing Q = 0 threshold circle
        
        let eulerVals = (polarMetric == .detQ) ? resolved.eulerDetQ : resolved.eulerLambdaMin
        let lagrVals  = (polarMetric == .detQ) ? resolved.lagrDetQ  : resolved.lagrLambdaMin
        let allVals   = eulerVals + lagrVals
        let maxPos    = max(1e-6, allVals.max() ?? 1.0)
        let maxNeg    = max(1e-6, -(allVals.min() ?? 0.0))
        
        func radiusForVal(_ v: Double) -> CGFloat {
            if v >= 0 {
                return zeroR + (maxR - zeroR) * CGFloat(min(1.0, v / maxPos))
            } else {
                return zeroR * (1.0 - 0.85 * CGFloat(min(1.0, (-v) / maxNeg)))
            }
        }
        
        // 1. Polar Grid Circles
        for frac in [0.25, 0.5, 0.75, 1.0] {
            let r = zeroR + (maxR - zeroR) * CGFloat(frac)
            let rect = CGRect(x: center.x - r, y: center.y - r, width: 2 * r, height: 2 * r)
            ctx.stroke(Path(ellipseIn: rect), with: .color(Color.primary.opacity(0.10)), lineWidth: 0.8)
        }
        
        // Unstable Core Circle (inside zeroR: det Q <= 0)
        let zeroRect = CGRect(x: center.x - zeroR, y: center.y - zeroR, width: 2 * zeroR, height: 2 * zeroR)
        ctx.fill(Path(ellipseIn: zeroRect), with: .color(Color.red.opacity(0.08)))
        ctx.stroke(Path(ellipseIn: zeroRect), with: .color(Color.yellow.opacity(0.75)),
                   style: StrokeStyle(lineWidth: 1.3, dash: [4, 3]))
        
        // 2. Angular Spokes & Degree Labels
        for deg in stride(from: 0, to: 360, by: 30) {
            let rad = Double(deg) * .pi / 180.0
            let pEnd = CGPoint(x: center.x + maxR * CGFloat(cos(rad)),
                               y: center.y - maxR * CGFloat(sin(rad)))
            var spoke = Path()
            spoke.move(to: center)
            spoke.addLine(to: pEnd)
            ctx.stroke(spoke, with: .color(Color.primary.opacity(deg % 90 == 0 ? 0.22 : 0.10)), lineWidth: 0.8)
            
            let lblR = maxR + 15
            let pLbl = CGPoint(x: center.x + lblR * CGFloat(cos(rad)),
                               y: center.y - lblR * CGFloat(sin(rad)))
            let txt = Text("\(deg)°")
                .font(.system(size: 9, design: .monospaced))
                .foregroundColor(.secondary)
            ctx.draw(txt, at: pLbl, anchor: .center)
        }
        
        // 3. Helper to draw a polar lobe curve
        func drawPolarProfile(vals: [Double], isDashed: Bool, baseColor: Color) {
            guard vals.count == resolved.eulerAngles.count, !vals.isEmpty else { return }
            var closedPath = Path()
            for i in 0..<vals.count {
                let rad = resolved.eulerAngles[i] * .pi / 180.0
                let r = radiusForVal(vals[i])
                let pt = CGPoint(x: center.x + r * CGFloat(cos(rad)),
                                 y: center.y - r * CGFloat(sin(rad)))
                if i == 0 { closedPath.move(to: pt) } else { closedPath.addLine(to: pt) }
            }
            closedPath.closeSubpath()
            if !isDashed {
                ctx.fill(closedPath, with: .color(baseColor.opacity(0.10)))
            }
            // Draw segment by segment colored by stability
            for i in 0..<vals.count {
                let j = (i + 1) % vals.count
                let rad0 = resolved.eulerAngles[i] * .pi / 180.0
                let rad1 = resolved.eulerAngles[j] * .pi / 180.0
                let r0 = radiusForVal(vals[i])
                let r1 = radiusForVal(vals[j])
                let p0 = CGPoint(x: center.x + r0 * CGFloat(cos(rad0)), y: center.y - r0 * CGFloat(sin(rad0)))
                let p1 = CGPoint(x: center.x + r1 * CGFloat(cos(rad1)), y: center.y - r1 * CGFloat(sin(rad1)))
                var seg = Path()
                seg.move(to: p0)
                seg.addLine(to: p1)
                let segStable = (vals[i] > 0 && vals[j] > 0)
                let col = isDashed ? Color.purple : (segStable ? Color.green : Color.red)
                ctx.stroke(seg, with: .color(col),
                           style: StrokeStyle(lineWidth: isDashed ? 1.8 : 2.3, dash: isDashed ? [5, 4] : []))
            }
        }
        
        if frameMode == .lagrangian || frameMode == .both {
            drawPolarProfile(vals: lagrVals, isDashed: frameMode == .both, baseColor: .purple)
        }
        if frameMode == .eulerian || frameMode == .both {
            drawPolarProfile(vals: eulerVals, isDashed: false, baseColor: .green)
        }
        
        // 4. Draw Critical Wave-Normal Axis ±n* (or ±N*) and Orthogonal Slip Plane ±t*
        let activeCritDeg = (frameMode == .lagrangian) ? resolved.lagrAngleDeg : resolved.eulerAngleDeg
        let critRad = activeCritDeg * .pi / 180.0
        let slipRad = critRad + .pi / 2.0
        
        // Slip Plane Axis ±t* (Cyan dashed)
        var slipLine = Path()
        slipLine.move(to: CGPoint(x: center.x - maxR * CGFloat(cos(slipRad)),
                                  y: center.y + maxR * CGFloat(sin(slipRad))))
        slipLine.addLine(to: CGPoint(x: center.x + maxR * CGFloat(cos(slipRad)),
                                     y: center.y - maxR * CGFloat(sin(slipRad))))
        ctx.stroke(slipLine, with: .color(Color.cyan.opacity(0.85)),
                   style: StrokeStyle(lineWidth: 1.5, dash: [6, 4]))
        
        // Critical Wave-Normal Axis ±n* (Orange/Red solid arrow)
        let nColor = (frameMode == .lagrangian) ? Color.purple : Color.orange
        let nEnd1 = CGPoint(x: center.x + maxR * CGFloat(cos(critRad)),
                            y: center.y - maxR * CGFloat(sin(critRad)))
        let nEnd2 = CGPoint(x: center.x - maxR * CGFloat(cos(critRad)),
                            y: center.y + maxR * CGFloat(sin(critRad)))
        var nLine = Path()
        nLine.move(to: nEnd2)
        nLine.addLine(to: nEnd1)
        ctx.stroke(nLine, with: .color(nColor), lineWidth: 2.0)
        drawArrowHead(ctx: &ctx, from: center, to: nEnd1, color: nColor, size: 8)
        
        // Also draw Lagrangian N* axis if in .both mode
        if frameMode == .both {
            let lRad = resolved.lagrAngleDeg * .pi / 180.0
            let lEnd1 = CGPoint(x: center.x + maxR * CGFloat(cos(lRad)),
                                y: center.y - maxR * CGFloat(sin(lRad)))
            let lEnd2 = CGPoint(x: center.x - maxR * CGFloat(cos(lRad)),
                                y: center.y + maxR * CGFloat(sin(lRad)))
            var lLine = Path()
            lLine.move(to: lEnd2)
            lLine.addLine(to: lEnd1)
            ctx.stroke(lLine, with: .color(Color.purple.opacity(0.85)),
                       style: StrokeStyle(lineWidth: 1.5, dash: [3, 3]))
        }
    }
    
    // MARK: - Right Canvas: 4-Atom Fundamental Unit Cell (Undeformed Lagrangian vs Deformed Eulerian) + (n*, m*) / (N*, M*)
    private func drawCrystalModeCell(ctx: inout GraphicsContext, size: CGSize, resolved: ResolvedModeData) {
        let center = CGPoint(x: size.width / 2, y: size.height / 2 + 8)
        let scale: CGFloat = min(size.width, size.height) * 0.25
        
        func toScreen(_ x: Double, _ y: Double) -> CGPoint {
            CGPoint(x: center.x + CGFloat(x) * scale,
                    y: center.y - CGFloat(y) * scale)
        }
        
        // 1. Subtle Cartesian Axes Crosshair at Origin
        let axisLen = 1.45
        var xAxis = Path()
        xAxis.move(to: toScreen(-axisLen, 0))
        xAxis.addLine(to: toScreen(axisLen, 0))
        ctx.stroke(xAxis, with: .color(Color.primary.opacity(0.08)), lineWidth: 0.8)
        var yAxis = Path()
        yAxis.move(to: toScreen(0, -axisLen))
        yAxis.addLine(to: toScreen(0, axisLen))
        ctx.stroke(yAxis, with: .color(Color.primary.opacity(0.08)), lineWidth: 0.8)
        
        // 2. Define the 4 Fundamental Unit-Cell Atoms in Undeformed Lagrangian Frame B_0 (centered at origin)
        // Atom 1: (0,0), Atom 2: E1, Atom 3: E1+E2, Atom 4: E2 (shifted by -(E1+E2)/2)
        let isTri = (manager.latticeSymmetry == .triangular)
        let rawAtoms: [(Double, Double)]
        if isTri {
            let a = pow(4.0 / 3.0, 0.25)
            let e1x = a, e1y = 0.0
            let e2x = 0.5 * a, e2y = 0.5 * sqrt(3.0) * a
            let cx = 0.5 * (e1x + e2x), cy = 0.5 * (e1y + e2y)
            rawAtoms = [
                (0.0 - cx, 0.0 - cy),
                (e1x - cx, e1y - cy),
                (e1x + e2x - cx, e1y + e2y - cy),
                (e2x - cx, e2y - cy)
            ]
        } else {
            rawAtoms = [
                (-0.65, -0.65),
                ( 0.65, -0.65),
                ( 0.65,  0.65),
                (-0.65,  0.65)
            ]
        }
        
        // Apply Rank-1 Mode Perturbation to the 4 atoms if slider > 0:
        // In Lagrangian: X_a' = X_a + A * (N* · X_a) * M*
        // In Eulerian:   x_a' = F*X_a + A * (n* · (F*X_a)) * m*
        let nRad = resolved.eulerAngleDeg * .pi / 180.0
        let nx = cos(nRad), ny = sin(nRad)
        let lRad = resolved.lagrAngleDeg * .pi / 180.0
        let NLx = cos(lRad), NLy = sin(lRad)
        
        let lagrAtomsBase: [(Double, Double)] = rawAtoms
        let lagrAtomsDisp: [(Double, Double)] = rawAtoms.map { X in
            let proj = NLx * X.0 + NLy * X.1
            return (X.0 + modeAmplitude * proj * resolved.lagrMx,
                    X.1 + modeAmplitude * proj * resolved.lagrMy)
        }
        
        let eulerAtomsBase: [(Double, Double)] = rawAtoms.map { X in
            (resolved.f11 * X.0 + resolved.f12 * X.1,
             resolved.f21 * X.0 + resolved.f22 * X.1)
        }
        let eulerAtomsDisp: [(Double, Double)] = eulerAtomsBase.map { x in
            let proj = nx * x.0 + ny * x.1
            return (x.0 + modeAmplitude * proj * resolved.mx,
                    x.1 + modeAmplitude * proj * resolved.my)
        }
        
        let isStable = resolved.minDetQ > 0
        let showLagrPrimary = (frameMode == .lagrangian || frameMode == .both)
        let showEulerPrimary = (frameMode == .eulerian || frameMode == .both)
        
        // Helper to draw a 4-atom unit cell (bonds + diagonal if triangular + 4 atom spheres)
        func drawFourAtomCell(basePts: [(Double, Double)],
                              dispPts: [(Double, Double)],
                              cellColor: Color,
                              atomColor: Color,
                              isGhost: Bool,
                              isDashed: Bool,
                              labelPrefix: String) {
            var poly = Path()
            for (i, pt) in dispPts.enumerated() {
                let p = toScreen(pt.0, pt.1)
                if i == 0 { poly.move(to: p) } else { poly.addLine(to: p) }
            }
            poly.closeSubpath()
            
            if !isGhost && !isDashed {
                ctx.fill(poly, with: .color(cellColor.opacity(0.08)))
            }
            ctx.stroke(poly,
                       with: .color(cellColor.opacity(isGhost ? 0.25 : 0.80)),
                       style: StrokeStyle(lineWidth: isGhost ? 1.2 : 2.2,
                                          dash: (isGhost || isDashed) ? [5, 4] : []))
            
            // Interior short diagonal bond for triangular Bravais rhombus (Atom 2 <-> Atom 4)
            if isTri && dispPts.count == 4 {
                var diag = Path()
                diag.move(to: toScreen(dispPts[1].0, dispPts[1].1))
                diag.addLine(to: toScreen(dispPts[3].0, dispPts[3].1))
                ctx.stroke(diag,
                           with: .color(cellColor.opacity(isGhost ? 0.15 : 0.40)),
                           style: StrokeStyle(lineWidth: 1.1, dash: [3, 3]))
            }
            
            // Draw the 4 corner atoms
            for i in 0..<dispPts.count {
                let pBase = toScreen(basePts[i].0, basePts[i].1)
                let pDisp = toScreen(dispPts[i].0, dispPts[i].1)
                
                if !isGhost && modeAmplitude > 0.01 {
                    var stem = Path()
                    stem.move(to: pBase)
                    stem.addLine(to: pDisp)
                    ctx.stroke(stem, with: .color(atomColor.opacity(0.75)), lineWidth: 1.5)
                    let ghostR: CGFloat = 3.5
                    ctx.stroke(Path(ellipseIn: CGRect(x: pBase.x - ghostR, y: pBase.y - ghostR, width: 2 * ghostR, height: 2 * ghostR)),
                               with: .color(atomColor.opacity(0.45)), lineWidth: 1.0)
                }
                
                let r: CGFloat = isGhost ? 5.0 : 8.5
                let rect = CGRect(x: pDisp.x - r, y: pDisp.y - r, width: 2 * r, height: 2 * r)
                if isGhost {
                    ctx.stroke(Path(ellipseIn: rect), with: .color(atomColor.opacity(0.35)),
                               style: StrokeStyle(lineWidth: 1.2, dash: [2, 2]))
                } else {
                    ctx.fill(Path(ellipseIn: rect), with: .color(atomColor))
                    ctx.stroke(Path(ellipseIn: rect), with: .color(Color.white.opacity(0.85)), lineWidth: 1.2)
                    ctx.draw(Text("\(i + 1)")
                                .font(.system(size: 9, weight: .bold, design: .monospaced))
                                .foregroundColor(.black),
                             at: pDisp, anchor: .center)
                }
            }
        }
        
        // 3. Render Undeformed (Lagrangian) and Deformed (Eulerian) 4-Atom Cells
        if frameMode == .eulerian {
            // Faint ghost of undeformed 4 atoms + bold deformed 4 atoms
            drawFourAtomCell(basePts: lagrAtomsBase, dispPts: lagrAtomsBase,
                             cellColor: .secondary, atomColor: .secondary,
                             isGhost: true, isDashed: true, labelPrefix: "X")
            drawFourAtomCell(basePts: eulerAtomsBase, dispPts: eulerAtomsDisp,
                             cellColor: isStable ? .green : .red,
                             atomColor: isStable ? .cyan : .orange,
                             isGhost: false, isDashed: false, labelPrefix: "x")
        } else if frameMode == .lagrangian {
            // Faint ghost of deformed 4 atoms + bold undeformed 4 atoms
            drawFourAtomCell(basePts: eulerAtomsBase, dispPts: eulerAtomsBase,
                             cellColor: .secondary, atomColor: .secondary,
                             isGhost: true, isDashed: true, labelPrefix: "x")
            drawFourAtomCell(basePts: lagrAtomsBase, dispPts: lagrAtomsDisp,
                             cellColor: .purple,
                             atomColor: .purple,
                             isGhost: false, isDashed: false, labelPrefix: "X")
        } else {
            // BOTH: Undeformed 4 atoms (dashed purple) + Deformed 4 atoms (solid)
            drawFourAtomCell(basePts: lagrAtomsBase, dispPts: lagrAtomsDisp,
                             cellColor: .purple,
                             atomColor: .purple,
                             isGhost: false, isDashed: true, labelPrefix: "X")
            drawFourAtomCell(basePts: eulerAtomsBase, dispPts: eulerAtomsDisp,
                             cellColor: isStable ? .green : .red,
                             atomColor: isStable ? .cyan : .orange,
                             isGhost: false, isDashed: false, labelPrefix: "x")
        }
        
        // 4. Overlay Slip Plane(s) and Wave-Normal & Polarization Vectors on top of the 4 atoms
        let vecLen = 0.92
        let bandLen = 1.28
        
        // A. Lagrangian Vectors (N*, M*, and slip plane T* ⊥ N*)
        if showLagrPrimary {
            let Tx = -NLy, Ty = NLx
            var lSlip = Path()
            lSlip.move(to: toScreen(-bandLen * Tx, -bandLen * Ty))
            lSlip.addLine(to: toScreen(bandLen * Tx, bandLen * Ty))
            ctx.stroke(lSlip, with: .color(Color.purple.opacity(frameMode == .both ? 0.55 : 0.85)),
                       style: StrokeStyle(lineWidth: 1.8, dash: [4, 4]))
            
            // N* vector (Purple)
            let NTip = toScreen(vecLen * NLx, vecLen * NLy)
            var NPath = Path()
            NPath.move(to: center)
            NPath.addLine(to: NTip)
            ctx.stroke(NPath, with: .color(Color.purple),
                       style: StrokeStyle(lineWidth: 2.5, dash: frameMode == .both ? [4, 3] : []))
            drawArrowHead(ctx: &ctx, from: center, to: NTip, color: .purple, size: 9)
            ctx.draw(Text(String(format: "N* (%.1f°)", resolved.lagrAngleDeg))
                        .font(.system(size: 10, weight: .bold, design: .monospaced))
                        .foregroundColor(.purple),
                     at: CGPoint(x: NTip.x + 22 * CGFloat(NLx), y: NTip.y - 14 * CGFloat(NLy)),
                     anchor: .center)
            
            // M* vector (Pink/Magenta)
            let MTip = toScreen(vecLen * resolved.lagrMx, vecLen * resolved.lagrMy)
            var MPath = Path()
            MPath.move(to: center)
            MPath.addLine(to: MTip)
            ctx.stroke(MPath, with: .color(Color.pink),
                       style: StrokeStyle(lineWidth: 2.5, dash: frameMode == .both ? [4, 3] : []))
            drawArrowHead(ctx: &ctx, from: center, to: MTip, color: .pink, size: 9)
            ctx.draw(Text(String(format: "M* (%.1f°)", resolved.lagrPolDeg))
                        .font(.system(size: 10, weight: .bold, design: .monospaced))
                        .foregroundColor(.pink),
                     at: CGPoint(x: MTip.x + 22 * CGFloat(resolved.lagrMx), y: MTip.y - 14 * CGFloat(resolved.lagrMy)),
                     anchor: .center)
        }
        
        // B. Eulerian Vectors (n*, m*, and slip plane t* ⊥ n*)
        if showEulerPrimary {
            let tx = -ny, ty = nx
            var eSlip = Path()
            eSlip.move(to: toScreen(-bandLen * tx, -bandLen * ty))
            eSlip.addLine(to: toScreen(bandLen * tx, bandLen * ty))
            ctx.stroke(eSlip, with: .color(Color.cyan),
                       style: StrokeStyle(lineWidth: 2.0, dash: [7, 4]))
            
            // n* vector (Orange)
            let nTip = toScreen(vecLen * nx, vecLen * ny)
            var nVecPath = Path()
            nVecPath.move(to: center)
            nVecPath.addLine(to: nTip)
            ctx.stroke(nVecPath, with: .color(Color.orange), lineWidth: 2.6)
            drawArrowHead(ctx: &ctx, from: center, to: nTip, color: .orange, size: 9)
            ctx.draw(Text(String(format: "n* (%.1f°)", resolved.eulerAngleDeg))
                        .font(.system(size: 10, weight: .bold, design: .monospaced))
                        .foregroundColor(.orange),
                     at: CGPoint(x: nTip.x + 22 * CGFloat(nx), y: nTip.y - 14 * CGFloat(ny)),
                     anchor: .center)
            
            // m* vector (Yellow)
            let mTip = toScreen(vecLen * resolved.mx, vecLen * resolved.my)
            var mVecPath = Path()
            mVecPath.move(to: center)
            mVecPath.addLine(to: mTip)
            ctx.stroke(mVecPath, with: .color(Color.yellow), lineWidth: 2.6)
            drawArrowHead(ctx: &ctx, from: center, to: mTip, color: .yellow, size: 9)
            ctx.draw(Text(String(format: "m* (%.1f°)", resolved.eulerPolDeg))
                        .font(.system(size: 10, weight: .bold, design: .monospaced))
                        .foregroundColor(.yellow),
                     at: CGPoint(x: mTip.x + 22 * CGFloat(resolved.mx), y: mTip.y - 14 * CGFloat(resolved.my)),
                     anchor: .center)
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
}
