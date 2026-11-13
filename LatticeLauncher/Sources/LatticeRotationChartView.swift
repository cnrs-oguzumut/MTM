import SwiftUI
import Charts

/// Shared grain colour palette (Geometry Studio canvas, rotation chart)
public enum GrainPalette {
    static let colors: [Color] = [
        Color(red: 0.56, green: 0.78, blue: 1.00),
        Color(red: 1.00, green: 0.76, blue: 0.42),
        Color(red: 0.62, green: 0.92, blue: 0.68),
        Color(red: 0.92, green: 0.62, blue: 0.95),
        Color(red: 1.00, green: 0.60, blue: 0.62),
        Color(red: 0.52, green: 0.92, blue: 0.90),
        Color(red: 0.95, green: 0.90, blue: 0.52),
        Color(red: 0.75, green: 0.68, blue: 1.00)
    ]

    public static func color(_ gid: Int) -> Color {
        colors[abs(gid) % colors.count]
    }
}

/// One row of lattice_rotation_log.csv (angles in degrees, counter-clockwise positive).
/// Rotation = lattice orientation minus the grain's starting angle, measured from current coordinates;
/// Net = Rotation minus the rigid rotation of the applied macroscopic deformation.
public struct LatticeRotationPoint: Identifiable, Equatable {
    public var id: Int { iteration }
    public let iteration: Int
    public let alpha: Double
    public let appliedRotation: Double
    public let meanRotation: Double
    public let meanRotationNet: Double
    public let stdRotation: Double
    public let maxAbsRotation: Double
    public let fracAbove2Deg: Double
    public let numElements: Int
    public let numBoundaryElements: Int
    public let numUnrecognizedElements: Int
    /// Per-grain mean rotation (nil when the grain had no measurable element at this step)
    public let grainRotation: [Double?]
    /// Per-grain current orientation φ = θ_start + rotation (continuous, not wrapped); empty for older logs
    public let grainOrientation: [Double?]
}

public enum LatticeRotationLog {
    /// Parses lattice_rotation_log.csv by column name (robust to added columns).
    public static func parse(_ content: String) -> [LatticeRotationPoint] {
        var points: [LatticeRotationPoint] = []
        var col: [String: Int] = [:]
        var grainCols: [Int] = []
        var orientCols: [Int] = []
        for raw in content.components(separatedBy: .newlines) {
            let line = raw.trimmingCharacters(in: .whitespaces)
            if line.isEmpty { continue }
            let parts = line.components(separatedBy: ",").map { $0.trimmingCharacters(in: .whitespaces) }
            if parts.first == "Iteration" {
                col = [:]
                for (i, name) in parts.enumerated() { col[name] = i }
                grainCols = parts.enumerated()
                    .filter { $0.element.hasPrefix("Grain") && $0.element.hasSuffix("RotationDeg") }
                    .map { $0.offset }
                orientCols = parts.enumerated()
                    .filter { $0.element.hasPrefix("Grain") && $0.element.hasSuffix("OrientationDeg") }
                    .map { $0.offset }
                continue
            }
            func d(_ name: String) -> Double? {
                guard let i = col[name], i < parts.count, let v = Double(parts[i]), v.isFinite else { return nil }
                return v
            }
            guard let iterIdx = col["Iteration"], iterIdx < parts.count, let iter = Int(parts[iterIdx]),
                  let alpha = d("Alpha"), let mean = d("MeanRotationDeg") else { continue }
            func optional(_ i: Int) -> Double? {
                guard i < parts.count, let v = Double(parts[i]), v.isFinite else { return nil }
                return v
            }
            let grains: [Double?] = grainCols.map(optional)
            let orientations: [Double?] = orientCols.map(optional)
            points.append(LatticeRotationPoint(
                iteration: iter,
                alpha: alpha,
                appliedRotation: d("AppliedRotationDeg") ?? 0,
                meanRotation: mean,
                meanRotationNet: d("MeanRotationNetDeg") ?? mean,
                stdRotation: d("StdRotationDeg") ?? 0,
                maxAbsRotation: d("MaxAbsRotationDeg") ?? 0,
                fracAbove2Deg: d("FracAbsRotationAbove2Deg") ?? 0,
                numElements: Int(d("NumElements") ?? 0),
                numBoundaryElements: Int(d("NumBoundaryElements") ?? 0),
                numUnrecognizedElements: Int(d("NumUnrecognizedElements") ?? 0),
                grainRotation: grains,
                grainOrientation: orientations
            ))
        }
        // A restarted run appends rows; keep the last row per iteration, in order
        var seen = Set<Int>()
        var unique: [LatticeRotationPoint] = []
        for p in points.reversed() where !seen.contains(p.iteration) {
            seen.insert(p.iteration)
            unique.append(p)
        }
        return unique.reversed()
    }
}

/// Lattice rotation ω(α): per-grain mean rotations, overall mean ± spread, and the applied rigid rotation.
struct LatticeRotationChartView: View {
    let points: [LatticeRotationPoint]

    private enum Quantity: String, CaseIterable, Identifiable {
        case rotation = "Rotation ω"
        case orientation = "Orientation φ"
        var id: String { rawValue }
    }

    @State private var quantity: Quantity

    init(points: [LatticeRotationPoint], startWithOrientation: Bool = false) {
        self.points = points
        _quantity = State(initialValue: startWithOrientation ? .orientation : .rotation)
    }
    @State private var subtractApplied: Bool = false
    @State private var showGrains: Bool = true
    @State private var hoverAlpha: Double? = nil

    private struct Series: Identifiable {
        let id: String
        let label: String
        let color: Color
        let dashed: Bool
        let values: [(alpha: Double, value: Double)]
    }

    private var grainCount: Int { points.map(\.grainRotation.count).max() ?? 0 }

    private var hasOrientation: Bool { points.contains { !$0.grainOrientation.isEmpty } }
    private var showsOrientation: Bool { quantity == .orientation && hasOrientation }

    private var series: [Series] {
        var out: [Series] = []
        if showsOrientation {
            // Current orientation of each grain (θ_start + rotation); a single line for a single crystal
            let n = points.map(\.grainOrientation.count).max() ?? 0
            for g in 0..<n {
                let vals = points.compactMap { p -> (alpha: Double, value: Double)? in
                    guard g < p.grainOrientation.count, let v = p.grainOrientation[g] else { return nil }
                    return (alpha: p.alpha, value: v)
                }
                out.append(Series(id: "o\(g)", label: "G\(g)", color: GrainPalette.color(g), dashed: false, values: vals))
            }
            return out
        }
        let shift: (LatticeRotationPoint) -> Double = { subtractApplied ? $0.appliedRotation : 0 }
        if showGrains && grainCount > 1 {
            for g in 0..<grainCount {
                let vals = points.compactMap { p -> (Double, Double)? in
                    guard g < p.grainRotation.count, let v = p.grainRotation[g] else { return nil }
                    return (p.alpha, v - shift(p))
                }
                out.append(Series(id: "g\(g)", label: "G\(g)", color: GrainPalette.color(g), dashed: false,
                                  values: vals.map { (alpha: $0.0, value: $0.1) }))
            }
        }
        out.append(Series(id: "mean", label: "Mean", color: .white, dashed: false,
                          values: points.map { (alpha: $0.alpha, value: $0.meanRotation - shift($0)) }))
        if !subtractApplied {
            out.append(Series(id: "applied", label: "Applied", color: .secondary, dashed: true,
                              values: points.map { (alpha: $0.alpha, value: $0.appliedRotation) }))
        }
        return out
    }

    private var alphaDomain: ClosedRange<Double> {
        let lo = points.map(\.alpha).min() ?? 0
        let hi = points.map(\.alpha).max() ?? 1
        return hi > lo ? lo...hi : (lo - 0.5)...(lo + 0.5)
    }

    private var hoverPoint: LatticeRotationPoint? {
        guard let a = hoverAlpha else { return nil }
        return points.min { abs($0.alpha - a) < abs($1.alpha - a) }
    }

    var body: some View {
        VStack(alignment: .leading, spacing: 6) {
            HStack(spacing: 10) {
                Text(showsOrientation ? "Grain Orientation φ(α)" : "Lattice Rotation ω(α)")
                    .font(.system(size: 12, weight: .bold))
                Picker("", selection: $quantity) {
                    ForEach(Quantity.allCases) { q in Text(q.rawValue).tag(q) }
                }
                .pickerStyle(.segmented)
                .labelsHidden()
                .controlSize(.small)
                .fixedSize()
                .help(hasOrientation
                      ? "Rotation ω: how far each grain's lattice has turned since the start. Orientation φ: each grain's current lattice angle (starting angle θ + ω)."
                      : "This run's log has no orientation columns (made before φ was logged): re-run to plot φ.")
                if showsOrientation {
                    EmptyView()
                } else if grainCount > 1 {
                    Toggle("Per grain", isOn: $showGrains)
                        .toggleStyle(.checkbox)
                        .font(.caption)
                        .help("Show each grain's area-weighted mean rotation (grain-boundary elements excluded)")
                }
                if !showsOrientation {
                    Toggle("Net of applied rotation", isOn: $subtractApplied)
                        .toggleStyle(.checkbox)
                        .font(.caption)
                        .help("Subtract the rigid rotation of the applied macroscopic deformation (−atan(γ/2) in simple shear), leaving the rotation caused by plasticity and grain boundaries")
                }
                Spacer()
                if quantity == .orientation && !hasOrientation {
                    Text("No φ in this log — re-run to record it")
                        .font(.system(size: 10))
                        .foregroundColor(.orange)
                } else if showsOrientation, let p = hoverPoint ?? points.last {
                    Text(String(format: "α = %.4f · φ = θ_start + ω (continuous, not wrapped)", p.alpha))
                        .font(.system(size: 10, design: .monospaced))
                        .foregroundColor(.secondary)
                } else if let p = hoverPoint ?? points.last {
                    Text(String(format: "α = %.4f · mean %.2f° · spread %.2f° · |ω|max %.1f° · %.0f%% > 2°",
                                p.alpha, subtractApplied ? p.meanRotationNet : p.meanRotation,
                                p.stdRotation, p.maxAbsRotation, 100 * p.fracAbove2Deg))
                        .font(.system(size: 10, design: .monospaced))
                        .foregroundColor(.secondary)
                }
            }

            Chart {
                if !showsOrientation && (!showGrains || grainCount <= 1) {
                    ForEach(points) { p in
                        let shift = subtractApplied ? p.appliedRotation : 0
                        AreaMark(
                            x: .value("α", p.alpha),
                            yStart: .value("ω − σ", p.meanRotation - shift - p.stdRotation),
                            yEnd: .value("ω + σ", p.meanRotation - shift + p.stdRotation)
                        )
                        .foregroundStyle(Color.white.opacity(0.12))
                    }
                }
                ForEach(series) { s in
                    ForEach(Array(s.values.enumerated()), id: \.offset) { _, v in
                        LineMark(
                            x: .value("α", v.alpha),
                            y: .value("ω (°)", v.value),
                            series: .value("Series", s.label)
                        )
                        .foregroundStyle(s.color)
                        .lineStyle(StrokeStyle(lineWidth: s.id == "mean" ? 2.4 : 1.6, dash: s.dashed ? [5, 3] : []))
                    }
                }
                if !showsOrientation {
                    RuleMark(y: .value("zero", 0))
                        .foregroundStyle(Color.secondary.opacity(0.35))
                }
                if let p = hoverPoint {
                    RuleMark(x: .value("α", p.alpha))
                        .foregroundStyle(Color.secondary.opacity(0.5))
                }
            }
            .chartXScale(domain: alphaDomain)
            .chartXAxisLabel("Load parameter α")
            .chartYAxisLabel(showsOrientation ? "Lattice orientation φ (°)"
                             : (subtractApplied ? "ω − ω_applied (°)" : "Lattice rotation ω (°)"))
            .chartOverlay { proxy in
                GeometryReader { geo in
                    Rectangle().fill(Color.clear).contentShape(Rectangle())
                        .onContinuousHover { phase in
                            switch phase {
                            case .active(let loc):
                                if let plot = proxy.plotFrame {
                                    let x = loc.x - geo[plot].origin.x
                                    hoverAlpha = proxy.value(atX: x, as: Double.self)
                                }
                            case .ended:
                                hoverAlpha = nil
                            }
                        }
                }
            }

            // Legend
            ScrollView(.horizontal, showsIndicators: false) {
                HStack(spacing: 12) {
                    ForEach(series) { s in
                        HStack(spacing: 4) {
                            Capsule().fill(s.color).frame(width: 14, height: 3)
                            Text(legendText(for: s))
                                .font(.system(size: 10, design: .monospaced))
                                .foregroundColor(.secondary)
                        }
                    }
                    if let p = hoverPoint ?? points.last, p.numUnrecognizedElements > 0 {
                        Text("· \(p.numUnrecognizedElements) core/transformed elements excluded")
                            .font(.system(size: 10))
                            .foregroundColor(.secondary)
                    }
                }
            }
        }
        .padding(12)
        .frame(maxWidth: .infinity, maxHeight: .infinity)
        .background(BlueHour.canvasBackground)
    }

    private func legendText(for s: Series) -> String {
        guard let p = hoverPoint ?? points.last else { return s.label }
        let a = p.alpha
        if let v = s.values.min(by: { abs($0.alpha - a) < abs($1.alpha - a) }) {
            return String(format: "%@ %.2f°", s.label, v.value)
        }
        return s.label
    }
}
