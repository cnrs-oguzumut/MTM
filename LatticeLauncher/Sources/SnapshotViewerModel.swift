import SwiftUI
import AppKit
import simd

// MARK: - Sources

enum SnapshotSource: String, CaseIterable, Identifiable {
    case configurations = "Snapshots"
    case preAvalanche = "Pre-avalanche"
    case defects = "Defects"
    case eigenModes = "Soft modes"
    case custom = "Opened file"
    
    var id: String { rawValue }
    
    var folder: String? {
        switch self {
        case .configurations, .preAvalanche: return "vtk_output"
        case .defects: return "vtk_defects"
        case .eigenModes: return "eigen_modes"
        case .custom: return nil
        }
    }
    
    func accepts(_ name: String) -> Bool {
        guard name.lowercased().hasSuffix(".vtk") else { return false }
        switch self {
        case .configurations: return name.hasPrefix("configuration_")
        case .preAvalanche: return name.hasPrefix("pre_avalanche_")
        default: return true
        }
    }
}

enum SnapshotRangeMode: String, CaseIterable, Identifiable {
    case percentile = "Auto (1–99 %)"
    case full = "Full min–max"
    case symmetric = "Symmetric ±max"
    case manual = "Manual"
    var id: String { rawValue }
}

enum QuiverStrideOption: Int, CaseIterable, Identifiable {
    case auto = 0
    case every1 = 1
    case every2 = 2
    case every4 = 4
    case every8 = 8
    
    var id: Int { rawValue }
    
    var label: String {
        switch self {
        case .auto: return "Auto (Zoom-Adaptive)"
        case .every1: return "All Atoms (1×)"
        case .every2: return "Every 2 a₀"
        case .every4: return "Every 4 a₀"
        case .every8: return "Every 8 a₀"
        }
    }
    
    var shortLabel: String {
        switch self {
        case .auto: return "Auto"
        case .every1: return "1×"
        case .every2: return "2×"
        case .every4: return "4×"
        case .every8: return "8×"
        }
    }
}

enum QuiverStylePreset: String, CaseIterable, Identifiable {
    case whiteHalo = "White + Dark Halo"
    case cyanHalo = "Cyan + Dark Halo"
    case goldHalo = "Gold + Dark Halo"
    case inkLightHalo = "Ink Black + Light Halo"
    
    var id: String { rawValue }
    
    var colors: (arrow: SIMD4<Float>, halo: SIMD4<Float>) {
        switch self {
        case .whiteHalo:
            return (SIMD4(0.99, 0.99, 1.0, 0.96), SIMD4(0.02, 0.04, 0.09, 0.82))
        case .cyanHalo:
            return (SIMD4(0.30, 0.95, 1.0, 0.96), SIMD4(0.02, 0.05, 0.11, 0.84))
        case .goldHalo:
            return (SIMD4(1.0, 0.86, 0.24, 0.96), SIMD4(0.05, 0.04, 0.02, 0.84))
        case .inkLightHalo:
            return (SIMD4(0.06, 0.08, 0.13, 0.95), SIMD4(1.0, 1.0, 1.0, 0.76))
        }
    }
}

struct FieldStats {
    var min: Float, max: Float, p01: Float, p99: Float, absMax: Float
}

struct HoverProbe: Equatable {
    let triangle: Int
    let value: Float?
    let vector: SIMD2<Float>?
}

struct LoadedFrame {
    let url: URL
    let dataset: VTKDataset
    let geometry: MeshGeometry
    let locator: TriangleLocator
    let latticeSpacing: Float
    /// Symmetry period of lattice angles in this snapshot (60 deg triangular, 90 deg square)
    var latticeAnglePeriod: Float = 60
}

// MARK: - Model

@MainActor
final class SnapshotViewerModel: ObservableObject {
    @Published var source: SnapshotSource = .configurations
    @Published private(set) var files: [URL] = []
    @Published private(set) var sourceCounts: [SnapshotSource: Int] = [:]
    @Published var index: Int = 0
    @Published private(set) var frame: LoadedFrame?
    @Published private(set) var isLoading = false
    @Published private(set) var errorMessage: String?
    
    @Published var fieldID: String? { didSet { autoSelectCyclicColormap() } }
    @Published var colormap: MeshColormap = .coolwarm {
        didSet { if oldValue.isCyclic != colormap.isCyclic { recomputeValues() } }
    }
    @Published var rangeMode: SnapshotRangeMode = .percentile
    @Published var manualMin: Double = 0
    @Published var manualMax: Double = 1
    @Published private(set) var values: MeshValues?
    @Published private(set) var stats: FieldStats?
    
    // Quiver / vector arrow state
    @Published var showArrows: Bool = false
    @Published var vectorFieldID: String?
    @Published var arrowScale: Double = 1.0
    @Published var arrowStrideOption: QuiverStrideOption = .auto
    @Published var arrowMinFraction: Double = 0.02
    @Published var arrowStyle: QuiverStylePreset = .whiteHalo
    @Published private(set) var quiver: MeshQuiver?
    @Published private(set) var effectiveArrowStride: Int = 1
    
    private(set) var customURL: URL?
    private(set) var runDirectory = ""
    private var cache: [URL: LoadedFrame] = [:]
    private var cacheOrder: [URL] = []
    private var loadTask: Task<Void, Never>?
    private var valuesTask: Task<Void, Never>?
    private var quiverTask: Task<Void, Never>?
    private var zoomDebounceTask: Task<Void, Never>?
    private var valueVersion = 0
    private var quiverVersion = 0
    private var currentZoomFactor: Float = 1.0
    private var autoEnableArrowsForCustomLoad = false
    private var colormapBeforeLatticeAngle: MeshColormap?
    
    /// Lattice orientation / rotation fields: angles only defined modulo the lattice symmetry
    static let latticeAngleFieldNames: Set<String> = [
        "ElementLatticeOrientationDeg", "ElementRotationAngleDeg", "NodalLatticeRotationDeg"
    ]
    
    var isLatticeAngleField: Bool {
        guard let f = currentField else { return false }
        return Self.latticeAngleFieldNames.contains(f.name)
    }
    
    /// Period used to fold lattice angles under the cyclic colour map (nil: plain linear mapping).
    /// -P/2 and +P/2 are the same lattice orientation, so they get the same colour.
    var cyclicPeriod: Float? {
        guard colormap.isCyclic, isLatticeAngleField, let f = frame else { return nil }
        return f.latticeAnglePeriod
    }
    
    /// Lattice angle fields switch to the cyclic map; leaving them restores the previous map
    private func autoSelectCyclicColormap() {
        if isLatticeAngleField {
            if !colormap.isCyclic {
                colormapBeforeLatticeAngle = colormap
                colormap = .cyclic
            }
        } else if let previous = colormapBeforeLatticeAngle {
            colormapBeforeLatticeAngle = nil
            if colormap.isCyclic { colormap = previous }
        }
    }
    
    /// 60 deg (triangular) or 90 deg (square): from the solver's VTK title ("... lattice triangular"),
    /// else from the cell shapes (triangular-lattice cells are near-equilateral, square-lattice cells
    /// are right triangles)
    nonisolated static func detectLatticeAnglePeriod(_ ds: VTKDataset) -> Float {
        let words = ds.title.lowercased().split(whereSeparator: { $0 == " " || $0 == "\t" })
        if let i = words.firstIndex(of: "lattice"), i + 1 < words.count {
            if words[i + 1] == "square" { return 90 }
            if words[i + 1] == "triangular" { return 60 }
        }
        let nT = ds.triangleCount
        guard nT > 0 else { return 60 }
        let stride = max(1, nT / 4000)
        var maxAngles: [Float] = []
        maxAngles.reserveCapacity(nT / stride + 1)
        var t = 0
        while t < nT {
            let i0 = Int(ds.triangles[3 * t]), i1 = Int(ds.triangles[3 * t + 1]), i2 = Int(ds.triangles[3 * t + 2])
            if i0 < ds.points.count, i1 < ds.points.count, i2 < ds.points.count {
                let p = [ds.points[i0], ds.points[i1], ds.points[i2]]
                var largest: Float = 0
                for k in 0..<3 {
                    let a = p[(k + 1) % 3] - p[k], b = p[(k + 2) % 3] - p[k]
                    let den = simd_length(a) * simd_length(b)
                    if den > 0 { largest = max(largest, acos(max(-1, min(1, simd_dot(a, b) / den))) * 180 / .pi) }
                }
                maxAngles.append(largest)
            }
            t += stride
        }
        guard !maxAngles.isEmpty else { return 60 }
        maxAngles.sort()
        return maxAngles[maxAngles.count / 2] > 78 ? 90 : 60
    }
    
    var currentURL: URL? {
        if source == .custom { return customURL }
        return files.indices.contains(index) ? files[index] : nil
    }
    
    var currentField: VTKField? {
        guard let id = fieldID else { return nil }
        return frame?.dataset.fields.first { $0.id == id }
    }
    
    var availableVectorFields: [VTKVectorField] {
        frame?.dataset.vectorFields ?? []
    }
    
    var currentVectorField: VTKVectorField? {
        guard let id = vectorFieldID else { return nil }
        return availableVectorFields.first { $0.id == id }
    }
    
    var displayRange: (Float, Float) {
        if let period = cyclicPeriod { return (-period / 2, period / 2) }
        guard let s = stats else { return (0, 1) }
        var lo: Float, hi: Float
        switch rangeMode {
        case .percentile: lo = s.p01; hi = s.p99
        case .full: lo = s.min; hi = s.max
        case .symmetric: lo = -s.absMax; hi = s.absMax
        case .manual: lo = Float(manualMin); hi = Float(manualMax)
        }
        if !(hi > lo) { hi = lo + max(1e-6, abs(lo) * 1e-3) }
        return (lo, hi)
    }
    
    // MARK: Files
    
    func setRunDirectory(_ dir: String) {
        guard dir != runDirectory else { return }
        runDirectory = dir
        cache.removeAll(); cacheOrder.removeAll()
        if source != .custom { rescan(jumpToLatest: true) }
    }
    
    func selectSource(_ s: SnapshotSource) {
        let wasEigen = (source == .eigenModes)
        source = s
        if s == .eigenModes {
            showArrows = true
        } else if wasEigen {
            showArrows = false
        }
        rescan(jumpToLatest: true)
    }
    
    func open(url: URL) {
        customURL = url
        source = .custom
        autoEnableArrowsForCustomLoad = true
        files = [url]
        index = 0
        loadCurrent()
    }
    
    func refreshSourceCounts() {
        guard !runDirectory.isEmpty else {
            sourceCounts = [:]
            return
        }
        var counts: [SnapshotSource: Int] = [:]
        let fm = FileManager.default
        for s in SnapshotSource.allCases where s != .custom {
            guard let folder = s.folder else { continue }
            let dir = URL(fileURLWithPath: runDirectory).appendingPathComponent(folder)
            let names = (try? fm.contentsOfDirectory(atPath: dir.path)) ?? []
            counts[s] = names.filter(s.accepts).count
        }
        if counts != sourceCounts {
            sourceCounts = counts
        }
    }

    /// Re-lists the folder for the current source. Returns true when the file list changed.
    @discardableResult
    func rescan(jumpToLatest: Bool = false) -> Bool {
        refreshSourceCounts()
        guard source != .custom, let folder = source.folder, !runDirectory.isEmpty else {
            if source != .custom { files = []; frame = nil; quiver = nil }
            return false
        }
        let dir = URL(fileURLWithPath: runDirectory).appendingPathComponent(folder)
        let names = (try? FileManager.default.contentsOfDirectory(atPath: dir.path)) ?? []
        let newFiles = names.filter(source.accepts).sorted().map { dir.appendingPathComponent($0) }
        let changed = newFiles != files
        let previous = currentURL
        files = newFiles
        if jumpToLatest {
            index = max(0, files.count - 1)
        } else if let prev = previous, let i = files.firstIndex(of: prev) {
            index = i
        } else {
            index = min(index, max(0, files.count - 1))
        }
        if files.isEmpty {
            frame = nil; values = nil; quiver = nil; stats = nil; errorMessage = nil
        } else if changed || jumpToLatest || frame?.url != currentURL {
            loadCurrent()
        }
        return changed
    }
    
    func step(_ delta: Int) {
        guard !files.isEmpty else { return }
        let n = min(files.count - 1, max(0, index + delta))
        if n != index { index = n; loadCurrent() }
    }
    
    func jump(to i: Int) {
        guard files.indices.contains(i), i != index else { return }
        index = i
        loadCurrent()
    }
    
    // MARK: Loading (background)
    
    func loadCurrent() {
        guard let url = currentURL else { return }
        if let hit = cache[url] {
            loadTask?.cancel()
            isLoading = false
            apply(hit)
            return
        }
        loadTask?.cancel()
        isLoading = true
        errorMessage = nil
        loadTask = Task.detached(priority: .userInitiated) { [weak self] in
            do {
                let ds = try VTKReader.load(url: url)
                if Task.isCancelled { return }
                let geo = MeshGeometry(dataset: ds)
                let loc = TriangleLocator(triPositions: geo.triPositions)
                let a0 = Self.estimateLatticeSpacing(dataset: ds)
                var loaded = LoadedFrame(url: url, dataset: ds, geometry: geo, locator: loc, latticeSpacing: a0)
                loaded.latticeAnglePeriod = SnapshotViewerModel.detectLatticeAnglePeriod(ds)
                await MainActor.run {
                    guard let self, self.currentURL == url else { return }
                    self.store(loaded)
                    self.isLoading = false
                    self.apply(loaded)
                }
            } catch {
                await MainActor.run {
                    guard let self, self.currentURL == url else { return }
                    self.isLoading = false
                    self.errorMessage = "\(url.lastPathComponent): \(error.localizedDescription)"
                }
            }
        }
    }
    
    private func store(_ f: LoadedFrame) {
        cache[f.url] = f
        cacheOrder.removeAll { $0 == f.url }
        cacheOrder.append(f.url)
        while cacheOrder.count > 3 { cache.removeValue(forKey: cacheOrder.removeFirst()) }
    }
    
    private func apply(_ f: LoadedFrame) {
        frame = f
        errorMessage = nil
        let fields = f.dataset.fields
        if fieldID == nil || !fields.contains(where: { $0.id == fieldID }) {
            // Keep the same quantity by name if the location changed, else pick a sensible default
            let preferred = ["ElementEnergy", "ElementVonMisesStress", "ElementBurgersMagnitude", "mode_1_amplitude", "BurgersMagnitude"]
            fieldID = preferred.lazy.compactMap { n in fields.first { $0.name == n }?.id }.first
                ?? fields.first(where: { $0.location == .cell })?.id
                ?? fields.first?.id
        }
        
        let vecFields = f.dataset.vectorFields
        if autoEnableArrowsForCustomLoad {
            autoEnableArrowsForCustomLoad = false
            showArrows = vecFields.contains { $0.name.hasPrefix("mode_") }
        }
        syncVectorFieldForScalarField(recompute: false)
        if vectorFieldID == nil || !vecFields.contains(where: { $0.id == vectorFieldID }) {
            vectorFieldID = vecFields.first(where: { $0.name == "mode_1" })?.id ?? vecFields.first?.id
        }
        
        recomputeValues()
        recomputeQuiver()
    }
    
    /// Keeps the active vector field aligned when switching between Mode 1 / Mode 2 scalar pills (or vector components).
    func syncVectorFieldForScalarField(recompute: Bool = true) {
        guard let f = frame, let fid = fieldID else { return }
        let vecFields = f.dataset.vectorFields
        guard !vecFields.isEmpty else { return }
        let parts = fid.split(separator: ":", maxSplits: 1)
        guard parts.count == 2 else { return }
        let locPrefix = String(parts[0])
        let rawName = String(parts[1])
        
        var candidateID: String? = nil
        if rawName.hasPrefix("mode_") {
            let toks = rawName.split(separator: "_")
            if toks.count >= 2 {
                candidateID = "\(locPrefix):mode_\(toks[1])"
            }
        } else {
            for suffix in ["_magnitude", "_amplitude", "_x", "_y"] where rawName.hasSuffix(suffix) {
                let base = String(rawName.dropLast(suffix.count))
                candidateID = "\(locPrefix):\(base)"
                break
            }
        }
        if let cid = candidateID, cid != vectorFieldID, vecFields.contains(where: { $0.id == cid }) {
            vectorFieldID = cid
            if recompute { recomputeQuiver() }
        }
    }
    
    // MARK: Field values & statistics (background)
    
    func recomputeValues() {
        valuesTask?.cancel()
        guard let f = frame, let field = currentField else {
            values = nil; stats = nil
            return
        }
        let tris = f.dataset.triangles
        let period = cyclicPeriod
        valuesTask = Task.detached(priority: .userInitiated) { [weak self] in
            // Cyclic map: fold into [-P/2, P/2) so equivalent lattice angles share a colour
            let fold: (Float) -> Float = { v in
                guard let P = period, v.isFinite else { return v }
                return v - P * (v / P + 0.5).rounded(.down)
            }
            var per = [Float](repeating: 0, count: tris.count)
            if field.location == .cell {
                let nT = min(field.values.count, tris.count / 3)
                for t in 0..<nT { let v = fold(field.values[t]); per[3 * t] = v; per[3 * t + 1] = v; per[3 * t + 2] = v }
            } else {
                let n = field.values.count
                for i in 0..<tris.count { let k = Int(tris[i]); per[i] = k < n ? fold(field.values[k]) : 0 }
            }
            let finite = field.values.filter { $0.isFinite }
            var st = FieldStats(min: 0, max: 1, p01: 0, p99: 1, absMax: 1)
            if !finite.isEmpty {
                let sorted = finite.sorted()
                let q = { (p: Double) -> Float in sorted[min(sorted.count - 1, max(0, Int(p * Double(sorted.count - 1))))] }
                st = FieldStats(min: sorted.first!, max: sorted.last!, p01: q(0.01), p99: q(0.99),
                                absMax: max(abs(sorted.first!), abs(sorted.last!)))
            }
            if Task.isCancelled { return }
            await MainActor.run {
                guard let self, self.fieldID == field.id else { return }
                self.valueVersion += 1
                self.values = MeshValues(version: self.valueVersion, perVertex: per)
                self.stats = st
                if self.rangeMode != .manual {
                    let r = self.displayRange
                    self.manualMin = Double(r.0); self.manualMax = Double(r.1)
                }
            }
        }
    }
    
    // MARK: Quiver arrows (background)
    
    private static func computeAutoStride(dataset ds: VTKDataset, latticeSpacing a0: Float, zoom: Float) -> Int {
        let b = ds.bounds
        let ext = b.max - b.min
        let spanInAtoms = max(ext.x, ext.y) / max(1e-6, a0)
        let visibleAtoms = spanInAtoms / max(1.0, zoom)
        let raw = visibleAtoms / 38.0
        return max(1, Int(raw.rounded()))
    }
    
    func updateZoomFactor(_ z: Float) {
        let clamped = max(0.25, z.isFinite ? z : 1.0)
        currentZoomFactor = clamped
        guard showArrows, arrowStrideOption == .auto, let f = frame else { return }
        let nextStride = Self.computeAutoStride(dataset: f.dataset, latticeSpacing: f.latticeSpacing, zoom: clamped)
        guard nextStride != effectiveArrowStride else {
            zoomDebounceTask?.cancel()
            return
        }
        zoomDebounceTask?.cancel()
        zoomDebounceTask = Task { @MainActor [weak self] in
            try? await Task.sleep(nanoseconds: 85_000_000)
            guard !Task.isCancelled, let self else { return }
            self.recomputeQuiver()
        }
    }
    
    func applyArrowStyle() {
        guard var q = quiver else { return }
        let cols = arrowStyle.colors
        q.arrowColor = cols.arrow
        q.haloColor = cols.halo
        quiver = q
    }
    
    func recomputeQuiver() {
        quiverTask?.cancel()
        guard showArrows, let f = frame, let vf = currentVectorField else {
            quiver = nil
            return
        }
        let strideSteps: Int = (arrowStrideOption == .auto)
            ? Self.computeAutoStride(dataset: f.dataset, latticeSpacing: f.latticeSpacing, zoom: currentZoomFactor)
            : max(1, arrowStrideOption.rawValue)
        effectiveArrowStride = strideSteps
        
        let ds = f.dataset
        let a0 = f.latticeSpacing
        let scale = Float(arrowScale)
        let minFrac = Float(arrowMinFraction)
        let cols = arrowStyle.colors
        let targetURL = f.url
        let targetVecID = vf.id
        quiverVersion += 1
        let nextVersion = quiverVersion
        
        quiverTask = Task.detached(priority: .userInitiated) { [weak self] in
            let built = Self.buildQuiver(
                dataset: ds,
                vectorField: vf,
                latticeSpacing: a0,
                effectiveStride: strideSteps,
                scale: scale,
                minFraction: minFrac,
                colors: cols,
                version: nextVersion
            )
            if Task.isCancelled { return }
            await MainActor.run {
                guard let self, self.showArrows, self.frame?.url == targetURL, self.vectorFieldID == targetVecID else { return }
                self.quiver = built
            }
        }
    }
    
    private nonisolated static func estimateLatticeSpacing(dataset ds: VTKDataset) -> Float {
        let nT = ds.triangles.count / 3
        let nP = ds.points.count
        if nT > 0 {
            var lengths: [Float] = []
            lengths.reserveCapacity(min(nT, 256) * 3)
            let step = max(1, nT / 256)
            for t in stride(from: 0, to: nT, by: step) {
                let i0 = Int(ds.triangles[3 * t])
                let i1 = Int(ds.triangles[3 * t + 1])
                let i2 = Int(ds.triangles[3 * t + 2])
                guard i0 < nP, i1 < nP, i2 < nP else { continue }
                let p0 = ds.points[i0], p1 = ds.points[i1], p2 = ds.points[i2]
                let l0 = simd_length(p1 - p0)
                let l1 = simd_length(p2 - p1)
                let l2 = simd_length(p0 - p2)
                if l0.isFinite && l0 > 1e-9 { lengths.append(l0) }
                if l1.isFinite && l1 > 1e-9 { lengths.append(l1) }
                if l2.isFinite && l2 > 1e-9 { lengths.append(l2) }
            }
            if !lengths.isEmpty {
                lengths.sort()
                // 33rd percentile gives exact nearest-neighbor spacing a₀ on both square (2 legs a₀, 1 diag √2 a₀)
                // and triangular/hexagonal (all 3 edges a₀) lattices.
                return lengths[lengths.count / 3]
            }
        }
        let b = ds.bounds
        let ext = simd_max(b.max - b.min, SIMD2(1e-5, 1e-5))
        return max(1e-5, ((ext.x * ext.y) / Float(max(1, nP))).squareRoot())
    }
    
    private nonisolated static func buildQuiver(
        dataset ds: VTKDataset,
        vectorField vf: VTKVectorField,
        latticeSpacing a0: Float,
        effectiveStride: Int,
        scale: Float,
        minFraction: Float,
        colors: (arrow: SIMD4<Float>, halo: SIMD4<Float>),
        version: Int
    ) -> MeshQuiver {
        let vecs = vf.vectors
        guard !vecs.isEmpty else {
            return MeshQuiver(version: version, arrowCount: 0, haloTriangles: [], arrowTriangles: [], arrowColor: colors.arrow, haloColor: colors.halo)
        }
        
        let anchors: [SIMD2<Float>]
        if vf.location == .point {
            anchors = ds.points
        } else {
            let nT = ds.triangles.count / 3
            let nP = ds.points.count
            var cents = [SIMD2<Float>]()
            cents.reserveCapacity(nT)
            for t in 0..<nT {
                let i0 = Int(ds.triangles[3 * t])
                let i1 = Int(ds.triangles[3 * t + 1])
                let i2 = Int(ds.triangles[3 * t + 2])
                if i0 < nP && i1 < nP && i2 < nP {
                    cents.append((ds.points[i0] + ds.points[i1] + ds.points[i2]) * (1.0 / 3.0))
                } else {
                    cents.append(.zero)
                }
            }
            anchors = cents
        }
        
        let count = min(vecs.count, anchors.count)
        guard count > 0 else {
            return MeshQuiver(version: version, arrowCount: 0, haloTriangles: [], arrowTriangles: [], arrowColor: colors.arrow, haloColor: colors.halo)
        }
        
        var maxMag: Float = 0
        for i in 0..<count {
            let m = simd_length(vecs[i])
            if m.isFinite && m > maxMag { maxMag = m }
        }
        guard maxMag > 1e-20 else {
            return MeshQuiver(version: version, arrowCount: 0, haloTriangles: [], arrowTriangles: [], arrowColor: colors.arrow, haloColor: colors.halo)
        }
        
        let cutoff = max(0, minFraction) * maxMag
        let strideSteps = max(1, effectiveStride)
        let hCell = max(1e-6, Float(strideSteps) * a0)
        
        // Spatial cell binning: pick the peak-amplitude site in each 2D cell so localized grain-boundary
        // or defect-core soft modes are always preserved at every stride density.
        var selected: [Int] = []
        if strideSteps <= 1 {
            selected.reserveCapacity(count)
            for i in 0..<count {
                let m = simd_length(vecs[i])
                if m.isFinite && m >= cutoff && m > 1e-20 {
                    selected.append(i)
                }
            }
        } else {
            let b = ds.bounds
            let ext = simd_max(b.max - b.min, SIMD2(hCell, hCell))
            let nx = max(1, min(2048, Int(ceil(ext.x / hCell))))
            let ny = max(1, min(2048, Int(ceil(ext.y / hCell))))
            var bestIdx = [Int32](repeating: -1, count: nx * ny)
            var bestMag = [Float](repeating: -1, count: nx * ny)
            for i in 0..<count {
                let m = simd_length(vecs[i])
                guard m.isFinite && m >= cutoff && m > 1e-20 else { continue }
                let p = anchors[i]
                let cx = min(nx - 1, max(0, Int((p.x - b.min.x) / hCell)))
                let cy = min(ny - 1, max(0, Int((p.y - b.min.y) / hCell)))
                let bin = cy * nx + cx
                if m > bestMag[bin] {
                    bestMag[bin] = m
                    bestIdx[bin] = Int32(i)
                }
            }
            selected.reserveCapacity(nx * ny / 2)
            for idx in bestIdx where idx >= 0 {
                selected.append(Int(idx))
            }
        }
        
        let lRef = max(1e-6, 1.15 * max(0.05, scale) * hCell)
        var arrowTris = [SIMD2<Float>]()
        var haloTris = [SIMD2<Float>]()
        arrowTris.reserveCapacity(selected.count * 9)
        haloTris.reserveCapacity(selected.count * 9)
        
        for i in selected {
            let p = anchors[i]
            let v = vecs[i]
            let m = simd_length(v)
            guard m > 1e-20 else { continue }
            let d = v / m
            let n = SIMD2(-d.y, d.x)
            let r = min(1.0, max(0.0, m / maxMag))
            
            let widthFactor: Float = 0.52 + 0.48 * r
            let wShaft = 0.052 * lRef * widthFactor
            let wHead = 2.65 * wShaft
            let len = max(r * lRef, 1.85 * wHead)
            let lHead = min(0.46 * len, 2.15 * wHead)
            
            let neck = p + (len - lHead) * d
            let shaftEnd = p + max(0, len - 0.90 * lHead) * d
            let tip = p + len * d
            
            // Contrast halo geometry
            let pad = 0.58 * wShaft
            let hTail = p - pad * d
            let hShaftW = wShaft + pad
            let hNeck = neck - 0.75 * pad * d
            let hHeadW = wHead + 1.25 * pad
            let hTip = tip + 1.55 * pad * d
            
            let hs0 = hTail - hShaftW * n
            let hs1 = hTail + hShaftW * n
            let hs2 = shaftEnd + hShaftW * n
            let hs3 = shaftEnd - hShaftW * n
            haloTris.append(hs0); haloTris.append(hs1); haloTris.append(hs2)
            haloTris.append(hs0); haloTris.append(hs2); haloTris.append(hs3)
            
            let hh0 = hNeck - hHeadW * n
            let hh1 = hNeck + hHeadW * n
            haloTris.append(hh0); haloTris.append(hh1); haloTris.append(hTip)
            
            // Foreground arrow geometry
            let s0 = p - wShaft * n
            let s1 = p + wShaft * n
            let s2 = shaftEnd + wShaft * n
            let s3 = shaftEnd - wShaft * n
            arrowTris.append(s0); arrowTris.append(s1); arrowTris.append(s2)
            arrowTris.append(s0); arrowTris.append(s2); arrowTris.append(s3)
            
            let h0 = neck - wHead * n
            let h1 = neck + wHead * n
            arrowTris.append(h0); arrowTris.append(h1); arrowTris.append(tip)
        }
        
        return MeshQuiver(
            version: version,
            arrowCount: selected.count,
            haloTriangles: haloTris,
            arrowTriangles: arrowTris,
            arrowColor: colors.arrow,
            haloColor: colors.halo
        )
    }
    
    // MARK: Hover probe
    
    func probe(_ p: SIMD2<Float>) -> HoverProbe? {
        guard let f = frame, let t = f.locator.triangle(at: p) else { return nil }
        let tri = f.dataset.triangles
        let triNodeIDs: [Int] = (3 * t + 2 < tri.count)
            ? [Int(tri[3 * t]), Int(tri[3 * t + 1]), Int(tri[3 * t + 2])]
            : []
        
        var scalarVal: Float? = nil
        if let field = currentField {
            if field.location == .cell {
                scalarVal = t < field.values.count ? field.values[t] : nil
            } else {
                let ids = triNodeIDs.filter { $0 < field.values.count }
                if !ids.isEmpty {
                    scalarVal = ids.map { field.values[$0] }.reduce(0, +) / Float(ids.count)
                }
            }
        }
        
        var vecVal: SIMD2<Float>? = nil
        if showArrows, let vf = currentVectorField {
            if vf.location == .cell {
                vecVal = t < vf.vectors.count ? vf.vectors[t] : nil
            } else {
                let pts = f.dataset.points
                if let nearest = triNodeIDs
                    .filter({ $0 < vf.vectors.count && $0 < pts.count })
                    .min(by: { simd_length_squared(pts[$0] - p) < simd_length_squared(pts[$1] - p) }) {
                    vecVal = vf.vectors[nearest]
                }
            }
        }
        
        return HoverProbe(triangle: t, value: scalarVal, vector: vecVal)
    }
}
