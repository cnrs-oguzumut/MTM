import SwiftUI
import AppKit

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

struct FieldStats {
    var min: Float, max: Float, p01: Float, p99: Float, absMax: Float
}

struct LoadedFrame {
    let url: URL
    let dataset: VTKDataset
    let geometry: MeshGeometry
    let locator: TriangleLocator
}

// MARK: - Model

@MainActor
final class SnapshotViewerModel: ObservableObject {
    @Published var source: SnapshotSource = .configurations
    @Published private(set) var files: [URL] = []
    @Published var index: Int = 0
    @Published private(set) var frame: LoadedFrame?
    @Published private(set) var isLoading = false
    @Published private(set) var errorMessage: String?
    
    @Published var fieldID: String?
    @Published var colormap: MeshColormap = .coolwarm
    @Published var rangeMode: SnapshotRangeMode = .percentile
    @Published var manualMin: Double = 0
    @Published var manualMax: Double = 1
    @Published private(set) var values: MeshValues?
    @Published private(set) var stats: FieldStats?
    
    private(set) var customURL: URL?
    private(set) var runDirectory = ""
    private var cache: [URL: LoadedFrame] = [:]
    private var cacheOrder: [URL] = []
    private var loadTask: Task<Void, Never>?
    private var valuesTask: Task<Void, Never>?
    private var valueVersion = 0
    
    var currentURL: URL? {
        if source == .custom { return customURL }
        return files.indices.contains(index) ? files[index] : nil
    }
    
    var currentField: VTKField? {
        guard let id = fieldID else { return nil }
        return frame?.dataset.fields.first { $0.id == id }
    }
    
    var displayRange: (Float, Float) {
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
        source = s
        rescan(jumpToLatest: true)
    }
    
    func open(url: URL) {
        customURL = url
        source = .custom
        files = [url]
        index = 0
        loadCurrent()
    }
    
    /// Re-lists the folder for the current source. Returns true when the file list changed.
    @discardableResult
    func rescan(jumpToLatest: Bool = false) -> Bool {
        guard source != .custom, let folder = source.folder, !runDirectory.isEmpty else {
            if source != .custom { files = []; frame = nil }
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
            frame = nil; values = nil; stats = nil; errorMessage = nil
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
                let loaded = LoadedFrame(url: url, dataset: ds, geometry: geo, locator: loc)
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
        recomputeValues()
    }
    
    // MARK: Field values & statistics (background)
    
    func recomputeValues() {
        valuesTask?.cancel()
        guard let f = frame, let field = currentField else {
            values = nil; stats = nil
            return
        }
        let tris = f.dataset.triangles
        valuesTask = Task.detached(priority: .userInitiated) { [weak self] in
            var per = [Float](repeating: 0, count: tris.count)
            if field.location == .cell {
                let nT = min(field.values.count, tris.count / 3)
                for t in 0..<nT { let v = field.values[t]; per[3 * t] = v; per[3 * t + 1] = v; per[3 * t + 2] = v }
            } else {
                let n = field.values.count
                for i in 0..<tris.count { let k = Int(tris[i]); per[i] = k < n ? field.values[k] : 0 }
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
    
    // MARK: Hover probe
    
    func probe(_ p: SIMD2<Float>) -> (triangle: Int, value: Float?)? {
        guard let f = frame, let t = f.locator.triangle(at: p) else { return nil }
        guard let field = currentField else { return (t, nil) }
        if field.location == .cell {
            return (t, t < field.values.count ? field.values[t] : nil)
        }
        let tri = f.dataset.triangles
        guard 3 * t + 2 < tri.count else { return (t, nil) }
        let n = field.values.count
        let ids = [Int(tri[3 * t]), Int(tri[3 * t + 1]), Int(tri[3 * t + 2])].filter { $0 < n }
        guard !ids.isEmpty else { return (t, nil) }
        return (t, ids.map { field.values[$0] }.reduce(0, +) / Float(ids.count))
    }
}
