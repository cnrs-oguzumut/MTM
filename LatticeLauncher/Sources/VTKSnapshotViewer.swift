import SwiftUI
import AppKit
import UniformTypeIdentifiers

/// Built-in viewer for the solver's `.vtk` snapshots (Results & Plots → Snapshots).
struct VTKSnapshotViewer: View {
    let runDirectory: String
    let isRunning: Bool
    
    @StateObject private var model = SnapshotViewerModel()
    @State private var showEdges = false
    @State private var showNodes = false
    @State private var approvedNodeCount: Int? = nil
    @State private var fitToken = 0
    @State private var hoverInfo: (triangle: Int, value: Float?)? = nil
    @State private var playing = false
    @State private var followLive = true
    @State private var paraViewMissing = false
    
    private let playTimer = Timer.publish(every: 0.3, on: .main, in: .common).autoconnect()
    private let liveTimer = Timer.publish(every: 2.0, on: .main, in: .common).autoconnect()
    
    private var nodeCount: Int { model.frame?.dataset.nodeCount ?? 0 }
    private var triangleCount: Int { model.frame?.dataset.triangleCount ?? 0 }
    private var nodesAutoHidden: Bool { LargeMeshPolicy.isLarge(nodeCount) && !showNodes && approvedNodeCount != nodeCount }
    
    private var appearance: MeshAppearance {
        var a = MeshAppearance()
        a.colormap = model.colormap
        let r = model.displayRange
        a.rangeMin = r.0
        a.rangeMax = r.1
        a.showEdges = showEdges
        a.showNodes = showNodes
        a.edgeColor = SIMD4(0, 0, 0, triangleCount > 40_000 ? 0.25 : 0.45)
        a.nodeSizePx = nodeCount > 20_000 ? 2.5 : 4
        return a
    }
    
    var body: some View {
        VStack(spacing: 0) {
            topBar
            Divider()
            controlBar
            Divider()
            canvas
            Divider()
            timeline
        }
        .onAppear { model.setRunDirectory(runDirectory) }
        .onChange(of: runDirectory) { model.setRunDirectory($0) }
        .onChange(of: model.frame?.url) { _ in
            hoverInfo = nil
            // Large crystal: hide nodes unless the user already approved this size
            if LargeMeshPolicy.isLarge(nodeCount) && approvedNodeCount != nodeCount { showNodes = false }
        }
        .onChange(of: triangleCount) { n in showEdges = n > 0 && n <= 8_000 }
        .onChange(of: showNodes) { on in if on { approvedNodeCount = nodeCount } }
        .onChange(of: model.fieldID) { _ in model.recomputeValues() }
        .onReceive(playTimer) { _ in
            guard playing else { return }
            if model.index >= model.files.count - 1 { playing = false } else { model.step(1) }
        }
        .onReceive(liveTimer) { _ in
            guard followLive, isRunning, model.source != .custom, !playing else { return }
            model.rescan(jumpToLatest: true)
        }
        .onReceive(NotificationCenter.default.publisher(for: .openVTKFilePanel)) { _ in openFile() }
        .alert("ParaView not found", isPresented: $paraViewMissing) {
            Button("Reveal in Finder") { if let u = model.currentURL { NSWorkspace.shared.activateFileViewerSelecting([u]) } }
            Button("OK", role: .cancel) {}
        } message: {
            Text("Install ParaView in /Applications to open snapshots there, or open the file from Finder.")
        }
    }
    
    // MARK: Top bar: source, file, open / ParaView
    
    private var topBar: some View {
        HStack(spacing: 10) {
            Picker("Source", selection: Binding(get: { model.source }, set: { model.selectSource($0) })) {
                ForEach(SnapshotSource.allCases.filter { $0 != .custom || model.source == .custom }) { s in
                    Text(s.rawValue).tag(s)
                }
            }
            .labelsHidden()
            .pickerStyle(.menu)
            .controlSize(.small)
            .frame(width: 130)
            
            if let url = model.currentURL {
                Text(url.lastPathComponent)
                    .font(.system(size: 11, weight: .semibold, design: .monospaced))
                    .foregroundColor(BlueHour.ice)
                    .lineLimit(1)
                    .truncationMode(.middle)
            }
            if let a = model.frame?.dataset.loadParameter {
                Text(String(format: "α = %.5g", a))
                    .font(.system(size: 11, design: .monospaced))
                    .foregroundColor(BlueHour.mist)
            }
            if model.frame != nil {
                Text("\(LargeMeshPolicy.formatted(nodeCount)) nodes · \(LargeMeshPolicy.formatted(triangleCount)) triangles")
                    .font(.system(size: 10))
                    .foregroundColor(BlueHour.mist)
                    .lineLimit(1)
            }
            if model.isLoading {
                ProgressView().scaleEffect(0.5).frame(width: 14, height: 14)
            }
            Spacer()
            Button(action: openFile) { Label("Open…", systemImage: "folder") }
                .controlSize(.small)
                .help("Open any .vtk file")
            Button(action: openInParaView) { Label("ParaView", systemImage: "arrow.up.forward.app") }
                .controlSize(.small)
                .disabled(model.currentURL == nil)
                .help("Open this snapshot in ParaView")
        }
        .padding(.horizontal, 12)
        .frame(height: 34)
        .background(BlueHour.controlSurface)
    }
    
    // MARK: Colour / display controls
    
    private var controlBar: some View {
        HStack(spacing: 10) {
            Text("Colour")
                .font(.system(size: 10.5, weight: .semibold))
                .foregroundColor(BlueHour.mist)
            fieldMenu
            Picker("Map", selection: $model.colormap) {
                ForEach(MeshColormap.allCases) { Text($0.rawValue).tag($0) }
            }
            .labelsHidden().pickerStyle(.menu).controlSize(.small).frame(width: 105)
            Picker("Range", selection: $model.rangeMode) {
                ForEach(SnapshotRangeMode.allCases) { Text($0.rawValue).tag($0) }
            }
            .labelsHidden().pickerStyle(.menu).controlSize(.small).frame(width: 130)
            if model.rangeMode == .manual {
                DoubleField(title: "min", value: $model.manualMin, width: 64)
                DoubleField(title: "max", value: $model.manualMax, width: 64)
            }
            Spacer()
            Toggle("Edges", isOn: $showEdges).toggleStyle(.checkbox).controlSize(.small)
            GuardedNodesToggle(isOn: $showNodes, nodeCount: nodeCount).toggleStyle(.checkbox).controlSize(.small)
            Button(action: { fitToken += 1 }) { Image(systemName: "arrow.up.left.and.down.right.magnifyingglass") }
                .controlSize(.small)
                .help("Fit to view (or double-click the mesh)")
        }
        .padding(.horizontal, 12)
        .frame(height: 32)
        .background(BlueHour.controlSurface.opacity(0.7))
    }
    
    private var fieldMenu: some View {
        let fields = model.frame?.dataset.fields ?? []
        let cells = fields.filter { $0.location == .cell }
        let points = fields.filter { $0.location == .point }
        return Menu {
            Button("Solid (no colouring)") { model.fieldID = nil }
            if !cells.isEmpty {
                Section("Per element") { ForEach(cells) { f in Button(f.displayName) { model.fieldID = f.id } } }
            }
            if !points.isEmpty {
                Section("Per node (smooth)") { ForEach(points) { f in Button(f.displayName) { model.fieldID = f.id } } }
            }
        } label: {
            Text(model.currentField?.displayName ?? "Solid")
                .font(.system(size: 11))
                .lineLimit(1)
        }
        .menuStyle(.borderedButton)
        .controlSize(.small)
        .frame(maxWidth: 240)
        .disabled(fields.isEmpty)
    }
    
    // MARK: Canvas
    
    private var canvas: some View {
        ZStack {
            MetalMeshView(
                geometry: model.frame?.geometry,
                values: model.currentField == nil ? nil : model.values,
                appearance: appearance,
                fitToken: fitToken,
                onHover: { p in
                    if let p = p { hoverInfo = model.probe(p) } else { hoverInfo = nil }
                }
            )
            
            if model.frame == nil {
                placeholder
            }
            
            VStack {
                if nodesAutoHidden && model.frame != nil {
                    NodesHiddenNote(nodeCount: nodeCount) {
                        approvedNodeCount = nodeCount
                        showNodes = true
                    }
                    .padding(.top, 8)
                }
                Spacer()
                HStack(alignment: .bottom) {
                    if let h = hoverInfo { hoverCard(h) }
                    Spacer()
                    if model.currentField != nil, model.stats != nil { colorbar }
                }
                .padding(10)
            }
        }
        .frame(maxWidth: .infinity, maxHeight: .infinity)
        .background(BlueHour.canvasBackground)
    }
    
    @ViewBuilder
    private var placeholder: some View {
        VStack(spacing: 10) {
            if model.isLoading {
                ProgressView()
                Text("Loading snapshot…").font(.caption).foregroundColor(BlueHour.mist)
            } else if let err = model.errorMessage {
                Image(systemName: "exclamationmark.triangle").font(.system(size: 30)).foregroundColor(.orange)
                Text(err).font(.caption).foregroundColor(BlueHour.mist).multilineTextAlignment(.center).frame(maxWidth: 420)
            } else {
                Image(systemName: "cube.transparent").font(.system(size: 36)).foregroundColor(BlueHour.mist)
                Text("No \(model.source.rawValue.lowercased()) found")
                    .font(.headline).foregroundColor(BlueHour.ice)
                Text(model.source.folder.map { "Looked in \($0)/ of the current output folder. Run a simulation, pick another source, or open a .vtk file." }
                     ?? "Open a .vtk file to view it.")
                    .font(.caption).foregroundColor(BlueHour.mist).multilineTextAlignment(.center).frame(maxWidth: 380)
                Button("Open .vtk File…", action: openFile).buttonStyle(.borderedProminent).controlSize(.small)
            }
        }
        .padding(20)
    }
    
    private func hoverCard(_ h: (triangle: Int, value: Float?)) -> some View {
        VStack(alignment: .leading, spacing: 3) {
            Text("Element \(h.triangle)")
                .font(.system(size: 10.5, weight: .bold, design: .monospaced))
            if let f = model.currentField, let v = h.value {
                Text("\(f.displayName): \(String(format: "%.5g", v))")
                    .font(.system(size: 10.5, design: .monospaced))
            }
        }
        .foregroundColor(BlueHour.ice)
        .padding(8)
        .background(.ultraThinMaterial, in: RoundedRectangle(cornerRadius: 7))
    }
    
    private var colorbar: some View {
        let r = model.displayRange
        return VStack(alignment: .leading, spacing: 4) {
            Text(model.currentField?.displayName ?? "")
                .font(.system(size: 10, weight: .bold))
                .lineLimit(1)
            HStack(spacing: 6) {
                Text(String(format: "%.3g", r.0)).font(.system(size: 9.5, design: .monospaced))
                LinearGradient(colors: model.colormap.gradientColors(), startPoint: .leading, endPoint: .trailing)
                    .frame(width: 150, height: 10)
                    .clipShape(RoundedRectangle(cornerRadius: 3))
                Text(String(format: "%.3g", r.1)).font(.system(size: 9.5, design: .monospaced))
            }
            if let s = model.stats, model.rangeMode == .percentile {
                Text(String(format: "data range %.3g … %.3g", s.min, s.max))
                    .font(.system(size: 9)).foregroundColor(BlueHour.mist)
            }
        }
        .foregroundColor(BlueHour.ice)
        .padding(8)
        .background(.ultraThinMaterial, in: RoundedRectangle(cornerRadius: 7))
    }
    
    // MARK: Timeline
    
    private var timeline: some View {
        HStack(spacing: 8) {
            Button(action: { playing = false; model.jump(to: 0) }) { Image(systemName: "backward.end.fill") }
            Button(action: { playing = false; model.step(-1) }) { Image(systemName: "chevron.left") }
            Button(action: {
                if model.index >= model.files.count - 1 { model.jump(to: 0) }
                playing.toggle()
            }) { Image(systemName: playing ? "pause.fill" : "play.fill") }
            Button(action: { playing = false; model.step(1) }) { Image(systemName: "chevron.right") }
            Button(action: { playing = false; model.jump(to: model.files.count - 1) }) { Image(systemName: "forward.end.fill") }
            
            if model.files.count > 1 {
                Slider(value: Binding(
                    get: { Double(model.index) },
                    set: { playing = false; model.jump(to: Int($0.rounded())) }
                ), in: 0...Double(model.files.count - 1), step: 1)
            } else {
                Spacer()
            }
            Text(model.files.isEmpty ? "–" : "\(model.index + 1) / \(model.files.count)")
                .font(.system(size: 10.5, design: .monospaced))
                .foregroundColor(BlueHour.mist)
                .frame(minWidth: 60)
            Toggle("Follow live run", isOn: $followLive)
                .toggleStyle(.checkbox)
                .controlSize(.small)
                .disabled(model.source == .custom)
                .help("While a simulation runs, jump to each new snapshot as it is written")
        }
        .buttonStyle(.borderless)
        .controlSize(.small)
        .disabled(model.files.isEmpty)
        .padding(.horizontal, 12)
        .frame(height: 34)
        .background(BlueHour.controlSurface)
    }
    
    // MARK: Actions
    
    private func openFile() {
        let panel = NSOpenPanel()
        panel.allowsMultipleSelection = false
        panel.canChooseDirectories = false
        if let vtk = UTType(filenameExtension: "vtk") { panel.allowedContentTypes = [vtk] }
        if !runDirectory.isEmpty { panel.directoryURL = URL(fileURLWithPath: runDirectory) }
        if panel.runModal() == .OK, let url = panel.url {
            playing = false
            model.open(url: url)
        }
    }
    
    private func openInParaView() {
        guard let url = model.currentURL else { return }
        let fm = FileManager.default
        var appURL: URL?
        for dir in ["/Applications", NSHomeDirectory() + "/Applications"] {
            if let names = try? fm.contentsOfDirectory(atPath: dir),
               let n = names.filter({ $0.lowercased().hasPrefix("paraview") && $0.hasSuffix(".app") }).sorted().last {
                appURL = URL(fileURLWithPath: dir).appendingPathComponent(n)
                break
            }
        }
        guard let app = appURL else { paraViewMissing = true; return }
        NSWorkspace.shared.open([url], withApplicationAt: app, configuration: NSWorkspace.OpenConfiguration())
    }
}
