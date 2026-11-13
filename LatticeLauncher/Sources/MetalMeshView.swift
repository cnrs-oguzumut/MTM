import SwiftUI
import AppKit
import MetalKit

// MARK: - SwiftUI wrapper

/// GPU mesh renderer. Geometry is uploaded once; pan / zoom / colour-range changes only update uniforms.
/// Redraws only when something changes (no continuous render loop).
/// Interaction: scroll / pinch / right-drag = zoom (around cursor), left-drag = pan, Shift+drag (or box mode) = zoom to box, double-click = fit.
public struct MetalMeshView: NSViewRepresentable {
    public var geometry: MeshGeometry?
    public var values: MeshValues?
    public var quiver: MeshQuiver? = nil
    public var appearance: MeshAppearance
    /// Increment to request "fit to view".
    public var fitToken: Int = 0
    /// Increment to step zoom in/out around the view centre.
    public var zoomInToken: Int = 0
    public var zoomOutToken: Int = 0
    /// When true, left-dragging draws a ParaView-style zoom box.
    public var boxZoomActive: Bool = false
    public var onHover: ((SIMD2<Float>?) -> Void)? = nil
    public var onZoomChange: ((Float) -> Void)? = nil
    /// Reports marquee rectangle in top-left (SwiftUI) view coordinates while dragging a zoom box.
    public var onBoxZoomChange: ((CGRect?) -> Void)? = nil
    public var onBoxZoomComplete: (() -> Void)? = nil
    
    public init(
        geometry: MeshGeometry?,
        values: MeshValues?,
        quiver: MeshQuiver? = nil,
        appearance: MeshAppearance,
        fitToken: Int = 0,
        zoomInToken: Int = 0,
        zoomOutToken: Int = 0,
        boxZoomActive: Bool = false,
        onHover: ((SIMD2<Float>?) -> Void)? = nil,
        onZoomChange: ((Float) -> Void)? = nil,
        onBoxZoomChange: ((CGRect?) -> Void)? = nil,
        onBoxZoomComplete: (() -> Void)? = nil
    ) {
        self.geometry = geometry
        self.values = values
        self.quiver = quiver
        self.appearance = appearance
        self.fitToken = fitToken
        self.zoomInToken = zoomInToken
        self.zoomOutToken = zoomOutToken
        self.boxZoomActive = boxZoomActive
        self.onHover = onHover
        self.onZoomChange = onZoomChange
        self.onBoxZoomChange = onBoxZoomChange
        self.onBoxZoomComplete = onBoxZoomComplete
    }
    
    public func makeCoordinator() -> MeshRenderer { MeshRenderer() }
    
    public func makeNSView(context: Context) -> NSView {
        guard let device = MTLCreateSystemDefaultDevice(), context.coordinator.setup(device: device) else {
            let label = NSTextField(labelWithString: "Metal is not available on this Mac – mesh view disabled.")
            label.textColor = .secondaryLabelColor
            return label
        }
        let view = InteractiveMTKView(frame: .zero, device: device)
        view.renderer = context.coordinator
        view.delegate = context.coordinator
        view.colorPixelFormat = .bgra8Unorm
        view.sampleCount = context.coordinator.sampleCount
        view.isPaused = true
        view.enableSetNeedsDisplay = true
        view.framebufferOnly = true
        context.coordinator.view = view
        return view
    }
    
    public func updateNSView(_ nsView: NSView, context: Context) {
        let r = context.coordinator
        r.onHover = onHover
        r.onZoomChange = onZoomChange
        r.onBoxZoomChange = onBoxZoomChange
        r.onBoxZoomComplete = onBoxZoomComplete
        if let v = nsView as? InteractiveMTKView {
            v.boxZoomActive = boxZoomActive
        }
        r.update(
            geometry: geometry,
            values: values,
            quiver: quiver,
            appearance: appearance,
            fitToken: fitToken,
            zoomInToken: zoomInToken,
            zoomOutToken: zoomOutToken
        )
    }
}

// MARK: - Interactive MTKView

final class InteractiveMTKView: MTKView {
    weak var renderer: MeshRenderer?
    var boxZoomActive: Bool = false
    private var lastDrag: CGPoint?
    private var rightDragAnchor: CGPoint?
    private var lastRightY: CGFloat?
    private var boxStart: CGPoint?
    private var tracking: NSTrackingArea?
    
    override var acceptsFirstResponder: Bool { true }
    override func acceptsFirstMouse(for event: NSEvent?) -> Bool { true }
    
    override func updateTrackingAreas() {
        super.updateTrackingAreas()
        if let t = tracking { removeTrackingArea(t) }
        let t = NSTrackingArea(rect: bounds, options: [.mouseMoved, .mouseEnteredAndExited, .activeInKeyWindow, .inVisibleRect], owner: self)
        addTrackingArea(t)
        tracking = t
    }
    
    private func local(_ e: NSEvent) -> CGPoint { convert(e.locationInWindow, from: nil) }
    
    private func swiftUIRect(from a: CGPoint, to b: CGPoint) -> CGRect {
        let minX = min(a.x, b.x)
        let maxX = max(a.x, b.x)
        let minY = min(a.y, b.y)
        let maxY = max(a.y, b.y)
        let h = bounds.height
        return CGRect(x: minX, y: h - maxY, width: maxX - minX, height: maxY - minY)
    }
    
    override func mouseDown(with e: NSEvent) {
        if e.clickCount == 2 {
            boxStart = nil
            renderer?.onBoxZoomChange?(nil)
            renderer?.fit()
            return
        }
        let p = local(e)
        if e.modifierFlags.contains(.control) {
            rightDragAnchor = p
            lastRightY = p.y
            return
        }
        if boxZoomActive || e.modifierFlags.contains(.shift) {
            boxStart = p
            renderer?.onBoxZoomChange?(CGRect(x: p.x, y: bounds.height - p.y, width: 0, height: 0))
            return
        }
        lastDrag = p
    }
    
    override func mouseDragged(with e: NSEvent) {
        let p = local(e)
        if let anchor = rightDragAnchor, let ly = lastRightY {
            let dy = Float(p.y - ly) * 0.015
            renderer?.zoom(by: exp(dy), around: anchor)
            lastRightY = p.y
            return
        }
        if let start = boxStart {
            renderer?.onBoxZoomChange?(swiftUIRect(from: start, to: p))
            return
        }
        if let l = lastDrag { renderer?.pan(dx: Float(p.x - l.x), dy: Float(p.y - l.y)) }
        lastDrag = p
    }
    
    override func mouseUp(with e: NSEvent) {
        let p = local(e)
        if let start = boxStart {
            boxStart = nil
            renderer?.onBoxZoomChange?(nil)
            let rect = CGRect(
                x: min(start.x, p.x),
                y: min(start.y, p.y),
                width: abs(p.x - start.x),
                height: abs(p.y - start.y)
            )
            if rect.width >= 6 && rect.height >= 6 {
                renderer?.zoomToViewRect(rect)
                renderer?.onBoxZoomComplete?()
            }
            return
        }
        lastDrag = nil
        rightDragAnchor = nil
        lastRightY = nil
    }
    
    override func rightMouseDown(with e: NSEvent) {
        let p = local(e)
        rightDragAnchor = p
        lastRightY = p.y
    }
    
    override func rightMouseDragged(with e: NSEvent) {
        let p = local(e)
        if let anchor = rightDragAnchor, let ly = lastRightY {
            let dy = Float(p.y - ly) * 0.015
            renderer?.zoom(by: exp(dy), around: anchor)
            lastRightY = p.y
        }
    }
    
    override func rightMouseUp(with e: NSEvent) {
        rightDragAnchor = nil
        lastRightY = nil
    }
    
    override func mouseMoved(with e: NSEvent) { renderer?.hover(at: local(e)) }
    override func mouseExited(with e: NSEvent) { renderer?.onHover?(nil) }
    
    override func scrollWheel(with e: NSEvent) {
        if e.modifierFlags.contains(.option) {
            // Option + scroll → pan
            renderer?.pan(dx: Float(e.scrollingDeltaX), dy: Float(-e.scrollingDeltaY))
        } else {
            // ParaView-style scroll wheel / two-finger trackpad scroll → zoom around cursor
            let dy = e.hasPreciseScrollingDeltas ? e.scrollingDeltaY * 0.012 : e.scrollingDeltaY * 0.12
            guard abs(dy) > 1e-6 else { return }
            renderer?.zoom(by: Float(exp(dy)), around: local(e))
        }
    }
    override func magnify(with e: NSEvent) {
        renderer?.zoom(by: Float(1 + e.magnification), around: local(e))
    }
}

// MARK: - Renderer

public final class MeshRenderer: NSObject, MTKViewDelegate {
    struct Uniforms {
        var center: SIMD2<Float>
        var scale: SIMD2<Float>
        var vmin: Float
        var vmax: Float
        var pointSize: Float
        var useValues: Int32
        var color: SIMD4<Float>
    }
    
    weak var view: MTKView?
    var onHover: ((SIMD2<Float>?) -> Void)?
    var onZoomChange: ((Float) -> Void)?
    var onBoxZoomChange: ((CGRect?) -> Void)?
    var onBoxZoomComplete: (() -> Void)?
    let sampleCount = 4
    
    private var device: MTLDevice!
    private var queue: MTLCommandQueue!
    private var fillPipeline: MTLRenderPipelineState!
    private var solidPipeline: MTLRenderPipelineState!
    private var pointPipeline: MTLRenderPipelineState!
    private var lutTexture: MTLTexture?
    private var lutColormap: MeshColormap?
    
    private var triBuffer: MTLBuffer?, triCount = 0
    private var edgeBuffer: MTLBuffer?, edgeCount = 0
    private var nodeBuffer: MTLBuffer?, nodeCount = 0
    private var nodeColorBuffer: MTLBuffer?
    private var outlineBuffer: MTLBuffer?, outlineCount = 0
    private var valueBuffer: MTLBuffer?
    private var quiverHaloBuffer: MTLBuffer?, quiverHaloCount = 0
    private var quiverArrowBuffer: MTLBuffer?, quiverArrowCount = 0
    private var quiverArrowColor = SIMD4<Float>(0.99, 0.99, 1.0, 0.96)
    private var quiverHaloColor = SIMD4<Float>(0.02, 0.04, 0.09, 0.80)
    
    private var geometryID: UUID?
    private var valuesVersion: Int?
    private var quiverVersion: Int?
    private var appearance = MeshAppearance()
    private var fitToken = -1
    private var zoomInToken = 0
    private var zoomOutToken = 0
    private var bounds: (min: SIMD2<Float>, max: SIMD2<Float>) = (.zero, SIMD2(1, 1))
    
    // View transform: world point shown at the view centre, and points-per-world-unit.
    private var center = SIMD2<Float>(0.5, 0.5)
    private var zoom: Float = 100
    private var lastReportedZoomPct: Int = -1
    private var needsFit = true
    
    func setup(device: MTLDevice) -> Bool {
        self.device = device
        guard let q = device.makeCommandQueue() else { return false }
        queue = q
        do {
            let lib = try device.makeLibrary(source: MeshRenderer.shaderSource, options: nil)
            func pipe(_ v: String, _ f: String) throws -> MTLRenderPipelineState {
                let d = MTLRenderPipelineDescriptor()
                d.vertexFunction = lib.makeFunction(name: v)
                d.fragmentFunction = lib.makeFunction(name: f)
                d.rasterSampleCount = sampleCount
                let ca = d.colorAttachments[0]!
                ca.pixelFormat = .bgra8Unorm
                ca.isBlendingEnabled = true
                ca.rgbBlendOperation = .add
                ca.alphaBlendOperation = .add
                ca.sourceRGBBlendFactor = .sourceAlpha
                ca.sourceAlphaBlendFactor = .one
                ca.destinationRGBBlendFactor = .oneMinusSourceAlpha
                ca.destinationAlphaBlendFactor = .oneMinusSourceAlpha
                return try device.makeRenderPipelineState(descriptor: d)
            }
            fillPipeline = try pipe("fill_vertex", "fill_fragment")
            solidPipeline = try pipe("solid_vertex", "solid_fragment")
            pointPipeline = try pipe("point_vertex", "point_fragment")
        } catch {
            NSLog("MeshRenderer shader error: \(error)")
            return false
        }
        return true
    }
    
    private func makeBuffer<T>(_ array: [T]) -> MTLBuffer? {
        guard !array.isEmpty else { return nil }
        return array.withUnsafeBytes { device.makeBuffer(bytes: $0.baseAddress!, length: $0.count, options: .storageModeShared) }
    }
    
    // MARK: Updates from SwiftUI
    
    func update(
        geometry: MeshGeometry?,
        values: MeshValues?,
        quiver: MeshQuiver? = nil,
        appearance newAppearance: MeshAppearance,
        fitToken newFit: Int,
        zoomInToken newZoomIn: Int = 0,
        zoomOutToken newZoomOut: Int = 0
    ) {
        var dirty = false
        if geometry?.id != geometryID {
            geometryID = geometry?.id
            upload(geometry)
            valuesVersion = nil
            dirty = true
        }
        if values?.version != valuesVersion {
            valuesVersion = values?.version
            if let v = values, v.perVertex.count == triCount {
                valueBuffer = makeBuffer(v.perVertex)
            } else {
                valueBuffer = nil
            }
            dirty = true
        }
        if quiver?.version != quiverVersion {
            quiverVersion = quiver?.version
            if let q = quiver {
                quiverHaloBuffer = makeBuffer(q.haloTriangles)
                quiverHaloCount = q.haloTriangles.count
                quiverArrowBuffer = makeBuffer(q.arrowTriangles)
                quiverArrowCount = q.arrowTriangles.count
                quiverArrowColor = q.arrowColor
                quiverHaloColor = q.haloColor
            } else {
                quiverHaloBuffer = nil; quiverHaloCount = 0
                quiverArrowBuffer = nil; quiverArrowCount = 0
            }
            dirty = true
        } else if let q = quiver, (q.arrowColor != quiverArrowColor || q.haloColor != quiverHaloColor) {
            quiverArrowColor = q.arrowColor
            quiverHaloColor = q.haloColor
            dirty = true
        }
        if newAppearance != appearance { appearance = newAppearance; dirty = true }
        if newFit != fitToken { fitToken = newFit; needsFit = true; dirty = true }
        if newZoomIn != zoomInToken {
            zoomInToken = newZoomIn
            let s = viewSize
            zoom(by: 1.35, around: CGPoint(x: CGFloat(s.x * 0.5), y: CGFloat(s.y * 0.5)))
        }
        if newZoomOut != zoomOutToken {
            zoomOutToken = newZoomOut
            let s = viewSize
            zoom(by: 1.0 / 1.35, around: CGPoint(x: CGFloat(s.x * 0.5), y: CGFloat(s.y * 0.5)))
        }
        if dirty { view?.needsDisplay = true }
    }
    
    private func upload(_ g: MeshGeometry?) {
        guard let g = g else {
            triBuffer = nil; edgeBuffer = nil; nodeBuffer = nil; nodeColorBuffer = nil; outlineBuffer = nil
            triCount = 0; edgeCount = 0; nodeCount = 0; outlineCount = 0
            return
        }
        triBuffer = makeBuffer(g.triPositions)
        triCount = g.triPositions.count
        edgeBuffer = makeBuffer(g.edgePositions)
        edgeCount = g.edgePositions.count
        nodeBuffer = makeBuffer(g.nodePositions)
        nodeCount = g.nodePositions.count
        nodeColorBuffer = (g.nodeColors?.count == g.nodePositions.count) ? makeBuffer(g.nodeColors!) : nil
        outlineBuffer = makeBuffer(g.outlineSegments)
        outlineCount = g.outlineSegments.count
        let newBounds = g.bounds
        // Re-fit only when the extent changes noticeably (keeps the view while scrubbing frames)
        let ext = newBounds.max - newBounds.min, oldExt = bounds.max - bounds.min
        if simd_length(ext - oldExt) > 0.05 * max(1e-6, simd_length(oldExt)) || simd_length(newBounds.min - bounds.min) > 0.05 * max(1e-6, simd_length(oldExt)) {
            needsFit = true
        }
        bounds = newBounds
    }
    
    // MARK: Interaction
    
    private var viewSize: SIMD2<Float> {
        guard let v = view else { return SIMD2(1, 1) }
        return SIMD2(Float(max(1, v.bounds.width)), Float(max(1, v.bounds.height)))
    }
    
    private var fitZoomValue: Float {
        let size = viewSize
        let ext = simd_max(bounds.max - bounds.min, SIMD2(1e-6, 1e-6))
        return max(1e-6, min(size.x / ext.x, size.y / ext.y) * 0.92)
    }
    
    private func notifyZoom() {
        let rel = zoom / fitZoomValue
        let pct = Int(round(rel * 100))
        guard pct != lastReportedZoomPct else { return }
        lastReportedZoomPct = pct
        DispatchQueue.main.async { [weak self] in
            self?.onZoomChange?(rel)
        }
    }
    
    func fit() {
        zoom = fitZoomValue
        center = (bounds.min + bounds.max) * 0.5
        needsFit = false
        view?.needsDisplay = true
        notifyZoom()
    }
    
    func pan(dx: Float, dy: Float) {
        center -= SIMD2(dx, dy) / zoom
        view?.needsDisplay = true
    }
    
    func zoom(by factor: Float, around p: CGPoint) {
        guard factor.isFinite, factor > 0 else { return }
        let size = viewSize
        let offset = SIMD2(Float(p.x), Float(p.y)) - size * 0.5
        let world = center + offset / zoom
        let fz = fitZoomValue
        let minZoom = fz * 0.1
        let maxZoom = fz * 250.0
        zoom = min(maxZoom, max(minZoom, zoom * factor))
        center = world - offset / zoom
        view?.needsDisplay = true
        notifyZoom()
    }
    
    func zoomToViewRect(_ r: CGRect) {
        let size = viewSize
        guard r.width >= 6, r.height >= 6 else { return }
        let mid = SIMD2(Float(r.midX), Float(r.midY))
        let offset = mid - size * 0.5
        let newCenter = center + offset / zoom
        let factor = min(size.x / Float(r.width), size.y / Float(r.height)) * 0.95
        let fz = fitZoomValue
        zoom = min(fz * 250.0, max(fz * 0.1, zoom * factor))
        center = newCenter
        needsFit = false
        view?.needsDisplay = true
        notifyZoom()
    }
    
    func hover(at p: CGPoint) {
        guard let cb = onHover else { return }
        let offset = SIMD2(Float(p.x), Float(p.y)) - viewSize * 0.5
        cb(center + offset / zoom)
    }
    
    // MARK: Drawing
    
    public func mtkView(_ view: MTKView, drawableSizeWillChange size: CGSize) {
        view.needsDisplay = true
    }
    
    public func draw(in view: MTKView) {
        if needsFit && view.bounds.width > 2 && triCount + nodeCount + outlineCount > 0 { fit() }
        guard let pass = view.currentRenderPassDescriptor, let drawable = view.currentDrawable,
              let cmd = queue.makeCommandBuffer() else { return }
        let bg = appearance.background
        pass.colorAttachments[0].clearColor = MTLClearColor(red: Double(bg.x), green: Double(bg.y), blue: Double(bg.z), alpha: Double(bg.w))
        pass.colorAttachments[0].loadAction = .clear
        guard let enc = cmd.makeRenderCommandEncoder(descriptor: pass) else { return }
        
        let size = viewSize
        let scale = SIMD2(2 * zoom / size.x, 2 * zoom / size.y)
        let backing = Float(view.window?.backingScaleFactor ?? 2)
        func uniforms(_ color: SIMD4<Float>, useValues: Bool = false, pointSize: Float = 1) -> Uniforms {
            Uniforms(center: center, scale: scale, vmin: appearance.rangeMin, vmax: appearance.rangeMax,
                     pointSize: pointSize * backing, useValues: useValues ? 1 : 0, color: color)
        }
        func setU(_ u: Uniforms) {
            var u = u
            enc.setVertexBytes(&u, length: MemoryLayout<Uniforms>.stride, index: 2)
            enc.setFragmentBytes(&u, length: MemoryLayout<Uniforms>.stride, index: 0)
        }
        
        // 1. Filled triangles
        if appearance.showFill, let tb = triBuffer, triCount > 0 {
            if let vb = valueBuffer {
                if lutColormap != appearance.colormap { makeLUT(appearance.colormap) }
                enc.setRenderPipelineState(fillPipeline)
                enc.setVertexBuffer(tb, offset: 0, index: 0)
                enc.setVertexBuffer(vb, offset: 0, index: 1)
                enc.setFragmentTexture(lutTexture, index: 0)
                setU(uniforms(SIMD4(1, 1, 1, 1), useValues: true))
            } else {
                enc.setRenderPipelineState(solidPipeline)
                enc.setVertexBuffer(tb, offset: 0, index: 0)
                setU(uniforms(appearance.solidFill))
            }
            enc.drawPrimitives(type: .triangle, vertexStart: 0, vertexCount: triCount)
        }
        // 2. Edges
        if appearance.showEdges, let eb = edgeBuffer, edgeCount > 0 {
            enc.setRenderPipelineState(solidPipeline)
            enc.setVertexBuffer(eb, offset: 0, index: 0)
            setU(uniforms(appearance.edgeColor))
            enc.drawPrimitives(type: .line, vertexStart: 0, vertexCount: edgeCount)
        }
        // 3. Outline (domain borders)
        if appearance.showOutline, let ob = outlineBuffer, outlineCount > 0 {
            enc.setRenderPipelineState(solidPipeline)
            enc.setVertexBuffer(ob, offset: 0, index: 0)
            setU(uniforms(appearance.outlineColor))
            enc.drawPrimitives(type: .line, vertexStart: 0, vertexCount: outlineCount)
        }
        // 4. Nodes
        if appearance.showNodes, let nb = nodeBuffer, nodeCount > 0 {
            enc.setRenderPipelineState(pointPipeline)
            enc.setVertexBuffer(nb, offset: 0, index: 0)
            enc.setVertexBuffer(nodeColorBuffer ?? nb, offset: 0, index: 1)
            var u = uniforms(appearance.nodeColor, useValues: nodeColorBuffer != nil, pointSize: appearance.nodeSizePx)
            enc.setVertexBytes(&u, length: MemoryLayout<Uniforms>.stride, index: 2)
            enc.setFragmentBytes(&u, length: MemoryLayout<Uniforms>.stride, index: 0)
            enc.drawPrimitives(type: .point, vertexStart: 0, vertexCount: nodeCount)
        }
        // 5. Quiver arrows (contrast halo pass followed by foreground arrow pass)
        if appearance.showArrows, let ab = quiverArrowBuffer, quiverArrowCount > 0 {
            enc.setRenderPipelineState(solidPipeline)
            if let hb = quiverHaloBuffer, quiverHaloCount > 0, quiverHaloColor.w > 0.001 {
                enc.setVertexBuffer(hb, offset: 0, index: 0)
                setU(uniforms(quiverHaloColor))
                enc.drawPrimitives(type: .triangle, vertexStart: 0, vertexCount: quiverHaloCount)
            }
            enc.setVertexBuffer(ab, offset: 0, index: 0)
            setU(uniforms(quiverArrowColor))
            enc.drawPrimitives(type: .triangle, vertexStart: 0, vertexCount: quiverArrowCount)
        }
        enc.endEncoding()
        cmd.present(drawable)
        cmd.commit()
    }
    
    private func makeLUT(_ cmap: MeshColormap) {
        let d = MTLTextureDescriptor.texture2DDescriptor(pixelFormat: .rgba8Unorm, width: 256, height: 1, mipmapped: false)
        d.usage = .shaderRead
        guard let tex = device.makeTexture(descriptor: d) else { return }
        let bytes = cmap.lut()
        bytes.withUnsafeBytes { tex.replace(region: MTLRegionMake2D(0, 0, 256, 1), mipmapLevel: 0, withBytes: $0.baseAddress!, bytesPerRow: 256 * 4) }
        lutTexture = tex
        lutColormap = cmap
    }
    
    // MARK: Shaders (compiled at runtime, so no Xcode build step is needed)
    
    static let shaderSource = """
    #include <metal_stdlib>
    using namespace metal;
    struct Uniforms { float2 center; float2 scale; float vmin; float vmax; float pointSize; int useValues; float4 color; };
    struct VOut { float4 pos [[position]]; float t; float4 color; float psize [[point_size]]; };
    
    vertex VOut fill_vertex(uint vid [[vertex_id]], const device float2* pos [[buffer(0)]],
                            const device float* val [[buffer(1)]], constant Uniforms& u [[buffer(2)]]) {
        VOut o;
        o.pos = float4((pos[vid] - u.center) * u.scale, 0.0, 1.0);
        float v = val[vid];
        o.t = isfinite(v) ? clamp((v - u.vmin) / max(u.vmax - u.vmin, 1e-20), 0.0, 1.0) : 0.0;
        o.color = u.color; o.psize = 1.0;
        return o;
    }
    fragment float4 fill_fragment(VOut in [[stage_in]], texture2d<float> lut [[texture(0)]]) {
        constexpr sampler s(filter::linear, address::clamp_to_edge);
        return lut.sample(s, float2(in.t, 0.5));
    }
    vertex VOut solid_vertex(uint vid [[vertex_id]], const device float2* pos [[buffer(0)]], constant Uniforms& u [[buffer(2)]]) {
        VOut o;
        o.pos = float4((pos[vid] - u.center) * u.scale, 0.0, 1.0);
        o.t = 0.0; o.color = u.color; o.psize = 1.0;
        return o;
    }
    fragment float4 solid_fragment(VOut in [[stage_in]]) { return in.color; }
    vertex VOut point_vertex(uint vid [[vertex_id]], const device float2* pos [[buffer(0)]],
                             const device float4* cols [[buffer(1)]], constant Uniforms& u [[buffer(2)]]) {
        VOut o;
        o.pos = float4((pos[vid] - u.center) * u.scale, 0.0, 1.0);
        o.t = 0.0; o.psize = u.pointSize;
        if (u.useValues != 0) { o.color = cols[vid]; } else { o.color = u.color; }
        return o;
    }
    fragment float4 point_fragment(VOut in [[stage_in]], float2 pc [[point_coord]]) {
        float d = length(pc - float2(0.5));
        if (d > 0.5) discard_fragment();
        float edge = smoothstep(0.38, 0.5, d);
        return float4(mix(in.color.rgb, float3(1.0), edge * 0.6), in.color.a);
    }
    """
}
