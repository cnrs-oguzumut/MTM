import SwiftUI

/// Shared rule for drawing individual nodes / lattice sites in large crystals.
public enum LargeMeshPolicy {
    /// Above this many nodes, nodes start hidden and turning them on asks for confirmation.
    public static let nodeLimit = 5_000
    /// Above this many nodes, node id labels are never drawn.
    public static let labelLimit = 2_000
    /// Above this many triangles, the FEM preview switches from Canvas to the GPU renderer.
    public static let gpuTriangleThreshold = 3_000
    
    public static func isLarge(_ nodeCount: Int) -> Bool { nodeCount > nodeLimit }
    
    public static func formatted(_ n: Int) -> String {
        let f = NumberFormatter()
        f.numberStyle = .decimal
        return f.string(from: NSNumber(value: n)) ?? "\(n)"
    }
}

/// A "Nodes" toggle that asks before turning nodes on in a large crystal.
public struct GuardedNodesToggle: View {
    let title: String
    @Binding var isOn: Bool
    let nodeCount: Int
    @State private var askConfirmation = false
    
    public init(_ title: String = "Nodes", isOn: Binding<Bool>, nodeCount: Int) {
        self.title = title
        self._isOn = isOn
        self.nodeCount = nodeCount
    }
    
    public var body: some View {
        Toggle(title, isOn: Binding(
            get: { isOn },
            set: { newValue in
                if newValue && LargeMeshPolicy.isLarge(nodeCount) {
                    askConfirmation = true
                } else {
                    isOn = newValue
                }
            }
        ))
        .largeNodeConfirmation(isPresented: $askConfirmation, nodeCount: nodeCount) { isOn = true }
    }
}

/// Small pill shown on top of a viewer when nodes were hidden automatically.
public struct NodesHiddenNote: View {
    let nodeCount: Int
    let noun: String
    let onShowAnyway: () -> Void
    @State private var askConfirmation = false
    
    public init(nodeCount: Int, noun: String = "nodes", onShowAnyway: @escaping () -> Void) {
        self.nodeCount = nodeCount
        self.noun = noun
        self.onShowAnyway = onShowAnyway
    }
    
    public var body: some View {
        HStack(spacing: 6) {
            Image(systemName: "eye.slash")
            Text("\(LargeMeshPolicy.formatted(nodeCount)) \(noun) hidden to keep the viewer fast")
            Button("Show anyway") { askConfirmation = true }
                .buttonStyle(.link)
        }
        .font(.system(size: 10.5))
        .foregroundColor(BlueHour.mist)
        .padding(.horizontal, 9)
        .padding(.vertical, 5)
        .background(.ultraThinMaterial, in: Capsule())
        .largeNodeConfirmation(isPresented: $askConfirmation, nodeCount: nodeCount, noun: noun, onConfirm: onShowAnyway)
    }
}

extension View {
    /// The standard "Show N nodes?" warning.
    func largeNodeConfirmation(isPresented: Binding<Bool>, nodeCount: Int, noun: String = "nodes", onConfirm: @escaping () -> Void) -> some View {
        alert("Show \(LargeMeshPolicy.formatted(nodeCount)) \(noun)?", isPresented: isPresented) {
            Button("Cancel", role: .cancel) {}
            Button("Show Anyway") { onConfirm() }
        } message: {
            Text("Drawing this many individual \(noun) can make the viewer slow and may briefly freeze the app. Triangles and colours are not affected."
                 + (nodeCount > LargeMeshPolicy.labelLimit ? " Id labels stay hidden at this size." : ""))
        }
    }
}
