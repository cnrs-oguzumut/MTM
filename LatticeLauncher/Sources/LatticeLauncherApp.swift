import SwiftUI
import AppKit

@main
struct LatticeLauncherApp: App {
    @StateObject private var config = SimulationConfig()
    @StateObject private var runner = ProcessRunner()
    
    var body: some Scene {
        WindowGroup {
            ContentView(config: config)
                .environmentObject(runner)
                .frame(minWidth: 1000, idealWidth: 1100, minHeight: 700, idealHeight: 750)
        }
        .windowStyle(.titleBar)
        .windowToolbarStyle(.unified)
    }
}
