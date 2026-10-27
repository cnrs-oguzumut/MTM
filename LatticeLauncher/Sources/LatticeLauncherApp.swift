import SwiftUI
import AppKit

@main
struct LatticeLauncherApp: App {
    @StateObject private var config = SimulationConfig()
    @StateObject private var runner = ProcessRunner()
    
    init() {
        // Apply dark appearance matching the Blue Hour palette
        NSApplication.shared.appearance = NSAppearance(named: .darkAqua)
        
        // Set Dock icon dynamically
        if let iconPath = Bundle.main.path(forResource: "AppIcon", ofType: "icns"),
           let iconImg = NSImage(contentsOfFile: iconPath) {
            NSApplication.shared.applicationIconImage = iconImg
        }
    }
    
    var body: some Scene {
        WindowGroup {
            ContentView(config: config)
                .environmentObject(runner)
                .environment(\.locale, Locale(identifier: "en_US_POSIX"))
                .preferredColorScheme(.dark)
                .tint(BlueHour.steel)
                .groupBoxStyle(BlueHourGroupBoxStyle())
                .background(BlueHour.windowBackground)
                .frame(minWidth: 1120, idealWidth: 1280, minHeight: 740, idealHeight: 840)
        }
        .windowStyle(.titleBar)
        .windowToolbarStyle(.unified)
    }
}
