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
        .commands {
            // 1. Experiment Menu
            CommandMenu("Experiment") {
                Button("Run Simulation") {
                    NotificationCenter.default.post(name: .runSimulationCommand, object: nil)
                }
                .keyboardShortcut("r", modifiers: [.command])
                .disabled(runner.status.isRunning)
                
                Button("Stop Simulation") {
                    NotificationCenter.default.post(name: .stopSimulationCommand, object: nil)
                }
                .keyboardShortcut(".", modifiers: [.command])
                .disabled(!runner.status.isRunning)
                
                Divider()
                
                Menu("Load Sample Experiment") {
                    ForEach(SampleExperimentPreset.allCases) { preset in
                        Button(preset.title) {
                            NotificationCenter.default.post(name: .applySampleExperiment, object: preset.rawValue)
                        }
                    }
                }
                
                Divider()
                
                Button("Open Output Folder in Finder") {
                    runner.openOutputFolder()
                }
                .keyboardShortcut("o", modifiers: [.command, .shift])
            }
            
            // 2. View Menu Shortcuts
            CommandGroup(after: .toolbar) {
                Divider()
                Button("Specimen Blueprint") {
                    NotificationCenter.default.post(name: .switchMainDisplayTab, object: "blueprint")
                }
                .keyboardShortcut("1", modifiers: [.command])
                
                Button("Specimen FEM Mesh") {
                    NotificationCenter.default.post(name: .switchMainDisplayTab, object: "femMesh")
                }
                .keyboardShortcut("2", modifiers: [.command])
                
                Button("Results & Plots") {
                    NotificationCenter.default.post(name: .switchMainDisplayTab, object: "results")
                }
                .keyboardShortcut("3", modifiers: [.command])
                
                Button("Split View (Specimen + Results)") {
                    NotificationCenter.default.post(name: .switchMainDisplayTab, object: "split")
                }
                .keyboardShortcut("4", modifiers: [.command])
                
                Divider()
                
                Button("Show / Hide Console") {
                    NotificationCenter.default.post(name: .toggleConsoleDrawer, object: nil)
                }
                .keyboardShortcut("j", modifiers: [.command])
                
                Button("Toggle Guided / Expert Hints") {
                    NotificationCenter.default.post(name: .toggleComplexityMode, object: nil)
                }
                .keyboardShortcut("e", modifiers: [.command, .option])
            }
            
            // 3. Help Menu (Mac App Store compliant Help & Onboarding)
            CommandGroup(replacing: .help) {
                Button("Welcome & Quick-Start Tour…") {
                    NotificationCenter.default.post(name: .showWelcomeTour, object: nil)
                }
                
                Button("Lattice Launcher Documentation & Physics Guide…") {
                    NotificationCenter.default.post(name: .showDocumentationCenter, object: DocumentationTopic.quickStart.rawValue)
                }
                .keyboardShortcut("?", modifiers: [.command])
                
                Divider()
                
                Button("Guide: How to Install Dislocations (2M Periodic & Volterra)…") {
                    NotificationCenter.default.post(name: .showDocumentationCenter, object: DocumentationTopic.installingDislocations.rawValue)
                }
                
                Button("Guide: Bicrystals, Voronoi & CSL Rationalizer…") {
                    NotificationCenter.default.post(name: .showDocumentationCenter, object: DocumentationTopic.grainBoundariesCSL.rawValue)
                }
                
                Button("Guide: Two-Stage Boundary Conditions & Nanoindenter…") {
                    NotificationCenter.default.post(name: .showDocumentationCenter, object: DocumentationTopic.boundaryConditionsIndenter.rawValue)
                }
                
                Button("Glossary of Physical Symbols & Parameters…") {
                    NotificationCenter.default.post(name: .showDocumentationCenter, object: DocumentationTopic.parameterGlossary.rawValue)
                }
            }
        }
    }
}
