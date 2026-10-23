import Foundation
import AppKit

public enum ExportDataHelper {
    /// Presents a native macOS NSSavePanel allowing the user to select the destination file path.
    public static func saveStringToFile(
        content: String,
        defaultFileName: String,
        allowedExtensions: [String] = ["csv", "txt", "json", "dat"]
    ) {
        DispatchQueue.main.async {
            let savePanel = NSSavePanel()
            savePanel.canCreateDirectories = true
            savePanel.showsTagField = false
            savePanel.nameFieldStringValue = defaultFileName
            
            // Set default directory to Downloads or Documents if available
            if let downloadsURL = FileManager.default.urls(for: .downloadsDirectory, in: .userDomainMask).first {
                savePanel.directoryURL = downloadsURL
            }
            
            savePanel.begin { response in
                guard response == .OK, let targetURL = savePanel.url else { return }
                do {
                    try content.write(to: targetURL, atomically: true, encoding: .utf8)
                    // Optionally show notification or log
                    print("Successfully exported data to \(targetURL.path)")
                } catch {
                    let alert = NSAlert()
                    alert.messageText = "Export Failed"
                    alert.informativeText = "Could not save file: \(error.localizedDescription)"
                    alert.alertStyle = .warning
                    alert.runModal()
                }
            }
        }
    }
    
    /// Copies text content directly to macOS clipboard
    public static func copyToClipboard(text: String) {
        let pasteboard = NSPasteboard.general
        pasteboard.clearContents()
        pasteboard.setString(text, forType: .string)
    }
}
