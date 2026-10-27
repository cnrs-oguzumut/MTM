import SwiftUI
import AppKit

// MARK: - "Blue Hour" Color Palette
// Swatches from darkest to lightest:
// 1. Abyss     #111622 / #0D1623 (RGB:  13,  22,  35)
// 2. Midnight  #1E293B / #1D293E (RGB:  29,  41,  62)
// 3. Slate     #31435F           (RGB:  49,  67,  95)
// 4. Steel     #586C89           (RGB:  88, 108, 137)
// 5. Mist      #97A8C0           (RGB: 151, 168, 192)
// 6. Ice       #D2E0ED           (RGB: 210, 224, 237)

public enum BlueHour {
    // Core 6-swatch palette
    public static let abyss    = Color(red: 13.0 / 255.0,  green: 22.0 / 255.0,  blue: 35.0 / 255.0)  // #0D1623
    public static let midnight = Color(red: 29.0 / 255.0,  green: 41.0 / 255.0,  blue: 62.0 / 255.0)  // #1D293E
    public static let slate    = Color(red: 49.0 / 255.0,  green: 67.0 / 255.0,  blue: 95.0 / 255.0)  // #31435F
    public static let steel    = Color(red: 88.0 / 255.0,  green: 108.0 / 255.0, blue: 137.0 / 255.0) // #586C89
    public static let mist     = Color(red: 151.0 / 255.0, green: 168.0 / 255.0, blue: 192.0 / 255.0) // #97A8C0
    public static let ice      = Color(red: 210.0 / 255.0, green: 224.0 / 255.0, blue: 237.0 / 255.0) // #D2E0ED

    // Semantic surface & UI tokens derived from the Blue Hour swatches
    public static let windowBackground  = abyss                                                       // #0D1623
    public static let canvasBackground  = Color(red: 10.0 / 255.0, green: 17.0 / 255.0, blue: 28.0 / 255.0) // #0A111C
    public static let sidebarBackground = Color(red: 19.0 / 255.0, green: 29.0 / 255.0, blue: 46.0 / 255.0) // #131D2E
    public static let controlSurface    = midnight                                                    // #1D293E
    public static let cardSurface       = Color(red: 34.0 / 255.0, green: 48.0 / 255.0, blue: 72.0 / 255.0) // #223048
    public static let elevatedSurface   = slate                                                       // #31435F
    public static let pillSurface       = slate.opacity(0.65)
    public static let border            = steel.opacity(0.45)
    public static let accent            = steel                                                       // #586C89
    public static let brightAccent      = mist                                                        // #97A8C0
    public static let textPrimary       = ice                                                         // #D2E0ED
    public static let textSecondary     = mist                                                        // #97A8C0

    // NSColor counterparts for AppKit window configuration
    public static let nsAbyss = NSColor(srgbRed: 13.0 / 255.0, green: 22.0 / 255.0, blue: 35.0 / 255.0, alpha: 1.0)
    public static let nsMidnight = NSColor(srgbRed: 29.0 / 255.0, green: 41.0 / 255.0, blue: 62.0 / 255.0, alpha: 1.0)
}

// MARK: - Custom GroupBox Style for Blue Hour Palette
public struct BlueHourGroupBoxStyle: GroupBoxStyle {
    public init() {}

    public func makeBody(configuration: Configuration) -> some View {
        VStack(alignment: .leading, spacing: 6) {
            configuration.label
                .font(.system(size: 12, weight: .semibold))
                .foregroundColor(BlueHour.ice)
                .tint(BlueHour.mist)

            configuration.content
                .padding(6)
                .frame(maxWidth: .infinity, alignment: .leading)
                .background(
                    RoundedRectangle(cornerRadius: 10, style: .continuous)
                        .fill(BlueHour.cardSurface.opacity(0.92))
                )
                .overlay(
                    RoundedRectangle(cornerRadius: 10, style: .continuous)
                        .strokeBorder(BlueHour.border, lineWidth: 1)
                )
        }
    }
}

