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

    // High-contrast interactive state tokens (eliminates low-contrast Dark Aqua grey pills)
    public static let activeAccent      = Color(red: 0.13, green: 0.76, blue: 0.68)                   // #21C2AD vibrant teal-cyan
    public static let tabTrack          = Color(red: 9.0 / 255.0, green: 15.0 / 255.0, blue: 25.0 / 255.0) // #090F19 recessed well
    public static let activeTabFill     = Color(red: 0.10, green: 0.50, blue: 0.45)                   // #1A8073 illuminated teal fill
    public static let activeTabBorder   = Color(red: 0.37, green: 0.92, blue: 0.83)                   // #5EEAD4 crisp luminous border

    // Semantic metric tab fills & borders (>= 5:1 white text contrast)
    public static let stressTabFill     = Color(red: 0.13, green: 0.39, blue: 0.82)                   // #2163D1 cobalt-sky
    public static let stressTabBorder   = Color(red: 0.45, green: 0.73, blue: 1.00)                   // #73BAFF
    public static let energyTabFill     = Color(red: 0.43, green: 0.22, blue: 0.78)                   // #6E38C7 royal violet
    public static let energyTabBorder   = Color(red: 0.74, green: 0.58, blue: 1.00)                   // #BD94FF
    public static let eigenTabFill      = Color(red: 0.06, green: 0.48, blue: 0.45)                   // #0F7A73 deep teal
    public static let eigenTabBorder    = Color(red: 0.35, green: 0.92, blue: 0.84)                   // #59EBD6
    public static let gbTabFill         = Color(red: 0.09, green: 0.50, blue: 0.27)                   // #178045 rich emerald
    public static let gbTabBorder       = Color(red: 0.43, green: 0.91, blue: 0.61)                   // #6EE89C

    // NSColor counterparts for AppKit window configuration
    public static let nsAbyss = NSColor(srgbRed: 13.0 / 255.0, green: 22.0 / 255.0, blue: 35.0 / 255.0, alpha: 1.0)
    public static let nsMidnight = NSColor(srgbRed: 29.0 / 255.0, green: 41.0 / 255.0, blue: 62.0 / 255.0, alpha: 1.0)
}

// MARK: - High-Contrast Segmented Picker & Toggle Chip Controls

public enum BlueHourControlScale {
    case micro      // Tight inline cells (e.g. Free | Fix | Drive)
    case compact    // Toolbars & sub-headers
    case regular    // Primary sidebar & chart tabs
    case prominent  // Global workspace switcher

    var fontSize: CGFloat {
        switch self {
        case .micro:     return 9.5
        case .compact:   return 10.5
        case .regular:   return 11.5
        case .prominent: return 12.0
        }
    }

    var iconSize: CGFloat {
        switch self {
        case .micro:     return 8.5
        case .compact:   return 9.5
        case .regular:   return 10.5
        case .prominent: return 11.5
        }
    }

    var hPadding: CGFloat {
        switch self {
        case .micro:     return 5
        case .compact:   return 8
        case .regular:   return 10
        case .prominent: return 13
        }
    }

    var vPadding: CGFloat {
        switch self {
        case .micro:     return 3
        case .compact:   return 4
        case .regular:   return 5
        case .prominent: return 5.5
        }
    }

    var cornerRadius: CGFloat {
        switch self {
        case .micro:     return 5
        case .compact:   return 6
        case .regular:   return 7
        case .prominent: return 8
        }
    }
}

public struct BlueHourTabItem<Value: Hashable>: Identifiable {
    public var id: Value { value }
    public let value: Value
    public let title: String
    public var icon: String?
    public var badge: String?
    public var hasData: Bool
    public var liveDotColor: Color?
    public var activeFill: Color
    public var activeBorder: Color
    public var helpText: String?

    public init(
        _ value: Value,
        title: String,
        icon: String? = nil,
        badge: String? = nil,
        hasData: Bool = true,
        liveDotColor: Color? = nil,
        activeFill: Color = BlueHour.activeTabFill,
        activeBorder: Color = BlueHour.activeTabBorder,
        helpText: String? = nil
    ) {
        self.value = value
        self.title = title
        self.icon = icon
        self.badge = badge
        self.hasData = hasData
        self.liveDotColor = liveDotColor
        self.activeFill = activeFill
        self.activeBorder = activeBorder
        self.helpText = helpText
    }
}

public struct BlueHourSegmentedPicker<Value: Hashable>: View {
    @Binding public var selection: Value
    public let items: [BlueHourTabItem<Value>]
    public var scale: BlueHourControlScale
    public var expandEqually: Bool

    @State private var hoveredValue: Value? = nil

    public init(
        selection: Binding<Value>,
        items: [BlueHourTabItem<Value>],
        scale: BlueHourControlScale = .regular,
        expandEqually: Bool = false
    ) {
        self._selection = selection
        self.items = items
        self.scale = scale
        self.expandEqually = expandEqually
    }

    public var body: some View {
        HStack(spacing: 2) {
            ForEach(items) { item in
                let isSelected = (selection == item.value)
                let isHovered = (hoveredValue == item.value)

                Button {
                    withAnimation(.easeOut(duration: 0.13)) {
                        selection = item.value
                    }
                } label: {
                    HStack(spacing: 4.5) {
                        if let dotColor = item.liveDotColor {
                            Circle()
                                .fill(isSelected ? Color.white : dotColor)
                                .frame(width: 6, height: 6)
                                .shadow(color: dotColor.opacity(0.65), radius: 2)
                        }
                        if let icon = item.icon {
                            Image(systemName: icon)
                                .font(.system(size: scale.iconSize, weight: isSelected ? .bold : .medium))
                        }
                        Text(item.title)
                            .font(.system(size: scale.fontSize, weight: isSelected ? .bold : .medium))
                            .lineLimit(1)
                            .fixedSize(horizontal: !expandEqually, vertical: true)
                        if let badge = item.badge, !badge.isEmpty {
                            Text(badge)
                                .font(.system(size: max(8.5, scale.fontSize - 2.0), weight: .bold, design: .monospaced))
                                .padding(.horizontal, 4.5)
                                .padding(.vertical, 1)
                                .background(
                                    Capsule()
                                        .fill(isSelected ? Color.black.opacity(0.28) : BlueHour.slate.opacity(0.65))
                                )
                                .foregroundColor(isSelected ? .white : BlueHour.ice)
                        }
                    }
                    .foregroundColor(
                        isSelected
                        ? .white
                        : (item.hasData
                           ? (isHovered ? BlueHour.ice : BlueHour.mist)
                           : (isHovered ? BlueHour.mist.opacity(0.78) : BlueHour.mist.opacity(0.45)))
                    )
                    .padding(.horizontal, scale.hPadding)
                    .padding(.vertical, scale.vPadding)
                    .frame(maxWidth: expandEqually ? .infinity : nil)
                    .contentShape(RoundedRectangle(cornerRadius: scale.cornerRadius, style: .continuous))
                    .background(
                        RoundedRectangle(cornerRadius: scale.cornerRadius, style: .continuous)
                            .fill(
                                isSelected
                                ? item.activeFill
                                : (isHovered ? BlueHour.slate.opacity(0.38) : Color.clear)
                            )
                    )
                    .overlay(
                        RoundedRectangle(cornerRadius: scale.cornerRadius, style: .continuous)
                            .strokeBorder(
                                isSelected
                                ? item.activeBorder
                                : (isHovered ? BlueHour.steel.opacity(0.45) : Color.clear),
                                lineWidth: isSelected ? 1.25 : 0.8
                            )
                    )
                    .shadow(color: isSelected ? item.activeFill.opacity(0.38) : .clear, radius: 3, y: 1)
                }
                .buttonStyle(.plain)
                .onHover { hovering in
                    if hovering {
                        hoveredValue = item.value
                    } else if hoveredValue == item.value {
                        hoveredValue = nil
                    }
                }
                .help(item.helpText ?? item.title)
            }
        }
        .padding(2.5)
        .background(
            RoundedRectangle(cornerRadius: scale.cornerRadius + 2.5, style: .continuous)
                .fill(BlueHour.tabTrack)
        )
        .overlay(
            RoundedRectangle(cornerRadius: scale.cornerRadius + 2.5, style: .continuous)
                .strokeBorder(BlueHour.steel.opacity(0.42), lineWidth: 1)
        )
    }
}

public struct BlueHourToggleChip: View {
    public let title: String
    public var icon: String?
    public var badge: String?
    @Binding public var isOn: Bool
    public var showsCheckmark: Bool
    public var activeFill: Color
    public var activeBorder: Color
    public var scale: BlueHourControlScale
    public var helpText: String?

    @State private var isHovered: Bool = false

    public init(
        _ title: String,
        icon: String? = nil,
        badge: String? = nil,
        isOn: Binding<Bool>,
        showsCheckmark: Bool = true,
        activeFill: Color = BlueHour.activeTabFill,
        activeBorder: Color = BlueHour.activeTabBorder,
        scale: BlueHourControlScale = .compact,
        helpText: String? = nil
    ) {
        self.title = title
        self.icon = icon
        self.badge = badge
        self._isOn = isOn
        self.showsCheckmark = showsCheckmark
        self.activeFill = activeFill
        self.activeBorder = activeBorder
        self.scale = scale
        self.helpText = helpText
    }

    public var body: some View {
        Button {
            withAnimation(.easeOut(duration: 0.12)) {
                isOn.toggle()
            }
        } label: {
            HStack(spacing: 4.5) {
                if showsCheckmark {
                    Image(systemName: isOn ? "checkmark.circle.fill" : "circle")
                        .font(.system(size: scale.iconSize + 0.5, weight: .semibold))
                        .foregroundColor(isOn ? activeBorder : BlueHour.steel)
                }
                if let icon = icon {
                    Image(systemName: icon)
                        .font(.system(size: scale.iconSize, weight: isOn ? .bold : .medium))
                }
                Text(title)
                    .font(.system(size: scale.fontSize, weight: isOn ? .bold : .medium))
                    .lineLimit(1)
                    .fixedSize(horizontal: true, vertical: true)
                if let badge = badge, !badge.isEmpty {
                    Text(badge)
                        .font(.system(size: max(8.5, scale.fontSize - 2.0), weight: .bold, design: .monospaced))
                        .padding(.horizontal, 4)
                        .padding(.vertical, 0.5)
                        .background(
                            Capsule()
                                .fill(isOn ? Color.black.opacity(0.28) : BlueHour.slate.opacity(0.55))
                        )
                        .foregroundColor(isOn ? .white : BlueHour.mist)
                }
            }
            .foregroundColor(isOn ? .white : (isHovered ? BlueHour.ice : BlueHour.mist))
            .padding(.horizontal, scale.hPadding)
            .padding(.vertical, scale.vPadding)
            .background(
                RoundedRectangle(cornerRadius: scale.cornerRadius + 1, style: .continuous)
                    .fill(isOn ? activeFill.opacity(0.88) : (isHovered ? BlueHour.slate.opacity(0.45) : BlueHour.tabTrack))
            )
            .overlay(
                RoundedRectangle(cornerRadius: scale.cornerRadius + 1, style: .continuous)
                    .strokeBorder(
                        isOn ? activeBorder : (isHovered ? BlueHour.steel.opacity(0.65) : BlueHour.steel.opacity(0.38)),
                        lineWidth: isOn ? 1.2 : 1.0
                    )
            )
        }
        .buttonStyle(.plain)
        .onHover { isHovered = $0 }
        .help(helpText ?? title)
    }
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

