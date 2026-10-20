import SwiftUI

// MARK: - Reliable Numeric Input Helpers
public struct DoubleField: View {
    public let title: String
    @Binding public var value: Double
    public var width: CGFloat
    public var minValue: Double?
    public var maxValue: Double?
    @State private var text: String = ""
    @FocusState private var isFocused: Bool
    
    public init(title: String, value: Binding<Double>, width: CGFloat = 80, minValue: Double? = nil, maxValue: Double? = nil) {
        self.title = title
        self._value = value
        self.width = width
        self.minValue = minValue
        self.maxValue = maxValue
    }
    
    private static func formatValue(_ v: Double) -> String {
        if v.isNaN || v.isInfinite { return "" }
        if v.truncatingRemainder(dividingBy: 1) == 0 && abs(v) < 1e7 {
            return String(format: "%.0f", locale: Locale(identifier: "en_US_POSIX"), v)
        }
        let s = String(format: "%.6g", locale: Locale(identifier: "en_US_POSIX"), v)
        return s.replacingOccurrences(of: ",", with: ".")
    }
    
    public var body: some View {
        TextField(title, text: $text)
            .textFieldStyle(.roundedBorder)
            .frame(width: width)
            .focused($isFocused)
            .onAppear {
                text = Self.formatValue(value)
            }
            .onChange(of: value) { _, newVal in
                let currentParsed = Double(text.replacingOccurrences(of: ",", with: "."))
                if currentParsed == nil || abs((currentParsed ?? 0) - newVal) > 1e-9 {
                    text = Self.formatValue(newVal)
                }
            }
            .onChange(of: text) { _, newText in
                // Automatically and immediately convert any comma ',' to dot '.'
                let clean = newText.replacingOccurrences(of: ",", with: ".")
                if clean != newText {
                    text = clean
                }
                if let parsed = Double(clean) {
                    if let minV = minValue, parsed < minV {
                        // User is still typing intermediate characters (e.g. 0.00 before typing 5)
                        return
                    }
                    if let maxV = maxValue, parsed > maxV {
                        return
                    }
                    value = parsed
                }
            }
            .onChange(of: isFocused) { _, focused in
                if !focused {
                    let clean = text.replacingOccurrences(of: ",", with: ".")
                    if let parsed = Double(clean) {
                        var clamped = parsed
                        if let minV = minValue { clamped = max(minV, clamped) }
                        if let maxV = maxValue { clamped = min(maxV, clamped) }
                        value = clamped
                        text = Self.formatValue(clamped)
                    } else {
                        text = Self.formatValue(value)
                    }
                }
            }
            .onSubmit {
                let clean = text.replacingOccurrences(of: ",", with: ".")
                if let parsed = Double(clean) {
                    var clamped = parsed
                    if let minV = minValue { clamped = max(minV, clamped) }
                    if let maxV = maxValue { clamped = min(maxV, clamped) }
                    value = clamped
                    text = Self.formatValue(clamped)
                } else {
                    text = Self.formatValue(value)
                }
            }
    }
}

public struct IntField: View {
    public let title: String
    @Binding public var value: Int
    public var width: CGFloat
    @State private var text: String = ""
    @FocusState private var isFocused: Bool
    
    public init(title: String, value: Binding<Int>, width: CGFloat = 70) {
        self.title = title
        self._value = value
        self.width = width
    }
    
    public var body: some View {
        TextField(title, text: $text)
            .textFieldStyle(.roundedBorder)
            .frame(width: width)
            .focused($isFocused)
            .onAppear {
                text = String(value)
            }
            .onChange(of: value) { _, newVal in
                if !isFocused {
                    text = String(newVal)
                }
            }
            .onChange(of: text) { _, newText in
                let filtered = newText.filter { "-0123456789".contains($0) }
                if filtered != newText {
                    text = filtered
                }
                if let parsed = Int(filtered) {
                    value = parsed
                }
            }
            .onChange(of: isFocused) { _, focused in
                if !focused {
                    if let parsed = Int(text) {
                        value = parsed
                        text = String(parsed)
                    } else {
                        text = String(value)
                    }
                }
            }
            .onSubmit {
                if let parsed = Int(text) {
                    value = parsed
                    text = String(parsed)
                } else {
                    text = String(value)
                }
            }
    }
}
