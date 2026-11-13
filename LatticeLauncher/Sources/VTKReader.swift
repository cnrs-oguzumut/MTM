import Foundation
import simd

// MARK: - Data Model

/// One colourable quantity stored in a VTK file (a scalar, or one component of a vector / tensor).
public struct VTKField: Identifiable, Hashable {
    public enum Location: String { case cell, point }
    
    public let name: String
    public let location: Location
    public let values: [Float]
    
    public var id: String { "\(location.rawValue):\(name)" }
    
    /// Human-friendly label, e.g. "ElementVonMisesStress" → "Von Mises Stress (element)".
    public var displayName: String {
        if name.hasPrefix("mode_") {
            let parts = name.split(separator: "_")
            if parts.count >= 3, let modeNum = Int(parts[1]) {
                let suffix = parts.dropFirst(2).joined(separator: "_")
                switch suffix {
                case "amplitude", "magnitude":
                    return "Mode \(modeNum) Amplitude |u| (node)"
                case "x":
                    return "Mode \(modeNum) Displacement uₓ (node)"
                case "y":
                    return "Mode \(modeNum) Displacement uᵧ (node)"
                default:
                    break
                }
            }
        }
        var n = name
        for prefix in ["Element", "Nodal", "Node"] where n.hasPrefix(prefix) && n.count > prefix.count {
            n = String(n.dropFirst(prefix.count))
            break
        }
        // Split CamelCase into words, keep suffixes like _xx untouched
        var out = ""
        var prevLower = false
        for ch in n {
            if ch.isUppercase && prevLower { out.append(" ") }
            out.append(ch)
            prevLower = ch.isLowercase || ch.isNumber
        }
        out = out.replacingOccurrences(of: "_", with: " ")
        return "\(out) (\(location == .cell ? "element" : "node"))"
    }
    
    public static func == (a: VTKField, b: VTKField) -> Bool { a.id == b.id }
    public func hash(into h: inout Hasher) { h.combine(id) }
}

/// A 2D vector field stored in a VTK `VECTORS` block (e.g. soft modes `mode_1`, `mode_2`, or displacements).
public struct VTKVectorField: Identifiable, Hashable {
    public let name: String
    public let location: VTKField.Location
    public let vectors: [SIMD2<Float>]
    
    public var id: String { "\(location.rawValue):\(name)" }
    
    public var displayName: String {
        if name.hasPrefix("mode_") {
            let parts = name.split(separator: "_")
            if parts.count >= 2, let modeNum = Int(parts[1]) {
                return "Mode \(modeNum) Vectors u (\(location == .cell ? "element" : "node"))"
            }
        }
        var n = name
        for prefix in ["Element", "Nodal", "Node"] where n.hasPrefix(prefix) && n.count > prefix.count {
            n = String(n.dropFirst(prefix.count))
            break
        }
        var out = ""
        var prevLower = false
        for ch in n {
            if ch.isUppercase && prevLower { out.append(" ") }
            out.append(ch)
            prevLower = ch.isLowercase || ch.isNumber
        }
        out = out.replacingOccurrences(of: "_", with: " ")
        return "\(out) (\(location == .cell ? "element" : "node"))"
    }
    
    public static func == (a: VTKVectorField, b: VTKVectorField) -> Bool { a.id == b.id }
    public func hash(into h: inout Hasher) { h.combine(id) }
}

public struct VTKDataset {
    public var title: String = ""
    public var points: [SIMD2<Float>] = []
    /// Three node indices per triangle (non-triangle cells are skipped).
    public var triangles: [UInt32] = []
    public var fields: [VTKField] = []
    public var vectorFields: [VTKVectorField] = []
    /// `FIELD FieldData` → `LoadParameter` (α), when present.
    public var loadParameter: Double?
    /// Eigenvalues parsed from soft-mode VTK title header (`... lambda <λ1> <λ2> ...`), when present.
    public var modeEigenvalues: [Double] = []
    
    public var nodeCount: Int { points.count }
    public var triangleCount: Int { triangles.count / 3 }
    
    public var bounds: (min: SIMD2<Float>, max: SIMD2<Float>) {
        guard var lo = points.first else { return (.zero, SIMD2(1, 1)) }
        var hi = lo
        for p in points {
            lo = simd_min(lo, p)
            hi = simd_max(hi, p)
        }
        return (lo, hi)
    }
}

public enum VTKReaderError: LocalizedError {
    case notVTK
    case unsupportedDataset(String)
    case truncated(String)
    case badNumber(String)
    
    public var errorDescription: String? {
        switch self {
        case .notVTK: return "Not a legacy VTK file (missing '# vtk' header)."
        case .unsupportedDataset(let s): return "Unsupported VTK dataset type '\(s)'. Only UNSTRUCTURED_GRID is supported."
        case .truncated(let s): return "The file ends unexpectedly while reading \(s)."
        case .badNumber(let s): return "Could not read a number while reading \(s)."
        }
    }
}

// MARK: - Reader

/// Reader for legacy `.vtk` files written by the MTM solver (ASCII or BINARY big-endian UNSTRUCTURED_GRID).
/// Call from a background task – it does not touch the UI.
public enum VTKReader {
    public static func load(url: URL) throws -> VTKDataset {
        var data = try Data(contentsOf: url, options: .mappedIfSafe)
        data.append(0) // sentinel so strtod never runs off the end
        return try data.withUnsafeBytes { raw -> VTKDataset in
            var parser = Parser(raw: raw)
            return try parser.parse()
        }
    }
}

private struct Parser {
    let base: UnsafeRawPointer
    let count: Int          // excluding the sentinel
    var pos: Int = 0
    var binary = false
    
    init(raw: UnsafeRawBufferPointer) {
        base = raw.baseAddress!
        count = raw.count - 1
    }
    
    @inline(__always) func byte(_ i: Int) -> UInt8 { base.load(fromByteOffset: i, as: UInt8.self) }
    @inline(__always) static func isSpace(_ b: UInt8) -> Bool { b == 32 || b == 10 || b == 13 || b == 9 }
    
    /// Reads up to (and consumes) the next '\n'. Returns nil at end of file.
    mutating func readLine() -> String? {
        guard pos < count else { return nil }
        let start = pos
        while pos < count && byte(pos) != 10 { pos += 1 }
        var end = pos
        if pos < count { pos += 1 } // consume \n
        if end > start && byte(end - 1) == 13 { end -= 1 }
        let buf = UnsafeRawBufferPointer(start: base + start, count: end - start)
        return String(decoding: buf, as: UTF8.self)
    }
    
    /// Next line that contains something other than whitespace, split into tokens.
    mutating func nextTokens() -> [String]? {
        while let line = readLine() {
            let toks = line.split(whereSeparator: { $0 == " " || $0 == "\t" }).map(String.init)
            if !toks.isEmpty { return toks }
        }
        return nil
    }
    
    // MARK: Numeric blocks
    
    static func byteSize(_ type: String) -> Int {
        switch type.lowercased() {
        case "double", "long", "unsigned_long", "vtktypeint64", "vtktypeuint64": return 8
        case "short", "unsigned_short": return 2
        case "char", "unsigned_char", "bit": return 1
        default: return 4 // float, int, unsigned_int
        }
    }
    
    mutating func readFloats(_ n: Int, type: String, what: String) throws -> [Float] {
        var out = [Float](repeating: 0, count: n)
        guard n > 0 else { return out }
        if binary {
            let t = type.lowercased()
            let sz = Parser.byteSize(t)
            guard pos + n * sz <= count else { throw VTKReaderError.truncated(what) }
            let p = base + pos
            switch t {
            case "float":
                for i in 0..<n { out[i] = Float(bitPattern: UInt32(bigEndian: p.loadUnaligned(fromByteOffset: i * 4, as: UInt32.self))) }
            case "double":
                for i in 0..<n { out[i] = Float(Double(bitPattern: UInt64(bigEndian: p.loadUnaligned(fromByteOffset: i * 8, as: UInt64.self)))) }
            case "int", "vtkidtype":
                for i in 0..<n { out[i] = Float(Int32(bigEndian: p.loadUnaligned(fromByteOffset: i * 4, as: Int32.self))) }
            case "unsigned_int":
                for i in 0..<n { out[i] = Float(UInt32(bigEndian: p.loadUnaligned(fromByteOffset: i * 4, as: UInt32.self))) }
            case "long", "vtktypeint64":
                for i in 0..<n { out[i] = Float(Int64(bigEndian: p.loadUnaligned(fromByteOffset: i * 8, as: Int64.self))) }
            case "short":
                for i in 0..<n { out[i] = Float(Int16(bigEndian: p.loadUnaligned(fromByteOffset: i * 2, as: Int16.self))) }
            case "unsigned_char", "char", "bit":
                for i in 0..<n { out[i] = Float(p.load(fromByteOffset: i, as: UInt8.self)) }
            default:
                for i in 0..<n { out[i] = Float(bitPattern: UInt32(bigEndian: p.loadUnaligned(fromByteOffset: i * 4, as: UInt32.self))) }
            }
            pos += n * sz
        } else {
            for i in 0..<n {
                out[i] = Float(try readASCIINumber(what: what))
            }
        }
        return out
    }
    
    mutating func readInts(_ n: Int, what: String) throws -> [Int32] {
        var out = [Int32](repeating: 0, count: n)
        guard n > 0 else { return out }
        if binary {
            guard pos + n * 4 <= count else { throw VTKReaderError.truncated(what) }
            let p = base + pos
            for i in 0..<n { out[i] = Int32(bigEndian: p.loadUnaligned(fromByteOffset: i * 4, as: Int32.self)) }
            pos += n * 4
        } else {
            for i in 0..<n { out[i] = Int32(try readASCIINumber(what: what)) }
        }
        return out
    }
    
    mutating func readASCIINumber(what: String) throws -> Double {
        while pos < count && Parser.isSpace(byte(pos)) { pos += 1 }
        guard pos < count else { throw VTKReaderError.truncated(what) }
        let start = (base + pos).assumingMemoryBound(to: CChar.self)
        var end: UnsafeMutablePointer<CChar>? = nil
        let v = strtod(start, &end)
        guard let e = end, e != UnsafeMutablePointer(mutating: start) else { throw VTKReaderError.badNumber(what) }
        pos += e - UnsafeMutablePointer(mutating: start)
        return v
    }
    
    // MARK: Main parse
    
    mutating func parse() throws -> VTKDataset {
        guard let header = readLine(), header.lowercased().hasPrefix("# vtk") else { throw VTKReaderError.notVTK }
        var ds = VTKDataset()
        ds.title = readLine() ?? ""
        let fmt = (readLine() ?? "").trimmingCharacters(in: .whitespaces).uppercased()
        binary = (fmt == "BINARY")
        
        var keepCell: [Bool] = []          // which cells became triangles
        var allTriangles = true
        var location: VTKField.Location = .point
        var tupleCount = 0
        
        func addField(_ name: String, _ vals: [Float], _ loc: VTKField.Location) {
            var v = vals
            if loc == .cell && !allTriangles && v.count == keepCell.count {
                v = zip(vals, keepCell).compactMap { $1 ? $0 : nil }
            }
            ds.fields.append(VTKField(name: name, location: loc, values: v))
        }
        
        parseLoop: while let toks = nextTokens() {
            let key = toks[0].uppercased()
            switch key {
            case "DATASET":
                let kind = toks.count > 1 ? toks[1].uppercased() : ""
                guard kind == "UNSTRUCTURED_GRID" else { throw VTKReaderError.unsupportedDataset(kind) }
                
            case "POINTS":
                guard toks.count >= 3, let n = Int(toks[1]) else { break parseLoop }
                let v = try readFloats(3 * n, type: toks[2], what: "POINTS")
                ds.points = (0..<n).map { SIMD2(v[3 * $0], v[3 * $0 + 1]) }
                
            case "CELLS":
                guard toks.count >= 3, let nCells = Int(toks[1]), let size = Int(toks[2]) else { break parseLoop }
                let ints = try readInts(size, what: "CELLS")
                keepCell = [Bool](repeating: false, count: nCells)
                var tris: [UInt32] = []
                tris.reserveCapacity(nCells * 3)
                var i = 0
                for c in 0..<nCells {
                    guard i < ints.count else { break }
                    let k = Int(ints[i])
                    if k == 3 && i + 3 < ints.count {
                        tris.append(UInt32(max(0, ints[i + 1])))
                        tris.append(UInt32(max(0, ints[i + 2])))
                        tris.append(UInt32(max(0, ints[i + 3])))
                        keepCell[c] = true
                    } else {
                        allTriangles = false
                    }
                    i += k + 1
                }
                ds.triangles = tris
                
            case "CELL_TYPES":
                guard toks.count >= 2, let n = Int(toks[1]) else { break parseLoop }
                _ = try readInts(n, what: "CELL_TYPES")
                
            case "FIELD":
                guard toks.count >= 3, let nArrays = Int(toks[2]) else { break parseLoop }
                for _ in 0..<nArrays {
                    guard let a = nextTokens(), a.count >= 4,
                          let nComp = Int(a[1]), let nTup = Int(a[2]) else { break parseLoop }
                    let vals = try readFloats(nComp * nTup, type: a[3], what: "FIELD \(a[0])")
                    if a[0] == "LoadParameter", let first = vals.first {
                        ds.loadParameter = Double(first)
                    }
                }
                
            case "CELL_DATA":
                location = .cell
                tupleCount = toks.count > 1 ? (Int(toks[1]) ?? 0) : 0
                
            case "POINT_DATA":
                location = .point
                tupleCount = toks.count > 1 ? (Int(toks[1]) ?? 0) : 0
                
            case "SCALARS":
                guard toks.count >= 3 else { break parseLoop }
                let name = toks[1]
                let nComp = toks.count >= 4 ? (Int(toks[3]) ?? 1) : 1
                // Optional LOOKUP_TABLE line
                let save = pos
                if let lt = nextTokens(), lt[0].uppercased() != "LOOKUP_TABLE" { pos = save }
                let vals = try readFloats(tupleCount * nComp, type: toks[2], what: "SCALARS \(name)")
                if nComp == 1 {
                    addField(name, vals, location)
                } else {
                    for c in 0..<nComp {
                        addField("\(name)_\(c)", stride(from: c, to: vals.count, by: nComp).map { vals[$0] }, location)
                    }
                }
                
            case "VECTORS", "NORMALS":
                guard toks.count >= 3 else { break parseLoop }
                let name = toks[1]
                let vals = try readFloats(tupleCount * 3, type: toks[2], what: "\(key) \(name)")
                if key == "VECTORS" {
                    var xs = [Float](repeating: 0, count: tupleCount)
                    var ys = xs, mag = xs
                    var vecs = [SIMD2<Float>]()
                    vecs.reserveCapacity(tupleCount)
                    for t in 0..<tupleCount {
                        let vx = vals[3 * t], vy = vals[3 * t + 1]
                        xs[t] = vx; ys[t] = vy
                        mag[t] = (vx * vx + vy * vy + vals[3 * t + 2] * vals[3 * t + 2]).squareRoot()
                        vecs.append(SIMD2(vx, vy))
                    }
                    if location == .cell && !allTriangles && vecs.count == keepCell.count {
                        vecs = zip(vecs, keepCell).compactMap { $1 ? $0 : nil }
                    }
                    ds.vectorFields.append(VTKVectorField(name: name, location: location, vectors: vecs))
                    addField("\(name)_magnitude", mag, location)
                    addField("\(name)_x", xs, location)
                    addField("\(name)_y", ys, location)
                }
                
            case "TENSORS":
                guard toks.count >= 3 else { break parseLoop }
                let name = toks[1]
                let vals = try readFloats(tupleCount * 9, type: toks[2], what: "TENSORS \(name)")
                let comps: [(String, Int)] = [("xx", 0), ("xy", 1), ("yx", 3), ("yy", 4)]
                for (suffix, k) in comps {
                    if suffix == "yx" {
                        // only add yx when the tensor is not symmetric
                        var asym = false
                        for t in stride(from: 0, to: tupleCount, by: max(1, tupleCount / 64)) where abs(vals[9 * t + 1] - vals[9 * t + 3]) > 1e-6 {
                            asym = true; break
                        }
                        if !asym { continue }
                    }
                    addField("\(name)_\(suffix)", (0..<tupleCount).map { vals[9 * $0 + k] }, location)
                }
                
            case "LOOKUP_TABLE":
                // Stand-alone colour table: LOOKUP_TABLE name size, then size*4 values
                if toks.count >= 3, let n = Int(toks[2]) {
                    _ = try readFloats(4 * n, type: binary ? "unsigned_char" : "float", what: "LOOKUP_TABLE")
                }
                
            case "METADATA":
                // Skip until an empty line
                while let line = readLine(), !line.trimmingCharacters(in: .whitespaces).isEmpty {}
                
            default:
                // Unknown section – keep what we have rather than failing.
                break parseLoop
            }
        }
        
        // Eigen-mode files carry α and λ_k in the title line ("... alpha 0.0425 lambda 1.01e-3 1.49e-3")
        let titleParts = ds.title.split(whereSeparator: { $0 == " " || $0 == "\t" || $0 == ":" })
        if ds.loadParameter == nil {
            if let i = titleParts.firstIndex(where: { $0.lowercased() == "alpha" }),
               i + 1 < titleParts.count,
               let a = Double(titleParts[i + 1]) {
                ds.loadParameter = a
            }
        }
        if let lamIdx = titleParts.firstIndex(where: { $0.lowercased() == "lambda" || $0.lowercased() == "lambda_min" }) {
            var parsedLambdas: [Double] = []
            var k = lamIdx + 1
            while k < titleParts.count, let v = Double(titleParts[k]) {
                parsedLambdas.append(v)
                k += 1
            }
            ds.modeEigenvalues = parsedLambdas
        }
        
        // Drop triangles referencing missing points (defensive)
        let np = UInt32(ds.points.count)
        if ds.triangles.contains(where: { $0 >= np }) {
            var clean: [UInt32] = []
            clean.reserveCapacity(ds.triangles.count)
            var t = 0
            while t + 2 < ds.triangles.count {
                let a = ds.triangles[t], b = ds.triangles[t + 1], c = ds.triangles[t + 2]
                if a < np && b < np && c < np { clean.append(contentsOf: [a, b, c]) }
                t += 3
            }
            ds.triangles = clean
        }
        return ds
    }
}
