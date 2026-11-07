import Foundation

public enum StressComponentType: String, CaseIterable, Identifiable {
    case primary = "Primary σ(α) / Load"
    case cauchyXY = "Cauchy Shear σ_xy"
    case cauchyXX = "Cauchy Normal σ_xx"
    case cauchyYY = "Cauchy Normal σ_yy"
    case hydrostaticPressure = "Hydrostatic Pressure p"
    case piola12 = "Projected Piola ⟨P_12⟩"
    case piola11 = "Piola Normal ⟨P_11⟩"
    case piola22 = "Piola Normal ⟨P_22⟩"
    case vonMisesAvg = "Von Mises ⟨σ_vM⟩ (Avg)"
    case vonMisesMax = "Von Mises max(σ_vM) (Peak)"
    case indenterLoad = "Indenter Load F_y(δ)"
    case meanContactPressure = "Mean Contact Pressure p_m"
    case contactAtoms = "Contact Atoms N_c"
    case effectiveStrain = "Conjugate Strain ε (force / stress control)"

    public var id: String { rawValue }

    public var shortLabel: String {
        switch self {
        case .primary: return "Primary σ"
        case .cauchyXY: return "σ_xy (Cauchy)"
        case .cauchyXX: return "σ_xx (Cauchy)"
        case .cauchyYY: return "σ_yy (Cauchy)"
        case .hydrostaticPressure: return "Pressure p"
        case .piola12: return "P_12 (Projected)"
        case .piola11: return "P_11 (Piola)"
        case .piola22: return "P_22 (Piola)"
        case .vonMisesAvg: return "⟨σ_vM⟩ (Avg)"
        case .vonMisesMax: return "max(σ_vM)"
        case .indenterLoad: return "Load F_y(δ)"
        case .meanContactPressure: return "Contact p_m"
        case .contactAtoms: return "Contact N_c"
        case .effectiveStrain: return "Strain ε"
        }
    }
}

public struct StressStrainPoint: Identifiable, Equatable {
    public var id: Int { iteration }
    public let iteration: Int
    public let alpha: Double
    public let preStress: Double
    public let postStress: Double
    public let preEnergy: Double
    public let postEnergy: Double
    public let shouldRemesh: Bool
    public let cauchyXX: Double
    public let cauchyXY: Double
    public let cauchyYY: Double
    public let hydrostaticPressure: Double
    public let piola11: Double
    public let piola12: Double
    public let piola22: Double
    public let vonMisesAvg: Double
    public let vonMisesMax: Double
    public let indenterLoad: Double
    public let meanContactPressure: Double
    public let contactAtoms: Double
    /// Work-conjugate strain of a dead load, sum f0 . u / A (0 when the run has no traction rules)
    public let effectiveStrain: Double
    /// Elastic energy E_el (force / stress control only; PostEnergy is then the enthalpy H = E_el - W_ext)
    public let elasticEnergy: Double?
    
    public init(
        iteration: Int,
        alpha: Double,
        preStress: Double,
        postStress: Double,
        preEnergy: Double,
        postEnergy: Double,
        shouldRemesh: Bool,
        cauchyXX: Double = 0.0,
        cauchyXY: Double? = nil,
        cauchyYY: Double = 0.0,
        hydrostaticPressure: Double = 0.0,
        piola11: Double = 0.0,
        piola12: Double = 0.0,
        piola22: Double = 0.0,
        vonMisesAvg: Double = 0.0,
        vonMisesMax: Double = 0.0,
        indenterLoad: Double = 0.0,
        meanContactPressure: Double = 0.0,
        contactAtoms: Double = 0.0,
        effectiveStrain: Double = 0.0,
        elasticEnergy: Double? = nil
    ) {
        self.iteration = iteration
        self.alpha = alpha
        self.preStress = preStress
        self.postStress = postStress
        self.preEnergy = preEnergy
        self.postEnergy = postEnergy
        self.shouldRemesh = shouldRemesh
        self.cauchyXX = cauchyXX
        self.cauchyXY = cauchyXY ?? postStress
        self.cauchyYY = cauchyYY
        self.hydrostaticPressure = hydrostaticPressure
        self.piola11 = piola11
        self.piola12 = piola12
        self.piola22 = piola22
        self.vonMisesAvg = vonMisesAvg
        self.vonMisesMax = vonMisesMax
        self.indenterLoad = indenterLoad
        self.meanContactPressure = meanContactPressure
        self.contactAtoms = contactAtoms
        self.effectiveStrain = effectiveStrain
        self.elasticEnergy = elasticEnergy
    }

    public func value(for component: StressComponentType) -> Double {
        switch component {
        case .primary: return postStress
        case .cauchyXY: return cauchyXY
        case .cauchyXX: return cauchyXX
        case .cauchyYY: return cauchyYY
        case .hydrostaticPressure: return hydrostaticPressure
        case .piola12: return piola12
        case .piola11: return piola11
        case .piola22: return piola22
        case .vonMisesAvg: return vonMisesAvg
        case .vonMisesMax: return vonMisesMax
        case .indenterLoad: return indenterLoad
        case .meanContactPressure: return meanContactPressure
        case .contactAtoms: return contactAtoms
        case .effectiveStrain: return effectiveStrain
        }
    }
}

public struct EigenPoint: Identifiable, Equatable {
    public let id: String
    public let iteration: Int
    public let alpha: Double
    public let trigger: String
    public let lambdaMin: Double
    public let lambda2: Double?
    public let failedShifts: Int
    
    public var isUnstable: Bool {
        lambdaMin <= 0.0 || failedShifts > 0
    }
    
    public init(
        iteration: Int,
        alpha: Double,
        trigger: String,
        lambdaMin: Double,
        lambda2: Double? = nil,
        failedShifts: Int = 0
    ) {
        self.id = "\(iteration)-\(trigger)-\(alpha)"
        self.iteration = iteration
        self.alpha = alpha
        self.trigger = trigger
        self.lambdaMin = lambdaMin
        self.lambda2 = lambda2
        self.failedShifts = failedShifts
    }
}

public struct GBEnergyPoint: Identifiable, Equatable {
    public var id: String { String(format: "%.4f", misorientationDeg) }
    public let misorientationDeg: Double
    public let halfAngleDeg: Double
    public let theta1Deg: Double
    public let theta2Deg: Double
    public let gbAngleDeg: Double
    public let nx: Int
    public let ny: Int
    public let lengthGBA0: Double
    public let lengthIntA0: Double
    public let totalEnergy: Double
    public let interiorEnergy: Double
    public let gbEnergyTotal: Double
    public let gbEnergyInterior: Double
    public let normalGapA0: Double
    public let seamDMinA0: Double
    public var sigma: Int
    public var m: Int
    public var n: Int
    public var relaxedVTKPath: String

    public init(
        misorientationDeg: Double,
        halfAngleDeg: Double,
        theta1Deg: Double,
        theta2Deg: Double,
        gbAngleDeg: Double,
        nx: Int,
        ny: Int,
        lengthGBA0: Double,
        lengthIntA0: Double,
        totalEnergy: Double,
        interiorEnergy: Double,
        gbEnergyTotal: Double,
        gbEnergyInterior: Double,
        normalGapA0: Double,
        seamDMinA0: Double,
        sigma: Int? = nil,
        m: Int? = nil,
        n: Int? = nil,
        relaxedVTKPath: String = ""
    ) {
        self.misorientationDeg = misorientationDeg
        self.halfAngleDeg = halfAngleDeg
        self.theta1Deg = theta1Deg
        self.theta2Deg = theta2Deg
        self.gbAngleDeg = gbAngleDeg
        self.nx = nx
        self.ny = ny
        self.lengthGBA0 = lengthGBA0
        self.lengthIntA0 = lengthIntA0
        self.totalEnergy = totalEnergy
        self.interiorEnergy = interiorEnergy
        self.gbEnergyTotal = gbEnergyTotal
        self.gbEnergyInterior = gbEnergyInterior
        self.normalGapA0 = normalGapA0
        self.seamDMinA0 = seamDMinA0
        let match = CSLRationalizer.closestMatch(forHalfAngleDeg: max(0.5, halfAngleDeg), maxIndex: 18)
        let isClose = abs(match.halfAngleDeg - halfAngleDeg) < 0.15
        self.sigma = sigma ?? (isClose ? match.sigma : 0)
        self.m = m ?? (isClose ? match.m : 0)
        self.n = n ?? (isClose ? match.n : 0)
        self.relaxedVTKPath = relaxedVTKPath
    }
}


