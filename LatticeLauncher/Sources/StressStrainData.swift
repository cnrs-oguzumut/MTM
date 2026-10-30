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
        contactAtoms: Double = 0.0
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

