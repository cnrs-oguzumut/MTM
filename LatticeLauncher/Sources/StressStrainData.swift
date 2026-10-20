import Foundation

public struct StressStrainPoint: Identifiable, Equatable {
    public var id: Int { iteration }
    public let iteration: Int
    public let alpha: Double
    public let preStress: Double
    public let postStress: Double
    public let preEnergy: Double
    public let postEnergy: Double
    public let shouldRemesh: Bool
    
    public init(
        iteration: Int,
        alpha: Double,
        preStress: Double,
        postStress: Double,
        preEnergy: Double,
        postEnergy: Double,
        shouldRemesh: Bool
    ) {
        self.iteration = iteration
        self.alpha = alpha
        self.preStress = preStress
        self.postStress = postStress
        self.preEnergy = preEnergy
        self.postEnergy = postEnergy
        self.shouldRemesh = shouldRemesh
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

