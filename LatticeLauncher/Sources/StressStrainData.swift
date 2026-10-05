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
