#include "acoustic_tensor.h"
#include <iostream>
#include <algorithm>
#include <cmath>

AcousticTensor::AcousticTensor(const Eigen::Matrix2d& F, const Eigen::Matrix2d& C, const Eigen::Matrix2d& Z)
    : F_(F), C_(C), Z_(Z), normalisation_(1.0) {
    dE_dC_.setZero(); // Initialize to zero
    // Initialize empty 4th order tensor for second derivatives
    auto dummy_i = itensor::Index(2, "i");
    auto dummy_j = itensor::Index(2, "j");
    auto dummy_k = itensor::Index(2, "k");
    auto dummy_l = itensor::Index(2, "l");
    dE2_dC2_ = itensor::ITensor(dummy_i, dummy_j, dummy_k, dummy_l);
    initializeIndices();
}

AcousticTensor::AcousticTensor(const Eigen::Matrix2d& F, const Eigen::Matrix2d& C, const Eigen::Matrix2d& Z,
                               const Eigen::Matrix2d& dE_dC, double normalisation)
    : F_(F), C_(C), Z_(Z), dE_dC_(dE_dC), normalisation_(normalisation) {
    // Initialize empty 4th order tensor for second derivatives
    auto dummy_i = itensor::Index(2, "i");
    auto dummy_j = itensor::Index(2, "j");
    auto dummy_k = itensor::Index(2, "k");
    auto dummy_l = itensor::Index(2, "l");    
    dE2_dC2_ = itensor::ITensor(dummy_i, dummy_j, dummy_k, dummy_l);
    initializeIndices();
}

// NEW: Constructor with both first and second derivatives
AcousticTensor::AcousticTensor(const Eigen::Matrix2d& F, const Eigen::Matrix2d& C, const Eigen::Matrix2d& Z,
                               const Eigen::Matrix2d& dE_dC, const itensor::ITensor& dE2_dC2, double normalisation)
    : F_(F), C_(C), Z_(Z), dE_dC_(dE_dC), dE2_dC2_(dE2_dC2), normalisation_(normalisation) {
    initializeIndices();
}

void AcousticTensor::setEnergyDerivative(const Eigen::Matrix2d& dE_dC, double normalisation) {
    dE_dC_ = dE_dC;
    normalisation_ = normalisation;
}

// NEW: Set both first and second derivatives
void AcousticTensor::setEnergyDerivatives(const Eigen::Matrix2d& dE_dC, const itensor::ITensor& dE2_dC2, double normalisation) {
    dE_dC_ = dE_dC;
    dE2_dC2_ = dE2_dC2;
    normalisation_ = normalisation;
}

// NEW (correct):
void AcousticTensor::computeEnergyDerivativeAnalytical(Strain_Energy_LatticeCalculator& calculator,
                                                     const std::function<Eigen::Matrix2d(const Eigen::Matrix2d&)>& dphi_func,
                                                     double normalisation) {    
    dE_dC_ = calculator.calculate_derivative_analytical(C_, dphi_func) / normalisation;
    normalisation_ = normalisation;
}

void AcousticTensor::initializeIndices() {
    i_ = itensor::Index(2, "i");    // Instead of Index("i", 2)
    j_ = itensor::Index(2, "j");    // Instead of Index("j", 2)
    k_ = itensor::Index(2, "k");    // Instead of Index("k", 2)
    l_ = itensor::Index(2, "l");    // Instead of Index("l", 2)
    m_ = itensor::Index(2, "m");    // And so on...
    n_ = itensor::Index(2, "n");
    r_ = itensor::Index(2, "r");
    s_ = itensor::Index(2, "s");
    a_ = itensor::Index(2, "a");
    b_ = itensor::Index(2, "b");
    w1_ = itensor::Index(2, "w1");
    w2_ = itensor::Index(2, "w2");
    w3_ = itensor::Index(2, "w3");
    w4_ = itensor::Index(2, "w4");}

void AcousticTensor::createT54(itensor::ITensor& T54) {
    auto a = itensor::Index(2);
    auto b = itensor::Index(2);
    auto kroneck = itensor::delta(a, b);
    
    for (int ii = 1; ii <= 2; ii++)
        for (int jj = 1; jj <= 2; jj++)
            for (int kk = 1; kk <= 2; kk++)
                for (int ll = 1; ll <= 2; ll++)
                    for (int mu = 1; mu <= 2; mu++)
                        for (int nn = 1; nn <= 2; nn++) {
                            T54.set(m_=ii, n_=jj, r_=kk, s_=ll, i_=mu, j_=nn,
                                kroneck.elt(a=ii, b=ll) * kroneck.elt(a=kk, b=mu) * kroneck.elt(a=jj, b=nn) +
                                kroneck.elt(a=kk, b=mu) * kroneck.elt(a=ii, b=nn) * kroneck.elt(a=jj, b=ll));
                        }
}

void AcousticTensor::createT53(itensor::ITensor& T53, itensor::ITensor& F_tensor) {
    auto a = itensor::Index(2);
    auto b = itensor::Index(2);
    auto kroneck = itensor::delta(a, b);
    
    // Fill F_tensor with deformation gradient values
    F_tensor.set(i_=1, j_=1, F_(0,0));
    F_tensor.set(i_=1, j_=2, F_(0,1));
    F_tensor.set(i_=2, j_=1, F_(1,0));
    F_tensor.set(i_=2, j_=2, F_(1,1));
    
    // Create T53 tensor
    for (int ll = 1; ll <= 2; ll++)
        for (int uu = 1; uu <= 2; uu++)
            for (int jj = 1; jj <= 2; jj++)
                for (int bb = 1; bb <= 2; bb++) {
                    T53.set(m_=ll, n_=uu, i_=jj, j_=bb,
                        kroneck.elt(a=ll, b=bb) * F_tensor.elt(i_=jj, j_=uu) +
                        kroneck.elt(a=uu, b=bb) * F_tensor.elt(i_=jj, j_=ll));
                }
}

// Fix computeJacobian to properly use calculated dE/dC
void AcousticTensor::computeJacobian(itensor::ITensor& dPhi) {
    // std::cout << "Computing Jacobian from calculated first derivative..." << std::endl;
  
        dPhi.set(k_=1, l_=1, dE_dC_(0,0));
        dPhi.set(k_=1, l_=2, dE_dC_(0,1)/2.);
        dPhi.set(k_=2, l_=1, dE_dC_(1,0)/2.);
        dPhi.set(k_=2, l_=2, dE_dC_(1,1));



    // if (has_first_derivative_ && !dE_dC_.isZero()) {
    //     // Use the computed first derivative directly
    //     dPhi.set(k_=1, l_=1, dE_dC_(0,0));
    //     dPhi.set(k_=1, l_=2, dE_dC_(0,1));
    //     dPhi.set(k_=2, l_=1, dE_dC_(1,0));
    //     dPhi.set(k_=2, l_=2, dE_dC_(1,1));
    //     std::cout << "✓ Using calculator-computed first derivative" << std::endl;
    // } else {
    //     std::cout << "Warning: Using fallback for Jacobian" << std::endl;
    //     // Your existing fallback code
    // }
}

// Fix computeHessian to properly map calculated d²E/dC²
void AcousticTensor::computeHessian(itensor::ITensor& ddPhi) {
    double div1 = 0.5;
    double div2 = 0.25;  // 0.5 * 0.5
    // div1=1;
    // div2=1;    
    
    
    for (int k = 1; k <= 2; k++) {
        for (int l = 1; l <= 2; l++) {
            for (int w3 = 1; w3 <= 2; w3++) {
                for (int w4 = 1; w4 <= 2; w4++) {
                    double hess_val = dE2_dC2_.elt(k, l, w3, w4);
                    
                    // Determine division factor based on off-diagonal components
                    double factor = 1.0;
                    
                    // Count how many indices are off-diagonal (i.e., k≠l or w3≠w4)
                    bool kl_offdiag = (k != l);   // True if (k,l) = (1,2) or (2,1)
                    bool w3w4_offdiag = (w3 != w4); // True if (w3,w4) = (1,2) or (2,1)
                    
                    if (kl_offdiag && w3w4_offdiag) {
                        // Both are off-diagonal (c12.c12 type): use div2
                        factor = div2;
                    } else if (kl_offdiag || w3w4_offdiag) {
                        // One is off-diagonal (mixed type): use div1
                        factor = div1;
                    }
                    // else: both diagonal (c11.c11, c22.c22, c11.c22 type): factor = 1.0
                    
                    ddPhi.set(k_=k, l_=l, w3_=w3, w4_=w4, factor * hess_val);
                }
            }
        }
    }
}


// Rest of the implementation remains the same...
// Replace the problematic analyzeAcousticTensor method in acoustic_tensor.cpp with this simpler version:

// MOST IMPORTANT: Replace analyzeAcousticTensor with original working physics
AcousticAnalysis AcousticTensor::analyzeAcousticTensor(bool lagrangian) {
    // std::cout << "Computing acoustic tensor using original working physics..." << std::endl;
    
    // T54 tensor using member indices
    auto T54 = itensor::ITensor(m_, n_, r_, s_, i_, j_);
    T54.fill(0.0);  // Explicit zero fill
    createT54(T54);
    
    // T53 tensor and F tensor using member indices  
    auto F = itensor::ITensor(i_, j_);
    auto T53 = itensor::ITensor(m_, n_, i_, j_);
    T53.fill(0.0);  // Explicit zero fill
    createT53(T53, F);
    
    // Jacobian - using computed first derivative
    auto dPhi = itensor::ITensor(k_, l_);
    dPhi.fill(0.0);  // Explicit zero fill
    computeJacobian(dPhi);
    
    
    // Create ddPhi with the right indices and copy values from dE2_dC2_
    auto ddPhi = itensor::ITensor(k_, l_, w3_, w4_);
    ddPhi.fill(0.0);  // Explicit zero fill
    
    computeHessian(ddPhi);
    
    // std::cout << "Mapped dE2_dC2_ to ddPhi. ddPhi order: " << itensor::order(ddPhi) << std::endl;
    
    // Material tensor using member indices
    auto mm = itensor::ITensor(i_, j_);
    mm.set(i_=1, j_=1, Z_(0,0));
    mm.set(i_=1, j_=2, Z_(0,1));
    mm.set(i_=2, j_=1, Z_(1,0));
    mm.set(i_=2, j_=2, Z_(1,1));
    
    // Material tensor contractions using member indices
    auto mm2 = mm * itensor::delta(i_, n_) * itensor::delta(j_, l_);
    auto mm1 = mm * itensor::delta(i_, m_) * itensor::delta(j_, k_);
    auto mm3 = mm * itensor::delta(i_, w1_) * itensor::delta(j_, w3_);
    auto mm4 = mm * itensor::delta(i_, w2_) * itensor::delta(j_, w4_);
    
    // EXACT ORIGINAL PHYSICS:
    auto f1 = mm1 * mm2 * dPhi * T54;
    
    auto T53_2 = T53 * itensor::delta(m_, w1_) * itensor::delta(n_, w2_) * itensor::delta(i_, r_) * itensor::delta(j_, s_);
    auto f2 = mm1 * mm2 * T53 * T53_2 * mm3 * mm4 * ddPhi;  // Now ddPhi has the right indices
    
    auto final_tensor = f1 + f2;
    AcousticAnalysis result;
    result.detAc = 100000.0;

    // Create additional indices for Eulerian transformation
    if(!lagrangian){
        auto p = itensor::Index(2, "p");
        auto q = itensor::Index(2, "q");
        
        // Create gradient tensors from deformation gradient F_
        auto gradF = itensor::ITensor(p, r_);
        auto gradFr = itensor::ITensor(q, i_);
        
        gradF.set(p=1, r_=1, F_(0,0));
        gradF.set(p=1, r_=2, F_(0,1));
        gradF.set(p=2, r_=1, F_(1,0));
        gradF.set(p=2, r_=2, F_(1,1));
        
        gradFr.set(q=1, i_=1, F_(0,0));
        gradFr.set(q=1, i_=2, F_(0,1));
        gradFr.set(q=2, i_=1, F_(1,0));
        gradFr.set(q=2, i_=2, F_(1,1));
        
        // Permute to correct index order
        final_tensor = itensor::permute(final_tensor, {r_, s_, i_, j_});
        
        // Transform to Eulerian frame
        auto final_tensor_euler = gradF * gradFr * final_tensor;
        final_tensor_euler = itensor::permute(final_tensor_euler, {p, s_, q, j_});
        
        // Direction vectors in Eulerian frame
        auto N1 = itensor::ITensor(p);
        auto N2 = itensor::ITensor(q);
        
        
        for (double ksi = 0; ksi < 2*M_PI; ksi += M_PI / 50.0) {
            N1.set(p=1, cos(ksi));
            N1.set(p=2, sin(ksi));
            N2.set(q=1, cos(ksi));
            N2.set(q=2, sin(ksi));
            
            auto acoustic = final_tensor_euler * N1 * N2;
            
            double detAc = acoustic.elt(s_=1, j_=1) * acoustic.elt(s_=2, j_=2) 
                        - acoustic.elt(s_=1, j_=2) * acoustic.elt(s_=2, j_=1);
            
            if (detAc <= -500) detAc = -500;
            if (detAc >= 800) detAc = 800;
            
            if (detAc < result.detAc) {
                result.detAc = detAc;
                result.xsi = toDegrees(ksi);
            }
        }
    }   
    // Direction vector and search using member indices
    else if(lagrangian){
        auto N = itensor::ITensor(s_);
        
        //for (double ksi = -M_PI; ksi < M_PI; ksi += 0.01) {
        for (double ksi = 0; ksi < 2*M_PI; ksi += M_PI / 60) {

            N.set(s_=1, cos(ksi));
            N.set(s_=2, sin(ksi));
            
            auto acoustic = final_tensor * N * N * itensor::delta(s_, j_);
            
            double detAc = acoustic.elt(r_=1, i_=1) * acoustic.elt(r_=2, i_=2) 
                        - acoustic.elt(r_=1, i_=2) * acoustic.elt(r_=2, i_=1);
            
            if (detAc <= -500) detAc = -500;
            if (detAc >= 800) detAc = 800;
            
            if (detAc < result.detAc) {
                result.detAc = detAc;
                result.xsi = toDegrees(ksi);
            }
        }
    }
        
    
    return result;
}

void AcousticTensor::printHessianComponents() {
    // Create ddPhi tensor with proper indices
    auto ddPhi = itensor::ITensor(k_, l_, w3_, w4_);
    ddPhi.fill(0.0);
    
    // Compute the Hessian
    computeHessian(ddPhi);
    
    // Print components with the same factors as in your example
    std::cout << "\n=== Hessian (ddPhi) Components ===" << std::endl;
    std::cout << "ddPhi(1,1,1,1) = " << 4*ddPhi.elt(k_=1, l_=1, w3_=1, w4_=1) << std::endl;
    std::cout << "ddPhi(2,2,2,2) = " << 4*ddPhi.elt(k_=2, l_=2, w3_=2, w4_=2) << std::endl;
    std::cout << "ddPhi(1,2,1,2) = " << 4*ddPhi.elt(k_=1, l_=2, w3_=1, w4_=2) << std::endl;
    std::cout << "ddPhi(1,1,2,2) = " << 4*ddPhi.elt(k_=1, l_=1, w3_=2, w4_=2) << std::endl;
    std::cout << "ddPhi(1,1,1,2) = " << 4*ddPhi.elt(k_=1, l_=1, w3_=1, w4_=2) << std::endl;
    std::cout << "ddPhi(2,2,1,2) = " << 4*ddPhi.elt(k_=2, l_=2, w3_=1, w4_=2) << std::endl;
    std::cout << "==================================\n" << std::endl;
}
// // Also add this simpler implementation for findMinDetAcousticTensor:

// AcousticAnalysis AcousticTensor::findMinDetAcousticTensor(double xsi) {
//     if (xsi >= 10.0) {
//         // Search mode - find global minimum
//         std::cout << "Performing global minimum search..." << std::endl;
//         return analyzeAcousticTensor(1.0, 0.0);  // The search is already built into analyzeAcousticTensor
//     } else {
//         // Specific angle mode
//         std::cout << "Computing at specific angle: " << xsi << " degrees" << std::endl;
//         double angle_rad = xsi * M_PI / 180.0;
        
//         // Use the same computation as analyzeAcousticTensor but for specific angle
//         Eigen::Vector2d n(cos(angle_rad), sin(angle_rad));
        
//         // Simplified acoustic tensor computation for this specific direction
//         Eigen::Matrix2d stress_term = F_ * dE_dC_;
//         Eigen::Matrix2d acoustic_matrix = Z_ * stress_term;
        
//         AcousticAnalysis result;
//         result.detAc = acoustic_matrix.determinant();
//         result.xsi = xsi;
        
//         return result;
//     }
// }
// Include other methods (computeAijkl, findMinDetAcousticTensor, etc.) from previous implementation...
// [I'll skip them here for brevity as they follow the same pattern]

void AcousticTensor::printInfo() const {
    std::cout << "Acoustic Tensor Information:" << std::endl;
    std::cout << "Deformation Gradient F:" << std::endl << F_ << std::endl;
    std::cout << "Cell Matrix C:" << std::endl << C_ << std::endl;
    std::cout << "Material Matrix Z:" << std::endl << Z_ << std::endl;
    if (!dE_dC_.isZero()) {
        std::cout << "First Energy Derivative dE/dC:" << std::endl << dE_dC_ << std::endl;
    }
    if (itensor::rank(dE2_dC2_) == 4 && itensor::norm(dE2_dC2_) > 1e-12) {
        std::cout << "Second Energy Derivative (Hessian) tensor available (4th order)" << std::endl;
        std::cout << "Hessian tensor norm: " << itensor::norm(dE2_dC2_) << std::endl;
    }
    std::cout << "Normalization: " << normalisation_ << std::endl;
}

// ... (include remaining methods from previous implementations)

itensor::ITensor AcousticTensor::eigenToITensor(const Eigen::Matrix2d& mat, const itensor::Index& idx1, const itensor::Index& idx2) const {
    auto tensor = itensor::ITensor(idx1, idx2);
    for (int ii = 1; ii <= 2; ii++) {
        for (int jj = 1; jj <= 2; jj++) {
            tensor.set(idx1=ii, idx2=jj, mat(ii-1, jj-1));
        }
    }
    return tensor;
}

double AcousticTensor::findMinValue(const std::vector<double>& vec) const {
    return *std::min_element(vec.begin(), vec.end());
}

size_t AcousticTensor::findMinIndex(const std::vector<double>& vec) const {
    return std::distance(vec.begin(), std::min_element(vec.begin(), vec.end()));
}

itensor::ITensor AcousticTensor::getAcousticTensor(bool lagrangian) {
    // Build acoustic tensor exactly as in analyzeAcousticTensor
    auto T54 = itensor::ITensor(m_, n_, r_, s_, i_, j_);
    T54.fill(0.0);
    createT54(T54);
    
    auto F = itensor::ITensor(i_, j_);
    auto T53 = itensor::ITensor(m_, n_, i_, j_);
    T53.fill(0.0);
    createT53(T53, F);
    
    auto dPhi = itensor::ITensor(k_, l_);
    dPhi.fill(0.0);
    computeJacobian(dPhi);
    
    auto ddPhi = itensor::ITensor(k_, l_, w3_, w4_);
    ddPhi.fill(0.0);
    computeHessian(ddPhi);
    
    auto mm = itensor::ITensor(i_, j_);
    mm.set(i_=1, j_=1, Z_(0,0));
    mm.set(i_=1, j_=2, Z_(0,1));
    mm.set(i_=2, j_=1, Z_(1,0));
    mm.set(i_=2, j_=2, Z_(1,1));
    
    auto mm2 = mm * itensor::delta(i_, n_) * itensor::delta(j_, l_);
    auto mm1 = mm * itensor::delta(i_, m_) * itensor::delta(j_, k_);
    auto mm3 = mm * itensor::delta(i_, w1_) * itensor::delta(j_, w3_);
    auto mm4 = mm * itensor::delta(i_, w2_) * itensor::delta(j_, w4_);
    
    auto f1 = mm1 * mm2 * dPhi * T54;
    auto T53_2 = T53 * itensor::delta(m_, w1_) * itensor::delta(n_, w2_) * 
                 itensor::delta(i_, r_) * itensor::delta(j_, s_);
    auto f2 = mm1 * mm2 * T53 * T53_2 * mm3 * mm4 * ddPhi;
    
    auto final_tensor = f1 + f2;
    
    if (!lagrangian) {
        auto p = itensor::Index(2, "p");
        auto q = itensor::Index(2, "q");
        auto gradF = itensor::ITensor(p, r_);
        auto gradFr = itensor::ITensor(q, i_);
        
        gradF.set(p=1, r_=1, F_(0,0));
        gradF.set(p=1, r_=2, F_(0,1));
        gradF.set(p=2, r_=1, F_(1,0));
        gradF.set(p=2, r_=2, F_(1,1));
        
        gradFr.set(q=1, i_=1, F_(0,0));
        gradFr.set(q=1, i_=2, F_(0,1));
        gradFr.set(q=2, i_=1, F_(1,0));
        gradFr.set(q=2, i_=2, F_(1,1));
        
        final_tensor = itensor::permute(final_tensor, {r_, s_, i_, j_});
        auto final_tensor_euler = gradF * gradFr * final_tensor;
        return itensor::permute(final_tensor_euler, {p, s_, q, j_});
    }
    
    return final_tensor;
}