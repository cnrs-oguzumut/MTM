#include "../include/mesh/ElementTriangle2D.h"

// Constructor implementation
ElementTriangle2D::ElementTriangle2D() : 
    nn{0, 0, 0}, 
    jacobian_det(0.0), 
    area(0.0),
    reference_area(0.0),
    use_external_deformation(false),
    shape_derivatives_calculated(false) {
    
    for (int i = 0; i < 3; i++) {
        trans[i].setZero();
    }
    
    dNdX.setZero();
    F.setIdentity();
    C.setIdentity();
    F_external.setIdentity();
}

// Set external deformation gradient and activate it
void ElementTriangle2D::setExternalDeformation(const Eigen::Matrix2d& F_ext) {
    F_external = F_ext;
    use_external_deformation = true;
}

// Disable external deformation (will use identity)
void ElementTriangle2D::disableExternalDeformation() {
    use_external_deformation = false;
}

// Get external deformation gradient
const Eigen::Matrix2d& ElementTriangle2D::getExternalDeformation() const {
    return F_external;
}

// Check if external deformation is active
bool ElementTriangle2D::isExternalDeformationActive() const {
    return use_external_deformation;
}

// Setters for node indices and translations
void ElementTriangle2D::setNodeIndex(int position, int index) {
    if (position >= 0 && position < 3) {
        nn[position] = index;
    }
}

void ElementTriangle2D::setTranslation(int position, const Eigen::Vector2d& translation) {
    if (position >= 0 && position < 3) {
        trans[position] = translation;
    }
}

// Getters
int ElementTriangle2D::getNodeIndex(int position) const {
    return (position >= 0 && position < 3) ? nn[position] : -1;
}

Eigen::Vector2d ElementTriangle2D::getTranslation(int position) const {
    return (position >= 0 && position < 3) ? trans[position] : Eigen::Vector2d::Zero();
}

// Get translation vector (deformed if external deformation is active)
Eigen::Vector2d ElementTriangle2D::getEffectiveTranslation(int position) const {
    if (position >= 0 && position < 3) {
        if (use_external_deformation) {
            return F_external * trans[position];
        } else {
            return trans[position];
        }
    }
    return Eigen::Vector2d::Zero();
}

double ElementTriangle2D::getArea() const { 
    return area; 
}

double ElementTriangle2D::getJacobianDet() const { 
    return jacobian_det; 
}

const Eigen::Matrix<double, 3, 2>& ElementTriangle2D::getDNdX() const { 
    return dNdX; 
}

const Eigen::Matrix2d& ElementTriangle2D::getDeformationGradient() const { 
    return F; 
}

const Eigen::Matrix2d& ElementTriangle2D::getMetricTensor() const { 
    return C; 
}

bool ElementTriangle2D::isInitialized() const { 
    return shape_derivatives_calculated; 
}

// Calculate shape function derivatives using reference configuration
double ElementTriangle2D::calculate_shape_derivatives(const std::vector<Point2D>& reference_points) {
    // Create matrix to hold reference coordinates
    Eigen::Matrix<double, 2, 3> X_e;
    
    // Fill the matrix with coordinates from reference configuration
    for (int i = 0; i < 3; i++) {
        // Use effective translation (deformed if active)
        Eigen::Vector2d effective_trans = getEffectiveTranslation(i);
           
        Eigen::Vector2d coord = reference_points[nn[i]].coord + effective_trans;
        X_e.col(i) = coord;  // Standard FEM notation - store as columns
    }
    
    // Define shape functions in natural coordinates
    Eigen::Matrix<double, 3, 2> dNdxi;
    dNdxi << -1.0, -1.0,     // dN1/d(xi,eta)
                1.0,  0.0,     // dN2/d(xi,eta)
                0.0,  1.0;     // dN3/d(xi,eta)
    
    // Calculate Jacobian matrix (standard FEM notation)
    Eigen::Matrix2d J = X_e * dNdxi;
    
    // Calculate determinant and area
    jacobian_det = J.determinant();
    area = std::abs(jacobian_det) / 2.0;
    
    
    // Calculate shape function derivatives in global coordinates
    if (std::abs(jacobian_det) > 1e-10) {
        Eigen::Matrix2d J_inv = J.inverse();
        if (reference_area <= 0.0){
            dNdX = dNdxi ;  // Standard FEM ordering
            // std::cout<<"dNdX"<<dNdX<<std::endl;
        }
        else{
            dNdX = dNdxi;  // Standard FEM ordering
            // std::cout<<"dNdX"<<dNdX<<std::endl;
        }

    } else {
        dNdX.setZero();
        std::cerr << "Warning: Near-singular Jacobian detected." << std::endl;
    }
    
    // Store reference area if not already set
    if (reference_area <= 0.0) {
        reference_area = area;
    }

    // Set flag that shape derivatives have been calculated
    shape_derivatives_calculated = true;
    
    return jacobian_det;
}

// Calculate deformation gradient using current configuration
void ElementTriangle2D::calculate_deformation_gradient(const std::vector<Point2D>& current_points) {
    // Check if shape derivatives have been calculated
    if (!shape_derivatives_calculated) {
        std::cerr << "Error: Shape derivatives not calculated. Call calculate_shape_derivatives first." << std::endl;
        return;
    }
    
    // Create matrix to hold current coordinates
    Eigen::Matrix<double, 2, 3> x_e;
    
    // Fill the matrix with coordinates from current configuration
    for (int i = 0; i < 3; i++) {
        // Use effective translation (deformed if active)
        Eigen::Vector2d effective_trans = getEffectiveTranslation(i);
        Eigen::Vector2d coord = current_points[nn[i]].coord + effective_trans;
        x_e.col(i) = coord;  // Store as columns for standard FEM notation
    }
    
    // Calculate deformation gradient (standard FEM notation)
    F = x_e * dNdX;
    
    // Calculate right Cauchy-Green deformation tensor
    C = F.transpose() * F;
}

// Get Green-Lagrange strain tensor
Eigen::Matrix2d ElementTriangle2D::get_green_lagrange_strain() const {
    Eigen::Matrix2d I = Eigen::Matrix2d::Identity();
    return 0.5 * (C - I);
}

// Calculate nodal forces from first Piola-Kirchhoff stress tensor
std::array<Eigen::Vector2d, 3> ElementTriangle2D::calculate_nodal_forces(const Eigen::Matrix2d& P) const {
    std::array<Eigen::Vector2d, 3> nodal_forces;
    
    // Check if shape derivatives have been calculated
    if (!shape_derivatives_calculated) {
        std::cerr << "Error: Shape derivatives not calculated. Cannot compute nodal forces." << std::endl;
        for (int i = 0; i < 3; i++) {
            nodal_forces[i].setZero();
        }
        return nodal_forces;
    }
    
    // Initialize all force vectors to zero
    for (int a = 0; a < 3; a++) {
        nodal_forces[a].setZero();
    }
    
    // Use reference area instead of current area
    double area_to_use = getReferenceArea();
    
    // For each node (a) and each spatial dimension (i)
    for (int a = 0; a < 3; a++) {
        for (int i = 0; i < 2; i++) {
            // Sum over reference dimensions (j)
            for (int j = 0; j < 2; j++) {
                // Note we add the contribution from each dimension
                nodal_forces[a](i) += P(i,j) * dNdX(a,j) * area_to_use;
            }
        }
    }            
    return nodal_forces;
}

// Overloaded version that takes a stress function (for variable stress within element)
std::array<Eigen::Vector2d, 3> ElementTriangle2D::calculate_nodal_forces(
    const std::function<Eigen::Matrix2d(const Eigen::Matrix2d&)>& stress_function) const {
    // Calculate stress tensor from current deformation gradient
    Eigen::Matrix2d P = stress_function(F);
    
    // Use the existing method to calculate nodal forces
    return calculate_nodal_forces(P);
}

// Implement both versions with the same code
void ElementTriangle2D::assemble_forces(const Eigen::Matrix2d& P, std::vector<Eigen::Vector2d>& global_forces) const {
    if (!shape_derivatives_calculated) {
        std::cerr << "Error: Shape derivatives not calculated. Cannot assemble forces." << std::endl;
        return;
    }
    
    std::array<Eigen::Vector2d, 3> element_forces = calculate_nodal_forces(P);
    
    for (int i = 0; i < 3; i++) {
        global_forces[nn[i]] += element_forces[i];
    }
}

// Implement the aligned version
using aligned_vector = std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>>;
void ElementTriangle2D::assemble_forces(const Eigen::Matrix2d& P, aligned_vector& global_forces) const {
    if (!shape_derivatives_calculated) {
        std::cerr << "Error: Shape derivatives not calculated. Cannot assemble forces." << std::endl;
        return;
    }
    
    std::array<Eigen::Vector2d, 3> element_forces = calculate_nodal_forces(P);
    
    for (int i = 0; i < 3; i++) {
        global_forces[nn[i]] += element_forces[i];
    }
}