#include "../include/mesh/ElementTriangle2D.h"

// Constructor implementation
ElementTriangle2D::ElementTriangle2D() : 
    nn{0, 0, 0}, 
    jacobian_det(0.0), 
    area(0.0),
    reference_area(0.0),
    use_external_deformation(false),
    shape_derivatives_calculated(false),
    reference_points_ptr(nullptr),
    dof_mapping(nullptr) {
    
    for (int i = 0; i < 3; i++) {
        trans[i].setZero();
    }
    
    dNdX.setZero();
    F.setIdentity();
    C.setIdentity();
    F_external.setIdentity();
}

// Reference configuration methods
void ElementTriangle2D::set_reference_mesh(const std::vector<Point2D>& points) {
    reference_points_ptr = &points;
}

void ElementTriangle2D::set_dof_mapping(const std::vector<std::pair<int, int>>& mapping) {
    dof_mapping = &mapping;
}

// External deformation methods
void ElementTriangle2D::setExternalDeformation(const Eigen::Matrix2d& F_ext) {
    F_external = F_ext;
    use_external_deformation = true;
}

void ElementTriangle2D::disableExternalDeformation() {
    use_external_deformation = false;
}

const Eigen::Matrix2d& ElementTriangle2D::getExternalDeformation() const {
    return F_external;
}

bool ElementTriangle2D::isExternalDeformationActive() const {
    return use_external_deformation;
}

// Node management
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

int ElementTriangle2D::getNodeIndex(int position) const {
    return (position >= 0 && position < 3) ? nn[position] : -1;
}

Eigen::Vector2d ElementTriangle2D::getTranslation(int position) const {
    return (position >= 0 && position < 3) ? trans[position] : Eigen::Vector2d::Zero();
}

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

// Property accessors
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

// Shape derivatives calculation
double ElementTriangle2D::calculate_shape_derivatives(const alglib::real_1d_array& free_dofs) {
    if (!dof_mapping || !reference_points_ptr) {
        std::cerr << "Error: DOF mapping or reference mesh not set." << std::endl;
        return 0.0;
    }

    Eigen::Matrix<double, 2, 3> X_e;
    const int n_free = free_dofs.length() / 2;
    for (int i = 0; i < 3; i++) {
        // Get effective translation
        Eigen::Vector2d effective_trans = getEffectiveTranslation(i);
        
        // Get the node index
        const int global_node_idx = nn[i];
        const auto& [original_idx, solver_idx] = (*dof_mapping)[global_node_idx];
        
        Eigen::Vector2d coord;
        if (solver_idx == -1) {
            // For fixed nodes, use reference position
            coord = (*reference_points_ptr)[original_idx].coord;
        } else {
            // For free nodes, use the position directly from free_dofs
            coord = Eigen::Vector2d(free_dofs[solver_idx], free_dofs[n_free + solver_idx]);
        }
        
        // Add effective translation
        coord += effective_trans;
        
        // Store as columns
        X_e.col(i) = coord;
    }    // for (int i = 0; i < 3; i++) {
    //     const int global_node_idx = nn[i];
    //     const auto& [original_idx, solver_idx] = (*dof_mapping)[global_node_idx];
    //     Eigen::Vector2d coord;
        
    //     if (solver_idx == -1) {
    //         // For fixed nodes, use reference position
    //         coord = (*reference_points_ptr)[original_idx].coord;
    //     } else {
    //         // For free nodes, directly use the position from free_dofs
    //         coord = Eigen::Vector2d(free_dofs[solver_idx], free_dofs[n_free + solver_idx]);
    //     }
        
    //     // Add effective translation
    //     coord += getEffectiveTranslation(i);
    //     X_e.col(i) = coord;
    // }
    Eigen::Matrix<double, 3, 2> dNdxi;
    dNdxi << -1.0, -1.0, 1.0, 0.0, 0.0, 1.0;
    
    Eigen::Matrix2d J = X_e * dNdxi;
    jacobian_det = J.determinant();
    area = std::abs(jacobian_det) / 2.0;

    if (std::abs(jacobian_det) > 1e-10) {
        dNdX = dNdxi * J.inverse();
        //dNdX = dNdxi ;
    } else {
        dNdX.setZero();
        std::cerr << "Warning: Near-singular Jacobian detected." << std::endl;
    }

    if (reference_area <= 0.0) {
        reference_area = area;
    }

    shape_derivatives_calculated = true;
    return jacobian_det;
}

double ElementTriangle2D::calculate_shape_derivatives(const std::vector<Point2D>& reference_points) {
    Eigen::Matrix<double, 2, 3> X_e;
    
    for (int i = 0; i < 3; i++) {
        Eigen::Vector2d effective_trans = getEffectiveTranslation(i);    
        std::cout << "Node " << i << " (global idx: " << nn[i] << ") effective translation: (" 
        << effective_trans.x() << ", " << effective_trans.y() << ")\n";
   
        Eigen::Vector2d coord = reference_points[nn[i]].coord + effective_trans;
        X_e.col(i) = coord;
    }
    
    Eigen::Matrix<double, 3, 2> dNdxi;
    dNdxi << -1.0, -1.0, 1.0, 0.0, 0.0, 1.0;
    
    Eigen::Matrix2d J = X_e * dNdxi;
    jacobian_det = J.determinant();
    area = std::abs(jacobian_det) / 2.0;

    if (std::abs(jacobian_det) > 1e-10) {
        dNdX = dNdxi * J.inverse();
        //dNdX = dNdxi ;
    } else {
        dNdX.setZero();
        std::cerr << "Warning: Near-singular Jacobian detected." << std::endl;
    }

    if (reference_area <= 0.0) {
        reference_area = area;
    }

    shape_derivatives_calculated = true;
    return jacobian_det;
}

// Deformation gradient calculation
void ElementTriangle2D::calculate_deformation_gradient(const alglib::real_1d_array& free_dofs) {
    if (!shape_derivatives_calculated) {
        std::cerr << "Error: Shape derivatives not calculated." << std::endl;
        return;
    }
    if (!dof_mapping || !reference_points_ptr) {
        std::cerr << "Error: DOF mapping or reference mesh not set." << std::endl;
        return;
    }

    Eigen::Matrix<double, 2, 3> x_e;
    const int n_free = free_dofs.length() / 2;

    for (int i = 0; i < 3; i++) {
        const int global_node_idx = nn[i];
        const auto& [original_idx, solver_idx] = (*dof_mapping)[global_node_idx];
        Eigen::Vector2d coord;
        
        if (solver_idx == -1) {
            // For fixed nodes, use reference position
            coord = (*reference_points_ptr)[original_idx].coord;
        } else {
            // For free nodes, use the position directly instead of treating as displacement
            coord = Eigen::Vector2d(free_dofs[solver_idx], free_dofs[n_free + solver_idx]);
        }
        
        coord += getEffectiveTranslation(i);
        x_e.col(i) = coord;
    }

    F = x_e * dNdX;
    C = F.transpose() * F;
}

void ElementTriangle2D::calculate_deformation_gradient(const std::vector<Point2D>& current_points) {
    if (!shape_derivatives_calculated) {
        std::cerr << "Error: Shape derivatives not calculated." << std::endl;
        return;
    }
    
    Eigen::Matrix<double, 2, 3> x_e;
    
    for (int i = 0; i < 3; i++) {
        Eigen::Vector2d effective_trans = getEffectiveTranslation(i);
        Eigen::Vector2d coord = current_points[nn[i]].coord + effective_trans;
        x_e.col(i) = coord;
    }
    
    F = x_e * dNdX;
    C = F.transpose() * F;
}

// Strain calculation
Eigen::Matrix2d ElementTriangle2D::get_green_lagrange_strain() const {
    Eigen::Matrix2d I = Eigen::Matrix2d::Identity();
    return 0.5 * (C - I);
}

// Force calculations
std::array<Eigen::Vector2d, 3> ElementTriangle2D::calculate_nodal_forces(const Eigen::Matrix2d& P) const {
    std::array<Eigen::Vector2d, 3> nodal_forces;
    
    if (!shape_derivatives_calculated) {
        std::cerr << "Error: Shape derivatives not calculated." << std::endl;
        for (auto& f : nodal_forces) f.setZero();
        return nodal_forces;
    }
    
    double area_to_use = reference_area > 0.0 ? reference_area : area;
    
    for (int a = 0; a < 3; a++) {
        nodal_forces[a] = (P * dNdX.row(a).transpose()) * area_to_use;
    }
            
    return nodal_forces;
}

std::array<Eigen::Vector2d, 3> ElementTriangle2D::calculate_nodal_forces(
    const std::function<Eigen::Matrix2d(const Eigen::Matrix2d&)>& stress_function) const {
    return calculate_nodal_forces(stress_function(F));
}

// Force assembly
void ElementTriangle2D::assemble_forces(const Eigen::Matrix2d& P, std::vector<Eigen::Vector2d>& global_forces) const {
    if (!shape_derivatives_calculated) {
        std::cerr << "Error: Shape derivatives not calculated." << std::endl;
        return;
    }
    
    auto element_forces = calculate_nodal_forces(P);
    for (int i = 0; i < 3; i++) {
        global_forces[nn[i]] += element_forces[i];
    }
}