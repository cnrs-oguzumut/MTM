#pragma once
#ifndef ELEMENT_TRIANGLE_2D_H
#define ELEMENT_TRIANGLE_2D_H
#include <array>
#include <vector>
#include <functional>
#include <iostream>
#include <Eigen/Dense>
#include "../geometry/Point2D.h"
#include "src/ap.h" // ALGLIB header

class ElementTriangle2D {
private:
    int nn[3]; // Node indices
    Eigen::Vector2d trans[3]; // Translations for each node
    int boundary_node_number;
    double jacobian_det; // Determinant of the Jacobian
    double area; // Element area
    double reference_area; // Reference (original) element area
    bool use_external_deformation; // Flag for using external deformation
    bool shape_derivatives_calculated; // Flag for shape derivatives calculation
    bool use_manual_shape_derivatives; // Flag for manual shape derivatives
    Eigen::Matrix<double, 3, 2> dNdX; // Shape function derivatives
    Eigen::Matrix2d F; // Deformation gradient
    Eigen::Matrix2d C; // Right Cauchy-Green deformation tensor
    Eigen::Matrix2d F_external; // External deformation gradient
    
    // Reduced DOF support
    const std::vector<Point2D>* reference_points_ptr; // Pointer to reference mesh
    const std::vector<std::pair<int, int>>* dof_mapping; // Pointer to DOF mapping

public:
    // Constructor
    ElementTriangle2D();
    
    // External deformation methods
    void setExternalDeformation(const Eigen::Matrix2d& F_ext);
    void disableExternalDeformation();
    const Eigen::Matrix2d& getExternalDeformation() const;
    bool isExternalDeformationActive() const;
    
    // Node index and translation methods
    void setNodeIndex(int position, int index);
    void setTranslation(int position, const Eigen::Vector2d& translation);
    void setBoundaryNodeNumber(int number);
    int getBoundaryNodeNumber();

    int getNodeIndex(int position) const;
    Eigen::Vector2d getTranslation(int position) const;
    Eigen::Vector2d getEffectiveTranslation(int position) const;
    
    // Reference configuration methods
    void setReferenceArea(double ref_area) { reference_area = ref_area; }
    double getReferenceArea() const { return reference_area > 0.0 ? reference_area : area; }
    
    // New area calculation methods
    double calculateReferenceArea(const std::vector<Point2D>& reference_points);
    double calculateCurrentArea(const std::vector<Point2D>& current_points);
    double calculateCurrentArea(const alglib::real_1d_array& current_dofs);

    
    // Reduced DOF methods
    void set_reference_mesh(const std::vector<Point2D>& points);
    void set_dof_mapping(const std::vector<std::pair<int, int>>& mapping);
    
    // Manual shape derivatives methods
    void set_shape_derivatives(const Eigen::Matrix<double, 3, 2>& dNdX_manual);
    bool has_manual_shape_derivatives() const;
    void reset_manual_shape_derivatives();
    
    // Getters
    double getArea() const;
    double getJacobianDet() const;
    const Eigen::Matrix<double, 3, 2>& getDNdX() const;
    const Eigen::Matrix2d& getDeformationGradient() const;
    const Eigen::Matrix2d& getMetricTensor() const;
    bool isInitialized() const;
    
    // Shape derivatives calculation
    double calculate_shape_derivatives(const std::vector<Point2D>& reference_points);
    double calculate_shape_derivatives(const alglib::real_1d_array& free_dofs);
    
    // Deformation gradient calculation
    void calculate_deformation_gradient(const std::vector<Point2D>& current_points);
    void calculate_deformation_gradient(const alglib::real_1d_array& current_dofs);
    
    // Strain calculation
    Eigen::Matrix2d get_green_lagrange_strain() const;
    
    // Force calculations
    std::array<Eigen::Vector2d, 3> calculate_nodal_forces(const Eigen::Matrix2d& P) const;
    std::array<Eigen::Vector2d, 3> calculate_nodal_forces(
        const std::function<Eigen::Matrix2d(const Eigen::Matrix2d&)>& stress_function) const;
    
    // Force assembly (standard vector version only)
    void assemble_forces(const Eigen::Matrix2d& P, std::vector<Eigen::Vector2d>& global_forces) const;
};

#endif // ELEMENT_TRIANGLE_2D_H