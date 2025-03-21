#pragma once
#ifndef ELEMENT_TRIANGLE_2D_H
#define ELEMENT_TRIANGLE_2D_H

#include <array>
#include <vector>
#include <functional>
#include <iostream>
#include <Eigen/Dense>
#include "../geometry/Point2D.h"

class ElementTriangle2D {
private:
    std::array<int, 3> nn;                        // Node indices
    std::array<Eigen::Vector2d, 3> trans;         // Translation vectors
    double jacobian_det;                          // Determinant of Jacobian
    double area;                                  // Element area
    Eigen::Matrix<double, 3, 2> dNdX;             // Shape function derivatives (standard FEM notation)
    Eigen::Matrix2d F;                            // Deformation gradient tensor
    Eigen::Matrix2d C;                            // Right Cauchy-Green tensor
    Eigen::Matrix2d F_external;                   // External deformation gradient
    bool use_external_deformation;                // Flag to activate external deformation
    bool shape_derivatives_calculated;            // Flag to track if shape derivatives have been calculated

public:
    // Constructor
    ElementTriangle2D();
    
    // Set external deformation gradient and activate it
    void setExternalDeformation(const Eigen::Matrix2d& F_ext);
    
    // Disable external deformation (will use identity)
    void disableExternalDeformation();
    
    // Get external deformation gradient
    const Eigen::Matrix2d& getExternalDeformation() const;
    
    // Check if external deformation is active
    bool isExternalDeformationActive() const;
    
    // Setters for node indices and translations
    void setNodeIndex(int position, int index);
    void setTranslation(int position, const Eigen::Vector2d& translation);
    
    // Getters
    int getNodeIndex(int position) const;
    Eigen::Vector2d getTranslation(int position) const;
    Eigen::Vector2d getEffectiveTranslation(int position) const;
    double getArea() const;
    double getJacobianDet() const;
    const Eigen::Matrix<double, 3, 2>& getDNdX() const;
    const Eigen::Matrix2d& getDeformationGradient() const;
    const Eigen::Matrix2d& getMetricTensor() const;
    bool isInitialized() const;
    
    // Calculate shape function derivatives using reference configuration
    double calculate_shape_derivatives(const std::vector<Point2D>& reference_points);
    
    // Calculate deformation gradient using current configuration
    void calculate_deformation_gradient(const std::vector<Point2D>& current_points);
    
    // Get Green-Lagrange strain tensor
    Eigen::Matrix2d get_green_lagrange_strain() const;
    
    // Calculate nodal forces from first Piola-Kirchhoff stress tensor
    std::array<Eigen::Vector2d, 3> calculate_nodal_forces(const Eigen::Matrix2d& P) const;
    
    // Overloaded version that takes a stress function (for variable stress within element)
    std::array<Eigen::Vector2d, 3> calculate_nodal_forces(
        const std::function<Eigen::Matrix2d(const Eigen::Matrix2d&)>& stress_function) const;
    
    // Calculate and assemble nodal forces into a global force vector
    void assemble_forces(const Eigen::Matrix2d& P, std::vector<Eigen::Vector2d>& global_forces) const;
};

#endif // ELEMENT_TRIANGLE_2D_H
