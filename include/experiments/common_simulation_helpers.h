#pragma once

#include <Eigen/Dense>
#include <array>
#include <cstddef>
#include <tuple>
#include <vector>

#include "../geometry/Point2D.h"
#include "../mesh/ElementTriangle2D.h"
#include "../optimization/LatticeOptimizer.h"

Eigen::Matrix<double, 3, 2>
calculateShapeDerivatives(const Eigen::Vector2d &p1, const Eigen::Vector2d &p2,
                          const Eigen::Vector2d &p3);

bool hasConnectivityChanged(const std::vector<ElementTriangle2D> &old_elements,
                            const std::vector<ElementTriangle2D> &new_elements,
                            const std::vector<size_t> &old_active,
                            const std::vector<size_t> &new_active);

std::tuple<double, Eigen::Matrix2d, int> perform_remeshing_loop_reduction(
    alglib::real_1d_array &x, const UserData *userData,
    const std::vector<int> &contact_atoms,
    const std::vector<int> &boundary_fixed_nodes, const Eigen::Matrix2d &F_ext,
    const Eigen::Matrix<double, 3, 2> &dndx,
    const std::array<double, 2> &offsets,
    const std::vector<int> &original_domain_map,
    const std::vector<std::tuple<double, double>> &translation_map,
    const Point2D &domain_dims_point, int &has_changes,
    int max_iterations, double reference_area, bool pbc,
    bool optimize_interior);

std::tuple<double, Eigen::Matrix2d, int> perform_remeshing_loop(
    alglib::real_1d_array &x, UserData *userData,
    const std::vector<int> &contact_atoms,
    const std::vector<int> &boundary_fixed_nodes, const Eigen::Matrix2d &F_ext,
    const Eigen::Matrix<double, 3, 2> &dndx,
    const std::array<double, 2> &offsets,
    const std::vector<int> &original_domain_map,
    const std::vector<std::tuple<double, double>> &translation_map,
    const Point2D &domain_dims_point, int max_iterations,
    double reference_area);

void writeSizesToFile(int Nx, int Ny);

size_t findMiddleAtom(const std::vector<Point2D> &points,
                      bool verbose = false);

std::vector<Point2D> scaleLattice(const std::vector<Point2D> &original_points,
                                  double scale_factor);

std::vector<Point2D>
scaleLatticeAroundPoint(const std::vector<Point2D> &original_points,
                        double scale_factor,
                        const Eigen::Vector2d &reference_point);
