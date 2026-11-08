#include <cassert>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "../include/experiments/acoustic_studies.h"
#include "../include/experiments/common_simulation_helpers.h"
#include "../include/utils/dislocation_utils.h"
#include "../include/geometry/DomainDimensions.h"
#include "../include/geometry/DomainInfo.h"
#include "../include/geometry/LatticeGenerator.h"
#include "../include/interatomic/inter_atomic.h"
#include "../include/lattice_energy/EnergyFunctions.h"
#include "../include/lattice_energy/ReferenceFrame.h"
#include "../include/mesh/ElementTriangle2D.h"
#include "../include/mesh/MeshGenerator.h"
#include "../include/mesh/Remesher.h"
#include "../include/optimization/LatticeOptimizer.h"
#include "../include/output/configuration_saver.h"

namespace {

void test_calculator_equivalence_and_gl2z_invariance() {
    const std::string lattice_type = "triangular";
    const double r_cutoff = 2.5;
    std::function<double(double)> pot = lennard_jones_energy_v2;
    std::function<double(double)> dpot = lennard_jones_energy_der_v2;
    std::function<double(double)> d2pot = lennard_jones_energy_sder_v2;
    const double a0 = find_optimal_lattice_parameter(pot, lattice_type, r_cutoff, dpot, d2pot);

    set_reference_frame(ReferenceFrame::SquareZ2);
    ReferenceSetup sq_setup = make_reference_setup(lattice_type, 1, a0, r_cutoff, pot);
    set_reference_frame(ReferenceFrame::Triangular);
    ReferenceSetup tri_setup = make_reference_setup(lattice_type, 1, a0, r_cutoff, pot);

    const double A_sq = sq_setup.calc->getUnitCellArea();
    const double A_tri = tri_setup.calc->getUnitCellArea();
    assert(std::abs(A_sq - A_tri) < 1e-12);
    assert(std::abs(sq_setup.element_area - tri_setup.element_area) < 1e-12);
    assert(std::abs(sq_setup.zero - tri_setup.zero) < 1e-10);

    // Test a wide set of deformation gradients F including large shears up to gamma ~ 5
    std::vector<Eigen::Matrix2d> test_Fs;
    for (double gamma : {0.0, 0.05, -0.15, 0.4, 1.2, -2.3, 3.7, 5.0}) {
        Eigen::Matrix2d F;
        F << 1.02, gamma,
             -0.03, 0.99;
        test_Fs.push_back(F);
    }
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> uni(-0.25, 0.25);
    std::uniform_real_distribution<double> shear_dist(-5.0, 5.0);
    for (int k = 0; k < 25; ++k) {
        Eigen::Matrix2d F;
        F << 1.0 + uni(rng), shear_dist(rng),
             0.2 * uni(rng), 1.0 + uni(rng);
        if (F.determinant() > 0.4) {
            test_Fs.push_back(F);
        }
    }

    const Eigen::Matrix2d F0 = sq_setup.F0;
    for (const Eigen::Matrix2d& F : test_Fs) {
        const Eigen::Matrix2d C_tri = F.transpose() * F;
        const double e_tri =
            tri_setup.calc->calculate_energy(C_tri, pot, tri_setup.zero) / A_tri;
        const Eigen::Matrix2d dE_tri =
            tri_setup.calc->calculate_derivative(C_tri, dpot) / A_tri;
        const Eigen::Matrix2d P_tri = 2.0 * F * dE_tri;
        const Eigen::Matrix2d sigma_tri = (1.0 / F.determinant()) * P_tri * F.transpose();

        const Eigen::Matrix2d F_sq = F * F0;
        const Eigen::Matrix2d C_sq = F_sq.transpose() * F_sq;
        const lagrange::Result red = lagrange::reduce(C_sq);
        const double e_sq =
            sq_setup.calc->calculate_energy(red.C_reduced, pot, sq_setup.zero) / A_sq;
        const Eigen::Matrix2d dE_sq =
            sq_setup.calc->calculate_derivative(red.C_reduced, dpot) / A_sq;
        const Eigen::Matrix2d P_sq =
            2.0 * F_sq * red.m_matrix * dE_sq * red.m_matrix.transpose();
        const Eigen::Matrix2d sigma_sq = (1.0 / F_sq.determinant()) * P_sq * F_sq.transpose();

        assert(std::abs(e_tri - e_sq) < 1e-10);
        assert((P_tri - P_sq * F0.transpose()).norm() < 1e-10);
        assert((sigma_tri - sigma_sq).norm() < 1e-10);
    }

    // GL(2, Z) invariance in the triangular frame:
    // E_tri(C) == E_tri(B^{-T} m^T B^T C B m B^{-1}) for unimodular integer m
    Eigen::Matrix2d B;
    B << 1.0, 0.5,
         0.0, std::sqrt(3.0) / 2.0;
    const Eigen::Matrix2d B_inv = B.inverse();

    std::vector<Eigen::Matrix2d> unimodular_ms;
    for (int s : {-5, -3, -1, 1, 2, 4, 5}) {
        Eigen::Matrix2d m1, m2, m3;
        m1 << 1, s, 0, 1;
        m2 << 1, 0, s, 1;
        m3 << 0, -1, 1, s;
        unimodular_ms.push_back(m1);
        unimodular_ms.push_back(m2);
        unimodular_ms.push_back(m3);
    }
    Eigen::Matrix2d F_base;
    F_base << 1.04, 0.23,
              -0.05, 0.97;
    const Eigen::Matrix2d C_base = F_base.transpose() * F_base;
    const double e_base = tri_setup.calc->calculate_energy(C_base, pot, tri_setup.zero);

    for (const Eigen::Matrix2d& m : unimodular_ms) {
        const Eigen::Matrix2d M_tri = B * m * B_inv;
        const Eigen::Matrix2d C_equiv = M_tri.transpose() * C_base * M_tri;
        const double e_equiv = tri_setup.calc->calculate_energy(C_equiv, pot, tri_setup.zero);
        assert(std::abs(e_base - e_equiv) < 1e-10);
    }

    std::cout << "[PASS] Triangular vs SquareZ2 equivalence (up to gamma=5) and GL(2,Z) invariance verified."
              << std::endl;
}

void test_pristine_mesh_and_piola_alignment() {
    const std::string lattice_type = "triangular";
    const double r_cutoff = 2.5;
    std::function<double(double)> pot = lennard_jones_energy_v2;
    std::function<double(double)> dpot = lennard_jones_energy_der_v2;
    std::function<double(double)> d2pot = lennard_jones_energy_sder_v2;
    const double a0 = find_optimal_lattice_parameter(pot, lattice_type, r_cutoff, dpot, d2pot);

    set_reference_frame(ReferenceFrame::Triangular);
    ReferenceSetup setup = make_reference_setup(lattice_type, 1, a0, r_cutoff, pot);

    const int nx = 8, ny = 8;
    std::vector<Point2D> points = LatticeGenerator::generate_2d_lattice(nx, ny, a0, lattice_type);
    const std::vector<Point2D> points_ref = points;
    DomainInfo dom_info = compute_domain_size(points);
    DomainDimensions domain_dims(dom_info.get_width(), dom_info.get_height());
    const std::array<double, 2> offsets = {a0 / 2.0, (std::sqrt(3.0) / 2.0) * a0};
    Point2D domain_dims_point(domain_dims.size_x, domain_dims.size_y);
    const bool pbc = true;

    auto [original_domain_map, translation_map] =
        MeshGenerator::create_domain_maps(static_cast<int>(points.size()), domain_dims, offsets);
    auto [interior_mapping, full_mapping] =
        create_dof_mapping_original(points, 1e-5, pbc);

    AdaptiveMesher mesher(domain_dims_point, offsets, original_domain_map,
                          translation_map, full_mapping, 1e-6, pbc);
    mesher.setUsePeriodicCopies(pbc);

    alglib::real_1d_array x;
    const int n_vars = static_cast<int>(interior_mapping.size());
    x.setlength(2 * n_vars);
    map_points_to_solver_array(x, points, interior_mapping, n_vars);

    auto [elements, active_elements] =
        mesher.createMesh(points, x, Eigen::Matrix2d::Identity(), &setup.dndx);
    assert(!active_elements.empty());

    const double norm_area = setup.calc->getUnitCellArea();
    for (size_t idx : active_elements) {
        auto& elem = elements[idx];
        elem.set_dof_mapping(full_mapping);
        elem.set_reference_mesh(points);
        elem.setExternalDeformation(Eigen::Matrix2d::Identity());
        elem.calculate_deformation_gradient(points);

        // 1. Pristine element metric C = I and energy = 0
        const Eigen::Matrix2d C = elem.getMetricTensor();
        assert((C - Eigen::Matrix2d::Identity()).norm() < 1e-12);
        const double e_el = setup.calc->calculate_energy(C, pot, setup.zero) / norm_area;
        assert(std::abs(e_el) < 1e-11);

        // 2. Aligned deformation gradient is I and polar rotation angle is 0 deg
        const Eigen::Matrix2d F_aligned =
            ConfigurationSaver::alignDeformationGradient(elem.getDeformationGradient(), 0.0);
        assert((F_aligned - Eigen::Matrix2d::Identity()).norm() < 1e-12);

        Eigen::JacobiSVD<Eigen::Matrix2d> svd(F_aligned, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Matrix2d Q = svd.matrixU() * svd.matrixV().transpose();
        const double rot_deg = std::atan2(Q(1, 0), Q(0, 0)) * 180.0 / M_PI;
        assert(std::abs(rot_deg) < 1e-10);
    }

    // Check averaged Piola = 0 on pristine mesh
    UserData ud0(points, elements, *setup.calc, pot, dpot, setup.zero, a0,
                 Eigen::Matrix2d::Identity(), interior_mapping, full_mapping, active_elements, false);
    double e0 = 0.0;
    Eigen::Matrix2d sig0 = Eigen::Matrix2d::Zero();
    ConfigurationSaver::calculateEnergyAndStress(&ud0, e0, sig0, true);
    assert(std::abs(e0) < 1e-9);
    assert(sig0.norm() < 1e-9);

    // 3. Small homogeneous shear: every element's aligned P must equal P_phys (no up/down cancellation)
    const double gamma = 0.05;
    Eigen::Matrix2d F_shear;
    F_shear << 1.0, gamma,
               0.0, 1.0;
    for (size_t i = 0; i < points.size(); ++i) {
        points[i].coord = F_shear * points_ref[i].coord;
    }

    const Eigen::Matrix2d C_phys = F_shear.transpose() * F_shear;
    const Eigen::Matrix2d dE_phys = setup.calc->calculate_derivative(C_phys, dpot) / norm_area;
    const Eigen::Matrix2d P_phys = 2.0 * F_shear * dE_phys;
    assert(std::abs(P_phys(0, 1)) > 1e-2);

    Eigen::Matrix2d sum_P = Eigen::Matrix2d::Zero();
    double sum_ref_area = 0.0;
    for (size_t idx : active_elements) {
        auto& elem = elements[idx];
        elem.setExternalDeformation(F_shear);
        elem.calculate_deformation_gradient(points);

        const Eigen::Matrix2d F_al =
            ConfigurationSaver::alignDeformationGradient(elem.getDeformationGradient(), 0.0);
        assert((F_al - F_shear).norm() < 1e-12);

        const Eigen::Matrix2d C_el = F_al.transpose() * F_al;
        const Eigen::Matrix2d dE_el = setup.calc->calculate_derivative(C_el, dpot) / norm_area;
        const Eigen::Matrix2d P_el = 2.0 * F_al * dE_el;
        assert((P_el - P_phys).norm() < 1e-11);

        sum_P += P_el * elem.getReferenceArea();
        sum_ref_area += elem.getReferenceArea();
    }
    const Eigen::Matrix2d P_avg = sum_P / sum_ref_area;
    assert((P_avg - P_phys).norm() < 1e-11);

    set_reference_frame(ReferenceFrame::SquareZ2);
    std::cout << "[PASS] Pristine mesh (C=I, F_aligned=I, angle=0, Piola=0) and sheared Piola alignment verified."
              << std::endl;
}

void test_guards_and_checkpoint_mismatch() {
    const auto pot = lennard_jones_energy_v2;

    // 1. Guards in make_reference_setup
    set_reference_frame(ReferenceFrame::Triangular);
    bool caught_square = false;
    try {
        (void)make_reference_setup("square", 1, 1.0, 2.5, pot);
    } catch (const std::invalid_argument&) {
        caught_square = true;
    }
    assert(caught_square);

    bool caught_zanzotto = false;
    try {
        (void)make_reference_setup("triangular", 3, 1.0, 2.5, pot);
    } catch (const std::invalid_argument&) {
        caught_zanzotto = true;
    }
    assert(caught_zanzotto);

    // 2. Checkpoint reference_frame mismatch rejection
    const std::string tmp_chk = "/tmp/mtm_test_triref_checkpoint.chk";
    ConfigurationSaver::CheckpointData chk_out;
    chk_out.nx = 4;
    chk_out.ny = 4;
    chk_out.iteration = 1;
    chk_out.points.emplace_back(0.0, 0.0);
    chk_out.points.emplace_back(1.0, 0.0);
    chk_out.points.emplace_back(0.5, std::sqrt(3.0) / 2.0);

    set_reference_frame(ReferenceFrame::Triangular);
    ConfigurationSaver::saveCheckpoint(tmp_chk, chk_out);

    Eigen::Matrix<double, 3, 2> dummy_dndx = Eigen::Matrix<double, 3, 2>::Zero();
    ConfigurationSaver::CheckpointData chk_in;
    // Loading in Triangular mode succeeds
    assert(ConfigurationSaver::loadCheckpoint(tmp_chk, chk_in, dummy_dndx));
    assert(chk_in.reference_frame_tag == "triangular");

    // Loading in SquareZ2 mode fails cleanly
    set_reference_frame(ReferenceFrame::SquareZ2);
    assert(!ConfigurationSaver::loadCheckpoint(tmp_chk, chk_in, dummy_dndx));

    std::filesystem::remove(tmp_chk);
    set_reference_frame(ReferenceFrame::SquareZ2);
    std::cout << "[PASS] ReferenceFrame guards and checkpoint reference_frame validation verified." << std::endl;
}

void test_elastic_constants_and_acoustic_consistency() {
    for (int mode : {1, 2, 4}) {
        AcousticStudyConfig cfg;
        cfg.mode = mode;
        cfg.lattice_type = "triangular";

        const AcousticElasticResults ec = computeAcousticElasticConstants(cfg);
        const DislocationSolutionParams dp = DislocationSolutionParams::fromPotentialMode(mode, "triangular");
        const ShearPathResult sp = computeSimpleShearPath(cfg, 0.0, 0.05, 0.05);
        assert(!sp.path.empty());
        const auto& pt0 = sp.path.front();
        assert(std::abs(pt0.alpha) < 1e-12);

        const double c11_path = pt0.q11_xx;
        const double c22_path = pt0.q22_yy;
        const double c66_path = pt0.q11_yy;
        const double c12_path = pt0.q12_xy - pt0.q11_yy;

        assert(std::abs(ec.C11 - dp.C11) < 1e-10);
        assert(std::abs(ec.C12 - dp.C12) < 1e-10);
        assert(std::abs(ec.C66 - dp.C66) < 1e-10);

        assert(std::abs(ec.C11 - c11_path) < 1e-9);
        assert(std::abs(ec.C22 - c22_path) < 1e-9);
        assert(std::abs(ec.C12 - c12_path) < 1e-9);
        assert(std::abs(ec.C66 - c66_path) < 1e-9);
    }
    set_reference_frame(ReferenceFrame::SquareZ2);
    std::cout << "[PASS] Elastic constants consistency (computeAcousticElasticConstants == fromPotentialMode == AcousticTensor at F0) verified."
              << std::endl;
}

} // namespace

void run_triangular_reference_tests() {
    test_calculator_equivalence_and_gl2z_invariance();
    test_pristine_mesh_and_piola_alignment();
    test_guards_and_checkpoint_mismatch();
    test_elastic_constants_and_acoustic_consistency();
}
