#!/usr/bin/env python3
"""
plot_shifted_dislocation_300x300.py

Comprehensive analysis of shifted dislocation cores (s = 0, 1, 2, 3, 4, 5):
1. Slip plane energy profiles (Full [-150, 150]h and Zoom [-25, 25]h)
2. Separate figures for Cauchy stress components:
   - sigma_xy (shear stress traction along slip plane: Full & Zoom)
   - sigma_xx (normal stress: Upper vs. Lower half-planes: Full & Zoom)
   - sigma_yy (normal stress: Upper vs. Lower half-planes: Full & Zoom)
3. Cumulative strain energy vs ln(R/h)
"""

import os
import struct
import numpy as np
import matplotlib.pyplot as plt

# Aesthetics
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['grid.color'] = '#e0e0e0'
plt.rcParams['grid.linestyle'] = '--'
plt.rcParams['grid.alpha'] = 0.6
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

def read_vtk_data(path, N=90000):
    with open(path, 'rb') as f:
        content = f.read()
    
    # 1. Coordinates
    idx_pts = content.find(f'POINTS {N} float'.encode())
    pts_start = content.find(b'\n', idx_pts) + 1
    pts = np.array(struct.unpack(f'>{3*N}f', content[pts_start : pts_start + 3*N*4])).reshape((N, 3))
    
    # 2. Nodal energy
    idx_e = content.find(b'SCALARS NodalEnergy float 1')
    h_e = content.find(b'\n', content.find(b'\n', idx_e) + 1) + 1
    e = np.array(struct.unpack(f'>{N}f', content[h_e : h_e + N*4]))
    
    # 3. Cells and Areas
    idx_cd = content.find(b'CELL_DATA ')
    idx_c = content.find(b'CELLS ')
    h_c = content.find(b'\n', idx_c)
    _, n_cells, n_ints = content[idx_c:h_c].decode('latin1').split()
    n_cells, n_ints = int(n_cells), int(n_ints)
    cell_data = np.array(struct.unpack(f'>{n_ints}i', content[h_c+1 : h_c+1 + n_ints*4])).reshape((n_cells, 4))
    cells = cell_data[:, 1:4]
    
    # 4. Element Cauchy stresses
    idx_es = content.find(b'TENSORS CauchyStress float', idx_cd)
    h_es = content.find(b'\n', idx_es) + 1
    elem_stresses = np.array(struct.unpack(f'>{9*n_cells}f', content[h_es : h_es + 9*n_cells*4])).reshape((n_cells, 3, 3))
    
    tri_centers = np.mean(pts[cells], axis=1)
    p0 = pts[cells[:, 0]]
    p1 = pts[cells[:, 1]]
    p2 = pts[cells[:, 2]]
    areas = 0.5 * np.abs((p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1]) - (p2[:, 0] - p0[:, 0]) * (p1[:, 1] - p0[:, 1]))
    
    return pts, e, cells, tri_centers, areas, elem_stresses

def extract_slip_plane_rows(pts, energies, y_center=149.5):
    mask_bot = np.abs(pts[:, 1] - (y_center - 0.5)) < 0.45
    mask_top = np.abs(pts[:, 1] - (y_center + 0.5)) < 0.45
    
    order_bot = np.argsort(pts[mask_bot, 0])
    x_bot = pts[mask_bot, 0][order_bot]
    e_bot = energies[mask_bot][order_bot]
    
    order_top = np.argsort(pts[mask_top, 0])
    x_top = pts[mask_top, 0][order_top]
    e_top = energies[mask_top][order_top]
    
    x_mid = 0.5 * (x_bot + x_top)
    e_mid = 0.5 * (e_bot + e_top)
    return x_mid, e_mid

def compute_slip_centroid(x, e):
    e_max = np.max(e)
    mask = e >= 0.5 * e_max
    return np.sum(x[mask] * e[mask]) / np.sum(e[mask])

def compute_stress_profiles(tri_centers, areas, elem_stresses, xc, yc=149.5, x_span=25.0, num_x=201, sigma=1.0):
    strip_mask = np.abs(tri_centers[:, 1] - yc) < 6.0 * sigma
    strip_centers = tri_centers[strip_mask]
    strip_areas = areas[strip_mask]
    strip_stresses = elem_stresses[strip_mask]
    
    x_eval = np.linspace(xc - x_span, xc + x_span, num_x)
    
    mask_up = strip_centers[:, 1] >= yc
    mask_dn = strip_centers[:, 1] <= yc
    
    up_c, up_a, up_s = strip_centers[mask_up], strip_areas[mask_up], strip_stresses[mask_up]
    dn_c, dn_a, dn_s = strip_centers[mask_dn], strip_areas[mask_dn], strip_stresses[mask_dn]
    
    s_xx_up, s_yy_up, s_xy_mid = [], [], []
    s_xx_dn, s_yy_dn = [], []
    
    for x0 in x_eval:
        d2_u = (up_c[:, 0] - x0)**2 + (up_c[:, 1] - (yc + 0.5))**2
        w_u = up_a * np.exp(-0.5 * d2_u / (sigma**2))
        wu_sum = np.sum(w_u)
        s_xx_up.append(np.sum(up_s[:, 0, 0] * w_u) / wu_sum)
        s_yy_up.append(np.sum(up_s[:, 1, 1] * w_u) / wu_sum)
        s_xy_u = np.sum(up_s[:, 0, 1] * w_u) / wu_sum
        
        d2_d = (dn_c[:, 0] - x0)**2 + (dn_c[:, 1] - (yc - 0.5))**2
        w_d = dn_a * np.exp(-0.5 * d2_d / (sigma**2))
        wd_sum = np.sum(w_d)
        s_xx_dn.append(np.sum(dn_s[:, 0, 0] * w_d) / wd_sum)
        s_yy_dn.append(np.sum(dn_s[:, 1, 1] * w_d) / wd_sum)
        s_xy_d = np.sum(dn_s[:, 0, 1] * w_d) / wd_sum
        
        s_xy_mid.append(0.5 * (s_xy_u + s_xy_d))
        
    return {
        'x_rel': x_eval - xc,
        'x_abs': x_eval,
        's_xy': np.array(s_xy_mid),
        's_xx_up': np.array(s_xx_up),
        's_xx_dn': np.array(s_xx_dn),
        's_yy_up': np.array(s_yy_up),
        's_yy_dn': np.array(s_yy_dn)
    }

def main():
    import sys
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if len(sys.argv) > 1 and os.path.isdir(sys.argv[1]):
        base_dir = os.path.abspath(sys.argv[1])
    else:
        base_dir = os.path.join(script_dir, 'shifted_dislocation_study_300x300')
    artifact_dir = '/Users/usalman/.gemini/antigravity-ide/brain/fac2e960-26a9-4cc5-af93-fee8a8b36021'
    shifts = [0, 1, 2, 3, 4, 5]
    
    # Elegant, distinguishable color palette
    colors = ['#1f4e79', '#008080', '#e05638', '#d9822b', '#7b2cbf', '#c1121f']
    markers = ['o', 's', '^', 'D', 'v', 'p']
    
    data = {}
    print("Parsing simulation data for all shifts...")
    for s in shifts:
        vtk_path = os.path.join(base_dir, f'shift_{s}/vtk_output/configuration_00001.vtk')
        if not os.path.exists(vtk_path):
            print(f"Warning: {vtk_path} does not exist.")
            continue
        pts, e, cells, tri_c, areas, stresses = read_vtk_data(vtk_path)
        x_slip, e_slip = extract_slip_plane_rows(pts, e)
        x_cm = compute_slip_centroid(x_slip, e_slip)
        
        print(f"Computing stress profiles for Shift s={s} (centroid x_cm={x_cm:.3f})...")
        # 1. Centroid-aligned stress profiles (relative x in [-25, 25])
        st_zoom = compute_stress_profiles(tri_c, areas, stresses, xc=x_cm, x_span=25.0, num_x=201)
        st_full = compute_stress_profiles(tri_c, areas, stresses, xc=x_cm, x_span=150.0, num_x=301, sigma=1.5)
        
        # 2. Unaligned stress profiles without x_c (fixed domain center xc=149.5)
        st_abs_zoom = compute_stress_profiles(tri_c, areas, stresses, xc=149.5, x_span=25.0, num_x=201)
        st_abs_full = compute_stress_profiles(tri_c, areas, stresses, xc=149.5, x_span=150.0, num_x=301, sigma=1.5)
        
        r_from_cm = np.linalg.norm(pts[:, :2] - np.array([x_cm, 149.5]), axis=1)
        r_from_fixed = np.linalg.norm(pts[:, :2] - np.array([149.5, 149.5]), axis=1)
        
        data[s] = {
            'pts': pts,
            'e': e,
            'x_slip': x_slip,
            'e_slip': e_slip,
            'x_cm': x_cm,
            'r_from_cm': r_from_cm,
            'r_from_fixed': r_from_fixed,
            'st_zoom': st_zoom,
            'st_full': st_full,
            'st_abs_zoom': st_abs_zoom,
            'st_abs_full': st_abs_full
        }
        print(f"Shift s={s}: centroid x_cm = {x_cm:.4f}, max energy = {np.max(e_slip):.5f}")

    # Helper function to save figure both locally and to artifact dir (both PNG and vectorial PDF)
    def save_fig(fig, filename):
        base_name, _ = os.path.splitext(filename)
        pdf_name = f"{base_name}.pdf"
        png_name = f"{base_name}.png"

        # 1. Base directory (shifted_dislocation_study_300x300)
        p1_png = os.path.join(base_dir, png_name)
        p1_pdf = os.path.join(base_dir, pdf_name)
        fig.savefig(p1_png, dpi=300)
        fig.savefig(p1_pdf)
        print(f"Saved: {p1_png} and {p1_pdf}")

        # 2. Artifact dir if it exists
        if os.path.exists(artifact_dir):
            p2_png = os.path.join(artifact_dir, png_name)
            p2_pdf = os.path.join(artifact_dir, pdf_name)
            fig.savefig(p2_png, dpi=300)
            fig.savefig(p2_pdf)

        # 3. Paper figures dir if relevant
        paper_fig_dir = '/Users/usalman/current_papers/Mesh_reconnection_MTM_ordered_submission/figures'
        if os.path.exists(paper_fig_dir):
            p3_png = os.path.join(paper_fig_dir, png_name)
            if os.path.exists(p3_png):
                p3_pdf = os.path.join(paper_fig_dir, pdf_name)
                fig.savefig(p3_pdf)
                print(f"Also saved to paper figures: {p3_pdf}")

        plt.close(fig)

    # ==============================================================================
    # FIGURE 1: Full Slip Plane Energy Profile (Across entire domain [-150, 150]h)
    # ==============================================================================
    fig1, ax1 = plt.subplots(figsize=(10, 5.5), dpi=300)
    for s, col in zip(shifts, colors):
        if s not in data: continue
        x_centered = data[s]['x_slip'] - data[s]['x_cm']
        ax1.plot(x_centered, data[s]['e_slip'], color=col, lw=2.0, label=f'Shift $s = {s}$')

    ax1.axvline(x=0.0, color='#333333', linestyle=':', lw=1.5, label='Dislocation Center')
    ax1.axvspan(-70.0, 70.0, color='#f0ad4e', alpha=0.12, label='Relaxation Zone')
    ax1.set_xlim(-150, 150)
    ax1.set_ylim(bottom=-0.002, top=0.075)
    ax1.set_xlabel(r'Relative Coordinate Along Slip Plane $(x - x_c) / h$', fontsize=13, fontweight='bold', labelpad=8)
    ax1.set_ylabel('Slip Plane Energy $E(x)$', fontsize=13, fontweight='bold', labelpad=8)
    ax1.set_title(r'Shifted Dislocation: Aligned Slip Plane Energy ($300 \times 300$, $R_{\rm free} = 70h$)', fontsize=14, fontweight='bold', pad=12)
    ax1.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.95, edgecolor='#cccccc', fontsize=10, ncol=2)
    ax1.grid(True)
    plt.tight_layout()
    save_fig(fig1, 'dislocation_300_shifted_slip_plane_full.png')

    # ==============================================================================
    # FIGURE 2: Zoomed Core Region Aligned at Centroid (x - x_cm) in [-25, 25]h
    # ==============================================================================
    fig2, ax2 = plt.subplots(figsize=(10, 5.5), dpi=300)
    for s, col, mark in zip(shifts, colors, markers):
        if s not in data: continue
        x_centered = data[s]['x_slip'] - data[s]['x_cm']
        mask_z = (x_centered >= -25.0) & (x_centered <= 25.0)
        ax2.plot(x_centered[mask_z], data[s]['e_slip'][mask_z], '-', color=col, lw=2.0, ms=4.5, 
                 marker=mark, label=f'Shift $s = {s}$')

    ax2.axvline(x=0.0, color='#333333', linestyle=':', lw=1.5, label='Dislocation Center')
    ax2.set_xlim(-25, 25)
    ax2.set_ylim(bottom=-0.002, top=0.075)
    ax2.set_xlabel(r'Relative Coordinate Along Slip Plane $(x - x_c) / h$', fontsize=13, fontweight='bold', labelpad=8)
    ax2.set_ylabel('Slip Plane Energy $E(x)$', fontsize=13, fontweight='bold', labelpad=8)
    ax2.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.95, edgecolor='#cccccc', fontsize=10, ncol=2)
    ax2.grid(True)
    plt.tight_layout()
    save_fig(fig2, 'dislocation_300_shifted_slip_plane_zoom.png')

    # ==============================================================================
    # Helper for Normal Stress Legend Construction
    # ==============================================================================
    from matplotlib.lines import Line2D
    def make_normal_stress_legend(ax, shifts_present, is_full=False):
        handles = []
        for s in shifts_present:
            handles.append(Line2D([0], [0], color=colors[shifts.index(s)], lw=2.0, label=f'Shift $s = {s}$'))
        handles.append(Line2D([0], [0], color='#444444', lw=2.0, ls='-', label='Upper ($y > y_c$)'))
        handles.append(Line2D([0], [0], color='#444444', lw=1.8, ls=':', label='Lower ($y < y_c$)'))
        handles.append(Line2D([0], [0], color='#333333', lw=1.5, ls=':', label='Dislocation Center'))
        if is_full:
            import matplotlib.patches as mpatches
            handles.append(mpatches.Patch(facecolor='#f0ad4e', alpha=0.2, edgecolor='none', label='Relaxation Zone'))
        ax.legend(handles=handles, loc='upper right', frameon=True, facecolor='white', 
                  framealpha=0.95, edgecolor='#cccccc', fontsize=9.0, ncol=2)

    # ==============================================================================
    # FIGURE 3: Zoomed Slip Plane Cauchy Shear Stress sigma_xy Around Core
    # ==============================================================================
    fig3, ax3 = plt.subplots(figsize=(10, 5.5), dpi=300)
    for s, col in zip(shifts, colors):
        if s not in data: continue
        st = data[s]['st_zoom']
        ax3.plot(st['x_rel'], st['s_xy'], color=col, lw=2.2, label=f'Shift $s = {s}$')

    ax3.axhline(y=0.0, color='#666666', linestyle='-', lw=0.8, alpha=0.6)
    ax3.axvline(x=0.0, color='#333333', linestyle=':', lw=1.5, label='Dislocation Center')
    ax3.set_xlim(-25, 25)
    ax3.set_ylim(-0.20, 0.20)
    ax3.set_xlabel(r'Relative Coordinate Along Slip Plane $(x - x_c) / h$', fontsize=13, fontweight='bold', labelpad=8)
    ax3.set_ylabel(r'Cauchy Shear Stress $\sigma_{xy}$', fontsize=13, fontweight='bold', labelpad=8)
    ax3.set_title(r'Shifted Dislocation: Slip Plane Cauchy Shear Stress $\sigma_{xy}$ ($(x - x_c) \in [-25, 25]h$)', fontsize=14, fontweight='bold', pad=12)
    ax3.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.95, edgecolor='#cccccc', fontsize=10, ncol=2)
    ax3.grid(True)
    plt.tight_layout()
    save_fig(fig3, 'dislocation_300_shifted_slip_plane_stress_xy_zoom.png')

    # ==============================================================================
    # FIGURE 4: Full Slip Plane Cauchy Shear Stress sigma_xy Across Domain
    # ==============================================================================
    fig4, ax4 = plt.subplots(figsize=(10, 5.5), dpi=300)
    for s, col in zip(shifts, colors):
        if s not in data: continue
        st = data[s]['st_full']
        ax4.plot(st['x_rel'], st['s_xy'], color=col, lw=2.0, label=f'Shift $s = {s}$')

    ax4.axhline(y=0.0, color='#666666', linestyle='-', lw=0.8, alpha=0.6)
    ax4.axvline(x=0.0, color='#333333', linestyle=':', lw=1.5, label='Dislocation Center')
    ax4.axvspan(-70.0, 70.0, color='#f0ad4e', alpha=0.12, label='Relaxation Zone')
    ax4.set_xlim(-150, 150)
    ax4.set_ylim(-0.20, 0.20)
    ax4.set_xlabel(r'Relative Coordinate Along Slip Plane $(x - x_c) / h$', fontsize=13, fontweight='bold', labelpad=8)
    ax4.set_ylabel(r'Cauchy Shear Stress $\sigma_{xy}$', fontsize=13, fontweight='bold', labelpad=8)
    ax4.set_title(r'Shifted Dislocation: Full Domain Shear Stress $\sigma_{xy}$ ($300 \times 300$, $R_{\rm free} = 70h$)', fontsize=14, fontweight='bold', pad=12)
    ax4.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.95, edgecolor='#cccccc', fontsize=10, ncol=2)
    ax4.grid(True)
    plt.tight_layout()
    save_fig(fig4, 'dislocation_300_shifted_slip_plane_stress_xy_full.png')

    # ==============================================================================
    # FIGURE 5: Zoomed Slip Plane Cauchy Normal Stress sigma_xx Around Core
    # ==============================================================================
    fig5, ax5 = plt.subplots(figsize=(10, 5.8), dpi=300)
    for s, col in zip(shifts, colors):
        if s not in data: continue
        st = data[s]['st_zoom']
        # Upper (compression) as solid line
        ax5.plot(st['x_rel'], st['s_xx_up'], color=col, lw=2.0)
        # Lower (tension) as dotted line
        ax5.plot(st['x_rel'], st['s_xx_dn'], color=col, lw=1.6, ls=':', alpha=0.85)

    ax5.axhline(0.0, color='#666666', linestyle='-', lw=0.8, alpha=0.6)
    ax5.axvline(0.0, color='#333333', linestyle=':', lw=1.5)
    ax5.set_xlim(-25, 25)
    ax5.set_ylim(-0.65, 0.65)
    ax5.set_xlabel(r'Relative Coordinate Along Slip Plane $(x - x_c) / h$', fontsize=13, fontweight='bold', labelpad=8)
    ax5.set_ylabel(r'Cauchy Normal Stress $\sigma_{xx}$', fontsize=13, fontweight='bold', labelpad=8)
    ax5.set_title(r'Shifted Dislocation: Core Zoom Normal Stress $\sigma_{xx}$ (Upper vs. Lower Half-Planes)', fontsize=14, fontweight='bold', pad=12)
    make_normal_stress_legend(ax5, [s for s in shifts if s in data], is_full=False)
    ax5.grid(True)
    plt.tight_layout()
    save_fig(fig5, 'dislocation_300_shifted_slip_plane_stress_xx_zoom.png')

    # ==============================================================================
    # FIGURE 6: Full Slip Plane Cauchy Normal Stress sigma_xx Across Domain
    # ==============================================================================
    fig6, ax6 = plt.subplots(figsize=(10, 5.8), dpi=300)
    for s, col in zip(shifts, colors):
        if s not in data: continue
        st = data[s]['st_full']
        ax6.plot(st['x_rel'], st['s_xx_up'], color=col, lw=2.0)
        ax6.plot(st['x_rel'], st['s_xx_dn'], color=col, lw=1.6, ls=':', alpha=0.85)

    ax6.axhline(0.0, color='#666666', linestyle='-', lw=0.8, alpha=0.6)
    ax6.axvline(0.0, color='#333333', linestyle=':', lw=1.5)
    ax6.axvspan(-70.0, 70.0, color='#f0ad4e', alpha=0.12)
    ax6.set_xlim(-150, 150)
    ax6.set_ylim(-0.65, 0.65)
    ax6.set_xlabel(r'Relative Coordinate Along Slip Plane $(x - x_c) / h$', fontsize=13, fontweight='bold', labelpad=8)
    ax6.set_ylabel(r'Cauchy Normal Stress $\sigma_{xx}$', fontsize=13, fontweight='bold', labelpad=8)
    ax6.set_title(r'Shifted Dislocation: Full Domain Normal Stress $\sigma_{xx}$ (Upper vs. Lower Half-Planes)', fontsize=14, fontweight='bold', pad=12)
    make_normal_stress_legend(ax6, [s for s in shifts if s in data], is_full=True)
    ax6.grid(True)
    plt.tight_layout()
    save_fig(fig6, 'dislocation_300_shifted_slip_plane_stress_xx_full.png')

    # ==============================================================================
    # FIGURE 7: Zoomed Slip Plane Cauchy Normal Stress sigma_yy Around Core
    # ==============================================================================
    fig7, ax7 = plt.subplots(figsize=(10, 5.8), dpi=300)
    for s, col in zip(shifts, colors):
        if s not in data: continue
        st = data[s]['st_zoom']
        ax7.plot(st['x_rel'], st['s_yy_up'], color=col, lw=2.0)
        ax7.plot(st['x_rel'], st['s_yy_dn'], color=col, lw=1.6, ls=':', alpha=0.85)

    ax7.axhline(0.0, color='#666666', linestyle='-', lw=0.8, alpha=0.6)
    ax7.axvline(0.0, color='#333333', linestyle=':', lw=1.5)
    ax7.set_xlim(-25, 25)
    ax7.set_ylim(-0.40, 0.40)
    ax7.set_xlabel(r'Relative Coordinate Along Slip Plane $(x - x_c) / h$', fontsize=13, fontweight='bold', labelpad=8)
    ax7.set_ylabel(r'Cauchy Normal Stress $\sigma_{yy}$', fontsize=13, fontweight='bold', labelpad=8)
    ax7.set_title(r'Shifted Dislocation: Core Zoom Normal Stress $\sigma_{yy}$ (Upper vs. Lower Half-Planes)', fontsize=14, fontweight='bold', pad=12)
    make_normal_stress_legend(ax7, [s for s in shifts if s in data], is_full=False)
    ax7.grid(True)
    plt.tight_layout()
    save_fig(fig7, 'dislocation_300_shifted_slip_plane_stress_yy_zoom.png')

    # ==============================================================================
    # FIGURE 8: Full Slip Plane Cauchy Normal Stress sigma_yy Across Domain
    # ==============================================================================
    fig8, ax8 = plt.subplots(figsize=(10, 5.8), dpi=300)
    for s, col in zip(shifts, colors):
        if s not in data: continue
        st = data[s]['st_full']
        ax8.plot(st['x_rel'], st['s_yy_up'], color=col, lw=2.0)
        ax8.plot(st['x_rel'], st['s_yy_dn'], color=col, lw=1.6, ls=':', alpha=0.85)

    ax8.axhline(0.0, color='#666666', linestyle='-', lw=0.8, alpha=0.6)
    ax8.axvline(0.0, color='#333333', linestyle=':', lw=1.5)
    ax8.axvspan(-70.0, 70.0, color='#f0ad4e', alpha=0.12)
    ax8.set_xlim(-150, 150)
    ax8.set_ylim(-0.40, 0.40)
    ax8.set_xlabel(r'Relative Coordinate Along Slip Plane $(x - x_c) / h$', fontsize=13, fontweight='bold', labelpad=8)
    ax8.set_ylabel(r'Cauchy Normal Stress $\sigma_{yy}$', fontsize=13, fontweight='bold', labelpad=8)
    ax8.set_title(r'Shifted Dislocation: Full Domain Normal Stress $\sigma_{yy}$ (Upper vs. Lower Half-Planes)', fontsize=14, fontweight='bold', pad=12)
    make_normal_stress_legend(ax8, [s for s in shifts if s in data], is_full=True)
    ax8.grid(True)
    plt.tight_layout()
    save_fig(fig8, 'dislocation_300_shifted_slip_plane_stress_yy_full.png')

    # ==============================================================================
    # FIGURE 9: Semi-Log Cumulative Energy E(R) vs ln(R/h)
    # ==============================================================================
    radii = np.linspace(0.5, 70.0, 100)
    fig9, ax9 = plt.subplots(figsize=(10, 5.5), dpi=300)

    for s, col, mark in zip(shifts, colors, markers):
        if s not in data: continue
        r = data[s]['r_from_cm']
        e = data[s]['e']
        cum_e = np.array([np.sum(e[r <= R]) for R in radii])
        ax9.plot(radii, cum_e, 'o-', color=col, lw=1.8, ms=3.0, alpha=0.9, label=f'Shift $s = {s}$')

    if 0 in data:
        r0_data = data[0]['r_from_cm']
        e0_data = data[0]['e']
        cum_e0 = np.array([np.sum(e0_data[r0_data <= R]) for R in radii])
        mask_fit = (radii >= 10.0) & (radii <= 60.0)
        fit_p = np.polyfit(np.log(radii[mask_fit]), cum_e0[mask_fit], 1)
        fit_r = np.linspace(5.0, 68.0, 100)
        ax9.plot(fit_r, fit_p[0] * np.log(fit_r) + fit_p[1], 'k-.', lw=2.0, label=f'Slope Fit ({fit_p[0]:.4f})')

    ax9.axvline(x=5.0, color='gray', ls=':', lw=1.5, label='Core Cutoff (5h)')
    ax9.set_xscale('log')
    ax9.set_xlim(0.5, 75.0)
    ax9.set_xlabel(r'Radius from Core $R / h$', fontsize=13, fontweight='bold', labelpad=8)
    ax9.set_ylabel(r'Cumulative Strain Energy $\sum_{r_i \leq R} E_i$', fontsize=13, fontweight='bold', labelpad=8)
    ax9.set_title(r'Semi-Log: Cumulative Strain Energy $E(R)$ vs. $\ln(R/h)$ ($300 \times 300$, $R_{\rm free} = 70h$)', fontsize=14, fontweight='bold', pad=12)
    ax9.legend(loc='lower right', frameon=True, facecolor='white', framealpha=0.95, edgecolor='#cccccc', fontsize=10, ncol=2)
    ax9.grid(True, which='both')
    plt.tight_layout()
    save_fig(fig9, 'dislocation_300_shifted_log_cumulative_energy.png')

    # ==============================================================================
    # FIGURES WITHOUT x_c SUBTRACTION (RAW / ABSOLUTE LATTICE COORDINATES)
    # ==============================================================================

    # FIGURE 10: Full Slip Plane Energy Profile without x_c (Absolute x in [0, 300]h)
    fig10, ax10 = plt.subplots(figsize=(10, 5.5), dpi=300)
    for s, col in zip(shifts, colors):
        if s not in data: continue
        ax10.plot(data[s]['x_slip'], data[s]['e_slip'], color=col, lw=2.0, label=f'Shift $s = {s}$')

    ax10.axvline(x=149.5, color='#333333', linestyle=':', lw=1.5, label='Crystal Center ($149.5h$)')
    ax10.axvspan(149.5 - 70.0, 149.5 + 70.0, color='#f0ad4e', alpha=0.12, label='Relaxation Zone')
    ax10.set_xlim(0, 300)
    ax10.set_ylim(bottom=-0.002, top=0.075)
    ax10.set_xlabel(r'Absolute Coordinate Along Slip Plane $x / h$', fontsize=13, fontweight='bold', labelpad=8)
    ax10.set_ylabel('Slip Plane Energy $E(x)$', fontsize=13, fontweight='bold', labelpad=8)
    ax10.set_title(r'Shifted Dislocation: Absolute Slip Plane Energy ($300 \times 300$, without $x_c$ subtraction)', fontsize=14, fontweight='bold', pad=12)
    ax10.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.95, edgecolor='#cccccc', fontsize=10, ncol=2)
    ax10.grid(True)
    plt.tight_layout()
    save_fig(fig10, 'dislocation_300_shifted_slip_plane_absolute_full.png')

    # FIGURE 11: Zoomed Core Slip Plane Energy without x_c (Absolute x in [135, 165]h)
    fig11, ax11 = plt.subplots(figsize=(10, 5.5), dpi=300)
    for s, col, mark in zip(shifts, colors, markers):
        if s not in data: continue
        x_s = data[s]['x_slip']
        e_s = data[s]['e_slip']
        mask_z = (x_s >= 135.0) & (x_s <= 165.0)
        ax11.plot(x_s[mask_z], e_s[mask_z], '-', color=col, lw=2.0, ms=4.5, 
                  marker=mark, label=f'Shift $s = {s}$ ($x_c = {data[s]["x_cm"]:.2f}h$)')

    ax11.axvline(x=149.5, color='#333333', linestyle=':', lw=1.5, label='Crystal Center ($149.5h$)')
    ax11.set_xlim(135, 165)
    ax11.set_ylim(bottom=-0.002, top=0.075)
    ax11.set_xlabel(r'Absolute Coordinate Along Slip Plane $x / h$', fontsize=13, fontweight='bold', labelpad=8)
    ax11.set_ylabel('Slip Plane Energy $E(x)$', fontsize=13, fontweight='bold', labelpad=8)
    ax11.set_title(r'Shifted Dislocation: Core Zoom without $x_c$ Subtraction (Gliding Along Slip Plane)', fontsize=14, fontweight='bold', pad=12)
    ax11.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.95, edgecolor='#cccccc', fontsize=9.5, ncol=2)
    ax11.grid(True)
    plt.tight_layout()
    save_fig(fig11, 'dislocation_300_shifted_slip_plane_absolute_zoom.png')

    # FIGURE 12: Zoomed Cauchy Shear Stress sigma_xy without x_c (Absolute x in [135, 165]h)
    fig12, ax12 = plt.subplots(figsize=(10, 5.5), dpi=300)
    for s, col in zip(shifts, colors):
        if s not in data: continue
        st = data[s]['st_abs_zoom']
        ax12.plot(st['x_abs'], st['s_xy'], color=col, lw=2.2, label=f'Shift $s = {s}$ ($x_c = {data[s]["x_cm"]:.2f}h$)')

    ax12.axhline(y=0.0, color='#666666', linestyle='-', lw=0.8, alpha=0.6)
    ax12.axvline(x=149.5, color='#333333', linestyle=':', lw=1.5, label='Crystal Center ($149.5h$)')
    ax12.set_xlim(135, 165)
    ax12.set_ylim(-0.20, 0.20)
    ax12.set_xlabel(r'Absolute Coordinate Along Slip Plane $x / h$', fontsize=13, fontweight='bold', labelpad=8)
    ax12.set_ylabel(r'Cauchy Shear Stress $\sigma_{xy}$', fontsize=13, fontweight='bold', labelpad=8)
    ax12.set_title(r'Shifted Dislocation: Core Zoom Shear Stress $\sigma_{xy}$ without $x_c$ Subtraction', fontsize=14, fontweight='bold', pad=12)
    ax12.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.95, edgecolor='#cccccc', fontsize=9.5, ncol=2)
    ax12.grid(True)
    plt.tight_layout()
    save_fig(fig12, 'dislocation_300_shifted_slip_plane_stress_xy_absolute_zoom.png')

    # FIGURE 13: Zoomed Cauchy Normal Stress sigma_xx without x_c (Absolute x in [135, 165]h)
    fig13, ax13 = plt.subplots(figsize=(10, 5.8), dpi=300)
    for s, col in zip(shifts, colors):
        if s not in data: continue
        st = data[s]['st_abs_zoom']
        ax13.plot(st['x_abs'], st['s_xx_up'], color=col, lw=2.0)
        ax13.plot(st['x_abs'], st['s_xx_dn'], color=col, lw=1.6, ls=':', alpha=0.85)

    ax13.axhline(0.0, color='#666666', linestyle='-', lw=0.8, alpha=0.6)
    ax13.axvline(149.5, color='#333333', linestyle=':', lw=1.5)
    ax13.set_xlim(135, 165)
    ax13.set_ylim(-0.65, 0.65)
    ax13.set_xlabel(r'Absolute Coordinate Along Slip Plane $x / h$', fontsize=13, fontweight='bold', labelpad=8)
    ax13.set_ylabel(r'Cauchy Normal Stress $\sigma_{xx}$', fontsize=13, fontweight='bold', labelpad=8)
    ax13.set_title(r'Shifted Dislocation: Core Zoom Normal Stress $\sigma_{xx}$ without $x_c$ Subtraction', fontsize=14, fontweight='bold', pad=12)
    make_normal_stress_legend(ax13, [s for s in shifts if s in data], is_full=False)
    ax13.grid(True)
    plt.tight_layout()
    save_fig(fig13, 'dislocation_300_shifted_slip_plane_stress_xx_absolute_zoom.png')

    # FIGURE 14: Zoomed Cauchy Normal Stress sigma_yy without x_c (Absolute x in [135, 165]h)
    fig14, ax14 = plt.subplots(figsize=(10, 5.8), dpi=300)
    for s, col in zip(shifts, colors):
        if s not in data: continue
        st = data[s]['st_abs_zoom']
        ax14.plot(st['x_abs'], st['s_yy_up'], color=col, lw=2.0)
        ax14.plot(st['x_abs'], st['s_yy_dn'], color=col, lw=1.6, ls=':', alpha=0.85)

    ax14.axhline(0.0, color='#666666', linestyle='-', lw=0.8, alpha=0.6)
    ax14.axvline(149.5, color='#333333', linestyle=':', lw=1.5)
    ax14.set_xlim(135, 165)
    ax14.set_ylim(-0.40, 0.40)
    ax14.set_xlabel(r'Absolute Coordinate Along Slip Plane $x / h$', fontsize=13, fontweight='bold', labelpad=8)
    ax14.set_ylabel(r'Cauchy Normal Stress $\sigma_{yy}$', fontsize=13, fontweight='bold', labelpad=8)
    ax14.set_title(r'Shifted Dislocation: Core Zoom Normal Stress $\sigma_{yy}$ without $x_c$ Subtraction', fontsize=14, fontweight='bold', pad=12)
    make_normal_stress_legend(ax14, [s for s in shifts if s in data], is_full=False)
    ax14.grid(True)
    plt.tight_layout()
    save_fig(fig14, 'dislocation_300_shifted_slip_plane_stress_yy_absolute_zoom.png')

    print("All shifted figures successfully generated!")

if __name__ == '__main__':
    main()
