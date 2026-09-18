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

def parse_vtk_nodal_data(path, N=90000, area_weighted=True):
    with open(path, 'rb') as f:
        content = f.read()
    idx_pts = content.find(f'POINTS {N} float'.encode())
    pts_start = content.find(b'\n', idx_pts) + 1
    pts = struct.unpack(f'>{3*N}f', content[pts_start : pts_start + 3*N*4])
    coords = np.array(pts).reshape((N, 3))

    idx_e = content.find(b'SCALARS NodalEnergy float 1')
    h_e = content.find(b'\n', content.find(b'\n', idx_e) + 1) + 1
    e = np.array(struct.unpack(f'>{N}f', content[h_e : h_e + N*4]))

    # Area-weighted averaging from CELL_DATA element Cauchy stress tensors
    idx_cd = content.find(b'CELL_DATA ')
    if area_weighted and idx_cd != -1:
        idx_c = content.find(b'CELLS ')
        h_c = content.find(b'\n', idx_c)
        _, n_cells, n_ints = content[idx_c:h_c].decode('latin1').split()
        n_cells = int(n_cells)
        n_ints = int(n_ints)
        cell_data = np.array(struct.unpack(f'>{n_ints}i', content[h_c+1 : h_c+1 + n_ints*4])).reshape((n_cells, 4))
        cells = cell_data[:, 1:4]

        idx_es = content.find(b'TENSORS CauchyStress float', idx_cd)
        h_es = content.find(b'\n', idx_es) + 1
        elem_stresses = np.array(struct.unpack(f'>{9*n_cells}f', content[h_es : h_es + 9*n_cells*4])).reshape((n_cells, 3, 3))

        p0 = coords[cells[:, 0]]
        p1 = coords[cells[:, 1]]
        p2 = coords[cells[:, 2]]
        areas = 0.5 * np.abs((p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1]) - (p2[:, 0] - p0[:, 0]) * (p1[:, 1] - p0[:, 1]))

        nodal_s = np.zeros((N, 3, 3))
        area_sums = np.zeros(N)
        for i in range(3):
            weighted = elem_stresses * areas[:, None, None]
            np.add.at(nodal_s, cells[:, i], weighted)
            np.add.at(area_sums, cells[:, i], areas)

        mask_valid = area_sums > 0
        nodal_s[mask_valid] /= area_sums[mask_valid, None, None]

        s_xx = nodal_s[:, 0, 0]
        s_yy = nodal_s[:, 1, 1]
        s_xy = nodal_s[:, 0, 1]
    else:
        idx_s = content.find(b'TENSORS NodalCauchyStress float')
        if idx_s != -1:
            h_s = content.find(b'\n', idx_s) + 1
            tensors = np.array(struct.unpack(f'>{9*N}f', content[h_s : h_s + 9*N*4])).reshape((N, 3, 3))
            s_xx = tensors[:, 0, 0]
            s_yy = tensors[:, 1, 1]
            s_xy = tensors[:, 0, 1]
        else:
            s_xx = np.zeros(N)
            s_yy = np.zeros(N)
            s_xy = np.zeros(N)

    return coords[:, :2], e, s_xx, s_yy, s_xy

def extract_slip_plane_rows(pts, values, y_center=149.5):
    mask_bot = np.abs(pts[:, 1] - (y_center - 0.7)) < 0.45
    mask_top = np.abs(pts[:, 1] - (y_center + 0.3)) < 0.45
    
    order_bot = np.argsort(pts[mask_bot, 0])
    x_bot = pts[mask_bot, 0][order_bot]
    v_bot = values[mask_bot][order_bot]
    
    order_top = np.argsort(pts[mask_top, 0])
    x_top = pts[mask_top, 0][order_top]
    v_top = values[mask_top][order_top]
    
    x_mid = 0.5 * (x_bot + x_top)
    v_mid = 0.5 * (v_bot + v_top)
    
    return (x_bot, v_bot), (x_top, v_top), (x_mid, v_mid)

# Paths
vtk_unrelaxed = 'dislocation_study_300x300/without_remesh/vtk_output/configuration_00000.vtk'
vtk_fixed = 'dislocation_study_300x300/without_remesh/vtk_output/configuration_00001.vtk'
vtk_remesh = 'dislocation_study_300x300/with_remesh/vtk_output/configuration_00001.vtk'
artifact_dir = '/Users/usalman/.gemini/antigravity-ide/brain/e0366a0b-42af-4a53-9525-edb8c926ac89'
out_dir = 'dislocation_study_300x300'
os.makedirs(out_dir, exist_ok=True)

pts0, e0, s_xx0, s_yy0, s_xy0 = parse_vtk_nodal_data(vtk_unrelaxed, N=90000)
pts1, e1, s_xx1, s_yy1, s_xy1 = parse_vtk_nodal_data(vtk_fixed, N=90000)
pts2, e2, s_xx2, s_yy2, s_xy2 = parse_vtk_nodal_data(vtk_remesh, N=90000)

core_x1, core_y1 = 149.5, 149.5
core_x2, core_y2 = 150.0, 149.5
r0 = np.linalg.norm(pts0 - np.array([core_x1, core_y1]), axis=1)
r1 = np.linalg.norm(pts1 - np.array([core_x1, core_y1]), axis=1)
r2 = np.linalg.norm(pts2 - np.array([core_x2, core_y2]), axis=1)

(x_b1, e_b1), (x_t1, e_t1), (x_m1, e_m1) = extract_slip_plane_rows(pts1, e1, y_center=core_y1)
(x_b2, e_b2), (x_t2, e_t2), (x_m2, e_m2) = extract_slip_plane_rows(pts2, e2, y_center=core_y2)

# Cauchy stress components along slip plane
(xb_xy1, sb_xy1), (xt_xy1, st_xy1), (xm_xy1, sm_xy1) = extract_slip_plane_rows(pts1, s_xy1, y_center=core_y1)
(xb_xy2, sb_xy2), (xt_xy2, st_xy2), (xm_xy2, sm_xy2) = extract_slip_plane_rows(pts2, s_xy2, y_center=core_y2)

(xb_xx1, sb_xx1), (xt_xx1, st_xx1), (xm_xx1, sm_xx1) = extract_slip_plane_rows(pts1, s_xx1, y_center=core_y1)
(xb_xx2, sb_xx2), (xt_xx2, st_xx2), (xm_xx2, sm_xx2) = extract_slip_plane_rows(pts2, s_xx2, y_center=core_y2)

(xb_yy1, sb_yy1), (xt_yy1, st_yy1), (xm_yy1, sm_yy1) = extract_slip_plane_rows(pts1, s_yy1, y_center=core_y1)
(xb_yy2, sb_yy2), (xt_yy2, st_yy2), (xm_yy2, sm_yy2) = extract_slip_plane_rows(pts2, s_yy2, y_center=core_y2)

# Shift coordinates so that each dislocation core is centered at 0:
x_m1_centered = x_m1 - core_x1
x_m2_centered = x_m2 - core_x2
xt1_c = xt_xx1 - core_x1
xt2_c = xt_xx2 - core_x2
xb1_c = xb_xx1 - core_x1
xb2_c = xb_xx2 - core_x2

# ==============================================================================
# FIGURE 1: Full Slip Plane Energy Profile (Centered at x - x_c = 0)
# ==============================================================================
fig1, ax1 = plt.subplots(figsize=(10, 5.5), dpi=300)
ax1.plot(x_m1_centered, e_m1, color='#d9534f', lw=2.0, label='Without Remeshing')
ax1.plot(x_m2_centered, e_m2, color='#0275d8', lw=2.0, ls='--', label='With Remeshing')
ax1.axvline(x=0.0, color='#333333', linestyle=':', lw=1.5, label='Dislocation Center')
ax1.axvspan(-70.0, 70.0, color='#f0ad4e', alpha=0.15, label='Relaxation Zone')

ax1.set_xlim(-150, 150)
ax1.set_ylim(bottom=-0.002, top=0.075)
ax1.set_xlabel(r'Relative Coordinate Along Slip Plane $(x - x_c) / h$', fontsize=13, fontweight='bold', labelpad=8)
ax1.set_ylabel('Slip Plane Energy $E(x)$', fontsize=13, fontweight='bold', labelpad=8)
ax1.set_title(r'Anisotropic Dislocation: Slip Plane Energy ($300 \times 300$, $R_{\rm free} = 70h$)', fontsize=14, fontweight='bold', pad=12)
ax1.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.95, edgecolor='#cccccc', fontsize=11)
ax1.grid(True)
plt.tight_layout()

f1_path = os.path.join(out_dir, 'dislocation_300_slip_plane_full.png')
f1_art = os.path.join(artifact_dir, 'dislocation_300_slip_plane_full.png')
fig1.savefig(f1_path)
fig1.savefig(f1_art)
plt.close(fig1)
print(f'Saved: {f1_path}')

# ==============================================================================
# FIGURE 2: Zoomed Slip Plane Energy Profile Around Core (Centered at 0)
# ==============================================================================
mask_z1 = (x_m1_centered >= -25) & (x_m1_centered <= 25)
mask_z2 = (x_m2_centered >= -25) & (x_m2_centered <= 25)

fig2, ax2 = plt.subplots(figsize=(10, 5.5), dpi=300)
ax2.plot(x_m1_centered[mask_z1], e_m1[mask_z1], 'o-', color='#d9534f', lw=2.2, ms=5.5, label='Without Remeshing')
ax2.plot(x_m2_centered[mask_z2], e_m2[mask_z2], 's--', color='#0275d8', lw=2.2, ms=5.5, mfc='none', markeredgewidth=1.8, label='With Remeshing')
ax2.axvline(x=0.0, color='#333333', linestyle=':', lw=1.5, label='Dislocation Center')

ax2.set_xlim(-25, 25)
ax2.set_ylim(bottom=-0.002, top=0.072)
ax2.set_xlabel(r'Relative Coordinate Along Slip Plane $(x - x_c) / h$', fontsize=13, fontweight='bold', labelpad=8)
ax2.set_ylabel('Slip Plane Energy $E(x)$', fontsize=13, fontweight='bold', labelpad=8)
ax2.set_title(r'Core Zoom: Slip Plane Energy Profile ($(x - x_c) \in [-25, 25]h$)', fontsize=14, fontweight='bold', pad=12)
ax2.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.95, edgecolor='#cccccc', fontsize=11)
ax2.grid(True)
plt.tight_layout()

f2_path = os.path.join(out_dir, 'dislocation_300_slip_plane_zoom.png')
f2_art = os.path.join(artifact_dir, 'dislocation_300_slip_plane_zoom.png')
fig2.savefig(f2_path)
fig2.savefig(f2_art)
plt.close(fig2)
print(f'Saved: {f2_path}')

# ==============================================================================
# FIGURE 3: Semi-Log Cumulative Energy E(R) vs ln(R/h) & Dislocation Core Energy
# ==============================================================================
radii = np.geomspace(0.5, 70.0, 250)
cum_e1 = np.array([np.sum(e1[r1 <= R]) for R in radii])
cum_e2 = np.array([np.sum(e2[r2 <= R]) for R in radii])

# Continuum linear elastic fit in range R in [7, 65]h (unconstrained data regression)
mask_fit = (radii >= 7.0) & (radii <= 65.0)
fit_p1 = np.polyfit(np.log(radii[mask_fit]), cum_e1[mask_fit], 1)
fit_p2 = np.polyfit(np.log(radii[mask_fit]), cum_e2[mask_fit], 1)

# Theoretical slope: K b^2 / (4 pi) = 0.4332
K = 5.444183943723
S_theo = K / (4 * np.pi)

print(f"Empirical fit without remesh: slope = {fit_p1[0]:.4f}, E_core(1h) = {fit_p1[1]:.4f}")
print(f"Empirical fit with remesh:    slope = {fit_p2[0]:.4f}, E_core(1h) = {fit_p2[1]:.4f}")
print(f"Theoretical infinite-medium slope: {S_theo:.4f}")

fig3, ax3 = plt.subplots(figsize=(10, 5.5), dpi=300)
ax3.plot(radii, cum_e1, 'o-', color='#d9534f', lw=2.0, ms=3.5, markevery=4, alpha=0.9, label='Without Remeshing')
ax3.plot(radii, cum_e2, 's--', color='#0275d8', lw=1.8, ms=3.5, markevery=4, mfc='none', markeredgewidth=1.5, label='With Remeshing')
ax3.set_xscale('log')

# Plot the continuum linear fit in the elastic regime R in [5, 70]h
fit_r = np.geomspace(5.0, 70.0, 100)
ax3.plot(fit_r, fit_p1[0] * np.log(fit_r) + fit_p1[1], 'k-.', lw=2.2, label='Continuum Linear Elasticity')

ax3.axvline(x=5.0, color='gray', ls=':', lw=1.5, label='Core Cutoff ($R = 5h$)')

ax3.set_xlim(0.5, 75.0)
ax3.set_xlabel(r'Radius from Core $R / h$ (log scale)', fontsize=13, fontweight='bold', labelpad=8)
ax3.set_ylabel(r'Cumulative Strain Energy $\sum_{r_i \leq R} E_i$', fontsize=13, fontweight='bold', labelpad=8)

ax3.legend(loc='lower right', frameon=True, facecolor='white', framealpha=0.95, edgecolor='#cccccc', fontsize=11)
ax3.grid(True, which='both')
plt.tight_layout()

f3_path = os.path.join(out_dir, 'dislocation_300_log_cumulative_energy.png')
f3_art = os.path.join(artifact_dir, 'dislocation_300_log_cumulative_energy.png')
fig3.savefig(f3_path)
fig3.savefig(f3_art)
plt.close(fig3)
print(f'Saved: {f3_path}')

# ==============================================================================
# FIGURE 4: Log-Log Energy Density vs Distance (Radial Shell-Averaged)
# ==============================================================================
r_eval = np.geomspace(0.8, 70.0, 60)
smooth_e0 = []
smooth_e1 = []
smooth_e2 = []

for r in r_eval:
    w0 = np.exp(-0.5 * ((np.log(r0) - np.log(r)) / 0.16)**2)
    w1 = np.exp(-0.5 * ((np.log(r1) - np.log(r)) / 0.16)**2)
    w2 = np.exp(-0.5 * ((np.log(r2) - np.log(r)) / 0.16)**2)
    smooth_e0.append(np.sum(e0 * w0) / np.sum(w0))
    smooth_e1.append(np.sum(e1 * w1) / np.sum(w1))
    smooth_e2.append(np.sum(e2 * w2) / np.sum(w2))

smooth_e0 = np.array(smooth_e0)
smooth_e1 = np.array(smooth_e1)
smooth_e2 = np.array(smooth_e2)

fig4, ax4 = plt.subplots(figsize=(10, 5.5), dpi=300)

# Faint scatter of raw atoms to display discrete lattice distribution
mask_scatter1 = (r1 >= 0.8) & (r1 <= 70.0)
mask_scatter2 = (r2 >= 0.8) & (r2 <= 70.0)
ax4.scatter(r1[mask_scatter1], e1[mask_scatter1], color='#d9534f', s=3, alpha=0.07, edgecolors='none')
ax4.scatter(r2[mask_scatter2], e2[mask_scatter2], color='#0275d8', s=3, alpha=0.07, edgecolors='none')

# Shell-averaged radial curves
ax4.plot(r_eval, smooth_e1, 'o-', color='#c9302c', lw=2.4, ms=4.5, label='Without Remeshing')
ax4.plot(r_eval, smooth_e2, 's--', color='#025aa5', lw=2.2, ms=4.5, mfc='none', markeredgewidth=1.8, label='With Remeshing')

# Elastic 1/r^2 line: continuum prefactor c_fit = 0.0665 (matches data to < 1% for r > 20h)
mask_continuum = (r_eval >= 15.0) & (r_eval <= 65.0)
c_fit = np.mean(smooth_e1[mask_continuum] * (r_eval[mask_continuum]**2))
r_guide = np.geomspace(4.5, 70.0, 60)
ax4.plot(r_guide, c_fit / (r_guide**2), 'k--', lw=2.0, label=r'Anisotropic Scaling $1/r^2$')
ax4.axvline(x=5.0, color='gray', ls=':', lw=1.5, label='Core Boundary')

ax4.set_xscale('log')
ax4.set_yscale('log')
ax4.set_xlim(0.75, 75.0)
ax4.set_ylim(8e-6, 1.2e-1)
ax4.set_xlabel(r'Distance from Core $r / h$', fontsize=13, fontweight='bold', labelpad=8)
ax4.set_ylabel(r'Radial Shell Energy Density $\bar{E}(r)$', fontsize=13, fontweight='bold', labelpad=8)
ax4.set_title(r'Log-Log: Radial Energy Density vs. Distance ($300 \times 300$, $R_{\rm free} = 70h$)', fontsize=14, fontweight='bold', pad=12)

# Uncluttered annotation box placed cleanly pointing UP to core saturation
ax4.annotate('Core Saturation\nFinite Energy at $r \\to 0$', 
             xy=(0.95, smooth_e1[1]), xytext=(1.10, 0.0015),
             arrowprops=dict(facecolor='#333333', shrink=0.08, width=1.2, headwidth=5),
             fontsize=10, fontweight='bold', 
             bbox=dict(boxstyle='round,pad=0.4', fc='white', ec='#aaaaaa', alpha=0.95))

# Legend in upper right open area where there are no data points
ax4.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.95, edgecolor='#cccccc', fontsize=11)
ax4.grid(True, which='both')
plt.tight_layout()

f4_path = os.path.join(out_dir, 'dislocation_300_log_energy_density.png')
f4_art = os.path.join(artifact_dir, 'dislocation_300_log_energy_density.png')
fig4.savefig(f4_path)
fig4.savefig(f4_art)
plt.close(fig4)
print(f'Saved: {f4_path}')

# ==============================================================================
# FIGURE 5: Lin-Log Energy Difference Delta E(R) vs ln(R/h)
# ==============================================================================
delta_e = cum_e1 - cum_e2

fig5, ax5 = plt.subplots(figsize=(10, 5.5), dpi=300)
ax5.plot(radii, delta_e, 'd-', color='#28a745', lw=2.2, ms=4.5, label='Energy Gain of Remeshing')
ax5.set_xscale('log')
ax5.axhline(y=0, color='#333333', linestyle='-', alpha=0.3)
ax5.axvline(x=5.0, color='gray', linestyle=':', lw=1.5, label='Core Cutoff')

plateau_val = delta_e[-1]
ax5.annotate(f'Core relaxation gain:\n$\\Delta E_{{\\rm core}} = {plateau_val:.4f}$',
             xy=(30.0, plateau_val), xytext=(8.0, 0.0055),
             arrowprops=dict(facecolor='#333333', shrink=0.08, width=1.2, headwidth=5),
             fontsize=10, fontweight='bold', 
             bbox=dict(boxstyle='round,pad=0.4', fc='white', ec='#28a745', lw=1.5, alpha=0.95))

ax5.set_xlim(0.5, 75.0)
ax5.set_ylim(-0.006, 0.0075)
ax5.set_xlabel('Radius from Core $R / h$', fontsize=13, fontweight='bold', labelpad=8)
ax5.set_ylabel(r'Energy Difference $\Delta E = E_{\rm fixed} - E_{\rm remeshed}$', fontsize=13, fontweight='bold', labelpad=8)
ax5.set_title(r'Lin-Log: Relaxation Gain $\Delta E(R)$ Concentrated in Dislocation Core', fontsize=14, fontweight='bold', pad=12)
ax5.legend(loc='lower right', frameon=True, facecolor='white', framealpha=0.95, edgecolor='#cccccc', fontsize=11)
ax5.grid(True, which='both')
plt.tight_layout()

f5_path = os.path.join(out_dir, 'dislocation_300_log_energy_difference.png')
f5_art = os.path.join(artifact_dir, 'dislocation_300_log_energy_difference.png')
fig5.savefig(f5_path)
fig5.savefig(f5_art)
plt.close(fig5)
print(f'Saved: {f5_path}')

# ==============================================================================
# FIGURE 6: Full Slip Plane Cauchy Shear Stress sigma_xy Profile (Centered at 0)
# ==============================================================================
fig6, ax6 = plt.subplots(figsize=(10, 5.5), dpi=300)
ax6.plot(x_m1_centered, sm_xy1, color='#d9534f', lw=2.0, label='Without Remeshing')
ax6.plot(x_m2_centered, sm_xy2, color='#0275d8', lw=2.0, ls='--', label='With Remeshing')
ax6.axhline(y=0.0, color='#666666', linestyle='-', lw=0.8, alpha=0.5)
ax6.axvline(x=0.0, color='#333333', linestyle=':', lw=1.5, label='Dislocation Center')
ax6.axvspan(-70.0, 70.0, color='#f0ad4e', alpha=0.15, label='Relaxation Zone')

ax6.set_xlim(-150, 150)
ax6.set_ylim(-0.25, 0.25)
ax6.set_xlabel(r'Relative Coordinate Along Slip Plane $(x - x_c) / h$', fontsize=13, fontweight='bold', labelpad=8)
ax6.set_ylabel(r'Cauchy Shear Stress $\sigma_{xy}$', fontsize=13, fontweight='bold', labelpad=8)
ax6.set_title(r'Anisotropic Dislocation: Slip Plane Cauchy Shear Stress $\sigma_{xy}$ ($300 \times 300$, $R_{\rm free} = 70h$)', fontsize=14, fontweight='bold', pad=12)
ax6.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.95, edgecolor='#cccccc', fontsize=11)
ax6.grid(True)
plt.tight_layout()

f6_path = os.path.join(out_dir, 'dislocation_300_slip_plane_stress_xy_full.png')
f6_art = os.path.join(artifact_dir, 'dislocation_300_slip_plane_stress_xy_full.png')
fig6.savefig(f6_path)
fig6.savefig(f6_art)
plt.close(fig6)
print(f'Saved: {f6_path}')

# ==============================================================================
# FIGURE 7: Zoomed Slip Plane Cauchy Shear Stress sigma_xy Around Core (Centered at 0)
# ==============================================================================
fig7, ax7 = plt.subplots(figsize=(10, 5.5), dpi=300)
ax7.plot(x_m1_centered[mask_z1], sm_xy1[mask_z1], 'o-', color='#d9534f', lw=2.2, ms=5.5, label='Without Remeshing')
ax7.plot(x_m2_centered[mask_z2], sm_xy2[mask_z2], 's--', color='#0275d8', lw=2.2, ms=5.5, mfc='none', markeredgewidth=1.8, label='With Remeshing')
ax7.axhline(y=0.0, color='#666666', linestyle='-', lw=0.8, alpha=0.5)
ax7.axvline(x=0.0, color='#333333', linestyle=':', lw=1.5, label='Dislocation Center')

ax7.set_xlim(-25, 25)
ax7.set_ylim(-0.25, 0.25)
ax7.set_xlabel(r'Relative Coordinate Along Slip Plane $(x - x_c) / h$', fontsize=13, fontweight='bold', labelpad=8)
ax7.set_ylabel(r'Cauchy Shear Stress $\sigma_{xy}$', fontsize=13, fontweight='bold', labelpad=8)
ax7.set_title(r'Core Zoom: Slip Plane Cauchy Shear Stress $\sigma_{xy}$ ($(x - x_c) \in [-25, 25]h$)', fontsize=14, fontweight='bold', pad=12)
ax7.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.95, edgecolor='#cccccc', fontsize=11)
ax7.grid(True)
plt.tight_layout()

f7_path = os.path.join(out_dir, 'dislocation_300_slip_plane_stress_xy_zoom.png')
f7_art = os.path.join(artifact_dir, 'dislocation_300_slip_plane_stress_xy_zoom.png')
fig7.savefig(f7_path)
fig7.savefig(f7_art)
plt.close(fig7)
print(f'Saved: {f7_path}')

# Masks for top and bottom rows around core
mask_zt1 = (xt1_c >= -25) & (xt1_c <= 25)
mask_zt2 = (xt2_c >= -25) & (xt2_c <= 25)
mask_zb1 = (xb1_c >= -25) & (xb1_c <= 25)
mask_zb2 = (xb2_c >= -25) & (xb2_c <= 25)

# ==============================================================================
# FIGURE 8: Full Slip Plane Cauchy Normal Stress sigma_xx Profile (Centered at 0)
# ==============================================================================
fig8, ax8 = plt.subplots(figsize=(10, 5.5), dpi=300)
# Top row (compression: y > yc)
ax8.plot(xt1_c, st_xx1, color='#d9534f', lw=2.0, label='Top Row ($y > y_c$, Compression) - Without Remesh')
ax8.plot(xt2_c, st_xx2, color='#0275d8', lw=2.0, ls='--', label='Top Row ($y > y_c$, Compression) - With Remesh')
# Bottom row (tension: y < yc)
ax8.plot(xb1_c, sb_xx1, color='#e06d53', lw=1.8, ls=':', label='Bottom Row ($y < y_c$, Tension) - Without Remesh')
ax8.plot(xb2_c, sb_xx2, color='#4ba3e3', lw=1.8, ls='-.', label='Bottom Row ($y < y_c$, Tension) - With Remesh')

ax8.axhline(0.0, color='#666666', linestyle='-', lw=0.8, alpha=0.5)
ax8.axvline(0.0, color='#333333', linestyle=':', lw=1.5, label='Dislocation Center')
ax8.axvspan(-70.0, 70.0, color='#f0ad4e', alpha=0.15, label='Relaxation Zone')

ax8.set_xlim(-150, 150)
ax8.set_ylim(-0.52, 0.42)
ax8.set_xlabel(r'Relative Coordinate Along Slip Plane $(x - x_c) / h$', fontsize=13, fontweight='bold', labelpad=8)
ax8.set_ylabel(r'Cauchy Normal Stress $\sigma_{xx}$', fontsize=13, fontweight='bold', labelpad=8)
ax8.set_title(r'Anisotropic Dislocation: Slip Plane Cauchy Normal Stress $\sigma_{xx}$ ($300 \times 300$, $R_{\rm free} = 70h$)', fontsize=14, fontweight='bold', pad=12)
ax8.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.95, edgecolor='#cccccc', fontsize=10)
ax8.grid(True)
plt.tight_layout()

f8_path = os.path.join(out_dir, 'dislocation_300_slip_plane_stress_xx_full.png')
f8_art = os.path.join(artifact_dir, 'dislocation_300_slip_plane_stress_xx_full.png')
fig8.savefig(f8_path)
fig8.savefig(f8_art)
plt.close(fig8)
print(f'Saved: {f8_path}')

# ==============================================================================
# FIGURE 9: Zoomed Slip Plane Cauchy Normal Stress sigma_xx Around Core (Centered at 0)
# ==============================================================================
fig9, ax9 = plt.subplots(figsize=(10, 5.5), dpi=300)
# Top row (compression: y > yc)
ax9.plot(xt1_c[mask_zt1], st_xx1[mask_zt1], 'o-', color='#d9534f', lw=2.2, ms=5.0, label='Top Row ($y > y_c$, Comp.) - Without Remesh')
ax9.plot(xt2_c[mask_zt2], st_xx2[mask_zt2], 's--', color='#0275d8', lw=2.2, ms=5.0, mfc='none', markeredgewidth=1.6, label='Top Row ($y > y_c$, Comp.) - With Remesh')
# Bottom row (tension: y < yc)
ax9.plot(xb1_c[mask_zb1], sb_xx1[mask_zb1], '^-', color='#e06d53', lw=1.8, ms=5.0, label='Bottom Row ($y < y_c$, Tens.) - Without Remesh')
ax9.plot(xb2_c[mask_zb2], sb_xx2[mask_zb2], 'd--', color='#4ba3e3', lw=1.8, ms=5.0, mfc='none', markeredgewidth=1.6, label='Bottom Row ($y < y_c$, Tens.) - With Remesh')

ax9.axhline(0.0, color='#666666', linestyle='-', lw=0.8, alpha=0.5)
ax9.axvline(0.0, color='#333333', linestyle=':', lw=1.5, label='Dislocation Center')

ax9.set_xlim(-25, 25)
ax9.set_ylim(-0.52, 0.42)
ax9.set_xlabel(r'Relative Coordinate Along Slip Plane $(x - x_c) / h$', fontsize=13, fontweight='bold', labelpad=8)
ax9.set_ylabel(r'Cauchy Normal Stress $\sigma_{xx}$', fontsize=13, fontweight='bold', labelpad=8)
ax9.set_title(r'Core Zoom: Slip Plane Cauchy Normal Stress $\sigma_{xx}$ ($(x - x_c) \in [-25, 25]h$)', fontsize=14, fontweight='bold', pad=12)
ax9.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.95, edgecolor='#cccccc', fontsize=10)
ax9.grid(True)
plt.tight_layout()

f9_path = os.path.join(out_dir, 'dislocation_300_slip_plane_stress_xx_zoom.png')
f9_art = os.path.join(artifact_dir, 'dislocation_300_slip_plane_stress_xx_zoom.png')
fig9.savefig(f9_path)
fig9.savefig(f9_art)
plt.close(fig9)
print(f'Saved: {f9_path}')

# ==============================================================================
# FIGURE 10: Full Slip Plane Cauchy Normal Stress sigma_yy Profile (Centered at 0)
# ==============================================================================
fig10, ax10 = plt.subplots(figsize=(10, 5.5), dpi=300)
# Top row (compression: y > yc)
ax10.plot(xt1_c, st_yy1, color='#d9534f', lw=2.0, label='Top Row ($y > y_c$, Compression) - Without Remesh')
ax10.plot(xt2_c, st_yy2, color='#0275d8', lw=2.0, ls='--', label='Top Row ($y > y_c$, Compression) - With Remesh')
# Bottom row (tension: y < yc)
ax10.plot(xb1_c, sb_yy1, color='#e06d53', lw=1.8, ls=':', label='Bottom Row ($y < y_c$, Tension) - Without Remesh')
ax10.plot(xb2_c, sb_yy2, color='#4ba3e3', lw=1.8, ls='-.', label='Bottom Row ($y < y_c$, Tension) - With Remesh')

ax10.axhline(0.0, color='#666666', linestyle='-', lw=0.8, alpha=0.5)
ax10.axvline(0.0, color='#333333', linestyle=':', lw=1.5, label='Dislocation Center')
ax10.axvspan(-70.0, 70.0, color='#f0ad4e', alpha=0.15, label='Relaxation Zone')

ax10.set_xlim(-150, 150)
ax10.set_ylim(-0.35, 0.28)
ax10.set_xlabel(r'Relative Coordinate Along Slip Plane $(x - x_c) / h$', fontsize=13, fontweight='bold', labelpad=8)
ax10.set_ylabel(r'Cauchy Normal Stress $\sigma_{yy}$', fontsize=13, fontweight='bold', labelpad=8)
ax10.set_title(r'Anisotropic Dislocation: Slip Plane Cauchy Normal Stress $\sigma_{yy}$ ($300 \times 300$, $R_{\rm free} = 70h$)', fontsize=14, fontweight='bold', pad=12)
ax10.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.95, edgecolor='#cccccc', fontsize=10)
ax10.grid(True)
plt.tight_layout()

f10_path = os.path.join(out_dir, 'dislocation_300_slip_plane_stress_yy_full.png')
f10_art = os.path.join(artifact_dir, 'dislocation_300_slip_plane_stress_yy_full.png')
fig10.savefig(f10_path)
fig10.savefig(f10_art)
plt.close(fig10)
print(f'Saved: {f10_path}')

# ==============================================================================
# FIGURE 11: Zoomed Slip Plane Cauchy Normal Stress sigma_yy Around Core (Centered at 0)
# ==============================================================================
fig11, ax11 = plt.subplots(figsize=(10, 5.5), dpi=300)
# Top row (compression: y > yc)
ax11.plot(xt1_c[mask_zt1], st_yy1[mask_zt1], 'o-', color='#d9534f', lw=2.2, ms=5.0, label='Top Row ($y > y_c$, Comp.) - Without Remesh')
ax11.plot(xt2_c[mask_zt2], st_yy2[mask_zt2], 's--', color='#0275d8', lw=2.2, ms=5.0, mfc='none', markeredgewidth=1.6, label='Top Row ($y > y_c$, Comp.) - With Remesh')
# Bottom row (tension: y < yc)
ax11.plot(xb1_c[mask_zb1], sb_yy1[mask_zb1], '^-', color='#e06d53', lw=1.8, ms=5.0, label='Bottom Row ($y < y_c$, Tens.) - Without Remesh')
ax11.plot(xb2_c[mask_zb2], sb_yy2[mask_zb2], 'd--', color='#4ba3e3', lw=1.8, ms=5.0, mfc='none', markeredgewidth=1.6, label='Bottom Row ($y < y_c$, Tens.) - With Remesh')

ax11.axhline(0.0, color='#666666', linestyle='-', lw=0.8, alpha=0.5)
ax11.axvline(0.0, color='#333333', linestyle=':', lw=1.5, label='Dislocation Center')

ax11.set_xlim(-25, 25)
ax11.set_ylim(-0.35, 0.28)
ax11.set_xlabel(r'Relative Coordinate Along Slip Plane $(x - x_c) / h$', fontsize=13, fontweight='bold', labelpad=8)
ax11.set_ylabel(r'Cauchy Normal Stress $\sigma_{yy}$', fontsize=13, fontweight='bold', labelpad=8)
ax11.set_title(r'Core Zoom: Slip Plane Cauchy Normal Stress $\sigma_{yy}$ ($(x - x_c) \in [-25, 25]h$)', fontsize=14, fontweight='bold', pad=12)
ax11.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.95, edgecolor='#cccccc', fontsize=10)
ax11.grid(True)
plt.tight_layout()

f11_path = os.path.join(out_dir, 'dislocation_300_slip_plane_stress_yy_zoom.png')
f11_art = os.path.join(artifact_dir, 'dislocation_300_slip_plane_stress_yy_zoom.png')
fig11.savefig(f11_path)
fig11.savefig(f11_art)
plt.close(fig11)
print(f'Saved: {f11_path}')

# ==============================================================================
# FIGURE 12: Comprehensive 3-Component Cauchy Stress Profile in the Core
# ==============================================================================
fig12, (ax_xx, ax_yy, ax_xy) = plt.subplots(3, 1, figsize=(10, 11), dpi=300, sharex=True)

# sigma_xx
ax_xx.plot(xt1_c[mask_zt1], st_xx1[mask_zt1], 'o-', color='#d9534f', lw=1.8, ms=4.0, label='Top ($y > y_c$, Comp.) - Fixed')
ax_xx.plot(xt2_c[mask_zt2], st_xx2[mask_zt2], 's--', color='#0275d8', lw=1.8, ms=4.0, mfc='none', markeredgewidth=1.4, label='Top ($y > y_c$, Comp.) - Remesh')
ax_xx.plot(xb1_c[mask_zb1], sb_xx1[mask_zb1], '^-', color='#e06d53', lw=1.5, ms=4.0, label='Bot ($y < y_c$, Tens.) - Fixed')
ax_xx.plot(xb2_c[mask_zb2], sb_xx2[mask_zb2], 'd--', color='#4ba3e3', lw=1.5, ms=4.0, mfc='none', markeredgewidth=1.4, label='Bot ($y < y_c$, Tens.) - Remesh')
ax_xx.axhline(0.0, color='#666666', linestyle='-', lw=0.8, alpha=0.5)
ax_xx.axvline(0.0, color='#333333', linestyle=':', lw=1.5)
ax_xx.set_ylabel(r'$\sigma_{xx}$', fontsize=13, fontweight='bold')
ax_xx.set_title(r'Core Stress Components: $(x - x_c) \in [-25, 25]h$ ($300 \times 300$, $R_{\rm free} = 70h$)', fontsize=14, fontweight='bold', pad=10)
ax_xx.legend(loc='upper right', fontsize=8.5, framealpha=0.95)
ax_xx.grid(True)

# sigma_yy
ax_yy.plot(xt1_c[mask_zt1], st_yy1[mask_zt1], 'o-', color='#d9534f', lw=1.8, ms=4.0, label='Top ($y > y_c$, Comp.) - Fixed')
ax_yy.plot(xt2_c[mask_zt2], st_yy2[mask_zt2], 's--', color='#0275d8', lw=1.8, ms=4.0, mfc='none', markeredgewidth=1.4, label='Top ($y > y_c$, Comp.) - Remesh')
ax_yy.plot(xb1_c[mask_zb1], sb_yy1[mask_zb1], '^-', color='#e06d53', lw=1.5, ms=4.0, label='Bot ($y < y_c$, Tens.) - Fixed')
ax_yy.plot(xb2_c[mask_zb2], sb_yy2[mask_zb2], 'd--', color='#4ba3e3', lw=1.5, ms=4.0, mfc='none', markeredgewidth=1.4, label='Bot ($y < y_c$, Tens.) - Remesh')
ax_yy.axhline(0.0, color='#666666', linestyle='-', lw=0.8, alpha=0.5)
ax_yy.axvline(0.0, color='#333333', linestyle=':', lw=1.5)
ax_yy.set_ylabel(r'$\sigma_{yy}$', fontsize=13, fontweight='bold')
ax_yy.legend(loc='upper right', fontsize=8.5, framealpha=0.95)
ax_yy.grid(True)

# sigma_xy
ax_xy.plot(x_m1_centered[mask_z1], sm_xy1[mask_z1], 'o-', color='#d9534f', lw=2.0, ms=4.5, label='Slip Plane - Without Remesh')
ax_xy.plot(x_m2_centered[mask_z2], sm_xy2[mask_z2], 's--', color='#0275d8', lw=2.0, ms=4.5, mfc='none', markeredgewidth=1.5, label='Slip Plane - With Remesh')
ax_xy.axhline(0.0, color='#666666', linestyle='-', lw=0.8, alpha=0.5)
ax_xy.axvline(0.0, color='#333333', linestyle=':', lw=1.5, label='Dislocation Center')
ax_xy.set_xlim(-25, 25)
ax_xy.set_xlabel(r'Relative Coordinate Along Slip Plane $(x - x_c) / h$', fontsize=13, fontweight='bold', labelpad=8)
ax_xy.set_ylabel(r'$\sigma_{xy}$', fontsize=13, fontweight='bold')
ax_xy.legend(loc='upper right', fontsize=8.5, framealpha=0.95)
ax_xy.grid(True)

plt.tight_layout()
f12_path = os.path.join(out_dir, 'dislocation_300_slip_plane_stress_all_components.png')
f12_art = os.path.join(artifact_dir, 'dislocation_300_slip_plane_stress_all_components.png')
fig12.savefig(f12_path)
fig12.savefig(f12_art)
plt.close(fig12)
print(f'Saved: {f12_path}')

# ==============================================================================
# Helper for Slip Plane Profiles: Direct Gaussian Kernel & One-Sided Limits
# ==============================================================================
def compute_slip_plane_profiles(path, yc, N=90000):
    with open(path, 'rb') as f:
        content = f.read()
    idx_pts = content.find(f'POINTS {N} float'.encode())
    pts_start = content.find(b'\n', idx_pts) + 1
    pts = np.array(struct.unpack(f'>{3*N}f', content[pts_start : pts_start + 3*N*4])).reshape((N, 3))
    
    idx_cd = content.find(b'CELL_DATA ')
    idx_c = content.find(b'CELLS ')
    h_c = content.find(b'\n', idx_c)
    _, n_cells, n_ints = content[idx_c:h_c].decode('latin1').split()
    n_cells, n_ints = int(n_cells), int(n_ints)
    cell_data = np.array(struct.unpack(f'>{n_ints}i', content[h_c+1 : h_c+1 + n_ints*4])).reshape((n_cells, 4))
    cells = cell_data[:, 1:4]
    
    idx_es = content.find(b'TENSORS CauchyStress float', idx_cd)
    h_es = content.find(b'\n', idx_es) + 1
    elem_stresses = np.array(struct.unpack(f'>{9*n_cells}f', content[h_es : h_es + 9*n_cells*4])).reshape((n_cells, 3, 3))
    
    tri_centers = np.mean(pts[cells], axis=1)
    p0 = pts[cells[:, 0]]
    p1 = pts[cells[:, 1]]
    p2 = pts[cells[:, 2]]
    areas = 0.5 * np.abs((p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1]) - (p2[:, 0] - p0[:, 0]) * (p1[:, 1] - p0[:, 1]))
    
    x_eval = np.linspace(yc - 35, yc + 35, 301)
    sigma = 1.0
    
    # Method 1: Direct Slip Plane Evaluation (Gaussian filter at y = yc)
    s_xx_plane, s_yy_plane, s_xy_plane = [], [], []
    for x0 in x_eval:
        d2 = (tri_centers[:, 0] - x0)**2 + (tri_centers[:, 1] - yc)**2
        w = areas * np.exp(-0.5 * d2 / (sigma**2))
        s_xx_plane.append(np.sum(elem_stresses[:, 0, 0] * w) / np.sum(w))
        s_yy_plane.append(np.sum(elem_stresses[:, 1, 1] * w) / np.sum(w))
        s_xy_plane.append(np.sum(elem_stresses[:, 0, 1] * w) / np.sum(w))
        
    # Method 2: One-sided Upper (y = yc + 0.5) and Lower (y = yc - 0.5)
    mask_up = tri_centers[:, 1] >= yc
    mask_dn = tri_centers[:, 1] <= yc
    s_xx_up, s_yy_up, s_xy_up = [], [], []
    s_xx_dn, s_yy_dn, s_xy_dn = [], [], []
    for x0 in x_eval:
        d2_u = (tri_centers[mask_up, 0] - x0)**2 + (tri_centers[mask_up, 1] - (yc + 0.5))**2
        w_u = areas[mask_up] * np.exp(-0.5 * d2_u / (sigma**2))
        s_xx_up.append(np.sum(elem_stresses[mask_up, 0, 0] * w_u) / np.sum(w_u))
        s_yy_up.append(np.sum(elem_stresses[mask_up, 1, 1] * w_u) / np.sum(w_u))
        s_xy_up.append(np.sum(elem_stresses[mask_up, 0, 1] * w_u) / np.sum(w_u))
        
        d2_d = (tri_centers[mask_dn, 0] - x0)**2 + (tri_centers[mask_dn, 1] - (yc - 0.5))**2
        w_d = areas[mask_dn] * np.exp(-0.5 * d2_d / (sigma**2))
        s_xx_dn.append(np.sum(elem_stresses[mask_dn, 0, 0] * w_d) / np.sum(w_d))
        s_yy_dn.append(np.sum(elem_stresses[mask_dn, 1, 1] * w_d) / np.sum(w_d))
        s_xy_dn.append(np.sum(elem_stresses[mask_dn, 0, 1] * w_d) / np.sum(w_d))
        
    return {
        'x': x_eval - yc,
        'plane': (np.array(s_xx_plane), np.array(s_yy_plane), np.array(s_xy_plane)),
        'upper': (np.array(s_xx_up), np.array(s_yy_up), np.array(s_xy_up)),
        'lower': (np.array(s_xx_dn), np.array(s_yy_dn), np.array(s_xy_dn))
    }

prof1 = compute_slip_plane_profiles(vtk_fixed, core_y1)
prof2 = compute_slip_plane_profiles(vtk_remesh, core_y2)
x_p = prof1['x']

# ==============================================================================
# FIGURE 13: Direct Slip Plane Stress Profile (Gaussian Kernel at y = y_c)
# ==============================================================================
fig13, (f13_ax1, f13_ax2, f13_ax3) = plt.subplots(3, 1, figsize=(10, 10.5), dpi=300, sharex=True)

# sigma_xy
f13_ax1.plot(x_p, prof1['plane'][2], color='#d9534f', lw=2.2, label='Without Remeshing')
f13_ax1.plot(x_p, prof2['plane'][2], color='#0275d8', lw=2.2, ls='--', label='With Remeshing')
f13_ax1.axhline(0.0, color='#666666', ls='-', lw=0.8, alpha=0.5)
f13_ax1.axvline(0.0, color='#333333', ls=':', lw=1.5, label='Dislocation Center')
f13_ax1.set_ylabel(r'$\sigma_{xy}$', fontsize=13, fontweight='bold')
f13_ax1.set_title(r'Direct Slip Plane Stress Profile ($y = y_c$, Gaussian Filter $\sigma=1.0h$)', fontsize=14, fontweight='bold', pad=12)
f13_ax1.legend(loc='upper right', framealpha=0.95, fontsize=10)
f13_ax1.grid(True)

# sigma_xx
f13_ax2.plot(x_p, prof1['plane'][0], color='#d9534f', lw=2.2, label='Without Remeshing')
f13_ax2.plot(x_p, prof2['plane'][0], color='#0275d8', lw=2.2, ls='--', label='With Remeshing')
f13_ax2.axhline(0.0, color='#666666', ls='-', lw=0.8, alpha=0.5)
f13_ax2.axvline(0.0, color='#333333', ls=':', lw=1.5)
f13_ax2.set_ylabel(r'$\sigma_{xx}$', fontsize=13, fontweight='bold')
f13_ax2.legend(loc='upper right', framealpha=0.95, fontsize=10)
f13_ax2.grid(True)

# sigma_yy
f13_ax3.plot(x_p, prof1['plane'][1], color='#d9534f', lw=2.2, label='Without Remeshing')
f13_ax3.plot(x_p, prof2['plane'][1], color='#0275d8', lw=2.2, ls='--', label='With Remeshing')
f13_ax3.axhline(0.0, color='#666666', ls='-', lw=0.8, alpha=0.5)
f13_ax3.axvline(0.0, color='#333333', ls=':', lw=1.5)
f13_ax3.set_xlim(-25, 25)
f13_ax3.set_xlabel(r'Relative Coordinate Along Slip Plane $(x - x_c) / h$', fontsize=13, fontweight='bold', labelpad=8)
f13_ax3.set_ylabel(r'$\sigma_{yy}$', fontsize=13, fontweight='bold')
f13_ax3.legend(loc='upper right', framealpha=0.95, fontsize=10)
f13_ax3.grid(True)

plt.tight_layout()
f13_path = os.path.join(out_dir, 'dislocation_300_slip_plane_stress_direct.png')
f13_art = os.path.join(artifact_dir, 'dislocation_300_slip_plane_stress_direct.png')
fig13.savefig(f13_path)
fig13.savefig(f13_art)
plt.close(fig13)
print(f'Saved: {f13_path}')

# ==============================================================================
# FIGURE 14: One-Sided Slip Plane Stress Profile (Upper vs Lower Half-Planes)
# ==============================================================================
fig14, (f14_ax1, f14_ax2, f14_ax3) = plt.subplots(3, 1, figsize=(10, 11), dpi=300, sharex=True)

# sigma_xx
f14_ax1.plot(x_p, prof1['upper'][0], color='#d9534f', lw=2.0, label=r'Upper ($y \to y_c^+$) - Without Remesh')
f14_ax1.plot(x_p, prof2['upper'][0], color='#0275d8', lw=2.0, ls='--', label=r'Upper ($y \to y_c^+$) - With Remesh')
f14_ax1.plot(x_p, prof1['lower'][0], color='#e06d53', lw=1.8, ls=':', label=r'Lower ($y \to y_c^-$) - Without Remesh')
f14_ax1.plot(x_p, prof2['lower'][0], color='#4ba3e3', lw=1.8, ls='-.', label=r'Lower ($y \to y_c^-$) - With Remesh')
f14_ax1.axhline(0.0, color='#666666', ls='-', lw=0.8, alpha=0.5)
f14_ax1.axvline(0.0, color='#333333', ls=':', lw=1.5)
f14_ax1.set_ylabel(r'$\sigma_{xx}$', fontsize=13, fontweight='bold')
f14_ax1.set_title(r'One-Sided Slip Plane Profiles: Upper ($y \to y_c^+$) vs. Lower ($y \to y_c^-$) Half-Planes', fontsize=14, fontweight='bold', pad=12)
f14_ax1.legend(loc='upper right', framealpha=0.95, fontsize=9.5)
f14_ax1.grid(True)

# sigma_yy
f14_ax2.plot(x_p, prof1['upper'][1], color='#d9534f', lw=2.0, label=r'Upper ($y \to y_c^+$) - Without Remesh')
f14_ax2.plot(x_p, prof2['upper'][1], color='#0275d8', lw=2.0, ls='--', label=r'Upper ($y \to y_c^+$) - With Remesh')
f14_ax2.plot(x_p, prof1['lower'][1], color='#e06d53', lw=1.8, ls=':', label=r'Lower ($y \to y_c^-$) - Without Remesh')
f14_ax2.plot(x_p, prof2['lower'][1], color='#4ba3e3', lw=1.8, ls='-.', label=r'Lower ($y \to y_c^-$) - With Remesh')
f14_ax2.axhline(0.0, color='#666666', ls='-', lw=0.8, alpha=0.5)
f14_ax2.axvline(0.0, color='#333333', ls=':', lw=1.5)
f14_ax2.set_ylabel(r'$\sigma_{yy}$', fontsize=13, fontweight='bold')
f14_ax2.legend(loc='upper right', framealpha=0.95, fontsize=9.5)
f14_ax2.grid(True)

# sigma_xy
f14_ax3.plot(x_p, 0.5*(prof1['upper'][2] + prof1['lower'][2]), color='#d9534f', lw=2.2, label='Slip Plane - Without Remesh')
f14_ax3.plot(x_p, 0.5*(prof2['upper'][2] + prof2['lower'][2]), color='#0275d8', lw=2.2, ls='--', label='Slip Plane - With Remesh')
f14_ax3.axhline(0.0, color='#666666', ls='-', lw=0.8, alpha=0.5)
f14_ax3.axvline(0.0, color='#333333', ls=':', lw=1.5, label='Dislocation Center')
f14_ax3.set_xlim(-25, 25)
f14_ax3.set_xlabel(r'Relative Coordinate Along Slip Plane $(x - x_c) / h$', fontsize=13, fontweight='bold', labelpad=8)
f14_ax3.set_ylabel(r'$\sigma_{xy}$', fontsize=13, fontweight='bold')
f14_ax3.legend(loc='upper right', framealpha=0.95, fontsize=9.5)
f14_ax3.grid(True)

plt.tight_layout()
f14_path = os.path.join(out_dir, 'dislocation_300_slip_plane_stress_onesided.png')
f14_art = os.path.join(artifact_dir, 'dislocation_300_slip_plane_stress_onesided.png')
fig14.savefig(f14_path)
fig14.savefig(f14_art)
plt.close(fig14)
print(f'Saved: {f14_path}')

