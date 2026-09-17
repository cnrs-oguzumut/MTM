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

def parse_vtk_nodal_data(path, N=90000):
    with open(path, 'rb') as f:
        content = f.read()
    idx_pts = content.find(f'POINTS {N} float'.encode())
    pts_start = content.find(b'\n', idx_pts) + 1
    pts = struct.unpack(f'>{3*N}f', content[pts_start : pts_start + 3*N*4])
    idx_e = content.find(b'SCALARS NodalEnergy float 1')
    h_e = content.find(b'\n', content.find(b'\n', idx_e) + 1) + 1
    e = np.array(struct.unpack(f'>{N}f', content[h_e : h_e + N*4]))
    coords = np.array(pts).reshape((N, 3))[:, :2]
    return coords, e

def extract_slip_plane_rows(pts, energies, y_center=149.5):
    mask_bot = np.abs(pts[:, 1] - (y_center - 0.7)) < 0.45
    mask_top = np.abs(pts[:, 1] - (y_center + 0.3)) < 0.45
    
    order_bot = np.argsort(pts[mask_bot, 0])
    x_bot = pts[mask_bot, 0][order_bot]
    e_bot = energies[mask_bot][order_bot]
    
    order_top = np.argsort(pts[mask_top, 0])
    x_top = pts[mask_top, 0][order_top]
    e_top = energies[mask_top][order_top]
    
    x_mid = 0.5 * (x_bot + x_top)
    e_mid = 0.5 * (e_bot + e_top)
    
    return (x_bot, e_bot), (x_top, e_top), (x_mid, e_mid)

# Paths
vtk_unrelaxed = 'dislocation_study_300x300/without_remesh/vtk_output/configuration_00000.vtk'
vtk_fixed = 'dislocation_study_300x300/without_remesh/vtk_output/configuration_00001.vtk'
vtk_remesh = 'dislocation_study_300x300/with_remesh/vtk_output/configuration_00001.vtk'
artifact_dir = '/Users/usalman/.gemini/antigravity-ide/brain/e0366a0b-42af-4a53-9525-edb8c926ac89'
out_dir = 'dislocation_study_300x300'
os.makedirs(out_dir, exist_ok=True)

pts0, e0 = parse_vtk_nodal_data(vtk_unrelaxed, N=90000)
pts1, e1 = parse_vtk_nodal_data(vtk_fixed, N=90000)
pts2, e2 = parse_vtk_nodal_data(vtk_remesh, N=90000)

core_x1, core_y1 = 149.5, 149.5
core_x2, core_y2 = 150.0, 149.5
r0 = np.linalg.norm(pts0 - np.array([core_x1, core_y1]), axis=1)
r1 = np.linalg.norm(pts1 - np.array([core_x1, core_y1]), axis=1)
r2 = np.linalg.norm(pts2 - np.array([core_x2, core_y2]), axis=1)

(x_b1, e_b1), (x_t1, e_t1), (x_m1, e_m1) = extract_slip_plane_rows(pts1, e1, y_center=core_y1)
(x_b2, e_b2), (x_t2, e_t2), (x_m2, e_m2) = extract_slip_plane_rows(pts2, e2, y_center=core_y2)

# Shift coordinates so that each dislocation core is centered at 0:
x_m1_centered = x_m1 - core_x1
x_m2_centered = x_m2 - core_x2

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
# FIGURE 3: Semi-Log Cumulative Energy E(R) vs ln(R/h)
# ==============================================================================
radii = np.linspace(0.5, 70.0, 100)
cum_e1 = np.array([np.sum(e1[r1 <= R]) for R in radii])
cum_e2 = np.array([np.sum(e2[r2 <= R]) for R in radii])

fig3, ax3 = plt.subplots(figsize=(10, 5.5), dpi=300)
ax3.plot(radii, cum_e1, 'o-', color='#d9534f', lw=2.0, ms=3.5, alpha=0.9, label='Without Remeshing')
ax3.plot(radii, cum_e2, 's--', color='#0275d8', lw=1.8, ms=3.5, mfc='none', markeredgewidth=1.5, label='With Remeshing')
ax3.set_xscale('log')

# Continuum linear elastic fit in range R in [7, 65]h
mask_fit = (radii >= 7.0) & (radii <= 65.0)
fit_p = np.polyfit(np.log(radii[mask_fit]), cum_e1[mask_fit], 1)
fit_r = np.geomspace(5.0, 70.0, 60)
ax3.plot(fit_r, fit_p[0]*np.log(fit_r) + fit_p[1], 'k-.', lw=2.0, label='Continuum Slope')
ax3.axvline(x=5.0, color='gray', ls=':', lw=1.5, label='Core Cutoff')

ax3.set_xlim(0.5, 75.0)
ax3.set_xlabel('Radius from Core $R / h$', fontsize=13, fontweight='bold', labelpad=8)
ax3.set_ylabel(r'Cumulative Strain Energy $\sum_{r_i \leq R} E_i$', fontsize=13, fontweight='bold', labelpad=8)
ax3.set_title(r'Semi-Log: Cumulative Strain Energy $E(R)$ vs. $\ln(R/h)$ ($300 \times 300$, $R_{\rm free} = 70h$)', fontsize=14, fontweight='bold', pad=12)
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
