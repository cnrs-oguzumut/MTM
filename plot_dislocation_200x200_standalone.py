import os
import struct
import numpy as np
import matplotlib.pyplot as plt

# Set aesthetic styling
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['grid.color'] = '#e0e0e0'
plt.rcParams['grid.linestyle'] = '--'
plt.rcParams['grid.alpha'] = 0.6

def parse_vtk_nodal_data(path, N=40000):
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

def extract_slip_plane_rows(pts, energies, y_center=99.5):
    # Slip plane is at y = 99.5.
    # Bottom row is at y ≈ 98.8, top row is at y ≈ 99.8
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
vtk_fixed = 'dislocation_study_200x200/without_remesh/vtk_output/configuration_00001.vtk'
vtk_remesh = 'dislocation_study_200x200/with_remesh/vtk_output/configuration_00001.vtk'
artifact_dir = '/Users/usalman/.gemini/antigravity-ide/brain/e0366a0b-42af-4a53-9525-edb8c926ac89'

pts1, e1 = parse_vtk_nodal_data(vtk_fixed)
pts2, e2 = parse_vtk_nodal_data(vtk_remesh)

(x_b1, e_b1), (x_t1, e_t1), (x_m1, e_m1) = extract_slip_plane_rows(pts1, e1)
(x_b2, e_b2), (x_t2, e_t2), (x_m2, e_m2) = extract_slip_plane_rows(pts2, e2)

core_x, core_y = 99.5, 99.5
r1 = np.linalg.norm(pts1 - np.array([core_x, core_y]), axis=1)
r2 = np.linalg.norm(pts2 - np.array([core_x, core_y]), axis=1)
r_ref = 0.5 * (r1 + r2)

out_dir = 'dislocation_study_200x200'

# ==============================================================================
# FIGURE 1: Full Slip Plane Energy Profile (x in [0, 200])
# ==============================================================================
fig1, ax1 = plt.subplots(figsize=(10, 5.5), dpi=300)
ax1.plot(x_m1, e_m1, color='#d9534f', lw=2.0, label='Without Remeshing')
ax1.plot(x_m2, e_m2, color='#0275d8', lw=2.0, ls='--', label='With Remeshing')
ax1.axvline(x=core_x, color='#333333', linestyle=':', lw=1.5, label='Dislocation Center')
ax1.axvspan(core_x - 30.0, core_x + 30.0, color='#f0ad4e', alpha=0.15, label='Relaxation Zone')

ax1.set_xlim(0, 200)
ax1.set_ylim(bottom=-0.002, top=0.075)
ax1.set_xlabel('Atomic Position Along Slip Plane $x / h$', fontsize=13, fontweight='bold', labelpad=8)
ax1.set_ylabel('Slip Plane Energy $E(x)$', fontsize=13, fontweight='bold', labelpad=8)
ax1.set_title(r'Dislocation Energy Profile Along Slip Plane ($200 \times 200$ Lattice)', fontsize=14, fontweight='bold', pad=12)
ax1.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.95, edgecolor='#cccccc', fontsize=11)
ax1.grid(True)
plt.tight_layout()

f1_path = os.path.join(out_dir, 'dislocation_200_slip_plane_full.png')
f1_art = os.path.join(artifact_dir, 'dislocation_200_slip_plane_full.png')
fig1.savefig(f1_path)
fig1.savefig(f1_art)
plt.close(fig1)
print(f'Saved: {f1_path}')

# ==============================================================================
# FIGURE 2: Zoomed Slip Plane Energy Profile Around Core (x in [80, 120])
# ==============================================================================
mask_z1 = (x_m1 >= 80) & (x_m1 <= 120)
mask_z2 = (x_m2 >= 80) & (x_m2 <= 120)

fig2, ax2 = plt.subplots(figsize=(10, 5.5), dpi=300)
ax2.plot(x_m1[mask_z1], e_m1[mask_z1], 'o-', color='#d9534f', lw=2.2, ms=5.5, label='Without Remeshing')
ax2.plot(x_m2[mask_z2], e_m2[mask_z2], 's--', color='#0275d8', lw=2.2, ms=5.5, mfc='none', markeredgewidth=1.8, label='With Remeshing')
ax2.axvline(x=core_x, color='#333333', linestyle=':', lw=1.5, label='Dislocation Center')

ax2.set_xlim(80, 120)
ax2.set_ylim(bottom=-0.002, top=0.072)
ax2.set_xlabel('Atomic Position Along Slip Plane $x / h$', fontsize=13, fontweight='bold', labelpad=8)
ax2.set_ylabel('Slip Plane Energy $E(x)$', fontsize=13, fontweight='bold', labelpad=8)
ax2.set_title(r'Core Zoom: Slip Plane Energy Profile ($x \in [80, 120]h$)', fontsize=14, fontweight='bold', pad=12)
ax2.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.95, edgecolor='#cccccc', fontsize=11)
ax2.grid(True)
plt.tight_layout()

f2_path = os.path.join(out_dir, 'dislocation_200_slip_plane_zoom.png')
f2_art = os.path.join(artifact_dir, 'dislocation_200_slip_plane_zoom.png')
fig2.savefig(f2_path)
fig2.savefig(f2_art)
plt.close(fig2)
print(f'Saved: {f2_path}')

# ==============================================================================
# FIGURE 3: Semi-Log Cumulative Energy E(R) vs ln(R/h)
# ==============================================================================
radii = np.linspace(0.5, 30.0, 60)
cum_e1 = np.array([np.sum(e1[r1 <= R]) for R in radii])
cum_e2 = np.array([np.sum(e2[r2 <= R]) for R in radii])

fig3, ax3 = plt.subplots(figsize=(10, 5.5), dpi=300)
ax3.plot(radii, cum_e1, 'o-', color='#d9534f', lw=2.0, ms=4, alpha=0.9, label='Without Remeshing')
ax3.plot(radii, cum_e2, 's--', color='#0275d8', lw=1.8, ms=4, mfc='none', markeredgewidth=1.6, label='With Remeshing')
ax3.set_xscale('log')

# Continuum linear elastic fit in range R in [6, 28]
mask_fit = (radii >= 6.0) & (radii <= 28.0)
fit_p = np.polyfit(np.log(radii[mask_fit]), cum_e1[mask_fit], 1)
fit_r = np.geomspace(5.0, 30.0, 50)
ax3.plot(fit_r, fit_p[0]*np.log(fit_r) + fit_p[1], 'k-.', lw=2.0, label='Continuum Slope')

ax3.axvline(x=5.0, color='gray', ls=':', lw=1.5, label='Core Cutoff')

ax3.set_xlim(0.5, 32.0)
ax3.set_xlabel('Radius from Core $R / h$', fontsize=13, fontweight='bold', labelpad=8)
ax3.set_ylabel(r'Cumulative Strain Energy $\sum_{r_i \leq R} E_i$', fontsize=13, fontweight='bold', labelpad=8)
ax3.set_title(r'Semi-Log: Cumulative Strain Energy $E(R)$ vs. $\ln(R/h)$', fontsize=14, fontweight='bold', pad=12)
ax3.legend(loc='lower right', frameon=True, facecolor='white', framealpha=0.95, edgecolor='#cccccc', fontsize=11)
ax3.grid(True, which='both')
plt.tight_layout()

f3_path = os.path.join(out_dir, 'dislocation_200_log_cumulative_energy.png')
f3_art = os.path.join(artifact_dir, 'dislocation_200_log_cumulative_energy.png')
fig3.savefig(f3_path)
fig3.savefig(f3_art)
plt.close(fig3)
print(f'Saved: {f3_path}')

# ==============================================================================
# FIGURE 4: Log-Log Energy Density vs Distance (Radial Shell-Averaged)
# ==============================================================================
r_eval = np.geomspace(0.7, 30.0, 45)
smooth_e1 = []
smooth_e2 = []

for r in r_eval:
    # Logarithmic radial kernel averaging over full 2pi azimuth:
    # filters out anisotropic angular variation cos^2(2theta) and discrete shell gaps
    w1 = np.exp(-0.5 * ((np.log(r1) - np.log(r)) / 0.14)**2)
    w2 = np.exp(-0.5 * ((np.log(r2) - np.log(r)) / 0.14)**2)
    smooth_e1.append(np.sum(e1 * w1) / np.sum(w1))
    smooth_e2.append(np.sum(e2 * w2) / np.sum(w2))

smooth_e1 = np.array(smooth_e1)
smooth_e2 = np.array(smooth_e2)

fig4, ax4 = plt.subplots(figsize=(10, 5.5), dpi=300)

# Faint scatter of raw atoms to illustrate the anisotropic angular fan
mask_scatter = (r1 >= 0.7) & (r1 <= 30.0)
ax4.scatter(r1[mask_scatter], e1[mask_scatter], color='#d9534f', s=4, alpha=0.10, edgecolors='none')
ax4.scatter(r2[mask_scatter], e2[mask_scatter], color='#0275d8', s=4, alpha=0.10, edgecolors='none')

# Shell-averaged radial curves
ax4.plot(r_eval, smooth_e1, 'o-', color='#c9302c', lw=2.4, ms=5, label='Without Remeshing')
ax4.plot(r_eval, smooth_e2, 's--', color='#025aa5', lw=2.2, ms=5, mfc='none', markeredgewidth=1.8, label='With Remeshing')

# Theoretical 1/r^2 continuum scaling
r_guide = np.geomspace(3.0, 30.0, 50)
c_guide = 0.082
ax4.plot(r_guide, c_guide / (r_guide**2), 'k--', lw=2.0, label=r'Elastic Scaling $1/r^2$')
ax4.axvline(x=5.0, color='gray', ls=':', lw=1.5, label='Core Boundary')

ax4.set_xscale('log')
ax4.set_yscale('log')
ax4.set_xlim(0.65, 32.0)
ax4.set_ylim(5e-5, 1.2e-1)
ax4.set_xlabel('Distance from Core $r / h$', fontsize=13, fontweight='bold', labelpad=8)
ax4.set_ylabel(r'Radial Shell Energy Density $\bar{E}(r)$', fontsize=13, fontweight='bold', labelpad=8)
ax4.set_title(r'Log-Log: Radial Energy Density vs. Distance ($200 \times 200$)', fontsize=14, fontweight='bold', pad=12)

# Clear annotation explaining core plateau vs 1/r^2
ax4.annotate('Core Saturation Plateau\n(Finite Energy at $r \\to 0$)', xy=(0.85, 0.06), xytext=(1.2, 0.015),
             arrowprops=dict(facecolor='#333333', shrink=0.08, width=1, headwidth=5),
             fontsize=10, fontweight='bold', bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#cccccc'))

ax4.legend(loc='lower left', frameon=True, facecolor='white', framealpha=0.95, edgecolor='#cccccc', fontsize=11)
ax4.grid(True, which='both')
plt.tight_layout()

f4_path = os.path.join(out_dir, 'dislocation_200_log_energy_density.png')
f4_art = os.path.join(artifact_dir, 'dislocation_200_log_energy_density.png')
fig4.savefig(f4_path)
fig4.savefig(f4_art)
plt.close(fig4)
print(f'Saved: {f4_path}')

# ==============================================================================
# FIGURE 5: Semi-Log Energy Difference Delta E(R) vs ln(R/h)
# ==============================================================================
delta_e = cum_e1 - cum_e2

fig5, ax5 = plt.subplots(figsize=(10, 5.5), dpi=300)
ax5.plot(radii, delta_e, 'd-', color='#28a745', lw=2.2, ms=5, label='Energy Gain of Remeshing')
ax5.set_xscale('log')
ax5.axhline(y=0, color='#333333', linestyle='-', alpha=0.3)
ax5.axvline(x=5.0, color='gray', linestyle=':', lw=1.5, label='Core Cutoff')

idx_5 = np.argmin(np.abs(radii - 5.0))
idx_plateau = np.argmin(np.abs(radii - 10.0))
ax5.annotate(f'Core relaxation gain:\n$\\Delta E(10h) = {delta_e[idx_plateau]:.4f}$',
             xy=(radii[idx_plateau], delta_e[idx_plateau]), xytext=(12.0, 0.0060),
             arrowprops=dict(facecolor='#333333', shrink=0.08, width=1, headwidth=5),
             fontsize=10, fontweight='bold', bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#28a745', lw=1.5))

ax5.set_xlim(0.5, 32.0)
ax5.set_ylim(-0.0045, 0.0075)
ax5.set_xlabel('Radius from Core $R / h$', fontsize=13, fontweight='bold', labelpad=8)
ax5.set_ylabel(r'Energy Difference $\Delta E = E_{\rm fixed} - E_{\rm remeshed}$', fontsize=13, fontweight='bold', labelpad=8)
ax5.set_title(r'Lin-Log: Relaxation Gain $\Delta E(R)$ Concentrated in Dislocation Core', fontsize=14, fontweight='bold', pad=12)
ax5.legend(loc='lower right', frameon=True, facecolor='white', framealpha=0.95, edgecolor='#cccccc', fontsize=11)
ax5.grid(True, which='both')
plt.tight_layout()

f5_path = os.path.join(out_dir, 'dislocation_200_log_energy_difference.png')
f5_art = os.path.join(artifact_dir, 'dislocation_200_log_energy_difference.png')
fig5.savefig(f5_path)
fig5.savefig(f5_art)
plt.close(fig5)
print(f'Saved: {f5_path}')
