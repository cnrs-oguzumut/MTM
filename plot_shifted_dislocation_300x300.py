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

base_dir = 'shifted_dislocation_study_300x300'
artifact_dir = '/Users/usalman/.gemini/antigravity-ide/brain/e0366a0b-42af-4a53-9525-edb8c926ac89'
shifts = [0, 1, 2, 3, 4, 5]
colors = ['#2b5c8f', '#2a9d8f', '#e76f51', '#f4a261', '#9b5de5', '#d62828']
markers = ['o', 's', '^', 'D', 'v', 'p']

data = {}
initial_core_x = 149.5

print("Parsing simulation data...")
for s in shifts:
    vtk_path = os.path.join(base_dir, f'shift_{s}/vtk_output/configuration_00001.vtk')
    if not os.path.exists(vtk_path):
        print(f"Warning: {vtk_path} does not exist yet.")
        continue
    pts, e = parse_vtk_nodal_data(vtk_path)
    x_slip, e_slip = extract_slip_plane_rows(pts, e)
    x_cm = compute_slip_centroid(x_slip, e_slip)
    
    # Distance from true centroid
    r_from_cm = np.linalg.norm(pts - np.array([x_cm, 149.5]), axis=1)
    
    data[s] = {
        'pts': pts,
        'e': e,
        'x_slip': x_slip,
        'e_slip': e_slip,
        'x_cm': x_cm,
        'r_from_cm': r_from_cm
    }
    print(f"Shift s={s}: centroid x_cm = {x_cm:.4f}, max energy = {np.max(e_slip):.5f}")

if not data:
    print("No data available yet to plot.")
    exit(0)

# ==============================================================================
# FIGURE 1: Full Slip Plane Energy Profile (Across entire domain [-150, 150]h)
# ==============================================================================
fig1, ax1 = plt.subplots(figsize=(10, 5.5), dpi=300)
for s, col, mark in zip(shifts, colors, markers):
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

f1_path = os.path.join(base_dir, 'dislocation_300_shifted_slip_plane_full.png')
f1_art = os.path.join(artifact_dir, 'dislocation_300_shifted_slip_plane_full.png')
fig1.savefig(f1_path)
fig1.savefig(f1_art)
plt.close(fig1)
print(f"Saved: {f1_path}")

# ==============================================================================
# FIGURE 2: Zoomed Core Region Aligned at Centroid (x - x_cm) in [-25, 25]h
# ==============================================================================
fig2, ax2 = plt.subplots(figsize=(10, 5.5), dpi=300)
for s, col, mark in zip(shifts, colors, markers):
    if s not in data: continue
    x_centered = data[s]['x_slip'] - data[s]['x_cm']
    mask_z = (x_centered >= -25.0) & (x_centered <= 25.0)
    ax2.plot(x_centered[mask_z], data[s]['e_slip'][mask_z], 'o-', color=col, lw=2.0, ms=4.5, label=f'Shift $s = {s}$')

ax2.axvline(x=0.0, color='#333333', linestyle=':', lw=1.5, label='Dislocation Center')
ax2.set_xlim(-25, 25)
ax2.set_ylim(bottom=-0.002, top=0.075)
ax2.set_xlabel(r'Relative Coordinate Along Slip Plane $(x - x_c) / h$', fontsize=13, fontweight='bold', labelpad=8)
ax2.set_ylabel('Slip Plane Energy $E(x)$', fontsize=13, fontweight='bold', labelpad=8)
ax2.set_title(r'Core Zoom: Aligned Slip Plane Profile ($(x - x_c) \in [-25, 25]h$)', fontsize=14, fontweight='bold', pad=12)
ax2.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.95, edgecolor='#cccccc', fontsize=10, ncol=2)
ax2.grid(True)
plt.tight_layout()

f2_path = os.path.join(base_dir, 'dislocation_300_shifted_slip_plane_zoom.png')
f2_art = os.path.join(artifact_dir, 'dislocation_300_shifted_slip_plane_zoom.png')
fig2.savefig(f2_path)
fig2.savefig(f2_art)
plt.close(fig2)
print(f"Saved: {f2_path}")

# ==============================================================================
# FIGURE 3: Semi-Log Cumulative Energy E(R) vs ln(R/h)
# ==============================================================================
radii = np.linspace(0.5, 70.0, 100)
fig3, ax3 = plt.subplots(figsize=(10, 5.5), dpi=300)

for s, col, mark in zip(shifts, colors, markers):
    if s not in data: continue
    r = data[s]['r_from_cm']
    e = data[s]['e']
    cum_e = np.array([np.sum(e[r <= R]) for R in radii])
    ax3.plot(radii, cum_e, 'o-', color=col, lw=1.8, ms=3.0, alpha=0.9, label=f'Shift $s = {s}$')

# Continuum slope reference
# Continuum slope prefactor: K b^2 / (4 pi) = 0.4332 (Foreman 1955)
c_slope = 0.433234
if 0 in data:
    r0_data = data[0]['r_from_cm']
    e0_data = data[0]['e']
    cum_e0 = np.array([np.sum(e0_data[r0_data <= R]) for R in radii])
    mask_fit = (radii >= 10.0) & (radii <= 60.0)
    fit_p = np.polyfit(np.log(radii[mask_fit]), cum_e0[mask_fit], 1)
    ax3.plot(fit_r, fit_p[0] * np.log(fit_r) + fit_p[1], 'k-.', lw=2.0, label='Continuum Slope')

ax3.axvline(x=5.0, color='gray', ls=':', lw=1.5, label='Core Cutoff')
ax3.set_xscale('log')
ax3.set_xlim(0.5, 75.0)
ax3.set_xlabel(r'Radius from Core $R / h$', fontsize=13, fontweight='bold', labelpad=8)
ax3.set_ylabel(r'Cumulative Strain Energy $\sum_{r_i \leq R} E_i$', fontsize=13, fontweight='bold', labelpad=8)
ax3.set_title(r'Semi-Log: Cumulative Strain Energy $E(R)$ vs. $\ln(R/h)$ ($300 \times 300$, $R_{\rm free} = 70h$)', fontsize=14, fontweight='bold', pad=12)
ax3.legend(loc='lower right', frameon=True, facecolor='white', framealpha=0.95, edgecolor='#cccccc', fontsize=10, ncol=2)
ax3.grid(True, which='both')
plt.tight_layout()

f3_path = os.path.join(base_dir, 'dislocation_300_shifted_log_cumulative_energy.png')
f3_art = os.path.join(artifact_dir, 'dislocation_300_shifted_log_cumulative_energy.png')
fig3.savefig(f3_path)
fig3.savefig(f3_art)
plt.close(fig3)
print(f"Saved: {f3_path}")
print("All 3 figures successfully generated.")
