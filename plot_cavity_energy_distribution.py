#!/usr/bin/env python3
"""
Generate Pre-Yield vs Post-Yield Cavity / Patch Topological Jump Distributions:
- Decomposes remeshing passes into spatially disjoint reconnected element cavities / patches p.
- For each cavity p:
    \\Delta E_cavity^(p) = sum_{T in added} E_T^(after) - sum_{T in removed} E_T^(before)
    \\Delta \\sigma_{xy, cavity}^(p) = \\bar{\\sigma}_{xy, after}^(p) - \\bar{\\sigma}_{xy, before}^(p)
- Plots:
    (a) Cavity energy jump upon remeshing |\\Delta E_cavity^(p)|
    (b) Cavity stress jump upon remeshing |\\Delta \\sigma_{xy, cavity}^(p)|
- Stepped histograms for pre-yield and post-yield with (a)/(b) labels and Min, Mean, Max stats.
"""

import csv
from pathlib import Path
from collections import defaultdict
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def setup_style():
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["DejaVu Serif", "Times New Roman", "Computer Modern Roman"],
        "mathtext.fontset": "cm",
        "font.size": 11.0,
        "axes.titlesize": 12.0,
        "axes.labelsize": 11.5,
        "xtick.labelsize": 10.5,
        "ytick.labelsize": 10.5,
        "legend.fontsize": 9.2,
        "axes.linewidth": 1.1,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.size": 4.5,
        "ytick.major.size": 4.5,
        "xtick.minor.size": 2.5,
        "ytick.minor.size": 2.5,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.08,
    })

def find_yield_alpha(sim_dir):
    log_file = sim_dir / "energy_stress_log.csv"
    with open(log_file, "r") as f:
        rows = list(csv.DictReader(f))
    max_row = max(rows, key=lambda r: float(r["PostStress"]))
    return float(max_row["Alpha"]), float(max_row["PostStress"])

def load_mesh(vtk_path):
    r = vtk.vtkUnstructuredGridReader()
    r.SetFileName(str(vtk_path))
    r.Update()
    grid = r.GetOutput()
    pts = vtk_to_numpy(grid.GetPoints().GetData())[:, :2]
    ee = vtk_to_numpy(grid.GetCellData().GetArray("ElementEnergy"))
    stress_arr = grid.GetCellData().GetArray("CauchyStress")
    if stress_arr is not None:
        s_xy = vtk_to_numpy(stress_arr)[:, 1]
    else:
        s_xy = np.zeros_like(ee)
    
    cells_data = grid.GetCells()
    tri = vtk_to_numpy(cells_data.GetConnectivityArray()).reshape(-1, 3)
    return pts, tri, ee, s_xy

def tri_key(pts, t):
    p = pts[t]
    coords = [(round(p[j, 0], 3), round(p[j, 1], 3)) for j in range(3)]
    coords.sort()
    return tuple(coords)

def tri_area(pts, t):
    p0, p1, p2 = pts[t[0]], pts[t[1]], pts[t[2]]
    return 0.5 * abs((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]))

def extract_cavity_data(sim_dir, alpha_yield):
    log_file = sim_dir / "energy_stress_log.csv"
    step_to_alpha = {}
    with open(log_file, "r") as f:
        for r in csv.DictReader(f):
            step_to_alpha[int(r["Iteration"])] = float(r["Alpha"])

    surgery_dir = sim_dir / "vtk_surgery"
    before_files = sorted(surgery_dir.glob("*_1_before_remesh.vtk"))
    print(f"Analyzing {len(before_files)} surgery passes for cavity decomposition...")

    pre_dE, post_dE = [], []
    pre_dS, post_dS = [], []

    for pass_idx, bf in enumerate(before_files):
        af = Path(str(bf).replace("_1_before_remesh.vtk", "_2_after_remesh.vtk"))
        if not af.exists():
            continue

        try:
            step_idx = int(bf.name.split("_")[1])
        except:
            step_idx = 0
        alpha_val = step_to_alpha.get(step_idx, 0.14)

        pts1, tri1, ee1, s1 = load_mesh(bf)
        pts2, tri2, ee2, s2 = load_mesh(af)

        # Precompute areas
        area1 = [tri_area(pts1, t) for t in tri1]
        area2 = [tri_area(pts2, t) for t in tri2]

        set1 = {tri_key(pts1, t): (ee1[i], s1[i], area1[i]) for i, t in enumerate(tri1)}
        set2 = {tri_key(pts2, t): (ee2[i], s2[i], area2[i]) for i, t in enumerate(tri2)}

        common = set(set1.keys()) & set(set2.keys())
        rem_keys = set(set1.keys()) - common
        add_keys = set(set2.keys()) - common

        if not rem_keys and not add_keys:
            continue

        # Disjoint-set union-find over shared vertices
        parent = {}
        def find(x):
            parent.setdefault(x, x)
            if parent[x] != x:
                parent[x] = find(parent[x])
            return parent[x]

        def union(x, y):
            rx, ry = find(x), find(y)
            if rx != ry:
                parent[rx] = ry

        for k in rem_keys:
            find(k)
            for pt in k:
                union(k, ("node", pt))

        for k in add_keys:
            find(k)
            for pt in k:
                union(k, ("node", pt))

        # Group elements by cavity
        cavities = defaultdict(lambda: {"E_rem": 0.0, "E_add": 0.0, "S_rem": 0.0, "S_add": 0.0, "A_rem": 0.0, "A_add": 0.0})
        for k in rem_keys:
            r = find(k)
            e, s, a = set1[k]
            cavities[r]["E_rem"] += e
            cavities[r]["S_rem"] += s * a
            cavities[r]["A_rem"] += a

        for k in add_keys:
            r = find(k)
            e, s, a = set2[k]
            cavities[r]["E_add"] += e
            cavities[r]["S_add"] += s * a
            cavities[r]["A_add"] += a

        for c in cavities.values():
            dE = abs(c["E_add"] - c["E_rem"])
            
            # Area-weighted stress average jump
            a_tot = max(c["A_add"], c["A_rem"], 1e-12)
            s_after = c["S_add"] / a_tot if c["A_add"] > 0 else 0.0
            s_before = c["S_rem"] / a_tot if c["A_rem"] > 0 else 0.0
            dS = abs(s_after - s_before)

            if alpha_val < alpha_yield:
                if dE > 1e-9:
                    pre_dE.append(dE)
                if dS > 1e-6:
                    pre_dS.append(dS)
            else:
                if dE > 1e-9:
                    post_dE.append(dE)
                if dS > 1e-6:
                    post_dS.append(dS)

        if (pass_idx + 1) % 50 == 0 or (pass_idx + 1) == len(before_files):
            print(f"  Processed {pass_idx + 1}/{len(before_files)} passes...")

    return {
        "pre_dE": np.array(pre_dE),
        "post_dE": np.array(post_dE),
        "pre_dS": np.array(pre_dS),
        "post_dS": np.array(post_dS),
    }

def main():
    setup_style()
    curr_dir = Path(__file__).resolve().parent
    if (curr_dir / "vtk_surgery").exists():
        sim_dir = curr_dir
        repo_root = curr_dir.parent
    elif (curr_dir / "study_100x100_positive" / "vtk_surgery").exists():
        sim_dir = curr_dir / "study_100x100_positive"
        repo_root = curr_dir
    else:
        raise FileNotFoundError("Could not find simulation directory containing 'vtk_surgery'")

    alpha_yield, peak_stress = find_yield_alpha(sim_dir)
    print("=================================================================")
    print(f"DYNAMIC YIELD DETECTION:")
    print(f"  alpha_yield = {alpha_yield:.5f}")
    print(f"  sigma_yield = {peak_stress:.5f}")
    print("=================================================================")

    data = extract_cavity_data(sim_dir, alpha_yield)

    pre_dE = data["pre_dE"]
    post_dE = data["post_dE"]
    pre_dS = data["pre_dS"]
    post_dS = data["post_dS"]

    stats = {
        "pre_E": {"n": len(pre_dE), "mean": np.mean(pre_dE), "min": np.min(pre_dE), "max": np.max(pre_dE)},
        "post_E": {"n": len(post_dE), "mean": np.mean(post_dE), "min": np.min(post_dE), "max": np.max(post_dE)},
        "pre_S": {"n": len(pre_dS), "mean": np.mean(pre_dS), "min": np.min(pre_dS), "max": np.max(pre_dS)},
        "post_S": {"n": len(post_dS), "mean": np.mean(post_dS), "min": np.min(post_dS), "max": np.max(post_dS)},
    }

    print("\n--- Summary Statistics (Cavity Jumps) ---")
    for k, v in stats.items():
        print(f"{k}: n={v['n']:,}, Min={v['min']:.6e}, Mean={v['mean']:.6e}, Max={v['max']:.6e}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.8), dpi=300)
    fig.subplots_adjust(wspace=0.28, left=0.08, right=0.97, top=0.92, bottom=0.14)

    color_pre = "#1f77b4"
    color_post = "#d95f02"

    # =========================================================================
    # LEFT PANEL: Cavity energy jump |\Delta E_cavity^(p)|
    # =========================================================================
    bins_E = np.logspace(-8, 1, 30)
    h_pre_E, _ = np.histogram(pre_dE, bins=bins_E)
    h_post_E, _ = np.histogram(post_dE, bins=bins_E)

    prob_pre_E = h_pre_E / len(pre_dE)
    prob_post_E = h_post_E / len(post_dE)

    ax1.stairs(prob_pre_E, bins_E, color=color_pre, lw=2.2, label=r"pre-yield")
    ax1.stairs(prob_post_E, bins_E, color=color_post, lw=2.2, label=r"post-yield")

    E_max_ref = 0.0481
    ax1.axvline(E_max_ref, color="#222222", ls="--", lw=1.3,
                label=r"single-element $E_{\mathrm{max}}$ ($\gamma = 0.5$)")

    ax1.set_xscale("log")
    ax1.set_xlim(1e-9, 1e1)
    ax1.set_ylim(0, max(np.max(prob_pre_E), np.max(prob_post_E)) * 1.25)
    ax1.set_xlabel(r"Cavity energy jump upon remeshing $|\Delta E_{\mathrm{cavity}}^{(p)}|$")
    ax1.set_ylabel(r"Probability per logarithmic bin")
    ax1.legend(frameon=True, facecolor="white", edgecolor="#d0d0d0",
               framealpha=0.92, loc="upper left", fontsize=9.2)

    ax1.text(0.04, 0.21, r"$\mathbf{(a)}$", transform=ax1.transAxes,
             fontsize=13.0, fontweight="bold", va="bottom", ha="left")

    stats_text_E = (
        f"Pre : Min = {stats['pre_E']['min']:.1e}   Mean = {stats['pre_E']['mean']:.2e}   Max = {stats['pre_E']['max']:.3f}\n"
        f"Post: Min = {stats['post_E']['min']:.1e}   Mean = {stats['post_E']['mean']:.2e}   Max = {stats['post_E']['max']:.3f}"
    )
    ax1.text(0.04, 0.05, stats_text_E, transform=ax1.transAxes,
             fontsize=8.5, family="monospace", va="bottom", ha="left",
             bbox=dict(boxstyle="round,pad=0.35", facecolor="#fafafa", edgecolor="#cccccc", lw=0.8))

    # =========================================================================
    # RIGHT PANEL: Cavity stress jump |\Delta \sigma_{xy, cavity}^(p)|
    # =========================================================================
    bins_S = np.logspace(-6, 0, 25)
    h_pre_S, _ = np.histogram(pre_dS, bins=bins_S)
    h_post_S, _ = np.histogram(post_dS, bins=bins_S)

    prob_pre_S = h_pre_S / len(pre_dS)
    prob_post_S = h_post_S / len(post_dS)

    ax2.stairs(prob_pre_S, bins_S, color=color_pre, lw=2.2, label=r"pre-yield")
    ax2.stairs(prob_post_S, bins_S, color=color_post, lw=2.2, label=r"post-yield")

    S_ellipticity_ref = 0.334
    ax2.axvline(S_ellipticity_ref, color="#222222", ls="--", lw=1.3,
                label=r"$\sigma_{xy}$ at loss of ellipticity ($\gamma = 0.1322$)")

    ax2.set_xscale("log")
    ax2.set_xlim(5e-7, 1.2e0)
    ax2.set_ylim(0, max(np.max(prob_pre_S), np.max(prob_post_S)) * 1.25)
    ax2.set_xlabel(r"Cavity stress jump upon remeshing $|\Delta \sigma_{xy, \mathrm{cavity}}^{(p)}|$")
    ax2.set_ylabel(r"Probability per logarithmic bin")
    ax2.legend(frameon=True, facecolor="white", edgecolor="#d0d0d0",
               framealpha=0.92, loc="upper left", fontsize=9.2)

    ax2.text(0.04, 0.21, r"$\mathbf{(b)}$", transform=ax2.transAxes,
             fontsize=13.0, fontweight="bold", va="bottom", ha="left")

    stats_text_S = (
        f"Pre : Min = {stats['pre_S']['min']:.1e}   Mean = {stats['pre_S']['mean']:.2e}   Max = {stats['pre_S']['max']:.3f}\n"
        f"Post: Min = {stats['post_S']['min']:.1e}   Mean = {stats['post_S']['mean']:.2e}   Max = {stats['post_S']['max']:.3f}"
    )
    ax2.text(0.04, 0.05, stats_text_S, transform=ax2.transAxes,
             fontsize=8.5, family="monospace", va="bottom", ha="left",
             bbox=dict(boxstyle="round,pad=0.35", facecolor="#fafafa", edgecolor="#cccccc", lw=0.8))

    out_dir = repo_root / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / "cavity_jump_distributions_comparison.png"
    out_pdf = out_dir / "cavity_jump_distributions_comparison.pdf"

    plt.savefig(out_pdf)
    plt.savefig(out_png)
    plt.close()
    print(f"\nFigure saved to:\n  {out_png}\n  {out_pdf}")

if __name__ == "__main__":
    main()
