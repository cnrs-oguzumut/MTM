#!/usr/bin/env python3
"""
Plot Internal Energy vs Staircase Translations (tx, ty)
Reconstructing the comparison between:
- Fixed connectivity, orientation I
- Fixed connectivity, orientation II
- Reconnection
"""

import argparse
import os
import shutil
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import csv
import numpy as np

def main():
    parser = argparse.ArgumentParser(description="Plot staircase shifting internal energy comparison.")
    parser.add_argument("--dir", type=str, default="final_tests", help="Base directory containing simulation runs")
    parser.add_argument("--output", type=str, default="shifting_internal_energy_comparison", help="Output filename base (without extension)")
    args = parser.parse_args()

    base_dir = Path(args.dir)
    if not base_dir.exists():
        alt = Path("final_simu") / "final_tests"
        if alt.exists():
            base_dir = alt
        else:
            raise FileNotFoundError(f"Could not find directory '{base_dir}' or '{alt}'")

    cases = [
        {
            "name": "Fixed connectivity, orientation I",
            "folder": base_dir / "simulation_01_left_bottom_no_remesh_amp2",
            "color": "black",
            "linestyle": "--",
            "linewidth": 1.8,
            "dashes": (5, 3),
            "zorder": 3,
        },
        {
            "name": "Fixed connectivity, orientation II",
            "folder": base_dir / "simulation_05_left_bottom_perturbed_no_remesh_amp2",
            "color": "red",
            "linestyle": "-.",
            "linewidth": 1.8,
            "dashes": (6, 2, 1, 2),
            "zorder": 4,
        },
        {
            "name": "Reconnection",
            "folder": base_dir / "simulation_02_left_bottom_remesh_amp2",
            "color": "blue",
            "linestyle": "-",
            "linewidth": 1.8,
            "dashes": None,
            "zorder": 5,
        },
    ]

    # Matplotlib styling for publication / gnuplot-like appearance
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 14,
        "axes.labelsize": 18,
        "axes.titlesize": 16,
        "xtick.labelsize": 13,
        "ytick.labelsize": 14,
        "legend.fontsize": 14,
        "mathtext.fontset": "cm",
        "lines.antialiased": True,
    })

    fig, ax = plt.subplots(figsize=(8.0, 5.8), dpi=300)

    for case in cases:
        csv_path = case["folder"] / "energy_stress_log.csv"
        if not csv_path.exists():
            print(f"Warning: {csv_path} does not exist yet. Skipping {case['name']}.")
            continue

        x_vals, y_vals = [], []
        with open(csv_path, mode="r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                x_vals.append(float(row["Alpha"]))
                y_vals.append(float(row["PostEnergy"]))
        x = np.array(x_vals)
        y = np.array(y_vals)

        plot_kwargs = {
            "label": case["name"],
            "color": case["color"],
            "linewidth": case["linewidth"],
            "zorder": case["zorder"],
        }
        if case["dashes"] is not None:
            plot_kwargs["dashes"] = case["dashes"]
        else:
            plot_kwargs["linestyle"] = case["linestyle"]

        ax.plot(x, y, **plot_kwargs)

    # Inward ticks on all 4 sides
    ax.tick_params(direction="in", top=True, right=True, which="both")
    ax.tick_params(which="major", length=6, width=0.9)
    ax.tick_params(which="minor", length=3, width=0.6)

    # X-axis setup
    ax.set_xlim(0.0, 4.0)
    ax.set_xticks([0.0, 1.0, 2.0, 3.0, 4.0])
    ax.set_xticklabels([
        r"$\mathrm{t}_x=0$" + "\n" + r"$\mathrm{t}_y=0$",
        r"$\mathrm{t}_x=1$" + "\n" + r"$\mathrm{t}_y=0$",
        r"$\mathrm{t}_x=1$" + "\n" + r"$\mathrm{t}_y=1$",
        r"$\mathrm{t}_x=2$" + "\n" + r"$\mathrm{t}_y=1$",
        r"$\mathrm{t}_x=2$" + "\n" + r"$\mathrm{t}_y=2$",
    ])
    # No x minor ticks in original gnuplot
    ax.xaxis.set_minor_locator(plt.NullLocator())

    # Y-axis setup
    ax.set_ylabel("Internal energy", labelpad=8)
    ax.set_ylim(0.0, 3.5)
    ax.set_yticks([0, 1, 2, 3])
    ax.yaxis.set_minor_locator(MultipleLocator(0.2))

    # Legend in upper left with crisp, compact border that doesn't overlap curves
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(0.015, 0.985),
        frameon=True,
        edgecolor="black",
        fancybox=False,
        framealpha=1.0,
        facecolor="white",
        borderpad=0.35,
        labelspacing=0.25,
        handlelength=1.8,
        handletextpad=0.5,
        fontsize=11.5,
    )

    fig.tight_layout()

    out_png = f"{args.output}.png"
    out_pdf = f"{args.output}.pdf"
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"Saved figure to {out_png} and {out_pdf}")

    # Also copy to eliaspaper if directory exists
    elias_dir = Path("eliaspaper")
    if elias_dir.exists():
        shutil.copy2(out_pdf, elias_dir / out_pdf)
        shutil.copy2(out_png, elias_dir / out_png)
        print(f"Copied figure to {elias_dir / out_png} and {elias_dir / out_pdf}")

if __name__ == "__main__":
    main()
