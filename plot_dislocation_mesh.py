#!/usr/bin/env python3
"""Plot the finite element triangulation mesh and atomic positions around a dislocation core, colored by nodal energy."""

import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri


def read_legacy_vtk(path):
    """Parse a legacy ASCII VTK unstructured grid, returning points, cells, and named nodal scalar fields."""
    with open(path, "rb") as f:
        v_line = f.readline().decode("latin1", errors="ignore").strip()
        t_line = f.readline().decode("latin1", errors="ignore").strip()
        fmt_line = f.readline().decode("latin1", errors="ignore").strip().upper()
        ds_line = f.readline().decode("latin1", errors="ignore").strip()

        is_binary = "BINARY" in fmt_line

        if not is_binary:
            f.seek(0)
            lines = f.read().decode("latin1", errors="ignore").splitlines()
            points = None
            cells = []
            scalars = {}
            mode = None
            n_point = None
            n_cell = None
            i = 0

            while i < len(lines):
                parts = lines[i].split()
                if not parts:
                    i += 1
                    continue
                key = parts[0]

                if key == "POINTS":
                    n = int(parts[1])
                    vals = []
                    i += 1
                    while len(vals) < 3 * n:
                        vals.extend(float(x) for x in lines[i].split())
                        i += 1
                    points = np.array(vals, dtype=float).reshape(n, 3)[:, :2]
                    continue

                if key == "CELLS":
                    n = int(parts[1])
                    i += 1
                    for _ in range(n):
                        row = [int(x) for x in lines[i].split()]
                        if row and row[0] == 3:
                            cells.append(row[1:4])
                        i += 1
                    continue

                if key == "CELL_DATA":
                    mode = "cell"
                    n_cell = int(parts[1])
                    i += 1
                    continue

                if key == "POINT_DATA":
                    mode = "point"
                    n_point = int(parts[1])
                    i += 1
                    continue

                if key == "SCALARS":
                    name = parts[1]
                    n_expected = n_point if mode == "point" else n_cell
                    i += 1
                    if i < len(lines) and lines[i].split()[:1] == ["LOOKUP_TABLE"]:
                        i += 1
                    vals = []
                    while len(vals) < n_expected:
                        vals.extend(float(x) for x in lines[i].split())
                        i += 1
                    if mode == "point":
                        scalars[name] = np.array(vals[:n_expected], dtype=float)
                    continue

                i += 1

            return {
                "points": np.asarray(points),
                "triangles": np.asarray(cells, dtype=int),
                "scalars": scalars,
            }

        # BINARY parsing
        points = None
        cells = []
        scalars = {}
        mode = None
        n_point = None
        n_cell = None

        while True:
            line_bytes = f.readline()
            if not line_bytes:
                break
            line = line_bytes.decode("latin1", errors="ignore").strip()
            if not line:
                continue
            parts = line.split()
            key = parts[0]

            if key == "POINTS":
                n = int(parts[1])
                raw = np.fromfile(f, dtype='>f4', count=3 * n)
                points = raw.reshape(n, 3)[:, :2].astype(np.float64)
                continue

            if key == "CELLS":
                n = int(parts[1])
                total_ints = int(parts[2])
                raw = np.fromfile(f, dtype='>i4', count=total_ints)
                if total_ints == 4 * n:
                    cells = raw.reshape(n, 4)[:, 1:4]
                else:
                    idx = 0
                    c_list = []
                    for _ in range(n):
                        c_len = raw[idx]
                        if c_len == 3:
                            c_list.append(raw[idx+1:idx+4])
                        idx += c_len + 1
                    cells = np.array(c_list, dtype=np.int32)
                continue

            if key == "CELL_DATA":
                mode = "cell"
                n_cell = int(parts[1])
                continue

            if key == "POINT_DATA":
                mode = "point"
                n_point = int(parts[1])
                continue

            if key == "SCALARS":
                name = parts[1]
                _ = f.readline()  # LOOKUP_TABLE
                n_expected = n_point if mode == "point" else n_cell
                arr = np.fromfile(f, dtype='>f4', count=n_expected).astype(np.float64)
                if mode == "point":
                    scalars[name] = arr
                continue

        return {
            "points": np.asarray(points),
            "triangles": np.asarray(cells, dtype=int),
            "scalars": scalars,
        }


def find_core(points, scalars):
    """Locate the dislocation core node using NodalEnergy."""
    if "NodalEnergy" in scalars:
        return int(np.argmax(scalars["NodalEnergy"]))
    centre = points.mean(axis=0)
    return int(np.argmin(((points - centre) ** 2).sum(axis=1)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vtk", type=Path, help="Input configuration VTK file.")
    parser.add_argument("--radius", type=float, default=15.0,
                        help="Radius of the window around the core to plot (default: 15.0).")
    parser.add_argument("--out-dir", type=Path, default=Path("plots"),
                        help="Output directory for plots (default: plots/).")
    parser.add_argument("--dpi", type=int, default=300,
                        help="Output resolution in dots per inch (default: 300).")
    parser.add_argument("--cmap", default="plasma",
                        help="Colormap for the energy coloring (default: plasma).")
    parser.add_argument("--shading", default="flat", choices=("flat", "gouraud"),
                        help="Shading style for triangulation: flat or gouraud (default: flat).")
    parser.add_argument("--format", default="png", choices=("png", "pdf", "both"),
                        help="Output plot format: png, pdf, or both (default: png).")
    parser.add_argument("--tag", default="",
                        help="Extra tag inserted into output filename.")
    args = parser.parse_args()

    # Load VTK data
    data = read_legacy_vtk(args.vtk)
    points = data["points"]
    triangles = data["triangles"]
    scalars = data["scalars"]

    if "NodalEnergy" not in scalars:
        raise SystemExit(
            f"NodalEnergy field not found in {args.vtk}. "
            f"Available scalar fields: {list(scalars.keys())}"
        )

    energy = scalars["NodalEnergy"]

    # Find core coordinates
    core_idx = find_core(points, scalars)
    core_x, core_y = points[core_idx]
    print(f"Auto-detected core at node {core_idx}, (x,y)=({core_x:.3f}, {core_y:.3f})")

    # Create Matplotlib Triangulation
    triang = mtri.Triangulation(points[:, 0], points[:, 1], triangles)

    # Plot setup
    fig, ax = plt.subplots(figsize=(6.0, 5.5), dpi=args.dpi)

    # Plot colored triangles
    tpc = ax.tripcolor(
        triang,
        energy,
        shading=args.shading,
        cmap=args.cmap
    )

    # Draw triangulation mesh lines
    ax.triplot(triang, color="black", linewidth=0.25, alpha=0.35)

    # Draw nodes/atoms
    ax.scatter(points[:, 0], points[:, 1], s=4.0, c="black", alpha=0.8, linewidths=0, zorder=3)

    # Mark the core center node with a distinct white star
    ax.scatter(core_x, core_y, s=140.0, color="white", edgecolor="black", marker="*",
               linewidths=0.8, zorder=5, label="Core Peak")

    # Crop viewport around the core
    ax.set_xlim(core_x - args.radius, core_x + args.radius)
    ax.set_ylim(core_y - args.radius, core_y + args.radius)
    ax.set_aspect("equal", adjustable="box")

    # Axis Labels
    ax.set_xlabel(r"$x$", fontsize=10)
    ax.set_ylabel(r"$y$", fontsize=10)
    ax.tick_params(labelsize=8)
    ax.grid(alpha=0.15, linestyle="--")

    # Add title and colorbar
    ax.set_title(f"Mesh Around Dislocation Core ({args.vtk.name})\nNodal Energy Coloring ({args.shading} shading)",
                 fontsize=9, pad=10)

    cbar = fig.colorbar(tpc, ax=ax, shrink=0.8)
    cbar.set_label("Nodal Energy", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    fig.tight_layout()

    # Save outputs
    args.out_dir.mkdir(parents=True, exist_ok=True)
    suffix = f"_{args.tag}" if args.tag else ""
    formats = ["png", "pdf"] if args.format == "both" else [args.format]
    for fmt in formats:
        out = args.out_dir / f"mesh_core_{args.vtk.stem}{suffix}.{fmt}"
        fig.savefig(out, dpi=args.dpi)
        print(f"Wrote {out}")

    plt.close(fig)


if __name__ == "__main__":
    main()
