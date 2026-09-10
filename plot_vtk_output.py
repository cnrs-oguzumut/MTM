#!/usr/bin/env python3
"""
plot_vtk_output.py - Flexible ParaView-style renderer for 2D VTK simulation outputs.

Features:
  - Fields:
      * 'nodal' / 'energy_nodal' (NodalEnergy)
      * 'element' / 'energy_element' (ElementEnergy)
      * 'stress_nodal', 'stress_element', 'coordination'
  - ParaView Representations:
      * 'colors' / 'surface': Smooth or flat colored surface without wireframe.
      * 'elements' / 'surface_with_edges': Colored surface with triangular element edges.
      * 'wireframe': Triangulation edges only.
      * 'nodes' / 'points': Point/atom scatter plot.
  - Flexible Overlays:
      * Combine any representation with --show-elements, --show-nodes, or --node-color field.
  - Batch Processing:
      * Multi-core rendering via multiprocessing (-j / --jobs).
      * Step filtering (--stride, --start, --end).
      * Optional movie generation (--movie / --fps).
  - High Customizability:
      * Colormaps (--cmap: turbo, viridis, plasma, jet, coolwarm, etc.).
      * Themes: 'white' (publication-ready) or 'dark' (ParaView default).
      * Zoom / Crop: --crop XMIN XMAX YMIN YMAX.

Usage Examples:
  # 1. Plot single file with nodal energy (surface only, like ParaView 'Surface'):
  python3 plot_vtk_output.py run_zanzotto_loading_150x150/vtk_output/configuration_00000.vtk --field nodal --mode colors

  # 2. Plot with element energy and mesh lines (like ParaView 'Surface With Edges'):
  python3 plot_vtk_output.py run_zanzotto_loading_150x150/vtk_output/configuration_00000.vtk --field element --mode elements

  # 3. Plot nodes/atoms:
  python3 plot_vtk_output.py run_zanzotto_loading_150x150/vtk_output/configuration_00000.vtk --field nodal --mode nodes

  # 4. Zoom in on a 40x40 window around a defect/avalanche:
  python3 plot_vtk_output.py run_zanzotto_loading_150x150/vtk_output/configuration_00000.vtk --crop 50 90 50 90 --show-elements

  # 5. Batch process an entire directory using 8 CPU cores:
  python3 plot_vtk_output.py run_zanzotto_loading_150x150/vtk_output/ --field nodal --mode colors -j 8

  # 6. Batch process and compile an MP4 movie:
  python3 plot_vtk_output.py run_zanzotto_loading_150x150/vtk_output/ --field nodal --mode colors --movie --fps 15
"""

import argparse
import sys
import shutil
import subprocess
from pathlib import Path
from multiprocessing import Pool, cpu_count
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri


def read_legacy_vtk(filepath):
    """
    Fast reader for legacy ASCII VTK Unstructured Grid files.
    Extracts points, triangular cells, cell scalars, point scalars, and field data.
    """
    path = Path(filepath)
    if not path.is_file():
        raise FileNotFoundError(f"File not found: {path}")

    lines = path.read_text().splitlines()
    points = None
    cells = []
    cell_scalars = {}
    point_scalars = {}
    field_data = {}
    mode = None
    current_count = 0
    i = 0
    n_lines = len(lines)

    while i < n_lines:
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        parts = line.split()
        tag = parts[0]

        if tag == "POINTS":
            n_pts = int(parts[1])
            i += 1
            pts_vals = []
            while len(pts_vals) < 3 * n_pts and i < n_lines:
                pts_vals.extend(float(x) for x in lines[i].split())
                i += 1
            points = np.array(pts_vals, dtype=np.float64).reshape(n_pts, 3)[:, :2]
            continue

        elif tag == "CELLS":
            n_cells = int(parts[1])
            i += 1
            for _ in range(n_cells):
                p = [int(x) for x in lines[i].split()]
                if len(p) >= 4 and p[0] == 3:
                    cells.append(p[1:4])
                i += 1
            continue

        elif tag == "CELL_DATA":
            mode = "cell"
            current_count = int(parts[1])
            i += 1
            continue

        elif tag == "POINT_DATA":
            mode = "point"
            current_count = int(parts[1])
            i += 1
            continue

        elif tag == "SCALARS":
            name = parts[1]
            i += 1
            if i < n_lines and lines[i].startswith("LOOKUP_TABLE"):
                i += 1
            vals = []
            while len(vals) < current_count and i < n_lines:
                vals.extend(float(x) for x in lines[i].split())
                i += 1
            arr = np.array(vals[:current_count], dtype=np.float64)
            if mode == "cell":
                cell_scalars[name] = arr
            elif mode == "point":
                point_scalars[name] = arr
            continue

        elif tag == "FIELD" and len(parts) >= 3:
            num_fields = int(parts[2])
            i += 1
            for _ in range(num_fields):
                f_meta = lines[i].split()
                f_name = f_meta[0]
                num_comp = int(f_meta[1])
                num_tuples = int(f_meta[2])
                i += 1
                f_vals = []
                while len(f_vals) < num_comp * num_tuples and i < n_lines:
                    f_vals.extend(float(x) for x in lines[i].split())
                    i += 1
                field_data[f_name] = f_vals[0] if len(f_vals) == 1 else f_vals
            continue

        i += 1

    return {
        "path": path,
        "points": points,
        "cells": np.array(cells, dtype=np.int32),
        "cell_scalars": cell_scalars,
        "point_scalars": point_scalars,
        "field_data": field_data,
    }


def resolve_field(data, requested_field):
    """
    Resolve requested field name to ('nodal'|'element', array, label).
    """
    req = requested_field.lower().strip()
    cell_scalars = data["cell_scalars"]
    point_scalars = data["point_scalars"]

    # Nodal Energy aliases
    if req in ("nodal", "energy_nodal", "nodalenergy", "nodal_energy"):
        if "NodalEnergy" in point_scalars:
            return "nodal", point_scalars["NodalEnergy"], "Nodal Energy"
        raise KeyError(f"'NodalEnergy' not found. Available point scalars: {list(point_scalars.keys())}")

    # Element Energy aliases
    if req in ("element", "cell", "energy_element", "elementenergy", "element_energy", "cell_energy"):
        if "ElementEnergy" in cell_scalars:
            return "element", cell_scalars["ElementEnergy"], "Element Energy"
        raise KeyError(f"'ElementEnergy' not found. Available cell scalars: {list(cell_scalars.keys())}")

    # Nodal Stress
    if req in ("stress_nodal", "nodalprojectedstress", "nodal_stress"):
        if "NodalProjectedStress" in point_scalars:
            return "nodal", point_scalars["NodalProjectedStress"], "Nodal Projected Stress"
        raise KeyError(f"'NodalProjectedStress' not found. Available: {list(point_scalars.keys())}")

    # Element Stress
    if req in ("stress_element", "elementprojectedstress", "element_stress"):
        if "ElementProjectedStress" in cell_scalars:
            return "element", cell_scalars["ElementProjectedStress"], "Element Projected Stress"
        raise KeyError(f"'ElementProjectedStress' not found. Available: {list(cell_scalars.keys())}")

    # Coordination
    if req in ("coordination", "referencecoordination"):
        if "ReferenceCoordination" in point_scalars:
            return "nodal", point_scalars["ReferenceCoordination"], "Reference Coordination"

    # Direct match in point scalars
    for k, v in point_scalars.items():
        if k.lower() == req:
            return "nodal", v, k

    # Direct match in cell scalars
    for k, v in cell_scalars.items():
        if k.lower() == req:
            return "element", v, k

    raise ValueError(
        f"Unknown field '{requested_field}'. Available point fields: {list(point_scalars.keys())}, "
        f"cell fields: {list(cell_scalars.keys())}"
    )


def render_vtk_to_png(vtk_file, out_file, args_dict):
    """
    Render a single VTK file to PNG according to specified style and options.
    """
    data = read_legacy_vtk(vtk_file)
    points = data["points"]
    cells = data["cells"]
    field_type, values, field_label = resolve_field(data, args_dict["field"])

    # Color limits
    vmin = args_dict.get("vmin")
    vmax = args_dict.get("vmax")
    if vmin is None or vmax is None:
        q_low, q_high = args_dict.get("quantiles", (0.001, 0.999))
        auto_vmin = np.nanquantile(values, q_low)
        auto_vmax = np.nanquantile(values, q_high)
        if vmin is None:
            vmin = auto_vmin
        if vmax is None:
            vmax = auto_vmax

    # Figure & Style setup
    theme = args_dict.get("theme", "white").lower()
    is_dark = theme == "dark"
    bg_color = "#1a1a1a" if is_dark else "white"
    fg_color = "#e0e0e0" if is_dark else "#222222"
    edge_default_color = "white" if is_dark else "black"

    fig, ax = plt.subplots(figsize=args_dict.get("figsize", (7.0, 6.5)), dpi=args_dict.get("dpi", 200))
    fig.patch.set_facecolor(bg_color)
    ax.set_facecolor(bg_color)

    triang = mtri.Triangulation(points[:, 0], points[:, 1], cells)
    mode = args_dict.get("mode", "colors").lower()

    # Determine what to draw
    draw_surface = mode in ("colors", "surface", "elements", "surface_with_edges")
    draw_edges = (mode in ("elements", "surface_with_edges", "wireframe")) or args_dict.get("show_elements", False)
    draw_nodes = (mode in ("nodes", "points")) or args_dict.get("show_nodes", False)

    cax_obj = None

    # 1. Surface coloring
    if draw_surface:
        cmap = args_dict.get("cmap", "turbo")
        if field_type == "nodal":
            shading = args_dict.get("shading", "gouraud")
            cax_obj = ax.tripcolor(
                triang,
                values,
                shading=shading,
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                zorder=1,
            )
        else:  # element field
            cax_obj = ax.tripcolor(
                triang,
                facecolors=values,
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                zorder=1,
            )

    # 2. Element mesh edges (Wireframe / Surface with Edges)
    if draw_edges:
        ec = args_dict.get("edge_color") or edge_default_color
        ew = args_dict.get("edge_width", 0.25)
        ea = args_dict.get("edge_alpha", 0.35)
        ax.triplot(triang, color=ec, linewidth=ew, alpha=ea, zorder=2)

    # 3. Nodes (Points / Atoms)
    if draw_nodes:
        ns = args_dict.get("node_size", 2.0)
        node_color_opt = args_dict.get("node_color", "auto")

        if node_color_opt == "field" and field_type == "nodal":
            cax_obj = ax.scatter(
                points[:, 0],
                points[:, 1],
                s=ns,
                c=values,
                cmap=args_dict.get("cmap", "turbo"),
                vmin=vmin,
                vmax=vmax,
                zorder=3,
                linewidths=0,
            )
        else:
            nc = fg_color if node_color_opt in ("auto", "default") else node_color_opt
            ax.scatter(
                points[:, 0],
                points[:, 1],
                s=ns,
                c=nc,
                zorder=3,
                linewidths=0,
                alpha=args_dict.get("node_alpha", 0.8),
            )

    # Aspect ratio & limits
    ax.set_aspect("equal")
    crop = args_dict.get("crop")
    if crop:
        xmin, xmax, ymin, ymax = crop
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
    else:
        if args_dict.get("xlim"):
            ax.set_xlim(args_dict["xlim"])
        if args_dict.get("ylim"):
            ax.set_ylim(args_dict["ylim"])

    # Colorbar
    if cax_obj is not None and not args_dict.get("no_colorbar", False):
        cbar = fig.colorbar(cax_obj, ax=ax, fraction=0.046, pad=0.04)
        cbar_label = args_dict.get("cbar_label") or field_label
        cbar.set_label(cbar_label, color=fg_color, fontsize=10)
        cbar.ax.tick_params(colors=fg_color, labelsize=9)
        cbar.outline.set_edgecolor(fg_color)

    # Title & Metadata
    if not args_dict.get("no_title", False):
        load = data["field_data"].get("LoadParameter")
        title_parts = [Path(vtk_file).name]
        if load is not None:
            title_parts.append(f"Load $\\alpha = {float(load):+.4f}$")
        title_parts.append(f"[{mode.upper()}]")
        ax.set_title("  |  ".join(title_parts), color=fg_color, fontsize=11, pad=8)

    # Clean axes option (ParaView viewport look)
    if args_dict.get("clean", False):
        ax.set_axis_off()
    else:
        ax.tick_params(colors=fg_color, labelsize=9)
        for spine in ax.spines.values():
            spine.set_color(fg_color)
            spine.set_linewidth(0.8)
        ax.set_xlabel("X", color=fg_color, fontsize=10)
        ax.set_ylabel("Y", color=fg_color, fontsize=10)

    plt.tight_layout()
    out_file = Path(out_file)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_file, dpi=args_dict.get("dpi", 200), facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)
    return out_file


def _worker(task):
    vtk_file, out_file, args_dict = task
    try:
        render_vtk_to_png(vtk_file, out_file, args_dict)
        return True, vtk_file.name, out_file.name, None
    except Exception as e:
        return False, vtk_file.name, out_file.name, str(e)


def make_movie(frame_files, output_movie_path, fps=15):
    """Stitch generated PNG frames into an MP4 video using ffmpeg."""
    if not shutil.which("ffmpeg"):
        print("Warning: ffmpeg not found in PATH; skipping movie generation.", file=sys.stderr)
        return False

    temp_list_file = output_movie_path.parent / "temp_ffmpeg_list.txt"
    try:
        with open(temp_list_file, "w") as f:
            for frame in frame_files:
                f.write(f"file '{frame.resolve()}'\n")
                f.write(f"duration {1.0 / fps:.5f}\n")
            # Repeat last frame
            if frame_files:
                f.write(f"file '{frame_files[-1].resolve()}'\n")

        cmd = [
            "ffmpeg", "-y",
            "-f", "concat", "-safe", "0",
            "-i", str(temp_list_file),
            "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
            "-c:v", "libx264",
            "-pix_fmt", "yuv420p",
            "-r", str(fps),
            str(output_movie_path)
        ]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        print(f"🎬 Movie successfully created: {output_movie_path}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error creating movie with ffmpeg: {e.stderr.decode()}", file=sys.stderr)
        return False
    finally:
        if temp_list_file.exists():
            temp_list_file.unlink()


def main():
    parser = argparse.ArgumentParser(
        description="Render 2D VTK simulation outputs to PNG (ParaView-style).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "input",
        type=Path,
        help="Input VTK file, or directory containing configuration_*.vtk files.",
    )
    parser.add_argument(
        "-o", "--out-dir",
        type=Path,
        default=None,
        help="Output directory for PNG images (default: <input_dir>/png_plots).",
    )
    parser.add_argument(
        "-f", "--field",
        default="nodal",
        choices=("nodal", "element", "energy_nodal", "energy_element", "stress_nodal", "stress_element"),
        help="Field to plot: 'nodal' (NodalEnergy) or 'element' (ElementEnergy) (default: nodal).",
    )
    parser.add_argument(
        "-m", "--mode",
        default="colors",
        choices=("colors", "surface", "elements", "surface_with_edges", "wireframe", "nodes", "points"),
        help="ParaView representation mode (default: colors / surface).",
    )
    parser.add_argument(
        "--show-elements", "--edges",
        action="store_true",
        help="Force drawing triangle element edges over the surface.",
    )
    parser.add_argument(
        "--show-nodes", "--nodes",
        action="store_true",
        help="Force drawing node points.",
    )
    parser.add_argument(
        "--crop",
        nargs=4,
        type=float,
        metavar=("XMIN", "XMAX", "YMIN", "YMAX"),
        help="Zoom in / crop to window: XMIN XMAX YMIN YMAX (e.g. --crop 50 90 50 90).",
    )
    parser.add_argument(
        "--xlim",
        nargs=2,
        type=float,
        metavar=("XMIN", "XMAX"),
        help="Fixed X axis limits for all frames (e.g. --xlim -10 310).",
    )
    parser.add_argument(
        "--ylim",
        nargs=2,
        type=float,
        metavar=("YMIN", "YMAX"),
        help="Fixed Y axis limits for all frames (e.g. --ylim -15 175).",
    )
    parser.add_argument(
        "--fixed-bounds",
        action="store_true",
        help="Automatically compute the global [Xmin, Xmax] and [Ymin, Ymax] envelope across all VTK files and fix the camera window.",
    )
    parser.add_argument(
        "--cmap",
        default="turbo",
        help="Colormap: turbo, viridis, plasma, inferno, coolwarm, jet, etc. (default: turbo).",
    )
    parser.add_argument(
        "--shading",
        default="gouraud",
        choices=("gouraud", "flat"),
        help="Shading interpolation for nodal field: 'gouraud' (smooth) or 'flat' (default: gouraud).",
    )
    parser.add_argument(
        "--theme",
        default="white",
        choices=("white", "dark"),
        help="Visual theme: 'white' (paper-ready) or 'dark' (ParaView default) (default: white).",
    )
    parser.add_argument(
        "--vmin",
        type=float,
        default=None,
        help="Manual lower color limit (default: auto quantile 0.001).",
    )
    parser.add_argument(
        "--vmax",
        type=float,
        default=None,
        help="Manual upper color limit (default: auto quantile 0.999).",
    )
    parser.add_argument(
        "--edge-color",
        default=None,
        help="Custom edge color for mesh lines (default: black in white theme, white in dark theme).",
    )
    parser.add_argument(
        "--edge-width",
        type=float,
        default=0.25,
        help="Mesh line width (default: 0.25).",
    )
    parser.add_argument(
        "--edge-alpha",
        type=float,
        default=0.35,
        help="Mesh line alpha transparency (default: 0.35).",
    )
    parser.add_argument(
        "--node-size",
        type=float,
        default=2.5,
        help="Node point size (default: 2.5).",
    )
    parser.add_argument(
        "--node-color",
        default="auto",
        help="Node color: 'auto', 'field' (color by scalar value), or color name like 'black'/'red'.",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Clean viewport style without axes, labels, and ticks (like ParaView screenshot).",
    )
    parser.add_argument(
        "--no-colorbar",
        action="store_true",
        help="Disable colorbar.",
    )
    parser.add_argument(
        "--no-title",
        action="store_true",
        help="Disable title header.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=200,
        help="Image resolution in DPI (default: 200).",
    )
    parser.add_argument(
        "--pattern",
        default="configuration_*.vtk",
        help="Filename pattern when input is a directory (default: configuration_*.vtk).",
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=1,
        help="Stride/step through files when processing a directory (e.g. --stride 5 plots every 5th file).",
    )
    parser.add_argument(
        "-j", "--jobs",
        type=int,
        default=cpu_count(),
        help=f"Number of parallel jobs for batch processing (default: all CPUs = {cpu_count()}).",
    )
    parser.add_argument(
        "--movie",
        action="store_true",
        help="Compile rendered frames into an MP4 video (requires ffmpeg).",
    )
    parser.add_argument(
        "--fps",
        type=int,
        default=15,
        help="Frames per second for movie output (default: 15).",
    )

    args = parser.parse_args()

    # Collect files
    input_path = args.input
    if input_path.is_file():
        files = [input_path]
        out_dir = args.out_dir or input_path.parent / "png_plots"
    elif input_path.is_dir():
        files = sorted(input_path.glob(args.pattern))
        if not files:
            print(f"Error: No files matching '{args.pattern}' found in {input_path}")
            sys.exit(1)
        if args.stride > 1:
            files = files[::args.stride]
        out_dir = args.out_dir or input_path / "png_plots"
    else:
        print(f"Error: Input path '{input_path}' does not exist.")
        sys.exit(1)

    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Processing {len(files)} file(s)...")
    print(f"  Field: {args.field}")
    print(f"  Representation: {args.mode}")
    print(f"  Output directory: {out_dir}")

    args_dict = dict(vars(args))

    # Auto-detect global bounds if requested
    if args.fixed_bounds and len(files) > 1:
        print("Computing global bounding box envelope across files...")
        glob_xmin, glob_xmax = float("inf"), float("-inf")
        glob_ymin, glob_ymax = float("inf"), float("-inf")
        sample_indices = sorted(list(set([0, len(files)//4, len(files)//2, 3*len(files)//4, len(files)-1])))
        for idx in sample_indices:
            d = read_legacy_vtk(files[idx])
            pts = d["points"]
            glob_xmin = min(glob_xmin, float(pts[:, 0].min()))
            glob_xmax = max(glob_xmax, float(pts[:, 0].max()))
            glob_ymin = min(glob_ymin, float(pts[:, 1].min()))
            glob_ymax = max(glob_ymax, float(pts[:, 1].max()))
        margin_x = (glob_xmax - glob_xmin) * 0.03
        margin_y = (glob_ymax - glob_ymin) * 0.03
        args_dict["xlim"] = [glob_xmin - margin_x, glob_xmax + margin_x]
        args_dict["ylim"] = [glob_ymin - margin_y, glob_ymax + margin_y]
        print(f"  Fixed camera bounds: X in [{args_dict['xlim'][0]:.1f}, {args_dict['xlim'][1]:.1f}], Y in [{args_dict['ylim'][0]:.1f}, {args_dict['ylim'][1]:.1f}]")

    # Build tasks
    field_tag = args.field.replace("energy_", "").replace("_energy", "")
    tasks = []
    ordered_out_files = []
    for f in files:
        stem = f.stem
        out_name = f"{stem}_{field_tag}_{args.mode}.png"
        out_file = out_dir / out_name
        tasks.append((f, out_file, args_dict))
        ordered_out_files.append(out_file)

    # Run
    if len(tasks) == 1:
        success, in_n, out_n, err = _worker(tasks[0])
        if success:
            print(f"✓ Saved: {out_dir / out_n}")
        else:
            print(f"✗ Failed {in_n}: {err}", file=sys.stderr)
            sys.exit(1)
    else:
        num_jobs = min(args.jobs, len(tasks))
        print(f"Rendering with {num_jobs} workers...")
        with Pool(processes=num_jobs) as pool:
            completed = 0
            for success, in_n, out_n, err in pool.imap_unordered(_worker, tasks):
                completed += 1
                if success:
                    if completed % 10 == 0 or completed == len(tasks):
                        print(f"  [{completed}/{len(tasks)}] Saved {out_n}")
                else:
                    print(f"  [ERROR] {in_n}: {err}", file=sys.stderr)

        print(f"✓ All {len(tasks)} images saved to: {out_dir}")

        if args.movie:
            movie_name = f"movie_{field_tag}_{args.mode}.mp4"
            movie_path = out_dir / movie_name
            make_movie(ordered_out_files, movie_path, fps=args.fps)


if __name__ == "__main__":
    main()
