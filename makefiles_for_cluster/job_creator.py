#!/usr/bin/env python3
"""
Automatic Job Setup Generator for Parallel Crystal Plasticity Simulations (Magi Cluster)
Creates SLURM script, parallel launcher, and individual run scripts
adapted for the lattice_triangulation interface:
    lattice_triangulation <nx> <ny> <mode> <seed> [remesh] [--precond=... --precond-from-step=N] [--alpha-end=... --alpha-start=... --step-size=...]
"""

import argparse
import os
import sys
from pathlib import Path


def get_input(prompt, default=None):
    """Get user input with optional default; uses default if non-interactive"""
    if not sys.stdin.isatty() and default is not None:
        return str(default)
    if default is not None:
        try:
            user_input = input(f"{prompt} [{default}]: ").strip()
            return user_input if user_input else str(default)
        except (EOFError, KeyboardInterrupt):
            return str(default)
    try:
        return input(f"{prompt}: ").strip()
    except (EOFError, KeyboardInterrupt):
        return ""


def resolve_val(cli_val, prompt, default, converter=str, auto_yes=False):
    """Use CLI value if specified; otherwise prompt interactively (or use default if auto_yes)"""
    if cli_val is not None:
        return converter(cli_val)
    if auto_yes:
        return converter(default)
    return converter(get_input(prompt, default))


def main():
    parser = argparse.ArgumentParser(
        description="Parallel Job Setup Generator for Crystal Plasticity Simulations (Magi Cluster)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 1. Hardware & Cores
    parser.add_argument("--n-jobs", type=int, default=None, help="Number of parallel jobs to run simultaneously")
    parser.add_argument("--threads", "--threads-per-job", dest="threads_per_job", type=int, default=None, help="OpenMP threads per job")

    # 2. Physics & System Parameters
    parser.add_argument("--nx", type=int, default=None, help="System size nx")
    parser.add_argument("--ny", type=int, default=None, help="System size ny")
    parser.add_argument("--mode", choices=["both", "positive", "negative"], default=None, help="Loading mode ('both', 'positive', 'negative')")
    parser.add_argument("--seed", "--start-seed", dest="start_seed", type=int, default=None, help="Starting random seed")
    parser.add_argument("--seed-strategy", choices=["unique", "paired"], default=None, help="Seed strategy for 'both' mode")
    parser.add_argument("--remesh", dest="remesh", action="store_true", default=None, help="Enable adaptive remeshing")
    parser.add_argument("--no-remesh", dest="remesh", action="store_false", help="Disable adaptive remeshing")

    # Loading schedule parameters
    parser.add_argument("--alpha-end", "--alpha-max", dest="alpha_end", type=float, default=None, help="Target shear strain magnitude |alpha_end| (e.g. 1.0 or 3.0)")
    parser.add_argument("--alpha-start", "--alpha-min", dest="alpha_start", type=float, default=None, help="Initial shear strain magnitude |alpha_start| (default: 0.14)")
    parser.add_argument("--step-size", type=float, default=None, help="Step size magnitude |d_alpha| (default: 6e-5)")

    # Data saving & checkpoint controls
    parser.add_argument("--checkpoint-interval", "--chk-interval", dest="checkpoint_interval", type=int, default=None, help="Periodic elastic checkpoint interval in steps (0 = disabled, default: 500)")
    parser.add_argument("--stress-drop-threshold", "--stress-drop", dest="stress_drop_threshold", type=float, default=None, help="Fractional stress drop threshold to trigger avalanche save (default: 0.10)")
    parser.add_argument("--save-triangle-data", dest="save_triangle_data", action="store_true", default=False, help="Enable legacy triangle_data output (default: disabled)")
    parser.add_argument("--no-perturbation", "--no-perturb", dest="no_perturbation", action="store_true", default=False, help="Disable triangulation perturbation for negative loading (keep same diagonal as positive)")
    parser.add_argument("--perturbation", type=float, default=None, help="Custom triangulation perturbation value (default: -1e-7 for negative loading)")

    # Preconditioner
    parser.add_argument("--precond", choices=["stiffness", "laplacian", "diag", "none"], default=None, help="L-BFGS preconditioner")
    parser.add_argument("--precond-from-step", type=int, default=None, help="Use preconditioned solver from step N")

    # Stability monitor
    parser.add_argument("--eig-every", type=int, default=None, help="Stability monitor: eigenvalues every N steps (0 = off)")
    parser.add_argument("--eig-at-avalanche", type=int, choices=[0, 1], default=None, help="Eigenvalues before/after avalanches (0 or 1)")

    # 3. Paths & Cluster
    parser.add_argument("--base-dir", default=None, help="Base directory for MTM on cluster")
    parser.add_argument("--exe-path", default=None, help="Executable path")
    parser.add_argument("--runs-dir", default=None, help="Runs output directory")

    # 4. SLURM Options
    parser.add_argument("--job-name", default=None, help="SLURM job name")
    parser.add_argument("--partition", default=None, help="SLURM partition")
    parser.add_argument("--email", default=None, help="Email for notifications")
    parser.add_argument("-y", "--yes", action="store_true", help="Auto-confirm script generation without interactive prompt")

    args = parser.parse_args()

    print("=" * 70)
    print("  Parallel Job Setup Generator (Lattice Triangulation - Magi Cluster)")
    print("=" * 70)
    print()

    # 1. Resource Configuration
    print("--- 1. Hardware & Core Allocation ---")
    n_jobs = resolve_val(args.n_jobs, "Number of parallel jobs to run simultaneously", "4", int, args.yes)
    threads_per_job = resolve_val(args.threads_per_job, "OpenMP threads per job", "32", int, args.yes)
    total_cores = n_jobs * threads_per_job

    print(f"  -> Total cores allocated: {total_cores} ({n_jobs} jobs x {threads_per_job} threads)")
    print()

    # 2. Physics & Simulation Parameters
    print("--- 2. Simulation Parameters ---")
    nx = resolve_val(args.nx, "System size nx", "150", int, args.yes)
    ny = resolve_val(args.ny, "System size ny", "150", int, args.yes)
    
    if args.mode is not None:
        mode_choice = args.mode.lower()
    elif args.yes:
        mode_choice = "both"
    else:
        print("Loading mode options:")
        print("  'both'     : Launch paired positive & negative runs with matching seeds")
        print("  'positive' : All jobs run positive loading with sequential seeds")
        print("  'negative' : All jobs run negative loading with sequential seeds")
        mode_choice = get_input("Select loading mode ('both', 'positive', 'negative')", "both").lower()
        if mode_choice not in ("both", "positive", "negative"):
            print(f"Warning: Unknown mode '{mode_choice}', defaulting to 'both'.")
            mode_choice = "both"

    start_seed = resolve_val(args.start_seed, "Starting random seed", "42", int, args.yes)
    seed_strategy = "unique"
    if mode_choice == "both":
        if args.seed_strategy is not None:
            seed_strategy = args.seed_strategy.lower()
        elif args.yes:
            seed_strategy = "unique"
        else:
            seed_strategy = get_input(
                "Seed strategy for 'both' mode: 'unique' (44, 45, 46, 47) or 'paired' (44, 44, 45, 45)",
                "unique"
            ).lower()

    # Remeshing configuration
    if args.remesh is not None:
        enable_remeshing = args.remesh
    elif args.yes:
        enable_remeshing = True
    else:
        remesh_choice = get_input("Enable adaptive remeshing? ('yes' / 'no')", "yes").lower()
        enable_remeshing = remesh_choice in ("y", "yes", "true", "1", "remesh")
    remesh_arg = "1" if enable_remeshing else "0"
    remesh_status = "enabled" if enable_remeshing else "disabled"
    suffix = "" if enable_remeshing else "_noremesh"
    print(f"  -> Remeshing: {remesh_status.upper()}")

    # Loading schedule parameters (alpha_start, alpha_end, step_size)
    print("Loading schedule parameters:")
    alpha_end = resolve_val(args.alpha_end, "Target shear strain magnitude |alpha_end|", "1.0", float, args.yes)
    alpha_start = resolve_val(args.alpha_start, "Initial shear strain magnitude |alpha_start|", "0.14", float, args.yes)
    step_size = resolve_val(args.step_size, "Step size magnitude |d_alpha|", "6e-5", float, args.yes)

    loading_flags_list = []
    if alpha_end != 1.0:
        loading_flags_list.append(f"--alpha-end={alpha_end:g}")
        suffix += f"_alpha{alpha_end:g}"
    if alpha_start != 0.14:
        loading_flags_list.append(f"--alpha-start={alpha_start:g}")
    if step_size != 6e-5:
        loading_flags_list.append(f"--step-size={step_size:g}")
    if args.checkpoint_interval is not None:
        loading_flags_list.append(f"--checkpoint-interval={args.checkpoint_interval}")
    if args.stress_drop_threshold is not None:
        loading_flags_list.append(f"--stress-drop-threshold={args.stress_drop_threshold:g}")
    if args.save_triangle_data:
        loading_flags_list.append("--save-triangle-data")
    if args.no_perturbation:
        loading_flags_list.append("--no-perturbation")
    elif args.perturbation is not None:
        loading_flags_list.append(f"--perturbation={args.perturbation:g}")
    loading_flags = " ".join(loading_flags_list)
    loading_status = f"|alpha| = {alpha_start:g} -> {alpha_end:g} (step {step_size:g})"
    print(f"  -> Loading schedule: {loading_status}")

    # L-BFGS preconditioner (see README "Preconditioned L-BFGS")
    if args.precond is not None:
        precond = args.precond.lower()
    elif args.yes:
        precond = "stiffness"
    else:
        print("L-BFGS preconditioner options:")
        print("  'stiffness' : analytic FEM stiffness, sparse Cholesky (fastest, ~10-15x on elastic steps)")
        print("  'laplacian' : reference Laplacian (~3x)")
        print("  'diag'      : ALGLIB diagonal preconditioner (little gain)")
        print("  'none'      : original plain L-BFGS")
        precond = get_input("Select preconditioner", "stiffness").lower()
    if precond in ("off", "no", "plain", "0", "false"):
        precond = "none"
    if precond not in ("stiffness", "laplacian", "diag", "none"):
        print(f"Warning: Unknown preconditioner '{precond}', defaulting to 'stiffness'.")
        precond = "stiffness"
    precond_flags = ""
    precond_status = "none (plain L-BFGS)"
    if precond != "none":
        from_step = resolve_val(
            args.precond_from_step,
            "Use plain L-BFGS for load steps before (1 = initial relaxation identical to plain runs)",
            "1",
            int,
            args.yes
        )
        precond_flags = f"--precond={precond} --precond-from-step={from_step}"
        precond_status = f"{precond} (from step {from_step})"
    # Preconditioned runs are the default; tag folders/jobs of the other choices
    precond_tag = "" if precond == "stiffness" else ("_plainlbfgs" if precond == "none" else f"_{precond}")
    suffix += precond_tag
    print(f"  -> Preconditioner: {precond_status}")

    # Stability monitor: lowest stiffness eigenvalues during the run
    eig_every = resolve_val(args.eig_every, "Stability monitor: eigenvalues every N load steps (0 = off)", "0", int, args.yes)
    eig_flags = ""
    eig_status = "off"
    if eig_every > 0:
        if args.eig_at_avalanche is not None:
            eig_at_avalanche = args.eig_at_avalanche
        elif args.yes:
            eig_at_avalanche = 1
        else:
            eig_aval = get_input("  Also before/after each avalanche, with soft modes? ('yes' / 'no')", "yes").lower()
            eig_at_avalanche = 1 if eig_aval in ("y", "yes", "true", "1") else 0
        eig_flags = f"--eig-every={eig_every} --eig-at-avalanche={eig_at_avalanche}"
        eig_status = f"every {eig_every} steps" + (", before/after avalanches" if eig_at_avalanche else "")
    print(f"  -> Stability monitor: {eig_status}")
    print()

    # 3. Paths and Directories
    print("--- 3. Cluster Paths ---")
    current_cwd = os.getcwd()
    default_base = current_cwd if "MTM" in current_cwd else "/home/dist/umut.salman/latest_MTM"
    base_dir = resolve_val(args.base_dir, "Base directory for MTM", default_base, str, args.yes)
    base_dir = os.path.abspath(os.path.expanduser(base_dir))
    exe_path = resolve_val(args.exe_path, "Executable path", f"{base_dir}/lattice_triangulation", str, args.yes)
    exe_path = os.path.abspath(os.path.expanduser(exe_path))
    
    default_runs_dir = f"{base_dir}/runs_{nx}x{ny}{suffix}"
    runs_dir = resolve_val(args.runs_dir, "Runs output directory", default_runs_dir, str, args.yes)

    # 4. SLURM Options
    print()
    print("--- 4. SLURM Cluster Options ---")
    default_job_name = f"shear_{nx}x{ny}{suffix}"
    job_name = resolve_val(args.job_name, "SLURM job name", default_job_name, str, args.yes)
    partition = resolve_val(args.partition, "SLURM partition", "COMPUTE2", str, args.yes)
    email = resolve_val(args.email, "Email for notifications", "umut.salman@lspm.cnrs.fr", str, args.yes)

    # Build per-job specifications: [(job_id, nx, ny, mode, seed, remesh_arg, run_folder)]
    job_specs = []
    for i in range(n_jobs):
        job_id = i + 1
        if mode_choice == "both":
            if seed_strategy.startswith("p"):
                seed = start_seed + (i // 2)
            else:
                seed = start_seed + i
            mode = "positive" if (i % 2 == 0) else "negative"
        elif mode_choice == "positive":
            seed = start_seed + i
            mode = "positive"
        else:  # negative
            seed = start_seed + i
            mode = "negative"
        
        folder_name = f"run_{job_id:02d}_{mode}_seed{seed}{suffix}"
        job_specs.append((job_id, nx, ny, mode, seed, remesh_arg, folder_name))

    print()
    print("=" * 70)
    print("Planned Job Configurations:")
    print(f"  Loading schedule: {loading_status}")
    print(f"  Preconditioner:   {precond_status}")
    print(f"  Stability monitor:{eig_status}")
    print("=" * 70)
    for j_id, j_nx, j_ny, j_mode, j_seed, j_remesh, j_dir in job_specs:
        start_cpu = (j_id - 1) * threads_per_job
        end_cpu = start_cpu + threads_per_job - 1
        remesh_desc = "remesh" if j_remesh == "1" else "no-remesh"
        print(f"  Job {j_id:2d}: CPUs {start_cpu:3d}-{end_cpu:3d} | {j_nx}x{j_ny} | mode={j_mode:8s} | seed={j_seed:<4d} | {remesh_desc:9s} | dir={j_dir}")
    print("=" * 70)
    print()

    if not args.yes:
        confirm = get_input("Generate scripts with these parameters? (y/n)", "y").lower()
        if confirm not in ("y", "yes"):
            print("Aborted.")
            sys.exit(0)

    print()
    print("Generating files...")

    # ==============================================================
    # 1. Generate run_one.sh
    # ==============================================================
    run_one_script = f"""#!/bin/bash
# Worker script for individual simulation run
# Usage: ./run_one.sh <RUN_ID> <NX> <NY> <MODE> <SEED> <REMESH> <RUN_DIR>

RUN_ID=$1
NX=$2
NY=$3
MODE=$4
SEED=$5
REMESH=$6
RUN_DIR=$7
EXE="{exe_path}"
PRECOND_FLAGS="{precond_flags}"
EIG_FLAGS="{eig_flags}"
LOADING_FLAGS="{loading_flags}"

# Create run directory and enter it
mkdir -p "$RUN_DIR"
cd "$RUN_DIR" || exit 1

# OpenMP thread settings
export OMP_NUM_THREADS={threads_per_job}
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

# Log header
echo "=========================================================="
echo "Run $RUN_ID starting at $(date)"
echo "  Working directory: $RUN_DIR"
echo "  Parameters:        size=${{NX}}x${{NY}}, mode=${{MODE}}, seed=${{SEED}}, remesh=${{REMESH}}"
echo "  CPU affinity:      $(taskset -cp $$ 2>/dev/null || echo 'N/A')"
echo "  OpenMP threads:    $OMP_NUM_THREADS"
echo "  Executable:        $EXE"
echo "  Preconditioner:    ${{PRECOND_FLAGS:-none (plain L-BFGS)}}"
echo "  Stability monitor: ${{EIG_FLAGS:-off}}"
echo "  Loading schedule:  ${{LOADING_FLAGS:-default (|alpha|=0.14->1.0, step 6e-5)}}"
echo "=========================================================="

# Execute simulation
"$EXE" "$NX" "$NY" "$MODE" "$SEED" "$REMESH" $PRECOND_FLAGS $EIG_FLAGS $LOADING_FLAGS > simulation.log 2>&1
EXIT_CODE=$?

echo ""
echo "=========================================================="
echo "Run $RUN_ID finished at $(date) with exit code: $EXIT_CODE"
echo "=========================================================="

exit $EXIT_CODE
"""
    run_one_file = f"{base_dir}/run_{job_name}_one.sh"
    with open(run_one_file, "w") as f:
        f.write(run_one_script)
    os.chmod(run_one_file, 0o755)
    print(f"✓ Created: {run_one_file}")

    # ==============================================================
    # 2. Generate parallel launcher script
    # ==============================================================
    parallel_script = f"""#!/bin/bash
# Parallel launcher script: dynamically pins each job to dedicated CPU cores
SCRIPT_DIR="{base_dir}"
RUNS_BASE="{runs_dir}"
N_JOBS={n_jobs}
THREADS_PER_JOB={threads_per_job}

mkdir -p "$RUNS_BASE"

echo "=========================================================="
echo "Starting {n_jobs} parallel jobs with dynamic CPU pinning..."
echo "Runs directory: $RUNS_BASE"
echo "=========================================================="
echo ""

# Dynamically determine the exact CPUs assigned to this SLURM job/process
CPU_ASSIGNMENTS=($(python3 -c "
import os

cpus = []
if hasattr(os, 'sched_getaffinity'):
    try:
        cpus = sorted(list(os.sched_getaffinity(0)))
    except Exception:
        pass

if not cpus:
    try:
        with open('/proc/self/status') as f:
            for line in f:
                if line.startswith('Cpus_allowed_list:'):
                    for part in line.split(':')[1].strip().split(','):
                        if '-' in part:
                            s, e = map(int, part.split('-'))
                            cpus.extend(range(s, e + 1))
                        elif part:
                            cpus.append(int(part))
                    cpus = sorted(cpus)
                    break
    except Exception:
        pass

cpus = sorted(list(dict.fromkeys(cpus)))
n_jobs = $N_JOBS
threads = $THREADS_PER_JOB

if len(cpus) >= n_jobs * threads:
    for i in range(n_jobs):
        job_cpus = cpus[i*threads : (i+1)*threads]
        min_c, max_c = min(job_cpus), max(job_cpus)
        if job_cpus == list(range(min_c, max_c + 1)):
            print(f'{{min_c}}-{{max_c}}')
        else:
            print(','.join(map(str, job_cpus)))
else:
    for i in range(n_jobs):
        print('auto')
" 2>/dev/null))

"""
    for j_id, j_nx, j_ny, j_mode, j_seed, j_remesh, j_dir in job_specs:
        idx = j_id - 1
        remesh_desc = "remesh" if j_remesh == "1" else "no-remesh"
        parallel_script += f"""
CPU_SPEC="${{CPU_ASSIGNMENTS[{idx}]}}"
if [ -n "$CPU_SPEC" ] && [ "$CPU_SPEC" != "auto" ]; then
    TASKSET_CMD="taskset -c $CPU_SPEC"
    echo "Launching Job {j_id}: mode={j_mode}, seed={j_seed}, remesh={remesh_desc} pinned to CPUs $CPU_SPEC..."
else
    TASKSET_CMD=""
    echo "Launching Job {j_id}: mode={j_mode}, seed={j_seed}, remesh={remesh_desc} (automatic affinity)..."
fi

$TASKSET_CMD "${{SCRIPT_DIR}}/run_{job_name}_one.sh" {j_id} {j_nx} {j_ny} {j_mode} {j_seed} {j_remesh} "${{RUNS_BASE}}/{j_dir}" &
"""

    parallel_script += """
# Wait for all parallel background jobs to finish
echo ""
echo "All jobs launched in background. Waiting for completion..."
wait

echo ""
echo "=========================================================="
echo "All parallel jobs completed at $(date)!"
echo "=========================================================="
"""
    parallel_file = f"{base_dir}/run_{job_name}_parallel.sh"
    with open(parallel_file, "w") as f:
        f.write(parallel_script)
    os.chmod(parallel_file, 0o755)
    print(f"✓ Created: {parallel_file}")

    # ==============================================================
    # 3. Generate SLURM submission script
    # ==============================================================
    slurm_script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={job_name}.out
#SBATCH --error={job_name}.err
#SBATCH --mail-type=end
#SBATCH --mail-user={email}
#SBATCH --partition={partition}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task={total_cores}
#SBATCH --hint=nomultithread

export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

echo "=========================================="
echo "Job:            {job_name}"
echo "Node:           $(hostname)"
echo "Date:           $(date)"
echo "Available CPUs: $(taskset -cp $$ 2>/dev/null || echo 'N/A')"
echo "=========================================="
echo ""
echo "Configuration:"
echo "  Parallel jobs:   {n_jobs}"
echo "  Threads per job: {threads_per_job}"
echo "  Total cores:     {total_cores}"
echo "  System size:     {nx}x{ny}"
echo "  Mode setting:    {mode_choice}"
echo "  Start seed:      {start_seed}"
echo "  Remeshing:       {remesh_status}"
echo "  Loading range:   {loading_status}"
echo "  Preconditioner:  {precond_status}"
echo "  Stability mon.:  {eig_status}"
echo "=========================================="
echo ""

# Launch parallel jobs
"{parallel_file}"

echo ""
echo "=========================================="
echo "SLURM job completed at $(date)"
echo "=========================================="
"""
    slurm_file = f"{base_dir}/run_{job_name}.sh"
    with open(slurm_file, "w") as f:
        f.write(slurm_script)
    os.chmod(slurm_file, 0o755)
    print(f"✓ Created: {slurm_file}")

    # ==============================================================
    # Summary & Submission Instructions
    # ==============================================================
    print()
    print("=" * 70)
    print("Setup Complete!")
    print("=" * 70)
    print()
    print("Files created:")
    print(f"  1. SLURM batch script:     {slurm_file}")
    print(f"  2. Parallel launcher:      {parallel_file}")
    print(f"  3. Single simulation run:  {run_one_file}")
    print()
    print("To submit this job to the cluster SLURM queue:")
    print(f"  cd {base_dir}")
    print(f"  sbatch run_{job_name}.sh")
    print()
    print("To monitor progress in real time:")
    print(f"  squeue -u $USER")
    print(f"  tail -f {job_name}.out")
    print(f"  tail -f {runs_dir}/*/simulation.log")
    print()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\nCancelled by user.")
        sys.exit(1)
    except Exception as e:
        print(f"\nError: {e}")
        sys.exit(1)
