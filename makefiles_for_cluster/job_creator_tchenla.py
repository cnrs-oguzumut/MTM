#!/usr/bin/env python3
"""
Automatic Job Setup Generator for Parallel Simulations on Tchenla Cluster
Paris 13 University (Université Sorbonne Paris Nord / LSPM)

Hardware on Tchenla:
  c[1-2]: 2x Intel Xeon-Gold 6430 (64 cores, 128 threads, 512GB RAM)
  c[3-4]: 2x Intel Xeon-Gold 6430 (64 cores, 128 threads, 256GB RAM)
  c[5-7]: 2x Intel Xeon E5-2630 v2 (12 cores, 24 threads, 125GB RAM)
  c8:     2x Intel Xeon E5-2650 v2 (16 cores, 32 threads, 64GB RAM)
"""

import os
import sys
from pathlib import Path


def get_input(prompt, default=None):
    """Get user input with optional default"""
    if default is not None:
        user_input = input(f"{prompt} [{default}]: ").strip()
        return user_input if user_input else str(default)
    return input(f"{prompt}: ").strip()


def main():
    print("=" * 70)
    print("  Tchenla Cluster - Parallel Job Generator (Lattice Triangulation)")
    print("=" * 70)
    print()

    # 1. Resource Configuration
    print("--- 1. Hardware & Core Allocation ---")
    print("Tchenla modern compute nodes c[1-4] have 64 physical cores (Intel Xeon-Gold 6430).")
    print("Typical allocations:")
    print("  - 1 job  x 64 threads = 64 cores (full node)")
    print("  - 2 jobs x 32 threads = 64 cores")
    print("  - 4 jobs x 16 threads = 64 cores")
    n_jobs = int(get_input("Number of parallel jobs to run simultaneously", "2"))
    threads_per_job = int(get_input("OpenMP threads per job", "32"))
    total_cores = n_jobs * threads_per_job

    print(f"  -> Total cores allocated: {total_cores} ({n_jobs} jobs x {threads_per_job} threads)")
    print()

    # 2. Simulation Parameters
    print("--- 2. Simulation Parameters ---")
    nx = int(get_input("System size nx", "150"))
    ny = int(get_input("System size ny", "150"))
    
    print("Loading mode options:")
    print("  'both'     : Launch paired positive & negative runs with matching seeds")
    print("  'positive' : All jobs run positive loading with sequential seeds")
    print("  'negative' : All jobs run negative loading with sequential seeds")
    mode_choice = get_input("Select loading mode ('both', 'positive', 'negative')", "both").lower()
    if mode_choice not in ("both", "positive", "negative"):
        print(f"Warning: Unknown mode '{mode_choice}', defaulting to 'both'.")
        mode_choice = "both"

    start_seed = int(get_input("Starting random seed", "42"))
    seed_strategy = "unique"
    if mode_choice == "both":
        seed_strategy = get_input(
            "Seed strategy for 'both' mode: 'unique' (42, 43, 44, 45) or 'paired' (42, 42, 43, 43)",
            "unique"
        ).lower()

    # Remeshing configuration
    remesh_choice = get_input("Enable adaptive remeshing? ('yes' / 'no')", "yes").lower()
    enable_remeshing = remesh_choice in ("y", "yes", "true", "1", "remesh")
    remesh_arg = "1" if enable_remeshing else "0"
    remesh_status = "enabled" if enable_remeshing else "disabled"
    suffix = "" if enable_remeshing else "_noremesh"
    print(f"  -> Remeshing: {remesh_status.upper()}")

    # L-BFGS preconditioner (see README "Preconditioned L-BFGS")
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
        from_step = int(get_input(
            "Use plain L-BFGS for load steps before (1 = initial relaxation identical to plain runs)",
            "1"))
        precond_flags = f"--precond={precond} --precond-from-step={from_step}"
        precond_status = f"{precond} (from step {from_step})"
    # Preconditioned runs are the default; tag folders/jobs of the other choices
    precond_tag = "" if precond == "stiffness" else ("_plainlbfgs" if precond == "none" else f"_{precond}")
    suffix += precond_tag
    print(f"  -> Preconditioner: {precond_status}")

    # Stability monitor: lowest stiffness eigenvalues during the run (eigen_log.csv,
    # eigen_modes/); about +7% run time with N = 5
    eig_every = int(get_input(
        "Stability monitor: eigenvalues every N load steps (0 = off)", "5"))
    eig_flags = ""
    eig_status = "off"
    if eig_every > 0:
        eig_aval = get_input(
            "  Also before/after each avalanche, with soft modes? ('yes' / 'no')", "yes").lower()
        eig_at_avalanche = 1 if eig_aval in ("y", "yes", "true", "1") else 0
        eig_flags = f"--eig-every={eig_every} --eig-at-avalanche={eig_at_avalanche}"
        eig_status = f"every {eig_every} steps" + (", before/after avalanches" if eig_at_avalanche else "")
    print(f"  -> Stability monitor: {eig_status}")
    print()

    # 3. Paths and Directories
    print("--- 3. Tchenla Paths ---")
    current_cwd = os.getcwd()
    default_base = current_cwd if "MTM" in current_cwd else "/home/usalman/MTM"
    base_dir = get_input("Base directory for MTM", default_base)
    exe_path = get_input("Executable path", f"{base_dir}/lattice_triangulation")
    
    default_runs_dir = f"{base_dir}/runs_{nx}x{ny}{suffix}"
    runs_dir = get_input("Runs output directory", default_runs_dir)

    # 4. SLURM Options
    print()
    print("--- 4. SLURM Options ---")
    default_job_name = f"shear_{nx}x{ny}{suffix}"
    job_name = get_input("SLURM job name", default_job_name)
    nodelist = get_input("Target specific nodes (e.g. 'c[1-4]' for 64-core nodes, or leave blank for any)", "c[1-4]").strip()

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
            
        remesh_tag = "remesh" if enable_remeshing else "noremesh"
        folder = f"{runs_dir}/run_{job_id:02d}_{mode}_seed{seed}_{remesh_tag}{precond_tag}"
        job_specs.append((job_id, nx, ny, mode, seed, remesh_arg, folder))

    # Review Configuration
    print()
    print("=" * 70)
    print("Configuration Summary:")
    print("=" * 70)
    print(f"  Base Directory:    {base_dir}")
    print(f"  Executable:        {exe_path}")
    print(f"  Runs Directory:    {runs_dir}")
    print(f"  System Size:       {nx}x{ny}")
    print(f"  Remeshing:         {remesh_status.upper()}")
    print(f"  Preconditioner:    {precond_status}")
    print(f"  Parallel Jobs:     {n_jobs}")
    print(f"  Threads Per Job:   {threads_per_job}")
    print(f"  Total Cores:       {total_cores}")
    if nodelist:
        print(f"  Target Node(s):    {nodelist}")
    print()
    print("Jobs to be launched:")
    for job_id, j_nx, j_ny, j_mode, j_seed, j_remesh, j_folder in job_specs:
        r_txt = "remesh" if j_remesh == "1" else "no-remesh"
        print(f"  Job {job_id:02d}: {j_mode:<8} | seed={j_seed} | {r_txt:<9} -> {Path(j_folder).name}")
    print("=" * 70)
    print()

    confirm = get_input("Proceed with file generation? (y/n)", "y").lower()
    if confirm not in ("y", "yes"):
        print("Aborted.")
        sys.exit(0)

    print()
    print("Generating files...")

    # ==============================================================
    # 1. Generate run_one.sh (Worker script)
    # ==============================================================
    run_one_script = f"""#!/bin/bash
# Worker script for individual simulation run on Tchenla
# Usage: ./run_{job_name}_one.sh <RUN_ID> <NX> <NY> <MODE> <SEED> <REMESH> <RUN_DIR>

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

# Ensure stack limit and environment modules are available
ulimit -s unlimited
if [ -f "{base_dir}/load_modules_tchenla.sh" ]; then
    source "{base_dir}/load_modules_tchenla.sh" > /dev/null 2>&1
fi

mkdir -p "$RUN_DIR"
cd "$RUN_DIR" || exit 1

# OpenMP thread settings
export OMP_NUM_THREADS={threads_per_job}
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export OMP_PROC_BIND=close
export OMP_PLACES=cores

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
echo "=========================================================="

# Execute simulation
"$EXE" "$NX" "$NY" "$MODE" "$SEED" "$REMESH" $PRECOND_FLAGS $EIG_FLAGS > simulation.log 2>&1
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
# Parallel launcher script for Tchenla: pins each job to dedicated CPU cores
SCRIPT_DIR="{base_dir}"
RUNS_BASE="{runs_dir}"
N_JOBS={n_jobs}
THREADS_PER_JOB={threads_per_job}

mkdir -p "$RUNS_BASE"

echo "=========================================================="
echo "Starting {n_jobs} parallel jobs with CPU pinning..."
echo "Runs directory: $RUNS_BASE"
echo "=========================================================="
echo ""

# Dynamically determine CPUs assigned to this SLURM job
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

if not cpus:
    cpus = list(range({total_cores}))

print(' '.join(map(str, cpus)))
"))

TOTAL_CPUS=${{#CPU_ASSIGNMENTS[@]}}
echo "Detected $TOTAL_CPUS assigned CPU cores."

PIDS=()
"""
    # Append job launches
    for i, (job_id, j_nx, j_ny, j_mode, j_seed, j_remesh, j_folder) in enumerate(job_specs):
        parallel_script += f"""
# Job {job_id:02d}: {j_mode} (seed {j_seed})
START_IDX=$(( {i} * THREADS_PER_JOB ))
END_IDX=$(( START_IDX + THREADS_PER_JOB - 1 ))

if [ $END_IDX -lt $TOTAL_CPUS ]; then
    CORE_LIST=$(IFS=,; echo "${{CPU_ASSIGNMENTS[*]:START_IDX:THREADS_PER_JOB}}")
    echo "Launching Job {job_id:02d} on cores: $CORE_LIST"
    taskset -c "$CORE_LIST" "{run_one_file}" \\
        {job_id} {j_nx} {j_ny} "{j_mode}" {j_seed} {j_remesh} "{j_folder}" &
    PIDS+=($!)
else
    echo "Warning: Not enough cores for Job {job_id:02d}. Launching without taskset affinity."
    "{run_one_file}" \\
        {job_id} {j_nx} {j_ny} "{j_mode}" {j_seed} {j_remesh} "{j_folder}" &
    PIDS+=($!)
fi
"""

    parallel_script += """
echo ""
echo "All jobs launched in background. Waiting for completion..."
echo ""

# Wait for all background jobs to complete
FAILURES=0
for pid in "${PIDS[@]}"; do
    wait "$pid" || FAILURES=$((FAILURES + 1))
done

echo ""
echo "=========================================================="
if [ $FAILURES -eq 0 ]; then
    echo "All parallel runs completed successfully!"
else
    echo "Completed with $FAILURES failed job(s)."
fi
echo "=========================================================="
exit $FAILURES
"""
    parallel_file = f"{base_dir}/launch_{job_name}_parallel.sh"
    with open(parallel_file, "w") as f:
        f.write(parallel_script)
    os.chmod(parallel_file, 0o755)
    print(f"✓ Created: {parallel_file}")

    # ==============================================================
    # 3. Generate SLURM submission script
    # ==============================================================
    nodelist_line = f"#SBATCH --nodelist={nodelist}" if nodelist else "# (no specific nodelist requested)"
    slurm_script = f"""#!/bin/bash
# ==============================================================================
# SLURM Job Script for Tchenla Cluster (Paris 13 University / LSPM)
# ==============================================================================
#SBATCH --job-name={job_name}
#SBATCH --output={job_name}.out
#SBATCH --error={job_name}.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={total_cores}
#SBATCH --time=48:00:00
{nodelist_line}

# 1. Environment & Stack Limit
ulimit -s unlimited

# 2. Load Tchenla Modules & Custom Libraries
if [ -f "{base_dir}/load_modules_tchenla.sh" ]; then
    source "{base_dir}/load_modules_tchenla.sh"
fi

export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

echo "=========================================================="
echo "Job:            {job_name} (SLURM ID: $SLURM_JOB_ID)"
echo "Node:           $(hostname)"
echo "Date:           $(date)"
echo "Available CPUs: $(taskset -cp $$ 2>/dev/null || echo 'N/A')"
echo "=========================================================="
echo ""
echo "Configuration:"
echo "  Parallel jobs:   {n_jobs}"
echo "  Threads per job: {threads_per_job}"
echo "  Total cores:     {total_cores}"
echo "  System size:     {nx}x{ny}"
echo "  Mode setting:    {mode_choice}"
echo "  Start seed:      {start_seed}"
echo "  Remeshing:       {remesh_status}"
echo "  Preconditioner:  {precond_status}"
echo "  Stability mon.:  {eig_status}"
echo "=========================================================="
echo ""

# Launch parallel jobs
"{parallel_file}"

echo ""
echo "=========================================================="
echo "SLURM job completed at $(date)"
echo "=========================================================="
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
    print("To submit this job to the Tchenla SLURM queue:")
    print(f"  cd {base_dir}")
    print(f"  sbatch run_{job_name}.sh")
    print()
    print("To monitor progress in real time:")
    print(f"  squeue -u $USER")
    print(f"  tail -f {job_name}.out")
    print(f"  tail -f {runs_dir}/*/simulation.log")
    print()


if __name__ == "__main__":
    main()
