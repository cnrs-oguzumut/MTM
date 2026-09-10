#!/usr/bin/env python3
"""
Automatic Job Setup Generator for Parallel Crystal Plasticity Simulations
Creates SLURM script, parallel launcher, and individual run scripts
adapted for the new lattice_triangulation interface:
    lattice_triangulation <nx> <ny> <mode> <seed>
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
    print("  Parallel Job Setup Generator (Lattice Triangulation)")
    print("=" * 70)
    print()

    # 1. Resource Configuration
    print("--- 1. Hardware & Core Allocation ---")
    n_jobs = int(get_input("Number of parallel jobs to run simultaneously", "4"))
    threads_per_job = int(get_input("OpenMP threads per job", "32"))
    total_cores = n_jobs * threads_per_job

    print(f"  -> Total cores allocated: {total_cores} ({n_jobs} jobs x {threads_per_job} threads)")
    print()

    # 2. Physics & Simulation Parameters
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
    print()

    # 3. Paths and Directories
    print("--- 3. Cluster Paths ---")
    base_dir = get_input("Base directory for MTM", "/home/dist/umut.salman/new/MTM")
    exe_path = get_input("Executable path", f"{base_dir}/lattice_triangulation")
    
    default_runs_dir = f"{base_dir}/runs_{nx}x{ny}"
    runs_dir = get_input("Runs output directory", default_runs_dir)

    # 4. SLURM Options
    print()
    print("--- 4. SLURM Cluster Options ---")
    job_name = get_input("SLURM job name", f"shear_{nx}x{ny}")
    partition = get_input("SLURM partition", "COMPUTE2")
    email = get_input("Email for notifications", "umut.salman@lspm.cnrs.fr")

    # Build per-job specifications: [(job_id, nx, ny, mode, seed, run_folder)]
    job_specs = []
    for i in range(n_jobs):
        job_id = i + 1
        if mode_choice == "both":
            seed = start_seed + (i // 2)
            mode = "positive" if (i % 2 == 0) else "negative"
        elif mode_choice == "positive":
            seed = start_seed + i
            mode = "positive"
        else:  # negative
            seed = start_seed + i
            mode = "negative"
        
        folder_name = f"run_{job_id:02d}_{mode}_seed{seed}"
        job_specs.append((job_id, nx, ny, mode, seed, folder_name))

    print()
    print("=" * 70)
    print("Planned Job Configurations:")
    print("=" * 70)
    for j_id, j_nx, j_ny, j_mode, j_seed, j_dir in job_specs:
        start_cpu = (j_id - 1) * threads_per_job
        end_cpu = start_cpu + threads_per_job - 1
        print(f"  Job {j_id:2d}: CPUs {start_cpu:3d}-{end_cpu:3d} | {j_nx}x{j_ny} | mode={j_mode:8s} | seed={j_seed:<4d} | dir={j_dir}")
    print("=" * 70)
    print()

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
# Usage: ./run_one.sh <RUN_ID> <NX> <NY> <MODE> <SEED> <RUN_DIR>

RUN_ID=$1
NX=$2
NY=$3
MODE=$4
SEED=$5
RUN_DIR=$6
EXE="{exe_path}"

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
echo "  Parameters:        size=${{NX}}x${{NY}}, mode=${{MODE}}, seed=${{SEED}}"
echo "  CPU affinity:      $(taskset -cp $$ 2>/dev/null || echo 'N/A')"
echo "  OpenMP threads:    $OMP_NUM_THREADS"
echo "  Executable:        $EXE"
echo "=========================================================="

# Execute simulation
"$EXE" "$NX" "$NY" "$MODE" "$SEED" > simulation.log 2>&1
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
# Parallel launcher script: pins each job to dedicated CPU cores
SCRIPT_DIR="{base_dir}"
RUNS_BASE="{runs_dir}"

mkdir -p "$RUNS_BASE"

echo "=========================================================="
echo "Starting {n_jobs} parallel jobs with CPU pinning..."
echo "Runs directory: $RUNS_BASE"
echo "=========================================================="
echo ""

"""
    for j_id, j_nx, j_ny, j_mode, j_seed, j_dir in job_specs:
        start_cpu = (j_id - 1) * threads_per_job
        end_cpu = start_cpu + threads_per_job - 1
        parallel_script += (
            f"echo \"Launching Job {j_id}: mode={j_mode}, seed={j_seed} on CPUs {start_cpu}-{end_cpu}...\"\n"
            f"taskset -c {start_cpu}-{end_cpu} \"${{SCRIPT_DIR}}/run_{job_name}_one.sh\" "
            f"{j_id} {j_nx} {j_ny} {j_mode} {j_seed} \"${{RUNS_BASE}}/{j_dir}\" &\n\n"
        )

    parallel_script += """
# Wait for all parallel background jobs to finish
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
