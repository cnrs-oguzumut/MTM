#!/bin/bash
# ==============================================================================
# SLURM Submission Script for Tchenla (LSPM CNRS / Université Sorbonne Paris Nord)
# Reference: https://wiki.lspm.cnrs.fr/doku.php/en/tchenla#submission_scripts
#
# Hardware on Tchenla:
#   c[1-2] : 2x Intel Xeon-Gold 6430 (64 physical cores, 128 threads, 512GB RAM)
#   c[3-4] : 2x Intel Xeon-Gold 6430 (64 physical cores, 128 threads, 256GB RAM)
#   c[5-7] : 2x Intel Xeon E5-2630 v2 (12 cores, 24 threads, 125GB RAM)
# ==============================================================================

#SBATCH --job-name="shear_sim"
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00

# Optional: Target specific compute node (uncomment if desired)
# #SBATCH --nodelist=c1

# ------------------------------------------------------------------------------
# 1. Stack Size Configuration (Recommended by Tchenla Wiki)
# ------------------------------------------------------------------------------
ulimit -s unlimited

# ------------------------------------------------------------------------------
# 2. Directory and Environment Discovery
# ------------------------------------------------------------------------------
# Use SLURM_SUBMIT_DIR (where sbatch was invoked) instead of spool directory
PROJECT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$PROJECT_DIR" || exit 1

if [ -f "${PROJECT_DIR}/load_modules_tchenla.sh" ]; then
    source "${PROJECT_DIR}/load_modules_tchenla.sh"
else
    echo "Warning: load_modules_tchenla.sh not found in ${PROJECT_DIR}"
fi

# Ensure Intel MKL runtime is in library path
if [ -d "/opt/intel/oneapi/mkl/2025.0/lib" ]; then
    export LD_LIBRARY_PATH="/opt/intel/oneapi/mkl/2025.0/lib:${LD_LIBRARY_PATH}"
fi

# ------------------------------------------------------------------------------
# 3. Threading Configuration (Wiki page 8: OpenMP Script)
# ------------------------------------------------------------------------------
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-32}
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export OMP_PROC_BIND=close
export OMP_PLACES=cores

# ------------------------------------------------------------------------------
# 4. Simulation Parameters
#    Defaults or passed as arguments: sbatch run_tchenla.sbatch [NX] [NY] [MODE] [SEED] [REMESH]
# ------------------------------------------------------------------------------
NX=${1:-150}
NY=${2:-150}
MODE=${3:-"positive"}       # "positive" or "negative"
SEED=${4:-42}
REMESH=${5:-1}              # 1 = enabled, 0 = disabled

REMESH_TAG=$([ "$REMESH" = "1" ] && echo "remesh" || echo "noremesh")
RUN_DIR="${PROJECT_DIR}/run_${MODE}_${NX}x${NY}_seed${SEED}_${REMESH_TAG}"

mkdir -p "$RUN_DIR"
cd "$RUN_DIR" || exit 1

echo "=================================================================="
echo "  Tchenla SLURM Simulation Job"
echo "=================================================================="
echo "  Job ID:            ${SLURM_JOB_ID:-N/A}"
echo "  Job Name:          ${SLURM_JOB_NAME:-shear_sim}"
echo "  Compute Node:      $(hostname)"
echo "  Start Time:        $(date)"
echo "  Working Directory: $(pwd)"
echo "  Allocated CPUs:    ${SLURM_CPUS_PER_TASK:-32}"
echo "  CPU Affinity:      $(taskset -cp $$ 2>/dev/null || echo 'N/A')"
echo "  OpenMP Threads:    $OMP_NUM_THREADS"
echo "  Parameters:        size=${NX}x${NY}, mode=${MODE}, seed=${SEED}, remesh=${REMESH} (${REMESH_TAG})"
echo "=================================================================="
echo ""

# ------------------------------------------------------------------------------
# 5. Launch Simulation
# ------------------------------------------------------------------------------
EXE="${PROJECT_DIR}/lattice_triangulation"
if [ ! -f "$EXE" ]; then
    echo "ERROR: Executable '$EXE' not found!"
    echo "Please compile first: make -f Makefile.tchenla -j 8"
    exit 1
fi

# Run simulation and mirror output to both slurm-%j.out and simulation.log
"$EXE" "$NX" "$NY" "$MODE" "$SEED" "$REMESH" 2>&1 | tee simulation.log
EXIT_CODE=${PIPESTATUS[0]}

echo ""
echo "=================================================================="
echo "Job finished at $(date) with exit code: $EXIT_CODE"
echo "Results saved in: $(pwd)"
echo "=================================================================="

exit $EXIT_CODE
