#!/bin/bash
# ==============================================================================
# SLURM Submission Script for Single Simulation on Tchenla Cluster
# Paris 13 University (Université Sorbonne Paris Nord / LSPM)
# Hardware:
#   c[1-2]: 2x Intel Xeon-Gold 6430 (64 cores, 128 threads, 512GB RAM)
#   c[3-4]: 2x Intel Xeon-Gold 6430 (64 cores, 128 threads, 256GB RAM)
# ==============================================================================
#SBATCH --job-name=shear_sim
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00

# 1. Environment & Stack Limit (Recommended by Tchenla wiki)
ulimit -s unlimited

# 2. Directory Discovery & Module Loading
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR" || exit 1

if [ -f "$SCRIPT_DIR/load_modules_tchenla.sh" ]; then
    source "$SCRIPT_DIR/load_modules_tchenla.sh"
else
    echo "Warning: load_modules_tchenla.sh not found in $SCRIPT_DIR"
fi

# 3. Thread and Thread-Binding Configuration
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-32}
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export OMP_PROC_BIND=close
export OMP_PLACES=cores

# 4. Simulation Parameters (Can be passed as CLI arguments or use defaults)
# Usage:
#   sbatch submit_single_tchenla.sh [NX] [NY] [MODE] [SEED] [REMESH]
# Examples:
#   sbatch submit_single_tchenla.sh 150 150 positive 42 1
#   sbatch submit_single_tchenla.sh 150 150 negative 43 0
NX=${1:-150}
NY=${2:-150}
MODE=${3:-"positive"}       # "positive" or "negative"
SEED=${4:-42}
REMESH=${5:-1}              # 1 = remeshing enabled, 0 = no remeshing

REMESH_TAG=$([ "$REMESH" = "1" ] && echo "remesh" || echo "noremesh")
RUN_DIR="${SCRIPT_DIR}/run_${MODE}_${NX}x${NY}_seed${SEED}_${REMESH_TAG}"

mkdir -p "$RUN_DIR"
cd "$RUN_DIR" || exit 1

echo "=================================================================="
echo "  Tchenla SLURM Simulation Job"
echo "=================================================================="
echo "  Job ID:            ${SLURM_JOB_ID:-N/A}"
echo "  Job Name:          ${SLURM_JOB_NAME:-shear_sim}"
echo "  Node:              $(hostname)"
echo "  Date:              $(date)"
echo "  Working Directory: $(pwd)"
echo "  Allocated CPUs:    $(taskset -cp $$ 2>/dev/null || echo 'N/A')"
echo "  OpenMP Threads:    $OMP_NUM_THREADS"
echo "  Parameters:        size=${NX}x${NY}, mode=${MODE}, seed=${SEED}, remesh=${REMESH} (${REMESH_TAG})"
echo "=================================================================="
echo ""

EXE="${SCRIPT_DIR}/lattice_triangulation"
if [ ! -f "$EXE" ]; then
    echo "ERROR: Executable '$EXE' does not exist!"
    echo "Please compile it first: 'make -f Makefile.tchenla -j 8'"
    exit 1
fi

# Execute simulation and stream to both console (slurm-%j.out) and simulation.log
"$EXE" "$NX" "$NY" "$MODE" "$SEED" "$REMESH" 2>&1 | tee simulation.log
EXIT_CODE=${PIPESTATUS[0]}

echo ""
echo "=================================================================="
echo "Simulation completed at $(date) with exit code: $EXIT_CODE"
echo "Output files saved in: $(pwd)"
echo "=================================================================="

exit $EXIT_CODE
