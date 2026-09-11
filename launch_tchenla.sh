#!/bin/bash
# ==============================================================================
# Interactive SLURM Job Launcher for Tchenla Cluster
# Paris 13 University (Université Sorbonne Paris Nord / LSPM)
# ==============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SBATCH_SCRIPT="${SCRIPT_DIR}/run_tchenla.sbatch"

if [ ! -f "$SBATCH_SCRIPT" ]; then
    echo "ERROR: $SBATCH_SCRIPT not found!"
    exit 1
fi

echo "=================================================================="
echo "  Tchenla SLURM Job Launcher (Lattice Triangulation)"
echo "=================================================================="
echo ""

# 1. Prompt for parameters with sensible defaults
read -p "System size NX [150]: " NX
NX=${NX:-150}

read -p "System size NY [150]: " NY
NY=${NY:-150}

echo "Loading mode options: 'positive', 'negative', 'both'"
read -p "Loading mode [positive]: " MODE
MODE=${MODE:-positive}

read -p "Random seed [42]: " SEED
SEED=${SEED:-42}

read -p "Enable adaptive remeshing? (yes/no) [yes]: " REMESH_IN
REMESH_IN=${REMESH_IN:-yes}
if [[ "$REMESH_IN" =~ ^(y|yes|1|true)$ ]]; then
    REMESH=1
    REMESH_STR="remesh"
else
    REMESH=0
    REMESH_STR="noremesh"
fi

read -p "Number of CPU cores (--cpus-per-task) [32]: " CORES
CORES=${CORES:-32}

read -p "Target compute nodes (--nodelist) [c[1-4]]: " NODELIST
NODELIST=${NODELIST:-c[1-4]}

echo ""
echo "------------------------------------------------------------------"
echo "Job Submission Summary:"
echo "  Nodes:     $NODELIST (Intel Xeon-Gold 6430)"
echo "  CPUs:      $CORES cores"
echo "  Size:      ${NX}x${NY}"
echo "  Mode:      $MODE"
echo "  Seed:      $SEED"
echo "  Remeshing: $REMESH_STR"
echo "------------------------------------------------------------------"
echo ""

# 2. Submit to SLURM
if [ "$MODE" = "both" ]; then
    echo "Submitting paired jobs (positive & negative)..."
    
    JOB_NAME_POS="shear_${NX}x${NY}_pos_s${SEED}_${REMESH_STR}"
    echo "Submitting Positive run..."
    sbatch --nodelist="$NODELIST" --cpus-per-task="$CORES" --job-name="$JOB_NAME_POS" \
        "$SBATCH_SCRIPT" "$NX" "$NY" "positive" "$SEED" "$REMESH"
        
    JOB_NAME_NEG="shear_${NX}x${NY}_neg_s${SEED}_${REMESH_STR}"
    echo "Submitting Negative run..."
    sbatch --nodelist="$NODELIST" --cpus-per-task="$CORES" --job-name="$JOB_NAME_NEG" \
        "$SBATCH_SCRIPT" "$NX" "$NY" "negative" "$SEED" "$REMESH"
else
    JOB_NAME="shear_${NX}x${NY}_${MODE}_s${SEED}_${REMESH_STR}"
    echo "Submitting SLURM job..."
    sbatch --nodelist="$NODELIST" --cpus-per-task="$CORES" --job-name="$JOB_NAME" \
        "$SBATCH_SCRIPT" "$NX" "$NY" "$MODE" "$SEED" "$REMESH"
fi

echo ""
echo "=================================================================="
echo "Job(s) submitted to SLURM queue!"
echo "Check status:  squeue -u \$USER"
echo "Check nodes:   sinfo"
echo "Live output:   tail -f slurm-*.out"
echo "=================================================================="
