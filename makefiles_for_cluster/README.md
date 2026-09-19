# Cluster Job Setup Guide (Magi & Tchenla)

This directory contains automated parallel job and SLURM submission script generators for the **Magi** and **Tchenla** high-performance computing clusters.

---

## 1. Quick Start Commands

### Magi Cluster

```bash
# Generate 2 parallel jobs (positive and negative shear to alpha = 3.0) on 32 threads
python3 makefiles_for_cluster/job_creator.py \
    --nx 300 --ny 300 \
    --alpha-end 3.0 \
    --checkpoint-interval 500 \
    --stress-drop-threshold 0.10 \
    --n-jobs 2 \
    --threads 32 \
    --mode both \
    -y

# Submit the generated SLURM job
sbatch run_shear_300x300_alpha3.sh
```

### Tchenla Cluster

```bash
python3 makefiles_for_cluster/job_creator_tchenla.py \
    --nx 300 --ny 300 \
    --alpha-end 3.0 \
    --checkpoint-interval 500 \
    --stress-drop-threshold 0.10 \
    --n-jobs 2 \
    --threads 32 \
    --mode both \
    -y

# Submit the generated SLURM job
sbatch run_shear_300x300_alpha3.sh
```

---

## 2. Command-Line Arguments Reference

| Flag | Type | Default | Description |
| :--- | :---: | :---: | :--- |
| `--nx`, `--ny` | Integer | Prompt / `150` | System dimensions in lattice cells (e.g. `300 300`). |
| `--mode` | Choice | `both` | Loading mode: `both` (positive & negative in parallel), `positive`, or `negative`. |
| `--alpha-end` / `--alpha-max` | Float | `1.0` | Target shear strain magnitude (e.g. `3.0` runs up to $\alpha = \pm 3.0$). |
| `--alpha-start` / `--alpha-min`| Float | `0.14` | Starting shear strain magnitude before incremental loading. |
| `--step-size` | Float | `6e-5` | Strain increment per step $|d\alpha|$. |
| `--checkpoint-interval` | Integer | `500` | Periodic elastic checkpoint interval in steps (`0` = disabled, saving only at avalanches). |
| `--stress-drop-threshold` | Float | `0.10` | Fractional stress drop $\Delta \sigma / \|\sigma\|$ threshold to trigger avalanche save (**larger saves less**, e.g. `0.20`). |
| `--save-triangle-data` | Flag | disabled | Enable legacy ASCII `triangle_data/` dump (disabled by default to save 15+ GB). |
| `--precond` | Choice | `stiffness` | L-BFGS preconditioner: `stiffness`, `laplacian`, `diag`, `none`. |
| `--precond-from-step` | Integer | `1` | Step from which preconditioning is activated (`1` keeps plain step 0 relaxation). |
| `--n-jobs` | Integer | `2` | Number of concurrent simulations running in parallel on the node. |
| `--threads` | Integer | `32` | OpenMP CPU threads allocated per job. |
| `-y`, `--yes` | Flag | false | Auto-confirm with defaults and skip interactive prompts. |

---

## 3. Useful Recipes

### Save Minimum Data (Maximum Speed & Quota Efficiency)
To save only major avalanches ($\ge 20\%$ drop) and disable periodic elastic checkpoints entirely:
```bash
python3 makefiles_for_cluster/job_creator.py \
    --nx 300 --ny 300 \
    --alpha-end 3.0 \
    --checkpoint-interval 0 \
    --stress-drop-threshold 0.20 \
    --n-jobs 2 \
    --threads 32 \
    --mode both \
    -y
```

### Resume an Interrupted Job
Every job script automatically detects `checkpoints/latest.chk` if the simulation was interrupted. Simply resubmit:
```bash
sbatch run_shear_300x300_alpha3.sh
```

---

## 4. Monitoring Jobs on Cluster

```bash
# Check queue status
squeue -u $USER

# Follow live simulation logs
tail -f runs_300x300_alpha3/run_positive_seed42/simulation.log
tail -f runs_300x300_alpha3/run_negative_seed43/simulation.log
```
