#!/bin/bash
#SBATCH --job-name=003_data-phenotype
#SBATCH --output=src/logs/003_data-phenotype.log
#SBATCH --error=src/logs/003_data-phenotype.err
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=truncate

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(dirname -- "${BASH_SOURCE[0]:-$0}")/../}"

Rscript src/003_data-phenotype.R
