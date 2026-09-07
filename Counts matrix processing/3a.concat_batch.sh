#!/bin/bash
# submit from inside this folder:  sbatch 3a.concat_batch.sh
#SBATCH --job-name=concat
#SBATCH --output=/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/3a.concat-old-forclustering/logs/concat_%j.out
#SBATCH --error=/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/3a.concat-old-forclustering/logs/concat_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --partition=engreitz

# Slurm copies the batch script to a spool dir, so resolve via the submit dir.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

source /home/groups/engreitz/Software/anaconda3/etc/profile.d/conda.sh
conda activate perturbprocessing

python "${SCRIPT_DIR}/3a.concatenate_mudata.py"
