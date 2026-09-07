#!/bin/bash
#submit with bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for SUBLIB in 1 2 3 4 5 6 7; do
    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=qc_filter_sublib${SUBLIB}
#SBATCH --output=/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/2.qcfiltering/float_logs/qc_filter_sublib${SUBLIB}_%j.out
#SBATCH --error=/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/2.qcfiltering/float_logs/qc_filter_sublib${SUBLIB}_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --partition=engreitz

source /home/groups/engreitz/Software/anaconda3/etc/profile.d/conda.sh
conda activate perturbprocessing

python "${SCRIPT_DIR}/2.qc_metrics_filter.py" --sublib ${SUBLIB}
EOF
done