#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for SUBLIB in 1 2 3 4 5 6 7; do
    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=cellranger_convert${SUBLIB}
#SBATCH --output=/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/3b.cellranger_format/logs/convert_sublib${SUBLIB}_%j.out
#SBATCH --error=/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/3b.cellranger_format/logs/convert_sublib${SUBLIB}_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --partition=engreitz

source /home/groups/engreitz/Software/anaconda3/etc/profile.d/conda.sh
conda activate perturbprocessing
python "${SCRIPT_DIR}/3b.converttocellranger.py" --sublib ${SUBLIB}
EOF
done