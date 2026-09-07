#!/bin/bash
#submit with bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

declare -A UMI_MIN
declare -A UMI_MAX

# hardcoded params for UMI filtering (min and max, determined based off of TZ pipeline BRI)
UMI_MIN[1]=1921;  UMI_MAX[1]=75000
UMI_MIN[2]=2414;  UMI_MAX[2]=75000
UMI_MIN[3]=2006;  UMI_MAX[3]=75000
UMI_MIN[4]=2345;  UMI_MAX[4]=75000
UMI_MIN[5]=3132;  UMI_MAX[5]=75000
UMI_MIN[6]=2925;  UMI_MAX[6]=75000
UMI_MIN[7]=1719;  UMI_MAX[7]=75000

for SUBLIB in 1 2 3 4 5 6 7; do
    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=umi_filter_sublib${SUBLIB}
#SBATCH --output=/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/1.umifiltering/logs/umi_filter_sublib${SUBLIB}_%j.out
#SBATCH --error=/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/1.umifiltering/logs/umi_filter_sublib${SUBLIB}_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=4
#SBATCH --partition=bigmem

source /home/groups/engreitz/Software/anaconda3/etc/profile.d/conda.sh
conda activate perturbprocessing
python "${SCRIPT_DIR}/1.umi_threshold_float.py" --sublib ${SUBLIB} --umi_threshold ${UMI_MIN[$SUBLIB]} --umi_max ${UMI_MAX[$SUBLIB]}
EOF
done