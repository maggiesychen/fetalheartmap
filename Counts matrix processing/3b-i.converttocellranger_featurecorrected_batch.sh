#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

base_output="/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/3b-i.cellranger_format_featurecorrected"
script="${SCRIPT_DIR}/3b-i.converttocellranger_featurecorrected.py"

# Step 1: run sublib1 first to build master gene list
MASTER_JOB=$(sbatch --partition=engreitz --mem=64G --time=2:00:00 \
    --job-name=master_gene_list \
    --output=${base_output}/logs/master_gene_list_%j.out \
    --error=${base_output}/logs/master_gene_list_%j.err \
    --wrap="source /home/groups/engreitz/Software/anaconda3/etc/profile.d/conda.sh && conda activate perturbprocessing && python ${script} --sublib 1" \
    | awk '{print $NF}')

echo "Sublib1/master gene list job: $MASTER_JOB"

# Step 2: submit sublibs 2-7 with dependency on sublib1 finishing
for SUBLIB in 2 3 4 5 6 7; do
    sbatch --dependency=afterok:${MASTER_JOB} \
           --partition=engreitz \
           --mem=64G \
           --time=12:00:00 \
           --job-name=cellranger_convert_sublib${SUBLIB} \
           --output=${base_output}/logs/convert_sublib${SUBLIB}_%j.out \
           --error=${base_output}/logs/convert_sublib${SUBLIB}_%j.err \
           --wrap="source /home/groups/engreitz/Software/anaconda3/etc/profile.d/conda.sh && conda activate perturbprocessing && python ${script} --sublib ${SUBLIB}"
done