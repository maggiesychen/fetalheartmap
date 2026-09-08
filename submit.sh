#!/usr/bin/env bash
# Submit a workflow stage to Slurm, one job per rule.
#
# Usage:
#   CONFIG=config/config.20260408_ipscvic_300k.yaml ./submit.sh [snakemake args]
#
# Examples:
#   CONFIG=config/config.20260408_ipscvic_300k.yaml ./submit.sh --dry-run
#   CONFIG=config/config.20260408_ipscvic_300k.yaml ./submit.sh --jobs 20
#   CONFIG=... WORKFLOW=workflows/00_counts_matrix_processing.smk ./submit.sh
#
# Per-rule memory, walltime, thread count and partition come from the
# `resources:` block of the config file. This uses the Snakemake 7 `--cluster`
# interface; see envs/snakemake.yml for the pinned version.

set -euo pipefail

CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW="${WORKFLOW:-$CODE_DIR/workflows/00_counts_matrix_processing.smk}"

if [[ -z "${CONFIG:-}" ]]; then
    echo "Error: CONFIG must be set." >&2
    echo "Usage: CONFIG=config/config.<analysis>.yaml ./submit.sh [args]" >&2
    exit 1
fi

if [[ ! -f "$CONFIG" ]]; then
    echo "Error: config file not found: $CONFIG" >&2
    exit 1
fi

if ! command -v snakemake >/dev/null 2>&1; then
    echo "Error: snakemake is not on PATH." >&2
    echo "Create it once with:  mamba env create -f envs/snakemake.yml" >&2
    echo "Then:                 conda activate fetalheartmap-snakemake" >&2
    exit 1
fi

# Cluster stdout/stderr go next to the workflow logs, under results_base.
CLUSTER_LOGS="$(
    python - "$CONFIG" <<'PY'
import os
import sys

import yaml

with open(sys.argv[1]) as handle:
    config = yaml.safe_load(handle)
print(os.path.join(config["results_base"], "logs", "cluster"))
PY
)"
mkdir -p "$CLUSTER_LOGS"

CLUSTER_CMD=(
    sbatch
    --parsable
    --partition="{resources.slurm_partition}"
    --mem="{resources.mem_mb}"
    --time="{resources.runtime}"
    --cpus-per-task="{threads}"
    --job-name="smk-{rule}"
    --output="$CLUSTER_LOGS/{rule}-%j.out"
    --error="$CLUSTER_LOGS/{rule}-%j.err"
)

ARGS=(
    --snakefile "$WORKFLOW"
    --configfile "$CONFIG"
    --cluster "${CLUSTER_CMD[*]}"
    --latency-wait 120
    --keep-going
    --rerun-incomplete
    --printshellcmds
    --rerun-triggers mtime
)

# Snakemake needs a --jobs cap when submitting to a cluster.
if [[ ! " $* " =~ " --jobs " ]] && [[ ! " $* " =~ " -j " ]]; then
    ARGS+=(--jobs 20)
fi

echo "Workflow    : $WORKFLOW"
echo "Config      : $CONFIG"
echo "Cluster logs: $CLUSTER_LOGS"
echo "EXECUTING: snakemake ${ARGS[*]} $*"
echo

snakemake "${ARGS[@]}" "$@"
