# Stage 02 -- per-gRNA knockdown of the target gene panel.
#
# SCEPTRE reports one fold change per target, pooling that target's gRNAs. This
# stage recomputes knockdown one gRNA at a time, by two independent estimators:
#
#   07a  Poisson GLM log2 fold change against sceptre's own null model
#   07b  model-free CP10K pseudobulk against a per-gene non-targeting pool
#
# They share the count data and the gRNA assignment but nothing else, so their
# agreement (reported by 07d) is evidence the knockdown is real rather than a
# modelling artefact.
#
# `per_guide_knockdown.sceptre_variant` names the stage-01 variant supplying the
# gRNA assignment and the pair sets; it defaults to trans_v2. gRNA assignment is
# identical across variants, so that choice only affects the pair sets.
#
# Run with:
#   CONFIG=config/config.20260408_ipscvic_300k.yaml \
#   WORKFLOW=workflows/02_per_guide_knockdown.smk ./submit.sh
#
# The DAG is wider than the script numbering suggests: 07a (targeting), 07a
# (non-targeting) and 07b are mutually independent and all three read only
# stage-01 outputs.

import sys
from pathlib import Path

CODE_DIR = Path(workflow.snakefile).resolve().parent.parent
sys.path.insert(0, str(CODE_DIR))

from scripts.pipeline_utils import resolve_config_paths
from scripts.snakemake_helpers import (
    get_base_shell_prefix,
    get_logs_path,
    get_python_env_prefix,
    get_r_env_prefix,
    get_resources,
    get_results_path,
    print_path_configuration,
)

resolve_config_paths(config, CODE_DIR)

RESULTS = get_results_path(config=config)
LOGS = get_logs_path(config=config)

SCEPTRE = config["sceptre"]
KD = config["per_guide_knockdown"]
VARIANT = KD["sceptre_variant"]

if VARIANT not in SCEPTRE["variants"]:
    raise ValueError(
        f"per_guide_knockdown.sceptre_variant = {VARIANT!r} is not a declared "
        f"sceptre variant: {sorted(SCEPTRE['variants'])}"
    )

ONDISC = f"{RESULTS}/{SCEPTRE['ondisc_dir']}"
SCEPTRE_OUT = f"{RESULTS}/sceptre/{VARIANT}/outputs"
OUT = f"{RESULTS}/per_guide_knockdown"

print_path_configuration(config)
print(f"knockdown source : sceptre variant {VARIANT}")

shell.prefix(get_base_shell_prefix(CODE_DIR))

PY_ENV = get_python_env_prefix(config)
R_ENV = get_r_env_prefix(config)


def resources_for(rule_name):
    return {
        key: get_resources(config, rule_name, key)
        for key in ("mem_mb", "runtime", "slurm_partition")
    }


def threads_for(rule_name):
    return get_resources(config, rule_name, "threads")


# Inputs shared by all three estimator rules.
SCEPTRE_INPUTS = dict(
    sceptre_object=f"{ONDISC}/sceptre_object_{VARIANT}.rds",
    gene_odm=f"{ONDISC}/gene.odm",
    grna_odm=f"{ONDISC}/grna.odm",
    assignment=f"{SCEPTRE_OUT}/grna_assignment_matrix.rds",
    guide_targets=SCEPTRE["guide_targets_file"],
    ens2symbol=SCEPTRE["ens2symbol_file"],
)

COMMON_ARGS = """\
            --sceptre-object {input.sceptre_object} \
            --response-odm {input.gene_odm} \
            --grna-odm {input.grna_odm} \
            --assignment-matrix {input.assignment} \
            --guide-targets {input.guide_targets} \
            --ens2symbol {input.ens2symbol} \
"""


rule all:
    input:
        f"{OUT}/fig_pct_knockdown_bar.svg",
        f"{OUT}/source_data_pct_knockdown_bar.csv",
        f"{OUT}/fig_tss_knockdown.svg",
        f"{OUT}/fig_nt_vs_tss_expression.svg",


# =============================================================================
# Step 7a -- Poisson per-gRNA fold changes (two guide sets, one script)
# =============================================================================
rule per_guide_fold_changes_targeting:
    """Per-gRNA log2 fold change for every targeting (target, response) pair."""
    input:
        **SCEPTRE_INPUTS,
        script=CODE_DIR / "scripts/07a_per_guide_fold_changes_poisson.R",
    output:
        tsv=f"{OUT}/per_guide_fold_changes_poisson.targeting.tsv",
    params:
        env=R_ENV,
        formula=KD["poisson"]["covariate_formula"],
        min_cells=KD["poisson"]["min_cells_per_guide"],
        nt_label=KD["nt_target_label"],
    log:
        f"{LOGS}/per_guide_knockdown/fold_changes_targeting.log",
    threads: threads_for("knockdown_poisson")
    resources:
        **resources_for("knockdown_poisson"),
    shell:
        """
        mkdir -p $(dirname {log})
        {params.env}
        Rscript {input.script} \
"""
        + COMMON_ARGS
        + """\
            --guide-set targeting \
            --nt-target-label "{params.nt_label}" \
            --covariate-formula "{params.formula}" \
            --min-cells-per-guide {params.min_cells} \
            --output {output.tsv} \
            &> {log}
        """


rule per_guide_fold_changes_non_targeting:
    """The same estimator on the non-targeting gRNAs -- the empirical null."""
    input:
        **SCEPTRE_INPUTS,
        gene_list=KD["chd_gene_table"],
        script=CODE_DIR / "scripts/07a_per_guide_fold_changes_poisson.R",
    output:
        tsv=f"{OUT}/per_guide_fold_changes_poisson.non_targeting.tsv",
    params:
        env=R_ENV,
        formula=KD["poisson"]["covariate_formula"],
        min_cells=KD["poisson"]["min_cells_per_guide"],
        nt_label=KD["nt_target_label"],
        gene_column=KD["chd_gene_column"],
    log:
        f"{LOGS}/per_guide_knockdown/fold_changes_non_targeting.log",
    threads: threads_for("knockdown_poisson")
    resources:
        **resources_for("knockdown_poisson"),
    shell:
        """
        mkdir -p $(dirname {log})
        {params.env}
        Rscript {input.script} \
"""
        + COMMON_ARGS
        + """\
            --guide-set non-targeting \
            --gene-list {input.gene_list} \
            --gene-list-column "{params.gene_column}" \
            --nt-target-label "{params.nt_label}" \
            --covariate-formula "{params.formula}" \
            --min-cells-per-guide {params.min_cells} \
            --output {output.tsv} \
            &> {log}
        """


# =============================================================================
# Step 7b -- model-free CP10K pseudobulk
# =============================================================================
rule target_gene_expression_per_guide:
    """Per-gRNA CP10K of each panel gene, plus a per-gene non-targeting pool."""
    input:
        **SCEPTRE_INPUTS,
        gene_list=KD["chd_gene_table"],
        script=CODE_DIR / "scripts/07b_target_gene_expression_per_guide.R",
    output:
        per_guide=f"{OUT}/target_gene_expression_per_guide.tsv",
        nt_pool=f"{OUT}/target_gene_expression_nt_pool.tsv",
    params:
        env=R_ENV,
        gene_column=KD["chd_gene_column"],
        nt_label=KD["nt_target_label"],
        tss_regex=KD["tss_target_regex"],
        cp10k_scale=KD["expression"]["cp10k_scale"],
        min_cells=KD["expression"]["min_cells_per_guide"],
    log:
        f"{LOGS}/per_guide_knockdown/target_gene_expression.log",
    threads: threads_for("knockdown_expression")
    resources:
        **resources_for("knockdown_expression"),
    shell:
        """
        mkdir -p $(dirname {log})
        {params.env}
        Rscript {input.script} \
"""
        + COMMON_ARGS
        + """\
            --gene-list {input.gene_list} \
            --gene-list-column "{params.gene_column}" \
            --nt-target-label "{params.nt_label}" \
            --tss-target-regex "{params.tss_regex}" \
            --cp10k-scale {params.cp10k_scale} \
            --min-cells-per-guide {params.min_cells} \
            --output-per-guide {output.per_guide} \
            --output-nt-pool {output.nt_pool} \
            &> {log}
        """


# =============================================================================
# Step 7c/7d -- tests and per-gene summaries
# =============================================================================
rule plot_tss_knockdown:
    """Mann-Whitney on log2 fold change, per gene; exports the per-gRNA dots."""
    input:
        targeting=f"{OUT}/per_guide_fold_changes_poisson.targeting.tsv",
        non_targeting=f"{OUT}/per_guide_fold_changes_poisson.non_targeting.tsv",
        gene_list=KD["chd_gene_table"],
        script=CODE_DIR / "scripts/07c_plot_tss_knockdown.py",
    output:
        per_guide=f"{OUT}/source_data_per_guide_tss_knockdown.csv",
        summary=f"{OUT}/source_data_tss_knockdown_summary.csv",
        figure=f"{OUT}/fig_tss_knockdown.svg",
    params:
        env=PY_ENV,
        config_file=lambda w: workflow.configfiles[0],
        stem=f"{OUT}/fig_tss_knockdown",
    log:
        f"{LOGS}/per_guide_knockdown/plot_tss_knockdown.log",
    threads: threads_for("knockdown_plot")
    resources:
        **resources_for("knockdown_plot"),
    shell:
        """
        mkdir -p $(dirname {log})
        {params.env}
        python -u {input.script} \
            --targeting {input.targeting} \
            --non-targeting {input.non_targeting} \
            --config {params.config_file} \
            --output-per-guide {output.per_guide} \
            --output-summary {output.summary} \
            --output-figure {params.stem} \
            &> {log}
        """


rule plot_nt_vs_tss:
    """Mann-Whitney on CP10K, per gene; exports the summary 07e orders by."""
    input:
        expression=f"{OUT}/target_gene_expression_per_guide.tsv",
        targeting=f"{OUT}/per_guide_fold_changes_poisson.targeting.tsv",
        gene_list=KD["chd_gene_table"],
        script=CODE_DIR / "scripts/07d_plot_nt_vs_tss.py",
    output:
        summary=f"{OUT}/source_data_nt_vs_tss_expression.csv",
        figure=f"{OUT}/fig_nt_vs_tss_expression.svg",
    params:
        env=PY_ENV,
        config_file=lambda w: workflow.configfiles[0],
        stem=f"{OUT}/fig_nt_vs_tss_expression",
    log:
        f"{LOGS}/per_guide_knockdown/plot_nt_vs_tss.log",
    threads: threads_for("knockdown_plot")
    resources:
        **resources_for("knockdown_plot"),
    shell:
        """
        mkdir -p $(dirname {log})
        {params.env}
        python -u {input.script} \
            --expression {input.expression} \
            --targeting {input.targeting} \
            --config {params.config_file} \
            --output-summary {output.summary} \
            --output-figure {params.stem} \
            &> {log}
        """


# =============================================================================
# Step 7e -- the combined figure
# =============================================================================
rule plot_pct_knockdown_bar:
    """Bars from the per-gRNA dots, asterisks and ordering from the CP10K test."""
    input:
        summary=f"{OUT}/source_data_nt_vs_tss_expression.csv",
        per_guide=f"{OUT}/source_data_per_guide_tss_knockdown.csv",
        script=CODE_DIR / "scripts/07e_plot_pct_knockdown_bar.py",
    output:
        figure=f"{OUT}/fig_pct_knockdown_bar.svg",
        source_data=f"{OUT}/source_data_pct_knockdown_bar.csv",
    params:
        env=PY_ENV,
        config_file=lambda w: workflow.configfiles[0],
        stem=f"{OUT}/fig_pct_knockdown_bar",
    log:
        f"{LOGS}/per_guide_knockdown/plot_pct_knockdown_bar.log",
    threads: threads_for("knockdown_plot")
    resources:
        **resources_for("knockdown_plot"),
    shell:
        """
        mkdir -p $(dirname {log})
        {params.env}
        python -u {input.script} \
            --summary {input.summary} \
            --per-guide {input.per_guide} \
            --config {params.config_file} \
            --output-figure {params.stem} \
            --output-source-data {output.source_data} \
            &> {log}
        """
