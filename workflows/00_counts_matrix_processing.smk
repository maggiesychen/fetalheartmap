# Stage 00 -- counts-matrix post-processing.
#
# Picks up the raw kb-python GEX and guide count matrices written by the
# upstream perturb-seq pipeline (github.com/tkzeng/perturb_pipeline) and
# produces:
#
#   * a clustered, annotated combined AnnData (input to cNMF), and
#   * per-sublibrary CellRanger feature-barcode matrices (input to SCEPTRE).
#
# Run with:
#   snakemake -s workflows/00_counts_matrix_processing.smk \
#       --configfile config/config.20260408_ipscvic_300k.yaml --cores 8
# or via ./submit.sh to spread the rules across Slurm.

import sys
from pathlib import Path

CODE_DIR = Path(workflow.snakefile).resolve().parent.parent
sys.path.insert(0, str(CODE_DIR))

from scripts.pipeline_utils import resolve_config_paths
from scripts.snakemake_helpers import (
    get_counts_h5ad,
    get_logs_path,
    get_resources,
    get_results_path,
    get_shell_prefix,
    get_sublibraries,
    get_umi_bounds,
    print_path_configuration,
)

resolve_config_paths(config, CODE_DIR)

RESULTS = get_results_path(config=config)
LOGS = get_logs_path(config=config)
SUBLIBRARIES = get_sublibraries(config)

print_path_configuration(config)

shell.prefix(get_shell_prefix(config, CODE_DIR))


wildcard_constraints:
    sublibrary="|".join(SUBLIBRARIES),


def resources_for(rule_name):
    """Slurm resources for one rule, read from config.resources."""
    return {
        key: get_resources(config, rule_name, key)
        for key in ("mem_mb", "runtime", "slurm_partition")
    }


def threads_for(rule_name):
    return get_resources(config, rule_name, "threads")


# The curated default target: the clustering output that feeds cNMF, plus the
# feature-corrected CellRanger export that feeds SCEPTRE and its validation
# report. The plain (non-feature-corrected) export is available as a separate
# target; see `snakemake ... cellranger_export_plain`.
rule all:
    input:
        f"{RESULTS}/clustering/combined_clustered.h5ad",
        expand(
            f"{RESULTS}/cellranger_format_featurecorrected/{{sublibrary}}/"
            "filtered_feature_bc_matrix/matrix.mtx.gz",
            sublibrary=SUBLIBRARIES,
        ),
        f"{RESULTS}/cellranger_format_featurecorrected/feature_consistency_report.txt",


rule cellranger_export_plain:
    """Optional target: the pre-feature-correction export (step 3b)."""
    input:
        expand(
            f"{RESULTS}/cellranger_format/{{sublibrary}}/"
            "filtered_feature_bc_matrix/matrix.mtx.gz",
            sublibrary=SUBLIBRARIES,
        ),


# =============================================================================
# Step 1 -- UMI thresholding
# =============================================================================
rule umi_filter:
    """Filter one sublibrary to its UMI window and pair GEX with guides.

    The lower bound is the barcode-rank inflection point reported by the
    upstream pipeline's cell calling; both bounds come from sample_info.
    """
    input:
        gex=lambda w: get_counts_h5ad(config, w.sublibrary, "gex"),
        guide=lambda w: get_counts_h5ad(config, w.sublibrary, "guide"),
        script=CODE_DIR / "scripts/01_umi_threshold_filter.py",
    output:
        h5mu=f"{RESULTS}/umi_filtered/{{sublibrary}}_umi_filtered_gex_and_guide.h5mu",
        plot=f"{RESULTS}/umi_filtered/plots/{{sublibrary}}_umi_rank_curve.png",
    params:
        umi_min=lambda w: get_umi_bounds(config, w.sublibrary)[0],
        umi_max=lambda w: get_umi_bounds(config, w.sublibrary)[1],
    log:
        f"{LOGS}/umi_filter/{{sublibrary}}.log",
    threads: threads_for("umi_filter")
    resources:
        **resources_for("umi_filter"),
    shell:
        """
        mkdir -p $(dirname {log})
        python -u {input.script} \
            --gex-input {input.gex} \
            --guide-input {input.guide} \
            --output {output.h5mu} \
            --umi-min {params.umi_min} \
            --umi-max {params.umi_max} \
            --rank-curve-plot {output.plot} \
            &> {log}
        """


# =============================================================================
# Step 2 -- QC filtering
# =============================================================================
rule qc_filter:
    """Per-cell MT/ribo/MAD filters, then per-gene minimum cell/count filters."""
    input:
        h5mu=f"{RESULTS}/umi_filtered/{{sublibrary}}_umi_filtered_gex_and_guide.h5mu",
        script=CODE_DIR / "scripts/02_qc_metrics_filter.py",
    output:
        h5mu=f"{RESULTS}/qc_filtered/{{sublibrary}}_qc_filtered_gex_and_guide.h5mu",
    params:
        config_file=lambda w: workflow.configfiles[0],
        plots_dir=f"{RESULTS}/qc_filtered/plots/{{sublibrary}}",
    log:
        f"{LOGS}/qc_filter/{{sublibrary}}.log",
    threads: threads_for("qc_filter")
    resources:
        **resources_for("qc_filter"),
    shell:
        """
        mkdir -p $(dirname {log})
        python -u {input.script} \
            --input {input.h5mu} \
            --output {output.h5mu} \
            --config {params.config_file} \
            --plots-dir {params.plots_dir} \
            &> {log}
        """


# =============================================================================
# Step 3a -- concatenate sublibraries
# =============================================================================
rule concatenate_sublibraries:
    """Outer-join all sublibraries into one MuData, barcodes suffixed by sublib."""
    input:
        h5mus=expand(
            f"{RESULTS}/qc_filtered/{{sublibrary}}_qc_filtered_gex_and_guide.h5mu",
            sublibrary=SUBLIBRARIES,
        ),
        script=CODE_DIR / "scripts/03a_concatenate_mudata.py",
    output:
        h5mu=f"{RESULTS}/combined/combined_all_sublibraries.h5mu",
    params:
        config_file=lambda w: workflow.configfiles[0],
        sublibraries=" ".join(SUBLIBRARIES),
    log:
        f"{LOGS}/concatenate/concatenate_sublibraries.log",
    threads: threads_for("concatenate")
    resources:
        **resources_for("concatenate"),
    shell:
        """
        mkdir -p $(dirname {log})
        python -u {input.script} \
            --inputs {input.h5mus} \
            --sublibraries {params.sublibraries} \
            --output {output.h5mu} \
            --config {params.config_file} \
            &> {log}
        """


# =============================================================================
# Step 3b -- CellRanger export, per-sublibrary gene sets
# =============================================================================
rule convert_to_cellranger:
    """Write one sublibrary as barcodes/features/matrix, genes above guides.

    Each sublibrary keeps whichever genes survived its own per-gene filters, so
    features.tsv.gz differs between sublibraries. Use the feature-corrected
    rules below for anything that reads the sublibraries as a set.
    """
    input:
        h5mu=f"{RESULTS}/qc_filtered/{{sublibrary}}_qc_filtered_gex_and_guide.h5mu",
        script=CODE_DIR / "scripts/03c_convert_to_cellranger.py",
    output:
        barcodes=f"{RESULTS}/cellranger_format/{{sublibrary}}/"
        "filtered_feature_bc_matrix/barcodes.tsv.gz",
        features=f"{RESULTS}/cellranger_format/{{sublibrary}}/"
        "filtered_feature_bc_matrix/features.tsv.gz",
        matrix=f"{RESULTS}/cellranger_format/{{sublibrary}}/"
        "filtered_feature_bc_matrix/matrix.mtx.gz",
    params:
        config_file=lambda w: workflow.configfiles[0],
        out_dir=f"{RESULTS}/cellranger_format/{{sublibrary}}/"
        "filtered_feature_bc_matrix",
    log:
        f"{LOGS}/convert_to_cellranger/{{sublibrary}}.log",
    threads: threads_for("cellranger_export")
    resources:
        **resources_for("cellranger_export"),
    shell:
        """
        mkdir -p $(dirname {log})
        python -u {input.script} \
            --input {input.h5mu} \
            --output-dir {params.out_dir} \
            --config {params.config_file} \
            --sublibrary {wildcards.sublibrary} \
            &> {log}
        """


# =============================================================================
# Step 3b-i -- CellRanger export onto a shared gene list
# =============================================================================
rule build_master_gene_list:
    """Union of gene IDs across sublibraries, so all exports share features."""
    input:
        h5mus=expand(
            f"{RESULTS}/qc_filtered/{{sublibrary}}_qc_filtered_gex_and_guide.h5mu",
            sublibrary=SUBLIBRARIES,
        ),
        script=CODE_DIR / "scripts/03b_build_master_gene_list.py",
    output:
        gene_list=f"{RESULTS}/cellranger_format_featurecorrected/master_gene_ids.txt",
    log:
        f"{LOGS}/build_master_gene_list/build_master_gene_list.log",
    threads: threads_for("master_gene_list")
    resources:
        **resources_for("master_gene_list"),
    shell:
        """
        mkdir -p $(dirname {log})
        python -u {input.script} \
            --inputs {input.h5mus} \
            --output {output.gene_list} \
            &> {log}
        """


rule convert_to_cellranger_featurecorrected:
    """Reindex one sublibrary onto the master gene list, zero-filling gaps."""
    input:
        h5mu=f"{RESULTS}/qc_filtered/{{sublibrary}}_qc_filtered_gex_and_guide.h5mu",
        gene_list=f"{RESULTS}/cellranger_format_featurecorrected/master_gene_ids.txt",
        script=CODE_DIR / "scripts/03c_convert_to_cellranger.py",
    output:
        barcodes=f"{RESULTS}/cellranger_format_featurecorrected/{{sublibrary}}/"
        "filtered_feature_bc_matrix/barcodes.tsv.gz",
        features=f"{RESULTS}/cellranger_format_featurecorrected/{{sublibrary}}/"
        "filtered_feature_bc_matrix/features.tsv.gz",
        matrix=f"{RESULTS}/cellranger_format_featurecorrected/{{sublibrary}}/"
        "filtered_feature_bc_matrix/matrix.mtx.gz",
    params:
        config_file=lambda w: workflow.configfiles[0],
        out_dir=f"{RESULTS}/cellranger_format_featurecorrected/{{sublibrary}}/"
        "filtered_feature_bc_matrix",
    log:
        f"{LOGS}/convert_to_cellranger_featurecorrected/{{sublibrary}}.log",
    threads: threads_for("cellranger_export")
    resources:
        **resources_for("cellranger_export"),
    shell:
        """
        mkdir -p $(dirname {log})
        python -u {input.script} \
            --input {input.h5mu} \
            --output-dir {params.out_dir} \
            --config {params.config_file} \
            --sublibrary {wildcards.sublibrary} \
            --master-gene-list {input.gene_list} \
            &> {log}
        """


rule validate_cellranger_features:
    """Fail loudly if the exported sublibraries disagree on features.tsv.gz.

    SCEPTRE indexes features positionally across directories, so a mismatch
    would corrupt every result rather than raise.
    """
    input:
        features=expand(
            f"{RESULTS}/cellranger_format_featurecorrected/{{sublibrary}}/"
            "filtered_feature_bc_matrix/features.tsv.gz",
            sublibrary=SUBLIBRARIES,
        ),
        script=CODE_DIR / "scripts/03d_validate_cellranger_features.py",
    output:
        report=f"{RESULTS}/cellranger_format_featurecorrected/"
        "feature_consistency_report.txt",
        reference=f"{RESULTS}/cellranger_format_featurecorrected/"
        "reference_features.txt",
    params:
        sublibraries=" ".join(SUBLIBRARIES),
    log:
        f"{LOGS}/validate_cellranger_features/validate_cellranger_features.log",
    threads: threads_for("master_gene_list")
    resources:
        **resources_for("master_gene_list"),
    shell:
        """
        mkdir -p $(dirname {log})
        python -u {input.script} \
            --features-files {input.features} \
            --sublibraries {params.sublibraries} \
            --output {output.report} \
            --reference-features {output.reference} \
            &> {log}
        """


# =============================================================================
# Step 4 -- clustering and annotation
# =============================================================================
rule cluster_and_annotate:
    """Normalise, score cell cycle, select HVGs, PCA/UMAP/Leiden, annotate."""
    input:
        h5mu=f"{RESULTS}/combined/combined_all_sublibraries.h5mu",
        s_genes=config["input_paths"]["cell_cycle_s_genes"],
        g2m_genes=config["input_paths"]["cell_cycle_g2m_genes"],
        script=CODE_DIR / "scripts/04_cluster_and_annotate.py",
    output:
        h5ad=f"{RESULTS}/clustering/combined_clustered.h5ad",
    params:
        config_file=lambda w: workflow.configfiles[0],
        plots_dir=f"{RESULTS}/clustering/plots",
    log:
        f"{LOGS}/cluster_and_annotate/cluster_and_annotate.log",
    threads: threads_for("clustering")
    resources:
        **resources_for("clustering"),
    shell:
        """
        mkdir -p $(dirname {log})
        python -u {input.script} \
            --input {input.h5mu} \
            --output {output.h5ad} \
            --config {params.config_file} \
            --plots-dir {params.plots_dir} \
            &> {log}
        """
