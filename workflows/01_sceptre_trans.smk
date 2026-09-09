# Stage 01 -- SCEPTRE differential-expression testing.
#
# Takes the feature-corrected CellRanger matrices from stage 00 and produces
# per-variant SCEPTRE results: a calibration check on negative-control pairs, a
# power check on TSS positive controls, and the discovery analysis.
#
# Variants are declared in config under sceptre.variants; `trans_v2` is the
# canonical run. The variant controls the pair set and the association
# covariates only -- gRNA-to-cell assignment is identical across variants,
# which is why one assignment underlies every result here.
#
# The heavy lifting is done by the Nextflow implementation of the SCEPTRE
# pipeline (github.com/timothy-barry/sceptre-pipeline), which this stage wraps
# rather than reimplements. That rule needs `nextflow` and `java` on PATH.
#
# Run with:
#   CONFIG=config/config.20260408_ipscvic_300k.yaml \
#   WORKFLOW=workflows/01_sceptre_trans.smk ./submit.sh
#
# Note on runtime: the trans discovery analysis submits thousands of Nextflow
# subjobs and took ~12 h wall clock for 2.75M pairs. The wrapper rule mostly
# sits idle waiting on them, so it asks for little memory and a long walltime.

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
    get_sublibraries,
    print_path_configuration,
)

resolve_config_paths(config, CODE_DIR)

RESULTS = get_results_path(config=config)
LOGS = get_logs_path(config=config)
SUBLIBRARIES = get_sublibraries(config)

SCEPTRE = config["sceptre"]
ONDISC = f"{RESULTS}/{SCEPTRE['ondisc_dir']}"
CELLRANGER = f"{RESULTS}/{SCEPTRE['cellranger_export']}"
VARIANTS = SCEPTRE["variants"]
RUN_VARIANTS = SCEPTRE["run_variants"]

unknown = set(RUN_VARIANTS) - set(VARIANTS)
if unknown:
    raise ValueError(
        f"sceptre.run_variants names undeclared variants: {sorted(unknown)}; "
        f"declared: {sorted(VARIANTS)}"
    )

print_path_configuration(config)
print(f"sceptre variants : {', '.join(RUN_VARIANTS)}")

shell.prefix(get_base_shell_prefix(CODE_DIR))

PY_ENV = get_python_env_prefix(config)
R_ENV = get_r_env_prefix(config)


wildcard_constraints:
    variant="|".join(VARIANTS),
    sublibrary="|".join(SUBLIBRARIES),


def resources_for(rule_name):
    return {
        key: get_resources(config, rule_name, key)
        for key in ("mem_mb", "runtime", "slurm_partition")
    }


def threads_for(rule_name):
    return get_resources(config, rule_name, "threads")


def variant_opt(variant, key, default=None):
    return VARIANTS[variant].get(key, default)


def barcode_files():
    return expand(
        f"{CELLRANGER}/{{sublibrary}}/filtered_feature_bc_matrix/barcodes.tsv.gz",
        sublibrary=SUBLIBRARIES,
    )


def cellranger_dirs():
    """Import order matters: it defines SCEPTRE's cell indexing."""
    return [
        f"{CELLRANGER}/{sublibrary}/filtered_feature_bc_matrix"
        for sublibrary in SUBLIBRARIES
    ]


rule all:
    input:
        f"{ONDISC}/sceptre_extra_covariates.csv",
        f"{RESULTS}/sceptre/validation/cell_ordering_report.txt",
        expand(
            f"{RESULTS}/sceptre/{{variant}}/sceptre_{{variant}}_power_check.tsv",
            variant=RUN_VARIANTS,
        ),
        expand(
            f"{RESULTS}/sceptre/{{variant}}/sceptre_{{variant}}_discovery_analysis.tsv",
            variant=RUN_VARIANTS,
        ),


# =============================================================================
# Step 6a -- import into SCEPTRE
# =============================================================================
rule import_into_sceptre:
    """Build the ondisc-backed sceptre object from the stage-00 export.

    Depends on stage 00's feature-consistency report: import_data_from_cellranger
    indexes features positionally across directories, so identical
    features.tsv.gz files are a correctness precondition, not a nicety.
    """
    input:
        matrices=expand(
            f"{CELLRANGER}/{{sublibrary}}/filtered_feature_bc_matrix/matrix.mtx.gz",
            sublibrary=SUBLIBRARIES,
        ),
        features_ok=f"{CELLRANGER}/feature_consistency_report.txt",
        guide_targets=SCEPTRE["guide_targets_file"],
        script=CODE_DIR / "scripts/06a_import_into_sceptre.R",
    output:
        sceptre_object=f"{ONDISC}/{SCEPTRE['base_object']}",
        gene_odm=f"{ONDISC}/gene.odm",
        grna_odm=f"{ONDISC}/grna.odm",
    params:
        env=R_ENV,
        dirs=",".join(cellranger_dirs()),
        moi=SCEPTRE["moi"],
        out_dir=ONDISC,
    log:
        f"{LOGS}/import_into_sceptre/import_into_sceptre.log",
    threads: threads_for("sceptre_import")
    resources:
        **resources_for("sceptre_import"),
    shell:
        """
        mkdir -p $(dirname {log})
        {params.env}
        Rscript {input.script} \
            --cellranger-dirs {params.dirs} \
            --guide-targets {input.guide_targets} \
            --moi {params.moi} \
            --output-dir {params.out_dir} \
            --output {output.sceptre_object} \
            &> {log}
        """


# =============================================================================
# Step 5 -- extra per-cell covariates
# =============================================================================
rule extract_sceptre_covariates:
    """Pull cell-state covariates out of the clustered AnnData, in import order.

    Row i of the CSV is SCEPTRE cell index i; the join downstream is positional.
    """
    input:
        barcodes=barcode_files(),
        script=CODE_DIR / "scripts/05_extract_sceptre_covariates.py",
    output:
        covariates=f"{ONDISC}/sceptre_extra_covariates.csv",
    params:
        env=PY_ENV,
        config_file=lambda w: workflow.configfiles[0],
        sublibraries=" ".join(SUBLIBRARIES),
    log:
        f"{LOGS}/extract_sceptre_covariates/extract_sceptre_covariates.log",
    threads: threads_for("sceptre_covariates")
    resources:
        **resources_for("sceptre_covariates"),
    shell:
        """
        mkdir -p $(dirname {log})
        {params.env}
        python -u {input.script} \
            --barcode-files {input.barcodes} \
            --sublibraries {params.sublibraries} \
            --output {output.covariates} \
            --config {params.config_file} \
            &> {log}
        """


rule validate_cell_ordering:
    """Check the positional cell contract the covariate join depends on."""
    input:
        covariates=f"{ONDISC}/sceptre_extra_covariates.csv",
        barcodes=barcode_files(),
        script=CODE_DIR / "scripts/06c_validate_cell_ordering.py",
    output:
        report=f"{RESULTS}/sceptre/validation/cell_ordering_report.txt",
    params:
        env=PY_ENV,
        config_file=lambda w: workflow.configfiles[0],
        sublibraries=" ".join(SUBLIBRARIES),
        checks="barcode,covariates",
    log:
        f"{LOGS}/validate_cell_ordering/validate_cell_ordering.log",
    threads: threads_for("sceptre_covariates")
    resources:
        **resources_for("sceptre_covariates"),
    shell:
        """
        mkdir -p $(dirname {log})
        {params.env}
        python -u {input.script} \
            --covariates {input.covariates} \
            --barcode-files {input.barcodes} \
            --sublibraries {params.sublibraries} \
            --checks {params.checks} \
            --output {output.report} \
            --config {params.config_file} \
            &> {log}
        """


# =============================================================================
# Step 6b -- analysis parameters, per variant
# =============================================================================
rule set_analysis_parameters:
    """Attach covariates, build the pair sets, and set the analysis parameters.

    The association formula is passed explicitly. Leaving it to sceptre's
    auto-construction would silently drop every continuous covariate with >= 15
    distinct values -- i.e. all four cell-state covariates.
    """
    input:
        sceptre_object=f"{ONDISC}/{SCEPTRE['base_object']}",
        gene_odm=f"{ONDISC}/gene.odm",
        grna_odm=f"{ONDISC}/grna.odm",
        covariates=f"{ONDISC}/sceptre_extra_covariates.csv",
        ordering_ok=f"{RESULTS}/sceptre/validation/cell_ordering_report.txt",
        script=CODE_DIR / "scripts/06b_set_analysis_parameters.R",
    output:
        sceptre_object=f"{ONDISC}/sceptre_object_{{variant}}.rds",
        summary=f"{RESULTS}/sceptre/{{variant}}/set_analysis_parameters_summary.txt",
    params:
        env=R_ENV,
        extra=lambda w: ",".join(variant_opt(w.variant, "extra_covariates", [])),
        categorical=lambda w: ",".join(
            variant_opt(w.variant, "categorical_covariates", [])
        ),
        pair_type=lambda w: variant_opt(w.variant, "pair_type", "trans"),
        distance=lambda w: variant_opt(w.variant, "distance_threshold", 1000000),
        protein_coding=lambda w: (
            "--protein-coding-only" if variant_opt(w.variant, "protein_coding_only")
            else ""
        ),
        gencode=SCEPTRE["gencode_gtf"],
        side=SCEPTRE["analysis"]["side"],
        resampling=SCEPTRE["analysis"]["resampling_mechanism"],
        control_group=SCEPTRE["analysis"]["control_group"],
        integration=SCEPTRE["analysis"]["grna_integration_strategy"],
        alpha=SCEPTRE["analysis"]["multiple_testing_alpha"],
    log:
        f"{LOGS}/set_analysis_parameters/{{variant}}.log",
    threads: threads_for("sceptre_setparams")
    resources:
        **resources_for("sceptre_setparams"),
    shell:
        """
        mkdir -p $(dirname {log})
        {params.env}
        Rscript {input.script} \
            --sceptre-object {input.sceptre_object} \
            --response-odm {input.gene_odm} \
            --grna-odm {input.grna_odm} \
            --covariates {input.covariates} \
            --extra-covariates "{params.extra}" \
            --categorical-covariates "{params.categorical}" \
            --pair-type {params.pair_type} \
            --distance-threshold {params.distance} \
            {params.protein_coding} \
            --gencode-gtf {params.gencode} \
            --side {params.side} \
            --resampling-mechanism {params.resampling} \
            --control-group {params.control_group} \
            --grna-integration-strategy {params.integration} \
            --multiple-testing-alpha {params.alpha} \
            --output {output.sceptre_object} \
            --summary {output.summary} \
            &> {log}
        """


# =============================================================================
# Step 6c -- run the SCEPTRE Nextflow pipeline
# =============================================================================
rule run_sceptre_pipeline:
    """Wrap `nextflow run timothy-barry/sceptre-pipeline`.

    gRNA assignment happens inside this pipeline and does NOT inherit the
    association formula set above -- assign_grnas() builds its own default
    (transcriptome depth + gRNA depth + batch). Set
    sceptre.pipeline.grna_assignment_formula to an .rds holding a formula
    object to override that; leaving it empty keeps sceptre's default, which is
    what every published run used.
    """
    input:
        sceptre_object=f"{ONDISC}/sceptre_object_{{variant}}.rds",
        gene_odm=f"{ONDISC}/gene.odm",
        grna_odm=f"{ONDISC}/grna.odm",
    output:
        summary=f"{RESULTS}/sceptre/{{variant}}/outputs/analysis_summary.txt",
        power=f"{RESULTS}/sceptre/{{variant}}/outputs/results_run_power_check.rds",
        discovery=f"{RESULTS}/sceptre/{{variant}}/outputs/results_run_discovery_analysis.rds",
        assignment=f"{RESULTS}/sceptre/{{variant}}/outputs/grna_assignment_matrix.rds",
    params:
        env=R_ENV,
        out_dir=f"{RESULTS}/sceptre/{{variant}}/outputs",
        run_dir=f"{RESULTS}/sceptre/{{variant}}/nextflow",
        work_dir=lambda w: f"{SCEPTRE['pipeline']['work_dir_base']}_{w.variant}",
        repo=SCEPTRE["pipeline"]["repo"],
        revision=SCEPTRE["pipeline"]["revision"],
        assignment_method=SCEPTRE["pipeline"]["grna_assignment_method"],
        assignment_formula=lambda w: (
            f"--grna_assignment_formula {SCEPTRE['pipeline']['grna_assignment_formula']}"
            if SCEPTRE["pipeline"].get("grna_assignment_formula")
            else ""
        ),
        pair_pod_size=SCEPTRE["pipeline"]["pair_pod_size"],
        grna_pod_size=SCEPTRE["pipeline"]["grna_pod_size"],
        n_calibration_pairs=SCEPTRE["pipeline"]["n_calibration_pairs"],
        nxf_opts=SCEPTRE["pipeline"]["nxf_opts"],
        java_module=SCEPTRE["pipeline"].get("java_module", ""),
    log:
        f"{LOGS}/run_sceptre_pipeline/{{variant}}.log",
    threads: threads_for("sceptre_pipeline")
    resources:
        **resources_for("sceptre_pipeline"),
    shell:
        """
        mkdir -p $(dirname {log}) {params.out_dir} {params.run_dir}
        {params.env}
        set +u
        if [ -n "{params.java_module}" ]; then module load {params.java_module}; fi
        set -u
        export NXF_OPTS="{params.nxf_opts}"

        # Nextflow writes .nextflow.log and its cache into the launch dir, so
        # give each variant its own to keep -resume caches separate.
        cd {params.run_dir}

        nextflow run {params.repo} -r {params.revision} -resume \
            -w {params.work_dir} \
            --sceptre_object_fp {input.sceptre_object} \
            --response_odm_fp {input.gene_odm} \
            --grna_odm_fp {input.grna_odm} \
            --output_directory {params.out_dir} \
            --grna_assignment_method {params.assignment_method} \
            {params.assignment_formula} \
            --pair_pod_size {params.pair_pod_size} \
            --grna_pod_size {params.grna_pod_size} \
            --n_calibration_pairs {params.n_calibration_pairs} \
            &> {log}
        """


# =============================================================================
# Step 6d -- export results as TSVs
# =============================================================================
rule export_sceptre_results:
    """Annotate the result RDS files with gene symbols and percent knockdown."""
    input:
        power=f"{RESULTS}/sceptre/{{variant}}/outputs/results_run_power_check.rds",
        discovery=f"{RESULTS}/sceptre/{{variant}}/outputs/results_run_discovery_analysis.rds",
        ens2symbol=SCEPTRE["ens2symbol_file"],
        script=CODE_DIR / "scripts/06d_export_sceptre_results.R",
    output:
        power=f"{RESULTS}/sceptre/{{variant}}/sceptre_{{variant}}_power_check.tsv",
        discovery=f"{RESULTS}/sceptre/{{variant}}/sceptre_{{variant}}_discovery_analysis.tsv",
    params:
        env=R_ENV,
        results_dir=f"{RESULTS}/sceptre/{{variant}}/outputs",
        out_dir=f"{RESULTS}/sceptre/{{variant}}",
        prefix=lambda w: f"sceptre_{w.variant}",
    log:
        f"{LOGS}/export_sceptre_results/{{variant}}.log",
    threads: threads_for("sceptre_export")
    resources:
        **resources_for("sceptre_export"),
    shell:
        """
        mkdir -p $(dirname {log})
        {params.env}
        Rscript {input.script} \
            --results-dir {params.results_dir} \
            --ens2symbol {input.ens2symbol} \
            --output-dir {params.out_dir} \
            --prefix {params.prefix} \
            &> {log}
        """
