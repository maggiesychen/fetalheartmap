# Fetal heart map — iPSC-VIC Perturb-seq analysis

Snakemake workflows for the iPSC-derived valve interstitial cell (VIC) high-MOI
CRISPRi Perturb-seq screen targeting congenital heart disease (CHD) genes.

**Contact:** Maggie Chen (msychen@stanford.edu), Engreitz Lab

**Stages so far:**

| Workflow | What it does |
|---|---|
| `00_counts_matrix_processing.smk` | kb counts → QC → clustering + CellRanger export |
| `01_sceptre_trans.smk` | SCEPTRE calibration / power / discovery, cis and trans |
| `02_per_guide_knockdown.smk` | per-gRNA knockdown of the CHD target panel |

---

## Stage 00 — counts-matrix processing

Picks up the raw per-sublibrary `adata.h5ad` files written by
[`tkzeng/perturb_pipeline`](https://github.com/tkzeng/perturb_pipeline) and
produces two things:

1. a clustered, annotated combined AnnData — the input to cNMF;
2. per-sublibrary CellRanger feature-barcode matrices — the input to SCEPTRE.

Scripts in `scripts/` are numbered in the order they run.

```
kb counts (GEX + guide, per sublibrary)
  │
  ├─ 01_umi_threshold_filter.py            rule: umi_filter
  │     Per-sublibrary UMI floor from the upstream barcode-rank inflection
  │     point, shared ceiling of 75,000. Intersects GEX ∩ guide barcodes.
  │     → umi_filtered/{sublibrary}_umi_filtered_gex_and_guide.h5mu
  │
  ├─ 02_qc_metrics_filter.py               rule: qc_filter
  │     MT ≤ 15%, ribo ≤ 4%, ±4 MAD outliers on log1p total counts and
  │     log1p n_genes, then per-gene filters (≥10 cells, ≥100 counts).
  │     → qc_filtered/{sublibrary}_qc_filtered_gex_and_guide.h5mu
  │
  ├─ 03a_concatenate_mudata.py             rule: concatenate_sublibraries
  │  │    Outer join across sublibraries, barcodes suffixed with the label.
  │  │    → combined/combined_all_sublibraries.h5mu
  │  │
  │  └─ 04_cluster_and_annotate.py         rule: cluster_and_annotate
  │        CP10k + log1p, cell-cycle scoring, 3,000 HVGs (seurat_v3,
  │        batch_key=sublib), 50 PCs / 30 used, k=15 neighbours, UMAP,
  │        Leiden at 0.1/0.25/0.5/1.0, marker scoring and dotplots.
  │        → clustering/combined_clustered.h5ad              [feeds cNMF]
  │
  └─ 03b_build_master_gene_list.py         rule: build_master_gene_list
     │    Union of gene IDs across all sublibraries.
     │    → cellranger_format_featurecorrected/master_gene_ids.txt
     │
     ├─ 03c_convert_to_cellranger.py       rules: convert_to_cellranger,
     │  │                                  convert_to_cellranger_featurecorrected
     │  │  Per-sublibrary barcodes/features/matrix triplets. Counts are
     │  │  rounded to int32 here. Genes stacked above guides
     │  │  ("Gene Expression" / "CRISPR Guide Capture").
     │  │  Without --master-gene-list → cellranger_format/{sublibrary}/…
     │  │  With    --master-gene-list → cellranger_format_featurecorrected/…
     │  │
     └─ 03d_validate_cellranger_features.py  rule: validate_cellranger_features
           Confirms every sublibrary emitted an identical features.tsv.gz.
           → cellranger_format_featurecorrected/…            [feeds SCEPTRE]
```

`03a`/`04` and the `03b`→`03c`→`03d` export branch are independent branches off
`02`, so they run in parallel. `03c` without `--master-gene-list` does not
depend on `03b` and can run at any point after `02`.

These map onto the numbered directories of the original analysis tree as:
`01` → `1.umifiltering`, `02` → `2.qcfiltering`, `03a` →
`3a.concat-old-forclustering`, `03c` → `3b.cellranger_format`, `03b`+`03c`+`03d`
→ `3b-i.cellranger_format_featurecorrected`, `04` → `4.clustering`.

### Why the feature-corrected export exists

`02_qc_metrics_filter.py` applies its per-gene filters one sublibrary at a time,
so the sublibraries end up with different gene sets — for the 2026-04-08 run,
25,703 / 25,266 / 17,354 / 24,522 / 25,451 / 22,960 / 25,234 genes, union
26,070. SCEPTRE's `import_data_from_cellranger` reads a list of matrix
directories and indexes features **positionally**, so mismatched
`features.tsv.gz` files corrupt every result rather than raising. The
`validate_cellranger_features` rule fails loudly if they ever diverge.

---

## Stage 01 — SCEPTRE trans differential expression

Tests every perturbation against every protein-coding gene genome-wide, with a
calibration check on negative-control pairs and a power check on the TSS
positive controls. All published results come from the trans analysis.

```
CellRanger matrices (stage 00, feature-corrected)
  │
  ├─ 05_extract_sceptre_covariates.py   rule: extract_sceptre_covariates
  │     S_score, G2M_score, pct_ribo, pct_mito read out of the clustered
  │     AnnData, in CellRanger import order.
  │     → sceptre_ondisc/sceptre_extra_covariates.csv
  │
  ├─ 06c_validate_cell_ordering.py      rule: validate_cell_ordering
  │     Confirms row i of that CSV is sceptre cell i. The covariate join is
  │     positional with no barcode key, so this check is load-bearing.
  │
  ├─ 06a_import_into_sceptre.R          rule: import_into_sceptre
  │     import_data_from_cellranger(moi="high"), ondisc-backed so the full
  │     240k x 26k matrix never has to be held in memory.
  │     → sceptre_ondisc/{sceptre_object.rds, gene.odm, grna.odm}
  │
  ├─ 06b_set_analysis_parameters.R      rule: set_analysis_parameters
  │     Attaches the covariates, builds the trans pair set, restricts responses
  │     to protein-coding genes, sets the association formula explicitly.
  │     → sceptre_ondisc/sceptre_object_trans_v2.rds
  │
  ├─ nextflow sceptre-pipeline          rule: run_sceptre_pipeline
  │     Wraps `nextflow run timothy-barry/sceptre-pipeline`, which does the
  │     gRNA assignment and all the association testing. ~12 h wall clock.
  │     → sceptre/trans_v2/outputs/results_run_*.rds
  │
  └─ 06d_export_sceptre_results.R       rule: export_sceptre_results
        Adds gene symbols and pct_knockdown = 100 * (1 - 2^log2FC).
        → sceptre/trans_v2/sceptre_trans_v2_*.tsv
```

Covariate sets and pair types are declared under `sceptre.variants`;
`sceptre.run_variants` picks which ones `rule all` builds. The default is
`trans_v2` — cell-state covariates, no Leiden cluster.

### Two things about the model that are easy to get wrong

**The association formula is passed explicitly.** Left to sceptre's
auto-construction, `auto_construct_formula_object()` silently drops any
continuous covariate with ≥ 15 distinct values — which is all four cell-state
covariates — with no warning and nothing in the printed object summary to show
it happened.

**gRNA assignment does not use those covariates, deliberately.**
`assign_grnas()` never reads the association formula; it builds its own default
adjusting for transcriptome depth, gRNA depth and batch. Cell-state covariates
belong in the association model because they describe the transcriptomic
response rather than gRNA capture efficiency, and are potentially downstream of
the perturbation itself. To override, point
`sceptre.pipeline.grna_assignment_formula` at an `.rds` holding a formula
object. Full write-up in the analysis directory at
`6.sceptre/readmes/20260908-methods.md`.

Also worth knowing: sceptre's built-in `response_p_mito` covariate is
identically zero here, because the feature-corrected export writes Ensembl IDs
into the feature *name* column so no `MT-` symbols are visible to sceptre. The
`pct_mito` supplied by `05_extract_sceptre_covariates.py` is the working
mitochondrial covariate.

---

## Stage 02 — per-gRNA knockdown

SCEPTRE reports one fold change per target, pooling that target's gRNAs. This
stage recomputes knockdown one gRNA at a time over the CHD gene panel, by two
independent estimators — they share the counts and the gRNA assignment but
nothing else, so their agreement is evidence the knockdown is real rather than
a modelling artefact.

```
sceptre object + gRNA assignment matrix (stage 01)
  │
  ├─ 07a_per_guide_fold_changes_poisson.R      x2, independent
  │     Poisson GLM log2FC against sceptre's own null model
  │     (sceptre book §10.4). One script, two guide sets:
  │       --guide-set targeting      → the TSS gRNAs
  │       --guide-set non-targeting  → the empirical null
  │
  ├─ 07b_target_gene_expression_per_guide.R
  │     Model-free CP10K pseudobulk vs a per-gene non-targeting pool
  │     (NT cells carrying none of that gene's TSS gRNAs).
  │
  ├─ 07c_plot_tss_knockdown.py     Mann-Whitney on log2FC → the per-gRNA dots
  ├─ 07d_plot_nt_vs_tss.py         Mann-Whitney on CP10K  → stars + ordering,
  │                                and the correlation between the two estimators
  │
  └─ 07e_plot_pct_knockdown_bar.py
        The combined figure. Bars are the mean of each gene's own dots.
        → per_guide_knockdown/fig_pct_knockdown_bar.svg
```

`07a` (targeting), `07a` (non-targeting) and `07b` are mutually independent
despite the numbering, and all three read only stage-01 outputs.

**Percent knockdown is `100 * (1 - 2^log2FC)`, a concave transform.** Do not
average percentages across gRNAs and compare that to the percentage of a mean
fold change — average log2 fold changes and convert once. `07e` asserts the
resulting Jensen inequality, along with four other guards (excluded genes must
exist, every plotted gene must have dots, each bar must equal the mean of its
dots to 1e-9, and the gene-label geometry must match the figure spec).

---

## Software

Snakemake 7.32.4 (the `--cluster` interface, not the 8+ `--executor slurm`
plugin) with Python 3.10.19 and:

| Package | Version | Package | Version |
|---|---:|---|---:|
| anndata | 0.11.4 | numpy | 2.2.6 |
| h5py | 3.15.1 | openpyxl | 3.1.5 |
| matplotlib | 3.10.8 | pandas | 2.3.3 |
| mudata | 0.3.3 | pynndescent | 0.6.0 |
| muon | 0.1.7 | python-igraph | 1.0.0 |
| natsort | 8.4.0 | PyYAML | 6.0.3 |
| numba | 0.64.0 | scanpy | 1.11.5 |
| scikit-learn | 1.7.2 | scikit-misc | 0.5.2 |
| scipy | 1.15.3 | seaborn | 0.13.2 |
| umap-learn | 0.5.11 | | |

Two environment files, so the orchestrator and the analysis stack solve
independently:

```bash
conda activate mamba-env          # see the lab manual, 03-software.md
mamba env create -f environment.yml        # -> fetalheartmap
mamba env create -f envs/snakemake.yml     # -> fetalheartmap-snakemake
```

`leidenalg` is deliberately absent: `04_cluster_and_annotate.py` uses
`sc.tl.leiden(flavor="igraph")`, which goes through `python-igraph`.

The versions above are transcribed from the environment the published run used
(`$OAK/Users/msychen/conda/envs/perturbprocessing`, which `config.conda.env`
points at). CI rebuilds `environment.yml` from scratch on every push, so a
solve that has drifted shows up as a red build rather than as a surprise months
later.

`envs/snakemake.yml` pins Python 3.11. An older Snakemake environment exists at
`~/.conda/envs/snakemake` on Python 3.8, which still parses these workflows but
emits `importlib.metadata has no attribute packages_distributions` on startup
and picks up `~/.local/lib/python3.8/site-packages`. Build the environment from
`envs/snakemake.yml` instead, or prefix commands with `PYTHONNOUSERSITE=1`.

---

## Layout

```
config/       one YAML per analysis; all paths and parameters live here
references/   small inputs that ship with the code (sample table, gene lists,
              gRNA targets, Ensembl->symbol map)
scripts/      step scripts, numbered in run order, plus shared helpers
              (pipeline_utils.py, snakemake_helpers.py, r_utils.R, figstyle.py)
workflows/    numbered Snakemake stages (*.smk)
envs/         the orchestrator environment
tests/        contract tests that run without the sequencing data
submit.sh     Slurm submission wrapper
```

Script numbers run across the whole analysis, not per stage: `01`–`04` are
stage 00, `05`–`06` stage 01, `07` stage 02.

R scripts take their arguments through `scripts/r_utils.R`, a ~100-line base-R
parser. That is deliberate: the shared R library these scripts run against has
neither `optparse` nor `getopt`, and its `argparse` is a reticulate wrapper
around Python's, which would make every R script depend on a working Python
inside R. The SCEPTRE pipeline's own scripts use plain `commandArgs()` for the
same reason.

Paths under `config/` and `references/` resolve relative to this code
directory. Runtime paths (`scratch_base`, `results_base`) resolve relative to
wherever Snakemake is launched; set them to absolute paths if you prefer.

---

## Required inputs

### 1. kb-python count matrices

`input_paths.kb_output_base` must point at the upstream pipeline's output tree.
For each sublibrary, the workflow reads:

```
{kb_output_base}/{counts_dir}/{counts_subdir}/adata.h5ad
```

with `counts_dir` and `counts_subdir` coming from the sample table. GEX and
guide libraries each get a row.

### 2. Sample table — `references/sample_info.<analysis>.tsv`

One row per library. Required columns:

| Column | Meaning |
|---|---|
| `sublibrary` | Sublibrary label; the unit of work, and the barcode suffix |
| `sample_id` | Unique library identifier |
| `sample_type` | `gex` or `guide` |
| `counts_dir` | Sublibrary directory under `kb_output_base` |
| `counts_subdir` | Path to the counts directory holding `adata.h5ad` |
| `min_umi_threshold` | UMI floor (gex rows); the barcode-rank inflection point |
| `max_umi_threshold` | UMI ceiling (gex rows); blank means no ceiling |
| `exclude` | `True` to drop a library without editing the workflow |

Every sublibrary needs exactly one `gex` and one `guide` row — `tests/` enforces
this.

### 3. Ensembl ↔ symbol dictionaries

`input_paths.ensembl_to_symbol` / `symbol_to_ensembl` are GENCODE v43 pickled
`dict` objects (`.npy`). They are too large to ship here and currently point at
a shared copy under `$OAK/Users/opushkar/genome/`. Set them to `""` to skip
symbol mapping; MT/ribo annotation and the marker panels then do not resolve.

### 4. Cell-cycle gene lists

Shipped in `references/` (Tirosh et al. 2016, 42 S and 54 G2M genes).

### 5. gRNA target table — `references/guide_targets.*.tsv`

`grna_id` / `grna_target` (plus `chr`/`start`/`end`), 3,026 gRNAs: 2,826
targeting across 179 targets, and 200 non-targeting. `grna_target` is an
Ensembl gene ID for a TSS target, an element name for an enhancer, or the
literal `non-targeting`. Shipped in `references/`.

### 6. Clustered AnnData, for the SCEPTRE covariates

`sceptre.covariates.clustered_h5ad` supplies `S_score`, `G2M_score`,
`pct_ribo` and `pct_mito`. Too large to ship (34-50 GB). The published runs
read the cell-cycle-regressed clustering; those four covariates are
byte-identical between it and the no-regression output (verified 2026-09-08),
so either reproduces them. Only `leiden_res0_25` differs between the two — 6
clusters vs 5 — and the default variant does not use it.

### 7. GENCODE annotation and the R toolchain

`sceptre.gencode_gtf` (GENCODE v43, 53 MB) restricts trans discovery pairs to
protein-coding responses. `sceptre.r_env` names the R modules and the shared
library providing `sceptre` 0.10.2, `ondisc` 1.2.0, `data.table` and
`rtracklayer`. `nextflow` and `java` must be on PATH for the pipeline rule.

Note that library has neither `optparse` nor `getopt`, which is why the R
scripts parse arguments through `scripts/r_utils.R` instead.

### 8. CHD gene panel — `references/chd_genes.suppt8.tsv`

26 genes with the flags stage 02 uses to pick the main-figure subset. Shipped
in `references/`.

---

## Running

```bash
conda activate fetalheartmap-snakemake
export CONFIG=config/config.20260408_ipscvic_300k.yaml

./submit.sh --dry-run        # check the DAG
./submit.sh --jobs 20        # one Slurm job per rule
```

Per-rule memory, walltime, threads and partition come from the `resources:`
block of the config; `submit.sh` maps them onto `sbatch`. Cluster stdout/stderr
land in `{results_base}/logs/cluster/`, per-rule script logs in
`{results_base}/logs/{rule}/`.

`submit.sh` defaults to stage 00. Select a later stage with `WORKFLOW`, and run
them in order — each depends on the one before:

```bash
export CONFIG=config/config.20260408_ipscvic_300k.yaml

WORKFLOW=workflows/00_counts_matrix_processing.smk ./submit.sh
WORKFLOW=workflows/01_sceptre_trans.smk            ./submit.sh
WORKFLOW=workflows/02_per_guide_knockdown.smk      ./submit.sh
```

Stage 01's Nextflow rule asks for a 3-day walltime and only 10 GB, because it
mostly waits on the subjobs it submits itself. Do not assume `-resume` will
save you a re-run: changing the analysis parameters changes that rule's input
hash, and everything downstream of gRNA assignment re-runs with it.

To run everything inside a single allocation instead:

```bash
conda activate fetalheartmap-snakemake
snakemake -s workflows/00_counts_matrix_processing.smk \
  --configfile config/config.20260408_ipscvic_300k.yaml --cores 8
```

The default target builds the clustering output and the feature-corrected
export. The pre-correction export is a separate target:

```bash
./submit.sh cellranger_export_plain
```

### Tests

```bash
conda activate fetalheartmap
pytest -q
```

These check the config/sample-table/workflow contract — required keys, that
`dotplot_resolution` is one of `leiden_resolutions`, that every rule declares a
log and resources, that every script is argparse-driven and referenced by a
workflow, and that no absolute path has crept back into the code. They need no
sequencing data, so CI runs them on every push.

---

## Adapting this to a new experiment

1. Copy `config/config.20260408_ipscvic_300k.yaml`, change `analysis_name`,
   `results_base` and `input_paths.kb_output_base`.
2. Copy `references/sample_info.20260408_ipscvic_300k.tsv`, fill in your
   sublibraries and their UMI bounds, and point the new config at it.
3. Set `qc_filtering.reference` to `mouse` if applicable, and replace
   `clustering.marker_panels` with panels for your cell types.
4. `pytest -q`, then `./submit.sh --dry-run`.

No code changes should be needed.

---

## Relationship to the other lab repos

| Repo | Overlap with this stage |
|---|---|
| [`tkzeng/perturb_pipeline`](https://github.com/tkzeng/perturb_pipeline) | Produces this stage's inputs. FASTQ → kallisto/bustools → cell calling → per-sublibrary `adata.h5ad`. Not vendored here. |
| [`EngreitzLab/telohaec-genomewide-perturb-seq-analysis`](https://github.com/EngreitzLab/telohaec-genomewide-perturb-seq-analysis) | Same conventions (numbered `.smk` stages, `config/`, `references/`, `scripts/`, `pipeline_utils` + `snakemake_helpers`), and this repo follows them. The *analyses* differ — see below. |


## Changes made when importing these scripts into the repo

The scripts came from `20260408-postprocessing/` as loose `sbatch` submitters
with hardcoded absolute paths. **No analysis logic, parameter or threshold was
changed.** What changed:

- **Snakemake replaces the `sbatch` loops.** The `for SUBLIB in 1 2 3 4 5 6 7`
  submitters and the `--dependency=afterok` chain that built the 3b-i master
  gene list are now rules. The chain became two rules
  (`build_master_gene_list`, then a fan-out), which removes the dependency on
  sublibrary 1 running first.
- **All paths moved into `config/` and `references/`.** Nothing in `scripts/` or
  `workflows/` contains an absolute path; `tests/` enforces that.
- **3b and 3b-i are one script.** `03c_convert_to_cellranger.py
  --master-gene-list` is the 3b-i behaviour; without the flag it is 3b. The
  duplicated file also drifted — 3b lacked 3b-i's feature-uniqueness assert,
  which is now applied in both modes.
- **The `reference_features.txt` cross-check became its own rule.** In the
  original it was written by whichever sublibrary ran first and compared against
  by the rest, which races once the jobs run in parallel.
  `validate_cellranger_features` compares all sublibraries after the fact and
  reports per-sublibrary feature counts.
- **`1.umi_threshold_int_layercheck.py` folded into
  `01_umi_threshold_filter.py`** as `--check-layer-sums` and `--round-counts`. It
  was a diagnostic variant (asserts `nascent + ambiguous + mature == X`, rounds
  to int32); nothing downstream consumed its output. Neither flag is set by the
  workflow, so the default path is unchanged.
- **Random seeds are now explicit.** `clustering.random_seed` is passed to PCA,
  the neighbour graph, UMAP, Leiden and gene scoring. The original relied on
  scanpy's defaults.
- **`use_highly_variable=True` → `mask_var="highly_variable"`** in `sc.tl.pca`,
  which is the non-deprecated spelling of the same thing.
- **`TEST_MODE` toggles removed** from the two CellRanger scripts.

