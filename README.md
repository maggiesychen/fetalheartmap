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

Still to come: cNMF and program-level DE, and the paper figures.

---

## What this stage does

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

### Re-running the 2026-04-08 analysis

The step-1 inputs are **no longer on disk**. All 14 `adata.h5ad` files under
`input_paths.kb_output_base` were removed by the `$SCRATCH` 90-day purge; the
directory skeleton survives but the files are gone. `--dry-run` therefore stops
at `umi_filter` with a `MissingInputException`, which is the workflow correctly
reporting missing inputs rather than a configuration error.

To reproduce from raw reads, re-run
[`tkzeng/perturb_pipeline`](https://github.com/tkzeng/perturb_pipeline) to
regenerate those matrices first. The step 1–4 outputs from the original run are
preserved on Oak under
`$OAK/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/`,
so downstream stages can start from there without re-running this one.

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

## Analysis decisions worth knowing before you reuse this

### Clustering does not regress out the cell cycle

`clustering.regress_cell_cycle` is `false`, and the published run's output is
what every downstream step consumes. A cell-cycle-**regressed** variant was
explored in the analysis directory
(`4.clustering/20260420-nocellcycle-clustering.ipynb`, whose filename is
misleading — it *does* call `regress_out(["S_score","G2M_score"])` and
`scale(max_value=10)`) and was not carried forward. Turning the flag on also
scales the data, which changes PCA and every clustering result downstream.

The published run's own sanity check confirmed that the clustered file is
compatible with the cNMF input:

```
Combined-clustered h5ad: 240504 cells x 25788 genes
cNMF input h5ad:         240504 cells x 25788 genes
obs_names identical after known barcode rename: True
var_names identical (order + content):          True
```

### Rounding, not truncating

Counts stay float32 through steps 1–3a and are rounded to int32 only at
CellRanger export. During the original run, step 2 was once re-run
**truncating** to int, which zeroed every sub-1.0 count and cost 52 cells and
1,855 genes in sublibrary 1 alone (45,938 × 25,703 → 45,886 × 23,848). That run
was abandoned. `03c_convert_to_cellranger.py` rounds, and nothing upstream of it
casts.

### `seurat_v3` sees non-integer input

Step 4 inherits float32 counts, so
`highly_variable_genes(flavor="seurat_v3")` logs `expects raw count data, but
non-integers were found`. Values are near-integer so HVG selection is
approximately right, but it is not what `seurat_v3` assumes.

### The ribosomal threshold is aggressive

`qc_filtering.ribo_threshold` is 4%, against 5 in the script this was adapted
from and 12 in the winter 2025 iPSC-VIC screen. It is the single biggest lever
on final cell count.

### MAD outliers are computed before the MT/ribo filters

…and applied after, so the cutoffs reflect the unfiltered population. Standard
practice, but the order matters if you change it.

### Notes that disagree with the code

`countsmatrices/filteringparameters.md` in the analysis directory records a UMI
ceiling of 100,000 for every sublibrary. The run used **75,000**, which is what
the Slurm logs confirm and what `references/sample_info.*.tsv` records. The
lower bounds agree. Treat the sample table as truth.

---

## Relationship to the other lab repos

| Repo | Overlap with this stage |
|---|---|
| [`tkzeng/perturb_pipeline`](https://github.com/tkzeng/perturb_pipeline) | Produces this stage's inputs. FASTQ → kallisto/bustools → cell calling → per-sublibrary `adata.h5ad`. Not vendored here. |
| [`EngreitzLab/telohaec-genomewide-perturb-seq-analysis`](https://github.com/EngreitzLab/telohaec-genomewide-perturb-seq-analysis) | Same conventions (numbered `.smk` stages, `config/`, `references/`, `scripts/`, `pipeline_utils` + `snakemake_helpers`), and this repo follows them. The *analyses* differ — see below. |

The TeloHAEC repo is a good structural template and a genuine methods reference
for the later program-level DE stage, but it is **not** a substitute for this
one:

- **Steps 1–4 (this stage).** TeloHAEC does the equivalent inside
  `00_fastq_to_sublibrary_h5ad.smk` and `01_combined_matrix_and_cnmf.smk`, but
  it calls cells *inside* the workflow (`BarcodeRanks_Inflection` on its own
  kallisto output) rather than taking a UMI floor from an upstream pipeline, and
  it has no CellRanger-format export at all, because it does not use SCEPTRE.
- **Guide assignment and DE (steps 5–7 of the analysis).** TeloHAEC contains no
  SCEPTRE code. It replaces gene-level SCEPTRE testing entirely with matched-cell
  program DE. There is nothing there to point at for this project's cis/trans
  SCEPTRE results.
- **cNMF (steps 8–9).** Both run cNMF, but TeloHAEC uses a modified CPU cNMF
  (`external/cnmf_modified.py`, K=50/150); this project uses a GPU `halsvar`
  torch-cNMF at K=30.
- **Program-level DE.** This is the real overlap. TeloHAEC's
  `scripts/run_matching_de_batch.R` and this project's `MatchedProgramDE` step
  are the same method and the same lineage of script — MatchIt propensity
  matching plus OLS with HC3 robust standard errors. TeloHAEC's
  `02_program_de.smk` is worth mirroring when that stage is added here.

---

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

### Two corrections to what the original scripts did

- **The cell-cycle fallback was not dead code — it ran.** Step 4 looked for
  `s_genes.txt` / `g2m_genes.txt` under `$OAK/Users/opushkar/common_sc`, but the
  files there are named `hs_cell_cycle_s_genes.txt` /
  `hs_cell_cycle_g2m_genes.txt`. The `os.path.exists` check therefore failed and
  the fallback branch executed: it downloaded `sc.datasets.pbmc3k_processed()`
  (23.5 MB, visible in the Slurm stderr) and used the hardcoded in-script Tirosh
  lists. The published run's `S genes found: 42, G2M genes found: 52` matches
  those hardcoded lists, not the files on disk. Those exact lists are now
  shipped as `references/cell_cycle_*.txt`, so the numbers are preserved and the
  compute-node download is gone.
- **Barcode intersection is now deterministic.** Steps 1 and 2 took shared
  barcodes via `list(set(...))`, so cell *order* varied between runs. Both
  modalities were subset with the same list, so GEX and guide stayed aligned and
  results were correct — but the step was not reproducible. It is now `sorted()`.
  A re-run will therefore not byte-match the April 2026 outputs; no re-run of
  the original could have either.
