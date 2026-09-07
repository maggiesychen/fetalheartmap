# Counts matrix processing

Post-processing of kallisto|bustools counts matrices for the
**20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq** experiment (iPSC-VIC high-MOI CRISPRi
screen, 300k cells, 7 sublibraries, day 9).

Picks up from the raw `adata.h5ad` files produced by the perturb-seq pipeline and ends
with (a) a clustered/annotated combined AnnData and (b) CellRanger-format matrices for
guide assignment. The upstream pipeline itself is not included here.

Guide assignment (step 5), SCEPTRE (steps 6–7) and cNMF (steps 8–10) live elsewhere.

## Pipeline

```
kb counts (GEX + guide, per sublibrary, on $SCRATCH)
  │
  ├─ 1. UMI filtering ──────────── 1.umi_threshold_batch.sh → 1.umi_threshold_float.py
  │     Per-sublib UMI floor from the TZ-pipeline barcode-rank inflection point,
  │     shared ceiling of 75,000. Intersects GEX ∩ guide barcodes, writes .h5mu.
  │     → 1.umifiltering/output/300k_sublib{N}_umi_filtered_gex_and_guide.h5mu
  │
  ├─ 2. QC filtering ───────────── 2.qc_metrics_batch.sh → 2.qc_metrics_filter.py
  │     MT ≤ 15%, ribo ≤ 4%, ±4 MAD outliers on log1p total counts and log1p n_genes,
  │     then per-gene filters (min 10 cells, min 100 counts).
  │     → 2.qcfiltering/float_output/qc_filtered_data/300k_sublib{N}_qc_filtered_*.h5mu
  │
  ├─ 3a. Concatenate ──────────── 3a.concat_batch.sh → 3a.concatenate_mudata.py
  │  │    Outer join across all 7 sublibs, barcodes suffixed with sublib label.
  │  │    → 3a.concat-old-forclustering/300k_combined_all_sublibs_newguidelist.h5mu
  │  │      (240,504 cells)
  │  │
  │  └─ 4. Cluster + annotate ─── 4.cluster_annotate_batch.sh → 4.cluster_and_annotate.py
  │        CP10k + log1p, cell-cycle scoring, 3,000 HVGs (seurat_v3, batch_key=sublib),
  │        50 PCs / 30 used, k=15 neighbours, UMAP, Leiden at 0.1/0.25/0.5/1.0.
  │        No cell-cycle regression, no scaling — see "Clustering variants" below.
  │        → 4.clustering/clustered_data/300k_combined_clustered.h5ad
  │          (240,504 cells x 25,788 genes) — this file feeds cNMF
  │
  ├─ 3b.  CellRanger format ───── 3b.converttocellranger_batch.sh
  │        Per-sublib barcodes/features/matrix triplets. Rounds counts to int32 here.
  │        Genes stacked above guides ("Gene Expression" / "CRISPR Guide Capture").
  │        → 3b.cellranger_format/sublib{N}/filtered_feature_bc_matrix/
  │
  └─ 3b-i. CellRanger format, feature-corrected ──────────────────────────────────
           3b-i.converttocellranger_featurecorrected_batch.sh
           Same as 3b, but first builds a master gene list (union across all 7
           sublibs, saved to disk) and reindexes every sublib onto it, padding
           missing genes with zeros. Needed because step 2 applies its per-gene
           filters per sublibrary, so the sublibs otherwise end up with different
           gene sets (25,703 / 25,451 / 25,266 / ... genes). Sublib 1 runs first to
           build the list; sublibs 2–7 follow via `--dependency=afterok`.
           → 3b-i.cellranger_format_featurecorrected/sublib{N}/filtered_feature_bc_matrix/
```

Steps 3a→4 and 3b/3b-i are independent branches off step 2.

## Running

Submit from inside this folder — the batch scripts resolve the Python scripts relative
to their own location, so a fresh clone works without editing paths.

```bash
bash   1.umi_threshold_batch.sh                             # loop submitter, 7 jobs
bash   2.qc_metrics_batch.sh                                # loop submitter, 7 jobs
sbatch 3a.concat_batch.sh                                   # single job
sbatch 4.cluster_annotate_batch.sh                          # single job
bash   3b.converttocellranger_batch.sh                      # loop submitter, 7 jobs
bash   3b-i.converttocellranger_featurecorrected_batch.sh   # chained, 1 then 2-7
```

All jobs expect the `perturbprocessing` conda env. Input/output paths are absolute and
still point at the original Oak locations — edit the `--- Parameters ---` block at the
top of each Python script to retarget.

## Clustering variants — read this before reusing step 4

There are two clustering paths in `4.clustering/`, and their names are misleading:

| File | Cell-cycle regression? | Used? |
|---|---|---|
| `4.cluster_and_annotate.py` (this repo) | **No** — no `regress_out`, no `scale` | **Yes** |
| `20260420-nocellcycle-clustering.ipynb` | **Yes** — `regress_out(["S_score","G2M_score"])` + `scale(max_value=10)` at cell 9 | No (superseded) |

Despite its filename, `20260420-nocellcycle-clustering.ipynb` is the cell-cycle-**regressed**
version; it produced `4.clustering/plots/clustering_ccregressed/` and was not carried forward.

`4.cluster_and_annotate.py` is the no-regression pipeline, and its output is what everything
downstream uses. `4.clustering/clustering_noccregression/` does **not** re-cluster — it loads
`300k_combined_clustered.h5ad`, reuses the stored PCA/neighbours/UMAP/Leiden, and layers on
marker genes, curated dotplots and annotation. Its own sanity check confirms that file is
byte-compatible with the cNMF input:

```
Combined-clustered h5ad: 240504 cells x 25788 genes
cNMF input h5ad:         240504 cells x 25788 genes
obs_names identical after known barcode rename: True
var_names identical (order + content):          True
```

## History: the float / int detour (2026-04-08 → 04-09)

Step 2 was run three times. Reconstructed from directory mtimes and Slurm logs, using
sublibrary 1 as the tracer:

| When | What happened | sublib1 result |
|---|---|---|
| Apr 8 23:50 | Step 1 (float) → `1.umifiltering/output/` | float32 |
| Apr 9 08:51 | Step 2, first run, into `2.qcfiltering/output/` | 45,938 × 25,703 |
| Apr 9 09:57 | Step 3a concat reads it | **45,938 × 25,703** |
| Apr 9 10:29 | Step 4 clustering | 240,504 × 25,788 |
| Apr 9 13:00 | Step 2 re-run **truncating** to int, overwriting `output/` in place | 45,886 × 23,848 |
| Apr 9 21:26 | `output/` renamed to `trunc_int_output/`, `logs/` to `trunc_int_logs/` | — |
| Apr 9 21:33 | Step 2 **float** re-run into a fresh `float_output/` | **45,938 × 25,703** |

Truncating instead of rounding zeroed every sub-1.0 count, costing 52 cells and 1,855
genes. That run was abandoned; the float re-run reproduced the original 08:51 dimensions
exactly. So the 3a→4 branch (built 09:57, from the original `output/`) and the 3b/3b-i
branch (built later, from `float_output/`) rest on step-2 output that agrees on every
filter dimension. Dimensions were verified, not values byte-for-byte.

Rounding to int32 now happens only in 3b / 3b-i, and only at write time.

The `-old-` in `3a.concat-old-forclustering/` is a leftover from this episode and does
**not** mean the concat is stale.

### Two renames happened after these runs

Both left the scripts pointing at directories that no longer existed:

- `2.qcfiltering/output` → `trunc_int_output` (Apr 9, 21:26)
- `3.concat-old-forclustering` → `3a.concat-old-forclustering` (Apr 10–11, when the
  3b / 3b-i split was introduced)

## Changes made when importing into this repo

Paths and shell bugs were fixed so the scripts run as-is. **No analysis logic, parameter
or threshold was changed**, so re-running reproduces the same numbers.

- `1.umi_threshold_batch.sh` called `1.umi_threshold.py`, which does not exist — now
  calls `1.umi_threshold_float.py` (the variant that produced the data actually used).
- `1.umi_threshold_float.py` wrote to `20260408-postprocessing/` — now
  `1.umifiltering/output/`, where its outputs really are.
- `1.umi_threshold_int_layercheck.py` wrote to `1.umifiltering/output/new` — now
  `.../output/20260409-intcorrected-output`.
- `3a.concatenate_mudata.py` read `2.qcfiltering/output/` — now `float_output/`; output
  dir `3.concat-old-forclustering` → `3a.`.
- `4.cluster_and_annotate.py` read `3.concat-old-forclustering` → `3a.`.
- `3a.concat_batch.sh` and `4.cluster_annotate_batch.sh` are plain `#SBATCH` scripts but
  still ended with a stray `EOF` / `done` from the loop-submitter pattern, making every
  job exit nonzero even on success. Removed. `3a` also called `3.concatenate_mudata.py`
  (the file is `3a.`), and its log paths pointed at the pre-rename directory.
- `3b.converttocellranger_batch.sh` had job name `cellranger_featurecorr` and wrote both
  logs into `3b-i.cellranger_format_featurecorrected/logs/` while running the **3b**
  script — so 3b's logs landed in 3b-i's folder. Repointed to `3b.cellranger_format/logs/`.
- All batch scripts now resolve their Python script via `SCRIPT_DIR` (submit-time `pwd`
  for the loop submitters, `$SLURM_SUBMIT_DIR` for the two direct `sbatch` scripts, since
  Slurm copies those to a spool directory) instead of hardcoding the old
  `perturbpipeline/processingperturboutputs/maggie_scripts/` path.

## Known quirks, left as-is

These reflect how the data was actually produced. Left untouched deliberately.

- **`filteringparameters.md` disagrees with what ran.** The notes record a UMI ceiling of
  100,000 for every sublibrary; the batch script uses **75,000**, and the Slurm logs
  confirm 75,000 is what executed. The lower bounds match. Treat the script as truth.
- **`seurat_v3` on non-integer input.** Step 4 inherits float32 counts, and
  `highly_variable_genes(flavor="seurat_v3")` logged
  `expects raw count data, but non-integers were found`. Values are near-integer so HVG
  selection is approximately right, but it is not what `seurat_v3` assumes.
- **Non-deterministic cell order.** `1.umi_threshold_float.py` and `2.qc_metrics_filter.py`
  take shared barcodes via `list(set(...))`, so row order varies between runs. Both
  modalities are subset with the same list, so GEX and guide stay aligned and results are
  correct; 3b / 3b-i re-sort anyway. `1.umi_threshold_int_layercheck.py` uses `sorted()`.
- **Ribosomal threshold is aggressive.** `ribo_threshold = 4` (%), against 5 in the
  reference script this was adapted from and 12 in the winter 2025 iPSC-VIC screen. This
  is the single biggest lever on final cell count.
- **MAD outliers are computed pre-MT/ribo-filter** and applied after, so the MADs reflect
  the unfiltered population. Standard practice, but order matters if you change it.
- **`3b.converttocellranger.py` has no feature-ID uniqueness assert** (3b-i does). Stripping
  Ensembl version suffixes could in principle collide; checked sublib 1 and found 28,729
  features, all unique, so it did not bite here.
- **Dead fallback in step 4.** The cell-cycle fallback branch calls
  `sc.datasets.pbmc3k_processed()` — a network download on a compute node — and never uses
  the result. Never triggered, because the real gene-list files exist.
- `1.umi_threshold_int_layercheck.py` is a diagnostic variant (verifies
  `nascent + ambiguous + mature == X` and rounds to int32). Its output landed in
  `1.umifiltering/output/20260409-intcorrected-output/` and **nothing downstream consumed
  it**. Kept for provenance.
