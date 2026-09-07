#!/usr/bin/env python
# %%
# Based on: https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html
import os

import matplotlib
import matplotlib.pyplot as plt
import muon as mu
import numpy as np
import scanpy as sc

matplotlib.rcParams["axes.spines.top"] = False
matplotlib.rcParams["axes.spines.right"] = False
sc.settings.set_figure_params(dpi=100, facecolor="white", frameon=False)

# %%
# --- Parameters ---
input_path = "/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/3a.concat-old-forclustering/300k_combined_all_sublibs_newguidelist.h5mu"
output_path = "/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/4.clustering"
output_name = "300k_combined_clustered"

ens2symbol_path = "/oak/stanford/groups/engreitz/Users/opushkar/genome/ensembl_to_symbol_dict_v43.npy"
symbol2ens_path = "/oak/stanford/groups/engreitz/Users/opushkar/genome/symbol_to_ensembl_dict_v43.npy"

# Cell cycle gene list directory (Regev lab format)
cc_genes_path = "/oak/stanford/groups/engreitz/Users/opushkar/common_sc"

# Normalization
target_sum = 1e4

# HVG selection
n_top_genes = 3000
hvg_batch_key = "sublib"   # set to None to ignore batch when selecting HVGs

# PCA / neighbors / UMAP
n_pcs = 50
n_pcs_use = 30   # PCs to use for neighbor graph (adjust after inspecting elbow plot)
n_neighbors = 15

# Leiden clustering resolutions to try
leiden_resolutions = [0.1, 0.25, 0.5, 1.0]

# %%
# --- Output directories ---
plots_dir = os.path.join(output_path, "plots", "clustering")
os.makedirs(plots_dir, exist_ok=True)
os.makedirs(os.path.join(output_path, "clustered_data"), exist_ok=True)


# %%
# --- Helper: save current figure ---
def savefig(name):
    path = os.path.join(plots_dir, name)
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {path}")


# %%
# =============================================================================
# 1. LOAD DATA
# =============================================================================
print("Loading MuData...")
mdata = mu.read(input_path)
adata = mdata["GEX"].copy()
adata.X = adata.X.astype("float32")
print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

# %%
# =============================================================================
# 2. GENE SYMBOL MAPPING
#    Resolve var_names to gene symbols so marker genes can be looked up by name.
# =============================================================================
var_names_are_ensembl = adata.var_names[0].startswith(("ENSG", "ENSMUSG"))

if var_names_are_ensembl and os.path.exists(ens2symbol_path):
    ens2symbol = np.load(ens2symbol_path, allow_pickle=True).item()
    symbol2ens = np.load(symbol2ens_path, allow_pickle=True).item()
    # store original ensembl IDs then switch var_names to symbols for easy lookup
    adata.var["gene_id"] = adata.var_names.str.split(".").str[0]
    adata.var["symbol"] = adata.var["gene_id"].map(ens2symbol)
    # drop genes that couldn't be mapped (NaN symbols)
    adata = adata[:, adata.var["symbol"].notna()].copy()
    adata.var_names = adata.var["symbol"].astype(str)
    adata.var_names_make_unique()
    print(f"After symbol mapping: {adata.n_obs} cells x {adata.n_vars} genes")
else:
    print("var_names appear to already be gene symbols — skipping mapping.")

# %%
# =============================================================================
# 3. QC METRICS  (re-calculated on the combined object)
# =============================================================================
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo"], inplace=True, log1p=True
)

# %%
# =============================================================================
# 4. NORMALIZATION  (CP10k + log1p)
# =============================================================================
adata.layers["counts"] = adata.X.copy()   # preserve raw integer counts

sc.pp.normalize_total(adata, target_sum=target_sum)
sc.pp.log1p(adata)
print("Normalization complete.")

# %%
# =============================================================================
# 5. CELL CYCLE SCORING
#    Requires s_genes.txt and g2m_genes.txt in cc_genes_path,
#    or falls back to scanpy's built-in PBMC gene list.
# =============================================================================
s_genes_file = os.path.join(cc_genes_path, "s_genes.txt")
g2m_genes_file = os.path.join(cc_genes_path, "g2m_genes.txt")

if os.path.exists(s_genes_file) and os.path.exists(g2m_genes_file):
    s_genes = [g.strip() for g in open(s_genes_file)]
    g2m_genes = [g.strip() for g in open(g2m_genes_file)]
else:
    # fallback: standard Tirosh et al. gene sets shipped with scanpy
    cell_cycle_genes = [
        x.strip() for x in sc.datasets.pbmc3k_processed().var_names
    ]
    # Use commonly cited S/G2M gene lists
    s_genes = [
        "MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6",
        "CDCA7","DTL","PRIM1","UHRF1","HELLS","RFC2","RPA2","NASP","RAD51AP1",
        "GMNN","WDR76","SLBP","CCNE2","UBR7","POLD3","MSH2","ATAD2","RAD51",
        "RRM2","CDC45","CDC6","EXO1","TIPIN","DSCC1","BLM","CASP8AP2","USP1",
        "CLSPN","POLA1","CHAF1B","BRIP1","E2F8",
    ]
    g2m_genes = [
        "HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2",
        "NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","FAM64A","SMC4","CCNB2",
        "CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1","KIF20B",
        "HJURP","CDCA3","HN1","CDC20","TTK","CDC25C","KIF2C","RANGAP1","NCAPD2",
        "DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR","AURKA","PSRC1","ANLN",
        "LBR","CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA",
    ]

s_genes = [g for g in s_genes if g in adata.var_names]
g2m_genes = [g for g in g2m_genes if g in adata.var_names]

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
print(f"Cell cycle scoring done. S genes found: {len(s_genes)}, G2M genes found: {len(g2m_genes)}")

# %%
# =============================================================================
# 6. HIGHLY VARIABLE GENES
# =============================================================================
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=n_top_genes,
    batch_key=hvg_batch_key,
    flavor="seurat_v3",
    layer="counts",
)
print(f"HVGs selected: {adata.var['highly_variable'].sum()}")

fig = sc.pl.highly_variable_genes(adata, show=False)
savefig("hvg_plot.png")

# %%
# =============================================================================
# 7. PCA
# =============================================================================
sc.tl.pca(adata, n_comps=n_pcs, use_highly_variable=True)

fig = sc.pl.pca_variance_ratio(adata, n_pcs=n_pcs, log=True, show=False)
savefig("pca_variance_ratio.png")

fig = sc.pl.pca(
    adata,
    color=["sublib", "pct_counts_mt", "total_counts", "phase"],
    ncols=2,
    size=3,
    show=False,
)
savefig("pca_qc_overview.png")

# %%
# =============================================================================
# 8. NEAREST NEIGHBOR GRAPH + UMAP
# =============================================================================
sc.pp.neighbors(adata, n_pcs=n_pcs_use, n_neighbors=n_neighbors)
sc.tl.umap(adata)
print("UMAP computed.")

# %%
# =============================================================================
# 9. CLUSTERING  (Leiden at multiple resolutions)
# =============================================================================
for res in leiden_resolutions:
    key = f"leiden_res{res:.2f}".replace(".", "_")
    sc.tl.leiden(adata, resolution=res, key_added=key, flavor="igraph", n_iterations=2)
    print(f"Leiden res={res}: {adata.obs[key].nunique()} clusters")

# Plot all resolutions side-by-side
fig = sc.pl.umap(
    adata,
    color=[f"leiden_res{r:.2f}".replace(".", "_") for r in leiden_resolutions],
    legend_loc="on data",
    ncols=2,
    show=False,
)
savefig("umap_leiden_all_resolutions.png")

# %%
# =============================================================================
# 10. UMAP — QC METRICS
# =============================================================================
fig = sc.pl.umap(
    adata,
    color=[
        "sublib",
        "total_counts",
        "log1p_total_counts",
        "n_genes_by_counts",
        "pct_counts_mt",
        "pct_counts_ribo",
    ],
    ncols=3,
    size=3,
    show=False,
)
savefig("umap_qc_metrics.png")

# %%
# =============================================================================
# 11. UMAP — CELL CYCLE
# =============================================================================
fig = sc.pl.umap(
    adata,
    color=["phase", "S_score", "G2M_score"],
    ncols=3,
    size=3,
    show=False,
)
savefig("umap_cell_cycle.png")

# %%
# =============================================================================
# 12. MARKER GENE SCORING + UMAP
#     sc.tl.score_genes computes an aggregate score per cell for each panel.
# =============================================================================
marker_panels = {
    "ipsc": [
        "SOX2", "NANOG", "POU5F1", "KLF4", "MYC", "LIN28A",
    ],
    "vic": [
        "VIM", "COL1A1", "COL3A1", "DCN", "LUM", "POSTN",
        "COL11A1", "PDGFRA", "ACTA2", "TAGLN",
    ],
    "ec": [
        "CDH5", "PECAM1", "CD34", "KDR", "ANGPT2",
        "SOX17", "CXCR4", "DLL4", "EFNB2",   # arterial
        "NR2F2", "NT5E", "FLRT2",              # venous
    ],
}

for panel_name, genes in marker_panels.items():
    present = [g for g in genes if g in adata.var_names]
    missing = [g for g in genes if g not in adata.var_names]
    if missing:
        print(f"[{panel_name}] genes not found in data: {missing}")
    if not present:
        print(f"[{panel_name}] No genes found — skipping score.")
        continue
    sc.tl.score_genes(adata, gene_list=present, score_name=f"{panel_name}_score")
    print(f"[{panel_name}] Scored {len(present)} genes.")

# Plot scores on UMAP
score_cols = [f"{p}_score" for p in marker_panels if f"{p}_score" in adata.obs.columns]
fig = sc.pl.umap(
    adata,
    color=score_cols,
    cmap="bwr",
    vcenter=0,
    ncols=3,
    size=3,
    show=False,
)
savefig("umap_marker_scores.png")

# Also plot individual marker genes directly
for panel_name, genes in marker_panels.items():
    present = [g for g in genes if g in adata.var_names]
    if not present:
        continue
    fig = sc.pl.umap(
        adata,
        color=present,
        ncols=4,
        size=3,
        show=False,
        vmin=0,
    )
    savefig(f"umap_{panel_name}_individual_genes.png")

# %%
# =============================================================================
# 13. DOTPLOT — marker genes by cluster
# =============================================================================
leiden_key = f"leiden_res{leiden_resolutions[1]:.2f}".replace(".", "_")  # use second-lowest res

for panel_name, genes in marker_panels.items():
    present = [g for g in genes if g in adata.var_names]
    if not present:
        continue
    fig = sc.pl.dotplot(
        adata,
        var_names=present,
        groupby=leiden_key,
        standard_scale="var",
        show=False,
    )
    savefig(f"dotplot_{panel_name}_by_{leiden_key}.png")

# %%
# =============================================================================
# 15. SAVE
# =============================================================================
out_file = os.path.join(output_path, "clustered_data", f"{output_name}.h5ad")
adata.var.index.name = "gene_symbol"
adata.write_h5ad(out_file)
print(f"\nSaved clustered AnnData to:\n  {out_file}")
print(f"Final: {adata.n_obs} cells x {adata.n_vars} genes")
# %%
