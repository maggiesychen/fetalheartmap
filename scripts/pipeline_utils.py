"""Config, sample-table, and AnnData helpers shared by the workflow scripts."""

import os

import numpy as np
import pandas as pd
import yaml

CODE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Paths that ship with this repo and are therefore resolved against CODE_DIR.
_BUNDLED_TOP_LEVEL = ("sample_info_file",)
_BUNDLED_INPUT_PATHS = (
    "cell_cycle_s_genes",
    "cell_cycle_g2m_genes",
)


def _resolve_bundle_path(value, bundle_dir):
    """Make a repo-relative path absolute; leave absolute paths alone."""
    if not value or os.path.isabs(value):
        return value
    return os.path.join(bundle_dir, value)


def resolve_config_paths(config, bundle_dir=CODE_DIR):
    """Resolve paths for files shipped with this code bundle.

    Runtime paths such as ``scratch_base`` and ``results_base`` stay relative to
    the directory Snakemake was launched from, unless the user made them
    absolute.
    """
    for key in _BUNDLED_TOP_LEVEL:
        if config.get(key):
            config[key] = _resolve_bundle_path(config[key], bundle_dir)

    input_paths = config.get("input_paths", {})
    for key in _BUNDLED_INPUT_PATHS:
        if input_paths.get(key):
            input_paths[key] = _resolve_bundle_path(input_paths[key], bundle_dir)

    return config


def load_config(config_file):
    """Load a YAML config and resolve its bundled input paths."""
    with open(config_file) as handle:
        config = yaml.safe_load(handle)
    return resolve_config_paths(config)


def read_sample_info_file(sample_info_file):
    """Read the sample table, accepting .tsv/.csv/.xlsx."""
    ext = os.path.splitext(sample_info_file)[1].lower()
    if ext in (".xlsx", ".xls"):
        return pd.read_excel(sample_info_file)
    sep = "," if ext == ".csv" else "\t"
    return pd.read_csv(sample_info_file, sep=sep)


def load_sample_info(config_or_path):
    """Return the sample table with excluded rows dropped.

    Accepts either a config dict or a path to the sample table.
    """
    if isinstance(config_or_path, dict):
        sample_info_file = config_or_path["sample_info_file"]
    else:
        sample_info_file = config_or_path

    df = read_sample_info_file(sample_info_file)

    required = {"sublibrary", "sample_id", "sample_type", "counts_dir"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"{sample_info_file} is missing required columns: {sorted(missing)}"
        )

    if "exclude" in df.columns:
        excluded = df["exclude"].astype(str).str.lower().isin(("true", "1", "yes"))
        df = df[~excluded]

    return df.reset_index(drop=True)


def get_sublibraries(config):
    """Sublibrary labels that have a GEX library, in table order."""
    df = load_sample_info(config)
    gex = df[df["sample_type"] == "gex"]
    return list(dict.fromkeys(gex["sublibrary"].astype(str)))


def get_sample_row(config, sublibrary, sample_type):
    """The single sample-table row for one sublibrary and modality."""
    df = load_sample_info(config)
    hits = df[
        (df["sublibrary"].astype(str) == str(sublibrary))
        & (df["sample_type"] == sample_type)
    ]
    if len(hits) != 1:
        raise ValueError(
            f"Expected exactly one {sample_type} row for sublibrary "
            f"{sublibrary}, found {len(hits)}"
        )
    return hits.iloc[0]


def get_counts_h5ad(config, sublibrary, sample_type):
    """Absolute path to one kb-python ``adata.h5ad`` for a sublibrary."""
    row = get_sample_row(config, sublibrary, sample_type)
    subdir = row.get("counts_subdir")
    parts = [config["input_paths"]["kb_output_base"], str(row["counts_dir"])]
    if isinstance(subdir, str) and subdir.strip():
        parts.append(subdir.strip())
    parts.append("adata.h5ad")
    return os.path.join(*parts)


def get_umi_bounds(config, sublibrary):
    """``(min_umi, max_umi)`` for a sublibrary; ``max_umi`` may be ``None``."""
    row = get_sample_row(config, sublibrary, "gex")
    umi_min = row.get("min_umi_threshold")
    umi_max = row.get("max_umi_threshold")
    if pd.isna(umi_min):
        raise ValueError(f"min_umi_threshold is missing for sublibrary {sublibrary}")
    umi_max = None if pd.isna(umi_max) else int(umi_max)
    return int(umi_min), umi_max


def mito_ribo_prefixes(reference):
    """MT and ribosomal gene-symbol prefixes for a species."""
    if reference == "human":
        return "MT-", ("RPS", "RPL")
    if reference == "mouse":
        return "Mt-", ("Rps", "Rpl")
    raise ValueError(f"Unsupported reference: {reference!r} (use human or mouse)")


def load_symbol_map(npy_path):
    """Load a pickled ``{ensembl_id: symbol}`` dict, or ``None`` if unavailable."""
    if not npy_path or not os.path.exists(npy_path):
        return None
    return np.load(npy_path, allow_pickle=True).item()


def annotate_mito_ribo(adata, reference, ensembl_to_symbol=None):
    """Add ``var['symbol']``, ``var['mt']`` and ``var['ribo']`` to ``adata``.

    Resolves gene symbols from ``var_names``, mapping through
    ``ensembl_to_symbol`` when ``var_names`` are Ensembl IDs.
    """
    mt_prefix, ribo_prefix = mito_ribo_prefixes(reference)

    var_names_are_ensembl = str(adata.var_names[0]).startswith(("ENSG", "ENSMUSG"))

    if var_names_are_ensembl:
        symbol_map = load_symbol_map(ensembl_to_symbol)
        if symbol_map is not None:
            symbols = (
                pd.Series(adata.var_names).str.split(".").str[0].map(symbol_map)
            )
        else:
            print(
                f"Warning: ensembl->symbol dict not found at {ensembl_to_symbol}. "
                "MT/ribo annotation may be incomplete."
            )
            symbols = pd.Series(adata.var_names)
    else:
        if "gene_id" in adata.var.columns:
            adata.var.reset_index(inplace=True)
            adata.var.index = adata.var["gene_id"]
        symbols = pd.Series(adata.var_names, index=adata.var_names)

    adata.var["symbol"] = symbols.values
    adata.var["mt"] = adata.var["symbol"].str.startswith(mt_prefix).fillna(False)
    adata.var["ribo"] = adata.var["symbol"].str.startswith(ribo_prefix).fillna(False)
    print(
        f"MT genes: {adata.var['mt'].sum()}  |  "
        f"Ribo genes: {adata.var['ribo'].sum()}"
    )
    return adata


def read_gene_list(path):
    """Read a one-gene-per-line text file, dropping blanks and comments."""
    with open(path) as handle:
        return [
            line.strip()
            for line in handle
            if line.strip() and not line.startswith("#")
        ]
