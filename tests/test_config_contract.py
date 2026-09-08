"""Checks that every shipped config and sample table satisfies the workflow's
contract. These run without the sequencing data, so CI can enforce them."""

import pandas as pd
import pytest

from scripts.pipeline_utils import (
    get_counts_h5ad,
    get_sublibraries,
    get_umi_bounds,
    load_sample_info,
    read_gene_list,
    resolve_config_paths,
)
from scripts.snakemake_helpers import (
    get_logs_path,
    get_resources,
    get_results_path,
    get_shell_prefix,
)

TOP_LEVEL_KEYS = (
    "analysis_name",
    "results_base",
    "sample_info_file",
    "input_paths",
    "qc_filtering",
    "concatenate",
    "cellranger_export",
    "clustering",
    "resources",
)

# Rules in workflows/00_counts_matrix_processing.smk that look up resources.
RESOURCE_RULES = (
    "umi_filter",
    "qc_filter",
    "concatenate",
    "cellranger_export",
    "master_gene_list",
    "clustering",
)

QC_KEYS = (
    "reference",
    "mt_threshold",
    "ribo_threshold",
    "n_mads_total_counts",
    "n_mads_n_genes",
    "filter_outliers",
    "filter_cells_by_min_genes",
    "min_genes_per_cell",
    "min_counts_per_cell",
    "min_cells_per_gene",
    "min_counts_per_gene",
)

CLUSTERING_KEYS = (
    "random_seed",
    "target_sum",
    "n_top_genes",
    "hvg_flavor",
    "hvg_batch_key",
    "n_pcs",
    "n_pcs_use",
    "n_neighbors",
    "leiden_resolutions",
    "dotplot_resolution",
    "regress_cell_cycle",
    "marker_panels",
)


@pytest.fixture
def config(raw_config, code_dir):
    """The config with repo-relative paths resolved, as the workflow sees it."""
    return resolve_config_paths(dict(raw_config), code_dir)


def test_top_level_keys_present(raw_config):
    missing = [key for key in TOP_LEVEL_KEYS if key not in raw_config]
    assert not missing, f"config is missing top-level keys: {missing}"


def test_qc_filtering_keys_present(raw_config):
    missing = [key for key in QC_KEYS if key not in raw_config["qc_filtering"]]
    assert not missing, f"qc_filtering is missing keys: {missing}"


def test_clustering_keys_present(raw_config):
    missing = [key for key in CLUSTERING_KEYS if key not in raw_config["clustering"]]
    assert not missing, f"clustering is missing keys: {missing}"


def test_qc_reference_is_supported(raw_config):
    assert raw_config["qc_filtering"]["reference"] in ("human", "mouse")


def test_dotplot_resolution_is_a_clustered_resolution(raw_config):
    params = raw_config["clustering"]
    assert params["dotplot_resolution"] in params["leiden_resolutions"], (
        "clustering.dotplot_resolution must be one of clustering.leiden_resolutions"
    )


def test_n_pcs_use_does_not_exceed_n_pcs(raw_config):
    params = raw_config["clustering"]
    assert params["n_pcs_use"] <= params["n_pcs"]


def test_marker_panels_are_non_empty_symbol_lists(raw_config):
    panels = raw_config["clustering"]["marker_panels"]
    assert panels, "clustering.marker_panels is empty"
    for name, genes in panels.items():
        assert genes, f"marker panel {name} is empty"
        assert len(set(genes)) == len(genes), f"marker panel {name} has duplicates"


def test_bundled_reference_files_exist(config):
    for key in ("cell_cycle_s_genes", "cell_cycle_g2m_genes"):
        path = config["input_paths"][key]
        assert path, f"input_paths.{key} is unset"
        assert path.endswith(".txt")
        assert read_gene_list(path), f"{path} is empty"


def test_cell_cycle_gene_lists_are_unique(config):
    s_genes = read_gene_list(config["input_paths"]["cell_cycle_s_genes"])
    g2m_genes = read_gene_list(config["input_paths"]["cell_cycle_g2m_genes"])
    assert len(set(s_genes)) == len(s_genes)
    assert len(set(g2m_genes)) == len(g2m_genes)
    assert not set(s_genes) & set(g2m_genes), "a gene appears in both phase lists"


def test_resources_cover_every_rule(raw_config):
    for rule_name in RESOURCE_RULES:
        for key in ("mem_mb", "runtime", "threads", "slurm_partition"):
            value = get_resources(raw_config, rule_name, key)
            assert value, f"resources.{rule_name}.{key} resolved to {value!r}"


def test_resources_fall_back_to_default(raw_config):
    assert get_resources(raw_config, "a_rule_with_no_entry", "mem_mb") == (
        raw_config["resources"]["default"]["mem_mb"]
    )


def test_missing_resource_raises(raw_config):
    with pytest.raises(KeyError):
        get_resources(raw_config, "umi_filter", "not_a_resource")


def test_path_helpers_are_rooted_in_results_base(raw_config):
    results = get_results_path(config=raw_config)
    assert get_results_path("clustering", config=raw_config).startswith(results)
    assert get_logs_path("umi_filter", config=raw_config).startswith(results)


def test_shell_prefix_exports_pythonpath(raw_config, code_dir):
    prefix = get_shell_prefix(raw_config, code_dir)
    assert str(code_dir) in prefix
    assert "export PYTHONPATH=" in prefix
    assert prefix.endswith("; ")
    if raw_config.get("conda", {}).get("env"):
        assert "conda activate" in prefix


def test_shell_prefix_braces_survive_snakemake_formatting(raw_config, code_dir):
    """Snakemake formats the prefix, so doubled braces must collapse to one."""
    prefix = get_shell_prefix(raw_config, code_dir)
    formatted = prefix.format()
    assert "${PYTHONPATH:+:$PYTHONPATH}" in formatted
    assert "{{" not in formatted and "}}" not in formatted


def test_shell_prefix_without_conda_env(raw_config, code_dir):
    config = dict(raw_config)
    config["conda"] = {}
    prefix = get_shell_prefix(config, code_dir)
    assert "conda activate" not in prefix
    assert "PYTHONPATH" in prefix


# --- sample table --------------------------------------------------------


def test_sample_info_file_loads(config):
    df = load_sample_info(config)
    assert len(df) > 0


def test_every_sublibrary_has_one_gex_and_one_guide(config):
    df = load_sample_info(config)
    counts = df.groupby(["sublibrary", "sample_type"]).size().unstack(fill_value=0)
    for sample_type in ("gex", "guide"):
        assert sample_type in counts.columns, f"no {sample_type} rows in sample table"
        offenders = counts.index[counts[sample_type] != 1].tolist()
        assert not offenders, (
            f"sublibraries without exactly one {sample_type} row: {offenders}"
        )


def test_sample_ids_are_unique(config):
    df = load_sample_info(config)
    assert df["sample_id"].is_unique, "sample_id values must be unique"


def test_sample_types_are_recognised(config):
    df = load_sample_info(config)
    unexpected = set(df["sample_type"]) - {"gex", "guide"}
    assert not unexpected, f"unrecognised sample_type values: {unexpected}"


def test_get_sublibraries_matches_gex_rows(config):
    df = load_sample_info(config)
    expected = list(
        dict.fromkeys(df[df["sample_type"] == "gex"]["sublibrary"].astype(str))
    )
    assert get_sublibraries(config) == expected


def test_umi_bounds_are_sane(config):
    for sublibrary in get_sublibraries(config):
        umi_min, umi_max = get_umi_bounds(config, sublibrary)
        assert umi_min > 0, f"{sublibrary}: min_umi_threshold must be positive"
        if umi_max is not None:
            assert umi_max > umi_min, (
                f"{sublibrary}: max_umi_threshold ({umi_max}) must exceed "
                f"min_umi_threshold ({umi_min})"
            )


def test_counts_paths_are_absolute_and_distinct(config):
    paths = []
    for sublibrary in get_sublibraries(config):
        for sample_type in ("gex", "guide"):
            path = get_counts_h5ad(config, sublibrary, sample_type)
            assert path.startswith("/"), f"{path} is not absolute"
            assert path.endswith("adata.h5ad")
            paths.append(path)
    assert len(set(paths)) == len(paths), "two samples map to the same adata.h5ad"


def test_excluded_rows_are_dropped(config, tmp_path):
    df = load_sample_info(config)
    df.loc[df["sample_type"] == "gex", "exclude"] = True
    path = tmp_path / "sample_info.tsv"
    df.to_csv(path, sep="\t", index=False)
    assert (load_sample_info(str(path))["sample_type"] == "guide").all()


def test_missing_required_column_raises(config, tmp_path):
    df = load_sample_info(config).drop(columns=["counts_dir"])
    path = tmp_path / "sample_info.tsv"
    df.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="counts_dir"):
        load_sample_info(str(path))


def test_min_umi_threshold_present_for_every_gex_row(config):
    df = load_sample_info(config)
    gex = df[df["sample_type"] == "gex"]
    assert not gex["min_umi_threshold"].apply(pd.isna).any(), (
        "every gex row needs a min_umi_threshold"
    )
