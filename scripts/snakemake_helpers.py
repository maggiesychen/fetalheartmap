"""Path and resource helpers for the Snakemake workflows."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipeline_utils import (  # noqa: E402
    get_counts_h5ad,
    get_sublibraries,
    get_umi_bounds,
    load_sample_info,
)

__all__ = [
    "get_counts_h5ad",
    "get_logs_path",
    "get_resources",
    "get_results_path",
    "get_scratch_path",
    "get_shell_prefix",
    "get_sublibraries",
    "get_umi_bounds",
    "load_sample_info",
    "print_path_configuration",
]


def _base(config, key):
    if not config:
        raise ValueError("config is required")
    if config.get(key):
        return config[key]
    if key == "scratch_base":
        return os.environ["SCRATCH"]
    raise KeyError(f"config is missing required key: {key}")


def get_scratch_path(*parts, config=None):
    """Build a path below the configured scratch directory."""
    scratch_dir = os.path.join(_base(config, "scratch_base"), config["analysis_name"])
    if parts and parts[0]:
        return os.path.join(scratch_dir, *[str(p) for p in parts])
    return scratch_dir


def get_results_path(*parts, config=None):
    """Build a path below the analysis results directory."""
    results_base = _base(config, "results_base")
    if parts and parts[0]:
        return os.path.join(results_base, *[str(p) for p in parts])
    return results_base


def get_logs_path(*parts, config=None):
    """Build a path below the results log directory."""
    logs_base = get_results_path("logs", config=config)
    if parts and parts[0]:
        return os.path.join(logs_base, *[str(p) for p in parts])
    return logs_base


def get_resources(config, rule_name, key):
    """Look up one Slurm resource for a rule, falling back to ``default``."""
    resources = config.get("resources", {})
    per_rule = resources.get(rule_name, {})
    if key in per_rule:
        return per_rule[key]
    default = resources.get("default", {})
    if key in default:
        return default[key]
    raise KeyError(
        f"resources.{rule_name}.{key} is not set and resources.default.{key} "
        "does not exist"
    )


def get_shell_prefix(config, code_dir):
    """Shell prefix that activates the analysis env and exports PYTHONPATH.

    Snakemake runs ``shell.prefix`` through ``str.format``, so any literal
    brace destined for the shell has to be doubled here.
    """
    parts = ["set -euo pipefail"]

    conda = config.get("conda", {})
    profile_script = conda.get("profile_script")
    env = conda.get("env")
    if env:
        if profile_script:
            # conda activate needs the hook sourced in a non-interactive shell.
            parts.append(f"set +u; source {profile_script}; set -u")
        parts.append(f"conda activate {env}")

    # ${{PYTHONPATH:+:$PYTHONPATH}} -> ${PYTHONPATH:+:$PYTHONPATH} after
    # Snakemake formats the prefix: append only if PYTHONPATH is already set.
    parts.append(
        f'export PYTHONPATH="{code_dir}${{{{PYTHONPATH:+:$PYTHONPATH}}}}"'
    )
    return "; ".join(parts) + "; "


def print_path_configuration(config):
    """Print the resolved paths and sublibraries at workflow startup."""
    print("=" * 72)
    print(f"analysis_name : {config['analysis_name']}")
    print(f"results       : {get_results_path(config=config)}")
    print(f"scratch       : {get_scratch_path(config=config)}")
    print(f"logs          : {get_logs_path(config=config)}")
    print(f"sample info   : {config['sample_info_file']}")
    print(f"kb output     : {config['input_paths']['kb_output_base']}")
    sublibraries = get_sublibraries(config)
    print(f"sublibraries  : {len(sublibraries)} -> {', '.join(sublibraries)}")
    print("=" * 72)
