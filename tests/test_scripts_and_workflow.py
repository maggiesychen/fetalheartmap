"""Static checks on the scripts and the Snakemake workflow.

These do not touch the sequencing data. They catch the failure modes that
otherwise only surface hours into a cluster run: a syntax error in a script, a
rule whose script path no longer exists, or a hardcoded absolute path creeping
back into the code.
"""

import ast
import py_compile
import re

import pytest

SCRIPT_DIR_NAME = "scripts"
WORKFLOW_DIR_NAME = "workflows"

# Absolute paths belong in config/ and references/, not in code. Any of these
# roots appearing in a .py or .smk file is a regression.
FORBIDDEN_PATH_PATTERN = re.compile(r"[\"'](/oak/|/scratch/|/home/|/share/)")

# Library modules, not argparse-driven step scripts.
HELPER_MODULES = {
    "pipeline_utils.py",
    "snakemake_helpers.py",
    "figstyle.py",
    "__init__.py",
    "r_utils.R",
}


def python_scripts(code_dir):
    return sorted((code_dir / SCRIPT_DIR_NAME).glob("*.py"))


def r_scripts(code_dir):
    return sorted((code_dir / SCRIPT_DIR_NAME).glob("*.R"))


def all_scripts(code_dir):
    return python_scripts(code_dir) + r_scripts(code_dir)


def workflow_files(code_dir):
    return sorted((code_dir / WORKFLOW_DIR_NAME).glob("*.smk"))


def test_scripts_directory_is_not_empty(code_dir):
    assert python_scripts(code_dir), "no scripts found"


def test_workflows_directory_is_not_empty(code_dir):
    assert workflow_files(code_dir), "no .smk workflows found"


def test_every_script_compiles(code_dir, tmp_path):
    for path in python_scripts(code_dir):
        # Write the bytecode into tmp_path so the repo stays clean.
        py_compile.compile(
            str(path), doraise=True, cfile=str(tmp_path / (path.stem + ".pyc"))
        )


def test_every_script_has_a_module_docstring(code_dir):
    missing = [
        path.name
        for path in python_scripts(code_dir)
        if not ast.get_docstring(ast.parse(path.read_text()))
    ]
    assert not missing, f"scripts without a module docstring: {missing}"


def test_step_scripts_are_argparse_driven(code_dir):
    """Every step script must take its inputs and outputs on the command line."""
    for path in python_scripts(code_dir):
        if path.name in HELPER_MODULES:
            continue
        source = path.read_text()
        assert "argparse.ArgumentParser" in source, f"{path.name} has no argparse"
        assert "--output" in source, f"{path.name} has no --output argument"


def test_r_scripts_are_cli_driven(code_dir):
    """R step scripts take their inputs on the command line, via r_utils.R.

    Deliberately not optparse: the shared R library has neither optparse nor
    getopt, so `scripts/r_utils.R` provides a base-R parser instead.
    """
    for path in r_scripts(code_dir):
        if path.name in HELPER_MODULES:
            continue
        source = path.read_text()
        assert "parse_cli_args" in source, f"{path.name} does not use parse_cli_args"
        assert 'source(file.path(this_script_dir_bootstrap(), "r_utils.R"))' in source, (
            f"{path.name} does not source r_utils.R"
        )
        assert "output" in source, f"{path.name} declares no output option"


def test_r_helper_has_no_package_dependencies(code_dir):
    """r_utils.R must stay dependency-free -- that is its whole purpose."""
    source = (code_dir / "scripts" / "r_utils.R").read_text()
    for line in source.splitlines():
        stripped = line.strip()
        if stripped.startswith("#"):
            continue
        assert not stripped.startswith(("library(", "require(")), (
            f"r_utils.R must not load packages: {stripped}"
        )


def test_r_scripts_have_a_header_comment(code_dir):
    """Every R script opens with a `##` block explaining what it does."""
    for path in r_scripts(code_dir):
        lines = [
            line for line in path.read_text().splitlines()
            if line.strip() and not line.startswith("#!")
        ]
        assert lines and lines[0].startswith("#"), (
            f"{path.name} has no header comment"
        )


def test_no_hardcoded_absolute_paths_in_code(code_dir):
    offenders = []
    for path in all_scripts(code_dir) + workflow_files(code_dir):
        for lineno, line in enumerate(path.read_text().splitlines(), start=1):
            stripped = line.strip()
            if stripped.startswith("#"):
                continue  # comments may cite real paths for provenance
            if FORBIDDEN_PATH_PATTERN.search(line):
                offenders.append(f"{path.name}:{lineno}: {stripped}")
    assert not offenders, (
        "absolute paths belong in config/, not in code:\n" + "\n".join(offenders)
    )


def test_workflow_script_references_exist(code_dir):
    """Every `CODE_DIR / "scripts/..."` in a workflow must resolve to a file."""
    pattern = re.compile(r'CODE_DIR\s*/\s*"([^"]+)"')
    referenced = set()
    for path in workflow_files(code_dir):
        referenced.update(pattern.findall(path.read_text()))
    assert referenced, "no script references found in the workflows"
    missing = [rel for rel in referenced if not (code_dir / rel).is_file()]
    assert not missing, f"workflow references missing files: {missing}"


def test_every_script_is_referenced_by_a_workflow(code_dir):
    """Guards against scripts that quietly fall out of the DAG."""
    helpers = HELPER_MODULES
    workflow_text = "\n".join(p.read_text() for p in workflow_files(code_dir))
    orphans = [
        path.name
        for path in all_scripts(code_dir)
        if path.name not in helpers and path.name not in workflow_text
    ]
    assert not orphans, f"scripts not referenced by any workflow: {orphans}"


def test_workflow_rules_declare_log_and_resources(code_dir):
    """Every rule that runs a shell command must log and declare resources."""
    for path in workflow_files(code_dir):
        blocks = re.split(r"^rule ", path.read_text(), flags=re.MULTILINE)[1:]
        for block in blocks:
            rule_name = block.split(":", 1)[0].strip()
            if "shell:" not in block:
                continue  # aggregation-only targets such as `all`
            assert "log:" in block, f"{path.name}: rule {rule_name} has no log"
            assert "resources:" in block, (
                f"{path.name}: rule {rule_name} declares no resources"
            )


def test_environment_files_pin_versions(code_dir):
    env_files = [code_dir / "environment.yml"] + sorted(
        (code_dir / "envs").glob("*.yml")
    )
    for path in env_files:
        assert path.is_file(), f"{path} is missing"
        pins = [
            line
            for line in path.read_text().splitlines()
            if re.match(r"\s*-\s*[A-Za-z].*[=]{1,2}\d", line)
        ]
        assert pins, f"{path.name} pins no package versions"


@pytest.mark.parametrize("required", ["README.md", "submit.sh", ".gitignore"])
def test_repo_scaffolding_present(code_dir, required):
    assert (code_dir / required).is_file(), f"{required} is missing"
