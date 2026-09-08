"""Shared fixtures. Adds the repo root to ``sys.path`` so ``scripts`` imports."""

import sys
from pathlib import Path

import pytest
import yaml

CODE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(CODE_DIR))

CONFIG_FILES = sorted((CODE_DIR / "config").glob("config.*.yaml"))


@pytest.fixture(scope="session")
def code_dir():
    return CODE_DIR


@pytest.fixture(params=[p.name for p in CONFIG_FILES], scope="session")
def config_path(request):
    """Every config file shipped in config/, one test invocation each."""
    return CODE_DIR / "config" / request.param


@pytest.fixture(scope="session")
def raw_config(config_path):
    with open(config_path) as handle:
        return yaml.safe_load(handle)
