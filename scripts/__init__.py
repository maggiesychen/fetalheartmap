"""Step scripts and shared helpers for the fetal heart map workflows.

Each step script is a standalone, argparse-driven entry point invoked by a
Snakemake rule; ``pipeline_utils`` and ``snakemake_helpers`` hold the code
shared between them. The repo root must be on ``PYTHONPATH`` for the
``scripts.*`` imports to resolve -- the workflows set this via
``shell.prefix``.
"""
