"""Reusable bioinformatics modules for the PGIRL pipeline.

Each module corresponds to one stage in the analysis plan. The goal is to
wrap validated external tools (Nextclade, MAFFT, IQ-TREE, etc.) in small,
pathogen-agnostic Python functions that can later be turned into Nextflow
processes.
"""
