"""Pathogen-specific bioinformatics workflows.

Each module implements a `run()` function that accepts the sample inputs and
returns paths to the per-sample outputs. The top-level pipeline dispatches to
the appropriate workflow based on the identified pathogen.
"""
