process CHECK_ENV {
    tag "check_env"

    output:
    stdout

    script:
    """
    echo "=== PGIRL environment check ==="
    echo "Python: \$(${params.python} --version)"
    echo "Nextclade: \$(nextclade --version 2>/dev/null || echo 'nextclade not found')"
    echo "Kraken2: \$(kraken2 --version 2>/dev/null || echo 'kraken2 not found')"
    echo "Kraken2 DB: \$(test -d '${params.kraken_db}' && echo '${params.kraken_db} (found)' || echo '${params.kraken_db} (NOT found)')"
    echo "seqkit: \$(seqkit --version 2>/dev/null || echo 'seqkit not found')"
    echo "minimap2: \$(minimap2 --version 2>/dev/null || echo 'minimap2 not found')"
    ${params.python} - <<'PY'
    import importlib
    pkgs = [
        "psycopg2",
        "Bio",
        "yaml",
        "requests",
        "pandas",
        "pyarrow",
        "pydantic",
        "scipy",
        "networkx",
        "dendropy",
        "matplotlib",
    ]
    for pkg in pkgs:
        try:
            m = importlib.import_module(pkg)
            print(f"[OK] {pkg} {getattr(m, '__version__', '')}")
        except Exception as exc:
            print(f"[MISSING] {pkg}: {exc}")
    PY
    echo "=== check complete ==="
    """
}
