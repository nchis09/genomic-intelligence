process CHECK_ENV {
    tag "check_env"

    output:
    stdout

    script:
    """
    echo "=== PGIRL environment check ==="
    echo "Python: \$(python3 --version)"
    echo "Nextclade: \$(nextclade --version 2>/dev/null || echo 'nextclade not found')"
    python3 - <<'PY'
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
