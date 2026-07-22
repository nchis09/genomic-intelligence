process SETUP_ENV {
    tag "setup_env"

    output:
    stdout

    script:
    """
    echo "=== PGIRL environment setup ==="
    ${params.python} - <<'PY'
    import importlib
    import subprocess
    import sys

    # Map package name (for pip) -> import name
    packages = {
        "psycopg2-binary": "psycopg2",
        "biopython":       "Bio",
        "pyyaml":          "yaml",
        "requests":        "requests",
        "pandas":          "pandas",
        "pyarrow":         "pyarrow",
        "pydantic>=2":     "pydantic",
        "scipy":           "scipy",
        "networkx":        "networkx",
        "dendropy":        "dendropy",
        "matplotlib":      "matplotlib",
        "typing_extensions": "typing_extensions",
    }

    missing = []
    for pkg_name, import_name in packages.items():
        try:
            importlib.import_module(import_name)
            print(f"[OK] {pkg_name}")
        except Exception:
            print(f"[INSTALL] {pkg_name}")
            missing.append(pkg_name)

    if missing:
        print(f"Installing {len(missing)} missing packages ...")
        cmd = [sys.executable, "-m", "pip", "install"] + missing
        subprocess.check_call(cmd)
        print("Installation complete.")
    else:
        print("All required Python packages are already available.")
    PY
    echo "=== PGIRL environment setup complete ==="
    """
}
