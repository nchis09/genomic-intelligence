#!/usr/bin/env python3
"""Generate a self-contained Mermaid HTML fallback when nf-metro fails."""
import argparse
import html
import pathlib
import sys


def main():
    p = argparse.ArgumentParser(
        description="Wrap a Mermaid MMD file in a static HTML viewer."
    )
    p.add_argument("--input", required=True, help="Input .mmd file")
    p.add_argument("--output", required=True, help="Output .html file")
    p.add_argument("--title", default="Workflow diagram", help="Page title")
    args = p.parse_args()

    mmd_text = pathlib.Path(args.input).read_text(encoding="utf-8")
    safe_mmd = html.escape(mmd_text)

    html_text = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>{html.escape(args.title)}</title>
</head>
<body>
  <h1>{html.escape(args.title)}</h1>
  <p>Mermaid workflow diagram (fallback because nf-metro render failed).</p>
  <pre class="mermaid">
{safe_mmd}
  </pre>
  <script type="module">
    import mermaid from 'https://cdn.jsdelivr.net/npm/mermaid@10/dist/mermaid.esm.min.mjs';
    mermaid.initialize({{ startOnLoad: true, securityLevel: 'loose' }});
  </script>
</body>
</html>
"""

    pathlib.Path(args.output).write_text(html_text, encoding="utf-8")
    print(f"Wrote Mermaid fallback HTML to {args.output}")


if __name__ == "__main__":
    main()
