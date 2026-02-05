#!/usr/bin/env bash
# Render manuscript.md and supplementary_materials.md to PDF via pandoc + xelatex
# Usage: bash render_manuscript.sh [--supplement-only]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Use RStudio's bundled pandoc if system pandoc is not available
if ! command -v pandoc &>/dev/null; then
    export PATH="/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/x86_64:$PATH"
fi

if ! command -v pandoc &>/dev/null; then
    echo "Error: pandoc not found. Install pandoc or RStudio." >&2
    exit 1
fi

PANDOC_OPTS=(
    --pdf-engine=xelatex
    --resource-path="$SCRIPT_DIR"
    -V geometry:margin=1in
    -V fontsize=11pt
    -V mainfont="Times New Roman"
    -V linestretch=1.5
    -V numbersections
    --standalone
)

if [[ "${1:-}" != "--supplement-only" ]]; then
    pandoc "$SCRIPT_DIR/manuscript.md" -o "$SCRIPT_DIR/manuscript.pdf" "${PANDOC_OPTS[@]}"
    echo "Rendered: manuscript.pdf"
fi

pandoc "$SCRIPT_DIR/supplementary_materials.md" -o "$SCRIPT_DIR/supplementary_materials.pdf" "${PANDOC_OPTS[@]}"
echo "Rendered: supplementary_materials.pdf"
