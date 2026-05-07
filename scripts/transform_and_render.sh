#!/usr/bin/env bash
# Transform every CellDesigner map in INPUT_DIR with pd2af and render
# each result to a diagram in OUTPUT_DIR using momapy visualize.

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "usage: $0 INPUT_DIR OUTPUT_DIR" >&2
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"

MODE="normal"
LAYOUT_MODE="auto"
RENDER_FORMAT="pdf"   # any momapy-supported format (pdf, svg, png, ...)

mkdir -p "$OUTPUT_DIR"

shopt -s nullglob
input_files=("$INPUT_DIR"/*.xml "$INPUT_DIR"/*.sbml)
shopt -u nullglob

if [ "${#input_files[@]}" -eq 0 ]; then
    echo "no .xml/.sbml files found in $INPUT_DIR" >&2
    exit 1
fi

for input_file in "${input_files[@]}"; do
    base_name="$(basename "${input_file%.*}")"
    output_file="$OUTPUT_DIR/${base_name}.${MODE}.${LAYOUT_MODE}.${RENDER_FORMAT}"
    echo "==> $input_file -> $output_file"
    pd2af -m "$MODE" -l "$LAYOUT_MODE" "$input_file" \
        | momapy visualize - -o "$output_file"
done
