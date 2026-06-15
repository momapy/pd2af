import argparse
import json
import os
import sys
import tempfile
import textwrap

import momapy.io.core

import pd2af
import pd2af.core
import pd2af.languages


_MODE_CHOICES = (
    "normal",
    "normal-no-complex",
    "keep-species",
    "keep-species-no-complex",
    "casq",
)

_LAYOUT_CHOICES = ("plain", "overlay", "auto")

_INFLUENCE_PAIRING_CHOICES = ("cross", "nearest")

# Human-facing descriptions of each transformation mode, keyed by the mode
# name. The PTM/complex behavior is presentation prose; the authoritative
# compatibility data is derived from pd2af.core / pd2af.languages, never
# duplicated here. `celldesigner_only` records a documented constraint that the
# core does not hard-enforce.
_TRANSFORMATION_MODE_INFO = {
    "normal": {
        "summary": (
            "merge forms of the same base species or entity pool into a "
            "single activity"
        ),
        "celldesigner_only": False,
    },
    "normal-no-complex": {
        "summary": (
            "merge forms of the same base species or entity pool into a "
            "single activity; additionally, replace complexes with their "
            "active subunits if any, promoting them to top-level activities"
        ),
        "celldesigner_only": False,
    },
    "keep-species": {
        "summary": "create one activity per distinct active species or entity pool",
        "celldesigner_only": False,
    },
    "keep-species-no-complex": {
        "summary": (
            "create one activity per distinct active species or entity pool; "
            "additionally, replace complexes with their active subunits if "
            "any, promoting them to top-level activities"
        ),
        "celldesigner_only": False,
    },
    "casq": {
        "summary": (
            "reproduce the CasQ transformation (Aghamiri et al., Bioinfo., 2020)"
        ),
        "celldesigner_only": True,
    },
}

_LAYOUT_MODE_INFO = {
    "plain": "reuse original positions",
    "overlay": "reuse full original layout with unmapped layout elements dimmed",
    "auto": "graphviz `dot` auto-layout (requires `dot` on PATH)",
}

# Human-facing display names for the input-language tokens (also the single
# source for the language order shown in the listing).
_INPUT_LANGUAGE_DISPLAY_NAMES = {
    pd2af.languages.CELLDESIGNER: "CellDesigner",
    pd2af.languages.SBGN_PD: "SBGN PD",
}


def _unsupported_layout_modes_for_language(language):
    """Concrete layout modes the given input language rejects.

    Derived from pd2af.core: SBGN-AF output (SBGN-PD input) does not support
    the `overlay` dimming; CellDesigner output supports every layout mode.
    """
    return tuple(
        layout_mode
        for layout_mode in pd2af.core.LAYOUT_MODES
        if language == pd2af.languages.SBGN_PD
        and layout_mode not in pd2af.core.SBGN_AF_LAYOUT_MODES
    )


def _compatible_language_names_for_mode(mode):
    """Display names of the input languages a transformation mode applies to.

    Every mode works on CellDesigner input; the `celldesigner_only` modes
    (currently `casq`) are restricted to it, the rest also accept SBGN-PD.
    """
    if _TRANSFORMATION_MODE_INFO[mode]["celldesigner_only"]:
        return [_INPUT_LANGUAGE_DISPLAY_NAMES[pd2af.languages.CELLDESIGNER]]
    return list(_INPUT_LANGUAGE_DISPLAY_NAMES.values())


def _compatible_language_names_for_layout_mode(layout_mode):
    """Display names of the input languages whose output supports a layout mode.

    Derived from pd2af.core: SBGN-AF output (SBGN-PD input) does not support
    the `overlay` dimming; CellDesigner output supports every layout mode.
    """
    return [
        display_name
        for language, display_name in _INPUT_LANGUAGE_DISPLAY_NAMES.items()
        if layout_mode not in _unsupported_layout_modes_for_language(language)
    ]

_EXTENSION_TO_WRITER = {
    ".xml": "celldesigner",
    ".sbml": "celldesigner",
    ".sbgn": "sbgnml",
    ".sbgnml": "sbgnml",
    ".pickle": "pickle",
    ".pkl": "pickle",
}


def _writer_for_output(output_path):
    ext = os.path.splitext(output_path)[1].lower()
    return _EXTENSION_TO_WRITER.get(ext, "pickle")


def _write_map_to_stdout(cd_map):
    with tempfile.NamedTemporaryFile(suffix=".pickle", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        momapy.io.core.write(cd_map, tmp_path, writer="pickle")
        with open(tmp_path, "rb") as f:
            sys.stdout.buffer.write(f.read())
    finally:
        os.unlink(tmp_path)


def _run(args):
    reader_result = momapy.io.core.read(args.input_file)
    input_map = reader_result.obj
    new_map = pd2af.transform(
        input_map,
        mode=args.transformation_mode,
        layout_mode=args.layout_mode,
        influence_pairing=args.influence_pairing,
    )
    if args.output is None:
        _write_map_to_stdout(new_map)
    else:
        writer = _writer_for_output(args.output)
        momapy.io.core.write(new_map, args.output, writer=writer)


def _build_modes_data():
    """Assemble the structured `list-modes` payload from the source-of-truth.

    The payload mirrors the rendered tables exactly: one entry per table row,
    one key per column. Compatibility is derived from pd2af.core, not
    duplicated, so the listing cannot drift from the validation the transform
    actually enforces.
    """
    transformation_modes = [
        {
            "transformation_mode": mode,
            "layout_modes": list(
                pd2af.core.get_compatible_layout_modes_for_transformation_mode(mode)
            ),
            "languages": _compatible_language_names_for_mode(mode),
            "description": info["summary"],
        }
        for mode, info in _TRANSFORMATION_MODE_INFO.items()
    ]
    layout_modes = [
        {
            "layout_mode": layout_mode,
            "languages": _compatible_language_names_for_layout_mode(layout_mode),
            "description": _LAYOUT_MODE_INFO[layout_mode],
        }
        for layout_mode in pd2af.core.LAYOUT_MODES
    ]
    return {
        "transformation_modes": transformation_modes,
        "layout_modes": layout_modes,
    }


def _render_table(headers, rows, max_widths=None):
    """Render an ASCII grid table, wrapping cells whose column has a max width.

    `max_widths` maps a column index to the width at which that column's cells
    are wrapped (others are never wrapped). Wrapped cells span several physical
    rows within the same logical row.
    """
    max_widths = max_widths or {}
    column_count = len(headers)
    wrapped_rows = []
    for row in rows:
        wrapped_cells = []
        for column_index, cell in enumerate(row):
            max_width = max_widths.get(column_index)
            if max_width is not None:
                wrapped_cells.append(textwrap.wrap(cell, max_width) or [""])
            else:
                wrapped_cells.append([cell])
        wrapped_rows.append(wrapped_cells)

    column_widths = []
    for column_index in range(column_count):
        width = len(headers[column_index])
        for wrapped_cells in wrapped_rows:
            for line in wrapped_cells[column_index]:
                width = max(width, len(line))
        column_widths.append(width)

    separator = "+" + "+".join("-" * (width + 2) for width in column_widths) + "+"

    def render_physical_row(cells):
        return (
            "| "
            + " | ".join(
                cell.ljust(column_widths[column_index])
                for column_index, cell in enumerate(cells)
            )
            + " |"
        )

    lines = [separator, render_physical_row(headers), separator]
    for wrapped_cells in wrapped_rows:
        height = max(len(cell_lines) for cell_lines in wrapped_cells)
        for line_index in range(height):
            lines.append(
                render_physical_row(
                    [
                        wrapped_cells[column_index][line_index]
                        if line_index < len(wrapped_cells[column_index])
                        else ""
                        for column_index in range(column_count)
                    ]
                )
            )
        lines.append(separator)
    return "\n".join(lines)


def _render_modes_table(title, columns, rows, max_widths=None):
    """Render one titled `list-modes` table from payload rows.

    `columns` is a list of (header, key) pairs naming each column and the row
    key it reads; list-valued cells are joined with ", " so the table is an
    exact view of the JSON payload.
    """
    headers = [header for header, _ in columns]
    table_rows = [
        [
            ", ".join(row[key]) if isinstance(row[key], list) else row[key]
            for _, key in columns
        ]
        for row in rows
    ]
    return title + "\n" + _render_table(headers, table_rows, max_widths=max_widths)


def _format_modes_tables(data):
    sections = [
        _render_modes_table(
            "Transformation modes (--transformation-mode, -m):",
            [
                ("transformation mode", "transformation_mode"),
                ("layout modes", "layout_modes"),
                ("languages", "languages"),
                ("description", "description"),
            ],
            data["transformation_modes"],
            max_widths={3: 48},
        ),
        _render_modes_table(
            "Layout modes (--layout-mode, -l):",
            [
                ("layout mode", "layout_mode"),
                ("languages", "languages"),
                ("description", "description"),
            ],
            data["layout_modes"],
        ),
    ]
    return "\n\n".join(sections)


def _list_modes(args):
    data = _build_modes_data()
    if args.json:
        print(json.dumps(data, indent=2))
    else:
        print(_format_modes_tables(data))


def _add_transform_parser(subparsers):
    parser = subparsers.add_parser(
        "transform",
        help="transform a process-description map into an activity-flow map",
        description=(
            "Transform a CellDesigner process-description map into an "
            "activity-flow map. Writes the map to stdout as a momapy pickle "
            "(preserves layout styling) so it can be piped into "
            "`momapy visualize`. With -o, the writer is chosen from the "
            "output file extension (.xml/.sbml -> CellDesigner XML, "
            ".pickle/.pkl -> pickle; defaults to pickle)."
        ),
    )
    parser.add_argument("input_file", help="input CellDesigner XML file")
    parser.add_argument(
        "-m",
        "--transformation-mode",
        choices=_MODE_CHOICES,
        default="normal",
        help=(
            "transformation mode (default: normal). 'normal' and "
            "'normal-no-complex' merge proteoforms of the same template and "
            "compartment into a single activity (true PD->AF transform) "
            "and require `--layout-mode auto`. 'keep-species' and "
            "'keep-species-no-complex' keep each PD species as its own "
            "activity. The '*-no-complex' variants drop complexes that "
            "have an active subunit, routing influences through the "
            "subunits. 'casq' emits one activity per surviving PD "
            "species after applying CASQ-style deletion rules "
            "(heterodimer simplification, name-preserving step pruning, "
            "transport collapse) with single-hop rewiring across "
            "deleted intermediates. Run `pd2af list-modes` for the full "
            "compatibility matrix."
        ),
    )
    parser.add_argument(
        "-l",
        "--layout-mode",
        choices=_LAYOUT_CHOICES,
        default="auto",
        help=(
            "layout strategy: auto (graphviz auto-layout, requires `dot`, "
            "default), plain (reuse original positions, model elements "
            "only), or overlay (reuse full original layout with non-model "
            "elements greyed out). 'normal' and 'normal-no-complex' transformation "
            "modes require `auto`."
        ),
    )
    parser.add_argument(
        "-p",
        "--influence-pairing",
        choices=_INFLUENCE_PAIRING_CHOICES,
        default="cross",
        help=(
            "how to draw an influence whose source or target maps to several "
            "layout glyphs: cross (default, one arc per source/target pair) or "
            "nearest (a single arc between the closest pair). 'nearest' only "
            "takes effect with `--layout-mode plain` or `overlay`, where glyph "
            "positions are real; in `auto` it is ignored."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="write output to this file instead of stdout",
    )
    parser.set_defaults(func=_run)


def _add_list_modes_parser(subparsers):
    parser = subparsers.add_parser(
        "list-modes",
        help="list transformation modes, layout modes, and compatibilities",
        description=(
            "List the available transformation modes, layout modes, their "
            "compatibilities, and the supported input languages."
        ),
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="emit the listing as JSON instead of human-readable tables",
    )
    parser.set_defaults(func=_list_modes)


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="pd2af",
        description=(
            "Transform a CellDesigner process-description map into an "
            "activity-flow map."
        ),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    _add_transform_parser(subparsers)
    _add_list_modes_parser(subparsers)
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
