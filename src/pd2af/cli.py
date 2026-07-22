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
import pd2af.rules


_MODE_CHOICES = (
    "normal",
    "normal-no-complex",
    "keep-species",
    "keep-species-no-complex",
    "casq",
)

_LAYOUT_CHOICES = ("plain", "overlay", "dot", "auto")

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
    "dot": "graphviz `dot` auto-layout (requires `dot` on PATH)",
    "auto": "pick automatically from the input (graphviz `dot` for a map)",
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


def _read_input_map(input_file):
    """Read the input map from a file path, or from stdin if `input_file` is None.

    `momapy.io.core.read` needs a file path (for content-based format
    auto-detection and for the reader), so stdin bytes are buffered to a
    temporary file before reading, mirroring `momapy visualize`.
    """
    if input_file is not None:
        return momapy.io.core.read(input_file)
    if sys.stdin.isatty():
        print("error: no input file and stdin is not a pipe", file=sys.stderr)
        sys.exit(1)
    data = sys.stdin.buffer.read()
    if not data:
        print("error: no input received on stdin", file=sys.stderr)
        sys.exit(1)
    with tempfile.NamedTemporaryFile(delete=False) as tmp:
        tmp_path = tmp.name
        tmp.write(data)
    try:
        return momapy.io.core.read(tmp_path)
    finally:
        os.unlink(tmp_path)


def _run(args):
    reader_result = _read_input_map(args.input_file)
    input_map = reader_result.obj
    transform_result = pd2af.transform(
        input_map,
        mode=args.transformation_mode,
        layout_mode=args.layout_mode,
        influence_pairing=args.influence_pairing,
        set_active=args.set_active,
        set_inactive=args.set_inactive,
        set_all_active=args.set_all_active,
        set_all_inactive=args.set_all_inactive,
        exclude_groups=tuple(args.exclude_groups or ()),
        exclude_rules=tuple(args.exclude_rules or ()),
    )
    new_map = transform_result.obj
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
        for layout_mode in _LAYOUT_CHOICES
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


def _build_groups_data():
    """Per-mode excludable / mandatory rule groups, from the live registry.

    Excludable groups are the dependency-graph leaves `--exclude-group` can
    drop cleanly; mandatory groups are depended-on by another included group.
    """
    return {
        mode: {
            "excludable": sorted(excludable),
            "mandatory": sorted(mandatory),
        }
        for mode in _MODE_CHOICES
        for excludable, mandatory in (
            pd2af.rules.get_excludable_groups(mode.replace("-", "_")),
        )
    }


def _format_groups_tables(data):
    sections = []
    for mode, groups in data.items():
        rows = [[group, "excludable"] for group in groups["excludable"]]
        rows += [[group, "mandatory"] for group in groups["mandatory"]]
        sections.append(
            f"Mode `{mode}`:\n"
            + _render_table(["group", "status"], rows)
        )
    return "\n\n".join(sections)


def _list_groups(args):
    data = _build_groups_data()
    if args.json:
        print(json.dumps(data, indent=2))
    else:
        print(_format_groups_tables(data))


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
    parser.add_argument(
        "input_file",
        nargs="?",
        default=None,
        help="input CellDesigner XML file (reads from stdin if omitted)",
    )
    parser.add_argument(
        "-m",
        "--transformation-mode",
        choices=_MODE_CHOICES,
        default="normal",
        help=(
            "transformation mode (default: normal). 'normal' and "
            "'normal-no-complex' merge proteoforms of the same template and "
            "compartment into a single activity (true PD->AF transform) "
            "and require `--layout-mode dot` (or `auto`). 'keep-species' and "
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
            "layout strategy: auto (pick automatically from the input, "
            "default; a map gets `dot`), dot (graphviz auto-layout, requires "
            "`dot`), plain (reuse original positions, model elements only), or "
            "overlay (reuse full original layout with non-model elements greyed "
            "out). 'normal' and 'normal-no-complex' transformation modes require "
            "`dot` (or `auto`)."
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
            "positions are real; in `dot` it is ignored."
        ),
    )
    parser.add_argument(
        "-a",
        "--set-active",
        action="append",
        default=None,
        metavar="ID",
        dest="set_active",
        help=(
            "mark the element with this id_ (a species or entity pool) as "
            "active, emitting hasActivity(..., isInputParameter). Repeatable: "
            "-a sa1 -a sa2. Wins over --set-all-inactive for these ids."
        ),
    )
    parser.add_argument(
        "-i",
        "--set-inactive",
        action="append",
        default=None,
        metavar="ID",
        dest="set_inactive",
        help=(
            "mark the element with this id_ (a species or entity pool) as "
            "NOT active, suppressing any hasActivity for it (overrides the "
            "automatic activity discovery: active flag/state, modulation "
            "source, reaction modifier, gate input, phenotype). Repeatable: "
            "-i sa1 -i sa2. Wins over --set-all-active for these ids. "
            "Erroring if an id is also passed to --set-active."
        ),
    )
    global_activity_group = parser.add_mutually_exclusive_group()
    global_activity_group.add_argument(
        "-A",
        "--set-all-active",
        action="store_true",
        dest="set_all_active",
        help=(
            "mark every top-level species / entity pool as active "
            "(subunits excluded). Per-id --set-inactive overrides this for "
            "the named ids. Mutually exclusive with --set-all-inactive."
        ),
    )
    global_activity_group.add_argument(
        "-I",
        "--set-all-inactive",
        action="store_true",
        dest="set_all_inactive",
        help=(
            "suppress activity for every element (including subunits). "
            "Per-id --set-active overrides this for the named ids. "
            "Mutually exclusive with --set-all-active."
        ),
    )
    parser.add_argument(
        "--exclude-group",
        action="append",
        default=None,
        metavar="GROUP",
        dest="exclude_groups",
        help=(
            "drop a whole rule group (e.g. `activity:phenotype` to stop "
            "treating phenotypes as activities, or `paths:chaining` to keep "
            "only single-hop influences). Repeatable. Excluding a group that "
            "another included group depends on is an error. Run "
            "`pd2af list-groups` to see the excludable groups per mode."
        ),
    )
    parser.add_argument(
        "--disable-rule",
        action="append",
        default=None,
        metavar="RULE",
        dest="exclude_rules",
        help=(
            "drop a single rule by identifier (the fine scalpel for the "
            "whole-with-scalpel table groups, e.g. "
            "`modulation_kind:celldesigner:catalysis`). Repeatable. Prefer "
            "`--exclude-group` for coherent behaviors."
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


def _add_list_groups_parser(subparsers):
    parser = subparsers.add_parser(
        "list-groups",
        help="list, per mode, the excludable and mandatory rule groups",
        description=(
            "List the rule groups each transformation mode uses, marked "
            "excludable (a dependency-graph leaf `--exclude-group` can drop) "
            "or mandatory (depended-on by another group, so not excludable)."
        ),
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="emit the listing as JSON instead of human-readable tables",
    )
    parser.set_defaults(func=_list_groups)


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    else:
        argv = list(argv)
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
    _add_list_groups_parser(subparsers)
    # Default to the transform subcommand when the first token isn't a known
    # subcommand or a help flag, so `pd2af map.xml` works like
    # `pd2af transform map.xml`.
    if not argv or (
        argv[0] not in subparsers.choices and argv[0] not in ("-h", "--help")
    ):
        argv = ["transform", *argv]
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
