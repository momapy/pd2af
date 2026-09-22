"""The ``pd2af`` command-line interface."""

import argparse
import importlib.metadata
import json
import os
import sys
import tempfile
import textwrap
import typing

import momapy.cli
import momapy.io.core

import pd2af
import pd2af.core
import pd2af.modes
import pd2af.asp.rules


_LAYOUT_CHOICES = tuple(pd2af.modes.LayoutMode) + (pd2af.modes.AUTO,)

_LAYOUT_MODE_INFO = pd2af.modes.LAYOUT_MODE_DESCRIPTIONS | {
    pd2af.modes.AUTO: pd2af.modes.AUTO_DESCRIPTION
}


def _compatible_language_names_for_mode(
    mode: pd2af.modes.TransformationMode,
) -> list[str]:
    """Display names of the input languages a transformation mode applies to."""
    return [
        properties["display_name"]
        for language, properties in pd2af.modes.LANGUAGES.items()
        if language in mode.compatible_languages
    ]


def _compatible_language_names_for_layout_mode(
    layout_mode: pd2af.modes.LayoutMode | str,
) -> list[str]:
    """Display names of the input languages whose output supports a layout mode.

    Derived from pd2af.modes: SBGN-AF output (SBGN-PD input) does not
    support the `overlay` dimming; CellDesigner output supports every layout
    mode. The `auto` meta value is not a concrete layout mode, so no language
    rejects it and every display name is returned.
    """
    return [
        properties["display_name"]
        for language, properties in pd2af.modes.LANGUAGES.items()
        if layout_mode == pd2af.modes.AUTO
        or layout_mode in pd2af.modes.LAYOUT_MODES_BY_LANGUAGE[language]
    ]


_EXTENSION_TO_WRITER = {
    ".xml": "celldesigner",
    ".sbml": "celldesigner",
    ".sbgn": "sbgnml",
    ".sbgnml": "sbgnml",
    ".pickle": "pickle",
    ".pkl": "pickle",
}


def _writer_for_output(output_path: str) -> str:
    ext = os.path.splitext(output_path)[1].lower()
    return _EXTENSION_TO_WRITER.get(ext, "pickle")


def _write_map_to_stdout(cd_map: typing.Any) -> None:
    with tempfile.NamedTemporaryFile(suffix=".pickle", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        momapy.io.core.write(cd_map, tmp_path, writer="pickle")
        with open(tmp_path, "rb") as f:
            sys.stdout.buffer.write(f.read())
    finally:
        os.unlink(tmp_path)


def _read_input_map(input_file: str | None) -> typing.Any:
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


def _run(args: argparse.Namespace) -> None:
    reader_result = _read_input_map(args.input_file)
    input_map = reader_result.obj
    transform_result = pd2af.transform(
        input_map,
        mode=args.transformation_mode,
        layout_mode=args.layout_mode,
        influence_pairing=args.influence_pairing,
        keep_species=args.keep_species,
        drop_compartments=args.drop_compartments,
        set_active=args.set_active,
        set_inactive=args.set_inactive,
        set_all_active=args.set_all_active,
        set_all_inactive=args.set_all_inactive,
        exclude_groups=tuple(args.exclude_groups or ()),
        exclude_rules=tuple(args.exclude_rules or ()),
        element_to_annotations=reader_result.element_to_annotations,
        element_to_notes=reader_result.element_to_notes,
    )
    new_map = transform_result.obj
    if args.output is not None:
        writer = _writer_for_output(args.output)
        momapy.io.core.write(
            new_map,
            args.output,
            writer=writer,
            element_to_annotations=transform_result.element_to_annotations,
            element_to_notes=transform_result.element_to_notes,
        )
    elif not args.visualize:
        # The stdout pickle serializes only the map object; annotations and
        # notes live in side-tables that a bare-map pickle cannot carry, so
        # they are dropped on this path (use -o file.xml / .sbgn to keep them).
        _write_map_to_stdout(new_map)
    if args.visualize:
        momapy.cli._visualize_map(new_map)


def _describe_option_default(name: str, option: dict) -> str:
    """An option's default, naming the modes whose `default_options` differ."""
    default = option["default"]
    exceptions = [
        mode.name
        for mode in pd2af.modes.get_transformation_modes().values()
        if mode.default_options.get(name, default) != default
    ]
    description = "on" if default else "off"
    if exceptions:
        description += f" ({'off' if default else 'on'} for {', '.join(exceptions)})"
    return description


def _build_modes_data() -> dict:
    """Assemble the structured `list-modes` payload from the source-of-truth.

    The payload mirrors the rendered tables exactly: one entry per table row,
    one key per column. Compatibility is read off each mode object, not
    duplicated, so the listing cannot drift from the validation the transform
    actually enforces.
    """
    transformation_modes = [
        {
            "transformation_mode": mode.name,
            "languages": _compatible_language_names_for_mode(mode),
            "description": mode.docs,
        }
        for mode in pd2af.modes.get_transformation_modes().values()
    ]
    transformation_options = [
        {
            "option": option["flag"],
            "default": _describe_option_default(name, option),
            "description": option["description"],
        }
        for name, option in pd2af.modes.TRANSFORMATION_OPTIONS.items()
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
        "transformation_options": transformation_options,
        "layout_modes": layout_modes,
    }


def _render_table(
    headers: list[str],
    rows: list[list[str]],
    max_widths: dict[int, int] | None = None,
) -> str:
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

    def render_physical_row(cells: list[str]) -> str:
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


def _render_modes_table(
    title: str,
    columns: list[tuple[str, str]],
    rows: list[dict],
    max_widths: dict[int, int] | None = None,
) -> str:
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


def _format_modes_tables(data: dict) -> str:
    sections = [
        _render_modes_table(
            "Transformation modes (--transformation-mode, -m):",
            [
                ("transformation mode", "transformation_mode"),
                ("languages", "languages"),
                ("description", "description"),
            ],
            data["transformation_modes"],
            max_widths={2: 48},
        )
        + "\nThe transformation options below change the layout modes on offer.",
        _render_modes_table(
            "Transformation options:",
            [
                ("option", "option"),
                ("default", "default"),
                ("description", "description"),
            ],
            data["transformation_options"],
            max_widths={2: 48},
        ),
        _render_modes_table(
            "Layout modes (--layout-mode, -l):",
            [
                ("layout mode", "layout_mode"),
                ("languages", "languages"),
                ("description", "description"),
            ],
            data["layout_modes"],
        )
        + "\n`plain` and `overlay` are available exactly when --keep-species is "
        "set without --drop-compartments.",
    ]
    return "\n\n".join(sections)


def _list_modes(args: argparse.Namespace) -> None:
    data = _build_modes_data()
    if args.json:
        print(json.dumps(data, indent=2))
    else:
        print(_format_modes_tables(data))


def _build_groups_data() -> dict:
    """Per-mode excludable / mandatory rule groups, from the live registry.

    Excludable groups are the dependency-graph leaves `--exclude-group` can
    drop cleanly; mandatory groups are depended-on by another included group.
    """
    data = {}
    for mode_name in pd2af.modes.get_transformation_modes():
        excludable, mandatory = pd2af.asp.rules.get_excludable_groups(mode_name)
        data[mode_name] = {
            "excludable": sorted(excludable),
            "mandatory": sorted(mandatory),
        }
    return data


def _format_groups_tables(data: dict) -> str:
    sections = []
    for mode, groups in data.items():
        rows = [[group, "excludable"] for group in groups["excludable"]]
        rows += [[group, "mandatory"] for group in groups["mandatory"]]
        sections.append(f"Mode `{mode}`:\n" + _render_table(["group", "status"], rows))
    return "\n\n".join(sections)


def _list_groups(args: argparse.Namespace) -> None:
    data = _build_groups_data()
    if args.json:
        print(json.dumps(data, indent=2))
    else:
        print(_format_groups_tables(data))


def _add_transform_parser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "transform",
        help="transform a process-description map into an activity-flow map",
        description=(
            "Transform a process-description map (CellDesigner or SBGN-PD) "
            "into an activity-flow map. Writes the map to stdout as a momapy "
            "pickle (preserves layout styling) so it can be piped into "
            "`momapy visualize`. With -o, the writer is chosen from the "
            "output file extension (.xml/.sbml -> CellDesigner XML, "
            ".sbgn/.sbgnml -> SBGN-ML, .pickle/.pkl -> pickle; defaults to "
            "pickle)."
        ),
    )
    parser.add_argument(
        "input_file",
        nargs="?",
        default=None,
        help=(
            "input process-description map, CellDesigner XML or SBGN-ML "
            "(reads from stdin if omitted)"
        ),
    )
    parser.add_argument(
        "-m",
        "--transformation-mode",
        choices=tuple(pd2af.modes.get_transformation_modes()),
        default="normal",
        help=(
            "transformation mode (default: normal). A mode decides what "
            "counts as an activity. 'normal' keeps a complex as an activity "
            "of its own, routing a subunit's influences to the complex it "
            "belongs to. 'no-complex' replaces a complex that has an active "
            "subunit with those subunits, promoting them to top-level "
            "activities. 'keep-reactions' keeps the PD topology itself: "
            "every species is an activity and every reaction becomes a "
            "positive influence from each of its reactants to each of "
            "its products, alongside the modulation arcs and reaction "
            "modifiers, with no multi-hop or consumption inference. "
            "Run `pd2af list-modes` for the full compatibility matrix."
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
            "out). `dot` is required whenever activities merge: unless "
            "`--keep-species` is set, and always under `--drop-compartments`."
        ),
    )
    parser.add_argument(
        "-p",
        "--influence-pairing",
        choices=tuple(pd2af.modes.InfluencePairingMode),
        default="cross",
        help=(
            "how to draw an influence whose source or target maps to several "
            "layout glyphs: cross (default, one arc per source/target pair) or "
            "nearest (a single arc between the closest pair). 'nearest' only "
            "takes effect with `--layout-mode plain` or `overlay`, where glyph "
            "positions are real; in `dot` it is ignored."
        ),
    )
    for name, option in pd2af.modes.TRANSFORMATION_OPTIONS.items():
        parser.add_argument(
            option["flag"],
            dest=name,
            action=argparse.BooleanOptionalAction,
            default=None,
            help=f"{option['description']}; default: "
            f"{_describe_option_default(name, option)}",
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
        "--exclude-rule",
        action="append",
        default=None,
        metavar="RULE",
        dest="exclude_rules",
        help=(
            "drop a single rule by identifier (the fine scalpel for the "
            "whole-with-scalpel table groups, e.g. "
            "`influences:kind:celldesigner:catalysis`). Repeatable. Prefer "
            "`--exclude-group` for coherent behaviors."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help=(
            "write output to this file instead of stdout. The writer comes "
            "from the extension (.xml/.sbml -> CellDesigner XML, "
            ".sbgn/.sbgnml -> SBGN-ML, .pickle/.pkl -> pickle); any other "
            "extension, including .svg, writes a pickle and never renders an "
            "image. Input RDF "
            "annotations and notes are carried onto the corresponding output "
            "elements only for file output (.xml/.sbml/.sbgn/.sbgnml); the "
            "stdout pickle cannot carry them and drops them."
        ),
    )
    parser.add_argument(
        "-V",
        "--visualize",
        action="store_true",
        help=(
            "open the output map in the momapy browser viewer. With -o the "
            "map is also written to the file; without -o nothing is written "
            "to stdout."
        ),
    )
    parser.set_defaults(func=_run)


def _add_list_modes_parser(subparsers: argparse._SubParsersAction) -> None:
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


def _add_list_groups_parser(subparsers: argparse._SubParsersAction) -> None:
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


def main(argv: list[str] | None = None) -> None:
    """Parse ``argv`` (``sys.argv`` by default) and run the named subcommand.

    Expected input, compatibility and file-access failures are reported as a
    short message on stderr with a nonzero exit status; any other exception
    keeps its traceback, and library callers of :func:`pd2af.transform` still
    see the exception itself.
    """
    if argv is None:
        argv = sys.argv[1:]
    else:
        argv = list(argv)
    parser = argparse.ArgumentParser(
        prog="pd2af",
        description=(
            "Transform a process-description map (CellDesigner or SBGN-PD) "
            "into an activity-flow map."
        ),
    )
    parser.add_argument(
        "--version",
        action="version",
        version=importlib.metadata.version("pd2af"),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    _add_transform_parser(subparsers)
    _add_list_modes_parser(subparsers)
    _add_list_groups_parser(subparsers)
    # Default to the transform subcommand when the first token isn't a known
    # subcommand or a help flag, so `pd2af map.xml` works like
    # `pd2af transform map.xml`.
    if not argv or (
        argv[0] not in subparsers.choices
        and argv[0] not in ("-h", "--help", "--version")
    ):
        argv = ["transform", *argv]
    args = parser.parse_args(argv)
    try:
        args.func(args)
    except (ValueError, OSError) as error:
        print(f"error: {error}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
