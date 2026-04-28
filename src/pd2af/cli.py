import argparse
import os
import sys
import tempfile

import momapy.io.core

import pd2af


_MODE_CHOICES = ("normal", "no-complex", "pure-af")

_LAYOUT_CHOICES = ("plain", "overlay", "auto")

_EXTENSION_TO_WRITER = {
    ".xml": "celldesigner",
    ".sbml": "celldesigner",
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
    cd_map = reader_result.obj
    new_map = pd2af.transform(
        cd_map,
        mode=args.mode,
        layout_mode=args.layout,
    )
    if args.output is None:
        _write_map_to_stdout(new_map)
    else:
        writer = _writer_for_output(args.output)
        momapy.io.core.write(new_map, args.output, writer=writer)


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="pd2af",
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
        "--mode",
        choices=_MODE_CHOICES,
        default="normal",
        help=(
            "transformation mode (default: normal). 'pure-af' merges "
            "proteoforms of the same template and compartment into a single "
            "activity and requires `--layout auto`."
        ),
    )
    parser.add_argument(
        "-l",
        "--layout",
        choices=_LAYOUT_CHOICES,
        default="plain",
        help=(
            "layout strategy: plain (reuse original positions, model "
            "elements only, default), overlay (reuse full original layout "
            "with non-model elements greyed out), or auto (graphviz "
            "auto-layout, requires `dot`)"
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="write output to this file instead of stdout",
    )
    args = parser.parse_args(argv)
    _run(args)


if __name__ == "__main__":
    main()
